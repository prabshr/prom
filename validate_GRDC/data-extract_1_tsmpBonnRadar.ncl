load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$TPS/functions/rootDist_clm_nclf.ncl"
external EXT1 "$TPS/external/pfread_dim.so"
external EXT2 "$TPS/external/pfread_var.so"
external EXT3 "$TPS/external/gwt_depth.so"
external EXT4 "$TPS/external/qoverland.so"
begin

 manningscoefConst = 5.52e-6 *3600. ;[s m-3]
 readclm           = True          ;
 readpfl           = True
 rootzone          = False
 fplot_root        = False
 diri_db           = "/p/project/chbn33/hbn331/database/bonnRadar/"
 diri_op           = "/p/largedata/hbn33/shrestha1/EA_PROM" + expt +"/HET2483/runE_" + inst + "_HET2483/"
 diri_db           = "/daten/daten04/pshrestha/shared-data/bonnRadar/" 
 diri_op           = "/daten/daten04/pshrestha/shared-data/bonnRadar_20160513/HET2483/runE_" + inst + "_HET2483/"
 ;expt
 ;year              = 2015;
 ;inst

 if (expt.eq.1) then
    year = 2015
 end if
 if (expt.eq.2) then
    year = 2016
 end if
 if (expt.eq.3) then
    year = 2017
 end if
;READ
;-------------------------------------------------------------------------
 ;
 fil_geo  = diri_db + "/cosmo/2014111500/laf2014111500.grb"
 fil_pft  = diri_db + "/clm/surfdata_HET"+year+"_0294x0294.nc"
 fil_sol  = diri_db + "/parflow/bonnRadarSoilInd.pfb"
 fil_xsl  = diri_db + "/parflow/xslope.pfb"
 fil_ysl  = diri_db + "/parflow/yslope.pfb"
 fil_por  = diri_db + "/parflow/rurlaf.out.porosity.pfb"
 fil_sps  = diri_db + "/parflow/rurlaf.out.specific_storage.pfb"

 ;
 diri_clm = diri_op + "/clmout/"
 diri_pfl = diri_op + "/pflout/"

;Latitutde-Longitude
;-----------------------------------
 geo_grb  = addfile(fil_geo,"r")
 lat2d    = geo_grb->g10_lat_0(3:296,3:296)
 lon2d    = geo_grb->g10_lon_1(3:296,3:296)

 if (readclm) then
 print("Reading CLM data ....")

;Land Surface States and Fluxes
;-----------------------------------
 fils_clm   = systemfunc("ls " + diri_clm + "/*.nc")
 ntim       = dimsizes(fils_clm)  ;number of history tapes
 ncdf_clm   = addfiles(fils_clm,"r")
 ListSetType (ncdf_clm, "cat") 
 time       = ncdf_clm[:]->time
 fshv       = ncdf_clm[:]->FSH_V
 fshg       = ncdf_clm[:]->FSH_G
 fctr       = ncdf_clm[:]->FCTR
 fgev       = ncdf_clm[:]->FGEV
 fcev       = ncdf_clm[:]->FCEV
 fira       = ncdf_clm[:]->FIRA    ;Net Longwave +ve upwards
 fsa        = ncdf_clm[:]->FSA     ;Net Shortwave
 fgr        = ncdf_clm[:]->FGR     ;Ground HEat Flux
 t2m        = ncdf_clm[:]->TSA
 h2osoi     = ncdf_clm[:]->H2OSOI
 dzsoi      = ncdf_clm[0]->DZSOI(:,0,0)
 
;Heterogenity of Plant Phenology, and weighted root zone Sw
;-----------------------------------
 dummy    = dzsoi*0.
 roota_par = (/-999.,7.,7.,7.,7.,7.,6.,6.,6.,7. ,7. ,7. ,11.,11.,11.,6.,6./)
 rootb_par = (/-999.,2.,2.,2.,1.,1.,2.,2.,2.,1.5,1.5,1.5,2. ,2. ,2. ,3.,3./)

 nz    = dimsizes(dzsoi)
 npft  = dimsizes(roota_par)
 rdist = new((/npft,nz/),"float")
 do ip = 1, dimsizes(roota_par)-1
   rdist(ip,:) = rootCLM(roota_par(ip), rootb_par(ip),dummy, dzsoi,False)
 end do

 if (fplot_root) then
   spft = (/1,7,13,15/)
   wks = gsn_open_wks("x11","debug_root")
   res = True
   res@tmXTOn   = False
   res@tmYROn   = False
   res@vpWidthF = 0.5
  
   res@xyMarkLineMode = "MarkLines"
   res@xyLineThicknessF = 3.0
   res@trYReverse = True
   plot = gsn_csm_xy(wks,rdist(spft,:),dzsoi,res)
   print("PFT    Sum Root Fractions")
   print(spft + "       " + dim_sum_n(rdist(spft,:),1))
 end if
;
;
 nc_pft     = addfile(fil_pft,"r")
 pft_pc     = nc_pft->PCT_PFT
 mlai       = tofloat(nc_pft->MONTHLY_LAI)
 ndim       = dimsizes(pft_pc)
 pft        = new((/ndim(1),ndim(2)/),"integer")
 lai        = new((/12,ndim(1),ndim(2)/),"float")
 swr        = new((/ntim,ndim(1),ndim(2)/),"float")
 sw5        = new((/ntim,ndim(1),ndim(2)/),"float")    ;Top 5 cm soil moisture

 if (rootzone) then
 do jj = 0, ndim(1)-1
 do ii = 0, ndim(2)-1
   pft(jj,ii)   = ind(pft_pc(:,jj,ii).eq.100)
   lai(:,jj,ii) = mlai(:,pft(jj,ii),jj,ii)
   wgt = h2osoi(:,0,jj,ii)*0.0
   dpdzsum = 0.
   stime = 38
   do k = 0, nz-1 
     ;print(rdist(pft(jj,ii),k)*h2osoi(stime,k,jj,ii))
     dp_dz = rdist(pft(jj,ii),k)/dzsoi(k)
     wgt   = wgt + h2osoi(:,k,jj,ii)*dp_dz*dzsoi(k)
     dpdzsum = dpdzsum + dp_dz*dzsoi(k)
   end do
   swr(:,jj,ii) = wgt/dpdzsum
   ;print("Root Zone Weighted Soil Moisture: " + swr(stime,jj,ii))
   sw5(:,jj,ii) = (h2osoi(:,0,jj,ii)*dzsoi(0) + h2osoi(:,1,jj,ii)*dzsoi(1))/(dzsoi(0) + dzsoi(1))  
 end do
 end do

 end if    ;root zone

 end if    ;readclm

 if (readpfl) then
 print("Reading ParFlow data ...")
;Groundwater level and unsaturated storage
;-----------------------------------
 fils_pflp  = systemfunc("ls " + diri_pfl + "/press/*.pfb")
 fils_pfls  = systemfunc("ls " + diri_pfl + "/satur/*.pfb")
 ntim       = dimsizes(fils_pflp)  ;number of history tapes

 ; Extract groundwater dimensions
 nx=new(1,integer)
 ny=new(1,integer)
 nz=new(1,integer)
 dx=new(1,double)
 dy=new(1,double)
 dz=new(1,double)
 EXT1::pfread_dim(fil_sol,nx,ny,nz,dx,dy,dz)

 ;Initialize
 soiltype    = new((/nz,ny,nx/),double)
 xslp        = new((/1,ny,nx/),double)
 yslp        = new((/1,ny,nx/),double)
 porosity    = new((/nz,ny,nx/),double)
 sp_storage  = new((/nz,ny,nx/),double)

 pressure = new((/nz,ny,nx/),double)
 satur    = new((/nz,ny,nx/),double)
 gwtd     = new((/ntim-1,ny,nx/),double)
 sss      = new((/ntim-1,ny,nx/),double)    ;Unsaturated Storage
 tws      = new((/ntim-1,ny,nx/),double)    ;Total Water Storage
 qo       = new((/ntim-1,2,ny,nx/),double)
 nm       = new((/ny,nx/),double)

 ;Estimate mass point depths
 var_dz          = new(nz,double)   ;m
 dz_3d           = new((/nz,ny,nx/),double)
 var_dz(0:9)     = (/0.02,0.03,0.05,0.07,0.13,0.20,0.30,0.50,0.70,1.00/)
 var_dz(10:nz-1) = 1.35
 zmassp          = new(nz,double)
 do k = 0, nz-1
   if (k.eq.0)
     zmassp(k) = 0.5*var_dz(k)
   else
     zmassp(k) = zmassp(k-1) + 0.5*(var_dz(k-1) + var_dz(k))
   end if
   dz_3d(k,:,:) = var_dz(nz-1 - k)
 end do

 zmassp    = zmassp(::-1)                ;reverse z_massp

 ;
 EXT2::pfread_var(fil_sol,nx,ny,nz,soiltype)
 EXT2::pfread_var(fil_xsl,nx,ny,1,xslp) 
 EXT2::pfread_var(fil_ysl,nx,ny,1,yslp) 
 EXT2::pfread_var(fil_por,nx,ny,nz,porosity)
 EXT2::pfread_var(fil_sps,nx,ny,nz,sp_storage)

 nm = manningscoefConst
 
 do fn = 1, ntim-1  ;skip the initial file to match with CLM
   fnamep = fils_pflp(fn)
   fnames = fils_pfls(fn)
   print(fn + " " + dx + " " + dy)
   EXT2::pfread_var(fnamep,nx,ny,nz,pressure)
   EXT2::pfread_var(fnames,nx,ny,nz,satur)
   EXT3::gwt_depth(nx,ny,nz,zmassp,pressure,gwtd(fn-1,:,:))
   print(xslp(0,50,50) + " " + yslp(0,50,50) + " " + max(pressure(nz-1,:,:)))
   EXT4::qoverland(nx,ny,dx,dy,nm,xslp(0,:,:),yslp(0,:,:),pressure(nz-1,:,:),qo(fn-1,:,:,:))

   ;TWS
   total_storage = dim_sum_n((satur*sp_storage*pressure + satur*porosity)*dx*dy*dz_3d,0)    ;m3
   tws(fn-1,:,:) = total_storage/(dx*dy)  ;m

   ;Unsaturated Storage
   ;Mask out saturated zone
   spress_1d  = ndtooned(pressure)
   satur_1d   = ndtooned(satur)
   bmask1d    = ind(spress_1d.gt.0)
   if (.not.all(ismissing(bmask1d))) then
     spress_1d(bmask1d) = 0.
     satur_1d(bmask1d)  = 0.
   end if
   pressure  = onedtond(spress_1d,(/nz,ny,nx/))
   satur     = onedtond(satur_1d,(/nz,ny,nx/))
   unsat_storage = dim_sum_n((satur*sp_storage*pressure + satur*porosity)*dx*dy*dz_3d,0)    ;m3
   sss(fn-1,:,:) = unsat_storage/(dx*dy)  ;m

   delete([/bmask1d/])
 end do 

 end if    ;readpfl
 
;Debug Output
;-----------------------------------
 foutnc = "debugTSMP_" + year + "_" + inst + ".nc"
 ;if (fileexists(foutnc)) then
   system("rm " + foutnc)
 ;end if
 fout = addfile(foutnc,"c")
;
 fout->lat2d = lat2d
 fout->lon2d = lon2d
 
 if (readclm) then
 rn_g   = (/fsa - fira - fgr/)                     ;available energy
 fsh    = (/fshv + fshg/)
 flh    = (/fctr+fcev+fgev/)
 t_et   = (/fctr/flh/)                             ;T/ET
 br0    = (/fsh/flh/)
 br0m   = where(rn_g.gt.0. ,  br0, br0@_FillValue)
 t_etm  = where(rn_g.gt.0. ,  t_et, t_et@_FillValue)

 fout->time  = time
 fout->PFT   = pft
 fout->LAI   = lai
 fout->Rn_G  = rn_g
 fout->T_ET  = t_etm
 fout->BR0   = br0m
 if (rootzone) then
   fout->SWR   = swr
   fout->SW5   = sw5 
 end if
 fout->H     = fsh
 fout->LE    = flh
 ;fout->fctr  = fctr
 ;fout->fcev  = fcev
 ;fout->fgev  = fgev
 ;fout->fshv  = fshv
 ;fout->fshg  = fshg

 end if     ;readclm

 if (readpfl) then
 fout->zsoi   = zmassp
 fout->STYP   = soiltype(::-1,:,:)
 fout->XSLP   = xslp(0,:,:)
 fout->YSLP   = yslp(0,:,:)
 fout->GWT    = gwtd
 fout->SSS    = sss
 fout->TWS    = tws
 fout->QO     = qo
 end if
end
