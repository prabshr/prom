;All variables preserve their native units
;scale/offset is only used for plotting
; Cumulative Frequency Altitude Diagram
; Based on Yuter and Houze 1995 (YH95)
; opt = True, generates CFADS for individual ensemble members
; opt = False,reads nc CFAD data for mulitple ensemble memebers to generate CFAD
; fcrop = True works with opt=False

;P. Shrestha
;Nov 18 2019
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$PROM/functions/cfads.ncl"
load "$PROM/functions/load_cfadColors.ncl"
load "$prom_root/functions/storm_tracker_nclf.ncl"
external EXT1 "$prom_root/functions/readdim.so"
external EXT2 "$prom_root/functions/readascii.so"

function read_cfad(cdiri,cens)
begin
 nc = dimsizes(cens)
 do nf = 0, nc-1
   fname  = cdiri + "/cfad_" + cens(nf) + ".nc"
   nc_fil = addfile(fname + ".nc","r")
   var    = nc_fil->Band1
   if (nf.eq.0) then
     cdim = dimsizes(var)
     vare = new((/nc,cdim(0),cdim(1),cdim(2),cdim(3),cdim(4)/),typeof(var))
   end if
   vare(nf,:,:,:,:,:) = var
 end do
 vare!0          = "ens"
 vare!1          = "time"
 vare!2          = "pv"
 vare!3          = "Z"
 vare!4          = "lat"
 vare!5          = "lon"
 vare&Z          = var&Z
 vare@_FillValue = -99
 vare@lat2d      = var@lat2d
 vare@lon2d      = var@lon2d

 return(vare)
end

function record_append(fns,fne,diri,vid,glon_idB,glat_idB,opt)
begin

  diri_output = diri + "/output"
  diri_input  = diri + "/input"
  fils_var    = systemfunc("ls " + diri_output + "/dbz-tm-1m_*")
  fils_elv    = systemfunc("ls " + diri_input + "/elevation_*")
  fil_lat     = systemfunc("ls " + diri_input + "/lat_*")
  fil_lon     = systemfunc("ls " + diri_input + "/lon_*")

  ;Read header information to extract nx,ny,nz
  nx=0 
  ny=0
  nz=0
  EXT1::readdim(nx,ny,nz,fils_elv(0))

  latPV = new((/1,1,ny,nx/),"float")
  lonPV = new((/1,1,ny,nx/),"float")
  hgtPV = new((/1,nz,ny,nx/),"float")
  EXT2::readascii(nx,ny,1,latPV,1,1,fil_lat)
  EXT2::readascii(nx,ny,1,lonPV,1,1,fil_lon)
  EXT2::readascii(nx,ny,nz,hgtPV,1,1,fils_elv(0))   ;Model mass point heights

  ;Extract lat-lon
  lat_c = latPV(0,0,:,:)
  lon_c = lonPV(0,0,:,:)

  ;Average height at iz
  ;Need to reverse the vertical levels
  mhgt  = hgtPV(0,::-1,ny/2,nx/2)
  mhgt  = (mhgt - mhgt(0))/1000.

  nt   = fne-fns + 1 
  varh =  new((/nt,4,nz,ny,nx/),"float")    ;Initialize variable
  varh@_FillValue  = -99
  varPV            = new((/8,nz,ny,nx/),"float")
  varPV@_FillValue = -99

  do nf = fns,fne
    it = nf-fns
    fname = fils_var(nf)
    print(nf + " " + fname)
    EXT2::readascii(nx,ny,nz,varPV,8,1,fname)

    ;Need to reverse the vertical levels
    varh(it,0,:,:,:) = varPV(1,::-1,:,:)
    varh(it,1,:,:,:) = varPV(2,::-1,:,:)
    varh(it,2,:,:,:) = varPV(3,::-1,:,:)
    varh(it,3,:,:,:) = varPV(4,::-1,:,:)

    Zh  = varh(it,0,:,:,:)
    Zv  = varh(it,1,:,:,:)
    KDP = varh(it,2,:,:,:)
    RHV = varh(it,3,:,:,:)
    varh(it,0,:,:,:) = where(Zh.lt.-5., varh@_FillValue, Zh)
    varh(it,1,:,:,:) = where(Zh.lt.-5., varh@_FillValue, Zh-Zv)
    varh(it,2,:,:,:) = where(Zh.lt.-5., varh@_FillValue, KDP)
    varh(it,3,:,:,:) = where(Zh.lt.-5., varh@_FillValue, RHV)

    if (opt@fstrack)
      wfac = opt@wfac   ;fraction of extent for search area
      wss  = opt@wss   ;SW fraction
      wse  = opt@wse ;NE fraction
      ext  = opt@ext
      ngridmax = opt@ngridmax

      ;Find location where W is maximum ....
      if (it.eq.0) then
        mpcorner    = (/glat_idB-ext(0),glon_idB-ext(1),glat_idB+ext(0),glon_idB+ext(1)/)
        glatlon = storm_tracker_nclf(Zh,lat_c,lon_c,wfac,wss,wse,ext,mpcorner)
      else
        glatlon = storm_tracker_nclf(Zh,lat_c,lon_c,wfac,wss,wse,ext,mpcornerw)
      end if
      glat_id    = glatlon(0)
      glon_id    = glatlon(1)

      mpcorner   = (/glat_id-wss*ext(0) ,glon_id-wss*ext(1) ,glat_id+wse*ext(0) ,glon_id+wse*ext(1)/)
      mpcornerw  = (/glat_id-wss*wfac*ext(0),glon_id-wss*wfac*ext(1),glat_id+wse*wfac*ext(0),glon_id+wse*wfac*ext(1)/)

      ji      = region_ind (lat_c,lon_c, mpcorner(0),mpcorner(2),mpcorner(1),mpcorner(3))
      ngrid   = ji(1)-ji(0)+1
      print("ny extent : " + ngrid)
      if (ngrid.gt.ngridmax) then
        ji(1) = ji(0) + ngridmax-1
        ngrid   = ji(1)-ji(0)+1
        print("ny extent reset : " + ngrid)
      end if
      ngrid   = ji(3)-ji(2)+1
      print("nx extent : " + ngrid)
      if (ngrid.gt.ngridmax) then
        ji(3) = ji(2) + ngridmax-1
        ngrid   = ji(3)-ji(2) + 1
        print("nx extent reset : " + ngrid)
      end if

      ;Cannot go beyond boundary
      if (ji(1).gt.ny-3) then
        ji(1) = ny-3
        ji(0) = ji(1) - ngrid + 1
        print("Y-Boundary reached ...")
      end if
      if (ji(3).gt.nx-3) then
        ji(3) = nx-3
        ji(2) = ji(3) - ngrid + 1
        print("X-Boundary reached ...")
      end if

      ;Wierd way of masking
      varh_c  = varh(it,:,:,:,:) 
      varh_c  = varh@_FillValue
      varh_c(:,:,ji(0):ji(1),ji(2):ji(3))  = varh(it,:,:,ji(0):ji(1),ji(2):ji(3))  
      varh(it,:,:,:,:) = varh_c
    end if
    delete([/Zh, Zv, KDP, RHV/])
    
 end do

 varh!0          = "time"
 varh!1          = "pv"
 varh!2          = "Z"
 varh!3          = "lat"
 varh!4          = "lon"
 varh&Z          = mhgt
 varh@_FillValue = -99
 varh@lat2d      = lat_c
 varh@lon2d      = lon_c

 fname = "cfad_" + opt@ens + ".nc"
 system("rm " + fname)
 fout = addfile(fname,"c")
 fout->Band1 = varh

 return(varh)
end


begin

;User Settings
;---------------------------------------------------------------------------
 ;pid       = 0 

 ;for PROM3 fopt=False only
 ;fnv          = (/0, 9,4 ,18,10,5/)    ;18 is 1530 UTC ;400 UTC
 ;fnsv         = (/0, 7,3 ,16, 8,4/)
 ;fnev         = (/3,10,6 ,19,11,7/)
 ;ensv         = (/7, 8,10,11, 12,19/)

 ydates    = (/"20150705","20160513","20170706"/)
 fnp       = (/11, 15, 18/)    
 fnsp      = (/ 9, 11,  16/)   
 fnep      = (/18, 20,  19/) 
 fcrop     = True ; spatial cropping
 ;
 fens      = (/6, 8, 11/)  

 ydate     = ydates(pid)
 fns       = fnsp(pid)
 fne       = fnep(pid)
 ens       = fens(pid)

 case      = "HET2483"
 diri      = "/daten04/pshrestha/shared-data/bonnRadar_" + ydate + "/Synth/" + case + "_" + ens 

 pvars     = (/"Z~B~H~N~","Z~B~DR~N~","K~B~DP~N~","~F8~r~B~~F10~hv~E~"/)
 punits    = (/"[dBZ]","[dB]","[deg/km]", "[-]"/)
 opt       = False ;- read netcdf

 wks_type     = "png"
 glon_idB    = 7.071663 ; Bonn Radar
 glat_idB    = 50.73052 ;

;Read netcdf files
;---------------------------------------------------------------------------
 if (opt) then

   opt@fstrack = False ;True
   opt@wfac = 0.15   ;fraction of extent for search area
   opt@wss  = 0.5    ;SW fraction
   opt@wse  = 1.+ opt@wss ;NE fraction
   opt@ext  = (/0.7,0.9/)
   opt@ngridmax = 180
   opt@ens   = ens
   vid = 0
   var =  record_append(fns,fne,diri,vid,glon_idB,glat_idB,opt)
 else
   if (pid.le.1) then
     nens = ens ;(/5,6,7,8,9/)
   else
     nens = (/11/)
     ;nens = (/7, 8,10,11, 19/)  ;dont change order, 3 indice crop
   end if
   cid  = pid + 1
   diri = "./scratch/data_prom" + cid + "/"
   vid  = 0
   var =  read_cfad(diri,nens)

   if (fcrop .and. pid.eq.0) then
      var(:,:,:,:,134:,0:100) = var@_FillValue
      var(:,:,:,:,0:65,:)    = var@_FillValue
   end if

   if (fcrop .and. pid.eq.1) then
      var(:,:,:,:,108:,:) = var@_FillValue
   end if
 
   if (fcrop .and. pid.eq.2) then
     var(:,:,:,:,0:47,:) = var@_FillValue
     var(:,:,:,:,116:,:) = var@_FillValue
   end if

 end if

; CFAD plot parameters
;--------------------------------------------------------
 copt          = True
 copt@YH95     = False  ;scale CFAD and remove data with low % of N
 copt@vid      = vid
 colors       = load_cfadColors(copt)
 minv         = colors@minv
 maxv         = colors@maxv
 binw         = colors@binw

 ndimv        = dimsizes(var)
 ndimv        = ndimv(::-1)
 nz           = ndimv(2) 
 ny           = ndimv(1) 
 nx           = ndimv(0) 
 print(nx + " " + ny + "" + nz)

 if (opt) then
   cvar = var(:,vid,:,:,:)
 else
   cvar = var(:,:,vid,:,:,:)
 end if
 zhh_cfad = cfads(var&Z,cvar,minv(vid),maxv(vid),binw(vid),pvars(vid),punits(vid),copt@YH95)

 vid      = 1
 copt@vid  = vid
 if (opt) then
   cvar = var(:,vid,:,:,:)
 else
   cvar = var(:,:,vid,:,:,:)
 end if
 zdr_cfad = cfads(var&Z,cvar,minv(vid),maxv(vid),binw(vid),pvars(vid),punits(vid),copt@YH95)

 vid      = 2
 copt@vid  = vid 
 if (opt) then
   cvar = var(:,vid,:,:,:)
 else
   cvar = var(:,:,vid,:,:,:)
 end if
 kdp_cfad = cfads(var&Z,cvar,minv(vid),maxv(vid),binw(vid),pvars(vid),punits(vid),copt@YH95)

 vid      = 3
 copt@vid  = vid
if (opt) then
   cvar = var(:,vid,:,:,:)
 else
   cvar = var(:,:,vid,:,:,:)
 end if
 rho_cfad = cfads(var&Z,cvar,minv(vid),maxv(vid),binw(vid),pvars(vid),punits(vid),copt@YH95)


;Plot
;---------------------------------------------------------------------------
 plot = new(4,"graphic")

 fname_out = "cfad_model_" + ydate + "_" + fns + "_" + fne + ".nc"

 wks = gsn_open_wks(wks_type,fname_out)
 colors =(/"white","black","white","gray30","gray40","gray50","gray61","gray75", \
          "gray80","gray90","paleturquoise4","paleturquoise3","orchid4","orchid3","orchid2","orchid1","orchid", \
         "royalblue4","darkturquoise","yellow"/)

 gsn_define_colormap(wks,colors)

 fncps_no = ":F21:"
 fncps_no = 21
 fncps_hgt = 0.025 ;2

 res                            = True           

 ;gsn Resources
 ;---------------------------------------------------
 res@gsnDraw                    = False
 res@gsnFrame                   = False
 res@gsnMaximize                = False
 res@gsnLeftStringFontHeightF   = fncps_hgt
 res@gsnRightStringFontHeightF  = fncps_hgt
 res@gsnLeftString              = "Model " ;+ timestr
; res@gsnRightString             = zhh_cfad@long_name

 ;stretched grid Resources
 ;---------------------------------------------------
 res@trGridType            = "TriangularMesh"
 res@trYMaxF               = 15.0 ;km
; res@sfXArray              = zhh_cfad&binval
; res@sfYArray              = zhh_cfad&height  ;km 

 ;tm Resources
 ;---------------------------------------------------
 res@tmXBLabelFont              = fncps_no
 res@tmYLLabelFont              = fncps_no
 res@tmYLLabelFontHeightF       = fncps_hgt
 res@tmXBLabelFontHeightF       = fncps_hgt
 res@tmYROn                     = False
 res@tmXTOn                     = False
 res@tmXBLabelFontColor         = "black"
 res@tmYLLabelFontColor         = "black"
 res@tmXBLabelFontThicknessF    = 2.0
 res@tmYLLabelFontThicknessF    = 2.0

 ;ti Resources
 ;---------------------------------------------------
 res@tiXAxisFont                = fncps_no
 res@tiYAxisFont                = fncps_no
 res@tiXAxisFontHeightF         = fncps_hgt
 res@tiYAxisFontHeightF         = fncps_hgt
 res@tiXAxisString              = zhh_cfad@units
 res@tiYAxisString              = "Height [km]"

 ;cn Resources
 ;---------------------------------------------------
 res@cnFillOn                   = True
 ;res@cnFillMode                 = "RasterFill"
 res@cnLinesOn                  = False
 res@cnInfoLabelOn              = False   ; Turn off annotations
 res@cnLineLabelsOn             = False
 res@cnLevelSelectionMode       = "ExplicitLevels"
 res@cnLevels                   = colors@cnLevels

 ;ln Resources
 ;---------------------------------------------------
 res@lbLabelBarOn               = False
 res@lbLabelAutoStride          = True
 res@lbOrientation              = "vertical"
 res@lbTitlePosition            = "Right"
 res@lbTitleString              = "%"
 res@lbTitleExtentF             = 0.15
 res@lbTitleAngleF              = 90
 res@lbLabelFontHeightF         = 0.02
 res@pmLabelBarWidthF           = 0.05
 ;res@pmLabelBarHeightF          = 0.35

 ;vp Resources
 ;---------------------------------------------------
 res@vpHeightF                = 0.6
 res@vpWidthF                 = 0.45
; res@vpYF                     = 0.9
; res@vpXF                     = 0.15

 ;plot specific
 ;---------------------------------------------------
 var_cfad                  = zhh_cfad
 res@gsnRightString        = var_cfad@long_name
 res@sfXArray              = var_cfad&binval
 res@sfYArray              = var_cfad&height  ;km
 res@tiXAxisString         = var_cfad@units
 plot(0)     = gsn_csm_contour(wks,var_cfad,res)
 res0        = res
 res0@cnFillOn                   = False
 res0@cnLinesOn                  = True
 res0@cnLineThicknessF           = 3.0
 delete(res0@cnLevels)
 res0@cnLevels                   = 5.
 dummy0   = gsn_csm_contour(wks,zhh_cfad,res0)
 overlay(plot(0),dummy0)

;----------------------------------------------
 res@gsnLeftString  = ""
 res@tiYAxisString  = ""
 res@tmYLLabelsOn   = False

 delete(var_cfad)
 delete(res@sfXArray)

 var_cfad                  = zdr_cfad
 res@gsnRightString        = var_cfad@long_name
 res@sfXArray              = var_cfad&binval
 res@sfYArray              = var_cfad&height  ;km
 res@tiXAxisString         = var_cfad@units
 plot(1)     = gsn_csm_contour(wks,var_cfad,res)

 delete(var_cfad)
 delete(res@sfXArray)
 var_cfad                  = kdp_cfad
 res@gsnRightString        = var_cfad@long_name
 res@sfXArray              = var_cfad&binval
 res@sfYArray              = var_cfad&height  ;km
 res@tiXAxisString         = var_cfad@units
 plot(2)     = gsn_csm_contour(wks,var_cfad,res)

 delete(var_cfad)
 delete(res@sfXArray)
 var_cfad                  = rho_cfad
 res@gsnRightString        = var_cfad@long_name
 res@sfXArray              = var_cfad&binval
 res@sfYArray              = var_cfad&height  ;km
 res@tiXAxisString         = var_cfad@units
 res@lbLabelBarOn          = True
 plot(3)     = gsn_csm_contour(wks,var_cfad,res)

; Add annotation
 txres               = True                      ; text mods desired
 txres@txFontHeightF = 0.015                    ; text font height
 txres@txJust        = "CenterLeft"              ; Default is "CenterCenter".

 gsn_text_ndc(wks,"b)",.01,.67,txres)


 resP = True
 resP@gsnMaximize    = True
 resP@gsnPanelXF = (/-1,0.28,0.51,0.74/)
; resP@gsnPanelDebug = True
 gsn_panel(wks,plot,(/1,4/), resP)

end
