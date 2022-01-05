;Uses maximum vertical velocity to track the strom propagation using
;a updated cropped region defined by "glat_id,glon_id and ext"
;Plots total hydrometeor, vertical velocity and wind vectors at "fpressure"
;      cross-section of hydrometeor, wind vectors and cold pool at "Ipangle"
;      vertical profile of hydrometor 

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

load "$PROM/functions/cosmo_time_nclf.ncl"
load "$PROM/functions/calcrho_cos_nclf.ncl"
load "$PROM/functions/storm_tracker_nclf.ncl"
load "$PROM/functions/radarline_nclf.ncl"
;Log: Nov. 1 2021, Addition of ART directory

begin

;-----------------------------------------------------------------------------------------
;#1 USER SETTINGS
;-----------------------------------------------------------------------------------------
 prom=4
 inst=6
 opt=1
 lnd = "HET"
 itype = "ART"
 itype_gscp  = lnd + itype

 if (prom.eq.1) then
   diri     = "/daten04/pshrestha/shared-data/bonnRadar_20150705"
   ;ixp      = 1  ;which hydrometeor to explore 
   fns      = 11 ;21 ; 16 ; 21 ;10
   fne      = 11 ;21
   dns      = 0  ; offset for dbz files
   xlon = (/50.7, 51.1/)
   offs = -4
 else if (prom.eq.2) then
   diri    = "/daten04/pshrestha/shared-data/bonnRadar_20160513"
   fns      = 488 ;19
   fne      = 48 ;55
 else if (prom.eq.3) then
   diri    = "//daten04/pshrestha/shared-data/bonnRadar_20170706"
   fns     = 0 ;12+12
   fne     = 0 ;48
 else if (prom.eq.4) then
   diri    = "/daten04/pshrestha/shared-data/bonnRadar_20150705"
   fns     = 16 
   fne     = 16 
   xlon = (/50.7, 51.1/)
   ;offs    = 0
 end if
 end if
 end if
 end if

 ext         = (/0.45,0.6/)* 1.0 ;2.0 ;3.0   ;*1.2
 ngridmax    = 90*3  ;110 ;90
 iopt        = True ;False    ;always True for accurate WRF interpolatoin
 fstrack     = True ;False

 glon_idB    = 7.071663 ; Bonn Radar
 glat_idB    = 50.73052 ;
 extB        = 0.91     ;extent in degrees

 Ipangles    = (/0.,0.,0.,27./) ;(/25.0,145.,170./)   ;consistent with plot_ppi
 qprom       = (/4.0, 3.0, 3.0,3.0/)    ;max for vertical profile

 hloc        = 70 ;70      ;max 99, height of cross-section for plotting
 uzmask      = False
 scale_uz    = 1.0
 hmax        = 15.0        ;km
 fpressure   = 550. ;350.   ;hPa, for pressure plot

 mpcorner    = (/glat_idB-ext(0),glon_idB-ext(1),glat_idB+ext(0),glon_idB+ext(1)/)
 
 hydro_colors     = (/"green","blue","skyblue","purple","grey40","orange"/)
 hydro_dash       = (/0, 0, 1, 1, 1, 1/)
 hydro_thickness  = hydro_dash + 4.

 subevent    = itype_gscp + "/runE_" + inst + "_" + itype_gscp
 diri_cos    = diri + "/" + "/" + subevent + "/cosout/"
 diri_cosI   = diri + "/" + "/" + subevent + "/cosout/ivr/"

;-----------------------------------------------------------------------------------------
;#2 Set variables to read and specify scales or offsets
;-------------------------------------------------------------------------------------
 hydrometeors_qx = (/"QC","QR","QI","QS","QG","QH"/)
 hydrometeors_qn = (/"NCCLOUD","NCRAIN","NCICE","NCSNOW","NCGRAUPEL","NCHAIL"/)

 if (opt.eq.0) then
   hydrometeors = hydrometeors_qn
   scale_varh   = 1/1.E3
   qcnlevs      = (/1.1, 1.E1, 1.E2, 1.E3, 1.E4, 1.E5, 1.E6, 1.E7/)
   hydunits     = "L~S~-1~N~"
   hydname      = "Number Concentration"
 else if (opt.eq.1) then
   hydrometeors = hydrometeors_qx
   scale_varh   = 1000.  ;kg/kg to g/kg
   qcnlevs      = (/0.01,0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0/)
   hydunits     = "g/m~S~3~N~"
   hydname      = "Sp. Hydrometeor Content"
 else if (opt.eq.2) then
   x3           = 1./3.
   hydrometeors = hydrometeors_qx
   a_geo        = (/0.124    , 0.124  , 0.835 , 2.4    , 0.142     , 0.1366/)  ;mkg-b
   b_geo        = (/x3       , x3     , 0.39  , 0.455  , 0.314     , x3/)      ;-
   xmin         = (/4.2E-15  , 2.6E-10, 1.E-12, 1.E-10, 1.E-9, 2.6E-9/)        ;kg
   xmax         = (/2.6E-10  , 3.0E-6 , 1.E-6 , 2.E-5 , 5.E-4, 5.E-4/)         ;g
   qmax         = 30.
   qmin         = 0.001 
   qcnlevs       = (/0.001, 0.1, 1.0, 5.0, 10., 15., 20., 25., 30./)
   scale_varh   = 1000.   ; m to  mm
   hydunits     = "mm"
   hydname      = "Max. Volume Diameter"
 else
   print("code not written for opt = " + opt)
   exit
 end if
 end if
 end if
;

;-----------------------------------------------------------------------------------------
;#3  READ NC FILES
;-------------------------------------------------------------------------------------
 nc_cosI = addfile(diri_cosI +"lfff00000000c.nc","r")
 lat = nc_cosI->lat
 lon = nc_cosI->lon
 topo = nc_cosI->HSURF(0,:,:)
; topo = smth9_Wrap(topo, 0.50, 0.25, False)

 hhl       = nc_cosI->HHL(0,::-1,:,:)
 ndimhhl   = dimsizes(hhl)
 nz        = ndimhhl(0)-1
 ny        = ndimhhl(1)
 nx        = ndimhhl(2)
 zmodel    = 0.5*(hhl(0:nz-1,:,:) + hhl(1:nz,:,:))/1000.  ;km
;
 fils_cos  = systemfunc("ls " + diri_cos + "lf*")
 nfils   = dimsizes(fils_cos)

 do fn  = fns, fne ;nfils -1      ;BIG LOOP

 nc_cos = addfile(fils_cos(fn),"r")
 time   = nc_cos->time
 ;
 str_time = cosmo_time_nclf(time,1)  ; 0 for date

 print("----------------------------------------------")
 print(fn + " " + str_time) 

;------------------------------------------------------------------------------
;#4 Prognostic state variables
;------------------------------------------------------------------------------
 tair            = nc_cos->T(0,::-1,:,:)     ;K
 pressure        = nc_cos->P(0,::-1,:,:)     ;Pa
 qv              = nc_cos->QV(0,::-1,:,:)    ;kg/kg
 ux              = nc_cos->U(0,::-1,:,:)
 uy              = nc_cos->V(0,::-1,:,:)
 uz              = nc_cos->W(0,::-1,:,:)
 uz_hhl          = uz(1:80,:,:)*0.
 do kk = 0,79
   uz_hhl(kk,:,:)  = 0.5*(uz(kk,:,:) + uz(kk+1,:,:))*scale_uz
 end do

; Fluxes
 tot_prec        = nc_cos->TOT_PREC(0,:,:)   ;kg/m2
 hail_gsp        = nc_cos->HAIL_GSP(0,:,:)   ;kg/m2

;Remove accumulated precipitation
 if (fn.eq.fns) then
   tot_prec0 = tot_prec
   hail_gsp0 = hail_gsp
 else
   tot_prec  = tot_prec - tot_prec0
   hail_gsp  = hail_gsp - hail_gsp0
 end if

 tot_prec0   = nc_cos->TOT_PREC(0,:,:)   ;kg/m2
 hail_gsp0   = nc_cos->HAIL_GSP(0,:,:)   ;kg/m2 

;-----------------------------------------------------------------------------------------
;#5 UPDATE STORM REGION 
;-------------------------------------------------------------------------------------
 if (fstrack) then

 wfac = 0.15   ;fraction of extent for search area
 wss  = 0.4    ;SW fraction
 wse  = 1.+wss ;NE fraction

 ;Find location where W is maximum ....
 if (fn.eq.fns) then
   glatlon = storm_tracker_nclf(uz_hhl,lat,lon,wfac,wss,wse,ext,mpcorner)
 else
   glatlon = storm_tracker_nclf(uz_hhl,lat,lon,wfac,wss,wse,ext,mpcornerw)
 end if
 glat_id    = glatlon(0)
 glon_id    = glatlon(1) 

 print("Cross-section through " + glon_id)

 print("Updating mpcorner and mpcornerw")
 mpcorner   = (/glat_id-wss*ext(0)     ,glon_id-wss*ext(1)     ,glat_id+wse*ext(0)     ,glon_id+wse*ext(1)/)
 mpcornerw  = (/glat_id-wss*wfac*ext(0),glon_id-wss*wfac*ext(1),glat_id+wse*wfac*ext(0),glon_id+wse*wfac*ext(1)/)

 else
 glat_id = glat_idB
 glon_id = glon_idB
 end if
;-----------------------------------------------------------------------------------------
;#6 DERIVE HYDROMETEORS 
;-------------------------------------------------------------------------------------

; Compute ambient air density 
 qc            = nc_cos->QC(0,::-1,:,:)    ;kg/kg
 qr            = nc_cos->QR(0,::-1,:,:)    ;kg/kg
 qi            = nc_cos->QI(0,::-1,:,:)    ;kg/kg
 qs            = nc_cos->QS(0,::-1,:,:)    ;kg/kg
 qg            = nc_cos->QG(0,::-1,:,:)    ;kg/kg
 qh            = nc_cos->QH(0,::-1,:,:)    ;kg/kg
 tlc           = qr + qi + qs + qg + qh    ;Precipitation loading
 rho           = calcrho_cos_nclf(tair,pressure,qv,qc,tlc) ;air density (kgm-3)
 nh           = dimsizes(hydrometeors)
 varh         = new((/nh,nz,ny,nx/),"float")    ;Initialize variable for plotting
 varh@units   = hydunits
 varh@long_name = hydname
 varh@_FillValue = -999.
 do ih = 0, nh-1
   if (opt.eq.0) then
     varh(ih,:,:,:) = (/nc_cos->$(hydrometeors_qn(ih))$(0,::-1,:,:)/)
     varh(ih,:,:,:) = (/scale_varh * varh(ih,:,:,:) * rho/)    ;(kg-1) to (L-1)
           ;Nt    = 1/1000. * Nd * rho            ; number density to number concentration
   else if (opt.eq.1) then
     varh(ih,:,:,:) = (/nc_cos->$(hydrometeors_qx(ih))$(0,::-1,:,:)/)
     varh(ih,:,:,:) = (/scale_varh * varh(ih,:,:,:) * rho/)    ;(kg/kg) to (g/m3)
           ;Qt    = 1000. * Qd * rho
   else if (opt.eq.2) then
     qtemp      = (/nc_cos->$(hydrometeors_qx(ih))$(0,::-1,:,:)/)
     Ntemp      = (/nc_cos->$(hydrometeors_qn(ih))$(0,::-1,:,:)/)
     qtemp      = where(qtemp.lt.0.0001, 0., qtemp)
     Ntemp      = where(Ntemp.eq.0, varh@_FillValue, Ntemp)        
     xmass      = (/qtemp/Ntemp/)   ; kg
     ;Clipping
     xmass      = where(xmass.gt.xmin(ih),xmass,0.)
     xmass      = where(xmass.gt.xmax(ih),xmax(ih),xmass)
     ;diameter mass relationship         
     varh(ih,:,:,:) = (/scale_varh* (a_geo(ih) * (xmass)^b_geo(ih))/)   ;mm
   end if
   end if
   end if
 end do


;-------------------------------------------------------------------------------
;#7 Crop data and vertical plane for interpolation
;-------------------------------------------------------------------------------
  ji      = region_ind (lat,lon, mpcorner(0),mpcorner(2),mpcorner(1),mpcorner(3))
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

 if (prom.eq.4 .and. fns.eq.16) then
   print("Hardwire setting")
  ;ART PROM setting
  glat_id = 50.8799
  glon_id = 7.14264 +0.09
  mpcorner = (/50.6999,6.90264,51.5099,7.98264/)
  ji = (/145,226,144,212/)
  end if

  print("Crop indices ")
  print(" " + ji)
 
  lat_c    = lat(ji(0):ji(1),ji(2):ji(3))
  lon_c    = lon(ji(0):ji(1),ji(2):ji(3))
  print("Cropped Extent")
  print("X " + min(lon_c) + " " + max(lon_c))
  print("Y " + min(lat_c) + " " + max(lat_c))
  zmodelc = zmodel(:,ji(0):ji(1),ji(2):ji(3))
  tmp_dim = dimsizes(lon_c)

  ; Obtain the radar line of view
  Ipangle = Ipangles(prom-1)
  plane = radarline_nclf(Ipangle,lat_c,lon_c,glat_id,glon_id,extB,iopt)
  ;;CPS
  plane(1) = plane(1) + offs

  glon_p = plane@X
  glat_p = plane@Y

  ;;0 : SN , 90 :WE
  lon_plane  = wrf_user_intrp2d(lon_c,plane,Ipangle,iopt)
  lat_plane  = wrf_user_intrp2d(lat_c,plane,Ipangle,iopt)
  hgt_plane  = wrf_user_intrp3d(zmodelc,zmodelc,"v",plane,Ipangle,iopt)

;
;crop
  topo_c      = topo(ji(0):ji(1),ji(2):ji(3))
  tk_c        = tair(:,ji(0):ji(1),ji(2):ji(3))
  qv_c        = qv(:,ji(0):ji(1),ji(2):ji(3))
  tlc_c       = tlc(:,ji(0):ji(1),ji(2):ji(3))
  ps_c        = pressure(:,ji(0):ji(1),ji(2):ji(3))
  ux_c        = ux(:,ji(0):ji(1),ji(2):ji(3))
  uy_c        = uy(:,ji(0):ji(1),ji(2):ji(3))
  uz_hhl_c    = uz_hhl(:,ji(0):ji(1),ji(2):ji(3))
  varh_c      = varh(:,:,ji(0):ji(1),ji(2):ji(3))
  prec_c      = tot_prec(ji(0):ji(1),ji(2):ji(3))
  hail_c      = hail_gsp(ji(0):ji(1),ji(2):ji(3))

  vx_c        = varh_c(ixp, :, :,:) ;dim_sum_n(varh_c,0)             ;total hydrometeor 
  ; Set precip to zero for fn=0
  if (fn.eq.fns) then
    prec_c = 0.
    hail_c = 0.
  end if

;Derived variables
  ;te_c        = wrf_eth ( qv_c, tk_c, ps_c )        ; theta-e, K
  th_c        = tk_c*(1.E5/ps_c)^(2./7.)             ; theta  , K
  te_c        = th_c*(1. + 0.608*qv_c - tlc_c)       ; density potential temp., K (DH2017)
;-------------------------------------------------------------------------------
;#8 Interpolate Cropped Data in Vertical Plane
;-------------------------------------------------------------------------------

  tc_plane    = (/wrf_user_intrp3d(tk_c,zmodelc,"v",plane,Ipangle,iopt) -273.15/)  ;C
  te_plane    = (/wrf_user_intrp3d(te_c,zmodelc,"v",plane,Ipangle,iopt) -273.15/)  ;C
  ux_plane    = (/wrf_user_intrp3d(ux_c,zmodelc,"v",plane,Ipangle,iopt)/)
  uy_plane    = (/wrf_user_intrp3d(uy_c,zmodelc,"v",plane,Ipangle,iopt)/)
  uz_plane    = (/wrf_user_intrp3d(uz_hhl_c,zmodelc,"v",plane,Ipangle,iopt)/)
  ndimp       = dimsizes(tc_plane)

  print(max(uz_plane) + "U Z MAX")
  varh_plane    = new((/nh,ndimp(0),ndimp(1)/),"float")
  do ih = 0, nh-1
    varh_plane(ih,:,:)    = (/wrf_user_intrp3d(varh_c(ih,:,:,:),zmodelc,"v",plane,Ipangle,iopt)/)
  end do

; Compute perturbation for THETAE to diagnose cold pool
  nplane     = dimsizes(te_plane)
  print("Masking density pot. temp.  1 km a.g.l.")
  hgt_plane_max = 1.0 ;
  do kk = 0, nplane(0)-1
    if (.not.all(ismissing(hgt_plane(kk,:)))) then
      te_plane(kk,:)   = te_plane(kk,:) - avg(te_plane(kk,:))
      hgt_avg_plane = avg(hgt_plane(kk,:))
      if (hgt_avg_plane.ge.hgt_plane_max) then
        te_plane(kk,:) = 0.
      end if
    end if
  end do
;
  x_plane               = lat_plane ;lon_plane

  tc_plane!0            = "height"
  tc_plane!1            = "longitude" 
  tc_plane&longitude    = (/x_plane/)
  tc_plane@long_name    = "T"
  tc_plane@units        = "C"

  te_plane!0            = "height"
  te_plane!1            = "longitude" 
  te_plane&longitude    = (/x_plane/)
  te_plane@long_name    = "THETAE'"
  te_plane@units        = "C"

  varh_plane!1            = "height"
  varh_plane!2            = "longitude" 
  varh_plane&longitude    = (/x_plane/)
  varh_plane@long_name    = hydname 
  varh_plane@units        = hydunits

  ux_plane!0         = "height"
  ux_plane!1         = "longitude" 
  ux_plane&longitude = (/x_plane/)

  uz_plane!0         = "height"
  uz_plane!1         = "longitude" 
  uz_plane&longitude = (/x_plane/)

;-------------------------------------------------------------------------------
;#9 Extract sounding data along cross-section
;-------------------------------------------------------------------------------
    vr_sum = dim_max_n(uz_plane,0)
    sloc = ind(vr_sum.eq.max(vr_sum))
    print("Sounding Location " + sloc)
    if (ismissing(sloc).or.sloc.eq.0.or.sloc.ge.dimsizes(lon_plane)-1) then
       sloc = tointeger(dimsizes(lon_plane)/2.)
    end if
    hloc = 80
    print("Extent of C-S " + hloc +  " Max: 99")
    ht_v = hgt_plane(0:hloc,sloc)     ;/1000.
    ht_v(0) = 0.
    sloc0 = sloc-1
    sloc1 = sloc+1
    tc_v = dim_avg_n(tc_plane(0:hloc,sloc0:sloc1),1)
    varh_v = dim_avg_n(varh_plane(:,0:hloc,sloc0:sloc1),2)
    ux_v = dim_avg_n(ux_plane(0:hloc,sloc0:sloc1),1)
    uy_v = dim_avg_n(uy_plane(0:hloc,sloc0:sloc1),1)
    uz_v = dim_avg_n(uz_plane(0:hloc,sloc0:sloc1),1)

    ;For correct profile plot
    do ih = 0, nh-1
      varh_v(ih,:) = mask(varh_v(ih,:),varh_v(ih,:).eq.0, False)
    end do 
 
;-------------------------------------------------------------------------------
;#10 Interpolate Cropped Data in horizontal at fixed pressure
;-------------------------------------------------------------------------------
;
   ps_c          = ps_c*0.01
   ux_p          = wrf_user_intrp3d(ux_c,ps_c,"h",fpressure,0.,False)
   uy_p          = wrf_user_intrp3d(uy_c,ps_c,"h",fpressure,0.,False)
   vx_p          = wrf_user_intrp3d(vx_c,ps_c,"h",fpressure,0.,False) 
   wx_p1         = wrf_user_intrp3d(uz_hhl_c,ps_c,"h",fpressure,0.,False)
   wx_p2         = wrf_user_intrp3d(uz_hhl_c,ps_c,"h",700.,0.,False)
   te_p          = wrf_user_intrp3d(te_c,ps_c,"h",850.,0.,False)
   qv_p          = wrf_user_intrp3d(qv_c,ps_c,"h",850.,0.,False)*1000.  

   ux_p@lat2d    = lat_c
   ux_p@lon2d    = lon_c
   uy_p@lat2d    = lat_c
   uy_p@lon2d    = lon_c

   wx_p1@lat2d    = lat_c
   wx_p1@lon2d    = lon_c
   wx_p2@lat2d    = lat_c
   wx_p2@lon2d    = lon_c

   vx_p@lat2d    = lat_c
   vx_p@lon2d    = lon_c

   topo_c@lat2d  = lat_c
   topo_c@lon2d  = lon_c

;------------------------------------------------------------------------------------------
;#11 Data Storage
;------------------------------------------------------------------------------------------
  if (fn.eq.fns) then
    latout = new((/nfils,tmp_dim(0),tmp_dim(1)/),"float")
    lonout = new((/nfils,tmp_dim(0),tmp_dim(1)/),"float")
    Band1 = new((/nfils,tmp_dim(0),tmp_dim(1)/),"float")
    Band2 = new((/nfils,tmp_dim(0),tmp_dim(1)/),"float")
    Band3 = new((/2,nfils,tmp_dim(0),tmp_dim(1)/),"float")
    Band4 = new((/nfils,tmp_dim(0),tmp_dim(1)/),"float")
    Band5 = new((/nfils,tmp_dim(0),tmp_dim(1)/),"float")
    Band6 = new((/nfils,2/),"float")
    Band1@_longname = "TOT_PREC"
    Band2@_longname = "HAIL_GSP"
    Band3@_longname = "UZ-500-700hPa"
    Band4@_longmame = "THETAV-850hPa"
    Band5@_longmame = "QV-850hPa"
    Band6@_longname = "UZ-track"
  end if

  latout(fn,:,:) = lat_c
  lonout(fn,:,:) = lon_c
  Band1(fn,:,:) = prec_c
  Band2(fn,:,:) = hail_c
  Band3(0,fn,:,:) = wx_p1
  Band3(1,fn,:,:) = wx_p2
  Band4(fn,:,:) = te_p
  Band5(fn,:,:) = qv_p
  Band6(fn,0)   = glat_id(0)
  Band6(fn,1)   = glon_id(0) 

  print("Rain " + " " + sum(prec_c) + "  Hail " + sum(hail_c) + " mm")
 
;-----------------------------------------------------------------------------------------
;#12  Plotting Resource
;-----------------------------------------------------------------------------------------
  plot = new(3,graphic)
  wks_type = "png"
  fname_out = "qx_prom_"+ prom + "_" +  itype_gscp + "_inst" + inst + "_opt_" + opt + "_" + sprinti("%0.4i",fn) + "_" + ixp
  wks = gsn_open_wks(wks_type,fname_out)

  gsn_define_colormap(wks,"precip2_17lev")

  fncps_no = ":F21:"
  fncps_no = 21
  fncps_hgt = 0.02

;-----------------------------------------------------------------------------------------
;#11a  Pressure plot resources
;-----------------------------------------------------------------------------------------

 res                            = True                 ; plot mods desired
 res@gsnDraw                    = False
 res@gsnFrame                   = False
 restopo                        = res
 resln                          = res
; res@gsnAddCyclic               = False                ; regional data
 res@gsnLeftString              = str_time
 res@gsnRightString             = sprintf("%0.0f",fpressure) + " hPa"
 res@gsnMaximize                = False
 res@gsnLeftStringFontHeightF   = fncps_hgt
 res@gsnRightStringFontHeightF  = fncps_hgt

 resw                           = res

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

 res@tiXAxisFont                = fncps_no
 res@tiYAxisFont                = fncps_no
 res@tiXAxisFontHeightF         = fncps_hgt
 res@tiYAxisFontHeightF         = fncps_hgt
 rescs                          = res
 resV                           = res

 res@pmTickMarkDisplayMode      = "Always"             ; turn on tick marks
 res@pmLabelBarWidthF           = 0.1
 res@pmLabelBarHeightF          = 0.5

;Map
 ;res@mpCenterLonF                 = lon@Longitude_of_southern_pole           ;lon of centre of map projection system
 ;res@mpCenterLatF                 = 90.0 + lat@Latitude_of_southern_pole         ;lat of centre of map projection system

 res@mpLimitMode                  = "Corners"           ; method to zoom
 res@mpLeftCornerLatF             = mpcorner(0) ;28; lat(0,0)
 res@mpLeftCornerLonF             = mpcorner(1) ;76 ;lon(0,0)
 res@mpRightCornerLatF            = mpcorner(2) ;32 ;lat(ndim(0)-1,ndim(1)-1)
 res@mpRightCornerLonF            = mpcorner(3) ;82 ;lon(ndim(0)-1,ndim(1)-1)
 res@mpOutlineOn                  = True
 res@mpOutlineBoundarySets        = "National"
 res@mpDataSetName                = "Earth..4"    ; database for non-USA divisions
 res@mpDataBaseVersion            = "MediumRes"   ; Medium resolution database
; res@mpOutlineSpecifiers          = "India:Uttaranchal"
 res@mpGridAndLimbOn              = True
 res@mpNationalLineThicknessF     = 3.0
 res@mpGeophysicalLineThicknessF  = 3.0
 res@mpProvincialLineThicknessF   = 2.0
 res@mpPerimOn                    = True
 res@mpFillOn                     = False

 ;Contours
 res@cnFillOn                   = True
 res@cnLinesOn                  = True ;False  ; False
 res@cnInfoLabelOn              = False   ; Turn off annotations
 res@cnLineLabelsOn             = False
 res@cnLevelSelectionMode       = "ExplicitLevels"
 rescs@cnLevels                      = qcnlevs
; res@cnLevels                   = min(qcnlevs) ;
 res@cnFillColors               = (/"white","slategrey"/)

 ;Labels
 res@lbLabelAutoStride          = True
 res@lbOrientation              = "vertical"
 res@lbTitleString              = hydunits 
 res@lbLabelFontHeightF         = 0.02
 res@pmLabelBarWidthF           = 0.05
 res@pmLabelBarHeightF          = 0.35

 ;Viewport
 res@vpHeightF                = 0.4
 res@vpWidthF                 = 0.4
 res@vpYF                     = 0.94
 res@vpXF                     = 0.15

 res@gsnScalarContour         = True

 res@tiXAxisString                 = ""
 res@tiYAxisString                 = ""

 ;Vectors
 res@vcRefMagnitudeF          = 25.0             ; define vector ref mag
 res@vcRefLengthF             = 0.045            ; define length of vec ref
 res@vcGlyphStyle             = "CurlyVector"    ; turn on curly vectors
 res@vcMinDistanceF           = 0.032            ; thin vectors
 res@vcRefAnnoOrthogonalPosF  = -.3             ; move ref vector down
 res@vcLineArrowColor         = "azure3" ;grey40"
 res@vcLineArrowThicknessF    = 2.
 res@vcGlyphOpacityF          = 0.7
 res@vcRefAnnoOn              = True

 plot(0)                      = gsn_csm_vector_scalar_map(wks,ux_p,uy_p,vx_p, res)

 ;Contours
 resw@cnFillOn                = False 
 resw@cnLinesOn               = True 
 resw@cnInfoLabelOn           = False   ; Turn off annotations
 resw@cnLineLabelsOn          = False
 resw@cnLineColor             = "red"
 resw@cnLineThicknessF        = 3.0
 resw@gsnLeftString          = ""
 resw@gsnRightString         = ""
 resw@cnLevelSelectionMode       = "ExplicitLevels"
 resw@cnLevels                   = ispan(5,40,5)*1. 
 resw@gsnContourNegLineDashPattern = 1
 print("500 hpa" + min(wx_p1) + " " + max(wx_p1))
 print("700 hpa" + min(wx_p2) + " " + max(wx_p2))
 
;CPS
; dumw                         = gsn_csm_contour(wks,wx_p1, resw)
; overlay(plot(0),dumw)

; resw@cnLineThicknessF        = 1.5
; resw@cnLineColor             = "orange"
; dumw1                        = gsn_csm_contour(wks,wx_p2, resw)
; overlay(plot(0),dumw1)
 delete(resw@cnLevels)

 resw@cnLineThicknessF        = 3.0
 resw@cnLineColor = "blue"
 resw@cnLevels  = ispan(2,10,1)*-1. 
;CPS
; dumw2                        = gsn_csm_contour(wks,wx_p1, resw)
; overlay(plot(0),dumw2) 

; resw@cnLineThicknessF        = 1.5
; resw@cnLineColor = "steelblue"
; dumw3                        = gsn_csm_contour(wks,wx_p2, resw)
; overlay(plot(0),dumw3)

 ;Topo Overlay

 ;Contour
 restopo@cnFillOn               = True
 restopo@cnLinesOn              = True
 restopo@cnInfoLabelOn          = False
 restopo@cnLineLabelsOn         = False
 restopo@cnLineThicknessF       = 2.
 restopo@cnLineDashPattern      = 8
 restopo@cnLevelSelectionMode   = "ExplicitLevels"
 restopo@cnLevels               = (/100., 300., 500./)
 restopo@cnMonoFillColor        = True
 restopo@cnMonoFillPattern      = False
 restopo@cnMonoFillScale        = False
 restopo@cnFillPatterns         = (/17,4,3,4,3,4,3,4,3,4/)
 restopo@cnFillScales           = (/10.,.45,.40,.35,.25,.20,.18,.15,.12,.10/)
 restopo@cnLineColor            = "black"
 restopo@cnFillDrawOrder        = "PreDraw"

 restopo@lbLabelBarOn           = False

 restopo@gsnLeftString          = ""
 restopo@gsnRightString         = ""
 ;dumtopo                        = gsn_csm_contour(wks,topo_c, restopo)
 ;overlay(plot(0),dumtopo)

 ;Cross-section Line
 plres = True
 plres@gsLineColor             = "green"
 plres@gsLineDashPattern        = 3
 plres@gsLineThicknessF         = 3.
 plres@tfPolyDrawOrder          = "PostDraw"
 dumpl = gsn_add_polyline(wks, plot(0),lon_plane, lat_plane, plres)

 plres@gsLineColor             = "red"
 dumxp = gsn_add_polyline(wks, plot(0),(/glon_idB,glon_p/), (/glat_idB,glat_p/), plres)

 txres               = True
 txres@txFont        = fncps_no
 txres@txFontHeightF = 0.02
 txres@txBackgroundFillColor = "white"
 dumtxt = gsn_add_text(wks,plot(0),"A",lon_plane(3),lat_plane(3),txres) ;
 dumtxt = gsn_add_text(wks,plot(0),"A'",lon_plane(dimsizes(lon_plane)-1),lat_plane(dimsizes(lat_plane)-1),txres) ;
 txres@txBackgroundFillColor = -1
 txres@txFontThicknessF      = 4.
 txres@txFontHeightF           = 0.025
 txres@txFontColor           = "black"
 dumtxt = gsn_add_text(wks,plot(0),"x",glon_id(0),glat_id(0),txres)

 pmres                 = True
 pmres@gsMarkerSizeF   = fncps_no
 pmres@gsMarkerThicknessF  = 3.
 pmres@gsMarkerColor   = "blue"
 pmres@gsMarkerIndex   = 7
 pmres@tfPolyDrawOrder  = "PostDraw"
 gsn_polymarker(wks,plot(0),lon_plane(sloc),lat_plane(sloc),pmres) ;

;Add annotation
 txres@txFontHeightF = 0.026
 txres@txBackgroundFillColor = "white"
 gsn_text_ndc(wks,"a)",0.07,0.95,txres) ;

;-----------------------------------------------------------------------------------------
;#11b  Cross-Section Plot Resources
;-----------------------------------------------------------------------------------------
 ;Viewport
 rescs@vpXF                          = 0.14
 rescs@vpYF                          = 0.48 ;0.45
 rescs@vpHeightF                     = 0.38   ;0.35
 rescs@vpWidthF                      = 0.33;   0.4

 rescsv = rescs

 rescs@gsnLeftString                 = " " ;  "c-s AA'"
 rescs@gsnRightString                = hydrometeors(ixp) + " ["+hydunits+"]"
 rescs@tiXAxisString                 = "Longitude" ;"Latitude"
 rescs@tiYAxisString                 = "Height (km)"

 rescs@trYMaxF                    = hloc ;80
 rescs@trXMaxF = glat_id + 0.4 ;0.3 ;30. ;xlon(1) ;51.1
 rescs@trXMinF = glat_id - 0.1 ;10. ;xlon(0) ;50.7

 loc_label                        = ispan(11,99,10)
 rescs@tmYLMode                   = "Explicit"
 rescs@tmYLValues                 = loc_label
 rescs@tmYLLabels                 = round(hgt_plane(loc_label,sloc),0)

 ;Label
 rescs@lbLabelAutoStride          = True
 rescs@lbOrientation              = "vertical"
 ;rescs@lbTitleString              = hydunits 
 rescs@lbLabelFontHeightF         = 0.015
 rescs@pmLabelBarWidthF           = 0.05
 rescs@pmLabelBarHeightF          = rescs@vpHeightF ;0.35


 ;Contour
 ;rescs@cnMissingValFillPattern       = 8
 ;rescs@cnMissingValFillColor         = "black"
 ;rescs@cnMissingValPerimOn           = True            ;
 ;rescs@cnMissingValPerimColor        = "black"         ;
 ;rescs@cnMissingValPerimDashPattern  = 1          ;
 ;rescs@cnMissingValPerimThicknessF   = 2.0
 rescs@cnLevelSelectionMode          = "ExplicitLevels"

; Cross-section contour
 rescs@cnLevels                      = qcnlevs
 rescs@cnFillOn                   = True
 rescs@cnFillMode                 = "RasterFill"
 rescs@cnLinesOn                  = False
 rescs@cnInfoLabelOn              = False   ; Turn off annotations
 rescs@cnLineLabelsOn             = False

 ;rescs@trXMinF       = min(lat_c)
 ;rescs@trXMaxF       = max(lat_c)
 rescs@tiXAxisString = "Latitude"
 plot(1)                          = gsn_csm_contour(wks, varh_plane(ixp,:,:), rescs)
 ;plot(1)                          = gsn_csm_contour(wks, uz_plane, rescs)

 if (ixp.eq.4) then
   rescs@cnFillOn = False
   rescs@cnLinesOn                  = True 
   rescs@cnLineThicknessF           = 2.0
   rescs@cnLineLabelsOn             = True
   dumhail = gsn_csm_contour(wks, varh_plane(5,:,:), rescs)
   overlay(plot(1),dumhail)
 end if

 rescs@cnLineThicknessF           = 2.0
 rescs@cnFillOn                   = False
 rescs@cnLineColor                = hydro_colors(0) 
 ;overlay(plot(1),dumqp)

 rescs@cnLineColor                = hydro_colors(2)
 ;overlay(plot(1),dumqp)

 rescs@cnLineColor                = hydro_colors(3) 
 ;overlay(plot(1),dumqp)

 rescs@cnLineColor                = hydro_colors(4) 
 ;overlay(plot(1),dumqp)

 rescs@cnLineColor                = hydro_colors(5) 
 ;overlay(plot(1),dumqp)

;;New Addition Oct 26
 spol          = dim_sum_n(varh_plane,0)         ;
 spol_val      = where(spol.gt.0.01,1.,0)
 spol_val!0    = "height"
 spol_val!1    = "longitude"
 spol_val&longitude = varh_plane&longitude

 ;For Overlaying plots
 ;--------------------
 reso = rescs
 delete(reso@cnLevels)
 reso@cnLevelSelectionMode       = "ExplicitLevels"
 reso@cnLinesOn                  = True
 reso@cnFillOn  = False
 reso@cnLineLabelsOn = False
 reso@cnInfoLabelOn = False
 reso@cnLineColor   = "grey30"
 reso@cnLevels = (/0.9,1.0,1.1/)
 reso@gsnLeftString =""
 reso@gsnRightString=""
 dumo     = gsn_csm_contour(wks,spol_val,reso)
 overlay(plot(1),dumo)

 ;;Vector overlay along the cross-section
 rescsv@vcRefAnnoFontHeightF     = 0.02
 rescsv@vcRefMagnitudeF          = 30.            ; define vector ref mag
 rescsv@vcRefLengthF             = 0.045            ; define length of vec ref
 rescsv@vcGlyphStyle             = "CurlyVector"    ; turn on curly vectors
 rescsv@vcMinDistanceF           = 0.019            ; thin vectors
 rescsv@vcRefAnnoOrthogonalPosF  = -1.11 ;-0.16 ;0.46               ; move ref vector down
 rescsv@vcLineArrowColor         = "azure3" ;"grey40"
 rescsv@vcLineArrowThicknessF    = 2.
 rescsv@vcMapDirection           = False

 rescsv@gsnLeftString            = ""
 rescsv@gsnRightString           = ""

 ;Make a copy for multiple plots
 ux_plane_c                      = ux_plane
 uz_plane_c                      = uz_plane

 if (uzmask) then
 ;Overlay horizontal vector for non hydrometeor regions
   ux_plane                        = mask(ux_plane,varh_plane(nhyd,:,:).ge.0.1,False)
   uz_plane                        = 0.
   rescsv@vcRefAnnoParallelPosF    = 0.1
   rescsv@vcRefAnnoString2         = "UX"
   dumvector                       = gsn_csm_vector(wks, ux_plane,uz_plane,rescsv)
   overlay(plot(1),dumvector)

 ;Overlay vertical vectors for hydromteor regions
   uz_plane_c                      = mask(uz_plane_c,varh_plane(nhyd,:,:).ge.0.1,True)
   ux_plane_c                      = 0.
   rescsv@vcRefMagnitudeF          = 10.
   rescsv@vcRefAnnoParallelPosF    = 0.2
   rescsv@vcRefAnnoString2         = "UZ"
   dumvector                       = gsn_csm_vector(wks, ux_plane_c,uz_plane_c,rescsv)
   overlay(plot(1),dumvector)
 else
   rescsv@vcRefAnnoString2         = "m/s"
   dumvector                       = gsn_csm_vector(wks, ux_plane,uz_plane,rescsv)
   overlay(plot(1),dumvector)
 end if
 ;;Temperature overlay
 rescs@cnFillOn                   = False
 rescs@cnLinesOn                  = True
 rescs@cnLineLabelsOn             = True
 rescs@cnLineLabelFontHeightF     = 0.015
 rescs@cnLineColor                = "black"
 rescs@cnLineThicknessF           = 3.
 delete(rescs@cnLevels)
 rescs@cnLevels                   = (/0./)
 dumtc                           = gsn_csm_contour(wks, tc_plane, rescs)
 overlay(plot(1),dumtc)

 delete(rescs@cnLevels)
 rescs@cnLevels                  = (/-4.,-3.0,-2.0,-1.0/)  
 dumte                            = gsn_csm_contour(wks, te_plane, rescs)
 overlay(plot(1),dumte)

 pmres@gsMarkerSizeF   = 0.045
 dumpm = gsn_add_polymarker(wks,plot(1),x_plane(sloc),0.,pmres) ;

 ;Add annotation
 txres@txFontHeightF = 0.026
 txres@txBackgroundFillColor = "white"
 ;gsn_text_ndc(wks,"c)",0.07,0.5,txres) ;

;-----------------------------------------------------------------------------------------
;#11c  Sounding Plot Resources
;-----------------------------------------------------------------------------------------
 ;Viewport
 resV@vpXF                          = 0.75
 resV@vpYF                          = 0.94
 resV@vpHeightF                     = 0.35
 resV@vpWidthF                      = 0.22

 resV@gsnLeftString                 = "Vertical Profile"
 resV@gsnRightString                = ""

 resV@xyLineThicknessF              = 3.

 resV@tiYAxisString                 = "Height (km)"
 resV@tiXAxisString                 = hydunits 

; resV@tmXBMode                      = "Explicit"
; resV@tmXBValues                    = ispan(0,9,1)*1.0
; resV@tmXBLabels                    = ispan(0,9,1)
 resV@tmLabelAutoStride             = True
 resV@tmXBAutoPrecision             = False
 resV@trXMinF                       = 0. ;min(qcnlevs) ;qmin ;0.0
 resV@trXMaxF                       = qprom(prom-1) ;max(qcnlevs) ;qmax 
 if (opt.eq.0 .or. opt.eq.2) then
   resV@trXLog                        = True
 end if
 resV@trYMaxF                       = max(ht_v)
 resV@xyLineColor                   = hydro_colors(ixp)
 plot(2)                            = gsn_csm_xy (wks,varh_v(ixp,:),ht_v,resV) ; create plot

;Add legend
 lgres                    = True
 lgres@lgLineColors       = hydro_colors
 lgres@lgLineThicknessF   = resV@xyLineThicknessF
 lgres@lgDashIndexes      = hydro_dash
 lgres@lgItemType         = "Lines"
 lgres@lgLabelFontHeightF = .1
 lgres@vpWidthF           = 0.1
 lgres@vpHeightF          = 0.18
 lgres@lgPerimOn          = False
 lbid = gsn_create_legend(wks,6,hydrometeors,lgres)

 amres = True
 amres@amParallelPosF   = 0.2
 amres@amOrthogonalPosF = -0.15
 annoid1 = gsn_add_annotation(plot(2),lbid,amres)

 ;Add annotation
 txres@txFontHeightF = 0.026
 txres@txBackgroundFillColor = -1 ;"white"
 gsn_text_ndc(wks,"b)",0.67,0.95,txres) ;

;Bonn Radar Extent
 ncirc    = 100
 circ_lat = new(ncirc,float)               ; Create arrays to hold circle.
 circ_lon = new(ncirc,float)
 raddeg   = extB ;
 nggcog(glat_idB,glon_idB,raddeg,circ_lat,circ_lon)

 lnres                  = True
 lnres@gsLineColor      = "black"
 lnres@gsLineThicknessF = 3.0
 lnres@gsLineDashPattern    = 1
 line_id1 = gsn_add_polyline(wks, plot(0), circ_lon, circ_lat, lnres)

 raddeg   = extB/2. ;
 nggcog(glat_idB,glon_idB,raddeg,circ_lat,circ_lon)

 line_id2 = gsn_add_polyline(wks, plot(0), circ_lon, circ_lat, lnres)

 raddeg   = extB/4. ;
 nggcog(glat_idB,glon_idB,raddeg,circ_lat,circ_lon)

 line_id3 = gsn_add_polyline(wks, plot(0), circ_lon, circ_lat, lnres)

; Draw cross hair
 dumch = new(2,"graphic")
 ply   = (/glat_idB ,glat_idB/)
 plx   = (/min(lon), max(lon)/)
 dumch(0) = gsn_add_polyline(wks, plot(0), plx, ply, lnres)

 plx = (/glon_idB ,glon_idB/)
 ply = (/min(lat), max(lat)/)
 dumch(1) = gsn_add_polyline(wks, plot(0), plx, ply, lnres)

;Draw panel plots
; psres = True
; maximize_output(wks,psres)

 draw(plot)
 frame(wks)


;---------------------------------------------------------------------------------------
;#13 Delete variables
;---------------------------------------------------------------------------------------
 delete([/plot,res, rescs, resw, rescsv, varh, varh_plane/])
 delete([/zmodelc, lon_plane,lat_plane,hgt_plane/])
 delete([/lat_c, lon_c, topo_c, prec_c, hail_c/])
 delete([/vx_c, tlc_c, te_c, th_c, tk_c, qv_c, ps_c, ux_c, uy_c, uz_hhl_c, varh_c/])
 delete([/tlc,te_p, qv_p,ux_p, uy_p, vx_p, wx_p1, wx_p2/])
 delete([/x_plane,te_plane,tc_plane,ux_plane,uy_plane,uz_plane/])
 delete([/vr_sum, ht_v/])
 delete([/ux_v, uy_v, uz_v/])
 delete([/ux_plane_c, uz_plane_c/])

 end do      ;BIG LOOP

;---------------------------------------------------------------------------------------
;#14 Data Output 
;---------------------------------------------------------------------------------------
 dataout = False
 if (dataout) then

 fnamenc = itype_gscp + "_" + inst + ".nc"
 system("rm " + fnamenc)

 fout = addfile(fnamenc,"c")
 fout->Band1 = Band1
 fout->Band2 = Band2
 fout->Band3 = Band3
 fout->Band4 = Band4
 fout->lat   = latout
 fout->lon   = lonout
 fout->Band5 = Band5
 fout->Band6 = (/Band6/)
 system("mkdir PNG")
 system("mv *.png PNG")
 system("mv PNG PNG_" + itype_gscp + "_" + inst)

 end if 
end
