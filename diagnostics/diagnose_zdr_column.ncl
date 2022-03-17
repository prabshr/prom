;ZDR column usually appears 15-20 minutes before hail reaches ground
;First examine the ground hail, then identify time when the ZDR column should appear
;Exporle mean rain diameter at 6 km height, it should show the location of ZDR column
;Now use this to explore model hydrometeor, aerosols and polarimetric radar data using
;c-s and vertical profile !
; Variables to plot are spol_** and cs_** at given height (zheight) and X-extent (XR)

;P. Shrestha
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$prom_root/functions/shade_terrain.ncl"
load "$prom_root/functions/cosmo_time_nclf.ncl"
load "$prom_root/functions/calcrho_cos_nclf.ncl"
load "$prom_root/functions/convol.ncl"
load "$prom_root/functions/generateXY.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$prom_root/functions/read_emvorado_pol.ncl"
load "$prom_root/functions/lndistART.ncl"
load "$PROM/functions/load_radarColors.ncl"

begin

;STEP 1 User Settings
;-----------------------------------------------------------------------------------------
 pname="EMPA2_PROM1";
 inst=6
 lnd = "HET"
 itype = "ART"
 itype_gscp  = lnd + itype

 ;diri   = "/daten/daten04/pshrestha/shared-data/bonnRadar_20150705"
 diri    = "/p/largedata/hbn33/shrestha1/"

 sfold       = "runE_" + inst + "_" + itype_gscp
 expname     = lnd + "_" + pname
 subevent    = expname + "/" + itype_gscp + "/" + sfold

 diri_cos    = diri + "/" + "/" + subevent + "/cosout/"
 diri_acos   = diri + "/" + "/" + subevent + "/a_cosout/"
 diri_cosI   = diri + "/" + "/" + subevent + "/cosout/ivr/"

; itype_gscp = "HETART_NH"
; sfold       = "runE_" + inst + "_" + itype_gscp
; expname     = lnd + "_" + pname
; subevent    = expname + "/" + itype_gscp + "/" + sfold
 diri_emvo   = "/p/largedata/hbn33/test-emvo/" + subevent

 ;Bonn Radar
 glon_idB    = 7.071663 ;
 glat_idB    = 50.73052 ;

 ;Time , Extent and C-S and variables to plot
 fns         = 12  ; For hydrometeor 5 minute and polarimetric data
 afn         = 33  ; For aerosol 1 hourly data
 zheight     = 6.0 ; km
 XR          = 10; 41 ; km
 YR          = 7;
 fcsx        = True ;  X or Y C-S
 fuv         = False ; Overlay wind vectors
 xycorner    = (/-5, -5, 40, 40/)

 readMODELHYDRO = True
   hopt         = 1   ;0-NX, 1-QX , 2-DX
   iqx          = 5 
 readEMVORADO   = True 
   popt         = True
   popt@profile = -1   ;return 3d data
   popt@nf      =  fns
   ipol         = 1 
 readART        = True 
   iaero        = 0 
   aopt         = True
   aopt@iy      = 90
   aopt@ix      = 150
   aopt@favg    = True
   aopt@fout    = 1  ;1-Ni(cm-3), 2-Mi(m3/m3), 3-Dp (nm)2
   aopt@nbdl     = 0

;STEP 2a Polarimetric Setup
;-------------------------------------------------------------------------------------
 pvars       = (/"Z~B~H~N~","Z~B~DR~N~","K~B~DP~N~","~F8~r~B~~F10~hv~E~"/)
 punits      = (/"dBZ","dB","deg/km", "-"/)
 plevels     = (/"zhh_levs","zdr_levs","kdp_levs","rho_levs"/)
 tmpcolors     = load_radarColors(0)
 pcolors     = tempcolors(2:)
 pcnlevs     = tmpcolors@$(plevels(ipol))$

;STEP 2b Aerosol Setup
;-----------------------------------------------------------------------------------------
; Geomeric standard deviations of lognormal modes
 sig = (/1.7, 2.0, 1.7, 2.0, 1.4, 2.5, 1.7, 1.6, 1.5, 1.9, 2.0, 1.7/)

;Obtain aggregate of 12 log-normal distribution modes
 avarname  = (/"Nx", "Mx", "Dp"/)
 aunits    = (/"cm~S~-3~N~", "m~S~3~N~/m~S~3~N~","nm"/)
 smodes = ispan(1,12,1)
 anames = (/"Nuc-Acc", "Soot-Coars", "Dust", "Sea"/)
 acnlevs= (/100,300,500,700,1000,1300,1500,1700,2000,2300,2500,2700,3000/)
 acolors = read_colormap_file("cmp_flux")
 acolors(0,:) = (/1., 1., 1., 1./)


;STEP 2b Hydrometeor Setup
;-------------------------------------------------------------------------------------
 hydrometeors_qx = (/"QC","QR","QI","QS","QG","QH"/)
 hydrometeors_qn = (/"NCCLOUD","NCRAIN","NCICE","NCSNOW","NCGRAUPEL","NCHAIL"/)

 qcolors = read_colormap_file("cmp_flux")
 qcolors(0,:) = (/1., 1., 1., 1./)

 if (hopt.eq.0) then
   hydrometeors = hydrometeors_qn
   if (iqx.eq.0) then
     scale_varh   = 1/1.E6
     qcnlevs      = ispan(0,5000,500) ;(/1.1, 1.E1, 1.E2, 1.E3, 1.E4, 1.E5, 1.E6, 1.E7/)
     qcnlevs(0)   = 100
     hydunits     = "cm~S~-3~N~"
   else
     scale_varh   = 1/1.E3
     qcnlevs      = (/1.1, 1.E1, 1.E2, 1.E3, 1.E4, 1.E5, 1.E6, 1.E7/)
     hydunits     = "L~S~-1~N~"
   end if
   hydname      = "Number Concentration"
 else if (hopt.eq.1) then
   hydrometeors = hydrometeors_qx
   scale_varh   = 1000.  ;kg/kg to g/kg
   qcnlevs      = (/0.01,0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0/)
   hydunits     = "g/m~S~3~N~"
   hydname      = "Sp. Hydrometeor Content"
 else if (hopt.eq.2) then
   x3           = 1./3.
   hydrometeors = hydrometeors_qx
   a_geo        = (/0.124    , 0.124  , 0.835 , 2.4    , 0.142     , 0.1366/)  ;mkg-b
   b_geo        = (/x3       , x3     , 0.39  , 0.455  , 0.314     , x3/)      ;-
   xmin         = (/4.2E-15  , 2.6E-10, 1.E-12, 1.E-10, 1.E-9, 2.6E-9/)        ;kg
   xmax         = (/2.6E-10  , 3.0E-6 , 1.E-6 , 2.E-5 , 5.E-4, 5.E-4/)         ;g
   qmax         = 30.
   qmin         = 0.001
   if (iqx.eq.0) then
     qcnlevs   = ispan(5,50,5)   ;um
     scale_varh   = 1.E6   ; m to  um
     hydunits     = "um"
   else
     qcnlevs      = (/100,300,500,1000,3000,5000,7000,10000,13000/)/1000.  ;m ;(/0.001, 0.1, 1.0, 5.0, 10., 15., 20., 25., 30./)
     scale_varh   = 1.E3   ; m to  mm
     hydunits     = "mm"
   end if
   hydname      = "Max. Volume Diameter"
 else
   print("code not written for hopt = " + hopt)
   exit
 end if
 end if
 end if
;

;STEP 3  READ COSMO NC FILES
;-------------------------------------------------------------------------------------
 nc_cosI   = addfile(diri_cosI +"lfff00000000c.nc","r")
 lat       = nc_cosI->lat
 lon       = nc_cosI->lon
 topo      = nc_cosI->HSURF(0,:,:)

 hhl       = nc_cosI->HHL(0,::-1,:,:)
 ndimhhl   = dimsizes(hhl)
 nz        = ndimhhl(0)-1
 ny        = ndimhhl(1)
 nx        = ndimhhl(2)
 zmodel    = 0.5*(hhl(0:nz-1,:,:) + hhl(1:nz,:,:))/1000.  ;km
;
 fils_cos  = systemfunc("ls " + diri_cos + "lf*")

 fname     = fils_cos(fns)
 print("Processing " + fname)

 nc_cos    = addfile(fname,"r")
 time      = nc_cos->time
 ;
 str_time  = cosmo_time_nclf(time,1)  ; 0 for date
 print("Processing ..." + fns + "  " + str_time)
 ;
 XYkm = generateXY(lat,lon,glat_idB,glon_idB)
 Ykm  = XYkm(1,:,:)
 Xkm  = XYkm(0,:,:)

 ixs = closest_val(XR,Xkm(150,:))
 iys = closest_val(YR,Ykm(:,150))
 print("Cross-section through " + XR + " km " + ixs)
 print("Cross-section through " + YR + " km " + iys)

;STEP 4 Prognostic state variables
;------------------------------------------------------------------------------
 tair            = nc_cos->T(0,::-1,:,:)     ;K
 pressure        = nc_cos->P(0,::-1,:,:)     ;Pa
 qv              = nc_cos->QV(0,::-1,:,:)    ;kg/kg
 ux              = nc_cos->U(0,::-1,:,:)
 uy              = nc_cos->V(0,::-1,:,:)
 uz              = nc_cos->W(0,::-1,:,:)
 uz_hhl          = uz(1:80,:,:)*0.
 do kk = 0,79
   uz_hhl(kk,:,:)  = 0.5*(uz(kk,:,:) + uz(kk+1,:,:))
 end do

 spol_uz          = wrf_user_intrp3d(uz_hhl,zmodel,"h",zheight,0.,False)
 spol_ux          = wrf_user_intrp3d(ux,zmodel,"h",zheight,0.,False)
 spol_uy          = wrf_user_intrp3d(uy,zmodel,"h",zheight,0.,False)

 spol_uz!0        = "Y"
 spol_uz!1        = "X"
 spol_uz&Y        = Ykm(:,150) ;
 spol_uz&X        = Xkm(150,:) ;

 copy_VarCoords(spol_uz,spol_ux)
 copy_VarCoords(spol_uz,spol_uy)

 if (fcsx) then
   cs_uz       = uz_hhl(:,:,ixs)
   sdims       = dimsizes(cs_uz)
   cs_uz!0     = "Z"
   cs_uz!1     = "Y"
   cs_uz&Y     = Ykm(:,ixs) ;
 else
   cs_uz       = uz_hhl(:,iys,:)
   sdims       = dimsizes(cs_uz)
   cs_uz!0     = "Z"
   cs_uz!1     = "X"
   cs_uz&X     = Xkm(iys,:) ;
 end if

 if (readMODELHYDRO) then
;STEP 5 Extract Hydrometeors
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
   if (hopt.eq.0) then
     varh(ih,:,:,:) = (/nc_cos->$(hydrometeors_qn(ih))$(0,::-1,:,:)/)
     varh(ih,:,:,:) = (/scale_varh * varh(ih,:,:,:) * rho/)    ;(kg-1) to (L-1)
           ;Nt    = 1/1000. * Nd * rho            ; number density to number concentration
   else if (hopt.eq.1) then
     varh(ih,:,:,:) = (/nc_cos->$(hydrometeors_qx(ih))$(0,::-1,:,:)/)
     varh(ih,:,:,:) = (/scale_varh * varh(ih,:,:,:) * rho/)    ;(kg/kg) to (g/m3)
           ;Qt    = 1000. * Qd * rho
   else if (hopt.eq.2) then
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

 spol_qx          = wrf_user_intrp3d(varh(iqx,:,:,:),zmodel,"h",zheight,0.,False)
 copy_VarCoords(spol_uz,spol_qx)

 ;Edge for C-S
 if (fcsx) then
   cs_spol          = dim_sum_n(varh(:,:,:,ixs),0)         ;
   cs_spol_val      = where(cs_spol.gt.0.01,1.,0)
   cs_spol_edge     = convol(cs_spol_val)

   xval2x        = conform_dims(sdims,Ykm(:,ixs),1)*1.
   zval2x        = zmodel(:,:,ixs)

   cs_tairC      = (/tair(:,:,ixs) - 273.16/)
   cs_qx         = varh(iqx,:,:,ixs)
 else
   cs_spol          = dim_sum_n(varh(:,:,iys,:),0)         ;
   cs_spol_val      = where(cs_spol.gt.0.01,1.,0)
   cs_spol_edge     = convol(cs_spol_val)

   xval2x        = conform_dims(sdims,Xkm(iys,:),1)*1.
   zval2x        = zmodel(:,iys,:)

   cs_tairC      = (/tair(:,iys,:) - 273.16/)
   cs_qx         = varh(iqx,:,iys,:)
 end if

 copy_VarCoords(cs_uz,cs_spol_val)
 copy_VarCoords(cs_uz,cs_tairC) 
 copy_VarCoords(cs_uz,cs_qx)  
 
; for terrain shading
 xypoly  = shade_terrain(xval2x,zval2x)
 xpoly   = xypoly(0,:)
 ypoly   = xypoly(1,:)

 end if   ;readMODELHYDRO


 if (readEMVORADO) then
;STEP 7 Read EMVORADO Outputs 
;------------------------------------------------------------------------------
 nh           = dimsizes(pvars)
 varp         = new((/nh,nz,ny,nx/),"float")    ;Initialize variable for plotting
 varp@_FillValue = -999.99

 popt@var  = "zh"
 Zh = read_emvorado_pol(diri_emvo,popt)
 popt@var  = "zv"
 Zv = read_emvorado_pol(diri_emvo,popt)
 popt@var  = "kdp"
 KDP = read_emvorado_pol(diri_emvo,popt)
 popt@var  = "rrhv"
 rRHV = read_emvorado_pol(diri_emvo,popt)
 popt@var  = "irhv"
 iRHV = read_emvorado_pol(diri_emvo,popt)

 RHV = ((rRHV^2+iRHV^2)^0.5) / ((Zh*Zv)^0.5)   ;[-]
 KDP = KDP*1000.                               ;[deg/km]
 Zh  = 10*log10(Zh/1.)                         ;[dBZ]
 Zv  = 10*log10(Zv/1.)                         ;[dBZ]

 varp(0,:,:,:) = where(Zh.lt.-10., varp@_FillValue, Zh)
 varp(1,:,:,:) = where(Zh.lt.-5., varp@_FillValue, Zh-Zv)
 varp(2,:,:,:) = where(Zh.lt.-5., varp@_FillValue, KDP)
 varp(3,:,:,:) = where(Zh.lt.-5., varp@_FillValue, RHV)
 delete([/Zh, Zv, KDP, RHV/])


 spol_pol      = wrf_user_intrp3d(varp(ipol,:,:,:),zmodel,"h",zheight,0.,False)
 spol_pol@_FillValue = varp@_FillValue
 copy_VarCoords(spol_uz,spol_pol)

 if (fcsx) then
   cs_pol    = (/varp(ipol,:,:,ixs)/)
 else
   cs_pol    = (/varp(ipol,:,iys,:)/)
 end if
 copy_VarCoords(cs_uz,cs_pol)

 end if  ;readEMVORADO

 if (readART) then

;STEP 7  Extract lognormal ASD for each mode
;-----------------------------------------------------------------------------------------

 afils_cos  = systemfunc("ls " + diri_acos + "lf*")
 afil_cos   = afils_cos(afn)

 varo  = new((/12,nz,ny,nx/),"float")
 do im = 0, dimsizes(smodes)-1
   aopt@smode = smodes(im)
   varo(im,:,:,:) = lndistART(afil_cos,aopt)
 end do

 varao  = new((/4,nz,ny,nx/),"float") 
 ;Extract sums at specific height
 varao(0,:,:,:) = dim_sum_n(varo(0:3,:,:,:),0)     ;nuc-acc mode ;small
 varao(1,:,:,:) = dim_sum_n(varo(4:5,:,:,:),0)     ;soot+coarse
 varao(2,:,:,:) = dim_sum_n(varo(6:8,:,:,:),0)     ;dust
 varao(3,:,:,:) = dim_sum_n(varo(9:11,:,:,:),0)    ;sea

 varaz  = new((/4,ny,nx/),"float")
 do ia = 0,3
   varaz(ia,:,:) = wrf_user_intrp3d(varao(ia,:,:,:),zmodel,"h",zheight,0., False)
 end do

 spol_aero    = varaz(iaero,:,:) 
 spol_aero@long_name = anames(iaero)
 copy_VarCoords(spol_uz,spol_aero)

 if (fcsx) then
   cs_aero    = varao(iaero,:,:,ixs)
 else
   cs_aero    = varao(iaero,:,iys,:)
 end if

 cs_aero@long_name = anames(iaero)
 copy_VarCoords(cs_uz,cs_aero)

 end if ;READART

;STEP 7 Plot
;------------------------------------------------------------------------------
 print("Plotting ...")

 fncps_no = ":F21:"
 fncps_no = 21
 fncps_hgt = 0.035

 plot     = new(2,"graphic")

 res      = True
 figname  = "Figure_QX_rhi_" + pname + "_" + itype_gscp + "_t" + fns
 wks      = gsn_open_wks("x11", figname)

;;gsn Resources
 ;---------------
 res@gsnDraw    = False
 res@gsnFrame   = False
 res@gsnLeftStringFontHeightF   = fncps_hgt
 res@gsnRightStringFontHeightF  = fncps_hgt
 res@gsnMaximize                = False
 res@gsnLeftString = ""
 res@gsnRightString = ""

 ;tm Resources
 ;---------------
 res@tmXBLabelFont              = fncps_no
 res@tmYLLabelFont              = fncps_no
 res@tmYLLabelFontHeightF       = fncps_hgt
 res@tmXBLabelFontHeightF       = fncps_hgt
 res@tmYROn                     = False
 res@tmXTOn                     = False
 res@tmXBLabelFontColor         = "black"
 res@tmYLLabelFontColor         = "black"

 ;ti Resources
 ;-----------------
 res@tiXAxisFont                = fncps_no
 res@tiYAxisFont                = fncps_no
 res@tiXAxisFontHeightF         = fncps_hgt
 res@tiYAxisFontHeightF         = fncps_hgt

 res@tiXAxisString = "X Range [km]" 
 res@tiYAxisString = "Y Range [km]"

 cs_res = res
 resu   = res

 ;pm resources
 ;--------------------
 res@pmTickMarkDisplayMode      = "Always"
 ;res@pmLabelBarOrthogonalPosF   = 0.3
 ;res@pmLabelBarParallelPosF     = 0.2

 ;lb resources
 ;--------------------
 res@lbLabelBarOn  = False
 res@lbOrientation = "Vertical"

 ;tr Resources
 ;-----------------------------
 res@trXMinF = xycorner(1)
 res@trYMinF = xycorner(0)
 res@trXMaxF = xycorner(3)
 res@trYMaxF = xycorner(2)

 ;cn resources
 ;--------------------
 res@cnFillOn   = True
 res@cnLineLabelsOn = False
 res@cnInfoLabelOn = False
 res@cnFillMode = "CellFill"
 res@cnLinesOn  = False
 res@cnLevelSelectionMode       = "ExplicitLevels"
 res@cnLevels       = qcnlevs 
 res@cnFillPalette  = qcolors ;qcolors cmap
 res@cnFillOpacityF = 0.8
 plot(0) = gsn_csm_contour(wks,spol_qx,res)

 pres = res
 pres@cnLinesOn = True
 pres@cnFillOn = False
 pres@cnLineThicknessF = 3.0
 pres@cnLineDashPattern = 1

 delete(pres@cnLevels)
 uzcnlevs = (/-7.,-5.,-3.,-1.,5.,10.,15.,20.,25.,30.,35.,40./)
 pres@cnLevels = uzcnlevs
 pres@cnLineDrawOrder = "PreDraw"
 dumw = gsn_csm_contour(wks,spol_uz,pres)
 dumw = ColorNegDashZeroPosContour(dumw,"blue","gray","red")
 overlay(plot(0),dumw)

 if (fuv) then
 ;-----------------------
 ;;Vector overlay along the cross-section
  resu@tmYLLabelsOn  = False
  resu@tmXBLabelsOn  = False

  resu@tiYAxisString = ""
  resu@vcRefAnnoFontHeightF     = 0.01
  resu@vcRefMagnitudeF          = 30.            ; define vector ref mag
  resu@vcRefLengthF             = 0.045            ; define length of vec ref
  resu@vcGlyphStyle             = "CurlyVector"    ; turn on curly vectors
  resu@vcMinDistanceF           = 0.020            ; thin vectors
  resu@vcRefAnnoOrthogonalPosF  = -0.16 ;0.46               ; move ref vector down
  resu@vcLineArrowColor         = "grey20" ;"azure3" ;"grey40"
  resu@vcLineArrowThicknessF    = 2.
  dumu     = gsn_csm_vector(wks,spol_ux,spol_uy,resu)
  overlay(plot(0),dumu)
 end if

 ; Trick to linearize the plot
 cs_res@trGridType = "TriangularMesh"
 cs_res@sfXArray   = xval2x
 cs_res@sfYArray   = zval2x

 ;tr resources
 ;--------------------
 cs_res@trXMinF    = xycorner(0)
 cs_res@trXMaxF    = xycorner(2)
 cs_res@trYMaxF    = 15.
 ;cs_res@trYMinF    = 0.

 cs_res@vpWidthF   = 0.5

 ;cn resources
 ;-------------------- 
 cs_res@cnLineLabelsOn        = False
 cs_res@cnInfoLabelOn         = False
 cs_res@cnFillOn              = True
 cs_res@cnFillMode            = "CellFill"
 cs_res@cnLinesOn             = False
 cs_res@cnLevelSelectionMode  = "ExplicitLevels"
 cs_res@cnFillPalette         = qcolors ;cmap
 cs_res@cnLevels              = qcnlevs ;qcnlevs
 cs_res@lbOrientation         = "Vertical"
 cs_res@tiYAxisString         = "Height [km]"
 if (fcsx) then
   cs_res@tiXAxisString         = "Y Range [km]"
 else
   cs_res@tiXAxisString         = "X Range [km]"
 end if
 cs_res@cnFillOpacityF = 0.8
 plot(1) = gsn_csm_contour(wks,cs_qx,cs_res)

 cs_res@cnLineDrawOrder = "PreDraw"
 cs_res@cnLinesOn           = True
 cs_res@cnFillOn            = False
 cs_res@cnLineThicknessF    = 3.0
 cs_res@cnLineDashPattern   = 1
 delete(cs_res@cnLevels)
 cs_res@cnLevels = uzcnlevs

 dumw_cs = gsn_csm_contour(wks,cs_uz,cs_res)
 dumw_cs = ColorNegDashZeroPosContour(dumw_cs,"blue","gray","red")
 overlay(plot(1),dumw_cs)

 ;For storm edge and air temperature
 ;--------------------
 delete(cs_res@cnLevels)
 cs_res@cnLineColor           = "grey30"
 cs_res@cnLevels              = (/0.9,1.0,1.1/)
 cs_res@cnLineDrawOrder = "PostDraw"
 dum_otl     = gsn_csm_contour(wks,cs_spol_val,cs_res)
 delete(cs_res@cnLevels)
 cs_res@cnLevels              = (/0./)
 cs_res@cnLineDashPattern     = 0
 cs_res@cnLineColor           = "black"
 cs_res@cnLineThicknessF      = 3.
 dum_tc = gsn_csm_contour(wks,cs_tairC,cs_res)
 ;overlay(plot(1),dum_otl) ; only use when qopt = 1
 overlay(plot(1),dum_tc)

 ;C-S overlay on plan view
 gres                  = True
 gres@gsLineColor      = "black"
 gres@gsLineDashPattern = 0
 gres@gsLineThicknessF = 3.0

 if (fcsx) then
   ypoly_cs = (/xycorner(0),xycorner(2)/)
   xpoly_cs = (/Xkm(0,ixs),Xkm(0,ixs)/)
 else
   xpoly_cs = (/xycorner(1),xycorner(3)/)
   ypoly_cs = (/Ykm(iys,0),Ykm(iys,0)/)
 end if
 dum_cs = gsn_add_polyline(wks,plot(0),xpoly_cs,ypoly_cs,gres)

 ;Panel Plot
 resP = True
 resP@gsnMaximize    = True
 resP@gsnPanelXF = (/-1,0.63/)
 resP@gsnPanelFigureStrings = (/"a)","b)"/)  
 resP@gsnPanelFigureStringsFontHeightF = fncps_hgt/2.
; resP@gsnPanelDebug = True
 gsn_panel(wks,plot,(/1,2/),resP)

end
