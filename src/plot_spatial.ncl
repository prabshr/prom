;Plot the plan veiw of hydrometeors at different heights or pressure levels 

;P. Shrestha
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$prom_root/functions/shade_terrain.ncl"
load "$prom_root/functions/cosmo_time_nclf.ncl"
load "$prom_root/functions/calcrho_cos_nclf.ncl"
load "$prom_root/functions/convol.ncl"

begin
;STEP 1 User Settings
;-----------------------------------------------------------------------------------------
 ;opt       = 2   ;0-NX, 1-QX , 2-DX
 ;qx = 4 
 iz          =47; 34 ;21(low level) ; 47(mid level) , 41? ;model level to plot
 fpress = 500.   ;hPa

 pname="EMPA2HC_PROM1";
 inst=6
 lnd = "HET"
 itype = "ART"
 itype_gscp  = lnd + itype

 diri    = "/daten/daten04/pshrestha/shared-data/bonnRadar_20150705"
 diri    = "/p/largedata/hbn33/shrestha1/"

 sfold       = "runE_" + inst + "_" + itype_gscp
 expname     = lnd + "_" + pname
 subevent    = expname + "/" + itype_gscp + "/" + sfold

 diri_cos    = diri + "/" + "/" + subevent + "/cosout/"
 diri_cosI   = diri + "/" + "/" + subevent + "/cosout/ivr/"

;STEP 2 Set variables to read and specify scales or offsets
;-------------------------------------------------------------------------------------
 hydrometeors_qx = (/"QC","QR","QI","QS","QG","QH"/)
 hydrometeors_qn = (/"NCCLOUD","NCRAIN","NCICE","NCSNOW","NCGRAUPEL","NCHAIL"/)

 if (opt.eq.0) then
   hydrometeors = hydrometeors_qn
   if (qx.eq.0) then
     scale_varh   = 1/1.E6
     qcnlevs      = ispan(100,5000,500) ;(/1.1, 1.E1, 1.E2, 1.E3, 1.E4, 1.E5, 1.E6, 1.E7/)
     hydunits     = "cm~S~-3~N~"
   else
     scale_varh   = 1/1.E3
     qcnlevs      = ispan(100,1000,100); (/1.1, 1.E1, 1.E2, 1.E3, 1.E4, 1.E5, 1.E6, 1.E7/)
     hydunits     = "L~S~-1~N~"
   end if
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
   if (qx.eq.0) then
     qcnlevs   = ispan(0,50,5)   ;um
     scale_varh   = 1.E6   ; m to  um
     hydunits     = "um"
   else
     qcnlevs      = ispan(100,2000,200)/1000.  ; (/0.001, 0.1, 1.0, 5.0, 10., 15., 20., 25., 30./)
     scale_varh   = 1000.   ; m to  mm
     hydunits     = "mm"
   end if
   hydname      = "Max. Volume Diameter"
 else
   print("code not written for opt = " + opt)
   exit
 end if
 end if
 end if
;

;STEP 3  READ NC FILES
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

 print("Processing " + fils_cos(fns))

 nc_cos    = addfile(fils_cos(fns),"r")
 time      = nc_cos->time
 ;
 str_time  = cosmo_time_nclf(time,1)  ; 0 for date
 print("Processing ..." + fns + "  " + str_time)

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

;STEP 6 Variable to plot
;------------------------------------------------------------------------------
 x1 = 120
 y1 = 120
 x2 = nx-1
 y2 = ny-1
 ;spol          = wrf_user_intrp3d(varh(qx,:,y1:y2,x1:x2),pressure(:,y1:y2,x1:x2),"h",fpress,0.,False)
 ;sux           = wrf_user_intrp3d(ux(:,y1:y2,x1:x2),pressure(:,y1:y2,x1:x2),"h",fpress,0.,False) 
 ;suy           = wrf_user_intrp3d(uy(:,y1:y2,x1:x2),pressure(:,y1:y2,x1:x2),"h",fpress,0.,False)

 spol    = varh(qx,iz,y1:y2,x1:x2)
 sux     = ux(iz,y1:y2,x1:x2)
 suy     = uy(iz,y1:y2,x1:x2)

;STEP 7 Plot
;------------------------------------------------------------------------------
 print("Plotting ...")

 fncps_no = ":F21:"
 fncps_no = 21
 fncps_hgt = 0.035

 res      = True
 wks      = gsn_open_wks("x11", "Figure_QX_"+opt)

 cmap = read_colormap_file("cmp_flux")
 cmap(0,:) = (/1., 1., 1., 1./)

;;gsn Resources
 ;---------------
 res@gsnDraw    = False
 res@gsnFrame   = False
 res@gsnLeftStringFontHeightF   = fncps_hgt
 res@gsnRightStringFontHeightF  = fncps_hgt
 res@gsnMaximize                = False
 res@gsnLeftString = hydrometeors(qx) + " " + hydunits
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
; res@tiXAxisString = "Latitude [~F34~0~F~N]"
; res@tiYAxisString = "Height [km]"

 resu  = res

 ;cn resources
 ;--------------------
 res@cnFillOn   = True
 res@cnFillPalette              = cmap
 res@cnFillMode = "CellFill"
 res@cnLinesOn  = False
 res@cnLevelSelectionMode       = "ExplicitLevels"
 res@cnLevels   = qcnlevs

 res@lbOrientation = "Vertical"
 ;res@vpWidthF = 0.4

 plot = gsn_csm_contour(wks,spol,res)

 ;-----------------------
 ;;Vector overlay 
 resu@tmYLLabelsOn  = False
 resu@tmXBLabelsOn  = False

 resu@tiYAxisString = ""
 resu@vcRefAnnoFontHeightF     = 0.01
 resu@vcRefMagnitudeF          = 30.            ; define vector ref mag
 resu@vcRefLengthF             = 0.045            ; define length of vec ref
 resu@vcGlyphStyle             = "CurlyVector"    ; turn on curly vectors
 resu@vcMinDistanceF           = 0.020            ; thin vectors
 resu@vcRefAnnoOrthogonalPosF  = -0.16 ;0.46               ; move ref vector down
 resu@vcLineArrowColor         = "black" ;"azure3" ;"grey40"
 resu@vcLineArrowThicknessF    = 2.
 resu@vcMapDirection           = False
 dumu     = gsn_csm_vector(wks,sux,suy,resu)
 overlay(plot,dumu)


 xval = (/xx,xx/)
 yval = (/0, 300/)
 gres = True
 duml = gsn_add_polyline(wks,plot,xval,yval,gres)
 draw(plot)
 frame(wks)

end
