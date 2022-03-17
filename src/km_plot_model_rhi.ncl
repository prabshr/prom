;Plot the c-s veiw of hydrometeors alond Y-axis
; specified model grid in X using linearized 

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

begin
;STEP 1 User Settings
;-----------------------------------------------------------------------------------------
 ;opt       = 2   ;0-NX, 1-QX , 2-DX
 iz          = 34 ;34 (2km a.g.l.) ;21(low level) ; 47(mid level) , 41? ;model level to plot

; fns         = 12  ;18
; XR          = 7; 41 ; km
; ixs         = 186 ;191 hail ZDR 
 xycorner    = (/-10, -10, 110, 110/)
 xycorner    = (/0, 0, 50, 50/)

 pname="EMPA2_PROM1";
 inst=6
 lnd = "HET"
 itype = "ART"
 itype_gscp  = lnd + itype

 diri    = "/daten/daten04/pshrestha/shared-data/bonnRadar_20150705"
 diri    = "/p/largedata/hbn33/shrestha1/"

 ;ydates    = (/"20150705","20160513","20170706"/)

 ;insts     = (/6, 8, 11/)
 ;fns       = (/11,34,18/) 
 ;ixs       = (/90,124,95/) + 65  ;Need to add offset due to clipping of radar

 ;;extract_ascii clips
 ;;Extent for clipping :58 237
 ;;Extent for clipping :65 244

 glon_idB    = 7.071663 ; Bonn Radar
 glat_idB    = 50.73052 ;

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
     scale_varh   = 1. ;1/1.E3
     qcnlevs      = ispan(100,1000,100)/100.; (/1.1, 1.E1, 1.E2, 1.E3, 1.E4, 1.E5, 1.E6, 1.E7/)
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
     qcnlevs   = ispan(5,50,5)   ;um
     scale_varh   = 1.E6   ; m to  um
     hydunits     = "um"
   else
     qcnlevs      = (/100,300,500,1000,3000,5000,7000,10000, 13000/)/1000.  ;m ;(/0.001, 0.1, 1.0, 5.0, 10., 15., 20., 25., 30./)
     scale_varh   = 1.E3   ; m to  mm
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
 ;
 XYkm = generateXY(lat,lon,glat_idB,glon_idB)
 Ykm  = XYkm(1,:,:)
 Xkm  = XYkm(0,:,:)

 ;
 ixs = closest_val(XR,Xkm(150,:))
 print("Cross-section through " + XR + " km " + ixs)

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
 spol          = dim_sum_n(varh(:,:,:,ixs),0)         ;
 spol_val      = where(spol.gt.0.01,1.,0)
 spol_edge     = convol(spol_val)

 sdims         = dimsizes(spol)
 spol!0        = "Z"
 spol!1        = "Y"
 spol&Y        = Ykm(:,ixs) ; 
 xval2x        = conform_dims(sdims,spol&Y,1)*1.
 zval2x        = zmodel(:,:,ixs)

 tairC         = tair(:,:,ixs) - 273.16
 uy_yp         = uy(:,:,ixs)
 uz_yp         = uz_hhl(:,:,ixs)

 copy_VarCoords(spol,tairC)   
 copy_VarCoords(spol, uy_yp)  
 copy_VarCoords(spol, uz_yp) 

 ;Need this for vector plot as stretching does not work
 uy_yp&Z      = zmodel(:,90,ixs)
 uz_yp&Z      = zmodel(:,90,ixs)
 
; for terrain shading
 xypoly  = shade_terrain(xval2x,zval2x)
 xpoly   = xypoly(0,:)
 ypoly   = xypoly(1,:)


;STEP 7 Plot
;------------------------------------------------------------------------------
 print("Plotting ...")

 fncps_no = ":F21:"
 fncps_no = 21
 fncps_hgt = 0.035

 plot     = new(4,"graphic")
 dum0     = new(4,"graphic")
 dum1     = new(4,"graphic")
 dumo     = new(4,"graphic")
 dumq     = new(4,"graphic")
 dumt     = new(4,"graphic")
 dumu     = new(4,"graphic")

 res      = True
 figname  = "Figure_QX_rhi_" + pname + "_" + itype_gscp + "_t" + fns + "_" +  opt 
 wks      = gsn_open_wks("x11", figname)

 cmap = read_colormap_file("cmp_flux")
 cmap(0,:) = (/1., 1., 1., 1./)

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
 ;res@tmXBLabelFontThicknessF    = 2.0
 ;res@tmXBMode   = "Explicit"
 ;res@tmXBValues = xval2x(0,ispan(0,sdims(1)-1,20))
 ;res@tmXBLabels = tointeger((res@tmXBValues)*100.)/100.
 ;res@tmXBLabelsOn  = False 

 ;ti Resources
 ;-----------------
 res@tiXAxisFont                = fncps_no
 res@tiYAxisFont                = fncps_no
 res@tiXAxisFontHeightF         = fncps_hgt
 res@tiYAxisFontHeightF         = fncps_hgt
 res@tiXAxisString = "Y Range [km]" 
 res@tiYAxisString = "Height [km]"

 pres = res

 ;pm resources
 ;--------------------
 res@pmTickMarkDisplayMode      = "Always"
 ;res@pmLabelBarOrthogonalPosF   = 0.3
 ;res@pmLabelBarParallelPosF     = 0.2

 ;
 res@lbOrientation = "Vertical"
 res@vpWidthF = 0.4

 ;tr resources
 ;--------------------
 res@trXMinF    = xycorner(0) 
 res@trXMaxF    = xycorner(2)
 res@trYMaxF    = 15.

 ;For vector plots
 resu  = res

 ; Trick to linearize the plot
 res@trGridType = "TriangularMesh"
 res@sfXArray   = xval2x
 res@sfYArray   = zval2x

 ;For other plots
 reso = res

 ;cn resources
 ;--------------------
 res@cnFillOn   = True
 res@cnFillPalette              = cmap
 res@cnFillMode = "CellFill"
 res@cnLinesOn  = False
 res@cnLevelSelectionMode       = "ExplicitLevels"
 res@cnLevels   = qcnlevs

 ;For Overlaying storm outline 
 ;--------------------
 reso@cnLevelSelectionMode       = "ExplicitLevels"
 reso@gsnLeftString = ""
 reso@gsnRightString = ""
 reso@cnFillOn  = False
 reso@cnLineLabelsOn = False
 reso@cnInfoLabelOn = False
 reso@cnLineColor   = "grey30"
 reso@cnLevels = (/0.9,1.0,1.1/)
 reso@tiYAxisString = ""
 reso@tmYLLabelsOn  = False
 reso@tmXBLabelsOn  = False
 do ic = 0,3
   dumo(ic)     = gsn_csm_contour(wks,spol_val,reso)
 end do
 ;;
 delete(reso@cnLevels)
 ;--------------------
 reso@cnLevels  = qcnlevs ;(/0.01,0.1,1.,2.,4./)
 reso@cnLineThicknessF = 2.0
 reso@cnLineColor = "grey20"
 ;reso@cnLineDashPattern = 15
 ;QC
 spol          = (/varh(0,:,:,ixs)/)
 dumq(0)       = gsn_csm_contour(wks,spol,reso)
 ;QI
 spol          = (/varh(2,:,:,ixs)/)
 dumq(3)       = gsn_csm_contour(wks,spol,reso)
 ;--------------------
 delete(reso@cnLevels)

 reso@cnLevels = (/0./)
 ;reso@tmXBLabelsOn  = True
 reso@cnLineColor = "black"
 reso@cnLineThicknessF = 3.
 do ic = 0,3
   dumt(ic)     = gsn_csm_contour(wks,tairC,reso)
 end do
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
 resu@vcLineArrowColor         = "azure3" ;"grey40"
 resu@vcLineArrowThicknessF    = 2.
 resu@vcMapDirection           = False
; do ic = 0,3
;   dumu(ic)     = gsn_csm_vector(wks,uy_yp,uz_yp,resu)
; end do

 delete(reso@cnLevels)
 reso@cnLineColor      = "grey20"
 reso@cnLineThicknessF = 1.
 reso@tmXBLabelsOn  = False
 reso@cnLevels   = (/5.,10.,20.,40./) ;;SPP Proposal (/10., 20., 30.,40./) ;(/5.,40./)
 reso@cnLineThicknessF  = 2.
 do ic = 0,3
   dumu(ic)     = gsn_csm_contour(wks,uz_yp,reso)
 end do

 fdebug = True 
 if (fdebug) then
   print("Generating debug plots ..")

   ;Read storm mask ------------------------------------------------------------------
   fsname   = "smask_" + pname + "_" + sfold + ".nc"
   diri_spath  = "/p/largedata/hbn33/storm_mask/" + expname + "/" + itype_gscp + "/nc/"

   nc_spath = addfile(diri_spath + "/" + fsname,"r")
   smsk     = nc_spath->Band1(fns,:,:)

   delete(smsk)
   smsk     = wrf_user_intrp3d(varh(qx,:,:,:),zmodel,"h",zheight,0.,False)

   smsk!0   = "Y"
   smsk!1   = "X"
   smsk&Y        = Ykm(:,150) ;
   smsk&X        = Xkm(150,:) ;

   res@gsnRightString = hydrometeors(qx) + " [" + hydunits + "]"
   delete(spol)
   ;zheight       = 5.0 ;km
   spol          = wrf_user_intrp3d(uz_hhl,zmodel,"h",zheight,0.,False)
   spol!0        = "Y"
   spol!1        = "X"
   spol&Y        = Ykm(:,150) ;
   spol&X        = Xkm(150,:) ;

   ;tr Resources
   ;-----------------------------
   pres@trXMinF = xycorner(1)
   pres@trYMinF = xycorner(0)
   pres@trXMaxF = xycorner(3)
   pres@trYMaxF = xycorner(2)
   pres@gsnLeftString = "Height : " + zheight + " km"

   pres@cnFillOn = True
   pres@cnFillPalette              = cmap
   pres@cnLevelSelectionMode       = "ExplicitLevels"
   pres@cnLevels = qcnlevs 
   pres@lbOrientation              = "vertical" 

   plot(0) = gsn_csm_contour(wks,smsk,pres)
   delete(pres@cnLevels)
   pres@cnFillOn = False
   pres@cnLineThicknessF = 3.0
   pres@cnLineDashPattern = 1
   pres@cnLevels = (/-3.,-2.,-1.,-0.5, 5.,10.,15.,20.,25.,30.,35.,40./)
   dumw = gsn_csm_contour(wks,spol,pres)
   dumw = ColorNegDashZeroPosContour(dumw,"blue","gray","red")
   overlay(plot(0),dumw)

   gres                  = True
   gres@gsLineColor      = "black"
   gres@gsLineDashPattern = 0
   gres@gsLineThicknessF = 3.0

   ypoly_cs = (/xycorner(0),xycorner(2)/)
   xpoly_cs = (/Xkm(0,ixs),Xkm(0,ixs)/)
   dum_cs = gsn_add_polyline(wks,plot(0),xpoly_cs,ypoly_cs,gres)

   draw(plot(0))
   frame(wks)
 else
 ;Plot and  Overlay
 do ic = 0,3
   if (ic.gt.0) then
     res@tiYAxisString = ""
     res@tmYLLabelsOn = False
   end if
   if (ic.eq.0) then
     res@gsnRightString = hydrometeors(ic+1) + " [" + hydunits + "]"
     spol          = (/varh(ic+1,:,:,ixs)/)
   else
     res@gsnRightString = hydrometeors(ic+2) + " [" + hydunits + "]"
     spol          = (/varh(ic+2,:,:,ixs)/)
   end if
   plot(ic) = gsn_csm_contour(wks,spol,res)
   overlay(plot(ic),dumo(ic))
   if (ic.eq.3) then
     overlay(plot(ic),dumq(ic))   ;QH overlay
   end if
   overlay(plot(ic),dumt(ic))
   if (ic.eq.1) then
     overlay(plot(ic),dumu(ic))
   end if
 end do

 pres                  = True
 pres@gsFillColor      = "gray"
 pres@gsLineThicknessF = 2.0          ; twice as thick
; pres@gsFillIndex      = 17           ; stipple, if desired

 do ip = 0,3
   dum0(ip) = gsn_add_polygon(wks,plot(ip),xpoly,ypoly,pres)
   dum1(ip) = gsn_add_polyline(wks,plot(ip),xpoly(0:sdims(1)-1),ypoly(0:sdims(1)-1),pres)
 end do

; Add legend
 lgres                    = True
 lgres@lgLineColors       = (/"black"/)
 lgres@lgLineThicknessF   = 2.
 lgres@lgDashIndexes      = 0
 lgres@lgItemType         = "Lines"
 lgres@lgLabelFontHeightF = 2.0
 lgres@vpWidthF           = 0.20
 lgres@vpHeightF          = 0.15
 lgres@lgPerimOn          = False
 lgres@lgLabelFontHeightF      = 0.15
 lbid1 = gsn_create_legend(wks,1,(/"W [ms~S~-1~N~]"/), lgres)
 lbid2 = gsn_create_legend(wks,1,(/hydrometeors(2) + " [" + hydunits + "]"/), lgres)

 amres = True
 amres@amParallelPosF   = 0.2
 amres@amOrthogonalPosF = -0.45
 annoid1 = gsn_add_annotation(plot(1),lbid1,amres)
 annoid3 = gsn_add_annotation(plot(3),lbid2,amres)

; Add annotation
 txres               = True                      ; text mods desired
 txres@txFontHeightF = 0.016                    ; text font height
 txres@txJust        = "CenterLeft"     
 gsn_text_ndc(wks,"d)",.001,.63,txres)

 ;Draw Panel
 resP = True
 resP@gsnMaximize    = True
 resP@gsnPanelXF = (/-1,0.26,0.46,.66/)  ;for 4 columns
; resP@gsnPanelDebug = True
 gsn_panel(wks,plot,(/1,4/), resP)

 end if

end
