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


begin
;STEP 1 User Settings
;-----------------------------------------------------------------------------------------
 ;pid       = 0 
 opt       = 1   ;0-NX, 1-QX , 2-DX

 ydates    = (/"20150705","20160513","20170706"/)

 insts     = (/6, 8, 11/)
 fns       = (/11,34,18/) - 11
 ixs       = (/90,124,95/) + 65  ;Need to add offset due to clipping of radar

 ;;extract_ascii clips
 ;;Extent for clipping :58 237
 ;;Extent for clipping :65 244

 lat_ext1  = (/51.2, 50.8, 51./)
 lat_ext0 = (/50.0+40./60.,50.+10./60.,50.42/)

 ;SPP proposal
 ;insts     = (/5, -1, -1/)
 ;fns       = (/21,-1,-1/)
 ;ixs       = (/184,-1,-1/)
 ;lat_ext1  = (/51.6, -1, -1/)
 ;----

 itype     = "HET2483"
 diri        = "/p/largedata/hbn33/shrestha1/HET_PROM1"
 ;diri      = "/daten04/pshrestha/shared-data/bonnRadar_" + ydates(pid)
 event  = "runE_" + insts(pid) + "_" + itype
 diri_cos  = diri + "/" + itype + "/" + event + "/cosout/"
 diri_cosI = diri + "/" + itype + "/" + event + "/cosout/ivr/"

;STEP 2 Set variables to read and specify scales or offsets
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

 print("Processing " + fils_cos(fns(pid)))

 nc_cos    = addfile(fils_cos(fns(pid)),"r")
 time      = nc_cos->time
 ;
 str_time  = cosmo_time_nclf(time,1)  ; 0 for date
 print("Processing ..." + fns(pid) + "  " + str_time)

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
 spol          = dim_sum_n(varh(:,:,:,ixs(pid)),0)         ;
 spol_val      = where(spol.gt.0.01,1.,0)
 spol_edge     = convol(spol_val)

 sdims         = dimsizes(spol)
 spol!0        = "Z"
 spol!1        = "Y"
 spol&Y        = lat(:,ixs(pid)) ; 
 xval2x        = conform_dims(sdims,spol&Y,1)*1.
 zval2x        = zmodel(:,:,ixs(pid))

 tairC         = tair(:,:,ixs(pid)) - 273.16
 uy_yp         = uy(:,:,ixs(pid))
 uz_yp         = uz_hhl(:,:,ixs(pid))

 copy_VarCoords(spol,tairC)   
 copy_VarCoords(spol, uy_yp)  
 copy_VarCoords(spol, uz_yp) 

 ;Need this for vector plot as stretching does not work
 uy_yp&Z      = zmodel(:,90,ixs(pid))
 uz_yp&Z      = zmodel(:,90,ixs(pid))
 
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
 wks      = gsn_open_wks("x11", "Figure_QX_rhi_"+pid)

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
 res@tmXBLabelFontThicknessF    = 2.0
 res@tmXBMode   = "Explicit"
 res@tmXBValues = xval2x(0,ispan(0,sdims(1)-1,20))
 res@tmXBLabels = tointeger((res@tmXBValues)*100.)/100.
 res@tmXBLabelsOn  = False 

 ;ti Resources
 ;-----------------
 res@tiXAxisFont                = fncps_no
 res@tiYAxisFont                = fncps_no
 res@tiXAxisFontHeightF         = fncps_hgt
 res@tiYAxisFontHeightF         = fncps_hgt
 res@tiXAxisString = "Latitude [~F34~0~F~N]"
 res@tiYAxisString = "Height [km]"

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
 res@trXMinF    = lat_ext0(pid)
 res@trXMaxF    = lat_ext1(pid)
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

 ;For Overlaying plots
 ;--------------------
 reso@cnLevelSelectionMode       = "ExplicitLevels"
 reso@cnFillOn  = False
 reso@cnLineLabelsOn = False
 reso@cnInfoLabelOn = False
 reso@cnLineColor   = "grey30"
 reso@cnLevels = (/0.9,1.0,1.1/)
 reso@tiYAxisString = ""
 reso@tmYLLabelsOn  = False
 do ic = 0,3
   dumo(ic)     = gsn_csm_contour(wks,spol_val,reso)
 end do
 ;;
 delete(reso@cnLevels)
 ;--------------------
 reso@cnLevels  = (/0.01,0.1,1.,2.,4./)
 reso@cnLineThicknessF = 2.0
 reso@cnLineColor = "black"
 ;QC
 spol          = (/varh(0,:,:,ixs(pid))/)
 dumq(0)       = gsn_csm_contour(wks,spol,reso)
 ;QH
 spol          = (/varh(5,:,:,ixs(pid))/)
 dumq(3)       = gsn_csm_contour(wks,spol,reso)
 ;--------------------
 delete(reso@cnLevels)
 reso@cnLevels = (/0./)
 reso@tmXBLabelsOn  = True
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
 reso@cnLineColor      = "black"
 reso@cnLineThicknessF = 1.
 reso@tmXBLabelsOn  = False
 reso@cnLevels   = (/5.,40./) ;;SPP Proposal (/10., 20., 30.,40./) ;(/5.,40./)
 reso@cnLineThicknessF  = 2.
 do ic = 0,3
   dumu(ic)     = gsn_csm_contour(wks,uz_yp,reso)
 end do
 

  
 ;QR
 res@gsnRightString = hydrometeors(1) + " [" + hydunits + "]" 
 spol          = (/varh(1,:,:,ixs(pid))/)  
 plot(0)  = gsn_csm_contour(wks,spol,res)

 ;QI
 res@tiYAxisString = ""
 res@tmYLLabelsOn  = False
 res@gsnRightString = hydrometeors(2) + " [" + hydunits + "]"
 spol          = (/varh(2,:,:,ixs(pid))/)   

 plot(1)  = gsn_csm_contour(wks,spol,res)

 ;QS
 res@gsnRightString = hydrometeors(3) + " [" + hydunits + "]"
 spol          = (/varh(3,:,:,ixs(pid))/)
 plot(2)  = gsn_csm_contour(wks,spol,res)

 ;QG
 res@gsnRightString = hydrometeors(4) + " [" + hydunits + "]"
 spol          = (/varh(4,:,:,ixs(pid))/)
 plot(3)  = gsn_csm_contour(wks,spol,res) 

 ;Now Overlay
 do ic = 0,3
   overlay(plot(ic),dumo(ic))
   if (ic.eq.3) then
     overlay(plot(ic),dumq(ic))
   end if
   overlay(plot(ic),dumt(ic))
   if (ic.eq.2) then
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
 lbid1 = gsn_create_legend(wks,1,(/"W [ms~S~-1~N~]"/), lgres)
 lbid2 = gsn_create_legend(wks,1,(/hydrometeors(5) + " [" + hydunits + "]"/), lgres)

 amres = True
 amres@amParallelPosF   = 0.2
 amres@amOrthogonalPosF = -0.45
; annoid1 = gsn_add_annotation(plot(2),lbid1,amres)
; annoid3 = gsn_add_annotation(plot(3),lbid2,amres)

; Add annotation
 txres               = True                      ; text mods desired
 txres@txFontHeightF = 0.016                    ; text font height
 txres@txJust        = "CenterLeft"     
 gsn_text_ndc(wks,"d)",.05,.65,txres)

 ;Draw Panel
 resP = True
 resP@gsnMaximize    = True
 ;resP@gsnPanelXF = (/-1,0.26,0.46,.66/)  ;for 4 columns
 ;resP@gsnPanelDebug = True
 gsn_panel(wks,plot,(/2,2/), resP)
end
