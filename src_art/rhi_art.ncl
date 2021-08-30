;Plot the c-s veiw of hydrometeors alond Y-axis
; specified model grid in X using linearized 

;P. Shrestha
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$PROM/functions/shade_terrain.ncl"
load "$PROM/functions/cosmo_time_nclf.ncl"
load "$PROM/functions/calc_rho_nclf.ncl"
load "$PROM/functions/convol.ncl"


begin
;STEP 1 User Settings
;-----------------------------------------------------------------------------------------
 pid       = 0
 ydates    = (/"20150705","20160513","20170706"/)

 insts     = (/6, 8, 11/)
 fns       = (/1,34,18/)
 ixs       = (/90,124,95/) + 65  ;Need to add offset due to clipping of radar

 lat_ext1  = (/51.2, 50.8, 51./)
 lat_ext0 = (/50.0+40./60.,50.+10./60.,50.42/)

 diri      = "/daten04/pshrestha/shared-data/ART/"

 diri_cos  = diri 
 diri_cosI = diri 

;STEP 2 Set variables to read and specify scales or offsets
;-------------------------------------------------------------------------------------
 tracers = (/"NO","SO2","ETH","VSOOT"/)

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

;STEP 5 Extract tracers 
;-------------------------------------------------------------------------------------

 nh           = dimsizes(tracers)
 varh         = new((/nh,nz,ny,nx/),"float")    ;Initialize variable for plotting
 varh@_FillValue = -999.
 do ih = 0, nh-1
   varh(ih,:,:,:) = (/nc_cos->$(tracers(ih))$(0,::-1,:,:)/)
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
 res@trYMaxF    = 5.

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

 ;For Overlaying plots
 ;--------------------
 reso@cnInfoLabelOn             = False
 reso@cnLevelSelectionMode       = "ExplicitLevels"
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

 ;
 tc = 0
 res@gsnRightString = tracers(tc) ;+ " [" + hydunits + "]" 
 spol          = (/varh(tc,:,:,ixs(pid))/)  
 plot(0)  = gsn_csm_contour(wks,spol,res)

 ;
 tc = 1
 res@tiYAxisString = ""
 res@tmYLLabelsOn  = False
 res@gsnRightString = tracers(tc) ; + " [" + hydunits + "]"
 spol          = (/varh(tc,:,:,ixs(pid))/)   

 plot(1)  = gsn_csm_contour(wks,spol,res)

 ;
 tc = 2
 res@gsnRightString = tracers(tc) ; + " [" + hydunits + "]"
 spol          = (/varh(tc,:,:,ixs(pid))/)
 plot(2)  = gsn_csm_contour(wks,spol,res)

 ;
 tc = 3
 res@gsnRightString = tracers(tc) ; + " [" + hydunits + "]"
 spol          = (/varh(tc,:,:,ixs(pid))/)
 plot(3)  = gsn_csm_contour(wks,spol,res) 

 ;Now Overlay
 do ic = 0,3
   overlay(plot(ic),dumt(ic))
   ;overlay(plot(ic),dumu(ic))
 end do

 pres                  = True
 pres@gsFillColor      = "gray"
 pres@gsLineThicknessF = 2.0          ; twice as thick
; pres@gsFillIndex      = 17           ; stipple, if desired

 do ip = 0,3
   dum0(ip) = gsn_add_polygon(wks,plot(ip),xpoly,ypoly,pres)
   dum1(ip) = gsn_add_polyline(wks,plot(ip),xpoly(0:sdims(1)-1),ypoly(0:sdims(1)-1),pres)
 end do

 ;Draw Panel
 resP = True
 resP@gsnMaximize    = True
 ;resP@gsnPanelXF = (/-1,0.26,0.46,.66/)  ;for 4 columns
 ;resP@gsnPanelDebug = True
 gsn_panel(wks,plot,(/2,2/), resP)
end
