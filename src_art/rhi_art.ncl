;Plot the c-s veiw of tracers along Y-axis
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

load "$prom_root/functions/lndistART.ncl"
load "$prom_root/functions/pm25ART.ncl"

begin
;STEP 1 User Settings
;-----------------------------------------------------------------------------------------
;-----------------------------------------------------------------------------------------
;#1 USER SETTINGS
;-----------------------------------------------------------------------------------------
  pname = "EMPA2_PROM1";
  inst  = 6
  lnd   = "HET"
  itype = "ART"
  diri  = "/p/largedata/hbn33/shrestha1/"

  nf          = 32
  glon_idB    = 7.071663 + 13*0.01 ; Bonn Radar
  glat_idB    = 50.73052

 ;--------------------------------------------------------------------------------
  itype_gscp  = lnd + itype
  sfold       = "runE_" + inst + "_" + itype_gscp
  expname     = lnd + "_" + pname
  subevent    = expname + "/" + itype_gscp + "/" + sfold
  diri_cos    = diri + "/" + subevent + "/a_cosout/"
  diri_cosI   = diri + "/" + subevent + "/a_cosout/ivr/"

;STEP 2 Set variables to read and specify scales or offsets
;-------------------------------------------------------------------------------------
 tracers = (/"NO2","SO2","VANTHA","VSOOT"/)
 fprocs  = True  ;tracers will be updated with processed Ncm3 and PM2.5

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

 nc_cos    = addfile(fils_cos(nf),"r")
 time      = nc_cos->time
 ;
 str_time  = cosmo_time_nclf(time,1)  ; 0 for date
 print(" " + fils_cos(nf))
 print("Processing ..." + nf + "  " + str_time)

 lat     = nc_cosI->lat
  lon     = nc_cosI->lon

 ;
 nm     = getind_latlon2d (lat,lon, glat_idB, glon_idB)
 iys    = nm(0,0)
 ixs    = nm(0,1) 


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

 if (fprocs) then
; Step 5b Retrieve Ncm3 and PM2.5
;--------------------------------------------------------------------------------
  ;Geomeric standard deviations of lognormal modes
  sig = (/1.7, 2.0, 1.7, 2.0, 1.4, 2.5, 1.7, 1.6, 1.5, 1.9, 2.0, 1.7/)

  ;Obtain aggregate of 12 log-normal distribution modes
  opt    = True
  opt@iy = -99 ;yloc
  opt@ix = -99 ;xloc
  opt@favg = True  ;Only for fout <1
  opt@fout = 1
  opt@nbdl = 0 

  smodes = ispan(1,12,1)


  nf_cos = dimsizes(fils_cos)
  fil_cos = fils_cos(nf)

  ;Ncm3
  do im = 0, dimsizes(smodes)-1
    opt@smode = smodes(im)
    if (im.eq.0) then
      temp = lndistART(fil_cos,opt)
      ndim = dimsizes(temp)
      varo = new((/12,ndim(0),ndim(1),ndim(2)/),"float")
      delete(ndim)
    end if
    varo(im,:,:,:) = lndistART(fil_cos,opt)
  end do

 ;Extract sums
  namemods   = (/"nuc./acc. (pure + mixed)", "soot", \
                 "dust", "PM2.5"/) ;"sea-salt"/)
  namemods   = (/"N~B~na~N~", "N~B~soot~N~","N~B~dust~N~","PM2.5"/)
  hydunits   = (/"cm~S~-3~N~","cm~S~-3~N~","cm~S~-3~N~","~F8~m~F21~g~N~m~S~-3~N~"/)

 ;Extract sums for cropped region
  varo1 = dim_sum_n(varo(0:3,:,:,:),0)     ;nuc-acc mode ;small
;  varo2 = dim_sum_n(varo(4:5,:,:,:),0)     ;soot+coarse
  varo2 = varo(4,:,:,:)    ;soot
  varo3 = dim_sum_n(varo(6:8,:,:,:),0)     ;dust
;  varo4 = dim_sum_n(varo(9:11,:,:,:),0)   ;sea

 ;PM2.5
  popt  = True
  popt@iy = iys
  popt@ix = ixs 
  popt@favg = False
  popt@fdry = True
  popt@nbld = 0
  popt@f3d  = True

  pm25 = pm25ART(fil_cos,popt)   ;Replace sea aerosol and concn very low

  ;;
  varh(0,:,:,:) = varo1
  varh(1,:,:,:) = varo2
  varh(2,:,:,:) = varo3
  varh(3,:,:,:) = pm25 

  tracers = namemods
 end if ;fprocs

;STEP 6 Variable to plot
;------------------------------------------------------------------------------
 spol          = dim_sum_n(varh(:,:,:,ixs),0)         ;
 spol_val      = where(spol.gt.0.01,1.,0)
 spol_edge     = convol(spol_val)

 sdims         = dimsizes(spol)
 spol!0        = "Z"
 spol!1        = "Y"
 spol&Y        = lat(:,ixs) ; 
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
 dumt     = new(4,"graphic")
 dumu     = new(4,"graphic")

 res      = True
 wks      = gsn_open_wks("png", "Figure_art_rhi_"+nf)

 cmap = read_colormap_file("ncl_default")
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
 res@vpHeightF = 0.55

 ;tr resources
 ;--------------------
 res@trXMinF    = glat_idB - .6
 res@trXMaxF    = glat_idB + 1.2 
 res@trYMaxF    = 3.

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
 if (fprocs) then
   res@trYMaxF = 15.
   res@cnLevelSelectionMode       = "ExplicitLevels"
   res@cnLevels = (/0.05,   \
                  0.1, 0.2, 0.4, 0.6, 0.8, \
                  1., 2., 4., 6., 8., \
                  10., 20., 40., 60., 80., \
                  100., 200., 400., 600., 800., \
                  1000., 2000., 4000., 6000., 8000. /)
 else
   res@cnLevelSelectionMode       = "ExplicitLevels"
   res@cnLevels = (/0.001,   \
                  0.01, 0.02, 0.04, 0.06, 0.08, \
                  0.1, 0.2, 0.4, 0.6, 0.8, \
                  1., 2., 4., 6., 8./)
 end if

 ;For Overlaying plots
 ;--------------------
 reso@cnInfoLabelOn             = False
 reso@cnLevelSelectionMode       = "ExplicitLevels"
 reso@cnLevels = (/0./)
 reso@tmXBLabelsOn  = False ;True
 reso@cnLineLabelFontHeightF = 0.02
 reso@cnLineThicknessF = 3.
 reso@tiXAxisString =""
 reso@tiYAxisString = ""

 do ic = 0,3
   dumt(ic)     = gsn_csm_contour(wks,tairC,reso)
 end do
 ;-----------------------
 ;;Vector overlay along the cross-section
 resu@tmYLLabelsOn  = False
 resu@tmXBLabelsOn  = False

 resu@tiYAxisString = ""
 resu@tiXAxisString = ""
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
 res@gsnRightString = tracers(tc) + " [" + hydunits(tc) + "]" 
 spol          = (/varh(tc,:,:,ixs)/)  
 res@tiXAxisString = " "
 res@tmXBLabelsOn  = False
 plot(0)  = gsn_csm_contour(wks,spol,res)

 ;
 tc = 1
 res@tiYAxisString = ""
 res@tmYLLabelsOn  = False
 res@gsnRightString = tracers(tc)  + " [" + hydunits(tc) + "]"
 spol          = (/varh(tc,:,:,ixs)/)   
 if (fprocs) then
  delete(res@cnLevels)
  res@cnLevels = (/0.05,   \
                  0.1, 0.2, 0.4, 0.6, 0.8, \
                  1., 2., 4., 6., 8., \
                  10., 20., 40., 60., 80./)
 end if
 plot(1)  = gsn_csm_contour(wks,spol,res)

 ;
 tc = 2
 res@gsnRightString = tracers(tc)  + " [" + hydunits(tc) + "]"
 spol          = (/varh(tc,:,:,ixs)/)
 if (fprocs) then
  delete(res@cnLevels)
  res@cnLevels = (/0.05,   \
                  0.1, 0.2, 0.4, 0.6, 0.8, 1.0/)
 end if
 res@tmXBLabelsOn  = True
 res@tiXAxisString = "Latitude [~F34~0~F~N]"
 res@tiYAxisString = "Height [km]"
 plot(2)  = gsn_csm_contour(wks,spol,res)

 ;
 tc = 3
 res@gsnRightString = tracers(tc)  + " [" + hydunits(tc) + "]"
 spol          = (/varh(tc,:,:,ixs)/)
 if (fprocs) then
  delete(res@cnLevels)
  cnlev = ispan(0,20,2)*1.0 ;
  cnlev(0) = 0.1
  res@cnLevels = cnlev
 end if
 res@tiYAxisString = ""
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
 resP@gsnPanelYF = (/-1,-1,0.47,.47/)  ;for 4 columns
 resP@gsnPanelDebug = True
 ;CPS resP@gsnPanelMainString = str_time
 gsn_panel(wks,plot,(/2,2/), resP)
end
