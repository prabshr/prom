;======================================================
; Hydrometeor profiles
;======================================================
; Prabhakar Shrestha
; E-mail: pshrestha@uni-bonn.de
;=======================================================

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"

begin

;-----------------------------------------------------------------------------------------
; USER SETTINGS
;-----------------------------------------------------------------------------------------

 diri        = "/p/scratch/chbn33/hbn331"
 event       = "bonnRadar_2014111500"
 case        = "HET2473"

 diri_cos    = diri + "/" + event + "/" + case + "/cosout/day2/"
 diri_cosI   = diri + "/" + event + "/" + case + "/cosout/ivr/"

 ;Bonn Radar
 glon_id = (/7.071663/) ;
 glat_id = (/50.73052/) ;

 gloc    = 0       ;Centre point of cross-section
 Ipangle = 298.6   ;128. ;117.   ;153.   ;w.r.t -Y-axis clockwise        ;angle of cross-section passing through gloc 
 hloc    = 90 ;70      ;max 99, height of cross-section for plotting
 hlev    = 49 
 pressure = 700 ;350.   ;hPa, for pressure plot
;177
 filindex      = ispan(0,144,1) ;ispan(200,201,1)
 ext           = 0.91 ;0.25 ;1.5 ;
 mpcorner      = (/glat_id(0)-ext,glon_id(0)-ext,glat_id(0)+ext,glon_id(0)+ext/)
 ;mpcorner      = (/50.5, 6.5, 51.1, 7.15/)
;-----------------------------------------------------------------------------------------
; READ INVARIANT FIELDS
;-------------------------------------------------------------------------------------
 nc_cosI = addfile(diri_cosI +"lfff00000000c.nc","r")
 lat = nc_cosI->lat
 lon = nc_cosI->lon
 topo = nc_cosI->HSURF(0,:,:)
 topo = smth9_Wrap(topo, 0.50, 0.25, False)

 
 hhl = nc_cosI->HHL(0,::-1,:,:)
 ndimhhl   = dimsizes(hhl)
 nz        = min((/hlev,ndimhhl(0)-1/))
 zmodel    = 0.5*(hhl(0:nz-1,:,:) + hhl(1:nz,:,:))
 ndim = dimsizes(lat)

;-----------------------------------------------------------------------------------------
; READ NC FILES 
;-----------------------------------------------------------------------------------------
 fils_cos  = systemfunc("ls " + diri_cos + "lf*")
 nf_cos   = dimsizes(fils_cos)
;Initialize empty arrays variables
 nf_ctr = dimsizes(filindex) 
;
 do fnctr = 0,nf_ctr-1   ;;BIG LOOP OVER FILES 
   
   fn = filindex(fnctr)
   nc_cos = addfile(fils_cos(fn),"r")
   time   = nc_cos->time
   dtime = calendar_decode2(time,-5)
   year  = dtime(0,0)
   month = dtime(0,1)
   day   = dtime(0,2)
   hh    = dtime(0,3)
   mm    = dtime(0,4)
   ss    = dtime(0,5)
 
   str_time  = sprinti("%0.4i",year) + "-" + sprinti("%0.2i",month) + \
               "-" + sprinti("%0.2i",day) + " " + sprinti("%0.2i",hh) + \
                sprinti("%0.2i",mm) ;+ " UTC" ;sprinti("%0.2i",ss) +  " UTC" 
   print("  " + fils_cos(fn))
   print(str_time + " ")
;------------------------------------------------------------------------------
;;Prognostic state variables 
;------------------------------------------------------------------------------
   ps          = nc_cos->P(0,:,:,:)
   tk          = nc_cos->T(0,:,:,:)
   qv          = nc_cos->QV(0,:,:,:)
   ux          = nc_cos->U(0,:,:,:)
   uy          = nc_cos->V(0,:,:,:)
   uz          = nc_cos->W(0,:,:,:)
   uz_hhl      = uz(1:nz,:,:)*0.
   do kk =  0,nz-1
     uz_hhl(kk,:,:)  = 0.5*(uz(kk,:,:) + uz(kk+1,:,:))
   end do
;-------------------------------------------------------------------------------
;;Prognostic hydrometeors, mixing ratios
;-------------------------------------------------------------------------------
   qc          = nc_cos->QC(0,:,:,:)
   qi          = nc_cos->QI(0,:,:,:)
   qg          = nc_cos->QG(0,:,:,:)
   qs          = nc_cos->QS(0,:,:,:)
   qh          = nc_cos->QH(0,:,:,:)
   qr          = nc_cos->QR(0,:,:,:)

;-------------------------------------------------------------------------------
;;Prognostic hydrometeors, number concentrationss
;-------------------------------------------------------------------------------
   nqc          = nc_cos->NCCLOUD(0,:,:,:)
   nqi          = nc_cos->NCICE(0,:,:,:)
   nqg          = nc_cos->NCGRAUPEL(0,:,:,:)
   nqs          = nc_cos->NCSNOW(0,:,:,:)
   nqh          = nc_cos->NCHAIL(0,:,:,:)
   nqr          = nc_cos->NCRAIN(0,:,:,:)

;-------------------------------------------------------------------------------
;;Crop data and vertical plane for interpolation
;-------------------------------------------------------------------------------
   if (fnctr.eq.0) then
     ji = region_ind (lat,lon, mpcorner(0),mpcorner(2),mpcorner(1),mpcorner(3))
    ; plane        = (/33,10/)
     latc = lat(ji(0):ji(1),ji(2):ji(3))
     lonc = lon(ji(0):ji(1),ji(2):ji(3))
     zmodelc = zmodel(0:nz-1,ji(0):ji(1),ji(2):ji(3))
     nm = getind_latlon2d (latc,lonc, glat_id(gloc), glon_id(gloc))
     plane = (/nm(0,1),nm(0,0)/)     ;x and y co-ordinates 
     ;;0 : SN , 90 :WE
     lon_plane  = wrf_user_intrp2d(lonc,plane,Ipangle,False)
     lat_plane  = wrf_user_intrp2d(latc,plane,Ipangle,False)
     ;;;hgt_plane  = wrf_user_intrp3d(zmodelc,zmodelc,"v",plane,Ipangle,False)
     hgt_coord     =  fspan(min(zmodelc),max(zmodelc),100)/1000.  ;[km] ;WRF3D interplates into 100 equal ht levels
     hgt_coord@units = "km"
     ;
     tmp_dim    = dimsizes(latc)
     qx_ptime   = new((/nf_ctr,tmp_dim(0),tmp_dim(1)/),"float")
   end if
 
   ;crop
    lat_c       = lat(ji(0):ji(1),ji(2):ji(3))
    lon_c       = lon(ji(0):ji(1),ji(2):ji(3))
    topo_c      = topo(ji(0):ji(1),ji(2):ji(3))
    tk_c        = tk(::-1,ji(0):ji(1),ji(2):ji(3))
    qv_c        = qv(::-1,ji(0):ji(1),ji(2):ji(3))
    ps_c        = ps(::-1,ji(0):ji(1),ji(2):ji(3))
    ux_c        = ux(::-1,ji(0):ji(1),ji(2):ji(3))
    uy_c        = uy(::-1,ji(0):ji(1),ji(2):ji(3))
    uz_hhl_c    = uz_hhl(::-1,ji(0):ji(1),ji(2):ji(3))

;kg/kg to g/kg
    qc_c        = qc(::-1,ji(0):ji(1),ji(2):ji(3))*1000. 
    qi_c        = qi(::-1,ji(0):ji(1),ji(2):ji(3))*1000.
    qg_c        = qg(::-1,ji(0):ji(1),ji(2):ji(3))*1000.
    qr_c        = qr(::-1,ji(0):ji(1),ji(2):ji(3))*1000.
    qs_c        = qs(::-1,ji(0):ji(1),ji(2):ji(3))*1000. 
    qh_c        = qh(::-1,ji(0):ji(1),ji(2):ji(3))*1000. 

;-------------------------------------------------------------------------------
;;Interpolate Cropped Data in Vertical Plane
;-------------------------------------------------------------------------------
;  We further crop the data in vertical to get high resolution interpolation...

    tc_plane    = wrf_user_intrp3d(tk_c(0:nz-1,:,:),zmodelc,"v",plane,Ipangle,False) -273.15  ;C 
    ux_plane    = wrf_user_intrp3d(ux_c(0:nz-1,:,:),zmodelc,"v",plane,Ipangle,False)
    uy_plane    = wrf_user_intrp3d(uy_c(0:nz-1,:,:),zmodelc,"v",plane,Ipangle,False)
    uz_plane    = wrf_user_intrp3d(uz_hhl_c(0:nz-1,:,:),zmodelc,"v",plane,Ipangle,False)
    qc_plane    = wrf_user_intrp3d(qc_c(0:nz-1,:,:),zmodelc,"v",plane,Ipangle,False)
    qi_plane    = wrf_user_intrp3d(qi_c(0:nz-1,:,:),zmodelc,"v",plane,Ipangle,False)
    qg_plane    = wrf_user_intrp3d(qg_c(0:nz-1,:,:),zmodelc,"v",plane,Ipangle,False)
    qs_plane    = wrf_user_intrp3d(qs_c(0:nz-1,:,:),zmodelc,"v",plane,Ipangle,False)
    qr_plane    = wrf_user_intrp3d(qr_c(0:nz-1,:,:),zmodelc,"v",plane,Ipangle,False)
    qh_plane    = wrf_user_intrp3d(qh_c(0:nz-1,:,:),zmodelc,"v",plane,Ipangle,False)

    qx_plane     = qc_plane + qi_plane + qg_plane + qs_plane + qh_plane + qr_plane
    print("Debug1 " + qc_plane@Orientation )
    print("Debug1.1 " + max(uz_plane) + "m/s  " + min(uz_plane) + "m/s")
    print("Debug1.1 " + max(qi_plane) + "g/kg  " + min(qi_plane) + "g/kg")
    print("Debug1.1 " + max(qs_plane) + "g/kg  " + min(qs_plane) + "g/kg")
    print("Debug1.1 " + max(qg_plane) + "g/kg  " + min(qg_plane) + "g/kg")
    print("Debug1.1 " + max(qh_plane) + "g/kg  " + min(qh_plane) + "g/kg")


    tc_plane!0            = "height"
    tc_plane!1            = "longitude"
    tc_plane&height       = (/hgt_coord/)
    tc_plane&longitude    = (/lon_plane/)
    tc_plane@long_name    = "T"
    tc_plane@units        = "C"

    qc_plane!0            = "height"
    qc_plane!1            = "longitude"
    qc_plane&height       = (/hgt_coord/)
    qc_plane&longitude    = (/lon_plane/)
    qc_plane@long_name    = "QC"
    qc_plane@units        = "g/kg"

    qi_plane!0            = "height"
    qi_plane!1            = "longitude"
    qi_plane&height       = (/hgt_coord/)
    qi_plane&longitude    = (/lon_plane/)
    qi_plane@long_name    = "QI"
    qi_plane@units        = "g/kg" 

    qg_plane!0            = "height"
    qg_plane!1            = "longitude"
    qg_plane&height       = (/hgt_coord/)
    qg_plane&longitude    = (/lon_plane/)
    qg_plane@long_name    = "QG"
    qg_plane@units        = "g/kg"

    qs_plane!0            = "height"
    qs_plane!1            = "longitude"
    qs_plane&height       = (/hgt_coord/)
    qs_plane&longitude    = (/lon_plane/)
    qs_plane@long_name    = "QS"
    qs_plane@units        = "g/kg"

    qh_plane!0            = "height"
    qh_plane!1            = "longitude"
    qh_plane&height       = (/hgt_coord/)
    qh_plane&longitude    = (/lon_plane/)
    qh_plane@long_name    = "QH"
    qh_plane@units        = "g/kg"

    qr_plane!0            = "height"
    qr_plane!1            = "longitude"
    qr_plane&height       = (/hgt_coord/)
    qr_plane&longitude    = (/lon_plane/)
    qr_plane@long_name    = "QR"
    qr_plane@units        = "g/kg"

    ux_plane!0         = "height"
    ux_plane!1         = "longitude"
    ux_plane&height    = (/hgt_coord/)
    ux_plane&longitude = (/lon_plane/)

    uz_plane!0         = "height"
    uz_plane!1         = "longitude"
    uz_plane&height    = (/hgt_coord/)
    uz_plane&longitude = (/lon_plane/)

;-------------------------------------------------------------------------------
;;Extract sounding data along cross-section
;-------------------------------------------------------------------------------
    qr_sum = dim_sum_n(qr_plane,0)
    sloc = ind(qr_sum.eq.max(qr_sum))
    print(sloc)
    if (ismissing(sloc).or.sloc.eq.0.or.sloc.ge.dimsizes(lon_plane)-1) then
       sloc = tointeger(dimsizes(lon_plane)/2.) 
    end if
;    ht_v = hgt_plane(0:hloc,sloc)/1000.
;    ht_v(0) = 0.
    ht_v  = hgt_coord(0:hloc) 
    tc_v = dim_avg_n(tc_plane(0:hloc,sloc-1:sloc+1),1)
    qc_v = dim_avg_n(qc_plane(0:hloc,sloc-1:sloc+1),1)
    qi_v = dim_avg_n(qi_plane(0:hloc,sloc-1:sloc+1),1)
    qg_v = dim_avg_n(qg_plane(0:hloc,sloc-1:sloc+1),1)
    qs_v = dim_avg_n(qs_plane(0:hloc,sloc-1:sloc+1),1)
    qh_v = dim_avg_n(qh_plane(0:hloc,sloc-1:sloc+1),1)
    qr_v = dim_avg_n(qr_plane(0:hloc,sloc-1:sloc+1),1)
    ux_v = dim_avg_n(ux_plane(0:hloc,sloc-1:sloc+1),1)
    uy_v = dim_avg_n(uy_plane(0:hloc,sloc-1:sloc+1),1)
    uz_v = dim_avg_n(uz_plane(0:hloc,sloc-1:sloc+1),1)

    ;;Give a minimum number for hydrometeor to make Log Plots
    qminval = 1.E-5
    do kk = 0, hloc
      qc_v(kk) = max((/qc_v(kk),qminval/))
      qi_v(kk) = max((/qi_v(kk),qminval/))
      qs_v(kk) = max((/qs_v(kk),qminval/))
      qg_v(kk) = max((/qg_v(kk),qminval/))
      qr_v(kk) = max((/qr_v(kk),qminval/))
      qh_v(kk) = max((/qh_v(kk),qminval/))
    end do
;-------------------------------------------------------------------------------
;;Interpolate Cropped Data in horizontal at fixed pressure 
;-------------------------------------------------------------------------------
;
   ps_c          = ps_c*0.01
   ux_p          = wrf_user_intrp3d(ux_c,ps_c,"h",pressure,0.,False)
   uy_p          = wrf_user_intrp3d(uy_c,ps_c,"h",pressure,0.,False)
   qx_p          = wrf_user_intrp3d(qs_c,ps_c,"h",pressure,0.,False)

   ux_p@lat2d    = lat_c
   ux_p@lon2d    = lon_c
   uy_p@lat2d    = lat_c
   uy_p@lon2d    = lon_c

   qx_p@lat2d    = lat_c
   qx_p@lon2d    = lon_c
    
   topo_c@lat2d  = lat_c
   topo_c@lon2d  = lon_c 

;-----------------------------------------------------------------------------------------
; Plotting Resource
;-----------------------------------------------------------------------------------------
  plot = new(3,graphic)
  wks_type = "png"
  fname_out = "Figure09_"+case + "_"+sprinti("%0.4i",fn)
  wks = gsn_open_wks(wks_type,fname_out)

  gsn_define_colormap(wks,"precip2_17lev")

  fncps_no = ":F21:"
  fncps_no = 21
  fncps_hgt = 0.02

;-----------------------------------------------------------------------------------------
; Pressure plot resources 
;-----------------------------------------------------------------------------------------

 res                            = True                 ; plot mods desired
 res@gsnDraw                    = False
 res@gsnFrame                   = False
 restopo                        = res
 resln                          = res
; res@gsnAddCyclic               = False                ; regional data
 res@gsnLeftString              = str_time
 res@gsnRightString             = sprintf("%0.0f",pressure) + " hPa"
 res@gsnMaximize                = False
 res@gsnLeftStringFontHeightF   = fncps_hgt
 res@gsnRightStringFontHeightF  = fncps_hgt
 
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
; res@mpCenterLonF                 = lon@Longitude_of_southern_pole           ;lon of centre of map projection system
; res@mpCenterLatF                 = 90.0 + lat@Latitude_of_southern_pole         ;lat of centre of map projection system

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
 res@cnLinesOn                  = False 
 res@cnInfoLabelOn              = False   ; Turn off annotations
 res@cnLineLabelsOn             = False
 res@cnLevelSelectionMode       = "ExplicitLevels"
 res@cnLevels                   = ispan(1,10, 1) * 0.1 ; g/kg 


 ;Labels
 res@lbLabelAutoStride          = True
 res@lbOrientation              = "vertical"
 res@lbTitleString              = "g/kg"
 res@lbLabelFontHeightF         = 0.02
 res@pmLabelBarWidthF           = 0.05
 res@pmLabelBarHeightF          = 0.35

 ;Vectors
 res@vcRefMagnitudeF          = 10.0 ;20.0             ; define vector ref mag
 res@vcRefLengthF             = 0.045            ; define length of vec ref
 res@vcGlyphStyle             = "CurlyVector"    ; turn on curly vectors
 res@vcMinDistanceF           = 0.032            ; thin vectors
 res@vcRefAnnoOrthogonalPosF  = -.3             ; move ref vector down
 res@vcLineArrowColor         = "grey40"
 res@vcLineArrowThicknessF    = 2.
 res@vcGlyphOpacityF          = 0.7
 res@vcRefAnnoOn              = True

 ;Viewport
 res@vpHeightF                = 0.35
 res@vpWidthF                 = 0.7
 res@vpYF                     = 0.94
 res@vpXF                     = 0.2

 res@gsnScalarContour         = True

 res@tiXAxisString                 = ""
 res@tiYAxisString                 = ""
 
 plot(0)                      = gsn_csm_vector_scalar_map(wks,ux_p,uy_p,qx_p, res)

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
 dumtopo                        = gsn_csm_contour(wks,topo_c, restopo)
 overlay(plot(0),dumtopo)
 
 ;Cross-section Line
 plres = True
 plres@gsLineDashPattern        = 3
 plres@gsLineThicknessF         = 3.
 plres@tfPolyDrawOrder          = "PostDraw"
 dumpl = gsn_add_polyline(wks, plot(0),lon_plane, lat_plane, plres)

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
; Cross-Section Plot Resources 
;-----------------------------------------------------------------------------------------
 ;Viewport
 rescs@vpXF                          = 0.39
 rescs@vpYF                          = 0.45
 rescs@vpHeightF                     = 0.35
 rescs@vpWidthF                      = 0.46

 rescsv = rescs

 rescs@gsnLeftString                 = "c-s AA'"
 rescs@gsnRightString                = ""
 rescs@tiXAxisString                 = "Longitude"
 rescs@tiYAxisString                 = ""

 rescs@trYMaxF                     = hgt_coord(hloc) ;hloc ;80
 rescs@trYMinF                     = 0.0
 rescs@tmYLLabelsOn                 = False
; loc_label                        = ispan(0,hloc,5)
; rescs@tmYLMode                   = "Explicit"
; rescs@tmYLValues                 = loc_label
; rescs@tmYLLabels                 = hgt_coord(loc_label)
 ;round(round(hgt_plane(loc_label,sloc)/10.,0)/100.,0)

 ;Label
 rescs@lbLabelAutoStride          = True
 rescs@lbOrientation              = "vertical"
 rescs@lbTitleString              = "g/kg"
 rescs@lbLabelFontHeightF         = 0.02
 rescs@pmLabelBarWidthF           = 0.05
 rescs@pmLabelBarHeightF          = 0.35

 ;Contour
 rescs@cnMissingValFillPattern       = 8        
 rescs@cnMissingValFillColor         = "black"
 rescs@cnMissingValPerimOn           = True            ; 
 rescs@cnMissingValPerimColor        = "black"         ;
 rescs@cnMissingValPerimDashPattern  = 1          ;
 rescs@cnMissingValPerimThicknessF   = 2.0
 rescs@cnLevelSelectionMode          = "ExplicitLevels"

; Cross-section contour
 hmax = 200 ;24
 hfac = 0.01
 hlev = 10 ;2   ;1
; rescs@cnLevels                      = ispan(1,hmax,hlev)*hfac
 rescs@cnLevels  = (/0.01, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0/)

; Cross-section contour

 rescs@cnFillOn                   = True 
 rescs@cnLinesOn                  = True ;False 
 rescs@cnInfoLabelOn              = False   ; Turn off annotations
 rescs@cnLineLabelsOn             = True; False

 plot(1)                          = gsn_csm_contour(wks, qr_plane, rescs)

 printVarSummary(qr_plane)

 rescs@cnLineThicknessF           = 2.0
 rescs@cnFillOn                   = False
 rescs@cnLineColor                = "purple" 
 dumqp                            = gsn_csm_contour(wks, qi_plane, rescs) 
; overlay(plot(1),dumqp)

 rescs@cnLineColor                = "red"
 dumqp                            = gsn_csm_contour(wks, qg_plane, rescs)
; overlay(plot(1),dumqp)

 rescs@cnLineColor                = "orange"
 dumqp                            = gsn_csm_contour(wks, qs_plane, rescs)
; overlay(plot(1),dumqp)

 rescs@cnLineColor                = "black"
 dumqp                            = gsn_csm_contour(wks, qc_plane, rescs)
; overlay(plot(1),dumqp)

 rescs@cnLineColor                = "magenta"
 dumqp                            = gsn_csm_contour(wks, qh_plane, rescs)
 overlay(plot(1),dumqp)
 
 ;;Vector overlay along the cross-section

 ;Vector 
 rescsv@vcRefAnnoFontHeightF     = 0.02
 rescsv@vcRefMagnitudeF          = 10. ;20.            ; define vector ref mag
 rescsv@vcRefLengthF             = 0.045            ; define length of vec ref
 rescsv@vcGlyphStyle             = "CurlyVector"    ; turn on curly vectors
 rescsv@vcMinDistanceF           = 0.019            ; thin vectors
 rescsv@vcRefAnnoOrthogonalPosF  = -0.16 ;0.46               ; move ref vector down
 rescsv@vcLineArrowColor         = "grey40"
 rescsv@vcLineArrowThicknessF    = 2.
 rescsv@gsnLeftString            = ""
 rescsv@gsnRightString           = ""

 ;Make a copy for multiple plots
 ux_plane_c                      = ux_plane
 uz_plane_c                      = uz_plane

 ;Overlay horizontal vector for non hydrometeor regions
 ux_plane                        = mask(ux_plane,qx_plane.ge.0.1,False) 
 uz_plane                        = 0.
 rescsv@vcRefAnnoParallelPosF    = 0.1
 rescsv@vcRefAnnoString2         = "UX"
 dumvector                       = gsn_csm_vector(wks, ux_plane,uz_plane,rescsv) 
 overlay(plot(1),dumvector)

 ;Overlay vertical vectors for hydromteor regions
 uz_plane_c                      = mask(uz_plane_c,qx_plane.ge.0.1,True)
 ux_plane_c                      = 0.
 rescsv@vcRefMagnitudeF          = 2. ;5.
 rescsv@vcRefAnnoParallelPosF    = 0.2
 rescsv@vcRefAnnoString2         = "UZ"
 dumvector                       = gsn_csm_vector(wks, ux_plane_c,uz_plane_c,rescsv) 
 overlay(plot(1),dumvector)

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

 pmres@gsMarkerSizeF   = 0.045
 dumpm = gsn_add_polymarker(wks,plot(1),lon_plane(sloc),hgt_coord(0),pmres) ;

 ;Add annotation
 txres@txFontHeightF = 0.026
 txres@txBackgroundFillColor = "white"
 gsn_text_ndc(wks,"c)",0.35,0.5,txres) ;

;-----------------------------------------------------------------------------------------
; Cross-Section Plot Resources
;-----------------------------------------------------------------------------------------
 ;Viewport
 resV@vpXF                          = 0.15 ;0.1
 resV@vpYF                          = 0.45
 resV@vpHeightF                     = 0.35
 resV@vpWidthF                      = 0.22

 resV@gsnLeftString                 = "Vertical Profile"
 resV@gsnRightString                = "" 

 resV@xyLineThicknessF              = 3.

 resV@tiYAxisString                 = "Height (km)"
 resV@tiXAxisString                 = "g/kg"

 resV@tmXBMode                      = "Explicit"
 resV@tmXBValues                    = ispan(0,10,5)*0.10
 resV@tmXBLabels                    = ispan(0,10,5)*0.10 
 resV@trXMinF                       = qminval ;0.0
 resV@trXMaxF                       = 1.; 5.;9. 
; resV@trXLog                        = True
 resV@trYMaxF                       = max(ht_v)
 plot(2)                            = gsn_csm_xy (wks,qc_v,ht_v,resV) ; create plot

 resV@xyLineColor                   = "purple"
 dumd                               = gsn_csm_xy (wks,qi_v,ht_v,resV) 
 overlay(plot(2),dumd)
  
 resV@xyLineColor                   = "red"
 dumd                               = gsn_csm_xy (wks,qg_v,ht_v,resV)
 overlay(plot(2),dumd)
 resV@xyLineColor                   = "orange"
 dumd                               = gsn_csm_xy (wks,qs_v,ht_v,resV)
 overlay(plot(2),dumd)

 resV@xyLineColor                   = "magenta"
 dumd                               = gsn_csm_xy (wks,qh_v,ht_v,resV)
 overlay(plot(2),dumd)
 
 resV@xyLineColor                   = "blue"
 dumd                               = gsn_csm_xy (wks,qr_v,ht_v,resV)
 overlay(plot(2),dumd)

;
;Add legend
 lgres                    = True
 lgres@lgLineColors       = (/"black","purple","red","orange","magenta","blue"/)
 lgres@lgLineThicknessF   = resV@xyLineThicknessF  
 lgres@lgDashIndexes      = (/0,0,0,0,0,0/)
 lgres@lgItemType         = "Lines"      
 lgres@lgLabelFontHeightF = .1 
 lgres@vpWidthF           = 0.1       
 lgres@vpHeightF          = 0.18          
 lgres@lgPerimOn          = False 
 lbid = gsn_create_legend(wks,6,(/" qc"," qi"," qg"," qs"," qh"," qr"/),lgres)

 amres = True
 amres@amParallelPosF   = 0.2
 amres@amOrthogonalPosF = -0.15
 annoid1 = gsn_add_annotation(plot(2),lbid,amres)

 ;Add annotation
 txres@txFontHeightF = 0.026
 txres@txBackgroundFillColor = -1 ;"white"
 gsn_text_ndc(wks,"b)",0.07,0.5,txres) ;

 delete([/qc_plane,qi_plane,qh_plane,qg_plane,qs_plane,qr_plane,tc_plane,ux_plane,uz_plane,resV,rescs,rescsv,txres/])

;Draw panel plots
 psres = True                                                               
 maximize_output(wks,psres)

 delete([/plot(0),plot(1),res/])

 qx_ptime(fnctr,:,:) = qx_p
 end do     ;;BIG LOOP OVER FILES

 ;;Write DEBUG OUTPUT FOR KML FILE
 system("rm test.nc")
 fout = addfile("test.nc","c")
 fout->QX = qx_ptime
 fout->lat = lat_c
 fout->lon = lon_c
 fout->topo = topo_c
end
