load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$prom_root/functions/phirot2phi_nclf.ncl"
load "$prom_root/functions/shapefile_utils.ncl"
load "$prom_root/functions/cosmo_time_nclf.ncl"
load "$prom_root/functions/extract_atmwater_nclf.ncl"
load "$prom_root/functions/vtc_nclf.ncl"

;Plot either AP tracers or vertically integrated water

begin

;STEP 1 User Settings
;-----------------------------------------------------------------------------------------

 pname = "EMPA2_PROM1";
 inst  = 6
 lnd   = "HET"
 itype = "ART"
 diri  = "/p/largedata/hbn33/shrestha1/"

 fns   = 6 
 fne   =fns 
 ixs   = 93 + 65
 nbdl  = 30    ;Cropping data extent

;Tracer
 varnames = (/"NO2","HCHO", "VSOILA0", "TAU"/)
 fvtcs    = (/True,  False, False, False/)
 htrops   = (/10000.0,-999,-999, -999/)  ;Tropospheric height -midlatitudes [m]
 uscales  = (/1E15,2.69E16,1.,1./) ; ; Dobson Units

 varname  = varnames(vid)
 fvtc     = fvtcs(vid)
 htrop    = htrops(vid)
 uscale   = uscales(vid)

;Vertically integrated water
 opt     = 2  ; Extract IWV(0),LWP(1),IWP(2)

 fstream   = False
 shp_fil   = "/daten01/z4/database/TestCases/bonnRadar/shape_map/vect_map.shp"
 fcountry  = False

 ;Bonn Radar
 lonc        = 7.071663; N
 latc        = 50.73052; E

 ;Bonn city
 latb = 50.733334
 lonb = 7.100000

; Read Variables
;-----------------------------------------------------------------------
 itype_gscp  = lnd + itype
 sfold       = "runE_" + inst + "_" + itype_gscp
 expname     = lnd + "_" + pname
 subevent    = expname + "/" + itype_gscp + "/" + sfold
 diri_cos    = diri + "/" + subevent + "/a_cosout/"
 diri_cosI   = diri + "/" + subevent + "/a_cosout/ivr/"


 fil_ivr    = diri_cosI + "lfff00000000c.nc"
 nc_cosI    = addfile(fil_ivr,"r")
 lon2d      = nc_cosI->lon
 lat2d      = nc_cosI->lat
 rotpole    = nc_cosI->rotated_pole
   pole_lon = rotpole@grid_north_pole_longitude
   pole_lat = rotpole@grid_north_pole_latitude

 topo       = nc_cosI->HSURF(0,:,:)
 hhl       = nc_cosI->HHL(0,::-1,:,:)

 ndimhhl   = dimsizes(hhl)
 nz        = ndimhhl(0)-1
 ny        = ndimhhl(1)
 nx        = ndimhhl(2)

 ;Crop topography,height and geo-reference
 hhlc        = hhl(:,nbdl:ny-nbdl-1,nbdl:nx-nbdl-1)
 topoc       = topo(nbdl:ny-nbdl-1,nbdl:nx-nbdl-1)
 topoc@lat2d = lat2d(nbdl:ny-nbdl-1,nbdl:nx-nbdl-1)
 topoc@lon2d = lon2d(nbdl:ny-nbdl-1,nbdl:nx-nbdl-1)
 delete(topoc@coordinates)
 
 fils_cos  = systemfunc("ls " + diri_cos + "lf*")

 do fn = fns,fne    ;BIG LOOP for animation

 nc_cos    = addfile(fils_cos(fn),"r")
 time      = nc_cos->time
 ;
 str_time  = cosmo_time_nclf(time,1)  ; 0 for date
 print("Processing ..." + fn + "  " + str_time)

 if (fvtc) then
   ;Extract vertical tropospheric column 
   tracer3d = nc_cos->$varname$(0,::-1,nbdl:ny-nbdl-1,nbdl:nx-nbdl-1)
   P        = nc_cos->P(0,::-1,nbdl:ny-nbdl-1,nbdl:nx-nbdl-1)
   T        = nc_cos->T(0,::-1,nbdl:ny-nbdl-1,nbdl:nx-nbdl-1)

   tracer   = vtc_nclf(tracer3d,P,T,hhlc,htrop,uscale)
   tracer@lat2d = lat2d(nbdl:ny-nbdl-1,nbdl:nx-nbdl-1)
   tracer@lon2d = lon2d(nbdl:ny-nbdl-1,nbdl:nx-nbdl-1)

 else 

   if (varname.eq."TAU") then
     tracer3d = nc_cos->TAU_AERO_A(0,::-1,nbdl:ny-nbdl-1,nbdl:nx-nbdl-1) + \
                nc_cos->TAU_DUST_A(0,::-1,nbdl:ny-nbdl-1,nbdl:nx-nbdl-1) + \
                nc_cos->TAU_SEAS_A(0,::-1,nbdl:ny-nbdl-1,nbdl:nx-nbdl-1)
     tracer   = dim_sum_n(tracer3d,0)
     tracer@long_name = "AOD at 550 nm"
     tracer@units     = "-"
   else
     tracer3d = nc_cos->$varname$(0,::-1,nbdl:ny-nbdl-1,nbdl:nx-nbdl-1)
     tracer   = (tracer3d(0,:,:))/1.E6 ; surface concentraion cm-3
     tracer3d@units = "cm~S~-3~N~"
   end if

   copy_VarAtts(tracer3d,tracer)
   ;tracer = extract_atmwater_nclf(fils_cos(fn),fil_ivr,nbdl,opt)
   ;Geo-reference
   tracer@lat2d = lat2d(nbdl:ny-nbdl-1,nbdl:nx-nbdl-1)
   tracer@lon2d = lon2d(nbdl:ny-nbdl-1,nbdl:nx-nbdl-1)

 end if

 vmin = min(tracer)
 vmax = max(tracer)

;Need to read shape file and convert rotated grids into geographical grids
;-----------------------------------------------------------------------
 if (fstream) then
 pollon= 180.0 + pole_lon
 pollat= 90. - pole_lat
 print_shapefile_info(shp_fil)
 shp_xy    = addfile(shp_fil,"r")
 geo_dims  = getfilevardimsizes(shp_xy,"x")
 seg_strm  = shp_xy->segments(:,0)
 x_strm    = shp_xy->x
 y_strm    = shp_xy->y

 lon_strm  = new(geo_dims(0),"float")
 lat_strm  = new(geo_dims(0),"float")

 do ist = 0, geo_dims(0)-1
   lon_strm(ist) = phirot2phi_nclf(y_strm(ist), x_strm(ist), pollat, pollon, 0.0, 1)
   lat_strm(ist) = phirot2phi_nclf(y_strm(ist), x_strm(ist), pollat, pollon, 0.0, 0)
 end do

 end if


;PLOT
;-----------------------------------------------------------------------
 fncps      = "~F21~"
 fncps_no   = 21
 fncps_hgt  = 0.021

 tfncps_no  = 22
 tfncps_hgt = fncps_hgt*1.2
 txbgColor  = "white"
 txfnColor  = "black"

 wks_type   = "x11"
 wks        = gsn_open_wks(wks_type,"Figure_" +varname+"_I" + inst + "_" + sprinti("%0.2i",fn))

 res = True

 ;gsn Resources
 ;--------------------
 res@gsnDraw                    = False
 res@gsnFrame                   = False
 res@gsnAddCyclic               = False
 res@gsnLeftString              = str_time
 if (tracer@long_name.eq."NO2") then
   tracer@long_name = "NO~B~2~N~"
 end if
 res@gsnRightString             = tracer@long_name + " [" + tracer@units + "]" 
 res@gsnRightStringFontHeightF= fncps_hgt
 res@gsnLeftStringFontHeightF = fncps_hgt

 ;vp Resources
 ;-----------------------
 res@vpWidthF                   = 0.57
 res@vpHeightF                  = 0.57
 res@vpXF                       = 0.15
 res@vpYF                       = 0.92;5

 mres                           = res
 ;tm Resources
 ;--------------------
 res@tmXBLabelFont              = fncps_no
 res@tmYLLabelFont              = fncps_no
 res@tmXTOn                     = False
 res@tmYROn                     = False
 res@tmYLLabelFontHeightF       = fncps_hgt
 res@tmXBLabelFontHeightF       = fncps_hgt

 ;ti Resources
 ;--------------------
 res@tiXAxisFontHeightF         = fncps_hgt
 res@tiYAxisFontHeightF         = fncps_hgt

 ;lg Resources
 ;--------------------
 res@lgLabelFont               = fncps_no

 ;lb Resources
 ;--------------------
 res@lbLabelAutoStride         = True
 res@lbLabelFont               = fncps_no
 res@lbOrientation             = "Vertical"
 res@lbLabelFontHeightF        = fncps_hgt
 res@lbLabelAutoStride         = True
 res@lbTitleString             = fncps + " "
 res@lbTitleFontHeightF        = 0.015
 ;res@lbAutoManage              = False

 ;tf Resources
 ;--------------------
 res@tfDoNDCOverlay            = True

 ;pm Resources
 ;--------------------
 res@pmTickMarkDisplayMode     = "conditional"
 res@pmLabelBarWidthF          = 0.03
; res@pmLabelBarOrthogonalPosF = .12


 cres = res
; cn Resources
 ;--------------------
 res@cnFillOn                  = False             ; turn on color
 res@cnLinesOn                 = True            ; no contour lines
 res@cnLineLabelsOn            = False            ; no contour labels
 res@cnInfoLabelOn             = False
 res@cnLineThicknessF          = 3.
 res@cnLineDashPattern         = 0
 res@cnLineLabelBackgroundColor  = "Transparent"
 res@cnLineColor               = "black"
 res@cnSmoothingOn             = True
 res@cnSmoothingDistanceF      = .0003
 res@cnLevelSelectionMode    = "ExplicitLevels"
 res@cnLevels                  = (/100./)*1.0

 ;Cropped corners
 mpcorner = (/lat2d(nbdl,nbdl),lon2d(nbdl,nbdl),lat2d(ny-nbdl-1,nx-nbdl-1),lon2d(ny-nbdl-1,nx-nbdl-1)/)

;mp Resources
;--------------------
 ;res@mpProjection           = "CylindricalEquidistant"
 res@mpCenterLonF           = 180.0 + pole_lon
 res@mpCenterLatF           = 90.0 - pole_lat
 
 res@mpLimitMode           = "LatLon"        ;Important
 res@mpMinLatF             = mpcorner(0) 
 res@mpMinLonF             = mpcorner(1)
 res@mpMaxLatF             = mpcorner(2)
 res@mpMaxLonF             = mpcorner(3)

 res@mpFillOn                    = True
 res@mpOceanFillColor            = "White"
 res@mpLandFillColor             = "transparent"
 res@mpGridAndLimbOn             = True
 res@mpDataSetName               = "Earth..4"
 res@mpDataBaseVersion           = "MediumRes"      ; use finer database
 res@mpOutlineOn                 = True
 res@mpOutlineBoundarySets       = "National"
 res@mpProvincialLineThicknessF  = 4.
 res@mpNationalLineDashPattern   = 3
 res@mpNationalLineThicknessF    = 8.
 res@mpGeophysicalLineThicknessF = 8.

 res@mpFillDrawOrder           = "PostDraw"

 plot = gsn_csm_contour_map(wks, topoc, res)

 ;cn Resources
 ;--------------------
 cres@cnFillOn                = True             ; turn on color
 cres@cnFillPalette           = "amwg" ;"sunshine_9lev"  ;"amwg"
 cres@cnLinesOn               = False            ; no contour lines
 cres@cnLineLabelsOn          = False            ; no contour labels
 cres@cnFillMode              = "CellFill"
 cres@cnLevelSelectionMode    = "ExplicitLevels"
 cres@cnFillOpacityF          = 0.5
 cres@gsnLeftString           = ""
 cres@gsnCenterString         = ""
 cres@cnInfoLabelFontColor    = "grey"
 if (vid.eq.0) then
   cres@cnLevels                = (/0.05,0.1,0.3,0.5,1.,3.,5.,10.,13.,15.,20./)

 else if (vid.eq.3) then
   cres@cnLevels                = ispan(1,10,1)*0.1
 else
   cres@cnLevels                = fspan(vmin,vmax,10) ;ispan(1,10,1)*0.1  ;fspan(vmin,vmax,20)
 end if
 end if

 cres@gsnRightString          = ""
 dum_tracer  = gsn_csm_contour(wks,tracer,cres)
 overlay(plot,dum_tracer)

;Add polymarkers and polylines
;-----------------------------------------------------------------------
 txres                        = True    ;City/Country Text
 pxres                        = True    ;City/Country PolyMarker
 lnres                        = True    ;Circle
 pres                         = True    ;ShapeFile

 pxres@gsMarkerColor          = "black" ;"firebrick2"
 pxres@gsMarkerSizeF          = 0.01
 pxres@gsMarkerThicknessF     = 6.0
 pxres@gsMarkerIndex          = 16

 txres@txFont                 = fncps_no
 txres@txFontHeightF          = tfncps_hgt*0.8
; txres@txFontThicknessF       = 10.0        ;creates shadow effect
; txres@txFontColor            = txfnColor
 ;txres@txBackgroundFillColor  = "white"
 ;txres@txPerimSpaceF          = 0.

 fbonn = gsn_add_text(wks,plot,"Bonn",lonb,latb,txres)

;
; Streams
;------------------
 if (fstream) then
 pres@gsLineColor      = "blue" ;steel blue"
 pres@gsLineThicknessF = 6.0
 pres@gsLineOpacityF   = 0.6
 pres@gsSegments       = seg_strm
 poly2 = gsn_add_polyline(wks,plot,lon_strm,lat_strm,pres)

 end if


 txres@txFontHeightF = tfncps_hgt;0.025
 txres@txFont = fncps_no
 txres@txFontColor = "black"
; gsn_text_ndc(wks,figno(figid)+")",0.15,0.97,txres) ;


 if (fcountry) then
 ;country
 deu = "/Users/pshrestha/data/shapefiles/gadm36_DEU_shp/gadm36_DEU_0.shp"
 bel = "/Users/pshrestha/data/shapefiles/gadm36_BEL_shp/gadm36_BEL_0.shp"
 fra  = "/Users/pshrestha/data/shapefiles/gadm36_FRA_shp/gadm36_FRA_0.shp"
 lux  = "/Users/pshrestha/data/shapefiles/gadm36_LUX_shp/gadm36_LUX_0.shp"
 nld = "/Users/pshrestha/data/shapefiles/gadm36_NLD_shp/gadm36_NLD_0.shp"

 ;states
 nlds  = "/Users/pshrestha/data/shapefiles/gadm36_DEU_shp/gadm36_DEU_1.shp"

 end if

;Bonn Radar Extent
;-----------------------------------------------------------------------
 ncirc    = 100
 circ_lat = new(ncirc,float)               ; Create arrays to hold circle.
 circ_lon = new(ncirc,float)

 nggcog(latc,lonc,0.91,circ_lat,circ_lon)

 ;print(circ_lat + " " + circ_lon)
 lnres@gsLineColor      = "Brown"
 lnres@gsLineThicknessF = 6.0
 ;line_id = gsn_add_polyline(wks, plot, circ_lon, circ_lat, lnres)

 ;;

 if (fcountry) then

 ppres             = True
 ;ppres@gsLineDashPattern = 16
 ppres@gsLineThicknessF = 6.0
 ppres@gsLineColor = "gray3"
 poly10 = gsn_add_shapefile_polylines(wks,plot,nlds,ppres)
 ppres@gsLineDashPattern = 3
 ppres@gsLineThicknessF = 10.0
 poly00 = gsn_add_shapefile_polylines(wks,plot,deu,ppres)
 poly01 = gsn_add_shapefile_polylines(wks,plot,bel,ppres)
 poly02 = gsn_add_shapefile_polylines(wks,plot,fra,ppres)
 poly03 = gsn_add_shapefile_polylines(wks,plot,lux,ppres)
 poly04 = gsn_add_shapefile_polylines(wks,plot,nld,ppres)

 end if

 ;Add cross-section

 ppres                    = True
 ppres@gsLineColor        = "gray3"
 ppres@gsLineDashPattern  = 3
 ppres@gsLineThicknessF   = 1.0
 xx   = (/lon2d(0,ixs),lon2d(ny-1,ixs)/)
 yy   = (/lat2d(0,ixs),lat2d(ny-1,ixs)/)
 ;csAA = gsn_add_polyline(wks,plot,xx,yy,ppres)

 if (vid.eq.3) then
   tres  = True
   tres@txFontHeightF = 0.015

   mres = True
   mres@gsMarkerIndex = 7
   mres@gsMarkerSizeF = 10.
   mres@gsMarkerThicknessF = 4.

   stname = "FJZ-JOYCE"
   mloc = (/50.908333,6.413056/)
   mid1 = gsn_add_polymarker(wks,plot,mloc(1),mloc(0),mres)
   lab1 = gsn_add_text(wks,plot,stname,mloc(1),mloc(0)+0.1,tres)
   stname = "MAINZ"
   mloc=(/49.999000,8.300000/)
   mid2 = gsn_add_polymarker(wks,plot,mloc(1),mloc(0),mres)
   lab2 = gsn_add_text(wks,plot,stname,mloc(1),mloc(0)+0.1,tres)
 end if
 ;
 ;Bonn Radar Extent
;-----------------------------------------------------------------------
 ncirc    = 100
 circ_lat = new(ncirc,float)               ; Create arrays to hold circle.
 circ_lon = new(ncirc,float)

 nggcog(latc,lonc,0.91,circ_lat,circ_lon)

 ;print(circ_lat + " " + circ_lon)
 lnres = True
 lnres@gsLineColor      = "Brown"
 lnres@gsLineThicknessF = 6.0
 line_id = gsn_add_polyline(wks, plot, circ_lon, circ_lat, lnres)

 draw(plot)
 frame(wks)
 delete([/res,cres/])
 end do ;BIG LOOP for animation
end
