load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$prom_root/functions/phirot2phi_nclf.ncl"
load "$prom_root/functions/shapefile_utils.ncl"
load "$prom_root/functions/cosmo_time_nclf.ncl"
load "$prom_root/functions/extract_atmwater_nclf.ncl"

;Plot either AP tracers or vertically integrated water

begin

;STEP 1 User Settings
;-----------------------------------------------------------------------------------------
 pid       = 0
 ydates    = (/"20150705","20160513","20170706"/)
 lndcov    = (/"HET","HET","HET"/)
 insts     = (/6, 8, 11/)
 fns       = (/36,34,18/) 
 ixs       = (/93,124,95/) + 65  ;Need to add offset due to clipping of radar
 ;ixs       = (/60,124,95/) + 65 ;Stack emission

 ftracer   = True 
;Tracer
 ilev      = 70
 varname   = "NO"

;Vertically integrated water
 opt     = 2  ; Extract IWV(0),LWP(1),IWP(2)
 nbdl    = 35

 fstream   = False
 shp_fil   = "/daten01/z4/database/TestCases/bonnRadar/shape_map/vect_map.shp"
 fcountry  = False


 diri      = "/p/largedata/hbn33/shrestha1/"
; diri      = "/daten04/pshrestha/shared-data/"

 diri_cos  = diri + "HET_ART_"+ydates(pid)+"/runART_" + insts(pid) +"_" +lndcov(pid)+"/cosout/"
 diri_cosI = diri_cos + "/ivr/"

 ;Bonn Radar
 lonc        = 7.071663; N
 latc        = 50.73052; E

; Read Variables
;-----------------------------------------------------------------------
 fil_ivr    = diri_cosI + "lfff00000000c.nc"
 nc_cosI    = addfile(fil_ivr,"r")
 lon2d      = nc_cosI->lon
 lat2d      = nc_cosI->lat
 rotpole    = nc_cosI->rotated_pole
 topo       = nc_cosI->HSURF(0,:,:)
 
 fils_cos  = systemfunc("ls " + diri_cos + "lf*")

 nc_cos    = addfile(fils_cos(fns(pid)),"r")
 time      = nc_cos->time
 ;
 str_time  = cosmo_time_nclf(time,1)  ; 0 for date
 print("Processing ..." + fns(pid) + "  " + str_time)

 if (ftracer) then
   ;Extract Tracers 
   tracer = nc_cos->$varname$(0,ilev,:,:)
   tracer@lat2d = lat2d
   tracer@lon2d = lon2d
 else
   tracer = extract_atmwater_nclf(fils_cos(fns(pid)),fil_ivr,nbdl,opt)
 end if

 vmin = min(tracer)
 vmax = max(tracer)

;Need to read shape file and convert rotated grids into geographical grids
;-----------------------------------------------------------------------
 if (fstream) then
 pollon= 180.0 + rotpole@grid_north_pole_longitude      
 pollat= 90. - rotpole@grid_north_pole_latitude
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


; Extract parameters
 ndim = dimsizes(lat2d)

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
 wks        = gsn_open_wks(wks_type,"Figure_tracers")

 res = True

 ;gsn Resources
 ;--------------------
 res@gsnDraw                    = False
 res@gsnFrame                   = False
 res@gsnSpreadColors            = True             ; use full color map
 res@gsnAddCyclic               = False
 res@gsnLeftString              = ""
 res@gsnRightString             = tracer@long_name + " [" + tracer@units + "]" 

 ;vp Resources
 ;-----------------------
 res@vpWidthF                   = 0.7
 res@vpHeightF                  = 0.7
 res@vpXF                       = 0.1
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


;mp Resources
;--------------------
 res@mpCenterLonF              = 180.0 + rotpole@grid_north_pole_longitude           ;lon of centre of map projection system
 res@mpCenterLatF              = 90.0 - rotpole@grid_north_pole_latitude    ;lat of centre of map projection system

 res@mpLimitMode               = "LatLon"
 res@mpMinLatF          = lat2d(0,0)
 res@mpMinLonF          = lon2d(0,0)
 res@mpMaxLatF         = lat2d(ndim(0)-1,ndim(1)-1)
 res@mpMaxLonF         = lon2d(ndim(0)-1,ndim(1)-1)

 res@mpFillOn                  = True
 res@mpOceanFillColor          = "White"
 res@mpLandFillColor           = "transparent"
 res@mpGridAndLimbOn           = True
 res@mpDataSetName             = "Earth..4"
 res@mpDataBaseVersion         = "MediumRes"      ; use finer database
 res@mpOutlineOn               = False ;True
 res@mpOutlineBoundarySets     = "National"
 res@mpProvincialLineThicknessF  = 4.
 res@mpNationalLineDashPattern   = 3
 res@mpNationalLineThicknessF    = 8.
 res@mpGeophysicalLineThicknessF = 8.
 res@mpFillDrawOrder           = "PostDraw"

 plot = gsn_csm_contour_map(wks, topo, res)

 ;cn Resources
 ;--------------------
 cres@cnFillOn                = True             ; turn on color
 cres@cnLinesOn               = False            ; no contour lines
 cres@cnLineLabelsOn          = False            ; no contour labels
 cres@cnFillMode              = "CellFill"
 cres@cnLevelSelectionMode    = "ExplicitLevels"
 cres@cnFillOpacityF          = 0.7
 cres@gsnLeftString           = ""
 cres@gsnCenterString         = ""
 cres@cnInfoLabelFontColor    = "grey"
 cres@cnLevels                = fspan(vmin,vmax,50)
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

 txres@txFont                 = tfncps_no
 txres@txFontHeightF          = tfncps_hgt
 txres@txFontThicknessF       = 10.0        ;creates shadow effect
 txres@txFontColor            = txfnColor
 ;txres@txBackgroundFillColor  = "white"
 ;txres@txPerimSpaceF          = 0.
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
 line_id = gsn_add_polyline(wks, plot, circ_lon, circ_lat, lnres)

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
 xx   = (/lon2d(0,ixs(pid)),lon2d(ndim(0)-1,ixs(pid))/)
 yy   = (/lat2d(0,ixs(pid)),lat2d(ndim(0)-1,ixs(pid))/)
 csAA = gsn_add_polyline(wks,plot,xx,yy,ppres)


 ;
 draw(plot)
 frame(wks)
end
