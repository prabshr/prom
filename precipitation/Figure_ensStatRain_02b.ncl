;======================================================
; Diagnose spatial rainfall
; 26.03.2012 Created 
; 20.04.2012 Added Flags readcos, readclm, readpfl
;======================================================
; Prabhakar Shrestha
; Transregional Collaborative Research Centre 32 (TR32)
; Z4: Model Development and Maintenance Services
; Meteorological Institute
; University of Bonn
; Auf dem Hügel 20
; D-53121 Bonn
; E-mail: pshrestha@uni-bonn.de
;=======================================================

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"

begin

;-----------------------------------------------------------------------------------------
; USER SETTINGS
;-----------------------------------------------------------------------------------------
 fplot       = True
 inset       = False
 nfs = 5 ;0
 nfe = 9 ;19
 plotRADOLAN = True ;False
 diriR       = "/daten01/z4/database/radola/"  ;RADOLAN data
 diriR       = "/p/home/jusers/shrestha1/juwels/Scripts/prom/precipitation/data_empa2_prom1/" 
 fstat       = 1
 diri        = "./data_empa2_prom" + expid
 prefx       = "accprecip_HETART"

 ;Bonn Radar
 glon_id = (/7.071663/) ;, 7.098308/) ;
 glat_id = (/50.73052/) ;, 50.737497/) ;
 gmarker = (/10,9/)
;;;

 if (fstat.eq.1) then
   cnLevels = (/1,5,10,15,20,25,30,40,50,60/)
 else if (fstat.eq.2) then
   cnLevels = ispan(1,20,2)
 end if
 end if

 do nf  = nfs,nfe ;0, nfils-1                   ;EPS Loop
   fname      = diri + "/" + prefx+ "_"+nf+".nc"
   fout       = addfile(fname,"r")
   if (nf.eq.nfs) then
     fgeo       = addfile(diri + "/georain.nc","r")
     topo       = (/fgeo->topo/)
     lat2d      = fgeo->lat2d
     lon2d      = fgeo->lon2d
     ndim       = dimsizes(topo)
     ny         = ndim(0)
     nx         = ndim(1)
     nt         = nfe-nfs+1 ;nfils
     pvar2d     = new((/nt,ny,nx/),"float")
   end if
   pvar2d(nf-nfs,:,:) = (/fout->ACCRAIN/)
   print(nf + " " + fname + " " + avg(pvar2d(nf-nfs,:,:)) + " " + stddev(pvar2d(nf-nfs,:,:)) )

 end do

 print("Ensemble Stat")
 print(avg(pvar2d) + " " + stddev(pvar2d) )

 if (plotRADOLAN) then
   ;foutname = diriR +  "/accprecip_prom"+expid +"_RADOLAN_0.nc"
   foutname = diriR +  "/accprecip_RADOLAN_0.nc"
   fout = addfile(foutname,"r")
   pvar1 = (/fout->ACCRAIN/)
   pvar1@lat2d = lat2d
   pvar1@lon2d = lon2d
   pvar1@units = "mm"
 end if

 if (fstat.eq.1) then
   pvar          = dim_avg_n(pvar2d(0:nfe-nfs,:,:),0)
   pvar@long_name = "Model" ;Ensemble Averaged Acc. Precip."
   pvar1@long_name = "RADOLAN"
 else if (fstat.eq.2) then
   pvar          = dim_stddev_n(pvar2d,0)
   pvar@long_name = "Ensemble Std. Dev. Acc. Precip."
 end if
 end if

 pvar@lat2d    = lat2d
 pvar@lon2d    = lon2d
 pvar@units     = "mm"
 
 topo@lat2d    = (/lat2d/)
 topo@lon2d    = (/lon2d/)
 mpcorner      = (/lat2d(0,0),lon2d(0,0),lat2d(ndim(0)-1,0),lon2d(0,ndim(1)-1)/)

 if (fplot) then
;-----------------------------------------------------------------------------------------
; Plotting Resource
;-----------------------------------------------------------------------------------------
 ;fname_out = "Figure_rain_"+ itype_gscp + "_" + fstat 
 fname_out = "rain_ensemble_" + fstat + "_" + expid
 fncps_no = ":F21:"
 fncps_no = 21
 fncps_hgt = 0.02 ;8

 wks_type                 = "x11"
 wks = gsn_open_wks(wks_type,fname_out)
 if (fstat.eq.1) then
   gsn_define_colormap(wks,"precip2_17lev")
 else if (fstat.eq.2) then
   gsn_define_colormap(wks,"perc2_9lev")
 end if
 end if

 plot = new(2,"graphic")

;Resources
 res                            = True                 ; plot mods desired
 res@gsnDraw                    = False
 res@gsnFrame                   = False
 ;res@gsnSpreadColor             = True
 restopo                        = res

 ;tm Resources
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
 
 ;ti Resources
 res@tiXAxisFont                = fncps_no
 res@tiYAxisFont                = fncps_no
 res@tiXAxisFontHeightF         = fncps_hgt
 res@tiYAxisFontHeightF         = fncps_hgt
 ;res@tiXAxisString              = " Longitude"
 ;res@tiYAxisString              = " Latitute"

 xyres                           = res                  ; for time-series plot

 ;gsn Resources
 res@gsnAddCyclic               = False
 res@gsnLeftStringFontHeightF   = fncps_hgt
 res@gsnRightStringFontHeightF  = fncps_hgt
 ;res@gsnLeftString              = "'
 ;res@gsnMaximize                = True


 res@pmTickMarkDisplayMode      = "Always"             ; turn on tick marks

 res@lbLabelAutoStride          = True
 res@lbOrientation              = "vertical"
 res@pmLabelBarWidthF           = 0.08
 res@pmLabelBarHeightF          = 0.35


 ; Mapping resource (pollat=39.0,   pollon=-173.0,)s
 res@mpCenterLonF                 = 180.-173.0
 res@mpCenterLatF                 = 90. -39.0 ;

 res@mpLimitMode                  = "Corners"           ; method to zoom
 res@mpLeftCornerLatF             = mpcorner(0)    
 res@mpLeftCornerLonF             = mpcorner(1)   
 res@mpRightCornerLatF            = mpcorner(2) 
 res@mpRightCornerLonF            = mpcorner(3) 

 res@mpOutlineOn                  = True
 res@mpOutlineBoundarySets        = "National"
 res@mpDataSetName                = "Earth..4"    ; database for non-USA divisions
 res@mpDataBaseVersion            = "MediumRes"   ; Medium resolution database
 res@mpGridAndLimbOn              = True
 res@mpNationalLineThicknessF     = 3.0
 res@mpGeophysicalLineThicknessF  = 3.0
 res@mpProvincialLineThicknessF   = 2.0
 res@mpPerimOn                    = True

 res@cnFillOn                   = True
 res@cnFillMode                 = "RasterFill"
 res@cnLinesOn                  = False 
 res@cnInfoLabelOn              = False   ; Turn off annotations
 res@cnLineLabelsOn             = False

 res@gsnRightString             = " mm"
 res@cnLevelSelectionMode       = "ExplicitLevels"
 res@cnLevels                   = cnLevels


 res@vpWidthF   = 0.35
 res@vpHeightF  = 0.35
 res@vpXF       = 0.15
 res@vpYF       = 0.9
 res@lbLabelBarOn = False
 plot(0)                         = gsn_csm_contour_map(wks, pvar, res)

 res@vpXF       = 0.56
 res@tmYLOn     = False
 res@lbLabelBarOn = True
 plot(1)                         = gsn_csm_contour_map(wks, pvar1, res)

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
 ;restopo@cnFillDrawOrder        = "PreDraw"

 restopo@lbLabelBarOn           = False

 restopo@gsnLeftString          = ""
 restopo@gsnRightString         = ""
 dumtopo                        = gsn_csm_contour(wks,topo, restopo)
 dumtopo1                        = gsn_csm_contour(wks,topo, restopo)
; overlay(plot(0),dumtopo)
; overlay(plot(1),dumtopo1)

 txres                       = True
 txres@gsMarkerSizeF = 12
 txres@gsMarkerThicknessF = 5.
 dumtxt = new(dimsizes(glat_id),"graphic")
 dumtxt1 = new(dimsizes(glat_id),"graphic")
 do ii= 0, dimsizes(glat_id)-1
 txres@gsMarkerIndex = gmarker(ii)
 dumtxt(ii) = gsn_add_polymarker(wks,plot(0),glon_id(ii),glat_id(ii),txres)
 dumtxt1(ii) = gsn_add_polymarker(wks,plot(1),glon_id(ii),glat_id(ii),txres)
 end do

 if (inset) then

 xyres@vpYF = 0.86
 xyres@vpXF = 0.57
 xyres@vpHeightF = 0.15
 xyres@vpWidthF = 0.3

 xyres@trXMinF = 0.
 xyres@trXMaxF = 24.
 xyres@trYMinF = min(cnLevels) ;(0) ; 0; -50.
 xyres@trYMaxF = max(cnLevels) ;
 xyres@tiYAxisString = " "
 xyres@tiXAxisString = ""

 xyres@xyMarkLineMode = "MarkLines"
 xyres@xyMarkerSizeF  = 12.
 xyres@xyMarkerThicknessF = 3.0
 xyres@xyLineThicknessF = 5.0
 xyres@xyLineColor = "white"
 do ii = 0, dimsizes(glat_id)-1
 xyres@xyMarker = gmarker(ii)
 if (ii.eq.0) then
   plot(2) = gsn_csm_xy(wks,xtime,varout(:,iy(ii),ix(ii)),xyres)
 else
  dum = gsn_csm_xy(wks,xtime,varout(:,iy(ii),ix(ii)),xyres)
  overlay(plot(2),dum)
 end if
 end do

 end if

 ;Bonn Radar Extent
;-----------------------------------------------------------------------
 ncirc    = 100
 circ_lat = new(ncirc,float)               ; Create arrays to hold circle.
 circ_lon = new(ncirc,float)
 raddeg   = 0.91 ;
 nggcog(glat_id(0),glon_id(0),raddeg,circ_lat,circ_lon)

 ;print(circ_lat + " " + circ_lon)
 lnres                  = True
 lnres@gsLineColor      = "Brown"
 lnres@gsLineThicknessF = 6.0
 ;line_id = gsn_add_polyline(wks, plot(0), circ_lon, circ_lat, lnres)
 ;line_id1 = gsn_add_polyline(wks, plot(1), circ_lon, circ_lat, lnres)


 draw(plot)
 frame(wks)

 delete(restopo)
 delete(res)
 end if     ; fplot
end
