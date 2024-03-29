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

 ;;MODEL OUTPUT
 pname = "EMPA2_PROM1";
 inst  = 6
 lnd   = "HET"
 itype = "ART"
 diri  = "/p/largedata/hbn33/shrestha1/"
 ;diri        = "/daten04/pshrestha/shared-data/bonnRadar_20150705/"

 itype_gscp  = lnd + itype
 sfold       = "runE_" + inst + "_" + itype_gscp
 expname     = lnd + "_" + pname
 subevent    = expname + "/" + itype_gscp + "/" + sfold

 diri_cos    = diri + "/" + subevent + "/cosout/"
 diri_cosI   = diri + "/" + subevent + "/cosout/ivr/"

 inset       = False
 
 ;Bonn Radar
 glon_id = (/7.071663, 7.00/) ;
 glat_id = (/50.73052, 50.6/) ;
 gmarker = (/10,9/)
 varname = "TOT_PREC" ;"HAIL_GSP" ;"TOT_PREC" ; 
 varminmax = (/0, 30, 5/) 

;-----------------------------------------------------------------------------------------
; SELECT FILES and EXTENTS for Plots 
;-----------------------------------------------------------------------------------------
 ;dext          = 1.5
 ;mpcorner      = (/glat_id(0)-dext,glon_id(0)-dext,glat_id(0)+dext,glon_id(0)+dext/)
 nbdl          = 3 + 17  ;exclude boundaries and relaxation zone

;-----------------------------------------------------------------------------------------
; READ INVARIANT FIELDS
;-------------------------------------------------------------------------------------
 nc_cosI = addfile(diri_cosI +"lfff00000000c.nc","r")
 lat = nc_cosI->lat
 lon = nc_cosI->lon
 ndim = dimsizes(lat)

 topo = nc_cosI->HSURF(0,nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1)
 lat2d = lat(nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1)
 lon2d = lon(nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1)
 topo@lat2d = lat2d
 topo@lon2d = lon2d

 ndimc      = dimsizes(lat2d)

 mpcorner      = (/lat2d(0,0),lon2d(0,0),lat2d(ndimc(0)-1,0),lon2d(0,ndimc(1)-1)/)
 mpcorner      = (/50.5, 6., 51.5, 8./)
;-----------------------------------------------------------------------------------------
; READ NC FILES 
;-----------------------------------------------------------------------------------------
 fils_cos  = systemfunc("ls " + diri_cos + "lf*")
 nf_cos   = dimsizes(fils_cos)

;Read specified variables
 filindex = ispan(1,nf_cos-1,1)
 nf_ctr = dimsizes(filindex)

 varout = new((/nf_ctr, ndimc(0),ndimc(1)/),float)
 str_time = new(nf_ctr, "string")
 xtime    = new(nf_ctr,"float")

 varout@lat2d = lat2d
 varout@lon2d = lon2d
;
 nm = getind_latlon2d (lat,lon, glat_id, glon_id)
 iy = new(dimsizes(glat_id),"integer")
 ix = new(dimsizes(glat_id),"integer")
 do ii = 0, dimsizes(glat_id)-1
  iy(ii) = nm(ii,0)
  ix(ii) = nm(ii,1)
 end do
;
 
 do fnctr = 0, nf_ctr-1   ;;BIG LOOP OVER FILES 
   
   fn = filindex(fnctr)
   nc_cos = addfile(fils_cos(fn),"r")
   nc_cos0 = addfile(fils_cos(fn-1),"r")
   print(" " + fils_cos(fn))
   time     = nc_cos->time
   dtime = calendar_decode2(time,-5)
   year  = dtime(0,0)
   month = dtime(0,1)
   day   = dtime(0,2)
   hh    = dtime(0,3)
   mm    = dtime(0,4)
   ss    = dtime(0,5)

   str_time(fnctr)  = sprinti("%0.4i",year) + "-" + sprinti("%0.2i",month) + \
               "-" + sprinti("%0.2i",day) + " " + sprinti("%0.2i",hh) + \
                sprinti("%0.2i",mm) + " UTC" ;sprinti("%0.2i",ss) +  " UTC"

   xtime(fnctr) = hh*1.0
;
   print(fils_cos(fn) + " " +str_time(fnctr))
   varout(fnctr,:,:)  = (nc_cos->$varname$(0,nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1)) \
                        - (nc_cos0->$varname$(0,nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1))
                
;-----------------------------------------------------------------------------------------
; Plotting Resource
;-----------------------------------------------------------------------------------------
 fname_out = "Figure_rain_"+ expname + "_"+sprinti("%0.2i", fn)
 fncps_no = ":F21:"
 fncps_no = 21
 fncps_hgt = 0.028

 wks_type                 = "png"
 wks = gsn_open_wks(wks_type,fname_out)
 gsn_define_colormap(wks,"precip2_17lev")

 plot = new(2,"graphic")

;Resources
 res                            = True                 ; plot mods desired
 res@gsnDraw                    = False
 res@gsnFrame                   = False
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
 res@tiXAxisString              = " Longitude"
 res@tiYAxisString              = " Latitute"

 xyres                           = res                  ; for time-series plot

 ;gsn Resources
 res@gsnAddCyclic               = False
 res@gsnLeftStringFontHeightF   = fncps_hgt
 res@gsnRightStringFontHeightF  = fncps_hgt
 res@gsnLeftString              = str_time(fnctr)
 res@gsnMaximize                = True


 res@pmTickMarkDisplayMode      = "Always"             ; turn on tick marks

 res@lbLabelAutoStride          = True
 res@lbOrientation              = "vertical"
 res@pmLabelBarWidthF           = 0.1
 res@pmLabelBarHeightF          = 0.5


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
;-----------------------------------------------------------------------------------------
; MAIN PLOT SPECIFIC RESOURCES 
;-----------------------------------------------------------------------------------------
 pvarout = varout(fnctr,:,:)
 
 pvarout = where(pvarout.ne.0.0,pvarout,pvarout@_FillValue)

 print("Min : " + min(pvarout) + " Max : " + max(pvarout))
 res@gsnRightString             = " mm"
 ;CPS res@cnLevelSelectionMode       = "ExplicitLevels"
 ;CPS res@cnLevels                   = ispan(varminmax(0),varminmax(1),varminmax(2))


 plot(0)                         = gsn_csm_contour_map(wks, pvarout, res)

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
 overlay(plot(0),dumtopo)

 txres                       = True
 txres@gsMarkerSizeF = 12
 txres@gsMarkerThicknessF = 5.
 dumtxt = new(dimsizes(glat_id),"graphic")
 do ii= 0, dimsizes(glat_id)-1
 txres@gsMarkerIndex = gmarker(ii)
 dumtxt(ii) = gsn_add_polymarker(wks,plot(0),glon_id(ii),glat_id(ii),txres)
 end do

 if (inset) then

 xyres@vpYF = 0.86
 xyres@vpXF = 0.57
 xyres@vpHeightF = 0.15
 xyres@vpWidthF = 0.3

 xyres@trXMinF = 0.
 xyres@trXMaxF = 24.
 xyres@trYMinF = varminmax(0) ; 0; -50.
 xyres@trYMaxF = varminmax(1) ;100. ;400.
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
   plot(1) = gsn_csm_xy(wks,xtime,varout(:,iy(ii),ix(ii)),xyres)
 else
  dum = gsn_csm_xy(wks,xtime,varout(:,iy(ii),ix(ii)),xyres)
  overlay(plot(1),dum)
 end if
 end do

 end if

 draw(plot)
 frame(wks)
 delete(restopo)
 delete(res)
 end do     ;;BIG LOOP OVER FILES

 system("rm precip.nc")
 fout = addfile("precip.nc","c")
 fout->PRECIP=varout
end
