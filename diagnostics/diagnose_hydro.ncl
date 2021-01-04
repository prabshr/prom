;======================================================
; Figure_03 
; Surface Fluxes and Soil States from the hydrological 
; component of the model 
; Average H and LE from 0600 to 1700 hrs
; See how ind_avg is computed ...

; 26.03.2012 Created 
; 20.04.2012 Added Flags readcos, readclm, readpfl
;======================================================
; Prabhakar Shrestha
; E-mail: pshrestha@uni-bonn.de
;=======================================================

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"

begin

;-----------------------------------------------------------------------------------------
; USER SETTINGS
;-----------------------------------------------------------------------------------------

;;MODEL OUTPUT
 itype_gscp  = "HET2483"
 diri        = "/daten04/pshrestha/shared-data/bonnRadar_20150705"
 subevent    = "runE_" + inst + "_" + itype_gscp

 diri_cosI    = diri + "/" + itype_gscp + "/" + subevent + "/cosout/ivr/"
 diri_clm    = diri + "/" + itype_gscp + "/" + subevent + "/clmout/"
 diri_pfl    = diri + "/" + itype_gscp + "/" + subevent + "/pflout/"

 ;Bonn Radar
 glon_id = (/7.071663, 7.00/) ;
 glat_id = (/50.73052, 50.6/) ;
 gmarker = (/10,9/)

 fplot   = True ;False
 ;fnctr   = 15
;-----------------------------------------------------------------------------------------
; READ COSMO GRIDS
;-------------------------------------------------------------------------------------
 nbdl     = 3
 nc_cosI  = addfile(diri_cosI +"lfff00000000c.nc","r")
 lat      = nc_cosI->lat
 lon      = nc_cosI->lon
 ndim     = dimsizes(lat)
 lat2d    = lat(nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1)
 lon2d    = lon(nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1)
 ndimc    = dimsizes(lat2d)

;-----------------------------------------------------------------------------------------
; READ NC FILES 
;-----------------------------------------------------------------------------------------
 fils_clm  = systemfunc("ls " + diri_clm + "clmoas*")
 ncdf_clm  = addfiles(fils_clm,"r")
 ListSetType (ncdf_clm, "cat")

 ; create timestamp for CLM
 tim_clm    = ncdf_clm[:]->time
 shf        = (ncdf_clm[:]->FSH_V) + (ncdf_clm[:]->FSH_G)
 lhf        = (ncdf_clm[:]->FCTR) + (ncdf_clm[:]->FGEV) + (ncdf_clm[:]->FCEV)
 fsds       = (ncdf_clm[:]->FSDS)
 flds       = (ncdf_clm[:]->FLDS)

 hrs        = tointeger(cd_string(tim_clm,"%H")) 
 days       = tointeger(cd_string(tim_clm,"%D")) 

 ind_avg    = ind(hrs.ge.6 .and. hrs.le.17 .and. days.eq.4)

 shf_avg    = dim_avg_n(shf(ind_avg,:,:),0)
 shf_avg@lat2d  = lat2d
 shf_avg@lon2d  = lon2d
 shf_avg@long_name  = "H"

 lhf_avg    = dim_avg_n(lhf(ind_avg,:,:),0)
 lhf_avg@lat2d  = lat2d
 lhf_avg@lon2d  = lon2d 
 lhf_avg@long_name = "LE"

 br_avg    = shf_avg/lhf_avg
 br_avg@lat2d  = lat2d
 br_avg@lon2d  = lon2d
 br_avg@long_name = "BR"

 br1   = avg(br_avg)
 br2   = stddev(br_avg)
 print("Bowen Ratio :" + br1 + " " + br2)

 ;Debug Output file
 fname = "flx_" + itype_gscp + "_" + inst + ".nc"
 system("rm " + fname)
 fout = addfile(fname,"c")
 fout->time = tim_clm
 fout->H    =  dim_avg_n(dim_avg_n(shf,2),1)
 fout->LE   =  dim_avg_n(dim_avg_n(lhf,2),1) 
 ;

 if (fplot) then

;-----------------------------------------------------------------------------------------
; Plotting Resource
;-----------------------------------------------------------------------------------------
;; fname_out = "Figure_pbl_"+case + "_"+sprinti("%0.4i",fn)
 fname_out = "Figure_pbl_"
 fncps_no = ":F21:"
 fncps_no = 21
 fncps_hgt = 0.02
 plot      = new(2,"graphic")

 wks_type                 = "x11"
 wks = gsn_open_wks(wks_type,fname_out)
 gsn_define_colormap(wks,"precip2_17lev")

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
 ;res@tiXAxisString              = " Longitude"
 ;res@tiYAxisString              = " Latitute"

 xyres                           = res                  ; for time-series plot

 ;gsn Resources
 res@gsnAddCyclic               = False
 res@gsnLeftStringFontHeightF   = fncps_hgt
 res@gsnRightStringFontHeightF  = fncps_hgt
 res@gsnLeftString              = " "
 ;res@gsnMaximize                = True

 res@pmTickMarkDisplayMode      = "Always"             ; turn on tick marks

 res@lbLabelAutoStride          = True
 res@lbOrientation              = "vertical"
 res@pmLabelBarWidthF           = 0.08
 res@pmLabelBarHeightF          = 0.35


 ; mp Resources Mapping resource (pollat=39.0,   pollon=-173.0,)s
 res@mpCenterLonF                 = 180.-173.0
 res@mpCenterLatF                 = 90. -39.0 ;
 res@mpLimitMode                  = "Corners"           ; method to zoom
 res@mpLeftCornerLatF             = lat2d(0,0)
 res@mpLeftCornerLonF             = lon2d(0,0)
 res@mpRightCornerLatF            = lat2d(ndimc(0)-1,0)
 res@mpRightCornerLonF            = lon2d(0,ndimc(1)-1)

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

 res@gsnRightString             = " "
 res@cnLevelSelectionMode       = "ExplicitLevels"
 res@cnLevels                   = ispan(0,30,2)*0.1 ;ispan(-50,400,50)*1.0 
 ;ispan(varminmax(0),varminmax(1),varminmax(2))
  ;0,100,10) ;ispan(0,500,20)*1.0 - 100.

 res@vpWidthF     = 0.35
 res@vpHeightF    = 0.35
 res@vpXF         = 0.15
 res@vpYF         = 0.9
 res@lbLabelBarOn = False

 plot(0)     = gsn_csm_contour_map(wks, br_avg, res)

 res@vpXF         = 0.55
 res@tmYLOn       = False
 res@lbLabelBarOn = True
 plot(1)     = gsn_csm_contour_map(wks, lhf_avg, res)

 txres                       = True
 txres@gsMarkerSizeF = 12
 txres@gsMarkerThicknessF = 5.
 dumtxt = new(dimsizes(glat_id),"graphic")
 do ii= 0, dimsizes(glat_id)-1
   txres@gsMarkerIndex = gmarker(ii)
   dumtxt(ii) = gsn_add_polymarker(wks,plot(0),glon_id(ii),glat_id(ii),txres)
 end do

 draw(plot)
 frame(wks)

 end if
end
