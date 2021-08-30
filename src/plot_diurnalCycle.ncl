load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"

begin
  itype_gscp = (/"HET2483","TGR2483","CUR2483"/)
  inst       = 5

  ndim       = dimsizes(itype_gscp) 
  do nf  = 0, ndim-1
    fname   = "./data/flx_" + itype_gscp(nf) + "_" + inst + ".nc"
    nc_fil = addfile(fname,"r")
    if (nf.eq.0) then
      time   = nc_fil->time
      data   = new((/2,ndim,dimsizes(time)/),"float")
    end if
    data(0,nf,:) = (/nc_fil->H/)
    data(1,nf,:) = (/nc_fil->LE/) 
  end do

;-----------------------------------------------------------------------------------------
; Plotting Resource
;-----------------------------------------------------------------------------------------
  wks_type = "x11"
  fname_out = "Figure_flx_" + inst
  wks = gsn_open_wks(wks_type,fname_out)

  fncps_no = ":F21:"
  fncps_no = 21
  fncps_hgt = 0.02

  res                            = True                 ; plot mods desired
 res@gsnDraw                    = False
 res@gsnFrame                   = False
 res@gsnLeftString              = " ";str_time
 res@gsnRightString             = " " ;varp@standard_name + " [" + varp@units+"]"
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

 res@xyLineThicknessF              = 3.
 res@xyDashPattern              = 0
 res@xyLineColors               = (/"orange","forestgreen","red"/)
 res@tiYAxisString                 = "Flux (W/m2)"
 res@tiXAxisString                 = "Time (UTC)" 

 plot                              = gsn_csm_xy (wks,time,data(1,:,:),res) ; create plot

 res@xyDashPattern              = 1

 dum                              = gsn_csm_xy (wks,time,data(0,:,:),res)

 overlay(plot,dum)

 draw(plot)
 frame(wks)



end
