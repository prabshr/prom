load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"

function data_extract(diri,nens,st,nt)
  begin
  do nf = 0, nens-1
    fname  = diri + "/CAF_HET2483_" + nf + ".nc"
    nc_fil = addfile(fname,"r")
    if (nf.eq.0) then
       temp = nc_fil->Band1
       ndim = dimsizes(temp)
       time = st + ispan(0,ndim-1,1)*nt/60.    ;5-15 minute interval data
       data  = new((/nens,ndim/),"float")
    end if
    ;print(" " + nf)
    data(nf,:) = nc_fil->Band1    ;CAF [-]
  end do
  time@units = "hours since 0000-01-01 00:00:00"
  data!1  = "time"
  data&time = time
  return(data)
end

begin

  diri1 = "$PROM/scratch/data_prom1/"
  diri2 = "$PROM/scratch/data_prom2/"
  diri3 = "$PROM/scratch/data_prom3/"
  nens = 20

; Read observations
  nc_obs1 = addfile(diri1 + "/CAF_BoXPol_20150705.nc","r")
  datao1  = nc_obs1->Band1
  timeo1  = datao1&time
  nc_obs2 = addfile(diri1 + "/CAF_JuXPol_20150705.nc","r")
  datao2  = nc_obs2->Band1
  timeo2  = datao2&time
 
; Read processed synthetic data
  data1 = data_extract(diri1,nens,14,5)
  data2 = data_extract(diri2,nens,9,5)
  data3 = data_extract(diri3,nens,11,15)

; Some staistics
  qStd1= 100.*dim_stddev_n(data1,0)/dim_avg_n(data1,0)
  qStd2= 100.*dim_stddev_n(data2,0)/dim_avg_n(data1,0)
  qStd3= 100.*dim_stddev_n(data3,0)/dim_avg_n(data1,0) 

  print(dim_median(qStd1) + " " + dim_median(qStd2) + " " + dim_median(qStd3))
;-----------------------------------------------------------------------------------------
; Plotting Resource
;-----------------------------------------------------------------------------------------
 plot = new(3,"graphic")
 wks_type = "x11"
 fname_out = "Figure_caf"
 wks = gsn_open_wks(wks_type,fname_out)

 fncps_no = ":F21:"
 fncps_no = 21
 fncps_hgt = 0.025

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

 res@trYMinF                    = -0.01
 res@trYMaxF                    = 0.40

 res@xyLineThicknessF           = 3.
 res@xyLineColor                = "gray"
; res@xyLineColors               = (/"aquamarine4","aquamarine3","aquamarine2","aquamarine1","aquamarine", \
;                                    "dodgerblue4","dodgerblue3","dodgerblue2","dodgerblue1","dodgerblue", \ 
;                                    "goldenrod4", "goldenrod3","goldenrod2","goldenrod1","goldenrod", \
;                                    "gray35", "gray29", "gray19", "gray10", "gray0"/)
 res@xyDashPattern              = 1 

 resd = res


 res@tiYAxisString              = "Convective Area Fraction [-]"
 res@tiXAxisString              = "Time (UTC)" 

 res@trXMaxF = max(data1&time)
 restick = True
 restick@ttmFormat = "%H:%m"
 restick@ttmAxis = "XB"
 time_axis_labels(data1&time,res,restick)

 plot(0)                        = gsn_csm_xy (wks,data1&time,data1,res) ; create plot

 resd@xyLineColor               = "red"
 dum0                           = gsn_csm_xy (wks,timeo1,datao1,resd)
 overlay(plot(0),dum0)
 resd@xyLineColor               = "blue"
 dum1                           = gsn_csm_xy (wks,timeo2,datao2,resd)
 overlay(plot(0),dum1)

;Overlay time of CFAD comparison
 fn           = 18  ;15:30 UTC  ;17
 fns          = fn  - 9 
 fne          = fn  + 2
 resd@xyLineColor               = "red"
 resd@xyDashPattern              = 0
 dum00                           = gsn_csm_xy (wks,timeo1(fns:fne),datao1(fns:fne),resd)
 overlay(plot(0),dum00)
 
;Overlay ensemble for cfads
 fn           = 18  ;15:30 UTC  ;17
 fns          = fn  - 9 
 fne          = fn  + 2 
 sens         = 6

; Overlay ensembles for CFAD comparison
 resd@xyLineColor  = (/"grey50"/)
 resd@xyDashPattern              = 1 
 dumc0 = gsn_csm_xy (wks,data1&time,data1(5:9,:),resd) ; create plot
 overlay(plot(0),dumc0)
 resd@xyDashPattern              = 0 
 dumc1 = gsn_csm_xy (wks,data1&time(fns:fne),data1(5:9,fns:fne),resd) ; create plot
 overlay(plot(0),dumc1)

;Overlay ensemble for polarimetric comparison
 resd@xyDashPattern              = 1 
 resd@xyLineColor    = (/"grey20"/)
 dump0 = gsn_csm_xy (wks,data1&time,data1(sens,:),resd) ; create plot
 overlay(plot(0),dump0)
 resd@xyDashPattern              = 0 
 dump1 = gsn_csm_xy (wks,data1&time(fns:fne),data1(sens,fns:fne),resd) ; create plot
 overlay(plot(0),dump1)
 


 res@trXMaxF = max(data2&time)
 time_axis_labels(data2&time,res,restick)
 plot(1)                        = gsn_csm_xy (wks,data2&time,data2,res)
 
 res@trXMaxF = max(data3&time)
 time_axis_labels(data3&time,res,restick)
 plot(2)                        = gsn_csm_xy (wks,data3&time,data3,res)
 resP = True
 resP@gsnPanelFigureStrings = (/"a)","b)","c)"/)
 resP@amJust   = "TopLeft"

;Add legend
 lgres                    = True
 lgres@lgLineColors       = (/"red", "blue", "gray"/)
 lgres@lgLineThicknessF   = res@xyLineThicknessF
 lgres@lgDashIndexes      = (/0,0,0/)
 lgres@lgItemType         = "Lines"
 lgres@lgLabelFontHeightF = 0.30
 lgres@vpWidthF           = 0.14
 lgres@vpHeightF          = 0.1
 lgres@lgPerimOn          = False
 lbid = gsn_create_legend(wks,3,(/"BoXPol","JuXPol","Model"/),lgres)

 amres = True
 amres@amParallelPosF   = 0.2
 amres@amOrthogonalPosF = -0.35
 annoid1 = gsn_add_annotation(plot(0),lbid,amres)

 gsn_panel(wks,plot,(/2,2/),resP)



end
