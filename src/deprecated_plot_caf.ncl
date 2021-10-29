load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"
;Flag - fenhance
;Oct 6 2021 - Figure revised restick set to %M from %m

function enhance_plot_model(wks,ploth,resd,fn,fns,fne,sens,ensv,time,data,opt)
begin
 print("Model fn = " + fn + " fns = " + fns + " fne = " + fne)

 resh    = resd
 ;Ensembles for CFAD comparison
 resh@xyLineColor    = opt@color1
 resh@xyDashPattern  = opt@dash1
 dumc0 = gsn_csm_xy (wks,time,data(ensv,:),resh) ; create plot
 overlay(ploth,dumc0)
 resh@xyDashPattern  = opt@dash2
 dumc1 = gsn_csm_xy (wks,time(fns:fne),data(ensv,fns:fne),resh) ; create plot
 overlay(ploth,dumc1)

 ;Member for polarimetric comparison
 resh@xyDashPattern  = opt@dash1
 resh@xyLineColor    = opt@color2
 dump0 = gsn_csm_xy (wks,time,data(sens,:),resh) ; create plot
 overlay(ploth,dump0)
 resh@xyDashPattern  = opt@dash2
 dump1 = gsn_csm_xy (wks,time(fns:fne),data(sens,fns:fne),resh) ; create plot
 overlay(ploth,dump1)
 ;
 if (opt@draw .eq.0) then
 thr                            = 0.02
 dum1                           = gsn_csm_xy (wks, (/time(fns),time(fns)/), (/data(sens,fns)-thr,data(sens,fns)+thr/),resh)
 dum2                           = gsn_csm_xy (wks, (/time(fne),time(fne)/), (/data(sens,fne)-thr,data(sens,fne)+thr/),resh)
 overlay(ploth,dum1)
 overlay(ploth,dum2)
 ;
 resh@xyMarkLineMode            = "Markers"
 resh@xyMarkers                 = opt@marker 
 resh@xyMarkerSizeF             = opt@msize
 resh@xyMarkerThicknessF        = opt@mthick
 resh@xyMarkerColor             = opt@color2
 dum3                           = gsn_csm_xy (wks, (/time(fn),time(fn+1)/), (/data(sens,fn),data@_FillValue/),resh)
 overlay(ploth,dum3)
 end if

 return(0)
end

function enhance_plot_obs(wks,ploth,resd,fn,fns,fne,time,data,opt)
begin
  print("Obs fn = " + fn + " fns = " + fns + " fne = " + fne)
  resh = resd
  ;CFAD comparison time
  resh@xyLineColor               = opt@color3
  resh@xyDashPattern             = opt@dash2
  dum0                           = gsn_csm_xy (wks,time(fns:fne),data(fns:fne),resh)
  overlay(ploth,dum0)
  ;
  thr                            = 0.02
  dum1                           = gsn_csm_xy (wks, (/time(fns),time(fns)/), (/data(fns)-thr,data(fns)+thr/),resh)
  dum2                           = gsn_csm_xy (wks, (/time(fne),time(fne)/), (/data(fne)-thr,data(fne)+thr/),resh)
  overlay(ploth,dum1)
  overlay(ploth,dum2)
  ;
  resh@xyMarkLineMode            = "Markers"
  resh@xyMarkers                 = opt@marker 
  resh@xyMarkerSizeF             = opt@msize
  resh@xyMarkerThicknessF        = opt@mthick
  resh@xyMarkerColor             = opt@color3
  dum3                           = gsn_csm_xy (wks, (/time(fn),time(fn+1)/), (/data(fn),data@_FillValue/),resh)
  overlay(ploth,dum3)
  delete(resh)
  return(0)
end

function data_extract(diri,nens,st,nt)
  begin
  print("Extracting "  + diri)
  do nf = 0, nens-1
    fname  = diri + "/CAF_HET2483_" + nf + ".nc"
    nc_fil = addfile(fname,"r")
    ;print(nf + " " + fname)
    if (nf.eq.0) then
       temp = nc_fil->Band1
       ndim = dimsizes(temp)
       time = st + ispan(0,ndim-1,1)*nt/60.    ;5-15 minute interval data
       data  = new((/nens,ndim/),"float")
    end if
    data(nf,:) = nc_fil->Band1    ;CAF [-]
  end do
  time@units = "hours since 0000-01-01 00:00:00"
  data!1  = "time"
  data&time = time
  return(data)
end

;--------------------------------------------------
begin

  diri1 = "$PROM/scratch/data_prom1/"
  diri2 = "$PROM/scratch/data_prom2/"
  diri3 = "$PROM/scratch/data_prom3/"
  nens = 20
  fenhance = True ;False


; Read observations PROM1
  nc_obsB1 = addfile(diri1 + "/CAF_BoXPol_20150705.nc","r")
  dataoB1  = nc_obsB1->Band1
  timeoB1  = dataoB1&time
  nc_obsJ1 = addfile(diri1 + "/CAF_JuXPol_20150705.nc","r")
  dataoJ1  = nc_obsJ1->Band1
  timeoJ1  = dataoJ1&time

; Read observations PROM2
  nc_obsB2 = addfile(diri2 + "/CAF_BoXPol_20160513.nc","r")
  dataoB2  = nc_obsB2->Band1
  timeoB2  = dataoB2&time

; Read observations PROM3
  nc_obsB3 = addfile(diri3 + "/CAF_BoXPol_20170706.nc","r")
  dataoB3 = nc_obsB3->Band1
  timeoB3  = dataoB3&time
  nc_obsJ3 = addfile(diri3 + "/CAF_JuXPol_20170706.nc","r")
  dataoJ3  = nc_obsJ3->Band1
  timeoJ3  = dataoJ3&time

   
; Read processed synthetic data PROM1,2,3, st= start times, 14, 9.59, 11 UTC, nt = 5 min interval
  data1 = data_extract(diri1 ,nens,14,5)
  data2 = data_extract(diri2 ,nens,9.59,5)
  data3 = data_extract(diri3 ,nens,11,15)

; Crop data at later stages of storm, where CAF might be affected due to border issue


; Some staistics
  qStd1= 100.*dim_stddev_n(data1,0)/dim_avg_n(data1,0)
  qStd2= 100.*dim_stddev_n(data2,0)/dim_avg_n(data1,0)
  qStd3= 100.*dim_stddev_n(data3,0)/dim_avg_n(data1,0) 

  print(dim_median(qStd1) + " " + dim_median(qStd2) + " " + dim_median(qStd3))
;-----------------------------------------------------------------------------------------
; Plotting Resource
;-----------------------------------------------------------------------------------------
 plot = new(3,"graphic")
 wks_type = "png"
 fname_out = "Figure_caf_revised"
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
 res@trYMaxF                    = 0.60

 res@xyLineThicknessF           = 3.
 res@xyLineColor                = "gray"
; res@xyLineColors               = (/"aquamarine4","aquamarine3","aquamarine2","aquamarine1","aquamarine", \
;                                    "dodgerblue4","dodgerblue3","dodgerblue2","dodgerblue1","dodgerblue", \ 
;                                    "goldenrod4", "goldenrod3","goldenrod2","goldenrod1","goldenrod", \
;                                    "gray35", "gray29", "gray19", "gray10", "gray0"/)
 res@xyDashPattern              = 15 

 resd = res


 res@tiYAxisString              = "Convective Area Fraction [-]"
 res@tiXAxisString              = "Time (UTC)" 

 restick = True
 restick@ttmFormat = "%H:%M"
 restick@ttmAxis = "XB"


;PROM1 ------------
 res@trXMinF = min(data1&time)
 res@trXMaxF = max(data1&time)
 time_axis_labels(data1&time,res,restick)

 plot(0)                        = gsn_csm_xy (wks,data1&time,data1,res) ; create plot

 resd@xyLineColor               = "red"
 dum0                           = gsn_csm_xy (wks,timeoB1,dataoB1,resd)
 overlay(plot(0),dum0)
 resd@xyLineColor               = "blue"
 dum1                           = gsn_csm_xy (wks,timeoJ1,dataoJ1,resd)
 overlay(plot(0),dum1)


;PROM2 ------------
 res@trXMinF = min(data2&time)
 res@trXMaxF = max(data2&time)
 time_axis_labels(data2&time,res,restick)
 plot(1)                        = gsn_csm_xy (wks,data2&time,data2,res)

 resd@xyLineColor               = "red"
 dumb0                           = gsn_csm_xy (wks,timeoB2,dataoB2,resd)
 overlay(plot(1),dumb0)

;PROM3 ------------
 res@trXMinF = min(data3&time) 
 res@trXMaxF = max(data3&time)
 time_axis_labels(data3&time,res,restick)
 plot(2)                        = gsn_csm_xy (wks,data3&time,data3,res)

 resd@xyLineColor               = "red"
 dumb3                           = gsn_csm_xy (wks,timeoB3,dataoB3,resd)
 overlay(plot(2),dumb3)
 resd@xyLineColor               = "blue"
 dumj3                           = gsn_csm_xy (wks,timeoJ3,dataoJ3,resd)
 overlay(plot(2),dumj3)

;STEP 3
 if (fenhance) then
;--------------Now we enhance the plots --------------------------------------------
 opt = True
 opt@draw   = 0
 opt@color1 = "grey50"
 opt@color2 = "grey20"
 opt@color3 = "red"
 opt@dash1  = 15
 opt@dash2  = 0
 opt@marker = 6 
 opt@msize  = 0.015
 opt@mthick = 2

;PROM1
;11,15,18/
;Overlay time of model ensemble used for cfads
 fn           = 11 ; 1455 ;18  ;15:30 UTC  ;17
 fns          = 9 ; fn  - 9
 fne          = 18 ; fn  + 0
 sens         = 6
 ensv         = (/5,6,7,8,9/)
 tmp1 = enhance_plot_model(wks,plot(0),resd,fn,fns,fne,sens,ensv,data1&time,data1,opt)

;Overlay time of CFAD comparison for observation
 fn           = 18  ;15:30 UTC  ;17
 fns          = fn  - 9
 fne          = fn  + 0
 tmp1 = enhance_plot_obs(wks,plot(0),resd,fn,fns,fne,timeoB1,dataoB1,opt)

;PROM2
;Overlay time of model ensemble used for cfads
 fn           = 15  ;10:50 UTC  ;17
 fns          = fn  - 4 
 fne          = fn  + 5 
 sens         = 8 
 ensv         = (/5,6,7,8,9/)
 tmp1 = enhance_plot_model(wks,plot(1),resd,fn,fns,fne,sens,ensv,data2&time,data2,opt)

; Obs start at 9:00 while model at 9:35
;Overlay time of CFAD comparison for observation
 fn           = 18  ;10:30 UTC  ;17
 fns          = fn  - 4 
 fne          = fn  + 5
 tmp1 = enhance_plot_obs(wks,plot(1),resd,fn,fns,fne,timeoB2,dataoB2,opt)


;PROM3
;Overlay time of model ensemble used for cfads
 fnv          = (/0, 9,4,18,10,5/)    ;18 is 1530 UTC ;400 UTC
 fnsv         = (/0, 7,3,16, 8,4/)   
 fnev         = (/3,10,6,19,11,7/)  
 sens         = 11 
 delete(ensv)
 ensv         = (/7, 8,10,11, 12,19/)
 topt         = opt
 topt@draw    = -1

 do i = 0, dimsizes(ensv)-1
   if (i.eq.3) then
     vopt = opt
   else
     vopt = topt
   end if
   tmp1 = enhance_plot_model(wks,plot(2),resd,fnv(i),fnsv(i),fnev(i),sens,ensv(i),data3&time,data3,vopt)
 end do

;Overlay time of CFAD comparison for observation
 fn           = 12 ;1400 UTC  ;17
 fns          = fn  - 2 
 fne          = fn  + 1 
 tmp1 = enhance_plot_obs(wks,plot(2),resd,fn,fns,fne,timeoB3,dataoB3,opt)

 end if ;fenhance
;---------Add annotation ----------------------------------------------------------------------

 resP = True
 resP@gsnPanelFigureStrings = (/"a) Case One","b) Case Two","c) Case Three"/)
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

 if (fenhance) then
   gsn_panel(wks,plot,(/2,2/),resP)
 else
   gsn_panel(wks,plot,(/1,3/),resP)
 end if
end
