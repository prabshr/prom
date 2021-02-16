load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"

function data_extract(diri,nens,nt)
  begin
  do nf = 0, nens-1
    fname  = diri + "/CAF_HET2483_" + nf + ".nc"
    nc_fil = addfile(fname,"r")
    if (nf.eq.0) then
       temp = nc_fil->Band1
       ndim = dimsizes(temp)
       time = ispan(0,ndim-1,1)*nt    ;5-15 minute interval data
       data  = new((/nens,ndim/),"float")
    end if
    print(" " + nf)
    data(nf,:) = nc_fil->Band1    ;CAF [-]
  end do
  data!1  = "time"
  data&time = time
  return(data)
end

begin

  diri1 = "$PROM/scratch/data_prom1/"
  diri2 = "$PROM/scratch/data_prom2/"
  diri3 = "$PROM/scratch/data_prom3/"
  nens = 20

; Read processed data
  data1 = data_extract(diri1,nens,5)
  data2 = data_extract(diri2,nens,5)
  data3 = data_extract(diri3,nens,15)

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
 res@xyLineColors               = (/"aquamarine4","aquamarine3","aquamarine2","aquamarine1","aquamarine", \
                                    "dodgerblue4","dodgerblue3","dodgerblue2","dodgerblue1","dodgerblue", \ 
                                    "goldenrod4", "goldenrod3","goldenrod2","goldenrod1","goldenrod", \
                                    "gray35", "gray29", "gray19", "gray10", "gray0"/)
 res@xyDashPattern              = 0
 res@tiYAxisString              = "Convective Area Fraction [-]"
 res@tiXAxisString              = "Time (minutes)" 

 res@trXMaxF = max(data1&time)
 plot(0)                        = gsn_csm_xy (wks,data1&time,data1,res) ; create plot

 res@trXMaxF = max(data2&time)
 plot(1)                        = gsn_csm_xy (wks,data2&time,data2,res)
 
 res@trXMaxF = max(data3&time)
 plot(2)                        = gsn_csm_xy (wks,data3&time,data3,res)
 resP = True
 resP@gsnPanelFigureStrings = (/"a)","b)","c)"/)
 resP@amJust   = "TopLeft"

 gsn_panel(wks,plot,(/2,2/),resP)



end
