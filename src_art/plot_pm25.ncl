;Plot chemical composition of PM2.5
;Parameters nf - file number
;Only species with percentage > thr_perc is used for plotting
;Pie charts are generated to zxy location
load "$prom_root/functions/pm25ART.ncl"

begin
;User Settings
;--------------------------------------------------------------------------------
;User Settings
;--------------------------------------------------------------------------------
  exps = (/"runART_6a_HET","runART_6b_HET","runART_6b_HET"/)

  heights     = (/10., 2000., 4000., 7000./)
  nf       = 18+15    ;midnighti + 15
  sid      = 0     ;used INT2lM-MOZART with ldust = True
  ;;
  diri = "/p/largedata/hbn33/shrestha1/HET_ART_20150705/" + exps(sid) + "/cosout/"
  fils_cos = systemfunc("ls " + diri + "*.nc")
  fil_cos  = fils_cos(nf)

  thr_perc = 0.1   ;Threshold percentage to plot
  zyx      = (/20, 150, 150/)   ;Model co-ordinate (Z,Y,X) 

;Retrieve PM25
;--------------------------------------------------------------------------------
  opt  = True
  opt@iy = zyx(1)
  opt@ix = zyx(2)
  opt@favg = False
  opt@fdry = True

  pm25_v = pm25ART(fil_cos,opt)
  mhgt   = pm25_v&Z   ;Model height a.g.l

; Compute percentage at specified level
  nlevs       = dimsizes(heights)
  izlevs      = new(nlevs,"integer")
  do iz       = 0, nlevs-1
    izlevs(iz) = closest_val(heights(iz),mhgt)
  end do


  do ip = 0, nlevs -1
    iz = izlevs(ip)
    pm25      = dim_sum_n(pm25_v(:,iz),0)
    pm25_str  = decimalPlaces(pm25 ,2,True) 

    percent = 100.*pm25_v(:,iz)/pm25

    names  = pm25_v@names
    nsp    = dimsizes(names)
    colors = new((/nsp,1/),"string")
    colors(:,0) = pm25_v@colors

    nind   = ind(percent.gt.thr_perc)
    names_f = names(nind)
    percent_f = percent(nind)
    colors_f  = colors(nind,:)

   ;Plot
   ;--------------------------------------------------------------------------------

   if (ip.eq.0) then
     pcRes   = True
     pcRes@gsnDraw = False
     pcRes@gsnFrame = False
     pcRes@gsnMaximize     = False ;True

     plot    = new(nlevs,"graphic")

     wks     = gsn_open_wks("png","Figure_pm25")               ; send graphics to PNG file
   end if
   pcRes@tiMainString    = heights(ip) + " m : " + pm25_str + " [" + pm25_v@units + "]"
   plot(ip)    = pie_chart(wks, percent_f, names_f, colors_f, pcRes)
   delete([/nind,percent_f,names_f,colors_f/])

 end do

 resP = True
 resP@gsnMaximize         = True
 resP@gsnPanelMainString  = pm25_v@long_name
 gsn_panel(wks,plot,(/2,2/),resP)
end
