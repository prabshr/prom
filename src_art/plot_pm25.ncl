;Plot chemical composition of PM2.5
;Parameters nf - file number
;Only species with percentage > thr_perc is used for plotting
;Pie charts are generated to zxy location
load "$prom_root/functions/pm25ART.ncl"

begin
;User Settings
;--------------------------------------------------------------------------------
  diri = "/p/largedata/hbn33/shrestha1/HET_ART_20150705/runART_6a_HET/cosout/"
;  diri = "/p/scratch/terrsysmp-art/runART_6_HET2483/cosmo_out/"

  fils_cos = systemfunc("ls " + diri + "*.nc")
  ;nf       = 0
  fil_cos  = fils_cos(nf)

  thr_perc = 0.1   ;Threshold percentage to plot
  zyx      = (/20, 150, 150/)   ;Model co-ordinate (Z,Y,X) 

;Retrieve PM25
;--------------------------------------------------------------------------------
  opt  = True
  opt@iy = zyx(1)
  opt@ix = zyx(2)

  pm25_v = pm25ART(fil_cos,opt)
  mhgt   = pm25_v&Z   ;Model height a.g.l

; Compute percentage at specified level
  iz     = zyx(0)

  do iz = 0, 50,10
  print("Height = " + mhgt(iz) + " a.g.l.")

  pm25   = dim_sum_n(pm25_v(:,iz),0)
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

  ;print(percent_f + " " + colors_f(:,0) + " " + names_f)

;Plot
;--------------------------------------------------------------------------------

 pcRes   = True
 pcRes@gsnMaximize     = True

 wks     = gsn_open_wks("x11","pie_chart")               ; send graphics to PNG file

 pcRes@tiMainString    = pm25_v@long_name + ": " + pm25_str + " [" + pm25_v@units + "]"
 plot    = pie_chart(wks, percent_f, names_f, colors_f, pcRes)

 delete([/nind,percent_f,names_f,colors_f/])

 end do

end
