;Plot chemical composition of PM2.5 for specified lat/long at different heights
;Parameters nf - file number
;Only species with percentage > thr_perc is used for plotting
load "$prom_root/functions/pm25ART.ncl"

begin
;-----------------------------------------------------------------------------------------
;#1 USER SETTINGS
;-----------------------------------------------------------------------------------------
  pname = "EMPA2_PROM1";
  inst  = 6
  lnd   = "HET"
  itype = "ART"
  diri  = "/p/largedata/hbn33/shrestha1/"

  heights     = (/2000./) ;(/10., 2000., 4000., 7000./)
  snf       = 33 ;18+15    ;midnighti + 15
  sid      = 0     ;used INT2lM-MOZART with ldust = True
  ;;
  thr_perc = 0.1 ;1.0  ;Threshold percentage to plot
  glon_idB = 7.071663 + 13*0.01 ; Bonn Radar
  glat_idB = 50.73052
  nbdl     = 0

;--------------------------------------------------------------------------------
  itype_gscp  = lnd + itype
  sfold       = "runE_" + inst + "_" + itype_gscp
  expname     = lnd + "_" + pname
  subevent    = expname + "/" + itype_gscp + "/" + sfold
  diri_cos    = diri + "/" + subevent + "/a_cosout/"
  diri_cosI   = diri + "/" + subevent + "/a_cosout/ivr/"

;;--------------------------------------------------------------------------------
  ;Identify geo-location
  nc_cosI = addfile(diri_cosI + "/lfff00000000c.nc","r")
  lat     = nc_cosI->lat
  lon     = nc_cosI->lon

  nm      = getind_latlon2d (lat,lon, glat_idB, glon_idB)
  yloc    = nm(0,0)
  xloc    = nm(0,1)

  print("Bonn Radar location " + yloc + " " + xloc)

;--------------------------------------------------------------------------------
  fils_cos = systemfunc("ls " + diri_cos + "*.nc")

  do nf = snf,snf ;0, dimsizes(fils_cos)-1 ;BIG LOOP
  fil_cos  = fils_cos(nf)

;Retrieve PM25
;--------------------------------------------------------------------------------
  opt  = True
  opt@iy = yloc 
  opt@ix = xloc
  opt@favg = False
  opt@fdry = True
  opt@nbdl = nbdl
  opt@f3d  = False

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
     pcRes@txFontHeightF   = 0.02 
     pcRes@gsMarkerSizeF   = 0.015;
     plot    = new(nlevs,"graphic")
     figname = "Figure_pm25_"+nf
     wks     = gsn_open_wks("x11", figname)               ; send graphics to PNG file
   end if
   pcRes@tiMainString    = "PM2.5 at " + heights(ip) + " m : " + pm25_str + " [" + pm25_v@units + "]"
   plot(ip)    = pie_chart(wks, percent_f, names_f, colors_f, pcRes)
   delete([/nind,percent_f,names_f,colors_f/])

 end do

 
 draw(plot)
 frame(wks)

 resP = True
 resP@gsnMaximize         = True
 resP@gsnPanelMainString  = pm25_v@long_name
 ;CPS gsn_panel(wks,plot,(/2,2/),resP)

 end do
end
