;Plot diurnal cycle of Nx or Mx of different aerosol modes
; at specified height in 4 different subplots
; Plots spatial average and standard deviation for 9x9 box centered at iy,ix in log scale

load "$prom_root/functions/lndistART.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"

function logtransform(xval)
begin
 if (xval.ge.0) then
   xt =  1. * log(log(abs(xval))+1.)
 else
   xt = -1. * log(log(abs(xval))+1.)
 end if
 return(xt)
end

begin
;-----------------------------------------------------------------------------------------
;#1 USER SETTINGS
;-----------------------------------------------------------------------------------------
  pname = "EMPA2_PROM1";
  einst = (/7/)
  lnd   = "HET"
  itype = "ART"
  diri  = "/p/largedata/hbn33/shrestha1/"

  colors      = (/"lightpink3","forestgreen","steelblue","gray20","gray"/)
  dashs       = (/0,0,0,0,0/)
  nbdl        = 0   ;Always 0 in this code, for iy,ix reference
  height      = 2000 ;m a.msl
  glon_idB    = 7.071663 + 13*0.01  ; Bonn Radar
  glat_idB    = 50.73052
  fprcs       = True ; Process or read processed data
  namemods   = (/"nuc./acc. (pure + mixed)", "soot + coarse", \
                 "dust", "sea-salt"/)

;--------------------------------------------------------------------------------

  if (fprcs) then       ;Process data extraction

  do ens = 0, dimsizes(einst)-1    ;ENSEMBLE LOOP

  itype_gscp  = lnd + itype
  sfold       = "runE_" + einst(ens) + "_" + itype_gscp
  expname     = lnd + "_" + pname
  subevent    = expname + "/" + itype_gscp + "/" + sfold

  diri_cos    = diri + "/" + subevent + "/aiop_cosout/"
  diri_cosI   = diri + "/" + subevent + "/aiop_cosout/ivr/"

  fils_cos = systemfunc("ls " + diri_cos + "*.nc")

;;--------------------------------------------------------------------------------
  ;Identify model level based on height
  nc_cosI = addfile(diri_cosI + "/lfff00000000c.nc","r")
  tmp     = nc_cosI->lat

  ndim    = dimsizes(tmp)
  hhl     = nc_cosI->HHL(0,:,nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1)
  lat     = nc_cosI->lat(nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1)
  lon     = nc_cosI->lon(nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1)

  ;
  nm      = getind_latlon2d (lat,lon, glat_idB, glon_idB)
  yloc    = nm(0,0)
  xloc    = nm(0,1)

  print("Location " + yloc + " " + xloc)

  delete(ndim)
  ndim    = dimsizes(hhl)
  Zw      = 0.5 * (hhl(0:ndim(0)-2,:,:) + hhl(1:ndim(0)-1,:,:))
  Z       = Zw(::-1,yloc-4:yloc+4,xloc-4:xloc+4)

;--------------------------------------------------------------------------------
  ;Geomeric standard deviations of lognormal modes
  sig = (/1.7, 2.0, 1.7, 2.0, 1.4, 2.5, 1.7, 1.6, 1.5, 1.9, 2.0, 1.7/)

  ;Obtain aggregate of 12 log-normal distribution modes
  opt    = True
  opt@iy = yloc
  opt@ix = xloc
  opt@favg = True
  opt@fout = 1  
  opt@nbdl = nbdl

  smodes = ispan(1,12,1) 

  nf_cos = dimsizes(fils_cos)

  time = ispan(0,nf_cos-1,1) 
  data = new((/4,nf_cos/),"float")
  sdata = data

  do nf = 0, nf_cos-1   ;BIG LOOP

  fil_cos = fils_cos(nf)

  do im = 0, dimsizes(smodes)-1
    opt@smode = smodes(im)
    if (nf.eq.0 .and. im.eq.0) then
      temp = lndistART(fil_cos,opt)
      ndim = dimsizes(temp)
      varo = new((/12,ndim(0),ndim(1),ndim(2)/),"float")
      delete(ndim)
    end if
    varo(im,:,:,:) = lndistART(fil_cos,opt)
  end do

 ;Extract sums at specific height
  varo1 = dim_sum_n(varo(0:3,:,yloc-4:yloc+4,xloc-4:xloc+4),0)     ;nuc-acc mode ;small
  varo2 = dim_sum_n(varo(4:5,:,yloc-4:yloc+4,xloc-4:xloc+4),0)     ;soot+coarse
  varo3 = dim_sum_n(varo(6:8,:,yloc-4:yloc+4,xloc-4:xloc+4),0)     ;dust
  varo4 = dim_sum_n(varo(9:11,:,yloc-4:yloc+4,xloc-4:xloc+4),0)   ;sea

  varoz1 = wrf_user_intrp3d(varo1,Z,"h",height,0., False)
  varoz2 = wrf_user_intrp3d(varo2,Z,"h",height,0., False)
  varoz3 = wrf_user_intrp3d(varo3,Z,"h",height,0., False)
  varoz4 = wrf_user_intrp3d(varo4,Z,"h",height,0., False)

  data(0,nf) = avg(varoz1)
  data(1,nf) = avg(varoz2)
  data(2,nf) = avg(varoz3)
  data(3,nf) = avg(varoz4)

  sdata(0,nf) = stddev(varoz1)
  sdata(1,nf) = stddev(varoz2)
  sdata(2,nf) = stddev(varoz3)
  sdata(3,nf) = stddev(varoz4)
 
 end do ;BIG LOOP

 if (ens.eq.0) then
    edata = new((/4,dimsizes(einst),nf_cos/),"float")
    esdata = edata
  end if

  edata(:,ens,:) = data
  esdata(:,ens,:) = sdata
  end do      ;ENSEMBLE LOOP

  foutname = "art-diurnal-profile_" + pname + "_" + height+".nc"
  system("rm " + foutname)
  fout = addfile(foutname,"c")
  fout->Band1 = (/edata/)
  fout->Band2 = (/esdata/)

  else
    foutname = "art-diurnal-profile_" + pname + "_" + height+".nc"
    fout = addfile("./data/" + foutname,"r")
    edata = fout->Band1 
   esdata = fout->Band2 
    ndim  = dimsizes(edata)
   nf_cos = ndim(2)
   time   = ispan(0,nf_cos-1,1) 
   time@units = "hours since 2015-07-04 06:00:00" 
  end if ;fprcs = True

;-----------------------------------------------------------------------------------------
; Plotting Resource
;-----------------------------------------------------------------------------------------

  plot = new(4,"graphic")
  dum  = new(4,"graphic")


  wks_type = "x11"
  fname_out = "Figure_aerosol_diurnal_" + height
  wks = gsn_open_wks(wks_type,fname_out)

  fncps_no = ":F21:"
  fncps_no = 21
  fncps_hgt = 0.02

 res                            = True                 ; plot mods desired
 res@gsnDraw                    = False
 res@gsnFrame                   = False
 res@gsnMaximize                = False
 res@gsnLeftStringFontHeightF   = fncps_hgt
 res@gsnRightStringFontHeightF  = fncps_hgt
 res@gsnRightString              = "Height: " + height + " m"

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

 lres = res

 res@tiXAxisFont                = fncps_no
 res@tiYAxisFont                = fncps_no
 res@tiXAxisFontHeightF         = fncps_hgt
 res@tiYAxisFontHeightF         = fncps_hgt
 
 res@xyMarkLineMode             = "MarkLines"
 res@xyMarker                   = 16 
 res@xyDashPattern              = 0
 res@xyLineThicknessF           = 4.
 res@xyLineColor                = "gray"
 res@tiXAxisString              = "Time (UTC)"

 res@vpWidthF       = 0.8
 res@vpHeightF      = 0.4
 res@vpXF           = 0.13

 res@trYLog         = True
 res@trYMaxF        = 3500.
 res@trYMinF        = 10 
 ;res@trXMaX

 gsres                   = True                        ; poly res
 gsres@tfPolyDrawOrder   = "Predraw"                   ; draw this first
 gsres@gsFillColor       = "SlateBlue"                 ; color chosen
 gsres@gsFillOpacityF    = 0.2


 restick = True
 restick@ttmFormat = "%h"
 restick@ttmAxis = "XB"
 restick@ttmNumTicks = 5
 restick@ttmValues = (/(/2015,7,4,0,0,0/), \
                       (/2015,7,4,6,0,0/), \
                       (/2015,7,4,12,0,0/), \
                       (/2015,7,4,18,0,0/), \
                       (/2015,7,5,0,0,0/), \
                       (/2015,7,5,6,0,0/), \
                       (/2015,7,5,12,0,0/), \
                       (/2015,7,5,18,0,0/), \
                       (/2015,7,6,0,0,0/)/)

 delete(restick@ttmValues)
 do ip = 0, 3
   res@gsnLeftString = namemods(ip)
   res@tiYAxisString = "N~B~na~N~ [cm~S~-3~N~]"
   time_axis_labels(time,res,restick)
   plot(ip) = gsn_csm_xy(wks,time,edata(ip,:,:),res)

 end do

 ip = 0
 ; Draw line
 lres@gsnRightString =""
 nt  = 32
 xx  = (/time(nt),time(nt)/)
 yy  = (/res@trYMinF, res@trYMaxF/)
 lres@xyLineColor = "blue"
 lres@xyLineThicknessF = 2.
 duml1 = gsn_csm_xy(wks,xx,yy,lres)
 overlay(plot(ip),duml1)

 nt  = 33
 xx  = (/time(nt),time(nt)/)
 yy  = (/res@trYMinF, res@trYMaxF/)
 lres@xyLineColor = "red"
 lres@xyLineThicknessF = 2.
 duml1 = gsn_csm_xy(wks,xx,yy,lres)
 overlay(plot(ip),duml1)

 draw(plot(ip))
 frame(wks)

 ;resP  = True
 ;resP@gsnPanelMainString = "Height = " + height + " m a. msl" ;nlDp@time 
 ;gsn_panel(wks,plot,(/2,2/), resP)


end
