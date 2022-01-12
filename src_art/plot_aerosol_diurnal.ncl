;Plot diurnal cycle of Nx or Mx of different aerosol modes
; at specified height in 4 different subplots
; Plots spatial average and standard deviation in log scale

load "$prom_root/functions/lndistART.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

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
  pname = "PROM1";
  inst  = 6
  lnd   = "HET"
  itype = "ART"
  diri  = "/p/largedata/hbn33/shrestha1/"

  colors      = (/"lightpink3","forestgreen","steelblue","gray20","gray"/)
  dashs       = (/0,0,0,0,0/)
  nbdl        = 20
  height      = 7000 ;m a.msl
;--------------------------------------------------------------------------------
  itype_gscp  = lnd + itype
  sfold       = "runART_" + inst + "_" + itype_gscp
  expname     = itype_gscp + "_" + pname
  subevent    = expname + "/" + sfold
  diri_cos    = diri + "/" + subevent + "/cosout/"
  diri_cosI   = diri + "/" + subevent + "/cosout/ivr/"

  fils_cos = systemfunc("ls " + diri_cos + "*.nc")

;;--------------------------------------------------------------------------------
  ;Identify model level based on height
  nc_cosI = addfile(diri_cosI + "/lfff00000000c.nc","r")
  tmp     = nc_cosI->lat

  ndim    = dimsizes(tmp)
  hhl     = nc_cosI->HHL(0,:,nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1)
  lat     = nc_cosI->lat(nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1)
  lon     = nc_cosI->lon(nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1)

  delete(ndim)
  ndim    = dimsizes(hhl)
  Z       = 0.5 * (hhl(0:ndim(0)-2,:,:) + hhl(1:ndim(0)-1,:,:))
  Z       = Z(::-1,:,:)

;--------------------------------------------------------------------------------
  ;Geomeric standard deviations of lognormal modes
  sig = (/1.7, 2.0, 1.7, 2.0, 1.4, 2.5, 1.7, 1.6, 1.5, 1.9, 2.0, 1.7/)

  ;Obtain aggregate of 12 log-normal distribution modes
  opt    = True
  opt@iy = 90
  opt@ix = 150
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
    end if
    varo(im,:,:,:) = lndistART(fil_cos,opt)
  end do

 ;Extract sums 
  namemods   = (/"nuc./acc. (pure + mixed)", "soot + coarse", \
                 "dust", "sea-salt"/)

 ;Extract sums at specific height
  varo1 = dim_sum_n(varo(0:3,:,:,:),0)     ;nuc-acc mode ;small
  varo2 = dim_sum_n(varo(4:5,:,:,:),0)     ;soot+coarse
  varo3 = dim_sum_n(varo(6:9,:,:,:),0)     ;dust
  varo4 = dim_sum_n(varo(10:11,:,:,:),0)   ;sea

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

;-----------------------------------------------------------------------------------------
; Plotting Resource
;-----------------------------------------------------------------------------------------

  plot = new(4,"graphic")
  dum  = new(4,"graphic")


  wks_type = "png"
  fname_out = "Figure_airmass_" + height
  wks = gsn_open_wks(wks_type,fname_out)

  fncps_no = ":F21:"
  fncps_no = 21
  fncps_hgt = 0.02

 res                            = True                 ; plot mods desired
 res@gsnDraw                    = False
 res@gsnFrame                   = False
 res@gsnRightString             = ""
 res@gsnMaximize                = False
 res@gsnLeftStringFontHeightF   = fncps_hgt
 res@gsnRightStringFontHeightF  = fncps_hgt
 res@gsnLeftString              = "" ;Z(iz) + " m"

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
 
 res@xyLineThicknessF           = 3.
 res@xyDashPattern              = 0
 res@tiXAxisString              = "Time (UTC)"

 res@vpWidthF       = 0.65
 res@vpHeightF      = 0.4

 res@trYLog         = True
 res@trYMaxF        = 5000.
 res@trYMinF        = 0.1

 gsres                   = True                        ; poly res
 gsres@tfPolyDrawOrder   = "Predraw"                   ; draw this first
 gsres@gsFillColor       = "SlateBlue"                 ; color chosen
 gsres@gsFillOpacityF    = 0.2

 xp = new((/nf_cos*2/),"float")
 yp = new((/nf_cos*2/),"float")
 do ip = 0,3
   do k = 0, nf_cos-1
    yp1           = data(ip,k) + sdata(ip,k)
    yp(k) = logtransform(yp1)

    xp(k)          = time(k)
    xp(2*nf_cos-1-k)  = time(k)
    yp2               = data(ip,k) - sdata(ip,k)
    yp(2*nf_cos-1-k)  = logtransform(yp2)
    ;CPS data(ip,k)  = logtransform(data(ip,k))
  end do
  res@gsnLeftString = namemods(ip)
  plot(ip) = gsn_csm_xy(wks,time,data(ip,:),res)
  ;CPS dum(ip) = gsn_add_polygon (wks,plot(ip),xp,yp,gsres)

 end do

 resP  = True
 resP@gsnPanelMainString = "Height = " + height + " m a. msl" ;nlDp@time 
 gsn_panel(wks,plot,(/2,2/), resP)


end
