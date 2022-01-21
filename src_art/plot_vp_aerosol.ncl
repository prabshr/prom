;Plot ensemble vertical profile of aerosols 
;Avgerage vertical profile plots for 9x9 box surrounding lat/lon

load "$prom_root/functions/lndistART.ncl"
load "$prom_root/functions/pm25ART.ncl"

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
  pname = "KNPROM1";
  einst = (/6/) ;(/5,6,7,8,9/)
  lnd   = "HET"
  itype = "ART"
  diri  = "/p/largedata/hbn33/shrestha1/"
  fprcs = True   ;Read data and process

  colors      = (/"lightpink3","forestgreen","steelblue","gray20","gray"/)
  dashs       = (/0,0,0,0,0/)
  nbdl        = 0    ;Keep 0 because of difference in lndist and pm25 reference for iy and ix
  glon_idB    = 7.071663 ; Bonn Radar
  glat_idB    = 50.73052
  ext         = 4        ; 4 + 4 + 1 ~ 10 grid cells averaging
  nf          = 32
;--------------------------------------------------------------------------------
  if (fprcs) then       ;Process data extraction

  do ens = 0, dimsizes(einst)-1    ;ENSEMBLE LOOP

  itype_gscp  = lnd + itype
  sfold       = "runART_" + einst(ens) + "_" + itype_gscp
  expname     = itype_gscp + "_" + pname
  subevent    = expname + "/" + sfold
  diri_cos    = diri + "/" + subevent + "/cosout/"
  diri_cosI   = diri + "/" + subevent + "/cosout/ivr/"

  fils_cos = systemfunc("ls " + diri_cos + "*.nc")

;;--------------------------------------------------------------------------------
  ;Identify geo-location and extract model height
  nc_cosI = addfile(diri_cosI + "/lfff00000000c.nc","r")
  tmp     = nc_cosI->lat

  ndim    = dimsizes(tmp)
  hhl     = nc_cosI->HHL(0,:,:,:)
  lat     = nc_cosI->lat
  lon     = nc_cosI->lon

  delete(ndim)
  ndim    = dimsizes(hhl)
  Z       = 0.5 * (hhl(0:ndim(0)-2,:,:) + hhl(1:ndim(0)-1,:,:))
  Z       = Z(::-1,:,:)
  nz      = ndim(0)-1
  delete(ndim)

  ;
  nm      = getind_latlon2d (lat,lon, glat_idB, glon_idB)
  yloc    = nm(0,0)
  xloc    = nm(0,1)
  height  = Z(:,yloc,xloc)/1000. ;m

  print("Location " + yloc + " " + xloc) 

;--------------------------------------------------------------------------------
  ;Geomeric standard deviations of lognormal modes
  sig = (/1.7, 2.0, 1.7, 2.0, 1.4, 2.5, 1.7, 1.6, 1.5, 1.9, 2.0, 1.7/)

  ;Obtain aggregate of 12 log-normal distribution modes
  opt    = True
  opt@iy = yloc
  opt@ix = xloc
  opt@favg = True  ;Only for fout <1
  opt@fout = 1  
  opt@nbdl = nbdl

  smodes = ispan(1,12,1) 

;Retrieve Ncm3 and PM2.5
;--------------------------------------------------------------------------------
  nf_cos = dimsizes(fils_cos)
  fil_cos = fils_cos(nf)

  ;Ncm3
  do im = 0, dimsizes(smodes)-1
    opt@smode = smodes(im)
    if (im.eq.0) then
      temp = lndistART(fil_cos,opt)
      ndim = dimsizes(temp)
      varo = new((/12,ndim(0),ndim(1),ndim(2)/),"float")
      delete(ndim)
    end if
    varo(im,:,:,:) = lndistART(fil_cos,opt)
  end do

 ;Extract sums 
  namemods   = (/"nuc./acc. (pure + mixed)", "soot + coarse", \
                 "dust", "sea-salt"/)

 ;Extract sums for cropped region
  varo1 = dim_sum_n(varo(0:3,:,yloc-ext:yloc+ext,xloc-ext:xloc+ext),0)     ;nuc-acc mode ;small
  varo2 = dim_sum_n(varo(4:5,:,yloc-ext:yloc+ext,xloc-ext:xloc+ext),0)     ;soot+coarse
  varo3 = dim_sum_n(varo(6:9,:,yloc-ext:yloc+ext,xloc-ext:xloc+ext),0)     ;dust
  varo4 = dim_sum_n(varo(10:11,:,yloc-ext:yloc+ext,xloc-ext:xloc+ext),0)   ;sea

 ;PM2.5
  popt  = True
  popt@iy = yloc
  popt@ix = xloc
  popt@favg = False
  popt@fdry = True
  popt@nbld = nbdl
  popt@f3d  = True

  pm25 = pm25ART(fil_cos,popt)   ;Replace sea aerosol and concn very low

  ;Extract sums for cropped region
  pm25_4 = pm25(:,yloc-ext:yloc+ext,xloc-ext:xloc+ext)

;--------------------------------------------------------------------------------
  data = new((/4,nz/),"float")
  sdata = data

  do iz = 0, nz-1
    data(0,iz) = avg(varo1(iz,:,:))
    data(1,iz) = avg(varo2(iz,:,:))
    data(2,iz) = avg(varo3(iz,:,:))
    ;data(3,iz) = avg(varo4(iz,:,:))
    data(3,iz) = avg(pm25_4(iz,:,:))    

    sdata(0,iz) = stddev(varo1(iz,:,:))
    sdata(1,iz) = stddev(varo2(iz,:,:))
    sdata(2,iz) = stddev(varo3(iz,:,:))
    ;sdata(3,iz) = stddev(varo4(iz,:,:))
    sdata(3,iz) = stddev(pm25_4(iz,:,:))
  end do

  if (ens.eq.0) then
    edata = new((/4,dimsizes(einst),nz/),"float")
    esdata = edata
  end if

  edata(:,ens,:) = data
  esdata(:,ens,:) = sdata
  end do      ;ENSEMBLE LOOP

  foutname = "art-vertical-profile.nc"
  system("rm " + foutname)
  fout = addfile(foutname,"c")
  fout->Band1 = (/edata/)
  fout->Band2 = (/esdata/)

  end if ;fprcs = True
;-----------------------------------------------------------------------------------------
; Plotting Resource
;-----------------------------------------------------------------------------------------

  plot = new(4,"graphic")
  dum  = new(4,"graphic")


  wks_type = "x11"
  fname_out = "Figure_vp_" + sprinti("%0.2i",nf) 
  wks = gsn_open_wks(wks_type,fname_out)

  fncps_no = ":F21:"
  fncps_no = 21
  fncps_hgt = 0.02 *1.8

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
 res@tiYAxisString              = "Height [km]"
 res@tiXAxisString              = "[cm~S~-3~N~]"

 res@vpWidthF       = 0.5
 res@vpHeightF      = 0.7

 res@trXLog         = True
 res@trXMaxF        = 1.E4
 res@trXMinF        = 0.01
 res@trYMaxF        = 16.0 ;km

 gsres                   = True                        ; poly res
 gsres@tfPolyDrawOrder   = "Predraw"                   ; draw this first
 gsres@gsFillColor       = "SlateBlue"                 ; color chosen
 gsres@gsFillOpacityF    = 0.2

 do ip = 0, 3
 res@gsnLeftString = namemods(ip)
 if (ip.eq.3) then
   res@gsnLeftString = "  PM2.5"
   res@trXLog = False
   res@trXMinF = 0.
   res@trXMaxF = 20.0
   res@tiXAxisString              = "[~F8~m~F21~g~N~m~S~-3~N~]"
 end if
 if (ip.gt.0) then
   res@tiYAxisString = ""
   res@tmYLLabelsOn = False
 end if

 plot(ip) = gsn_csm_xy(wks,edata(ip,:,:),height,res)

 end do

 resP  = True
 resP@gsnPanelMainString = varo@time
 resP@gsnPanelDebug = True
 resP@gsnPanelXF = (/0.05,0.26,0.47,0.68/) + 0.02

 gsn_panel(wks,plot,(/1,4/), resP)


end
