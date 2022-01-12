;Plot Nx or Mx of different aerosol modes at specified height
;to explore patterns of different air mass
; Data are ouput at specified heights in 4 different subplots
; for nuc+acc, dust, sea, soot+coarse
; Data are interpolated to specified height

load "$prom_root/functions/lndistART.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

begin

;-----------------------------------------------------------------------------------------
;#1 USER SETTINGS
;-----------------------------------------------------------------------------------------
 pname = "PROM1";
 inst  = 6
 lnd   = "HET"
 itype = "ART"
 diri  = "/p/largedata/hbn33/shrestha1/"

 height  = 7000. ;4000. ; 2000. ; m
 colors  = (/"lightpink3","forestgreen","steelblue","gray20","gray"/)
 dashs   = (/0,0,0,0,0/)
; nf      = 0 ;18+15    ;midnight+15 hrs
 nbdl    = 20
;--------------------------------------------------------------------------------

  itype_gscp  = lnd + itype
  sfold       = "runART_" + inst + "_" + itype_gscp
  expname     = itype_gscp + "_" + pname
  subevent    = expname + "/" + sfold
  diri_cos    = diri + "/" + subevent + "/cosout/"
  diri_cosI   = diri + "/" + subevent + "/cosout/ivr/"

  fils_cos = systemfunc("ls " + diri_cos + "*.nc")

;-----------------------------------------------------------------------------------------
  ;Geomeric standard deviations of lognormal modes
  sig = (/1.7, 2.0, 1.7, 2.0, 1.4, 2.5, 1.7, 1.6, 1.5, 1.9, 2.0, 1.7/)

  ;Obtain aggregate of 12 log-normal distribution modes
  opt    = True
  opt@iy = 90
  opt@ix = 150
  opt@favg = True
  opt@fout = 1  ;1-Ni(cm-3), 2-Mi(m3/m3), 3-Dp (nm)2  
  varname  = (/"Nx (cm-3)", "Mx (m3/m3)", "Dp (nm)"/)
  opt@nbdl = nbdl

  smodes = ispan(1,12,1) 

;-----------------------------------------------------------------------------------------
;#2a Extract GEO 
;-----------------------------------------------------------------------------------------
  nc_cosI = addfile(diri_cosI + "lfff00000000c.nc","r")
  tmp     = nc_cosI->lat
  ndim    = dimsizes(tmp)

  hhl     = nc_cosI->HHL(0,:,nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1)
  lat     = nc_cosI->lat(nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1)
  lon     = nc_cosI->lon(nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1)

  delete(ndim)
  ndim    = dimsizes(hhl)
  Z       = 0.5 * (hhl(0:ndim(0)-2,:,:) + hhl(1:ndim(0)-1,:,:))
  Z       = Z(::-1,:,:)

;-----------------------------------------------------------------------------------------
;#2b Extract lognormal ASD for each mode
;-----------------------------------------------------------------------------------------
  do nf = 0, dimsizes(fils_cos)-1

  fil_cos  = fils_cos(nf)

  do im = 0, dimsizes(smodes)-1
    opt@smode = smodes(im)
    if (im.eq.0) then
      temp = lndistART(fil_cos,opt)
      ndim = dimsizes(temp)
      varo = new((/12,ndim(0),ndim(1),ndim(2)/),"float")
    end if
    varo(im,:,:,:) = lndistART(fil_cos,opt)
  end do

  printVarSummary(varo)
  printVarSummary(Z)

 ;Extract sums at specific height
  varo1 = dim_sum_n(varo(0:3,:,:,:),0)     ;nuc-acc mode ;small
  varo2 = dim_sum_n(varo(4:5,:,:,:),0)     ;soot+coarse
  varo3 = dim_sum_n(varo(6:9,:,:,:),0)     ;dust
  varo4 = dim_sum_n(varo(10:11,:,:,:),0)   ;sea

  varoz1 = wrf_user_intrp3d(varo1,Z,"h",height,0., False)
  varoz2 = wrf_user_intrp3d(varo2,Z,"h",height,0., False)
  varoz3 = wrf_user_intrp3d(varo3,Z,"h",height,0., False)
  varoz4 = wrf_user_intrp3d(varo4,Z,"h",height,0., False) 

  varoz1@long_name = "Nuc-Acc"
  varoz2@long_name = "Soot-Coars"
  varoz3@long_name = "Dust"
  varoz4@long_name = "Sea"

;-----------------------------------------------------------------------------------------
; Plotting Resource
;-----------------------------------------------------------------------------------------

  plot = new(4,"graphic")

  wks_type = "png"
  fname_out = "Figure_airmass" + sprinti("%0.2i",nf)
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

 res@cnFillOn = True
 res@cnLevelSelectionMode    = "ExplicitLevels"
 res@cnFillOpacityF          = 0.5
 res@cnInfoLabelFontColor    = "grey"
 res@cnLevels                = (/0.05,0.1,0.5,1.0,5.0,10.,50.,100.,500.,1000.,1500.,2000./)
 res@cnLineLabelsOn           = False
 res@cnInfoLabelOn           = False
 res@lbLabelBarOn            = False


 plot(0) = gsn_csm_contour(wks,varoz1,res)
 plot(1) = gsn_csm_contour(wks,varoz2,res)
 plot(2) = gsn_csm_contour(wks,varoz3,res)
 plot(3) = gsn_csm_contour(wks,varoz4,res)

 resP  = True
 resP@gsnPanelLabelBar = True
 resP@lbOrientation    = "Vertical"
 resP@gsnPanelFigureStrings = (/varoz1@long_name,varoz2@long_name,varoz3@long_name,varoz4@long_name/)
 resP@gsnPanelMainString = varname(opt@fout-1) + " Height: " + height + " m a.m.s.l." + "(" + varo@time + ")"
 gsn_panel(wks,plot,(/2,2/), resP)

 end do
end
