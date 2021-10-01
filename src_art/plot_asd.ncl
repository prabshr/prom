;Plot aerosol size distribution for select modes

load "$prom_root/functions/lndistART.ncl"

begin
;User Settings
;--------------------------------------------------------------------------------
;  diri = "/p/largedata/hbn33/shrestha1/HET_ART_20150705/runART_6_HET/cosout/"
  diri = "/p/scratch/terrsysmp-art/runART_6_HET2483/cosmo_out/"

  fils_cos = systemfunc("ls " + diri + "*.nc")
  ;nf       = 0
  fil_cos  = fils_cos(nf)

;--------------------------------------------------------------------------------
  ;Geomeric standard deviations of lognormal modes
  sig = (/1.7, 2.0, 1.7, 2.0, 1.4, 2.5, 1.7, 1.6, 1.5, 1.9, 2.0, 1.7/)

  ;Obtain aggregate of 12 log-normal distribution modes
  opt    = True
  opt@iy = 90
  opt@ix = 150
  smodes = (/1,2,3,4,5,6,8,9,10,11,12/)
  do im = 0, dimsizes(smodes)-1
    opt@smode = smodes(im)
    if (im.eq.0) then
      nlDp  = lndistART(fil_cos,opt)
      Dp   = nlDp&Dp
      Z    = nlDp&Z
    else
      nlDp2 = lndistART(fil_cos,opt)
      nlDp  = nlDp + nlDp2
    end if
  end do

  opt@smode = 7
  nlDp_dust = lndistART(fil_cos,opt)

  do iz = 40,40 ;0, dimsizes(Z)-1 
;-----------------------------------------------------------------------------------------
; Plotting Resource
;-----------------------------------------------------------------------------------------
  wks_type = "x11"
  fname_out = "Figure_ASD_" + sprinti("%0.2i",iz)
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
 res@gsnLeftString              = Z(iz) + " m"

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

 res@xyLineThicknessF              = 3.
 res@xyDashPattern              = 0
 res@tiXAxisString              = "Dp [nm]"
 res@tiYAxisString             = "dN/dlnDp [cm-3]" 

 res@trYMaxF = max(nlDp_dust(:,iz))
 res@trYLog  = True

 res@trXMinF = 1.
 res@trXLog = True
 plot  =  gsn_csm_xy (wks,Dp,nlDp(:,iz),res) ; create plot

 res@xyLineColor = "red"
 dum   = gsn_csm_xy (wks,Dp,nlDp_dust(:,iz),res) ; create plot
 overlay(plot,dum)

 draw(plot)
 frame(wks)

 end do ;BIG LOOP     

end
