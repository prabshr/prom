begin
  filename = "discharge_2017.nc"
  fil_GRDC = "$PROJECT_chbn33/hbn331/database/bonnRadar/validation/NC_GRDC/GRDC_level_3.txt"

  ;;
  nc_fil   = addfile(filename,"r")
  data     = nc_fil->Band1
  obs      = data(0,0,:)        ;GRDC data
  model    = data(1,:,:)        ;TSMP model
  
  ;;
  stid     = asciiread(fil_GRDC,-1,"string")
  ndim     = dimsizes(stid)
  river    = new(ndim,"string")
  xloc     = ispan(0,ndim-1,1)
  
  ;;
 ;Open Plot resources
 fncps_no     = ":F21:"
 fncps_no     = 21
 fncps_hgt    = 0.02

  wks          = gsn_open_wks("x11","Figure_GRDC_validation")
  res          = True

 ;gsn Resources
 ;-------------------------------------------
  res@gsnDraw                    = False
  res@gsnFrame                   = False
 ;res@gsnMaximize                = False
 ;res@gsnLeftStringFontHeightF   = fncps_hgt
 ;res@gsnRightStringFontHeightF  = fncps_hgt

 ;tm Resources
 ;-------------------------------------------
  res@tmXBLabelFont              = fncps_no
  res@tmYLLabelFont              = fncps_no
  res@tmYLLabelFontHeightF       = fncps_hgt
  res@tmXBLabelFontHeightF       = fncps_hgt
  res@tmYROn                     = False
  res@tmXBOn                     = False
  res@tmXTOn                     = False
  res@tmXBLabelFontColor         = "black"
  res@tmYLLabelFontColor         = "black"
  res@tmXBLabelFontThicknessF    = 2.0
  res@tmYLLabelFontThicknessF    = 2.0
  res@tmXBMode                   = "Explicit"
  res@tmXBValues                 = xloc

  res@tiXAxisFontHeightF         = fncps_hgt
  res@tiYAxisFontHeightF         = fncps_hgt

  res@gsnXYBarChart       = True 

  ;res@gsnXYBarChartColors = colors_rgb

  res@trXMinF = -2 
  res@trXMaxF = ndim+ 1
  res@trYMaxF = 40.

  res@tiYAxisString  = "Mean daily discharge (m/s)"
  res@tiXAxisString  = "Rivers"

  plot    = gsn_csm_xy(wks,xloc, obs, res)

  res@gsnXYBarChart       = False

  res@xyMarkLineMode      = "Markers"
  res@xyMarkerColor       = "red"
  res@xyMarker            = 3

  dum     = gsn_csm_xy(wks,xloc, model, res)
  overlay(plot,dum)

  ;add river name
  txres = True
  txres@txAngleF = 90.
  txres@txFontHeightF = 0.01
  do nf    = 0, ndim-1
    temp   = str_split(stid(nf),";")
    gridno = temp(0)
    river(nf)  = temp(1)
    jy     = tointeger(temp(6))
    ix     = tointeger(temp(7))
    yoff   = max((/obs(nf),max(model(:,nf))/)) + 5
    text   = gsn_add_text(wks,plot,river(nf),nf, yoff , txres)
  end do

  draw(plot)
  frame(wks)

end
