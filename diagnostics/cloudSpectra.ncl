load "$PROM/functions/dropSizeDistribution.ncl"
load "$PROM/functions/gfct.ncl"
;gfct and gamma function give same results
begin

  qx = 1.d0          ;g/m3 
  Nx = 300.d0        ;cm-3

  opt = True
  opt@ag = 0.124*1.d0
  opt@bg = 1./3.*1.d0
  opt@mu = 1./3.*1.d0
  opt@nu = 2./3.*1.d0

  ;
;  cnames = (/"cloudfastautoconv","cloudstandard","cloud_nue1mue1","cloudslowautoconv", \
;            "cloudheike","cloudsepa"/)
;  colors = (/"black", "red", "blue", "green","orange","cyan"/)
;  dpats  = (/0, 1, 2, 3, 4, 5 /)
;  cmu   = (/0.   , 1./3., 1., 1./3., 6., 6./)*1.d0
;  cnu   = (/1./3., 2./3., 1., 1.0,  1./3, 1./)*1.d0


  cnames   = (/"~F5~m~F21~=0,~F5~n~F21~=1/3","~F5~m~F21~=6,~F5~n~F21~=1"/)
  colors   = (/"black","red"/)
  dpats    = (/0, 1/)
  cmu      = (/0., 6./) *1.d0
  cnu      = (/1./3., 1./) * 1.d0
 
  ndim     = dimsizes(cnames)
  fdd      = dropSizeDistribution(qx,Nx,opt)
  D        = fdd&D
  fD       = new((/ndim,dimsizes(D)/),typeof(fdd))

  do ic    = 0,ndim-1
    opt@mu   = cmu(ic)
    opt@nu   = cnu(ic)
    fD(ic,:) = dropSizeDistribution(qx,Nx,opt)
  end do

  
  ;Plot
  fncps_no = 21
  fncps_hgt = 0.028
  wks   = gsn_open_wks ("png","Figure_hydrometeor_spectra")               ; send graphics to PNG file

  res                  = True 
  res@gsnDraw          = False
  res@gsnFrame         = False
  res@gsnLeftStringFontHeightF   = fncps_hgt
  res@gsnRightStringFontHeightF  = fncps_hgt

  res@gsnLeftString    = "QC = " + qx + " gm~S~-3~N~,  NC = " + Nx + " cm~S~-3~N~"
  res@trYLog           = True
  res@trYMinF          = 1.E-3
  ;res@trYMaxF          = 1.E2                   

  ;tm Resources
  res@tmXBLabelFont              = fncps_no
  res@tmYLLabelFont              = fncps_no
  res@tmYLLabelFontHeightF       = fncps_hgt
  res@tmXBLabelFontHeightF       = fncps_hgt
  res@tmXTOn           = False
  res@tmYROn           = False

  ;ti Resources
  res@tiXAxisFont                = fncps_no
  res@tiYAxisFont                = fncps_no
  res@tiXAxisFontHeightF         = fncps_hgt
  res@tiYAxisFontHeightF         = fncps_hgt
  res@tiYAxisString    = "f(D) " + fD@units
  res@tiXAxisString    = "D~B~p~N~ " + "[~F5~m~F21~m]"

  res@xyLineThicknessF = 3.0
  res@xyLineColors     = colors
  res@xyDashPatterns   = dpats

  res@vpXF  = 0.25 

  plot  = gsn_csm_xy (wks,D,fD,res) ; create plot

  lgres                    = True
  lgres@lgLineThicknessF   = 3.0
  lgres@lgLineColors       = colors    ; colors for legend lines
  lgres@vpWidthF           = 0.18                   ; width of legend (NDC)
  lgres@vpHeightF          = 0.14                  ; height of legend (NDC)
  lgres@lgPerimOn          = False                 ; turn off perimeter
  lgres@lgDashIndexes      = dpats              ; dash indexes
  lgres@lgLabelFontHeightF = 0.8

  lbid = gsn_create_legend(wks,ndim,cnames,lgres)
  amres = True
  amres@amParallelPosF   = 0.35		    ; move legend to the right
  amres@amOrthogonalPosF = -0.3             ; move the legend down
  annoid1 = gsn_add_annotation(plot,lbid,amres)   ; add legend to plot

  ; Add annotation
  txres               = True                      ; text mods desired
  txres@txFontHeightF = 0.032                    ; text font height
  txres@txJust        = "CenterLeft"              ; Default is "CenterCenter".

  gsn_text_ndc(wks,"a)",.005,.83,txres)

  draw(plot)
  frame(wks)

end
