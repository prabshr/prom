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
  cnames = (/"cloudfastautoconv","cloudstandard","cloud_nue1mue1","cloudslowautoconv", \
            "cloudheike"/)
  colors = (/"black", "red", "blue", "green","orange"/)
  dpats  = (/0, 1, 2, 3, 4 /)
  cmu   = (/0.   , 1./3., 1., 1./3., 6./)*1.d0
  cnu   = (/1./3., 2./3., 1., 1.0, 1./3./)*1.d0
 
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
  wks   = gsn_open_wks ("x11","Figure_hydrometeor_spectra")               ; send graphics to PNG file

  res                  = True 
  res@gsnDraw          = False
  res@gsnFrame         = False
 
  res@gsnLeftString    = "qx = " + qx + " gm~S~-3~N~,  Nx = " + Nx + " cm~S~-3~N~"
  res@trYLog           = True
  res@trYMinF          = 1.E-3
  res@trYMaxF          = 1.E2                   

  res@tmXTOn           = False
  res@tmYROn           = False

  res@tiYAxisString    = "f(D) " + fD@units
  res@tiXAxisString    = "Dp " + "[~F5~m~F21~m]"

  res@xyLineThicknessF = 3.0
  res@xyLineColors     = colors
  res@xyDashPatterns   = dpats

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

  draw(plot)
  frame(wks)

end
