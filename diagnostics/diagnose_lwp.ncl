;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;Diagnostic liquid/ice water paths or cloud top heights
;based on 2D time-series of masked storms
 

function data_extract(diri,opt)
  begin
  fname  = diri + "/" + opt@fname + "/" + opt@fname + "_" + \
               opt@pname + "_runE_" + opt@inst + "_" + opt@lnd + opt@itype + ".nc"
  print("Extracting "  +fname)
  nc_fil = addfile(fname,"r")
  data3d = (/nc_fil->$opt@varname$/)  ;/100.0  ;scaled to 100 kg/m2
  data = dim_avg_n(dim_avg_n(data3d,2),1)
  ndim = dimsizes(data)
  time = opt@sthr + ispan(0,ndim-1,1)*opt@nt/60.    ;5-15 minute interval data
  time@units = "hours since 0000-01-01 00:00:00"
  data!0  = "time"
  data&time = time
  return(data)
end

; Main Code

begin
  diri = "$PROM/scratch/"

  ; Read processed synthetic data PROM1,2,3, st= start times, 14, 9.59, 11 UTC, nt = 5 min interval
  opt         = True
  opt@inst    = 6
  opt@nt      = 5    ;data interval , mins
  opt@sthr    = 14   ;start hour
  opt@lnd     = "HET"
  opt@itype   = "ART"
  opt@pname   = "EMPA2_PROM1"
  opt@varname = "DX"
  opt@fname   = "QX"

  temp = data_extract(diri ,opt)
  ndim = dimsizes(temp)

;Compare different experiments
;  itypes = (/2483,2463,2493,2583,2563,2593,2683,2663,2693/) ;, \
;            ;  2483,2463,2493,2583,2563,2593,2663,2663, \
 ;;           ;  2483,2463,2493,2583,2563,2593,2663,2693/)
 ; colors = (/"red","red","red","red","red","red","red","red","red"/) ;, \
 ;            ;    "blue","blue","blue","blue","blue","blue","blue","blue","blue", \
 ;            ;    "green","green","green","green","green","green","green","green","green"/)
;  colors = (/"orchid3","yellow3","grey20","orchid3","yellow3","grey20","orchid3","yellow3","grey20"/)
;
;  dpats  = (/0,0,0,0,0,0,0,0,0/) ;, \
 ;           ; 0,0,0,0,0,0,0,0,0, \
  ;          ; 0,0,0,0,0,0,0,0,0/)
  dpats  = (/0,0,0,2,2,2,3,3,3/) 
;  lnds   = (/"HET","HET","HET","HET","HET","HET","HET","HET","HET"/) ;, \
;            ;   "RHD","RHD","RHD","RHD","RHD","RHD","RHD","RHD","RHD", \
 ;            ;     "EHD","EHD","EHD","EHD","EHD","EHD","EHD","EHD","EHD"/)
;
  ;einst = (/5,6,7,8,9/)
  pnames = (/"EMPA2_PROM1","EMPA2HC_PROM1"/)
  colors = (/"blue","red"/)
;;;;

  nprf  = dimsizes(pnames)
  data = new((/nprf,ndim/),"float")

  do ip = 0, nprf-1
    ;opt@itype = itypes(ip)
    ;opt@lnd = lnds(ip)
    ;opt@inst = einst(ip)
    opt@pname = pnames(ip)
    data(ip,:) = data_extract(diri ,opt)
  end do

 ;Plot

  wks_type = "x11"
  fname_out = "Figure_caf_revised"
  wks = gsn_open_wks(wks_type,fname_out)

  res = True
  res@gsnDraw = False
  res@gsnFrame = False

  res@xyLineThicknessF           = 3.
  res@xyLineColors  = colors
  res@xyDashPatterns = dpats

  res@tiYAxisString              = opt@varname ;+ " [100 * kgm~S~-2~N~]"
  res@tiXAxisString              = "Time (UTC)"

  res@vpWidthF  = 0.7
  res@vpHeightF = 0.5

  res@tmXTOn    = False
  res@tmYROn    = False

  res@trXMaxF = 16.5

  plot                        = gsn_csm_xy (wks,data&time,data,res) ; create plot

 ; Add legend
 lgres                    = True
 ;lgres@lgLineColors       = colors(0:ndim-1)
 lgres@lgLineThicknessF   = 3.
 lgres@lgDashIndexes      = (/0,0,0/)
 lgres@lgItemType         = "Lines"
 lgres@lgLabelFontHeightF = 0.10
 lgres@vpWidthF           = 0.15
 lgres@vpHeightF          = 0.1
 lgres@lgPerimOn          = False

 lbid1 = gsn_create_legend(wks,3,(/"con. CN","mar. CN","poll. CN"/), lgres)
 lgres@lgLineColors       = (/"black","black","black"/)
 lgres@lgDashIndexes      = (/0,2,3/)
 lbid2 = gsn_create_legend(wks,3,(/"def. IN","low IN","high IN"/), lgres)

 amres = True
 amres@amParallelPosF   = 0.2
 amres@amOrthogonalPosF = -0.35
; annoid1 = gsn_add_annotation(plot,lbid1,amres)
 amres@amParallelPosF   = -0.2
; annoid3 = gsn_add_annotation(plot,lbid2,amres)


 draw(plot)
 frame(wks)

end
