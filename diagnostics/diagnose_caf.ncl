;Compares model ensemble members with BoxPol and JuXPol

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"

;Diagnostic plots of ensemble CAF

function data_extract(diri,opt)
  begin
  print("Extracting "  + diri)
  if (opt@ver.eq.1) then
    fname  = diri + "/CAF_" + opt@lnd + opt@itype + "_" + opt@inst + ".nc"
  else
    fname = diri + "/CAF_" + opt@pname + "_runE_" + opt@inst + "_" + opt@lnd + opt@itype + ".nc"
  end if
  nc_fil = addfile(fname,"r")
    ;print(nf + " " + fname)
  data = nc_fil->Band1    ;CAF [-]
  ndim = dimsizes(data)
  time = opt@sthr + ispan(0,ndim-1,1)*opt@nt/60.    ;5-15 minute interval data
  time@units = "hours since 0000-01-01 00:00:00"
  data!0  = "time"
  data&time = time
  return(data)
end

;Main Code

begin
  diri1 = "$PROM/scratch/data_prom1/"
  diri2 = "$PROM/scratch/CAF/"
 
  ; Read processed synthetic data PROM1,2,3, st= start times, 14, 9.59, 11 UTC, nt = 5 min interval
  opt       = True
  opt@inst  = 6      ;ens no
  opt@nt    = 5    ;data interval , mins
  opt@sthr  = 14   ;start hour
  opt@lnd   = "HET"
  opt@itype = "ART_SP4"
  opt@pname = "EMPA2_PROM1"

  opt@ver = 2 
  temp = data_extract(diri2 ,opt)
  ndim = dimsizes(temp)

  ;Read observations
  fobs1 = addfile(diri2 + "CAF_BoXPol_20150705.nc","r")
  data1 = fobs1->Band1
  tobs1  = fobs1->time

  fobs2 = addfile(diri2 + "CAF_JuXPol_20150705.nc","r")
  data2 = fobs2->Band1
  tobs2  = fobs2->time

;  itypes = (/2483,2463,2493,2583,2563,2593,2683,2663,2693/) ;, \
;            ;  2483,2463,2493,2583,2563,2593,2663,2663, \
;            ;  2483,2463,2493,2583,2563,2593,2663,2693/)
;  colors = (/"red","red","red","red","red","red","red","red","red"/) ;, \
;             ;    "blue","blue","blue","blue","blue","blue","blue","blue","blue", \
;             ;    "green","green","green","green","green","green","green","green","green"/)
;  colors = (/"orchid3","yellow3","grey20","orchid3","yellow3","grey20","orchid3","yellow3","grey20"/)
;
;  dpats  = (/0,0,0,0,0,0,0,0,0/) ;, \
;            ; 0,0,0,0,0,0,0,0,0, \
;            ; 0,0,0,0,0,0,0,0,0/)
;  dpats  = (/0,0,0,2,2,2,3,3,3/)
;  lnds   = (/"HET","HET","HET","HET","HET","HET","HET","HET","HET"/) ;, \
;            ;   "RHD","RHD","RHD","RHD","RHD","RHD","RHD","RHD","RHD", \
;             ;     "EHD","EHD","EHD","EHD","EHD","EHD","EHD","EHD","EHD"/)
;

  einst   = (/5,6,7,8,9/)
  colors  = (/"lightpink3","lightpink3","lightpink3","lightpink3","lightpink3"/)
  dpats   = (/1,0,1,1,1/)
  nprf    = dimsizes(einst) 
  data    = new((/nprf,ndim/),"float")

  do ip = 0, nprf-1
    ;opt@itype = itypes(ip)
    ;opt@lnd = lnds(ip)
    opt@inst = einst(ip)
    ;opt@pname = pnames(ip)
    data(ip,:) = data_extract(diri2 ,opt)
  end do

; Plot
  fncps_no = 21
  fncps_hgt = 0.028
  wks_type = "x11"
  fname_out = "Figure_caf_" + opt@lnd + "_" + opt@itype + "_" + opt@pname
  wks = gsn_open_wks(wks_type,fname_out)

  res = True

  ;gsn Resources
  res@gsnDraw = False
  res@gsnFrame = False

  ;tm Resources
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

  ;ti Resources
  res@tiXAxisFont                = fncps_no
  res@tiYAxisFont                = fncps_no
  res@tiXAxisFontHeightF         = fncps_hgt
  res@tiYAxisFontHeightF         = fncps_hgt
  res@tiYAxisString              = "Convective Area Fraction [-]"
  res@tiXAxisString              = "Time (UTC)"

  ;xy Resources
  res@xyDashPattern = 0
  res@xyLineThicknessF           = 3.
  res@xyLineColors  = colors
  res@xyDashPatterns = dpats

  ;vp Resources
  res@vpWidthF = 0.7
  res@vpHeightF = 0.5
  res@vpWidthF = 0.7

  ;tr Resources
  res@trYMaxF  = 0.6
  res@trXMaxF  = 16.5
  
  ;Time Labels 
  restick = True
  restick@ttmFormat = "%H%M"
  time_axis_labels(data&time,res,restick) ;

  plot                        = gsn_csm_xy (wks,data&time,data,res) ; create plot

; Observations
 delete(res@xyLineColors)
 res@xyLineColor                = "black"
 res@xyDashPatterns = dpats*0
 dum1      = gsn_csm_xy (wks,tobs1,data1,res) ; create plot
 res@xyLineColor                = "black"
 res@xyDashPatterns = dpats*0 + 1 
 dum2      = gsn_csm_xy (wks,tobs2,data2,res) ; create plot

 overlay(plot,dum1)
 overlay(plot,dum2)

; Bars for CFADS 
  xval = (/15.5,15.5/)
  yval = (/0, 1.0/)
  gres = True
  gres@gsLineColor = "grey20"
  duml1 = gsn_add_polyline(wks,plot,xval,yval,gres)
  xval = (/14.+45./60.,14. + 45./60./)
  yval = (/0, 1.0/)
  gres = True
  duml2 = gsn_add_polyline(wks,plot,xval,yval,gres)

;Add legend
 lgres                    = True
 lgres@lgLineColors       = (/"black","black","lightpink3"/) 
 lgres@lgLineThicknessF   = (/3.,3.,3./)
 lgres@lgDashIndexes      = (/0,1,0/)
 lgres@lgItemType         = "Lines"
 lgres@lgLabelFontHeightF = .3
 lgres@vpWidthF           = 0.13
 lgres@vpHeightF          = 0.13
 lgres@lgPerimOn          = False
 lbid = gsn_create_legend(wks,3,(/"BoXPol","JuXPol","Model"/),lgres)

 amres = True
 amres@amParallelPosF   = 0.25
 amres@amOrthogonalPosF = -0.35
 annoid1 = gsn_add_annotation(plot,lbid,amres)

 draw(plot)
 frame(wks)

end
