;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

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

begin
  diri1 = "$PROM/scratch/data_prom1/"
  diri2 = "$PROM/scratch/CAF/"
 
  ; Read processed synthetic data PROM1,2,3, st= start times, 14, 9.59, 11 UTC, nt = 5 min interval
  opt = True
  opt@inst = 6      ;ens no
  opt@nt   = 5    ;data interval , mins
  opt@sthr = 14   ;start hour
  opt@lnd  = "HET"
  opt@itype = "ART"
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
  colors = (/"orchid3","yellow3","grey20","orchid3","yellow3","grey20","orchid3","yellow3","grey20"/)
;
;  dpats  = (/0,0,0,0,0,0,0,0,0/) ;, \
;            ; 0,0,0,0,0,0,0,0,0, \
;            ; 0,0,0,0,0,0,0,0,0/)
  dpats  = (/0,0,0,2,2,2,3,3,3/)
;  lnds   = (/"HET","HET","HET","HET","HET","HET","HET","HET","HET"/) ;, \
;            ;   "RHD","RHD","RHD","RHD","RHD","RHD","RHD","RHD","RHD", \
;             ;     "EHD","EHD","EHD","EHD","EHD","EHD","EHD","EHD","EHD"/)
;

   einst = (/5,6,7,8,9/)

;  nprf  = dimsizes(itypes)
  nprf = dimsizes(einst) 
  data = new((/nprf,ndim/),"float")

  do ip = 0, nprf-1
   ; opt@itype = itypes(ip)
   ; opt@lnd = lnds(ip)
    opt@inst = einst(ip)
    data(ip,:) = data_extract(diri2 ,opt)
  end do
  ;opt@ver = 1
  ;data2 = data_extract(diri1 ,opt)


  wks_type = "png"
  fname_out = "Figure_caf_revised"
  wks = gsn_open_wks(wks_type,fname_out)

  res = True
  res@gsnDraw = False
  res@gsnFrame = False

  res@xyLineThicknessF           = 3.
  ;res@xyLineColors  = colors
  ;res@xyDashPatterns = dpats

  res@xyDashPattern = 0
  res@xyLineColor = "lightpink3"
  res@tiYAxisString              = "Convective Area Fraction [-]"
  res@tiXAxisString              = "Time (UTC)"

  res@vpWidthF = 0.7
  res@vpHeightF = 0.5

  res@tmXTOn  = False
  res@tmYROn  = False

  res@trYMaxF  = 0.6
  res@trXMaxF  = 16.5
  res@vpWidthF = 0.7
;PROM1 ------------
  plot                        = gsn_csm_xy (wks,data&time,data,res) ; create plot

; Obs
 delete(res@xyLineColors)
   ; Obs
  delete(res@xyLineColors)
  res@xyLineColor                = "black"
  dum1      = gsn_csm_xy (wks,tobs1,data1,res) ; create plot


  res@xyLineColor                = "black"
  res@xyDashPattern = 1
  dum2      = gsn_csm_xy (wks,tobs2,data2,res) ; create plot

  overlay(plot,dum1)
  overlay(plot,dum2)

  draw(plot)
  frame(wks)

end
