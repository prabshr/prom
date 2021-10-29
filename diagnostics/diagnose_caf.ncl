;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;Diagnostic plots of ensemble CAF

function data_extract(diri,opt)
  begin
  print("Extracting "  + diri)
  if (opt@ver.eq.1) then
    fname  = diri + "/CAF_" + opt@lnd + opt@itype + "_" + opt@nf + ".nc"
  else
    fname = diri + "/CAF_PROM1_runE_" + opt@nf + "_" + opt@lnd + opt@itype + ".nc"
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
  diri2 = "$PROM/scratch/data_v2/"
 
  ; Read processed synthetic data PROM1,2,3, st= start times, 14, 9.59, 11 UTC, nt = 5 min interval
  opt = True
  opt@nf = 6      ;ens no
  opt@nt   = 5    ;data interval , mins
  opt@sthr = 14   ;start hour
  opt@lnd  = "HET"
  opt@itype = 2483

  opt@ver = 2 
  temp = data_extract(diri2 ,opt)
  ndim = dimsizes(temp)
  itypes = (/2483,2463,2493,2583,2563,2593,2663,2693/)
  colors = (/"red","blue","green", \
                 "red","blue","green", \
                  "blue","green"/)
  dpats  = (/0,0,0,1,1,1,2,2,2/)
  nprf  = dimsizes(itypes)
  data = new((/nprf,ndim/),"float")

  do ip = 0, nprf-1
    opt@itype = itypes(ip)
    data(ip,:) = data_extract(diri2 ,opt)
  end do
  ;opt@ver = 1
  ;data2 = data_extract(diri1 ,opt)


  wks_type = "x11"
  fname_out = "Figure_caf_revised"
  wks = gsn_open_wks(wks_type,fname_out)

  res = True
  res@gsnDraw = False
  res@gsnFrame = False

  res@xyLineThicknessF           = 3.
  res@xyLineColors  = colors
  res@xyDashPatterns = dpats

  res@tiYAxisString              = "Convective Area Fraction [-]"
  res@tiXAxisString              = "Time (UTC)"

;PROM1 ------------
  plot                        = gsn_csm_xy (wks,data&time,data,res) ; create plot

  ;res@xyLineColor                = "black"
  ;dum      = gsn_csm_xy (wks,data1&time,data1,res) ; create plot
  ;overlay(plot,dum)
  draw(plot)
  frame(wks)

end
