;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"

;Diagnostic plots of H and LE

function data_extract(diri,opt)
  begin
  print("Extracting "  + diri)
  fname  = diri + "/flx_" + opt@lnd + opt@itype + "_" + opt@nf + ".nc"
  nc_fil = addfile(fname,"r")

  time = nc_fil->time
  data = new((/2,dimsizes(time)/),"float")
  data(0,:) = (/nc_fil->H/)     ;W/m2 [-]
  data(1,:) = (/nc_fil->LE/)

  data!1    = "time"
  data&time = time
  return(data)
end

begin
  diri = "./data_0/"
 
  ; Read processed synthetic data PROM1,2,3, 
  opt = True
  opt@nf = 6      ;ens no
  opt@lnd  = "HET"
  opt@itype = "2483"

  temp = data_extract(diri ,opt)
  ndim = dimsizes(temp)
  itypes = (/2483,2463,2493,2583,2563,2593,2663,2693, \
              2483,2463,2493,2583,2563,2593,2663,2663, \
                2483,2463,2493,2583,2563,2593,2663,2693/)
  colors = (/"red","red","red","red","red","red","red","red","red", \
                 "blue","blue","blue","blue","blue","blue","blue","blue","blue", \
                  "forestgreen","forestgreen","forestgreen","forestgreen","forestgreen", \
                  "forestgreen","forestgreen","forestgreen","forestgreen"/)
  dpats  = (/0,0,0,0,0,0,0,0,0, \
             0,0,0,0,0,0,0,0,0, \
             0,0,0,0,0,0,0,0,0/)
  lnds   = (/"EHD","EHD","EHD","EHD","EHD","EHD","EHD","EHD","EHD", \
             "HET","HET","HET","HET","HET","HET","HET","HET","HET", \
               "RHD","RHD","RHD","RHD","RHD","RHD","RHD","RHD","RHD"/)

  nprf  = dimsizes(itypes)
  data = new((/nprf,ndim(0),ndim(1)/),"float")

  do ip = 0, nprf-1
    opt@itype = itypes(ip)
    opt@lnd = lnds(ip)
    data(ip,:,:) = data_extract(diri ,opt)
  end do
  time = data&time

  ;time = tofloat(3. + (data&time)*24.)
  ;time@units = "hours since 2015-07-04 00:00:00"
  ;printVarSummary(time)

  plot = new(2,"graphic")

  wks_type = "png"
  fname_out = "Figure_surfaceflux"
  wks = gsn_open_wks(wks_type,fname_out)

  res = True
  res@gsnDraw = False
  res@gsnFrame = False

  res@tmXTOn   = False
  res@tmXBOn   = False

  res@xyLineThicknessF           = 3.
  res@xyLineColors  = colors
  res@xyDashPatterns = dpats

  res@tiXAxisString              = "Time (UTC)"

  res@trYMaxF = 350.
  res@vpHeightF = 0.4
  res@vpWidthF  = 0.7
 ;
  restick = True
  restick@ttmFormat = "%H"
  time_axis_labels(time,res,restick) ;

 ;PROM1 ------------
  res@tiYAxisString = "H [Wm-2]"
  plot(0)                        = gsn_csm_xy (wks,time,data(:,0,:),res) ; create plot
  res@tiYAxisString = "LE [Wm-2]"
  plot(1)                        = gsn_csm_xy (wks,time,data(:,1,:),res) ;

  ;Add legend
 lgres                    = True
 lgres@lgLineColors       = (/"blue","red","forestgreen"/)
 lgres@lgLineThicknessF   = (/3.,3.,3./)
 lgres@lgDashIndexes      = (/0, 0, 0/)
 lgres@lgItemType         = "Lines"
 lgres@lgLabelFontHeightF = 0.10
 lgres@vpWidthF           = 0.15
 lgres@vpHeightF          = 0.1
 lgres@lgPerimOn          = False
 lbid = gsn_create_legend(wks,3,(/"HET","EHD","RHD"/),lgres)

 amres = True
 amres@amParallelPosF   = 0.25
 amres@amOrthogonalPosF = -0.35
 annoid1 = gsn_add_annotation(plot(0),lbid,amres)

  resP = True
  gsn_panel(wks,plot,(/2,1/),resP)

end
