load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "./functions/cosmo_time_nclf.ncl"

;Explore ensemble statistics of melting layer
begin

 diri="/p/scratch/chbn33/hbn331/runHET2473/"
 ensemble_size = 40
 ntime  = 43 

 ;Read instances
 fil   = new(ensemble_size,"string")
 do nf = 0, ensemble_size-1
   if (nf.eq.0) then
     fil_ivr = addfile(diri + "tsmp_instance_" + nf + "/cosmo_out/" + "lfff00000000c.nc","r")
     hhl     = fil_ivr->HHL
   end if
   fils = systemfunc("ls " + diri + "tsmp_instance_" + nf + "/cosmo_out/" + "lfff*.nc")
   fil(nf)  = fils(ntime)
 end do
 str_time = cosmo_time_nclf(fil(0))

 nc_data = addfiles(fil,"r")
 T         = nc_data[:]->T
 QI        = nc_data[:]->QI
 QS        = nc_data[:]->QS
 QR        = nc_data[:]->QR
 QG        = nc_data[:]->QG
 QC        = nc_data[:]->QC
 QH        = nc_data[:]->QH

 ndim      = dimsizes(hhl)
 nz        = ndim(1)-1
 iy        = tointeger(ndim(2)/2)
 ix        = tointeger(ndim(3)/2)
 zmodel    = 0.5*(hhl(0,0:nz-1,:,:) + hhl(0,1:nz,:,:))

; 
 var1d_T     = T(:,nz-55:nz-1,iy,ix)  - 273.14 ; degree C
 var1d_QS    = QS(:,nz-55:nz-1,iy,ix)*1000.    ;g/kg
 var1d_QR    = QR(:,nz-55:nz-1,iy,ix)*1000.    ;g/kg
 var1d_QG    = QG(:,nz-55:nz-1,iy,ix)*1000.    ;g/kg
 var1d_QI    = QI(:,nz-55:nz-1,iy,ix)*1000.    ;g/kg
 var1d_QC    = QC(:,nz-55:nz-1,iy,ix)*1000.    ;g/kg
 var1d_QH    = QH(:,nz-55:nz-1,iy,ix)*1000.    ;g/kg

 ht_v        = zmodel(nz-55:nz-1,iy,ix)/1000.  ; km
; Plotting
 wks_type = "png"
 fname_out = "Figure_ensemble"
 wks = gsn_open_wks(wks_type,fname_out)

 fncps_no = ":F21:"
 fncps_no = 21
 fncps_hgt = 0.02

 res                            = True                 ; plot mods desired
 res@gsnDraw                    = False
 res@gsnFrame                   = False
 res@gsnMaximize                = False
 res@gsnLeftStringFontHeightF   = fncps_hgt
 res@gsnRightStringFontHeightF  = fncps_hgt
 res@gsnLeftString              = str_time

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


 res@trYMaxF                    = 5.0
 res@trXMaxF                    = 0.50
; res@vpXF                          = 0.15 ;0.1
; res@vpYF                          = 0.45
; res@vpHeightF                     = 0.35
; res@vpWidthF                      = 0.22

 res@xyLineThicknessF              = 3.

 res@tiYAxisString                 = "Height (km)"
 res@tiXAxisString                 = "g/kg"

 res@xyLineColor                   = "orange"
 plot = gsn_csm_xy (wks,var1d_QS,ht_v,res) ; create plot

 res@xyLineOpacityF                = 0.2
 res@xyLineColor                   = "purple"
 dumd                               = gsn_csm_xy (wks,var1d_QI,ht_v,res)
 overlay(plot,dumd)

 res@xyLineOpacityF                = 0.2
 res@xyLineColor                   = "red"
 dumd                               = gsn_csm_xy (wks,var1d_QG,ht_v,res)
 overlay(plot,dumd)

 res@xyLineOpacityF                = 0.2
 res@xyLineColor                   = "black"
 dumd                               = gsn_csm_xy (wks,var1d_QC,ht_v,res)
 overlay(plot,dumd)

 res@xyLineOpacityF                = 0.2
 res@xyLineColor                   = "magenta"
 dumd                               = gsn_csm_xy (wks,var1d_QH,ht_v,res)

 res@xyLineOpacityF                = 0.7
 res@xyLineColor                   = "blue"
 dumd                               = gsn_csm_xy (wks,var1d_QR,ht_v,res)
 overlay(plot,dumd)

 ;Add legend
 lgres                    = True
 lgres@lgLineColors       = (/"black","purple","red","orange","magenta","blue"/)
 lgres@lgLineThicknessF   = res@xyLineThicknessF
 lgres@lgDashIndexes      = (/0,0,0,0,0,0/)
 lgres@lgItemType         = "Lines"
 lgres@lgLabelFontHeightF = .1
 lgres@vpWidthF           = 0.1
 lgres@vpHeightF          = 0.18
 lgres@lgPerimOn          = False
 lbid = gsn_create_legend(wks,6,(/" qc"," qi"," qg"," qs"," qh"," qr"/),lgres)

 amres = True
 amres@amParallelPosF   = 0.4
 amres@amOrthogonalPosF = -0.35
 annoid1 = gsn_add_annotation(plot,lbid,amres)

 draw(plot)
 frame(wks)
end

