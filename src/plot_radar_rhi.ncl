;Plot the plan veiw of polarimetric variables at low and mid levels
; and cross-section (frhi = True) replaces ppi at high levels 

;P. Shrestha
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$PROM/functions/shade_terrain.ncl"
load "$PROM/functions/read_radarnc.ncl"
load "$PROM/functions/load_radarColors.ncl"
begin

;User Settings
;---------------------------------------------------------------------------
 ;pid       = 0 
 ydates    = (/"20150705","20160513","20170706"/)
 fnss      = (/18, 18, 12/)    ;(/15:30 UTC ,10:30 UTC, UTC/)  
 trXMin    = (/-25,-30,-40/) 
 trXMax    = (/60, 10, 50/)
 xdists    = (/-28, 40, -28/)  ;km; ix           = 150   ;132
 updates   = (/1, 0, 1/)  ;

 frhi        = True

 glon_idB    = 7.071663 ; Bonn Radar
 glat_idB    = 50.73052 ;
 radar     = "BoXPol" 

 fn        = fnss(pid)  
 ydate     = ydates(pid)
 xdist     = xdists(pid) 
 diri      = "/daten01/z4/database/PROM/" + ydate + "/" + radar + "/GRIDS/"
;Read netcdf files 
;---------------------------------------------------------------------------

 fils_var     = systemfunc("ls " + diri + "/*.nc")

 ndimf        = dimsizes(fils_var)

 fns = fn
 fne = fn

 do nf= fns,fne ;0, ndimf-1
   ;Read data
   fname = fils_var(nf)
   if (updates(pid).eq.1 .and. pid.eq.2) then
     fname = "/daten01/z4/database/PROM/updated/GRIDS//grid_volume_bon_eac_201707061400.nc"
   end if
   if (updates(pid).eq.1 .and. pid.eq.0) then
     fname = "/daten01/z4/database/PROM/updated/GRIDS//grid_volume_bon_eac_201507051530.nc"
   end if
   fname = "./grid_volume_bon_new_201605131030.nc" 
   print("Extracting  " + fname)
   Zh  = read_radarnc(fname,0, glon_idB,glat_idB) 
   Zdr = read_radarnc(fname,1, glon_idB,glat_idB)
   Kdp = read_radarnc(fname,2, glon_idB,glat_idB)
   Rho = read_radarnc(fname,3, glon_idB,glat_idB)
 end do

 x = Zh&X
 ix = closest_val(xdist,x)
 print("c-s at ix = " + ix + " X:" + x(ix) + " km" )

 Zh_cs  = Zh(:,:,ix)
 Zdr_cs = Zdr(:,:,ix)
 Kdp_cs = Kdp(:,:,ix)
 Rho_cs = Rho(:,:,ix)

 sdims  = dimsizes(Zh_cs)
 xval2x = conform_dims (sdims,Zh_cs&Y,1)
 zval2x = conform_dims (sdims,Zh_cs&Z,0) 
 
 print("Plotting ...")
 ;------------------------------------------------------------------------
 fncps_no = ":F21:"
 fncps_no = 21
 fncps_hgt = 0.035

 npol = 4 
 plot     = new(npol,"graphic")
 if (frhi) then
   dum0   = new(npol,"graphic")
   dum1   = new(npol,"graphic")
   dum2   = new(npol,"graphic")
   dum3   = new(npol,"graphic")
   dum4   = new(npol,"graphic")
   dum_cs = new(npol,"graphic")
 end if


 res      = True
 wks      = gsn_open_wks("x11", "Figure_" + radar +  "_t" + fn + "_ix" + ix)

 colors = load_radarColors(0)
 colors(2) = "white"
 gsn_define_colormap(wks,colors)

 res      = True
 ;gsn Resources
 ;---------------
 res@gsnDraw                    = False
 res@gsnFrame                   = False
 res@gsnLeftStringFontHeightF   = fncps_hgt
 res@gsnRightStringFontHeightF  = fncps_hgt
 res@gsnMaximize                = False
 ; makes no diff res@gsnSpreadColors            = True
 res2 = res   ;for Zh contour

 ;tm Resources
 ;--------------- 
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
 ;-----------------
 res@tiXAxisFont                = fncps_no
 res@tiYAxisFont                = fncps_no
 res@tiXAxisFontHeightF         = fncps_hgt
 res@tiYAxisFontHeightF         = fncps_hgt

 res3 = res   ;for frhi = True

 ; for terrain shading
 xypoly  = shade_terrain(xval2x,zval2x)
 xpoly   = xypoly(0,:)
 ypoly   = xypoly(1,:)

 ;Trick to linearize the plot
 res3@cnFillOn   = True
 res3@cnFillMode = "CellFill"
 res3@cnLinesOn  = False
 res3@cnLevelSelectionMode       = "ExplicitLevels"
 res3@cnInfoLabelOn    = False

 res3@trGridType = "TriangularMesh"
 res3@sfXArray   = xval2x
 res3@sfYArray   = zval2x
 res3@gsnLeftString = ""
 res3@vpWidthF = 0.52

 res3@trXMinF    = trXMin(pid) ;-25 ;min(xval2x)
 res3@trXMaxF    = trXMax(pid) ; 60 ;max(xval2x) 
 res3@trYMaxF    =  15

 res3@tmXBMode   = "Explicit"
 res3@tmXBValues = xval2x(0,ispan(0,sdims(1)-1,50)) 
 res3@tmXBLabels = tointeger(res3@tmXBValues*100.)/100.

 res3@lbOrientation = "Vertical"
 res3@lbLabelFontHeightF = 0.025

 res3@tmYLLabelsOn = True
 res3@tiYAxisString = "Height [km]"
 res3@tiXAxisString = "Y Range [km]"
 res3@gsnLeftString = radar 
 res3@gsnRightString = Zh_cs@long_name + " " + Zh_cs@units
 res3@cnLevels                   =  colors@zhh_levs
 plot(0)  = gsn_csm_contour(wks,Zh_cs,res3)

 res3@tiYAxisString = ""
 res3@gsnLeftString = ""
 res3@tmYLLabelsOn = False
 res3@gsnRightString = Zdr_cs@long_name + " " + Zdr_cs@units
 res3@cnLevels                   =  colors@zdr_levs
 plot(1)  = gsn_csm_contour(wks,Zdr_cs,res3)

 ;res3@tmYLLabelsOn = True
 ;res3@tiYAxisString = "Height [km]" 
 res3@gsnRightString = Kdp_cs@long_name + " " + Kdp_cs@units
 res3@cnLevels                   =  colors@kdp_levs
 plot(2)  = gsn_csm_contour(wks,Kdp_cs,res3)

 res3@tiYAxisString = ""
 res3@tmYLLabelsOn = False
 res3@gsnRightString = Rho_cs@long_name + " " + Rho_cs@units
 res3@cnLevels                   =  colors@rho_levs
 plot(3)  = gsn_csm_contour(wks,Rho_cs,res3)

 if (npol.gt.10) then
   delete(res3@cnLevels)
   res3@cnLineThicknessF = 2.0
   res3@cnFillOn = False
   res3@cnLinesOn = True
   res3@cnLineLabelsOn = False
   res3@cnInfoLabelOn = False
   res3@cnLevels                   = (/20.,30.,40.,50./) ; ipol_1levs(1:)
   dum = gsn_csm_contour(wks,Zh_cs,res3)
   overlay(plot(1),dum)
   delete(res3@cnLevels)
 end if

 pres                  = True
 pres@gsFillColor      = "gray"
 ;pres@gsLineThicknessF = 2.0          ; twice as thick
 do ipol = 0, npol-1
   dum0(ipol) = gsn_add_polygon(wks,plot(ipol),xpoly,ypoly,pres)
   dum1(ipol) = gsn_add_polyline(wks,plot(ipol),xpoly(0:sdims(1)-1),ypoly(0:sdims(1)-1),pres)
 end do

 txres = True
 txres@txFontHeightF = 0.012
; gsn_text_ndc(wks,"b)",.012,.62,txres)

 resP = True
 resP@gsnMaximize    = True
 ;resP@gsnPanelDebug = True
 resP@gsnPanelXF = (/-1,0.29,0.53,0.77/)
 gsn_panel(wks,plot,(/1,4/), resP)

 fout = addfile("temp.nc","c")
 fout->ZDR = Zdr_cs
end
