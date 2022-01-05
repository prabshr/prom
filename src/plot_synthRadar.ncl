;Plot the plan veiw of polarimetric variables at low and mid levels
; and cross-section (frhi = True) replaces ppi at high levels 

;P. Shrestha
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$PROM/functions/shade_terrain.ncl"
load "$PROM/functions/load_radarColors.ncl"
external EXT1 "$prom_root/functions/readdim.so"
external EXT2 "$prom_root/functions/readascii.so"

begin
 ;pid       = 1 
 ydates    = (/"20150705","20160513","20170706"/)
 insts     = (/6, 8, 11/)
 fns       = (/11,15,18/)
 ixs       = (/90,124,95/)
 lat_exts  = (/51.2, 50.8, 51./)
 lat_exts0 = (/40.0, 40.0, 50.42/)

;SPP proposal
; insts     = (/5, -1, -1/)
; fns       = (/21,-1,-1/)
; ixs       = (/184,-1,-1/) - 65
; lat_exts  = (/51.6, -1, -1/)
 ;----

 lat_ext0  = lat_exts0(pid)
 lat_ext   = lat_exts(pid)
 fn        = fns(pid)
 ix        = ixs(pid)
 inst      = insts(pid)
 ydate     = ydates(pid)

 diri_dBz    = "/daten04/pshrestha/shared-data/bonnRadar_" + ydate + "/Synth/HET2483_" + inst +"/"
 pvars       = (/"Z~B~H~N~","Z~B~DR~N~","K~B~DP~N~","~F8~r~B~~F10~hv~E~"/)
 punits      = (/"[dBZ]","[dB]","[deg/km]", "[-]"/)

 frhi        = True ;replace 3rd level ppi with rhi
 glon_idB    = 7.071663 ; Bonn Radar
 glat_idB    = 50.73052 ;
 dlat        = 1.0 ;0.5
 dlon        = 1.0 ; 0.5

 if (pid.eq.0) then
   mpc  = (/50. + 40./60., 6. + 50./60., 51. + 15./60., 7. + 30./60./)
   izlevs      = (/25, 44, 54/)
 else if (pid.eq.1) then
   mpc  = (/50. + 10./60., 7. + 0./60., 51. + 0./60., 8. + 0./60./) 
   izlevs      = (/25, 41, 54/)
 else if (pid.eq.2) then
   mpc  = (/50. + 5./60., 6. + 40./60., 51. + 0./60., 7. + 50./60./)
   izlevs      = (/25, 44, 54/)
 else
   mpc  = (/glat_idB - dlat, glon_idB - dlon, glat_idB + dlat, glon_idB + dlon/)
 end if
 end if
 end if

;;;;;;;;;;

 diri_input  = diri_dBz + "/input"
 diri_output = diri_dBz + "/output"

 fils_var    = systemfunc("ls " + diri_output + "/dbz-tm-1m_*")
 fils_elv    = systemfunc("ls " + diri_input + "/elevation_*")
 fil_lat     = systemfunc("ls " + diri_input + "/lat_*")
 fil_lon     = systemfunc("ls " + diri_input + "/lon_*")

 ;Read header information to extract nx,ny,nz
 array = new(3,integer)
 EXT1::readdim(array(0),array(1),array(2),fils_elv(fn))

 latPV = new((/1,1,array(1),array(0)/),"float")
 lonPV = new((/1,1,array(1),array(0)/),"float")
 hgtPV = new((/1,array(2),array(1),array(0)/),"float")
 EXT2::readascii(array(0),array(1),1,latPV,1,1,fil_lat)
 EXT2::readascii(array(0),array(1),1,lonPV,1,1,fil_lon)
 EXT2::readascii(array(0),array(1),array(2),hgtPV,1,1,fils_elv(fn))   ;Model mass point heights

 ;Extract lat-lon
 lat_c = latPV(0,0,:,:)
 lon_c = lonPV(0,0,:,:)
 
 ;Average height at iz
 ;Need to reverse the vertical levels
 mhgt  = hgtPV(0,::-1,:,:) 
 phgt  = mhgt(0,:,:)
 phgt@lat2d  = lat_c
 phgt@lon2d  = lon_c

 ;
 npol          = dimsizes(pvars) 
 nz            = dimsizes(izlevs)
 ;np            = (npol-1)*nz           ;skip Zh
 np            = npol*nz

 varh         = new((/npol,array(2),array(1),array(0)/),"float")    ;Initialize variable for plotting
 varh@_FillValue = -99
 varPV        = new((/8,array(2),array(1),array(0)/),"float")
 varPV@_FillValue = -99
 fname        = fils_var(fn)
 print(" " + fname)

 EXT2::readascii(array(0),array(1),array(2),varPV,8,1,fname)
 ;Need to reverse the vertical levels
 varh(0,:,:,:) = varPV(1,::-1,:,:)
 varh(1,:,:,:) = varPV(2,::-1,:,:)
 varh(2,:,:,:) = varPV(3,::-1,:,:)
 varh(3,:,:,:) = varPV(4,::-1,:,:)

 Zh  = varh(0,:,:,:)
 Zv  = varh(1,:,:,:)
 KDP = varh(2,:,:,:)
 RHV = varh(3,:,:,:)
 varh(0,:,:,:) = where(Zh.lt.-5., varh@_FillValue, Zh)
 varh(1,:,:,:) = where(Zh.lt.-5., varh@_FillValue, Zh-Zv)
 varh(2,:,:,:) = where(Zh.lt.-5., varh@_FillValue, KDP)
 varh(3,:,:,:) = where(Zh.lt.-5., varh@_FillValue, RHV)
 delete([/Zh, Zv, KDP, RHV/])

 ;Variable to plot
 varh@lat2d    = lat_c
 varh@lon2d    = lon_c

 print("Plotting ...")

 fncps_no = ":F21:"
 fncps_no = 21
 fncps_hgt = 0.035

 plot     = new(np,"graphic")
 if (frhi) then
   dum0   = new(np,"graphic")
   dum1   = new(np,"graphic")
   dum2   = new(np,"graphic")
   dum3   = new(np,"graphic")
   dum4   = new(np,"graphic")
   dum_cs = new(np,"graphic")
 end if


 res      = True
 ppid     = pid + 1
 wks      = gsn_open_wks("png", "Figure_case" + ppid + "_i" + inst + "_t" + fn + "_ix" + ix)

 colors   = load_radarColors(0)

 gsn_define_colormap(wks,colors)

 ;CPS
 res@vpXF = 0.4
 res@gsnMaximize = False

 ;gsn Resources
 ;---------------
 res@gsnDraw                    = False
 res@gsnFrame                   = False
 res@gsnLeftStringFontHeightF   = fncps_hgt
 res@gsnRightStringFontHeightF  = fncps_hgt
 res@gsnMaximize                = False
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

 ;ti Resources
 ;-----------------
 res@tiXAxisFont                = fncps_no
 res@tiYAxisFont                = fncps_no
 res@tiXAxisFontHeightF         = fncps_hgt
 res@tiYAxisFontHeightF         = fncps_hgt

 res3 = res   ;for frhi = True

 ;pm resources
 ;--------------------
 res@pmTickMarkDisplayMode      = "Always"
 res@pmLabelBarOrthogonalPosF   = 0.3
 ;res@pmLabelBarParallelPosF     = 0.2

 ;mp Resources
 ;----------------------

 res@mpLimitMode                  = "Corners"           ; method to zoom
 res@mpLeftCornerLatF             = max((/mpc(0),min(lat_c)/)) 
 res@mpLeftCornerLonF             = max((/mpc(1),min(lon_c)/))
 res@mpRightCornerLatF            = min((/mpc(2),max(lat_c)/))
 res@mpRightCornerLonF            = min((/mpc(3),max(lon_c)/))
 res@mpOutlineOn                  = True
 res@mpOutlineBoundarySets        = "National"
 res@mpDataSetName                = "Earth..4"    ; database for non-USA divisions
 res@mpDataBaseVersion            = "MediumRes"   ; Medium resolution database
; res@mpOutlineSpecifiers          = "India:Uttaranchal"
 res@mpGridAndLimbOn              = True
 res@mpNationalLineThicknessF     = 3.0
 res@mpGeophysicalLineThicknessF  = 3.0
 res@mpProvincialLineThicknessF   = 2.0
 res@mpPerimOn                    = True
 res@mpFillOn                     = False


 ;cn Resources
 ;-------------------------
 res@cnFillOn   = True
 res@cnFillMode = "CellFill"
 res@cnLinesOn  = False
 res@cnLevelSelectionMode       = "ExplicitLevels"
 res@cnInfoLabelOn    = False

 ;tr Resources
 ;-----------------------------
 ;res@trXMinF = 60
 ;res@trYMinF = 60
 ;res@trXMaxF = 140
 ;res@trYMaxF = 140

 delete(res2@gsnLeftStringFontHeightF)
 delete(res2@gsnRightStringFontHeightF)
 res2@cnLineThicknessF = 2.0
 res2@cnFillOn = False
 res2@cnLinesOn = True
 res2@cnLineLabelsOn = False
 ;res2@cnLineColor    = "Gray"
 res2@cnInfoLabelOn = False
 res2@cnLevelSelectionMode       = "ExplicitLevels"

 ip = -1
 ; BIG LOOP

 do ik = 0, nz-1

 iz = izlevs(ik)
 temp = mhgt(iz,:,:) - phgt ;.ag.l

 print(iz + " " + avg(temp))
 avghgt = round((avg(temp)/100.),3)*100
 print("Average model mass point height at iz =" +iz + " is " + avghgt + " m a.m.s.l")

 if (frhi) then
   nzz = nz-2
 else
   nzz = nz-1
 end if

 if (ik.eq.nzz) then
   res@lbLabelBarOn = True
   res@tmXBLabelsOn = True
 else
   res@lbLabelBarOn = False
   res@tmXBLabelsOn = False
 end if

 ;Remove redundant label bar
 if (frhi) then
   res@lbLabelBarOn = False
 end if

 do ipol = 0, npol-1   

   if (ipol.eq.0) then
     ipol_levs = colors@zhh_levs
     res@tmYLLabelsOn = True
     res@gsnRightString = ""
   else
     res@tmYLLabelsOn = False
   end if

   if (ipol.eq.1) then
     ipol_levs = colors@zdr_levs
     res@gsnRightString = ""
   end if
   if (ipol.eq.2) then
     ipol_levs = colors@kdp_levs
     res@gsnRightString = ""
   end if
   if (ipol.eq.3) then
     ipol_levs = colors@rho_levs
     res@gsnRightString = avghgt + " m a.g.l."
   end if

   ip = ip + 1

   if (ik.eq.0) then
     res@gsnLeftString              = pvars(ipol) + " " + punits(ipol) 
   else
     res@gsnLeftString   = ""
   end if

   res@cnLevels                   =  ipol_levs
   plot(ip)  = gsn_csm_contour_map(wks,varh(ipol,iz,:,:) ,res)

   if (ipol.eq.1) then
     res2@cnLevels                   = (/20.,30.,40.,50./) ; ipol_1levs(1:)
     dum      = gsn_csm_contour(wks,varh(0,iz,:,:),res2)
     ;overlay(plot(ip),dum)
   end if

   if (frhi) then
   if (ik.lt.nz-1) then
     pres                  = True
     pres@gsLineColor      = "gray"
     ypoly_cs = (/mpc(0),mpc(2)/)
     xpoly_cs = (/avg(lon_c(:,ix)),avg(lon_c(:,ix))/)
     dum_cs(ip) = gsn_add_polyline(wks,plot(ip),xpoly_cs,ypoly_cs,pres)
     delete(pres)
   end if
   end if 

   txres               = True
   txres@txFont        = fncps_no
   txres@txBackgroundFillColor = "white"
   ;dumtxt = gsn_add_text(wks,plot(0),"A",lon_plane(3),lat_plane(3),txres) ;
   ;dumtxt = gsn_add_text(wks,plot(0),"A'",lon_plane(dimsizes(lon_plane)-1),lat_plane(dimsizes(lat_plane)-1),txres) ;
   txres@txBackgroundFillColor = -1
   txres@txFontThicknessF      = 5.
   txres@txFontHeightF           = 0.045
   txres@txFontColor           = "black"
   dumtxt = gsn_add_text(wks,plot(ip),"x",glon_idB,glat_idB,txres)

   ;dum      = gsn_csm_contour(wks,phgt,res2)
   ;overlay(plot(ip),dum)
   ;print(ip + " "  +ik + " " + ipol +" "+ res@tmXBLabelsOn)

   ;Here we replace 3rd level ppi with rhi
   ;----------------------------------------------------------
   if (frhi) then
   if (ik.eq.nz-1 ) then
      spol          = (/varh(ipol,:,:,ix)/)         ;
      sdims         = dimsizes(spol)
      spol!0        = "Z"
      spol!1        = "Y"
      spol&Y        = lat_c(:,ix) ; degrees 
      spol0          = (/varh(0,:,:,ix)/)    
      spol0!0        = "Z"
      spol0!1        = "Y"
      spol0&Y        = lat_c(:,ix) ; degrees

      xval2x        = conform_dims(sdims,spol&Y,1)*1.
      zval2x        = mhgt(:,:,ix)/1000.  ;km

      ; for terrain shading
      xypoly  = shade_terrain(xval2x,zval2x)
      xpoly   = xypoly(0,:)
      ypoly   = xypoly(1,:)
      xpoly1   = (/mpc(0),mpc(2)/)
      ypoly1  = (/avg(mhgt(izlevs(0),:,:)),avg(mhgt(izlevs(0),:,:))/)
      ypoly2  = (/avg(mhgt(izlevs(1),:,:)),avg(mhgt(izlevs(1),:,:))/) 

      ;Trick to linearize the plot
      res3@cnFillOn   = True
      res3@cnFillMode = "CellFill"
      res3@cnLinesOn  = False
      res3@cnLevelSelectionMode       = "ExplicitLevels"
      res3@cnInfoLabelOn    = False
      res3@cnLevels                   =  ipol_levs

      ;
      res3@trGridType = "TriangularMesh"
      res3@sfXArray   = xval2x
      res3@sfYArray   = zval2x
      res3@gsnLeftString = ""
      res3@gsnRightString = ""

      res3@vpWidthF = 0.4

      res3@trXMinF    = max((/lat_ext0,mpc(0)/))
      res3@trXMaxF    = lat_ext ;min((/lat_ext,mpc(2)/))
      res3@trYMaxF    = 15.

      res3@tmXBMode   = "Explicit"
      res3@tmXBValues = xval2x(0,ispan(0,sdims(1)-1,20)) 
      res3@tmXBLabels = tointeger((res3@tmXBValues)*100)/100.
 ;CPS     res3@tmXBLabelsOn = False

      res3@lbOrientation = "Vertical"

      res3@tiXAxisString = " " ;"Latitude [~F34~0~F~N]"
      if (ipol.eq.0) then
        res3@tmYLLabelsOn = True
        res3@tiYAxisString = "Height [km]"
      else
;CPS        res3@tiYAxisString = ""
;CPS        res3@tmYLLabelsOn = False
      end if

      plot(ip)  = gsn_csm_contour(wks,spol,res3)
      if (ipol.eq.1) then
        delete(res3@cnLevels)
        res3@cnLineThicknessF = 2.0
        res3@cnFillOn = False
        res3@cnLinesOn = True
        res3@cnLineLabelsOn = False
        res3@cnInfoLabelOn = False
        res3@cnLevels                   = (/20.,30.,40.,50./) ; ipol_1levs(1:)
        dum = gsn_csm_contour(wks,spol0,res3)
        ;overlay(plot(ip),dum)
        delete(res3@cnLevels)
      end if

      pres                  = True
      pres@gsFillColor      = "gray"
      pres@gsLineThicknessF = 2.0          ; twice as thick
      dum0(ip) = gsn_add_polygon(wks,plot(ip),xpoly,ypoly,pres)
      dum1(ip) = gsn_add_polyline(wks,plot(ip),xpoly(0:sdims(1)-1),ypoly(0:sdims(1)-1),pres)
      pres@gsLineColor = "gray"
      dum2(ip) = gsn_add_polyline(wks,plot(ip),xpoly1,ypoly1,pres)
      dum3(ip) = gsn_add_polyline(wks,plot(ip),xpoly1,ypoly2,pres)
      delete(pres)
   end if 
   end if   ;frhi
 end do
 end do

; Add annotation
 txres               = True                      ; text mods desired
 txres@txFontHeightF = 0.016                    ; text font height
 txres@txJust        = "CenterLeft"              ; Default is "CenterCenter".

;CPS  gsn_text_ndc(wks,"a)",.05,.77,txres)
;CPS  gsn_text_ndc(wks,"b)",.05,.58,txres)
;CPS  gsn_text_ndc(wks,"c)",.05,.40,txres)  ;0.39

 
; just draw ZDR c-s oknly
 draw(plot(9))
 frame(wks)

 resP = True
 resP@gsnMaximize    = True 
 ;resP@gsnPanelXF = (/-1,0.36,0.64,-1,0.36,0.64,-1,0.36,0.64/)
 resP@gsnPanelXF = (/-1,0.26,0.46,.66,-1,0.26,0.46,.66,-1,0.26,0.46,.66/)  ;for 4 columns  
 ;resP@gsnPanelDebug = True
; gsn_panel(wks,plot,(/nz,npol/), resP)

end
