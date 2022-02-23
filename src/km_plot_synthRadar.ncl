;Plot the plan veiw of polarimetric variables at low and mid levels
; and cross-section (frhi = True) replaces ppi at high levels 

;P. Shrestha
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$PROM/functions/shade_terrain.ncl"
load "$PROM/functions/load_radarColors.ncl"
load "$prom_root/functions/read_emvorado_pol.ncl"
load "$prom_root/functions/generateXY.ncl"
external EXT1 "$prom_root/functions/readdim.so"
external EXT2 "$prom_root/functions/readascii.so"

begin
;-----------------------------------------------------------------------------------------
;#1 USER SETTINGS
;-----------------------------------------------------------------------------------------
 prom=4
 pname="EMPA2_PROM1";
 inst=6
 readBPRO = False ;True    ;
 lnd = "HET"
 itype = "ART"
 itype_gscp  = lnd + itype

 ;
 fn          = 18
 ix          = 191
 izlevs      = (/25, 44, 54/)
                ;y1 x1 y2   x2
 xycorner    = (/-10, -10, 110, 110/)
 ;

 diri        = "/p/largedata/hbn33/shrestha1/"
 sfold       = "runE_" + inst + "_" + itype_gscp
 expname     = lnd + "_" + pname
 subevent    = expname + "/" + itype_gscp + "/" + sfold

 diri_cos    = diri + "/" + "/" + subevent + "/cosout/"
 diri_cosI   = diri + "/" + "/" + subevent + "/cosout/ivr/"
 ;diri_dBz    = "/daten04/pshrestha/shared-data/bonnRadar_" + ydate + "/Synth/HET2483_" + inst +"/"
 diri_dBz    = diri + "/Synth/" + itype_gscp +"_" + inst+"/"
 diri_emvo   = "/daten/daten04/pshrestha/shared-data/test-emvo/" + pname + "/" + sfold
 diri_emvo   = "/p/largedata/hbn33/test-emvo/" + subevent
 dataFillVal    = -999        ;Polarimetric flags

 pvars       = (/"Z~B~H~N~","Z~B~DR~N~","K~B~DP~N~","~F8~r~B~~F10~hv~E~"/)
 punits      = (/"[dBZ]","[dB]","[deg/km]", "[-]"/)

 frhi        = True ;replace 3rd level ppi with rhi
 glon_idB    = 7.071663 ; Bonn Radar
 glat_idB    = 50.73052 ;
 dxy         = 1.132  ;km

;;;;;;;;;;
;-----------------------------------------------------------------------------------------
;#3  READ NC AND ASCII FILES
;-------------------------------------------------------------------------------------
 nc_cosI = addfile(diri_cosI +"lfff00000000c.nc","r")
 lat_c = nc_cosI->lat
 lon_c = nc_cosI->lon
 phgt = nc_cosI->HSURF(0,:,:)

 hhl       = nc_cosI->HHL(0,::-1,:,:)
 ndimhhl   = dimsizes(hhl)
 nz        = ndimhhl(0)-1
 ny        = ndimhhl(1)
 nx        = ndimhhl(2)
 mhgt      = 0.5*(hhl(0:nz-1,:,:) + hhl(1:nz,:,:)) 
;
 fils_cos  = systemfunc("ls " + diri_cos + "lf*")
 ;nfils   = dimsizes(fils_cos)
 ;nfils    = fne - fns + 1

 ;Generate X-Y coordinates in km
 XYkm = generateXY(lat_c,lon_c,glat_idB,glon_idB)
 Ykm  = XYkm(1,:,:)
 Xkm  = XYkm(0,:,:)

 ;----------------------------------------------------
 npol          = dimsizes(pvars)
 nzp            = dimsizes(izlevs)
 ;np            = (npol-1)*nzp           ;skip Zh
 np            = npol*nzp


 if (readBPRO) then

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

 else
  ;Read EMVORADO outputs ------------------------------------------------------------
   nh           = dimsizes(pvars)
   varh         = new((/nh,nz,ny,nx/),"float")    ;Initialize variable for plotting
   varh@_FillValue = -999.99
   opt      = True
   opt@profile = -1   ;return 3d data
   opt@nf      =  fn


   opt@var  = "zh"
   Zh = read_emvorado_pol(diri_emvo,opt)
   opt@var  = "zv"
   Zv = read_emvorado_pol(diri_emvo,opt)
   opt@var  = "kdp"
   KDP = read_emvorado_pol(diri_emvo,opt)
   opt@var  = "rrhv"
   rRHV = read_emvorado_pol(diri_emvo,opt)
   opt@var  = "irhv"
   iRHV = read_emvorado_pol(diri_emvo,opt)

   RHV = ((rRHV^2+iRHV^2)^0.5) / ((Zh*Zv)^0.5)   ;[-]
   KDP = KDP*1000.                               ;[deg/km]
   Zh  = 10*log10(Zh/1.)                         ;[dBZ]
   Zv  = 10*log10(Zv/1.)                         ;[dBZ]
 end if

 varh(0,:,:,:) = where(Zh.lt.-10., varh@_FillValue, Zh)
 varh(1,:,:,:) = where(Zh.lt.-5., varh@_FillValue, Zh-Zv)
 varh(2,:,:,:) = where(Zh.lt.-5., varh@_FillValue, KDP)
 varh(3,:,:,:) = where(Zh.lt.-5., varh@_FillValue, RHV)
 delete([/Zh, Zv, KDP, RHV/])

 varh!2 = "Y"
 varh!3 = "X"
 varh&Y = Ykm(:,0)
 varh&X = Xkm(0,:)

 
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
   dum_hy = new(np,"graphic")
   dum_hx = new(np,"graphic")
 end if


 res      = True
 wks      = gsn_open_wks("x11", "Figure_case_i" + inst + "_t" + fn + "_ix" + ix)

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

 res@tiYAxisString              = "Y Range [km]"
 res@tiXAxisString              = "X Range [km]"


 res3 = res   ;for frhi = True

 ;pm resources
 ;--------------------
 res@pmTickMarkDisplayMode      = "Always"
 res@pmLabelBarOrthogonalPosF   = 0.3
 ;res@pmLabelBarParallelPosF     = 0.2

 ;cn Resources
 ;-------------------------
 res@cnFillOn   = True
 res@cnFillMode = "CellFill"
 res@cnLinesOn  = False
 res@cnLevelSelectionMode       = "ExplicitLevels"
 res@cnInfoLabelOn    = False

 ;tr Resources
 ;-----------------------------
 res@trXMinF = xycorner(1)  
 res@trYMinF = xycorner(0) 
 res@trXMaxF = xycorner(3)
 res@trYMaxF = xycorner(2)

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

 do ik = 0, nzp-1

 iz = izlevs(ik)

 temp = mhgt(iz,:,:) - phgt ;.ag.l

 print(iz + " " + avg(temp))
 avghgt = round((avg(temp)/100.),3)*100
 print("Average model mass point height at iz =" +iz + " is " + avghgt + " m a.m.s.l")

 if (frhi) then
   nzz = nzp-2
 else
   nzz = nzp-1
 end if
 ;
 if (ik.eq.1) then
   res@tiXAxisString   = "X Range [km]"
 else
   res@tiXAxisString   = ""
 end if
 ;
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
     res@tiYAxisString = "Y Range [km]"
   else
     res@tiYAxisString = ""
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
   plot(ip)  = gsn_csm_contour(wks,varh(ipol,iz,:,:) ,res)

   if (ipol.eq.1) then
     res2@cnLevels                   = (/20.,30.,40.,50./) ; ipol_1levs(1:)
     dum      = gsn_csm_contour(wks,varh(0,iz,:,:),res2)
     ;overlay(plot(ip),dum)
   end if

   if (frhi) then
   if (ik.lt.nzp-1) then
     pres                  = True
     pres@gsLineColor      = "black"
     pres@gsLineDashPattern = 0
     ypoly_cs = (/xycorner(0),xycorner(2)/)
     xpoly_cs = (/Xkm(0,ix),Xkm(0,ix)/)
     dum_cs(ip) = gsn_add_polyline(wks,plot(ip),xpoly_cs,ypoly_cs,pres)
     pres@gsLineDashPattern = 2
     xpoly_cs = (/-150,150/)
     ypoly_cs = (/ 0, 0/)
     dum_hy(ip) = gsn_add_polyline(wks,plot(ip),ypoly_cs,xpoly_cs,pres)
     dum_hx(ip) = gsn_add_polyline(wks,plot(ip),xpoly_cs,ypoly_cs,pres)
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
   dumtxt = gsn_add_text(wks,plot(ip),"x",0.,0.,txres)

   ;dum      = gsn_csm_contour(wks,phgt,res2)
   ;overlay(plot(ip),dum)
   ;print(ip + " "  +ik + " " + ipol +" "+ res@tmXBLabelsOn)

   ;Here we replace 3rd level ppi with rhi
   ;----------------------------------------------------------
   if (frhi) then
   if (ik.eq.nzp-1 ) then
      spol          = (/varh(ipol,:,:,ix)/)         ;
      sdims         = dimsizes(spol)
      spol!0        = "Z"
      spol!1        = "Y"
      spol&Y        = Ykm(:,ix) ; km 
      spol0          = (/varh(0,:,:,ix)/)    
      spol0!0        = "Z"
      spol0!1        = "Y"
      spol0&Y        = Ykm(:,ix) ; km 

      xval2x        = conform_dims(sdims,spol&Y,1)*1.
      zval2x        = mhgt(:,:,ix)/1000.  ;km

      ; for terrain shading
      xypoly  = shade_terrain(xval2x,zval2x)
      xpoly   = xypoly(0,:)
      ypoly   = xypoly(1,:)
      xpoly1   = (/xycorner(0),xycorner(2)/)
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

      res3@trXMinF    = xycorner(0) 
      res3@trXMaxF    = xycorner(2) 
      res3@trYMaxF    = 15.

 ;CPS     res3@tmXBLabelsOn = False

      res3@lbOrientation = "Vertical"

      res3@tiXAxisString = "Y Range [km]" ;"Latitude [~F34~0~F~N]"
      if (ipol.eq.0) then
        res3@tmYLLabelsOn = True
        res3@tiYAxisString = "Height [km]"
      else
        res3@tiYAxisString = ""
        res3@tmYLLabelsOn = False
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

 gsn_text_ndc(wks,"a)",.005,.83,txres)
 gsn_text_ndc(wks,"b)",.005,0.62,txres)
 gsn_text_ndc(wks,"c)",.005,.36,txres)  ;0.39

 
 resP = True
 resP@gsnMaximize    = True 
 ;resP@gsnPanelXF = (/-1,0.36,0.64,-1,0.36,0.64,-1,0.36,0.64/)
 resP@gsnPanelXF = (/-1,0.26,0.46,.66,-1,0.26,0.46,.66,-1,0.26,0.46,.66/)  ;for 4 columns  
 resP@gsnPanelYF = (/-1,-1,-1,-1,0.61,0.61,0.61,0.61,0.35,0.35,0.35,0.35/)

 resP@gsnPanelDebug = True
 gsn_panel(wks,plot,(/nzp,npol/), resP)

end
