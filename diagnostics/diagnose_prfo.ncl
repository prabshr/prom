;Plot the plan veiw of polarimetric variables at 
;a) Model Level
;b) Interpolated data at average model mass point height (a.m.s.l) for that model level
; Also overlaid in plot a) is the spatial pattern of model mass point height at that level

;P. Shrestha
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$PROM/functions/shade_terrain.ncl"

external EXT1 "$prom_root/functions/readdim.so"
external EXT2 "$prom_root/functions/readascii.so"

begin
 diri_dBz    = "/daten04/pshrestha/shared-data/bonnRadar_20150705/Synth/"
 itype_gscp  = "HET2483"
 inst        = 5 
 pvars       = (/"Zh","Zdr","KDP","RHV"/)
 punits      = (/"dBZ","dB","deg/km ", "-"/)
 ipol        =  0     ;polarimetric moment to plot
; iz          = 21 ;21(low level) ; 47(mid level) , 41? ;model level to plot
; ix         = 95   ;location for C-S
 fn          = 10;20    ;file to read

 ipol_1levs  = ispan(5,65,5)*1. -10.
 ipol_2levs  = (/-1.,0,0.1,0.2,0.3,0.4,.5,0.6,0.8,1.0,2.,3.0,4./)
 ipol_3levs  = (/-0.5,-0.1, 0, 0.05, 0.1, 0.2, 0.3, 0.4,0.6, 0.8, 1.0, 2.0, 3.0/)
 ipol_4levs  = (/0.8, 0.85, 0.9, 0.92, 0.94, 0.95, 0.96,0.97, 0.98, 0.985, 0.99, 0.995, 1.0/)
 
;;;;;;;;;;

 diri_input  = diri_dBz + itype_gscp + "_" + inst +  "/input"
 diri_output = diri_dBz + itype_gscp + "_" + inst +  "/output"

 fils_var    = systemfunc("ls " + diri_output + "/dbz-tm-1m_*")
 fils_elv    = systemfunc("ls " + diri_input + "/elevation_*")
 fil_lat     = systemfunc("ls " + diri_input + "/lat_*")
 fil_lon     = systemfunc("ls " + diri_input + "/lon_*")

 ; Add a loop
 ;do fn = 10,15

 ;Read header information to extract nx,ny,nz
 array = new(3,integer)
 EXT1::readdim(array(0),array(1),array(2),fils_elv(fn))

 latPV = new((/1,1,array(1),array(0)/),"float")
 lonPV = new((/1,1,array(1),array(0)/),"float")
 hgtPV = new((/1,array(2),array(1),array(0)/),"float")
 EXT2::readascii(array(0),array(1),1,latPV,1,1,fil_lat)
 EXT2::readascii(array(0),array(1),1,lonPV,1,1,fil_lon)
 EXT2::readascii(array(0),array(1),array(2),hgtPV,1,1,fils_elv(fn))   ;Model mass point heights
 
 ;Average height at iz
 ;Need to reverse the vertical levels
 hgtPV  = hgtPV(:,::-1,:,:) 
 mhgt   = hgtPV(0,iz,:,:)
 avghgt = avg(mhgt)
 print("Average model mass point height at iz =" +iz + " is " + avghgt + " m a.g.l")

 nh           = dimsizes(pvars)
 varh         = new((/nh,array(2),array(1),array(0)/),"float")    ;Initialize variable for plotting
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
 varh(1,:,:,:) = where(Zh.lt.-5., varh@_FillValue, Zh-Zv)
 varh(2,:,:,:) = where(Zh.lt.-5., varh@_FillValue, KDP)
 varh(3,:,:,:) = where(Zh.lt.-5., varh@_FillValue, RHV)
 delete([/Zh, Zv, KDP, RHV/])

 ;Variable to plot
 spol          = varh(ipol,iz,:,:)         ;
 spol!0        = "Y"
 spol!1        = "X"
 spol&Y        = ispan(0,array(1)-1,1)  ;~ km
 spol&X        = ispan(0,array(0)-1,1)  ;~ km

 spol_ip       = varh(ipol,:,:,ix)           ;
 sdims         = dimsizes(spol_ip)
 spol_ip!0     = "Z"
 spol_ip!1     = "Y"
 spol_ip&Y     = ispan(0,array(1)-1,1) ; ~ km
 xval2x        = conform_dims(sdims,spol_ip&Y,1)*1.
 zval2x        = hgtPV(0,:,:,ix)
; for terrain shading
 xypoly  = shade_terrain(xval2x,zval2x)
 xpoly   = xypoly(0,:)
 ypoly   = xypoly(1,:)


 print("Plotting ...")

 plot     = new(2,"graphic")
 res      = True
 wks      = gsn_open_wks("x11", "Figure")

 res@gsnDraw    = False
 res@gsnFrame   = False
 res@cnFillOn   = True
 res@cnFillMode = "CellFill"
 res@cnLinesOn  = False
 res@cnLevelSelectionMode       = "ExplicitLevels"
 if (ipol.eq.0) then
   ipol_levs = ipol_1levs
 end if
 if (ipol.eq.1) then
   ipol_levs = ipol_2levs
 end if
 if (ipol.eq.2) then
   ipol_levs = ipol_3levs
 end if
 if (ipol.eq.3) then
   ipol_levs = ipol_4levs
 end if

 res@cnLevels                   =  ipol_levs

 res@trXMinF = 60
 res@trYMinF = 60
 res@trXMaxF = 140
 res@trYMaxF = 140

 res@tiMainString = "Model level " + iz
 plot(0)  = gsn_csm_contour(wks,spol,res)

 ;Overlay model height contours
 sres     = res
 delete(sres@cnLevelSelectionMode)
 delete(sres@cnLevels)
 sres@cnFillOn = False
 sres@cnLinesOn = True
 sres@cnLineLabelsOn = False
 sres@cnLineColor    = "Gray"
 dum      = gsn_csm_contour(wks,mhgt,sres)
 overlay(plot(0),dum)

 ; Trick to linearize the plot
 res@trGridType = "TriangularMesh"
 res@sfXArray   = xval2x
 res@sfYArray   = zval2x
 res@tiMainString = "C-S at ix " + ix
 delete(res@trYMinF)
 res@trYMinF    = 0.
 res@trYMaxF    = 15000
 
 plot(1)  = gsn_csm_contour(wks,spol_ip,res)

; Add shading of topography
 pres                  = True
 pres@gsFillColor      = "gray"
 pres@gsLineThicknessF = 2.0          ; twice as thick
; pres@gsFillIndex      = 17           ; stipple, if desired
 dum0 = gsn_add_polygon(wks,plot(1),xpoly,ypoly,pres)
 dum1 = gsn_add_polyline(wks,plot(1),xpoly(0:sdims(1)-1),ypoly(0:sdims(1)-1),pres)

; Add c-s lines
 xppi = (/ix,ix/)
 yppi = (/0, array(1)/)
 dum2 = gsn_add_polyline(wks,plot(0),xppi,yppi,pres) 

 delete([/xppi,yppi/])
 xppi = ispan(0,array(1)-1,1)
 yppi = hgtPV(0,iz,:,ix)
 ;might looks confusing with melting layer T, so maybe remove it
 dum3 = gsn_add_polyline(wks,plot(1),xppi,yppi,pres)
 delete([/xppi,yppi/])

 gsn_panel(wks,plot,(/1,2/),False)

 delete([/res,spol,spol_ip/])

; end do

end
