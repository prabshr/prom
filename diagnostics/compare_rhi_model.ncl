;Plot the c-s veiw of polarimetric variables along Y-axis 
;a) specified model grid in X using linearized plot
;b) Interpolated data at specified model grid in X

;P. Shrestha
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$PROM/functions/shade_terrain.ncl"

external EXT1 "$prom_root/functions/readdim.so"
external EXT2 "$prom_root/functions/readascii.so"

begin
 diri_dBz    = "/daten04/pshrestha/shared-data/bonnRadar_20150705/Synth/HET2483_8/"
 pvars       = (/"Zh","Zdr","KDP","RHV"/)
 punits      = (/"dBZ","dB","deg/km ", "-"/)
 ipol        = 1     ;polarimetric moment to plot
 ;ix          = 95   ;location for C-S
 fn          = 10;20    ;file to read

 ipol_1levs  = ispan(5,65,5)*1. -10.
 ipol_2levs  = (/-1.,0,0.1,0.2,0.3,0.4,.5,0.6,0.8,1.0,2.,3.0,4./)
 ipol_3levs  = (/-0.5,-0.1, 0, 0.05, 0.1, 0.2, 0.3, 0.4,0.6, 0.8, 1.0, 2.0, 3.0/)
 ipol_4levs  = (/0.8, 0.85, 0.9, 0.92, 0.94, 0.95, 0.96,0.97, 0.98, 0.985, 0.99, 0.995, 1.0/)
 
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
 
 ;Average height at iz
 ;Need to reverse the vertical levels
 hgtPV  = hgtPV(:,::-1,:,:) 

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
 spol          = varh(ipol,:,:,ix)         ;
 sdims         = dimsizes(spol)
 spol!0        = "Z"
 spol!1        = "Y"
 spol&Y        = ispan(0,array(1)-1,1) ; ~ km
 xval2x        = conform_dims(sdims,spol&Y,1)*1.
 zval2x        = hgtPV(0,:,:,ix)

; for terrain shading
 xypoly  = shade_terrain(xval2x,zval2x)
 xpoly   = xypoly(0,:)
 ypoly   = xypoly(1,:)

 ;Now, interpolate also the variable to average model mass point height
 print("Interpolating " + pvars(ipol) + " data at ix " + ix )
 plane = (/ ix, 0/) 
 spol_ip       = wrf_user_intrp3d(varh(ipol,:,:,:),hgtPV(0,:,:,:),"v",plane,0.,False)  ;N-S c-s at ix
 spol_ip@_FillValue = varh@_FillValue
 ndimz            = dimsizes(spol_ip)
 maxhgt           = max(hgtPV(0,:,:,:))
 minhgt           = min(hgtPV(0,:,:,:))
 spol_ip!0        = "Z"
 spol_ip&Z        = fspan(minhgt,maxhgt,ndimz(0))
 spol_ip!1        = "Y"
 spol_ip&Y        = ispan(1,array(1)-1,1) ; ~km  0 is missing bug in wrf_interp

 print("Plotting ...")

;Differences
 ;print("RMSD " + rms)

 plot     = new(2,"graphic")
 res      = True
 wks      = gsn_open_wks("x11", "Figure")

 res@gsnDraw    = False
 res@gsnFrame   = False
 res@cnFillOn   = True
 res@cnFillMode = "CellFill"
 res@cnLinesOn  = False
 res@cnLevelSelectionMode       = "ExplicitLevels"

 res@trXMinF    = 60
 res@trXMaxF    = 140
 res@trYMaxF    = 15000

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

 res@tiMainString = "Interpolated C-S "
 plot(1)  = gsn_csm_contour(wks,spol_ip,res)

 ; Trick to linearize the plot
 res@trGridType = "TriangularMesh"
 res@sfXArray   = xval2x
 res@sfYArray   = zval2x
 res@tiMainString = "C-S at ix " + ix
 plot(0)  = gsn_csm_contour(wks,spol,res)

 pres                  = True
 pres@gsFillColor      = "gray"
 pres@gsLineThicknessF = 2.0          ; twice as thick
; pres@gsFillIndex      = 17           ; stipple, if desired

 dum0 = gsn_add_polygon(wks,plot(0),xpoly,ypoly,pres)
 dum1 = gsn_add_polyline(wks,plot(0),xpoly(0:sdims(1)-1),ypoly(0:sdims(1)-1),pres)



 gsn_panel(wks,plot,(/1,2/),False)


end
