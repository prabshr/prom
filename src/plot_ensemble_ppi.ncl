;Plot the plan veiw of one polarimetric variables for all ensemble members
; at Model Level

;P. Shrestha
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
external EXT1 "$prom_root/functions/readdim.so"
external EXT2 "$prom_root/functions/readascii.so"

begin
 ydates    = (/"20150705","20160513","20170706"/)
 diri_dBz    = "/daten04/pshrestha/shared-data/bonnRadar_" + ydates(pid) + "/Synth/"
 itype_gscp  = "HET2483"
 pvars       = (/"Z~B~H~N~","Z~B~DR~N~","K~B~DP~N~","~F8~r~B~~F10~hv~E~"/)
 punits      = (/"[dBZ]","[dB]","[deg/km]", "[-]"/)
 izlev       = 44 ; (/25, 44, 54/) 
 ipol        = 0   ;select pvars
 nens        = 20
; fn          = 0  ;time
 ipol_1levs  = ispan(5,65,5)*1. -10.
 ipol_2levs  = (/-1.,0,0.1,0.2,0.3,0.4,.5,0.6,0.8,1.0,2.,3.0,4./)
 ipol_3levs  = (/-0.5,-0.1, 0, 0.05, 0.1, 0.2, 0.3, 0.4,0.6, 0.8, 1.0, 2.0, 3.0/)
 ipol_4levs  = (/0.8, 0.85, 0.9, 0.92, 0.94, 0.95, 0.96,0.97, 0.98, 0.985, 0.99, 0.995, 1.0/)
 
;;;;;;;;;;

 do inst = 0, nens-1   ;BIG LOOP
   diri_input  = diri_dBz + itype_gscp + "_" + inst + "/input"
   diri_output = diri_dBz + itype_gscp + "_" + inst + "/output"

   if (inst.eq.0) then
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
     varPV        = new((/8,array(2),array(1),array(0)/),"float")
     varPV@_FillValue = -99

     ;
     varh         = new((/4,array(2),array(1),array(0)/),"float")    ;Initialize 
     varh@_FillValue = -99
     ;
     var_ens     = new((/nens,array(1),array(0)/),"float")    ;Initialize   
     var_ens@_FillValue = -99

   end if

   fils_var    = systemfunc("ls " + diri_output + "/dbz-tm-1m_*")
   fname       = fils_var(fn)
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

   var_ens(inst,:,:) = varh(ipol,izlev,:,:)
 
 end do      ;BIG LOOP inst

 ;Variable to plot
 var_ens@lat2d    = lat_c
 var_ens@lon2d    = lon_c

 print("Plotting ...")

 fncps_no = ":F21:"
 fncps_no = 21
 fncps_hgt = 0.035

 plot     = new(nens,"graphic")
 res      = True
 wks      = gsn_open_wks("png", "Figure_" + sprinti("%0.2i",fn))

 res      = True
 ;gsn Resources
 ;---------------
 res@gsnDraw                    = False
 res@gsnFrame                   = False
 res@gsnMaximize                = False

 ;pm resources
 ;--------------------
 res@pmTickMarkDisplayMode      = "Always" 
 res@pmLabelBarOrthogonalPosF   = 0.3
 ;res@pmLabelBarParallelPosF     = 0.2

 ;lb resources
 ;------------------
 res@lbLabelBarOn               = False

 ;tm Resources
 ;--------------- 
 res@tmXBLabelFont              = fncps_no
 res@tmYLLabelFont              = fncps_no
 res@tmYLLabelFontHeightF       = fncps_hgt
 res@tmXBLabelFontHeightF       = fncps_hgt
 res@tmXBOn = False
 res@tmYLOn = False
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

 ;mp Resources
 ;----------------------
 res@mpLimitMode                  = "Corners"           ; method to zoom
 res@mpLeftCornerLatF             = min(lat_c) ;mpcorner(0) ;28; lat(0,0)
 res@mpLeftCornerLonF             = min(lon_c) ;mpcorner(1) ;76 ;lon(0,0)
 res@mpRightCornerLatF            = max(lat_c) ;mpcorner(2) ;32 ;lat(ndim(0)-1,ndim(1)-1)
 res@mpRightCornerLonF            = max(lon_c) ;mpcorner(3) ;82 ;lon(ndim(0)-1,ndim(1)-1)
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

 do ip = 0, nens-1
   plot(ip)  = gsn_csm_contour_map(wks,var_ens(ip,:,:) ,res)
 end do

 resP = True
 resP@gsnMaximize    = True 
 resP@gsnPanelLabelBar    = True                ; add common colorbar
 resP@lbLabelFontHeightF  = 0.007    
 resP@gsnPanelMainString  = "Time " + fn + " | " + pvars(ipol) + "  " + punits(ipol)
 ;resP@gsnPanelDebug = True
 gsn_panel(wks,plot,(/4,5/), resP)

end
