;======================================================
; Diagnose vertical profiles
; COSMO outputs M1,i specific hydrometeor content (kg/kg)
;               M0,i specific hydrometeor number  (kg-1)
; Here we plot,     first moment as [g/kg], opt = 0
;                    zero moment as [m-3] , opt = 1
;  D = ax^b mean volume diameter as [m]   , opt = 2
; for opt = 0, air density is also computed
;
;======================================================
; Prabhakar Shrestha
; E-mail: pshrestha@uni-bonn.de
;=======================================================

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$PROM/functions/calc_rho_nclf.ncl"
load "$PROM/functions/melt_layer_height_nclf.ncl"
load "$PROM/functions/cosmo_time_nclf.ncl"

begin

;-----------------------------------------------------------------------------------------
; USER SETTINGS
;-----------------------------------------------------------------------------------------
 ;opt         = 0  ; [0 - zero moment, 1- first moment]
 diri        = "/p/scratch/chbn33/hbn331"
 event       = "bonnRadar_2014111500"
 case        = "graupelhail2test4"
 subevent    = "HET2473"

 diri_cos    = diri + "/" + event + "/" + case + "/" + subevent + "/cosout/"
 diri_cosI   = diri + "/" + event + "/" + case + "/" + subevent + "/ivr/"

 ;
 hydro_colors = (/"green","blue","skyblue","purple","black","orange"/)
 hydro_dash   = (/0, 0, 1, 1, 1, 1/)
 hydro_thickness = hydro_dash + 4. 
 hmax         = 10000   ;m     ;for vertical extent of plot

 ;Bonn Radar
 glon_id = (/7.071663/) ;
 glat_id = (/50.73052/) ;

 gloc          = 0       ;Location for vertical Profile
 filindex      = ispan(0,24,1)   ;files to plot

;-----------------------------------------------------------------------------------------
; Set variables to read and specify scales or offsets
;-------------------------------------------------------------------------------------
 hydrometeors_qx = (/"QC","QR","QI","QS","QG","QH"/)
 hydrometeors_qn = (/"NCCLOUD","NCRAIN","NCICE","NCSNOW","NCGRAUPEL","NCHAIL"/)

 if (opt.eq.0) then
   hydrometeors = hydrometeors_qn
   scale_varh   = 1/1000.
   qmax         = 1E7
   hydunits     = "L~S~-1~N~"
   hydname      = "Number Concentration"
else if (opt.eq.1) then
   hydrometeors = hydrometeors_qx
   scale_varh   = 1000.  ;kg/kg to g/kg
   qmax         = 0.8     ;g/kg  ;for horizontal extent of plot
   hydunits     = "g/m~S~3~N~"
   hydname      = "Sp. Hydrometeor Content"
else if (opt.eq.2) then
   x3           = 1./3.
   hydrometeors = hydrometeors_qx
   a_geo        = (/0.124    , 0.124  , 0.835 , 2.4    , 0.142     , 0.1366/)  ;mkg-b
   b_geo        = (/x3       , x3     , 0.39  , 0.455  , 0.314     , x3/)      ;-
   qmax         = 15.
   scale_varh   = 1000.
   hydunits     = "mm"
   hydname      = "Mean Vol. Hydrometeor Diameter"
 else
   print("code not written for opt = " + opt)
   exit
 end if
 end if
 end if
;
;-----------------------------------------------------------------------------------------
; READ INVARIANT FIELDS
;-------------------------------------------------------------------------------------

 nc_cosI = addfile(diri_cosI +"lfff00000000c.nc","r")
 lat = nc_cosI->lat
 lon = nc_cosI->lon
 hsurf = nc_cosI->HSURF(0,:,:)

 nm = getind_latlon2d (lat,lon, glat_id(gloc), glon_id(gloc))
 iy = nm(0,0)
 ix = nm(0,1)
 print(glat_id + " " + glon_id)
 print(iy + " " + ix)
 
 hhl = nc_cosI->HHL(0,::-1,:,:)
 ndimhhl   = dimsizes(hhl)
 nz        = ndimhhl(0)-1 ;
 zmodel    = 0.5*(hhl(0:nz-1,iy,ix) + hhl(1:nz,iy,ix)) ;- hsurf(iy,ix)
 ndim = dimsizes(lat)

;-----------------------------------------------------------------------------------------
; READ NC FILES 
;-----------------------------------------------------------------------------------------
 fils_cos  = systemfunc("ls " + diri_cos + "lf*")
 nf_cos   = dimsizes(fils_cos)
;Initialize empty arrays variables
; nf_ctr = dimsizes(filindex) 
 nf_ctr = dimsizes(fils_cos)
;
 do fnctr = 0,nf_ctr-1   ;;BIG LOOP OVER FILES 
   
   fn = fnctr ;filindex(fnctr)
   nc_cos = addfile(fils_cos(fn),"r")
   time   = nc_cos->time
   ;
   str_time = cosmo_time_nclf(time)
 
;------------------------------------------------------------------------------
;;Prognostic state variables
;------------------------------------------------------------------------------
   tair             = nc_cos->T(0,::-1,iy,ix)     ;K
   if (opt.ne.2) then
     pressure         = nc_cos->P(0,::-1,iy,ix)     ;Pa
     qv               = nc_cos->QV(0,::-1,iy,ix)    ;kg/kg
     qc               = nc_cos->QC(0,::-1,iy,ix)    ;kg/kg
     qr               = nc_cos->QR(0,::-1,iy,ix)    ;kg/kg
     qi               = nc_cos->QI(0,::-1,iy,ix)    ;kg/kg
     qs               = nc_cos->QS(0,::-1,iy,ix)    ;kg/kg
     qg               = nc_cos->QG(0,::-1,iy,ix)    ;kg/kg
     qh               = nc_cos->QH(0,::-1,iy,ix)    ;kg/kg
   end if
   if (opt.ne.2) then
     qrs              = qr + qi + qs + qg + qh      ;Precipitation loading
     rho              = calc_rho_nclf(tair,pressure,qv,qc,qrs) ;air density (kgm-3)
   end if

   nh           = dimsizes(hydrometeors)
   varh         = new((/nh,nz/),"float")    ;Initialize variable for plotting
   varh@units   = hydunits
   varh@long_name = hydname

   do ih = 0, nh-1
     if (opt.eq.0) then
       varh(ih,:) = (/nc_cos->$(hydrometeors_qn(ih))$(0,::-1,iy,ix)/)
       varh(ih,:) = (/scale_varh * varh(ih,:) * rho/)    ;(kg-1) to (L-1)
           ;Nt    = 1/1000. * Nd * rho            ; number density to number concentration
     else if (opt.eq.1) then
       varh(ih,:) = (/nc_cos->$(hydrometeors_qx(ih))$(0,::-1,iy,ix)/)
       varh(ih,:) = (/scale_varh * varh(ih,:) * rho/)    ;(kg/kg) to (g/m3)
           ;Qt    = 1000. * Qd * rho 
     else if (opt.eq.2) then
       qtemp      = (/nc_cos->$(hydrometeors_qx(ih))$(0,::-1,iy,ix)/)
       Ntemp      = (/nc_cos->$(hydrometeors_qn(ih))$(0,::-1,iy,ix)/)
       Ntemp(ind(Ntemp.eq.0)) = varh@_FillValue
       varh(ih,:) = (/qtemp/Ntemp/)
       varh(ih,:) = (/scale_varh* (a_geo(ih) * (varh(ih,:))^b_geo(ih))/)
     end if
     end if
     end if
   end do 

; Set 0 as missing value
   varh1d  = ndtooned(varh)
   varh0   = ind(varh1d.eq.0)
   if (.not.all(ismissing(varh0))) then
     varh1d(varh0) = varh@_FillValue
   end if
   varh    = onedtond(varh1d,(/nh,nz/))
   delete(varh0)
;

;---Melt layer height -----
   zmelt       = melt_layer_height_nclf(tair,zmodel)
   ;For melting layer plot
   qx = fspan(1E-6,qmax,10)     ;g/kt
   freezing = qx*0.0 + zmelt ; m 

;-----------------------------------------------------------------------------------------
; Plotting Resource
;-----------------------------------------------------------------------------------------
  wks_type = "x11"
  fname_out = "Figure_diagnostic_"+ opt+"_"+ case + "_"+sprinti("%0.4i",fn)
  wks = gsn_open_wks(wks_type,fname_out)

  fncps_no = ":F21:"
  fncps_no = 21
  fncps_hgt = 0.02

;-----------------------------------------------------------------------------------------
; Pressure plot resources 
;-----------------------------------------------------------------------------------------

 res                            = True                 ; plot mods desired

 ;gsn Resources
 ;-------------------------------------------
 res@gsnDraw                    = False
 res@gsnFrame                   = False
 res@gsnLeftString              = str_time
 res@gsnRightString             = "Profiles at " + lat(iy,ix) + " N ," + lon(iy,ix) + " E" 
 res@gsnMaximize                = False
 res@gsnLeftStringFontHeightF   = fncps_hgt
 res@gsnRightStringFontHeightF  = fncps_hgt

 ;tm Resources
 ;------------------------------------------- 
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
 ;-------------------------------------------
 res@tiXAxisFont                = fncps_no
 res@tiYAxisFont                = fncps_no
 res@tiXAxisFontHeightF         = fncps_hgt
 res@tiYAxisFontHeightF         = fncps_hgt

 ;vp Resources
 ;-------------------------------------------
; res@vpXF                          = 0.15 ;0.1
 res@vpYF                          = 0.95
 res@vpHeightF                     = 0.75
 res@vpWidthF                      = 0.72

 ;xy Resources
 ;-------------------------------------------
 res@xyLineThicknesses             = hydro_thickness 
 res@xyMarkLineMode                = "MarkLines"
 res@xyLineColors                = hydro_colors
 res@xyDashPatterns            = hydro_dash
 res@xyMarker                  = 14
 res@xyMarkerColors            = hydro_colors 
 res@xyMarkerSizeF            = 7. ;hydro_colors

 ;ti Resources
 ;-------------------------------------------
 res@tiYAxisString                 = "Height a.m.s.l. (km)"
 res@tiXAxisString                 = varh@long_name + " (" + varh@units + ")"

 ;tr Resources
 ;-------------------------------------------
 res@trYMaxF                       = hmax
 res@trXMaxF                       = qmax
 if (opt.eq.0) then
   res@trXMinF                     = 1.0
   res@trXLog                      = True
 end if

 plot                              = gsn_csm_xy (wks,varh,zmodel,res) ; create plot

 resV                              = True
 resV@gsnDraw                      = False
 resV@gsnFrame                     = False
 resV@xyLineColor                  = "gray"
 resV@xyMarkLineMode               = "Lines"
 resV@xyDashPattern                = 1
 dum                               = gsn_csm_xy(wks,qx,freezing,resV)
 overlay(plot,dum)

 ;Add legend
 lgres                    = True
 lgres@lgLineColors       = hydro_colors
 lgres@lgLineThicknessF   = hydro_thickness
 lgres@lgDashIndexes      = hydro_dash
 lgres@lgItemType         = "Lines"
 lgres@lgLabelFontHeightF = .1
 lgres@vpWidthF           = 0.1
 lgres@vpHeightF          = 0.18
 lgres@lgPerimOn          = False
 lbid = gsn_create_legend(wks,nh,hydrometeors,lgres)

 amres = True
 amres@amParallelPosF   = 0.25
 amres@amOrthogonalPosF = -0.35
 annoid1 = gsn_add_annotation(plot,lbid,amres)

 draw(plot)
 frame(wks)

 delete([/varh,res,resV/])

 end do
 system("rm -r qxDiagnostics")
 system("mkdir qxDiagnostics")
 system("mv Figure_diagnostic_*.png qxDiagnostics")
 system("convert -delay 100 qxDiagnostics/Figure*.png "+"Figure_diagnostic_"+subevent+"_"+opt +".gif")

end
