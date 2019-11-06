;======================================================
; Diagnose vertical profiles
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
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"

begin

;-----------------------------------------------------------------------------------------
; USER SETTINGS
;-----------------------------------------------------------------------------------------

 diri        = "/p/scratch/chbn33/hbn331"
 event       = "bonnRadar_2014111500"
 case        = "graupelhail2test4"
 subevent    = "/"

 diri_cos    = diri + "/" + event + "/" + case + "/" + subevent + "/cosout/"
 diri_cosI   = diri + "/" + event + "/" + case + "/" + subevent + "/ivr/"

 ;Bonn Radar
 glon_id = (/7.071663/) ;
 glat_id = (/50.73052/) ;

 gloc          = 0       ;Location for vertical Profile
 filindex      = ispan(0,24,1)   ;Files to select for plotting

;-----------------------------------------------------------------------------------------
; READ INVARIANT FIELDS
;-------------------------------------------------------------------------------------
 nc_cosI = addfile(diri_cosI +"lfff00000000c.nc","r")
 lat = nc_cosI->lat
 lon = nc_cosI->lon
 hsurf = nc_cosI->HSURF(0,:,:)

 nm = getind_latlon2d (lat,lon, glat_id(gloc), glon_id(gloc))
 ;iy = nm(0,1)
 ;ix = nm(0,0)
 
 hhl = nc_cosI->HHL(0,::-1,:,:)
 ndimhhl   = dimsizes(hhl)
 nz        = ndimhhl(0)-1 ;
 zmodel    = 0.5*(hhl(0:nz-1,iy,ix) + hhl(1:nz,iy,ix)) - hsurf(iy,ix)
 ;zmodel = hhl(:,iy,ix) - hsurf(iy,ix)
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
   dtime = calendar_decode2(time,-5)
   year  = dtime(0,0)
   month = dtime(0,1)
   day   = dtime(0,2)
   hh    = dtime(0,3)
   mm    = dtime(0,4)
   ss    = dtime(0,5)

   str_time  = sprinti("%0.4i",year) + "-" + sprinti("%0.2i",month) + \
               "-" + sprinti("%0.2i",day) + " " + sprinti("%0.2i",hh) + \
                sprinti("%0.2i",mm) ;+ " UTC" ;sprinti("%0.2i",ss) +  " UTC" 
   print("  " + fils_cos(fn))
   print(str_time + " ")
;------------------------------------------------------------------------------
;;Prognostic state variables 
;------------------------------------------------------------------------------
   varp         = nc_cos->$varname$(0,::-1,iy,ix)
   
;-----------------------------------------------------------------------------------------
; Plotting Resource
;-----------------------------------------------------------------------------------------
  wks_type = "x11"
  fname_out = "Figure_diagnostic_"+case + "_"+sprinti("%0.4i",fn)
  wks = gsn_open_wks(wks_type,fname_out)

  fncps_no = ":F21:"
  fncps_no = 21
  fncps_hgt = 0.02

;-----------------------------------------------------------------------------------------
; Pressure plot resources 
;-----------------------------------------------------------------------------------------

 res                            = True                 ; plot mods desired
 res@gsnDraw                    = False
 res@gsnFrame                   = False
 res@gsnLeftString              = str_time
 res@gsnRightString             = varp@standard_name + " [" + varp@units+"]"
 res@gsnMaximize                = False
 res@gsnLeftStringFontHeightF   = fncps_hgt
 res@gsnRightStringFontHeightF  = fncps_hgt
 
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

 ;Viewport
; res@vpXF                          = 0.15 ;0.1
 res@vpYF                          = 0.95
 res@vpHeightF                     = 0.75
 res@vpWidthF                      = 0.72

 res@xyLineThicknessF              = 3.
 res@xyMarkLineMode                = "MarkLines"
 res@xyMarker                      = 2

 res@tiYAxisString                 = "Height (km)"
 res@tiXAxisString                 = varp@units 

 if (varname.eq."T") then
 res@trXMinF                       = 260 
 res@trXMaxF                       = 290 
 end if
 res@trYMaxF                       = 10000 ;max(zmodel)
 plot                               = gsn_csm_xy (wks,varp,zmodel,res) ; create plot

 draw(plot)
 frame(wks)

 end do
end
