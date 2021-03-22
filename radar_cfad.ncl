;All variables preserve their native units
;scale/offset is only used for plotting
; Cumulative Frequency Altitude Diagram
; Based on Yuter and Houze 1995 (YH95)
;P. Shrestha
;Nov 18 2019
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$PROM/functions/cfads.ncl"
load "$PROM/functions/read_radarnc.ncl"
load "$PROM/functions/load_cfadColors.ncl"

function record_append(fns,fne,fils_var,vid,glon_idB,glat_idB)
 begin
 do nf = fns,fne
   fname = fils_var(nf)
   vara  = read_radarnc(fname,vid, glon_idB,glat_idB)
   if (nf.eq.fns) then
     ndim = dimsizes(vara)
     var  = new((/fne-fns+1,ndim(0),ndim(1),ndim(2)/),typeof(vara))
     copy_VarAtts(vara,var)
   end if
   var(nf-fns,:,:,:) = vara
 end do
 print("Cropping data................")
 varc = var(:,:,:,0:199)
 delete(var)
 var  = varc
 return(var)
end


begin

;User Settings
;---------------------------------------------------------------------------
 ydate        = "20150705"
 diri         = "/daten01/z4/database/PROM/" + ydate + "/GRIDS/"
 fn           = 18  ;15:30 UTC  ;17
 fns          = fn  - 6
 fne          = fn  + 2

 wks_type     = "png"

 glon_idB    = 7.071663 ; Bonn Radar
 glat_idB    = 50.73052 ;

;Read netcdf files
;---------------------------------------------------------------------------

 fils_var     = systemfunc("ls " + diri + "/*.nc")
 ndimf        = dimsizes(fils_var)

 vid = 0
 var =  record_append(fns,fne,fils_var,vid,glon_idB,glat_idB)

; CFAD plot parameters
;--------------------------------------------------------
 opt          = True
 opt@YH95     = False  ;scale CFAD and remove data with low % of N
 opt@vid      = vid
 colors       = load_cfadColors(opt)
 minv         = colors@minv
 maxv         = colors@maxv
 binw         = colors@binw
 nz           = dimsizes(var&Z)
 ny           = dimsizes(var&Y)
 nx           = dimsizes(var&X)

 zhh_cfad = cfads(var&Z,var,minv(vid),maxv(vid),binw(vid),var@long_name,var@units,opt@YH95)

 vid      = 1
 opt@vid  = vid
 var =  record_append(fns,fne,fils_var,opt@vid,glon_idB,glat_idB)
 zdr_cfad = cfads(var&Z,var,minv(vid),maxv(vid),binw(vid),var@long_name,var@units,opt@YH95)

 vid      = 2
 opt@vid  = vid 
 var =  record_append(fns,fne,fils_var,opt@vid,glon_idB,glat_idB)
 kdp_cfad = cfads(var&Z,var,minv(vid),maxv(vid),binw(vid),var@long_name,var@units,opt@YH95)

 vid      = 3
 opt@vid  = vid
 var =  record_append(fns,fne,fils_var,opt@vid,glon_idB,glat_idB)
 rho_cfad = cfads(var&Z,var,minv(vid),maxv(vid),binw(vid),var@long_name,var@units,opt@YH95)


;Plot
;---------------------------------------------------------------------------
 plot = new(4,"graphic")

 fname_out = "cfad_radar_" + ydate + "_" + fn +  ".nc"

 wks = gsn_open_wks(wks_type,fname_out)
 colors =(/"white","black","white","gray30","gray40","gray50","gray61","gray75", \
          "gray80","gray90","paleturquoise4","paleturquoise3","orchid4","orchid3","orchid2","orchid1","orchid", \
         "royalblue4","darkturquoise","yellow"/)

 gsn_define_colormap(wks,colors)

 fncps_no = ":F21:"
 fncps_no = 21
 fncps_hgt = 0.025 ;2

 res                            = True           

 ;gsn Resources
 ;---------------------------------------------------
 res@gsnDraw                    = False
 res@gsnFrame                   = False
 res@gsnMaximize                = False
 res@gsnLeftStringFontHeightF   = fncps_hgt
 res@gsnRightStringFontHeightF  = fncps_hgt
 res@gsnLeftString              = "BoxPol " ;+ timestr
; res@gsnRightString             = zhh_cfad@long_name

 ;stretched grid Resources
 ;---------------------------------------------------
 res@trGridType            = "TriangularMesh"
 res@trYMaxF               = 15.0 ;km
; res@sfXArray              = zhh_cfad&binval
; res@sfYArray              = zhh_cfad&height  ;km 

 ;tm Resources
 ;---------------------------------------------------
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
 ;---------------------------------------------------
 res@tiXAxisFont                = fncps_no
 res@tiYAxisFont                = fncps_no
 res@tiXAxisFontHeightF         = fncps_hgt
 res@tiYAxisFontHeightF         = fncps_hgt
 res@tiXAxisString              = zhh_cfad@units
 res@tiYAxisString              = "Height [km]"

 ;cn Resources
 ;---------------------------------------------------
 res@cnFillOn                   = True
 ;res@cnFillMode                 = "RasterFill"
 res@cnLinesOn                  = False
 res@cnInfoLabelOn              = False   ; Turn off annotations
 res@cnLineLabelsOn             = False
 res@cnLevelSelectionMode       = "ExplicitLevels"
 res@cnLevels                   = colors@cnLevels

 ;ln Resources
 ;---------------------------------------------------
 res@lbLabelBarOn               = False
 res@lbLabelAutoStride          = True
 res@lbOrientation              = "vertical"
 res@lbTitlePosition            = "Right"
 res@lbTitleString              = "%"
 res@lbTitleExtentF             = 0.15
 res@lbTitleAngleF              = 90
 res@lbLabelFontHeightF         = 0.02
 res@pmLabelBarWidthF           = 0.05
 ;res@pmLabelBarHeightF          = 0.35

 ;vp Resources
 ;---------------------------------------------------
 res@vpHeightF                = 0.6
 res@vpWidthF                 = 0.45
; res@vpYF                     = 0.9
; res@vpXF                     = 0.15

 ;plot specific
 ;---------------------------------------------------
 var_cfad                  = zhh_cfad
 res@gsnRightString        = var_cfad@long_name
 res@sfXArray              = var_cfad&binval
 res@sfYArray              = var_cfad&height  ;km
 res@tiXAxisString         = var_cfad@units
 plot(0)     = gsn_csm_contour(wks,var_cfad,res)
 res0        = res
 res0@cnFillOn                   = False
 res0@cnLinesOn                  = True
 res0@cnLineThicknessF           = 3.0
 delete(res0@cnLevels)
 res0@cnLevels                   = 5.
 dummy0   = gsn_csm_contour(wks,zhh_cfad,res0)
 overlay(plot(0),dummy0)

;----------------------------------------------
 res@gsnLeftString  = ""
 res@tiYAxisString  = ""
 res@tmYLLabelsOn   = False

 delete(var_cfad)
 delete(res@sfXArray)

 var_cfad                  = zdr_cfad
 res@gsnRightString        = var_cfad@long_name
 res@sfXArray              = var_cfad&binval
 res@sfYArray              = var_cfad&height  ;km
 res@tiXAxisString         = var_cfad@units
 plot(1)     = gsn_csm_contour(wks,var_cfad,res)

 delete(var_cfad)
 delete(res@sfXArray)
 var_cfad                  = kdp_cfad
 res@gsnRightString        = var_cfad@long_name
 res@sfXArray              = var_cfad&binval
 res@sfYArray              = var_cfad&height  ;km
 res@tiXAxisString         = var_cfad@units
 plot(2)     = gsn_csm_contour(wks,var_cfad,res)

 delete(var_cfad)
 delete(res@sfXArray)
 var_cfad                  = rho_cfad
 res@gsnRightString        = var_cfad@long_name
 res@sfXArray              = var_cfad&binval
 res@sfYArray              = var_cfad&height  ;km
 res@tiXAxisString         = var_cfad@units
 res@lbLabelBarOn          = True
 plot(3)     = gsn_csm_contour(wks,var_cfad,res)

 resP = True
 resP@gsnMaximize    = True
 resP@gsnPanelXF = (/-1,0.28,0.51,0.74/)
 resP@gsnPanelDebug = True
 gsn_panel(wks,plot,(/1,4/), resP)

end
