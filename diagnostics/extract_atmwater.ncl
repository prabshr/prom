load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$TPS/functions/calcrho_cos_nclf.ncl"

;Computes Precipitatble Water, LWP and IWP

begin

;-----------------------------------------------------------------------------------------
; USER SETTINGS
;-----------------------------------------------------------------------------------------
 pid         = 3 
 ;inst        = 6
 lnd         = "HET"
 itype       = "ART" 
 casenames   = (/"bonnRadar_20150705"/)
 pnames      = (/"PROM1","PROM2","PROM3","EMPA2_PROM1"/)
;;;
 ;casename    = casenames(pid)
 pname       = pnames(pid)
 itype_gscp  = lnd + itype
 sfold       = "runE_" + inst + "_" + itype_gscp
 expname     = lnd + "_" + pname
 subevent    = expname + "/" + itype_gscp + "/" + sfold

 ;CLUMA2
 ;diri_cosmo  = "/daten04/pshrestha/shared-data/" + casename + "/" + expname + "/" + sfold
 diri_spath  = "/daten/daten04/pshrestha/shared-data/storm_mask/nc"

 ;JUWELS
 diri_spath  = "/p/largedata/hbn33/storm_mask/" + expname + "/" + itype_gscp + "/nc/"
 diri_cosmo  = "/p/largedata/hbn33/shrestha1/" + subevent

 ;Bonn Radar
 glon_id = (/7.071663/) ;, 7.098308/) ;
 glat_id = (/50.73052/) ;, 50.737497/) ;
 gmarker = (/10,9/)
 pressure    = (/500.,700./) ; hPa for wind vectors
 fplot       = False

;-----------------------------------------------------------------------------------------
; SELECT FILES and EXTENTS for Plots
;-----------------------------------------------------------------------------------------
 filindex      = ispan(0,1,1)
 ;dext          = 1.5
 ;mpcorner      = (/glat_id(0)-dext,glon_id(0)-dext,glat_id(0)+dext,glon_id(0)+dext/)
 nbdl          = 3 + 17  ;exclude boundaries and relaxation zone

;-----------------------------------------------------------------------------------------
; READ INVARIANT FIELDS
;-------------------------------------------------------------------------------------
 fil_ivr = diri_cosmo + "/cosout/ivr/lfff00000000c.nc"
 nc_cosI = addfile(fil_ivr,"r")
 lat     = nc_cosI->lat
 lon     = nc_cosI->lon
 ndim    = dimsizes(lat)

 topo    = nc_cosI->HSURF(0,nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1)
 lat2d   = lat(nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1)
 lon2d   = lon(nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1)

 hhl     = nc_cosI->HHL(0,::-1,nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1)
 ndimc   = dimsizes(hhl)

 zmodel     = 0.5*(hhl(0:ndimc(0)-2,:,:) + hhl(1:ndimc(0)-1,:,:))
 topo@lat2d = lat2d
 topo@lon2d = lon2d

 ndimc      = dimsizes(hhl)

 mpcorner      = (/lat2d(0,0),lon2d(0,0),lat2d(ndimc(1)-1,0),lon2d(0,ndimc(2)-1)/)

;-----------------------------------------------------------------------------------------
;Read storm mask
;-----------------------------------------------------------------------------------------
 fsname   = "smask_" + pname + "_" + sfold + ".nc"
 nc_spath = addfile(diri_spath + "/" + fsname,"r")
 smsk     = nc_spath->Band1
 st       = smsk@st
 et       = smsk@et
 ntime    = et - st + 1

 print("st : " + st + " et: " + et)

;-----------------------------------------------------------------------------------------
; READ NC FILES
;-----------------------------------------------------------------------------------------
 fils_cos  = systemfunc("ls " + diri_cosmo + "/cosout/lf*")
 nf_cos    = dimsizes(fils_cos)

 do nf   = st, et

 if (nf.eq.st) then
   ;Initialize variables
   iwv           = new((/nf_cos,ndimc(1),ndimc(2)/),"float")
   iwv@lat2d     = lat2d
   iwv@lon2d     = lon2d
   iwv@long_name = "Precipitable Water"   ;or water vapor
   iwv@units     = "mm"
   iwv@_FillValue = -999.
   iwv           = iwv@_FillValue

   lwp           = iwv
   lwp@long_name = "Liquid Water Path"

   iwp           = iwv
   iwp@long_name = "Ice Water Path" 
   lwp       = 0.
   iwp       = 0.
   iwv       = 0.
 end if

 filC    = addfile(fils_cos(nf),"r")
 qv      = filC->QV(0,::-1,nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1)
 t       = filC->T (0,::-1,nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1)
 p       = filC->P (0,::-1,nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1)
 qc      = filC->QC(0,::-1,nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1)
 qr      = filC->QR(0,::-1,nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1)
 qi      = filC->QI(0,::-1,nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1)
 qs      = filC->QS(0,::-1,nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1)
 qg      = filC->QG(0,::-1,nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1)
 qh      = filC->QH(0,::-1,nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1)
 u       = filC->U (0,::-1,nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1)
 v       = filC->V (0,::-1,nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1)

 qliq    = qc + qr
 qice    = qi + qs + qg + qh
 qrs     = qliq + qice
 qcd     = qc*0.  ; dummy, 
 rhoa    = calcrho_cos_nclf(t,p,qv,qcd,qrs)   ;air density kg/m3

 do iz = 0, ndimc(0)-2
   iwv(nf,:,:) = iwv(nf,:,:) + rhoa(iz,:,:)*qv  (iz,:,:)*(hhl(iz+1,:,:)-hhl(iz,:,:))    ;kg/m2 or mm  
   lwp(nf,:,:) = lwp(nf,:,:) + rhoa(iz,:,:)*qliq(iz,:,:)*(hhl(iz+1,:,:)-hhl(iz,:,:))    ;kg/m2 or mm
   iwp(nf,:,:) = iwp(nf,:,:) + rhoa(iz,:,:)*qice(iz,:,:)*(hhl(iz+1,:,:)-hhl(iz,:,:))    ;kg/m2 or mm
 end do

 print("Processing " + nf + "..." )

 if (fplot) then

 ;-----------------------------------------------------------------------------------------
 ; PLOT
 ;-----------------------------------------------------------------------------------------

 p    = p*0.01 ;hPa
 z_p  = wrf_user_intrp3d(zmodel,p,"h",pressure(0),0.,False)
 u_p  = wrf_user_intrp3d(u,p,"h",pressure(1),0.,False)
 v_p  = wrf_user_intrp3d(v,p,"h",pressure(1),0.,False)
 u_p@lat2d = lat2d
 u_p@lon2d = lon2d
 v_p@lat2d = lat2d
 v_p@lon2d = lon2d
 z_p@lat2d = lat2d
 z_p@lon2d = lon2d

 fncps_no = 21
 fncps_hgt = 0.028
 fname_out = "cps"
 wks_type                 = "png"
 wks = gsn_open_wks(wks_type,fname_out)
 ;gsn_define_colormap(wks,"precip2_17lev")

;Resources
 res                            = True                 ; plot mods desired
 ;gsn Resources
 res@gsnDraw                    = False
 res@gsnFrame                   = False
 restopo                        = res
 resgpm                         = res
 res@gsnScalarContour           = True
 res@gsnAddCyclic               = False
 res@gsnLeftStringFontHeightF   = fncps_hgt
 res@gsnRightStringFontHeightF  = fncps_hgt
 res@gsnLeftString              = "Precipitable Water"
 res@gsnRightString             = " mm"
 res@gsnMaximize                = True

 ;tm Resources
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
 res@tiXAxisFont                = fncps_no
 res@tiYAxisFont                = fncps_no
 res@tiXAxisFontHeightF         = fncps_hgt
 res@tiYAxisFontHeightF         = fncps_hgt
 res@tiXAxisString              = " Longitude"
 res@tiYAxisString              = " Latitute"

 ;lb Resources
 res@lbLabelAutoStride          = True
 res@lbOrientation              = "vertical"

 ;pm Resournces
 res@pmTickMarkDisplayMode      = "Always"   
 res@pmLabelBarWidthF           = 0.1
 res@pmLabelBarHeightF          = 0.5

 ; mp Resources Mapping resource (pollat=39.0,   pollon=-173.0,)s
 res@mpCenterLonF                 = 180.-173.0
 res@mpCenterLatF                 = 90. -39.0 ;
 res@mpLimitMode                  = "Corners"           ; method to zoom
 res@mpLeftCornerLatF             = mpcorner(0)
 res@mpLeftCornerLonF             = mpcorner(1)
 res@mpRightCornerLatF            = mpcorner(2)
 res@mpRightCornerLonF            = mpcorner(3)
 res@mpOutlineOn                  = True
 res@mpOutlineBoundarySets        = "National"
 res@mpDataSetName                = "Earth..4"    ; database for non-USA divisions
 res@mpDataBaseVersion            = "MediumRes"   ; Medium resolution database
 res@mpGridAndLimbOn              = True
 res@mpNationalLineThicknessF     = 3.0
 res@mpGeophysicalLineThicknessF  = 3.0
 res@mpProvincialLineThicknessF   = 2.0
 res@mpPerimOn                    = True
 res@mpFillOn                     = False

 ;cn Resources
 res@cnFillOn                   = True
 res@cnFillMode                 = "RasterFill"
 res@cnLinesOn                  = False
 res@cnInfoLabelOn              = False   ; Turn off annotations
 res@cnLineLabelsOn             = False
 res@cnLevelSelectionMode       = "ExplicitLevels"
 res@cnLevels                   = cnLevels

; vp Resournces
 res@vpXF                       = 0.150

; vc Resources
 res@vcRefMagnitudeF          = 20.0             ; define vector ref mag
 res@vcRefLengthF             = 0.045            ; define length of vec ref
 res@vcGlyphStyle             = "CurlyVector"    ; turn on curly vectors
 res@vcMinDistanceF           = 0.032            ; thin vectors
 res@vcRefAnnoOrthogonalPosF  = -.3             ; move ref vector down
 res@vcLineArrowColor         = "white" ;"grey40"
 res@vcLineArrowThicknessF    = 2.
 res@vcGlyphOpacityF          = 0.7
 res@vcRefAnnoOn              = True

 plot                          = gsn_csm_vector_scalar_map(wks,u_p,u_p,iwv, res)


 ;cn Resources
 restopo@cnFillOn               = True
 restopo@cnLinesOn              = True
 restopo@cnInfoLabelOn          = False
 restopo@cnLineLabelsOn         = False
 restopo@cnLineThicknessF       = 2.
 restopo@cnLineDashPattern      = 8
 restopo@cnLevelSelectionMode   = "ExplicitLevels"
 restopo@cnLevels               = (/100., 300., 500./)
 restopo@cnMonoFillColor        = True
 restopo@cnMonoFillPattern      = False
 restopo@cnMonoFillScale        = False
 restopo@cnFillPatterns         = (/17,4,3,4,3,4,3,4,3,4/)
 restopo@cnFillScales           = (/10.,.45,.40,.35,.25,.20,.18,.15,.12,.10/)
 restopo@cnLineColor            = "black"
 restopo@cnFillDrawOrder        = "PreDraw"

 ;lb Resources
 restopo@lbLabelBarOn           = False

 ;gsn Resources
 restopo@gsnLeftString          = ""
 restopo@gsnRightString         = ""
 dumtopo                        = gsn_csm_contour(wks,topo, restopo)
 overlay(plot,dumtopo)

 resgpm@cnFillOn               = False
 resgpm@cnLinesOn              = True
 resgpm@cnLineLabelsOn         = True
 resgpm@cnInfoLabelOn         = False
 resgpm@cnLineThicknessF       = 5.0
 resgpm@cnLineColor            = "black"

 resgpm@cnLevelSelectionMode   = "ExplicitLevels"
 resgpm@cnLevels               = ispan(5850,5920,10) ;

 dumgpm                         = gsn_csm_contour(wks,z_p, resgpm)
 overlay(plot,dumgpm)

 txres                       = True
 txres@gsMarkerSizeF = 12
 txres@gsMarkerThicknessF = 5.
 dumtxt = new(dimsizes(glat_id),"graphic")
 do ii= 0, dimsizes(glat_id)-1
   txres@gsMarkerIndex = gmarker(ii)
   dumtxt(ii) = gsn_add_polymarker(wks,plot,glon_id(ii),glat_id(ii),txres)
 end do

 draw(plot)
 frame(wks)
 delete([/res,restopo/])

 end if

 end do

 ;Crop smask
 smskc = smsk(:, nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1)
 lwpx   = where(smskc.eq.1, lwp,lwp@_FillValue)
 iwpx   = where(smskc.eq.1, iwp,iwp@_FillValue)

 fname = "LWP_" + pname + "_" + sfold + ".nc"
 system("rm " + fname)
 fout = addfile(fname,"c")
; fout->IWV = iwv
; fout->LWP = lwp
; fout->IWP = iwp
 fout->LWP = lwpx
 fout->IWP = iwpx
; fout->smsk = smskc

end
