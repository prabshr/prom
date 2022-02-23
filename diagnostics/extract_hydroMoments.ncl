load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$TPS/functions/calcrho_cos_nclf.ncl"

;Extracts mass density, number density, particle diameters

begin

;-----------------------------------------------------------------------------------------
; USER SETTINGS
;-----------------------------------------------------------------------------------------
 ic          = 0 ;select hydrometeor

 pid         = 3 
 ;inst        = 6
 lnd         = "HET"
 itype       = "ART" 
 casenames   = (/"bonnRadar_20150705"/)
 pnames      = (/"PROM1","PROM2","PROM3","EMPA2HC_PROM1"/)
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
 hydrometeors_qx = (/"QC","QR","QI","QS","QG","QH"/)
 hydrometeors_qn = (/"NCCLOUD","NCRAIN","NCICE","NCSNOW","NCGRAUPEL","NCHAIL"/)
 a_geo        = (/0.124    , 0.124  , 0.835 , 2.4    , 0.142     , 0.1366/)  ;mkg-b
 x3           = 1./3.
 b_geo        = (/x3       , x3     , 0.39  , 0.455  , 0.314     , x3/)      ;-
 xmin         = (/4.2E-15  , 2.6E-10, 1.E-12, 1.E-10, 1.E-9, 2.6E-9/)        ;kg
 xmax         = (/2.6E-10  , 3.0E-6 , 1.E-6 , 2.E-5 , 5.E-4, 5.E-4/)         ;g

 ;
 fils_cos  = systemfunc("ls " + diri_cosmo + "/cosout/lf*")
 nf_cos    = dimsizes(fils_cos)

 do nf   = st, et

 if (nf.eq.st) then
   ;Initialize variables
   hval           = new((/nf_cos,ndimc(1),ndimc(2)/),"float")
   hval@lat2d     = lat2d
   hval@lon2d     = lon2d
   hval@_FillValue = -999.
   hval = 0. 
   hval1 = hval
   hval2 = hval
 end if

 filC    = addfile(fils_cos(nf),"r")
 qx      = filC->$(hydrometeors_qx(ic))$(0,::-1,nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1)
 qn      = filC->$(hydrometeors_qn(ic))$(0,::-1,nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1)

 ;density
 t       = filC->T (0,::-1,nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1)
 p       = filC->P (0,::-1,nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1)
 qv      = filC->QV(0,::-1,nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1)
 qc      = filC->QC(0,::-1,nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1)
 qr      = filC->QR(0,::-1,nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1)
 qi      = filC->QI(0,::-1,nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1)
 qs      = filC->QS(0,::-1,nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1)
 qg      = filC->QG(0,::-1,nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1)
 qh      = filC->QH(0,::-1,nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1)

 qliq    = qc + qr
 qice    = qi + qs + qg + qh
 qrs     = qliq + qice
 qcd     = qc*0.  ; dummy,
 rhoa    = calcrho_cos_nclf(t,p,qv,qcd,qrs)   ;air density kg/m3

 ;********************
 qx      = qx*rhoa    ;kg/m3
 qn      = qn*rhoa    ;m-3

 ;********************
 qn      = where(qn.eq.0, hval@_FillValue, qn)
 qn@_FillValue = hval@_FillValue

 xmass   = (/(qx/qn)/)   ; kg
 ;Clipping
 xmass   = where(xmass.gt.xmin(ic),xmass,0.)
 xmass   = where(xmass.gt.xmax(ic),xmax(ic),xmass)
 ;diameter mass relationship
 dx      = (/(a_geo(ic) * (xmass)^b_geo(ic))/)   
 dx@_FillValue = hval@_FillValue
 dx      = where(dx.eq.0, hval@_FillValue, dx)

 qx@_FillValue = hval@_FillValue
 qx      = where(qx.eq.0, hval@_FillValue, qx)
 

 ;********************

 hval(nf,:,:) = dim_avg_n(dx,0)*1.E6         ;um
 hval1(nf,:,:) = dim_avg_n(qx,0)*1.E3   ;g/m3
 hval2(nf,:,:) = dim_avg_n(qn,0)/1.E6   ;cm-3

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
 hvalx   = where(smskc.eq.1, hval,hval@_FillValue)
 hval1x  = where(smskc.eq.1, hval1,hval@_FillValue)
 hval2x  = where(smskc.eq.1, hval2,hval@_FillValue)


 hvalx@long_name = "Average Diameter of " + hydrometeors_qx(ic)
 hvalx@units     = "um"

 hval1x@long_name = "Mass density of " + hydrometeors_qx(ic)
 hval1x@units   =  "g/m3"

 hval2x@long_name = "Number density of " + hydrometeors_qx(ic)
 hval2x@units   =  "cm-3"

 fname = "QX_" + pname + "_" + sfold + ".nc"
 system("rm " + fname)
 fout = addfile(fname,"c")
 fout->DX = hvalx
 fout->QX = hval1x
 fout->QN = hval2x
end
