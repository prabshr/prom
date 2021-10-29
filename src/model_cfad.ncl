;New version, uses storm mask and EMVORADO ouputs
;Log: Oct 22 2021

;All variables preserve their native units
;scale/offset is only used for plotting
; Cumulative Frequency Altitude Diagram
; Based on Yuter and Houze 1995 (YH95)
; opt = True, generates CFADS for individual ensemble members
; opt = False,reads nc CFAD data for mulitple ensemble memebers to generate CFAD
; fcrop = True works with opt=False

;P. Shrestha
;Nov 18 2019
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$PROM/functions/cfads.ncl"
load "$PROM/functions/load_cfadColors.ncl"
load "$prom_root/functions/read_emvorado_pol.ncl"

function read_cfad(cdiri,cens)
begin
 nc = dimsizes(cens)
 do nf = 0, nc-1
   fname  = cdiri + "/cfad_" + cens(nf) + ".nc"
   nc_fil = addfile(fname + ".nc","r")
   var    = nc_fil->Band1
   if (nf.eq.0) then
     cdim = dimsizes(var)
     vare = new((/nc,cdim(0),cdim(1),cdim(2),cdim(3),cdim(4)/),typeof(var))
   end if
   vare(nf,:,:,:,:,:) = var
 end do
 vare!0          = "ens"
 vare!1          = "time"
 vare!2          = "pv"
 vare!3          = "Z"
 vare!4          = "lat"
 vare!5          = "lon"
 vare&Z          = var&Z
 vare@_FillValue = -99
 vare@lat2d      = var@lat2d
 vare@lon2d      = var@lon2d

 return(vare)
end

begin

 ;User Settings ---------------------------------------------------------------
 pid         = 0
 lnd         = "HET"
 ;itype  = "2463"
 inst        = 6 
 opt         = True
 wks_type     = "x11"
;;
 casenames   = (/"bonnRadar_20150705"/)
 pnames      = (/"PROM1","PROM2","PROM3"/)
 pvars     = (/"Z~B~H~N~","Z~B~DR~N~","K~B~DP~N~","~F8~r~B~~F10~hv~E~"/)
 punits    = (/"[dBZ]","[dB]","[deg/km]", "[-]"/)
;;;
 npol        = dimsizes(pvars)
 casename    = casenames(pid)
 pname       = pnames(pid)
 expname     = lnd + itype
 sfold       = "runE_" + inst + "_" + expname
 diri_cosmo  = "/daten04/pshrestha/shared-data/" + casename + "/" + expname + "/" + sfold
 diri_emvo   = "/daten/daten04/pshrestha/shared-data/test-emvo/PROM1/" + sfold
 diri_spath  = "/daten/daten04/pshrestha/shared-data/storm_mask/nc"

;Read netcdf files
;---------------------------------------------------------------------------
 if (opt) then
   ;Read Grid data ---------------------------------------------------------------
   fil_ivr = diri_cosmo + "/cosout/ivr/lfff00000000c.nc"
   ivr_nc  = addfile(fil_ivr,"r")
   lon2d   = ivr_nc->lon
   lat2d   = ivr_nc->lat
   hhl     = ivr_nc->HHL
   hsurf   = ivr_nc->HSURF

   ndim    = dimsizes(hhl)
   nz      = ndim(1)-1
   ny      = ndim(2)
   nx      = ndim(3)

   ;Extract model level height from half-level heights
   mhgt  = (hhl(0,0:nz-1,:,:) + hhl(0,1:nz,:,:) ) *0.5
   mhgt  = mhgt(::-1,:,:)


   ;Read storm mask ------------------------------------------------------------------
   fsname   = "smask_" + pname + "_" + sfold + ".nc"
   nc_spath = addfile(diri_spath + "/" + fsname,"r")
   smsk     = nc_spath->Band1
   st       = smsk@st
   et       = smsk@et

   ;Initialize for cropping extent while plotting
   x = ispan(0,nx-1,1)
   y = ispan(0,ny-1,1)
   xmat  = conform_dims((/ny,nx/),x,1)
   ymat  = conform_dims((/ny,nx/),y,0)

  ;BIG LOOP
  st = 12
  et = 24
  print("Manually changing st and et " + st + " " + et)
  print(" ********* ")

  ntime        = et - st + 1
  varh         = new((/ntime,npol,nz,ny,nx/),"float") 
  varh@_FillValue = -999.99
  varh!2       = "Z"
  varh&Z       = (mhgt(:,ny/2,nx/2) - hsurf(0,ny/2,nx/2))/1000.  ;height in km a.g.l.

  do nf    = st, et

    ;Read EMVORADO outputs ------------------------------------------------------------
    opt      = True
    opt@profile = -1   ;return 3d data
    opt@nf      =  nf

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

    iens=0
    varh(nf-st,0,:,:,:) = where(Zh.lt.-5., varh@_FillValue, Zh)
    varh(nf-st,1,:,:,:) = where(Zh.lt.-5., varh@_FillValue, Zh-Zv)
    varh(nf-st,2,:,:,:) = where(Zh.lt.-5., varh@_FillValue, KDP)
    varh(nf-st,3,:,:,:) = where(Zh.lt.-5., varh@_FillValue, RHV)
    delete([/Zh, Zv, KDP, RHV/])

    ;Mask the output --------------------------------------------------------
    varh(nf-st,:,:,:,:)  = mask(varh(nf-st,:,:,:,:),smsk(nf,:,:),1)

  end do
  
  ;foutname = "cfad_" + pname + "_" + expname + "_" + inst + ".nc"
  ;system("rm " + foutname)
  ;fout = addfile(foutname,"c")
  ;fout->Band1 = varh

  var = varh
  delete(varh)

 else

   varh =  read_cfad(diri,nens)

 end if

; CFAD plot parameters
;--------------------------------------------------------
 copt          = True
 copt@YH95     = False  ;scale CFAD and remove data with low % of N
 colors       = load_cfadColors(copt)
 minv         = colors@minv
 maxv         = colors@maxv
 binw         = colors@binw

 ndimv        = dimsizes(var)
 ndimv        = ndimv(::-1)
 nz           = ndimv(2) 
 ny           = ndimv(1) 
 nx           = ndimv(0) 
 print(nx + " " + ny + "" + nz)

 vid = 0
 copt@vid      = vid
 if (opt) then
   cvar = var(:,vid,:,:,:)
 else
   cvar = var(:,:,vid,:,:,:)
 end if
 zhh_cfad = cfads(var&Z,cvar,minv(vid),maxv(vid),binw(vid),pvars(vid),punits(vid),copt@YH95)

 vid      = 1
 copt@vid  = vid
 if (opt) then
   cvar = var(:,vid,:,:,:)
 else
   cvar = var(:,:,vid,:,:,:)
 end if
 zdr_cfad = cfads(var&Z,cvar,minv(vid),maxv(vid),binw(vid),pvars(vid),punits(vid),copt@YH95)

 vid      = 2
 copt@vid  = vid 
 if (opt) then
   cvar = var(:,vid,:,:,:)
 else
   cvar = var(:,:,vid,:,:,:)
 end if
 kdp_cfad = cfads(var&Z,cvar,minv(vid),maxv(vid),binw(vid),pvars(vid),punits(vid),copt@YH95)

 vid      = 3
 copt@vid  = vid
 if (opt) then
   cvar = var(:,vid,:,:,:)
 else
   cvar = var(:,:,vid,:,:,:)
 end if
 rho_cfad = cfads(var&Z,cvar,minv(vid),maxv(vid),binw(vid),pvars(vid),punits(vid),copt@YH95)

 ;Debug outputpfile
 foutname = "CFAD_ZH_" + pname + "_" + sfold + ".nc"
 system("rm " + foutname)
 fout = addfile(foutname,"c")
 fout->Band1 = zhh_cfad
 foutname = "CFAD_ZDR_" + pname + "_" + sfold + ".nc"
 system("rm " + foutname)
 fout = addfile(foutname,"c")
 fout->Band1 = zdr_cfad
 foutname = "CFAD_KDP_" + pname + "_" + sfold + ".nc"
 system("rm " + foutname)
 fout = addfile(foutname,"c")
 fout->Band1 = kdp_cfad
 foutname = "CFAD_RHO_" + pname + "_" + sfold + ".nc"
 system("rm " + foutname)
 fout = addfile(foutname,"c")
 fout->Band1 = rho_cfad

;Plot
;---------------------------------------------------------------------------
 plot = new(4,"graphic")

 fname_out = "cfad_model_" + pname + "_" + sfold + "_" + inst

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
 res@gsnLeftString              = "Model " ;+ timestr
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

; Add annotation
 txres               = True                      ; text mods desired
 txres@txFontHeightF = 0.015                    ; text font height
 txres@txJust        = "CenterLeft"              ; Default is "CenterCenter".

 gsn_text_ndc(wks,"b)",.01,.67,txres)


 resP = True
 resP@gsnMaximize    = True
 resP@gsnPanelXF = (/-1,0.28,0.51,0.74/)
; resP@gsnPanelDebug = True
 gsn_panel(wks,plot,(/1,4/), resP)

end
