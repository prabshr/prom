;Plot the plan veiw of polarimetric variables at model level
;Examine enhanced column of ZDR and KDP above melting layer 

;P. Shrestha
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$prom_root/functions/read_emvorado_pol.ncl"
load "$PROM/functions/load_radarColors.ncl"
external EXT1 "$prom_root/functions/readdim.so"
external EXT2 "$prom_root/functions/readascii.so"

begin
;User Settings ---------------------------------------------------------------
 pid         = 0
 lnd         = "HET"
 itype_gscp  = "2483"
 inst        = 7 
 iz          = 47 ;21(low level) ; 47(mid level) , 41? ;model level to plot
;;
 casenames   = (/"bonnRadar_20150705"/)
 pnames      = (/"PROM1","PROM2","PROM3"/)
 pvars       = (/"Zh","Zdr","KDP","RHV"/)
 punits      = (/"dBZ","dB","deg/km ", "-"/)
;;;
 casename    = casenames(pid)
 pname       = pnames(pid)
 expname     = lnd + itype_gscp
 sfold       = "runE_" + inst + "_" + expname 
 diri_cosmo  = "/daten04/pshrestha/shared-data/" + casename + "/" + expname + "/" + sfold
 diri_emvo   = "/daten/daten04/pshrestha/shared-data/test-emvo/PROM1/" + sfold 
 diri_spath  = "/daten/daten04/pshrestha/shared-data/storm_mask/nc"

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
 do nf    = st,et

;Read EMVORADO outputs ------------------------------------------------------------
 ;varh         = new((/npol,nz,ny,nx/),"float")    ;Initialize variable for plotting
 varh         = new((/2,nz,ny,nx/),"float")
 varh@_FillValue = -999.99
 opt      = True
 opt@profile = -1   ;return 3d data
 opt@nf      =  nf 

 opt@var  = "zh"
 Zh = read_emvorado_pol(diri_emvo,opt)
 opt@var  = "zv"
 Zv = read_emvorado_pol(diri_emvo,opt)
 opt@var  = "kdp"
 KDP = read_emvorado_pol(diri_emvo,opt)
; opt@var  = "rrhv"
; rRHV = read_emvorado_pol(diri_emvo,opt)
; opt@var  = "irhv"
; iRHV = read_emvorado_pol(diri_emvo,opt)

; RHV = ((rRHV^2+iRHV^2)^0.5) / ((Zh*Zv)^0.5)   ;[-]
 KDP = KDP*1000.                               ;[deg/km]
 Zh  = 10*log10(Zh/1.)                         ;[dBZ]
 Zv  = 10*log10(Zv/1.)                         ;[dBZ]

; varh(0,:,:,:) = where(Zh.lt.-5., varh@_FillValue, Zh)
 varh(0,:,:,:) = where(Zh.lt.-5., varh@_FillValue, Zh-Zv)
 varh(1,:,:,:) = where(Zh.lt.-5., varh@_FillValue, KDP)
; varh(3,:,:,:) = where(Zh.lt.-5., varh@_FillValue, RHV)
; delete([/Zh, Zv, KDP, RHV/])
 delete([/Zh, Zv/])

;Mask the output --------------------------------------------------------

 varh = mask(varh,smsk(opt@nf,:,:),1)

;Plotting --------- ------------------------------------------------------
 ;Variable to plot
 spol1          = (/varh(0,iz,:,:)/)         ;
 spol1!0        = "Y"
 spol1!1        = "X"
 spol1&Y        = ispan(0,ny-1,1)  ;~ km
 spol1&X        = ispan(0,nx-1,1)  ;~ km

 spol2          = (/varh(1,iz,:,:)/)         ;
 spol2!0        = "Y"
 spol2!1        = "X"
 spol2&Y        = ispan(0,ny-1,1)  ;~ km
 spol2&X        = ispan(0,nx-1,1)  ;~ km

 plot     = new(2,"graphic")
 res      = True
 wks      = gsn_open_wks("png", "Figure_" + pname + "_" + sfold + "_" + opt@nf)
 colors   = load_radarColors(0)
 gsn_define_colormap(wks,colors)

 res@gsnDraw    = False
 res@gsnFrame   = False
 res@cnFillOn   = True
 res@cnFillMode = "CellFill"
 res@cnLinesOn  = False
 res@cnLevelSelectionMode       = "ExplicitLevels"

 ;Crop extents for plot based on mask
 smat  = where(smsk(opt@nf,:,:).eq.1,1,smsk@_FillValue)
 smsky = smat + ymat
 smskx = smat + xmat 

 res@trXMinF = max((/min(smskx) - 2,0/))
 res@trXMaxF = min((/max(smskx) + 2,nx-1/))
 res@trYMinF = max((/min(smsky) - 2,0/))
 res@trYMaxF = min((/max(smsky) + 2,ny-1/))

 res@tiMainString = "Model level: " + iz  ;;;;+ "| Avg. hgt.: " + avghgt + " m"

 res@cnLevels   = colors@zdr_levs
 res@gsnLeftString = pvars(1)
 res@gsnRightString = punits(1)
 plot(0)  = gsn_csm_contour(wks,spol1,res)

 res@gsnLeftString = pvars(2)
 res@gsnRightString = punits(2)
 res@cnLevels   = colors@kdp_levs
 plot(1)  = gsn_csm_contour(wks,spol2,res)


 delete(res@cnLevelSelectionMode)
 delete(res@cnLevels)
 res@gsnLeftString = ""
 res@gsnRightString = ""
 res@cnFillOn = False
 res@cnLinesOn = True
 res@cnLineLabelsOn = False
 res@cnLineColor    = "Gray"
 dum0      = gsn_csm_contour(wks,mhgt(iz,:,:),res)
 dum1      = gsn_csm_contour(wks,mhgt(iz,:,:),res)
 overlay(plot(0),dum0)
 overlay(plot(1),dum1)
 
 resP = True
 gsn_panel(wks,plot,(/1,2/),resP)

 end do    ;BIG LOOP

 foutdir = pname + "_" + expname + "_inst" + inst
 system("mkdir " + foutdir)
 system("mv *.png " + foutdir)

end
