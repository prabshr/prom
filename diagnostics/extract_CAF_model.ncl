;Plot the time-series of convective area fraction for all ensemble members
;at specified  model level
; CAF is the fractional area of storm where the ZH > 40 dBZ
;V2.0
;Uses storm tracking and masking including EMVORADO outputs
;Oct 25 2021
;P. Shrestha
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$prom_root/functions/read_emvorado_pol.ncl"
load "$PROM/functions/load_radarColors.ncl"

begin
;User Settings ---------------------------------------------------------------
 pid         = 3 
 lnd         = "HET"
 itype       = "ART_SP4"
 ;inst        = 6 
 iz          = 34 ;21(low level) ; 47(mid level) , 41? ;model level to plot
 dbz_th      = 40.  ;dBZ
 dbz_min     =  0.  ; -5.  ;dBZ

;;
 casenames   = (/"bonnRadar_20150705"/)
 pnames      = (/"PROM1","PROM2","PROM3","EMPA2_PROM1"/)
 pvars       = (/"Zh","Zdr","KDP","RHV"/)
 punits      = (/"dBZ","dB","deg/km ", "-"/)
;;;
 ;casename    = casenames(pid)
 pname       = pnames(pid)
 itype_gscp  = lnd + itype
 sfold       = "runE_" + inst + "_" + itype_gscp
 expname     = lnd + "_" + pname
 subevent    = expname + "/" + itype_gscp + "/" + sfold

 ;CLUMA2
 ;diri_cosmo  = "/daten04/pshrestha/shared-data/" + casename + "/" + expname + "/" + sfold
 diri_emvo   = "/daten/daten04/pshrestha/shared-data/test-emvo/" + pname + "/" + lnd + "/" + sfold
 diri_spath  = "/daten/daten04/pshrestha/shared-data/storm_mask/nc"

 ;JUWELS
 diri_emvo   = "/p/largedata/hbn33/test-emvo/" + subevent
 diri_spath  = "/p/largedata/hbn33/storm_mask/" + expname + "/" + itype_gscp + "/nc/"
 diri_cosmo  = "/p/largedata/hbn33/shrestha1/" + subevent
 
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
 avghgt = avg((/mhgt(iz,:,:) - hsurf(0,:,:)/))  ;
 print("Avg. height = " + avghgt + " a.g.l.") 

;Read storm mask ------------------------------------------------------------------
 fsname   = "smask_" + pname + "_" + sfold + ".nc"
 nc_spath = addfile(diri_spath + "/" + fsname,"r")
 smsk     = nc_spath->Band1
 st       = smsk@st 
 et       = smsk@et
 ntime    = et - st + 1

;Initialize variables
 caf     = new(37,"float")
 tvar    = new(37,"float")

;BIG LOOP
 do nf    = st,et

;Read EMVORADO outputs ------------------------------------------------------------
 varh         = new((/nz,ny,nx/),"float")
 varh@_FillValue = -999.99
 opt      = True
 opt@profile = -1   ;return 3d data
 opt@nf      =  nf

 opt@var  = "zh"
 Zh = read_emvorado_pol(diri_emvo,opt)
 Zh  = 10*log10(Zh/1.)                         ;[dBZ]

 varh = where(Zh.lt.-5., varh@_FillValue, Zh)
 delete([/Zh/])

;Mask the output --------------------------------------------------------
 varh = mask(varh,smsk(opt@nf,:,:),1)

;Estimate CAF
 val1d = ndtooned(varh(iz,:,:))
 carea = tofloat(num(val1d.ge.dbz_th))
 tarea = tofloat(num(val1d.ge.dbz_min))
 if (tarea.gt.0) then    ;Avoid division by zero
    caf(nf) = carea/tarea
 end if
 tvar(nf) = tarea
 delete([/val1d/]) 
 print("Convective area fraction " + carea + " " + tarea)
 end do


 ;CAF output
 caf@model_level = iz
 caf@model_height = avghgt

 foname = "CAF_" + pname + "_" + sfold + ".nc"
 system("rm " + foname)
 fout = addfile(foname,"c")
 fout->Band1 = caf 
 fout->Band2 = tvar
end
