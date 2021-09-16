;To get comparable storm statistics, need to have the mask for the storm
;as it propagates
;For deep convective storms, ZH is integrated from 9 to 12 km height
;else, this is done by generating a ZH map by integrating ZH from ML-1km to ML+2 km
;Iterative smoothing is also applied to the above integrated values

;P. Shrestha
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$PROM/functions/shade_terrain.ncl"
load "$PROM/functions/load_radarColors.ncl"
load "$prom_root/functions/read_emvorado_pol.ncl"

begin
 ;User Settings 
 ;---------------------------------------------------------------------------------
 pid       = 0            ;Experiment ID
 expname   = "HET2483"    ;Perturbation ID
 Zh_thr    = -5.0   ;dBz
 fgauss    = True         ;Default
 fdeep     = True         ;Deep Convective Storm or Shallow
 fsmooth   = False        ;Not applied

 ;inst      = 5            ;Ensemble number
 ;
 cnames    = (/"bonnRadar_20150705","bonnRadar_20160513","bonnRadar_20170706"/)
 pnames    = (/"PROM1","PROM2","PROM3"/)

 ;Pre-process
 ;---------------------------------------------------------------------------------
 casename    = cnames(pid)
 pname       = pnames(pid)
 sfold       = "runE_" + inst + "_" + expname
 diri_emvo   = "/daten/daten04/pshrestha/shared-data/test-emvo/" + pname + "/" + sfold
 diri_cosmo  = "/daten04/pshrestha/shared-data/" + casename + "/" + expname + "/" + sfold

;First extract geographical data
 ;---------------------------------------------------------------------------------
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
 mhgt  = (hhl(0,0:nz-1,ny/2,nx/2) + hhl(0,1:nz,ny/2,nx/2) ) *0.5 
 mhgt  = mhgt(::-1)/1000.

 ;Check number of COSMO outputs and extract height of ML
 fils_cos = systemfunc("ls " + diri_cosmo + "/cosout/*.nc")
 nfils = dimsizes(fils_cos)

 nc_cos  = addfile(fils_cos(0),"r")
 tair    = nc_cos->T(0,:,ny/2,nx/2)
 tair    = abs(tair(::-1) - 273.14)
 mind    = ind(tair.eq.min(tair))
 hmkm    = mhgt(mind)

 ;Now find indiecs of upper and lower limits from ML 
 if (fdeep) then
   hu5    = closest_val(10.0,mhgt)
   hl5    = closest_val(9.0,mhgt)
 else
   hu5    = closest_val(hmkm+2.0,mhgt)
   hl5    = closest_val(hmkm-1.0,mhgt)
 end if

 print("Melting layer height " + hmkm + " -  model layer " + mind)
 print("Integrating from model layers " + hl5 + " to " + hu5)
;; Read Zh data
 ;---------------------------------------------------------------------------------

 Zht   = new((/nfils,ny,nx/),"float")
 Zht@_FillValue = -999.99
 Zht@lat2d = lat2d
 Zht@lon2d = lon2d
 Zht@long_name = "Zh vertical sum"
 do nf = 0, nfils-1   ;BIG LOOP

   opt         = True
   ;opt@profile = 1   ;return verticl profles at ix, iy
   ;opt@iy      = iy 
   ;opt@ix      = ix
   opt@profile = 0   ;return 3D data
   opt@nf      = nf 

   opt@var  = "zh"
   Zh = read_emvorado_pol(diri_emvo,opt) 
   Zh  = 10*log10(Zh/1.)                         ;[dBZ]

   Zhsum = where(Zh.lt.Zh_thr, Zht@_FillValue, Zh)
   Zhsum@_FillValue = Zht@_FillValue
   Zhint  = dim_sum_n(Zhsum(hl5:hu5,:,:),0)        ;Integrate over ML-1 to ML+2 km height
   Zht(nf,:,:)  = Zhint

   if (fsmooth) then
    Zhsmth = Zht(nf,:,:)
    ;Smooth data (5-point weighte average)
    opt = True
    opt@cyclic = False
    Zhs_thr = 10.0 ;dBZ (integrated) threhsold
    do iter = 0,5
      wgt     = where(Zhsmth.lt.Zhs_thr,0.01,0.5)       ;Zhs_thr = 10.0 ;dBZ (integrated)
      Zhsmth  = wgt_area_smooth(Zhsmth,wgt,opt)
    end do
    Zht(nf,:,:)  = Zhsmth
    delete([/wgt,opt,Zhsmth/])
   end if
 
   if (fgauss) then
    Zhsmth = Zht(nf,:,:)
    ;Low pass filter using Gaussian Kernel
    nwt   = 3  ;nwt = 2*sigma + 1 
    sig   = 0.4 
    niter = 1 
    ;
    wgt   = filwgts_normal (nwt, sig , 0)
    do iter = 0, niter-1  
    do iy = 0, ny-1
      Zhsmth(iy,:) = wgt_runave(Zhsmth(iy,:), wgt, 0)
    end do
    do ix = 0, nx-1
      Zhsmth(:,ix) = wgt_runave(Zhsmth(:,ix), wgt, 0)
    end do
    end do
    Zht(nf,:,:)  = Zhsmth
    delete([/wgt,Zhsmth/])
   end if

   delete([/Zhsum, Zh, Zhint/])
 end do

;Write Output
 fname = pname + "_" + sfold + ".nc"
 system("rm " + fname)
 fout = addfile(fname,"c")
 fout->Band1 = Zht
end
