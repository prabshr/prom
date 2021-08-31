load "$prom_root/functions/calcrho_cos_nclf.ncl"
function extract_atmwater_nclf(fil_cos,fil_ivr,nbdl,opt)
;Computes Precipitatble Water, LWP and IWP

begin


 filC    = addfile(fil_cos,"r")
 filI    = addfile(fil_ivr,"r")

 lat     = filI->lat
 lon     = filI->lon

 ndim    = dimsizes(lat)
 lat2d   = lat(nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1)
 lon2d   = lon(nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1)
 hhl     = filI->HHL(0,::-1,nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1)
 ndimc   = dimsizes(hhl)

;Initialize variables
 iwv           = new((/ndimc(1),ndimc(2)/),"float")
 iwv           = 0.
 iwv@lat2d     = lat2d
 iwv@lon2d     = lon2d
 iwv@long_name = "Precipitable Water"   ;or water vapor
 iwv@units     = "mm"

 lwp           = iwv
 lwp@long_name = "Liquid Water Path"

 iwp           = iwv
 iwp@long_name = "Ice Water Path"
 
;Extract
 qv      = filC->QV(0,::-1,nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1)
 t       = filC->T (0,::-1,nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1)
 p       = filC->P (0,::-1,nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1)
 qc      = filC->QC(0,::-1,nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1)
 qr      = filC->QR(0,::-1,nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1)
 qi      = filC->QI(0,::-1,nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1)
 qs      = filC->QS(0,::-1,nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1)
 qg      = filC->QG(0,::-1,nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1)
 qh      = filC->QH(0,::-1,nbdl:ndim(0)-nbdl-1,nbdl:ndim(1)-nbdl-1)

 qliq    = qc + qr
 qice    = qi + qs + qg + qh
 qrs     = qliq + qice
 rhoa    = calcrho_cos_nclf(t,p,qv,qc,qrs)   ;air density kg/m3

 do iz = 0, ndimc(0)-2
    dz = (hhl(iz+1,:,:)-hhl(iz,:,:))
   iwv =  iwv + rhoa(iz,:,:)*qv(iz,:,:)*dz      ;kg/m2 or mm  
   lwp =  lwp + rhoa(iz,:,:)*qliq(iz,:,:)*dz    ;kg/m2 or mm
   iwp =  iwp + rhoa(iz,:,:)*qice(iz,:,:)*dz    ;kg/m2 or mm
 end do

 if (opt.eq.0) then
   return(iwv)
 else if (opt.eq.1) then
   return(lwp)
 else if (opt.eq.2) then
   return(iwp)
 end if
 end if
 end if

end
