undef("pdfRootMoisture")

; Volumetric Soil Moisture weighted by root distribution n
; OPTIONS

function pdfRootMoisture(fils_clm,nbdl,rootfr)
begin

;Read State Variables
;---------------------------------------------------------------
 ;CLM tapes and forcing files , leap year not supported

 ncdf_clm  = addfiles(fils_clm,"r")
 ListSetType (ncdf_clm, "cat")

 dzsoi_full = ncdf_clm[0]->DZSOI
 ndim  = dimsizes(dzsoi_full)
 nz    = ndim(0)

 ; Read CLM history tape
 Svol   = ncdf_clm[:]->H2OSOI(:,:,nbdl:ndim(1)-nbdl-1,nbdl:ndim(2)-nbdl-1)  ;Volumetric (m3/m3) 
 dzsoi  = ncdf_clm[0]->DZSOI(:,nbdl:ndim(1)-nbdl-1,nbdl:ndim(2)-nbdl-1)

; Computation
;---------------------------------------------------------------
 Sv_scale = Svol(:,0,:,:)*0.0

 do it = 0, dimsizes(fils_clm)-1
   wgt = Svol(0,0,:,:)*0.0 
   dpdzsum = wgt 
   do k = 0, nz-1
     dp_dz = rootfr(k,:,:)/dzsoi(k,:,:)
     wgt   = wgt + Svol(it,k,:,:)*dp_dz*dzsoi(k,:,:)
     dpdzsum = dpdzsum + dp_dz*dzsoi(k,:,:)
   end do
   Sv_scale(it,:,:) = wgt/dpdzsum      ;Root Distribution Scaled Volumetric Soil Moisture 
 end do

 varout = Sv_scale
 varout@_FillValue = -999.99
 return(varout)
end
