undef("rootCLM")
function rootCLM(fsurfil,fclmfil,nbdl)
begin
 print("Using .." + fsurfil)

 ;CLM PFT root parameters
 roota  = (/-999.,7.,7.,7.,7.,7.,6.,6.,6.,7. ,7. ,7. ,11.,11.,11.,6.,7./)
 rootb  = (/-999.,2.,2.,2.,1.,1.,2.,2.,2.,1.5,1.5,1.5, 2., 2., 2.,3.,2./)

 ;First extract 2D root parameters 
 nc_surf   = addfile(fsurfil,"r")
 pft_3d    = tofloat(nc_surf->PCT_PFT)
 ndim      = dimsizes(pft_3d)

 ;Initialize bare soil
 roota_par    = pft_3d(0,nbdl:ndim(1)-nbdl-1,nbdl:ndim(2)-nbdl-1)*0./100.
 rootb_par    = pft_3d(0,nbdl:ndim(1)-nbdl-1,nbdl:ndim(2)-nbdl-1)*0./100.

 do ipft =  1, 16
   roota_par  = roota_par + pft_3d(ipft,nbdl:ndim(1)-nbdl-1,nbdl:ndim(2)-nbdl-1)*roota(ipft)/100.
   rootb_par  = rootb_par + pft_3d(ipft,nbdl:ndim(1)-nbdl-1,nbdl:ndim(2)-nbdl-1)*rootb(ipft)/100.
 end do

 ; Compute root fraction based on CLM PFT parameters 

 nc_clm     = addfile(fclmfil,"r")
 zsoi       = nc_clm->ZSOI(:,nbdl:ndim(1)-nbdl-1,nbdl:ndim(2)-nbdl-1)
 dzsoi      = nc_clm->DZSOI(:,nbdl:ndim(1)-nbdl-1,nbdl:ndim(2)-nbdl-1)

 nsoi       = dimsizes(zsoi)
 nlevsoi    = nsoi(0)
 ny         = nsoi(1)
 nx         = nsoi(2)

 zi         = new((/nlevsoi+1,ny,nx/),"float")
 rootfr     = new((/nlevsoi,ny,nx/),"float")
 zi(0,:,:)  = 0.

 do lev = 1, nlevsoi-1
   zi(lev,:,:)   = zi(lev-1,:,:) + dzsoi(lev-1,:,:)
   rootfr(lev-1,:,:) = .5*( exp(-roota_par * zi(lev-1,:,:))  \
                       + exp(-rootb_par * zi(lev-1,:,:))  \
                       - exp(-roota_par * zi(lev,:,:  ))  \
                       - exp(-rootb_par * zi(lev,:,:  )) )
 end do
;
  zi(nlevsoi,:,:)     = zi(nlevsoi-1,:,:) + dzsoi(nlevsoi-1,:,:)
  rootfr(nlevsoi-1,:,:) = .5*( exp(-roota_par * zi(nlevsoi-1,:,:))  \
                     + exp(-rootb_par * zi(nlevsoi-1,:,:)) )

; For bare-soil, set first layer root fraction as 1.0
  rootfr(0,:,:)         = where(roota_par.eq.0,1.0,rootfr(0,:,:))
  rootfr(nlevsoi-1,:,:) = where(roota_par.eq.0,0.0,rootfr(nlevsoi-1,:,:))

  return(rootfr)

end
