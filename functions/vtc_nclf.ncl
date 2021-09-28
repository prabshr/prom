undef("vtc_nclf")
function vtc_nclf(data3d,P,T,hhl,htrop,uscale)
;extract vertical tropospheric columns for gases (molecules per cm2)
;from model outputs
;hhl in m, data3d in ppmv (i.e. parts per million by volume) for gaseous tracers
;P in Pascals, T in K
;htrop in m
;uscale to scale the processed data

begin
  R      = 8.314 ; J/K/mol
  Na     = 6.022 * 1E23 ; molecules/mol

 ; ppmv = 1 [m3 gas]   * Na * dz [m]          Na * dz  [molecules]
 ;        -------------------------------  = -------------------
 ;         1E6 [m3 air] * V [m3 gas] *1E4      V * 1E10 [cm2 air] 

  scale  = 1E-10 ;

 ;
  ndim   = dimsizes(hhl)
  zind   = closest_val(htrop,hhl(:,ndim(1)/2,ndim(2)/2))
  ;print(zind)

  ;Compute volume per mole 
  V = R * T / P             ;m3

  ;Compute model dz
  dz = hhl(1:ndim(0)-1,:,:) - hhl(0:ndim(0)-2,:,:)

  ;Scale and aggregate over htrop
  data2d = new((/ndim(1),ndim(2)/),"float")
  data2d = 0.
  do iz = 0, zind-1
    data2d = data2d + data3d(iz,:,:) * dz(iz,:,:)  * scale * (Na/V(iz,:,:))
  end do

  data2d = data2d / uscale 
  data2d@long_name = data3d@long_name
  data2d@units     = uscale + "  molecules cm-2"
  return(data2d)
end
