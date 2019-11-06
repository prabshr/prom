undef("melt_layer_height_nclf")
function melt_layer_height_nclf(tair,zmodel)
 begin
  ;tair       ;vertical profile of air temperature in K
  ;zmodel     ;model main level height
  ;zmelt      ;melting layer height in m
  

  vart        = tair - 273.15   ;
  index       = closest_val(0.0,vart)
  if (vart(index).eq.0) then
    zmelt = zmodel(index)
  else
    nzm = 10
    tair1d = new(nzm,"float")
    zm1d   = new(nzm,"float")
    do iz = 0, nzm-1
      tair1d(iz) = vart(index+nzm/2-iz)
      zm1d(iz)   = zmodel(index+nzm/2-iz)
    end do
    tairxo       = tair1d
    tairxo(iz/2) = 0.
    zmxo = linint1 (tair1d,zm1d,False, tairxo, 0)
    zmelt = zmxo(iz/2)
   ;print(zm1d + " " + tair1d )
 end if
 print("Melting layer near :" + zmodel(index) + " Tair =" + vart(index))
 print("Interpolated height :" + zmelt)
 return(zmelt)
end

