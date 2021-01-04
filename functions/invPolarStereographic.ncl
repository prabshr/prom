function invPolarStereographic(x,y,nps)
;Based on https://docs.wradlib.org/en/stable/notebooks/radolan/radolan_grid.html
begin
  ndimx = dimsizes(x)
  ndimy = dimsizes(y)

  R2    = (nps@earth_radius_of_projection)^2
  lamda0= nps@straight_vertical_longitude_from_pole
  phi0  = nps@standard_parallel

  latlon   = new((/2,ndimy,ndimx/),"double")
  d2r      = 4*atan(1)/180.

  do jj= 0, ndimy-1
  do ii= 0, ndimx-1
    latlon(0,jj,ii) = atan(-x(ii)/y(jj))*1./d2r + lamda0
    numer = R2*(1 + sin(d2r*phi0))^2  - (x(ii)*x(ii) + y(jj)*y(jj))
    denom = R2*(1 + sin(d2r*phi0))^2  + (x(ii)*x(ii) + y(jj)*y(jj))
    latlon(1,jj,ii) = asin(numer/denom)*1./d2r
  end do
  end do
 
  return(latlon)

end
