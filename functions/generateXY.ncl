function generateXY(lat,lon,glat,glon)
begin
 ;Generate X-Y coordinates in km
  dxy     = 1.132 ;km
  ndim    = dimsizes(lat)
  ny      = ndim(0)
  nx      = ndim(1)

  nm      = getind_latlon2d (lat,lon, glat, glon)
  yc      = nm(0,0)
  xc      = nm(0,1)
  XYkm    = new((/2,ny,nx/),"float")
  do ji = 0, ny-1
    XYkm(1,ji,:) = (ji - yc)*dxy
  end do
  do ii = 0, nx-1
    XYkm(0,:,ii) = (ii - xc)*dxy
  end do
  return(XYkm)
end
