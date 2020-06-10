undef("storm_tracker_nclf")
;Identifies the location of maximum vertical wind speed

function storm_tracker_nclf(uz_hhl,lat,lon,wfac,wss,wse,ext,mpcorner)
begin
 uz_fill = dim_sum_n(uz_hhl,0)    ; sum of vertical velocity

 ;Mask regions outside to zero
 nm = getind_latlon2d (lat,lon, mpcorner(0), mpcorner(1))
 uz_fill(0:nm(0,0),:) = 0.
 uz_fill(:,0:nm(0,1)) = 0.

 nm = getind_latlon2d (lat,lon, mpcorner(2), mpcorner(3))
 uz_fill(nm(0,0):,:) = 0.
 uz_fill(:,nm(0,1):) = 0.

 ;Compute maximum vertical velocity
 wmax = max(uz_fill)

 uz_fill_1d    = ndtooned(uz_fill)
 dsizes        = dimsizes(uz_fill)
 indices       = ind_resolve(ind(uz_fill_1d.eq.wmax),dsizes)

 glatlon       = new(2,"float")
 glatlon(0)    = lat(indices(0,0),indices(0,1))
 glatlon(1)    = lon(indices(0,0),indices(0,1))

 return(glatlon)
end
