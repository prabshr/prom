function read_radarnc(fname,varid,clon,clat)
begin
  varname  = (/"ZH_grid","ZDR_grid","KDP_grid","RHO_grid"/)
  pvars    = (/"Z~B~H~N~","Z~B~DR~N~","K~B~DP~N~","~F8~r~B~~F10~hv~E~"/)
  punits   = (/"[dBZ]","[dB]","[deg/km]", "[-]"/)

  nc_fil = addfile(fname,"r")
  x       = nc_fil->x
  y       = nc_fil->y
  z       = nc_fil->z
  var     = nc_fil->$varname(varid)$
  ;
  tstr    = str_split(nc_fil@ScanStart,":")
  time    = tofloat(tstr(0)) + tofloat(tstr(1))/60. 

  ;Estimate lat-lon
  dxy = 111.0  ;km
  d2r = 4*atan(1)/180.
  ;lat  = lat(100) + (y(ix,:,0)*cos(glat_idB*d2r))/dxy

  ;Attributes
  var@_FillValue = 999
  var!0   = "Z"
  var!1   = "Y"
  var!2   = "X"
  var&Z   = z
  var&Y   = y
  var&X   = x
  var@long_name     = pvars(varid)
  var@units  =  punits(varid)
  var@time      = time
  print("Extracting " + varname(varid) + " " + nc_fil@ScanStart)
  return(var)
end

