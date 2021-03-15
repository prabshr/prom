function water_path_extract_nclf(var_in,var_out,lat,lon,opt)
begin
  ndim = dimsizes(var_in)
  glat_idB = opt@glat_idB
  glon_idB = opt@glon_idB
  ext      = opt@ext
  wfac     = opt@wfac
  wss      = opt@wss
  wse      = opt@wse
  offset   = opt@offset

  mpc      = 1.*(/glat_idB-ext(0),glon_idB-ext(1),glat_idB+ext(0),glon_idB+ext(1)/)
  varo = new((/ndim(0),ndim(1),ndim(2)/),"float")
  data1d = new(ndim(0),"float")
  do fn = 0, ndim(0)-1
    uz_fill = var_in(fn,:,:)

    ;Find location where W is maximum ....
    if (fn.eq.0) then
      mpc_find = mpc
    else
      mpc_find = mpcw
    end if

    ;Mask regions outside to zero
    nm = getind_latlon2d (lat,lon, mpc_find(0), mpc_find(1))
    uz_fill(0:nm(0,0),:) = 0.
    uz_fill(:,0:nm(0,1)) = 0.

    nm = getind_latlon2d (lat,lon, mpc_find(2), mpc_find(3))
    uz_fill(nm(0,0):,:) = 0.
    uz_fill(:,nm(0,1):) = 0.

   ;Compute maximum
    wmax = max(uz_fill)

    uz_fill_1d    = ndtooned(uz_fill)
    dsizes        = dimsizes(uz_fill)
    indices       = ind_resolve(ind(uz_fill_1d.eq.wmax),dsizes)

    glatlon      = new(2,"float")
    glat_id      = lat(indices(0,0),indices(0,1))
    glon_id      = lon(indices(0,0),indices(0,1))

    print("Updating mpcorner and mpcornerw")
    mpcw  = tointeger((/glat_id-wss*wfac*ext(0),glon_id-wss*wfac*ext(1),glat_id+wse*wfac*ext(0),glon_id+wse*wfac*ext(1)/))
    mpcw(0) = max((/mpcw(0)-offset,0/))
    mpcw(1) = max((/mpcw(1)-offset,0/))
    mpcw(2) = min((/mpcw(2)+offset,ndim(1)-1/))
    mpcw(3) = min((/mpcw(3)+offset,ndim(2)-1/))

    varo(fn,mpcw(0):mpcw(2),mpcw(1):mpcw(3)) = var_out(fn,mpcw(0):mpcw(2),mpcw(1):mpcw(3))

    data1d(fn) = sum(varo(fn,:,:))
  end do


 ;Debug Output
  system("rm " + opt@fname + ".nc")
  fout = addfile(opt@fname + ".nc","c")
  fout->Band1 = varo

  return(data1d)
end

