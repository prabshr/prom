undef("readGRDC")
function readGRDC(fil_GRDC,year,opt)
begin
  nc_GRDC   = addfile(fil_GRDC,"r")
  qo_GRDC   = nc_GRDC->Band1
  yr_GRDC   = qo_GRDC(0,:)
  fstatus   = 0

  print("Extracting data for " + qo_GRDC@river + " " + qo_GRDC@lat + " " + qo_GRDC@lon)
  crpind = ind(yr_GRDC.eq.year)
  if (dimsizes(crpind).gt.1) then
    qoc_GRDC = qo_GRDC(:,crpind)
    x_GRDC   = new(dimsizes(crpind),"integer")
    do it = 0, dimsizes(crpind)-1
      yyyy     = tointeger(qoc_GRDC(0,it))
      mm       = tointeger(qoc_GRDC(1,it))
      dd       = tointeger(qoc_GRDC(2,it))
      x_GRDC(it) = day_of_year(yyyy,mm,dd)
    end do
  else
    print("Only few Data present for year " + year)
    fstatus = 1
  end if

  if (fstatus.ne.1) then
    if (opt) then
      ;ParFlow first file is skipped 
      mtime     = ispan(6,366,opt@freq)  
      ndim      = dimsizes(mtime)
      qoverland = new(ndim,"float")
      do it     = 0, ndim-1
        indm    = ind(x_GRDC.eq.mtime(it))
        if (.not.ismissing(indm)) then
          qoverland(it) = qoc_GRDC(3,indm)
        end if
      end do
      time      = mtime   ;-1
    else
      qoverland           = qoc_GRDC(3,:)
      time                = x_GRDC
    end if
      qoverland!0         = "time"
      qoverland&time      = time 
      qoverland@long_name = year + "Discharge "
      qoverland@status    = fstatus
      qoverland@river     = qo_GRDC@river
  else
   qoverland =  -999.0
   qoverland@status     = fstatus
  end if
  return(qoverland)
end

