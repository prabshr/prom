function cfads(height,var,minv,maxv,wbin,lname,unit,YH95)
begin
   ;Initialize
   nbins    = tointeger((maxv - minv)/wbin)
   print("Using " + nbins)
   var_cfad = new((/dimsizes(height),nbins/),"float")
   var_cfad!0 = "height"
   var_cfad!1 = "binval"

;;YH95
   nMax = new(dimsizes(height),"integer")
   do k = 0, dimsizes(height)-1
     opt         = True
     opt@bin_min = minv
     opt@bin_max = maxv
     if (dimsizes(dimsizes(var)).eq.5) then
       zpdf        = pdfx(var(:,:,k,:,:), nbins, opt)
     else if (dimsizes(dimsizes(var)).eq.4) then
       zpdf        = pdfx(var(:,k,:,:), nbins, opt)
     else
       zpdf        = pdfx(var(k,:,:), nbins, opt)
     end if
     end if

     if (all(ismissing(zpdf))) then
       nMax(k)        = nMax@_FillValue
     else
       nMax(k)          = zpdf@nMax
       bin_center       = zpdf@bin_center
     end if
     if (YH95) then
       var_cfad(k,:)    = tofloat(zpdf)/dzV(k)/binwidth(varid)    ;% YH95 per km-1 per dbZ
     else
       var_cfad(k,:)    = tofloat(zpdf)                            ;%
     end if
     print(k + " " + sum(zpdf))
     delete(zpdf)
   end do

  ;Mask CFAD if nMax(k) .lt. 20% of max(nMax)
  thrmax = 0.20
  do k = 0, dimsizes(height)-1
    if (ismissing(nMax(k)) .or. nMax(k).lt.thrmax*max(nMax)) then
      print("Masking at " + k + " " + nMax(k) + " < 20 % of max nMax")
      var_cfad(k,:) = var_cfad@_FillValue
    end if
  end do

  if (YH95) then
    print("YH95 units % perkm-1 per dbZ")
  else
    print("     units %")
  end if

  ;Assign co-ordinates
  var_cfad&height = (/height/)
  var_cfad&binval = tofloat(bin_center)
  var_cfad@long_name = lname
  var_cfad@units     = unit
  ;
  ;Write netcdf files for quick plotting
  foutname = "CFAD_" + lname + ".nc"
  system("rm " + foutname)
  fout = addfile(foutname,"c")
  fout->cfad = var_cfad

  return(var_cfad)
end

