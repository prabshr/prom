function radarline_nclf(Ipangle,lat_c,lon_c,glat_idB,glon_idB,extB,iopt)

;iopt = False
;Obtains the second co-ordintate and index passing throuhg a c-s at specified angle
;iopt = True
;Obtains two co-ordintates and indexes passing throuhg a c-s at specified angle
begin
  ;Line of Radar View and point on the Box
  print("radarline_nclf inputs "+ " " + Ipangle + " degree " + glat_idB + " N " + glon_idB + " E")

  pi     = 4.0*atan(1.0)
  if (Ipangle.le.90.) then
    Sangle = 90. - Ipangle
  else
    Sangle = 360. + 90. - Ipangle
  end if
  mslp   = tan(pi*(Sangle)/180.)
  cosf   = abs(cos(pi*glat_idB/180.))
  cinc   = glat_idB - mslp*glon_idB*cosf
  ;;print("Radar line " + cinc + " " + mslp)

  ndim   = dimsizes(lon_c)

  if (iopt) then
    ;#Search for left boundary
    bpx1   = min(lon_c(:,0))
    bpy1   = mslp * bpx1 * cosf + cinc
    if (bpy1 .lt. min(lat_c) .or. bpy1.gt.max(lat_c)) then
      if ((Ipangle.lt.90.) .or. ((Ipangle.gt.180.) .and. (Ipangle.lt.270.))) then
        bpy1 = min(lat_c)
      else
        bpy1 = max(lat_c)
      end if
      bpx1 = (bpy1 - cinc)/mslp/cosf
    end if
    ;#Search for  right boundary
    bpx2   = max(lon_c(:,ndim(1)-1))
    bpy2   = mslp * bpx2 * cosf + cinc
    if (bpy2 .lt. min(lat_c) .or. bpy2.gt. max(lat_c)) then
      if ((Ipangle.lt.90.) .or. ((Ipangle.gt.180.) .and. (Ipangle.lt.270.))) then
        bpy2 = max(lat_c)
      else
        bpy2 = min(lat_c)
      end if
      bpx2 = (bpy2 - cinc)/mslp/cosf
    end if
    nm1     = getind_latlon2d (lat_c,lon_c, bpy1, bpx1)
    nm2     = getind_latlon2d (lat_c,lon_c, bpy2, bpx2)
    plane   = (/nm1(0,1),nm2(0,1),nm1(0,0),nm2(0,0)/)*1.0 + 1.
   
    plane   = (/nm1(0,1),nm1(0,0),nm2(0,1),nm2(0,0)/)*1.0  ;NCL Version 6.6.2

    if (Ipangle.le.180) then
      glon_p  = glon_idB + extB* abs(cos(pi*(90.-Ipangle)/180.))
    else
      glon_p  = glon_idB - extB* abs(cos(pi*(90.-Ipangle)/180.))
    end if
    glat_p  = mslp*glon_p*cosf + cinc
  else
    ;Find a point inside box passing through the radar line
    ;distB  = sqrt((avg(lat_c)-glat_idB)^2. + (cos(pi*glat_idB/180.)*(avg(lon_c)-glon_idB))^2.)
    distB  = gc_latlon(glat_idB,glon_idB,avg(lat_c),avg(lon_c),2,2)
    if (Ipangle.le.180) then
      glon_p = glon_idB + distB*cos(pi*(90.-Ipangle)/180.)  ;guess of x on the box along radar line
      glat_p = mslp*glon_p*cosf + cinc
    else
      glon_p = glon_idB - distB*cos(pi*(90.-Ipangle)/180.)  ;guess of x on the box along radar line
      glat_p = mslp*glon_p*cosf + cinc
    end if

    nm      = getind_latlon2d (lat_c,lon_c, glat_p, glon_p)
    plane   = (/nm(0,1),nm(0,0)/)*1. + 1.      ;x and y co-ordinates
  end if

  plane@X = glon_p
  plane@Y = glat_p
  return(plane)
end
