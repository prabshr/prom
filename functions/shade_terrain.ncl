function shade_terrain(range,depth)
begin

 dimsr  = dimsizes(range)
 ndepth = dimsr(0)
 nrange = dimsr(1)
;
 min_range = min(range)
 max_range = max(range)
 min_depth = min(depth) 
 max_depth = max(depth)

;
; Create arrays to hold the points that will outline the
; terrain area to be shaded in gray.
;
  npoly = nrange + 3    ; The terrain line plus 3 points to 
                        ; close the polygon

  xpoly = new(npoly,typeof(range))
  ypoly = new(npoly,typeof(range))

  xpoly(0:nrange-1) = (/range(0,:)/)    ; The terrain
  ypoly(0:nrange-1) = (/depth(0,:)/)

  xpoly(nrange)     = max_range         ; Bottom right corner of plot
  ypoly(nrange)     = min_depth

  xpoly(nrange+1)   = min_range         ; Bottom left corner of plot
  ypoly(nrange+1)   = min_depth

  xpoly(nrange+2)   = range(0,0)        ; Start of terrain. This closes
  ypoly(nrange+2)   = depth(0,0)        ; the polygon.

  xypoly = new((/2,npoly/),typeof(range))
  xypoly(0,:) = xpoly
  xypoly(1,:) = ypoly

  return(xypoly)

end
