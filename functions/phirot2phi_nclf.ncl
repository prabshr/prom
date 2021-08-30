undef("phirot2phi_nclf")
;convertis phi from on rotated system to another
;------------------------------------------------------
function phirot2phi_nclf(phirot, rlarot, polphi, pollam, polgam, opt)
begin
  zrpi18 = 57.2957795d                      
  zpir18 = 0.0174532925d

  zsinpol    = sin (zpir18 * polphi)
  zcospol    = cos (zpir18 * polphi)

  if (opt.eq.0) then ; phi

    zphis      = zpir18 * phirot
    if (rlarot .gt. 180.0d) then
      zrlas = rlarot - 360.0d
    else
      zrlas = rlarot
    end if
    zrlas       = zpir18 * zrlas

    if (polgam .ne. 0.0d) then
      zgam  = zpir18 * polgam
      zarg  = zsinpol*sin (zphis) +    \ 
          zcospol*cos(zphis) * ( cos(zrlas)*cos(zgam) - sin(zgam)*sin(zrlas) )
    else 
      zarg  = zcospol * cos (zphis) * cos (zrlas) + zsinpol * sin (zphis)
    end if

    rot2geo  = zrpi18 * asin (zarg)

  else if (opt.eq.1) then   ;lambda

    zlampol = zpir18 * pollam
    zphis   = zpir18 * phirot
    if (rlarot .gt. 180.0d) then
      zrlas = rlarot - 360.0d
    else
      zrlas = rlarot
    end if
    zrlas   = zpir18 * zrlas

    if (polgam .ne. 0.0d) then 
      zgam    = zpir18 * polgam
      zarg1   = sin (zlampol) *                                               \ 
        (- zsinpol*cos(zphis) * (cos(zrlas)*cos(zgam) - sin(zrlas)*sin(zgam))  \
         + zcospol * sin(zphis))                                               \
      - cos (zlampol)*cos(zphis) * (sin(zrlas)*cos(zgam) + cos(zrlas)*sin(zgam))

      zarg2   = cos (zlampol) *                                                \
        (- zsinpol*cos(zphis) * (cos(zrlas)*cos(zgam) - sin(zrlas)*sin(zgam))  \
         + zcospol * sin(zphis))                                               \ 
         +  sin (zlampol)*cos(zphis) * (sin(zrlas)*cos(zgam) + cos(zrlas)*sin(zgam))
    else 
      zarg1   = sin (zlampol) * (-zsinpol * cos(zrlas) * cos(zphis)  +   \ 
                                zcospol *              sin(zphis)) -     \
              cos (zlampol) *             sin(zrlas) * cos(zphis)
      zarg2   = cos (zlampol) * (-zsinpol * cos(zrlas) * cos(zphis)  +   \
                                zcospol *              sin(zphis)) +     \
              sin (zlampol) *             sin(zrlas) * cos(zphis)
    end if

    if (zarg2 .eq. 0.0d) then
      zarg2 = 1.0E-20
    end if
    rot2geo = zrpi18 * atan2(zarg1,zarg2)

  end if
  end if

  return(tofloat(rot2geo))

end
;-----------------------------------------------

