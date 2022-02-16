function gfct(x)
;Obtained from SB2M COSMOv5.1
;P. Shrestha
begin

;!*******************************************************************************
;!                                                                              *
;!       Gammafunktion  Numerical Recipes (F77)  SB2M                           *
;!                                                                              *
;!*******************************************************************************

  cof = (/76.18009173d0,-86.50532033d0,24.01409822d0,  \ 
            -1.231739516d0,.120858003d-2,-.536382d-5/)
  stp = 2.50662827465d0
  half = 0.5d0
  one  = 1.0d0
  fpf  = 5.5d0

  xx  = x  - one
  tmp = xx + fpf
  tmp = (xx + half) * log(tmp) - tmp
  ser = one
  do j = 1,6
    xx  = xx  + one
    ser = ser + cof(j-1) / xx
  end do
  gamma1 = tmp + log(stp*ser)
  gamma1 = exp(gamma1)

  return(gamma1) 
end 

