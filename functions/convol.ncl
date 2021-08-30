function convol(data)
;Roberts cross convolution
begin
  print("Roberts cross convolution to find edges ....")
  ;Kernel
  kx = (/(/1,0/),(/0,-1/)/)
  ky = (/(/0,1/),(/-1,0/)/)
  ;
  ndim = dimsizes(data)
  ny   = ndim(0)
  nx   = ndim(1)
  yy   = data^(0.5)
  GradI = new((/2,ny,nx/),"float")
  do j = 0, ny-2
    jj = j+1
    do i = 0, nx-2
      ii = i + 1
      Gx   = (kx(0,0)*yy(j,i) + kx(1,1)*yy(jj,ii)) \
              + (kx(1,0)*yy(jj,i) + kx(0,1)*yy(j,ii))
      Gy   = (ky(0,0)*yy(j,i) + ky(1,1)*yy(jj,ii)) \
              + (ky(1,0)*yy(jj,i) + ky(0,1)*yy(j,ii))
      GradI(0,j,i) = (Gx^2 + Gy^2)^(0.5)
      if (Gx.ne.0) then
        GradI(1,j,i) = atan(Gy/Gx) * 180./3.14159
      end if
    end do
  end do
  GradI = where(GradI.lt.1,0,1)
  return(GradI(0,:,:))
end
