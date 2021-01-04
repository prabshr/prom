 SUBROUTINE readdim(nx,ny,nz,fname)

 IMPLICIT NONE
!
 CHARACTER*(*),   INTENT(IN)        :: fname
 INTEGER(KIND=4), INTENT(OUT)       :: nx,ny,nz
 CHARACTER(len=250)                 :: header
 INTEGER(KIND=4)                    :: iunit
! 
! 
 iunit = 100

 OPEN(iunit,file=trim(fname),status='old') 

 header = REPEAT(' ', LEN(header))
 READ (iunit, '(A)') header
 WRITE (*,*) TRIM(header)
 READ (iunit, *) nx, ny, nz 

 CLOSE(iunit)
 RETURN
 END SUBROUTINE readdim
