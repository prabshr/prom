 SUBROUTINE readascii(nx,ny,nz,avar,ncols,scol,fname)

 IMPLICIT NONE
!
 CHARACTER*(*),   INTENT(IN)        :: fname
 INTEGER(KIND=4), INTENT(IN)        :: nx,ny,nz,ncols,scol
 REAL,            INTENT(OUT)       :: avar(nx,ny,nz,ncols)
 CHARACTER(len=250)                 :: header
 INTEGER(KIND=4)                    :: i,j,k,iunit,ie,je,ke,nc
 REAL                               :: nvar(ncols)
! 
! 
 iunit = 100

 OPEN(iunit,file=trim(fname),status='old') 

 header = REPEAT(' ', LEN(header))
 READ (iunit, '(A)') header
 READ (iunit, *) ie, je, ke 

 !read whole array READ (iunit, *) avar 
 DO k=1,nz
   DO j=1,ny
     DO i=1,nx
       IF (ncols == 1) THEN
         READ(iunit,*) avar(i,j,k,ncols)
       ELSEIF (ncols == 8) THEN
         READ(iunit,*) (nvar(nc),nc=1,ncols)
         avar(i,j,k,:) = nvar !;(scol)
       ELSE
         WRITE(*,*) "CODE NOT WRITTEN FOR ncols=",ncols
         CALL ABORT
       ENDIF
     END DO
   END DO
 END DO

 CLOSE(iunit)
 RETURN
 END SUBROUTINE readascii
