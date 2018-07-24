SUBROUTINE advance_field(nx, ny, f, g)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nx, ny
  REAL(KIND=8), DIMENSION(nx,ny), INTENT(IN)    :: g
  REAL(KIND=8), DIMENSION(nx,ny), INTENT(INOUT) :: f

  INTEGER :: i,j

  DO j=2,ny-1
    DO i=2,nx-1
      f(i,j) = 0.49D0*(g(i+1,j) + g(i-1,j) + g(i,j+1) + g(i,j-1) - 4*g(i,j)) &
               + 2*g(i,j) - f(i,j)
    END DO
  END DO
END SUBROUTINE




PROGRAM wave2d
  IMPLICIT NONE
  INTEGER, PARAMETER :: nx=2000, ny=2000
  INTEGER, PARAMETER :: tmax=500
  REAL(KIND=8), DIMENSION(nx,ny) :: f,g
  INTEGER :: i,j,tstep


  !----------------------------------------------------------------------------
  ! initialize coefficient and fields
  !----------------------------------------------------------------------------
  DO j=1,ny
    DO i=1,nx
      f(i,j) = 0D0
      g(i,j) = 0D0
    END DO
  END DO


  !----------------------------------------------------------------------------
  ! main loop for the time evolution
  !----------------------------------------------------------------------------
  DO tstep=1,tmax
    !------------------------------------
    ! point source
    !------------------------------------
    g(nx/2,ny/2) = SIN(0.02*tstep)


    CALL advance_field(nx, ny, f, g)
    CALL advance_field(nx, ny, g, f)
  END DO


  !----------------------------------------------------------------------------
  ! gather fields and save as binary files
  !----------------------------------------------------------------------------
  OPEN(UNIT=11, FILE='field.bin', ACCESS='DIRECT', FORM='UNFORMATTED', &
      RECL=nx*ny*8, STATUS='NEW')
  WRITE(11,rec=1) f
  CLOSE(11)
END PROGRAM
