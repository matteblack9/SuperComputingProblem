SUBROUTINE advance_field(nx, ny, f, g, c)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nx, ny
  REAL(KIND=8), DIMENSION(nx,ny), INTENT(IN)    :: g, c
  REAL(KIND=8), DIMENSION(nx,ny), INTENT(INOUT) :: f

  INTEGER :: i,j

  DO j=2,ny-1
    DO i=2,nx-1
      f(i,j) = c(i,j)*(g(i+1,j) + g(i-1,j) + g(i,j+1) + g(i,j-1) - 4*g(i,j)) &
               + 2*g(i,j) - f(i,j)
    END DO
  END DO
END SUBROUTINE




SUBROUTINE periodic_x(nx, ny, f)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nx, ny
  REAL(KIND=8), DIMENSION(nx,ny), INTENT(INOUT) :: f

  f(nx,:) = f(2,:)
  f(1,:) = f(nx-1,:)
END SUBROUTINE




PROGRAM wave2d
  IMPLICIT NONE
  INTEGER, PARAMETER :: nx=2000, ny=2000
  INTEGER, PARAMETER :: tmax=1500
  INTEGER, PARAMETER :: width=120, thick=60, gap=1200   ! slit parameters
  INTEGER, PARAMETER :: distance=800
  REAL(KIND=8), PARAMETER :: c0=0.49D0
  REAL(KIND=8), DIMENSION(nx,ny) :: f, g, c
  INTEGER :: i,j,tstep


  !----------------------------------------------------------------------------
  ! initialize coefficient and fields
  !----------------------------------------------------------------------------
  DO j=1,ny
    DO i=1,nx
      f(i,j) = 0D0
      g(i,j) = 0D0
      c(i,j) = c0
    END DO
  END DO


  !----------------------------------------------------------------------------
  ! slit structure
  !----------------------------------------------------------------------------
  c(:,ny/2:ny/2+thick) = 0D0
  c(nx/2-(gap+width)/2:nx/2+(gap+width)/2,ny/2:ny/2+thick) = c0
  c(nx/2-(gap-width)/2:nx/2+(gap-width)/2,ny/2:ny/2+thick) = 0D0
  

  !----------------------------------------------------------------------------
  ! main loop for the time evolution
  !----------------------------------------------------------------------------
  DO tstep=1,tmax
    !------------------------------------
    ! point source
    !------------------------------------
    ! IF (myrank == 0) THEN
    !   g(200,ny-200) = SIN(0.02*tstep)
    ! END IF


    !------------------------------------
    ! line source
    !------------------------------------
    g(:,ny/2-distance) = SIN(0.02*tstep)


    CALL advance_field(nx, ny, f, g, c)
    CALL periodic_x(nx, ny, f)

    CALL advance_field(nx, ny, g, f, c)
    CALL periodic_x(nx, ny, g)
  END DO


  !----------------------------------------------------------------------------
  ! save as binary files
  !----------------------------------------------------------------------------
  OPEN(UNIT=11, FILE='field.bin', ACCESS='DIRECT', FORM='UNFORMATTED', &
      RECL=nx*ny*8, STATUS='NEW')
  WRITE(11,rec=1) f
  CLOSE(11)
END PROGRAM
