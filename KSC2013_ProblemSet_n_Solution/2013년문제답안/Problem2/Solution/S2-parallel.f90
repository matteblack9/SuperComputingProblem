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




SUBROUTINE exchange_boundary(nx, ny, f)
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  INTEGER, INTENT(IN) :: nx, ny
  REAL(KIND=8), DIMENSION(nx,ny), INTENT(INOUT) :: f
  INTEGER :: ierr, nprocs, myrank
  INTEGER :: ireq1, ireq2, ireq3, ireq4, status(MPI_STATUS_SIZE)

  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)


  IF (myrank == 0) THEN
    CALL MPI_ISEND(f(1,ny-1), nx, MPI_REAL8, 1, 10, MPI_COMM_WORLD, ireq1, ierr)
    CALL MPI_IRECV(f(1,ny), nx, MPI_REAL8, 1, 20, MPI_COMM_WORLD, ireq2, ierr)

  ELSE IF (myrank == nprocs-1) THEN
    CALL MPI_ISEND(f(1,2), nx, MPI_REAL8, nprocs-2, 20, MPI_COMM_WORLD, ireq1, ierr)
    CALL MPI_IRECV(f(1,1), nx, MPI_REAL8, nprocs-2, 10, MPI_COMM_WORLD, ireq2, ierr)

  ELSE
    CALL MPI_ISEND(f(1,ny-1), nx, MPI_REAL8, myrank+1, 10, MPI_COMM_WORLD, ireq1, ierr)
    CALL MPI_ISEND(f(1,2), nx, MPI_REAL8, myrank-1, 20, MPI_COMM_WORLD, ireq2, ierr)
    CALL MPI_IRECV(f(1,ny), nx, MPI_REAL8, myrank+1, 20, MPI_COMM_WORLD, ireq3, ierr)
    CALL MPI_IRECV(f(1,1), nx, MPI_REAL8, myrank-1, 10, MPI_COMM_WORLD, ireq4, ierr)

    CALL MPI_WAIT(ireq3, status, ierr)
    CALL MPI_WAIT(ireq4, status, ierr)
  END IF

  CALL MPI_WAIT(ireq1, status, ierr)
  CALL MPI_WAIT(ireq2, status, ierr)
END SUBROUTINE



PROGRAM wave2d
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  INTEGER, PARAMETER :: tnx=8000, tny=8000
  INTEGER, PARAMETER :: tmax=1000
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: f, g
  REAL(KIND=8), DIMENSION(tnx,tny) :: fout
  INTEGER :: nx, ny
  INTEGER :: i,j,tstep
  INTEGER :: ierr, nprocs, myrank, rank


  !----------------------------------------------------------------------------
  ! initialize the MPI environmnet
  !----------------------------------------------------------------------------
  CALL MPI_INIT(ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)


  !----------------------------------------------------------------------------
  ! allocate the field arrays
  !----------------------------------------------------------------------------
  nx = tnx
  ny = tny/nprocs + 2
  ALLOCATE(f(nx,ny))
  ALLOCATE(g(nx,ny))


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
    IF (MOD(nprocs,2) == 0) THEN
      IF (myrank == nprocs/2) THEN
        g(nx/2,2) = SIN(0.02*tstep)
      END IF
    ELSE
      IF (myrank == nprocs/2) THEN
        g(nx/2,ny/2) = SIN(0.02*tstep)
      END IF
    END IF


    CALL advance_field(nx, ny, f, g)
    CALL exchange_boundary(nx, ny, f)

    CALL advance_field(nx, ny, g, f)
    CALL exchange_boundary(nx, ny, g)
  END DO


  !----------------------------------------------------------------------------
  ! gather fields and save as binary files
  !----------------------------------------------------------------------------
  CALL MPI_GATHER(f(1,2), nx*(ny-2), MPI_REAL8, fout, nx*(ny-2), MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

  IF (myrank == 0) THEN
    OPEN(UNIT=11, FILE='field.bin', ACCESS='DIRECT', FORM='UNFORMATTED', &
        RECL=tnx*tny*8, STATUS='NEW')
    WRITE(11,rec=1) fout
    CLOSE(11)
  END IF


  !----------------------------------------------------------------------------
  ! finalize the MPI environmnet
  !----------------------------------------------------------------------------
  CALL MPI_FINALIZE(ierr)
END PROGRAM
