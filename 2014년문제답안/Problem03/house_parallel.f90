PROGRAM houses
  USE MPI
  IMPLICIT NONE
  INTEGER, PARAMETER :: nx=17
  INTEGER, PARAMETER :: GRANULARITY=5   ! it should be less than nx

  TYPE :: JOB
      SEQUENCE
      LOGICAL :: working        ! true: working, false: quit
      INTEGER, DIMENSION(GRANULARITY) :: site
  END TYPE

  TYPE :: JOB_MSG
      SEQUENCE
      INTEGER :: solutions_found
      INTEGER :: origin
  END TYPE

  INTEGER, DIMENSION(nx) :: site
  INTEGER :: total_num_sol
  CHARACTER(10) :: date, time, zone     ! measure time
  INTEGER, DIMENSION(8) :: t1, t2
  REAL :: elapsed
  INTEGER :: nprocs, myrank, ierr, stat(MPI_STATUS_SIZE)      ! MPI variables


  CALL MPI_INIT(ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)


  IF (myrank == 0) THEN
    PRINT *, 'nx=', nx
    PRINT *, 'nprocs=', nprocs

    CALL DATE_AND_TIME(date, time, zone, t1)
    total_num_sol = master_build_house(1)
    total_num_sol = total_num_sol + wait_remaining_results()
    CALL DATE_AND_TIME(date, time, zone, t2)
    elapsed = REAL(t2(5)-t1(5))*3600 + REAL(t2(6)-t1(6))*60 + &
              REAL(t2(7)-t1(7))      + REAL(t2(8)-t1(8))*1e-3   ! sec

    PRINT *, 'num_sol=', total_num_sol
    PRINT *, 'elapsed_time=', elapsed, 'sec'

  ELSE
    CALL worker()
  END IF

  CALL MPI_FINALIZE(ierr)


CONTAINS
  RECURSIVE FUNCTION master_build_house(column) RESULT(num_sol)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: column
    LOGICAL :: is_sol
    INTEGER :: num_sol, i, j
    
    num_sol = 0

    DO i=1,nx
      ! Try to build a house in each line of column
      site(column) = i

      ! Check if this placement is still a solution
      is_sol = .TRUE.
      DO j=column-1,1,-1
        IF (site(column) == site(j) .OR. &
            site(column) == site(j)-(column-j) .OR. &
            site(column) == site(j)+(column-j)) THEN
          is_sol = .FALSE.
          EXIT
        END IF
      END DO

      IF (is_sol) THEN
        IF (column == GRANULARITY) THEN
          ! If this is the last level (granularity of the job),
          ! send a next job to a worker
          num_sol = num_sol + send_job_worker()
        ELSE   
          ! The placement is not complete.
          ! Try to place the house on the next column
          num_sol = num_sol + master_build_house(column+1)
        END IF
      END IF
    END DO
  END FUNCTION



  RECURSIVE FUNCTION worker_build_house(column, sub_site) RESULT(num_sol)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: column
    INTEGER, DIMENSION(GRANULARITY), INTENT(OUT) :: sub_site
    LOGICAL :: is_sol
    INTEGER :: num_sol, i, j

    num_sol = 0

    DO i=1,nx
      ! Try to build a house in each line of column
      sub_site(column) = i

      ! Check if this placement is still a solution
      is_sol = .TRUE.
      DO j=column-1,1,-1
        IF (sub_site(column) == sub_site(j) .OR. &
            sub_site(column) == sub_site(j)-(column-j) .OR. &
            sub_site(column) == sub_site(j)+(column-j)) THEN
          is_sol = .FALSE.
          EXIT
        END IF
      END DO

      IF (is_sol) THEN
        IF (column == nx) THEN
          ! If this is the last column, printout the solution
          num_sol = num_sol + 1
        ELSE   
          ! The placement is not complete.
          ! Try to place the house on the next column
          num_sol = num_sol + worker_build_house(column+1, sub_site)
        END IF
      END IF
    END DO
  END FUNCTION



  FUNCTION send_job_worker() RESULT(num_sol)
    IMPLICIT NONE
    LOGICAL :: is_sol
    INTEGER :: num_sol, i, j
    TYPE(JOB) :: todo
    TYPE(JOB_MSG) :: msg

    num_sol = 0     ! The number of solutions found meanwhile

    ! Set the job
    todo%working = .TRUE.

    DO i=1,GRANULARITY
      todo%site(i) = site(i)
    END DO

    ! Recieve the last result from a worker
    CALL MPI_RECV(msg, SIZEOF(msg), MPI_BYTE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, stat, ierr)

    num_sol = msg%solutions_found

    ! Send the new job to the worker
    CALL MPI_SEND(todo, SIZEOF(todo), MPI_BYTE, msg%origin, 1, MPI_COMM_WORLD, ierr)
  END FUNCTION



  FUNCTION wait_remaining_results() RESULT(num_sol)
    ! Wait for remaining results, sending a quit whenever a new result arrives
    IMPLICIT NONE
    INTEGER :: num_sol, n_workers
    TYPE(JOB) :: todo
    TYPE(JOB_MSG) :: msg

    num_sol = 0
    n_workers = nprocs-1

    todo%working = .FALSE.

    DO WHILE (n_workers > 0)
      ! Receive a message from a worker
      CALL MPI_RECV(msg, SIZEOF(msg), MPI_BYTE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, stat, ierr)
      num_sol = num_sol + msg%solutions_found

      CALL MPI_SEND(todo, SIZEOF(todo), MPI_BYTE, msg%origin, 1, MPI_COMM_WORLD, ierr)

      n_workers = n_workers - 1;
    END DO
  END FUNCTION



  SUBROUTINE worker()
    ! There is a default message which lets a worker request a 
    ! job reporting the number of solutions found in the last iteration
    IMPLICIT NONE
    INTEGER :: num_sol, i
    TYPE(JOB) :: todo
    TYPE(JOB_MSG) :: msg

    msg%solutions_found = 0
    msg%origin = myrank

    ! Request initial job
    CALL MPI_SEND(msg, SIZEOF(msg), MPI_BYTE, 0, 0, MPI_COMM_WORLD, ierr)

    DO WHILE (.TRUE.)
      ! Wait for a job or a quit message
      CALL MPI_RECV(todo, SIZEOF(todo), MPI_BYTE, 0, 1, MPI_COMM_WORLD, stat, ierr)

      IF (.NOT. todo%working) EXIT

      num_sol = worker_build_house(GRANULARITY+1, todo%site)

      ! Ask for more work
      msg%solutions_found = num_sol
      msg%origin = myrank
      CALL MPI_SEND(msg, SIZEOF(msg), MPI_BYTE, 0, 0, MPI_COMM_WORLD, ierr)
    END DO
  END SUBROUTINE
END PROGRAM
