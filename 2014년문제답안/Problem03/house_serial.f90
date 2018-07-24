PROGRAM houses
  IMPLICIT NONE
  INTEGER, PARAMETER :: nx=17
  INTEGER :: num_sol
  INTEGER, DIMENSION(nx) :: site
  CHARACTER(10) :: date, time, zone
  INTEGER, DIMENSION(8) :: t1, t2
  REAL :: elapsed

  PRINT *, 'nx=', nx

  CALL DATE_AND_TIME(date, time, zone, t1)
  num_sol = build_house(1)
  CALL DATE_AND_TIME(date, time, zone, t2)
  elapsed = REAL(t2(5)-t1(5))*3600 + REAL(t2(6)-t1(6))*60 + &
            REAL(t2(7)-t1(7))      + REAL(t2(8)-t1(8))*1e-3   ! sec

  PRINT *, 'num_sol=', num_sol
  PRINT *, 'elapsed_time=', elapsed, 'sec'


CONTAINS
  RECURSIVE FUNCTION build_house(column) RESULT(num_sol)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: column
    INTEGER :: num_sol
    INTEGER :: i,j
    LOGICAL :: is_sol

    num_sol = 0

    ! Try to build a house in each line of column
    DO i=1,nx
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
        IF (column == nx) THEN
          ! If this is the last column, printout the solution
          num_sol = num_sol + 1
        ELSE
          ! The placement is not complete.
          ! Try to place the house on the next column
          num_sol = num_sol + build_house(column+1)
        END IF
      END IF
    END DO
  END FUNCTION
END PROGRAM




