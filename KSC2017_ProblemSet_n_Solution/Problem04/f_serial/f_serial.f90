
PROGRAM MAIN 
    ! Step 1 : Generate covariance matrix & standard normal random number files 

    integer :: m_size = 3200
    integer :: d_size = 30
    integer :: debug_step = 0

    real(8), allocatable,dimension(:) :: cov_numbers
    real(8), allocatable,dimension(:) :: z_numbers
    real(8), allocatable,dimension(:) :: c2, c3

!Step 1 : Generate covariance matrix & standard normal random number files
    CALL cov_gen(cov_numbers, m_size)
    print*, "Step",debug_step,"completed - generation of a covariance matrix"
    CALL z_gen(z_numbers, m_size, d_size)
    print*, "Step",debug_step+1,"completed - generation of a Z matrix"

!Step 2 : Cholesky Decomposition of Sigma-matrix 
    CALL cholesky(c2, cov_numbers, m_size)
    print*, "Step",debug_step+2,"completed - Cholesky decomposition"

!Step 3 : Correlated Random Number Generator: Matrix Multiplication
    CALL mat_mul(c3, c2, z_numbers, m_size, d_size)
    print*, "Step",debug_step+3,"completed - matrix multiplication"

!Step 4 : Write result files*/
    CALL system("rm -rf result")
    CALL system("mkdir result")

    CALL file_write("./result/result.txt", m_size, d_size, c3)
    print*, "Step",debug_step+4,"completed - generation of results"

    deallocate(cov_numbers)
    deallocate(z_numbers)
    deallocate(c2)
    deallocate(c3)

    print*, "Step",debug_step+5,"completed - memory release"

CONTAINS

!File write
SUBROUTINE file_write(f_name, m_size, d_size, matrix)
    integer :: m_size, d_size
    real(8), dimension(0:m_size*d_size-1) :: matrix
    character(*):: f_name

    OPEN(UNIT=30,FILE=TRIM(f_name),ACTION="WRITE")

    do i=0,m_size*d_size-1
       write(30,10) matrix(i)
    enddo

 10 format(F10.6)
   
    close(30)
 
END SUBROUTINE file_write

! Covariance matrix generate 
SUBROUTINE cov_gen(L,n)
    integer :: ierr
    integer :: i, j, n
    real(8), allocatable, dimension(:) :: L

    allocate(L(0:n*n-1),STAT=ierr)
    if (ierr /= 0) then
        print*, "ALLOCATION Failed:cov_gen"
        STOP
    endif

    DO i=0,n-1
        L(i*n + i) = log(dble(n)) * log(dble(n)) - cos(dble(i * n))
    enddo
      
    do j=1,n-1
        L((j-1)*n+j) = sin(dble(j * n))
        L(j*n+j-1) = L((j-1)*n+j)
    enddo

END SUBROUTINE cov_gen

! Multivariate number generate 
SUBROUTINE z_gen(L, n, d)
    integer :: i, j, n, d, ierr
    real(8), allocatable, dimension(:) :: L

    allocate(L(0:n*d-1),STAT=ierr)
    if (ierr /= 0) then
        print*, "ALLOCATION Failed:z_gen"
        STOP
    endif

    do i=0, n-1
       do j=0, d-1
          L(i * d + j) = exp( - cos(dble((i+1) * j)) * sin(dble((i+1)* (j+1))) )
       enddo
    enddo

END SUBROUTINE z_gen

!Cholesky Decomposition
SUBROUTINE cholesky(L, A, n) 
    integer :: i, j, k, n, ierr
    real(8) :: s
    real(8), dimension(0:n*n-1) :: A
    real(8), allocatable, dimension(:) :: L

    allocate(L(0:n*n-1),STAT=ierr)
    if (ierr /= 0) then
        print*, "ALLOCATION Failed:cholesky"
        STOP
    endif

    do i=0,n-1
       do j=0,i
          s=0
          do k=0,j-1
             s = s + L(i * n + k) * L(j * n + k)
          enddo
          if(i==j) then
             L(i*n+j) = SQRT(A(i*n+i)-s)
          else 
             L(i*n+j) = (1.0 / L(j * n + j) * (A(i * n + j) - s))
          endif 
       enddo
    enddo

END SUBROUTINE cholesky

!Matrix multiplication 
SUBROUTINE mat_mul(L, A, B, n, m)
    integer :: i, j, k, n, m, ierr
    real(8), dimension(0:n*n-1) :: A
    real(8), dimension(0:n*m-1) :: B
    real(8), allocatable, dimension(:) :: L

    allocate(L(0:n*m-1),STAT=ierr)
    if (ierr /= 0) then
        print*, "ALLOCATION Failed:mat_mul"
        STOP
    endif

    do i=0,n-1
      do j=0,m-1
        do k=0,n-1
           L(i * m + j) = L(i * m +j) + A(i * n + k) * B(k * m + j)
        enddo
      enddo
    enddo

END SUBROUTINE mat_mul

SUBROUTINE show_matrix(A, n, m) 
    integer :: i, n, m
    
    real(8), dimension(0:n*m-1) :: A

    do i=0,n-1
       print*, (A(i * m + j),j=0,m-1)
    enddo

END SUBROUTINE show_matrix

END PROGRAM MAIN



