PROGRAM MAIN 

include 'mpif.h'

    integer :: m_size = 3200
    integer :: d_size = 30
    integer :: debug_step = 0

    real(8), allocatable,dimension(:) :: cov_numbers
    real(8), allocatable,dimension(:) :: z_numbers
    real(8), allocatable,dimension(:) :: temp, c2, c3

    integer :: nprocs, myrank, ierr
    integer :: i,j,k
    integer :: start_index, last_index, my_size
    real(8) :: elapsed_time = 0.0
    integer :: status(MPI_STATUS_SIZE)

    CALL mpi_init(ierr)
    CALL mpi_comm_size(MPI_COMM_WORLD, nprocs, ierr)
    CALL mpi_comm_rank(MPI_COMM_WORLD, myrank, ierr)

    if( MOD(m_size,nprocs) /= 0) then
       if(myrank==0) print*, "m_size must be an integer multiple of nprocs!"
       CALL MPI_Finalize(ierr)
       STOP
    endif

    elapsed_time = elapsed_time - MPI_Wtime()
    start_index = m_size*myrank/nprocs + 1
    last_index = m_size*(myrank+1)/nprocs
    my_size = last_index - start_index + 1

!Step 1 : Generate covariance matrix & standard normal random number files
    CALL cov_gen(cov_numbers, start_index, my_size, m_size, myrank, nprocs)
    if(myrank == 0) print*, "Step",debug_step,"completed - generation of a covariance matrix"

    CALL z_gen(z_numbers, start_index, my_size, d_size);
    if(myrank == 0) print*, "Step",debug_step+1,"completed - generation of a Z matrix"


!Step 2 : Cholesky Decomposition of Sigma-matrix 
    if(myrank /= 0) allocate(temp(0:(start_index-1)*m_size-1))

    if(nprocs > 1) then
        if(myrank == 0) then
            CALL cholesky(c2, cov_numbers, start_index, last_index, m_size, myrank, temp)
            do ii=1, nprocs-1 
                CALL MPI_Send(c2, my_size*m_size, MPI_REAL8, ii, 1002, MPI_COMM_WORLD, ierr)
            enddo
        else
            CALL MPI_Recv(temp, my_size*m_size, MPI_REAL8, 0, 1002, MPI_COMM_WORLD, status, ierr)
        endif
 
        do jj=1, nprocs-2
            if(myrank == jj) then
                CALL cholesky(c2, cov_numbers, start_index, last_index, m_size, myrank, temp)
                do ii=myrank+1, nprocs-1
                   CALL MPI_Send(c2, my_size*m_size, MPI_REAL8, ii, 1002, MPI_COMM_WORLD, ierr)
                enddo
            else if(myrank > jj) then
                CALL MPI_Recv(temp(my_size*jj*m_size), my_size*m_size, MPI_REAL8, jj, 1002, &
                              MPI_COMM_WORLD, status, ierr)
            endif

            CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
        enddo 
    endif 

    if(myrank == nprocs-1) then
        CALL cholesky(c2, cov_numbers, start_index, last_index, m_size, myrank, temp)
    endif

    if(myrank /= 0) deallocate(temp)

    if(myrank == 0) print*, "Step",debug_step+2,"completed - Cholesky decomposition"

!Step 3 : Correlated Random Number Generator: Matrix Multiplication
    CALL mat_mul(c3, c2, z_numbers, m_size, d_size, myrank, nprocs, my_size)
    if(myrank==0) print*, "Step",debug_step+3,"completed - matrix multiplication"

!Step 4 : Write result files*/
    if(myrank==0) then
      CALL system("rm -rf result")
      CALL system("mkdir result")
    endif

    CALL file_write("./result/result.txt", my_size, d_size, c3, myrank, nprocs)
    if(myrank==0) print*, "Step",debug_step+4,"completed - generation of results"

    deallocate(cov_numbers)
    deallocate(z_numbers)
    deallocate(c2)
    deallocate(c3)

    if(myrank==0) print*, "Step",debug_step+5,"completed - memory release"

    elapsed_time = elapsed_time + MPI_Wtime()

    if(myrank==0) print*, "Total elapsed time = ", elapsed_time, " (secs)"
   
    CALL MPI_Finalize(ierr)

CONTAINS

!File write
SUBROUTINE file_write(f_name, m_size, d_size, matrix, myrank, nprocs)
    integer :: m_size, d_size
    real(8), dimension(0:m_size*d_size-1) :: matrix
    character(*):: f_name
    integer :: myrank, nprocs, i, j, ierr
   
!    OPEN(UNIT=30,FILE=TRIM(f_name),STATUS="NEW",ACTION="WRITE")

    do i=0, nprocs-1
       if(myrank==i) then
          if(i==0) then
            OPEN(UNIT=30,FILE=TRIM(f_name),STATUS="NEW",ACTION="WRITE")
          else
            OPEN(UNIT=30,FILE=TRIM(f_name),STATUS="OLD",ACTION="WRITE",POSITION="APPEND")
          endif

          do j=0,m_size*d_size-1
             write(30,10) matrix(j)
          enddo

          close(30)
       endif
       CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
    enddo

 10 format(F10.6)
   
END SUBROUTINE file_write

! Covariance matrix generate 
SUBROUTINE cov_gen(L, start_index, my_n, n, myrank, nprocs)
    integer :: ierr
    integer :: i, j, my_n, n, row_offset
    integer :: start_index, myrank, nprocs
    real(8), allocatable, dimension(:) :: L !output

    allocate(L(0:my_n*n-1),STAT=ierr);
    if (ierr /= 0) then
        print*, "ALLOCATION Failed:cov_gen"
        STOP
    endif

    do i=0, my_n*n-1
       L(i) = 0.0
    enddo

    row_offset = start_index - 1

    do i=0,my_n-1
        L(i*n + i + row_offset) = log(dble(n)) * log(dble(n)) - cos(dble((i + row_offset) * n))
        if(.NOT.(i==0 .AND. myrank == 0)) then
            L(i*n + i + row_offset - 1) = sin(dble((i+row_offset) * n))
        endif
        if(.NOT.(i==my_n-1 .AND. myrank == nprocs-1)) then
            L(i*n + i + row_offset + 1) = sin(dble((i+row_offset+1) * n))
        endif
    enddo

END SUBROUTINE cov_gen

! Multivariate number generate 
SUBROUTINE z_gen(L, start_index, my_n, d)
    integer :: i, j, row_offset
    integer :: my_n, d, start_index, ierr
    real(8), allocatable, dimension(:) :: L

    allocate(L(0:my_n*d-1),STAT=ierr);
    if (ierr /= 0) then
        print*, "ALLOCATION Failed:z_gen"
        STOP
    endif

    do i=0, my_n*d-1
       L(i) = 0.0
    enddo

    row_offset = start_index - 1

    do i=0, my_n-1
       do j=0, d-1
          L(i * d + j) = exp( - cos(dble((row_offset+i+1) * j)) * sin(dble((row_offset+i+1)* (j+1))) )
       enddo
    enddo

END SUBROUTINE z_gen

!Cholesky Decomposition
SUBROUTINE cholesky(rslt, A, start_index, last_index, n, myrank, prslt) 
    integer :: i, j, k, n, ierr(2), row_offset
    integer :: start_index, last_index, myrank
    real(8) :: s
    real(8), dimension(0:(last_index-start_index+1)*n-1) :: A !cov_numbers
    real(8), allocatable, dimension(:) :: L, rslt
    real(8) :: prslt(0:(start_index-1)*n-1)

    allocate(L(0:last_index*n-1),STAT=ierr(1))
    allocate(rslt(0:(last_index-start_index+1)*n-1),STAT=ierr(2))
    if (ierr(1) /= 0 .OR. ierr(2) /= 0) then
        print*, "ALLOCATION Failed:cholesky"
        CALL MPI_Finalize(ierr(1))
        STOP
    endif

    do i = 0, last_index*n-1 
       L(i) = 0.0
    enddo
    
    do i = 0, (last_index-start_index+1)*n-1-1
        rslt(i) = 0.0
    enddo
        
    if(myrank /= 0) L(0:) = prslt(0:(start_index-1)*n-1)
  
    row_offset = start_index - 1
 
    do i = row_offset, last_index-1
       do j=0,i
          s=0.
          do k=0,j-1
             s = s + L(i * n + k) * L(j * n + k)
          enddo
          if(i==j) then
             L(i*n+j) = SQRT(A((i-row_offset)*n+i)-s)
          else 
             L(i*n+j) = (1.0 / L(j * n + j) * (A((i-row_offset) * n + j) - s))
          endif
       enddo
    enddo

    rslt(0:) = L(row_offset*n:row_offset*n+(last_index-start_index+1)*n - 1)

    deallocate(L)
END SUBROUTINE cholesky

!Matrix multiplication 
SUBROUTINE mat_mul(L, A, B, n, m, myrank, nprocs, mysize)
    integer :: i, j, k, n, m, ierr(3)
    real(8), dimension(0:mysize*n-1) :: A
    real(8), dimension(0:mysize*m-1) :: B
    real(8), allocatable, dimension(:) :: L, Bl, Br
    integer :: myrank, nprocs, mysize
    integer :: leftrank, rightrank
    integer :: req(2), status(MPI_STATUS_SIZE)
  
    allocate(L(0:mysize*m-1),STAT=ierr(1))
    allocate(Bl(0:mysize*m-1),STAT=ierr(2))
    allocate(Br(0:mysize*m-1),STAT=ierr(3))
    if (ierr(1) /= 0 .OR. ierr(2) /= 0 .OR. ierr(3) /= 0) then
        print*, "ALLOCATION Failed:mat_mul"
        CALL MPI_Finalize(ierr(1))
        STOP
    endif

    do i=0, mysize*m-1
       L(i) = 0.
       Bl(i) = 0.
       Br(i) = 0.
    enddo

    leftrank = MOD((myrank-1+nprocs),nprocs)
    rightrank = MOD((myrank+1),nprocs)

    if(nprocs > 1) then
       CALL MPI_Irecv(Bl, mysize*m, MPI_REAL8, rightrank, 1002, MPI_COMM_WORLD, req(1), ierr(1))
       CALL MPI_Isend(B, mysize*m, MPI_REAL8, leftrank, 1002, MPI_COMM_WORLD, req(2), ierr(2))
    endif

    do i=0,mysize-1
      do j=0,m-1
        do k=0,mysize-1
           L(i * m + j) = L(i * m +j) + A(i * n + k + mysize*myrank) * B(k * m + j)
        enddo
      enddo
    enddo

    if(nprocs > 1) CALL MPI_Waitall(2, req, status,ierr(1))

    if(nprocs > 2) then
       CALL MPI_Irecv(Br, mysize*m, MPI_REAL8, leftrank, 1002, MPI_COMM_WORLD, req(1), ierr(1))
       CALL MPI_Isend(B, mysize*m, MPI_REAL8, rightrank, 1002, MPI_COMM_WORLD, req(2), ierr(2))
    endif

    if(nprocs > 1) then
       do i=0,mysize-1
         do j=0,m-1
           do k=0,mysize-1
              L(i * m + j) = L(i * m +j) + A(i * n + k + mysize*rightrank) * Bl(k * m + j)
           enddo
         enddo
       enddo
    endif 

    if(nprocs > 2) CALL MPI_Waitall(2, req, status,ierr(1))

    if(nprocs > 2) then
       do i=0,mysize-1
         do j=0,m-1
           do k=0,mysize-1
              L(i * m + j) = L(i * m +j) + A(i * n + k + mysize*leftrank) * Br(k * m + j)
           enddo
         enddo
       enddo
    endif

    deallocate(Bl)
    deallocate(Br)

END SUBROUTINE mat_mul

SUBROUTINE show_matrix(A, n, m) 
    integer :: i, n, m
    
    real(8), dimension(0:n*m-1) :: A

    do i=0,n-1
       write(500,10) (A(i * m + j),j=0,m-1)
    enddo
 10 format(F10.6)

END SUBROUTINE show_matrix

END PROGRAM MAIN



