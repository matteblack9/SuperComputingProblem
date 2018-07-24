program main
    
    use mpi
    use multigrid
    implicit none
    
    real(8), parameter :: PI = 3.141592653589793
    integer(4), parameter :: nlevel = 15
    integer(4), parameter :: maxiteration = 20
    real(8), parameter :: length_x = 1.0
    real(8), parameter :: length_y = 1.0

    integer(4) :: ngrid
    integer(4) :: mlevel, gsize
    integer(4) :: numgrid_x(nlevel), numgrid_y(nlevel)
    integer(4) :: matrixDOF(nlevel), index_level(nlevel), sum_matrixDOF, index
    integer(4) :: count, ii, jj, kk, i, iter
    real(8)    :: gridsize(nlevel)
    real(8), allocatable :: rhs(:), solution(:)
    real(8)    :: coord_x, coord_y
    real(8) :: time_start, time_end, time_elapsed

    integer(4) :: mpisize, my_rank, rank
    integer(4) :: prev_rank, next_rank, distance
    integer(4) :: status(MPI_STATUS_SIZE),ierr
   
    character(256) :: myfilename

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,mpisize,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)

    time_start = MPI_WTIME()

    ngrid = 2**(nlevel-1)

    if(my_rank.eq.0) then
        prev_rank = MPI_PROC_NULL
    else
        prev_rank = my_rank - 1
    endif

    if(my_rank.eq.mpisize-1) then
        next_rank = MPI_PROC_NULL
    else
        next_rank = my_rank + 1
    endif

    numgrid_x(nlevel) = ngrid
    numgrid_y(nlevel) = ngrid/mpisize
    gridsize(nlevel) = length_x/ngrid
    matrixDOF(nlevel) = numgrid_x(nlevel)*numgrid_y(nlevel)

    sum_matrixDOF = matrixDOF(nlevel)

    gsize=mpisize
    mlevel=1

    do i=1,nlevel
        gsize=gsize/2
        mlevel=mlevel+1
        if(gsize.eq.1) exit
    enddo

    do i=nlevel,mlevel+1,-1
        numgrid_x(i-1) = numgrid_x(i)/2
        numgrid_y(i-1) = numgrid_y(i)/2
        gridsize(i-1) = length_x/numgrid_x(i-1)
        matrixDOF(i-1) = numgrid_x(i-1)*numgrid_y(i-1)
        sum_matrixDOF = sum_matrixDOF + matrixDOF(i-1)
    enddo
    
    do i=mlevel,2,-1
        numgrid_x(i-1) = numgrid_x(i)/2
        numgrid_y(i-1) = 1
        gridsize(i-1) = length_x/numgrid_x(i-1)
        matrixDOF(i-1) = numgrid_x(i-1)*numgrid_y(i-1)
        sum_matrixDOF = sum_matrixDOF + matrixDOF(i-1)
    enddo

    if(my_rank.eq.0) print*,'[Multigrid] Geometry and matrix size initialized.'

    allocate( rhs(sum_matrixDOF) )
    allocate( solution(sum_matrixDOF) )
    
    index_level(1)=1
    do i=2,nlevel
        index_level(i)=index_level(i-1)+matrixDOF(i-1)    
    enddo
    
    count = 0

    index = 2*index_level(nlevel)-1
    do jj=1,numgrid_y(nlevel)
        coord_y = length_y/dble(mpisize)*dble(my_rank)+jj*gridsize(nlevel)-0.5*gridsize(nlevel)
        do ii=1,numgrid_x(nlevel)
            coord_x = ii*gridsize(nlevel)-0.5*gridsize(nlevel)
            rhs(index_level(nlevel)+count) = sin(coord_x/length_x*PI)* &
                                             sin(coord_y/length_y*PI);
            count=count+1
        enddo 
    enddo
    if(my_rank.eq.0) print*,'[Multigrid] Geometry and rhs constructed.'
    if(my_rank.eq.0) print*,'[Multigrid] Start solving equations.'

    do iter=1,20
        if(my_rank.eq.0) write(*,*) '[Multigrid] Iteration = ',iter
        do i=nlevel,mlevel+1,-1 
            call restriction(rhs(index_level(i)),rhs(index_level(i-1)),numgrid_x(i-1),numgrid_y(i-1))
            if(my_rank.eq.0) print 102, i, i-1
        enddo
        distance=1
        do i=mlevel,2,-1 
            call restriction_mlevel(rhs(index_level(i)),rhs(index_level(i-1)),numgrid_x(i-1),numgrid_y(i-1),my_rank,distance)
            distance = distance * 2
            if(my_rank.eq.0) print 102, i, i-1
        enddo

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
            solution(index_level(1)) = -rhs(index_level(1))/4.0
            write(*,*) '[Multigrid] Solution at the coarsest level = ',solution(index_level(1))
        endif

        do i=2,mlevel
            distance = distance / 2
            call interpolation_mlevel(solution(index_level(i-1)),solution(index_level(i)), &
                                  numgrid_x(i),numgrid_y(i),my_rank,mpisize,distance)
            if(my_rank.eq.0) print 103, i-1, i
        enddo
        do i=mlevel+1,nlevel
            call interpolation_mpi(solution(index_level(i-1)),solution(index_level(i)), &
                               numgrid_x(i),numgrid_y(i),my_rank,prev_rank,next_rank)
            if(my_rank.eq.0) print 103, i-1, i
        enddo
        do i=1,matrixDOF(nlevel)
            rhs(index_level(nlevel)+i-1) = pi*pi*pi*pi*solution(index_level(nlevel)+i-1)
        enddo
    enddo

    if(my_rank.eq.0) print*, '[Multigrid] Final solution obtained.'

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    i=nlevel-5
    index = 2*index_level(i)-1

    if(my_rank.eq.(mpisize/2-1)) then 
        write(myfilename,100)
        open(40,file=myfilename,status='new', form='formatted')
        do jj=-7,0
            coord_y = length_y/dble(mpisize)*dble(my_rank)+(jj+numgrid_y(i))*gridsize(i)-0.5*gridsize(i)
            do ii=-7,8
                count = (jj+numgrid_y(i)-1)*numgrid_x(i) + ii -1 + numgrid_x(i)/2
                coord_x = (ii+numgrid_x(i)/2)*gridsize(i)-0.5*gridsize(i)
                write(40,101) coord_x,coord_y, &
                              solution(index_level(i)+count)
            enddo
        enddo
        close(40)
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    if(my_rank.eq.(mpisize/2)) then 
        write(myfilename,100)
        open(40,file=myfilename,status='old', position='append', form='formatted')
        do jj=1,8
            coord_y = length_y/dble(mpisize)*dble(my_rank)+(jj)*gridsize(i)-0.5*gridsize(i)
            do ii=-7,8
                count = (jj-1)*numgrid_x(i) + ii-1 + numgrid_x(i)/2
                coord_x = (ii+numgrid_x(i)/2)*gridsize(i)-0.5*gridsize(i)
                write(40,101) coord_x,coord_y, &
                              solution(index_level(i)+count)
            enddo
        enddo
        close(40)
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    if(my_rank.eq.0) print*, '[Multigrid] Final solution printed.'

    deallocate(rhs)
    deallocate(solution)
    if(my_rank.eq.0) print*, '[Multigrid] Memory deallocated.'
    time_end = MPI_WTIME()

    if(my_rank.eq.0) write(*,105), time_end-time_start
    

    call MPI_FINALIZE(ierr)
    
100 format('./solution.dat')
101 format(3(F12.6,x))
102 format(' [Multigrid] Restriction of RHS done (level ',I2.1,' -> level ',I2.1').')
103 format(' [Multigrid] Interpolation of solution done (level ',I2.1,' -> level ',I2.1').')
105 format(' [Multigrid] Mission completed in ',F10.6,' (secs).')
end
