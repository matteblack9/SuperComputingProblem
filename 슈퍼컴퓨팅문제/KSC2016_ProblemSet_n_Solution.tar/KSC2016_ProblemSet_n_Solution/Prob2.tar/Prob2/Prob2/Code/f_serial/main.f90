program main
    
    use multigrid
    implicit none
    
    real(8), parameter :: PI = 3.141592653589793
    integer(4), parameter :: nlevel = 15
    integer(4), parameter :: maxiteration = 20
    real(8), parameter :: length_x = 1.0
    real(8), parameter :: length_y = 1.0

    integer(4) :: ngrid
    integer(4) :: numgrid_x(nlevel), numgrid_y(nlevel)
    integer(4) :: matrixDOF(nlevel), index_level(nlevel), sum_matrixDOF
    integer(4) :: count, ii, jj, kk, i, iter
    real(8)    :: gridsize(nlevel)
    real(8), allocatable :: rhs(:), solution(:)
    real(8)    :: coord_x, coord_y
    real(8) :: time_start, time_end

   
    character(256) :: myfilename
    call CPU_TIME(time_start)

    ngrid = 2**(nlevel-1)

    numgrid_x(nlevel) = ngrid
    numgrid_y(nlevel) = ngrid
    gridsize(nlevel) = length_x/ngrid
    matrixDOF(nlevel) = numgrid_x(nlevel)*numgrid_y(nlevel)

    sum_matrixDOF = matrixDOF(nlevel)

    do i=nlevel,2,-1
        numgrid_x(i-1) = numgrid_x(i)/2
        numgrid_y(i-1) = numgrid_y(i)/2
        gridsize(i-1) = length_x/numgrid_x(i-1)
        matrixDOF(i-1) = numgrid_x(i-1)*numgrid_y(i-1)
        sum_matrixDOF = sum_matrixDOF + matrixDOF(i-1)
    enddo
    
    print*,'[Multigrid] Geometry and matrix size initialized.'

    allocate( rhs(sum_matrixDOF) )
    allocate( solution(sum_matrixDOF) )
    
    index_level(1)=1
    do i=2,nlevel
        index_level(i)=index_level(i-1)+matrixDOF(i-1)    
    enddo
    
    count = 0

    do jj=1,numgrid_y(nlevel)
        coord_y = jj*gridsize(nlevel)-0.5*gridsize(nlevel)
        do ii=1,numgrid_x(nlevel)
            coord_x = ii*gridsize(nlevel)-0.5*gridsize(nlevel)
            rhs(index_level(nlevel)+count) = sin(coord_x/length_x*PI)* &
                                             sin(coord_y/length_y*PI)
            count=count+1
        enddo 
    enddo
    print*,'[Multigrid] Geometry and rhs constructed.'
    print*,'[Multigrid] Start solving equations.'

    do iter = 1, maxiteration
        write(*,*) '[Multigrid] Iteration = ', iter
        do i=nlevel,2,-1 
            call restriction(rhs(index_level(i)),rhs(index_level(i-1)),numgrid_x(i-1))
            print 102, i, i-1
        enddo

        solution(index_level(1)) = -rhs(index_level(1))/4.0
        write(*,*) '[Multigrid] Solution at the coarsest level = ',solution(index_level(1))

        do i=2,nlevel
            call interpolation(solution(index_level(i-1)),solution(index_level(i)),numgrid_x(i))
            print 103, i-1, i
        enddo
        do i=1,matrixDOF(nlevel)
            rhs(index_level(nlevel)+i-1) = pi*pi*pi*pi*solution(index_level(nlevel)+i-1)
        enddo

    enddo

    print*, '[Multigrid] Final solution obtained.'

    i=nlevel-5
    write(myfilename,100)
    open(40,file=myfilename,status='unknown', form='formatted')
    do jj=-7,8
        coord_y = (jj+numgrid_y(i)/2)*gridsize(i)-0.5*gridsize(i)
        do ii=-7,8
            count = (jj-1+numgrid_y(i)/2)*numgrid_x(i) + ii-1 + numgrid_x(i)/2
            coord_x = (ii+numgrid_x(i)/2)*gridsize(i)-0.5*gridsize(i)
            write(40,101) coord_x,coord_y, &
                          solution(index_level(i)+count)
        enddo
    enddo
    
    close(40)

    print*, '[Multigrid] Final solution printed.'

    deallocate(rhs)
    deallocate(solution)
    print*, '[Multigrid] Memory deallocated.'

    call CPU_TIME(time_end)
    write(*,105) time_end-time_start
    
100 format('./solution.dat')
101 format(3(F12.6,x))
102 format(' [Multigrid] Restriction of RHS done (level ',I2.1,' -> level ',I2.1').')
103 format(' [Multigrid] Interpolation of solution done (level ',I2.1,' -> level ',I2.1').')
105 format(' [Multigrid] Mission completed in ',F10.6,' (secs).')

end
