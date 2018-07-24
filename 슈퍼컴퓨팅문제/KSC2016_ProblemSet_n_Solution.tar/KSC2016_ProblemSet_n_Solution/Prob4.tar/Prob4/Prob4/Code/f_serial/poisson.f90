program poisson
    
    use multigrid
    use m_cgsolver
    use matrixconstructor
    
    implicit none
    
    real(8), parameter :: PI = 3.141592653589793
    integer(4) :: nlevel, initial_ngrid

    integer(4), allocatable :: numgrid(:), matrixDOF(:), index_level(:), index2_level(:)
    integer(4) :: sum_matrixDOF, sum2_matrixDOF, index
    integer(4) :: count, ii, jj, kk, i, maxiteration
    real(8)    ::length, tolerance
    real(8), allocatable :: gridsize(:), poissonmatrix(:),rhs(:), solution(:), coordinate(:)
    real(8) :: time_start, time_end
   
    character(LEN=256) :: myfilename

    nlevel = 3
    initial_ngrid = 128
    maxiteration = 20000
    tolerance = 1.0e-10
    length = 1.0

    call CPU_TIME(time_start)

    allocate( numgrid(nlevel) )
    allocate( matrixDOF(nlevel) )
    allocate( gridsize(nlevel) )
    allocate( index_level(nlevel) )
    allocate( index2_level(nlevel) )

    numgrid(nlevel) = initial_ngrid
    gridsize(nlevel) = length/initial_ngrid
    matrixDOF(nlevel) = numgrid(nlevel)*numgrid(nlevel)

    sum_matrixDOF = matrixDOF(nlevel)
    sum2_matrixDOF = matrixDOF(nlevel)*matrixDOF(nlevel)

    do i=nlevel,2,-1
        numgrid(i-1) = numgrid(i)/2
        gridsize(i-1) = length/numgrid(i-1)
        matrixDOF(i-1) = numgrid(i-1)*numgrid(i-1)
        sum_matrixDOF = sum_matrixDOF + matrixDOF(i-1)
        sum2_matrixDOF = sum2_matrixDOF + matrixDOF(i-1)*matrixDOF(i-1)
    enddo
    
    print*,'[Poisson] Geometry and matrix size initialized.'

    allocate( rhs(sum_matrixDOF) )
    allocate( solution(sum_matrixDOF) )
    allocate( coordinate(sum_matrixDOF*2) )
    allocate( poissonmatrix(sum2_matrixDOF) )
    
    index_level(1)=1
    index2_level(1)=1
    do i=2,nlevel
        index_level(i)=index_level(i-1)+matrixDOF(i-1)
        index2_level(i)=index2_level(i-1)+matrixDOF(i-1)*matrixDOF(i-1)
    enddo
    
    do i=1,nlevel
        count = 0
        index = 2*index_level(i)-1
        do jj=1,numgrid(i)
            do ii=1,numgrid(i)
                coordinate(index+2*count) = ii*gridsize(i)-0.5*gridsize(i);
                coordinate(index+2*count+1) = jj*gridsize(i)-0.5*gridsize(i);
               	count=count+1
            enddo 
        enddo
    enddo

    count = 0

    index = 2*index_level(nlevel)-1
    do jj=1,numgrid(nlevel)
        do ii=1,numgrid(nlevel)
            rhs(index_level(nlevel)+count) = sin(coordinate(index+2*count)/length*PI)* &
                   sin(coordinate(index+2*count+1)/length*PI)*gridsize(nlevel)*gridsize(nlevel)
            count=count+1
        enddo 
    enddo
    print*,'[Poisson] Geometry and rhs constructed.'
    
    do i=1,nlevel
        do ii=1,matrixDOF(i)*matrixDOF(i)
            poissonmatrix(index2_level(i)+ii-1) = 0.0
        enddo
    enddo
    
    do i=1,nlevel
        call construct_poissonmatrix(1, 1, numgrid(i), numgrid(i), numgrid(i), poissonmatrix(index2_level(i)))
    enddo

    print*,'[Poisson] Poisson matrix constructed.'
    
    do i=1,matrixDOF(nlevel)
        solution(index_level(nlevel)+i-1) = 1.0
    enddo
    
    print*,'[Poisson] Start solving equations.'
    
    do i=nlevel,2,-1 
        call restriction(rhs(index_level(i)),rhs(index_level(i-1)),numgrid(i-1))
        print 102, i, i-1
        do ii=1,matrixDOF(i-1)
            solution(index_level(i-1)+ii-1) = 1.0
        enddo
    enddo

    call cgsolver(matrixDOF(1), poissonmatrix(index2_level(1)), rhs(index_level(1)), solution(index_level(1)), &
                  maxiteration, tolerance)

    do i=2,nlevel
        call interpolation(solution(index_level(i-1)),solution(index_level(i)),numgrid(i))
        print 103, i-1, i
        call cgsolver(matrixDOF(i), poissonmatrix(index2_level(i)), rhs(index_level(i)), solution(index_level(i)), &
             maxiteration, tolerance)
        print 104, i

    enddo

    print*,'[Poisson] Final solution obtained.'

    call system("rm -rf result")
    call system("mkdir result")

    do i=1,nlevel
        index = 2*index_level(i)-1
        write(myfilename,100) i
        open(40,file=myfilename,status='unknown', form='formatted')
        count=0
        do ii=1,numgrid(i)
            do jj=1,numgrid(i)
                write(40,101) coordinate(index+2*count),coordinate(index+2*count+1), &
                              solution(index_level(i)+count)
                count=count+1
            enddo
        enddo
        
        close(40)
    enddo

    print*,'[Poisson] Final solution printed.'

    deallocate(rhs)
    deallocate(solution)
    deallocate(coordinate)
    deallocate(poissonmatrix)
    
    deallocate( numgrid )
    deallocate( matrixDOF )
    deallocate( gridsize )
    deallocate( index_level )
    deallocate( index2_level )

    print*,'[Poisson] Memory deallocated.'

    call CPU_TIME(time_end)

    write(*,105) time_end-time_start
    
100 format('./result/solution_level',I1,'.dat')
101 format(3(F12.6,x))
102 format(' [Multigrid] Restriction of RHS done (level ',I2.1,' -> level ',I2.1').')
103 format(' [Multigrid] Interpolation of solution done (level ',I2.1,' -> level ',I2.1').')
104 format(' [Poisson] Solution in level ',I2.1,' obtained.')
105 format(' [Poisson] Mission completed in ',F10.6,' (secs).')

end
