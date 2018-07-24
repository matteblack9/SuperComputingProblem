program poisson
    
    use multigrid
    use m_cgsolver
    use matrixconstructor
    use mpi
    
    implicit none
    
    real(8), parameter :: PI = 3.141592653589793
    integer(4) :: nlevel, initial_ngrid
    integer(4), allocatable :: numgrid(:), globalnumgrid(:), matrixDOF(:), globalmatrixDOF(:)
    integer(4), allocatable :: index_level(:), index2_level(:)
    integer(4), allocatable :: globalindex_level(:), globalindex2_level(:)
    integer(4), allocatable :: firstgrid_y(:), lastgrid_y(:), firstrow(:), lastrow(:)
    integer(4) :: sum_matrixDOF, sum2_matrixDOF, index
    integer(4) :: globalsum_matrixDOF, globalsum2_matrixDOF, globalindex
    integer(4) :: count, ii, jj, kk, ll, i, maxiteration, myrank, ncpus, rank
    integer(4) :: ierr
    real(8)    ::length, tolerance
    real(8), allocatable :: gridsize(:), poissonmatrix(:),rhs(:), solution(:), coordinate(:)
    real(8) :: time_start, time_end
   
    character(LEN=256) :: myfilename

    nlevel = 3
    initial_ngrid = 128
    maxiteration = 20000
    tolerance = 1.0e-10
    length = 1.0

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,ncpus,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

    time_start = MPI_WTIME()

    allocate( numgrid(nlevel) )
    allocate( globalnumgrid(nlevel) )
    allocate( matrixDOF(nlevel) )
    allocate( globalmatrixDOF(nlevel) )
    allocate( firstgrid_y(nlevel) )
    allocate( lastgrid_y(nlevel) )
    allocate( firstrow(nlevel) )
    allocate( lastrow(nlevel) )
    allocate( gridsize(nlevel) )
    allocate( index_level(nlevel) )
    allocate( index2_level(nlevel) )
    allocate( globalindex_level(nlevel) )
    allocate( globalindex2_level(nlevel) )

    globalnumgrid(nlevel) = initial_ngrid
    gridsize(nlevel) = length/initial_ngrid
    globalmatrixDOF(nlevel) = globalnumgrid(nlevel)*globalnumgrid(nlevel)

    firstgrid_y(nlevel) = (globalnumgrid(nlevel)*myrank/ncpus)+1
    lastgrid_y(nlevel) = globalnumgrid(nlevel)*(myrank+1)/ncpus
    numgrid(nlevel) = lastgrid_y(nlevel)-firstgrid_y(nlevel)+1
    matrixDOF(nlevel) = numgrid(nlevel)*globalnumgrid(nlevel)
    firstrow(nlevel) = (firstgrid_y(nlevel)-1)*globalnumgrid(nlevel)+1
    lastrow(nlevel) = lastgrid_y(nlevel)*globalnumgrid(nlevel)

    sum_matrixDOF = matrixDOF(nlevel)
    sum2_matrixDOF = matrixDOF(nlevel)*globalmatrixDOF(nlevel)
    globalsum_matrixDOF = globalmatrixDOF(nlevel)
    globalsum2_matrixDOF = globalmatrixDOF(nlevel)*globalmatrixDOF(nlevel)

    do i=nlevel,2,-1
        globalnumgrid(i-1) = globalnumgrid(i)/2
        gridsize(i-1) = length/globalnumgrid(i-1)
        globalmatrixDOF(i-1) = globalnumgrid(i-1)*globalnumgrid(i-1)
        globalsum_matrixDOF = globalsum_matrixDOF + globalmatrixDOF(i-1)
        globalsum2_matrixDOF = globalsum2_matrixDOF + globalmatrixDOF(i-1)*globalmatrixDOF(i-1)

        firstgrid_y(i-1) = (globalnumgrid(i-1)*myrank/ncpus)+1
        lastgrid_y(i-1) = globalnumgrid(i-1)*(myrank+1)/ncpus
        numgrid(i-1) = lastgrid_y(i-1)-firstgrid_y(i-1)+1
        matrixDOF(i-1) = numgrid(i-1)*globalnumgrid(i-1)
        firstrow(i-1) = (firstgrid_y(i-1)-1)*globalnumgrid(i-1)+1
        lastrow(i-1) = lastgrid_y(i-1)*globalnumgrid(i-1)
        sum_matrixDOF = sum_matrixDOF + matrixDOF(i-1)
        sum2_matrixDOF = sum2_matrixDOF + globalmatrixDOF(i-1)*matrixDOF(i-1)
    enddo

    if((globalnumgrid(1)<ncpus).or.((nlevel.gt.1).and.(mod(globalnumgrid(1),ncpus).ne.0))) then
        if((globalnumgrid(1).lt.ncpus).and.(myrank.eq.0)) then
            print*,'[Poisson] # of MPI processes must be equal or larger than # of grids in the (smallest) domain!'
        endif
        if((nlevel.gt.1).and.(mod(globalnumgrid(1),ncpus).ne.0).and.myrank.eq.0) then
            print*,'[Poisson/Multigrid] # of grids in the smallest domain must be a multiple of # of MPI processes!'
        endif

            deallocate(gridsize)
            deallocate(globalnumgrid)
            deallocate(globalmatrixDOF)
            deallocate(matrixDOF)
            deallocate(numgrid)
            deallocate(firstrow)
            deallocate(lastrow)
            deallocate(firstgrid_y)
            deallocate(lastgrid_y)
            deallocate( index_level )
            deallocate( index2_level )
            deallocate( globalindex_level )
            deallocate( globalindex2_level )

            call MPI_Finalize(ierr)

           if(myrank.eq.0) print*,'[Poisson] MPI finalized.'
           stop

    endif

    
    if(myrank.eq.0) print*,'[Poisson] Geometry and matrix size initialized.'

    allocate( rhs(sum_matrixDOF) )
    allocate( solution(sum_matrixDOF) )
    allocate( coordinate(sum_matrixDOF*2) )
    allocate( poissonmatrix(sum2_matrixDOF) )
    
    index_level(1)=1
    index2_level(1)=1
    do i=2,nlevel
        index_level(i)=index_level(i-1)+matrixDOF(i-1)
        index2_level(i)=index2_level(i-1)+globalmatrixDOF(i-1)*matrixDOF(i-1)
    enddo
    
    do i=1,nlevel
        count = 0
        index = 2*index_level(i)-1
        do jj=firstgrid_y(i),lastgrid_y(i)
            do ii=1,globalnumgrid(i)
                coordinate(index+2*count) = ii*gridsize(i)-0.5*gridsize(i);
                coordinate(index+2*count+1) = jj*gridsize(i)-0.5*gridsize(i);
               	count=count+1
            enddo 
        enddo
!        write(*,*) i, myrank, numgrid(i), coordinate(index+1), coordinate(index+2*count-2+1), &
!                   firstrow(i), lastrow(i), index_level(i)
    enddo

    count = 0

    index = 2*index_level(nlevel)-1
    do jj=1,numgrid(nlevel)
        do ii=1,globalnumgrid(nlevel)
            rhs(index_level(nlevel)+count) = sin(coordinate(index+2*count)/length*PI)* &
                   sin(coordinate(index+2*count+1)/length*PI)*gridsize(nlevel)*gridsize(nlevel)
            count=count+1
        enddo 
    enddo
    if(myrank.eq.0) print*,'[Poisson] Geometry and rhs constructed.'
    
    do i=1,nlevel
        do ii=1,matrixDOF(i)*globalmatrixDOF(i)
            poissonmatrix(index2_level(i)+ii-1) = 0.0
        enddo
    enddo
    
    do i=1,nlevel
        call construct_poissonmatrix(firstrow(i), firstgrid_y(i), lastgrid_y(i), &
                                     globalnumgrid(i), globalnumgrid(i), poissonmatrix(index2_level(i)))
    enddo

!    do i=1,matrixDOF(index_level(nlevel))*globalmatrixDOF(index_level(nlevel))
!        write(*,*) myrank, i+(firstgrid_y(nlevel)-1)*matrixDOF(nlevel)*globalmatrixDOF(nlevel) ,&
!                   poissonmatrix(index2_level(nlevel)+i-1)
!    enddo
    
    if(myrank.eq.0) print*,'[Poisson] Poisson matrix constructed.'
    
    do i=1,matrixDOF(nlevel)
        solution(index_level(nlevel)+i-1) = 1.0
    enddo
    
    if(myrank.eq.0) print*,'[Poisson] Start solving equations.'
    
    do i=nlevel,2,-1 
        call restriction(rhs(index_level(i)),rhs(index_level(i-1)),globalnumgrid(i-1),numgrid(i-1))
        if(myrank.eq.0) print 102, i, i-1
        do ii=1,matrixDOF(i-1)
            solution(index_level(i-1)+ii-1) = 1.0
        enddo
    enddo

    call cgsolver(matrixDOF(1), globalmatrixDOF(1), myrank, ncpus, poissonmatrix(index2_level(1)), rhs(index_level(1)), &
                  solution(index_level(1)), maxiteration, tolerance)

    do i=2,nlevel
        call interpolation(solution(index_level(i-1)),solution(index_level(i)),globalnumgrid(i),numgrid(i),myrank,ncpus)
        if(myrank.eq.0) print 103, i-1, i
        call cgsolver(matrixDOF(i), globalmatrixDOF(i), myrank, ncpus, poissonmatrix(index2_level(i)), rhs(index_level(i)), &
                      solution(index_level(i)), maxiteration, tolerance)
        if(myrank.eq.0) print 104, i

    enddo

    if(myrank.eq.0) print*,'[Poisson] Final solution obtained.'

    if(myrank.eq.0) call system("rm -rf result")
    if(myrank.eq.0) call system("mkdir result")

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    do i=1,nlevel
        index = 2*index_level(i)-1
        write(myfilename,100) i

        do rank=0,ncpus-1
            if(myrank.eq.rank) then
                if(myrank.eq.0) then
                    open(40,file=myfilename,status="new", form="formatted", action="write")
                else
                    open(40,file=myfilename,status="old", position="append", form="formatted", action="write")
                endif
                    count=0
                    do jj=1,numgrid(i)
                        do ii=1,globalnumgrid(i)
                            write(40,101) coordinate(index+2*count),coordinate(index+2*count+1), &
                                              solution(index_level(i)+count)
                            count=count+1
                        enddo
                    enddo
                close(40)
            endif
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        enddo
    enddo

    if(myrank.eq.0) print*,'[Poisson] Final solution printed.'

    deallocate(rhs)
    deallocate(solution)
    deallocate(coordinate)
    deallocate(poissonmatrix)
    
    deallocate( numgrid )
    deallocate( globalnumgrid )
    deallocate( matrixDOF )
    deallocate( globalmatrixDOF )
    deallocate( gridsize )
    deallocate( index_level )
    deallocate( index2_level )
    deallocate( firstgrid_y )
    deallocate( lastgrid_y )
    deallocate( firstrow )
    deallocate( lastrow )
    deallocate( globalindex_level )
    deallocate( globalindex2_level )

    if(myrank.eq.0) print*,'[Poisson] Memory deallocated.'

    time_end=MPI_WTIME()

    if(myrank.eq.0) write(*,105) time_end-time_start
    call MPI_Finalize(ierr)
    if(myrank.eq.0) print*,'[Poisson] MPI finalized.'
    
100 format('./result/solution_level',I1,'.dat')
101 format(3(F12.6,x))
102 format(' [Multigrid] Restriction of RHS done (level ',I2.1,' -> level ',I2.1').')
103 format(' [Multigrid] Interpolation of solution done (level ',I2.1,' -> level ',I2.1').')
104 format(' [Poisson] Solution in level ',I2.1,' obtained.')
105 format(' [Poisson] Mission completed in ',F10.6,' (secs).')

end
