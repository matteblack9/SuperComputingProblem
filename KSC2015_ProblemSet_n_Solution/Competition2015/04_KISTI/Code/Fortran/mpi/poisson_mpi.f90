program main

      use mpi

      implicit none

      real(8), parameter :: PI = 3.1415926535897932384626
      integer(4) :: numgrid_x, numgrid_y, firstgrid_x, lastgrid_x
      integer(4) :: firstrow, lastrow, matrixDOF, globalmatrixDOF
      integer(4) :: count, ii, jj, kk, maxiteration, ierr, myrank, ncpus
      real(8) :: length_x, length_y, gridsize, tolerance
      real(8), allocatable :: poissonmatrix(:), rhs(:), solution(:), coordinate(:)
      real(8) :: elapsed_time
      character(256) ::  myfilename

      call MPI_Init(ierr)
      call MPI_Comm_size(MPI_COMM_WORLD,ncpus,ierr)
      call MPI_Comm_rank(MPI_COMM_WORLD,myrank,ierr)

      if(myrank.EQ.0) then
          write(*,*) '[Poisson] MPI initialized.'
      endif

      elapsed_time = -MPI_Wtime()

      maxiteration = 20000
      tolerance = 1.0e-10
      gridsize = 0.01
      length_x = 1.0
      length_y = 1.0

      numgrid_x = nint(length_x/gridsize) + 1
      numgrid_y = nint(length_y/gridsize) + 1
      globalmatrixDOF = numgrid_x * numgrid_y

      firstgrid_x = (numgrid_x * myrank / ncpus) + 1
      lastgrid_x = numgrid_x * (myrank + 1) / ncpus
      firstrow = (firstgrid_x - 1) * numgrid_y + 1
      lastrow = lastgrid_x * numgrid_y
      matrixDOF =  lastrow - firstrow + 1

      if(myrank.EQ.0) then
          write(*,*) '[Poisson] Geometry and matrix size initialized.'
      endif

      allocate( coordinate(matrixDOF*2) ) 
      allocate( rhs(matrixDOF) )

      count = 1

      do ii = firstgrid_x, lastgrid_x
          do jj = 1, numgrid_y
              coordinate(2*count-1) = (ii-1)*gridsize;
              coordinate(2*count) = (jj-1)*gridsize;
              rhs(count) = sin(coordinate(2*count-1)/length_x*PI)*sin(coordinate(2*count)/length_y*PI)*gridsize*gridsize
              count = count+1
          enddo
      enddo

      if(myrank.EQ.0) then
          write(*,*) '[Poisson] Geometry and rhs constructed.'
      endif

      allocate( poissonmatrix(matrixDOF * globalmatrixDOF) )

      do ii=1, matrixDOF * globalmatrixDOF
          poissonmatrix(ii) = 0.0
      enddo

      call construct_poissonmatrix(firstrow, firstgrid_x, lastgrid_x, numgrid_x, numgrid_y, poissonmatrix)

      if(myrank.EQ.0) then
          write(*,*) '[Poisson] Poisson matrix constructed.'
      endif

      allocate( solution(matrixDOF) )

      do ii=1, matrixDOF
          solution(ii) = 1.0;
      enddo

      if(myrank.EQ.0) then
          write(*,*) '[Poisson] Start solving equations.'
      endif

      call cgsolver(matrixDOF, myrank, ncpus, poissonmatrix, rhs, solution, maxiteration, tolerance)

      if(myrank.EQ.0) then
          write(*,*) '[Poisson] Solution obtained.'
          call system('rm -rf result')
          call system('mkdir result')
      endif

      call MPI_Barrier(MPI_COMM_WORLD,ierr)
      
      do ii = 0, ncpus-1
          if(myrank.EQ.ii) then
              open(11, file='./result/solution.dat', form='formatted',position='append')
              count = 1

              do jj=firstgrid_x, lastgrid_x
                  do kk=1, numgrid_y
                      write(11,101) coordinate(2*count-1), coordinate(2*count), solution(count)
                      count = count + 1
                  enddo
              enddo
              close(11)
          endif
          call MPI_Barrier(MPI_COMM_WORLD,ierr)
      enddo
101   format(2(f9.5,x),f12.6)

      if(myrank.EQ.0) then
          write(*,*) '[Poisson] Solution printed.'
      endif

      deallocate(solution)
      deallocate(poissonmatrix)
      deallocate(rhs)
      deallocate(coordinate)

      if(myrank.EQ.0) then
          write(*,*) '[Poisson] Memory deallocated.'
      endif

      elapsed_time = elapsed_time + MPI_Wtime()

102   format('[Poisson] Finalizing MPI - Computing finished in ',e15.9,' (secs).')
      if(myrank.EQ.0) then
          write(*,102) elapsed_time
      endif

      call MPI_Finalize(ierr)

end program

