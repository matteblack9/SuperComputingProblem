program main

      implicit none

      real(8), parameter :: PI = 3.1415926535897932384626
      integer(4) :: numgrid_x, numgrid_y
      integer(4) :: matrixDOF
      integer(4) :: count, ii, jj, kk, maxiteration
      real(8) :: length_x, length_y, gridsize, tolerance
      real(8), allocatable :: poissonmatrix(:), rhs(:), solution(:), coordinate(:)
      character(256) ::  myfilename
!FILE *fp;

      maxiteration = 20000
      tolerance = 1.0e-10
      gridsize = 0.01
      length_x = 1.0
      length_y = 1.0

      numgrid_x = nint(length_x/gridsize) + 1
      numgrid_y = nint(length_y/gridsize) + 1
      matrixDOF = numgrid_x * numgrid_y

      write(*,*) '[Poisson] Geometry and matrix size initialized.'

      allocate( coordinate(matrixDOF*2) ) 
      allocate( rhs(matrixDOF) )

      count = 1

      do ii = 1, numgrid_x
          do jj = 1, numgrid_y
              coordinate(2*count-1) = (ii-1)*gridsize;
              coordinate(2*count) = (jj-1)*gridsize;
              rhs(count) = sin(coordinate(2*count-1)/length_x*PI)*sin(coordinate(2*count)/length_y*PI)*gridsize*gridsize
              count = count+1
          enddo
      enddo

      write(*,*) '[Poisson] Geometry and rhs constructed.'

      allocate( poissonmatrix(matrixDOF*matrixDOF) )

      do ii=1, matrixDOF*matrixDOF
          poissonmatrix(ii) = 0.0
      enddo

      call construct_poissonmatrix(1, 1, numgrid_x, numgrid_x, numgrid_y, poissonmatrix)

      write(*,*) '[Poisson] Poisson matrix constructed.'

      allocate( solution(matrixDOF) )

      do ii=1, matrixDOF
          solution(ii) = 1.0;
      enddo

      write(*,*) '[Poisson] Start solving equations.'

      call cgsolver(matrixDOF, poissonmatrix, rhs, solution, maxiteration, tolerance)

      write(*,*) '[Poisson] Solution obtained.'

      call system('rm -rf result')
      call system('mkdir result')
      open(11, file='./result/solution.dat', form='formatted')
      
      count = 1
      do ii=1, numgrid_x
          do jj=1, numgrid_y
              write(11,101) coordinate(2*count-1), coordinate(2*count), solution(count)
              count = count + 1
          enddo
      enddo
101   format(2(f9.5,x),f12.6)

      close(11)

      write(*,*) '[Poisson] Solution printed.'

      deallocate(solution)
      deallocate(poissonmatrix)
      deallocate(rhs)
      deallocate(coordinate)

      write(*,*) '[Poisson] Memory deallocated.'

end program

