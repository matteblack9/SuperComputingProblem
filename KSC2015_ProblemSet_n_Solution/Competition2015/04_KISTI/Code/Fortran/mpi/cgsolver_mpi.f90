subroutine cgsolver(nsize, myrank, ncpus, matrix, rhs, solution, maxiteration, tolerance)

      implicit none

      include 'mpif.h'

      integer(4) :: nsize, myrank, ncpus, maxiteration
      real(8) :: matrix(nsize*nsize)
      real(8) :: rhs(nsize)
      real(8) :: solution(nsize)
      real(8) :: tolerance

      integer(4) :: ii, jj, kk, rank(3), maxnsize
      integer(4),allocatable :: load(:), sindex(:)
      real(8) :: alpha=0.0, beta=0.0, temp1, temp2, res0tol=0.0
      real(8), allocatable :: res(:), p(:), Ax(:), Ap(:), xtl(:), xtr(:)
      integer(4) :: ierr

      allocate(   load(0:ncpus-1) )
      allocate( sindex(0:ncpus-1) )
      call MPI_Allgather(nsize, 1, MPI_INTEGER, load(0), 1, MPI_INTEGER, MPI_COMM_WORLD,ierr)

      allocate( res(nsize) )
      allocate(   p(nsize) )
      allocate(  Ax(nsize) )
      allocate(  Ap(nsize) )

      rank(1) = mod( (myrank - 1 + ncpus), ncpus)
      rank(2) = myrank
      rank(3) = mod( (myrank + 1 + ncpus), ncpus)

      maxnsize = 0
      sindex(0) = 0
      do ii = 0, ncpus-1
          if(maxnsize.LT.load(ii)) then
              maxnsize = load(ii)
          endif
          if(ii.NE.0) then
              sindex(ii) = sindex(ii-1)+load(ii-1)
          endif
      enddo

      allocate( xtl(maxnsize) )
      allocate( xtr(maxnsize) )

      call multiply(rank, load, sindex, ncpus, matrix, solution, Ax, xtl, xtr, maxnsize)

      do ii=1, nsize
          res(ii) = rhs(ii) - Ax(ii);
          p(ii)   = res(ii)
      enddo

      call innerproduct(res, res, nsize,res0tol)

      if(myrank.EQ.0) then
          write(*,*) '[CG] Conjugate gradient is started.'
      endif

101   format('[CG] mse ',e12.6,' with a tolerance criteria of ',e12.6,' at ',i5,' iterations.')
      do ii=0, maxiteration-1
          if ( (myrank.EQ.0) .AND. (mod(ii,20).EQ.0) .AND. (ii.NE.0)) then
              write(*,101) sqrt(temp2/res0tol), tolerance, ii
          endif

          call innerproduct(res, res, nsize, temp1)
          call multiply(rank, load, sindex, ncpus, matrix, p, Ap, xtl, xtr, maxnsize)
          call innerproduct(Ap, p, nsize, temp2)

          alpha=temp1/temp2

          do jj=1, nsize
              solution(jj) = solution(jj) + alpha * p(jj)
              res(jj) = res(jj) - alpha * Ap(jj)
          enddo

          call innerproduct(res, res, nsize, temp2);

          if (sqrt(temp2/res0tol) .LT. tolerance) exit 

          beta = temp2/temp1

          do jj=1, nsize
              p(jj)= res(jj) + beta * p(jj)
          enddo
      enddo

102   format('[CG] Finished with total iteration = ',i5,', mse = ',e15.9)
      if(myrank.EQ.0) then
          write(*,102) ii, sqrt(temp2/res0tol)
      endif

      deallocate(xtl)
      deallocate(xtr)
      deallocate(res)
      deallocate(p)
      deallocate(Ax)
      deallocate(Ap)
      deallocate(sindex) 
      deallocate(load) 

end subroutine cgsolver

subroutine innerproduct(x, y, nsize, inprd)

      implicit none

      include 'mpif.h'

      real(8) :: x(nsize)
      real(8) :: y(nsize)
      integer(4) :: nsize
      real(8) :: inprd

      integer(4) :: ii, ierr
      real(8) :: result, globalresult

      result = 0.0

      do ii=1, nsize
          result = result + x(ii) * y(ii)
      enddo

      call MPI_Allreduce(result, globalresult, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

      inprd = globalresult

end subroutine innerproduct

subroutine multiply(rank, load, sindex, ncpus, matrix, x, y, xtl, xtr, maxnsize)

      implicit none

      include 'mpif.h'

      integer(4) :: rank(3), load(0:ncpus-1), sindex(0:ncpus-1), ncpus, maxnsize
      real(8) :: matrix(load(rank(2)) * load(rank(2)))
      real(8) :: x(load(rank(2)))
      real(8) :: y(load(rank(2)))
      real(8) :: xtl(maxnsize)
      real(8) :: xtr(maxnsize)

      integer(4) :: ii, jj, leftrank, myrank, rightrank, globalmatrixDOF
      integer(4) :: req(2), status(mpi_status_size,2), ierr

      leftrank  = rank(1)
      myrank    = rank(2)
      rightrank = rank(3)
      globalmatrixDOF = sindex(ncpus-1) + load(ncpus-1)

      do ii=1, load(myrank)       ! initialize y
          y(ii) = 0.0
      enddo

      if(ncpus.GT.1) then
          call MPI_Irecv(xtl, load(rightrank), MPI_DOUBLE_PRECISION, rightrank, 1002, MPI_COMM_WORLD, req(1), ierr)
          call MPI_Isend(x,   load(myrank),    MPI_DOUBLE_PRECISION,  leftrank, 1002, MPI_COMM_WORLD, req(2), ierr)
      endif

      do ii=1, load(myrank)
          do jj=1, load(myrank)
              y(ii) = y(ii) + matrix( (ii-1) * globalmatrixDOF + sindex(myrank) + jj) * x(jj)
          enddo
      enddo

      if(ncpus.GT.1) then
          call MPI_Waitall(2,req, status, ierr)
      endif

      if(ncpus.GT.2) then
          call MPI_Irecv(xtr, load(leftrank), MPI_DOUBLE_PRECISION,  leftrank, 1002, MPI_COMM_WORLD, req(1), ierr)
          call MPI_Isend(x,   load(myrank),   MPI_DOUBLE_PRECISION, rightrank, 1002, MPI_COMM_WORLD, req(2), ierr)
      endif

      if(ncpus.GT.1) then
          do ii=1, load(myrank)
              do jj=1, load(rightrank)
                  y(ii) = y(ii) + matrix( (ii-1) * globalmatrixDOF + sindex(rightrank) + jj) * xtl(jj)
              enddo
          enddo
      endif

      if(ncpus.GT.2) then
          call MPI_Waitall(2, req, status, ierr)
      endif

      if(ncpus.GT.2) then
          do ii=1, load(myrank)
              do jj=1, load(leftrank)
                  y(ii) = y(ii) + matrix( (ii-1) * globalmatrixDOF + sindex(leftrank) + jj) * xtr(jj)
              enddo
          enddo
      endif

end subroutine multiply
