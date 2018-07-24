subroutine cgsolver(nsize, matrix, rhs, solution, maxiteration, tolerance)

      implicit none

      integer(4) :: nsize, maxiteration
      real(8) :: matrix(nsize*nsize)
      real(8) :: rhs(nsize)
      real(8) :: solution(nsize)
      real(8) :: tolerance

      integer(4) :: ii, jj, kk
      real(8) :: alpha=0.0, beta=0.0, temp1, temp2, res0tol=0.0
      real(8), allocatable :: res(:), p(:), Ax(:), Ap(:)

      allocate( res(nsize) )
      allocate(   p(nsize) )
      allocate(  Ax(nsize) )
      allocate(  Ap(nsize) )

      call multiply(nsize, matrix, solution, Ax)

      do ii=1, nsize
          res(ii) = rhs(ii) - Ax(ii);
          p(ii)   = res(ii)
      enddo

      call innerproduct(res, res, nsize,res0tol)

      write(*,*) '[CG] Conjugate gradient is started.'
101   format('[CG] mse ',e12.6,' with a tolerance criteria of ',e12.6,' at ',i5,' iterations.')

      do ii=0, maxiteration-1
          if ((mod(ii,20).EQ.0).AND.(ii.NE.0)) then
              write(*,101) sqrt(temp2/res0tol), tolerance, ii
          endif

          call innerproduct(res, res, nsize, temp1)
          call multiply(nsize, matrix, p, Ap)
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
      write(*,102) ii, sqrt(temp2/res0tol)

      deallocate(res)
      deallocate(p)
      deallocate(Ax)
      deallocate(Ap)

end subroutine cgsolver

subroutine innerproduct(x, y, nsize, inprd)

      implicit none

      real(8) :: x(nsize)
      real(8) :: y(nsize)
      integer(4) :: nsize
      real(8) :: inprd

      integer(4) :: ii
      real(8) :: result

      result = 0.0

      do ii=1, nsize
          result = result + x(ii) * y(ii)
      enddo

      inprd = result

end subroutine innerproduct

subroutine multiply(nsize, matrix, x, y)
      implicit none

      integer(4) :: nsize
      real(8) :: matrix(nsize * nsize)
      real(8) :: x(nsize)
      real(8) :: y(nsize)

      integer(4) :: ii, jj

      do ii=1, nsize       ! initialize y
          y(ii) = 0.0
      enddo

      do ii=1, nsize
          do jj=1, nsize
              y(ii) = y(ii) + matrix( (ii-1) * nsize + jj) * x(jj)
          enddo
      enddo

end subroutine multiply
