module m_cgsolver

    implicit none
    
contains

subroutine cgsolver(size, matrix, rhs, solution, maxiteration, tolerance)

    implicit none
    
    integer(4) :: size, maxiteration
    real(8) :: matrix(size*size), rhs(size), solution(size)
    real(8) :: tolerance

    integer(4) :: i,j,k
    real(8) :: alpha=0.0, beta=0.0, temp1, temp2, res0tol=0.0
    real(8), allocatable :: p(:), Ax(:), Ap(:), res(:)

    allocate(res(size))
    allocate(p(size))
    allocate(Ax(size))
    allocate(Ap(size))

    call multiply(size, matrix, solution, Ax)

    do i=1, size
        res(i) = rhs(i) - Ax(i)
        p(i) = res(i)
    enddo 
    
    res0tol = innerproduct(res, res, size)

    print*,'[CG] Conjugate gradient is started.'

    do i=0,maxiteration-1
        if((mod(i,20).EQ.0).AND.(i.NE.0)) then
            print 100, sqrt(temp2/res0tol), tolerance, i
        endif
        
        temp1 = innerproduct(res, res, size)
        call multiply(size, matrix, p, Ap)
        temp2 = innerproduct(Ap, p, size)
        
        alpha = temp1/temp2
        
        do j=1,size
            solution(j) = solution(j) + alpha*p(j)
            res(j) = res(j) - alpha*Ap(j)
        enddo
        
        temp2 = innerproduct(res, res, size)
        
        if( sqrt(temp2/res0tol).LT.tolerance) then
            exit
        endif
        
        beta = temp2/temp1
        
        do j=1,size
            p(j) = res(j) + beta * p(j)
        enddo
    enddo
            

    print 101, i+1, sqrt(temp2/res0tol)
100 format (' [CG] mse ',E15.6,' with a tolerance criteria of ',E15.6,' at ',I5,' iterations.')
101 format (' [CG] Finished with total iteration = ',I5,', mse = ',E15.6)

    deallocate(res)
    deallocate(p)
    deallocate(Ax)
    deallocate(Ap)
    
end subroutine cgsolver

subroutine multiply(size, matrix, x, y)

    implicit none

    integer(4) :: size
    real(8) :: matrix(size*size)
    real(8) :: x(size), y(size)
    integer(4) :: i,j
    
    do i=1, size
        y(i)=0.0
    enddo
    
    do i=1, size
        do j=1, size
            y(i) = y(i) + matrix((i-1)*size+j) * x(j)
        enddo
    enddo

end subroutine multiply

function innerproduct(x, y, size)

    implicit none
    real(8) :: innerproduct
    integer(4) :: size
    real(8) :: x(size), y(size)
    integer(4) :: i
    
    innerproduct = 0.0
    do i=1,size
        innerproduct = innerproduct + x(i) * y(i)
    enddo
        
end function innerproduct
    
end module m_cgsolver


