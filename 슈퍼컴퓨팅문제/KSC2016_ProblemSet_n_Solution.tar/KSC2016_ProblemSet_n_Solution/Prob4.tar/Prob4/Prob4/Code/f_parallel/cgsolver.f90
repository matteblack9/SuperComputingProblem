module m_cgsolver

    use mpi
    implicit none
    
contains

subroutine cgsolver(size, globalsize, myrank, ncpus, matrix, rhs, solution, maxiteration, tolerance)

    implicit none
    
    integer(4) :: size, globalsize, myrank, ncpus, maxiteration
    real(8) :: matrix(size*globalsize), rhs(size), solution(size)
    real(8) :: tolerance

    integer(4) :: i,j,k,rank(3),ierr
    integer(4) :: load(ncpus), sindex(ncpus), maxsize
    real(8) :: alpha=0.0, beta=0.0, temp1, temp2, res0tol=0.0
    real(8), allocatable :: p(:), Ax(:), Ap(:), res(:), xtp(:), xtn(:)

    call MPI_ALLGATHER(size,1,MPI_INT,load,1,MPI_INT,MPI_COMM_WORLD,ierr)

    allocate(res(size))
    allocate(p(size))
    allocate(Ax(size))
    allocate(Ap(size))

    rank(1) = mod(myrank-1+ncpus,ncpus)
    rank(2) = myrank
    rank(3) = mod(myrank+1,ncpus)

    maxsize = 0
    sindex(1) = 0

    do i = 1,ncpus
        if(maxsize.lt.load(i)) then
            maxsize = load(i)
        endif
        if(i.ne.1) then
            sindex(i) = sindex(i-1)+load(i-1)
        endif
    enddo

    allocate(xtp(maxsize))
    allocate(xtn(maxsize))

    call multiply(size,globalsize,maxsize,rank,load,sindex,ncpus, matrix, solution, Ax,xtp,xtn)

    do i=1, size
        res(i) = rhs(i) - Ax(i)
        p(i) = res(i)
    enddo 
    
    res0tol = innerproduct(res, res, size)

    if(myrank.eq.0) print*,'[CG] Conjugate gradient is started.'

    do i=0,maxiteration-1
        if((mod(i,20).EQ.0).AND.(i.NE.0)) then
            if(myrank.eq.0) print 100, sqrt(temp2/res0tol), tolerance, i
        endif

        temp1 = innerproduct(res, res, size)
        call multiply(size, globalsize, maxsize, rank,load,sindex,ncpus, matrix, p, Ap, xtp,xtn)
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
            

    if(myrank.eq.0) print 101, i+1, sqrt(temp2/res0tol)
100 format (' [CG] mse ',E15.6,' with a tolerance criteria of ',E15.6,' at ',I5,' iterations.')
101 format (' [CG] Finished with total iteration = ',I5,', mse = ',E15.6)

    deallocate(res)
    deallocate(p)
    deallocate(Ax)
    deallocate(Ap)
    deallocate(xtn)
    deallocate(xtp)
    
end subroutine cgsolver

subroutine multiply(size,globalsize,maxsize,rank,load,sindex,ncpus, matrix, x, y,xtp,xtn)

    use mpi
    implicit none

    integer(4) :: size,globalsize,maxsize,rank(3),load(ncpus),sindex(ncpus),ncpus
    real(8) :: matrix(size*globalsize)
    real(8) :: x(size), y(size), xtp(maxsize), xtn(maxsize)
    integer(4) :: i,j, prevrank,nextrank,myrank,globalmatrixDOF
    integer(4) :: req(2), status(MPI_STATUS_SIZE), ierr
    
    prevrank = rank(1)
    myrank = rank(2)
    nextrank = rank(3)
    globalmatrixDOF = sindex(ncpus) + load(ncpus)

    do i=1, load(myrank+1)
        y(i)=0.0
    enddo

    call MPI_Irecv(xtp, load(nextrank+1), MPI_REAL8, nextrank, 1002, MPI_COMM_WORLD, req(1),ierr)
    call MPI_Isend(x, load(myrank+1), MPI_REAL8, prevrank, 1002, MPI_COMM_WORLD, req(2),ierr)
    
    do j=1, load(myrank+1)
        do i=1, load(myrank+1)
            y(j) = y(j) + matrix((j-1)*globalsize+sindex(myrank+1)+i) * x(i)
        enddo
    enddo
    call MPI_WAITALL(2,req,status,ierr)

    call MPI_Irecv(xtn, load(prevrank+1), MPI_REAL8, prevrank, 1002, MPI_COMM_WORLD, req(1),ierr)
    call MPI_Isend(x, load(myrank+1), MPI_REAL8, nextrank, 1002, MPI_COMM_WORLD, req(2),ierr)

    do j=1, load(myrank+1)
        do i=1, load(nextrank+1)
            y(j) = y(j) + matrix((j-1)*globalsize+sindex(nextrank+1)+i) * xtp(i)
        enddo
    enddo
    call MPI_WAITALL(2,req,status,ierr)

    do j=1, load(myrank+1)
        do i=1, load(prevrank+1)
            y(j) = y(j) + matrix((j-1)*globalsize+sindex(prevrank+1)+i) * xtn(i)
        enddo
    enddo

end subroutine multiply

function innerproduct(x, y, size)

    use mpi
    implicit none
    real(8) :: innerproduct, innerproduct_local
    integer(4) :: size
    real(8) :: x(size), y(size)
    integer(4) :: i,ierr
    
    innerproduct_local = 0.0
    innerproduct = 0.0
    do i=1,size
        innerproduct_local = innerproduct_local + x(i) * y(i)
    enddo
    call MPI_ALLREDUCE(innerproduct_local,innerproduct,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
        
end function innerproduct
    
end module m_cgsolver
