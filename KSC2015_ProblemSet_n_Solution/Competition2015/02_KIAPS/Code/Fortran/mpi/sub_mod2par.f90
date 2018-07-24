module sub_mod2par

use sub_mod1

implicit none

contains

subroutine equal_load(n1,n2,nprocs,myid,istart,ifinish)

implicit none

integer, intent(in) :: nprocs, myid, n1, n2
integer, intent(out) :: istart, ifinish
integer :: iw1, iw2

iw1=(n2-n1+1)/nprocs
iw2=mod(n2-n1+1,nprocs)
istart=myid*iw1+n1+min(myid,iw2)
ifinish=istart+iw1-1 
if(iw2 > myid) ifinish=ifinish+1
if(n2 < istart) ifinish=istart-1

end subroutine equal_load


subroutine interface_flux(qq, fstar, ib, speed, np, ne, nprocs, myrank)

use mpi

implicit none

integer, intent(in) :: np, ne, nprocs, myrank
real(8), intent(in) :: speed
real(8), dimension(np, ne), intent(in) :: qq
real(8), dimension(2, ne), intent(out) :: fstar
real(8), intent(in) :: ib(np, 2)
integer :: ee, el, er, l, r
real(8) :: nl, nr, ql, qr, h
real(8) :: qb(2, ne), qb_start, qb_end
integer :: ierr, status(mpi_status_size), ireq1, ireq2, ireq3, ireq4

fstar = 0D0
l = 1 ! left
r = 2 ! right

do ee = 1, ne
    ! qb(l,ee) is the left edge interpolated value of qq at element ee
    ! qb(r,ee) is the right edge interpolated value of qq at element ee
    qb(l,ee) = dot_product(ib(:,1), qq(:,ee))
    qb(r,ee) = dot_product(ib(:,2), qq(:,ee))
enddo

! calculating numerical flux (fstar) with periodic boundary condition
if (myrank == 0) then
    call mpi_isend(qb(r,ne), 1, mpi_real8, myrank+1, 10, mpi_comm_world, ireq1, ierr) 
    call mpi_irecv(qb_start, 1, mpi_real8, nprocs-1, 10, mpi_comm_world, ireq2, ierr)
else if (myrank == nprocs-1) then
    call mpi_isend(qb(r,ne), 1, mpi_real8, 0,        10, mpi_comm_world, ireq1, ierr)
    call mpi_irecv(qb_start, 1, mpi_real8, myrank-1, 10, mpi_comm_world, ireq2, ierr)
else
    call mpi_isend(qb(r,ne), 1, mpi_real8, myrank+1, 10, mpi_comm_world, ireq1, ierr) 
    call mpi_irecv(qb_start, 1, mpi_real8, myrank-1, 10, mpi_comm_world, ireq2, ierr)
endif

call mpi_wait(ireq1, status, ierr)
call mpi_wait(ireq2, status, ierr)

ee = 1
fstar(l,ee) = -((qb_start+qb(l,ee))/2D0*speed + abs(speed)*(qb_start-qb(l,ee))/2D0)

if (myrank == 0) then
    call mpi_isend(-fstar(l,ee), 1, mpi_real8, nprocs-1, 20, mpi_comm_world, ireq3, ierr)
else
    call mpi_isend(-fstar(l,ee), 1, mpi_real8, myrank-1, 20, mpi_comm_world, ireq3, ierr)
endif

if (myrank == nprocs-1) then
    call mpi_irecv(fstar(r,ne), 1, mpi_real8, 0, 20, mpi_comm_world, ireq4, ierr)
else
    call mpi_irecv(fstar(r,ne), 1, mpi_real8, myrank+1, 20, mpi_comm_world, ireq4, ierr)
endif

do ee = 2, ne
    fstar(l,ee) = -((qb(r,ee-1)+qb(l,ee))/2D0*speed + abs(speed)*(qb(r,ee-1)-qb(l,ee))/2D0)
    fstar(r,ee-1) = -fstar(l,ee)
enddo

call mpi_wait(ireq3, status, ierr)
call mpi_wait(ireq4, status, ierr)

end subroutine interface_flux 


subroutine rhs(qq, rr, dv, df, ib, speed, ne, np, nprocs, myrank)

implicit none

integer, intent(in) :: ne, np, nprocs, myrank
real(8), intent(in) :: qq(np, ne), dv(np, np, ne),df(np, 2, ne)
real(8), intent(in) :: ib(np, 2)
real(8), intent(in) :: speed
real(8), intent(out) :: rr(np, ne)
real(8) :: fstar(2, ne)
integer :: ee, ii, jj

rr = 0D0

do ee = 1, ne
    do ii = 1, np
        do jj = 1, np
            rr(ii,ee) = rr(ii,ee) + speed*dv(ii,jj,ee)*qq(jj,ee)
        enddo
    enddo
enddo

call interface_flux(qq, fstar, ib, speed, np, ne, nprocs, myrank)

do ee = 1, ne
    do ii = 1, np
        do jj = 1, 2
            rr(ii,ee) = rr(ii,ee) - df(ii,jj,ee)*fstar(jj,ee)
        enddo
    enddo
enddo

end subroutine rhs

end module sub_mod2par
