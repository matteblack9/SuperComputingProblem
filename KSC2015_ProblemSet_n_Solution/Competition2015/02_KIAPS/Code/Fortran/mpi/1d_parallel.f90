program parallel_transport

use sub_mod1
use sub_mod2par
use mpi

implicit none

integer, parameter :: np = 30, tne = 512
! integer, parameter :: np = 20, tne = 16 
! integer, parameter :: np = 4, tne = 5 
real(8), parameter :: tend = 1D1, speed = 1D0
! real(8), parameter :: tend = 1D-1, speed = 1D0
character(*), parameter :: init_type = 'mixed2'
character(10) :: c1, c2
real(8) :: roots(np), weights(np), ll(np), dl(np), xmin, xmax, &
    mf(np,2), ib(np, 2), deltax, jac, xr, xl, cfl, dt, rtime
integer :: err, ii, jj, kk, ee, eres
integer(8) :: nstep
real(8), dimension(np, np) :: smat!, minv

real(8), allocatable :: dx(:), mesh(:), df(:,:,:), xx(:,:), &
    qq(:,:), k1(:,:), k2(:,:), k3(:,:), k4(:,:), mmat(:,:,:), & 
    dv(:,:,:), gxx(:,:), gqq(:,:), buff1(:,:), buff2(:,:), minv_vec(:,:)

integer :: ierr, nprocs, myrank, istart, ifinish, ne, ireq1, ireq2, ireq3, &
    ireq4, status(mpi_status_size)
real(8) :: lxmax, lxmin

! initialize the MPI environmnet
call mpi_init(ierr)
call mpi_comm_size(mpi_comm_world, nprocs, ierr)
call mpi_comm_rank(mpi_comm_world, myrank, ierr)

call equal_load(1, tne, nprocs, myrank, istart, ifinish)

ne = ifinish - istart + 1

allocate(dx(ne), mesh(ne+1), df(np, 2, ne), xx(np, ne), &
    qq(np, ne), k1(np, ne), k2(np, ne), k3(np, ne), k4(np, ne), &
    mmat(np, np, ne), dv(np, np, ne), minv_vec(np, ne))

! setup mesh
! global xmin, xmax
xmin = 0D0
xmax = 1D1
deltax = (xmax-xmin)/real(tne,8)

! local lxmin, lxmax
lxmin = xmin+(istart-1)*deltax
lxmax = xmin+(ifinish)*deltax
mesh(1) = lxmin
mesh(ne+1) = lxmax
do ee = 2, ne
    mesh(ee) = lxmin+(ee-1)*deltax
enddo

! gauss lobatto quadrature point, weight setup
call gausslobatto_quadrature(np, roots, weights, err)

! coordinates and element size
do ee = 1, ne
    xl = mesh(ee)
    xr = mesh(ee+1)
    dx(ee) = xr-xl
    do ii = 1, np
        xx(ii, ee) = xl + 5D-1*(1D0+roots(ii))*dx(ee)
    enddo
enddo

! mass matrix
mmat = 0D0
do ee = 1, ne
    jac = abs(dx(ee))/2D0
    do kk = 1, np
        call lagrange(roots(kk), ll, roots, np)
        do jj = 1, np
            do ii = 1, np
                mmat(ii,jj,ee) = mmat(ii,jj,ee) + jac*weights(kk)*ll(ii)*ll(jj)
            enddo
        enddo
    enddo
enddo

! stiffness matrix
smat = 0D0
do kk = 1, np
    call lagrange(roots(kk), ll, roots, np)
    call lagrange_deriv(roots(kk), dl, roots, np)
    do jj = 1, np
        do ii = 1, np
            smat(ii, jj) = smat(ii, jj) + weights(kk)*ll(jj)*dl(ii)
        enddo
    enddo
enddo

! face integration
mf = 0D0
call lagrange(-1D0, mf(:,1), roots, np)
call lagrange( 1D0, mf(:,2), roots, np)

! boundary interpolation
ib = 0D0
call lagrange(-1D0,ib(:,1), roots, np)
call lagrange( 1D0,ib(:,2), roots, np)

! divergence operators
dv = 0D0
df = 0D0

! do ee = 1, ne
!     call inverse(mmat(:,:,ee), minv, np)
!     dv(:,:,ee) = matmul(minv, smat)
!     df(:,:,ee) = matmul(minv, mf)
! enddo

! mass matrix is diagonal. for the inversion of the mass matrix, we extract 
! the inverse of the diagonal components only
do ee = 1, ne
    do ii = 1, np
        minv_vec(ii,ee) = 1D0/mmat(ii,ii,ee)
        dv(ii,:,ee) = minv_vec(ii,ee)*smat(ii,:)
        df(ii,:,ee) = minv_vec(ii,ee)*mf(ii,:)
    enddo
enddo

call initialize(qq, xx, ne, np, xmax, xmin, init_type)
cfl = 1D0 / np**2
dt = cfl * minval(dx) / abs(speed)
rtime = 0D0
nstep = 0

if (myrank == 0) then
    write(*,20) "Start Time Integration"
end if

do while (rtime < tend)
    dt = min(dt, tend-rtime)
    call rhs(qq,            k1, dv, df, ib, speed, ne, np, nprocs, myrank)
    call rhs(qq+5D-1*dt*k1, k2, dv, df, ib, speed, ne, np, nprocs, myrank)
    call rhs(qq+5D-1*dt*k2, k3, dv, df, ib, speed, ne, np, nprocs, myrank)
    call rhs(qq+     dt*k3, k4, dv, df, ib, speed, ne, np, nprocs, myrank)
    qq = qq + 1D0/6D0*dt*(k1 + 2D0*k2 + 2D0*k3 + k4)
    rtime = rtime + dt
    nstep = nstep + 1

    if (myrank == 0 .and. mod(nstep,10000) == 0) then
        write(*,24) 'nstep = ',nstep, ', ',rtime/tend*1D2, '% complete'
    endif
enddo

if (myrank == 0) then

    if (tne > 200) then
        eres = 2
    elseif (tne > 60) then
        eres = 3
    elseif (tne > 30) then
        eres = 6
    else
        eres = 10
    endif

    write(*,20) "Integration complete"
    write(*,20) '-----------------------------------------------'
    write(*,20) 'code type   : fortran parallel'
    write(*,22) 'Final time  : ', rtime
    write(*,22) 'CFL         : ', cfl
    write(*,21) 'DOF         : ', tne*np
    write(*,21) 'No. of Elem : ', tne
    write(*,21) 'Order       : ', np
    write(*,21) 'eres        : ', eres 
    write(*,21) 'time steps  : ', nstep
    write(*,21) 'nprocs      : ', nprocs
    write(*,20) '-----------------------------------------------'

    allocate(gxx(np,tne), gqq(np,tne), buff1(np,tne/nprocs+1), buff2(np,tne/nprocs+1))
endif


if (myrank /= 0) then
    call mpi_isend(xx, np*ne, mpi_real8, 0, 10, mpi_comm_world, ireq1, ierr)
    call mpi_isend(qq, np*ne, mpi_real8, 0, 20, mpi_comm_world, ireq2, ierr)
    call mpi_isend(istart, 1, mpi_integer, 0, 30, mpi_comm_world, ireq3, ierr)
    call mpi_isend(ifinish, 1, mpi_integer, 0, 40, mpi_comm_world, ireq4, ierr)
    call mpi_wait(ireq3, status, ierr)
    call mpi_wait(ireq4, status, ierr)
    call mpi_wait(ireq1, status, ierr)
    call mpi_wait(ireq2, status, ierr)
endif

if (myrank == 0) then
    do jj=istart, ifinish
        gxx(:,jj) = xx(:,jj)
        gqq(:,jj) = qq(:,jj)
    enddo
    
    do ii=1, nprocs-1
        call mpi_irecv(istart, 1, mpi_integer, ii, 30, mpi_comm_world, ireq3, ierr)
        call mpi_irecv(ifinish, 1, mpi_integer, ii, 40, mpi_comm_world, ireq4, ierr)

        call mpi_wait(ireq3, status, ierr)
        call mpi_wait(ireq4, status, ierr)

        ne = ifinish-istart+1

        call mpi_irecv(buff1, np*ne, mpi_real8, ii, 10, mpi_comm_world, ireq1, ierr)
        call mpi_irecv(buff2, np*ne, mpi_real8, ii, 20, mpi_comm_world, ireq2, ierr)

        call mpi_wait(ireq1, status, ierr)
        call mpi_wait(ireq2, status, ierr)

        kk = 1
        do jj=istart, ifinish
            gxx(:,jj) = buff1(:,kk)
            gqq(:,jj) = buff2(:,kk)
            kk = kk + 1
        enddo
    enddo

    call save_field(gxx, gqq, tne, np, roots, eres)

    deallocate(gxx, gqq, buff1, buff2)
endif

deallocate(dx, mesh, df, xx, qq, k1, k2, k3, k4, mmat, dv, minv_vec)

20 format (A)
21 format (A, I13)
22 format (A, E13.5)
23 format (F5.1, A)
24 format (A, I10, A,F5.1, A)
25 format (A, F13.3)

call mpi_finalize(ierr)

end program


