program serial_transport

use sub_mod1
use sub_mod2

implicit none

! integer, parameter :: np = 30, ne = 512
integer, parameter :: np = 5, ne = 10
                                       ! np: number of nodes per element
                                       ! ne: number of elements
                                       ! caution: np 20, ne 256 will take approx 10 minutes on a serial job
! real(8), parameter :: tend = 1D4, speed = 1D0
real(8), parameter :: tend = 1D2, speed = 1D0
! real(8), parameter :: tend = 1D-1, speed = 1D-0
character(*), parameter :: init_type = 'mixed2'
character(10) :: c1, c2
real(8) :: roots(np), weights(np), ll(np), dl(np), xmin, xmax, &
    mf(np,2), ib(np, 2), deltax, jac, xr, xl, cfl, dt, rtime
integer :: err, ii, jj, kk, ee, eres
integer(8) :: nstep
real(8), dimension(np, np) :: smat, minv

real(8) :: dx(ne), mesh(ne+1), df(np, 2, ne), fstar(2, ne)
real(8), dimension(np, ne) :: xx, qq, k1, k2, k3, k4, minv_vec
real(8), dimension(np, np, ne) :: mmat, dv

! setup mesh
xmin = 0D0
xmax = 1D1
deltax = (xmax-xmin)/real(ne,8)
mesh(1) = xmin
mesh(ne+1) = xmax
do ee = 2, ne
    mesh(ee) = xmin+(ee-1)*deltax
enddo

! gauss lobatto quadrature point, weight setup
call gausslobatto_quadrature(np, roots, weights, err)

! coordinates and element size
do ee = 1, ne
    xl = mesh(ee)
    xr = mesh(ee+1)
    dx(ee) = xr-xl ! size of each element
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

! initialize qq field
call initialize(qq, xx, ne, np, xmax, xmin, init_type)
cfl = 1D0 / np**2
dt = cfl * minval(dx) / abs(speed)
rtime = 0D0
nstep = 0

write(*,20) "Start Time Integration"

! Runge-Kutta 4th order Time Integration
do while (rtime < tend)
    dt = min(dt, tend-rtime)
    call rhs(qq,            k1, dv, df, ib, speed, ne, np)
    call rhs(qq+5D-1*dt*k1, k2, dv, df, ib, speed, ne, np)
    call rhs(qq+5D-1*dt*k2, k3, dv, df, ib, speed, ne, np)
    call rhs(qq+     dt*k3, k4, dv, df, ib, speed, ne, np)
    qq = qq + 1D0/6D0*dt*(k1 + 2D0*k2 + 2D0*k3 + k4)

    rtime = rtime + dt
    nstep = nstep + 1
    if (mod(nstep,10000) == 0) then
        write(*,24) 'nstep = ',nstep, ', ',rtime/tend*1D2, '% complete'
    endif
enddo

write(*,20) "Integration complete"

if (ne > 200) then
    eres = 2
elseif (ne > 60) then
    eres = 3
elseif (ne > 30) then
    eres = 6
else
    eres = 10
endif

! final report
write(*,20) '-----------------------------------------------'
write(*,20) 'code type   : fortran serial'
write(*,22) 'Final time  : ', rtime
write(*,22) 'CFL         : ', cfl
write(*,21) 'DOF         : ', ne*np
write(*,21) 'No. of Elem : ', ne
write(*,21) 'Order       : ', np
write(*,21) 'eres        : ', eres 
write(*,21) 'time steps  : ', nstep
write(*,20) '-----------------------------------------------'

call save_field(xx, qq, ne, np, roots, eres)

20 format (A)
21 format (A, I13)
22 format (A, E13.5)
23 format (F5.1, A)
24 format (A, I10, A,F5.1, A)
25 format (A, F13.3)

end program

