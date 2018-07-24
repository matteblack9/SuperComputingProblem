module sub_mod1

implicit none

contains

subroutine lagrange(xx, ll, pts, lsize)

implicit none
integer, intent(in) :: lsize
real(8), intent(in) :: xx, pts(lsize)
real(8), intent(out) :: ll(lsize)
integer :: ii, jj

do ii=1, lsize
    ll(ii) = 1D0
    do jj=1, lsize
        if (jj .ne. ii) then
            ll(ii) = ll(ii)*(xx-pts(jj))/(pts(ii)-pts(jj))
        endif
    enddo
enddo

end subroutine lagrange

subroutine lagrange_deriv(xx, dl, pts, dlsize)

implicit none
integer, intent(in) :: dlsize
real(8), intent(in) :: xx, pts(dlsize)
real(8), intent(out) :: dl(dlsize)
integer :: ii, jj, kk
real(8) :: term

do ii=1, dlsize
    dl(ii) = 0D0
    do jj=1, dlsize
        if (jj .ne. ii) then
            term = 1D0 / (pts(ii)-pts(jj))
            do kk=1, dlsize
                if ((kk .ne. ii) .and. (kk .ne. jj)) then
                    term = term*(xx-pts(kk))/(pts(ii)-pts(kk))
                endif
            enddo
            dl(ii) = dl(ii) + term
        endif
    enddo
enddo
end subroutine lagrange_deriv

recursive subroutine legendre(n, x, xsize, output)

implicit none
integer, intent(in) :: n, xsize
real(8), intent(in) :: x(xsize)
real(8), intent(out) :: output(xsize)
real(8) :: temp1(xsize), temp2(xsize)

if (n .eq. 0) then
    output = x*0D0+1.D0
else if (n .eq. 1) then
    output = x
else
    call legendre(n-1, x, xsize, temp1)
    call legendre(n-2, x, xsize, temp2)
    output = ((2D0*n-1D0)*x*temp1-(n-1)*temp2)/n
endif

end subroutine legendre

recursive subroutine legendre_scalar(n, x, output)

implicit none
integer, intent(in) :: n
real(8), intent(in) :: x
real(8), intent(out) :: output
real(8) :: temp1, temp2

if (n .eq. 0) then
    output = x*0D0+1.D0
else if (n .eq. 1) then
    output = x
else
    call legendre_scalar(n-1, x, temp1)
    call legendre_scalar(n-2, x, temp2)
    output = ((2D0*n-1D0)*x*temp1-(n-1)*temp2)/n
endif

end subroutine legendre_scalar

subroutine dlegendre_scalar(n, x, output)

implicit none
integer, intent(in) :: n
real(8), intent(in) :: x
real(8), intent(out) :: output
real(8) :: temp1, temp2

call legendre_scalar(n-1,x, temp1)
call legendre_scalar(n, x, temp2)

output = n*(temp1-x*temp2)/(1D0-x**2)

end subroutine dlegendre_scalar

subroutine ddlegendre_scalar(n, x, output)

implicit none
integer, intent(in) :: n
real(8), intent(in) :: x
real(8), intent(out) :: output
real(8) :: temp1, temp2

call dlegendre_scalar(n,x, temp1)
call legendre_scalar(n, x, temp2)

output = (2D0*x*temp1-n*(n+1)*temp2)/(1D0-x**2)

end subroutine ddlegendre_scalar

subroutine dddlegendre_scalar(n, x, output)

implicit none
integer, intent(in) :: n
real(8), intent(in) :: x
real(8), intent(out) :: output
real(8) :: temp1, temp2

call ddlegendre_scalar(n,x, temp1)
call dlegendre_scalar(n, x, temp2)

output = (2D0*x*temp1-(n*(n+1)-2D0)*temp2)/(1D0-x**2)

end subroutine dddlegendre_scalar

subroutine dlegendre_roots(polyorder, output, err, tol)

implicit none
integer, intent(in) :: polyorder
real(8), intent(out) :: output(polyorder)
integer, intent(out) :: err
real(8), intent(in), optional :: tol
real(8) :: tolerance, x, pi, error, dx, y, dy, ddy
real(8), allocatable :: temproot(:)
integer :: ii, iters, n

n = polyorder

pi = acos(0D0)*2D0

if (present(tol)) then
    tolerance = tol
else
    tolerance = 1D-16
endif

output = 1D0

if (polyorder < 2) then
    err = 1
else
    allocate(temproot(1:polyorder/2))
    do ii=1,polyorder/2-1
        x = (1D0 - (3D0*(n-2D0)/(8D0*(n-1D0)**3)))*cos((4D0*(ii+1D0)-3D0)/(4D0*(n-1D0)+1D0)*pi)
        error = 10*tolerance
        iters = 0
        do while ((error > tolerance) .and. (iters < 1000))
            call dlegendre_scalar(polyorder-1, x, y)
            call ddlegendre_scalar(polyorder-1, x, dy)
            call dddlegendre_scalar(polyorder-1, x, ddy)
            dx =  -2D0*y*dy/(2D0*dy**2-y*ddy)
            x = x + dx
            iters = iters + 1
            error = abs(dx)
        enddo
        output(ii+1) = x
    enddo
    temproot(:) = output(1:polyorder/2)
    if (mod(polyorder, 2) .eq. 0) then
        output(:polyorder/2) = -temproot(:)
        output(polyorder/2+1:) = temproot(polyorder/2:1:-1)
    else
        output(:polyorder/2) = -temproot(:)
        output(polyorder/2+1) = 0D0
        output(polyorder/2+2:) = temproot(polyorder/2:1:-1)
    endif
    err = 0
    deallocate(temproot)
endif

end subroutine dlegendre_roots

subroutine gausslobatto_quadrature(polyorder, roots, weights, err)

implicit none
integer, intent(in) :: polyorder
real(8), intent(out) :: roots(polyorder), weights(polyorder)
integer, intent(out) :: err
integer :: n
real(8) :: temp(polyorder)

if (polyorder .eq. 1) then
    roots = 0D0
    weights = 2D0
    err = 0
else
    n = polyorder
    ! call dlegendre_roots(n, roots, err, tol)
    call dlegendre_roots(n, roots, err)
    if (err .eq. 0) then
        call legendre(n-1, roots, n, temp)
        weights = 2D0/(n*(n-1D0)*(temp**2))
        err=0
    else
        err=1
    endif
endif

end subroutine gausslobatto_quadrature

subroutine inverse(inmat,outmat,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! inmat(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! outmat(n,n) - inverse matrix of A
! source from
! http://ww2.odu.edu/~agodunov/computing/programs/book2/Ch06/Inverse.f90
implicit none 
integer, intent(in) :: n
real(8), dimension(n,n), intent(in) :: inmat
real(8), dimension(n,n), intent(out) :: outmat 
real(8), dimension(n,n) :: a, L, U
real(8), dimension(n) :: b, d, x
real(8) :: coeff
integer :: i, j, k

a = inmat
L=0.0
U=0.0
b=0.0

do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

do i=1,n
  L(i,i) = 1.0
end do
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

do k=1,n
  b(k)=1.0
  d(1) = b(1)
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
  do i=1,n
    outmat(i,k) = x(i)
  end do
  b(k)=0.0
end do
end subroutine inverse

subroutine save_field(xx, qq, ne, np, roots, eres)

implicit none

integer, intent(in) :: ne, np, eres
real(8), intent(in) :: qq(np, ne), xx(np, ne), roots(np)
integer :: ii, ee 
real(8) :: plot_coords(eres), ll(np), xcoord, ycoord, dx


dx = 2D0/real(eres-1,8)
plot_coords(1) = -1D0
plot_coords(eres) = 1D0
do ii = 2, eres-1
    plot_coords(ii) = -1D0 + real(ii-1,8)*dx
enddo

! save as ascii file
open(unit=11, file='data.txt',status='replace')
do ee = 1, ne
    do ii = 1, eres
        call lagrange(plot_coords(ii), ll, roots, np)
        xcoord = dot_product(ll, xx(:,ee))
        ycoord = dot_product(ll, qq(:,ee))
        write(11, '(F21.6, 3X, F21.6)') xcoord, ycoord
    enddo
enddo

close(11)

! ! save as binary files
! open(unit=11, file='xcoord.bin', access='direct', form='unformatted', &
!   recl=eres*ne*8, status='replace')
! open(unit=12, file='ycoord.bin', access='direct', form='unformatted', &
!   recl=eres*ne*8, status='replace')
! write(11,rec=1) xcoord
! write(12,rec=1) ycoord
!
! close(11)
! close(12)

end subroutine save_field

subroutine initialize(qq, xx, ne, np, xmax, xmin, init_type)

implicit none

integer, intent(in) :: ne, np
character(*) :: init_type
real(8), intent(in) :: xmax, xmin, xx(np, ne)
real(8), intent(out) :: qq(np, ne)
integer :: ii, ee
real(8) :: pi, mu, sig

pi = 4D0 * atan(1D0)

if (init_type == 'sin') then
    qq = sin(2D0 * pi * xx/(xmax-xmin))
elseif (init_type == 'gaussian') then
    mu = 5D-1*(xmin+xmax)
    sig = 1D-1*(xmax-xmin)
    qq = exp( -(xx-mu)**2 / (2D0*sig**2))
elseif (init_type == 'box') then
    mu = 5D-1*(xmin+xmax)
    sig = 1D-1*(xmax-xmin)
    qq = 0D0
    do ee= 1, ne
        do ii= 1, np
            if (abs(mu-xx(ii,ee))<sig) then
                qq(ii,ee) = 1D0
            endif
        enddo
    enddo
elseif (init_type == 'mixed1') then
    qq = 2D-1*sin(2D0 * pi * xx/(xmax-xmin))
    mu = 4D-1*(xmin+xmax)
    sig = 5D-2*(xmax-xmin)
    qq = qq +2D0*exp( -(xx-mu)**2 / (2D0*sig**2))
elseif (init_type == 'mixed2') then
    mu = 3D-1*(xmin+xmax)
    sig = 1D-1*(xmax-xmin)
    qq = 0D0
    do ee= 1, ne
        do ii= 1, np
            if (abs(mu-xx(ii,ee))<sig) then
                qq(ii,ee) = 1D0
            endif
        enddo
    enddo
    qq = qq + 2D-1*sin(2D0 * pi * xx/(xmax-xmin))
    mu = 6D-1*(xmin+xmax)
    sig = 1D-1*(xmax-xmin)
    qq = qq + exp( -(xx-mu)**2 / (2D0*sig**2))
else
    stop "Wrong init_type"
endif

end subroutine initialize

end module sub_mod1


