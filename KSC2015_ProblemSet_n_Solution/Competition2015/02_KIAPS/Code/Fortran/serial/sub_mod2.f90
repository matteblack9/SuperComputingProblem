module sub_mod2

use sub_mod1

implicit none

contains

subroutine interface_flux(qq, fstar, ib, speed, np, ne)

implicit none

integer, intent(in) :: np, ne
real(8), intent(in) :: speed
real(8), dimension(np, ne), intent(in) :: qq
real(8), dimension(2, ne), intent(out) :: fstar
real(8), intent(in) :: ib(np, 2)
integer :: ee, el, er, l, r
real(8) :: nl, nr, ql, qr, h, qb(2,ne)

fstar = 0D0
l = 1 ! left
r = 2 ! right

do ee = 1, ne
    ! qb(l,ee) is the left edge interpolated value of qq at element ee
    ! qb(r,ee) is the right edge interpolated value of qq at element ee
    qb(l,ee) = dot_product(ib(:,1), qq(:,ee))
    qb(r,ee) = dot_product(ib(:,2), qq(:,ee))
enddo

ee = 1
! calculating numerical flux (fstar) with periodic boundary condition
fstar(l,ee) = -((qb(r,ne)+qb(l,ee))/2D0*speed + abs(speed)*(qb(r,ne)-qb(l,ee))/2D0)
fstar(r,ne) = -fstar(l,1)
do ee = 2, ne
    fstar(l,ee) = -((qb(r,ee-1)+qb(l,ee))/2D0*speed + abs(speed)*(qb(r,ee-1)-qb(l,ee))/2D0)
    fstar(r,ee-1) = -fstar(l,ee)
enddo

end subroutine interface_flux 

subroutine rhs(qq, rr, dv, df, ib, speed, ne, np)

implicit none

integer, intent(in) :: ne, np
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

call interface_flux(qq, fstar, ib, speed, np, ne)

do ee = 1, ne
    do ii = 1, np
        do jj = 1, 2
            rr(ii,ee) = rr(ii,ee) - df(ii,jj,ee)*fstar(jj,ee)
        enddo
    enddo
enddo


    
end subroutine rhs

end module sub_mod2


