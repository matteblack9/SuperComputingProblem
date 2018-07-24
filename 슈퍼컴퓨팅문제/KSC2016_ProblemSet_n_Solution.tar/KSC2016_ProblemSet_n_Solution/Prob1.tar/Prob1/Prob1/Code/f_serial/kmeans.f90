module kmeans_mod
    use util_mod
    implicit none

contains

subroutine assignment_step(kmeans, pt, kindex)
    implicit none
    type(POINT), intent(in) :: kmeans(N_K), pt(N_PT)
    integer, intent(out) :: kindex(N_PT)
    integer :: i,j
    real(8) :: d1, d2

    do i=1, N_PT
        kindex(i) = 1
        d1 = distance(pt(i), kmeans(1))
        do j=2,N_K
            d2 = distance(pt(i), kmeans(j))
            if( d1 > d2 ) then
                d1 = d2
                kindex(i) = j
            end if
        end do
    end do
end subroutine

subroutine update_step(kmeans, pt, kindex)
    implicit none
    integer, intent(in) :: kindex(N_PT)
    type(POINT), intent(in) :: pt(N_PT)
    type(POINT), intent(out) :: kmeans(N_K)
    integer :: i, idx, num_pt(N_K), ierr

    kmeans%x = 0.0
    kmeans%y = 0.0
    num_pt = 0

    do i=1, N_PT
        idx = kindex(i)
        kmeans(idx)%x = kmeans(idx)%x + pt(i)%x
        kmeans(idx)%y = kmeans(idx)%y + pt(i)%y
        num_pt(idx) = num_pt(idx) + 1
    end do
	
    kmeans%x = kmeans%x / num_pt
    kmeans%y = kmeans%y / num_pt
end subroutine

end module
