module kmeans_mod
    use util_mod
    implicit none

contains

subroutine assignment_step(kmeans, pt, istart, iend, kindex)
    implicit none
    integer, intent(in) :: istart, iend
    type(POINT), intent(in) :: kmeans(N_K), pt(N_PT)
    integer, intent(out) :: kindex(N_PT)
    integer :: i,j
    real(8) :: d1, d2

    do i=istart, iend
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

subroutine update_step(kmeans, pt, istart, iend, kindex)
    implicit none
	include 'mpif.h'
    integer, intent(in) :: istart, iend, kindex(N_PT)
    type(POINT), intent(in) :: pt(N_PT)
    type(POINT), intent(out) :: kmeans(N_K)
    type(POINT)  :: local_kmeans(N_K)
    integer :: i, idx, num_pt(N_K), local_num_pt(N_K), ierr

    kmeans%x = 0.0
    kmeans%y = 0.0
    local_kmeans%x = 0.0
    local_kmeans%y = 0.0
    num_pt = 0
    local_num_pt = 0

    do i=istart, iend
        idx = kindex(i)
        local_kmeans(idx)%x = local_kmeans(idx)%x + pt(i)%x
        local_kmeans(idx)%y = local_kmeans(idx)%y + pt(i)%y
        local_num_pt(idx) = local_num_pt(idx) + 1
    end do
	
	call MPI_Allreduce(local_kmeans, kmeans, N_K*2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
	call MPI_Allreduce(local_num_pt, num_pt, N_K, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)

    kmeans%x = kmeans%x / num_pt
    kmeans%y = kmeans%y / num_pt
end subroutine

end module
