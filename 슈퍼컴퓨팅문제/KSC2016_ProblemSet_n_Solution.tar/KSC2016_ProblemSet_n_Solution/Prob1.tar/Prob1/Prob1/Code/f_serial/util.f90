module util_mod
    implicit none
    integer, parameter :: N_K = 16
    integer, parameter :: N_PT = 20000000

	type POINT
		real(8) :: x, y
	end type

contains

real*8 function distance(pt1, pt2)
    implicit none
    type(POINT), intent(in) :: pt1, pt2
    distance = sqrt((pt2%x-pt1%x)**2+(pt2%y-pt1%y)**2)
    return
end function distance

integer function check_diff(kmeans, kmeans_old)
    implicit none
    type(POINT), intent(in) :: kmeans(N_K), kmeans_old(N_K)
    integer i

    check_diff = 0
    do i=1,N_K
        if( distance(kmeans(i), kmeans_old(i)) > 0 ) then
            check_diff = 1
            return
        end if
    end do
    return
end function check_diff

subroutine initialize_data(kmeans, pt, kmeans_old, kindex)
    implicit none
    type(POINT), dimension(:), allocatable, intent(out) :: kmeans, pt, kmeans_old
    integer, dimension(:), allocatable, intent(out) :: kindex
    
    allocate(pt(N_PT), kindex(N_PT), kmeans(N_K), kmeans_old(N_K))

    open(unit=1, file='input.dat', form='unformatted', access='direct', recl=N_PT*16)
    read(1, rec=1) pt
    read(1, rec=2) kmeans
    close(1)
end subroutine

subroutine check_result(kmeans, kmeans_old)
    implicit none
    type(POINT), intent(in) :: kmeans(N_K)
	type(POINT), intent(inout) :: kmeans_old(N_K)

    integer :: i
    real(8) :: err

    open(unit=1, file='result.dat', form='unformatted', access='direct', recl=N_K*16)
    read(1, rec=1) kmeans_old
    close(1)

    err = 0.0
    do i=1, N_K
        err = err + distance(kmeans_old(i), kmeans(i))
    end do

    if( err < 1.0E-8 ) then 
        print *, "Result: Pass!"
    else
        print *, "Result: Fail!"
    end if
end subroutine

subroutine release_data(kmeans, pt, kmeans_old, kindex)
    implicit none
    type(POINT), dimension(:), allocatable, intent(inout) :: kmeans, pt, kmeans_old
    integer, dimension(:), allocatable, intent(inout) :: kindex

    deallocate(kmeans, pt, kmeans_old, kindex)
end subroutine

end module
