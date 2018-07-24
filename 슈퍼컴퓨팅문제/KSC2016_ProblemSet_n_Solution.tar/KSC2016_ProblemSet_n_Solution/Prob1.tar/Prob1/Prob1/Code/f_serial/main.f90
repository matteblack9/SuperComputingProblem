program kmeans_program
    use kmeans_mod
    implicit none
	include 'mpif.h'

    integer :: iter
    type(POINT), dimension(:), allocatable :: pt, kmeans, kmeans_old
    integer, dimension(:), allocatable :: kindex
	integer :: nprocs, rank, ierr
    real(8) time1, time2
    
	call MPI_Init(ierr)
	call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
	call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  1. read point data and centers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call initialize_data(kmeans, pt, kmeans_old, kindex)

    time1 = MPI_Wtime()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  2. execute k-means clustering
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do iter = 0, 99
        call assignment_step(kmeans, pt, kindex)
        kmeans_old = kmeans
        call update_step(kmeans, pt, kindex)
        if( check_diff(kmeans, kmeans_old) == 0 ) exit
    end do

    time2 = MPI_Wtime()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  3. report and verify the results
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if( rank == 0 ) then
	    print *, 'Iteration: ', iter, ' Execution Time: ', time2-time1
		call check_result(kmeans, kmeans_old);
	end if
    call release_data(kmeans, pt, kmeans_old, kindex)

	call MPI_Finalize(ierr)
end program
