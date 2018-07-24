module multigrid

    use mpi
    implicit none
    
contains

    subroutine restriction(uf,uc,ncx,ncy)
    
    implicit none
    
    real(8) :: uf(2*ncx,2*ncy), uc(ncx,ncy)
    integer(4) :: ncx, ncy
    integer(4) :: ic, jc
    
    do jc=1, ncy
        do ic=1,ncx
            uc(ic,jc) = 1.0/4.0*(uf(2*ic-1,2*jc-1)+uf(2*ic,2*jc-1)+uf(2*ic-1,2*jc)+uf(2*ic,2*jc))
        enddo
    enddo
    
    end subroutine restriction
    
    subroutine restriction_mlevel(uf,uc,ncx,ncy,my_rank,distance)
    
    implicit none
    
    real(8) :: uf(2*ncx,2*ncy), uc(ncx,ncy)
    integer(4) :: ncx, ncy, my_rank, distance
    integer(4) :: ic, jc
    real(8), allocatable :: uf_next(:,:)
    integer(4) :: status(MPI_STATUS_SIZE),ierr

    allocate(uf_next(2*ncx,ncy))

    if(mod(my_rank,2*distance).eq.0) then
        call MPI_RECV(uf_next, ncx*2, MPI_REAL8, my_rank+distance, distance, MPI_COMM_WORLD, status, ierr) 
        jc=1
        do ic=1,ncx
            uc(ic,jc) = 1.0/4.0*(uf(2*ic-1,1)+uf(2*ic,1)+uf_next(2*ic-1,1)+uf_next(2*ic,1))
        enddo
    else if(mod(my_rank,2*distance).eq.distance) then
        call MPI_SEND(uf, ncx*2, MPI_REAL8, my_rank-distance, distance, MPI_COMM_WORLD, ierr)
    endif

    deallocate(uf_next)
    
    end subroutine restriction_mlevel
    
    subroutine interpolation_mlevel(uc,uf,nfx,nfy,my_rank,mpisize,distance)
    
    implicit none

    real(8) :: uf(nfx,nfy), uc(nfx/2,nfy)
    integer(4) :: nfx,nfy,my_rank,mpisize,distance
    integer(4) :: i,j,ic,jc,ncx,ncy,count=0
    real(8),allocatable :: uc_ghost(:,:),uf_double(:,:)
    integer(4) :: status(MPI_STATUS_SIZE),ierr

    ncx = nfx/2
    ncy = nfy
    allocate(uc_ghost(0:ncx+1,0:ncy+1))
    allocate(uf_double(nfx,2*nfy))

    if(mod(my_rank,2*distance).eq.0) then
        do jc=1,ncy
            do ic=1,ncx
                uc_ghost(ic,jc) = uc(ic,jc)
            enddo
        enddo

        do jc=1,ncy
            uc_ghost(0,jc) = -uc_ghost(1,jc)
            uc_ghost(ncx+1,jc) = -uc_ghost(ncx,jc)
        enddo

        if(my_rank.ne.0) then
            call MPI_RECV(uc_ghost(0,0),ncx+2,MPI_REAL8,my_rank-2*distance, 0, MPI_COMM_WORLD, status, ierr)
            call MPI_SEND(uc_ghost(0,1),ncx+2,MPI_REAL8,my_rank-2*distance, 0, MPI_COMM_WORLD, ierr)
        endif

        if(my_rank+2*distance.ne.mpisize) then
            call MPI_SEND(uc_ghost(0,ncy  ),ncx+2,MPI_REAL8,my_rank+2*distance, 0, MPI_COMM_WORLD, ierr)
            call MPI_RECV(uc_ghost(0,ncy+1),ncx+2,MPI_REAL8,my_rank+2*distance, 0, MPI_COMM_WORLD, status, ierr)
        endif

        if(my_rank.eq.0) then
            do ic=1,ncx
                uc_ghost(ic,0) = -uc_ghost(ic,1)
            enddo
            uc_ghost(0,0) = -uc_ghost(1,1);
            uc_ghost(ncx+1,0) = -uc_ghost(ncx,1);
        endif

        if(my_rank+2*distance.eq.mpisize) then
            do ic=1,ncx
                uc_ghost(ic,ncy+1) = -uc_ghost(ic,ncy)
            enddo
            uc_ghost(0,ncy+1) = -uc_ghost(1,ncy);
            uc_ghost(ncx+1,ncy+1) = -uc_ghost(ncx,ncy);
        endif

        do jc=1,ncy
            do ic=1,ncx
                uf_double(2*ic-1,2*jc-1) = (9.0*uc_ghost(ic,jc) + 3.0*uc_ghost(ic-1,jc) &
                                          + 3.0*uc_ghost(ic,jc-1) + 1.0*uc_ghost(ic-1,jc-1))/16.0
                uf_double(2*ic-1,2*jc)   = (9.0*uc_ghost(ic,jc) + 3.0*uc_ghost(ic-1,jc) &
                                          + 3.0*uc_ghost(ic,jc+1) + 1.0*uc_ghost(ic-1,jc+1))/16.0
                uf_double(2*ic,  2*jc-1) = (9.0*uc_ghost(ic,jc) + 3.0*uc_ghost(ic+1,jc) &
                                          + 3.0*uc_ghost(ic,jc-1) + 1.0*uc_ghost(ic+1,jc-1))/16.0
                uf_double(2*ic,  2*jc)   = (9.0*uc_ghost(ic,jc) + 3.0*uc_ghost(ic+1,jc) &
                                          + 3.0*uc_ghost(ic,jc+1) + 1.0*uc_ghost(ic+1,jc+1))/16.0
            enddo
        enddo
        do j=1,nfy
            do i=1,nfx
                uf(i,j) = uf_double(i,j)
            enddo
        enddo
        call MPI_SEND(uf_double(1,2), nfx, MPI_REAL8, my_rank+distance, distance, MPI_COMM_WORLD, ierr)
    else if(mod(my_rank,2*distance).eq.distance) then
        call MPI_RECV(uf(1,1), nfx, MPI_REAL8, my_rank-distance, distance, MPI_COMM_WORLD, status, ierr)
    endif
    
    deallocate(uf_double);
    deallocate(uc_ghost);
    
    end subroutine interpolation_mlevel
    
    subroutine interpolation_mpi(uc,uf,nfx,nfy,my_rank,prev_rank,next_rank)
    
    implicit none

    real(8) :: uf(nfx,nfy), uc(nfx/2,nfy/2)
    integer(4) :: nfx,nfy,my_rank,prev_rank,next_rank
    integer(4) :: i,j,ic,jc,ncx,ncy,count=0
    real(8),allocatable :: uc_ghost(:,:)
    integer(4) :: status(MPI_STATUS_SIZE), ierr

    ncx = nfx/2
    ncy = nfy/2
    allocate(uc_ghost(0:ncx+1,0:ncx+1))

    do jc=1,ncy
        do ic=1,ncx
            uc_ghost(ic,jc) = uc(ic,jc)
        enddo
    enddo

    do jc=1,ncy
        uc_ghost(0,jc) = -uc_ghost(1,jc)
        uc_ghost(ncx+1,jc) = -uc_ghost(ncx,jc)
    enddo

    if(prev_rank.ne.MPI_PROC_NULL) then
        call MPI_RECV(uc_ghost(0,0),ncx+2,MPI_REAL8,prev_rank, 0, MPI_COMM_WORLD, status, ierr)
        call MPI_SEND(uc_ghost(0,1),ncx+2,MPI_REAL8,prev_rank, 0, MPI_COMM_WORLD, ierr)
    endif

    if(next_rank.ne.MPI_PROC_NULL) then
        call MPI_SEND(uc_ghost(0,ncy  ),ncx+2,MPI_REAL8,next_rank, 0, MPI_COMM_WORLD, ierr)
        call MPI_RECV(uc_ghost(0,ncy+1),ncx+2,MPI_REAL8,next_rank, 0, MPI_COMM_WORLD, status, ierr)
    endif

    if(prev_rank.eq.MPI_PROC_NULL) then
        do ic=1,ncx
            uc_ghost(ic,0) = -uc_ghost(ic,1)
        enddo
        uc_ghost(0,0) = -uc_ghost(1,1);
        uc_ghost(ncx+1,0) = -uc_ghost(ncx,1);
    endif

    if(next_rank.eq.MPI_PROC_NULL) then
        do ic=1,ncx
            uc_ghost(ic,ncy+1) = -uc_ghost(ic,ncy)
        enddo
        uc_ghost(0,ncy+1) = -uc_ghost(1,ncy);
        uc_ghost(ncx+1,ncy+1) = -uc_ghost(ncx,ncy);
    endif



    do jc=1,ncy
        do ic=1,ncx
            uf(2*ic-1,2*jc-1) = (9.0*uc_ghost(ic,jc) + 3.0*uc_ghost(ic-1,jc) + 3.0*uc_ghost(ic,jc-1) + 1.0*uc_ghost(ic-1,jc-1))/16.0
            uf(2*ic-1,2*jc)   = (9.0*uc_ghost(ic,jc) + 3.0*uc_ghost(ic-1,jc) + 3.0*uc_ghost(ic,jc+1) + 1.0*uc_ghost(ic-1,jc+1))/16.0
            uf(2*ic,  2*jc-1) = (9.0*uc_ghost(ic,jc) + 3.0*uc_ghost(ic+1,jc) + 3.0*uc_ghost(ic,jc-1) + 1.0*uc_ghost(ic+1,jc-1))/16.0
            uf(2*ic,  2*jc)   = (9.0*uc_ghost(ic,jc) + 3.0*uc_ghost(ic+1,jc) + 3.0*uc_ghost(ic,jc+1) + 1.0*uc_ghost(ic+1,jc+1))/16.0
        enddo
    enddo
    
    deallocate(uc_ghost);

    
    end subroutine interpolation_mpi
    
end module multigrid
    
