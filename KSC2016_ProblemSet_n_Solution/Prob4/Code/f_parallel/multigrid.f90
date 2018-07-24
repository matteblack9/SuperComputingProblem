module multigrid

    use mpi

    implicit none
    
contains

    subroutine restriction(uf,uc,nx,ny)
    
    implicit none
    
    real(8) :: uf(2*nx,2*ny), uc(nx,ny)
    integer(4) :: nx,ny
    integer(4) :: ic, jc, ierr
    
    do jc=1, ny
        do ic=1,nx
            uc(ic,jc) = 1.0/4.0*(uf(2*ic-1,2*jc-1)+uf(2*ic,2*jc-1)+uf(2*ic-1,2*jc)+uf(2*ic,2*jc))
        enddo
    enddo

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    end subroutine restriction
    
    subroutine interpolation(uc,uf,nx,ny,myrank,ncpus)
    
    implicit none

    real(8) :: uf(nx,ny), uc(nx/2,ny/2)
    integer(4) :: nx,ny,myrank,ncpus
    integer(4) :: i,j,ic,jc,ncx,ncy,count=0
    integer(4) :: prevrank, nextrank
    integer(4) :: ierr, request(2), status(MPI_STATUS_SIZE)
    real(8),allocatable :: uc_ghost(:,:)

    ncx = nx/2
    ncy = ny/2

    prevrank = mod(myrank-1+ncpus,ncpus)
    nextrank = mod(myrank+1,ncpus)

    allocate(uc_ghost(0:ncx+1,0:ncy+1))

    do jc=1,ncy
        do ic=1,ncx
            uc_ghost(ic,jc) = uc(ic,jc)
        enddo
    enddo

    do jc=1,ncy
        uc_ghost(0,jc) = -uc_ghost(1,jc)
        uc_ghost(ncx+1,jc) = -uc_ghost(ncx,jc)
    enddo


    if(ncpus.gt.1) then
        call MPI_IRECV(uc_ghost(0,ncy+1),ncx+2,MPI_REAL8,nextrank,1002,MPI_COMM_WORLD,request(1),ierr)
        call MPI_ISEND(uc_ghost(0,1),ncx+2,MPI_REAL8,prevrank,1002,MPI_COMM_WORLD,request(2),ierr)
        call MPI_WAITALL(2,request,status,ierr)

        call MPI_IRECV(uc_ghost(0,0),ncx+2,MPI_REAL8,prevrank,1002,MPI_COMM_WORLD,request(1),ierr)
        call MPI_ISEND(uc_ghost(0,ncy),ncx+2,MPI_REAL8,nextrank,1002,MPI_COMM_WORLD,request(2),ierr)
        call MPI_WAITALL(2,request,status,ierr)
    endif

    do ic=1,ncx
        if(myrank.eq.0) then
            uc_ghost(ic,0) = -uc_ghost(ic,1)
        endif
        if(myrank.eq.ncpus-1) then
            uc_ghost(ic,ncy+1) = -uc_ghost(ic,ncy)
        endif
    enddo

    if(myrank.eq.0) then
        uc_ghost(0,0) = -uc_ghost(1,1);
        uc_ghost(ncx+1,0) = -uc_ghost(ncx,1);
    endif
    if(myrank.eq.ncpus-1) then
        uc_ghost(0,ncy+1) = -uc_ghost(1,ncy);
        uc_ghost(ncx+1,ncy+1) = -uc_ghost(ncx,ncy);
    endif

    do jc=1,ncy
        do ic=1,ncx
            uf(2*ic-1,2*jc-1) = (9.0*uc_ghost(ic,jc) + 3.0*uc_ghost(ic-1,jc) + 3.0*uc_ghost(ic,jc-1) + 1.0*uc_ghost(ic-1,jc-1))/16.0
            uf(2*ic,  2*jc-1) = (9.0*uc_ghost(ic,jc) + 3.0*uc_ghost(ic+1,jc) + 3.0*uc_ghost(ic,jc-1) + 1.0*uc_ghost(ic+1,jc-1))/16.0
            uf(2*ic-1,2*jc)   = (9.0*uc_ghost(ic,jc) + 3.0*uc_ghost(ic-1,jc) + 3.0*uc_ghost(ic,jc+1) + 1.0*uc_ghost(ic-1,jc+1))/16.0
            uf(2*ic,  2*jc)   = (9.0*uc_ghost(ic,jc) + 3.0*uc_ghost(ic+1,jc) + 3.0*uc_ghost(ic,jc+1) + 1.0*uc_ghost(ic+1,jc+1))/16.0
        enddo
    enddo
    
    deallocate(uc_ghost);

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    end subroutine interpolation
    
end module multigrid
    
