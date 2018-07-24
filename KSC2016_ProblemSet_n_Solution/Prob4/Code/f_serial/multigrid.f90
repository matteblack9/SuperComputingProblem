module multigrid

    implicit none
    
contains

    subroutine restriction(uf,uc,nc)
    
    implicit none
    
    real(8) :: uf(2*nc,2*nc), uc(nc,nc)
    integer(4) :: nc
    integer(4) :: ic, jc
    
    do ic=1,nc
        do jc=1, nc
            uc(ic,jc) = 1.0/4.0*(uf(2*ic-1,2*jc-1)+uf(2*ic,2*jc-1)+uf(2*ic-1,2*jc)+uf(2*ic,2*jc))
        enddo
    enddo
    
    end subroutine restriction
    
    subroutine interpolation(uc,uf,nf)
    
    implicit none

    real(8) :: uf(nf,nf), uc(nf/2,nf/2)
    integer(4) :: nf
    integer(4) :: i,j,ic,jc,nc,count=0
    real(8),allocatable :: uc_ghost(:,:)

    nc = nf/2
    allocate(uc_ghost(0:nc+1,0:nc+1))

    do jc=1,nc
        do ic=1,nc
            uc_ghost(ic,jc) = uc(ic,jc)
        enddo
    enddo

    do ic=1,nc
        uc_ghost(ic,0) = -uc_ghost(ic,1)
        uc_ghost(ic,nc+1) = -uc_ghost(ic,nc)
    enddo

    do jc=1,nc
        uc_ghost(0,jc) = -uc_ghost(1,jc)
        uc_ghost(nc+1,jc) = -uc_ghost(nc,jc)
    enddo

    uc_ghost(0,0) = -uc_ghost(1,1);
    uc_ghost(0,nc+1) = -uc_ghost(1,nc);
    uc_ghost(nc+1,0) = -uc_ghost(nc,1);
    uc_ghost(nc+1,nc+1) = -uc_ghost(nc,nc);

    do ic=1,nc
        do jc=1,nc
            uf(2*ic-1,2*jc-1) = (9.0*uc_ghost(ic,jc) + 3.0*uc_ghost(ic-1,jc) + 3.0*uc_ghost(ic,jc-1) + 1.0*uc_ghost(ic-1,jc-1))/16.0
            uf(2*ic,  2*jc-1) = (9.0*uc_ghost(ic,jc) + 3.0*uc_ghost(ic+1,jc) + 3.0*uc_ghost(ic,jc-1) + 1.0*uc_ghost(ic+1,jc-1))/16.0
            uf(2*ic-1,2*jc)   = (9.0*uc_ghost(ic,jc) + 3.0*uc_ghost(ic-1,jc) + 3.0*uc_ghost(ic,jc+1) + 1.0*uc_ghost(ic-1,jc+1))/16.0
            uf(2*ic,  2*jc)   = (9.0*uc_ghost(ic,jc) + 3.0*uc_ghost(ic+1,jc) + 3.0*uc_ghost(ic,jc+1) + 1.0*uc_ghost(ic+1,jc+1))/16.0
        enddo
    enddo
    
    deallocate(uc_ghost);

    
    end subroutine interpolation
    
end module multigrid
    
