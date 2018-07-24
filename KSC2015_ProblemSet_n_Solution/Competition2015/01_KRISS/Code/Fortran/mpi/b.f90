!234567890
       program problem_fd
       implicit none
       include 'mpif.h'
       integer i,n1,n2,j,jsta,jend
       integer iter,niter
       integer istatus(MPI_STATUS_SIZE)
       integer ierr,myid,nprocs
       integer iprev,inext,ista,iend
       integer isd1,isd2,irv1,irv2,itag,iroot
       real*8 xi,xf,dx
       real*8 tmr
       real*8, allocatable:: ar(:),br(:)
       real*8, external :: genvv
       real*8 ptmr
       real*8 tic,toc

!  do not change -----{
       n1=1
       n2=100
       n2=100000000
       niter=3
!  do not change -----}

       allocate ( ar(n2), br(n2) )

       xi=0.d0
       xf=1.d0
       dx=(xf-xi)/dble(n2-n1)
       do i=n1,n2
       br(i)=xi+dble(i-n1)*dx
       enddo
       CALL MPI_INIT(ierr)
       tic=MPI_WTIME()
       CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
       CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
       CALL para_range(n1,n2,nprocs,myid,ista,iend)

       write(*,*) 'rank :', myid, 'ista =', ista, 'iend =', iend

       jsta=ista
       jend=iend
       if(myid == 0) jsta=n1+1
       if(myid == nprocs-1) jend=n2-1
       inext=myid+1
       iprev=myid-1
       if(myid == nprocs-1) inext=MPI_PROC_NULL
       if(myid == 0) iprev=MPI_PROC_NULL
       do i=ista,iend
       br(i)=xi+dble(i-n1)*dx
       enddo
       do iter=1,niter
       itag=101
       CALL MPI_ISEND(br(iend), 1,MPI_REAL8, inext,itag,MPI_COMM_WORLD,isd1,ierr)
       CALL MPI_ISEND(br(ista), 1,MPI_REAL8, iprev,itag,MPI_COMM_WORLD,isd2,ierr)
       CALL MPI_IRECV(br(ista-1), 1,MPI_REAL8, iprev,itag,MPI_COMM_WORLD,irv1,ierr)
       CALL MPI_IRECV(br(iend+1), 1,MPI_REAL8, inext,itag,MPI_COMM_WORLD,irv2,ierr)
       CALL MPI_WAIT(isd1,istatus,ierr)
       CALL MPI_WAIT(isd2,istatus,ierr)
       CALL MPI_WAIT(irv1,istatus,ierr)
       CALL MPI_WAIT(irv2,istatus,ierr)
       do j=jsta,jend
!  do not change -----{
       ar(j)=(br(j-1)+br(j+1))/4.d0+br(j)/2.d0+1.d0/genvv(br(j))
!  do not change -----}
       enddo
       do i=ista,iend
!  do not change -----{
       br(i)=ar(i)
!  do not change -----}
       enddo
       enddo
      
       ptmr=0.d0
       do j=jsta,jend
       ptmr=ptmr+ar(j)
       enddo
       iroot=0
       CALL MPI_REDUCE(ptmr, tmr, 1,MPI_REAL8, MPI_SUM,iroot,MPI_COMM_WORLD,ierr)
       if(myid == 0) print *,'tmr = ',tmr
       toc=MPI_WTIME()
       if(myid == 0) print*, toc-tic,' sec'   
   
       deallocate( ar, br )

       CALL MPI_FINALIZE(ierr)
       stop
       end program problem_fd
!234567890
       real*8 function genvv(x)
       implicit none
       real*8 x
       
       genvv=x**2+x**4+x**6+exp(-x**2)+cos(x)+sin(x)+tan(x)
       return
       end
!234567890
       subroutine para_range(n1,n2,nprocs,myid,ista,iend)
       implicit none
       integer n1,n2,nprocs,myid,ista,iend
       integer iwork1,iwork2

       iwork1 = (n2 - n1 + 1) / nprocs
       iwork2 = mod(n2 - n1 + 1, nprocs)
       ista = myid * iwork1 + n1 + min(myid, iwork2)
       iend = ista + iwork1 - 1
       if(iwork2 > myid) iend = iend + 1
       return
       end
