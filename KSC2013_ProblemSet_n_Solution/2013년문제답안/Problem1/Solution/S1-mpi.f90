       module intdata
       implicit none  
       private
       save
       integer ncall
       public :: ncall

       end module intdata
!234567890
       program simpson
       USE intdata, ONLY : ncall
       implicit none  
       include 'mpif.h'
       real*8, external :: func 
       real*8 rslt,aa,bb
       real*8 rslt0
       integer n,ncall0
       integer ierr,kount,iroot
       integer myid,nproc
       real*8 time1,time2,time_start,time_end

       call MPI_INIT( ierr )
       call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
       call MPI_COMM_SIZE( MPI_COMM_WORLD, nproc, ierr )
       if(myid == 0 .and. nproc > 1) print *,  nproc," processes are alive"
       if(myid == 0 .and. nproc ==1) print *,  nproc," process is alive"

       time_start=MPI_WTIME()

       n=20
       n=6
       n=200000
       n=2000000000

       aa=0.0d0
       bb=1.0d0
       ncall=0
       ncall=0
       call simpson2(func,n,aa,bb,rslt,myid,nproc)

       iroot=0 ; kount=1
       call MPI_REDUCE(rslt,rslt0,kount,MPI_DOUBLE_PRECISION,MPI_SUM,iroot,MPI_COMM_WORLD,ierr)
       if(myid == 0) write(6,*) n,rslt0,' n,rslt0'
       iroot=0 ; kount=1
       call MPI_REDUCE(ncall,ncall0,kount,MPI_INTEGER,MPI_SUM,iroot,MPI_COMM_WORLD,ierr)
       if(myid == 0) write(6,*) ncall0,' ncall0'

       time_end=MPI_WTIME()
       if(myid == 0) then     ! -------=== { process id =0
       write(6,'(4(f14.5,1x,a))') (time_end-time_start),'s', (time_end-time_start)/60.d0,'m', (time_end-time_start)/3600.d0,'h', (time_end-time_start)/3600.d0/24.d0,'d'
                     endif    ! -------=== } process id =0
       call MPI_FINALIZE(ierr)
       stop
       end program simpson
!234567890
       subroutine simpson2(func,n,aa,bb,rslt,myid,nproc)
       implicit none
       integer myid,nproc
       real*8 func
       integer n
       real*8 rslt,aa,bb
       real*8 h,xx
       integer j
       integer n1,n2,istart,ifinish

       rslt=0.d0
       if(mod(n,2) /= 0)then
       print*,' input error, n must be even number',n
       stop
                        endif
       n1=1 ; n2=n-1
       call equal_load(n1,n2,nproc,myid,istart,ifinish)

       h=(bb-aa)/dble(n)
       if(myid == 0) rslt=(func(aa)+func(bb))
!!     do j=1,n-1
       do j=istart,ifinish
       xx=aa+h*dble(j)
       if(mod(j,2) == 1) then
        rslt=rslt+4.0d0*func(xx)
                         else
        rslt=rslt+2.0d0*func(xx)
                         endif
       enddo

       rslt=rslt*h/3.0d0
       return
       end
       subroutine equal_load(n1,n2,nproc,myid,istart,ifinish)
       implicit none
       integer nproc,myid,istart,ifinish,n1,n2
       integer iw1,iw2
       iw1=(n2-n1+1)/nproc ; iw2=mod(n2-n1+1,nproc)
       istart=myid*iw1+n1+min(myid,iw2)
       ifinish=istart+iw1-1 ; if(iw2 > myid) ifinish=ifinish+1
!      print*, n1,n2,myid,nproc,istart,ifinish
       if(n2 < istart) ifinish=istart-1
       return
       end
!234567890
       real*8 function func(x)
       USE intdata, ONLY : ncall
       implicit none
       real*8 x
!      func=x*x
       func=4.d0/(x*x+1.d0)
       ncall=ncall+1
       return
       end
!234567890
