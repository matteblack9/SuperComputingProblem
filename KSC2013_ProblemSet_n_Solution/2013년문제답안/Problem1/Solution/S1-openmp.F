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
       USE omp_lib
       implicit none  
       real*8, external :: func 
       real*8 rslt,aa,bb
       real*8 rslt0
       integer n,ncall0
       integer ierr,kount,iroot
       integer myid,nproc
       real*8 time_start,time_end
       real*4 tarray(2),rate
       integer c1,c2,cr,cm

       call system_clock(count_rate=cr)
       rate = REAL(cr)
       call system_clock(c1)
       time_start = c1/rate
       time_start = omp_get_wtime()


       n=20
       n=6
       n=200000
       n=2000000000

       aa=0.0d0
       bb=1.0d0
       ncall=0
       ncall=0
       call simpson2(func,n,aa,bb,rslt,myid,nproc)

       call system_clock(c2)
       time_end = c2/rate
       time_end = omp_get_wtime()
       write(6,*) ncall, rslt
       write(6,'(4(f14.5,1x,a))') (time_end-time_start),'s', 
     &  (time_end-time_start)/60.d0,'m', (time_end-time_start)/3600.d0, 
     &   'h', (time_end-time_start)/3600.d0/24.d0,'d'
       stop
       end program simpson
!234567890
       subroutine simpson2(func,n,aa,bb,rslt,myid,nproc)
       USE intdata, ONLY : ncall
       implicit none
       integer myid,nproc
       real*8 func
       integer n
       real*8 rslt,aa,bb
       real*8 h,xx
       integer j
       integer n1,n2,istart,ifinish
       integer omp_get_num_threads,omp_get_thread_num

       ncall = 0
       rslt=0.d0
       if(mod(n,2) /= 0)then
       print*,' input error, n must be even number',n
       stop
                        endif
       n1=1 ; n2=n-1
       h=(bb-aa)/dble(n)

#ifdef OLD
!$omp parallel
       nproc = omp_get_num_threads()
       myid = omp_get_thread_num()
       if(myid.eq.0) print *,"Total number of threads is ",nproc
!$omp end parallel
!$omp parallel do private (xx,j) reduction(+:rslt) reduction(+:ncall)
       do j=1,n-1
       xx=aa+h*dble(j)
      if(mod(j,2) == 1) then
        rslt=rslt+4.0d0*func(xx)
       else
        rslt=rslt+2.0d0*func(xx)
       endif
       enddo
!$omp end parallel do
#else
!$omp parallel private(nproc,myid,istart,ifinish,xx) 
!$omp+ reduction(+:rslt,ncall)
       nproc = omp_get_num_threads()
       myid = omp_get_thread_num()
       if(myid.eq.0) print *,"Total number of threads is ",nproc
       call equal_load(n1,n2,nproc,myid,istart,ifinish)
       do j = istart,ifinish
       xx=aa+h*dble(j)
       if(mod(j,2) == 1) then
        rslt=rslt+4.0d0*func(xx)
       else
        rslt=rslt+2.0d0*func(xx)
       endif
       ncall = ncall + 1
       enddo
!$omp end parallel
#endif

       rslt = rslt + (func(aa)+func(bb))
       ncall = ncall + 2

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
       return
       end
!234567890
