[문제1번 포트란 기본 코드]
      module intdata
       implicit none  
       private
       save
       integer ncall
       public :: ncall
       end module intdata

      program simpson
       USE intdata, ONLY : ncall
       implicit none  
       real*8, external :: func 
       real*8 rslt,aa,bb
       integer n
	   integer c1,c2,cr,cm
	   real rate,time_start,time_end
	   call system_clock(count_rate = cr)
	   rate = REAL(cr)
	   call system_clock(c1)
	   time_start = c1/rate
       n=2000000000
       aa=0.0d0
       bb=1.0d0
       ncall=0
       call simpson0(func,n,aa,bb,rslt)
	   call system_clock(c2)
	   time_end = c2/rate
      write(6,*) n,rslt,' n,rslt',ncall
      write(6,*) 'Wallclock time',time_end-time_start
       stop
       end program simpson

       subroutine simpson0(func,n,aa,bb,rslt)
       implicit none
       real*8 func
       integer n
       real*8 rslt,aa,bb
       real*8 h,xx
       integer j
       logical lodd
       rslt=0.d0
       if(mod(n,2) /= 0)then
       print*,' input error, n must be even number',n
       stop
       endif

       h=(bb-aa)/dble(n)
       rslt=(func(aa)+func(bb))
       lodd=.true.
       do j=1,n-1
       xx=aa+h*dble(j)
       if(lodd)then
       rslt=rslt+4.0d0*func(xx)
       else
       rslt=rslt+2.0d0*func(xx)
       endif
       lodd=(.not. lodd)
       enddo

       rslt=rslt*h/3.0d0
       return
       end

       real*8 function func(x)
       USE intdata, ONLY : ncall
       implicit none
       real*8 x
       func=4.d0/(x*x+1.d0)
       ncall=ncall+1
       return
       end
