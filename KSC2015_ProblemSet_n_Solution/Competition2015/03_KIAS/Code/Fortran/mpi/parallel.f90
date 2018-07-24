module data
   real,allocatable,dimension(:) :: base
end module data
!----
      recursive subroutine domaindecomp(nmem, valmin, valmax, Comm)
      use data
      include 'mpif.h'
!      real,allocatable,dimension(:) :: base
!      common /data/ base
      real, allocatable, target:: reallocbase(:)
      real, allocatable, target :: sbase(:)
      integer*8 nowmem,nsend,nrecv
      real valmin,valmax
      integer Comm
      integer myid,nid
      integer*8, intent(in out):: nmem
      integer status(MPI_STATUS_SIZE)
      integer*8 left,right
      real swaptmp
      real halfval 
      integer subgroupid, nsubgroup, subgroupsize
      integer dest,src
      integer key
      integer newcom
      real newvalmin,newvalmax

      halfval = (valmax+valmin)*0.5

      call MPI_Comm_size(Comm,nid,ierror)
      if(nid .eq.1) return
      call MPI_Comm_rank(Comm,myid,ierror)

       left =1
       right = nmem+1
       if(myid .lt. nid/2) then
        do while(left .lt. right)
           if(base(left) .ge. halfval) then
              right = right - 1
              swaptmp = base(left)
              base(left) = base(right)
              base(right) = swaptmp
           else
              left = left + 1
           endif
        enddo
       else
       do while(left .lt. right)
          if(base(left) .lt. halfval) then
             right = right - 1
             swaptmp = base(left)
             base(left) = base(right)
             base(right) = swaptmp
          else
             left = left + 1
          endif
       enddo
       endif

      nsend = nmem - left +1

      nsubgroup = 2
      subgroupsize = nid/nsubgroup
      subgroupid = myid/subgroupsize

      dest = mod(myid+subgroupsize+nid,nid)
      src  = mod(myid-subgroupsize+nid,nid)
      call MPI_Sendrecv(nsend,1, MPI_INTEGER8,dest,0, nrecv,1, MPI_INTEGER8, src, 0, Comm, status,ierror)

      allocate(sbase(1:nsend))
      do i = 1, nsend
         sbase(i) = base(nmem-nsend+i)
      enddo
      call MPI_Sendrecv(sbase(1), nsend, MPI_REAL, dest,0,base(nmem-nsend+1),nrecv,MPI_REAL, src, 0, Comm,status,ierror)
      nmem = nrecv + nmem - nsend
      deallocate(sbase)

     key = mod(myid,subgroupsize)
     newvalmin = (valmax-valmin)/nsubgroup * subgroupid + valmin
     newvalmax = (valmax-valmin)/nsubgroup * (subgroupid+1) + valmin
     if(subgroupid .eq.1) newvalmax = valmax
     call MPI_Comm_split(Comm,subgroupid,key, newcom,ierror)
     call domaindecomp(nmem, newvalmin, newvalmax, newcom)
     call MPI_Comm_free(newcom, ierror)


      return
      end


      program main
      use data
      include 'mpif.h'
      parameter(maxnp= 100000000)
      integer i,j
      real timearray(2),lvalmin,lvalmax
      real*8 time1,time2,walltime
      integer iseed
      character*80 arg
!      real, allocatable, dimension( : ) ::  base
!      common /data/ base
      integer*8 np
      real step,valmin,valmax
      integer myid,nid



      call MPI_Init(ierror)
      call MPI_Comm_size(MPI_COMM_WORLD,nid,ierror)
      call MPI_Comm_rank(MPI_COMM_WORLD,myid,ierror)

      time1 = MPI_Wtime()
      step = (1.-0.)/nid
      allocate(base(1:maxnp*3))
      iseed = -1*(myid+1284)
      do j = 1, maxnp
         base(j) = ran2(iseed)
      enddo
      np = maxnp

      valmin = 0
      valmax = 1
       call domaindecomp(np,valmin,valmax,MPI_COMM_WORLD)

      print *, 'p', myid, 'has ',np,'members'
      call MPI_Barrier(MPI_COMM_WORLD,ierror)
      time2 =  MPI_Wtime()
      walltime = time2-time1
      if(myid.eq.0) print *,'Wallclock time=',walltime
      lvalmin = 1./nid *myid
      lvalmax = 1./nid *(myid+1)
      call Check(base,np,lvalmin, lvalmax)
      call MPI_Finalize(ierror)
      stop
      end

      subroutine Check(base,np,lvalmin,lvalmax)
      integer*8 np,i
      real base(np),lvalmin,lvalmax
      do i = 1, np
         if(base(i).lt. lvalmin .or. base(i).ge.lvalmax)then
            print *,'Error in DD ', base(i),i,lvalmin,lvalmax
            stop
         endif
      enddo

      return
      end
      FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, &
      IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791, &
      NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END
!  (C) Copr. 1986-92 Numerical Recipes Software +)-*1a311.
