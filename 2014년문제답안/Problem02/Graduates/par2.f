c234567
      function potential(r, np, jstart, jfinal)
      double precision potential
      real r(3,np)
      integer np, jstart, jfinal
      potential = 0
      do j=jstart+1, jfinal
         do i = j+1, np
            x = r(1,i) - r(1,j)
            y = r(2,i) - r(2,j)
            z = r(3,i) - r(3,j)
            dist = sqrt(x*x+y*y+z*z)
            if(dist .gt.0) potential = potential + 1.d0/dist
         enddo
      enddo
      jjstart = min(np-jfinal,np/2)
      jjfinal = min(np/2, np-jstart)
      do j = jjstart + 1, jjfinal 
         do i = j +1, np
            x = r(1,i) - r(1,j)
            y = r(2,i) - r(2,j)
            z = r(3,i) - r(3,j)
            dist = sqrt(x*x+y*y+z*z)
            if(dist .gt.0) potential = potential + 1.d0/dist
         enddo
      enddo
      return
      end

      function getrandomnp(istep, niter)
      integer getrandomnp
      integer iseed
      common /seed/ iseed
      if(mod(istep-1,3).eq.0) then
         getrandomnp = (100000*ran2(iseed))
      else
         getrandomnp = (500*ran2(iseed))
      endif
      return
      end
      program mpi_c
      include 'mpif.h'
      parameter(maxnp=5000000)
	  integer iREADY,iWRITING,iFINISH,NP_TAG,iR_TAG, iPOT_TAG
      parameter(iREADY=0,iWRITING=-1,iFINISH=-991,NP_TAG=99,iR_TAG=990)
      parameter(iPOT_TAG=123)
      integer i,j,k
      integer np,niter
      real r(3,maxnp)
      double precision totpotent,potent,potential
      integer finish,ready
      integer iseed
      common /seed/ iseed
      real timearray(2), time1,time2,ran2
	  external ran2,potential
	  integer*8 nlocalsize, nwork, nsplit, njump

      integer mstatus(MPI_STATUS_SIZE), request,ierror
	  external real etime


      finish = iFINISH

      call MPI_INIT(ierror)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nid, ierror)
      call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierror)

      niter = 50*3
      iseed = -9

      xx = ran2(iseed)

      time1 = etime(timearray)
      nlocalsize = 1000*1000

      if(myid.eq.0) then
         nsend = 0
		 nrecv = 0
         do i = 1, niter
            np = getrandomnp(i,niter)
            do j = 1,  np
               r(1,j) = 2*ran2(iseed) -1
               r(2,j) = 2*ran2(iseed) -1
               r(3,j) = 2*ran2(iseed) -1
            enddo
            nwork = np*(np-1d0)/2d0
            nsplit = (nwork-1)/nlocalsize + 1
            if(nsplit .le. 1) then
               njump = np/2+1
            else
               njump = (np/2-1)/nsplit + 1
            endif
            do j = np, np/2, -njump
               jstart = max(np/2,j-njump)
               jfinal = j
20             continue
                  call MPI_PROBE(MPI_ANY_SOURCE, iREADY, MPI_COMM_WORLD,
     &                       mstatus, ierror)
                  call MPI_RECV(ready, 1, MPI_INTEGER, 
     &                       mstatus(MPI_SOURCE),iREADY,MPI_COMM_WORLD,
     &                       mstatus, ierror)
                  if(ready .eq. iREADY) then
                     nsend = nsend + 1
                     call MPI_SEND(np,1,MPI_INTEGER, mstatus(MPI_SOURCE)
     &                    , NP_TAG,MPI_COMM_WORLD,IERROR)
                     call MPI_SEND(r,np*3,MPI_REAL, mstatus(MPI_SOURCE),
     &                    iR_TAG,MPI_COMM_WORLD,IERROR)
                     call MPI_SEND(jstart,1,MPI_INTEGER, 
     &                    mstatus(MPI_SOURCE), iR_TAG,
     &                    MPI_COMM_WORLD,IERROR)
                     call MPI_SEND(jfinal,1,MPI_INTEGER, 
     &                 mstatus(MPI_SOURCE),iR_TAG,MPI_COMM_WORLD,IERROR)
                  else
                     nrecv = nrecv + 1
                     call MPI_RECV(potent, 1, MPI_DOUBLE_PRECISION, 
     &                 mstatus(MPI_SOURCE), iPOT_TAG,MPI_COMM_WORLD,
     &                 mstatus, IERROR)
                     totpotent = totpotent + potent
                  endif
               if(ready .ne. iREADY) goto 20
            enddo
         enddo
         j = 0
         i = nrecv
         do while(i .lt. nsend)
            call MPI_PROBE(MPI_ANY_SOURCE, iREADY, MPI_COMM_WORLD, 
     &                     mstatus, IERROR)
            call MPI_RECV(ready, 1, MPI_INTEGER, mstatus(MPI_SOURCE), 
     &                     iREADY, MPI_COMM_WORLD, mstatus, IERROR)
            if(ready .eq. iWRITING) then
               call MPI_RECV(potent, 1, MPI_DOUBLE_PRECISION, 
     &                  mstatus(MPI_SOURCE), iPOT_TAG, MPI_COMM_WORLD,
     &                    mstatus, IERROR)
               totpotent = totpotent + potent
               call MPI_SEND(finish, 1, MPI_INTEGER, 
     &              mstatus(MPI_SOURCE), NP_TAG, MPI_COMM_WORLD, IERROR)
               i = i + 1
            else
               call MPI_SEND(finish, 1, MPI_INTEGER, mstatus(MPI_SOURCE)
     &               , NP_TAG, MPI_COMM_WORLD, IERROR)
               j = j + 1
            endif
         enddo
         do i = 1, nid-j-1
            call MPI_RECV(ready, 1, MPI_INTEGER, MPI_ANY_SOURCE, iREADY,
     &               MPI_COMM_WORLD, mstatus, IERROR)
            call MPI_SEND(finish, 1, MPI_INTEGER, mstatus(MPI_SOURCE), 
     &           NP_TAG, MPI_COMM_WORLD, IERROR)
         enddo
      else
         ready = iREADY
         call MPI_SEND(ready, 1, MPI_INTEGER, 0, iREADY, MPI_COMM_WORLD,
     &         IERROR)
         call MPI_RECV(np, 1, MPI_INTEGER, 0, NP_TAG, MPI_COMM_WORLD, 
     &         mstatus, IERROR)
         do while (np .ne. iFINISH)
            call MPI_RECV(r, np*3, MPI_REAL, 0, iR_TAG, MPI_COMM_WORLD, 
     &         mstatus, IERROR)
            call MPI_RECV(jstart, 1, MPI_INTEGER, 0, iR_TAG, 
     &           MPI_COMM_WORLD, mstatus, IERROR)
            call MPI_RECV(jfinal, 1, MPI_INTEGER, 0, iR_TAG, 
     &               MPI_COMM_WORLD, mstatus, IERROR)
            potent = potential(r, np, jstart, jfinal)
            ready = iWRITING
            call MPI_ISSEND(ready,1,MPI_INTEGER,0,iREADY,MPI_COMM_WORLD,
     &                    request, IERROR)
            call MPI_WAIT(request, mstatus,IERROR)
            call MPI_SEND(potent, 1, MPI_DOUBLE_PRECISION, 0, iPOT_TAG, 
     &              MPI_COMM_WORLD, IERROR)

            ready = iREADY
            call MPI_SEND(ready, 1, MPI_INTEGER, 0, iREADY, 
     &            MPI_COMM_WORLD,IERROR)
            call MPI_RECV(np, 1, MPI_INTEGER, 0, NP_TAG, MPI_COMM_WORLD,
     &            mstatus, IERROR)
         enddo
      endif

      time2 = etime(timearray)
      if(myid.eq.0) print *, '[F] np=',nid,'TOP is', totpotent, 'in ', 
     &            (time2-time1),' seconds'
      call MPI_FINALIZE(IERROR)

	  stop
	  end

      FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
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
C  (C) Copr. 1986-92 Numerical Recipes Software +)-*1a311.
