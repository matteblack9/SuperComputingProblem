c234567
      function potential(x,y,z, np) 
      real*8 potential
      real*4 x(np),y(np),z(np),dist,tmpx,tmpy,tmpz
      potential = 0
      do j = 1, np
      do i = j+1, np
         tmpx = x(i)-x(j)
         tmpy = y(i)-y(j)
         tmpz = z(i)-z(j)
         dist = sqrt(tmpx*tmpx + tmpy*tmpy + tmpz*tmpz)
		 if(dist .gt.0) potential = potential + 1./dist
      enddo
      enddo
      return 
      end

      function getrandomnp(istep, niter)
      integer getrandomnp,istep,niter
      if(istep .lt. niter/10) then
         getrandomnp = (3000+5000*rand())
      else
         getrandomnp = (100*rand())
      endif
      return
      end


      program main
	  parameter(maxnp= 5000000)
      integer i,j,np,niter
      real*8 totpotent,potential
      integer iseed
      real x(maxnp),y(maxnp),z(maxnp)
      real timearray(2),time1,time2,walltime
	  integer getrandomnp
	  external getrandomnp



	  totpotent = 0
	  iseed = 100
      niter = 4000
      call srand(iseed)


      call etime(timearray, time1)
      do i = 1, niter
         np = getrandomnp(i,niter)
         do j = 1, np
            x(j) = 2.*rand()-1.
            y(j) = 2.*rand()-1.
            z(j) = 2.*rand()-1.
         enddo
         totpotent = totpotent + potential(x,y,z,np)
      enddo
      call etime(timearray, time2)
      walltime = time2-time1
      print *,'total potential =',totpotent,' and wallclock time=',
     &          walltime
      stop
      end
