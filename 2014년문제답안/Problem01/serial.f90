!234567890
       program two_dim
       implicit none
       integer m,n
       real*8, allocatable :: uu(:,:),vv(:,:),ww(:,:)
       real*8, allocatable :: xg(:),yg(:)
       integer i,j,miter,iter
       real*8 pi,test0,hh,hh2

       m=1000
       n=m
       allocate(xg(m),yg(n))
       allocate(uu(m,n)) ; allocate(vv(m,n),ww(m,n))
       do i=1,m
       xg(i)=0.d0+(1.d0-0.d0)/float(m-1)*float(i-1)
       enddo
       do j=1,n
       yg(j)=0.d0+(1.d0-0.d0)/float(n-1)*float(j-1)
       enddo
       pi=4.d0*datan(1.d0)
       hh=xg(2)-xg(1) ; hh2=hh*hh

       uu=0.0d0
       call vnbd(m,n,uu)
       vv=uu

       miter=1000000000
       miter=10000
       do iter=1,miter
       call vnbd(m,n,uu)
       call axmb(m,n,uu,ww,xg,yg)
       call vnbd(m,n,ww)
       do i=2,m-1  
       do j=2,n-1  
       uu(i,j)=vv(i,j)+ww(i,j)*(hh2/4.d0) *0.9
       enddo
       enddo
       test0=maxval(abs(vv-uu))/(maxval(abs(uu))+1.d-8)
       vv=uu
       if(mod(iter,1000) == 1 .or. iter < 1000) print*, iter,test0
       if(test0 < 1.d-6) exit
       enddo

       print*, uu(1,1),uu(m,1)
       print*, uu(1,n),uu(m,n)
       open(11,file='fort.11',form='formatted')
       do i=1,m
       do j=1,n
       write(11,*) xg(i),yg(j), uu(i,j)
       enddo
       write(11,*)
       enddo
       close(11)
!      open(12,file='fort.12',form='formatted')
!      do i=1,m
!      do j=1,n
!      write(12,*) xg(i),yg(j), uu(i,j)-cos(pi*xg(i))*cos(pi*yg(j))
!      enddo
!      write(12,*)
!      enddo
!      close(12)
!
!      set pm3d
!      set palette rgbformulae 33,13,10
!      splot 'fort.11' with pm3d
!
       deallocate(xg,yg)
       deallocate(uu) ; deallocate(vv,ww)
       stop
       end program two_dim

       subroutine axmb(m,n,uu,ww,xg,yg)
       implicit none
       integer m,n
       real*8 uu(m,n),ww(m,n),xg(m),yg(n)
       integer i,j
       real*8 pi,hh,hh2

       hh=xg(2)-xg(1) ; hh2=hh*hh
       ww=0.d0
       do i=2,m-1
       do j=2,n-1
       ww(i,j)=ww(i,j)+(-2.d0*uu(i,j)+uu(i-1,j)+uu(i+1,j))/hh2
       enddo
       enddo
       hh=yg(2)-yg(1) ; hh2=hh*hh
       do i=2,m-1
       do j=2,n-1
       ww(i,j)=ww(i,j)+(-2.d0*uu(i,j)+uu(i,j-1)+uu(i,j+1))/hh2
       enddo
       enddo
       pi=4.d0*datan(1.d0)
       do i=1,m
       do j=1,n
       ww(i,j)=ww(i,j)-(-2.d0*pi*pi*cos(pi*xg(i))*cos(pi*yg(j)))
       enddo
       enddo
       end subroutine axmb

       subroutine vnbd(m,n,uu)
       implicit none
       integer m,n
       real*8 uu(m,n)

       uu(1,:)=uu(2,:)
       uu(m,:)=uu(m-1,:)
       uu(:,1)=uu(:,2)
       uu(:,n)=uu(:,n-1)
       end subroutine vnbd

