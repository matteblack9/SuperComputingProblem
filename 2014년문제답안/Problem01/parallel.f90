!234567890
       program abcd
       implicit none
       include 'mpif.h'
       integer m,n,ns1,ns2
       real*8, allocatable :: uu(:,:),vv(:,:),ww(:,:)
       real*8, allocatable :: xg(:),yg(:),uurow(:)
       integer jsta,jend,jsta2,jend1,iprev,inext,i,j,isend1,isend2,irecv1,irecv2
       integer miter,iter,iplot
       integer ierr,myrank,nprocs
       integer istatus(MPI_STATUS_SIZE)
       real*8 pi,hh,hh2,test0,test1
       real*8 wtime1,wtime2
       logical lexit

       call MPI_INIT(ierr)
       call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
       call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
       wtime1=MPI_WTIME()
       
       m=1000
       n=m
       allocate(xg(m),yg(n))
       do i=1,m
       xg(i)=0.d0+(1.d0-0.d0)/float(m-1)*float(i-1)
       enddo
       do j=1,n
       yg(j)=0.d0+(1.d0-0.d0)/float(n-1)*float(j-1)
       enddo
       hh=xg(2)-xg(1)
       hh2=hh*hh
       pi=4.d0*atan(1.d0)

       call para_range(1,n,nprocs,myrank,jsta,jend)
       jsta2=jsta
       jend1=jend
       if(myrank == 0) jsta2=2
       if(myrank == nprocs-1) jend1=n-1

       ns1=max(1,jsta-1)
       ns2=min(n,jend+1)
       allocate(uu(1:m,ns1:ns2))
       allocate(vv(1:m,ns1:ns2))
       allocate(ww(1:m,ns1:ns2))
       inext=myrank+1
       iprev=myrank-1
       if(myrank == nprocs-1) inext=MPI_PROC_NULL
       if(myrank == 0) iprev=MPI_PROC_NULL

       do j=jsta,jend
       do i=1,m
       uu(i,j)=0.d0
       enddo
       enddo
       call vnbd(m,n,ns1,ns2,uu)
       vv=uu
       miter=1000000000
       do iter=1,miter
       call vnbd(m,n,ns1,ns2,uu)

       call MPI_ISEND(uu(1,jend),   m,MPI_DOUBLE_PRECISION,inext,1,MPI_COMM_WORLD,isend1,ierr)
       call MPI_ISEND(uu(1,jsta),   m,MPI_DOUBLE_PRECISION,iprev,1,MPI_COMM_WORLD,isend2,ierr)
       call MPI_IRECV(uu(1,jsta-1), m,MPI_DOUBLE_PRECISION,iprev,1,MPI_COMM_WORLD,irecv1,ierr)
       call MPI_IRECV(uu(1,jend+1), m,MPI_DOUBLE_PRECISION,inext,1,MPI_COMM_WORLD,irecv2,ierr)
       call MPI_WAIT(isend1,istatus,ierr)
       call MPI_WAIT(isend2,istatus,ierr)
       call MPI_WAIT(irecv1,istatus,ierr)
       call MPI_WAIT(irecv2,istatus,ierr)

       call axmb(m,n,ns1,ns2,uu,ww,xg,yg,jsta2,jend1)
       call vnbd(m,n,ns1,ns2,ww)

       do j=jsta,jend
       do i=1,m
       uu(i,j)=vv(i,j)+ww(i,j)*(hh2/4.d0) *0.9
       enddo
       enddo
       test1=maxval(abs(vv-uu))/(maxval(abs(uu))+1.d-8)
       vv=uu
       call MPI_REDUCE(test1,test0,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
       if(myrank ==0)then
       if(mod(iter,1000) == 1 .or. iter < 1000) print*, iter,test0
       lexit=.false.
       if(test0 < 1.d-5) lexit=.true.
                    endif
       call MPI_BCAST(lexit,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
       if(lexit) exit
       enddo

       if(myrank == 0)then
       print*, uu(1,1),uu(m,1)
                      endif
       if(myrank == nprocs-1)then
       print*, uu(1,n),uu(m,n)
                             endif
!      set pm3d
!      set palette rgbformulae 33,13,10
!      splot 'fort.11' with pm3d
!
       open(13,file='sol',form='unformatted',access='direct',recl=8*(2+m))
       do j=jsta,jend
       write(13,rec=j) (uu(i,j),i=1,m)
       enddo
       close(13)
       wtime2=MPI_WTIME()
       if(myrank ==0) print*,wtime2-wtime1,' (sec)'
       iplot=0
       iplot=1
       if(iplot ==1)then
       if(myrank == 0)then
       allocate(uurow(m))
       open(11,file='fort.11',form='formatted')
       open(13,file='sol',form='unformatted',access='direct',recl=8*(2+m))
       do j=1,n
       read(13,rec=j) (uurow(i),i=1,m)
       do i=1,m
       write(11,*) xg(i),yg(j),uurow(i)
       enddo
       write(11,*)
       enddo
       close(13)
       close(11)
       deallocate(uurow)
                      endif
                    endif
      
       deallocate(xg,yg)
       deallocate(uu,vv,ww)
       call MPI_FINALIZE(ierr)
       end program abcd
!
       subroutine vnbd(m,n,ns1,ns2,uu)
       implicit none
       integer m,n,ns1,ns2
       real*8 uu(m,ns1:ns2)

       uu(1,:)=uu(2,:)
       uu(m,:)=uu(m-1,:)
       if(ns1 == 1) uu(:,1)=uu(:,2)
       if(ns2 == n) uu(:,n)=uu(:,n-1)
       end subroutine vnbd

       subroutine axmb(m,n,ns1,ns2,uu,ww,xg,yg,jsta2,jend1)
       implicit none
       integer m,n,ns1,ns2,jsta2,jend1
       real*8 uu(m,ns1:ns2),ww(m,ns1:ns2),xg(m),yg(n)
       integer i,j
       real*8 pi,hh,hh2

       hh=xg(2)-xg(1)
       hh2=hh*hh
       ww=0.d0
       do i=2,m-1
       ww(i,:)=ww(i,:)+(-2.d0*uu(i,:)+uu(i-1,:)+uu(i+1,:))/hh2
       enddo
       hh=yg(2)-yg(1)
       hh2=hh*hh
       do j=ns1+1,ns2-1
       ww(:,j)=ww(:,j)+(-2.d0*uu(:,j)+uu(:,j-1)+uu(:,j+1))/hh2
       enddo
       pi=4.d0*atan(1.d0)
       do i=1,m
       do j=ns1,ns2
       ww(i,j)=ww(i,j)-(-2.d0*pi*pi*cos(pi*xg(i))*cos(pi*yg(j)))
       enddo
       enddo
       end subroutine axmb

       subroutine para_range(n1,n2,nprocs,myrank,ista,iend)
       implicit none
       integer n1,n2,nprocs,myrank,ista,iend
       integer iwork1,iwork2

       iwork1=(n2-n1+1)/nprocs
       iwork2=mod(n2-n1+1,nprocs)
       ista=myrank*iwork1+n1+min(myrank,iwork2)
       iend=ista+iwork1-1
       if(iwork2 > myrank) iend=iend+1
       return
       end
       subroutine para_range1(n1,n2,nprocs,myrank,ista,iend)
       implicit none
       integer n1,n2,nprocs,myrank,ista,iend
       integer iwork

       iwork=(n2-n1)/nprocs+1
       ista=min(myrank*iwork+n1,n2+1)
       iend=min(ista+iwork-1,n2)
       return
       end

