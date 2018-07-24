!234567890
       program problem_fd
       implicit none
       integer i,n1,n2,j,jsta,jend
       integer iter,niter
       real*8 xi,xf,dx
       real*8 tmr
       real*8, allocatable:: ar(:),br(:)
       real*8, external :: genvv

! do not change ----{
       n1=1
       n2=100
       n2=100000000
       niter=3
! do not change ----}

       allocate ( ar(n2), br(n2) )

       jsta=n1
       jend=n2
       jsta=n1+1
       jend=n2-1
       xi=0.d0
       xf=1.d0
       dx=(xf-xi)/dble(n2-n1)
       do i=n1,n2
       br(i)=xi+dble(i-n1)*dx
       enddo
       do iter=1,niter
       do j=jsta,jend
! do not change ----{
       ar(j)=(br(j-1)+br(j+1))/4.d0+br(j)/2.d0+1.d0/genvv(br(j))
! do not change ----}
       enddo
       do i=n1,n2
! do not change ----{
       br(i)=ar(i)
! do not change ----}
       enddo
       enddo
       tmr=0.d0
       do j=jsta,jend
       tmr=tmr+ar(j)
       enddo
       print *,'tmr = ',tmr
   
       deallocate( ar, br );
       stop
       end program problem_fd
!234567890
       real*8 function genvv(x)
       implicit none
       real*8 x
       
       genvv=x**2+x**4+x**6+exp(-x**2)+cos(x)+sin(x)+tan(x)
       return
       end

