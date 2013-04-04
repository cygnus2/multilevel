        subroutine staple
!
!   Calculates the sum of staples for the time-like link (site) 
!
      use lattsize
      use precdef
      use filecom
      use stapcom
      use umastcom
      use neibcom

      implicit none

      integer i,j,k,l1,l2,n1,n2,n3,m1,m2,m3
      complex (kind=double), dimension(3,3) :: un1,un2,um1,um2

      staples=0 ! initializes staples to zero

      do i=1,nsite
      do j=1,4
!
      do k=1,4
         if (k.eq.j) go to 10
!
! positive neighbours
!
         n1=neib(j,i)
         n2=neib(k,i)
         n3=i
!
! negative neighbours
!
         m3=neib(-k,i)
         m2=m3
         m1=neib(j,m3)
!
!  Extract links and multiply
!
         call uextract(un1,k,n1)
         call udagextract(un2,j,n2)
!
         call umult(un1,un2)
!
         call udagextract(un2,k,n3)
         call umult(un1,un2)
!
!
         call udagextract(um1,k,m1)
         call udagextract(um2,j,m2)
!
         call umult(um1,um2)
!
         call uextract(um2,k,m3)
         call umult(um1,um2)
!
!  create staples
!
         staples(1,1,j,i)=staples(1,1,j,i)+un1(1,1)+um1(1,1)
         staples(1,2,j,i)=staples(1,2,j,i)+un1(1,2)+um1(1,2)
         staples(1,3,j,i)=staples(1,3,j,i)+un1(1,3)+um1(1,3)
         staples(2,1,j,i)=staples(2,1,j,i)+un1(2,1)+um1(2,1)
         staples(2,2,j,i)=staples(2,2,j,i)+un1(2,2)+um1(2,2)
         staples(2,3,j,i)=staples(2,3,j,i)+un1(2,3)+um1(2,3)
         staples(3,1,j,i)=staples(3,1,j,i)+un1(3,1)+um1(3,1)
         staples(3,2,j,i)=staples(3,2,j,i)+un1(3,2)+um1(3,2)
         staples(3,3,j,i)=staples(3,3,j,i)+un1(3,3)+um1(3,3)
!
10       continue
         enddo
         enddo
         enddo
!
         return
         end         
