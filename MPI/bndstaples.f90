        subroutine boundstaples(istart,base)
!
!   Calculates the sum of staples for the time-like link (site) 
!
      use lattsize
      use precdef
      use filecom
      use stapcom
      use umastcom
      use neibcom
      use paramod

      implicit none

      integer :: i,j,k,l1,l2,n1,n2,n3,m1,m2,m3,istart,iti
      integer :: icount,base
      complex (kind=double), dimension(3,3) :: un1,un2,um1,um2
!
      icount=0
      istart=istart*ncu
      do i=istart+1,istart+ncu
      icount=icount+1
      do j=1,3
!
! negative t neighbours
!
         m1=neib(-4,i)
         m1=m1-base
         m2=neib(j,m1)
!         write(10+rank,*) base,m1,m2
!
!  Extract links and multiply
!
         call udagextract(um1,4,m2)
         call udagextract(um2,j,m1)
!
         call umult(um1,um2)
!
         call uextract(um2,4,m1)
         call umult(um1,um2)
!
!  create staples
!
         bndstp(:,:,j,icount)=um1(:,:)
!
         enddo
         enddo
!
         return
         end         
