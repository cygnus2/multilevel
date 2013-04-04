        subroutine boundstaples(iti)
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

      integer :: i,j,k,l1,l2,n1,n2,n3,m1,m2,m3,iti,istart,iend
      integer :: icount,tcount
      integer :: tsloc
      complex (kind=double), dimension(3,3) :: un1,un2,um1,um2
!
!      bndstp=0 ! initializes staples to zero
!
      icount=0
      tsloc=2
      tcount=iti/itdim + 1
!
      istart=(iti+tsloc-2)*ncu
      iend=istart+ncu
      do i=istart+1,iend
      icount = icount + 1
      do j=1,3
!
      k=4
!
! Calculate -t dir staples
!
         m3=neib(j,i)
!
!  Extract links and multiply
!
         call udagextract(um1,4,m3)
         call udagextract(um2,j,i)
!
         call umult(um1,um2)
!
         call uextract(um2,4,i)
         call umult(um1,um2)
!
!  create staples
!
         bndstp(:,:,j,icount,tcount)=um1(:,:)
!
!10       continue
         enddo
         enddo
!
         return
         end         
