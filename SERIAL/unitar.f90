      subroutine unitar
                                                                       
!     Re-unitarise SU(3) matrices (using the 1st two columns).          
                                                                       
      use lattsize
      use precdef
      use filecom
      use umastcom

      implicit none

      integer i,j
      real (kind=double) :: wn
      complex (kind=double) :: scp 

      do i = 1,nsite
      do j = 1,4        

!
      wn=dsqrt((zabs(umast(1,1,j,i)))**2+(zabs(umast(1,2,j,i)))**2 &
       +(zabs(umast(1,3,j,i)))**2)
         umast(1,1,j,i)=umast(1,1,j,i)/wn
         umast(1,2,j,i)=umast(1,2,j,i)/wn
         umast(1,3,j,i)=umast(1,3,j,i)/wn
!
      scp=dconjg(umast(1,1,j,i))*umast(2,1,j,i)+dconjg(umast(1,2,j,i)) &
       *umast(2,2,j,i)+dconjg(umast(1,3,j,i))*umast(2,3,j,i)
         umast(2,1,j,i)=umast(2,1,j,i)-scp*umast(1,1,j,i)
         umast(2,2,j,i)=umast(2,2,j,i)-scp*umast(1,2,j,i)
         umast(2,3,j,i)=umast(2,3,j,i)-scp*umast(1,3,j,i)
!
      wn=dsqrt(zabs(umast(2,1,j,i))**2+zabs(umast(2,2,j,i))**2 &
       +zabs(umast(2,3,j,i))**2)
         umast(2,1,j,i)=umast(2,1,j,i)/wn
         umast(2,2,j,i)=umast(2,2,j,i)/wn
         umast(2,3,j,i)=umast(2,3,j,i)/wn
!
      umast(3,1,j,i)=dconjg(umast(1,2,j,i)*umast(2,3,j,i) &
       -umast(1,3,j,i)*umast(2,2,j,i))
      umast(3,2,j,i)=dconjg(umast(1,3,j,i)*umast(2,1,j,i) &
       -umast(1,1,j,i)*umast(2,3,j,i))
      umast(3,3,j,i)=dconjg(umast(1,1,j,i)*umast(2,2,j,i) &
       -umast(1,2,j,i)*umast(2,1,j,i))
!
      enddo
      enddo
      return
      end
