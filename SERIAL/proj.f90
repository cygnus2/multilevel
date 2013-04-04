      subroutine proj
                                                                       
!     Re-unitarise SU(3) matrices (using the 1st two columns).          
                                                                       
      use lattsize
      use precdef
      use filecom
      use hitcom

      implicit none

      integer i
      real (kind=double) :: wn
      complex (kind=double) :: scp 

      do i = 1,nsite
!
      wn=dsqrt((zabs(hit(1,1,i)))**2+(zabs(hit(1,2,i)))**2 &
       +(zabs(hit(1,3,i)))**2)
         hit(1,1,i)=hit(1,1,i)/wn
         hit(1,2,i)=hit(1,2,i)/wn
         hit(1,3,i)=hit(1,3,i)/wn
!
      scp=dconjg(hit(1,1,i))*hit(2,1,i)+dconjg(hit(1,2,i)) &
       *hit(2,2,i)+dconjg(hit(1,3,i))*hit(2,3,i)
         hit(2,1,i)=hit(2,1,i)-scp*hit(1,1,i)
         hit(2,2,i)=hit(2,2,i)-scp*hit(1,2,i)
         hit(2,3,i)=hit(2,3,i)-scp*hit(1,3,i)
!
      wn=dsqrt(zabs(hit(2,1,i))**2+zabs(hit(2,2,i))**2 &
       +zabs(hit(2,3,i))**2)
         hit(2,1,i)=hit(2,1,i)/wn
         hit(2,2,i)=hit(2,2,i)/wn
         hit(2,3,i)=hit(2,3,i)/wn
!
      hit(3,1,i)=dconjg(hit(1,2,i)*hit(2,3,i) &
       -hit(1,3,i)*hit(2,2,i))
      hit(3,2,i)=dconjg(hit(1,3,i)*hit(2,1,i) &
       -hit(1,1,i)*hit(2,3,i))
      hit(3,3,i)=dconjg(hit(1,1,j,i)*hit(2,2,i) &
       -hit(1,2,i)*hit(2,1,i))
!
      enddo
      return
      end
