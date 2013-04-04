       subroutine multihit(beta)
! 
!  updates the whole lattice using 3 su(2) subgroups
! 
      use lattsize
      use precdef
      use filecom
      use umastcom
      use neibcom
      use stapcom
      use hitcom

      implicit none

      integer :: i,j,j1,j2

      complex (kind=double), dimension(3,3) :: un,u1,smsta

      real :: beta
! 
        call init_arrays
!
        do i=1,nsite
        j=4
!      
!
!
         smsta(1,1)=staples(1,1,4,i)
         smsta(1,2)=staples(1,2,4,i)
         smsta(1,3)=staples(1,3,4,i)
         smsta(2,1)=staples(2,1,4,i)
         smsta(2,2)=staples(2,2,4,i)
         smsta(2,3)=staples(2,3,4,i)
         smsta(3,1)=staples(3,1,4,i)
         smsta(3,2)=staples(3,2,4,i)
         smsta(3,3)=staples(3,3,4,i)
!
! Cabibbo-Marinari begins. smsta contains sum of staples.
!
        do j1=1,3
        do j2=1,3
        u1(j1,j2)=conjg(smsta(j2,j1))
        enddo
        enddo
        call gauss_integ(beta,u1,un)
!
        hit(:,:,i)=un(:,:)
!
!	End of multihit    
!
        enddo
!
        return
        end
