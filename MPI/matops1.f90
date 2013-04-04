
      subroutine uextract1(uu,mu,x)

!     Extracts the (x,mu) link variable from umast.

      use lattsize
      use precdef
      use umastcom

      implicit none

      integer :: x,mu,i,j
      complex (kind=double), dimension(3,3):: uu

      do j = 1,3
      do i = 1,3

      uu(i,j) = umastloc1(i,j,mu,x)

      end do
      end do

      return
      end

!     ----------------------------------------------------------

      subroutine udagextract1(uu,mu,x)

!     Extracts the (x,mu) link variable from umast and daggers it.

      use lattsize
      use precdef
      use umastcom

      implicit none

      integer :: x,mu,i,j
      complex (kind=double), dimension(3,3):: uu

      do j = 1,3
      do i = 1,3

      uu(j,i) = conjg(umastloc1(i,j,mu,x))

      end do
      end do

      return
      end
