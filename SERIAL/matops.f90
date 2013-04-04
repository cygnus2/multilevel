
      subroutine uextract(uu,mu,x)

!     Extracts the (x,mu) link variable from umast.

      use lattsize
      use precdef
      use umastcom

      implicit none

      integer :: x,mu,i,j
      complex (kind=double), dimension(3,3):: uu

      do j = 1,3
      do i = 1,3

      uu(i,j) = umast(i,j,mu,x)

      end do
      end do

      return
      end

!     ----------------------------------------------------------

      subroutine udagextract(uu,mu,x)

!     Extracts the (x,mu) link variable from umast and daggers it.

      use lattsize
      use precdef
      use umastcom

      implicit none

      integer :: x,mu,i,j
      complex (kind=double), dimension(3,3):: uu

      do j = 1,3
      do i = 1,3

      uu(j,i) = conjg(umast(i,j,mu,x))

      end do
      end do

      return
      end

!     ----------------------------------------------------------

      subroutine umult(ua,ub)

!     ua <-- ua*ub

      use precdef

      implicit none


      integer i,j
      complex (kind=double), dimension(3,3) ::  ua,ub,uc

!    The k loop is written explicitly - uc doesn't need to
!    set to zero.

      do i = 1,3
      do j = 1,3

      uc(j,i) =   ua(j,1)*ub(1,i)&
               + ua(j,2)*ub(2,i)&
               + ua(j,3)*ub(3,i)

      end do
      end do

      do i = 1,3
      do j = 1,3

      ua(j,i) = uc(j,i)

      end do
      end do

      return
      end
!------------------------------------------------------------
        subroutine mult(A,B,C)
!
        use precdef
!
        complex (kind=double), dimension(3,3) :: A,B,C
!
        C(1,1)=A(1,1)*B(1,1)+A(1,2)*B(2,1)+A(1,3)*B(3,1)
        C(1,2)=A(1,1)*B(1,2)+A(1,2)*B(2,2)+A(1,3)*B(3,2)
        C(1,3)=A(1,1)*B(1,3)+A(1,2)*B(2,3)+A(1,3)*B(3,3)
!
        C(2,1)=A(2,1)*B(1,1)+A(2,2)*B(2,1)+A(2,3)*B(3,1)
        C(2,2)=A(2,1)*B(1,2)+A(2,2)*B(2,2)+A(2,3)*B(3,2)
        C(2,3)=A(2,1)*B(1,3)+A(2,2)*B(2,3)+A(2,3)*B(3,3)
!
        C(3,1)=A(3,1)*B(1,1)+A(3,2)*B(2,1)+A(3,3)*B(3,1)
        C(3,2)=A(3,1)*B(1,2)+A(3,2)*B(2,2)+A(3,3)*B(3,2)
        C(3,3)=A(3,1)*B(1,3)+A(3,2)*B(2,3)+A(3,3)*B(3,3)
!
        return
        end
