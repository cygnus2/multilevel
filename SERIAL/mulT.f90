        SUBROUTINE mulT2(T,TT)

        use precdef
        use lattsize
!
        implicit none
!
        complex(kind=double),dimension(3,3,3,3) :: T1,T2,T3
        complex(kind=double),dimension(3,3,3,3,ncu,leng4,ir) :: T
        complex(kind=double),dimension(3,3,3,3,ncu,leng4/2,ir) :: TT

        integer :: i3,i4,j,k1,k2,k3,k4,l1,l2
!
        do j=1,ir
        do i4=1,leng4,2
        do i3=1,ncu
!
        do k1=1,3
        do k2=1,3
        do k3=1,3
        do k4=1,3
          T1(k1,k2,k3,k4)=T(k1,k2,k3,k4,i3,i4,j)
          T2(k1,k2,k3,k4)=T(k1,k2,k3,k4,i3,i4+1,j)
        end do
        end do
        end do
        end do
!
        do k1=1,3
        do k2=1,3
        do k3=1,3
        do k4=1,3
           T3(k1,k2,k3,k4)=cmplx(0d0,0d0)
           do L1=1,3
           do L2=1,3
              T3(k1,k2,k3,k4)=T3(k1,k2,k3,k4)+T1(k1,l1,k3,l2)*T2(l1,k2,l2,k4)
           end do
           end do
           TT(k1,k2,k3,k4,i3,(i4+1)/2,j)=T3(k1,k2,k3,k4)
        end do
        end do
        end do
        end do
        end do
        end do
        end do
        return
        stop
        end
