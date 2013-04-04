         subroutine makePOL
!
      use lattsize
      use precdef
      use umastcom
      use neibcom
      use polcom
      use hitcom
!
      implicit none
!
      complex (kind=double), dimension(3,3) :: usm1,usm2,usm3
      integer :: i,i1,i2,i3,i4,k,k1,k2,k3,k4,j,l,lx,ly,lz,l1,l2,m,kk
!
!
        m=0
        do i4=1,leng4,itdim
        m=m+1
        k=1
        do i3=1,leng3
        do i2=1,leng2
        do i1=1,leng1
!
        l=(i4-1)*ncu+(i3-1)*ns+(i2-1)*ni+i1

        usm3=0
        usm3(1,1)=cmplx(1d0,0d0)
        usm3(2,2)=cmplx(1d0,0d0)
        usm3(3,3)=cmplx(1d0,0d0)

        do kk=1,itdim
        usm1(:,:)=hit(:,:,l)
!         call uextract(usm1,4,l)
!
         l=neib(4,l)
         usm3=matmul(usm3,usm1)
        enddo
!
        POL(:,:,k,m)=POL(:,:,k,m)+usm3(:,:)
!
        k=k+1
!
        end do
        end do
        end do
        end do
!
        return
        stop
        end
