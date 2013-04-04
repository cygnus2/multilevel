         subroutine makeTP
!
      use lattsize
      use precdef
      use umastcom
      use neibcom
      use prodcom
      use hitcom
!
      implicit none
!
      complex (kind=double), dimension(3,3) :: usm1,usm2,usm3,usm4,usm5,usm
      integer :: i,i1,i2,i3,k,j,jj,l,l1,l2,ll2
      integer, dimension(leng4/2) :: flag
!
        do jj=1,dir-1
        do i1=1,ncu
!
!  Calculation in the first sublattice
!
        
!  Make the first bend: common to all ops
        l=i1
        l1=neib(-jj,l)
!        write(*,*) "Starting at ",i1," Direction ",jj
!        write(*,*) l,l1
        call udagextract(usm1,jj,l1)
        call uextract(usm2,4,l1)
        usm1=matmul(usm1,usm2)
        l2=neib(4,l1)
!        write(*,*) l2
        call uextract(usm3,4,l)
        call udagextract(usm4,jj,l2)
        usm3=matmul(usm3,usm4)
        usm(:,:)=usm1(:,:)-usm3(:,:)
        ll2=neib(4,i1)

!  Complete the calculation for the remaining sublattice
        flag=0
        do i=1,leng4/2
        l2=ll2

        usm3=0
        usm3(1,1)=cmplx(1d0,0d0)
        usm3(2,2)=cmplx(1d0,0d0)
        usm3(3,3)=cmplx(1d0,0d0)

        do k=2,itdim
           if(i.EQ.(k-1))then
            flag(i)=1
            l=neib(-jj,l2)
!            write(*,*) l
            call uextract(usm1,jj,l)
            call uextract(usm2,4,l2)
            usm1=matmul(usm1,usm2)
            call uextract(usm4,4,l)
            l=neib(4,l)
!            write(*,*) l
            call uextract(usm5,jj,l)
            usm4=matmul(usm4,usm5)
            usm1(:,:)=usm1(:,:)-usm4(:,:)
!            write(*,*) l2
           else if(flag(i).EQ.1)then
!             call uextract(usm1,4,l2)
              usm1(:,:)=hit(:,:,l2)
            else
             l1=neib(-jj,l2)
!             write(*,*) l1
!             call uextract(usm1,4,l1)
              usm1(:,:)=hit(:,:,l1)
           endif
           usm3=matmul(usm3,usm1)
           l2=neib(4,l2)
!           write(*,*) l2
         enddo

         usm1=matmul(usm,usm3)
         TP(:,:,i1,jj,1,i)=TP(:,:,i1,jj,1,i)+usm1(:,:)
!
         enddo
         ll2=l2
!
!  Now for the other sublattices; all ops
!
        do i=1,leng4/2
        l2=ll2
        do j=itdim,leng4-1,itdim

        usm3=0
        usm3(1,1)=cmplx(1d0,0d0)
        usm3(2,2)=cmplx(1d0,0d0)
        usm3(3,3)=cmplx(1d0,0d0)

        do k=1,itdim
           if(i.EQ.(j-1+k))then
            flag(i)=1
            l=neib(-jj,l2)
!            write(*,*) l
            call uextract(usm1,jj,l)
            call uextract(usm2,4,l2)
            usm1=matmul(usm1,usm2)
            call uextract(usm4,4,l)
            l=neib(4,l)
!            write(*,*) l
            call uextract(usm5,jj,l)
            usm4=matmul(usm4,usm5)
            usm1(:,:)=usm1(:,:)-usm4(:,:)
!            write(*,*) l2
           else if(flag(i).EQ.1)then
!             call uextract(usm1,4,l2)
              usm1(:,:)=hit(:,:,l2)
            else
             l1=neib(-jj,l2)
!             write(*,*) l1
!             call uextract(usm1,4,l1)
              usm1(:,:)=hit(:,:,l1)
           endif
           usm3=matmul(usm3,usm1)
           l2=neib(4,l2)
!           write(*,*) l2
         enddo

         TP(:,:,i1,jj,j/itdim+1,i)=TP(:,:,i1,jj,j/itdim+1,i)+usm3(:,:)

         enddo
!         write(*,*) "........."
         enddo
!         write(*,*) "---------------------"
!
        end do
        end do
         
!
        return
        stop
        end

