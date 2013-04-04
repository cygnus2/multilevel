       subroutine subover(iti,base)    
! 
!  updates the whole lattice using 3 su(2) subgroups
! 
      use lattsize
      use precdef
      use filecom
      use umastcom
      use neibcom

      implicit none

      integer :: i,iti,j,k,l1,l2,n1,n2,n3,m1,m2,m3,base
      integer :: maxlimit

      complex (kind=double), dimension(3,3) :: ua,un,u1,u2,u3,w,smsta
      complex (kind=double), dimension(3,3) :: un1,un2,um1,um2

      real (kind=double), dimension(4) :: vv,ww
! 
!
! set maxlimit
        maxlimit=ncu*(tsloc+1)
        do i=iti*ncu+1,(iti+itdim-1)*ncu
        do j=1,4
!      
!     initialize sum of staples to zero
!
      do l1=1,3
      do l2=1,3
      smsta(l1,l2)=dcmplx(0d0,0d0)
      enddo
      enddo
!
      do k=1,4
         if (k.eq.j) go to 10
!
! positive neighbours
!
         n1=neib(j,i)
         if(n1 > maxlimit) n1=n1-base
         n2=neib(k,i)
         if(n2 > maxlimit) n2=n2-base
         n3=i-base
!
! negative neighbours
!
         m3=neib(-k,i)
         if(m3 > maxlimit) m3=m3-base
         m1=neib(j,m3)
!         m3=m3-base
         m2=m3
!         m1=m1-base
!         write(10,*) m1,m2,m3
!
!  Extract links and multiply
!
         call uextract(un1,k,n1)
         call udagextract(un2,j,n2)
!
         call umult(un1,un2)
!
         call udagextract(un2,k,n3)
         call umult(un1,un2)
!
!
         call udagextract(um1,k,m1)
         call udagextract(um2,j,m2)
!
         call umult(um1,um2)
!
         call uextract(um2,k,m3)
         call umult(um1,um2)
!
!  create sum of staples
!
         smsta(1,1)=smsta(1,1)+un1(1,1)+um1(1,1)
         smsta(1,2)=smsta(1,2)+un1(1,2)+um1(1,2)
         smsta(1,3)=smsta(1,3)+un1(1,3)+um1(1,3)
         smsta(2,1)=smsta(2,1)+un1(2,1)+um1(2,1)
         smsta(2,2)=smsta(2,2)+un1(2,2)+um1(2,2)
         smsta(2,3)=smsta(2,3)+un1(2,3)+um1(2,3)
         smsta(3,1)=smsta(3,1)+un1(3,1)+um1(3,1)
         smsta(3,2)=smsta(3,2)+un1(3,2)+um1(3,2)
         smsta(3,3)=smsta(3,3)+un1(3,3)+um1(3,3)
!
10       continue
         enddo
!
! Cabibbo-Marinari begins. smsta contains sum of staples.
!
!      link to be updated
!
        call uextract(ua,j,n3)
!
        call mult(ua,smsta,w)
!
!	su2 sub-group (1,2)
!
        vv(1)=REAL(w(1,1)+w(2,2))
        vv(2)=-IMAG(w(1,2)+w(2,1))
        vv(3)=-REAL(w(1,2)-w(2,1))
        vv(4)=-IMAG(w(1,1)-w(2,2))
!
        call su2relax(vv,ww)
!
        u1(1,1)=CMPLX(ww(1),ww(4))
        u1(1,2)=CMPLX(ww(3),ww(2))
        u1(1,3)=CMPLX(0,0)
        u1(2,1)=CMPLX(-ww(3),ww(2))
        u1(2,2)=CMPLX(ww(1),-ww(4))
        u1(2,3)=CMPLX(0,0)
        u1(3,1)=CMPLX(0,0)
        u1(3,2)=CMPLX(0,0)
        u1(3,3)=CMPLX(1,0)
!
        call mult(u1,ua,un)
        call mult(un,smsta,w)
!
!       su2 sub-group (1,3)
!
        vv(1)=REAL(w(1,1)+w(3,3))
        vv(2)=-IMAG(w(1,3)+w(3,1))
        vv(3)=-REAL(w(1,3)-w(3,1))
        vv(4)=-IMAG(w(1,1)-w(3,3))
!
        call su2relax(vv,ww)
!
        u3(1,1)=CMPLX(ww(1),ww(4))
        u3(1,3)=CMPLX(ww(3),ww(2))
        u3(1,2)=CMPLX(0,0)
        u3(3,1)=CMPLX(-ww(3),ww(2))
        u3(3,3)=CMPLX(ww(1),-ww(4))
        u3(3,2)=CMPLX(0,0)
        u3(2,1)=CMPLX(0,0)
        u3(2,3)=CMPLX(0,0)
        u3(2,2)=CMPLX(1,0)
!
        call mult(u3,un,u1)
        call mult(u1,smsta,w)
!
!       su2 sub-group (2,3)
!
        vv(1)=REAL(w(2,2)+w(3,3))
        vv(2)=-IMAG(w(2,3)+w(3,2))
        vv(3)=-REAL(w(2,3)-w(3,2))
        vv(4)=-IMAG(w(2,2)-w(3,3))
!
        call su2relax(vv,ww)
!
        u2(2,2)=CMPLX(ww(1),ww(4))
        u2(2,3)=CMPLX(ww(3),ww(2))
        u2(2,1)=CMPLX(0,0)
        u2(3,2)=CMPLX(-ww(3),ww(2))
        u2(3,3)=CMPLX(ww(1),-ww(4))
        u2(3,1)=CMPLX(0,0)
        u2(1,2)=CMPLX(0,0)
        u2(1,3)=CMPLX(0,0)
        u2(1,1)=CMPLX(1,0)
!
        call mult(u2,u1,un)
!
        umastloc(:,:,j,i-base)=un(:,:)
!
!	End of link update
!
        enddo
        enddo
!
        return
        end