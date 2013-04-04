       subroutine bs_update(beta,iti,base)
! 
!  updates the boundary lattice using 3 su(2) subgroups
! 
      use lattsize
      use precdef
      use filecom
      use umastcom
      use neibcom
      use stapcom
      use paramod

      implicit none

      integer :: i,iti,j,k,l1,l2,n1,n2,n3,m1,m2,m3,icount,base

      complex (kind=double), dimension(3,3) :: ua,un,u1,u2,u3,w,smsta
      complex (kind=double), dimension(3,3) :: un1,un2,um1,um2

      real (kind=double), dimension(4) :: vv,ww
      real :: beta,sbeta
      real :: plaq
      integer :: outf
! 
!
        icount=0
        do i=(iti-1)*ncu+1,iti*ncu
        icount=icount+1
        do j=1,3
!      
!     initialize sum of staples 
!
      smsta(:,:)=bndstp_r(:,:,j,icount)
!
      do k=1,3 ! both fwd and bwd staples need to be calculated in x,y,z
         if (k.eq.j) go to 10
!
! positive neighbours
!
         n1=neib(j,i)-base
         n2=neib(k,i)-base
         n3=i-base
!
! negative neighbours
!
         m3=neib(-k,i)-base
         m1=neib(j,m3)
         m2=m3
!         write(10+rank,*) base,"fwd to ",i,":",n1,n2,n3
!         write(10+rank,*) base,"bwd to ",i,":",m1,m2,m3
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
         smsta(:,:)=smsta(:,:)+un1(:,:)+um1(:,:)
!
10       continue
         enddo
! positive t-staples   
         n1=neib(4,i)-base
         n2=neib(j,i)-base
         n3=i-base
!
!  Extract links and multiply
!
         call uextract(un1,4,n2)
         call udagextract(un2,j,n1)
!
         call umult(un1,un2)
!
         call udagextract(un2,4,n3)
         call umult(un1,un2)
!
!  create sum of staples
!
         smsta(:,:)=smsta(:,:)+un1(:,:)
!
! Cabibbo-Marinari begins. smsta contains sum of staples.
!
        sbeta=beta/3.0
!
!      link to be updated
!
        call uextract(ua,j,n3)
!
        call mult(ua,smsta,w)
!       call umult(un,smsta)
        plaq = plaq + real(w(1,1)+w(2,2)+w(3,3))
!
!	su2 sub-group (1,2)
!
        vv(1)=REAL(w(1,1)+w(2,2))
        vv(2)=-IMAG(w(1,2)+w(2,1))
        vv(3)=-REAL(w(1,2)-w(2,1))
        vv(4)=-IMAG(w(1,1)-w(2,2))
!
        call su2(vv,sbeta,ww)
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
        call su2(vv,sbeta,ww)
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
        call su2(vv,sbeta,ww)
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
        !call umult(un,smsta)
        !plaq = plaq + real(un(1,1)+un(2,2)+un(3,3))
!
!	End of link update
!
        enddo
        enddo
        outf=base/ncu
        plaq=plaq/(36.d0*ncu)
        write(32+outf,*) plaq
!
        return
        end
