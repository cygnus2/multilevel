       subroutine cabmar(beta)
! 
!  updates the whole lattice using 3 su(2) subgroups
! 
      use lattsize
      use precdef
      use filecom
      use umastcom
      use neibcom
      use ranlxd_generator
      use omp_lib

      implicit none

      integer :: i,j,k
      integer :: eot,eoz,eoy,eox,x,y,z,t
      integer :: l1,l2,n1,n2,n3,m1,m2,m3
      complex (kind=double), dimension(3,3) :: ua,un,u1,u2,u3,w,smsta
      complex (kind=double), dimension(3,3) :: un1,un2,um1,um2

      real (kind=double), dimension(4) :: vv,ww
      real :: beta,sbeta
!
       sbeta=beta/3.0
!
!
!   Update all links in x-dir
!
      do eot=0,1
      do eoz=0,1
      do eoy=0,1
      do t=eot+1,leng4,2
      do z=eoz+1,leng3,2
      do y=eoy+1,leng2,2
      do x=1,leng1
      
      i = x + ni*(y-1) + ns*(z-1) + ncu*(t-1)
      j = 1
! Update link in the x-dir belonging to site i
!      write(*,*) "i=",i,"j=",j
!      call subcabmar(beta,j,i)

      smsta=0
!
      do k=1,4
         if (k.eq.j) go to 10
!
! positive neighbours
!
         n1=neib(j,i)
         n2=neib(k,i)
         n3=i
!
! negative neighbours
!
         m3=neib(-k,i)
         m2=m3
         m1=neib(j,m3)
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
!
!      link to be updated
!
        call uextract(ua,j,i)
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
        umast(1,1,j,i)=un(1,1)
        umast(1,2,j,i)=un(1,2)
        umast(1,3,j,i)=un(1,3)
        umast(2,1,j,i)=un(2,1)
        umast(2,2,j,i)=un(2,2)
        umast(2,3,j,i)=un(2,3)
        umast(3,1,j,i)=un(3,1)
        umast(3,2,j,i)=un(3,2)
        umast(3,3,j,i)=un(3,3)
!
!	End of link update
!
!
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
!
!   Update all links in y-dir
!
      do eot=0,1
      do eoz=0,1
      do eox=0,1
      do t=eot+1,leng4,2
      do z=eoz+1,leng3,2
      do y=1,leng2
      do x=eox+1,leng1,2
      
      i = x + ni*(y-1) + ns*(z-1) + ncu*(t-1)
      j = 2
! Update link in the y-dir belonging to site i
!      write(*,*) "i=",i,"j=",j
!      call subcabmar(beta,j,i)

      smsta=0
!
      do k=1,4
         if (k.eq.j) go to 20
!
! positive neighbours
!
         n1=neib(j,i)
         n2=neib(k,i)
         n3=i
!
! negative neighbours
!
         m3=neib(-k,i)
         m2=m3
         m1=neib(j,m3)
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
20       continue
         enddo
!
! Cabibbo-Marinari begins. smsta contains sum of staples.
!
!
!      link to be updated
!
        call uextract(ua,j,i)
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
        umast(1,1,j,i)=un(1,1)
        umast(1,2,j,i)=un(1,2)
        umast(1,3,j,i)=un(1,3)
        umast(2,1,j,i)=un(2,1)
        umast(2,2,j,i)=un(2,2)
        umast(2,3,j,i)=un(2,3)
        umast(3,1,j,i)=un(3,1)
        umast(3,2,j,i)=un(3,2)
        umast(3,3,j,i)=un(3,3)
!
!	End of link update
!
!
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
!
!   Update all links in z-dir
!
      do eot=0,1
      do eoy=0,1
      do eox=0,1
      do t=eot+1,leng4,2
      do z=1,leng3
      do y=eoy+1,leng2,2
      do x=eox+1,leng1,2
      
      i = x + ni*(y-1) + ns*(z-1) + ncu*(t-1)
      j = 3
! Update link in the z-dir belonging to site i
!      write(*,*) "i=",i,"j=",j
!      call subcabmar(beta,j,i)

      smsta=0
!
      do k=1,4
         if (k.eq.j) go to 30
!
! positive neighbours
!
         n1=neib(j,i)
         n2=neib(k,i)
         n3=i
!
! negative neighbours
!
         m3=neib(-k,i)
         m2=m3
         m1=neib(j,m3)
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
30       continue
         enddo
!
! Cabibbo-Marinari begins. smsta contains sum of staples.
!
!
!      link to be updated
!
        call uextract(ua,j,i)
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
        umast(1,1,j,i)=un(1,1)
        umast(1,2,j,i)=un(1,2)
        umast(1,3,j,i)=un(1,3)
        umast(2,1,j,i)=un(2,1)
        umast(2,2,j,i)=un(2,2)
        umast(2,3,j,i)=un(2,3)
        umast(3,1,j,i)=un(3,1)
        umast(3,2,j,i)=un(3,2)
        umast(3,3,j,i)=un(3,3)
!
!	End of link update
!
!
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
!
!   Update all links in t-dir
!
      do eoz=0,1
      do eoy=0,1
      do eox=0,1
      do t=1,leng4
      do z=eoz+1,leng3,2
      do y=eoy+1,leng2,2
      do x=eox+1,leng1,2
      
      i = x + ni*(y-1) + ns*(z-1) + ncu*(t-1)
      j = 4
! Update link in the t-dir belonging to site i
!      write(*,*) "i=",i,"j=",j
!      call subcabmar(beta,j,i)

      smsta=0
!
      do k=1,4
         if (k.eq.j) go to 40
!
! positive neighbours
!
         n1=neib(j,i)
         n2=neib(k,i)
         n3=i
!
! negative neighbours
!
         m3=neib(-k,i)
         m2=m3
         m1=neib(j,m3)
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
40       continue
         enddo
!
! Cabibbo-Marinari begins. smsta contains sum of staples.
!
!
!      link to be updated
!
        call uextract(ua,j,i)
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
        umast(1,1,j,i)=un(1,1)
        umast(1,2,j,i)=un(1,2)
        umast(1,3,j,i)=un(1,3)
        umast(2,1,j,i)=un(2,1)
        umast(2,2,j,i)=un(2,2)
        umast(2,3,j,i)=un(2,3)
        umast(3,1,j,i)=un(3,1)
        umast(3,2,j,i)=un(3,2)
        umast(3,3,j,i)=un(3,3)
!
!	End of link update
!
!
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo

!
      return
      end
