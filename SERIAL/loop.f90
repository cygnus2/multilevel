	subroutine loopxy
!
! COMPUTES X-Y WILSON LOOP
!                                                                       
      use lattsize
      use precdef
      use umastcom
      use neibcom
!
      implicit none
!
      complex (kind=double), dimension(3,3) :: x,r,v,w
      complex (kind=double), dimension(3,3,nsite) :: q1
      integer :: i,jd,ll,jt,jz,jy,jx,l,ii,ij,ji,jj
      integer :: mm,nn,nar
      real :: sum,ave
!
!----------------------------------
	do 300 jd=mjd,mjd
	ll=1
!
	r(1,1)=dcmplx(1d0,0d0)
        r(1,2)=dcmplx(0d0,0d0)
        r(1,3)=dcmplx(0d0,0d0)
        r(2,1)=dcmplx(0d0,0d0)
        r(2,2)=dcmplx(1d0,0d0)
        r(2,3)=dcmplx(0d0,0d0)
        r(3,1)=dcmplx(0d0,0d0)
        r(3,2)=dcmplx(0d0,0d0)
        r(3,3)=dcmplx(1d0,0d0)	
!	
	do 175 jt=1,leng4
	do 170 jz=1,leng3
	do 160 jy=1,leng2
	do 110 jx=1,leng1
        l=(jt-1)*ncu+(jz-1)*ns+(jy-1)*ni+jx
        do 120 ii=1,jd
!
        call uextract(v,1,l)
	call mult(r,v,x)
!
	l=neib(1,l)
!
        do mm=1,3
        do nn=1,3
        r(mm,nn)=x(mm,nn)
        enddo
        enddo
!
120	continue
	do 130 ij=1,JD
!
        call uextract(v,2,l)
        call mult(r,v,x)
!
	l=neib(2,l)
!
        do mm=1,3
        do nn=1,3
        r(mm,nn)=x(mm,nn)
        enddo
        enddo
!
130	continue
	do 140 ji=1,jd
!
        l=neib(-1,l)
        call udagextract(v,1,l)
        call mult(r,v,x)
!
        do mm=1,3
        do nn=1,3
        r(mm,nn)=x(mm,nn)
        enddo
        enddo
!
140	continue
	do 150 jj=1,jd
!
        l=neib(-2,l)
        call udagextract(v,2,l)
        call mult(r,v,x)
!
        do mm=1,3
        do nn=1,3
        r(mm,nn)=x(mm,nn)
        enddo
        enddo
!
150	continue
        nar=(jt-1)*ncu+(jz-1)*ns+(jy-1)*ni+jx
!
        do mm=1,3
        do nn=1,3
        q1(mm,nn,nar)=r(mm,nn)
        enddo
        enddo
!
	ll=ll+1
!
        r(1,1)=dcmplx(1d0,0d0)
        r(1,2)=dcmplx(0d0,0d0)
        r(1,3)=dcmplx(0d0,0d0)
        r(2,1)=dcmplx(0d0,0d0)
        r(2,2)=dcmplx(1d0,0d0)
        r(2,3)=dcmplx(0d0,0d0)
        r(3,1)=dcmplx(0d0,0d0)
        r(3,2)=dcmplx(0d0,0d0)
        r(3,3)=dcmplx(1d0,0d0)
!
110	continue
160	continue
170	continue
175	continue

!	GENERATES AVERAGE OF WILSON LOOPS
!-------------------------------------------------------------
	sum=0
	do i=1,nsite
	sum=sum+real(q1(1,1,i)+q1(2,2,i)+q1(3,3,i))
	enddo
	ave=sum/(3*nsite)
	write(31,*) ave
300 	continue
	return
	end



	subroutine loopxz
!
! COMPUTES X-Z WILSON LOOP
!                                                                       
      use lattsize
      use precdef
      use umastcom
      use neibcom
!
      implicit none
!
      complex (kind=double), dimension(3,3) :: x,r,v,w
      complex (kind=double), dimension(3,3,nsite) :: q1
      integer :: i,jd,ll,jt,jz,jy,jx,l,ii,ij,ji,jj
      integer :: mm,nn,nar
      real :: sum,ave
!
!----------------------------------
	do 300 jd=mjd,mjd
	ll=1
!
	r(1,1)=dcmplx(1d0,0d0)
        r(1,2)=dcmplx(0d0,0d0)
        r(1,3)=dcmplx(0d0,0d0)
        r(2,1)=dcmplx(0d0,0d0)
        r(2,2)=dcmplx(1d0,0d0)
        r(2,3)=dcmplx(0d0,0d0)
        r(3,1)=dcmplx(0d0,0d0)
        r(3,2)=dcmplx(0d0,0d0)
        r(3,3)=dcmplx(1d0,0d0)	
!	
	do 175 jt=1,leng4
	do 170 jz=1,leng3
	do 160 jy=1,leng2
	do 110 jx=1,leng1
        l=(jt-1)*ncu+(jz-1)*ns+(jy-1)*ni+jx
        do 120 ii=1,jd
!
        call uextract(v,1,l)
	call mult(r,v,x)
!
	l=neib(1,l)
!
        do mm=1,3
        do nn=1,3
        r(mm,nn)=x(mm,nn)
        enddo
        enddo
!
120	continue
	do 130 ij=1,JD
!
        call uextract(v,3,l)
        call mult(r,v,x)
!
	l=neib(3,l)
!
        do mm=1,3
        do nn=1,3
        r(mm,nn)=x(mm,nn)
        enddo
        enddo
!
130	continue
	do 140 ji=1,jd
!
        l=neib(-1,l)
        call udagextract(v,1,l)
        call mult(r,v,x)
!
        do mm=1,3
        do nn=1,3
        r(mm,nn)=x(mm,nn)
        enddo
        enddo
!
140	continue
	do 150 jj=1,jd
!
        l=neib(-3,l)
        call udagextract(v,3,l)
        call mult(r,v,x)
!
        do mm=1,3
        do nn=1,3
        r(mm,nn)=x(mm,nn)
        enddo
        enddo
!
150	continue
        nar=(jt-1)*ncu+(jz-1)*ns+(jy-1)*ni+jx
!
        do mm=1,3
        do nn=1,3
        q1(mm,nn,nar)=r(mm,nn)
        enddo
        enddo
!
	ll=ll+1
!
        r(1,1)=dcmplx(1d0,0d0)
        r(1,2)=dcmplx(0d0,0d0)
        r(1,3)=dcmplx(0d0,0d0)
        r(2,1)=dcmplx(0d0,0d0)
        r(2,2)=dcmplx(1d0,0d0)
        r(2,3)=dcmplx(0d0,0d0)
        r(3,1)=dcmplx(0d0,0d0)
        r(3,2)=dcmplx(0d0,0d0)
        r(3,3)=dcmplx(1d0,0d0)
!
110	continue
160	continue
170	continue
175	continue

!	GENERATES AVERAGE OF WILSON LOOPS
!-------------------------------------------------------------
	sum=0
	do i=1,nsite
	sum=sum+real(q1(1,1,i)+q1(2,2,i)+q1(3,3,i))
	enddo
	ave=sum/(3*nsite)
	write(32,*) ave
300 	continue
	return
	end



	subroutine loopyz
!
! COMPUTES Y-Z WILSON LOOP
!                                                                       
      use lattsize
      use precdef
      use umastcom
      use neibcom
!
      implicit none
!
      complex (kind=double), dimension(3,3) :: x,r,v,w
      complex (kind=double), dimension(3,3,nsite) :: q1
      integer :: i,jd,ll,jt,jz,jy,jx,l,ii,ij,ji,jj
      integer :: mm,nn,nar
      real :: sum,ave
!
!----------------------------------
	do 300 jd=mjd,mjd
	ll=1
!
	r(1,1)=dcmplx(1d0,0d0)
        r(1,2)=dcmplx(0d0,0d0)
        r(1,3)=dcmplx(0d0,0d0)
        r(2,1)=dcmplx(0d0,0d0)
        r(2,2)=dcmplx(1d0,0d0)
        r(2,3)=dcmplx(0d0,0d0)
        r(3,1)=dcmplx(0d0,0d0)
        r(3,2)=dcmplx(0d0,0d0)
        r(3,3)=dcmplx(1d0,0d0)	
!	
	do 175 jt=1,leng4
	do 170 jz=1,leng3
	do 160 jy=1,leng2
	do 110 jx=1,leng1
        l=(jt-1)*ncu+(jz-1)*ns+(jy-1)*ni+jx
        do 120 ii=1,jd
!
        call uextract(v,2,l)
	call mult(r,v,x)
!
	l=neib(2,l)
!
        do mm=1,3
        do nn=1,3
        r(mm,nn)=x(mm,nn)
        enddo
        enddo
!
120	continue
	do 130 ij=1,JD
!
        call uextract(v,3,l)
        call mult(r,v,x)
!
	l=neib(3,l)
!
        do mm=1,3
        do nn=1,3
        r(mm,nn)=x(mm,nn)
        enddo
        enddo
!
130	continue
	do 140 ji=1,jd
!
        l=neib(-2,l)
        call udagextract(v,2,l)
        call mult(r,v,x)
!
        do mm=1,3
        do nn=1,3
        r(mm,nn)=x(mm,nn)
        enddo
        enddo
!
140	continue
	do 150 jj=1,jd
!
        l=neib(-3,l)
        call udagextract(v,3,l)
        call mult(r,v,x)
!
        do mm=1,3
        do nn=1,3
        r(mm,nn)=x(mm,nn)
        enddo
        enddo
!
150	continue
        nar=(jt-1)*ncu+(jz-1)*ns+(jy-1)*ni+jx
!
        do mm=1,3
        do nn=1,3
        q1(mm,nn,nar)=r(mm,nn)
        enddo
        enddo
!
	ll=ll+1
!
        r(1,1)=dcmplx(1d0,0d0)
        r(1,2)=dcmplx(0d0,0d0)
        r(1,3)=dcmplx(0d0,0d0)
        r(2,1)=dcmplx(0d0,0d0)
        r(2,2)=dcmplx(1d0,0d0)
        r(2,3)=dcmplx(0d0,0d0)
        r(3,1)=dcmplx(0d0,0d0)
        r(3,2)=dcmplx(0d0,0d0)
        r(3,3)=dcmplx(1d0,0d0)
!
110	continue
160	continue
170	continue
175	continue

!	GENERATES AVERAGE OF WILSON LOOPS
!-------------------------------------------------------------
	sum=0
	do i=1,nsite
	sum=sum+real(q1(1,1,i)+q1(2,2,i)+q1(3,3,i))
	enddo
	ave=sum/(3*nsite)
	write(33,*) ave
300 	continue
	return
	end





	subroutine loopxt
!
! COMPUTES X-T WILSON LOOP
!                                                                       
      use lattsize
      use precdef
      use umastcom
      use neibcom
!
      implicit none
!
      complex (kind=double), dimension(3,3) :: x,r,v,w
      complex (kind=double), dimension(3,3,nsite) :: q1
      integer :: i,jd,ll,jt,jz,jy,jx,l,ii,ij,ji,jj
      integer :: mm,nn,nar
      real :: sum,ave
!
!----------------------------------
	do 300 jd=mjd,mjd
	ll=1
!
	r(1,1)=dcmplx(1d0,0d0)
        r(1,2)=dcmplx(0d0,0d0)
        r(1,3)=dcmplx(0d0,0d0)
        r(2,1)=dcmplx(0d0,0d0)
        r(2,2)=dcmplx(1d0,0d0)
        r(2,3)=dcmplx(0d0,0d0)
        r(3,1)=dcmplx(0d0,0d0)
        r(3,2)=dcmplx(0d0,0d0)
        r(3,3)=dcmplx(1d0,0d0)	
!	
	do 175 jt=1,leng4
	do 170 jz=1,leng3
	do 160 jy=1,leng2
	do 110 jx=1,leng1
        l=(jt-1)*ncu+(jz-1)*ns+(jy-1)*ni+jx
        do 120 ii=1,jd
!
        call uextract(v,1,l)
	call mult(r,v,x)
!
	l=neib(1,l)
!
        do mm=1,3
        do nn=1,3
        r(mm,nn)=x(mm,nn)
        enddo
        enddo
!
120	continue
	do 130 ij=1,JD
!
        call uextract(v,4,l)
        call mult(r,v,x)
!
	l=neib(4,l)
!
        do mm=1,3
        do nn=1,3
        r(mm,nn)=x(mm,nn)
        enddo
        enddo
!
130	continue
	do 140 ji=1,jd
!
        l=neib(-1,l)
        call udagextract(v,1,l)
        call mult(r,v,x)
!
        do mm=1,3
        do nn=1,3
        r(mm,nn)=x(mm,nn)
        enddo
        enddo
!
140	continue
	do 150 jj=1,jd
!
        l=neib(-4,l)
        call udagextract(v,4,l)
        call mult(r,v,x)
!
        do mm=1,3
        do nn=1,3
        r(mm,nn)=x(mm,nn)
        enddo
        enddo
!
150	continue
        nar=(jt-1)*ncu+(jz-1)*ns+(jy-1)*ni+jx
!
        do mm=1,3
        do nn=1,3
        q1(mm,nn,nar)=r(mm,nn)
        enddo
        enddo
!
	ll=ll+1
!
        r(1,1)=dcmplx(1d0,0d0)
        r(1,2)=dcmplx(0d0,0d0)
        r(1,3)=dcmplx(0d0,0d0)
        r(2,1)=dcmplx(0d0,0d0)
        r(2,2)=dcmplx(1d0,0d0)
        r(2,3)=dcmplx(0d0,0d0)
        r(3,1)=dcmplx(0d0,0d0)
        r(3,2)=dcmplx(0d0,0d0)
        r(3,3)=dcmplx(1d0,0d0)
!
110	continue
160	continue
170	continue
175	continue

!	GENERATES AVERAGE OF WILSON LOOPS
!-------------------------------------------------------------
	sum=0
	do i=1,nsite
	sum=sum+real(q1(1,1,i)+q1(2,2,i)+q1(3,3,i))
	enddo
	ave=sum/(3*nsite)
	write(34,*) ave
300 	continue
	return
	end




	subroutine loopyt
!
! COMPUTES Y-T WILSON LOOP
!                                                                       
      use lattsize
      use precdef
      use umastcom
      use neibcom
!
      implicit none
!
      complex (kind=double), dimension(3,3) :: x,r,v,w
      complex (kind=double), dimension(3,3,nsite) :: q1
      integer :: i,jd,ll,jt,jz,jy,jx,l,ii,ij,ji,jj
      integer :: mm,nn,nar
      real :: sum,ave
!
!----------------------------------
	do 300 jd=mjd,mjd
	ll=1
!
	r(1,1)=dcmplx(1d0,0d0)
        r(1,2)=dcmplx(0d0,0d0)
        r(1,3)=dcmplx(0d0,0d0)
        r(2,1)=dcmplx(0d0,0d0)
        r(2,2)=dcmplx(1d0,0d0)
        r(2,3)=dcmplx(0d0,0d0)
        r(3,1)=dcmplx(0d0,0d0)
        r(3,2)=dcmplx(0d0,0d0)
        r(3,3)=dcmplx(1d0,0d0)	
!	
	do 175 jt=1,leng4
	do 170 jz=1,leng3
	do 160 jy=1,leng2
	do 110 jx=1,leng1
        l=(jt-1)*ncu+(jz-1)*ns+(jy-1)*ni+jx
        do 120 ii=1,jd
!
        call uextract(v,2,l)
	call mult(r,v,x)
!
	l=neib(2,l)
!
        do mm=1,3
        do nn=1,3
        r(mm,nn)=x(mm,nn)
        enddo
        enddo
!
120	continue
	do 130 ij=1,JD
!
        call uextract(v,4,l)
        call mult(r,v,x)
!
	l=neib(4,l)
!
        do mm=1,3
        do nn=1,3
        r(mm,nn)=x(mm,nn)
        enddo
        enddo
!
130	continue
	do 140 ji=1,jd
!
        l=neib(-2,l)
        call udagextract(v,2,l)
        call mult(r,v,x)
!
        do mm=1,3
        do nn=1,3
        r(mm,nn)=x(mm,nn)
        enddo
        enddo
!
140	continue
	do 150 jj=1,jd
!
        l=neib(-4,l)
        call udagextract(v,4,l)
        call mult(r,v,x)
!
        do mm=1,3
        do nn=1,3
        r(mm,nn)=x(mm,nn)
        enddo
        enddo
!
150	continue
        nar=(jt-1)*ncu+(jz-1)*ns+(jy-1)*ni+jx
!
        do mm=1,3
        do nn=1,3
        q1(mm,nn,nar)=r(mm,nn)
        enddo
        enddo
!
	ll=ll+1
!
        r(1,1)=dcmplx(1d0,0d0)
        r(1,2)=dcmplx(0d0,0d0)
        r(1,3)=dcmplx(0d0,0d0)
        r(2,1)=dcmplx(0d0,0d0)
        r(2,2)=dcmplx(1d0,0d0)
        r(2,3)=dcmplx(0d0,0d0)
        r(3,1)=dcmplx(0d0,0d0)
        r(3,2)=dcmplx(0d0,0d0)
        r(3,3)=dcmplx(1d0,0d0)
!
110	continue
160	continue
170	continue
175	continue

!	GENERATES AVERAGE OF WILSON LOOPS
!-------------------------------------------------------------
	sum=0
	do i=1,nsite
	sum=sum+real(q1(1,1,i)+q1(2,2,i)+q1(3,3,i))
	enddo
	ave=sum/(3*nsite)
	write(35,*) ave
300 	continue
	return
	end





	subroutine loopzt
!
! COMPUTES Z-T WILSON LOOP
!                                                                       
      use lattsize
      use precdef
      use umastcom
      use neibcom
!
      implicit none
!
      complex (kind=double), dimension(3,3) :: x,r,v,w
      complex (kind=double), dimension(3,3,nsite) :: q1
      integer :: i,jd,ll,jt,jz,jy,jx,l,ii,ij,ji,jj
      integer :: mm,nn,nar
      real :: sum,ave
!
!----------------------------------
	do 300 jd=mjd,mjd
	ll=1
!
	r(1,1)=dcmplx(1d0,0d0)
        r(1,2)=dcmplx(0d0,0d0)
        r(1,3)=dcmplx(0d0,0d0)
        r(2,1)=dcmplx(0d0,0d0)
        r(2,2)=dcmplx(1d0,0d0)
        r(2,3)=dcmplx(0d0,0d0)
        r(3,1)=dcmplx(0d0,0d0)
        r(3,2)=dcmplx(0d0,0d0)
        r(3,3)=dcmplx(1d0,0d0)	
!	
	do 175 jt=1,leng4
	do 170 jz=1,leng3
	do 160 jy=1,leng2
	do 110 jx=1,leng1
        l=(jt-1)*ncu+(jz-1)*ns+(jy-1)*ni+jx
        do 120 ii=1,jd
!
        call uextract(v,3,l)
	call mult(r,v,x)
!
	l=neib(3,l)
!
        do mm=1,3
        do nn=1,3
        r(mm,nn)=x(mm,nn)
        enddo
        enddo
!
120	continue
	do 130 ij=1,JD
!
        call uextract(v,4,l)
        call mult(r,v,x)
!
	l=neib(4,l)
!
        do mm=1,3
        do nn=1,3
        r(mm,nn)=x(mm,nn)
        enddo
        enddo
!
130	continue
	do 140 ji=1,jd
!
        l=neib(-3,l)
        call udagextract(v,3,l)
        call mult(r,v,x)
!
        do mm=1,3
        do nn=1,3
        r(mm,nn)=x(mm,nn)
        enddo
        enddo
!
140	continue
	do 150 jj=1,jd
!
        l=neib(-4,l)
        call udagextract(v,4,l)
        call mult(r,v,x)
!
        do mm=1,3
        do nn=1,3
        r(mm,nn)=x(mm,nn)
        enddo
        enddo
!
150	continue
        nar=(jt-1)*ncu+(jz-1)*ns+(jy-1)*ni+jx
!
        do mm=1,3
        do nn=1,3
        q1(mm,nn,nar)=r(mm,nn)
        enddo
        enddo
!
	ll=ll+1
!
        r(1,1)=dcmplx(1d0,0d0)
        r(1,2)=dcmplx(0d0,0d0)
        r(1,3)=dcmplx(0d0,0d0)
        r(2,1)=dcmplx(0d0,0d0)
        r(2,2)=dcmplx(1d0,0d0)
        r(2,3)=dcmplx(0d0,0d0)
        r(3,1)=dcmplx(0d0,0d0)
        r(3,2)=dcmplx(0d0,0d0)
        r(3,3)=dcmplx(1d0,0d0)
!
110	continue
160	continue
170	continue
175	continue

!	GENERATES AVERAGE OF WILSON LOOPS
!-------------------------------------------------------------
	sum=0
	do i=1,nsite
	sum=sum+real(q1(1,1,i)+q1(2,2,i)+q1(3,3,i))
	enddo
	ave=sum/(3*nsite)
	write(36,*) ave
300 	continue
	return
	end
