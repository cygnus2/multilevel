       subroutine su2(smsta,sbeta,un)
!
        use precdef
        use ranlxd_generator
        implicit none
!
        integer :: ii
        real (kind=double), dimension (4) :: smsta,un,ua,s
        real (kind=double), dimension (6) :: r
        real (kind=double) :: pi,det,ss,yl,y,xi,yy,arho,rho
        real :: sbeta
!
!	Heat bath starts
!
        pi=2d0*dasin(1d0)
        det=0d0
        do ii=1,4
        det=det+smsta(ii)*smsta(ii)
        enddo
!
        arho=dsqrt(DABS(det))
        rho=sbeta*arho
        s(1)=smsta(1)/arho
        s(2)=smsta(2)/arho
        s(3)=smsta(3)/arho
        s(4)=smsta(4)/arho
!
10      continue
        call ranlxd(r)
        ss=dsin(pi*r(1)/2.)
        yl= -dlog(1.-r(2))
        y=yl*ss*ss - dlog(1.-r(3))
        ua(1)=1.-(y/rho)
        if ((2*r(6)*r(6)).LE.(1+ua(1))) then
        xi=dsqrt(1-ua(1)*ua(1))
        yy=1.-2.*r(4)
        ua(2)=yy*xi
        ua(3)=(dsqrt(1-yy*yy))*dcos(2.*pi*r(5))*xi
        ua(4)=(dsqrt(1-yy*yy))*dsin(2.*pi*r(5))*xi
!
        un(1)=ua(1)*s(1)-ua(2)*s(2)-ua(3)*s(3)-ua(4)*s(4)
        un(2)=ua(1)*s(2)+ua(2)*s(1)-ua(3)*s(4)+ua(4)*s(3)
        un(3)=ua(1)*s(3)+ua(2)*s(4)+ua(3)*s(1)-ua(4)*s(2)
        un(4)=ua(1)*s(4)-ua(2)*s(3)+ua(3)*s(2)+ua(4)*s(1)
         else
          go to 10
        endif
!	
!	end of heat bath
!
       return
       end
