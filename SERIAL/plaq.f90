     subroutine plaqmeas

!     As a cross check I evaluate the expectation value of the
!     plaquette for the background gauge configuration.

      use lattsize
      use precdef
      use neibcom

      implicit none

      integer :: x,mu,nu,xmu,xnu
      complex (kind=double) :: up(3,3),u2(3,3),u3(3,3),u4(3,3)
      real (kind=double) :: plaqs,plaqt


      plaqs = 0.0d0
      plaqt = 0.0d0

      do x = 1,nsite
      do mu = 1,3
      do nu = mu+1,4

      xmu = neib(mu,x)
      xnu = neib(nu,x)

      call uextract(up,mu,x)
      call uextract(u2,nu,xmu)
      call udagextract(u3,mu,xnu)
      call udagextract(u4,nu,x)

      call umult(up,u2)
      call umult(up,u3)
      call umult(up,u4)

      if(nu == 4) then
        plaqt = plaqt + real(up(1,1) + up(2,2) + up(3,3))
      else
        plaqs = plaqs + real(up(1,1) + up(2,2) + up(3,3))
      endif

      enddo
      enddo
      enddo

      plaqs = plaqs/(9.d0*float(nsite))
      plaqt = plaqt/(9.d0*float(nsite))

      write(20,*) plaqs,plaqt

      return
      end
