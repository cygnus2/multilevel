     subroutine plaqmeas(rnk)

!     As a cross check I evaluate the expectation value of the
!     plaquette for the background gauge configuration.

      use lattsize
      use precdef
      use neibcom
      use paramod

      implicit none

      integer  :: x,mu,nu,xmu,xnu,rnk
      complex (kind=double) :: up(3,3),u2(3,3),u3(3,3),u4(3,3)
      real (kind=double) :: plaqs,plaqt

      plaqs = 0.0d0
      plaqt = 0.0d0

      do x = ncu,itdim*ncu
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

      plaqs = plaqs/(9.d0*float(ncu*(itdim-1)))
      plaqt = plaqt/(9.d0*float(ncu*(itdim-1)))
      write(20+rnk,*) plaqs,plaqt

      return
      end
