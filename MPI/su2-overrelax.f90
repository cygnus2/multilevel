       subroutine su2relax(smsta,un)
!
        use precdef
        implicit none
!
        integer :: ii
        real (kind=double), dimension (4) :: smsta,un
        real (kind=double) :: factor, det
!
!	Overrelaxation starts
!
        det=0d0
        do ii=1,4
        det=det+smsta(ii)*smsta(ii)
        enddo
!
        factor=2d0*smsta(1)/(DABS(det))
!
        un(1)=factor*smsta(1)-1d0
        un(2)=factor*smsta(2)
        un(3)=factor*smsta(3)
        un(4)=factor*smsta(4)
!	
!	end of overrelaxation 
!
       return
       end
