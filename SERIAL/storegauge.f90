      subroutine storegauge(it)
!    ********************************

      use lattsize
      use precdef
      use filecom
      use umastcom

      implicit none
!    -------------------------------
      integer i,j,it
      open(11, file='fort.11', status='replace')

      do i = 1,nsite
      do j = 1,4

      write(11,*) umast(1,1,j,i)
      write(11,*) umast(1,2,j,i)
      write(11,*) umast(1,3,j,i)
      write(11,*) umast(2,1,j,i)
      write(11,*) umast(2,2,j,i)
      write(11,*) umast(2,3,j,i)
      write(11,*) umast(3,1,j,i)
      write(11,*) umast(3,2,j,i)
      write(11,*) umast(3,3,j,i)

      end do
      end do

      close(11)
      write(*,*) 'Gauge config written at iteration ',it
    
      return
      end
