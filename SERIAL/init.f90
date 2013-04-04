!     ******************************************************************
      subroutine init(gauge_init)
!     ******************************************************************
!     ** indirac            **   Program by C.B. Lang, 27/02/00       **  
!     **                    **   last modified:2/2/04 (P.Majumdar)    **  
!     ******************************************************************
!     ------------------------------------------------------------------

      use lattsize
      use precdef
      use filecom
      
      implicit none
!     ------------------------------------------------------------------
!     The files cfile, ufile, dfile are defined externally
!     (dfile is not used for Wilson)
!     ------------------------------------------------------------------
!     Initialize the gauge configuration:
!
!     Inputflags:      gauge_init <0    config. read from file  
!			   	  =0    cold start          
!				  >0    hot start
!
!     ------------------------------------------------------------------
      integer :: gauge_init
!     ------------------------------------------------------------------

      write(*,*) 
      write(*,100)leng1,leng2,leng3,leng4
      write(*,*)

100   format('  Lattice: ',i4,3(' x',i4))

      if(gauge_init.lt.0) then
 
!       assuming here that ufile contains already the 
!       name of gauge configuration file

        call readgauge
	else
            if(gauge_init.eq.0) then
               call orderedgauge
            else
               call disorderedgauge
            endif

      endif

      end
!     ----------------------------------------------------------

      subroutine readgauge

!     Reads the gauge field configuration from file.


      use lattsize
      use precdef
      use filecom
      use umastcom

      implicit none

      integer i,j

!      open(11, file=ufile, status='old')
       write(*,*) 'Going to open fort.11'
       open(11, file='fort.11', status='old')
       write(*,*) 'Opened fort.11,to read data'

!      old format
!      do i = 1,22
!      do i = 1,28
!      read(11,*)
!      end do

      do i = 1,nsite
      do j = 1,4

      read(11,*) umast(1,1,j,i)
      read(11,*) umast(1,2,j,i)
      read(11,*) umast(1,3,j,i)
      read(11,*) umast(2,1,j,i)
      read(11,*) umast(2,2,j,i)
      read(11,*) umast(2,3,j,i)
      read(11,*) umast(3,1,j,i)
      read(11,*) umast(3,2,j,i)
      read(11,*) umast(3,3,j,i)

      end do
      end do

      close(11)

      write(*,*)
      write(*,*) 'Gauge field read from file '
      write(*,*)

!      call plaqmeas

      return
      end
!     ----------------------------------------------------------

      subroutine orderedgauge

! Initializes to ordered gauge field configuration.


      use lattsize
      use precdef
      use filecom
      use umastcom

      implicit none

      integer i,j

      do i = 1,nsite
      do j = 1,4

      umast(1,1,j,i)=dcmplx(1d0,0d0)
      umast(1,2,j,i)=dcmplx(0d0,0d0)
      umast(1,3,j,i)=dcmplx(0d0,0d0)
      umast(2,1,j,i)=dcmplx(0d0,0d0)
      umast(2,2,j,i)=dcmplx(1d0,0d0)
      umast(2,3,j,i)=dcmplx(0d0,0d0)
      umast(3,1,j,i)=dcmplx(0d0,0d0)
      umast(3,2,j,i)=dcmplx(0d0,0d0)
      umast(3,3,j,i)=dcmplx(1d0,0d0)

      end do
      end do

      return
      end
!     ----------------------------------------------------------

      subroutine disorderedgauge

! Initializes to disordered gauge field configuration.


      use lattsize
      use precdef
      use filecom
      use umastcom
      use ranlxd_generator

      implicit none

      integer i,j
      real (kind=double), dimension(12) :: v

      do i = 1,nsite
      do j = 1,4

      call ranlxd(v)

      umast(1,1,j,i)=dcmplx(v(1),v(2))
      umast(1,2,j,i)=dcmplx(v(3),v(4))
      umast(1,3,j,i)=dcmplx(v(5),v(6))
      umast(2,1,j,i)=dcmplx(v(7),v(8))
      umast(2,2,j,i)=dcmplx(v(9),v(10))
      umast(2,3,j,i)=dcmplx(v(11),v(12))
      umast(3,1,j,i)=dcmplx(0d0,0d0)
      umast(3,2,j,i)=dcmplx(0d0,0d0)
      umast(3,3,j,i)=dcmplx(0d0,0d0)

      end do
      end do

      call unitar
      return
      end
