  module precdef
      
      integer, parameter :: double=selected_real_kind(p=14,r=30)
      integer, parameter :: single=selected_real_kind(6)
      
  end module precdef
!!
  module lattsize

      use precdef

      implicit none
      save

!     lattice dimensions

      integer, parameter :: leng1=12, leng2=12, leng3=12, leng4=4,ni=leng1
      integer, parameter :: ns=ni*leng2,ncu=ns*leng3,nsite=ncu*leng4
      integer, parameter :: dir=4, itdim=2
  end module lattsize
!!
  module neibcom

      use lattsize

      implicit none
      integer, save,  dimension(-4:4,nsite) :: neib


  end module neibcom
!!
  module umastcom

      use precdef
      use lattsize
      implicit none
      save
      complex (kind=double), dimension(3,3,4,nsite) :: umast

  end module umastcom
!!
  module stapcom 

      use precdef
      use lattsize
      implicit none
      save
      complex (kind=double), dimension(3,3,4,nsite) :: staples

  end module stapcom 
!!
  module polcom

       use precdef
       use lattsize
       implicit none
       save
       complex(kind=double),dimension(3,3,ncu,leng4/itdim):: POL

  end module polcom
!!
  module prodcom

       use precdef
       use lattsize
       implicit none
       save
       complex(kind=double),dimension(3,3,ncu,dir-1,leng4/itdim,leng4/2):: TP

  end module prodcom
!!
  module prodncom

       use precdef
       use lattsize
       implicit none
       save
       complex(kind=double),dimension(3,3,ncu,dir-1,leng4/itdim,leng4/2):: TNP

  end module prodncom
!!
  module hitcom

       use precdef
       use lattsize
       implicit none
       save
       complex (kind=double),dimension(3,3,nsite):: hit

 end module hitcom
!!
  module filecom

      implicit none
      save
      character (len=64):: cfile,ufile,dfile
      
  end module filecom
