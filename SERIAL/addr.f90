!    ******************************************************************
      subroutine mk_index
!     ******************************************************************
!     **                    **                                        **
!     **  MK_INDEX          **    Program by C.B. Lang                **
!     **                    **                                        **
!     ******************************************************************
!     ------------------------------------------------------------------

      use lattsize
      use precdef
      use neibcom

      implicit none
!     ------------------------------------------------------------------
      
!     ------------------------------------------------------------------
!     Program for construction of d=4 index vectors
!     ------------------------------------------------------------------
!     Part of package:  
!     ------------------------------------------------------------------
!     Output variables:   neib(-mu:mu,1:nsite)
!     ------------------------------------------------------------------
!     Remarks: neib (1:mu,x) = Index(x+mu)
!              neib (-mu:-1,x) = Index(x-mu)
!              siteindex(i1,i2,i3,i4)=(((i4*leng3+i3)*leng2+i2)*leng1+i1    
!    ------------------------------------------------------------------
!     ******************************************************************
!     ------------------------------------------------------------------
      integer idir,n,i1,i2,i3,i4,i1p,i2p,i3p,i4p,i1m,i2m,i3m,i4m,&
     &        is(4,4),j
      data is /1,0,0,0,  0,1,0,0,  0,0,1,0,  0,0,0,1/
!     ------------------------------------------------------------------
!     first the n.n. index arrays
!     ------------------------------------------------------------------
      do idir=1,4
        do i4=0,leng4-1
          i4p=mod(i4+is(4,idir),leng4)
          i4m=mod(i4-is(4,idir)+leng4,leng4)
          do i3=0,leng3-1
            i3p=mod(i3+is(3,idir),leng3)
            i3m=mod(i3-is(3,idir)+leng3,leng3)
            do i2=0,leng2-1
              i2p=mod(i2+is(2,idir),leng2)
              i2m=mod(i2-is(2,idir)+leng2,leng2)
              do i1=0,leng1-1
                i1p=mod(i1+is(1,idir),leng1)
                i1m=mod(i1-is(1,idir)+leng1,leng1)
                n=1+((i4*leng3+i3)*leng2+i2)*leng1+i1
                neib(0,n)=n
                neib(idir,n) =1+((i4p*leng3+i3p)*leng2+i2p)*leng1+i1p
                neib(-idir,n)=1+((i4m*leng3+i3m)*leng2+i2m)*leng1+i1m
              enddo
            enddo
          enddo
        enddo
      enddo
      
      write(*,*) 
      write(*,*)' Neigbour indices computed.'
      write(*,*)
      return
      end
