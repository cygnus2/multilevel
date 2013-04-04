         program main
!
      use lattsize
      use precdef
      use filecom
      use stapcom
      use umastcom
      use polcom
      use prodcom
      use prodncom
      use ranlxd_generator
!
      implicit none
!
      real :: beta
      real(kind=double) :: xnorm
      complex(kind=double) :: pave,pav
      integer :: i,j,k,istart,mm,nn,jj,iupd,iti
      integer :: iseed
      integer :: l1,l2,l3,l4,l5,l6,l7,l8,ll1,ll2,ll3,ll4,ll5,ll6
      integer :: m1,m2,m3,m4
      integer :: itherm,imes
      complex(kind=double),dimension(3,3,ncu,3,leng4/2) :: TNF,TF
      complex(kind=double),dimension(3,3,leng4/itdim) :: TN1,T1
      complex(kind=double),dimension(3,3,ncu) :: PF
      complex(kind=double),dimension(3,3,1+leng4/itdim) :: PF1
!
       read(*,*) istart,beta,iseed,iupd,itherm,imes
       write(*,*) istart,beta,iseed,iupd,itherm,imes
!
!  Setting up address variables
!
        call mk_index
!
        call rlxd_init(1,iseed)
        call init(istart)
!
        write(*,*)"Starting with itdim=",itdim," beta=",beta
        write(*,*)"Going to do ",imes," measurements"
!   Thermalisation
        if(istart.ge.0) then
         do j=1,itherm
!         do i=1,10
         call cabmar(beta) 
         do k=1,3
         call overrelax
         enddo
!         enddo
!         call unitar
        call plaqmeas
         enddo
        endif
!   Measurements
!!        write(*,*) 'Started with measurements'
!!        do j=1,imes
!        if(mod(j,50).eq.0) then
!          call storegauge(j)
!        endif
!!        do i=1,10
!!        call cabmar(beta)
!!        do k=1,3
!!        call overrelax
!!        enddo 
!!        enddo
!
! Going into measurement routine
!
!!       TP=0
!!       TNP=0
!!       POL=0
!
!!       do jj=1,iupd
!
!!          call staple         ! computes the staples
!!          call multihit(beta)      ! analytic SU(3) integral (de Forcrand)
!
!!          call makePOL
!!          call makeTP
!!          call makeTNP
!
!!          do iti=1,leng4,itdim
!!             call b_update(beta,iti)
!!             call bover(iti)
!!             call sub_update(beta,iti)
!!             call subover(iti)
!!          enddo
!!       enddo
!
!!        POL=POL/float(iupd)
!!        TP=TP/float(iupd)
!!        TNP=TNP/float(iupd)
!
!   Calculate the Polyakov loop
!
!!        do l1=1,ncu
!
!!             PF1=0
!!             PF1(1,1,1)=cmplx(1d0,0d0)
!!             PF1(2,2,1)=cmplx(1d0,0d0)
!!             PF1(3,3,1)=cmplx(1d0,0d0)
!
!!        do l8=1,leng4/itdim
!!           do l3=1,3
!!           do l4=1,3
!!              PF(l3,l4,l1)=cmplx(0d0,0d0)
!!              do ll1=1,3
!!                 PF(l3,l4,l1)=PF(l3,l4,l1) + PF1(l3,ll1,l8)*POL(ll1,l4,l1,l8)
!!              end do
!!              PF1(l3,l4,l8+1)=PF(l3,l4,l1)
!!           end do
!!           end do
!!        end do
!
!!        end do
!
!   Polyakov loop measurement
!
!!           pave=0
!!           do m1=1,ncu
!!           do m3=1,3
!!              pave=pave+PF(m3,m3,m1)
!!           end do
!!           end do
!!           pav=pave/float(3*ncu)
!!        write(40,*) REAL(pav),AIMAG(pav)
!
!   Calculate the correlation functions
!
!!        do l2=1,dir-1
!!        do l1=1,ncu
!
!!             T1=0
!!             TN1=0
!
!!        do l5=1,leng4/2
!!        T1(:,:,1)=TP(:,:,l1,l2,1,l5)
!!        TN1(:,:,1)=TNP(:,:,l1,l2,1,l5)
!!        do l8=2,leng4/itdim
!!           do l3=1,3
!!           do l4=1,3
!!              TF(l3,l4,l1,l2,l5)=cmplx(0d0,0d0)
!!              TNF(l3,l4,l1,l2,l5)=cmplx(0d0,0d0)
!!              do ll1=1,3
!!                 TF(l3,l4,l1,l2,l5)=TF(l3,l4,l1,l2,l5) + T1(l3,ll1,l8-1)*TP(ll1,l4,l1,l2,l8,l5)
!!                 TNF(l3,l4,l1,l2,l5)=TNF(l3,l4,l1,l2,l5) + TN1(l3,ll1,l8-1)*TNP(ll1,l4,l1,l2,l8,l5)
!!              end do
!!              T1(l3,l4,l8)=TF(l3,l4,l1,l2,l5)
!!              TN1(l3,l4,l8)=TNF(l3,l4,l1,l2,l5)
!!            end do
!!            end do
!!        end do
!
!!        end do
!
!!        end do
!!        end do
!
!       Electric Field Correlator Measurements
!
!!           do m4=1,leng4/2
!!           pave=0
!!           do m2=1,dir-1
!!           do m1=1,ncu
!!           do m3=1,3
!!              pave=pave+TNF(m3,m3,m1,m2,m4)+TF(m3,m3,m1,m2,m4)
!!           end do
!!           end do
!!           end do
!!           pav=pave/float(2*3*3*ncu)
!!           write(40+m4,*) REAL(pav),AIMAG(pav)
!!           end do
!        call plaqmeas
!
! Measurement complete
!
!!        call unitar
!!        end do
!10       continue
!        call storegauge(j)
         end
