      subroutine init_arrays

      use precdef

      implicit none

      real (kind=double), dimension(16) :: rabscis, weight
      real (kind=double) :: radius, pi, twopi
      integer :: k, iflag
      
      complex(kind=double), dimension(16,3) :: x
      complex(kind=double), dimension(16) :: weight_x

      common / init_cons / x, weight_x

      iflag = 1
      pi = dacos(-1d0)
      twopi = 2. * pi

      radius = 0.01d0              !! Radius der Kontour des Integrals
  
      rabscis( 1) = 0.99726386d0
      rabscis( 2) = 0.98561151d0
      rabscis( 3) = 0.96476226d0 
      rabscis( 4) = 0.93490608d0
      rabscis( 5) = 0.89632116d0
      rabscis( 6) = 0.84936761d0
      rabscis( 7) = 0.79448380d0
      rabscis( 8) = 0.73218212d0
      rabscis( 9) = 0.66304427d0
      rabscis(10) = 0.58771576d0
      rabscis(11) = 0.50689991d0
      rabscis(12) = 0.42135128d0
      rabscis(13) = 0.33186860d0
      rabscis(14) = 0.23928736d0
      rabscis(15) = 0.14447196d0
      rabscis(16) = 0.04830767d0
  
      weight( 1) = 0.00701861d0
      weight( 2) = 0.01627439d0
      weight( 3) = 0.02539207d0
      weight( 4) = 0.03427386d0
      weight( 5) = 0.04283590d0
      weight( 6) = 0.05099806d0
      weight( 7) = 0.05868409d0
      weight( 8) = 0.06582222d0
      weight( 9) = 0.07234579d0
      weight(10) = 0.07819390d0
      weight(11) = 0.08331192d0
      weight(12) = 0.08765209d0
      weight(13) = 0.09117388d0
      weight(14) = 0.09384440d0
      weight(15) = 0.09563872d0
      weight(16) = 0.09654009d0
  
      do k = 1, 16
        x(k,1) = radius * cdexp(dcmplx(0.,pi * rabscis(k)))
        x(k,2) = x(k,1)**2
        x(k,3) = x(k,2) * x(k,1)
        weight_x(k) = x(k,1) * weight(k)
      end do
    
      return
      end

      subroutine adjlink(A,B)

      use precdef

      complex(kind=double), dimension(3,3) :: A, B
      integer :: i,j
      
      do i=1,3
      do j=1,3
      B(i,j)=conjg(A(j,i))
      enddo
      enddo
 
      return
      end

      subroutine Det(u,c)

      use precdef

      implicit none

      complex(kind=double), dimension(3,3) :: u
      complex(kind=double) :: c

       c = u(1,1) * u(2,2) * u(3,3) &
         + u(2,1) * u(3,2) * u(1,3) &
         + u(3,1) * u(1,2) * u(2,3) &
         - u(1,3) * u(2,2) * u(3,1) &
         - u(1,2) * u(2,1) * u(3,3) &
         - u(1,1) * u(2,3) * u(3,2)

      return
      end



!! Original, scalar version
      subroutine gauss_integ(beta,staple,int_link)
!! ==============================================================
!! Integriert die Staple der links in Raum- & Zeit-Richtung, 
!! um ein besseres Signal fuer die Wilsonloops zu erhalten
!! Semi-analytische Methode nach  Ph. de Forcrand
!! ==============================================================
      
      use precdef
      
      implicit none

      complex(kind=double),dimension(0:2,0:2) :: int_link,staple,st,st_a  
      complex(kind=double),dimension(0:2,0:2) :: stst_a,st1, st3, st4, tmp1
      complex(kind=double) :: udet1, udet2, udet3,p,px,xpsqr,tpz,z_,zi
      complex(kind=double) :: c_tmp,bessel1, bessel2
      complex(kind=double), dimension(16,3) :: x
      complex(kind=double), dimension(16) :: weight_x
      real(kind=double) :: a,b,c,q,rpx,int0r,int1r,int2r,int3r,int4r
      real(kind=double) :: pi, twopi
      integer k,j1
      real :: beta
  
      common / init_cons / x, weight_x

!!      call init_arrays
      pi = dacos(-1d0)
      twopi = 2d0 * pi

!!
!! zuerst werden die Komponenten des Integranden berechnet
!! 
      st = beta / 6.0 * staple   !! staple und staple^adjungiert
      call adjlink(st,st_a)      !! muessen vertauscht sein !
!!                q = Det(st) + Det(st~) = 2*Re( Det(st) )
      call Det(st,q)
      q=2d0 * q
                                 !! stst_a ist symetrisch und reell
!!  stst_a = st * st_a           !! stst_a.[0,0]=[1,1]=[2,2]
                                 !! so sind auch a,b,c reell
      call mult(st,st_a,stst_a)
!!                a = Tr st*st~
      a = stst_a(0,0)+stst_a(1,1)+stst_a(2,2)

!!                b = 1/2 ( (Tr st*st~)^2 - Tr (st*st~)^2 )
      call mult(stst_a,stst_a,tmp1)
      b = 0.5 * (a**2 - (tmp1(0,0)+tmp1(1,1)+tmp1(2,2)))

!!                c = Det( st*st~ )
      call Det(stst_a,c) 

!!                st1 = d q / d st~

      st1(0,0) = st_a(1,1) * st_a(2,2) - st_a(1,2) * st_a(2,1)
      st1(0,1) = st_a(0,2) * st_a(2,1) - st_a(0,1) * st_a(2,2)
      st1(0,2) = st_a(0,1) * st_a(1,2) - st_a(0,2) * st_a(1,1)
      st1(1,0) = st_a(1,2) * st_a(2,0) - st_a(1,0) * st_a(2,2)
      st1(1,1) = st_a(0,0) * st_a(2,2) - st_a(0,2) * st_a(2,0)
      st1(1,2) = st_a(0,2) * st_a(1,0) - st_a(0,0) * st_a(1,2)
      st1(2,0) = st_a(1,0) * st_a(2,1) - st_a(1,1) * st_a(2,0)
      st1(2,1) = st_a(0,1) * st_a(2,0) - st_a(0,0) * st_a(2,1)
      st1(2,2) = st_a(0,0) * st_a(1,1) - st_a(0,1) * st_a(1,0)
!!
!!             st2 = d Tr st*st~ / d st~ = st  (keine Rechnung notwendig)
!!
!!             1/2  d (Tr st*st~)^2 / d st~ = (Tr st*st~) * st
!!             1/2  d Tr (st*st~)^2 / d st~ = st*st~ * st
!!             st3 =  (Tr st*st~) * st - st*st~ * st
      call mult(stst_a,st,tmp1)
      st3(:,:) = a * st(:,:) - tmp1(:,:)
!!
!!             st4 = d Det(st*st~) / d st~
!!

!! Entwicklung nach 1.Spalte
      udet1 = stst_a(1,1) * stst_a(2,2) - stst_a(1,2) * stst_a(2,1) 
      udet2 = stst_a(2,1) * stst_a(0,2) - stst_a(2,2) * stst_a(0,1)
      udet3 = stst_a(0,1) * stst_a(1,2) - stst_a(0,2) * stst_a(1,1)

      do j1 = 0, 2
        st4(0,j1) = st(0,j1)*udet1 + st(1,j1)*udet2 + st(2,j1)*udet3
      end do
!!
!! Entwicklung nach 2.Spalte
      udet1 = stst_a(1,2) * stst_a(2,0) - stst_a(1,0) * stst_a(2,2) 
      udet2 = stst_a(2,2) * stst_a(0,0) - stst_a(2,0) * stst_a(0,2)
      udet3 = stst_a(0,2) * stst_a(1,0) - stst_a(0,0) * stst_a(1,2)

      do j1 = 0, 2
        st4(1,j1) = st(0,j1)*udet1 + st(1,j1)*udet2 + st(2,j1)*udet3
      end do
!!
!! Entwicklung nach 3.Spalte
      udet1 = stst_a(1,0) * stst_a(2,1) - stst_a(1,1) * stst_a(2,0) 
      udet2 = stst_a(2,0) * stst_a(0,1) - stst_a(2,1) * stst_a(0,0)
      udet3 = stst_a(0,0) * stst_a(1,1) - stst_a(0,1) * stst_a(1,0)

      do j1 = 0, 2
        st4(2,j1) = st(0,j1)*udet1 + st(1,j1)*udet2 + st(2,j1)*udet3
      end do
  
!! ------------------------------------------------------------------
!! Hier beginnt die Gauss'sche-Integration mit 32 Stuetzstellen auf
!! einem Kreis um den Ursprung. Das Integral der oberen Kontour ist
!! das komplex konjugierte der unteren Kontour => nur das halbe
!! Integral muss berechnet werden.
!! ------------------------------------------------------------------

      int0r = 0.
      int1r = 0.
      int2r = 0.
      int3r = 0.
      int4r = 0.
      do k = 1, 16
        p  = dcmplx(1.0,0.0) + a * x(k,1) + b * x(k,2) + c * x(k,3)
        px = p / x(k,1)

        rpx = dsqrt(real(px)**2 + aimag(px)**2)
        z_  = dcmplx(dsqrt(2.*(rpx+real(px)))  &
                  ,dsqrt(2.*(rpx-real(px)))) !!z=2*(P(x)/x)^1/2
        if (aimag(px) .lt. 0.) z_ = conjg(z_)       !! damit |arg(z)| < Pi/2
    
        rpx   = 4.0 * rpx                           !! nun rpx = |z|
        zi    = dconjg(z_) / rpx 
        xpsqr = 2.0 * zi                            !! xpsqr = 2.0 / z_
    
        tpz   = cdsqrt(twopi * z_ )             !! = (2Pi*z)^1/2

        bessel1 = dcmplx(1.0,0.0) + zi * (-0.375           &
                                + zi * (-0.1171875         &
                                + zi * (-0.1025390625      &
                                + zi * (-0.1441955566      &
                                - zi *   0.2775764465       ))))

        bessel2 = dcmplx(1.0,0.0) + zi * (-1.875           &
                                + zi * ( 0.8203125         &
                                + zi * ( 0.3076171875      &
                                + zi * ( 0.3172302246      & 
                                + zi *   0.5154991150       ))))    

        c_tmp = cdexp(x(k,1) * q + z_) / tpz * weight_x(k)

        bessel1 = bessel1 * c_tmp * xpsqr
        bessel2 = bessel2 * c_tmp / p    

        int0r = int0r +           bessel1  !! Z
        int1r = int1r + x(k,1) *  bessel1  !! 1.Teil von dZ/dst~, Fak. dQ/dst~
        int2r = int2r + x(k,1) *  bessel2  !! 2.Teil von dZ/dst~
        int3r = int3r + x(k,2) *  bessel2  !! Ableitungen der 3 Teile von P
        int4r = int4r + x(k,3) *  bessel2  !! sie sind unabh. von z

      end do
     
      int1r = int1r / int0r   !! log. Ableitung von Z nach staple
      int2r = int2r / int0r   !! U = 1/Z *  dZ/dst~
      int3r = int3r / int0r   
      int4r = int4r / int0r   

      int_link = int1r*st1 + int2r*st + int3r*st3 + int4r*st4     !! st2=st
     
      return
      end
