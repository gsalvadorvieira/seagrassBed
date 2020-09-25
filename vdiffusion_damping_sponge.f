      subroutine vdiffusion(sdif,Tdif,udif,vdif,wdif,n)
c     ------------------------------------------------
c     Wind Stress can be specified here
c     Kz*du/dz = -Tx/rho,   Kz*dv/dz = -Ty/rho
c     use level n (explicit computation)
c     computes d2s/dz2 at the cell centers.
c
      implicit logical (a-z)
      include 'header.f'
      integer i,j,k,it,n

c     The viscosity is Kz*Kzmax
c     Use the same value for mixing tracers and momentum
c     CdT and CdN are the normal and tangential drag coeffs Cd

      double precision dia(NKg),dsdzfc(0:NK),
     &     dudzfc(0:NK),dvdzfc(0:NK),dwdzfc(0:NK),dTdzfc(ntr,0:NK),
     &     sdif(NI,NJ,NK),Tdif(ntr,NI,NJ,NK),udif(NI,NJ,NK),
     &     vdif(NI,NJ,NK),wdif(NI,NJ,NK),Kdudzb,Kdvdzb,Kdudzt,
C      &     Kdvdzt,wzkth,fact,CgNL,CsDIA,func,unorml,utangt,
     &     Kdvdzt,wzkth,fact,CsDIA,func,unorml,utangt,
     &     udrag,vdrag,wdrag,facz
      double precision RR,facb,fac,dfactor,Kzmax,facy,rhoinv,taux,tauy
      double precision crowdfac,dxc,dxcinv,CgNLfac,wfac,rinv,
     &     denom,wdampfac(NI),outdamp(NI)

      double precision CgNL(1:NI)

      double precision grassFrac
      
C     AUTOMATICALLY IMPORTED PARAMETER
      parameter(grassFrac=210.d0)

c     negative taux implies upshore wind, neg tauy implies onshore
c=      parameter (tauy= 0.025, taux=0.0)
      parameter (taux= 0.0, tauy=0.0)

C     Kzmax= 1.d-3 (Chapman and Lentz value)  Gawarkiewicz-Chap use 2.d-3
c=== const kz
c     RR is the bottom friction parameter in m/s - why is it necessary?
      parameter (RR= 0.0005, Kzmax= 1.d-5)
c      parameter (RR= 0.d0, Kzmax= 1.d-5)
c     no bottom fric => RR = 0
c     Kv*du/dz = r*u  at the bottom
c
c     rho must be evaluated to compute the bouyancy freq N, if
c     we are using the KPP scheme (subroutine viscosity).

c     Not needed for constant viscosity case
c     call evalrho(rho,n)

c     characteristic CS area (1cm2) of stems
c     parameter (CsA=1.d-4 ) set in header.f
      CsDIA = dsqrt(4.d0*CsA/PI)

      facb = RR*DL/UL
      fac= 1.d0/(UL*DL*delta)
c     fac*Kz = 1/ReV
      fact = DL/UL

c     define zNL, remove/sparsify grass after NI*grassFrac
      do i=1,NI

C       LINEAR
C         if (i.gt.(NI*dampGrass)) then
C C           zNL(i) = 100.d0*(1.d0 - dble(i)/(NI))
C           zNL(i) = 100.d0 - 100.d0/((1.d0-grassFrac)*dble(NI))*
C      &          (dble(i)-dble(NI)*grassFrac)     
C C           zNL(i) = 0.d0
C         else
C           zNL(i) = 100.d0
C         end if

C       TANH STARTING AT FRAC AND STEEPEST BETWEEN THERE AND OUTLET
        zNL(i) = 20.d0 + 80.d0*(1.d0 - 0.5d0*
     &         ( 1.d0 + TANH(15.d0/288.d0*(dble(i)-grassFrac)) ) ) 

      end do

c     Cd= 0.3 (in header), zNL= 100.0, CgNL= 30.0
      do i=1,NI
        CgNL(i) = Cd*zNL(i)
      end do

      do j=1,NJ
         do i=1,NI
c     u,v= 0 at lower bndry
            u(i,j,0,n)= -u(i,j,1,n)
            v(i,j,0,n)= -v(i,j,1,n)

            do k=1,NK-1
               wzkth= wzk(i,j,k)
               dsdzfc(k)= wzkth*(s(i,j,k+1,n)-s(i,j,k,n))
               dudzfc(k)= wzkth*(u(i,j,k+1,n)-u(i,j,k,n))
               dvdzfc(k)= wzkth*(v(i,j,k+1,n)-v(i,j,k,n))
               dwdzfc(k)= wzkth*(w(i,j,k+1,n)-w(i,j,k,n))
               do it=1,ntr
                  dTdzfc(it,k)= wzkth*(T(it,i,j,k+1,n)-T(it,i,j,k,n))
               end do
            end do
            k=0
            dsdzfc(0)= 0.d0
            do it=1,ntr
               dTdzfc(it,0)= 0.d0
            end do
            dudzfc(0)= wzk(i,j,k)*(u(i,j,k+1,n)-u(i,j,k,n))
            dvdzfc(0)= wzk(i,j,k)*(v(i,j,k+1,n)-v(i,j,k,n))
            dwdzfc(0)= wzk(i,j,k)*(w(i,j,k+1,n)-w(i,j,k,n))
c
c changed on Aug 7,2002, for wind stress. If no wind, dudzfc(NK)=0
            dsdzfc(NK)= 0.0
            dudzfc(NK)= 0.0
            dvdzfc(NK)= 0.0
            dwdzfc(NK)= 0.0
            do it=1,ntr
               dTdzfc(it,NK)= 0.0
            end do

c= const visc. kz=1
            do k=0,NK
             Kz(i,j,k)=1.d0

C     ! Modified to add viscosity closer to the end of the tunnel 
C              increase Kz linearly from treshold
C               if (i .gt. NI*3/4) then
C                 Kz(i,j,k) = 1.d0 + (((i-NI*3/4)/(NI/4))**1)*199.d0
C               else
C                 Kz(i,j,k) = 1.d0
C               end if

C       TANH STARTING AT FRAC AND STEEPEST BETWEEN THERE AND OUTLET
             Kz(i,j,k) = 1.d0 + 9.d0 * 0.5d0 *
     &         ( 1.d0 + TANH(15.d0/288.d0*(dble(i)-grassFrac) ) ) 

            end do

            Kdudzb= facb*u(i,j,1,n)
            Kdvdzb= facb*v(i,j,1,n)

            rhoinv = 1.d0/(rho(i,j,NK) + R0)
            Kdudzt= taux*rhoinv*fact
            Kdvdzt= tauy*rhoinv*fact

c     fac= 1.d0/(UL*DL*delta)
c     fac*Kz = 1/ReV
            do k=1,NK-1
               dfactor=  fac*Kzmax*Jac(i,j,k)*wz(i,j,k)
               sdif(i,j,k)= dfactor
     &              *(dsdzfc(k)*Kz(i,j,k)- dsdzfc(k-1)*Kz(i,j,k-1) )
               do it=1,ntr
                  Tdif(it,i,j,k)= dfactor
     &                 *(dTdzfc(it,k)*Kz(i,j,k)-
     &                 dTdzfc(it,k-1)*Kz(i,j,k-1) )
               end do
               udif(i,j,k)= dfactor
     &              *(dudzfc(k)*Kz(i,j,k)- dudzfc(k-1)*Kz(i,j,k-1) )
               vdif(i,j,k)= dfactor
     &              *(dvdzfc(k)*Kz(i,j,k)- dvdzfc(k-1)*Kz(i,j,k-1) )
               wdif(i,j,k)= dfactor
     &              *(dwdzfc(k)*Kz(i,j,k)- dwdzfc(k-1)*Kz(i,j,k-1) )

            end do

            k=NK
c           ----
            dfactor=  fac*Jac(i,j,k)*wz(i,j,k)
            sdif(i,j,k)= dfactor*Kzmax
     &           *(dsdzfc(k)*Kz(i,j,k)- dsdzfc(k-1)*Kz(i,j,k-1) )
            do it=1,ntr
               Tdif(it,i,j,k)= dfactor*Kzmax
     &             *(dTdzfc(it,k)*Kz(i,j,k)-dTdzfc(it,k-1)*Kz(i,j,k-1))
            end do
            udif(i,j,k)= dfactor
     &           *(Kdudzt - dudzfc(k-1)*Kz(i,j,k-1)*Kzmax )
            vdif(i,j,k)= dfactor
     &           *(Kdvdzt - dvdzfc(k-1)*Kz(i,j,k-1)*Kzmax )
            wdif(i,j,k)= dfactor
     &           *(dwdzfc(k)*Kz(i,j,k)- dwdzfc(k-1)*Kz(i,j,k-1))*Kzmax 

         end do
      end do

      goto 502

c     MODEL 1 : Drag proportional to crowding
c     =========
      do k=1,NKg
         do i=1,NI
            dxc= xc(i) -xc(i-1)
            dxcinv = 1.d0/dxc
            do j=1,NJ
               crowdfac = (dxc + xg(i,j,k,n) -xg(i-1,j,k,n))*dxcinv
c     crowdfac should vary between 0 (when touching) and 2 (when distant)
               if (crowdfac.gt.2.d0) crowdfac= 2.d0
               if (crowdfac.lt.0.d0) crowdfac= 0.d0
c     Hence CgNLfac varies between 0 (separation >= 2*orig sep.) and
c     2*CgNL (separation = 0 or overlap)
               CgNLfac = CgNL(i) + (1.d0 -crowdfac)*CgNL(i)
               fac= CgNLfac*Jac(i,j,k)
               udif(i,j,k)= udif(i,j,k)-fac*u(i,j,k,n)*
     &              dabs(u(i,j,k,n))
               vdif(i,j,k)= vdif(i,j,k)-fac*v(i,j,k,n)*
     &              dabs(v(i,j,k,n))
               wdif(i,j,k)= wdif(i,j,k)-fac*w(i,j,k,n)*
     &              dabs(w(i,j,k,n))
            end do
         end do
      end do


c     MODEL 2 : Drag is inversely proportional to deflection, i.e.
c     ========
c     horiz. drag is maximum when grass is vertical, less when grass
c     deflected along flow.

      do k=1,NKg
         do i=1,NI
            dxc= xc(i) -xc(i-1)
            dxcinv = 1.d0/dxc
            do j=1,NJ
               crowdfac = xg(i,j,k,n)*dxcinv
c     crowdfac varies between 0 and 2
               crowdfac= dabs(crowdfac)
               if (crowdfac.gt.2.d0) crowdfac= 2.d0
c     Hence CgNLfac varies between 0 and CgNL
c     CgNL fac is less when plants are flat
c
               CgNLfac = 0.5*(2.d0 -crowdfac)*CgNL(i)
               fac= CgNLfac*Jac(i,j,k)
               udif(i,j,k)= udif(i,j,k) -fac*u(i,j,k,n)
     &              *dabs(u(i,j,k,n))
               vdif(i,j,k)= vdif(i,j,k) -fac*v(i,j,k,n)
     &              *dabs(v(i,j,k,n))
               wdif(i,j,k)= wdif(i,j,k) -fac*w(i,j,k,n)
     &              *dabs(w(i,j,k,n))
            end do
         end do
      end do

c     MODEL 3
c     ======
 501  do k=1,NKg
         func = 1.d0
c         func = (dble(NKg-k)+0.5)/dble(NKg)
c         func = (dble(NKg-k)+0.5)/dble(NKg)*0.8d0 +0.2d0
         dia(k) = CsDIA*func
      end do
      do k=1,NKg
         do j=1,NJ
            do i=1,NI
c     CgNLfac is less when plants are flat, varies bet 0 and CgNL, as
c     the square of the cosine of the angle with the vertical.
               CgNLfac = costhe(i,j,k)*costhe(i,j,k)*CgNL(i)
c               fac= CgNLfac*Jac(i,j,k)
               fac= CgNLfac*Jac(i,j,k)*dia(k)*LEN
c     wfac is prop to the square of sin theta.
c               wfac = (CgNL - CgNLfac)*Jac(i,j,k)
               wfac = (CgNL(i) - CgNLfac)*Jac(i,j,k)*dia(k)*LEN*delta
               udif(i,j,k)= udif(i,j,k) -fac*u(i,j,k,n)*
     &              dabs(u(i,j,k,n))
               vdif(i,j,k)= vdif(i,j,k) -fac*v(i,j,k,n)*
     &              dabs(v(i,j,k,n))
               wdif(i,j,k)= wdif(i,j,k) -wfac*w(i,j,k,n)*
     &              dabs(w(i,j,k,n))
            end do
         end do
      end do

      return

c     MODEL 4
c     ======
c     CdN and CdT set in main.f
c     Cg = 0.3  in header
c     CdT and CdN are the normal and tangential drag coeffs Cd

C     Can modify grass diameter here if wanted
 502  do k=1,NKg
         func = 1.d0
c         func = (dble(NKg-k)+0.5)/dble(NKg)*0.8d0 +0.2d0
c         func = (dble(NKg-k)+0.5)/dble(NKg)
         dia(k) = CsDIA*func
      end do
c      delta = DL/LEN
c      delinv = LEN/DL
c     WL/UL = delta
c     UL/WL = delinv

      do k=1,NKg
         do j=1,NJ
            do i=1,NI
               facz = dia(k)*LEN*zNL(i)
C               if (i.gt.(NI-15)) then
C                  outdamp(i) = dexp(-dble((NI-i)*(NI-i))/16.d0)
C               else
C                  outdamp(i) = 0.d0
C               end if
c               write(6,*) 'outdampfac',i,outdamp(i)

c            if (i.le.20) then
c               wdampfac(i) = 1.d0 -dexp(-dble(i*i)/36.d0)
c            else
               wdampfac(i) = 1.d0
c            end if
c            write(6,*) 'wdampfac',i,wdampfac

c     to account for crowding.....
c     ........................
c     xc is dimensional (in m), xg is also dim (in m)
c.               denom = xc(i+1)+xg(i+1,j,k,n)-(xc(i)+xg(i,j,k))
c.               if (denom.gt.0.d0) then
c.                  rinv = (xc(i+1)-xc(i))/denom
c.               else
c.                  rinv = 10.d0
c.               end if
c     in order that zNL is not reduced by anything greater than a factor of 10.
c.               rinv = min(rinv, 10.d0)
c.
c     zNL is to be multiplied by rinv
c.               fac= rinv*facz*Jac(i,j,k)

C                if (i.gt.(NI-15)) then
C                   fac=0;
C                else
C                   fac= facz*Jac(i,j,k)
C                end if

               fac = facz*Jac(i,j,k)

               unorml= -u(i,j,k,n)*sinthe(i,j,k) +
     &              delta*w(i,j,k,n)*costhe(i,j,k)
               utangt= u(i,j,k,n)*costhe(i,j,k) +
     &              delta*w(i,j,k,n)*sinthe(i,j,k)
c     unorml and utangt are normalized by UL*UL
               udrag= -CdN*unorml*dabs(unorml)*sinthe(i,j,k)
     &              + CdT*utangt*dabs(utangt)*costhe(i,j,k)
               wdrag= CdT*utangt*dabs(utangt)*sinthe(i,j,k)
     &              + CdN*unorml*dabs(unorml)*costhe(i,j,k)
c     for now set vdrag is simplified
               vdrag= CdN*v(i,j,k,n)*dabs(v(i,j,k,n))
c
               udif(i,j,k)= udif(i,j,k) -fac*udrag
               vdif(i,j,k)= vdif(i,j,k) -fac*vdrag
c     &              - 0.05*outdamp(i)*v(i,j,k,n)
               wdif(i,j,k)= wdif(i,j,k) -fac*wdrag*delinv*wdampfac(i)
c     &              - 0.05*outdamp(i)*w(i,j,k,n)

            end do
         end do
      end do

C Deal with the actual damping
C       do i=1,NI
C          if (i.gt.(NI-15)) then
C             outdamp(i) = 0.5d0 * dexp(-dble((NI-i)*(NI-i))/16.d0)
C          else
C             outdamp(i) = 0.d0
C          end if
C       end do

C       do k=1,NK
C          do j=1,NJ
C             do i=1,NI

C                wdampfac(i) = 1.d0

C                udif(i,j,k)= udif(i,j,k) -
C      &                        outdamp(i)*(udif(i,j,k)-ufbcw(j,k))
C                vdif(i,j,k)= vdif(i,j,k) - outdamp(i)*(vdif(i,j,k)-0.d0)
C                wdif(i,j,k)= wdif(i,j,k) - outdamp(i)*(wdif(i,j,k)-0.d0)
C c     &              - 0.05*outdamp(i)*w(i,j,k,n)

C             end do
C          end do
C       end do

      return

      end
