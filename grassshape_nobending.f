      subroutine grassshape(n)
c----------------------------------------------------
c     We calculate the coordinates xg,zg (x,z coordinates) of a "single"
c     plant stem which extends from s=0 to s=l.
c     Modified for open bc
c     No bending stiffness -- buoyancy only
c     We assume that the grass density is inconsequential to its bending.
c     ---------------------------------------------------------------

      implicit double precision (a-h,o-z)
      include 'header.f'

      integer i,j,k,n

      double precision dia(NKg),csarea(NKg),
     &     Tdzds(0:NKg),Tdxds(0:NKg),intTdzds,intTdxds,WL
     &     alphax,alphaz,alpha,unorml,utangt,udrag,vdrag,wdrag

c     drho = density diff between water and grass (grass in buoyant)
c     CsA = charac. cross sectional area of grass stem (1cm2)
c     B0= charac buoyancy force per unit length (kg/s^2)
c     NKg = number of cells in vertical
c     UgL = charac velocity in grass = 0.1*UL
c     Grasslen = grass length non-dim by DL, dsprime= Grasslen/NKg
c     CdT and CdN are the normal and tangential drag coeffs Cd

C     AUTOMATICALLY IMPORTED PARAMETER
      parameter(drho=50.d0)
C       parameter (drho=25.d0)
C       parameter (drho=990.d0)

c     delta = DL/LEN
c     delta = WL/UL
      CsDIA = dsqrt(4.d0*CsA/PI)
      Ug = 0.1d0*UL
      dsprime= Grasslen/dble(NKg)

c     for inflex grass      goto 2000

c     This can be changed for varying thickness of grass stems
c     -------------------------------------------------------
      do k=1,NKg
         dia(k)= 1.d0
c         dia(k) = ((dble(NKg-k)+0.5)/dble(NKg) )*0.8d0 + 0.2d0
c         dia(k) = ((dble(NKg-k)+0.5)/dble(NKg) )

         csarea(k)= dia(k)*dia(k)
      end do

      B0 = drho*GL*gpr*CsA
c     alpha = R0*Ug*Ug*CsDIA* Cd / B0  ~1
c     Therefore  Cd approx =  alpha*B0 /(R0*Ug*Ug*CsDIA) = 0.03

c     Cd = 0.3 (same as Cd in header)
c     CdT and CdN are the normal and tangential drag coeffs Cd set in main.f

      alpha = R0*Ug*Ug*CsDIA*CdN/B0

      if (dabs(alpha-1.d0).gt.30.d0) then
         write(6,*) 'In grassbend, alpha is much larger than 1'
         write(6,*) 'alpha = ', alpha
      endif

c     use UL and WL since these alphax and alphaz multiply u*u and w*w

      alphax = R0*UL*UL*CsDIA/B0
      WL= delta*UL
      alphaz = R0*WL*WL*CsDIA/B0

      do iter=1,3
c      write(6,*) 'in grassbend,  iter =  ',iter


      do j=1,NJ
         do i=1,NI

c     Tdzds and Tdxds are at cell faces

c     Since the tension, T is 0 at the stem tip
            Tdzds(NKg) = 0.d0
            Tdxds(NKg) = 0.d0

            do k=NKg,1,-1

c     rhoprime = (R0 + rho(i,j,k))/R0
c     for constant density,
               rhoprime = 1.d0
               unorml= -u(i,j,k,n)*sinthe(i,j,k) +
     &              delta*w(i,j,k,n)*costhe(i,j,k)
               utangt= u(i,j,k,n)*costhe(i,j,k) +
     &              delta*w(i,j,k,n)*sinthe(i,j,k)
c     unorml and utangt are normalized by UL
               udrag= -CdN*unorml*dabs(unorml)*sinthe(i,j,k)
     &              + CdT*utangt*dabs(utangt)*costhe(i,j,k)
               wdrag= CdT*utangt*dabs(utangt)*sinthe(i,j,k)
     &              + CdN*unorml*dabs(unorml)*costhe(i,j,k)
c     for now, vdrag is not used
               vdrag= 0.d0

c     bprime is the buoyancy force per unit length normalized by B0
               bprime= csarea(k)

c     alphax = R0*UL*UL*CsDIA/B0
               Tdzds(k-1)= Tdzds(k) +dsprime *(
     &              alphax*wdrag*rhoprime*dia(k)
     &              + bprime )

               Tdxds(k-1)= Tdxds(k) +dsprime *(
     &              alphax*udrag*rhoprime*dia(k) )

            end do

            do k=1,NKg

               fac = dsqrt ( Tdzds(k-1)*Tdzds(k-1) + Tdxds(k-1)*
     &              Tdxds(k-1) )
               dzdskm1 = Tdzds(k-1) /fac
               dxdskm1 = Tdxds(k-1) /fac

               costhe(i,j,k)= dxdskm1
               sinthe(i,j,k)= dzdskm1
               theta(i,j,k)= datan2(sinthe(i,j,k),costhe(i,j,k))
               if (k.eq.1) then
                  zg(i,j,k,n)= dzdskm1*0.5*dsprime
                  xg(i,j,k,n)= dxdskm1*0.5*dsprime
               else
                  xg(i,j,k,n)= xg(i,j,k-1,n) + dxdskm1*dsprime
                  zg(i,j,k,n)= zg(i,j,k-1,n) + dzdskm1*dsprime

               endif
c               write(6,*) 'xg,zg',k,xg(i,j,k,n),zg(i,j,k,n)
c	       if (i.eq.1) then
c		write(6,*) 'i,j,k,theta',i,j,k,theta(i,j,k)
c	       endif
            end do
c            stop
            k= NKg
            intTdzds= 0.25*Tdzds(k-1) + 0.75*Tdzds(k)
            intTdxds= 0.25*Tdxds(k-1) + 0.75*Tdxds(k)
            fac = dsqrt ( intTdzds*intTdzds + intTdxds*
     &              intTdxds )
            dzdskm1 = intTdzds /fac
            dxdskm1 = intTdxds /fac

            ght(i,j) = zg(i,j,NKg,n) + dzdskm1*0.5d0*dsprime
            gdef(i,j) = xg(i,j,NKg,n) + dxdskm1*0.5d0*dsprime

c            write(6,*) 'i = ', i , 'j = ', j
c            write(6,*) 'Tdzds = ', (Tdzds(k),k=0,NKg)
c            write(6,*) 'Tdxds = ', (Tdxds(k),k=0,NKg)
c            write(6,*) 'ght,gdef',ght(i,j),gdef(i,j)
c            stop
         end do
      end do
      end do

 2000 do i=1,NI
         ght(i,0)= ght(i,1)
         ght(i,NJ+1)= ght(i,NJ)
         gdef(i,0)= gdef(i,1)
         gdef(i,NJ+1)= gdef(i,NJ)
         do k=1,NKg
            xg(i,0,k,n)= xg(i,1,k,n)
            xg(i,NJ+1,k,n)= xg(i,NJ,k,n)
            zg(i,0,k,n)= zg(i,1,k,n)
            zg(i,NJ+1,k,n)= zg(i,NJ,k,n)
         end do
      end do


c     periodic-ew boundaries
      if (periodicew) then
      do j=0,NJ+1
         ght(0,j)= ght(NI,j)
         ght(NI+1,j)= ght(1,j)
         gdef(0,j)= gdef(NI,j)
         gdef(NI+1,j)= gdef(1,j)
         do k=1,NKg
            xg(0,j,k,n)= xg(NI,j,k,n)
            xg(NI+1,j,k,n)= xg(1,j,k,n)
            zg(0,j,k,n)= zg(NI,j,k,n)
            zg(NI+1,j,k,n)= zg(1,j,k,n)
         end do
      end do
      else
      do j=0,NJ+1
         ght(0,j)= ght(1,j)
         ght(NI+1,j)= ght(NI,j)
         gdef(0,j)= gdef(1,j)
         gdef(NI+1,j)= gdef(NI,j)
         do k=1,NKg
            xg(0,j,k,n)= xg(1,j,k,n)
            xg(NI+1,j,k,n)= xg(NI,j,k,n)
            zg(0,j,k,n)= zg(1,j,k,n)
            zg(NI+1,j,k,n)= zg(NI,j,k,n)
            costhe(NI+1,j,k)= costhe(NI,j,k)
            costhe(0,j,k)= costhe(1,j,k)
            sinthe(NI+1,j,k)= sinthe(NI,j,k)
            sinthe(0,j,k)= sinthe(1,j,k)
            theta(NI+1,j,k)= theta(NI,j,k)
            theta(0,j,k)= theta(1,j,k)
         end do
      end do
      endif
c
      do j=0,NJ+1
         do i=0,NI+1
            Dg(i,j)= D(i,j) - ght(i,j)
c     ght is non-dim by DL and positive
c     D(i,j) is positive and non-dim by DL .
c     Dg(i,j) (depth of grass) is positive and non-dim by DL

         end do
      end do

      return

c     CONVERT xg,zg to dimensional form (meters) by xlying by DL
      do k=1,NKg
         do j=0,NJ+1
            do i=0,NI+1
               xg(i,j,k,n)= xg(i,j,k,n)*DL
               zg(i,j,k,n)= zg(i,j,k,n)*DL
            end do
         end do
      end do

      return
      end
