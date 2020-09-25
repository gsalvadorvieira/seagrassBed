      subroutine grassshape(n)
c----------------------------------------------------
c     In this version, i.e. grassshape_inflexible, we prevent the grass
c     from bending.
c     theta is the angle with the horizontal, so theta =90deg, sinthe=1,
c     costhe=0.d0
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

c     delta = DL/LEN (specified in main.f)
      CsDIA = dsqrt(4.d0*CsA/PI)
      dsprime= Grasslen/dble(NKg)

c     This can be changed for varying thickness of grass stems
c     -------------------------------------------------------
      do k=1,NKg
         dia(k)= 1.d0
         csarea(k)= dia(k)*dia(k)
      end do

      B0 = drho*GL*gpr*CsA

c     Cd = 0.3 (same as Cd in header)
c     CdT and CdN are the normal and tangential drag coeffs Cd set in main.f

c     use UL and WL since these alphax and alphaz multiply u*u and w*w

      alphax = R0*UL*UL*CsDIA/B0

c=inflex      do iter=1,3
c      write(6,*) 'in grassbend,  iter =  ',iter


      do j=1,NJ
         do i=1,NI
            do k=1,NKg

               theta(i,j,k)= pi/2.d0
               costhe(i,j,k)= 0.d0
               sinthe(i,j,k)= 1.d0
               dzdskm1= 1.d0
               dxdskm1= 0.d0
               if (k.eq.1) then
                  zg(i,j,k,n)= dzdskm1*0.5*dsprime
                  xg(i,j,k,n)= dxdskm1*0.5*dsprime
               else
                  xg(i,j,k,n)= xg(i,j,k-1,n) + dxdskm1*dsprime
                  zg(i,j,k,n)= zg(i,j,k-1,n) + dzdskm1*dsprime
               endif
c               write(6,*) 'xg,zg',k,xg(i,j,k,n),zg(i,j,k,n)

            end do
            
            dzdskm1= 1.d0
            dxdskm1= 0.d0
            ght(i,j) = zg(i,j,NKg,n) + dzdskm1*0.5d0*dsprime
            gdef(i,j) = xg(i,j,NKg,n) + dxdskm1*0.5d0*dsprime

         end do

      end do

      do i=1,NI
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
