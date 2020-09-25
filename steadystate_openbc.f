      subroutine steadystate(Niter)
c     ============================

c     This routine calculates the steady state velocity profile uss(z)
c     as a function of the horizontal pressure gradient, viscosity,
c     grass drag, grass height, and fluid height

      implicit double precision (a-h,o-z)
      include 'header.f'
      double precision dia(NKg),temp,unormal,utangt,udrag,uxi,dH
      double precision drpx(NI,NJ,NK),drpy(NI,NJ,NK),
     &     grpifc(0:NI,NJ,NK),grpjfc(NI,0:NJ,NK)
      common/rpgrads/drpx,drpy,grpifc,grpjfc

c     NK= number of cells
c     NKg = number of cells in grass
c     DL = charac vertical ht
c     LEN = charac horizontal distance (m)
c     UL = charac horiz velocity   (m/s)
c     WL = charac vertical velocity (m/s)  = DL*UL/LEN = delta*UL
c     vmu = viscosity (in vertical)  (m2/s)
c     REv = vertical Reynolds number = WL*DL/vmu = delta*UL*DL/vmu
c     dH = delta free-surface (over distance LEN)
c     Cd = drag coeff of grass
c     zNL = number of grass blades/clusters per distance LEN
c     CgNL is non-dim = Cd*NL
c     dz = vertical grid spacing (non-dim by DL)
c     GL*gpr = accn due to gravity
c     CsDIA = charac stem dia

      integer Niter,iter,i,j,k,itmp,jtmp,it
      double precision tol,maxres,res,rlx,ustar

      double precision CgNL(0:NI+1)


c     Important that Cg*zNL approx =  dH*GL/(UL*UL)

      double precision cpf(3,NK),wzs(NK),wzks(0:NK),rhs(NK)
c     usual vmu=1.d-4
C     originally Amala had 1.d-4 here
      vmu= 1.d-4
C      vmu= 3.d-5
c      parameter (dH = 0.5d-2 )
      dH= dHinit

      tol= 1.d-12

      i=NI/2
      j=NJ/2

      CsDIA = dsqrt(4.d0*CsA/PI)
      do k=1,NKg
c==         dia(k) = CsDIA*(((dble(NKg-k)+0.5)/dble(NKg) )*0.8d0 +0.2d0)
         dia(k) = CsDIA
c=         dia(k) = CsDIA*((dble(NKg-k)+0.5)/dble(NKg) )
      end do

c     define zNL and CgNL
      do i=0,NI+1
C         zNL(i) = 100.d0*(1.d0 - dble(i)/(NI+1))
        zNL(i) = 100.d0
        CgNL(i) = Cd*zNL(i)
      end do  

c     dH is per LEN, i.e.dh=1.d-6 is 0.001 mm ht per 1 m, i.e. 1/million slope,
c     if LEN=1m

c      write(6,*) 'DL',DL,'LEN',LEN,'delta',delta
      WL = delta*UL
      REv = delta*UL*DL/vmu
      REVinv = vmu/(delta*UL*DL)

c      write(6,*) 'delta,UL,DL,vmu',delta,UL,DL,vmu
c      write(6,*) 'Re_V = ', REv ,'  CgNL = ',CgNL

c     wzks and wzs are constantk in x,y

      do k=1,NK
c         wzs(k) = 1.d0/dz
         wzs(k) = wz(NI/2,NJ/2,k)
      end do
      do k=0,NK
         wzks(k) = wzk(NI/2,NJ/2,k)
      end do

c     Initialize uss(k) to zero
      do k=1,NK
         uss(k)= u(NI/2,NJ/2,k,0)
      end do

c     Find the coefficients of the elliptic operator and the RHS

c      if (Niter.gt.20) then
c         write(6,*) 'Niter must be 20 in steadystate'
c         stop
c      end if
c     RHS = dp/dx = G*dH*HL/(UL*UL)
      const = -GL*gpr*dH*HL/(UL*UL)
c     HL is just UL*UL/GL, so const is just = gpr*dH

      do k=1,NK
         rhs(k) = const
      end do

c     Coef cpf(1,k) is for k-1 cell, cpf(2,k) for k cell, cpf(3,k) for k+1 cell

      do k=2,NKg
         cpf(2,k) = -(wzks(k) + wzks(k-1))*wzs(k)*REVinv -CgNL(1)*uss(k)
     &        *dia(k)
         cpf(1,k) = wzks(k-1)*wzs(k)*REVinv
         cpf(3,k) = wzks(k)*wzs(k)*REVinv
      end do

      do k=NKg+1,NK-1
         cpf(2,k) = -(wzks(k) + wzks(k-1))*wzs(k)*REVinv
         cpf(1,k) = wzks(k-1)*wzs(k)*REVinv
         cpf(3,k) = wzks(k)*wzs(k)*REVinv
      end do

c     Boundary conditions
c     Lower bndry, u=0, ie, uss(0)= -uss(1)
      k=1
      cpf(2,k)= -(wzks(k) +2.d0*wzks(k-1))*wzs(k)*REVinv -CgNL(1)*uss(k)
     &     *dia(k)
      cpf(1,k) = 0.d0
      cpf(3,k) = wzks(k)*wzs(k)*REVinv

c     Upper bndry, du/dz =0
      k= NK
      cpf(2,k) = -(wzks(k-1)*wzs(k))*REVinv
      cpf(1,k) = wzks(k-1)*wzs(k)*REVinv
      cpf(3,k) = 0.d0

c     Use SOR to solve for uss(k)

      rlx= 1.72
      itmp= NI/2
      jtmp= NJ/2
      do iter=1,1000

c     cpf(2,k) has to be recalculated at each iteration
c
         do k=2,NKg
c     calculate drag as in vdiffusion (no w vel)
            unorml= -uss(k)*sinthe(itmp,jtmp,k)
            utangt= uss(k)*costhe(itmp,jtmp,k)
c     unorml and utangt are normalized by UL*UL
            udrag= -CdN*unorml*dabs(unorml)*sinthe(itmp,jtmp,k)
     &           + CdT*utangt*dabs(utangt)*costhe(itmp,jtmp,k)
c            cpf(2,k)= -(wzks(k) +wzks(k-1))*wzs(k)*REVinv -CgNL*uss(k)
c     &           *dia(k)
            cpf(2,k)= -(wzks(k) +wzks(k-1))*wzs(k)*REVinv -zNL(1)*LEN*
     &           udrag*dia(k)
         end do
         do k=NKg+1,NK-1
            cpf(2,k) = -(wzks(k) + wzks(k-1))*wzs(k)*REVinv
         end do

c     Boundary conditions
c     Lower bndry, u=0, ie, uss(0)= -uss(1)
         k=1
         unorml= -uss(k)*sinthe(itmp,jtmp,k)
         utangt= uss(k)*costhe(itmp,jtmp,k)
c     unorml and utangt are normalized by UL*UL
         udrag= -CdN*unorml*dabs(unorml)*sinthe(itmp,jtmp,k)
     &        + CdT*utangt*dabs(utangt)*costhe(itmp,jtmp,k)
         cpf(2,k)= -(wzks(k) +2.d0*wzks(k-1))*wzs(k)*REVinv
     &        -zNL(1)*LEN*udrag*dia(k)
c     &        -CgNL*uss(k)*dia(k)

c     Upper bndry, du/dz =0 need not recalc cpf(2,k)

         do k=1,NK
            if (dabs(cpf(2,k)).le.1.d-15) then
               write(6,*) 'cpf(2) should be non-zero'
               write(6,*) cpf(2,k),k
               stop
            end if
         end do


         k=1
         ustar= ( cpf(3,k)*uss(k+1)
     &        - rhs(k) )/(-cpf(2,k))
         uss(k)= (1.d0 -rlx)*uss(k) +rlx*ustar


         do k=2,NK-1
            ustar= ( cpf(3,k)*uss(k+1)
     &           +cpf(1,k)*uss(k-1)
     &           - rhs(k) )/(-cpf(2,k))
            uss(k)= (1.d0 -rlx)*uss(k) +rlx*ustar



         end do
         k=NK
         ustar= ( cpf(1,k)*uss(k-1)
     &           - rhs(k) )/(-cpf(2,k))


         uss(k)= (1.d0 -rlx)*uss(k) +rlx*ustar


c     Compute residual
c     res = b - A x'
         maxres= 0.d0

         k=1
         res= rhs(k)-
     &        (cpf(2,k)*uss(k)
     &        +cpf(3,k)*uss(k+1))

         if (dabs(res).gt.maxres) then
            maxres= dabs(res)
         end if

         do k=2,NK-1
            res= rhs(k)-
     &           (cpf(1,k)*uss(k-1)
     &           +cpf(2,k)*uss(k)
     &           +cpf(3,k)*uss(k+1))
            if (dabs(res).gt.maxres) then
               maxres= dabs(res)
            end if
         end do

         k=NK
         res= rhs(k)-
     &        (cpf(1,k)*uss(k-1)
     &        +cpf(2,k)*uss(k) )
         if (dabs(res).gt.maxres) then
            maxres= dabs(res)
         end if

         if (mod(iter,100).eq.0) then
            write(6, '(A, I6, E11.4)') 'iter,maxres',iter,maxres
         end if
         if (maxres.gt.5000.d0) then
            do k=1,NKg
               write(6,*) 'sin,cos,theta', sinthe(i,j,k),
     &              costhe(i,j,k),theta(i,j,k)
            end do

            write(6,*) 'STOP in steadystate. Res too big, i,j,maxres=',
     &           i,j,maxres
            stop
         end if
         if (maxres.lt.tol) goto 13
      end do

 13   write(6,*) 'steadstate converged with iter=',iter,maxres
c=      write(6,*) 'write velocities'

c     Lower boundary uss(z=0)= 0
c      uss(0) = -uss(1)
c     Upper boundary du/dz = 0
c      uss(NK+1)= uss(NK)

      do k=1,NK
         ufss(k)= Jac(NI/2,NJ/2,k)*ux(NI/2,NJ/2)*uss(k)
         do j=1,NJ
            ufbcw(j,k)= ufss(k)
            ufbce(j,k)= ufss(k)   ! This will be modified with time.
            do i=0,NI
               uf(i,j,k) = ufss(k)
            end do
         end do
      end do

c     Initialize u(0)
c     -------------------
      do k=1,NK
         do j=0,NJ+1
            do i=0,NI+1
               u(i,j,k,0)= uss(k)
            end do
         end do
      end do

c      do k=1,NK
c         write(6,*) k,zc(NI/2,NJ/2,k)*DL,uss(k)*UL
c      end do
c      stop
c     Compute drpx and grpifc which stay constant throughout (this is the
c     background pressure gradient
c     ---------------------------
      const= -gpr*dH
      do k=1,NK
         do j=1,NJ
            do i=1,NI
               drpx(i,j,k)= const
               drpy(i,j,k)= 0.d0
               si(i,j,k)= drpx(i,j,k)
               sj(i,j,k)= drpy(i,j,k)
            end do
         end do
         do j=1,NJ
            do i=0,NI
               uxi= 0.5*(ux(i,j) + ux(i+1,j))
               grpifc(i,j,k)= Jifc(i,j,k)*uxi*const
c               if ((j.eq.4).and.(k.eq.7)) write(6,*)
c     &              Jifc(i,j,k),uxi,grpifc(i,j,k)
            end do
         end do
         do j=0,NJ
            do i=1,NI
               grpjfc(i,j,k)= 0.d0
            end do
         end do
      end do


      goto 201

c     Sin wave to u   or random perturb

      dum = ran3(iseed)
      halfL = 0.5*( 0.5*(xc(NI+1)+xc(NI)) - 0.5*(xc(1)+xc(0)) )
      halfLinv = 1.d0/halfL
      xcenter = 0.5*(xc(NI/2) +xc(NI/2+1))
      do i=1,NI
         pulse = dSin (0.5d0*PI*(1.d0 - dabs(xc(i)-xcenter)*halfLinv) )
c         write(6,*) 'pulse',pulse
         temp= 0.0001*UL*ran3(iseed)
         do k=1,NK
            do j=0,NJ+1
c               u(i,j,k,0)= uss(k) + 0.25*pulse*uss(k)
c               u(i,j,k,0)= uss(k) + 0.025*pulse*uss(k)
               u(i,j,k,0)= uss(k) + temp
            end do
         end do
      end do

c     Initialize outflow  (not needed)
c 201  do k=1,NK
c         uout(k)= uss(k)
c         vout(k) = 0.d0
c         wout(k) = 0.d0
c         sout(k) = s(NI,NJ/2,k,0)
c         do it=1,ntr
c            Tout(it,k) = T(it,NI,NJ/2,k,0)
c         end do
c      end do

 201  return
c      do i=0,NI
c         uxi= 0.5*(ux(i,j) + ux(i+1,j))
c         write(6,*) grpifc(i,4,7),Jifc(i,4,7),uxi
c      end do
      end
