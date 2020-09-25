      subroutine updatedh(step)
c     ==================================
c     Not doing anything as for now (dHfinal-dHinit)=0

      implicit double precision (a-h,o-z)
      include 'header.f'
      integer step,nsteps,ilocate
      double precision const,dH

      integer Niter,iter,i,j,k,itmp,jtmp,it
      double precision tol,maxres,res,rlx,ustar
      double precision cpf(3,NK),wzs(NK),wzks(0:NK),rhs(NK)
      double precision dia(NKg),temp,unormal,utangt,udrag,uxi,CsDIA
      double precision drpx(NI,NJ,NK),drpy(NI,NJ,NK),
     &     grpifc(0:NI,NJ,NK),grpjfc(NI,0:NJ,NK)
      common/rpgrads/drpx,drpy,grpifc,grpjfc

      double precision CgNL(0:NI+1)

      ilocate= NI/2
C      ilocate= NI/4
C      nsteps= 12000
      nsteps = 6000

      dH= dHinit + (dHfinal-dHinit)*dble(step)/dble(nsteps)
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


c=      if (step.eq.1) then
c     NOW UPDATE ufbcw - upstream flux
c     ===================================
      dH= (dHfinal-dHinit)/dble(nsteps)
c     This dH represents the change in height (from the previous time step)

      tol= 1.d-18

      j=NJ/2

      CsDIA = dsqrt(4.d0*CsA/PI)
      do k=1,NKg
c==         dia(k) = CsDIA*(((dble(NKg-k)+0.5)/dble(NKg) )*0.8d0 +0.2d0)
         dia(k) = CsDIA
c=         dia(k) = CsDIA*((dble(NKg-k)+0.5)/dble(NKg) )
      end do


c     dH is per LEN, i.e.dh=1.d-6 is 0.001 mm ht per 1 m, i.e. 1/million slope,
c     if LEN=1m

c     define zNL
      do i=0,NI+1
C         zNL(i) = 100.d0*(1.d0 - dble(i)/(NI+1))
        zNL(i) = 100.d0
      end do  

      WL = delta*UL
      REv = delta*UL*DL/vmu
      REVinv = vmu/(delta*UL*DL)
      CgNL = Cd*zNL
c      write(6,*) 'delta,UL,DL,vmu',delta,UL,DL,vmu
c      write(6,*) 'Re_V = ', REv ,'  CgNL = ',CgNL

c     wzks and wzs are constantk in x,y
      i=ilocate

      do k=1,NK
c         wzs(k) = 1.d0/dz
         wzs(k) = wz(ilocate,NJ/2,k)
      end do
      do k=0,NK
         wzks(k) = wzk(ilocate,NJ/2,k)
      end do


c     Initialize uprof(k)  to zero since this is the change in u
      do k=1,NK
c         uprof(k)= u(ilocate,NJ/2,k,0)
         uprof(k)= 0.d0
      end do

c     Find the coefficients of the elliptic operator and the RHS
c     RHS = dp/dx = G*dH*HL/(UL*UL)
      const = -GL*gpr*dH*HL/(UL*UL)
c     HL is just UL*UL/GL, so const is just = gpr*dH

      do k=1,NK
         rhs(k) = const
      end do

c     Coef cpf(1,k) is for k-1 cell, cpf(2,k) for k cell, cpf(3,k) for k+1 cell

      do k=2,NKg
         cpf(2,k) = -(wzks(k) + wzks(k-1))*wzs(k)*REVinv
     &         -CgNL(1)*uprof(k)*dia(k) ! using first stem for steady state
         cpf(1,k) = wzks(k-1)*wzs(k)*REVinv
         cpf(3,k) = wzks(k)*wzs(k)*REVinv
      end do

      do k=NKg,NK-1
         cpf(2,k) = -(wzks(k) + wzks(k-1))*wzs(k)*REVinv
         cpf(1,k) = wzks(k-1)*wzs(k)*REVinv
         cpf(3,k) = wzks(k)*wzs(k)*REVinv
      end do

c     Boundary conditions
c     Lower bndry, u=0, ie, uprof(0)= -uprof(1)
      k=1
      cpf(2,k)= -(wzks(k) +2.d0*wzks(k-1))*wzs(k)*REVinv
     &     -CgNL(1)*uprof(k)*dia(k)
      cpf(1,k) = 0.d0
      cpf(3,k) = wzks(k)*wzs(k)*REVinv

c     Upper bndry, du/dz =0
      k= NK
      cpf(2,k) = -(wzks(k-1)*wzs(k))*REVinv
      cpf(1,k) = wzks(k-1)*wzs(k)*REVinv
      cpf(3,k) = 0.d0

c     Use SOR to solve for uprof(k)

      rlx= 1.72
      itmp= ilocate
      jtmp= NJ/2
      do iter=1,1000

c     cpf(2,k) has to be recalculated at each iteration
c
         do k=2,NKg
c     calculate drag as in vdiffusion (no w vel)
            unorml= -uprof(k)*sinthe(itmp,jtmp,k)
            utangt= uprof(k)*costhe(itmp,jtmp,k)
c     unorml and utangt are normalized by UL*UL
            udrag= -CdN*unorml*dabs(unorml)*sinthe(itmp,jtmp,k)
     &           + CdT*utangt*dabs(utangt)*costhe(itmp,jtmp,k)
c            cpf(2,k)= -(wzks(k) +wzks(k-1))*wzs(k)*REVinv -CgNL*uprof(k)
c     &           *dia(k)
            cpf(2,k)= -(wzks(k) +wzks(k-1))*wzs(k)*REVinv -zNL(1)*LEN*
     &           udrag*dia(k)
         end do
         do k=NKg,NK-1
            cpf(2,k) = -(wzks(k) + wzks(k-1))*wzs(k)*REVinv
         end do

c     Boundary conditions
c     Lower bndry, u=0, ie, uprof(0)= -uprof(1)
         k=1
         unorml= -uprof(k)*sinthe(itmp,jtmp,k)
         utangt= uprof(k)*costhe(itmp,jtmp,k)
c     unorml and utangt are normalized by UL*UL
         udrag= -CdN*unorml*dabs(unorml)*sinthe(itmp,jtmp,k)
     &        + CdT*utangt*dabs(utangt)*costhe(itmp,jtmp,k)
         cpf(2,k)= -(wzks(k) +2.d0*wzks(k-1))*wzs(k)*REVinv
     &        -zNL(1)*LEN*udrag*dia(k)
c     &        -CgNL*uprof(k)*dia(k)

c     Upper bndry, du/dz =0 need not recalc cpf(2,k)

         do k=1,NK
            if (dabs(cpf(2,k)).le.1.d-15) then
               write(6,*) 'cpf(2) should be non-zero'
               write(6,*) cpf(2,k),k
               stop
            end if
         end do


         k=1
         ustar= ( cpf(3,k)*uprof(k+1)
     &        - rhs(k) )/(-cpf(2,k))
         uprof(k)= (1.d0 -rlx)*uprof(k) +rlx*ustar


         do k=2,NK-1
            ustar= ( cpf(3,k)*uprof(k+1)
     &           +cpf(1,k)*uprof(k-1)
     &           - rhs(k) )/(-cpf(2,k))
            uprof(k)= (1.d0 -rlx)*uprof(k) +rlx*ustar



         end do
         k=NK
         ustar= ( cpf(1,k)*uprof(k-1)
     &           - rhs(k) )/(-cpf(2,k))


         uprof(k)= (1.d0 -rlx)*uprof(k) +rlx*ustar


c     Compute residual
c     res = b - A x'
         maxres= 0.d0

         k=1
         res= rhs(k)-
     &        (cpf(2,k)*uprof(k)
     &        +cpf(3,k)*uprof(k+1))

         if (dabs(res).gt.maxres) then
            maxres= dabs(res)
         end if

         do k=2,NK-1
            res= rhs(k)-
     &           (cpf(1,k)*uprof(k-1)
     &           +cpf(2,k)*uprof(k)
     &           +cpf(3,k)*uprof(k+1))
            if (dabs(res).gt.maxres) then
               maxres= dabs(res)
            end if
         end do

         k=NK
         res= rhs(k)-
     &        (cpf(1,k)*uprof(k-1)
     &        +cpf(2,k)*uprof(k) )
         if (dabs(res).gt.maxres) then
            maxres= dabs(res)
         end if

c         if (mod(iter,100).eq.0) then
c            write(6,*) 'iter,maxres',iter,maxres
c         end if
         if (maxres.gt.5000.d0) then
            do k=1,NKg
               write(6,*) 'sin,cos,theta', sinthe(i,j,k),
     &              costhe(i,j,k),theta(i,j,k)
            end do

            write(6,*) 'STOP in updateh. Res too big, i,j,maxres=',
     &           i,j,maxres
            stop
         end if
         if (maxres.lt.tol) goto 13
      end do

 13   continue

      do k=1,NK
         ufprof(k)= Jac(ilocate,NJ/2,k)*ux(ilocate,NJ/2)*uprof(k)
      end do

c=      endif

      do k=1,NK
         ufss(k) = ufss(k) + ufprof(k)
         uss(k) = uss(k) + uprof(k)
         do j=1,NJ
            ufbcw(j,k)= ufbcw(j,k) + ufprof(k)
            uf(0,j,k) = ufbcw(j,k)
            do i=1,NI
               uf(i,j,k) = uf(i,j,k) + ufprof(k)
c            ufbce(j,k) will be set in the beg. of momentum
               u(i,j,k,0)= u(i,j,k,0)+ uprof(k)
            end do
         end do
      end do
c     Note that uss(k) is used in advec, but ufss is not used again.


      return
      end
