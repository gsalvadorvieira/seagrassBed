      program main
c---------------------------------------------------
c     test model
      implicit logical (a-z)
c
      include 'header.f'
c     parameters:
c     NI,NJ,NK are the no. of points in each of the 3 co-ordintate directions.
c     L, DL are the characteristic length and depth scales.
c     AH and AV are the eddy diffusivities in the horizontal and vertical dir
c     lexp is a constant in the equations denoting the length scale

c     dtime is the time step
c     S0 and T0 are the mean salinity and temp. sreal= S0+s,Treal=T0+T
c     dtf is the(dimensionless) time step
c     delta is the ratio of vertical to horizontal length scale = D/L
c     delinv= 1/delta

c     kappah is the implicitness parameter in the Crank-Nicolson scheme for
c     the h eqn. kaphinv is its inverse.
c     set kappah=1/2.  kaphinv=1/kappah

      integer step,nsteps,n,i,j,k,frame_int,ngraph,iter,niter
      double precision fdiv,pcorr(maxout),ctrdiv,hsum,hmean
      double precision pgrad(0:NI+1,0:NJ+1,0:NK+1),
     &     hgrad(0:NI+1,0:NJ+1,0:NK+1),bgrad(0:NI+1,0:NJ+1,0:NK+1)
      real dtime, tarray(2), tim
      open(unit=31, file='gulf.in',status='old')
      read(31,*)
      read(31,*) nsteps,dtf,ngraph,frame_int
      read(31,*)
      read(31,*) dx,dy
      close(31)
      write(6,*) 'nsteps= ',nsteps
c     define the parameters
c      frame_int= 10

      delta= DL/LEN
      delinv= LEN/DL
      qpr= 1.d0
c
c     Kappah is the implicitness parameter for h (1=fully implicit, 0=explicit)
      kappah= 0.65d0
      kaphinv= 1.d0/kappah
c
c     UL= 0.1 in header
c     UL in m/s
      TL= LEN/UL
c     EPS= 1.d0  (in header) EPS = (UL*UL)/(GL*HL) = inverse Froude number
c     HL is such that the Fr number GL*HL/(UL*UL) = 1
      HL=UL*UL/GL
c     HL=1.d-3 for UL=0.1
c     HDL= 1.d-3
      HDL= HL/DL
c
c     grass: drag coeff in normal and tgt directions
c     -----
      CdN = Cd
      CdT = 0.d0
C      CdT = 0.25d0*Cd

      write(6,*) 'Parameters'
      write(6,*) 'delta=',delta,'dtime=',dtf,
     &     'LEN ', LEN, ' DL ', DL
      write(6,*) 'R0, UL, P1, HL, HL/DL', R0,UL,P1,HL,HDL
c
      call init(pcorr)
      write(6,*) 'init called'
c      call meanh(NI,NJ,h,hmean)
c      write(6,*) 'initial  hmean',hmean

      call sigma

c     Important to call sigma before steady state

      call steadystate(1)
c=      call grassbend(0)
      call findzall
      call sigma
      do iter = 1,10
         call steadystate(iter)   ! Finds 1-D flow profile and initializes flow
c     diff grassshape solicit diff solutions particulary when EI= 10^-2
c     call grassshape after steadystate
         call grassshape(0) ! computes grass deformation and ght
         write(6,*) 'iter ',iter, 'write out grassshape'
C         do k=1,NKg
C            write(6,*) 'xg,zg',k,xg(1,1,k,0),zg(1,1,k,0),theta(1,1,k)
C         end do
         call findzall ! computes z position from new grass height ght
         call sigma
      end do
      write(6,*) 'steady profiles computed'

      call hsave ! saves the hydrostatic pressure gradients at the initial time step
      call facediv(EPS,fdiv) ! checks divergence
      call cdiv(EPS,ctrdiv,0) ! checks divergence using u,v,w at cell centers
C       write(6,*) ' fdiv  ',fdiv,' cdiv  ',ctrdiv
      write(6,'(2(A,ES13.6))') 'fdiv  ',fdiv,', cdiv  ',ctrdiv
c     ===========
c     write out initial values
      call vort(0)
      call outcdf(xc,yc,zc,p,h,s,T,u,v,w,vor,uf,vf,wf,theta,xg,
     &     zg,ght,gdef,0)

      step= 0
      call vort(0) ! double call?

      call writeksurf(frame_int,step,24,s,T,rho,u,v,w,vor,xc,yc,zc)
      call writeslice(frame_int,step,s,T,rho,u,v,w,vor,yc,zc)
      call calcpgrad(pgrad,hgrad,bgrad)
      call writeyslice(frame_int,step,s,T,rho,u,v,w,p,vor,h,ght,gdef,
     &     xg,zg,xc,zc,pgrad,hgrad,bgrad,theta)
c
      write(6,*) 'To start momentum ...'
      write(6,*)
c=      tim= dtime(tarray)
      call checks ! Checks and writes a few things, but nothing crucial

c     ===========
c     start time steps

 50   step=step+1 ! MAIN LOOP IS HERE: goto 50 until step.eq.nsteps

C     Force v=0 (2D) - 7/27/2020
      do i=0,NI+1
         do j=0,NJ+1
            do k=0,NK+1
               v(i,j,k,0)=0.d0
            end do
         end do
      end do

C      call updatedh(step) ! No need to call because dh is time invariant
      call momentum(pcorr,step) ! Most of the code is here. Solves the eqs using a RK-3 scheme
c     compute the KE
      call energy(step)
      call meanh(NI,NJ,h,hmean)

      if (mod(step,frame_int).eq.0) then
         call vort(0)
c         call writeframe(frame_int,step,h,s,T,rho,u,v,w,vor,xc,yc,zc)
         call writeksurf(frame_int,step,24,s,T,rho,u,v,w,vor,xc,yc,zc)
         call writeslice(frame_int,step,s,T,rho,u,v,w,vor,yc,zc)
         call calcpgrad(pgrad,hgrad,bgrad)
         call writeyslice(frame_int,step,s,T,rho,u,v,w,p,vor,h,ght,
     &     gdef,xg,zg,xc,zc,pgrad,hgrad,bgrad,theta)
      endif

      write(6,'(A,I8,A,ES13.6)') '--- STEPS: ',step, ', hmean = ',hmean
      write(6,*)

c      if (mod(step,100).eq.0) then
c         tim= dtime(tarray)
c         write(6,*) 'cpu time for 100 steps=', tim
c      endif

      if (step.eq.nsteps) goto 100
      if (mod(step,ngraph).eq.0) goto 100
      if (step.eq.50) goto 100
      if (step.eq.100) goto 100
      if (step.eq.250) goto 100
      if (step.eq.500) goto 100
c      sigma is called before next time steps to re-evaluate metrics.
c      call sigma  - called in momentum
      goto 50

c     ===========
c     write out when appropriate
 100  n=0
      call vort(0)
      call outcdf(xc,yc,zc,p,h,s,T,u,v,w,vor,uf,vf,wf,theta,xg,zg,ght,
     &     gdef,step)

C     finish run
      if (step.lt.nsteps) goto 50
         stop
      end
