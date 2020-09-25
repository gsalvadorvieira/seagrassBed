      subroutine mgrid(p,dtime,edt,cfcdiv)
      implicit logical (a-z)
c     --------------
c     multigrid solver for the elliptic equation Del2 p= f
c     Run pgrogram "preprocess" to get the value of maxout
c     input :
c     -----
c     call the subroutine mgrid with p_initial on fine grid
c     maxout : max dimension on the one-dim array - includes outer points
c     maxint: max dimension of the one-dim array for variables at interior
c              points of grid
c     int1 :  number of internal points on the  fine grid
c     ngrid : number of grid levels
c     nx(m),ny(m),nz(m) = number of int grid points at level m. m=1...ngrid
c     m=1 is the finest grid, m=ngrid refers to the coarsest grid.
c     ntout(m) is the total number of grid points (storage locations) at
c               grid level m - including the outer ficticious points
c     ntint(m) is the total number of interior grid points (storage locations)
c               at grid level m
c     loco(m) is the storage location of a scalar variable in a one-dim array
c             of dimension  Summation(m=1,ngrid) {ntout(m)}
c     loci(m) is the storage location of a scalar specified only at interior
c             grid points
c     loccp(m) is the storage location for cp on the m-th level
c     res(i)  temporary stores the residual before it interpolated onto
c             the coarse grid
c     dxf,dyf,dzf  are the grid spacings on the fine grid

      integer ngrid,m,maxout,maxint,int1,l,ncycle,ll,
     &     nu1,nu2,noc

C     AUTOMATED UPDATE
      parameter(ngrid=3,maxout=257168,maxint=189216,int1=165888)

c     96x8x48 (grassbed)
C       parameter(ngrid=3,maxout=58256,maxint=42048,int1=36864)

c     192x8x48 (grassbed)
C      parameter(ngrid=3,maxout=115088,maxint=84096,int1=73728)

c     288x8x48 (grassbed)
C      parameter(ngrid=3,maxout=171920,maxint=126144,int1=110592)

cc     48x32x16 grid
cc      parameter(ngrid=4,maxout=36312,maxint=28080,int1=24576)
c     48x24x32 grid
c      parameter(ngrid=4,maxout=13272,maxint=9360,int1=8192)
c-      parameter(ngrid=5,maxout=46416,maxint=37448,int1=32768)
c1k     96x96x24
c      parameter(ngrid=4,maxout=291092,maxint=252720,int1=221184)
c1ktest   16x96x24
c      parameter(ngrid=4,maxout=54392,maxint=42120,int1=36864)
c.5k  192x96x32
c      parameter(ngrid=5,maxout=750240,maxint=674064,int1=589824)
c1k   48x96x24
c      parameter(ngrid=4,maxout=149072,maxint=126360,int1=110592)
c.25k  192x384x24  (insuff. memory on vayu)
c      parameter(ngrid=4,maxout=2258852,maxint=2021760,int1=1769472)
c125m  96x192x24
c       parameter(ngrid=4,maxout=575132,maxint=505440,int1=442368)
ctest125m  16x96x24
c      parameter(ngrid=4,maxout=54392,maxint=42120,int1=36864)
c.5k  192x96x16
c      parameter(ngrid=4,maxout=400472,maxint=336960,int1=294912)
c-     48x48x24
csq      parameter(ngrid=4,maxout=76352,maxint=63180,int1=55296)
c     32x32x16
c      parameter(ngrid=4,maxout=24792,maxint=18720,int1=16384)
c=32x64x16      parameter(ngrid=4,maxout=47832,maxint=37440,int1=32768)
c=      parameter(ngrid=4,maxout=92312,maxint=74880,int1=65536)
c      parameter(ngrid=4,maxout=92312,maxint=74880,int1=65536)

c=     96x8x72 (deep grassbed)
c=      parameter(ngrid=3,maxout=86000,maxint=63072,int1=55296)
c     96x8x36 (shallow water)
c      parameter(ngrid=3,maxout=44384,maxint=31536,int1=27648)
c=     96x8x60 (deeper)
c      parameter(ngrid=3,maxout=72128,maxint=52560,int1=46080)
c     96x8x56 (deeper)
c      parameter(ngrid=3,maxout=67504,maxint=49056,int1=43008)
c     96x8x52 (deeper)
c      parameter(ngrid=3,maxout=62880,maxint=45552,int1=39936)
c     96x48x48
c      parameter(ngrid=5,maxout=284992,maxint=252774,int1=221184)
c     192x8x48 (long grassbed)
c      parameter(ngrid=3,maxout=115088,maxint=84096,int1=73728)

      integer ntout(ngrid),ntint(ngrid),nx(ngrid),ny(ngrid),nz(ngrid),
     &     loco(ngrid),loci(ngrid),loccp(ngrid)
      double precision cp(19*maxint),rhs(maxint),p(maxout),
     &     res(int1),dtime,tol,maxres,edt,ratio,oldres,cfcdiv      


c     edt= EPS/dtime. If the tolerance is on (u_x +v_y +eps*w_z)
c     then  the tolerance on the residual r= edt*(u_x+ v_y + eps*w_z)
c     is equal to edt*(the specified tolerance)
c
      open (unit=90, file='mg.in')
      read(90,*) nx(1),ny(1),nz(1),noc,nu1,nu2,tol
      close(90)

      if (cfcdiv.le.tol) then
         do m=1,maxout
            p(m)= 0.d0
         end do
         maxres = edt*cfcdiv
         ncycle = 0
         goto 101
      end if

c     redefine the tolerance as tol*edt
      tol= edt*tol
      ntint(1)= nx(1)*ny(1)*nz(1)
      ntout(1)= (nx(1)+2)*(ny(1)+2)*(nz(1)+2)
c
      loco(1)= 1
      loci(1)= 1
      loccp(1)= 1
      do m=2,ngrid
         if (mod(nx(m-1),2).ne.0) goto 20
         if (mod(ny(m-1),2).ne.0) goto 20
         if (mod(nz(m-1),2).ne.0) goto 20
         nx(m)= nx(m-1)/2
         ny(m)= ny(m-1)/2
         nz(m)= nz(m-1)/2
         ntint(m)= nx(m)*ny(m)*nz(m)
         ntout(m)= (nx(m)+2)*(ny(m)+2)*(nz(m)+2)
         loco(m)= loco(m-1) + ntout(m-1)
         loci(m)= loci(m-1) + ntint(m-1)
c     for the 19 point stencil
         loccp(m)= loccp(m-1) + 19*ntint(m-1)
      end do
c
c     compute fine grid coefficients and rhs on fine grid
      call cpfine(dtime,cp,rhs)
c      call checkcp(1,nx(1),ny(1),nz(1),cp(loccp(1)))
c      write(6,*) 'cpfine checked'
c     **********************************
      do m=2,ngrid
         call cpcors(nx(m),ny(m),nz(m),cp(loccp(m-1)),cp(loccp(m)))
c         call checkcp(m,nx(m),ny(m),nz(m),cp(loccp(m)))
      end do
c
      do 1000 ncycle=1,noc
c+         if (mod(ncycle,3).eq.0) then
c+c     then we update the rhs
c+            call intermq(dtime)
c+            call remcor
c+            call newsrc
c+            call exitp
c+            call newqfn(dtime,rhs)
c+         endif
cw         write(6,*) 'cycle=',ncycle
c     initialize the coarse grid values of p
         do l=loco(2),maxout
            p(l)= 0.d0
         end do
c     DOWN CYCLE
c     for m=1
      m=1
c     -----------------------
cwc     write out starting resid
cw      call resid(m,nx(m),ny(m),nz(m),cp(loccp(m)),
cw     &        p(loco(m)),rhs(loci(m)),res,maxres)
cw      write(6,*) 'starting resid=',maxres
c     ---------------------------------
      do 11 l=1,nu1
         call linerelax(nx(m),ny(m),nz(m),cp(loccp(m)),p(loco(m)),
     &        rhs(loci(m)))
c     *******try and put mgpfill here ******
         call mgpfill(dtime,p)
c         call mgpfill(dtime,p)
 11   continue
cc      call mgpfill(dtime,p)
      call resid(m,nx(m),ny(m),nz(m),cp(loccp(m)),
     &        p(loco(m)),rhs(loci(m)),res,maxres)
cw      if (ncycle.ne.1) then
cw         ratio = maxres/oldres
cw         write(6,*) ncycle,'maxres',maxres,'conv ratio=',ratio
cw      endif
cw      oldres= maxres
cw      write(6,*) 'm= ',m, 'maxres=',maxres
      if (maxres.lt.tol) goto 101
      call restrict(nx(m),ny(m),nz(m),res,rhs(loci(m+1)))
c      do 100 m=1,ngrid-1
      do 100 m=2,ngrid-1
         do 21 l=1,nu1
            call linerelax(nx(m),ny(m),nz(m),cp(loccp(m)),p(loco(m)),
     &           rhs(loci(m)))
 21      continue
         call resid(m,nx(m),ny(m),nz(m),cp(loccp(m)),
     &        p(loco(m)),rhs(loci(m)),res,maxres)
cw         write(6,*) 'm= ',m, 'maxres=',maxres
         call restrict(nx(m),ny(m),nz(m),res,rhs(loci(m+1)))
 100  continue
c
c     UP CYCLE
      do 200 m= ngrid,2,-1
         if (m.eq.ngrid) then
 25         call sor(nx(m),ny(m),nz(m),cp(loccp(m)),p(loco(m)),
     &           rhs(loci(m)))
         else
            do 41 l=1,nu2
               call linerelax(nx(m),ny(m),nz(m),cp(loccp(m)),
     &              p(loco(m)),rhs(loci(m)))
 41         continue
         endif
c     this call to resid is not necessary - temporary for res check
cw         call resid(m,nx(m),ny(m),nz(m),cp(loccp(m)),
cw     &        p(loco(m)),rhs(loci(m)),res,maxres)
cw         write(6,*) 'm= ',m, 'maxres=',maxres
         call efill(nx(m),ny(m),nz(m),p(loco(m)))
         call prolong(nx(m),ny(m),nz(m),p(loco(m)),p(loco(m-1)) )
 200  continue
c      do l=1,3
c         call  mgpfill(dtime,p)
c      end do
 1000 continue
c*      write(6,*) '100 mg V cycles, not converged'
c     added Nov 17,1993
      m=1
      do 111 l=1,nu1
         call linerelax(nx(m),ny(m),nz(m),cp(loccp(m)),p(loco(m)),
     &        rhs(loci(m)))
c     *******try and put mgpfill here ******
         call  mgpfill(dtime,p)
 111  continue
c
      m=1
      call resid(m,nx(m),ny(m),nz(m),cp(loccp(m)),
     &     p(loco(m)),rhs(loci(m)),res,maxres)
cw      write(6,*) 'max-res=',maxres
c     normliz can be called only for the case of dirichlet bcs.
 101  continue
c 101  call normliz(nx(1),ny(1),nz(1),p)
c      do l=1,3
c           call  mgpfill(dtime,p)
c      end do
c     maxres/edt is the value of (u_x +v_y +ep*w_Z)
      maxres= maxres/edt
      write(6,'(A,I8,A,ES13.6)') 'mg-cycle ', ncycle,',  maxres=',maxres
C       write(6,*) 'mg-cycle ', ncycle,',  maxres=',maxres
c     this call to resid is not necessary - temporary for res check
c      call resid(m,nx(m),ny(m),nz(m),cp(loccp(m)),
c     &     p(loco(m)),rhs(loci(m)),res,maxres)
      goto 201
 20   write(6,*) 'cannot coarsen this grid as specified'

 201  return
      end
