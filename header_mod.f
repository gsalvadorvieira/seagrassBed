c     header for subroutines
c----------------------------------------------------------
      include 'dims.f'
      logical rect, periodicew
c
C       double precision LEN,DL,UL,EPS,vmu,Cd,GL,gpr,zNL,dz,dx,dy
      double precision LEN,DL,UL,EPS,vmu,Cd,GL,gpr,dz,dx,dy
      double precision delta,qpr,PI,Grasslen,CsA,CdN,CdT,dHinit,
     &     dHfinal
c     dH is the free surface slope (or charac free-surface ht diff per LEN)
c     dHinit is initially, dHfinal is finally
c     mu is the viscosity
c     zNL is the number of grass stems per unit area (i.e.per m^2)
c     CsA = charac. cross sectional area of grass stem (1cm2)

      parameter (rect = .true., periodicew = .false.)
      parameter (DL = 1.d0, LEN = 1.d0, UL = 0.1d0, EPS=1.d0,
c     &     vmu= 1.d-4, Cd= 3.d-1, GL= 10.d0, gpr=0.981d0,
     &     Cd=0.6d0, GL=10.d0, gpr=0.981d0,
C      &     zNL= 100.d0, dz=0.05d0, PI=3.141592653589793, CsA=1.d-4,
     &     PI=3.141592653589793, CsA=1.d-4,
C    AUTOMATICALLY MODIFYING NEXT LINES (23,24)  
     &     dHinit=3.d-3, dHfinal=3.d-3,
     &     dz=0.05d0)

      double precision zNL(0:NI+1)

      double precision xc(0:NI+1),yc(0:NJ+1),
     &     zc(0:NI+1,0:NJ+1,0:NK+1),zf(0:NI+1,0:NJ+1,-1:NK+1)
      double precision uss(NK),ufss(NK),uprof(NK),ufprof(NK),
     &     tinfl(0:NK+1),
     &     udifout(NJ,NK),vdifout(NJ,NK),wdifout(NJ,NK)
      double precision cufu(NJ,NK),cufv(NJ,NK),cufw(NJ,NK),
     &     cufs(NJ,NK),cufT(ntr,NJ,NK),ufbold(NJ,NK)
      double precision
     &     uf(0:NI,NJ,NK),u(0:NI+1,0:NJ+1,0:NK+1,0:1),
     &     v(0:NI+1,0:NJ+1,0:NK+1,0:1),
     &     vf(NI,0:NJ,NK),
     &     w(0:NI+1,0:NJ+1,0:NK+1,0:1),wf(NI,NJ,0:NK),
     &     p(0:NI+1,0:NJ+1,0:NK+1),
     &     rho(0:NI+1,0:NJ+1,0:NK+1),vor(0:NI+1,0:NJ+1,0:NK+1),
     &     gradhn(0:NI+1,0:NJ+1,2),
     &     hxn(0:NI,NJ,NK),hyn(NI,0:NJ,NK),
     &     si(0:NI+1,0:NJ+1,0:NK+1),sj(0:NI+1,0:NJ+1,0:NK+1),
     &     sk(0:NI+1,0:NJ+1,0:NK+1),
c     &     oldh(0:NI+1,0:NJ+1),
     &     sifc(0:NI,NJ,NK),sjfc(NI,0:NJ,NK),skfc(0:NI+1,0:NJ+1,0:NK)
      double precision ght(0:NI+1,0:NJ+1),gdef(0:NI+1,0:NJ+1),
     &     xg(0:NI+1,0:NJ+1,NKg,0:1),zg(0:NI+1,0:NJ+1,NKg,0:1),
     &     costhe(0:NI+1,0:NJ+1,NKg),sinthe(0:NI+1,0:NJ+1,NKg),
     &     theta(0:NI+1,0:NJ+1,NKg)
      double precision
     &     s(0:NI+1,0:NJ+1,0:NK+1,0:1),
     &     T(ntr,0:NI+1,0:NJ+1,0:NK+1,0:1),
     &     cx(0:NI+1,0:NJ+1,0:NK+1),cy(0:NI+1,0:NJ+1,0:NK+1),
     &     cz(0:NI+1,0:NJ+1,0:NK+1),cxf(0:NI,NJ,NK),cyf(NI,0:NJ,NK),
     &     czf(NI,NJ,0:NK)
c     &     h(0:NI+1,0:NJ+1)
      double precision
     &     ufbce(NJ,NK),ufbcw(NJ,NK),vfbcn(NI,NK),vfbcs(NI,NK),
     &     wfbcb(NI,NJ),ueast(0:NJ+1,NK),uwest(0:NJ+1,NK),
     &     vnorth(NI,NK),vsouth(NI,NK),ssouth(NI,NK),sbackgrnd

      double precision D(0:NI+1,0:NJ+1),Ddx(0:NI+1,0:NJ+1),
     &     Ddy(0:NI+1,0:NJ+1),Dg(0:NI+1,0:NJ+1),
     &     J2d(0:NI+1,0:NJ+1),Jac(0:NI+1,0:NJ+1,0:NK+1),
     &     g11(0:NI+1,0:NJ+1),g22(0:NI+1,0:NJ+1),
     &     g12(0:NI+1,0:NJ+1),gi(0:NI,NJ,NK,2),gj(NI,0:NJ,NK,2),
     &     gi3(0:NI,NJ,NK),gj3(NI,0:NJ,NK),gqi(0:NI,NJ,NK,2),
     &     gqj(NI,0:NJ,NK,2),gqi3(0:NI,NJ,NK),
     &     gqj3(NI,0:NJ,NK),gqk(0:NI+1,0:NJ+1,0:NK,3),
     &     Jifc(0:NI,NJ,NK),Jjfc(NI,0:NJ,NK)
      double precision dtf,P1,delinv,HL,HDL,TL,kappah,kaphinv
      double precision ux(0:NI+1,0:NJ+1),uy(0:NI+1,0:NJ+1),
     &     vx(0:NI+1,0:NJ+1),vy(0:NI+1,0:NJ+1),
     &     wx(0:NI+1,0:NJ+1,0:NK+1),
     &     wy(0:NI+1,0:NJ+1,0:NK+1),wz(0:NI+1,0:NJ+1,0:NK+1),
     &     wzk(0:NI+1,0:NJ+1,0:NK)

      double precision uvis(NI,NJ,NK),vvis(NI,NJ,NK),wvis(NI,NJ,NK)
     &     ,Kz(NI,NJ,0:NK)
      common/param/delta,qpr,Grasslen,CdN,CdT,vmu
      common/vis/uvis,vvis,wvis,Kz
      common/grid/xc,yc,zc,zf,D,Ddx,Ddy,Dg,ux,vx,uy,vy,wx,wy,wz,wzk
      common/metrics/J2d,Jac,g11,g22,g12,gi,gj,
     &     gi3,gj3,gqi,gqj,gqi3,gqj3,gqk,Jifc,Jjfc
c      common/std/xs,ys,xt,yt,bs,cs,ds,bt,ct,dt
      common/other/ dtf,P1,delinv,HL,HDL,TL,kappah,kaphinv,
     &     dx,dy
      common/orlanski/uss,ufss,tinfl,udifout,vdifout,wdifout,
     &     uprof,ufprof
      common/variable/u,uf,v,vf,w,wf,p,rho,vor,gradhn,
     &     hxn,hyn,si,sj,sk,sifc,sjfc,skfc,s,T,cx,cy,cz,cxf,cyf,
     &     czf
c    &     oldh, h
      common/grass/ght,gdef,xg,zg,costhe,sinthe,theta
      common/bcs/ufbce,ufbcw,
     &     vfbcn,vfbcs,wfbcb,ueast,uwest,vnorth,vsouth,ssouth,sbackgrnd,
     &     cufu,cufv,cufw,cufs,cufT,ufbold

