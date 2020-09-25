      subroutine chfine(dtime,chf,fn)
c     -----------------------------------------------
c     call chfine(dtime,n,cp(1),rhs(1))
c     computes the coefficients in the Laplace operator
c     This routine is called right after intpol when cxf,cyf are U*,V*.
c     uses the 9 point stencil (in 2d).
c     
      implicit logical (a-z)
      include 'header.f'
c      double precision v1(0:NI+1,0:NJ+1,0:NK+1),
c     &     v2(0:NI+1,0:NJ+1,0:NK+1),v3(0:NI+1,0:NJ+1,0:NK+1),
c     &     v4(0:NI+1,0:NJ+1,0:NK+1),v5(0:NI+1,0:NJ+1,0:NK+1)
c     
c     cpim1 is the coeff of h(i-1,j), cpip1 is the coeff of h(i+1,j),etc
c     fn is the contribution to the source terms on RHS.
c     Whenever possible we substitue the the velocities themselves
c     (at the boundaries). We cannot substitute the u,w vels at the entrance.
c     Our eqn is
c     cpim1*p(i-1) +cpip1*p(i+1) +cpjm1*p(j-1) ...- cp0*p(i,j,k)= fn
c     
c     chf(1,..)=cpkm1(..), chf(2,..)=cpjp1(..), chf(3,..)=cpip1(..),
c     chf(4,..)=cpjm1(..), chf(5,..)=cpim1(..), chf(6,..)=cp0(..),
c     chf(7,..)=cpkp1(..)
      integer i,j,k,l,im1
      double precision chf(9,NI,NJ),fn(NI,NJ),edt,dtime,
     &     edtg,gpkinv,dtinv,sumsif(NI,NJ),sumsjf(NI,NJ),
     &     sumuf(NI,NJ),sumvf(NI,NJ),
     &     sumhxn(NI,NJ),sumhyn(NI,NJ),sumgi(NI,NJ,2),sumgj(NI,NJ,2),
     &     sumcxf(NI,NJ),sumcyf(NI,NJ),const
c     
      edt= EPS/dtime
      edtg= edt/(gpr*kappah)
      gpkinv= 1.d0/(gpr*kappah)
      dtinv= HDL/(dtime*kappah)
c      dnk= dfloat(NK)
      const= kaphinv -1.d0
c      write(6,*) 'const=',const
c      write(6,*) 'in chfine,edt,eg,gpkinv,dtinv,dnk',
c     &     edt,edtg,gpkinv,dtinv,dnk 
c     
c     
      call hbc(chf,fn,dtime)
c     
c     For i=2...NI-1, j=2...NJ-1, k=2,NK-1  it is straightforward
      do 10 i=1,NI-1  
c      do 10 i=1,NI  ! for periodicew bc
         do 20 j=1,NJ-1
            sumsif(i,j)= 0.d0
            sumsjf(i,j)= 0.d0
            sumcxf(i,j)= 0.d0
            sumcyf(i,j)= 0.d0
            sumuf(i,j)= 0.d0
            sumvf(i,j)= 0.d0
            sumhxn(i,j)= 0.d0
            sumhyn(i,j)= 0.d0
            do l=1,2
               sumgi(i,j,l)= 0.d0
               sumgj(i,j,l)= 0.d0
            end do
            do 30 k=1,NK
               sumsif(i,j)= sumsif(i,j) +sifc(i,j,k)
               sumsjf(i,j)= sumsjf(i,j) +sjfc(i,j,k)
               sumcxf(i,j)= sumcxf(i,j) +cxf(i,j,k)
               sumcyf(i,j)= sumcyf(i,j) +cyf(i,j,k)
               sumuf(i,j)= sumuf(i,j) +uf(i,j,k)
               sumvf(i,j)= sumvf(i,j) +vf(i,j,k)
               sumhxn(i,j)= sumhxn(i,j) +hxn(i,j,k)
               sumhyn(i,j)= sumhyn(i,j) +hyn(i,j,k)
               do l=1,2
                  sumgi(i,j,l)= sumgi(i,j,l) + gi(i,j,k,l)
                  sumgj(i,j,l)= sumgj(i,j,l) + gj(i,j,k,l)
               end do
 30         continue
 20      continue
 10   continue
c     
c
cc      do 40 i=2,NI for periodicew
      do 40 i=2,NI-1
         do 50 j=2,NJ-1
            fn(i,j)= edtg*( kaphinv*wfbcb(i,j)
     &           -J2d(i,j)*oldh(i,j)*dtinv
     &           + sumcxf(i,j) - sumcxf(i-1,j)
     &           + sumcyf(i,j) - sumcyf(i,j-1)
     &           + const*(sumuf(i,j) - sumuf(i-1,j)
     &           + sumvf(i,j) - sumvf(i,j-1) ) )
     &           -gpkinv*( sumsif(i,j) - sumsif(i-1,j) + 
     &           sumsjf(i,j) - sumsjf(i,j-1) )
     &           - const*( sumhxn(i,j) - sumhxn(i-1,j) + 
     &           sumhyn(i,j) - sumhyn(i,j-1) )
c     
            chf(1,i,j)= -edtg*J2d(i,j)*dtinv - (sumgi(i,j,1) +
     &           sumgi(i-1,j,1)+sumgj(i,j,2) +sumgj(i,j-1,2))
            chf(2,i,j)= (sumgi(i,j,1) +0.25*(sumgj(i,j,1) 
     &           -sumgj(i,j-1,1)))
            chf(3,i,j)= (sumgi(i-1,j,1) +0.25*(-sumgj(i,j,1)
     &           +sumgj(i,j-1,1)) )
            chf(4,i,j)= (sumgj(i,j,2) +0.25*(-sumgi(i-1,j,2)
     &           +sumgi(i,j,2)) )
            chf(5,i,j)= (sumgj(i,j-1,2) +0.25*(-sumgi(i,j,2)+
     &           sumgi(i-1,j,2)) )
            chf(6,i,j)= 0.25*(sumgi(i,j,2)+sumgj(i,j,1))
            chf(7,i,j)= -0.25*(sumgi(i-1,j,2)+sumgj(i,j,1))
            chf(8,i,j)= -0.25*(sumgi(i,j,2) +sumgj(i,j-1,1))
            chf(9,i,j)= 0.25*(sumgi(i-1,j,2) +sumgj(i,j-1,1))
 50      continue
 40   continue
c     
c      do k=1,NK
c         do j=2,NJ
c            do i=2,NI
c               v1(i,j,k)= sumcxf(i,j) - sumcxf(im1,j)
c               v2(i,j,k)= sumcyf(i,j) - sumcyf(i,j-1)
cc               v3(i,j,k)= grpjfc(i,j,k)
c               v4(i,j,k)= sjfc(i,j,k)
c               v5(i,j,k)= sumsjf(i,j)
c            end do
c         end do
c      end do
c      call outarray(v1,v2,v3,v4,v5)
c      write(6,*) 'stopping in chfine'
c      stop
c      call writout(chf,fn)
      return
      end
      
