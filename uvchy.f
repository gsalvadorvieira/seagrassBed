      subroutine uvchy(dtime)
c----------------------------------------------------
c     compute the tilde velocities at cell centers and store them in n=1
      implicit logical (a-z)
      include 'header.f'
      integer i,j,k
      double precision dte,dtime,hxi,heta,hx,hy,kaph1,
     &     Ub,Vb,Wb,d2,c1,c2,epsinv,temp1,temp2,wfk
c     operation:  u(n)= u(0) + dtime* f(u(m))
c
      double precision pxi,peta,psig,px,py
      double precision drpx(NI,NJ,NK),drpy(NI,NJ,NK),
     &     grpifc(0:NI,NJ,NK),grpjfc(NI,0:NJ,NK)
      common/rpgrads/drpx,drpy,grpifc,grpjfc
c
      dte= dtime/EPS
      epsinv= 1.d0/EPS
      kaph1= 1.d0 - kappah
c
      do 10 j=1,NJ
         do 20 i=1,NI
            hxi= 0.5d0*( h(i+1,j)-h(i-1,j) )
            heta= 0.5d0*( h(i,j+1)-h(i,j-1) )
            hx= ux(i,j)*hxi +vx(i,j)*heta
            hy= uy(i,j)*hxi +vy(i,j)*heta
            do 30 k=1,NK

               pxi= 0.5d0*(p(i+1,j,k)-p(i-1,j,k))
               peta= 0.5d0*(p(i,j+1,k)-p(i,j-1,k))
               psig= 0.5d0*(p(i,j,k+1)-p(i,j,k-1))
               px= ux(i,j)*pxi +vx(i,j)*peta +wx(i,j,k)*psig
               py= uy(i,j)*pxi +vy(i,j)*peta +wy(i,j,k)*psig

c     cx and cy contain the convective terms.
               cx(i,j,k)=  cx(i,j,k) -dte*( gpr*(kappah*hx
     &              + kaph1*gradhn(i,j,1)) + si(i,j,k) 
     &              + qpr *px )
               cy(i,j,k)=  cy(i,j,k) -dte*( gpr*(kappah*hy
     &              + kaph1*gradhn(i,j,2)) + sj(i,j,k) 
     &              + qpr *py )
c               cz(i,j,k)=  cz(i,j,k)
 30         continue
            gradhn(i,j,1)= hx
            gradhn(i,j,2)= hy
 20      continue
 10   continue
c
c     czf is already computed
      do 11 j=1,NJ
         do 21 i=1,NI
            do 31 k=1,NK
               wfk= 0.5d0*(czf(i,j,k) +czf(i,j,k-1))
               cz(i,j,k)= (wfk/Jac(i,j,k) -cx(i,j,k)*wx(i,j,k)
     &              -cy(i,j,k)*wy(i,j,k) )/(EPS*wz(i,j,k))
 31         continue
 21      continue
 11   continue
      return
c     ======
c     Boundaries (these are not probably needed)
      do j=1,NJ
         do i=1,NI
            cx(i,j,0)= cx(i,j,1)
            cx(i,j,NK+1)= cx(i,j,NK)
            cy(i,j,0)= cy(i,j,1)
            cy(i,j,NK+1)= cy(i,j,NK)
         end do
      end do
      do k=0,NK+1
         do i=1,NI
            cx(i,0,k)= cx(i,1,k)
            cx(i,NJ+1,k)= cx(i,NJ,k)
            cy(i,0,k)= cy(i,1,k)
            cy(i,NJ+1,k)= cy(i,NJ,k)
         end do
      end do
      do k=0,NK+1
         do j=0,NJ+1
            cx(0,j,k)= cx(NI,j,k)
            cx(NI+1,j,k)= cx(1,j,k)
            cy(0,j,k)= cy(NI,j,k)
            cy(NI+1,j,k)= cy(1,j,k)
         end do
      end do
c
cc     Boundaries
c      i=NI+1
c      do 35 j=jex1,jex2
c         hxi= h(i,j)-h(i-1,j)
c         if (j.eq.jex1) then
c            heta= h(i,j+1)-h(i,j)
c         else if (j.eq.jex2) then
c            heta= h(i,j)-h(i,j-1)
c         else
c            heta= 0.5*(h(i,j+1)-h(i,j-1))
c         endif
c         hx= ux(i,j)*hxi +vx(i,j)*heta
c         hy= uy(i,j)*hxi +vy(i,j)*heta
c         do 36  k=1,NK
c            cx(i,j,k)=  cx(i,j,k) -dte*( gpr*(kappah*hx
c     &           + kaph1*gradhn(i,j,1)) + si(i,j,k) )
c            cy(i,j,k)=  cy(i,j,k) -dte*( gpr*(kappah*hy
c     &           + kaph1*gradhn(i,j,2)) + sj(i,j,k) )
cc            cz(i,j,k)=  cz(i,j,k)
c 36      continue
c         gradhn(i,j,1)= hx
c         gradhn(i,j,2)= hy
c 35   continue
c
cc     Need not fill the ghost points. The ghost points must not be used
cc     in newcor.
cc     *******************************
cc
cc      do k=1,NK
cc         do j=1,NJ
cc            if (ufbcw(j,k).lt.0.d0) then
cc               cx(0,j,k)= uwest(j,k)
cc               cy(0,j,k)= cy(1,j,k)
cc               cz(0,j,k)= cz(1,j,k)
cc            else
ccc     incoming at west boundary
cc               cx(0,j,k)= uwest(j,k)
cc               cy(0,j,k)= vwest(j,k)
cc               cz(0,j,k)= wwest(j,k)
cc            end if
cc            if (ufbce(j,k).gt.0.d0) then
cc               cx(NI+1,j,k)= ueast(j,k)
cc               cy(NI+1,j,k)= cy(NI,j,k)
cc               cz(NI+1,j,k)= cz(NI,j,k)
cc            else
ccc     incoming at west boundary
cc               cx(NI+1,j,k)= ueast(j,k)
cc               cy(NI+1,j,k)= veast(j,k)
cc               cz(NI+1,j,k)= weast(j,k)
cc            end if
cc         end do
cc      end do
cc      do k=1,NK
cc         do i=1,NI
cc            if (vfbcs(i,k).lt.0.d0) then
cc               cy(i,0,k)= vsouth(i,k)
cc               cx(i,0,k)= cx(i,1,k)
cc               cz(i,0,k)= cz(i,1,k)
cc            else
ccc     incoming at south boundary
cc               cy(i,0,k)= vsouth(i,k)
cc               cx(i,0,k)= usouth(i,k)
cc               cz(i,0,k)= wsouth(i,k)
cc            end if
cc            if (vfbcn(i,k).gt.0.d0) then
cc               cy(i,NJ+1,k)= vnorth(i,k)
cc               cx(i,NJ+1,k)= cx(i,NJ,k)
cc               cz(i,NJ+1,k)= cz(i,NJ,k)
cc            else
ccc     incoming at north boundary
cc               cy(i,NJ+1,k)= vnorth(i,k)
cc               cx(i,NJ+1,k)= unorth(i,k)
cc               cz(i,NJ+1,k)= wnorth(i,k)
cc            end if
cc         end do
cc      end do
ccc     vertical edges  (use east and west boundary arrays)
cc      do k=1,NK
cc         cx(0,0,k)= uwest(0,k)
cc         cx(NI+1,0,k)= ueast(0,k)
cc         cx(0,NJ+1,k)= uwest(NJ+1,k)
cc         cx(NI+1,NJ+1,k)= ueast(NJ+1,k)
cc         cy(0,0,k)= vwest(0,k)
cc         cy(NI+1,0,k)= veast(0,k)
cc         cy(0,NJ+1,k)= vwest(NJ+1,k)
cc         cy(NI+1,NJ+1,k)= veast(NJ+1,k)
cc         cz(0,0,k)= wwest(0,k)
cc         cz(NI+1,0,k)= weast(0,k)
cc         cz(0,NJ+1,k)= wwest(NJ+1,k)
cc         cz(NI+1,NJ+1,k)= weast(NJ+1,k)
cc      end do
cc      do j=0,NJ+1
cc         do i=0,NI+1
cc            cx(i,j,NK+1)= cx(i,j,NK)
cc            cy(i,j,NK+1)= cy(i,j,NK)
cc            cz(i,j,NK+1)= cz(i,j,NK)
cc            cx(i,j,0)= cx(i,j,1)
cc            cy(i,j,0)= cy(i,j,1)
cc            cz(i,j,0)= -cz(i,j,1)
cc         end do
cc      end do
c




      return
      end
      

