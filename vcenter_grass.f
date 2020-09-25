      subroutine vcenter(pf,dtime,n)
c----------------------------------------------------
c     compute the final vel (n+1 th step) at the cell centers
      implicit logical (a-z)
      include 'header.f'
      integer i,j,k,n
      double precision pxi,peta,psig,dte,dtime,dJtinv,
     &     px,py,pz,fac,pf(0:NI+1,0:NJ+1,0:NK+1),wfk
      double precision ubold(NJ,NK),vbold(NJ,NK),wbold(NJ,NK),du
c     operation:  u(n)= u(0) + dtime* f(u(m))
c
c     Save old values of u,v,w to calculate cufu,etc.
c     ------------------------------------------------
      i=NI
      do k=1,NK
         do j=1,NJ
            ubold(j,k) =u(i,j,k,0)
            vbold(j,k)= v(i,j,k,0)
            wbold(j,k)= w(i,j,k,0)
         end do
      end do
c     -------------------------

      dte= dtime/EPS
c
      do 10 j=1,NJ
         do 20 i=1,NI
            do 30 k=1,NK
               pxi= 0.5d0*(pf(i+1,j,k)-pf(i-1,j,k))
               peta= 0.5d0*(pf(i,j+1,k)-pf(i,j-1,k)) 
               psig= 0.5d0*(pf(i,j,k+1)-pf(i,j,k-1)) 
               px= ux(i,j)*pxi +vx(i,j)*peta +wx(i,j,k)*psig
               py= uy(i,j)*pxi +vy(i,j)*peta +wy(i,j,k)*psig
               pz= wz(i,j,k)*psig
               u(i,j,k,n)=  cx(i,j,k) -dte*(qpr*px )
               v(i,j,k,n)=  cy(i,j,k) -dte*(qpr*py )
cc
c=               w(i,j,k,n)=  cz(i,j,k) -dtbeta*( pz )
 30         continue
 20      continue
 10   continue
c
c     vface is already called - compute w from wf
      do 11 j=1,NJ
         do 21 i=1,NI
            do 31 k=1,NK
               wfk= 0.5d0*(wf(i,j,k) +wf(i,j,k-1))
               w(i,j,k,n)= (wfk/Jac(i,j,k) -u(i,j,k,n)*wx(i,j,k)
     &              -v(i,j,k,n)*wy(i,j,k) )/(EPS*wz(i,j,k))
 31         continue
 21      continue
 11   continue
c
c     Now calculate cufu,cufv,cufw
c     ----------------------------
c     cufu,cufv,and cufw are the Uf velocities for u,v,w

      i=NI
      do k=1,NK
         do j=1,NJ
            dJtinv=Jac(i,j,k)/dtime
            du = (u(i,j,k,n) -u(i-1,j,k,n))
            if (dabs(du).gt.1.d-16) then
               cufu(j,k) = -(u(i,j,k,n)-ubold(j,k))*dJtinv/du
            else 
               cufu(j,k)= 0.d0
            end if
            du = (v(i,j,k,n) -v(i-1,j,k,n))
            if (dabs(du).gt.1.d-16) then
               cufv(j,k) = -(v(i,j,k,n)-vbold(j,k))*dJtinv/du
            else
               cufv(j,k) = 0.d0
            end if
            du = (w(i,j,k,n) -w(i-1,j,k,n))
            if (dabs(du).gt.1.d-16) then
               cufw(j,k) = -(w(i,j,k,n)-wbold(j,k))*dJtinv /du
            else 
               cufw(j,k) = 0.d0
            end if
         end do
      end do

      return
      end
      


