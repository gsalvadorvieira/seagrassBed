      subroutine calcpgrad(pgrad,hgrad,bgrad)
c----------------------------------------------------
c     compute the tilde velocities at cell centers and store them in n=1
      implicit logical (a-z)
      include 'header.f'
      integer i,j,k
      double precision dte,dtime,hxi,heta,hx,hy,kaph1,epsinv
c     operation:  u(n)= u(0) + dtime* f(u(m))
c
      double precision pxi,peta,psig,px,py
      double precision pgrad(0:NI+1,0:NJ+1,0:NK+1),
     &     hgrad(0:NI+1,0:NJ+1,0:NK+1),bgrad(0:NI+1,0:NJ+1,0:NK+1)
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
               hgrad(i,j,k)=   gpr*(kappah*hx
     &              + kaph1*gradhn(i,j,1))
               bgrad(i,j,k)= si(i,j,k)
               pgrad(i,j,k)= qpr *px
 30         continue
 20      continue
 10   continue

      return
      end
