      subroutine vort(n)
c     ---------------------------------------------------
c     computes the vorticity (for now just the vertical component)
      include 'header.f'
      integer i,j,k,n
      double precision wdx,udz,LENinv,DLinv,WL

c     j=0,NJ
      do i=1,NI
         do k=1,NK
            u(i,0,k,n)= 2.0*u(i,1,k,n) -u(i,2,k,n)
            v(i,0,k,n)= 2.0*v(i,1,k,n) -v(i,2,k,n)
            w(i,0,k,n)= 2.0*w(i,1,k,n) -w(i,2,k,n)
c
            u(i,NJ+1,k,n)= 2.0*u(i,NJ,k,n) -u(i,NJ-1,k,n)
            v(i,NJ+1,k,n)= 2.0*v(i,NJ,k,n) -v(i,NJ-1,k,n)
            w(i,NJ+1,k,n)= 2.0*w(i,NJ,k,n) -w(i,NJ-1,k,n)
         end do
      end do
      do j=0,NJ+1
         do i=1,NI
            u(i,j,0,n)= 2.0*u(i,j,1,n) -u(i,j,2,n)
            v(i,j,0,n)= 2.0*v(i,j,1,n) -v(i,j,2,n)
            w(i,j,0,n)= 2.0*w(i,j,1,n) -w(i,j,2,n)
c
            u(i,j,NK+1,n)= 2.0*u(i,j,NK,n) -u(i,j,NK-1,n)
            v(i,j,NK+1,n)= 2.0*v(i,j,NK,n) -v(i,j,NK-1,n)
            w(i,j,NK+1,n)= 2.0*w(i,j,NK,n) -w(i,j,NK-1,n)
         end do
      end do

      if (periodicew) then
c     periodic-ew boundaries
      do k=0,NK+1
         do j=0,NJ+1
            u(0,j,k,n)= u(NI,j,k,n)
            v(0,j,k,n)= v(NI,j,k,n)
            w(0,j,k,n)= w(NI,j,k,n)
c
            u(NI+1,j,k,n)= u(1,j,k,n)
            v(NI+1,j,k,n)= v(1,j,k,n)
            w(NI+1,j,k,n)= w(1,j,k,n)
         end do
      end do
      endif

      DLinv= 1.0/DL
      LENinv=1.0/LEN
      WL = EPS*delta *UL
      do j=1,NJ
         do i=1,NI
            do k=1,NK
               wdx= 0.5*((w(i+1,j,k,0)-w(i-1,j,k,0))*ux(i,j) +
     &              (w(i,j+1,k,0) -w(i,j-1,k,0))*vx(i,j) +
     &              (w(i,j,k+1,0) -w(i,j,k-1,0))*wx(i,j,k) )*WL*LENinv
               udz= 0.5*(u(i,j,k+1,0)-u(i,j,k-1,0))*wz(i,j,k)*UL*DLinv
               vor(i,j,k)= udz - wdx
               
c=               udz = 0.5*( (u(i,j,k+1,0) -u(i,j,k-1,0))*wz(i,j,k) )
c=               wdx = 0.5*((w(i+1,j,k,0)-w(i-1,j,k,0))*ux(i,j) +
c=     &              (w(i,j+1,k,0) -w(i,j-1,k,0))*vx(i,j) +
c=     &              (w(i,j,k+1,0) -w(i,j,k-1,0))*wx(i,j,k) )
c=               z2(i,j,k)= udz*UL*DLinv -wdx*WL*LENinv
c=
c=               wdy= 0.5*((w(i+1,j,k,0)-w(i-1,j,k,0))*uy(i,j) +
c=     &              (w(i,j+1,k,0) -w(i,j-1,k,0))*vy(i,j) +
c=     &              (w(i,j,k+1,0) -w(i,j,k-1,0))*wy(i,j,k) )
c=               vdz= 0.5*( (v(i,j,k+1,0) -v(i,j,k-1,0))*wz(i,j,k) )
c=               z1(i,j,k)= wdy*WL*LENinv -vdz*UL*DLinv
            end do
         end do
      end do

c     boundary points
c     j=0,NJ
      do i=1,NI
         do k=1,NK
            vor(i,0,k)= 2.0*vor(i,1,k) -vor(i,2,k)
            vor(i,NJ+1,k)= 2.0*vor(i,NJ,k) -vor(i,NJ-1,k)
         end do
      end do
      do j=0,NJ+1
         do i=1,NI
            vor(i,j,0)= 2.0*vor(i,j,1) -vor(i,j,2)
            vor(i,j,NK+1)= 2.0*vor(i,j,NK) -vor(i,j,NK-1)
         end do
      end do
c     periodic-ew boundaries
      if (periodicew) then
      do k=0,NK+1
         do j=0,NJ+1
            vor(0,j,k)= vor(NI,j,k)
            vor(NI+1,j,k)= vor(1,j,k)
         end do
      end do
      else 
      do k=0,NK+1
         do j=0,NJ+1
            vor(0,j,k)= 2.0*vor(1,j,k) -vor(2,j,k)
            vor(NI+1,j,k)= 2.0*vor(NI,j,k) -vor(NI-1,j,k)
         end do
      end do
      end if
      return
      end




