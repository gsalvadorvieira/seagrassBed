      subroutine energy(step)
c---------------------------------------------------
      implicit logical (a-z)
      include 'header.f'
c
      integer i,j,k,step
      double precision const,enkin,enin,enout,enbal
c
      const= EPS*EPS*delta*delta
      enkin= 0.d0
      do k=1,NK
         do j=1,NJ
            do i=1,NI
               enkin= Jac(i,j,k)*(u(i,j,k,0)*u(i,j,k,0) +v(i,j,k,0)*
     &              v(i,j,k,0) + const*w(i,j,k,0)*w(i,j,k,0))
     &              +enkin
            end do
         end do
      end do
c
cc     flow in
c      enin= 0.d0
c      do k=1,NK
c         do i=ient1,ient2
c            enin= enin + vfbcs(i,k)*dtf*(uent(i,k)*uent(i,k) +
c     &           vent(i,k)*vent(i,k) )
c
c         end do
c      end do
cc
cc     flow out
c      enout = 0.d0
c      do k=1,NK
c         do j=jex1,jex2
c            enout= enout +ufbce(j,k)*dtf*(u(NI,j,k,0)*u(NI,j,k,0) +
c     &           v(NI,j,k,0)*v(NI,j,k,0) +const*w(NI,j,k,0)*
c     &           w(NI,j,k,0) )
c         end do
c      end do
cc
c      enbal= enin + enout + enkin
c
c      write(901,11) step,enin,enout,enkin,enbal
      write(91,*) step,enkin
C 11   format(I5,5x,4(e10.4,4x))
      return
      end
