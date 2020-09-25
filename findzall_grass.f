      subroutine findzall
c     -------------------------------------------------------
c     finds the value of z (non-dim by DL) given the value of sigma.

c     D and Dg are  positive and non-dim by DL.

      include 'header.f'
      integer i,j,k
      double precision sigc,sigf,dNKg,dnkmnginv,Dgph,dnkginv,Dmdg,
     &     DLinv
c
      dNKg= dble(NKg)
      dnkginv= 1.d0/dble(NKg)
      dnkmnginv = 1.d0/dble(NK -NKg)
      DLinv = 1.d0/DL

      do j=0,NJ+1
         do i=0,NI+1
            Dg(i,j)= D(i,j) - ght(i,j)*DLinv
            Dgph = Dg(i,j) + HDL*h(i,j)
            Dmdg = D(i,j) -Dg(i,j)
            do k=0,NKg
               sigc = dble(k) - 0.5d0
               sigf = dble(k)
               zc(i,j,k)= (sigc)*dnkginv*Dmdg -D(i,j)
               zf(i,j,k)= (sigf)*dnkginv*Dmdg -D(i,j)
            end do
c
            do k=NKg+1,NK+1
               sigc= dble(k)-0.5
               zc(i,j,k)= (sigc -dNKg)*dnkmnginv*Dgph -Dg(i,j)
            end do
            do k=NKg+1,NK
               sigf= dble(k)
               zf(i,j,k)= (sigf -dNKg)*dnkmnginv*Dgph -Dg(i,j)
            end do
         end do
      end do
c
      return
      end
