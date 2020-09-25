      subroutine topog
c     ----------------
c     creates the bottom topography

c     D (total depth) is positive and non-dim by DL
c     hg is positive, (ht of grass) and non-dim by DL
c     Dg (depth of grass canopy) is positive and non-dim by DL

      include 'header.f' ! includes variable dz
      integer i,j,k
      double precision DLinv

      DLinv= 1.d0/DL

C     Modified to keep elements about the same size,
C     grass length is a function of NKg

      Grasslen = DLinv*dble(NKg)*dz
C       Grasslen = 1.2d0*DLinv
c      Grasslen = 0.6d0*DLinv

      do j=0,NJ+1
         do i=0,NI+1

            D(i,j) = DLinv*dble(NK)*dz
c            D(i,j) = 2.8d0*DLinv
C             D(i,j) = 2.4d0*DLinv

c     initially, grass is vertical
            ght(i,j)= Grasslen
            Dg(i,j)= D(i,j) - ght(i,j)
            Ddx(i,j) = 0.d0
            Ddy(i,j) = 0.d0
         end do
      end do

      return
      end

      
