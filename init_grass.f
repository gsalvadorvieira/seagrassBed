      subroutine init(pcorr)
c     ------------------------------------------------
c     For the curvilinear grid.
c     sigma levels are evenly spaced.  Hence wz is a function of
c     time but not of z.   Here u,v,w refer to the xi,eta,sigma direcitons.
c     Those metric quantities that do not change with time are evaluated here.
c     The rest are evaluated in "sigma.f" which will be called at every time
c     step.
c     Also it is absolutely essential that the integral of the flux
c     over the entrance section is equal to that over the exit section.
c     -No longer necessary with the free surface.
c     This may require some special attention when the free surface is
c     moving at the exit. At the entrance vf is held fixed.
c
      implicit double precision (a-h,o-z)
      include 'header.f'
      double precision xdu(0:NI+1,0:NJ+1),ydu(0:NI+1,0:NJ+1),
     &     xdv(0:NI+1,0:NJ+1),ydv(0:NI+1,0:NJ+1),
     &     pcorr(maxout),ds
c
      integer i,j,k,it

      write(6,*) 'in init'
c
c     dx,dy are the grid dims (in m) read in in main.f

      do 30 j=0,NJ+1
         do 30 i=0,NI+1
            xdu(i,j)= dx/LEN
            ydu(i,j)= 0.d0
            xdv(i,j)= 0.d0
            ydv(i,j)= dy/LEN
 30   continue
c
c     xc,yc are the grid coords (m)
      yc(0) = -0.5*dy
      do j=1,NJ+1
         yc(j)= yc(j-1) + dy
      end do
      xc(0) = -0.5*dx
      do i=1,NI+1
         xc(i)= xc(i-1) + dx
      end do
c
      do 50 j=0,NJ+1
         do 60 i=0,NI+1
            J2d(i,j)= xdu(i,j)*ydv(i,j) -xdv(i,j)*ydu(i,j)
c            write(300,13) xdu(i,j),ydv(i,j),xdv(i,j),ydu(i,j),J2d(i,j)
c 13         format(5(e8.2,2x))
            ux(i,j)= ydv(i,j)/J2d(i,j)
            vx(i,j)= -ydu(i,j)/J2d(i,j)
            uy(i,j)= -xdv(i,j)/J2d(i,j)
            vy(i,j)= xdu(i,j)/J2d(i,j)
            g11(i,j)= ux(i,j)*ux(i,j) +uy(i,j)*uy(i,j)
            g12(i,j)= ux(i,j)*vx(i,j) +uy(i,j)*vy(i,j)
            g22(i,j)= vx(i,j)*vx(i,j) +vy(i,j)*vy(i,j)
 60      continue
 50   continue
c
      call topog
 101  continue
c
      do 55 j=0,NJ+1
         do 65 i=0,NI+1
            do 75 k=0,NK+1
                  u(i,j,k,0)= 0.d0
                  v(i,j,k,0)= 0.d0
                  w(i,j,k,0)= 0.d0
c     s is the density variation, which is constant
                  s(i,j,k,0)= 0.d0
                  do it=1,ntr
                     T(it,i,j,k,0)= 0.d0
                  end do
                  rho(i,j,k)= 0.d0
 75         continue
 65      continue
 55   continue
      do k=1,maxout
         pcorr(k)= 0.d0
      end do
c
      do k=1,NK
         do j=1,NJ
            do i=1,NI
               uvis(i,j,k)= 0.d0
               vvis(i,j,k)= 0.d0
               wvis(i,j,k)= 0.d0
            end do
         end do
      end do
c
c     specify the initial fields
c     read in the cell-centered pressure and u,v velocities
      do 81 j=0,NJ+1
         do 82 i=0,NI+1
            h(i,j)= 0.d0
 82      continue
 81   continue
c
      do 202 k=1,NK
         do 202 j=1,NJ
            do 202 i=0,NI
               uf(i,j,k)= 0.d0
 202  continue
c
      do 203 k=1,NK
         do 203 j=0,NJ
            do 203 i=1,NI
               vf(i,j,k)= 0.d0
 203  continue
c
      do k=0,NK
         do j=1,NJ
            do i=1,NI
               wf(i,j,k)= 0.d0
            end do
         end do
      end do
c
      do k=1,NK
         do j=1,NJ
            cufu(j,k)= 0.d0
            cufv(j,k)= 0.d0
            cufw(j,k)= 0.d0
            cufw(j,k)= 0.d0
            do it=1,ntr
               cufT(it,j,k)= 0.d0
            end do
         end do
      end do
c
c     Initialize s,T
      call findzall
      call tracerinit
      call evalrho(rho,0)

c     write out z-grid
c     -----------------
      open (unit=60,file='zgrid.out')
      write(60,*) 'vertical grid'
      do k=0,NK+1
         write(60,*) zc(NI/2,NJ/2,k)*DL
      end do

      write(60,*) 'face values'
      do k=-1,NK+1
         write(60,*) zf(NI/2,NJ/2,k)*DL
      end do
      close(60)
c     -----------

      do k=0,NK+1
         do j=0,NJ+1
            do i=0,NI+1
               si(i,j,k)= 0.d0
               sj(i,j,k)= 0.d0
               sk(i,j,k)= 0.d0
            end do
         end do
      end do
      do k=1,NK
         do j=1,NJ
            do i=0,NI
               sifc(i,j,k) = 0.0
            end do
         end do
      end do

      do k=1,NK
         do j=0,NJ
            do i=1,NI
               sjfc(i,j,k)= 0.0
            end do
         end do
      end do

      do k=0,NK
         do j=1,NJ
            do i=1,NI
               skfc(i,j,k)= 0.0
            end do
         end do
      end do

      do k=1,NK
         do i=1,NI
            vfbcn(i,k)= 0.0
            vfbcs(i,k)= 0.0
         end do
      end do
c     ---------------------
c     Initially grass is vertical, theta,angle with the vertical=0
c
      do k=1,NKg
         do j=1,NJ
            do i=1,NI
               costhe(i,j,k)= 0.d0
               sinthe(i,j,k)= 1.d0
               xg(i,j,k,0)= 0.d0
c               theta(i,j,k)= datan(sinthe(i,j,k)/costhe(i,j,k))
               theta(i,j,k)= datan2(sinthe(i,j,k),costhe(i,j,k))
c
c               if ((i.eq.NI/2).and.(j.eq.NJ/2))
c     &              write(6,*) sinthe(i,j,k),costhe(i,j,k),theta(i,j,k)

            end do
         end do
      end do
      ds= Grasslen*DL/dble(NKg)
      do j=1,NJ
         do i=1,NI
            ght(i,j)= Grasslen
            gdef(i,j)= 0.d0
            zg(i,j,1,0)= 0.5*ds
            do k=2,NKg
               zg(i,j,k,0)= zg(i,j,k-1,0)+ ds
            end do
         end do
      end do


      return
      end
