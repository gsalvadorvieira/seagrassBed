      subroutine sigma
c     ------------------------
c     for openbc
c     This subroutine updates all those quantities that are a function
c     of time and  (time and depth). It is called  every time step.
c
c     h = ht of free surface
c     D = depth of water
c     Dg= depth of grass-water interface (changes with time)
c     
      implicit logical (a-z)
      include 'header.f'
      integer i,j,k

      double precision hpdg,hpdginv,hu,hv,hx,hy,z,temp,
     &     be2,wxk,wyk,d2,dnkg,dnkmng,sig,Dgu,Dgv,Dgx,Dgy,dmdg,
     &     hxpdgx,hypdgy,sigc,
     &     dmdginv,dmdginv2,g13(0:NI+1,0:NJ+1,0:NK+1),
     &     g23(0:NI+1,0:NJ+1,0:NK+1)
c     Ddx and Ddy are known from init
c
      dnkmng= dble(NK-NKg)
      dnkg= dble(NKg)
      be2 = 1.d0

      do j=0,NJ+1
         do i=0,NI+1
c     Done in grassben      Dg(i,j)= D(i,j) - ght(i,j)*DLinv
c     ght is non-dim by DL and positive
c     D(i,j) is positive and non-dim by DL . 
c     Dg(i,j) (depth of grass) is positive and non-dim by DL

            hpdg= h(i,j)*HDL +Dg(i,j)
            hpdginv= 1.d0/hpdg
            dmdg = D(i,j) - Dg(i,j)
            dmdginv = 1.d0/dmdg
            dmdginv2 = 1.d0/(dmdg*dmdg)
            if (i.eq.0) then
               hu= HDL*(h(i+1,j)-h(i,j))
               Dgu = Dg(i+1,j)-Dg(i,j)
            else if  (i.eq.NI+1) then
               hu= HDL*(h(i,j)-h(i-1,j))
               Dgu= Dg(i,j)-Dg(i-1,j)
            else
               hu= 0.5d0*HDL*(h(i+1,j)-h(i-1,j))
               Dgu= 0.5d0*(Dg(i+1,j)-Dg(i-1,j))
            endif
            if (j.eq.0) then
               hv= HDL*(h(i,j+1)-h(i,j))
               Dgv= Dg(i,j+1)-Dg(i,j)
            else if (j.eq.NJ+1) then
               hv= HDL*(h(i,j)-h(i,j-1))
               Dgv= Dg(i,j)-Dg(i,j-1)
            else
               hv= 0.5d0*HDL*(h(i,j+1)-h(i,j-1))
               Dgv= 0.5d0*(Dg(i,j+1)-Dg(i,j-1))
            endif
            hx= hu*ux(i,j) + hv*vx(i,j)
            hy= hu*uy(i,j) + hv*vy(i,j)
            Dgx= Dgu*ux(i,j) + Dgv*vx(i,j)
            Dgy= Dgu*uy(i,j) + Dgv*vy(i,j)
            hxpdgx= hx + Ddx(i,j)
            hypdgy= hy + Ddy(i,j)
c     wz is not a function of depth when the sigma lines are equally spaced.
c     Hence wz is wz(i,j,time). For a stretched grid wz would be w(i,j,k,time).
c     Then Jac which is now Jac(i,j,time) would  become  Jac(i,j,k,time)
c
c            write(100,*) J2d(i,j)
c            write(200,*) Jac(i,j)
c
            do  k=0,NKg
               wz(i,j,k)= dmdginv*dnkg
               Jac(i,j,k)= J2d(i,j)/wz(i,j,k)
               sigc = dble(k) -0.5
               wx(i,j,k)= (dnkg**Ddx(i,j)-sigc*(Ddx(i,j)-Dgx))*dmdginv
               wy(i,j,k)= (dnkg**Ddy(i,j)-sigc*(Ddy(i,j)-Dgy))*dmdginv

               g13(i,j,k)= ux(i,j)*wx(i,j,k) +uy(i,j)*wy(i,j,k)
               g23(i,j,k)= vx(i,j)*wx(i,j,k) +vy(i,j)*wy(i,j,k)
c               g33(i,j,k)= wx(i,j,k)*wx(i,j,k) +wy(i,j,k)*wy(i,j,k) +
c     &              be2*wz(i,j)*wz(i,j)
            end do
            do k= NKg+1,NK+1
               wz(i,j,k)= hpdginv*dnkmng
               Jac(i,j,k)= J2d(i,j)/wz(i,j,k)
               sigc = dble(k) -0.5
               wx(i,j,k)= (dnkmng*Dgx - (sigc-dnkg)*hxpdgx )*hpdginv
               wy(i,j,k)= (dnkmng*Dgy - (sigc-dnkg)*hypdgy )*hpdginv
c
               g13(i,j,k)= ux(i,j)*wx(i,j,k) +uy(i,j)*wy(i,j,k)
               g23(i,j,k)= vx(i,j)*wx(i,j,k) +vy(i,j)*wy(i,j,k)
            end do
            do k=0,NKg
               wzk(i,j,k)= dmdginv*dnkg
            end do
            do k=NKg+1,NK
               wzk(i,j,k)= hpdginv*dnkmng
            end do

c     
c     We evaluate gk at i=0,NI+1,j=0,NJ+1 only because these are used
c     at k=0,NK to fill the horizontal edges in mgpfill or hnbc.
            do k=0,NKg
               sigc= dble(k)
               wxk= (dnkg**Ddx(i,j)-sigc*(Ddx(i,j)-Dgx))*dmdginv
               wyk= (dnkg**Ddy(i,j)-sigc*(Ddy(i,j)-Dgy))*dmdginv
               gqk(i,j,k,1)= qpr*Jac(i,j,k)*(ux(i,j)*wxk +uy(i,j)*wyk)
               gqk(i,j,k,2)= qpr*Jac(i,j,k)*(vx(i,j)*wxk +vy(i,j)*wyk)
               gqk(i,j,k,3)= Jac(i,j,k)*(qpr*(wxk*wxk +wyk*wyk) +
     &              be2*wzk(i,j,k)*wzk(i,j,k))
            end do
            do k=NKg+1,NK
               sigc= dble(k)
               wxk= (dnkmng*Dgx - (sigc-dnkg)*hxpdgx )*hpdginv
               wyk= (dnkmng*Dgy - (sigc-dnkg)*hypdgy )*hpdginv
               gqk(i,j,k,1)= qpr*Jac(i,j,k)*(ux(i,j)*wxk +uy(i,j)*wyk)
               gqk(i,j,k,2)= qpr*Jac(i,j,k)*(vx(i,j)*wxk +vy(i,j)*wyk)
               gqk(i,j,k,3)= Jac(i,j,k)*(qpr*(wxk*wxk +wyk*wyk) +
     &              be2*wzk(i,j,k)*wzk(i,j,k))
            end do
         end do
      end do
c
      do 19 k=1,NK
      do 21 i=0,NI
         do 31 j=1,NJ
            Jifc(i,j,k)= 0.5d0*(Jac(i,j,k)+ Jac(i+1,j,k))
            gi(i,j,k,1)= 0.5d0*(g11(i,j) +g11(i+1,j))*Jifc(i,j,k)
            gi(i,j,k,2)= 0.5d0*(g12(i,j) +g12(i+1,j))*Jifc(i,j,k)
            gqi(i,j,k,1)= qpr*gi(i,j,k,1) 
            gqi(i,j,k,2)= qpr*gi(i,j,k,2) 
 31      continue
 21   continue
 19   continue
      do 28 k=1,NK
      do 22 i=1,NI
         do 32 j=0,NJ
            Jjfc(i,j,k)= 0.5d0*(Jac(i,j,k)+ Jac(i,j+1,k))
            gj(i,j,k,1)= 0.5d0*(g12(i,j) +g12(i,j+1))*Jjfc(i,j,k)
            gj(i,j,k,2)= 0.5d0*(g22(i,j) +g22(i,j+1))*Jjfc(i,j,k)
            gqj(i,j,k,1)= qpr*gj(i,j,k,1) 
            gqj(i,j,k,2)= qpr*gj(i,j,k,2) 
 32      continue
 22   continue
 28   continue
c
      do 30 j=1,NJ
         do 30 i=0,NI
            do 30 k=1,NK
               gi3(i,j,k)= 0.5d0*(g13(i,j,k) +g13(i+1,j,k))*Jifc(i,j,k)
               gqi3(i,j,k)= qpr*gi3(i,j,k) 
 30   continue
      do 40 j=0,NJ
         do 40 i=1,NI
            do 40 k=1,NK
               gj3(i,j,k)= 0.5d0*(g23(i,j,k) +g23(i,j+1,k))*Jjfc(i,j,k)
               gqj3(i,j,k)= qpr*gj3(i,j,k) 
 40   continue
c
      return
      end
      


