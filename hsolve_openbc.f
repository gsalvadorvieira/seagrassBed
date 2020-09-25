      subroutine hsolve(h,oldh,dtime)
c     --------------
c     openbc
c     ----------------------------------
      implicit logical (a-z)
      include 'header_mod.f'
      integer iter,i,j,kount,l, imax, jmax
      double precision ch(9,NI,NJ),rhs(NI,NJ),
     &     oldh(0:NI+1,0:NJ+1),h(0:NI+1,0:NJ+1),res,
     &     dti,dtime,tol,maxres,hstar,rlx
c
cc      tol= 1.d-8
c     tol of 1.d-7 doesn't seem to change the solution. 1.d-6 would be ok
      tol=1.d-12
      dti= 1.d0/dtime
      kount=0
c
c     maxres= 0.d0
c     res= 0.d0
c     iter=0
c     write(6,*) 'in hsolve',tol,dti,maxres,res,iter
c
c     save old value of h
c
c     compute fine grid coefficients and rhs on fine grid
      call chfine(dtime,ch,rhs)
c      call checkch(1,NI,NJ,ch)
c     **********************************
c      write(6,*) 'enter rlx'
c      read(5,*) rlx
c      open (unit=90,file='test.out')
c      do j=1,NJ
c         write(90,*) 'j=',j
c         do i=1,NI
c            write(90,*) (ch(l,i,j),l=1,9),rhs(i,j)
c         end do
c      end do
c      close(90)
c      stop

      rlx= 1.72
C 999  do 1000 iter=1,100
 999  do 1000 iter=1,100
C$DOACROSSSHARE(ch,h,rlx,rhs),
C$&LOCAL(i,j,hstar)
         do 120 j=1,NJ
            do 130 i=1,NI
               hstar= ( ch(2,i,j)*h(i+1,j)
     &              +ch(3,i,j)*h(i-1,j)
     &              +ch(4,i,j)*h(i,j+1)
     &              +ch(5,i,j)*h(i,j-1)
     &              +ch(6,i,j)*h(i+1,j+1)
     &              +ch(7,i,j)*h(i-1,j+1)
     &              +ch(8,i,j)*h(i+1,j-1)
     &              +ch(9,i,j)*h(i-1,j-1)
     &              - rhs(i,j) )/(-ch(1,i,j))
               h(i,j)= (1.d0 -rlx)*h(i,j) +rlx*hstar
 130        continue
 120     continue
c         write(6,*) 'writing h'
c         do j=2,2
c            do i=1,10
c               write(6,*) h(i,j)
c            end do
c         end do
c         write(6,*) 'stopping in h'
c         stop
c
         do 140 l=1,3
            call hfill(dtime,h)
 140     continue

         maxres= 0.d0
C$DOACROSSSHARE(ch,h,rhs,maxres),
C$&LOCAL(i,j,res)
         do 220 j=1,NJ
            do 230 i=1,NI
c     res = b - A x'
               res = rhs(i,j)-
     &              (ch(1,i,j)*h(i,j)
     &              +ch(2,i,j)*h(i+1,j)
     &              +ch(3,i,j)*h(i-1,j)
     &              +ch(4,i,j)*h(i,j+1)
     &              +ch(5,i,j)*h(i,j-1)
     &              +ch(6,i,j)*h(i+1,j+1)
     &              +ch(7,i,j)*h(i-1,j+1)
     &              +ch(8,i,j)*h(i+1,j-1)
     &              +ch(9,i,j)*h(i-1,j-1))
               if (dabs(res).gt.maxres) then
                  maxres= dabs(res)
                  imax = i
                  jmax = j
               end if
 230        continue
 220     continue

         if (maxres.gt.3000.d0) then
            write(6,'(A, I8, I8, ES13.6)')
     &    'STOP in hsolve. res too large. imax,jmax,maxres=',
     &    imax,jmax,maxres

            call outcdf(xc,yc,zc,p,h,s,T,u,v,w,vor,uf,vf,wf,theta,xg,zg,
     &          ght,gdef,step)
            
            stop
         end if
c     may have to call hfill a few times if g12 is non-zero
c     if g12 is zero, we may call hfill after convergence is
c     reached or may not call it at all.
c+         do 240 l=1,3
c+            call  hfill(dtime,h)
c+ 240     continue
c         call mprove(h,ch,rhs,dtime)
         if (maxres.lt.tol) goto 13

 1000 continue

 13   write(6,'(A,I8,A,ES13.6)') 'iter = ',iter,', maxres = ',maxres
c      call  hfill(dtime,h)
c+      kount=kount+1
c
c     Bring the Coriolis terms up to the (n+1) time level. New sifc, sjfc
c+      if (kount.ge.3) goto 301
c+      call uvinterm(h,dtime)
c+      call coriolis(1)
c+      call srcface(1)
c+      call newrhs(dtime,rhs)
c+      goto 999
c
cc     hdt is computed in vhydro for i=1,NI,j=1,NJ... but for outer points
cc     hdt is actually (dh'/dt')*(H/D)
c 301  do 310 j=0,NJ+1
c         hdt(0,j)= dti*( h(0,j) -oldh(0,j) )*HDL
c         hdt(NI+1,j)= dti*( h(NI+1,j) -oldh(NI+1,j) )*HDL
c 310  continue
c      do 320 i=1,NI
c         hdt(i,0)= dti*( h(i,0) -oldh(i,0))*HDL
c         hdt(i,NJ+1)= dti*( h(i,NJ+1) -oldh(i,NJ+1))*HDL
c 320  continue
c
      do l=1,1
         call mprove(h,ch,rhs,dtime)
      end do

      return

      end
