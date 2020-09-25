      subroutine vhydro(dtime)
c     -----------------------------------------------
c     openbc
c     computes, u,v,w -tilde and writes these over cxf,cyf,czf.
c     czf is unchanged.
      implicit logical (a-z)
      integer i,j,k
      include 'header.f'
c      double precision v1(0:NI+1,0:NJ+1,0:NK+1),
c     &     v2(0:NI+1,0:NJ+1,0:NK+1),v3(0:NI+1,0:NJ+1,0:NK+1),
c     &     v4(0:NI+1,0:NJ+1,0:NK+1),v5(0:NI+1,0:NJ+1,0:NK+1)
      double precision hxi,heta,dte,dtime,gradh,kaph1,temp,pz,hdt,
     &     sumuf(0:NI,NJ),sumvf(NI,0:NJ),hx,hy,Jack,wfbct,hmean
      double precision drpx(NI,NJ,NK),drpy(NI,NJ,NK),
     &     grpifc(0:NI,NJ,NK),grpjfc(NI,0:NJ,NK)
      common/rpgrads/drpx,drpy,grpifc,grpjfc
      dte= dtime/EPS
      kaph1= 1.d0 -kappah
c
      do 10 j=1,NJ
         do 20 i=0,NI
            hxi= h(i+1,j)- h(i,j)
            heta= 0.25*(h(i+1,j+1) +h(i,j+1) -h(i+1,j-1)
     &            -h(i,j-1))
            sumuf(i,j)= 0.d0
            do 30 k=1,NK
               hx= gi(i,j,k,1)*hxi +gi(i,j,k,2)*heta
               gradh= gpr*(kappah*hx + kaph1*hxn(i,j,k))
               cxf(i,j,k)= cxf(i,j,k) -dte*(gradh +sifc(i,j,k))
               sumuf(i,j)= sumuf(i,j) + kappah*cxf(i,j,k)
     &              + kaph1*uf(i,j,k)
               hxn(i,j,k)= hx
 30         continue
 20      continue
 10   continue
c     ************added on Mar 21,94 ***********
      goto 201 
      do 15 j=1,NJ
c         sumuf(0,j)= 0.d0
         do  16 k=1,NK
            cxf(0,j,k)= ufbcw(j,k)
c            sumuf(0,j)= sumuf(0,j) + kappah*cxf(0,j,k)
c     &           + kaph1*uf(0,j,k)
 16      continue
c         sumuf(NI,j)= 0.d0
         do 17 k=1,NK
            cxf(NI,j,k)= ufbce(j,k)
c            sumuf(NI,j)= sumuf(NI,j) + kappah*cxf(NI,j,k)
c     &           + kaph1*uf(NI,j,k)
 17      continue
 15   continue
c     **************************************
c
 201  do 40 i=1,NI
         do 50 j=0,NJ
            hxi= 0.25*(h(i+1,j+1) +h(i+1,j) -h(i-1,j+1) -h(i-1,j))
            heta= h(i,j+1) -h(i,j)
            sumvf(i,j)= 0.d0
            do 60 k=1,NK
               hy= gj(i,j,k,1)*hxi +gj(i,j,k,2)*heta
               gradh= gpr*(kappah*hy + kaph1*hyn(i,j,k))
               cyf(i,j,k)= cyf(i,j,k) -dte*(gradh +sjfc(i,j,k))
               sumvf(i,j)= sumvf(i,j) + kappah*cyf(i,j,k)
     &              + kaph1*vf(i,j,k)
               hyn(i,j,k)= hy
 60         continue
 50      continue
 40   continue
c
c     ************added on Mar 21,94 ***********
      goto 202
      do 25 i=1,NI
c=         sumvf(i,NJ)= 0.d0
         do  26 k=1,NK
            cyf(i,NJ,k)= vfbcn(i,k)
c=            sumvf(i,NJ)= sumvf(i,NJ) + kappah*cyf(i,NJ,k)
c=     &           + kaph1*vf(i,NJ,k)
 26      continue
c=         sumvf(i,0)= 0.d0
         do 27 k=1,NK
            cyf(i,0,k)= vfbcs(i,k)
c=            sumvf(i,0)= sumvf(i,0) + kappah*cyf(i,0,k)
c=     &           + kaph1*vf(i,0,k)
 27      continue
 25   continue
c     **** added July 1, 2002
c= 202  goto 401 remove this to reinsert re-eval of h to prevent drift 7/2/04
c     **************************************
c     Evaluate b.c. on wf(at k=NK) (excludes wt). Bc on wf(k=0)= 0
c     This would be the Bc if the residual of the h-eqn were zero
c     For the edges (i=0,NI+1,j=0,NJ+1) hdt is computed in hsolve.
 202  temp= dtime/HDL
      hmean= 0.d0
      do 310 j=1,NJ
         do 310 i=1,NI
            wfbct= -(sumuf(i,j) -sumuf(i-1,j) +sumvf(i,j)
     &           -sumvf(i,j-1)) +wfbcb(i,j)
c     reinsert this to prevent drifting
            hdt= wfbct/J2d(i,j)
c            hdt(i,j)= wfbct/J2d(i,j)
c     hdt is actually (dh'/dt')*(H/D)
c     coorect h(i,j), so that the residual does not accumulate
c     reinsert to prevent h drifting
            h(i,j)= oldh(i,j) + temp*hdt
c==            h(i,j)= oldh(i,j) + temp*hdt(i,j)
            hmean= hmean +h(i,j)
 310  continue
      hmean  = hmean/dble(NI*NJ)
c     remove hmean from h (added 7/13/04 to remove drift)
      do j=1,NJ
         do i=1,NI
            h(i,j)= h(i,j) - hmean
         end do
      end do

c     added 12/99
      call hfill(dtime,h)
c
c      do 330 i=1,NI
c         hdt(i,0)= hdt(i,1)
c         hdt(i,NJ+1)= hdt(i,NJ)
c 330  continue
c      do 320 j=0,NJ+1
c         hdt(0,j)= hdt(NI,j)
c         hdt(NI+1,j)= hdt(1,j)
c 320  continue
c
c-      call findztop (no longer use this - now move the whole grid)
c===
c     no longer call this here - remain explicit in wz,wx,wy,Jac -
c     these are called only before the next time step in momentum.
c===      call findzall
c===      call sigma
c
c     Compute czf 
c     If HY, then pz=0, skfc=0
c
 401  do j=1,NJ
         do i=1,NI
c czf must be 0 only in intpol, not here.  czf(i,j,0)= wfbcb(i,j)
c czf will become 0 if HY. Else, it may be non-zero at this stage,but 0 finally
c            if (i.eq.16)          write(200,*) 'j = ',j
            do k=0,NK
               pz= (p(i,j,k+1) -p(i,j,k))*gqk(i,j,k,3) +
     &              0.25*(p(i+1,j,k+1)
     &              +p(i+1,j,k)-p(i-1,j,k+1)-p(i-1,j,k))*gqk(i,j,k,1)
     &              +0.25*(p(i,j+1,k+1)+p(i,j+1,k)-p(i,j-1,k+1)
     &              -p(i,j-1,k))*gqk(i,j,k,2)
c==
c               if ((k.eq.NK).and.(j.eq.40)) then 
c                  pz= (p(i,j,k+1) -p(i,j,k))*gqk(i,j,k,3) +
c     &                 0.25*(p(i+1,j,k+1)+p(i+1,j,k)
c     &                 -p(i-1,j,k+1)-p(i-1,j,k))*gqk(i,j,k,1)
c     &                 +0.25*(p(i,j+1,k+1)+p(i,j+1,k)-p(i,j-1,k+1)
c     &                 -p(i,j-1,k))*gqk(i,j,k,2)
c
c                  write(6,*)
c     &              i,czf(i,j,k),pz,skfc(i,j,k),
c     &              czf(i,j,k)-dte*(pz+skfc(i,j,k))
c               end if
c==
               czf(i,j,k)= czf(i,j,k) -dte*(pz +skfc(i,j,k))
c               ufdu= (cxf(i,j,k)-cxf(i-1,j,k))
c               vfdv= (cyf(i,j,k)-cyf(i,j-1,k))
c               czf(i,j,k)= czf(i,j,k-1) - ufdu -vfdv
c               if (i.eq.16)  then
c                  write(200,*) pz,skfc(i,j,k),czf(i,j,k)
c               endif
            end do
         end do
      end do
c      write(6,*) 'writing done in vhydro'

c     call uvchy after czf is computed
      call uvchy(dtime)  ! updates cell-centered velocities with hydrostatic

c8**************
c      do k=1,NK
c         do j=1,NJ
c            do i=1,NI
c               v1(i,j,k)= cxf(i,j,k) - cxf(i-1,j,k)
c               v2(i,j,k)= cyf(i,j,k)- cyf(i,j-1,k)
c               v3(i,j,k)= czf(i,j,k)
c               v4(i,j,k)= hxn(i,j,k)
c               v5(i,j,k)= hyn(i,j,k)- hyn(i,j-1,k)
c            end do
c         end do
c      end do
c
c      call outarray(v1,v2,v3,v4,v5)
c      write(6,*) 'stopping in vhydro'
c      stop
c********************
      return
 101  continue

c     U-tilde and V-tilde are computed from the equation in vhydro.f
c     W-tilde is found from interpolation (and stored in czf)
      do 210 k=1,NK-1
         do 220 j=1,NJ
            do 230 i=1,NI
               Jack= 0.5d0*(Jac(i,j,k) +Jac(i,j,k+1))
               czf(i,j,k)= 0.5*Jack*( wx(i,j,k)*cx(i,j,k) +
     &              wx(i,j,k+1)*cx(i,j,k+1) +
     &              wy(i,j,k)*cy(i,j,k) +wy(i,j,k+1)*cy(i,j,k+1) +
     &              EPS*(wz(i,j,k)*cz(i,j,k) 
     &              +wz(i,j,k+1)*cz(i,j,k+1)) )
 230        continue
 220     continue
 210  continue
c
c     boundaries
      k= NK
      do 240 j=1,NJ
         do 250 i=1,NI
            czf(i,j,0)= wfbcb(i,j)
c*            czf(i,j,NK)= wfbct(i,j)
c     Use cx,cy,cz from k=NK, but use mean of wx,wy,wz,Jac in case 
c     grid is sretched.
               Jack= 0.5d0*(Jac(i,j,k) +Jac(i,j,k+1))
               czf(i,j,k)= 0.5*Jack*( (wx(i,j,k) +
     &              wx(i,j,k+1))*cx(i,j,k) +
     &              (wy(i,j,k) +wy(i,j,k+1))*cy(i,j,k) +
     &              EPS*((wz(i,j,k) +wz(i,j,k+1))*cz(i,j,k) ) )
c-            czf(i,j,NK)= Jac(i,j,NK)*( wx(i,j,NK)*cx(i,j,NK) +
c-     &           wy(i,j,NK)*cy(i,j,NK) +
c-     &           EPS*wz(i,j,NK)*cz(i,j,NK) )
 250     continue
 240  continue
c
      return
      end
      
            





