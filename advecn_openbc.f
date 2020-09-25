      subroutine advecn(m,n,dtime,step)
c     ---------------------------------------------
c     convecn(m,n,dtime) means use level m and write s,T to level n
c     computes Cx,Cy,Cz at the cell centers.
c     uflx,vflx,wflx,sflx,Tflx are the fluxes defined at
c     cell faces (using QUICK).
c
      implicit logical (a-z)
      include 'header.f' ! uss is loaded here
      integer i,j,k,n,m,step,it,counter,i0
c     ntr is the number of tracers in variable T (Defined in header.f)
      double precision Eighth,Half,absvel,left,right,ctr,upos,uneg,
     &     SevenEighths

      double precision dtime,uflx(0:NI,0:NJ,0:NK),
     &     vflx(0:NI,0:NJ,0:NK),wflx(0:NI,0:NJ,0:NK),
     &     sflx(0:NI,0:NJ,0:NK),Tflx(ntr,0:NI,0:NJ,0:NK),
     &     usx(0:NI+1,0:NJ+1,0:NK+1),uTx(ntr,0:NI+1,0:NJ+1,0:NK+1),
     &     uux(0:NI+1,0:NJ+1,0:NK+1),uvx(0:NI+1,0:NJ+1,0:NK+1),
     &     uwx(0:NI+1,0:NJ+1,0:NK+1),
     &     sdif(NI,NJ,NK),Tdif(NI,NJ,NK),udif(NI,NJ,NK),
     &     vdif(NI,NJ,NK),wdif(NI,NJ,NK),dtJ
      double precision dTx(ntr,0:NI+1),dsx(0:NI+1),dux(0:NI+1),
     &     dvx(0:NI+1),dwx(0:NI+1),dTy(ntr,0:NJ),dsy(0:NJ),
     &     duy(0:NJ),dvy(0:NJ),dwy(0:NJ),dTz(ntr,0:NK+1),dsz(0:NK+1),
     &     duz(0:NK+1),dvz(0:NK+1),dwz(0:NK+1)
      double precision s_restore(0:NJ+1,0:NK+1),alph_restore(0:NJ+1)
      common/restore/s_restore,alph_restore
      double precision courinv,fc
      double precision sbold(NJ,NK),Tbold(ntr,NJ,NK),dJtinv,du
      double precision dampAmp, dampFrac, ussmean, uoutmean, flux_ratio

      parameter (Eighth=0.125d0, Half=0.5d0, courinv=2.d0)
C     AUTOMATIC ENTRY ON LINE 35
      parameter(dampAmp=5.d-2, dampFrac=60.d0)

      double precision restore_rate(0:NI+1), u_outlet(NK), w_outlet(NK)

      do i=1,NI

C       STEP FUNCTION AT FRACTION
C         if (i.gt.NI*dampFrac) then
C           restore_rate(i) = dampAmp
C         else 
C           restore_rate(i) = 0.d0
C         end if

C       LINEAR FROM FRACTION
C         if (i.gt.NI*dampFrac) then
C           restore_rate(i) = dampAmp/((1.d0-dampFrac)*dble(NI))*
C      &          (dble(i)-dble(NI)*dampFrac)     
C         else 
C           restore_rate(i) = 0.d0
C         end if

C       TANH STARTING AT FRAC AND STEEPEST BETWEEN THERE AND OUTLET
        restore_rate(i) = dampAmp*0.5d0*(1.d0 + 
     &     TANH(15.d0/288.d0/2.d0*(dble(i) - (dble(NI)-dampFrac))))

      end do

      do k=1,NK
         do j=1,NJ

c     for Inflow/Outflow boundaries
c     dTx,dsx,.... at faces...
            i=0
c     s,v,w are zero at inflow, u=uss and T=tinfl at inflow
            do it=1,ntr
               dTx(it,i)= T(it,1,j,k,m) - tinfl(k)
            end do
            dsx(i)= s(1,j,k,m) - 0.d0
            dux(i)= u(1,j,k,m) - uss(k)
            dvx(i)= v(1,j,k,m) - 0.d0
            dwx(i)= w(1,j,k,m) - 0.d0
c     outflow boundaries  : do i=1,NI  (instead of i=1,NI-1)
c     Call to openbc puts data at NI+1
            do i=1,NI
               do it=1,ntr
                  dTx(it,i)= T(it,i+1,j,k,m) - T(it,i,j,k,m)
               end do
               dsx(i)= s(i+1,j,k,m) - s(i,j,k,m)
               dux(i)= u(i+1,j,k,m) - u(i,j,k,m)
               dvx(i)= v(i+1,j,k,m) - v(i,j,k,m)
               dwx(i)= w(i+1,j,k,m) - w(i,j,k,m)
            end do

c            do i=1,NI ! periodic boundaries
c     do i=1,NI-1  (when bcs specified) (going up to NI gives NAN)
            do i=1,NI-1
               absvel= dabs(uf(i,j,k))
               upos= Half*(uf(i,j,k) +absvel)
               uneg= Half*(uf(i,j,k) -absvel)

               do it=1,ntr
               left= Eighth*(dTx(it,i) -dTx(it,i-1))
               right= Eighth*(dTx(it,i+1) -dTx(it,i))
               ctr= Half* (T(it,i+1,j,k,m) +T(it,i,j,k,m))
               Tflx(it,i,j,k) = uf(i,j,k)*ctr - upos*left -uneg*right
               end do

               left= Eighth*(dsx(i) -dsx(i-1))
               right= Eighth*(dsx(i+1) -dsx(i))
               ctr= Half* (s(i+1,j,k,m) +s(i,j,k,m))
               sflx(i,j,k) = uf(i,j,k)*ctr - upos*left -uneg*right

               left= Eighth*(dux(i) -dux(i-1))
               right= Eighth*(dux(i+1) -dux(i))
               ctr= Half* (u(i+1,j,k,m) +u(i,j,k,m))
               uflx(i,j,k) = uf(i,j,k)*ctr - upos*left -uneg*right

               left= Eighth*(dvx(i) -dvx(i-1))
               right= Eighth*(dvx(i+1) -dvx(i))
               ctr= Half* (v(i+1,j,k,m) +v(i,j,k,m))
               vflx(i,j,k) = uf(i,j,k)*ctr - upos*left -uneg*right

               left= Eighth*(dwx(i) -dwx(i-1))
               right= Eighth*(dwx(i+1) -dwx(i))
               ctr= Half* (w(i+1,j,k,m) +w(i,j,k,m))
               wflx(i,j,k) = uf(i,j,k)*ctr - upos*left -uneg*right

            end do
c     Inflow boundary conditions
            i=0
            sflx(i,j,k) = 0.d0
            do it=1,ntr
               Tflx(it,i,j,k) = tinfl(k)*uf(i,j,k)
            end do
            uflx(i,j,k) = uss(k)*uf(i,j,k)
            vflx(i,j,k) = 0.d0
            wflx(i,j,k) = 0.d0
c     Outflow boundary
            i=NI
            sflx(i,j,k) = s(i,j,k,m)*uf(i,j,k)
            do it=1,ntr
               Tflx(it,i,j,k) = T(it,i,j,k,m)*uf(i,j,k)
            end do
            uflx(i,j,k) = u(i,j,k,m)*uf(i,j,k)
            vflx(i,j,k) = v(i,j,k,m)*uf(i,j,k)
            wflx(i,j,k) = w(i,j,k,m)*uf(i,j,k)
         end do
      end do

      do k=1,NK
         do j=1,NJ
            do i=1,NI
               uux(i,j,k)=  uflx(i,j,k) -uflx(i-1,j,k)
               uvx(i,j,k)=  vflx(i,j,k) -vflx(i-1,j,k)
               uwx(i,j,k)=  wflx(i,j,k) -wflx(i-1,j,k)
               usx(i,j,k)=  sflx(i,j,k) -sflx(i-1,j,k)
               do it=1,ntr
                  uTx(it,i,j,k)=  Tflx(it,i,j,k) -Tflx(it,i-1,j,k)
               end do
            end do
         end do
      end do

c     y-direction
c     -----------
      do k=1,NK
         do i=1,NI
c     dTy,dsy,.... at faces...
            do j=0,NJ,NJ
               do it=1,ntr
                  dTy(it,j)= 0.d0
               end do
               dsy(j)= 0.d0
               duy(j)= 0.d0
               dvy(j)= 0.d0
               dwy(j)= 0.d0
            end do
            do j=1,NJ-1
               do it=1,ntr
                  dTy(it,j)= T(it,i,j+1,k,m) - T(it,i,j,k,m)
               end do
               dsy(j)= s(i,j+1,k,m) - s(i,j,k,m)
               duy(j)= u(i,j+1,k,m) - u(i,j,k,m)
               dvy(j)= v(i,j+1,k,m) - v(i,j,k,m)
               dwy(j)= w(i,j+1,k,m) - w(i,j,k,m)
            end do
            s(i,0,k,m)= s(i,1,k,m)
            do it=1,ntr
               T(it,i,0,k,m)= T(it,i,1,k,m)
            end do
            u(i,0,k,m)= u(i,1,k,m)
            v(i,0,k,m)= v(i,1,k,m)
            w(i,0,k,m)= w(i,1,k,m)

            s(i,NJ+1,k,m)= s(i,NJ,k,m)
            do it=1,ntr
               T(it,i,NJ+1,k,m)= T(it,i,NJ,k,m)
            end do
            u(i,NJ+1,k,m)= u(i,NJ,k,m)
            v(i,NJ+1,k,m)= v(i,NJ,k,m)
            w(i,NJ+1,k,m)= w(i,NJ,k,m)

            do j=1,NJ-1

               if (vf(i,j,k).ge.0.d0) then

                  do it=1,ntr
                     left= Eighth*(dTy(it,j) -dTy(it,j-1))
                     ctr= Half* (T(it,i,j+1,k,m) +T(it,i,j,k,m))
                     fc= ctr -left
                     call ultim(T(it,i,j-1,k,m),T(it,i,j,k,m),
     &                    T(it,i,j+1,k,m),dTy(it,j),dTy(it,j-1),
     &                    courinv,fc)
                     Tflx(it,i,j,k)= vf(i,j,k)*fc
                  end do

                  left= Eighth*(dsy(j) -dsy(j-1))
                  ctr= Half* (s(i,j+1,k,m) +s(i,j,k,m))
                  fc= ctr -left
                  call ultim(s(i,j-1,k,m),s(i,j,k,m),
     &                 s(i,j+1,k,m),dsy(j),dsy(j-1),courinv,fc)
                  sflx(i,j,k)= vf(i,j,k)*fc

                  left= Eighth*(duy(j) -duy(j-1))
                  ctr= Half* (u(i,j+1,k,m) +u(i,j,k,m))
                  fc= ctr -left
                  call ultim(u(i,j-1,k,m),u(i,j,k,m),
     &                 u(i,j+1,k,m),duy(j),duy(j-1),courinv,fc)
                  uflx(i,j,k)= vf(i,j,k)*fc

                  left= Eighth*(dvy(j) -dvy(j-1))
                  ctr= Half* (v(i,j+1,k,m) +v(i,j,k,m))
                  fc= ctr -left
                  call ultim(v(i,j-1,k,m),v(i,j,k,m),
     &                 v(i,j+1,k,m),dvy(j),dvy(j-1),courinv,fc)
                  vflx(i,j,k)= vf(i,j,k)*fc

                  left= Eighth*(dwy(j) -dwy(j-1))
                  ctr= Half* (w(i,j+1,k,m) +w(i,j,k,m))
                  fc= ctr -left
                  call ultim(w(i,j-1,k,m),w(i,j,k,m),
     &                 w(i,j+1,k,m),dwy(j),dwy(j-1),courinv,fc)
                  wflx(i,j,k)= vf(i,j,k)*fc

               else

                  do it=1,ntr
                     right= Eighth*(dTy(it,j+1) -dTy(it,j))
                     ctr= Half* (T(it,i,j+1,k,m) +T(it,i,j,k,m))
                     fc= ctr - right
                     call ultim(T(it,i,j+2,k,m),T(it,i,j+1,k,m),
     &                    T(it,i,j,k,m),dTy(it,j),dTy(it,j+1),
     &                    courinv,fc)
                     Tflx(it,i,j,k)= vf(i,j,k)*fc
                  end do

                  right= Eighth*(dsy(j+1) -dsy(j))
                  ctr= Half* (s(i,j+1,k,m) +s(i,j,k,m))
                  fc= ctr - right
                  call ultim(s(i,j+2,k,m),s(i,j+1,k,m),
     &                 s(i,j,k,m),dsy(j),dsy(j+1),courinv,fc)
                  sflx(i,j,k)= vf(i,j,k)*fc

                  right= Eighth*(duy(j+1) -duy(j))
                  ctr= Half* (u(i,j+1,k,m) +u(i,j,k,m))
                  fc= ctr - right
                  call ultim(u(i,j+2,k,m),u(i,j+1,k,m),
     &                 u(i,j,k,m),duy(j),duy(j+1),courinv,fc)
                  uflx(i,j,k)= vf(i,j,k)*fc

                  right= Eighth*(dvy(j+1) -dvy(j))
                  ctr= Half* (v(i,j+1,k,m) +v(i,j,k,m))
                  fc= ctr - right
                  call ultim(v(i,j+2,k,m),v(i,j+1,k,m),
     &                 v(i,j,k,m),dvy(j),dvy(j+1),courinv,fc)
                  vflx(i,j,k)= vf(i,j,k)*fc

                  right= Eighth*(dwy(j+1) -dwy(j))
                  ctr= Half* (w(i,j+1,k,m) +w(i,j,k,m))
                  fc= ctr - right
                  call ultim(w(i,j+2,k,m),w(i,j+1,k,m),
     &                 w(i,j,k,m),dwy(j),dwy(j+1),courinv,fc)
                  wflx(i,j,k)= vf(i,j,k)*fc

               end if

c               absvel= dabs(vf(i,j,k))
c               upos= Half*(vf(i,j,k) +absvel)
c               uneg= Half*(vf(i,j,k) -absvel)
c
c               left= Eighth*(dTy(j) -dTy(j-1))
c               right= Eighth*(dTy(j+1) -dTy(j))
c               ctr= Half* (T(i,j+1,k,m) +T(i,j,k,m))
c               Tflx(i,j,k) = vf(i,j,k)*ctr - upos*left -uneg*right
c
c               left= Eighth*(dsy(j) -dsy(j-1))
c               right= Eighth*(dsy(j+1) -dsy(j))
c               ctr= Half* (s(i,j+1,k,m) +s(i,j,k,m))
c               sflx(i,j,k) = vf(i,j,k)*ctr - upos*left -uneg*right
c
c               left= Eighth*(duy(j) -duy(j-1))
c               right= Eighth*(duy(j+1) -duy(j))
c               ctr= Half* (u(i,j+1,k,m) +u(i,j,k,m))
c               uflx(i,j,k) = vf(i,j,k)*ctr - upos*left -uneg*right
c
c               left= Eighth*(dvy(j) -dvy(j-1))
c               right= Eighth*(dvy(j+1) -dvy(j))
c               ctr= Half* (v(i,j+1,k,m) +v(i,j,k,m))
c               vflx(i,j,k) = vf(i,j,k)*ctr - upos*left -uneg*right
c
c               left= Eighth*(dwy(j) -dwy(j-1))
c               right= Eighth*(dwy(j+1) -dwy(j))
c               ctr= Half* (w(i,j+1,k,m) +w(i,j,k,m))
c               wflx(i,j,k) = vf(i,j,k)*ctr - upos*left -uneg*right

            end do
c     Solid Boundaries
            do j=0,NJ,NJ
               do it=1,ntr
                  Tflx(it,i,j,k)= 0.d0
               end do
               sflx(i,j,k)= 0.d0
               uflx(i,j,k)= 0.d0
               vflx(i,j,k)= 0.d0
               wflx(i,j,k)= 0.d0
            end do
         end do
      end do

      do k=1,NK
         do i=1,NI
            do j=1,NJ
               uux(i,j,k)= (uflx(i,j,k) -uflx(i,j-1,k) ) + uux(i,j,k)
               uvx(i,j,k)= (vflx(i,j,k) -vflx(i,j-1,k) ) + uvx(i,j,k)
               uwx(i,j,k)= (wflx(i,j,k) -wflx(i,j-1,k) ) + uwx(i,j,k)
               usx(i,j,k)= (sflx(i,j,k) -sflx(i,j-1,k) ) + usx(i,j,k)
               do it=1,ntr
                  uTx(it,i,j,k)= (Tflx(it,i,j,k) -Tflx(it,i,j-1,k) )
     &                 + uTx(it,i,j,k)
               end do
            end do
         end do
      end do

c     z-direction
c     -----------
      do j=1,NJ
         do i=1,NI
c     dTz,dsz,.... at faces...
c            do k=0,NK,NK
            k=0
            do it=1,ntr
               dTz(it,k)= 0.d0
            end do
            dsz(k)= 0.d0
            duz(k)= 0.d0
            dvz(k)= 0.d0
            dwz(k)= 0.d0
c
            s(i,j,k,m)= s(i,j,k+1,m)
            do it=1,ntr
               T(it,i,j,k,m)= T(it,i,j,k+1,m)
            end do
            u(i,j,k,m)= u(i,j,k+1,m)
            v(i,j,k,m)= v(i,j,k+1,m)
            w(i,j,k,m)= w(i,j,k+1,m)
c
            do k=1,NK-1
               do it=1,ntr
                  dTz(it,k)= T(it,i,j,k+1,m) - T(it,i,j,k,m)
               end do
               dsz(k)= s(i,j,k+1,m) - s(i,j,k,m)
               duz(k)= u(i,j,k+1,m) - u(i,j,k,m)
               dvz(k)= v(i,j,k+1,m) - v(i,j,k,m)
               dwz(k)= w(i,j,k+1,m) - w(i,j,k,m)
            end do
c-c     Top boundary - linear extrapolation, constant gradient
c-            do k=NK,NK+1
c-               dTz(k)= dtz(k-1)
c-               dsz(k)= dsz(k-1)
c-               duz(k)= duz(k-1)
c-               dvz(k)= dvz(k-1)
c-               dwz(k)= dwz(k-1)
c-            end do
c-c     Linear extrapolation
c-            s(i,j,NK+1,m)= 2.d0*s(i,j,NK,m) - s(i,j,NK-1,m)
c-            T(i,j,NK+1,m)= 2.d0*T(i,j,NK,m) - T(i,j,NK-1,m)
c-            u(i,j,NK+1,m)= 2.d0*u(i,j,NK,m) - u(i,j,NK-1,m)
c-            v(i,j,NK+1,m)= 2.d0*v(i,j,NK,m) - v(i,j,NK-1,m)
c-            w(i,j,NK+1,m)= 2.d0*w(i,j,NK,m) - w(i,j,NK-1,m)
c
c     Top boundary - Zero gradient
            do k=NK,NK+1
               do it=1,ntr
                  dTz(it,k)= 0.d0
               end do
               dsz(k)= 0.d0
               duz(k)= 0.d0
               dvz(k)= 0.d0
               dwz(k)= 0.d0
            end do
c     Linear extrapolation
            s(i,j,NK+1,m)= s(i,j,NK,m)
            do it=1,ntr
               T(it,i,j,NK+1,m)= T(it,i,j,NK,m)
            end do
            u(i,j,NK+1,m)= u(i,j,NK,m)
            v(i,j,NK+1,m)= v(i,j,NK,m)
            w(i,j,NK+1,m)= w(i,j,NK,m)

            do k=1,NK

               if (wf(i,j,k).ge.0.d0) then

                  do it=1,ntr
                     left= Eighth*(dTz(it,k) -dTz(it,k-1))
                     ctr= Half* (T(it,i,j,k+1,m) +T(it,i,j,k,m))
                     fc= ctr -left
                     call ultim(T(it,i,j,k-1,m),T(it,i,j,k,m),
     &                    T(it,i,j,k+1,m),dTz(it,k),dTz(it,k-1),
     &                    courinv,fc)
                     Tflx(it,i,j,k)= wf(i,j,k)*fc
                  end do

                  left= Eighth*(dsz(k) -dsz(k-1))
                  ctr= Half* (s(i,j,k+1,m) +s(i,j,k,m))
                  fc= ctr -left
                  call ultim(s(i,j,k-1,m),s(i,j,k,m),
     &                 s(i,j,k+1,m),dsz(k),dsz(k-1),courinv,fc)
                  sflx(i,j,k)= wf(i,j,k)*fc

                  left= Eighth*(duz(k) -duz(k-1))
                  ctr= Half* (u(i,j,k+1,m) +u(i,j,k,m))
                  fc= ctr -left
                  call ultim(u(i,j,k-1,m),u(i,j,k,m),
     &                 u(i,j,k+1,m),duz(k),duz(k-1),courinv,fc)
                  uflx(i,j,k)= wf(i,j,k)*fc

                  left= Eighth*(dvz(k) -dvz(k-1))
                  ctr= Half* (v(i,j,k+1,m) +v(i,j,k,m))
                  fc= ctr -left
                  call ultim(v(i,j,k-1,m),v(i,j,k,m),
     &                 v(i,j,k+1,m),dvz(k),dvz(k-1),courinv,fc)
                  vflx(i,j,k)= wf(i,j,k)*fc

                  left= Eighth*(dwz(k) -dwz(k-1))
                  ctr= Half* (w(i,j,k+1,m) +w(i,j,k,m))
                  fc= ctr -left
                  call ultim(w(i,j,k-1,m),w(i,j,k,m),
     &                 w(i,j,k+1,m),dwz(k),dwz(k-1),courinv,fc)
                  wflx(i,j,k)= wf(i,j,k)*fc

               else

                  if (k.eq.NK) then
                     do it=1,ntr
                        Tflx(it,i,j,k)= wf(i,j,k)*Half*
     &                    (T(it,i,j,k+1,m) +T(it,i,j,k,m))
                     end do
                     sflx(i,j,k)= wf(i,j,k)*Half*
     &                    (s(i,j,k+1,m) +s(i,j,k,m))
                     uflx(i,j,k)= wf(i,j,k)*Half*
     &                    (u(i,j,k+1,m) +u(i,j,k,m))
                     vflx(i,j,k)= wf(i,j,k)*Half*
     &                    (v(i,j,k+1,m) +v(i,j,k,m))
                     wflx(i,j,k)= wf(i,j,k)*Half*
     &                    (w(i,j,k+1,m) +w(i,j,k,m))
                     goto 209
                  endif

                  do it=1,ntr
                     right= Eighth*(dTz(it,k+1) -dTz(it,k))
                     ctr= Half* (T(it,i,j,k+1,m) +T(it,i,j,k,m))
                     fc= ctr - right
                     call ultim(T(it,i,j,k+2,m),T(it,i,j,k+1,m),
     &                    T(it,i,j,k,m),dTz(it,k),dTz(it,k+1),
     &                    courinv,fc)
                     Tflx(it,i,j,k)= wf(i,j,k)*fc
                  end do

                  right= Eighth*(dsz(k+1) -dsz(k))
                  ctr= Half* (s(i,j,k+1,m) +s(i,j,k,m))
                  fc= ctr - right
                  call ultim(s(i,j,k+2,m),s(i,j,k+1,m),
     &                 s(i,j,k,m),dsz(k),dsz(k+1),courinv,fc)
                  sflx(i,j,k)= wf(i,j,k)*fc

                  right= Eighth*(duz(k+1) -duz(k))
                  ctr= Half* (u(i,j,k+1,m) +u(i,j,k,m))
                  fc= ctr - right
                  call ultim(u(i,j,k+2,m),u(i,j,k+1,m),
     &                 u(i,j,k,m),duz(k),duz(k+1),courinv,fc)
                  uflx(i,j,k)= wf(i,j,k)*fc

                  right= Eighth*(dvz(k+1) -dvz(k))
                  ctr= Half* (v(i,j,k+1,m) +v(i,j,k,m))
                  fc= ctr - right
                  call ultim(v(i,j,k+2,m),v(i,j,k+1,m),
     &                 v(i,j,k,m),dvz(k),dvz(k+1),courinv,fc)
                  vflx(i,j,k)= wf(i,j,k)*fc

                  right= Eighth*(dwz(k+1) -dwz(k))
                  ctr= Half* (w(i,j,k+1,m) +w(i,j,k,m))
                  fc= ctr - right
                  call ultim(w(i,j,k+2,m),w(i,j,k+1,m),
     &                 w(i,j,k,m),dwz(k),dwz(k+1),courinv,fc)
                  wflx(i,j,k)= wf(i,j,k)*fc

               end if

c               absvel= dabs(wf(i,j,k))
c               upos= Half*(wf(i,j,k) +absvel)
c               uneg= Half*(wf(i,j,k) -absvel)
c
c               left= Eighth*(dTz(k) -dTz(k-1))
c               right= Eighth*(dTz(k+1) -dTz(k))
c               ctr= Half* (T(i,j,k+1,m) +T(i,j,k,m))
c               Tflx(i,j,k) = wf(i,j,k)*ctr - upos*left -uneg*right
c
c               left= Eighth*(dsz(k) -dsz(k-1))
c               right= Eighth*(dsz(k+1) -dsz(k))
c               ctr= Half* (s(i,j,k+1,m) +s(i,j,k,m))
c               sflx(i,j,k) = wf(i,j,k)*ctr - upos*left -uneg*right
c
c               left= Eighth*(duz(k) -duz(k-1))
c               right= Eighth*(duz(k+1) -duz(k))
c               ctr= Half* (u(i,j,k+1,m) +u(i,j,k,m))
c               uflx(i,j,k) = wf(i,j,k)*ctr - upos*left -uneg*right
c
c               left= Eighth*(dvz(k) -dvz(k-1))
c               right= Eighth*(dvz(k+1) -dvz(k))
c               ctr= Half* (v(i,j,k+1,m) +v(i,j,k,m))
c               vflx(i,j,k) = wf(i,j,k)*ctr - upos*left -uneg*right
c
c               left= Eighth*(dwz(k) -dwz(k-1))
c               right= Eighth*(dwz(k+1) -dwz(k))
c               ctr= Half* (w(i,j,k+1,m) +w(i,j,k,m))
c               wflx(i,j,k) = wf(i,j,k)*ctr - upos*left -uneg*right

 209        end do
c     Solid Boundary
            k=0
            do it=1,ntr
               Tflx(it,i,j,k)= 0.d0
            end do
            sflx(i,j,k)= 0.d0
            uflx(i,j,k)= 0.d0
            vflx(i,j,k)= 0.d0
            wflx(i,j,k)= 0.d0
c     Top surface (done already)
c            k=NK
c            Tflx(i,j,k)= T(i,j,k,m)*wf(i,j,k)
c            sflx(i,j,k)= s(i,j,k,m)*wf(i,j,k)
c            uflx(i,j,k)= u(i,j,k,m)*wf(i,j,k)
c            vflx(i,j,k)= v(i,j,k,m)*wf(i,j,k)
c            wflx(i,j,k)= w(i,j,k,m)*wf(i,j,k)
         end do
      end do

      call vdiffusion(sdif,Tdif,udif,vdif,wdif,m)
c=      call biharmonic(udif,vdif,m,step)
c     Save outflow diffusion for use in openbc
c     ----
      do k=1,NK
         do j=1,NJ
            udifout(j,k)= udif(NI,j,k)
            vdifout(j,k)= vdif(NI,j,k)
            wdifout(j,k)= wdif(NI,j,k)
         end do
      end do
c     -------
c
      do 250 k=1,NK
         do 251 j=1,NJ
            do 252 i=1,NI
               uux(i,j,k)= uux(i,j,k) + (uflx(i,j,k) -uflx(i,j,k-1) )
     &              - udif(i,j,k)
               uvx(i,j,k)= uvx(i,j,k) + (vflx(i,j,k) -vflx(i,j,k-1) )
     &              - vdif(i,j,k)
               uwx(i,j,k)= uwx(i,j,k) + (wflx(i,j,k) -wflx(i,j,k-1) )
     &              - wdif(i,j,k)
               usx(i,j,k)= usx(i,j,k) + (sflx(i,j,k) -sflx(i,j,k-1) )
c=     &              - sdif(i,j,k)
               do it=1,ntr
                  uTx(it,i,j,k)= uTx(it,i,j,k) + (Tflx(it,i,j,k)
     &                 -Tflx(it,i,j,k-1) )
cc     &              - Tdif(i,j,k)
               end do
 252        continue
 251     continue
 250  continue
c
      go to 202
c     Horizontal eddy viscosity
      call viscous(udif,vdif,wdif,m)
      do k=1,NK
         do j=1,NJ
            do i=1,NI
               uux(i,j,k)= uux(i,j,k) - udif(i,j,k)
               uvx(i,j,k)= uvx(i,j,k) - vdif(i,j,k)
               uwx(i,j,k)= uwx(i,j,k) - wdif(i,j,k)
c               usx(i,j,k)= usx(i,j,k) - sdif(i,j,k)
c               do it=1,ntr
c                  uTx(it,i,j,k)= uTx(it,i,j,k) - Tdif(i,j,k)
c               end do
            end do
         end do
      end do
c
c     Save old values of s,T to calculate cufs,cufT
c     -----------------------------------------------
      i=NI
      do k=1,NK
         do j=1,NJ
            sbold(j,k) = s(i,j,k,0)
            do it=1,ntr
               Tbold(it,j,k)= T(it,i,j,k,0)
            end do
         end do
      end do
c     -------------------------

C     *$*  ASSERT CONCURRENT CALL

C     Compute mean value of uss (steady state) and averaged outlet profile
 202  ussmean = 0.d0
      uoutmean = 0.d0
      j = NJ/2
      i0 = NI - 2*INT(dampFrac) + 1
      do k=1,NK
         ussmean = ussmean + uss(k)
         counter = 0
         do i = i0,NI
            u_outlet(k) = u_outlet(k) + u(i,j,k,0)
            w_outlet(k) = w_outlet(k) + w(i,j,k,0)
            counter = counter + 1
         end do
         u_outlet(k) = u_outlet(k)/counter
         w_outlet(k) = w_outlet(k)/counter
         uoutmean = uoutmean + u_outlet(k)
      end do
      ussmean = ussmean/NK
      uoutmean = uoutmean/NK
      flux_ratio = ussmean/uoutmean

      do 300 j=1,NJ
         do 301 i=1,NI
            do 302 k=1,NK
               dtJ= dtime/Jac(i,j,k)
C                cx(i,j,k)= u(i,j,k,0) -dtJ*uux(i,j,k)
C                cx(i,j,k)= u(i,j,k,0) -dtJ*uux(i,j,k)
C      &                - restore_rate(i) * (u(i,j,k,0) - uf(1,j,k)) ! wrong, uf is flux
               cx(i,j,k)= u(i,j,k,0) -dtJ*uux(i,j,k)
     &                - restore_rate(i) *
     &                (u(i,j,k,0) - flux_ratio*u_outlet(k)) ! check if works (seems to)

C                cy(i,j,k)= v(i,j,k,0) -dtJ*uvx(i,j,k)
               cy(i,j,k)= v(i,j,k,0) -dtJ*uvx(i,j,k)
     &               - restore_rate(i) * v(i,j,k,0)

C                cz(i,j,k)= w(i,j,k,0) -dtJ*uwx(i,j,k)
               cz(i,j,k)= w(i,j,k,0) -dtJ*uwx(i,j,k)
     &               - restore_rate(i) *
     &               (w(i,j,k,0) - w_outlet(k))

               s(i,j,k,n)= s(i,j,k,0) -dtJ*usx(i,j,k)
c==     &              - alph_restore(j)*(s(i,j,k,m)-s_restore(j))
               do it=1,ntr
                  T(it,i,j,k,n)= T(it,i,j,k,0) -dtJ*uTx(it,i,j,k)
               end do
c     T_restore= 0.d0
 302        continue
 301     continue
 300  continue

c     Now calculate cufs,cufT
c     ----------------------------
c     cufs and cufT are the Uf velocities for s and T
      i=NI
      do k=1,NK
         do j=1,NJ
            dJtinv=Jac(i,j,k)/dtime
            du = (s(i,j,k,n) -s(i-1,j,k,n))
            if (dabs(du).gt.1.d-16) then
               cufs(j,k) = -(s(i,j,k,n)-sbold(j,k))*dJtinv /du
            else
               cufs(j,k) = 0.d0
            endif

            do it=1,ntr
               du = (T(it,i,j,k,n) -T(it,i-1,j,k,n))
               if (dabs(du).gt.1.d-16) then
                  cufT(it,j,k)=-(T(it,i,j,k,n)-Tbold(it,j,k))*dJtinv/du
               else
                  cufT(it,j,k) = 0.d0
               endif
            end do
         end do
      end do
c     ---------------------------------------------------

cc     RESTORE s for 1000 time steps with e-folding time 100 time steps
c      if (step.le.1000) then
c         do k=1,NK
c            do j=1,NJ
c               do i=1,NI
c                  s(i,j,k,n)= s(i,j,k,n)
c     &                 - 0.01*(s(i,j,k,m)-s_restore(j,k))
c               end do
c            end do
c         end do
c      end if
c
c+      call remineralize(n) ( called from momentum)
c     Surface flux
c=      call surfaceflux(dtime,n)
c
c     Boundary conditions
c     call sTbc
c==      call sTbc_periodicew(n)
c
      return
      end


      subroutine ultim(up,ct,dn,gdn,gup,courinv,fc)
c     ----------------------------------------------
c     Implementation of the ULTIMATE limiter based on upstrm,dnstrm,ctr
c     points and d(phi)/dx at the up- and dnstrm locations
c     Returns the "ulitimate limited" face value of the variable
c     that should be multiplied by uf,vf or wf to give the flux.

      double precision up,ct,dn,gup,gdn,courinv,fc,ref,del,adel,acurv

      del= dn -up
      adel= dabs(del)
      acurv= dabs(gdn -gup)
      if (acurv.ge.adel)  then
         fc= ct
         return
      else
         ref= up + (ct -up)*courinv
         if (del.gt.0) then
            fc= dmax1(ct,fc)
            ref= dmin1(ref,dn)
            fc= dmin1(fc,ref)
         else
            fc= dmin1(ct,fc)
            ref= dmax1(ref,dn)
            fc= dmax1(ref,fc)
         endif
      end if
      return
      end
