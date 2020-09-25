      subroutine momentum(pcorr,step)
c----------------------------------------------------
      implicit logical (a-z)
      include 'header.f'
      integer step,i,j,k,icount
      double precision dtim,fdiv,cfcdiv,ctrdiv,edt,pcorr(maxout)
      real dtime, tarray(2),tim
      double precision ufav
c
c     --------------------------------------------------------
c     initialized with u(t=0) and uf(t=0), ie. u(n=0),uf(n=0)
c     ------------------------------------------------------
c     get ready for the 1st step
c     evalrp and srcface were callled at the finish or the previous step
c     or in init
c
c     save old value of h
      do 10 j=0,NJ+1
         do 20 i=0,NI+1
            oldh(i,j)= h(i,j)
c            write(50,*) oldh(i,j)
 20      continue
 10   continue


C     Related to the outflow bc issue

c     save old value of uf at boundary
C      do k=1,NK
C         do j=1,NJ
C            ufbold(j,k)= uf(NI,j,k)
C         end do
C      end do
C     Added Mar 27
C      do k=1,NK
C         do j=1,NJ
C            ufav= 0.d0
C            icount= 0
C            do i=NI/4,9*NI/10
C               ufav= ufav + uf(i,j,k)
C               icount= icount +1
C            end do
C            ufbce(j,k) = ufav/dble(icount)
C         end do
C      end do

C     -----------------
c     1st R-K step
c     -----------------
      dtim= dtf/3.d0 ! (Phi^* = Phi_n + dt/3*F_n)

c     compute new zc,zf and wz,wx,wy and Jac (since grid has moved)
c=      call grassbend(0)
cc      call movegrass(dtim,0,1)
c      write(6,*) 'in mom, call grassshape'
c      if (step.eq.2) then
c         call grassshape2(0)
c      else
      call grassshape(0)
c      end if
c      write(6,*) 'in mom, finish grassshape'
c      do k=1,NKg
c         write(6,*) xg(NI/2,NJ/2,k,0),zg(NI/2,NJ/2,k,0),theta(1,1,k)
c      end do
c      stop
      call findzall
      call sigma
c     advecn(m,n,dtim) means use level m and write s,T to level n
      call advecn(0,1,dtim,step) ! Solves ustar equation (not until u1)
c=      call tracersource(1,dtim)
c=      call ecol(0,1,dtim)
c     call intpol with the time level at which uf is known
c     besides interpolating the convective terms, intpol computes u*
      ! Interpolate onto faces
      call intpol
c      call outcdf(xc,yc,zc,p,h,s,T,cx,cy,cz,vor,cxf,cyf,czf,
c     &     theta,xg(0,0,1,1),zg(0,0,1,1),ght,gdef,1)
c      stop
c     we call the following 3 after advecn, in which case
c     evalrp would use the new value of s,T(n+1), but old value of h(n)
c=      call rpevalgrad(1)
c     evalrp must be called before srcface or coriolis
c=      call coriolis(0)
      call srcface(0,step) !Source terms at the faces
      call openbc(0,1,dtim)
c      call openbc(0,1,dtim,1)   !openbc_prev
c
c      open(unit=70,file='hbef.out')
c      do j=0,NJ+1
c         write(70,*) h(NI/2,j)
c      end do
c      close(70)
      call hsolve(h,oldh,dtim) !compute h^(n+1)
c=      call calcskfc
      call vhydro(dtim) ! Calculate Uf tilde, Vf tilde
      call cfdiv(cfcdiv) ! run diagnostic
c=      call openbc2(0,1,dtim)
c      open(unit=70,file='haft.out')
c      do j=0,NJ+1
c         write(70,*) h(NI/2,j)
c      end do
c      close(70)
c      stop
c      call outcdf(xc,yc,zc,p,h,s,T,cx,cy,cz,vor,cxf,cyf,czf,
c     &     theta,xg(0,0,1,1),zg(0,0,1,1),ght,gdef,1)
c      stop
c     v,w, present beacuse cannot bal both v and vf in geostrophic perfectly
c     in calling newcor, n should be the same as in coriolis
cc      call evalrp(1)
c=      call newcor(dtim,0)
c=      call newsrc
c     compute q (correction to NH press p)
c     cpuf must be called before mgrid which calls vface and overwrites uf
cc      call exitp(p)
      edt= EPS/dtim
c     Calculate dynamic pressure
      call mgrid(pcorr,dtim,edt,cfcdiv) ! Solves for deltaP (non-hydrostatic part)
cc      call linesolve(p,dtim,step,edt)
c      call psolve(p,dtim,step,edt)
cc      call pnsolve(p,dtim,step,edt)
c      call outcdf(xc,yc,zc,pcorr,h,s(0,0,0,1),T,cx,cy,cz,vor,cxf,
c     &     cyf,czf,theta,xg(0,0,1,1),zg(0,0,1,1),ght,gdef,1)
c      stop
      call vface(pcorr,dtim) ! Compute Uf^(n+1) and so on
cc      call vnface(p,dtim)
c     fill f(n=1), use f(n=0)
c     call vcenter with the value of n that we wish to fill
      call vcenter(pcorr,dtim,1) ! Compute Uf^(n+1) and so on
c     computes the final vel (at n+1 -th time step) at the cell centers
c     Now correct the NH pressure
      call pcorrect(pcorr)
c     call openbc(m,n,nnew,dtim,isubstep)
c=      call openbc(0,1,1,dtim,1)  openbc_atend (didn't work well)
      call facediv(dtim,fdiv)
      call cdiv(dtim,ctrdiv,1) ! small but not guaranteed to go to zero

      write(6,'(3(A,ES13.6))')
     &      'rk-1:  cfcdiv ',cfcdiv,' fdiv ',fdiv,' cdiv ',ctrdiv
C       write(6,*)  'rk-1 cfcdiv ',cfcdiv,' fdiv ',fdiv,' cdiv ',ctrdiv
c      call outcdf(xc,yc,zc,pcorr,h,s(0,0,0,1),T,cx,cy,cz,vor,cxf,
c     &     cyf,czf,theta,xg(0,0,1,1),zg(0,0,1,1),ght,gdef,1)
c      stop
c      call outcdf(xc,yc,zc,pcorr,h,s(0,0,0,1),T,u(0,0,0,1),v(0,0,0,1),
c     &     w(0,0,0,1),vor,uf,vf,wf,theta,xg(0,0,1,1),zg(0,0,1,1),ght,
c     &     gdef,1)
c      write(6,*) 'stopping in momentum.f'
c      stop

c     -----------------
c     2nd R-K step
c     -----------------
      dtim= 0.5d0*dtf  ! (Phi^** = Phi_n + dt/2*F^*)
c     compute new zc,zf and wz,wx,wy and Jac (since grid has moved)
c=      call grassbend(1)
c=       call movegrass(dtim,1,1)
c      call outcdf(xc,yc,zc,pcorr,h,s(0,0,0,1),T,u(0,0,0,1),v(0,0,0,1),
c     &     w(0,0,0,1),vor,uf,vf,wf,theta,xg(0,0,1,1),zg(0,0,1,1),ght,
c     &     gdef,1)
c      stop
      call findzall
      call sigma
c      call outflowbc(1,dtim)
c     advecn(m,n,dtim) means use level m and write s,T to level n
      call advecn(1,1,dtim,step)
c=      call tracersource(1,dtim)
c=      call ecol(1,1,dtim)
cc      call correctbc
c     call intpol with the time level at which uf is known
c     besides interpolating the convective terms, intpol computes u*
      call intpol
c      call outcdf(xc,yc,zc,p,h,s,T,cx,cy,cz,vor,cxf,cyf,czf,
c     &     theta,xg(0,0,1,1),zg(0,0,1,1),ght,gdef,1)
c      stop
c     we call the following 3 after advecn, in which case
c     evalrp would use the new value of s,T(n+1), but old value of h(n)
c      tim= dtime(tarray)
c=      call rpevalgrad(1)
c      call rpflatbtm(1)
c      tim= dtime(tarray)
c      write(6,*) 'time for rpevalgrad',tim
cc      call evalrp(1)
c     evalrp must be called before srcface or coriolis
c      call viscous(1)
cc      call gradrp
c-      call vertmix(1)
c=      call coriolis(1)
      call srcface(1,step)
      call openbc(1,1,dtim)
c      call openbc(1,1,dtim,2)  !openbc_prev
ccc      call bcmake(dtim)
      call hsolve(h,oldh,dtim)
c=      call calcskfc
cc      call hnsolve(h,oldh,hdt,dtim)
      call vhydro(dtim)
c=      call openbc2(1,1,dtim)
c=      call cfdiv(cfcdiv)
c      call outcdf(xc,yc,zc,pcorr,h,s(0,0,0,1),T,cx,cy,cz,vor,cxf,cyf,czf,
c     &     theta,xg(0,0,1,1),zg(0,0,1,1),ght,
c     &     gdef,1)
c      stop
c      call ufcorrector
c     in calling newcor, n should be the same as in coriolis.
cc      call evalrp(1)
c=      call newcor(dtim,1)
c=      call newsrc
cc      call cpsi
c     compute q
c      call cpuf
cc      call exitp(p)
      edt= EPS/dtim
      call mgrid(pcorr,dtim,edt,cfcdiv)
cc      call linesolve(p,dtim,step,edt)
c-      call psolve(p,dtim,step,edt)
cc      call pnsolve(p,dtim,step,edt)
c     ********************
      call vface(pcorr,dtim)
cc      call vnface(p,dtim)
c     computes the vel at the cell faces
c     fill f(n=1), use f(n=0)
c     call vcenter with the value of n that we wish to fill
      call vcenter(pcorr,dtim,1)
c     computes the final vel (at n+1 -th time step) at the cell centers
      call pcorrect(pcorr)
c=      call openbc(1,1,0,dtim,2)  !openbc_atend
      call facediv(dtim,fdiv)
c      if (step.eq.27) write(6,*) 'rk-2, step=27 to call cdiv,fdiv',fdiv
      call cdiv(dtim,ctrdiv,1)

      write(6,'(3(A,ES13.6))')
     &      'rk-2:  cfcdiv ',cfcdiv,' fdiv ',fdiv,' cdiv ',ctrdiv
C       write(6,*)  'rk-2 cfcdiv ',cfcdiv,' fdiv ',fdiv,' cdiv ',ctrdiv
c      write(6,*)  'rk-2  fdiv  ',fdiv,' cdiv  ',ctrdiv
c
c     3rd R-K step
c     -----------------
      dtim= dtf   ! (Phi_{n+1} = Phi_n + dt*F^**)

c     compute new zc,zf and wz,wx,wy and Jac (since grid has moved)
c=      call grassbend(1)
c=      call movegrass(dtim,1,0)
      call findzall
      call sigma
c      call outflowbc(1,dtim)
c     advecn(m,n,dtim) means use level m and write s,T to level n
      call advecn(1,0,dtim,step)
c=      call tracersource(0,dtim)
c=      call ecol(1,0,dtim)
cc      call correctbc
c     call intpol with the time level at which uf is known
c     besides interpolating the convective terms, intpol computes u*
      call intpol
c=      call rpevalgrad(0)
c      call rpflatbtm(0)
cc      call evalrp(0)
c     evalrp must be called before srcface or coriolis
c      call viscous(1)
cc      call gradrp
c-      call vertmix(1)
c=      call coriolis(1)
      call srcface(1,step)
      call openbc(1,0,dtim)
c      call openbc(1,0,dtim,3)  !openbc_prev
ccc      call bcmake(dtim)
      call hsolve(h,oldh,dtim)
c=      call calcskfc
cc      call hnsolve(h,oldh,hdt,dtim)
      call vhydro(dtim)
c=      call openbc2(1,0,dtim)
c=      call cfdiv(cfcdiv)
c      call ufcorrector
c     in calling newcor, n should be the same as in coriolis.
cc      call evalrp(0)
c=      call newcor(dtim,1)
c=      call newsrc
cc      call cpsi
c     compute q
c      call cpuf
cc      call exitp(p)
      edt= EPS/dtim
      call mgrid(pcorr,dtim,edt,cfcdiv)
cc      call linesolve(p,dtim,step,edt)
c-      call psolve(p,dtim,step,edt)
cc      call pnsolve(p,dtim,step,edt)
c     *********************
      call vface(pcorr,dtim)
cc      call vnface(p,dtim)
c     computes the vel at the cell faces
c     fill f(n=1), use f(n=0)
c     call vcenter with the value of n that we wish to fill
      call vcenter(pcorr,dtim,0)
c     computes the final vel (at n+1 -th time step) at the cell centers
      call pcorrect(pcorr)
c=      call openbc(1,0,1,dtim,3)  ! openbc_atend didn't work well
      call facediv(dtim,fdiv)
      call cdiv(dtim,ctrdiv,0)

      write(6,'(3(A,ES13.6))')
     &      'rk-3:  cfcdiv ',cfcdiv,' fdiv ',fdiv,' cdiv ',ctrdiv
C       write(6,*)  'rk-3 cfcdiv ',cfcdiv,' fdiv ',fdiv,' cdiv ',ctrdiv   
c      write(6,*)  'rk-3  fdiv  ',fdiv,' cdiv  ',ctrdiv
c
c     Convective adjustment carried out every 9 time steps.
c-      if (mod(step,9).eq.0) then
c-         call conadjust(0)
c-      end if
c
c=      call resettr(step)
c
cc      call openbc(1,0,dtim)

c      open(unit=79,file='vtimeser.out',access='append')
c      open(unit=78,file='utimeser.out',access='append')
c      open(unit=76,file='htimeser.out',access='append')
c      open(unit=77,file='wtimeser.out',access='append')
c      i=NI/2
c      k=NK
c      write(76,*) h(i,1),h(i,NJ/2),h(i,NJ)
c      write(77,*) w(i,1,k,0),w(i,NJ/2,k,0),w(i,NJ,k,0)
c      write(78,*) u(i,1,k,0),u(i,NJ/2,k,0),u(i,NJ,k,0)
c      write(79,*) v(i,1,k,0),v(i,NJ/2,k,0),v(i,NJ,k,0)
c      close(78 79)
c      close(77 76)
c
cc     compute and write out extreme values
cc=      call extremes(step)
c
      return
      end
