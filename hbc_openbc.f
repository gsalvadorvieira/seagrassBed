      subroutine hbc(chf,fn,dtime)
c-------------------------------------------------------------------c
c     evaluates the boundaries in psolve c fills in
c     ch0,chim1,chip1,chjm1,chjp1,chkm1,chkp1,fn for the boundaries
      implicit logical (a-z)
      include 'header.f'
c     3-D Poisson solver - uses SOR
c     Written for the control volume formulation
c     Uses the 19 point stencil (in 3d).
c     
c     fn is the contribution to the source terms on the RHS.
c     Whenever possible we substitue the the velocities themselves
c     (at the boundaries). We cannot substitute the u,w vels at the entrance.
c     Our eqn is
c     chim1*p(i-1) +chip1*p(i+1) +chjm1p(j-1) ... -ch0*p(i,j,k)= f
      integer i,j,k,l
      double precision dtime,edt,chf(9,NI,NJ),
c      oldh(0:NI+1,0:NJ+1),
     &     fn(NI,NJ),dtinv,sumuf,sumvf,sumsif,sumsjf,
     &     sumcxf,sumcyf,gprinv,eg,const,kaph1,
     &     sumhxn,sumgi(2),sumhyn,sumgj(2)
c     
      kaph1= 1.d0 -kappah
      edt= EPS/dtime
      eg= edt/(gpr*kappah)
      gprinv= 1.d0/(gpr*kappah)
      dtinv= HDL/(dtime*kappah)
c      dnk= dfloat(NK)
      const= kaphinv -1.d0
c      write(6,*) 'const=',const,'kappah=',kappah
c      write(6,*) 'in hbc,edt,eg,gprinv,dtinv,dnk',
c     &     edt,eg,gprinv,dtinv,dnk
c
c         
c     For i=2...NI-1, j=2...NJ-1, k=2,NK-1  it is straightforward
      do 10 i=1,NI
         do 20 j=1,NJ
            if ((i.eq.1).or.(i.eq.NI).or.(j.eq.1).or.
     &           (j.eq.NJ)) then
               fn(i,j)= eg*( kaphinv*wfbcb(i,j)
     &              -J2d(i,j)*oldh(i,j)*dtinv )
               chf(1,i,j)= -eg*J2d(i,j)*dtinv 
               do l=2,9
                  chf(l,i,j)= 0.d0
               end do
               if (i.eq.NI) then
c     Take sum U-tilde(n+1) equal to sum uf(n+1) 
                  sumuf= 0.d0
                  do 30 k=1,NK
                     sumuf= sumuf +ufbce(j,k)
 30               continue
c                  write(6,*) 'fn',i,j,fn(i,j),sumuf,eg,kaphinv
                  fn(i,j)= fn(i,j) + eg*sumuf*kaphinv
c                  write(6,*) 'fn',i,j,eg*sumuf*kaphinv,fn(i,j)
c     Can set ufbce equal to uf in vface.
c     use uf from previous time step:  just here because of the exit
c                     sumuf= sumuf +uf(i,j,k)
c     try sum U* as bc on  sum U-tilde
c     doesn't work   sumuf= sumuf +cxf(i,j,k)
c     here cxf is uses as bc for  u-tilde
c                     sumuf= sumuf +kappah*cxf(i,j,k) +kaph1*ufbce(j,k)
               else
                  sumsif= 0.d0
                  sumcxf= 0.d0
                  sumuf= 0.d0
                  sumhxn= 0.d0
                  do l=1,2
                     sumgi(l)= 0.d0
                  end do
                  do 35 k=1,NK
                     sumsif= sumsif +sifc(i,j,k)
                     sumcxf= sumcxf +cxf(i,j,k)
                     sumuf= sumuf +uf(i,j,k)
                     sumhxn= sumhxn +hxn(i,j,k)
                     do l=1,2
                        sumgi(l)= sumgi(l) + gi(i,j,k,l)
                     end do
 35               continue
                  fn(i,j)= fn(i,j) +eg*(sumcxf +const*sumuf)
     &                 -gprinv*sumsif -const*sumhxn
                  chf(1,i,j)= chf(1,i,j) -sumgi(1)
                  chf(2,i,j)= chf(2,i,j) +sumgi(1)
                  chf(4,i,j)= chf(4,i,j) +0.25*sumgi(2)
                  chf(5,i,j)= chf(5,i,j) -0.25*sumgi(2)
                  chf(6,i,j)= chf(6,i,j) +0.25*sumgi(2)
                  chf(8,i,j)= chf(8,i,j) -0.25*sumgi(2)
               endif
               if (i.eq.1) then
                  sumuf= 0.d0
                  do 40 k=1,NK
                     sumuf= sumuf + ufbcw(j,k)
c                     sumuf= sumuf +kappah*cxf(0,j,k) +kaph1*ufbcw(j,k)
 40               continue
c                  write(6,*) 'fn',i,j,fn(i,j),sumuf,eg,kaphinv
                  fn(i,j)= fn(i,j) - eg*sumuf*kaphinv
c                  write(6,*) 'fn2',i,j,fn(i,j),sumuf
               else
                  sumsif= 0.d0
                  sumcxf= 0.d0
                  sumuf= 0.d0
                  sumhxn= 0.d0
                  do l=1,2
                     sumgi(l)= 0.d0
                  end do
                  do 45 k=1,NK
                     sumsif= sumsif +sifc(i-1,j,k)
                     sumcxf= sumcxf +cxf(i-1,j,k)
                     sumuf= sumuf +uf(i-1,j,k)
                     sumhxn= sumhxn +hxn(i-1,j,k)
                     do l=1,2
                        sumgi(l)= sumgi(l) +gi(i-1,j,k,l)
                     end do
 45               continue
                  fn(i,j)= fn(i,j) -eg*(sumcxf +const*sumuf)
     &                 + gprinv*sumsif + const*sumhxn
c                  if (i.eq.NI) write(6,*) j,eg*(sumcxf +const*sumuf),
c     &                 gprinv*sumsif + const*sumhxn, fn(i,j)
                  chf(1,i,j)= chf(1,i,j) -sumgi(1)
                  chf(3,i,j)= chf(3,i,j) +sumgi(1)
                  chf(4,i,j)= chf(4,i,j) -0.25*sumgi(2) 
                  chf(5,i,j)= chf(5,i,j) +0.25*sumgi(2) 
                  chf(7,i,j)= chf(7,i,j) -0.25*sumgi(2) 
                  chf(9,i,j)= chf(9,i,j) +0.25*sumgi(2)
               endif
c               if ((i.eq.1).or.(i.eq.NI)) write(6,*) 'fn',i,j,fn(i,j)
               if (j.eq.NJ) then
                  sumvf= 0.d0
                  do 50 k=1,NK
                     sumvf= sumvf +vfbcn(i,k)
c                     sumvf= sumvf +kappah*cyf(i,j,k) +kaph1*vfbcn(i,k)
 50               continue
                  fn(i,j)= fn(i,j) +eg*sumvf*kaphinv
               else 
                  sumsjf= 0.d0
                  sumcyf= 0.d0
                  sumvf= 0.d0
                  sumhyn= 0.d0
                  do l=1,2
                     sumgj(l)= 0.d0
                  end do
                  do 55 k=1,NK
                     sumsjf= sumsjf +sjfc(i,j,k)
                     sumcyf= sumcyf +cyf(i,j,k)
                     sumvf= sumvf +vf(i,j,k)
                     sumhyn= sumhyn +hyn(i,j,k)
                     do l=1,2
                        sumgj(l)= sumgj(l) +gj(i,j,k,l)
                     end do
 55               continue
                  fn(i,j)= fn(i,j) +eg*(sumcyf +const*sumvf)
     &                 - gprinv*sumsjf - const*sumhyn
                  chf(1,i,j)= chf(1,i,j) -sumgj(2)
                  chf(2,i,j)= chf(2,i,j) +0.25*sumgj(1)
                  chf(3,i,j)= chf(3,i,j) -0.25*sumgj(1)
                  chf(4,i,j)= chf(4,i,j) +sumgj(2)
                  chf(6,i,j)= chf(6,i,j) +0.25*sumgj(1)
                  chf(7,i,j)= chf(7,i,j) -0.25*sumgj(1)
               endif
               if (j.eq.1) then
                  sumvf= 0.d0
                  do 60 k=1,NK
                     sumvf= sumvf +vfbcs(i,k)
c                     sumvf= sumvf +kappah*cyf(i,0,k) +kaph1*vfbcs(i,k)
 60               continue
                  fn(i,j)= fn(i,j) -eg*sumvf*kaphinv
               else
                  sumsjf= 0.d0
                  sumcyf= 0.d0
                  sumvf= 0.d0
                  sumhyn= 0.d0
                  do l=1,2
                     sumgj(l)= 0.d0
                  end do
                  do 65 k=1,NK
                     sumsjf= sumsjf +sjfc(i,j-1,k)
                     sumcyf= sumcyf +cyf(i,j-1,k)
                     sumvf= sumvf +vf(i,j-1,k)
                     sumhyn= sumhyn +hyn(i,j-1,k)
                     do l=1,2
                        sumgj(l)= sumgj(l) +gj(i,j-1,k,l)
                     end do
 65               continue
                  fn(i,j)= fn(i,j) -eg*(sumcyf +const*sumvf)
     &                 + gprinv*sumsjf + const*sumhyn
                  chf(1,i,j)= chf(1,i,j) -sumgj(2)
                  chf(2,i,j)= chf(2,i,j) -0.25*sumgj(1)
                  chf(3,i,j)= chf(3,i,j) +0.25*sumgj(1)
                  chf(5,i,j)= chf(5,i,j) +sumgj(2)
                  chf(8,i,j)= chf(8,i,j) -0.25*sumgj(1)
                  chf(9,i,j)= chf(9,i,j) +0.25*sumgj(1)
               endif
            endif
 20      continue
 10   continue
c     
      return
      end
      
