      subroutine hfill(dtime,hf)
c     -----------------------------------
c     call hfill(dtime,h)
c     computes h at points outside the boundaries
c     by using central difference to specify the h gradient
c     at the boundary which is known from the condition v_normal=0.
c     We use uf,vf,wf at n=0  as the BC on uf,vf,wf.  This is ok for
c     the case where the BCs are not changing witth time.
c     We make the approx that tilde-u_nonrmal = uf_normal at the boundary.
c     We should ideally use one-sided differencing near the edges and corners
c     so as to not use these values which we are filling just by averaging.
c     Then we should also modify cpfine so as to not access these values.
c     Right now we will not bother to do this -  but will call this routine
c     a few times (alternating with hsolve) to iterate when g12 is non-zero.
c     If the BCs are changing with time we must modify this subroutine.
      implicit logical (a-z)
      include 'header.f'
      integer i,j,k,l,step,rk
      double precision dtime,edt,sumcxf,sumuf,sumsif,sumcyf,sumvf,
     &     sumsjf,hxi,heta,ginv,hf(0:NI+1,0:NJ+1),cons,dnk,gkn,gkninv,
     &     gradh,sumhxn,sumhyn,sumgi(2),sumgj(2)
c

      edt= EPS/dtime
      dnk= dfloat(NK)
c      cons= (1.d0 -kappah)*dnk*gpr
c      gkn= gpr*kappah*dnk
      cons= (1.d0 -kappah)*gpr
      gkn= gpr*kappah
      gkninv= 1.d0/gkn
c      
c     faces i=0,i=NI
c     --------------
      do 10 j=1,NJ
c         ginv= 1.d0/gi(0,j,1)
c         if (j.eq.1) then
c            heta= 0.5d0*(hf(0,j+1)+hf(1,j+1) -hf(0,j)-hf(1,j))
c         else if (j.eq.NJ) then
c            heta= 0.5d0*(hf(0,j)+hf(1,j) -hf(0,j-1)-hf(1,j-1))
c         else
            heta= 0.25*(hf(0,j+1)+hf(1,j+1) -hf(0,j-1)-hf(1,j-1))
c         endif
         sumcxf= 0.d0
         sumuf= 0.d0
         sumsif= 0.d0
         sumhxn= 0.d0
         do l=1,2
            sumgi(l)= 0.d0
         end do
         do 15 k=1,NK
            sumcxf= sumcxf +cxf(0,j,k)
            sumuf= sumuf +ufbcw(j,k)
            sumsif= sumsif +sifc(0,j,k)
            sumhxn= sumhxn +hxn(0,j,k)
            do l=1,2
               sumgi(l)= sumgi(l) +gi(0,j,k,l)
            end do
 15      continue
         gradh= (edt*(sumcxf -sumuf) -sumsif -cons*sumhxn)*gkninv
         hf(0,j)= hf(1,j) + (sumgi(2)*heta -gradh)/sumgi(1)
 10   continue
c
      do 20 j=1,NJ
c         ginv= 1.d0/gi(NI,j,1)
c         if (j.eq.1) then
c            heta= 0.5d0*(hf(NI,j+1)+hf(NI+1,j+1) -hf(NI,j)-hf(NI+1,j))
c         else if (j.eq.NJ) then
c            heta= 0.5d0*(hf(NI,j)+hf(NI+1,j) -hf(NI,j-1)-hf(NI+1,j-1))
c         else
            heta=0.25*(hf(NI,j+1)+hf(NI+1,j+1)-hf(NI,j-1)-hf(NI+1,j-1))
c         end if
         sumcxf= 0.d0
         sumuf= 0.d0
         sumsif= 0.d0
         sumhxn= 0.d0
         do l=1,2
            sumgi(l)= 0.d0
         end do
         do 25 k=1,NK
            sumcxf= sumcxf +cxf(NI,j,k)
            sumuf= sumuf +ufbce(j,k)
            sumsif= sumsif +sifc(NI,j,k)
            sumhxn= sumhxn +hxn(NI,j,k)
            do l=1,2
               sumgi(l)= sumgi(l) +gi(NI,j,k,l)
            end do
 25      continue
         gradh= (edt*(sumcxf -sumuf) -sumsif -cons*sumhxn)*gkninv
         hf(NI+1,j)= hf(NI,j) + (-sumgi(2)*heta +gradh)/sumgi(1)
 20   continue
c
c     faces j=0,j=NJ
c     --------------
      do 30 i=1,NI
c         ginv= 1.d0/gj(i,0,2)
c         if (i.eq.1) then
c            hxi= 0.5d0*(hf(i+1,0)+hf(i+1,1) -hf(i,0)-hf(i,1))
c         else if (i.eq.NI) then
c            hxi= 0.5d0*(hf(i,0)+hf(i,1) -hf(i-1,0)-hf(i-1,1))
c         else
            hxi= 0.25*(hf(i+1,0)+hf(i+1,1) -hf(i-1,0)-hf(i-1,1))
c         end if
         sumcyf= 0.d0
         sumvf= 0.d0
         sumsjf= 0.d0
         sumhyn= 0.d0
         do l=1,2
            sumgj(l)= 0.d0
         end do
         do 35 k=1,NK
            sumcyf= sumcyf +cyf(i,0,k)
            sumvf= sumvf +vfbcs(i,k)
            sumsjf= sumsjf +sjfc(i,0,k)
            sumhyn= sumhyn +hyn(i,0,k)
            do l=1,2
               sumgj(l)= sumgj(l) + gj(i,0,k,l)
            end do
 35      continue
         gradh= (edt*(sumcyf -sumvf) -sumsjf -cons*sumhyn)*gkninv
         hf(i,0)= hf(i,1) + (sumgj(1)*hxi -gradh)/sumgj(2)
 30   continue
c
      do 40 i=1,NI
c         ginv= 1.d0/gj(i,NJ,2)
c         if (i.eq.1) then
c            hxi= 0.5d0*(hf(i+1,NJ)+hf(i+1,NJ+1) -hf(i,NJ)-hf(i,NJ+1))
c         else if (i.eq.NI) then
c            hxi= 0.5d0*(hf(i,NJ)+hf(i,NJ+1) -hf(i-1,NJ)-hf(i-1,NJ+1))
c         else 
            hxi= 0.25*(hf(i+1,NJ)+hf(i+1,NJ+1)-hf(i-1,NJ)-hf(i-1,NJ+1))
c         end if
         sumcyf= 0.d0
         sumvf= 0.d0
         sumsjf= 0.d0
         sumhyn= 0.d0
         do l=1,2
            sumgj(l)= 0.d0
         end do
         do 45 k=1,NK
            sumcyf= sumcyf +cyf(i,NJ,k)
            sumvf= sumvf +vfbcn(i,k)
            sumsjf= sumsjf +sjfc(i,NJ,k)
            sumhyn= sumhyn +hyn(i,NJ,k)
            do l=1,2
               sumgj(l)= sumgj(l) +gj(i,NJ,k,l)
            end do
 45      continue
         gradh= (edt*(sumcyf -sumvf) -sumsjf -cons*sumhyn)*gkninv
         hf(i,NJ+1)= hf(i,NJ) + (-sumgj(1)*hxi +gradh)/sumgj(2)
 40   continue
c
c     4 corners
      hf(0,0)= 0.5d0*(hf(1,0)+hf(0,1))
      hf(NI+1,0)= 0.5d0*(hf(NI,0)+hf(NI+1,1))
      hf(0,NJ+1)= 0.5d0*(hf(0,NJ)+hf(1,NJ+1))
      hf(NI+1,NJ+1)= 0.5d0*(hf(NI,NJ+1)+hf(NI+1,NJ))
c
      return
      end
      

