      subroutine mgpfill(dtime,pf)
c     -----------------------------------
c     for openbc
c     --------------------------
c     In this version of mgpfill, the dpdsig term is left out
c     at the vertical faces because g13 g23 are larger than
c     g11 and g12.
c     call mgpfill(dtime,p(loco(m)))
c     computes the pressure at points outside the impermeable boundaries
c     by using central difference to specify the pressure gradient
c     at the boundary which is known from the condition v_normal=0.
c     We use uf,vf,wf at n=0  as the BC on uf,vf,wf.  This is ok for
c     the case where the BCs are not changing witth time.
c     We should ideally use one-sided differencing near the edges and corners
c     so as to not use these values which we are filling just by averaging.
c     Then we should also modify cpfine so as to not access these values.
c     Right now we will not bother to do this. 
c     If the BCs are changing with time we must modify this subroutine.
c     Remove sifc  in the version with grass
      implicit logical (a-z)
      include 'header.f'
      integer i,j,k
      double precision dtime,edt,dpdxi,dpdeta,dpdsig,pxi,peta,psig,
     &     pf(0:NI+1,0:NJ+1,0:NK+1),gin
c
      if (dtime.gt.1.d-16) then
         edt= EPS/dtime
      else
c         write(6,*) 'dtime 0, set edt to 0 in mgpfill'
         edt = 0.0
      end if
      

c      delinv= 1/delta, delta= D/L
c      
      if (qpr.ne.0) then

cc     face k=0
cc     --------
c      k=0
c      do 50 i=1,NI
c         do  50 j=1,NJ
c            gin= 1.d0/gqk(i,j,k,3)
cc     note that quantities are not defined at i=0,NI+1,j=0,NJ+1, and
cc     so we revert to filling pf inside the the domain only.
c            dpdxi= 0.25*( pf(i+1,j,k+1) +pf(i+1,j,k)
c     &           -pf(i-1,j,k+1) -pf(i-1,j,k) )
c            dpdeta= 0.25*( pf(i,j+1,k+1) +pf(i,j+1,k)
c     &           -pf(i,j-1,k+1) -pf(i,j-1,k) )
c            psig= gin*(-skfc(i,j,k)
c     &           - gqk(i,j,k,1)*dpdxi - gqk(i,j,k,2)*dpdeta +
c     &           edt*(czf(i,j,k)-wfbcb(i,j)) )
cc     since  : wf(i,j,k)= czf(i,j,k) -dte*delta*(pz +skfc(i,j,k))
c            pf(i,j,0)= pf(i,j,1) - psig
c 50   continue
c
cc     bottom edges (k=0)
c      j=0
c      do i=1,NI
c         gin= 1.d0/gqk(i,j,0,3)
c         dpdxi= 0.25*( pf(i+1,j,k+1) +pf(i+1,j,k)
c     &        -pf(i-1,j,k+1) -pf(i-1,j,k) )
c         dpdeta= 0.5*( pf(i,j+1,k+1) +pf(i,j+1,k)
c     &        -pf(i,j,k+1) -pf(i,j,k) )
c         psig= gin*(-skfc(i,j,k)
c     &        - gqk(i,j,k,1)*dpdxi - gqk(i,j,k,2)*dpdeta )
c         pf(i,j,0)= pf(i,j,1) - psig
c      end do
c      j=NJ+1
c      do i=1,NI
c         gin= 1.d0/gqk(i,j,0,3)
c         dpdxi= 0.25*( pf(i+1,j,k+1) +pf(i+1,j,k)
c     &        -pf(i-1,j,k+1) -pf(i-1,j,k) )
c         dpdeta= 0.5*( pf(i,j,k+1) +pf(i,j,k)
c     &        -pf(i,j-1,k+1) -pf(i,j-1,k) )
c         psig= gin*(-skfc(i,j,k)
c     &        - gqk(i,j,k,1)*dpdxi - gqk(i,j,k,2)*dpdeta )
c         pf(i,j,0)= pf(i,j,1) - psig
c      end do
cc
cc     face k=NK
cc     --------
c      k=NK
c      do i=1,NI
c         do  60 j=0,NJ+1
c            pf(i,j,NK+1)= -pf(i,j,NK)
cc            pf(i,j,NK+1)= 0.
c 60   continue
c
c      pf(0,0,0)= (pf(1,0,0)+pf(0,1,0)+pf(0,0,1))/3.d0
c      pf(NI+1,0,0)= (pf(NI,0,0)+pf(NI+1,1,0)+pf(NI+1,0,1))/3.d0
c      pf(NI+1,NJ+1,0)=(pf(NI,NJ+1,0)+pf(NI+1,NJ,0)+pf(NI+1,NJ+1,1))/3.d0
c      pf(0,NJ+1,0)= (pf(1,NJ+1,0)+pf(0,NJ,0)+pf(0,NJ+1,1))/3.d0

c     faces i=0,i=NI
c     --------------
      i=0
      do 10 j=1,NJ
c         if (ufbcw(j,NK).gt.0.d0) then
         do 11 k=1,NK
c-         if (ufbcw(j,k).ge.0.d0) then
            gin= 1.d0/gqi(i,j,k,1)
            dpdeta= 0.25*( pf(i,j+1,k) +pf(i+1,j+1,k)
     &           -pf(i,j-1,k) -pf(i+1,j-1,k) )
            dpdsig= 0.25*( pf(i,j,k+1) +pf(i+1,j,k+1)
     &           -pf(i,j,k-1) -pf(i+1,j,k-1) )
c     but cxf(which is u-tilde)= uf at the boundaries
            pxi= gin*(
cc     &           - gqi(i,j,2)*dpdeta - gqi3(i,j,k)*dpdsig
     &           - gqi(i,j,k,2)*dpdeta 
     &           +(edt*(cxf(i,j,k)-ufbcw(j,k)) ) )
c     &           -sifc(i,j,k)) )
            pf(0,j,k)= pf(1,j,k) - pxi
c-         end if
 11      continue
c         endif
 10   continue
c
c
      i=NI
      do 20 j=1,NJ
c         if (ufbce(j,NK).lt.0.d0) then
         do 21 k=1,NK
c-         if (ufbce(j,k).le.0.d0) then
            gin= 1.d0/gqi(i,j,k,1)
            dpdeta= 0.25*( pf(i,j+1,k) +pf(i+1,j+1,k)
     &           -pf(i,j-1,k) -pf(i+1,j-1,k) )
            dpdsig= 0.25*( pf(i,j,k+1) +pf(i+1,j,k+1)
     &           -pf(i,j,k-1) -pf(i+1,j,k-1) )
            pxi= gin*(
cc     &           - gqi(i,j,2)*dpdeta - gqi3(i,j,k)*dpdsig
     &           - gqi(i,j,k,2)*dpdeta 
     &           +(edt*(cxf(i,j,k)-ufbce(j,k)) ) )
c     &           -sifc(i,j,k)) )
            pf(NI+1,j,k)= pf(NI,j,k) + pxi
c-            end if
 21      continue
c         endif
 20   continue
c


c     faces j=0,j=NJ
c     --------------
      j=0
      do 30 i=1,NI
         do 31 k=1,NK
            gin= 1.d0/gqj(i,j,k,2)
            dpdxi= 0.25*( pf(i+1,j+1,k) +pf(i+1,j,k)
     &           -pf(i-1,j+1,k) -pf(i-1,j,k) )
            dpdsig= 0.25*( pf(i,j+1,k+1) +pf(i,j,k+1)
     &           -pf(i,j+1,k-1) -pf(i,j,k-1) )
            peta= gin*(
cc     &           - gqj(i,j,1)*dpdxi - gqj3(i,j,k)*dpdsig
     &           - gqj(i,j,k,1)*dpdxi 
     &           +(edt*(cyf(i,j,k)-vfbcs(i,k))  ) )
c     &           -sjfc(i,j,k)) )
            pf(i,0,k)= pf(i,1,k) - peta
 31      continue
 30   continue
c
      j=NJ
      do 40 i=1,NI
         do 41 k=1,NK
            gin= 1.d0/gqj(i,j,k,2)
            dpdxi= 0.25*( pf(i+1,j+1,k) +pf(i+1,j,k)
     &           -pf(i-1,j+1,k) -pf(i-1,j,k) )
            dpdsig= 0.25*( pf(i,j+1,k+1) +pf(i,j,k+1)
     &           -pf(i,j+1,k-1) -pf(i,j,k-1) )
c            peta= gin*(edt*(cyf(i,j,k)-cyf(i,j,k,0)) 
            peta= gin*(
cc     &           - gqj(i,j,1)*dpdxi - gqj3(i,j,k)*dpdsig
     &           - gqj(i,j,k,1)*dpdxi 
     &           +(edt*(cyf(i,j,k)-vfbcn(i,k)) ) )
c     &-sjfc(i,j,k)) )
            pf(i,NJ+1,k)= pf(i,NJ,k) + peta
 41      continue
 40   continue
c
c     Edges : columns
      do k=1,NK
         pf(NI+1,NJ+1,k)= 0.25*(pf(NI+1,NJ,k)+ pf(NI,NJ+1,k)) 
     &        + 0.5*pf(NI,NJ,k)
         pf(0,0,k)= 0.25*(pf(0,1,k)+ pf(1,0,k))
     &        + 0.5*pf(1,1,k)
         pf(NI+1,0,k)= 0.25*(pf(NI+1,1,k)+ pf(NI,0,k)) 
     &        + 0.5*pf(NI,1,k)
         pf(0,NJ+1,k)= 0.25*(pf(0,NJ,k)+ pf(1,NJ+1,k)) 
     &        + 0.5*pf(1,NJ,k)
      end do

c     face k=0
c     --------
      k=0
      do i=1,NI
         do  j=1,NJ
            gin= 1.d0/gqk(i,j,k,3)
c     note that quantities are not defined at i=0,NI+1,j=0,NJ+1, and
c     so we revert to filling pf inside the the domain only.
            dpdxi= 0.25*( pf(i+1,j,k+1) +pf(i+1,j,k)
     &           -pf(i-1,j,k+1) -pf(i-1,j,k) )
            dpdeta= 0.25*( pf(i,j+1,k+1) +pf(i,j+1,k)
     &           -pf(i,j-1,k+1) -pf(i,j-1,k) )
            psig= gin*(-skfc(i,j,k)
     &           - gqk(i,j,k,1)*dpdxi - gqk(i,j,k,2)*dpdeta +
     &           edt*(czf(i,j,k)-wfbcb(i,j)) )
c     since  : wf(i,j,k)= czf(i,j,k) -dte*delta*(pz +skfc(i,j,k))
            pf(i,j,0)= pf(i,j,1) - psig
         end do
      end do


c     Bottom edges (k=0)
      do i=1,NI
         pf(i,0,0)= 0.25*(pf(i,1,0)+pf(i,0,1))+0.5*pf(i,1,1)
         pf(i,NJ+1,0)= 0.25*(pf(i,NJ,0)+pf(i,NJ+1,1))+0.5*pf(i,NJ,1)
      end do
      do j=0,NJ+1
         pf(0,j,0)= 0.25*(pf(1,j,0)+pf(0,j,1))+0.5*pf(1,j,1)
         pf(NI+1,j,0)= 0.25*(pf(NI,j,0)+pf(NI+1,j,1))+0.5*pf(NI,j,1)
      end do
      goto 99 ! the following is not needed : I'm averaging instead
c     bottom edges (k=0)
      j=0
      do i=1,NI
         gin= 1.d0/gqk(i,j,0,3)
         dpdxi= 0.25*( pf(i+1,j,k+1) +pf(i+1,j,k)
     &        -pf(i-1,j,k+1) -pf(i-1,j,k) )
         dpdeta= 0.5*( pf(i,j+1,k+1) +pf(i,j+1,k)
     &        -pf(i,j,k+1) -pf(i,j,k) )
         psig= gin*(-skfc(i,j,k)
     &        - gqk(i,j,k,1)*dpdxi - gqk(i,j,k,2)*dpdeta )
         pf(i,j,0)= pf(i,j,1) - psig
      end do
      j=NJ+1
      do i=1,NI
         gin= 1.d0/gqk(i,j,0,3)
         dpdxi= 0.25*( pf(i+1,j,k+1) +pf(i+1,j,k)
     &        -pf(i-1,j,k+1) -pf(i-1,j,k) )
         dpdeta= 0.5*( pf(i,j,k+1) +pf(i,j,k)
     &        -pf(i,j-1,k+1) -pf(i,j-1,k) )
         psig= gin*(-skfc(i,j,k)
     &        - gqk(i,j,k,1)*dpdxi - gqk(i,j,k,2)*dpdeta )
         pf(i,j,0)= pf(i,j,1) - psig
      end do
c
c     face k=NK
c     --------
 99   k=NK
      do i=0,NI+1
         do j=0,NJ+1
            pf(i,j,NK+1)= -pf(i,j,NK)
c            pf(i,j,NK+1)= 0.
         end do
      end do

      else
c     face k=NK
c     --------
      k=NK
      do i=1,NI
         do j=0,NJ+1
            pf(i,j,NK+1)= -pf(i,j,NK)
         end do
      end do
      end if
c     the above is done only if qpr is non-zero
c
c
      return
      end

