      subroutine restrict(nxm,nym,nzm,fin,cor)
c     call restrict(nx(m),ny(m),nz(m),res,rhs(loci(m+1)))
c     -----------------------------------------------------------
      implicit logical (a-z)
c     interpolate a value from grid level "lev" down to the next coarser grid
c     at level "lev+1".
c     The array at lev is specified as a 3-d array
c     fin(nx(lev),ny(lev),nz(lev)). The array at lev+1 is placed in
c     the correct position in the one-dim array cor(maxint)
      integer nxm,nym,nzm,i,j,k,ii,jj,kk
      double precision fin(nxm,nym,nzm),cor(nxm/2,nym/2,nzm/2),coeff
c
      coeff= 0.125*4.d0
c     0.125 is due to the averaging. The factor 4 is here because we
c     have doubled the grid spacing and we are multiplying the rhs by 4
c     to account for this. Alternately we could have divied the coefficents
c     by 4. 
      do 10 k=1,nzm/2
         kk= 2*k-1
         do 20 j=1,nym/2
            jj= 2*j-1
            do 30 i=1,nxm/2
               ii= 2*i-1
               cor(i,j,k)= coeff*(fin(ii,jj,kk)+fin(ii,jj+1,kk)
     &              +fin(ii+1,jj,kk)+fin(ii+1,jj+1,kk) +
     &              fin(ii,jj,kk+1) + fin(ii,jj+1,kk+1)
     &              + fin(ii+1,jj,kk+1) + fin(ii+1,jj+1,kk+1) )
 30         continue
 20      continue
 10   continue
      return
      end
      
