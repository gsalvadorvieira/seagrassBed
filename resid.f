      subroutine resid(m,nxm,nym,nzm,cp,pp,fn,res,maxres)
c     --------------------------------------------------
c     call resid(m,nx(m),ny(m),nz(m),cp(loccp(m)),p(loco(m)),
c    &       rhs(loci(m)),res )
c
c     check the residual
      implicit logical (a-z)
      include 'header.f'
      integer i,j,k,nxm,nym,nzm,m, imax,jmax,kmax,n,step
      double precision pp(0:nxm+1,0:nym+1,0:nzm+1), 
     &     fn(nxm,nym,nzm),cp(19,nxm,nym,nzm),
     &     res(nxm,nym,nzm),maxres
c
c
      maxres= 0.d0
C$DOACROSSSHARE(cp,pp,fn,maxres,res),
C$&   LOCAL(i,j,k)
      do 110 k=1,nzm
         do 120 j=1,nym
            do 130 i=1,nxm
c     cpim1*pp(i-1) +cpip1*pp(i+1) +cpjm1p(j-1) ...+ fn = cp0*pp(i,j,k)
c     res = b - A x'     
               res(i,j,k)= fn(i,j,k)- 
     &              (cp(1,i,j,k)*pp(i,j,k)
     &              +cp(2,i,j,k)*pp(i+1,j,k)
     &              +cp(3,i,j,k)*pp(i,j+1,k)
     &              +cp(4,i,j,k)*pp(i-1,j,k)
     &              +cp(5,i,j,k)*pp(i,j-1,k)
     &              +cp(6,i,j,k)*pp(i,j,k+1)
     &              +cp(7,i,j,k)*pp(i,j,k-1)
     &              +cp(8,i,j,k)*pp(i-1,j+1,k)
     &              +cp(9,i,j,k)*pp(i-1,j-1,k)
     &              +cp(10,i,j,k)*pp(i+1,j-1,k)
     &              +cp(11,i,j,k)*pp(i+1,j+1,k)
     &              +cp(12,i,j,k)*pp(i-1,j,k-1)
     &              +cp(13,i,j,k)*pp(i+1,j,k-1)
     &              +cp(14,i,j,k)*pp(i+1,j,k+1)
     &              +cp(15,i,j,k)*pp(i-1,j,k+1)
     &              +cp(16,i,j,k)*pp(i,j-1,k-1)
     &              +cp(17,i,j,k)*pp(i,j+1,k-1)
     &              +cp(18,i,j,k)*pp(i,j+1,k+1)
     &              +cp(19,i,j,k)*pp(i,j-1,k+1))
               if (dabs(res(i,j,k)).gt.maxres) then
                  maxres= dabs(res(i,j,k))
                  imax = i
                  jmax = j
                  kmax = k
               end if
 130        continue
 120     continue
 110  continue
c
      if (maxres.gt.3000.d0) then
C          write(6,*) 'STOP in resid. res too large, i,j,k,maxres=',
C      &        maxres
            write(6,'(A, 3(I6), ES16.6)')
     &    'STOP in resid. res too large. imax,jmax,kmax,maxres=',
     &     imax, jmax, kmax, maxres

            step=999999
            write(6,*) 'failed at step=',step
            call outcdf(xc,yc,zc,p,h,s,T,u,v,w,vor,uf,vf,wf,theta,xg,zg,
     &          ght,gdef,step)      
           stop
      end if
      write(6,*) 'level',m,'     maxres',maxres
 202  return
      end

