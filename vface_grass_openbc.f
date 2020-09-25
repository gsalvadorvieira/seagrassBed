      subroutine vface(pf,dtime)
c----------------------------------------------------
c     openbc
c     -------------------------- 
c     compute the final vel (n+1 th step) at the cell faces
c     applies the boundary conditions
      implicit logical (a-z)
      include 'header.f'
      integer i,j,k,nexit,iflag
      double precision dte,px,py,pz,dtime,
     &     pf(0:NI+1,0:NJ+1,0:NK+1)
c
      dte= dtime/EPS
c      kaph1= 1.d0 -kappah
c
      do 10 j=1,NJ
         do 15 i=1,NI
            do 16 k=1,NK
               px= (pf(i+1,j,k) -pf(i,j,k))*gqi(i,j,k,1) +0.25*
     &           (pf(i+1,j+1,k)+pf(i,j+1,k)
     &              -pf(i+1,j-1,k)-pf(i,j-1,k))
     &            *gqi(i,j,k,2)+0.25*(pf(i+1,j,k+1)+pf(i,j,k+1)-
     &            pf(i+1,j,k-1)-pf(i,j,k-1))*gqi3(i,j,k)
               uf(i,j,k)= cxf(i,j,k) -dte*( px )

 16         continue
 15      continue
 10   continue

c     ufbc found from previous call to openbc
c     -----------------
c     Done in openbc

c      do k=1,NK
c         do j=1,NJ
c            uf(0,j,k)=ufbcw(j,k)
c            uf(NI,j,k)=ufbce(j,k)
c         end do
c      end do

c     y-direction
c     -----------
      do 20 i=1,NI
         do 25 j=0,NJ
            do 26 k=1,NK
               if (j.eq.0) then
                  vf(i,j,k)= vfbcs(i,k)
               else if (j.eq.NJ) then
                  vf(i,j,k)= vfbcn(i,k)
               else
               py= (pf(i,j+1,k) 
     &                 -pf(i,j,k))*gqj(i,j,k,2) +0.25*(pf(i+1,j+1,k)
     &          +pf(i+1,j,k)-pf(i-1,j+1,k)-pf(i-1,j,k))*gqj(i,j,k,1)
     &          +0.25*(pf(i,j+1,k+1)+pf(i,j,k+1)-pf(i,j+1,k-1)
     &          -pf(i,j,k-1))*gqj3(i,j,k)
               vf(i,j,k)= cyf(i,j,k) -dte*( py )
               endif
 26         continue
 25      continue
 20   continue
c         
      do 30 j=1,NJ
         do 30 i=1,NI
            wf(i,j,0)= wfbcb(i,j)
c     wf(i,j,NK)= wfbct(i,j)
c     do 35 k=1,NK-1
c+    const= dtime/(J2d(i,j)*HDL)
            do 35 k=1,NK
               pz= (pf(i,j,k+1) -pf(i,j,k))*gqk(i,j,k,3) +
     &           0.25*(pf(i+1,j,k+1)
     &           +pf(i+1,j,k)-pf(i-1,j,k+1)-pf(i-1,j,k))*gqk(i,j,k,1)
     &           +0.25*(pf(i,j+1,k+1)+pf(i,j+1,k)-pf(i,j-1,k+1)
     &           -pf(i,j-1,k))*gqk(i,j,k,2)
               wf(i,j,k)= czf(i,j,k) -dte*(pz )
 35         continue
 30   continue
c
c      j=17
c      i=15
c      do k=1,NK
c         write(6,*) uf(i,j,k),vf(i,j,k)
c      end do
c      write(6,*) 'czf'
c      do k=0,NK
c         write(6,*) czf(i,j,k),wf(i,j,k)
c      end do
c      stop
      return
      end
      

