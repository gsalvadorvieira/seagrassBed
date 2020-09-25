      subroutine srcface(n,step)
c----------------------------------------------------
c     FOR PERIODICEW boundaries
c     --------------------------
      implicit double precision (a-z)
      include 'header.f'
      integer i,j,k,n,in,step,iter

      double precision drpx(NI,NJ,NK),drpy(NI,NJ,NK),
     &     grpifc(0:NI,NJ,NK),grpjfc(NI,0:NJ,NK)
      common/rpgrads/drpx,drpy,grpifc,grpjfc

c     si,sj are already set to drpx,drpy in steadystate and not
c     changed with time. Since si,sj are not used in hsolve, the
c     px and py terms can be added later in uvchy.

c     for openbc.  mgpfill should have been called

      do k=1,NK
            do j=1,NJ
               do i=0,NI
                  px= (p(i+1,j,k) -p(i,j,k))*gqi(i,j,k,1) +0.25*
     &              (p(i+1,j+1,k)+p(i,j+1,k)-p(i+1,j-1,k)-p(i,j-1,k))
     &              *gqi(i,j,k,2)+0.25*(p(i+1,j,k+1)+p(i,j,k+1)-
     &              p(i+1,j,k-1)-p(i,j,k-1))*gqi3(i,j,k)
                  sifc(i,j,k)= grpifc(i,j,k)  + px
            end do
         end do
      end do
c
c========= no y gradient
c      return
c=============

c
      do k=1,NK
         do i=1,NI
            do j=0,NJ
               py= (p(i,j+1,k)
     &              -p(i,j,k))*gqj(i,j,k,2) +0.25*(p(i+1,j+1,k)
     &              +p(i+1,j,k)-p(i-1,j+1,k)-p(i-1,j,k))*gqj(i,j,k,1)
     &              +0.25*(p(i,j+1,k+1)+p(i,j,k+1)-p(i,j+1,k-1)
     &              -p(i,j,k-1))*gqj3(i,j,k)
               sjfc(i,j,k)= grpjfc(i,j,k) + py
            end do
c     mgpfill should have been called before
         end do
      end do
c
      return

      end
