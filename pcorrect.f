      subroutine pcorrect(pcorr)
c----------------------------------------------------
c     Correct the NH presure (add the correction pcorr to p)
c     -------------------------- 
c     
      implicit logical (a-z)
      include 'header.f'
      integer i,j,k
      double precision pcorr(0:NI+1,0:NJ+1,0:NK+1)

      do k=0,NK+1
         do j=0,NJ+1
            do i=0,NI+1
               p(i,j,k) = p(i,j,k) + pcorr(i,j,k)
            end do
         end do
      end do

      return
      end
