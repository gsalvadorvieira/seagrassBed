      subroutine openbc(m,n,dtime)
c     ------------------------------------------------
c     use level m, write to level n
c     sets the face velocities at outflow and the
c     variables at NI+1.
c     sets the face velocities at inflow and the
c     variables at i=0.
c
      implicit logical (a-z)
      include 'header.f'
      integer i,j,k,m,n,it,counter
C      added this declaration for dtime and dtJ
      double precision dtime, dtJ
C      double precision ufmean

C     Compute mean velocity to impose at outlet
C       counter = 0
C       ufmean = 0.d0
C       do k=1,NK
C          do j=1,NJ
C            ufmean = ufmean + uf(1,j,k)
C            counter = counter + 1
C          end do
C       end do
C       ufmean = ufmean/counter

      do k=1,NK
         do j=1,NJ
C             ufbce(j,k)= uf(NI/2,j,k)
C            ufbce(j,k)= uf(1,j,k)
C            ufbce(j,k)= ufmean
C            ufbce(j,k)= uf(NI*15/16,j,k)
           ufbce(j,k)= uf(NI-9,j,k)
C            ufbce(j,k)= 2.d0*uf(NI-1,j,k)-uf(NI-2,j,k)
         end do
      end do
      
      do k=1,NK
         do j=1,NJ
C            s(NI+1,j,k,n)= s(NI+1,j,k,0) - merge(1.d0, 0.d0, dtJ)*
C     &           u(NI,j,k,n)*(s(NI+1,j,k,m)

C           I am unsure if it should be NI or NI+1 here, but likely NI from previous studies
C             dtJ = dtime/Jac(NI,j,k)
C             dtJ = dtime/Jac(NI+1,j,k)

            s(NI+1,j,k,n)= s(NI,j,k,n)
            do it=1,ntr
               T(it,NI+1,j,k,n)= T(it,NI,j,k,n)
            end do
         end do
      end do

      return
      end
