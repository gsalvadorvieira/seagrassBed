      subroutine tracerinit
c     --------------------
c     initializes tracer fields
c
      include 'header.f'
      integer  i,j,k,l
c
      do j=0,NJ+1
         do i=0,NI+1

C           Step function
C            do k=0,NKg
C               T(1,i,j,k,0)= 1.d0
C            end do
C            do k=NKg+1,NK+1
C               T(1,i,j,k,0)= 0.d0
C            end do

C           Linear profile
            do k=0,NK+1
               T(1,i,j,k,0)= 1.d0 - 1.d0*k/(NK+1)
            end do

            if ((i.eq.1).and.(j.eq.1)) then
               do k=0,NK+1
                  tinfl(k)= T(1,i,j,k,0)
               end do
            end if

         end do
      end do
      return

      end
