      subroutine mprove(h,ch,oldrhs,dtime)
c     --------------
      implicit logical (a-z)
      include 'dims.f'
      integer iter,i,j,l
      double precision ch(9,NI,NJ),rhs(NI,NJ),oldrhs(NI,NJ),
     &     he(0:NI+1,0:NJ+1),h(0:NI+1,0:NJ+1),res,
     &     dtime,tol,maxres,hstar,rlx

      tol= 1.d-20
      do j=0,NJ+1
         do i=0,NI+1
            he(i,j)= 0.d0
         end do
      end do

c     Form the new RHS : A*(h +dh) - b
      do j=1,NJ
         i=1
         rhs(i,j) = 
     &           (ch(1,i,j)*h(i,j)
     &           +ch(2,i,j)*h(i+1,j)
     &           +ch(3,i,j)*h(NI,j)
     &           +ch(4,i,j)*h(i,j+1)
     &           +ch(5,i,j)*h(i,j-1)
     &           +ch(6,i,j)*h(i+1,j+1)
     &           +ch(7,i,j)*h(NI,j+1)
     &           +ch(8,i,j)*h(i+1,j-1)
     &           +ch(9,i,j)*h(NI,j-1))
     &           -oldrhs(i,j)
         do i=2,NI-1
            rhs(i,j) =
     &           (ch(1,i,j)*h(i,j)
     &           +ch(2,i,j)*h(i+1,j)
     &           +ch(3,i,j)*h(i-1,j)
     &           +ch(4,i,j)*h(i,j+1)
     &           +ch(5,i,j)*h(i,j-1)
     &           +ch(6,i,j)*h(i+1,j+1)
     &           +ch(7,i,j)*h(i-1,j+1)
     &           +ch(8,i,j)*h(i+1,j-1)
     &           +ch(9,i,j)*h(i-1,j-1))
     &           -oldrhs(i,j)
         end do
         i=NI
         rhs(i,j) = 
     &        (ch(1,i,j)*h(i,j)
     &        +ch(2,i,j)*h(1,j)
     &        +ch(3,i,j)*h(i-1,j)
     &        +ch(4,i,j)*h(i,j+1)
     &        +ch(5,i,j)*h(i,j-1)
     &        +ch(6,i,j)*h(1,j+1)
     &        +ch(7,i,j)*h(i-1,j+1)
     &        +ch(8,i,j)*h(1,j-1)
     &        +ch(9,i,j)*h(i-1,j-1))
     &        -oldrhs(i,j)
      end do

c     Now solve for the error : he
c     -------------------------


      rlx= 1.72
 999  do 1000 iter=1,1000
C$DOACROSSSHARE(ch,he,rlx,rhs),
C$&LOCAL(i,j,hstar)
         do 120 j=1,NJ
            i=1
            hstar= ( ch(2,i,j)*he(i+1,j)
     &           +ch(3,i,j)*he(NI,j)
     &           +ch(4,i,j)*he(i,j+1)
     &           +ch(5,i,j)*he(i,j-1)
     &           +ch(6,i,j)*he(i+1,j+1)
     &           +ch(7,i,j)*he(NI,j+1)
     &           +ch(8,i,j)*he(i+1,j-1)
     &           +ch(9,i,j)*he(NI,j-1)
     &           - rhs(i,j) )/(-ch(1,i,j))
            he(i,j)= (1.d0 -rlx)*he(i,j) +rlx*hstar
            do 130 i=2,NI-1
               hstar= ( ch(2,i,j)*he(i+1,j)
     &              +ch(3,i,j)*he(i-1,j)
     &              +ch(4,i,j)*he(i,j+1)
     &              +ch(5,i,j)*he(i,j-1)
     &              +ch(6,i,j)*he(i+1,j+1)
     &              +ch(7,i,j)*he(i-1,j+1)
     &              +ch(8,i,j)*he(i+1,j-1)
     &              +ch(9,i,j)*he(i-1,j-1)
     &              - rhs(i,j) )/(-ch(1,i,j))
               he(i,j)= (1.d0 -rlx)*he(i,j) +rlx*hstar
 130        continue
            i=NI
            hstar= ( ch(2,i,j)*he(1,j)
     &           +ch(3,i,j)*he(i-1,j)
     &           +ch(4,i,j)*he(i,j+1)
     &           +ch(5,i,j)*he(i,j-1)
     &           +ch(6,i,j)*he(1,j+1)
     &           +ch(7,i,j)*he(i-1,j+1)
     &           +ch(8,i,j)*he(1,j-1)
     &           +ch(9,i,j)*he(i-1,j-1)
     &           - rhs(i,j) )/(-ch(1,i,j))
            he(i,j)= (1.d0 -rlx)*he(i,j) +rlx*hstar
 120     continue
c     
         
c     Boundary condition :
c     Since the old bc dh/dn is done perfectly, d(herror)/dn = 0
c
         do i=1,NI
            he(i,0) = he(i,1)
            he(i,NJ+1) = he(i,NJ)
         end do
         do j=0,NJ+1
            he(0,j)= he(1,j)
            he(NI+1,j)= he(NI,j)
         end do
            
         maxres= 0.d0
C$DOACROSSSHARE(ch,he,rhs,maxres),
C$&LOCAL(i,j,res)
         do 220 j=1,NJ
            do 230 i=1,NI
c     res = b - A x'     
               res= rhs(i,j)- 
     &              (ch(1,i,j)*he(i,j)
     &              +ch(2,i,j)*he(i+1,j)
     &              +ch(3,i,j)*he(i-1,j)
     &              +ch(4,i,j)*he(i,j+1)
     &              +ch(5,i,j)*he(i,j-1)
     &              +ch(6,i,j)*he(i+1,j+1)
     &              +ch(7,i,j)*he(i-1,j+1)
     &              +ch(8,i,j)*he(i+1,j-1)
     &              +ch(9,i,j)*he(i-1,j-1))
               if (dabs(res).gt.maxres) then
                  maxres= dabs(res)
               end if
 230        continue
 220     continue
c     
c         write(6,*) 'iter,maxres',iter,maxres
         if (maxres.gt.3000.d0) then
            write(6,*) 'STOP in mprove. res too large, i,j,maxres=',
     &           i,j,maxres
            stop
         end if
         if (maxres.lt.tol) goto 13
 1000 continue
 13   continue
c     write(6,*) 'iter=',iter,maxres
      
      do j=0,NJ+1
         do i=0,NI+1
            h(i,j)= h(i,j) - he(i,j)
         end do
      end do

      return
      end
