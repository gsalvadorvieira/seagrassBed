      subroutine evalrho(rhonew,n)
c     ---------------------------------------------
c     rho is the same as s(i,j,k,n)

      implicit logical (a-z)
      include 'header.f'
      integer i,j,k,n
      double precision 
     &     rhonew(0:NI+1,0:NJ+1,0:NK+1)

      do 10 i=0,NI+1
         do 20 j=0,NJ+1
            do 30 k=0,NK+1
               rhonew(i,j,k)= s(i,j,k,n)
 30         continue
 20      continue
 10   continue
c
      return
      end



