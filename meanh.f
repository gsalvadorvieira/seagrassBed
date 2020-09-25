      subroutine meanh(NI,NJ,h,hmean)
c     -------------------------------
c     computes the mean free-surface elevation within the domain
      integer NI,NJ,i,j
      double precision h(0:NI+1,0:NJ+1),hsum,hmean
c
      hsum=0.d0
      do 45 i=1,NI
         do 55 j=1,NJ
            hsum= hsum +h(i,j)
 55      continue
 45   continue
      hmean= hsum/dble(NI*NJ)
      if (dabs(hmean).lt.100.d0) then
         continue
      else
c     if "hmean" is nan, the program will stop
c         call output(p,step,0)
c         call finalout(p,step,0)
         write(6,*) 'error, hmean=',hmean
         stop
      end if
c
      return
      end
