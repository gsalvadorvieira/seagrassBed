      subroutine checks
c     ------------------
      include 'header.f'
      integer i,j,k
c     Checks and writes a few things
      

      do j=1,NJ
         do i=0,NI
            if (dabs(D(i+1,j)-D(i,j)).gt.1.d-4) then
               write(6,*) 'topog in both directions'
               goto 202
c               write(6,*) i,j,D(i+1,j),D(i,j)
c               write(6,*) 'need to modify rpevalgrad (rgradients) for
C      &              slope in x direcn: grpifc needs to be modified'
c               write(6,*) 'need to modify velbc_periodicew'
c               stop
            end if
         end do
      end do

 202  do j=0,NJ+1
         do i=0,NI+1
            if ((uy(i,j).ne.0.d0).or.(vx(i,j).ne.0.d0)) then
               write(6,*) 'Need to modify biharmonic for curvi grid'
               stop
            end if
         end do
      end do

      if (rect.eqv.(.false.)) then
c        need to modify grpifc,grpjfc to contain cross diff terms
         write(6,*) 'modify grpifc,grpjfc, stop in checks'
         stop
      end if
c
c      open(unit=88,file='umax.out',status='new')
c      open(unit=89,file='umin.out',status='new')
c      close(88 89)
c      open(unit=79,file='vtimeser.out',status='new')
c      open(unit=78,file='utimeser.out',status='new')
c      close(78 79)

      return
      end 
