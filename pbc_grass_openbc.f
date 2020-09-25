      subroutine pbc(cpf,fn,dtime)
c-------------------------------------------------------------------c
c     evaluates the boundaries in psolve c fills in
c     cp0,cpim1,cpip1,cpjm1,cpjp1,cpkm1,cpkp1,fn for the boundaries
      implicit logical (a-z)
      include 'header.f'
c     3-D Poisson solver - uses SOR
c     Written for the control volume formulation
c     Uses the 19 point stencil (in 3d).
c     
c     fn is the contribution to the source terms on the RHS.
c     Whenever possible we substitue the the velocities themselves
c     (at the boundaries). We cannot substitute the u,w vels at the entrance.
c     Our eqn is
c     cpim1*p(i-1) +cpip1*p(i+1) +cpjm1p(j-1) ... -cp0*p(i,j,k)= f
      integer i,j,k,l
      double precision dtime,edt,cpf(19,NI,NJ,NK),
     &     fn(NI,NJ,NK)
c     
c     delinv= 1/delta, delta= D/L
      edt= EPS/dtime
c     
c     *****************************
c     do the edges separately
      do 40 k=1,NK
         do 50 j=1,NJ
            do 60 i=1,NI
c     periodic_ew boundaries
c     this would be:   if ((j.eq.1).or.(j.eq.NJ).or.(k.eq.1).or.(k.eq.NK)) then

               if ((i.eq.1).or.(i.eq.NI).or.(j.eq.1).or. 
     &              (j.eq.NJ).or.(k.eq.1).or.(k.eq.NK)) then
                  fn(i,j,k)= 0.d0
                  do l=1,19
                     cpf(l,i,j,k)= 0.d0
                  end do
                  if (i.eq.NI) then
                     fn(i,j,k)= fn(i,j,k) + edt*ufbce(j,k)
                  else
c     add values for the east face if (i.ne.NI) 
                     fn(i,j,k)= fn(i,j,k) + edt*cxf(i,j,k)
c     
                     cpf(1,i,j,k)= cpf(1,i,j,k) - gqi(i,j,k,1)
                     cpf(2,i,j,k)= cpf(2,i,j,k) +gqi(i,j,k,1)
                     cpf(3,i,j,k)= cpf(3,i,j,k) +0.25*gqi(i,j,k,2)
                     cpf(5,i,j,k)= cpf(5,i,j,k) -0.25*gqi(i,j,k,2)
                     cpf(6,i,j,k)= cpf(6,i,j,k) +0.25*gqi3(i,j,k)
                     cpf(7,i,j,k)= cpf(7,i,j,k) -0.25*gqi3(i,j,k)
                     cpf(10,i,j,k)= cpf(10,i,j,k) -0.25*gqi(i,j,k,2)
                     cpf(11,i,j,k)= cpf(11,i,j,k) +0.25*gqi(i,j,k,2)
                     cpf(13,i,j,k)= cpf(13,i,j,k) -0.25*gqi3(i,j,k)
                     cpf(14,i,j,k)= cpf(14,i,j,k) +0.25*gqi3(i,j,k)
                  end if
                  if (i.eq.1) then
                     fn(i,j,k)= fn(i,j,k) - edt*ufbcw(j,k)
                  else
c     add values for the west face if (i.ne.1)
                     fn(i,j,k)=  fn(i,j,k) -edt*cxf(i-1,j,k)
c     
                     cpf(1,i,j,k)= cpf(1,i,j,k) -gqi(i-1,j,k,1)
                     cpf(3,i,j,k)= cpf(3,i,j,k) -0.25*gqi(i-1,j,k,2)
                     cpf(4,i,j,k)= cpf(4,i,j,k) +gqi(i-1,j,k,1)
                     cpf(5,i,j,k)= cpf(5,i,j,k) +0.25*gqi(i-1,j,k,2)
                     cpf(6,i,j,k)= cpf(6,i,j,k) -0.25*gqi3(i-1,j,k)
                     cpf(7,i,j,k)= cpf(7,i,j,k) +0.25*gqi3(i-1,j,k)
                     cpf(8,i,j,k)= cpf(8,i,j,k) -0.25*gqi(i-1,j,k,2)
                     cpf(9,i,j,k)= cpf(9,i,j,k) +0.25*gqi(i-1,j,k,2)
                     cpf(12,i,j,k)= cpf(12,i,j,k) +0.25*gqi3(i-1,j,k)
                     cpf(15,i,j,k)= cpf(15,i,j,k) -0.25*gqi3(i-1,j,k)
                  end if
                  if (j.eq.NJ) then
                     fn(i,j,k)= fn(i,j,k) + edt*vfbcn(i,k)
                  else
c     add values for the north face if (j.ne.NJ)
                     fn(i,j,k)=  fn(i,j,k) +edt*cyf(i,j,k)
c
                     cpf(1,i,j,k)= cpf(1,i,j,k) -gqj(i,j,k,2)
                     cpf(2,i,j,k)= cpf(2,i,j,k) +0.25*gqj(i,j,k,1)
                     cpf(3,i,j,k)= cpf(3,i,j,k) +gqj(i,j,k,2)
                     cpf(4,i,j,k)= cpf(4,i,j,k) -0.25*gqj(i,j,k,1)
                     cpf(6,i,j,k)= cpf(6,i,j,k) +0.25*gqj3(i,j,k)
                     cpf(7,i,j,k)= cpf(7,i,j,k) -0.25*gqj3(i,j,k)
                     cpf(8,i,j,k)= cpf(8,i,j,k) -0.25*gqj(i,j,k,1)
                     cpf(11,i,j,k)= cpf(11,i,j,k) +0.25*gqj(i,j,k,1)
                     cpf(17,i,j,k)= cpf(17,i,j,k) -0.25*gqj3(i,j,k)
                     cpf(18,i,j,k)= cpf(18,i,j,k) +0.25*gqj3(i,j,k)
                  endif
                  if (j.eq.1) then
                     fn(i,j,k)= fn(i,j,k) - edt*vfbcs(i,k)
                  else
c     add values for the south face if (j.ne.1)
                     fn(i,j,k)=  fn(i,j,k) -edt*cyf(i,j-1,k)
c
                     cpf(1,i,j,k)= cpf(1,i,j,k) -gqj(i,j-1,k,2)
                     cpf(2,i,j,k)= cpf(2,i,j,k) -0.25*gqj(i,j-1,k,1)
                     cpf(4,i,j,k)= cpf(4,i,j,k) +0.25*gqj(i,j-1,k,1)
                     cpf(5,i,j,k)= cpf(5,i,j,k) +gqj(i,j-1,k,2)
                     cpf(6,i,j,k)= cpf(6,i,j,k) -0.25*gqj3(i,j-1,k)
                     cpf(7,i,j,k)= cpf(7,i,j,k) +0.25*gqj3(i,j-1,k)
                     cpf(9,i,j,k)= cpf(9,i,j,k) +0.25*gqj(i,j-1,k,1)
                     cpf(10,i,j,k)= cpf(10,i,j,k) -0.25*gqj(i,j-1,k,1)
                     cpf(16,i,j,k)= cpf(16,i,j,k) +0.25*gqj3(i,j-1,k)
                     cpf(19,i,j,k)= cpf(19,i,j,k) -0.25*gqj3(i,j-1,k)
                  end if                     
                  if  (k.eq.NK) then
                     fn(i,j,k)= fn(i,j,k) + edt*czf(i,j,k) 
                     cpf(1,i,j,k)= cpf(1,i,j,k) -2.d0*gqk(i,j,k,3)
                  else
c     add values for the top face if (k.ne.NK)
cc                     temp= fn(i,j,k)
                     fn(i,j,k)=  fn(i,j,k) +edt*czf(i,j,k)
c
                     cpf(1,i,j,k)= cpf(1,i,j,k) -gqk(i,j,k,3)
                     cpf(2,i,j,k)= cpf(2,i,j,k) +0.25*gqk(i,j,k,1)
                     cpf(3,i,j,k)= cpf(3,i,j,k) +0.25*gqk(i,j,k,2)
                     cpf(4,i,j,k)= cpf(4,i,j,k) -0.25*gqk(i,j,k,1)
                     cpf(5,i,j,k)= cpf(5,i,j,k) -0.25*gqk(i,j,k,2)
                     cpf(6,i,j,k)= cpf(6,i,j,k) +gqk(i,j,k,3)
                     cpf(14,i,j,k)= cpf(14,i,j,k) +0.25*gqk(i,j,k,1)
                     cpf(15,i,j,k)= cpf(15,i,j,k) -0.25*gqk(i,j,k,1)
                     cpf(18,i,j,k)= cpf(18,i,j,k) +0.25*gqk(i,j,k,2)
                     cpf(19,i,j,k)= cpf(19,i,j,k) -0.25*gqk(i,j,k,2)
                  endif
                  if (k.eq.1) then
                     fn(i,j,k)= fn(i,j,k) - edt*wfbcb(i,j)
cc    since wf(k=0) is zero, we need not write this line
                  else
c                  if (k.ne.1) then
c     add values for the bottom face if (k.ne.1)
cc                     temp= fn(i,j,k)
                     fn(i,j,k)=  fn(i,j,k) -edt*czf(i,j,k-1)
c
cc                     if (k.eq.NK) write(400,*) temp,fn(i,j,k)
                     cpf(1,i,j,k)= cpf(1,i,j,k) - gqk(i,j,k-1,3) 
                     cpf(2,i,j,k)= cpf(2,i,j,k) -0.25*gqk(i,j,k-1,1)
                     cpf(3,i,j,k)= cpf(3,i,j,k) -0.25*gqk(i,j,k-1,2)
                     cpf(4,i,j,k)= cpf(4,i,j,k) +0.25*gqk(i,j,k-1,1)
                     cpf(5,i,j,k)= cpf(5,i,j,k) +0.25*gqk(i,j,k-1,2)
                     cpf(7,i,j,k)= cpf(7,i,j,k) +gqk(i,j,k-1,3)
                     cpf(12,i,j,k)= cpf(12,i,j,k) +0.25*gqk(i,j,k-1,1)
                     cpf(13,i,j,k)= cpf(13,i,j,k) -0.25*gqk(i,j,k-1,1)
                     cpf(16,i,j,k)= cpf(16,i,j,k) +0.25*gqk(i,j,k-1,2)
                     cpf(17,i,j,k)= cpf(17,i,j,k) -0.25*gqk(i,j,k-1,2)
                  endif
               endif
 60         continue
 50      continue
 40   continue
c
      return
      end
      
      
