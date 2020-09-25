      subroutine cpfine(dtime,cpf,fn)
c     -----------------------------------------------
c     call cpfine(dtime,n,cp(1),rhs(1))
c     computes the coefficients in the Laplace operator
c     uses the 19 point stencil (in 3d).
c
      implicit logical (a-z)
      include 'header.f'
c
c     cp0 is the coeff of p(i,j,k)
c     cpim1 is the coeff of p(i-1,j,k), cpip1 is the coeff of p(i+1,j,k),etc
c     fn is the contribution to the source terms on RHS.
c     Whenever possible we substitue the the velocities themselves
c     (at the boundaries). We cannot substitute the u,w vels at the entrance.
c     Our eqn is
c     cpim1*p(i-1) +cpip1*p(i+1) +cpjm1*p(j-1) ...- cp0*p(i,j,k)= fn
c
c     cpf(1,..)=cpkm1(..), cpf(2,..)=cpjp1(..), cpf(3,..)=cpip1(..),
c     cpf(4,..)=cpjm1(..), cpf(5,..)=cpim1(..), cpf(6,..)=cp0(..),
c     cpf(7,..)=cpkp1(..)
      integer i,j,k
      double precision cpf(19,NI,NJ,NK),fn(NI,NJ,NK),edt,dtime
c
      edt= EPS/dtime
c      delinv= 1/delta
c
      call pbc(cpf,fn,dtime)
c      call pbc(cpf,fn,dtime,n)
c         
c     For i=2...NI-1, j=2...NJ-1, k=2,NK-1  it is straightforward
      do 10 i=2,NI-1 
c      do 10 i=1,NI ! for periodicew bc
         do 20 j=2,NJ-1
            do 30 k=2,NK-1
               fn(i,j,k)=  edt*( cxf(i,j,k)-cxf(i-1,j,k)
     &              + cyf(i,j,k) - cyf(i,j-1,k)
     &              + czf(i,j,k) - czf(i,j,k-1) )
c
             cpf(1,i,j,k)= -(gqi(i-1,j,k,1)+gqi(i,j,k,1) 
     &              +gqj(i,j-1,k,2)+
     &              gqj(i,j,k,2) +gqk(i,j,k-1,3) +gqk(i,j,k,3))
             cpf(2,i,j,k)= gqi(i,j,k,1) +0.25*(gqj(i,j,k,1) 
     &            -gqj(i,j-1,k,1) +gqk(i,j,k,1) -gqk(i,j,k-1,1))
             cpf(3,i,j,k)= gqj(i,j,k,2) +0.25*(gqi(i,j,k,2) 
     &           -gqi(i-1,j,k,2) +gqk(i,j,k,2) -gqk(i,j,k-1,2))
             cpf(4,i,j,k)= gqi(i-1,j,k,1) +0.25*(-gqj(i,j,k,1) 
     &            +gqj(i,j-1,k,1) -gqk(i,j,k,1) +gqk(i,j,k-1,1))
             cpf(5,i,j,k)= gqj(i,j-1,k,2) +0.25*(-gqi(i,j,k,2)
     &            +gqi(i-1,j,k,2) -gqk(i,j,k,2) +gqk(i,j,k-1,2))
             cpf(6,i,j,k)= gqk(i,j,k,3) +0.25*(gqi3(i,j,k)
     &            -gqi3(i-1,j,k) +gqj3(i,j,k) -gqj3(i,j-1,k))
             cpf(7,i,j,k)= gqk(i,j,k-1,3) +0.25*(-gqi3(i,j,k)
     &            +gqi3(i-1,j,k) -gqj3(i,j,k) +gqj3(i,j-1,k))
c
             cpf(8,i,j,k)= 0.25*(-gqi(i-1,j,k,2) -gqj(i,j,k,1))
             cpf(9,i,j,k)= 0.25*(gqi(i-1,j,k,2) +gqj(i,j-1,k,1))
             cpf(10,i,j,k)= 0.25*(-gqi(i,j,k,2) -gqj(i,j-1,k,1))
             cpf(11,i,j,k)= 0.25*(gqi(i,j,k,2) +gqj(i,j,k,1))
             cpf(12,i,j,k)= 0.25*(gqi3(i-1,j,k) +gqk(i,j,k-1,1))
             cpf(13,i,j,k)= 0.25*(-gqi3(i,j,k) -gqk(i,j,k-1,1))
             cpf(14,i,j,k)= 0.25*(gqi3(i,j,k) +gqk(i,j,k,1))
             cpf(15,i,j,k)= 0.25*(-gqi3(i-1,j,k) -gqk(i,j,k,1))
             cpf(16,i,j,k)= 0.25*(gqj3(i,j-1,k) +gqk(i,j,k-1,2))
             cpf(17,i,j,k)= 0.25*(-gqj3(i,j,k) -gqk(i,j,k-1,2))
             cpf(18,i,j,k)= 0.25*(gqj3(i,j,k) +gqk(i,j,k,2))
             cpf(19,i,j,k)= 0.25*(-gqj3(i,j-1,k) -gqk(i,j,k,2))
 30         continue
 20      continue
 10   continue

      return
      end
      

