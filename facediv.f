      subroutine facediv(dtime,maxdiv)
c     ------------------
c     checks divergence
c     wf does not contain  J*wt and is not exactly the contravariant vel W
      include 'header.f'
      integer i,j,k,imax,jmax,kmax
      double precision ufdx,vfdy,wfdz,div,maxdiv,dtime
c     edd is a scaling factor that makes this divergence
c     the same as the residual from the pressure-Poisson equation. 
c     This is a check on the code.
c
c      double precision dtheta,dphi,phi0,cosphi0,sinphi0,C,dx,dy,dz,
c     & lon0,lat0
c      C= PI/180.d0
cc
c      lon0= 260.d0*C
c      lat0= 20.d0*C
c      dtheta= C/8.d0
c      dphi= C/8.d0
c      phi0= lat0 + 0.5d0*(dphi*dfloat(NJ))
c      cosphi0= dCos(phi0)
c      sinphi0= dSin(phi0)
c      write(6,*) dphi,dtheta,phi0,cosphi0,sinphi0
cc     for the rectilinear case, dx,dy,dz are constant
c      dx= R*cosphi0*dtheta /LEN
c      dy= R*dphi/LEN
c      dz= 1000.d0/(dfloat(NK)*DL)
c      write(6,*) 'dx dy',dx,dy
c
c      edd= EPS/(dtime*delta*kappa)
c     We have now changed the defn of tol so that it is just (u_x +v_y +ep*w_z)
c     We do not need edd any more.
c      edd= EPS*delinv/dtime
      maxdiv= 0.d0
      do 10 k=1,NK
         do 20 j=1,NJ
            do 30 i=1,NI
               ufdx= (uf(i,j,k)-uf(i-1,j,k))
               vfdy= (vf(i,j,k)-vf(i,j-1,k))
               wfdz= (wf(i,j,k)-wf(i,j,k-1))
c               div= dabs((ufdx+ vfdy + wfdz)*edd)
               div= dabs(ufdx+ vfdy + wfdz)
               if (div.gt.maxdiv) then
                  maxdiv=div
                  imax=i
                  jmax=j
                  kmax=k
               end if
 30         continue
 20      continue
 10   continue
c      write(6,*) 'in facediv, i,j,k',imax,jmax,kmax
      return
      end
      
