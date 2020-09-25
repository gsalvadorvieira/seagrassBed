      subroutine writeframe(frame_int,step,h,s,T,rho,u,v,w,vor,
     &     xc,yc,zc)
c     ---------------------------------------------------
c     reads in a netcdf file and interpolates the data from the sigma grid
c     onto a z level

C      INCLUDE '/usr/local/include/netcdf.inc'

      INCLUDE '/home/salvadorvieira.g/a/include/netcdf.inc'
C      INCLUDE '../netcdf-fortran-4.5.2_install/include/netcdf.inc'

c      INCLUDE '/gaeaopt/netcdf-3.4atd/include/netcdf.inc'
c      INCLUDE
c     &     '/home/scratch/vayu/public/netcdf-3.5.0/include/netcdf.inc'
      include 'dims.f'
      integer NI2,NJ2,NK2,i,j,k,step,frame_int,itr,ntr
      parameter ( NI2=NI+2,NJ2=NJ+2,NK2=NK+2)
      double precision DL,zmax
      parameter (DL= 1000.d0)
      double precision svert(NK),Tvert(NK,ntr),rvert(NK),uvert(NK),
     &     vvert(NK),wvert(NK),zvert(NK),vorvert(NK),seval

      double precision h(0:NI+1,0:NJ+1),
     &     rho(0:NI+1,0:NJ+1,0:NK+1),u(0:NI+1,0:NJ+1,0:NK+1),
     &     v(0:NI+1,0:NJ+1,0:NK+1),w(0:NI+1,0:NJ+1,0:NK+1),
     &     vor(0:NI+1,0:NJ+1,0:NK+1),
     &     s(0:NI+1,0:NJ+1,0:NK+1),T(ntr,0:NI+1,0:NJ+1,0:NK+1),
     &     zc(0:NI+1,0:NJ+1,0:NK+1),xc(0:NI+1),
     &     yc(0:NJ+1),bs1(NK),cs1(NK),ds1(NK)
      real*4 xsurf(NI),ysurf(NJ),zsurf(NI,NJ),ssurf(NI,NJ),
     &     Tsurf(NI,NJ,ntr),rsurf(NI,NJ),usurf(NI,NJ),vsurf(NI,NJ),
     &     wsurf(NI,NJ),hsurf(NI,NJ),vorsurf(NI,NJ)
c

      character*80 outname

      integer start(3),count(3),dims(3),start2d(2),count2d(2),
     &     dims4(4),start4(4),count4(4)
      zmax= 0.d0
c      zmax= -1.d16
cc     find zmax
c      do j=1,NJ
c         do i=1,NI
c            if (zc(i,j,NK).gt.zmax) zmax= zc(i,j,NK)
c         end do
c      end do
c     write(6,*) 'interpolate data to ' , zmax*DL , 'm'

      do i=1,NI
         xsurf(i)= xc(i)
      end do
      do j=1,NJ
         ysurf(j)= yc(j)
      end do
      do j=1,NJ
         do i=1,NI
            zsurf(i,j)= zmax*DL
            hsurf(i,j)= h(i,j)
            do k=1,NK
               svert(k)= s(i,j,k) +S0
c               Tvert(k)= T(i,j,k) +T0
               rvert(k)= rho(i,j,k) +R0 -1000.
               uvert(k)= u(i,j,k)
               vvert(k)= v(i,j,k)
               wvert(k)= w(i,j,k)
               vorvert(k)= vor(i,j,k)
               zvert(k)= zc(i,j,k)
               do itr=1,ntr
                  Tvert(k,itr)= T(itr,i,j,k)
               end do
            end do
c=            call spline (NK,zvert,svert,bs1,cs1,ds1)
c=            ssurf(i,j)= seval(NK,zmax,zvert,svert,bs1,cs1,ds1)

c=            call spline (NK,zvert,Tvert,bs1,cs1,ds1)
c=            Tsurf(i,j)= seval(NK,zmax,zvert,Tvert,bs1,cs1,ds1)

            do itr=1,ntr
               call spline (NK,zvert,Tvert(:,itr),bs1,cs1,ds1)
               Tsurf(i,j,itr)= seval(NK,zmax,zvert,Tvert(:,itr),bs1,
     &              cs1,ds1)
            end do

            call spline (NK,zvert,rvert,bs1,cs1,ds1)
            rsurf(i,j)= seval(NK,zmax,zvert,rvert,bs1,cs1,ds1)

            call spline (NK,zvert,uvert,bs1,cs1,ds1)
            usurf(i,j)= seval(NK,zmax,zvert,uvert,bs1,cs1,ds1)

            call spline (NK,zvert,vvert,bs1,cs1,ds1)
            vsurf(i,j)= seval(NK,zmax,zvert,vvert,bs1,cs1,ds1)

            call spline (NK,zvert,wvert,bs1,cs1,ds1)
            wsurf(i,j)= seval(NK,zmax,zvert,wvert,bs1,cs1,ds1)

            call spline (NK,zvert,vorvert,bs1,cs1,ds1)
            vorsurf(i,j)= seval(NK,zmax,zvert,vorvert,bs1,cs1,ds1)

         end do
      end do

c     write to a netcdf file

      outname= 'surfmovie.cdf'

      if (step.eq.0) then
         idDatFile =  nccre(outname,NCCLOB,rcode)

c         dims2d(1) = ncddef(idDatFile,'x',NI,rcode)
c         dims2d(2) = ncddef(idDatFile,'y',NJ,rcode)

         dims(1) = ncddef(idDatFile,'x',NI,rcode)
         dims(2) = ncddef(idDatFile,'y',NJ,rcode)
         dims(3) = ncddef(idDatFile,'time',NCUNLIM,rcode)

         dims4(1) = dims(1)
         dims4(2) = dims(2)
         dims4(3) = ncddef(idDatFile,'ntr',ntr,rcode)
         dims4(4) = dims(3)

         idvx = ncvdef(idDatFile,'xc',NCFLOAT,1,dims(1),rcode)
         idvy = ncvdef(idDatFile,'yc',NCFLOAT,1,dims(2),rcode)
         idvz = ncvdef(idDatFile,'zc',NCFLOAT,2,dims,rcode)
         idvh = ncvdef(idDatFile,'h',NCFlOAT,3,dims,rcode)
c=         idvs = ncvdef(idDatFile,'s',NCFLOAT,3,dims,rcode)
c=         idvt = ncvdef(idDatFile,'temp',NCFLOAT,3,dims,rcode)
         idvt = ncvdef(idDatFile,'tr',NCFLOAT,4,dims4,rcode)
         idvrho = ncvdef(idDatFile,'rho',NCFLOAT,3,dims,rcode)
         idvu = ncvdef(idDatFile,'u',NCFLOAT,3,dims,rcode)
         idvv = ncvdef(idDatFile,'v',NCFLOAT,3,dims,rcode)
         idvw = ncvdef(idDatFile,'w',NCFLOAT,3,dims,rcode)
         idvvor = ncvdef(idDatFile,'vor',NCFLOAT,3,dims,rcode)

         CALL ncendf(idDatFile,rcode)

      else
         idDatFile = ncopn(outname, NCWRITE,rcode)
c         call ncredf(idDatFile,rcode)
c         idvx = NCVID(idDatFile, 'xc', RCODE)
c         idvy = NCVID(idDatFile, 'yc', RCODE)
c         idvz = NCVID(idDatFile, 'zc', RCODE)
         idvh = NCVID(idDatFile, 'h', RCODE)
c=         idvs = NCVID(idDatFile, 's', RCODE)
c=         idvt = NCVID(idDatFile, 'temp', RCODE)
         idvt = NCVID(idDatFile, 'tr', RCODE)
         idvrho = NCVID(idDatFile, 'rho', RCODE)
         idvu = NCVID(idDatFile, 'u', RCODE)
         idvv = NCVID(idDatFile, 'v', RCODE)
         idvw = NCVID(idDatFile, 'w', RCODE)
         idvvor = NCVID(idDatFile, 'vor', RCODE)
      endif

      count2d(1)= NI
      count2d(2)= NJ

      count(1)= NI
      count(2)= NJ
      count(3)= 1

      count4(1)= NI
      count4(2)= NJ
      count4(3)= ntr
      count4(4)= 1

      start2d(1)= 1
      start2d(2)= 1

      start(1)= 1
      start(2)= 1
      start(3)= step/frame_int +1

      start4(1)= 1
      start4(2)= 1
      start4(3)= 1
      start4(4)= step/frame_int +1

      if (step.eq.0) then
         CALL ncvpt(idDatFile,idvx, start(1), count(1), xsurf, rcode)
         CALL ncvpt(idDatFile,idvy, start(2), count(2), ysurf, rcode)
         CALL ncvpt(idDatFile,idvz, start2d, count2d, zsurf, rcode)
      endif
      CALL ncvpt(idDatFile,idvh, start, count, hsurf, rcode)
c=      CALL ncvpt(idDatFile,idvs, start, count, ssurf, rcode)
c=      CALL ncvpt(idDatFile,idvT, start, count, Tsurf, rcode)
      CALL ncvpt(idDatFile,idvt, start4, count4, Tsurf, rcode)
      CALL ncvpt(idDatFile,idvrho, start, count, rsurf, rcode)
      CALL ncvpt(idDatFile,idvu, start, count, usurf, rcode)
      CALL ncvpt(idDatFile,idvv, start, count, vsurf, rcode)
      CALL ncvpt(idDatFile,idvw, start, count, wsurf, rcode)
      CALL ncvpt(idDatFile,idvvor, start, count, vorsurf, rcode)

      CALL ncclos(idDatFile,rcode)

      end
