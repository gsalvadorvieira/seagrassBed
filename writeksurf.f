      subroutine writeksurf
     &     (frame_int,step,k,s,T,rho,u,v,w,vor,xc,yc,zc)
c     ------------------------------------------------------------------
c     periodically writes out the data on the k-th surface

c      INCLUDE '/usr/local/include/netcdf.inc'

      INCLUDE '/home/salvadorvieira.g/a/include/netcdf.inc'
C      INCLUDE '../netcdf-fortran-4.5.2_install/include/netcdf.inc'

c      INCLUDE '/gaeaopt/netcdf-3.4atd/include/netcdf.inc'
c      INCLUDE
c     &     '/home/scratch/vayu/public/netcdf-3.5.0/include/netcdf.inc'
      include 'dims.f'
      integer NI2,NJ2,NK2,i,j,k,step,frame_int
      parameter ( NI2=NI+2,NJ2=NJ+2,NK2=NK+2)
      double precision DL,zmax
      parameter (DL= 1000.d0)

      double precision h(0:NI+1,0:NJ+1),
     &     rho(0:NI+1,0:NJ+1,0:NK+1),u(0:NI+1,0:NJ+1,0:NK+1),
     &     v(0:NI+1,0:NJ+1,0:NK+1),w(0:NI+1,0:NJ+1,0:NK+1),
     &     vor(0:NI+1,0:NJ+1,0:NK+1),
     &     s(0:NI+1,0:NJ+1,0:NK+1),T(0:NI+1,0:NJ+1,0:NK+1),
     &     zc(0:NI+1,0:NJ+1,0:NK+1),xc(0:NI+1),
     &     yc(0:NJ+1)
      real*4 xsurf(NI),ysurf(NJ),zsurf(NI,NJ),ssurf(NI,NJ),
     &     Tsurf(NI,NJ),rsurf(NI,NJ),usurf(NI,NJ),vsurf(NI,NJ),
     &     wsurf(NI,NJ),hsurf(NI,NJ),vorsurf(NI,NJ)
c

      character*80 outname

      integer start(3),count(3),dims(3),start2d(2),count2d(2),dims2d(2)

      do j=1,NJ
         do i=1,NI
            zsurf(i,j)= zc(i,j,k)*DL
c            xsurf(i,j)= xc(i)
c            ysurf(i,j)= yc(j)
c            hsurf(i,j)= h(i,j)
            rsurf(i,j)= rho(i,j,k)
            usurf(i,j)= u(i,j,k)
            vsurf(i,j)= v(i,j,k)
            wsurf(i,j)= w(i,j,k)
            vorsurf(i,j)= vor(i,j,k)
         end do
      end do
      do i=1,NI
         xsurf(i)= xc(i)
      end do
      do j=1,NJ
         ysurf(j)= yc(j)
      end do

c     write to a netcdf file

      outname= 'movk24surf.cdf'

      if (step.eq.0) then
         idDatFile =  nccre(outname,NCCLOB,rcode)

c         dims2d(1) = ncddef(idDatFile,'xi',NI,rcode)
c         dims2d(2) = ncddef(idDatFile,'eta',NJ,rcode)

         dims(1) = ncddef(idDatFile,'x',NI,rcode)
         dims(2) = ncddef(idDatFile,'y',NJ,rcode)
         dims(3) = ncddef(idDatFile,'time',NCUNLIM,rcode)

         idvx = ncvdef(idDatFile,'xc',NCFLOAT,1,dims(1),rcode)
         idvy = ncvdef(idDatFile,'yc',NCFLOAT,1,dims(2),rcode)
         idvz = ncvdef(idDatFile,'zc',NCFLOAT,2,dims,rcode)
c         idvh = ncvdef(idDatFile,'h',NCFlOAT,3,dims,rcode)
c=         idvs = ncvdef(idDatFile,'s',NCFLOAT,3,dims,rcode)
c=         idvt = ncvdef(idDatFile,'temp',NCFLOAT,3,dims,rcode)
         idvrho = ncvdef(idDatFile,'rho',NCFLOAT,3,dims,rcode)
         idvu = ncvdef(idDatFile,'u',NCFLOAT,3,dims,rcode)
         idvv = ncvdef(idDatFile,'v',NCFLOAT,3,dims,rcode)
         idvw = ncvdef(idDatFile,'w',NCFLOAT,3,dims,rcode)
         idvvor = ncvdef(idDatFile,'vor',NCFLOAT,3,dims,rcode)

c         write(6,*) 'defined 2d vars'
         CALL ncendf(idDatFile,rcode)

      else
         idDatFile = ncopn(outname, NCWRITE,rcode)
c         call ncredf(idDatFile,rcode)
c         idvx = NCVID(idDatFile, 'xc', RCODE)
c         idvy = NCVID(idDatFile, 'yc', RCODE)
c         idvz = NCVID(idDatFile, 'zc', RCODE)
c         idvh = NCVID(idDatFile, 'h', RCODE)
c=         idvs = NCVID(idDatFile, 's', RCODE)
c=         idvt = NCVID(idDatFile, 'temp', RCODE)
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

      start2d(1)= 1
      start2d(2)= 1

      start(1)= 1
      start(2)= 1
      start(3)= step/frame_int +1

      if (step.eq.0) then
         CALL ncvpt(idDatFile,idvx, start(1), count(1), xsurf, rcode)
         CALL ncvpt(idDatFile,idvy, start(2), count(2), ysurf, rcode)
         CALL ncvpt(idDatFile,idvz, start2d, count2d, zsurf, rcode)
      endif
c      CALL ncvpt(idDatFile,idvh, start, count, hsurf, rcode)
c=      CALL ncvpt(idDatFile,idvs, start, count, ssurf, rcode)
c=      CALL ncvpt(idDatFile,idvT, start, count, Tsurf, rcode)
      CALL ncvpt(idDatFile,idvrho, start, count, rsurf, rcode)
      CALL ncvpt(idDatFile,idvu, start, count, usurf, rcode)
      CALL ncvpt(idDatFile,idvv, start, count, vsurf, rcode)
      CALL ncvpt(idDatFile,idvw, start, count, wsurf, rcode)
      CALL ncvpt(idDatFile,idvvor, start, count, vorsurf, rcode)

      CALL ncclos(idDatFile,rcode)

      end
