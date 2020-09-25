      subroutine writeslice(frame_int,step,s,T,rho,u,v,w,vor,yc,zc)
c     ------------------------------------------------------------
c     reads in a netcdf file and interpolates the data from the sigma grid
c     onto a z level

C      character*100 aux
C      parameter(aux="/home/salvadorvieira.g/PSOM/netcdf-fortran"//
C     & "-4.5.2_install/include/netcdf.inc")
C      parameter(aux="/home/salvadorvieira.g/a/"//
C     & "include/netcdf.inc")
C      parameter (aux="/home/salvadorvieira.g/a/include/netcdf.inc")
C      INCLUDE : aux

      INCLUDE '/home/salvadorvieira.g/a/include/netcdf.inc'
C      INCLUDE '../netcdf-fortran-4.5.2_install/include/netcdf.inc'

C      INCLUDE '/home/salvadorvieira.g/PSOM/netcdf-fortran-4.5.2_install/include/netcdf.inc'


c      INCLUDE '/usr/local/include/netcdf.inc'

c      INCLUDE '/gaeaopt/netcdf-3.4atd/include/netcdf.inc'
c      INCLUDE
c     &     '/home/scratch/vayu/public/netcdf-3.5.0/include/netcdf.inc'
      include 'dims.f'
      integer NI2,NJ2,NK2,i,j,k,step,frame_int,itr
      parameter ( NI2=NI+2,NJ2=NJ+2,NK2=NK+2)
      double precision DL
      parameter (DL= 1000.d0)
      double precision svert(NK),Tvert(NK),rvert(NK),uvert(NK),
     &     vvert(NK),wvert(NK),zvert(NK),seval

      double precision
     &     rho(0:NI+1,0:NJ+1,0:NK+1),u(0:NI+1,0:NJ+1,0:NK+1),
     &     v(0:NI+1,0:NJ+1,0:NK+1),w(0:NI+1,0:NJ+1,0:NK+1),
     &     s(0:NI+1,0:NJ+1,0:NK+1),T(ntr,0:NI+1,0:NJ+1,0:NK+1),
     &     vor(0:NI+1,0:NJ+1,0:NK+1),
     &     zc(0:NI+1,0:NJ+1,0:NK+1),
     &     yc(0:NJ+1),bs1(NK),cs1(NK),ds1(NK)
      real*4 yslice(NJ),zslice(NJ,0:NK+1),sslice(NJ,0:NK+1),
     &     Tslice(ntr,NJ,0:NK+1),rslice(NJ,0:NK+1),uslice(NJ,0:NK+1),
     &     vslice(NJ,0:NK+1),wslice(NJ,0:NK+1),vorslice(NJ,0:NK+1)
c

      character*80 outname

      integer start(3),count(3),dims(3),start2d(2),dims4(4),start4(4),
     &     count4(4)


      i= NI/2
      do k=1,NK
         do j=1,NJ
            zslice(j,k)= zc(i,j,k)*DL
            yslice(j)= yc(j)

c=            sslice(j,k)= s(i,j,k) +S0
c=            Tslice(j,k)= T(i,j,k) +T0
            rslice(j,k)= rho(i,j,k) +R0 -1000.d0
            uslice(j,k)= u(i,j,k)
            vslice(j,k)= v(i,j,k)
            wslice(j,k)= w(i,j,k)
            vorslice(j,k)= vor(i,j,k)
            do itr=1,ntr
               Tslice(itr,j,k)= T(itr,i,j,k)
            end do
         end do
      end do
c     FOR THE SAKE OF BETTER PLOTS
      do j=1,NJ
         zslice(j,NK+1)= 0.d0
         zslice(j,0)= 0.5*DL*(zc(i,j,0)+zc(i,j,1))
         rslice(j,NK+1)= rslice(j,NK)
         rslice(j,0)= rslice(j,1)
         uslice(j,NK+1)= uslice(j,NK)
         uslice(j,0)= uslice(j,1)
         vslice(j,NK+1)= vslice(j,NK)
         vslice(j,0)= vslice(j,1)
         wslice(j,NK+1)= wslice(j,NK)
         wslice(j,0)= wslice(j,1)
         vorslice(j,NK+1)= vorslice(j,NK)
         vorslice(j,0)= vorslice(j,1)
         do itr=1,ntr
            Tslice(itr,j,NK+1)= Tslice(itr,j,NK)
            Tslice(itr,j,0)= Tslice(itr,j,1)
         end do
      end do

c     write to a netcdf file

      outname= 'slicemovie.cdf'

      if (step.eq.0) then
         idDatFile =  nccre(outname,NCCLOB,rcode)

         dims(1) = ncddef(idDatFile,'y',NJ,rcode)
         dims(2) = ncddef(idDatFile,'z',NK+2,rcode)
         dims(3) = ncddef(idDatFile,'time',NCUNLIM,rcode)

         dims4(1) = ncddef(idDatFile,'ntr',ntr,rcode)
         dims4(2) = dims(1)
         dims4(3) = dims(2)
         dims4(4) = dims(3)

         idvy = ncvdef(idDatFile,'yc',NCFLOAT,1,dims(1),rcode)
         idvz = ncvdef(idDatFile,'zc',NCFLOAT,2,dims,rcode)
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
c=         idvs = NCVID(idDatFile, 's', RCODE)
c=         idvt = NCVID(idDatFile, 'temp', RCODE)
         idvt = NCVID(idDatFile, 'tr', RCODE)
         idvrho = NCVID(idDatFile, 'rho', RCODE)
         idvu = NCVID(idDatFile, 'u', RCODE)
         idvv = NCVID(idDatFile, 'v', RCODE)
         idvw = NCVID(idDatFile, 'w', RCODE)
         idvvor = NCVID(idDatFile, 'vor', RCODE)
      endif

      count(1)= NJ
      count(2)= NK+2
      count(3)= 1

      count4(1)= ntr
      count4(2)= NJ
      count4(3)= NK+2
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
         CALL ncvpt(idDatFile,idvy, start(1), count(1), yslice, rcode)
         CALL ncvpt(idDatFile,idvz, start, count, zslice, rcode)
      endif
c=      CALL ncvpt(idDatFile,idvs, start, count, sslice, rcode)
      CALL ncvpt(idDatFile,idvt, start4, count4, Tslice, rcode)
      CALL ncvpt(idDatFile,idvrho, start, count, rslice, rcode)
      CALL ncvpt(idDatFile,idvu, start, count, uslice, rcode)
      CALL ncvpt(idDatFile,idvv, start, count, vslice, rcode)
      CALL ncvpt(idDatFile,idvw, start, count, wslice, rcode)
      CALL ncvpt(idDatFile,idvvor, start, count, vorslice, rcode)

      CALL ncclos(idDatFile,rcode)

      end
