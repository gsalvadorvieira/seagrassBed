      subroutine writeyslice(frame_int,step,s,T,rho,u,v,w,p,vor,h,
     &     ght,gdef,xg,zg,xc,zc,pgrad,hgrad,bgrad,theta)
c     ------------------------------------------------------------
c     reads in a netcdf file and interpolates the data from the sigma grid
c     onto a z level

C       character*100 fname
c      INCLUDE '/usr/local/include/netcdf.inc'

      INCLUDE "/home/salvadorvieira.g/a/include/netcdf.inc"

C      INCLUDE '../netcdf-fortran-4.5.2_install/include/netcdfinc'

C       parameter (foldername = "/home/salvadorvieira.g/" //
C      &   "netcdf-fortran-4.5.2_install/include/netcdf.inc")
C       INCLUDE : foldername

C       parameter(fname='/home/salvadorvieira.g/a/include/netcdf.inc')
C       include: fname

c      INCLUDE '/gaeaopt/netcdf-3.4atd/include/netcdf.inc'
c      INCLUDE
c     &     '/home/scratch/vayu/public/netcdf-3.5.0/include/netcdf.inc'

      include 'dims.f'
      integer NI2,NJ2,NK2,i,j,k,step,frame_int
      parameter (NI2=NI+2,NJ2=NJ+2,NK2=NK+2)
      double precision PI,const
      parameter (PI=3.141592653589793)
      double precision DL
      parameter (DL= 1.d0,LEN=1.d0)
      double precision svert(NK),Tvert(NK),rvert(NK),uvert(NK),
     &     vvert(NK),wvert(NK),zvert(NK),seval
      double precision ght(0:NI+1,0:NJ+1),gdef(0:NI+1,0:NJ+1),
     &     xg(0:NI+1,0:NJ+1,NKg),zg(0:NI+1,0:NJ+1,NKg),
     &     h(0:NI+1,0:NJ+1)
      double precision pgrad(0:NI+1,0:NJ+1,0:NK+1),
     &     hgrad(0:NI+1,0:NJ+1,0:NK+1),
     &     bgrad(0:NI+1,0:NJ+1,0:NK+1),
     &     theta(0:NI+1,0:NJ+1,NKg)

      double precision
     &     rho(0:NI+1,0:NJ+1,0:NK+1,0:1),u(0:NI+1,0:NJ+1,0:NK+1,0:1),
     &     v(0:NI+1,0:NJ+1,0:NK+1,0:1),w(0:NI+1,0:NJ+1,0:NK+1,0:1),
     &     s(0:NI+1,0:NJ+1,0:NK+1,0:1),T(1,0:NI+1,0:NJ+1,0:NK+1,0:1),
     &     vor(0:NI+1,0:NJ+1,0:NK+1),p(0:NI+1,0:NJ+1,0:NK+1,0:1),
     &     zc(0:NI+1,0:NJ+1,0:NK+1),
C     &     xc(0:NJ+1),bs1(NK),cs1(NK),ds1(NK)
     &     xc(0:NI+1),bs1(NK),cs1(NK),ds1(NK)

      real*4 xslice(NI),zslice(NI,0:NK+1),sslice(NI,0:NK+1),
     &     Tslice(NI,0:NK+1),rslice(NI,0:NK+1),uslice(NI,0:NK+1),
     &     vslice(NI,0:NK+1),wslice(NI,0:NK+1),vorslice(NI,0:NK+1),
     &     pslice(NI,0:NK+1),pgslice(NI,0:NK+1),hgslice(NI,0:NK+1),
     &     bgslice(NI,0:NK+1),
     &     xgslice(NI,NKg),zgslice(NI,NKg),ghtslice(NI),gdefslice(NI),
     &     hslice(NI),theslice(NI,NKg),dsslice(NI,NKg)
c

      character*80 outname

      integer start(3),count(3),countg(3),count1d(2),start1d(2),
     &     dims(3),dimg(3),dim1d(2),dim2d(2)


      j= NJ/2 +1
      do k=1,NK
         do i=1,NI
            zslice(i,k)= zc(i,j,k)*DL
            xslice(i)= xc(i)

c=           sslice(j,k)= s(i,j,k) +S0
C            Tslice(j,k)= T(i,j,k) +T0

            Tslice(i,k) = T(1,i,j,k,0) ! using only first tracer
            rslice(i,k)= rho(i,j,k,0) + R0 - 1000.d0
            uslice(i,k)= u(i,j,k,0)
            vslice(i,k)= v(i,j,k,0)
            wslice(i,k)= w(i,j,k,0)
            pslice(i,k)= p(i,j,k,0)
            pgslice(i,k)= pgrad(i,j,k)
            hgslice(i,k)= hgrad(i,j,k)
            bgslice(i,k)= bgrad(i,j,k)
            vorslice(i,k)= vor(i,j,k)

         end do
      end do

      const= 180.d0/PI
      do i=1,NI
         ghtslice(i)= ght(i,j)
         gdefslice(i)= gdef(i,j)
         hslice(i)= h(i,j)
         do k=1,NKg
            xgslice(i,k)= xg(i,j,k)*LEN
            zgslice(i,k)= zg(i,j,k)*DL
            theslice(i,k)= theta(i,j,k)*const
C	     write(6,*) 'step,i,k,xg,zg,theta',step,i,k,
C     &    	xgslice(i,k),zgslice(i,k),theslice(i,k)
	       end do
         dsslice(i,1)= 0.0
         do k=2,NKg
            dsslice(i,k)= sqrt((xgslice(i,k)-xgslice(i,k-1))**2
     &           + (zgslice(i,k)-zgslice(i,k-1))**2 )
         end do
      end do

c     FOR THE SAKE OF BETTER PLOTS
      do i=1,NI
         zslice(i,NK+1)= 0.d0
         zslice(i,0)= 0.5*DL*(zc(i,j,0)+zc(i,j,1))
         rslice(i,NK+1)= rslice(i,NK)
         rslice(i,0)= rslice(i,1)
         uslice(i,NK+1)= uslice(i,NK)
         uslice(i,0)= 2.0*uslice(i,1) -uslice(i,2)
         vslice(i,NK+1)= 2.0*vslice(i,NK) -vslice(i,NK-1)
         vslice(i,0)= 2.0*vslice(i,1) -vslice(i,2)
         wslice(i,NK+1)= 2.0*wslice(i,NK) -wslice(i,NK-1)
         wslice(i,0)= 0.0
         pgslice(i,NK+1)= 2.0*pgslice(i,NK) -pgslice(i,NK-1)
         pgslice(i,0)= 2.0*pgslice(i,1) -pgslice(i,2)
         hgslice(i,NK+1)= 2.0*hgslice(i,NK) -hgslice(i,NK-1)
         hgslice(i,0)= 2.0*hgslice(i,1) -hgslice(i,2)
         bgslice(i,NK+1)= 2.0*bgslice(i,NK) -bgslice(i,NK-1)
         bgslice(i,0)= 2.0*bgslice(i,1) -bgslice(i,2)
         pslice(i,NK+1)= -pslice(i,NK)
         pslice(i,0)= 2.0*pslice(i,1) -pslice(i,2)
         vorslice(i,NK+1)= 2.0*vorslice(i,NK) -vorslice(i,NK-1)
         vorslice(i,0)= 2.0*vorslice(i,1) -vorslice(i,2)

         Tslice(i,NK+1)= Tslice(i,NK)
         Tslice(i,0)= Tslice(i,1)

      end do

c     write to a netcdf file

      outname= 'yslicemovie.cdf'

      if (step.eq.0) then
         idDatFile =  nccre(outname,NCCLOB,rcode)

         dims(1) = ncddef(idDatFile,'x',NI,rcode)
         dims(2) = ncddef(idDatFile,'z',NK+2,rcode)
         dims(3) = ncddef(idDatFile,'time',NCUNLIM,rcode)

         dimg(1)= dims(1)
         dimg(2)= ncddef(idDatFile,'sig',NKg,rcode)
         dimg(3) = dims(3)

         dim1d(1)= dims(1)
         dim1d(2) = dims(3)

         dim2d(1) = dims(1)
         dim2d(2) = dims(2)

         idvy = ncvdef(idDatFile,'xc',NCFLOAT,1,dims(1),rcode)
C         idvz = ncvdef(idDatFile,'zc',NCFLOAT,2,dims,rcode)
         idvz = ncvdef(idDatFile,'zc',NCFLOAT,3,dims,rcode)
c=         idvs = ncvdef(idDatFile,'s',NCFLOAT,3,dims,rcode)
         idvt = ncvdef(idDatFile,'tr',NCFLOAT,3,dims,rcode)
         idvh = ncvdef(idDatFile,'h',NCFLOAT,2,dim1d,rcode)
         idvght = ncvdef(idDatFile,'ght',NCFLOAT,2,dim1d,rcode)
         idvgdef = ncvdef(idDatFile,'gdef',NCFLOAT,2,dim1d,rcode)
         idvrho = ncvdef(idDatFile,'rho',NCFLOAT,3,dims,rcode)
         idvu = ncvdef(idDatFile,'u',NCFLOAT,3,dims,rcode)
         idvv = ncvdef(idDatFile,'v',NCFLOAT,3,dims,rcode)
         idvw = ncvdef(idDatFile,'w',NCFLOAT,3,dims,rcode)
         idvp = ncvdef(idDatFile,'p',NCFLOAT,3,dims,rcode)
         idvpg = ncvdef(idDatFile,'pgrad',NCFLOAT,3,dims,rcode)
         idvhg = ncvdef(idDatFile,'hgrad',NCFLOAT,3,dims,rcode)
         idvbg = ncvdef(idDatFile,'bgrad',NCFLOAT,3,dims,rcode)
         idvvor = ncvdef(idDatFile,'vor',NCFLOAT,3,dims,rcode)
         idvxg = ncvdef(idDatFile,'xg',NCFLOAT,3,dimg,rcode)
         idvzg = ncvdef(idDatFile,'zg',NCFLOAT,3,dimg,rcode)
         idvthe = ncvdef(idDatFile,'theta',NCFLOAT,3,dimg,rcode)
         idvds = ncvdef(idDatFile,'ds',NCFLOAT,3,dimg,rcode)

         CALL ncendf(idDatFile,rcode)

      else

         idDatFile = ncopn(outname, NCWRITE,rcode)
c         call ncredf(idDatFile,rcode)
c=         idvs = NCVID(idDatFile, 's', RCODE)

         idvt = NCVID(idDatFile, 'tr', RCODE)

C        added this line
         idvz = NCVID(idDatFile,'zc', RCODE)

         idvh = NCVID(idDatFile, 'h', RCODE)
         idvght = NCVID(idDatFile, 'ght', RCODE)
         idvgdef = NCVID(idDatFile, 'gdef', RCODE)
         idvrho = NCVID(idDatFile, 'rho', RCODE)
         idvu = NCVID(idDatFile, 'u', RCODE)
         idvv = NCVID(idDatFile, 'v', RCODE)
         idvw = NCVID(idDatFile, 'w', RCODE)
         idvp = NCVID(idDatFile, 'p', RCODE)
         idvpg = NCVID(idDatFile, 'pgrad', RCODE)
         idvhg = NCVID(idDatFile, 'hgrad', RCODE)
         idvbg = NCVID(idDatFile, 'bgrad', RCODE)
         idvvor = NCVID(idDatFile, 'vor', RCODE)
         idvxg = NCVID(idDatFile, 'xg', RCODE)
         idvzg = NCVID(idDatFile, 'zg', RCODE)
         idvthe = NCVID(idDatFile, 'theta', RCODE)
         idvds = NCVID(idDatFile, 'ds', RCODE)


      endif

      count(1)= NI
      count(2)= NK+2
      count(3)= 1

      countg(1)= NI
      countg(2)= NKg
      countg(3)= 1

      count1d(1)= NI
      count1d(2)= 1

      start1d(1)= 1
      start1d(2)= step/frame_int +1

      start(1)= 1
      start(2)= 1
      start(3)= step/frame_int +1

      if (step.eq.0) then
         CALL ncvpt(idDatFile,idvy, start(1), count(1), xslice, rcode)
         CALL ncvpt(idDatFile,idvz, start, count, zslice, rcode)
      endif

C      Modified to account for deformation in z grid (zc)
      CALL ncvpt(idDatFile,idvz, start, count, zslice, rcode)

c=      CALL ncvpt(idDatFile,idvs, start, count, sslice, rcode)
      CALL ncvpt(idDatFile,idvt, start, count, Tslice, rcode)
      CALL ncvpt(idDatFile,idvh, start1d, count1d, hslice, rcode)
      CALL ncvpt(idDatFile,idvght, start1d, count1d, ghtslice, rcode)
      CALL ncvpt(idDatFile,idvgdef, start1d, count1d, gdefslice, rcode)
      CALL ncvpt(idDatFile,idvrho, start, count, rslice, rcode)
      CALL ncvpt(idDatFile,idvu, start, count, uslice, rcode)
      CALL ncvpt(idDatFile,idvv, start, count, vslice, rcode)
      CALL ncvpt(idDatFile,idvw, start, count, wslice, rcode)
      CALL ncvpt(idDatFile,idvp, start, count, pslice, rcode)
      CALL ncvpt(idDatFile,idvpg, start, count, pgslice, rcode)
      CALL ncvpt(idDatFile,idvhg, start, count, hgslice, rcode)
      CALL ncvpt(idDatFile,idvbg, start, count, bgslice, rcode)
      CALL ncvpt(idDatFile,idvvor, start, count, vorslice, rcode)
      CALL ncvpt(idDatFile,idvxg, start, countg, xgslice, rcode)
      CALL ncvpt(idDatFile,idvzg, start, countg, zgslice, rcode)
      CALL ncvpt(idDatFile,idvthe, start, countg, theslice, rcode)
      CALL ncvpt(idDatFile,idvds, start, countg, dsslice, rcode)

      CALL ncclos(idDatFile,rcode)

      end
