      subroutine outcdf(xc,yc,zc,p,h,s,T,u,v,w,vor,uf,vf,wf,
     &     theta,xg,zg,ght,gdef,step)
c----------------------------------------------------------
c     outputs the solution as a netCDF file
c      implicit logical (a-z)
c      include 'header.f'
      include 'dims.f'
      integer i,j,k,step,nstp,idigit,ipos,it
      double precision rho(0:NI+1,0:NJ+1,0:NK+1),
     &     p(0:NI+1,0:NJ+1,0:NK+1),u(0:NI+1,0:NJ+1,0:NK+1),
     &     v(0:NI+1,0:NJ+1,0:NK+1),w(0:NI+1,0:NJ+1,0:NK+1),
     &     vor(0:NI+1,0:NJ+1,0:NK+1),
     &     s(0:NI+1,0:NJ+1,0:NK+1),T(ntr,0:NI+1,0:NJ+1,0:NK+1),
     &     Tr(0:NI+1,0:NJ+1,0:NK+1,ntr),
     &     h(0:NI+1,0:NJ+1),uf(0:NI,NJ,NK),
     &     vf(NI,0:NJ,NK),wf(NI,NJ,0:NK)
      double precision ght(0:NI+1,0:NJ+1),gdef(0:NI+1,0:NJ+1),
     &     xg(0:NI+1,0:NJ+1,NKg),zg(0:NI+1,0:NJ+1,NKg),
     &     theta(0:NI+1,0:NJ+1,NKg),thedeg(0:NI+1,0:NJ+1,NKg)
      double precision xc(0:NI+1),yc(0:NJ+1),
     &     zc(0:NI+1,0:NJ+1,0:NK+1)
c      real*8 x(0:NI+1,0:NJ+1,0:NK+1),y(0:NI+1,0:NJ+1,0:NK+1),
c      real*8 x(0:NI+1,0:NJ+1),y(0:NI+1,0:NJ+1),
       real*8 z(0:NI+1,0:NJ+1,0:NK+1)
      real*8 smax,Tmax,umax,vmax,wmax,smin,Tmin,umin,vmin,wmin
      double precision PI,const
      character  outname*12, facename*12
c      include
c     &     '/home/scratch/vayu/public/netcdf-3.5.0/include/netcdf.inc'
c      INCLUDE '/usr/local/include/netcdf.inc'

      INCLUDE '/home/salvadorvieira.g/a/include/netcdf.inc'
C      INCLUDE '../netcdf-fortran-4.5.2_install/include/netcdf.inc'

c     INCLUDE '/gaeaopt/netcdf-3.4atd/include/netcdf.inc'

      integer start(3), count(3), start2d(2), count2d(2),count2dtr(3),
     &     countuf(3), countvf(3), countwf(3),count4(4),start4(4),
     &     start2dtr(3)
      integer dims(3),dims2d(2),dims2dtr(3),dimuf(3), dimvf(3),
     &     dimwf(3),dims4(4),dimg(3),countg(3)

      parameter (PI=3.141592653589793)
      DATA start /1, 1, 1/
      DATA start4 /1, 1, 1, 1/
      DATA start2d /1, 1/
      DATA start2dtr /1, 1, 1/

      const= 180.d0/PI
      do k=1,NKg
      do j=0,NJ+1
         do i=0,NI+1
            thedeg(i,j,k)= theta(i,j,k)*const
c            x(i,j)= xc(i)
c            y(i,j)= yc(j)
         end do
      end do
      end do
      do k=0,NK+1
         do j=0,NJ+1
            do i=0,NI+1
c               x(i,j,k)= xc(i)
c               y(i,j,k)= yc(j)
               z(i,j,k)= zc(i,j,k)
               do it=1,ntr
                  Tr(i,j,k,it)= T(it,i,j,k)
               end do
            end do
         end do
      end do

c     FOR THE SAKE OF BETTER PLOTS
      do j=0,NJ+1
         do i=0,NI+1
            z(i,j,NK+1)= 0.d0
            z(i,j,0)= 0.5*(z(i,j,0)+z(i,j,1))
            s(i,j,NK+1)= s(i,j,NK)
            s(i,j,0)= s(i,j,1)
            u(i,j,NK+1)= u(i,j,NK)
            u(i,j,0)= u(i,j,1)
            v(i,j,NK+1)= v(i,j,NK)
            v(i,j,0)= v(i,j,1)
            w(i,j,NK+1)= w(i,j,NK)
            w(i,j,0)= w(i,j,1)
            vor(i,j,NK+1)= vor(i,j,NK)
            vor(i,j,0)= vor(i,j,1)
            do it=1,ntr
               Tr(i,j,NK+1,it)= Tr(i,j,NK,it)
               Tr(i,j,0,it)= Tr(i,j,1,it)
            end do
         end do
      end do

      count(1)= NI+2
      count(2)= NJ+2
      count(3)= NK+2

      countg(1)= NI+2
      countg(2)= NJ+2
      countg(3)= NKg

      count4(1)= NI+2
      count4(2)= NJ+2
      count4(3)= NK+2
      count4(4)= ntr

      count2d(1)= NI+2
      count2d(2)= NJ+2

      count2dtr(1)= NI+2
      count2dtr(2)= NJ+2
      count2dtr(3)= ntr

      countuf(1)= NI+1
      countuf(2)= NJ
      countuf(3)= NK

      countvf(1)= NI
      countvf(2)= NJ+1
      countvf(3)= NK

      countwf(1)= NI
      countwf(2)= NJ
      countwf(3)= NK+1

      call evalrho(rho,0)
      do k=0,NK+1
         do j=0,NJ+1
            do i=0,NI+1
               rho(i,j,k)= rho(i,j,k) +R0 -1000.d0
            end do
         end do
      end do

c     name the output file
c
      outname = 'outxxxxx.cdf'
c-      facename = 'vfcxxxxx.cdf'
      facename = 'vfc.cdf'
      nstp = step
      do 10 ipos = 8, 4, -1
	     idigit = mod(nstp,10)
	     outname(ipos:ipos) = char(ichar('0') + idigit)
c-	      facename(ipos:ipos) = char(ichar('0') + idigit)
	     nstp = nstp / 10
 10   continue
c
c      write(6,*) 'in outcdf, outname',outname
c      write(6,*) 'in outcdf, facename',facename

      idDatFile =  nccre(outname,NCCLOB,rcode)

      dims(1) = ncddef(idDatFile,'x',NI+2,rcode)
      dims(2) = ncddef(idDatFile,'y',NJ+2,rcode)
      dims(3) = ncddef(idDatFile,'sigma',NK+2,rcode)

      dims4(1) = dims(1)
      dims4(2) = dims(2)
      dims4(3) = dims(3)
      dims4(4) = ncddef(idDatFile,'ntr',ntr,rcode)
c
      dims2d(1)= dims(1)
      dims2d(2)= dims(2)
c
      dims2dtr(1)= dims(1)
      dims2dtr(2)= dims(2)
      dims2dtr(3)= dims4(4)
c
      dimg(1) = dims(1)
      dimg(2) = dims(2)
      dimg(3) = ncddef(idDatFile,'sig',NKg,rcode)

      idvx = ncvdef(idDatFile,'xc',NCDOUBLE,1,dims(1),rcode)
      idvy = ncvdef(idDatFile,'yc',NCDOUBLE,1,dims(2),rcode)
      idvz = ncvdef(idDatFile,'zc',NCDOUBLE,3,dims,rcode)
      idvh = ncvdef(idDatFile,'h',NCDOUBLE,2,dims2d,rcode)
      idvght = ncvdef(idDatFile,'ght',NCDOUBLE,2,dims2d,rcode)
      idvgdef = ncvdef(idDatFile,'gdef',NCDOUBLE,2,dims2d,rcode)
c      idvhdt = ncvdef(idDatFile,'consump',NCDOUBLE,3,dims2dtr,rcode)
      idvs = ncvdef(idDatFile,'s',NCDOUBLE,3,dims,rcode)
      idvt = ncvdef(idDatFile,'tr',NCDOUBLE,4,dims4,rcode)
      idvrho = ncvdef(idDatFile,'rho',NCDOUBLE,3,dims,rcode)
      idvp = ncvdef(idDatFile,'p',NCDOUBLE,3,dims,rcode)
      idvu = ncvdef(idDatFile,'u',NCDOUBLE,3,dims,rcode)
      idvv = ncvdef(idDatFile,'v',NCDOUBLE,3,dims,rcode)
      idvw = ncvdef(idDatFile,'w',NCDOUBLE,3,dims,rcode)
      idvz3 = ncvdef(idDatFile,'vor',NCDOUBLE,3,dims,rcode)
      idvth = ncvdef(idDatFile,'theta',NCDOUBLE,3,dimg,rcode)
      idvxg = ncvdef(idDatFile,'xg',NCDOUBLE,3,dimg,rcode)
      idvzg = ncvdef(idDatFile,'zg',NCDOUBLE,3,dimg,rcode)

      CALL ncendf(idDatFile,rcode)

      CALL ncvpt(idDatFile,idvx, start(1), count(1), xc, rcode)
      CALL ncvpt(idDatFile,idvy, start(2), count(2), yc, rcode)
      CALL ncvpt(idDatFile,idvz, start, count, z, rcode)
      CALL ncvpt(idDatFile,idvh, start2d, count2d, h, rcode)
      CALL ncvpt(idDatFile,idvght,start2d,count2d,ght,rcode)
      CALL ncvpt(idDatFile,idvgdef,start2d,count2d,gdef,rcode)
c      CALL ncvpt(idDatFile,idvhdt,start2dtr,count2dtr,consump,rcode)
      CALL ncvpt(idDatFile,idvs, start, count, s, rcode)
      CALL ncvpt(idDatFile,idvt, start4, count4, Tr, rcode)
      CALL ncvpt(idDatFile,idvrho, start, count, rho, rcode)
      CALL ncvpt(idDatFile,idvp, start, count, p, rcode)
      CALL ncvpt(idDatFile,idvu, start, count, u, rcode)
      CALL ncvpt(idDatFile,idvv, start, count, v, rcode)
      CALL ncvpt(idDatFile,idvw, start, count, w, rcode)
      CALL ncvpt(idDatFile,idvz3, start, count, vor, rcode)
      CALL ncvpt(idDatFile,idvth, start, countg, thedeg, rcode)
      CALL ncvpt(idDatFile,idvxg, start, countg, xg, rcode)
      CALL ncvpt(idDatFile,idvzg, start, countg, zg, rcode)
c
      CALL ncclos(idDatFile,rcode)
c
c     ----------------------------------------------------------------
c     write face velocities, uf,vf,wf
c     -------------------------------
c      if (mod(step,3000).ne.0) return
c
      idFaceFile =  nccre(facename,NCCLOB,rcode)

      dimuf(1) = ncddef(idFaceFile,'xi-ew',NI+1,rcode)
      dimuf(2) = ncddef(idFaceFile,'eta-ew',NJ,rcode)
      dimuf(3) = ncddef(idFaceFile,'sigma-ew',NK,rcode)
c
      dimvf(1) = ncddef(idFaceFile,'xi-ns',NI,rcode)
      dimvf(2) = ncddef(idFaceFile,'eta-ns',NJ+1,rcode)
      dimvf(3) = ncddef(idFaceFile,'sigma-ns',NK,rcode)
c
      dimwf(1) = ncddef(idFaceFile,'xi-tb',NI,rcode)
      dimwf(2) = ncddef(idFaceFile,'eta-tb',NJ,rcode)
      dimwf(3) = ncddef(idFaceFile,'sigma-tb',NK+1,rcode)
c
      idvuf = ncvdef(idFaceFile,'uf',NCDOUBLE,3,dimuf,rcode)
      idvvf = ncvdef(idFaceFile,'vf',NCDOUBLE,3,dimvf,rcode)
      idvwf = ncvdef(idFaceFile,'wf',NCDOUBLE,3,dimwf,rcode)
c      idvKf = ncvdef(idFaceFile,'Kz',NCDOUBLE,3,dimwf,rcode)

      CALL ncendf(idFaceFile,rcode)
c
      CALL ncvpt(idFaceFile,idvuf, start, countuf, uf, rcode)
      CALL ncvpt(idFaceFile,idvvf, start, countvf, vf, rcode)
      CALL ncvpt(idFaceFile,idvwf, start, countwf, wf, rcode)
c      CALL ncvpt(idFaceFile,idvKf, start, countwf, Kz, rcode)
c
      CALL ncclos(idFaceFile,rcode)

      return
      end
