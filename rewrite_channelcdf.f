      program rewritecdf
c
c     reads in a netCDF file (on my grid) and rewrites it
c     using depth(i,j,k) instead of  sigma(xi,eta,sigma)


      INCLUDE '/usr/local/include/netcdf.inc'

      integer NI2,NJ2,NK2
      include 'dims.f'
c      parameter ( NI=48,NJ=24,NK=16)
      parameter ( NI2=NI+2,NJ2=NJ+2,NK2=NK+2)
c
      character*80 filename,outname
c     
      integer i,j,k
      double precision  x(0:NI+1),y(0:NJ+1),
     &     z(0:NK+1),depth,delz
      double precision s(0:NI+1,0:NJ+1,0:NK+1),
     &     T(0:NI+1,0:NJ+1,0:NK+1),rho(0:NI+1,0:NJ+1,0:NK+1),
     &     u(0:NI+1,0:NJ+1,0:NK+1),p(0:NI+1,0:NJ+1,0:NK+1),
     &     v(0:NI+1,0:NJ+1,0:NK+1),w(0:NI+1,0:NJ+1,0:NK+1),
     &     h(0:NI+1,0:NJ+1)

      
      integer start(3), count(3), start2d(2), count2d(2)
      integer dims(3), dims2d(2)

      DATA start /1, 1, 1/
      DATA start2d /1, 1/
      DATA count /NI2, NJ2, NK2/
      DATA count2d /NI2, NJ2/


      open(unit=12, file='xyD.dat')
      do j=0,NJ+1
         do i=0,NI+1
            read(12,*) x(i),y(j),depth
         end do
      end do
      depth=depth*1000.
      delz=depth/dble(NK)
      do k=0,NK+1
         z(k)= -depth + (dble(k)-0.5) *delz
      end do

c     read netCDF file
      write(6,*) 'enter (within single quotes) the name of cdf file 
     &     to be converted:'
      read(5,*) filename
      write(6,*) filename
      idT = ncopn(filename, NCNOWRIT,rcode)
c

c      idvx = ncvid(idT,'x',rcode)
c      idvy = ncvid(idT,'y',rcode)
c      idvz = ncvid(idT,'z',rcode)
c      write(6,*) 'xyz ids defned'
c
      idvh = ncvid(idT,'h',rcode)
c      idvhdt = ncvid(idT,'hdt',rcode)
c      idvq2 = ncvid(idT,'q2deatn',rcode)
c      idvml = ncvid(idT,'dml',rcode)
      idvs = ncvid(idT,'s',rcode)
      idvt = ncvid(idT,'T',rcode)
      idvrho = ncvid(idT,'rho',rcode)
c-      idvq = ncvid(idT,'q',rcode)
      idvp = ncvid(idT,'p',rcode)
      idvu = ncvid(idT,'u',rcode)
      idvv = ncvid(idT,'v',rcode)
      idvw = ncvid(idT,'w',rcode)
c     
      write(6,*) 'idv defined'
c
c      call ncvgt( idT, idvx, start(1), count(1), x, rcode )
c      call ncvgt( idT, idvy, start(2), count(2), y, rcode )
c      call ncvgt( idT, idvz, start, count, z, rcode )
c      write(6,*) 'grid read'
c
      call ncvgt( idT, idvh, start2d, count2d, h, rcode )
c      call ncvgt( idT, idvq2, start2d, count2d, q2deatn, rcode )
c      call ncvgt( idT, idvml, start2d, count2d, dml, rcode )
      call ncvgt( idT, idvs, start, count, s, rcode )
      call ncvgt( idT, idvt, start, count, T, rcode )
      call ncvgt( idT, idvrho, start, count, rho, rcode )
c      call ncvgt( idT, idvq, start, count, q, rcode )
      call ncvgt( idT, idvp, start, count, p, rcode )
      call ncvgt( idT, idvu, start, count, u, rcode )
      call ncvgt( idT, idvv, start, count, v, rcode )
      call ncvgt( idT, idvw, start, count, w, rcode )
      write(6,*) 'data read'
      
      call ncclos(idT, rcode)

      write(6,*) 'old cdf file read'
c
c------------------------------------------------------------
c     write new cdf
      outname = 'newsemt02000.cdf'
      write(6,*) 'now create file',outname

      idDatFile =  nccre(outname,NCCLOB,rcode)      

      dims(1) = ncddef(idDatFile,'long',NI+2,rcode)
      dims(2) = ncddef(idDatFile,'lat',NJ+2,rcode)
      dims(3) = ncddef(idDatFile,'k',NK+2,rcode)
c
      dims2d(1)= dims(1)
      dims2d(2)= dims(2)
c
      idvx = ncvdef(idDatFile,'long',NCDOUBLE,1,dims(1),rcode)
      idvy = ncvdef(idDatFile,'lat',NCDOUBLE,1,dims(2),rcode)
      idvz = ncvdef(idDatFile,'depth',NCDOUBLE,1,dims(3),rcode)
      idvh = ncvdef(idDatFile,'h',NCDOUBLE,2,dims2d,rcode)
c      idvhdt = ncvdef(idDatFile,'hdt',NCDOUBLE,2,dims2d,rcode)
c-      idvq2 = ncvdef(idDatFile,'q2deatn',NCDOUBLE,2,dims2d,rcode)
c-      idvml = ncvdef(idDatFile,'dml',NCDOUBLE,2,dims2d,rcode)
      idvs = ncvdef(idDatFile,'sal',NCDOUBLE,3,dims,rcode)
      idvt = ncvdef(idDatFile,'temp',NCDOUBLE,3,dims,rcode)
      idvrho = ncvdef(idDatFile,'rho',NCDOUBLE,3,dims,rcode)
c      idvq = ncvdef(idDatFile,'q',NCDOUBLE,3,dims,rcode)
c      idvp = ncvdef(idDatFile,'p',NCDOUBLE,3,dims,rcode)
      idvu = ncvdef(idDatFile,'u',NCDOUBLE,3,dims,rcode)
      idvv = ncvdef(idDatFile,'v',NCDOUBLE,3,dims,rcode)
      idvw = ncvdef(idDatFile,'w',NCDOUBLE,3,dims,rcode)

c      call ncapt(idDatFile,idvx,'deg E',NCCHAR,5,'example',rcode)
c      call ncapt(idDatFile,idvy,'deg N',NCCHAR,5,'example',rcode)
      CALL ncendf(idDatFile,rcode)

      CALL ncvpt(idDatFile,idvx, start(1), count(1), x, rcode)
      CALL ncvpt(idDatFile,idvy, start(2), count(2), y, rcode)
      CALL ncvpt(idDatFile,idvz, start(3), count(3), z, rcode)
      CALL ncvpt(idDatFile,idvh, start2d, count2d, h, rcode)
c      CALL ncvpt(idDatFile,idvhdt, start2d, count2d, hdt, rcode)
c-      CALL ncvpt(idDatFile,idvq2, start2d, count2d, q2deatn, rcode)
c-      CALL ncvpt(idDatFile,idvml, start2d, count2d, dml, rcode)
      CALL ncvpt(idDatFile,idvs, start, count, s, rcode)
      CALL ncvpt(idDatFile,idvt, start, count, T, rcode)
      CALL ncvpt(idDatFile,idvrho, start, count, rho, rcode)
c      CALL ncvpt(idDatFile,idvq, start, count, q, rcode)
c      CALL ncvpt(idDatFile,idvp, start, count, p, rcode)
      CALL ncvpt(idDatFile,idvu, start, count, u, rcode)
      CALL ncvpt(idDatFile,idvv, start, count, v, rcode)
      CALL ncvpt(idDatFile,idvw, start, count, w, rcode)
c
      CALL ncclos(idDatFile,rcode)
c
c
      end


