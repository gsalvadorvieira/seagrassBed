      integer NI,NJ,NK,NKg,maxout,ntr

C      NKg = number of points in the grassy bed

C     AUTOMATED VALUES
      parameter(NI=432,NJ=8,NK=48,NKg=24,ntr=1)
      parameter(maxout=257168)

C      parameter(NI=288,NJ=8,NK=48,NKg=24,ntr=1)
C      parameter(maxout=171920)

c      parameter(NI=96,NJ=8,NK=56,NKg=24,ntr=1)
c      parameter( maxout=67504)
c      parameter(NI=96,NJ=8,NK=52,NKg=24,ntr=1)
c      parameter( maxout=62880)
c      parameter(NI=96,NJ=8,NK=60,NKg=24,ntr=1)
c      parameter( maxout=72128)
cshall      parameter(NI=96,NJ=8,NK=36,NKg=24,ntr=1)
cshall      parameter( maxout=44384)
c      parameter(NI=96,NJ=8,NK=72,NKg=24,ntr=1)
c      parameter( maxout=86000)

cwide      parameter(NI=96,NJ=48,NK=48,NKg=24,ntr=1)
cwide      parameter( maxout=284992)
c=      parameter(NI=96,NJ=8,NK=48,NKg=36,ntr=1)
clong      parameter(NI=192,NJ=8,NK=48,NKg=24,ntr=1)
clong      parameter( maxout=115088)

      real*8 S0,T0,R0
      parameter (S0=35.7d0, T0=18.d0 ,R0=1000.d0)
c      parameter (S0=35.7d0, T0=15.d0 ,R0=1027.d0)
