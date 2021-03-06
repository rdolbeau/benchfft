      SUBROUTINE MFFTP4(EXPTAB,SPETAB,FACTAB,N,N1)
*
*     THIS SUBROUTINE BUILDS THE TWIDDLE FACTOR TABLES FOR USE
*     OF MFFT?S ROUTINES (SPECIAL OPTIMIZATION FOR SMALL
*     DATA MATRICES); IT MUST BE USED BY MFFTP ONLY.
*
*     PARAMETERS:
*     EXPTAB: TWIDDLE FACTOR TABLE
*     SPETAB: SPECIAL TWIDDLE FACTOR TABLE
*     FACTAB: FACTORIZATION OF N
*     N: ORDER OF THE TRANSFORM
*     N1: FIRST DIMENSION OF THE ARRAY TO BE TRANSFORMED (.GE.N)
*
      INTEGER MAXFAC
      PARAMETER(MAXFAC=3)
      INTEGER FACTAB(14)
      COMPLEX EXPTAB(0:*),SPETAB(0:*)
      INTEGER FACTOR(MAXFAC)
      DATA FACTOR/2,3,5/
 
      MX=N
      LX=1
      J=0
      DO 50 IFACT=FACTAB(14),1,-1
        NX=FACTOR(IFACT)
        DO 40 IPOW=1,FACTAB(IFACT)
          MX=MX/NX
          DO 30 NU=1,NX-1
            DO 20 MU=0,MX-1
              DO 10 I1=0,N1-1
                SPETAB(J)=CONJG(EXPTAB(MU*LX*NU))
                J=J+1
10            CONTINUE
20          CONTINUE
30        CONTINUE
          LX=LX*NX
40      CONTINUE
50    CONTINUE
      MX=N
      LX=1
      J=N*N1
      DO 100 IFACT=1,FACTAB(14)
        NX=FACTOR(IFACT)
        DO 90 IPOW=1,FACTAB(IFACT)
          MX=MX/NX
          DO 80 NU=1,NX-1
            DO 70 LAMBDA=0,LX-1
              DO 60 I1=0,N1-1
                SPETAB(J)=EXPTAB(MX*LAMBDA*NU)
                J=J+1
60            CONTINUE
70          CONTINUE
80        CONTINUE
          LX=LX*NX
90      CONTINUE
100   CONTINUE
      END
