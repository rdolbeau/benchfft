*=======================================================================
* Single Precision Complex Fast Fourier Transform
*
*  A subroutine to compute the discrete Fourier transform by the fastest
* available algorithm for any length input series.
*
* Reference:
*        Ferguson, W., (1979),   A Simple Derivation of Glassmans's
*          General N Fast Fourier Transform, MRC Tech. Summ. Rept. 2029,
*          Math. Res. Cent. U. of Wisconsin, Madison, Wis.
*
*  REFERENCES
*  ----------
*
* Routines called:
* SPCPFT
*
* Functions called:
* MOD
* FLOAT
*
* VAX extensions:
* DO WHILE
*=======================================================================

      SUBROUTINE SPCFFT(U,N,ISIGN,WORK,INTERP)

* VARIABLES
* ---------

      IMPLICIT NONE

      LOGICAL*1
     |  INU      
     | ,SCALE    

      INTEGER*4
     |  A        
     | ,B        
     | ,C        
     | ,N        
     | ,I        
     | ,ISIGN    

      REAL*4
     |  INTERP   

      COMPLEX*8
     |  U(*)            
     | ,WORK(*)         

*     Initialize parameters.

      A = 1
      B = N
      C = 1

      INU = .TRUE.

      IF (ISIGN.EQ.1) THEN

         SCALE = .TRUE.

      ELSE IF (ISIGN.EQ.-1) THEN

         SCALE = .FALSE.

      END IF

* Calculate Fourier transform by means of Glassman's algorithm

      DO WHILE ( B .GT. 1 )

         A = C * A

* Start of routine to get factors of N

         C = 2

* Increment C until it is an integer factor of B

         DO WHILE ( MOD(B,C) .NE. 0 )

                  C = C + 1

         END DO

* Calculate largest factor of B

         B = B / C


* Call Glassman's Fast Fourier Transform routine

         IF ( INU ) THEN

            CALL SPCPFT (A,B,C,U,WORK,ISIGN)

          ELSE

            CALL SPCPFT (A,B,C,WORK,U,ISIGN)

         END IF

* Set flag to swap input & output arrays to SPCPFT

         INU = ( .NOT. INU )

      END DO

* If odd number of calls to SPCPFT swap WORK into U

      IF ( .NOT. INU ) THEN

         DO I = 1, N
            U(I) = WORK(I)
         END DO

      END IF

* Scale the output for inverse Fourier transforms.

      IF ( SCALE ) THEN

         DO I = 1, N
            U(I) = U(I) / (N/INTERP)
         END DO

      END IF


* TERMINATION
* -----------

      RETURN
      END


*=======================================================================
* Single Precision Complex Prime Factor Transform
*
*  REFERENCES
*  ----------
*
* Calling routines:
* SPCFFT
*
* Subroutines called:
* -none-
*
* Functions called:
* CMLPX
* COS
* SIN
* FLOAT
*=======================================================================

      SUBROUTINE SPCPFT( A, B, C, UIN, UOUT, ISIGN )

* VARIABLES
* ---------

      IMPLICIT NONE

      INTEGER*4
     |  ISIGN           
     | ,A               
     | ,B               
     | ,C               
     | ,IA              
     | ,IB              
     | ,IC              
     | ,JCR             
     | ,JC              

      REAL*8
     |  ANGLE

      COMPLEX*8
     |  UIN(B,C,A)      
     | ,UOUT(B,A,C)     
     | ,DELTA           
     | ,OMEGA           
     | ,SUM             

* ALGORITHM
* ---------

* Initialize run time parameters.


      ANGLE =6.28318530717958 / FLOAT( A * C )
      OMEGA = CMPLX( 1.0, 0.0 )

* Check the ISIGN of the transform.

      IF( ISIGN .EQ. 1 ) THEN

         DELTA = CMPLX( DCOS(ANGLE), DSIN(ANGLE) )

      ELSE

         DELTA = CMPLX( DCOS(ANGLE), -DSIN(ANGLE) )

      END IF



* Do the computations.

      DO IC = 1, C

         DO IA = 1, A

            DO IB = 1, B

               SUM = UIN( IB, C, IA )

               DO JCR = 2, C

                  JC = C + 1 - JCR
                  SUM = UIN( IB, JC, IA ) + OMEGA * SUM

               END DO

               UOUT( IB, IA, IC ) = SUM

            END DO

            OMEGA = DELTA * OMEGA

         END DO

      END DO

* TERMINATION
* -----------

      RETURN
      END


