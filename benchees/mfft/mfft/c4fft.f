        SUBROUTINE C4FFT(C,ID1,ID2,N1,N2,N3,N4,W1,W2,W3,W4,ISIG,IORD,
     $                   IWORK,IERR)
*
*PURPOSE:
*       THIS ROUTINE PERFORMS A 4-DIMENSIONAL FOURIER TRANSFORM,
*       OF ORDER N1*N2*N3*N4
****USAGE:
*       THE USER IS EXPECTED TO PROVIDE THE DATA IN A 4-DIMENSIONAL
*       COMPLEX ARRAY C, DIMENSIONED IN THE CALLING PROGRAM
*       DIMENSION C(ID1,ID2,N3,N4)
*       FOR OPTIMAL PERFORMANCE, ID1 AND ID2 SHOULD BE ODD, AND EQUAL
*       RESPECTIVELY TO N1 OR N1+1, N2 OR N2+1, DEPENDING IF N1 AND N2
*       ARE ODD OR EVEN. THE ROUTINE IS INTENDED FOR REPEATED USAGE, THUS
*       SEPARATE SET-UP AND OPERATING CALLS ARE AVAILABLE: THE USER SHOULD
*       ALLWAYS PERFORM A SET-UP CALL ( ISIG=0 ) PASSING THE
*       CORRECT VALUES OF THE ARGUMENTS, BEFORE PERFORMING THE
*       ACTUAL TRANSFORM ( ISIG= +1 OR -1 ); THE USER CAN CHOOSE
*       WHETHER TO OBTAIN THE RESULT OF A DIRECT TRANSFORM IN
*       NATURAL ORDER (ISIG=-1,IORD=1) OR IN "BIT-REVERSED" ORDER
*       (ISIG=-1,IORD=0); THIS CHOICE SAVES SOME COMPUTER TIME
*       IN CASES DISCUSSED IN THE LONG WRITE-UP. ANALOGOUSLY, THE
*       INVERSE TRANSFORM ACCEPTS INPUT IN NATURAL ORDER (ISIG=1,
*       IORD=1) OR IN "BIT-REVERSED" ORDER ( ISIG=1,IORD=-1).
****ARGUMENTS
*       C:
*          DECLARED COMPLEX C(ID1,ID2,N3,N4) IN THE CALLING PROGRAM
*          INPUT: ARRAY TO BE TRANSFORMED;
*          OUTPUT: TRANSFORMED ARRAY IF ISIG.NE.0; UNCHANGED IF ISIG.EQ.0
*       ID1:
*          INTEGER;
*          INPUT: FIRST DIMENSION OF C, AS DECLARED IN THE CALLING PROGRAM
*       ID2:
*          INTEGER;
*          INPUT: SECOND DIMENSION OF C, AS DECLARED IN THE CALLING
*               PROGRAM
*       N1:
*          INTEGER;
*          INPUT: ORDER OF THE TRANSFORM ALONG THE FIRST DIMENSION;
*                N1.LE.ID1;  N1 MUST BE A PRODUCT OF POWERS OF 2,3,5;
*       N2:
*          INTEGER;
*          INPUT: ORDER OF THE TRANSFORM ALONG THE SECOND DIMENSION;
*                 N2.LE.ID2 ; N2 MUST BE A PRODUCT OF POWERS OF 2,3,5;
*       N3:
*          INTEGER;
*          INPUT: THIRD DIMENSION OF C, AS DECLARED IN THE CALLING PROGRAM,
*                 AND ORDER OF THE TRANSFORM ALONG THE THIRD DIMENSION
*                 N3 MUST BE A PRODUCT OF POWERS OF 2,3,5;
*       N4:
*         INTEGER;
*         INPUT: ORDER OF THE TRANSFORM ALONG THE FOURTH DIMENSION
*                 N1 MUST BE A PRODUCT OF POWERS OF 2,3,5;
*       W1,W2,W3,W4:
*         THEY MUST BE DECLARED IN THE CALLING PROGRAM AS
*         COMPLEX W1(2*N1+7),W2(2*N2+7),W3(2*N3+7),W4(2*N4+7)
*         IF ANY OF THE N'S ARE EQUAL, THE CORRESPONDING W'S
*         DON'T NEED TO BE DISTINCT.
*         INPUT: IF ISIG.NE.0, THE W'S MUST CONTAIN THE TABLES PREPARED
*                BY A PREVIOUS CALL WITH THE SAME ARGUMENTS BUT ISIG.EQ.0.
*         OUTPUT: IF ISIG.EQ.0, THE W ARE FILLED WITH THE TABLES REQUIRED
*                 FOR THE TRANSFORM. IF ISIG.NE.0, UNCHANGED.
*       IORD:
*         INTEGER;
*         INPUT: IF IORD.EQ.0, NO REORDERING IS PERFORMED AFTER A DIRECT
*                TRANSFORM, SO THAT THE RESULTS ARE LEFT IN "BIT REVERSED"
*                ORDER; AND NO PERMUTATION IS PERFORMED BEFORE THE INVERSE
*                TRANSFORM, SO THAT THE INPUT MATRIX IS REQUIRED TO BE
*               ALREADY IN "BIT-REVERSED" ORDER; IF IORD.NE.0, OUTPUT OF
*              DIRECT  TRANSFORM AND INPUT TO INVERSE TRANSFORM ARE IN
*            NATURAL ORDER
*       ISIG:
*         INTEGER;
*         INPUT: IF ISIG.LT.0  , DIRECT TRANSFORM ( C(J)EXP(-2*PI*I*J*K/N))
*                IF ISIG.GT.0  ,INVERSE TRANSFORM ( C(J)EXP(2*PI*I*J*K/N))
*                IF ISIG.EQ.0  , SETUP RUN ( FILLING OF THE W'S )
*       IWORK:
*         INTEGER ARRAY, OF SIZE MAX(N1,N2,N3,N4)
*         WORK AREA;
*       IERR:
*         INTEGER:
*         OUTPUT: ERROR CODE:
*                 0 :NO ERRORS;
*                 1 :ID1.LT.N1 OR ID2.LT.N2
*                 2 :ONE OF THE N'S CONTAINS A FACTOR DIFFERENT FROM 2,3,5
*                 3 :IN A CALL WITH ISIG.NE.0, SOME OF THE W'S CONTAINED
*                      ILLEGAL VALUES
*
        COMPLEX C(0:*)
        INTEGER W1(-14:*),W2(-14:*),W3(-14:*),W4(-14:*)
        INTEGER IWORK(*)
        INTEGER IDERR,FACERR,TBERR
        PARAMETER(IDERR=1,FACERR=2,TBERR=3)
 
        IF(ID1.LT.N1.OR.ID2.LT.N2)THEN
                IERR=IDERR
                RETURN
        ENDIF
 
        IERR=0
        IF(ISIG.EQ.0)THEN
 
                CALL MFFTP(N1,W1,0,IERR)
                IF(IERR.NE.0)RETURN
                IF(N2.EQ.N1)THEN
                        CALL MFFTZ0(W1,1,4*N1+14,W2,1)
                ELSE
                        CALL MFFTP(N2,W2,0,IERR)
                        IF(IERR.NE.0)RETURN
                ENDIF
                IF(N3.EQ.N1)THEN
                        CALL MFFTZ0(W1,1,4*N1+14,W3,1)
                ELSE IF(N3.EQ.N2)THEN
                        CALL MFFTZ0(W2,1,4*N2+14,W3,1)
                ELSE
                        CALL MFFTP(N3,W3,0,IERR)
                        IF(IERR.NE.0)RETURN
                ENDIF
                IF(N4.EQ.N1)THEN
                        CALL MFFTZ0(W1,1,4*N1+14,W4,1)
                ELSE IF(N4.EQ.N2)THEN
                        CALL MFFTZ0(W2,1,4*N2+14,W4,1)
                ELSE IF (N4.EQ.N3)THEN
                        CALL MFFTZ0(W3,1,4*N3+14,W4,1)
                ELSE
                        CALL MFFTP(N4,W4,0,IERR)
                        IF(IERR.NE.0)RETURN
                ENDIF
                RETURN
 
        ELSE  IF(ISIG.GT.0)THEN
 
                IF(IORD.NE.0)THEN
                        CALL MFFTOV(C,1,ID1,N1,ID2*N3*N4,W1(3*N1),IWORK)
                        CALL MFFTOM(C,ID1,1,ID1*ID2,N2,N1,N3*N4,W2(3*N2)
     $                             ,IWORK)
                        CALL MFFTOM(C,ID1*ID2,ID1*ID2*N3,1,N3,N4,ID1*N2,
     $                              W3(3*N3),IWORK)
                        CALL MFFTOV(C,ID1*ID2*N3,1,N4,ID1*ID2*N3,
     $                              W4(3*N4),IWORK)
                ENDIF
 
                CALL MFFTIV(C,1,ID1,N1,ID2*N3*N4,W1,IERR)
                IF(IERR.NE.0)RETURN
                CALL MFFTIM(C,ID1,1,ID1*ID2,N2,N1,N3*N4,W2,IERR)
                IF(IERR.NE.0)RETURN
                CALL MFFTIM(C,ID1*ID2,ID1*ID2*N3,1,N3,N4,ID1*N2,W3,IERR)
                IF(IERR.NE.0)RETURN
                CALL MFFTIV(C,ID1*ID2*N3,1,N4,ID1*ID2*N3,W4,IERR)
 
        ELSE
 
                CALL MFFTDV(C,1,ID1,N1,ID2*N3*N4,W1,IERR)
                IF(IERR.NE.0)RETURN
                CALL MFFTDM(C,ID1,1,ID1*ID2,N2,N1,N3*N4,W2,IERR)
                IF(IERR.NE.0)RETURN
                CALL MFFTDM(C,ID1*ID2,ID1*ID2*N3,1,N3,N4,ID1*N2,W3,IERR)
                IF(IERR.NE.0)RETURN
                CALL MFFTDV(C,ID1*ID2*N3,1,N4,ID1*ID2*N3,W4,IERR)
 
                IF(IORD.NE.0)THEN
                        CALL MFFTOV(C,1,ID1,N1,ID2*N3*N4,W1(2*N1),IWORK)
                        CALL MFFTOM(C,ID1,1,ID1*ID2,N2,N1,N3*N4,W2(2*N2)
     $                             ,IWORK)
                        CALL MFFTOM(C,ID1*ID2,ID1*ID2*N3,1,N3,N4,ID1*N2,
     $                              W3(2*N3),IWORK)
                        CALL MFFTOV(C,ID1*ID2*N3,1,N4,ID1*ID2*N3,
     $                              W4(2*N4),IWORK)
                ENDIF
        ENDIF
        END
