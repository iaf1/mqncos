      PROGRAM MAIN
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION XR(1000), FREV(1000), FREO(1000), FRED(1000), XMU(1000),
     + FREN(1000), DEN(1000), U(1000)
C     AA=10000. NUMBER OF ATOMS
C     N1 NUMBER OF INTEGRATION STEPS IN R-GRID
C     STEP, INTEGRATION STEP IN R-SPACE
C     A0 SCATTERING LENGTH IN H.O. UNITS
C     ALPHA, THE SCATERING PARAMETER FOR THE H.O. FUNCTION
C     TIME, THE STEP IN TIME
C     ITER, NUMBER OF ITERATIONS

      PRINT*, "INPUT: A0, N1, STEP, AA, TIME, ALPHA, ITER"
      READ(5,*) A0, N1, STEP, AA, TIME, ALPHA, ITER
 340  FORMAT(2D15.5)
      WRITE(6,1000) AA, N1, STEP, A0, ALPHA, TIME, ITER
 1000 FORMAT(/5X,'*',F10.0,' BOSONS IN A SPHERICAL TRAP*',/
     + 5X,'* R-GRID IN N1=',I4,' POINTS','R-STEP',F8.3/
     + 5X,'* A0=',E12.4, 2X, ' ALPHA=',E12.4/
     + 5X,'* TIME=',E12.4,' NUMBER-ITER =',I8)
      PI = 4.D0*DATAN(1.D0)
      PIIN = 1.D0/(4.D0*PI)
      PI2IN = DSQRT(PIIN)
      ALPHA2 = ALPHA*ALPHA
      CVAR = 2.D0*DSQRT(ALPHA)**3/DSQRT(DSQRT(PI)) !NORMAL CTANT
      
C     BUILDING THE STARTING WAVE FUNCTION R(r). Phi(r)=R(r)/r*Y00
      DO I=1,N1
            XR(I) = STEP*DFLOAT(I-1)
            XR2 = XR(I)*XR(I)
            FREV(I) = CVAR*XR(I)*DEXP(-0.5D0*ALPHA2*XR2)
            FREO(I) = FREV(I)
      ENDDO

C     STARTING THE CONVERGENCE PROCESS
C*************************************
C TO REDUCE TO THE ARMONIC OSCILLATOR CEQU=0
C*************************************
      CEQU = A0*AA
C     IF YOU PUT CEQU=0 YOU RECOVER THE H.O.
      AS3N = AA*A0*A0*A0
C     CEQU = 0.D0
      ITW = 0
      DO 2340 IT = 1, ITER
            ITW = ITW+1
 500        XNORM = 0.D0
            ENE0 = 0.D0

            DO I = 2, N1-1
                FRED(I) = (FREO(I-1)+FREO(I+1)-2.D0*FREO(I))/(STEP*STEP)
            ENDDO
            FRED(N1) = (FREO(N1-1)-2.D0*FREO(I))/(STEP*STEP)
            DO I = 1, N1
                  XR2 = XR(I)*XR(I)
                  IF (I .EQ. 1) THEN
                        XMU(I) = 0.D0
                  ELSE
                        ENE0 = ENE0-FREO(I)*FRED(I)*0.5D0
     +                        +0.5D0*XR2*FREO(I)*FREO(I)
     +                        +0.5D0*CEQU*XR2*(FREO(I)/XR(I))**4
                        XMU(I) = -0.5D0*FRED(I)/FREO(I)
     +                        +0.5D0*XR2+CEQU*(FREO(I)/XR(I))**2
                  ENDIF
                  FREN(I) = FREO(I) - TIME*XMU(I)*FREO(I)
                  XNORM = XNORM+FREN(I)*FREN(I)
            ENDDO
            XNORM = DSQRT(XNORM*STEP)
            ENE0 = ENE0*STEP
            IF (ITW .EQ. 200) THEN
                  WRITE(6,*) 'ENE0 =', ENE0
                  ITW = 0
            ENDIF
C           I DEFINE THE NEW WF.
            DO I = 1,N1
                  FREO(I) = FREN(I)/XNORM
            ENDDO
            IF (IT .EQ. ITER) THEN
                  WRITE(10,340) (XR(I), XMU(I), I=2,N1)
            ELSE
            ENDIF
 2340 CONTINUE
C     CALCULATION OF THE RADIOUS, POTENTIAL AND KINETIC ENERGY
C     SINGLE PARTICLE POTENTIAL
      DO I = 2, N1-1
            FRED(I) = (FREO(I-1)+FREO(I+1)-2.D0*FREO(I))/(STEP*STEP)
      ENDDO
      FRED(N1) = (FREO(N1-1)-2.D0*FREO(I))/(STEP*STEP)
      
      RADIOUS = 0.D0
      XKIN = 0.D0
      POTHO = 0.D0
      POTSELF = 0.D0
      CHEM = 0.D0
      XAVER = 0.D0
      XNORMDEN = 0.D0

      DO I = 2,N1
            XR2 = XR(I)*XR(I)
            RADIOUS = RADIOUS + XR2*FREO(I)*FREO(I)
            XKIN = XKIN + FREO(I)*FRED(I)
            POTHO = POTHO + XR2*FREO(I)*FREO(I)
            POTSELF = POTSELF + XR2*(FREO(I)/XR(I))**4
            CHEM = CHEM + XMU(I)*FREO(I)*FREO(I)
            U(I) = 0.5D0*XR2 + CEQU*(FREO(I)/XR(I))**2
            DEN(I) = (FREO(I)/XR(I))**2
            XNORMDEN = XNORMDEN + DEN(I)*XR2
            XAVER = XAVER + FREO(I)*FREO(I)*AS3N*DEN(I)
      ENDDO

      RADIOUS2 = RADIOUS*STEP
      RADIOUS = DSQRT(RADIOUS*STEP)
      XAVER = XAVER*STEP
      CHEM = CHEM*STEP
      XKIN = -XKIN*STEP*0.5D0
      POTHO = 0.5D0*POTHO*STEP
      POTSELF = POTSELF*STEP*CEQU*0.5D0
      POT = POTSELF + POTHO
      XNORMDEN = XNORMDEN*STEP
      
      WRITE(6,*) 'XNORMDEN =', XNORMDEN
      WRITE(6,740) ENE0,CHEM,XKIN,POT,POTHO,POTSELF,RADIOUS,RADIOUS2
 740  FORMAT(/5X,'* ENER', E12.5, '    AVERAGE CHEMICAL=', E12.5/
     + 5X, '* KIN-ENER=', E12.5,'    TOTAL-POT =', E12.5,/
     + 5X, '* POTHO =', E12.4, 2X, '    POTINT =', E12.4/
     + 5X, '* RADIOUS =', E12.4, 5X, 'RADIOUS2 =', E15.7, /)
      DO I = 2, N1
            WRITE (9, '(2E15.5)') XR(I), DEN(I)
      ENDDO
 345  CONTINUE
      STOP
      END

