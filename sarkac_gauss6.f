      PROGRAM SARKAC
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 L
      COMMON GENLIK
      DATA A,G,L/0.D0,9.8D0,1.D0/
      PI=DACOS(-1.D0)
      B=PI/2
      T_YAKLASIK=2*PI*DSQRT(L/G)

      DO 10 I=0,12
         GENLIK=5*I*PI/180.D0

         CALL GAUSS6(A,B,6,S)

         T_TAM=4*DSQRT(L/G)*S
         WRITE(*,20) GENLIK*180D0/PI,T_YAKLASIK,T_TAM, T_TAM-T_YAKLASIK
10    CONTINUE
20    FORMAT(2X,F6.1,2X,3F12.6)
      END


      FUNCTION F(X)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON GENLIK
      F=1.D0/DSQRT(1.D0-(DSIN(GENLIK/2)*DSIN(X))**2)
      RETURN
      END
      
      SUBROUTINE GAUSS6(A,B,N,S)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(6),W(6)
      IF(N.NE.6) STOP 'N SAYISI 6 OLMALI!'
            X(1) = -0.932469514203152
            X(2) = -0.661209386466265
            X(3) = -0.238619186083197
            X(4) = -X(3)
            X(5) = -X(2)
            X(6) = -X(1)
            W(1) = 0.171324492379170
            W(2) = 0.360761573048139
            W(3) = 0.467913934572691
            W(4) = W(3)
            W(5) = W(2)
            W(6) = W(1)
            S=0.D0
            DO 10 I=1,N
                  X(I)=0.5D0*((B-A)*X(I)+B+A)
                  S=S+W(I)*F(X(I))
10          CONTINUE
            S=0.5D0*(B-A)*S
      RETURN
      END