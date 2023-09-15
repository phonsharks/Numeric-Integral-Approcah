      PROGRAM SARKAC
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 L
      COMMON GENLIK
      DATA A,G,L/0.D0,9.8D0,1.D0/
      PI=DACOS(-1.D0)
      B=PI/2
      T_YAKLASIK=2*PI*DSQRT(L/G)
      N=200

      DO 10 I=0,12
         GENLIK=5*I*PI/180.D0
         CALL TRAPEZ(A,B,N,S)
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
      
      SUBROUTINE TRAPEZ(A,B,N,S)
      IMPLICIT REAL*8(A-H,O-Z)
      IF(N.LE.1.OR.A.GT.B) STOP 'HATALI VERI!'
      H=(B-A)/N
      S=0.5*(F(A)+F(B))
      DO 10 I=1,N-1
         X=A+I*H
         S=S+F(X)
10    CONTINUE
      S=H*S
      RETURN
      END
