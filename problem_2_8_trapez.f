      PROGRAM INTEGRAL
      IMPLICIT REAL*8(A-H,O-Z)

      A=2.0
      B=1.0
      s1=0.0
      s2=A
      PI=DACOS(-1.D0)
      STAM=PI*A*B
         CALL TRAPEZ(s1,s2,60,STRAPEZ)
         WRITE(*,20) 4*STRAPEZ,STAM,4*STRAPEZ-STAM

20    FORMAT(3(3X,F16.8))
      END
      
      FUNCTION F(X)
      IMPLICIT REAL*8(A-H,O-Z)
      A=2.0
	B=1.0
      F=DSQRT((1-X**2/A**2)*B**2)
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
