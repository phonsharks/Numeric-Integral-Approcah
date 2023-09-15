      PROGRAM INTEGRAL
      IMPLICIT REAL*8(A-H,O-Z)

      A=2.0
      B=1.0
      s1=0.0
      s2=A
      PI=DACOS(-1.D0)
      STAM=PI*A*B
         CALL GAUSS6(s1,s2,6,SGAUSS)
         WRITE(*,20) 4*SGAUSS,STAM,4*SGAUSS-STAM

20    FORMAT(3(3X,F16.8))
      END
      
      FUNCTION F(X)
      IMPLICIT REAL*8(A-H,O-Z)
      A=2.0
	B=1.0
      F=DSQRT((1-X**2/A**2)*B**2)
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
