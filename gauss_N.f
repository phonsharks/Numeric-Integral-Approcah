      PROGRAM GAUSS_N
      IMPLICIT REAL*8(A-H,O-Z)
      
      A=0.D0
      B=3.D0
	STAM=1.D0-DEXP(-3.D0)
      
      CALL GAUSS2(A,B,2,SGAUSS2)
	CALL GAUSS3(A,B,3,SGAUSS3)
	CALL GAUSS4(A,B,4,SGAUSS4)
	CALL GAUSS6(A,B,6,SGAUSS6)



      WRITE(*,20) SGAUSS2-STAM,SGAUSS3-STAM, SGAUSS4-STAM, SGAUSS6-STAM 
20    FORMAT(4(2X,F10.8))
      END

      FUNCTION F(X)
      IMPLICIT REAL*8(A-H,O-Z)
      F=DEXP(-X**2)
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
      
      SUBROUTINE GAUSS2(A,B,N,S)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(2),W(2)
            X(1) = -0.57735027
            X(2) = -X(1)
            W(1) = 1.
            W(2) = W(1)
            S=0.D0
            DO 10 I=1,N
                  X(I)=0.5D0*((B-A)*X(I)+B+A)
                  S=S+W(I)*F(X(I))
10          CONTINUE
            S=0.5D0*(B-A)*S
      RETURN
      END

	SUBROUTINE GAUSS3(A,B,N,S)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(3),W(3)
            X(1) = -0.77459667
            X(2) = -X(1)
		  X(3)=0.D0
            W(1) = 0.55555556
            W(2) = W(1)
		  W(3) = 0.88888889
            S=0.D0
            DO 10 I=1,N
                  X(I)=0.5D0*((B-A)*X(I)+B+A)
                  S=S+W(I)*F(X(I))
10          CONTINUE
            S=0.5D0*(B-A)*S
      RETURN
      END

	SUBROUTINE GAUSS4(A,B,N,S)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(4),W(4)
            X(1) = -0.86113631
            X(2) = -X(1)
		  X(3) =-0.33998104
		  X(4) =-X(3)
            W(1) = 0.34785485
            W(2) = W(1)
		  W(3) = 0.65214515
		  W(4) = W(3)
            S=0.D0
            DO 10 I=1,N
                  X(I)=0.5D0*((B-A)*X(I)+B+A)
                  S=S+W(I)*F(X(I))
10          CONTINUE
            S=0.5D0*(B-A)*S
      RETURN
      END

