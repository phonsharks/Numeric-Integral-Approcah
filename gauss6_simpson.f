      PROGRAM GAUSS_VE_SIMPSON
      IMPLICIT REAL*8(A-H,O-Z)
      N=6
      A=0.D0
      B=3.D0
	STAM=1.D0-DEXP(-3.D0)
      CALL SIMPSON(A,B,N,SSIMPSON)
      CALL GAUSS6(A,B,N,SGAUSS)
      WRITE(*,20) SSIMPSON-STAM, SGAUSS-STAM
20    FORMAT(2X,'SIMPSON =',F10.8,6X,'GAUSS =',F10.8)
      END

      FUNCTION F(X)
      IMPLICIT REAL*8(A-H,O-Z)
      F=DEXP(-X)
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
      
      SUBROUTINE SIMPSON(A,B,N,S)
      IMPLICIT REAL*8(A-H,O-Z)
      IF(N.LE.1.OR.A.GT.B) STOP 'HATALI VERI!'
      IF(MOD(N,2).EQ.1) STOP 'N SAYISI CIFT DEGIL.'
      H=(B-A)/N
      S=F(A)+F(B)
      DO 10 I=1,N-1
           X=A+I*H
           KATSAYI=2*(MOD(I,2)+1)
           S=S+KATSAYI*F(X)
10    CONTINUE
      S=S*H/3.0
      RETURN
      END
