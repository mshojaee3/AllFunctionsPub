      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)

      INCLUDE 'ABA_PARAM.INC'

      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4 JSTEP(4)

      INTEGER I,J,K
      DOUBLE PRECISION young,nu,mu,lambda
      DOUBLE PRECISION F(3,3),Ft(3,3),C(3,3),E(3,3),S(3,3)
      DOUBLE PRECISION sig(3,3),tmp(3,3),Id(3,3)
      DOUBLE PRECISION trE,Jdet
      DOUBLE PRECISION h
      DOUBLE PRECISION sig0(6),sigp(6)
      DOUBLE PRECISION dEp(3,3),Fp(3,3),A(3,3)

C----- material
      young  = PROPS(1)
      nu     = PROPS(2)
      mu     = young/(2.D0*(1.D0+nu))
      lambda = young*nu/((1.D0+nu)*(1.D0-2.D0*nu))

C----- copy F
      DO I=1,3
        DO J=1,3
          F(I,J)=DFGRD1(I,J)
        END DO
      END DO

C----- identity
      DO I=1,3
        DO J=1,3
          Id(I,J)=0.D0
        END DO
        Id(I,I)=1.D0
      END DO

C----- compute stress
      CALL SVK_STRESS(F,Id,mu,lambda,sig)

C----- export stress
      CALL TENS_TO_VOIGT(sig,STRESS,NTENS)

C----- numerical tangent (robust)
      DO I=1,NTENS
        DO J=1,NTENS
          DDSDDE(I,J)=0.D0
        END DO
      END DO

      h = 1.D-8
      CALL TENS_TO_VOIGT6(sig,sig0)

      DO K=1,6
        CALL VOIGT_BASIS_SYM(K,dEp)
        DO I=1,3
          DO J=1,3
            dEp(I,J)=h*dEp(I,J)
          END DO
        END DO

C       Fp = (I + dEp) * F
        DO I=1,3
          DO J=1,3
            A(I,J)=Id(I,J)+dEp(I,J)
          END DO
        END DO
        CALL MATMUL33(A,F,Fp)

        CALL SVK_STRESS(Fp,Id,mu,lambda,sig)
        CALL TENS_TO_VOIGT6(sig,sigp)

        DO I=1,6
          DDSDDE(I,K)=(sigp(I)-sig0(I))/h
        END DO
      END DO

      SSE=0.D0
      SPD=0.D0
      SCD=0.D0
      RPL=0.D0
      RETURN
      END

C===========================================================
      SUBROUTINE SVK_STRESS(F,Id,mu,lambda,sig)
      IMPLICIT NONE
      DOUBLE PRECISION F(3,3),Id(3,3),sig(3,3)
      DOUBLE PRECISION mu,lambda
      DOUBLE PRECISION Ft(3,3),C(3,3),E(3,3),S(3,3),tmp(3,3)
      DOUBLE PRECISION trE,Jdet
      INTEGER I,J

      CALL TRANSPOSE33(F,Ft)
      CALL MATMUL33(Ft,F,C)

      DO I=1,3
        DO J=1,3
          E(I,J)=0.5D0*(C(I,J)-Id(I,J))
        END DO
      END DO

      trE = E(1,1)+E(2,2)+E(3,3)

      DO I=1,3
        DO J=1,3
          S(I,J)=2.D0*mu*E(I,J)
        END DO
      END DO
      DO I=1,3
        S(I,I)=S(I,I)+lambda*trE
      END DO

      CALL DET33(F,Jdet)
      CALL MATMUL33(F,S,tmp)
      CALL MATMUL33(tmp,Ft,sig)
      DO I=1,3
        DO J=1,3
          sig(I,J)=sig(I,J)/Jdet
        END DO
      END DO

      RETURN
      END

C===========================================================
      SUBROUTINE MATMUL33(A,B,C)
      IMPLICIT NONE
      DOUBLE PRECISION A(3,3),B(3,3),C(3,3)
      INTEGER I,J,K
      DO I=1,3
        DO J=1,3
          C(I,J)=0.D0
          DO K=1,3
            C(I,J)=C(I,J)+A(I,K)*B(K,J)
          END DO
        END DO
      END DO
      END

      SUBROUTINE TRANSPOSE33(A,AT)
      IMPLICIT NONE
      DOUBLE PRECISION A(3,3),AT(3,3)
      INTEGER I,J
      DO I=1,3
        DO J=1,3
          AT(I,J)=A(J,I)
        END DO
      END DO
      END

      SUBROUTINE DET33(A,detA)
      IMPLICIT NONE
      DOUBLE PRECISION A(3,3),detA
      detA = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))
     &     - A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1))
     &     + A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
      END

      SUBROUTINE VOIGT_BASIS_SYM(k,B)
      IMPLICIT NONE
      INTEGER k
      DOUBLE PRECISION B(3,3)
      INTEGER I,J
      DO I=1,3
        DO J=1,3
          B(I,J)=0.D0
        END DO
      END DO
      IF (k.EQ.1) THEN
        B(1,1)=1.D0
      ELSEIF (k.EQ.2) THEN
        B(2,2)=1.D0
      ELSEIF (k.EQ.3) THEN
        B(3,3)=1.D0
      ELSEIF (k.EQ.4) THEN
        B(1,2)=1.D0
        B(2,1)=1.D0
      ELSEIF (k.EQ.5) THEN
        B(1,3)=1.D0
        B(3,1)=1.D0
      ELSEIF (k.EQ.6) THEN
        B(2,3)=1.D0
        B(3,2)=1.D0
      END IF
      END

      SUBROUTINE TENS_TO_VOIGT6(A,v)
      IMPLICIT NONE
      DOUBLE PRECISION A(3,3),v(6)
      v(1)=A(1,1)
      v(2)=A(2,2)
      v(3)=A(3,3)
      v(4)=A(1,2)
      v(5)=A(1,3)
      v(6)=A(2,3)
      END

      SUBROUTINE TENS_TO_VOIGT(A,STRESS,NTENS)
      IMPLICIT NONE
      INTEGER NTENS
      DOUBLE PRECISION A(3,3),STRESS(NTENS)
      DOUBLE PRECISION v(6)
      CALL TENS_TO_VOIGT6(A,v)
      STRESS(1)=v(1)
      STRESS(2)=v(2)
      IF (NTENS.GE.3) STRESS(3)=v(3)
      IF (NTENS.GE.4) STRESS(4)=v(4)
      IF (NTENS.GE.5) STRESS(5)=v(5)
      IF (NTENS.GE.6) STRESS(6)=v(6)
      END
