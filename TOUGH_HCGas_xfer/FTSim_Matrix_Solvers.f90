!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>                                                                     >
!>                                                                     >
!>                ITERATIVE AND DIRECT MATRIX SOLVERS                  >
!>                                                                     >
!>                                                                     >
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
      MODULE Matrix_Solvers
!
!
!***********************************************************************
!*                                                                     *
!*                 Direct Solvers (LU Decomposition)                   *
!*                                                                     *
!***********************************************************************
!
!
         PUBLIC :: DGBTRF,DGBTRS
!
!
!
         CONTAINS
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!      D I R E C T    S O L V E R S
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE DGBTRF(M,N,KL,KU,AB,LDAB,IPIV,INFO)
!
! ----------
! ... Statement of implicitness
! ----------
!
      IMPLICIT NONE
!
! ----------
! ... Integer parameters
! ----------
!
      INTEGER, PARAMETER :: NBMAX = 32, LDWORK = NBMAX+1
!
! ----------
! ... Integer variables
! ----------
!
      INTEGER :: M, N, KL, KU, LDAB, INFO
      INTEGER :: i, j, jp, JM, KV, IP
      INTEGER :: NB, JU, JB, I2, I3, KM, NW, JJ, J2, J3, K2, II
!
! ----------
! ... Integer arrays
! ----------
!
      INTEGER, DIMENSION(N) :: IPIV
      INTEGER, DIMENSION(1) :: MaxL
!
! ----------
! ... Double precision variables
! ----------
!
      REAL(KIND = 8) :: TEMP
!
! ----------
! ... Double precision arrays
! ----------
!
      REAL(KIND = 8), DIMENSION(LDAB,N)       :: AB
      REAL(KIND = 8), DIMENSION(LDWORK,NBMAX) :: WORK13,WORK31
!
!
!=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of DGBTRF
!
!
      KV = KU + KL
      IF(KU <= 64) THEN
         NB = 1
      ELSE
         NB = 32
      END IF
!
! >>>>>>>
!
      NB = MIN(NB,NBMAX)
!
      IF_Option: IF(NB <= 1 .OR. NB > KL) THEN
!
         CALL DGBTF2(M,N,KL,KU,AB,LDAB,IPIV,INFO)
!
      ELSE
!
         FORALL (j = 1:NB, i = 1:NB, i /= j ) WORK13(i,j) = 0.0d0
!
         JU = 1
!
!
         DO_BigLoop: DO j = 1, MIN( M, N ), NB
!
            JB = MIN( NB, MIN( M, N )-j+1 )
            I2 = MIN( KL-JB, M-j-JB+1 )
            I3 = MIN( JB, M-j-KL+1 )
!
!
            DO_MLoop1: DO JJ = j, j + JB - 1
!
               KM        = MIN( KL, M-JJ )
               MaxL(1:1) = MAXLOC(ABS(AB(KV+1:KV+KM+1,JJ))) - KV
!
               jp         = MaxL(1)
               IPIV( JJ ) = JP + JJ - j
!
!
               IF_Pos1: IF( AB(KV+JP,JJ) /= 0.0d0 ) THEN
!
                  JU = MAX(JU,MIN(JJ+KU+JP-1,N))
                  IF(JP /= 1) THEN
                     IF( (JP+JJ-1) < (J+KL) ) THEN
                        CALL DSWAP(JB,AB(KV+1+JJ-J,J),LDAB-1,AB(KV+JP+JJ-J,J),LDAB-1)
                     ELSE
                        CALL DSWAP(JJ-J,AB(KV+1+JJ-J,J),LDAB-1,WORK31(JP+JJ-J-KL,1),LDWORK)
                        CALL DSWAP(J+JB-JJ,AB(KV+1,JJ),LDAB-1,AB(KV+JP,JJ),LDAB-1)
                     END IF
                  END IF
!
                  AB(KV+2:KV+1+km,JJ) = AB(KV+2:KV+1+km,JJ)/AB(KV+1,JJ)
!
                  JM = MIN(JU, J+JB-1)
                  IF(JM > JJ) CALL DGER(KM,JM-JJ,-1.0d0,AB(KV+2,JJ),1, AB(KV,JJ+1),LDAB-1,AB(KV+1,JJ+1),LDAB-1)
!
               ELSE
!
                  IF(INFO == 0) INFO = JJ
!
               END IF IF_Pos1
!
               nw = MIN( JJ-J+1, I3 )
               IF(NW > 0) work31(1:nw,jj-j+1) = ab(kv+kl-jj+j+1:kv+kl-jj+j+nw,jj)
!
            END DO DO_MLoop1
!
!
!
            IF_LessN: IF( J+JB <= N ) THEN
!
               J2 = MIN( JU-J+1, KV ) - JB
               J3 = MAX( 0, JU-J-KV+1 )
!
               CALL DLASWP(J2,AB(KV+1-JB,J+JB),LDAB-1,1,JB,IPIV(J))
!
               DO I = J, J + JB - 1
                  IPIV(I) = IPIV(I) + J - 1
               END DO
!
               K2 = J - 1 + JB + J2
               DO_SLoop1: DO i = 1, J3
                  JJ = K2 + i
                  DO II = j+i-1, j+JB-1
                     IP = IPIV( II )
                     IF( IP /= II ) THEN
                        TEMP                 = AB( KV+1+II-JJ, JJ )
                        AB( KV+1+II-JJ, JJ ) = AB( KV+1+IP-JJ, JJ )
                        AB( KV+1+IP-JJ, JJ ) = TEMP
                     END IF
                  END DO
               END DO DO_SLoop1
!
!
               IF_J2: IF(J2 > 0) THEN
!
                  CALL DTRSM( JB, J2, AB( KV+1, J ), LDAB-1, AB( KV+1-JB, J+JB ), LDAB-1 )
!
                  IF(I2 > 0) CALL DGEMM(I2,J2,JB, -1.0d0, AB(KV+1+JB,J), LDAB-1, AB(KV+1-JB, J+JB), LDAB-1, AB(KV+1, J+JB), LDAB-1)
!
                  IF(I3 > 0) CALL DGEMM(I3,J2,JB, -1.0d0, WORK31, LDWORK, AB(KV+1-JB, J+JB), LDAB-1, AB(KV+KL+1-JB, J+JB), LDAB-1)
!
               END IF IF_J2
!
!
               IF_J3: IF(J3 > 0) THEN
!
                  DO JJ = 1, J3
                     WORK13( JJ:JB, JJ ) = AB( 1:JB-JJ+1, JJ+j+KV-1 )
                  END DO
!
                  CALL DTRSM( JB, J3, AB( KV+1, J ), LDAB-1, WORK13, LDWORK )
!
                  IF(I2 > 0) CALL DGEMM(I2, J3, JB, -1.0d0, AB( KV+1+JB, J ), LDAB-1, WORK13, LDWORK, AB( 1+JB, J+KV ), LDAB-1)
!
                  IF(I3 > 0) CALL DGEMM(I3, J3, JB, -1.0d0, WORK31, LDWORK, WORK13, LDWORK, AB( 1+KL, J+KV ), LDAB-1)
!
                  DO JJ = 1, J3
                     AB( 1:JB-JJ+1, JJ+j+KV-1 ) = WORK13( JJ:JB, JJ )
                  END DO
               END IF  IF_J3
!
            ELSE
!
               DO i = j, j + JB - 1
                  IPIV(i) = IPIV(i) + j - 1
               END DO
!
            END IF IF_LessN
!
!
            DO_JJLoop: DO JJ = J + JB - 1, J, -1
!
               JP = IPIV(JJ) - JJ + 1
!
               IF( JP /= 1 ) THEN
                  IF( (JP+JJ-1) < (J+KL) ) THEN
                     CALL DSWAP(JJ-J,AB(KV+1+JJ-J,J),LDAB-1,AB(KV+JP+JJ-J,J),LDAB-1)
                  ELSE
                     CALL DSWAP(JJ-J,AB(KV+1+JJ-J,J),LDAB-1,WORK31(JP+JJ-J-KL,1),LDWORK)
                  END IF
               END IF
!
               NW = MIN( I3, JJ-J+1 )
               IF(NW > 0) ab(kv+kl-jj+j+1:kv+kl-jj+j+nw,jj) = work31(1:nw,jj-j+1)
!
            END DO DO_JJLoop
!
         END DO DO_BigLoop
!
      END IF IF_Option
!
!
!=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of DGBTRF
!
!
      RETURN
      END SUBROUTINE DGBTRF
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE DGBTF2(M,N,KL,KU,AB,LDAB,IPIV,INFO)
!
! ----------
! ... Statement of implicitness
! ----------
!
      IMPLICIT NONE
!
! ----------
! ... Integer variables
! ----------
!
      INTEGER :: M, N, KL, KU, LDAB, INFO
      INTEGER :: j, i6, jp, KV, JU, KM
!
! ----------
! ... Integer arrays
! ----------
!
      INTEGER, DIMENSION(N) :: IPIV
      INTEGER, DIMENSION(1) :: MaxL
!
! ----------
! ... Double precision variables
! ----------
!
      REAL(KIND = 8) :: ddttmm
!
! ----------
! ... Double precision arrays
! ----------
!
      REAL(KIND = 8), DIMENSION(LDAB,N) :: AB
!
!
!=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of DGBTF2
!
!
      KV = KU + KL
      JU = 1
!
! >>>>>
!
      DO_BigLoop: DO j = 1, MIN( M, N )
!
         KM        = MIN(KL,M-j)
         MaxL(1:1) = MAXLOC(ABS(AB(KV+1:KV+KM+1,j)))
!
         jp      = MaxL(1)
         IPIV(j) = JP + j - 1
!
         IF(AB(KV+JP,j) /= 0.0d0) THEN
!
            JU = MAX(JU,MIN(j+KU+JP-1,N))
!
            IF(JP /= 1) THEN
               DO i6=1,(JU-j+1)*(LDAB-1),LDAB-1
                  ddttmm           = ab(KV+jp-1+i6,j)
                  ab(KV+jp-1+i6,j) = ab(KV+i6,j)
                  ab(KV+i6,J)      = ddttmm
               END DO
            END IF
!
            IF( KM > 0 ) THEN
               ab(KV+2:KV+1+km,j) = ab(KV+2:KV+1+km,j)/AB(KV+1,j)
               IF( JU > J ) CALL DGER(KM,JU-J,-1.0d0,AB(KV+2,J),1,AB(KV,J+1),LDAB-1,AB(KV+1,J+1),LDAB-1)
            END IF
!
         ELSE
!
            IF( INFO == 0 ) INFO = J
!
         END IF
!
      END DO DO_BigLoop
!
!
!=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of DGBTF2
!
!
      RETURN
      END SUBROUTINE DGBTF2
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE DGBTRS(N,KL,KU,AB,LDAB,IPIV,B,LDB)
!
! ----------
! ... Statement of implicitness
! ----------
!
      IMPLICIT NONE
!
! ----------
! ... Integer variables
! ----------
!
      INTEGER :: N, KL, KU, LDAB, LDB
      INTEGER :: j, L, KD, LM
!
! ----------
! ... Integer arrays
! ----------
!
      INTEGER, DIMENSION(N) :: IPIV
!
! ----------
! ... Double precision arrays
! ----------
!
      REAL(KIND = 8), DIMENSION(LDAB,N) :: AB
      REAL(KIND = 8), DIMENSION(LDB,1)  :: B
!
!
!=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of DGBTRS
!
!
      KD = KU + KL + 1
      DO j = 1,N-1
         LM = MIN(KL,N-J)
         L  = IPIV(j)
         IF(L /= j) CALL DSWAP(1,B(L,1),LDB,B(j,1),LDB)
         CALL DGER(LM,1,-1.0d0,AB(KD+1,j),1,B(j,1),LDB,B(j+1,1),LDB)
      END DO
!
! >>>>>
!
      CALL DTBSV(N,KL+KU,AB,LDAB,B(1,1))
!
!
!=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of DGBTRS
!
!
      RETURN
      END SUBROUTINE DGBTRS
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE DLASWP(N,A,LDA,K1,K2,IPIV)
!
! ----------
! ... Statement of implicitness
! ----------
!
      IMPLICIT NONE
!
! ----------
! ... Integer variables
! ----------
!
      INTEGER :: N, LDA, K1, K2
      INTEGER :: i, IP
!
! ----------
! ... Integer arrays
! ----------
!
      INTEGER, DIMENSION(K1:K2) :: IPIV
!
! ----------
! ... Double precision arrays
! ----------
!
      REAL(KIND = 8), DIMENSION(LDA,*) :: A
!
!
!=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of DLASWP
!
!
      DO i = K1, K2
         IP = IPIV(i)
         IF(IP /= i) CALL DSWAP(N,A(i,1),LDA,A(IP,1),LDA)
      END DO
!
!
!=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of DLASWP
!
!
      RETURN
      END SUBROUTINE DLASWP
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE DGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
!
! ----------
! ... Statement of implicitness
! ----------
!
      IMPLICIT NONE
!
! ----------
! ... Integer variables
! ----------
!
      INTEGER :: M, N, INCX, INCY, LDA
      INTEGER :: j, JY, KX
!
! ----------
! ... Double precision variables
! ----------
!
      REAL(KIND = 8) :: ALPHA
!
! ----------
! ... Double precision arrays
! ----------
!
      REAL(KIND = 8), DIMENSION(*)     :: X, Y
      REAL(KIND = 8), DIMENSION(LDA,*) :: A
!
!
!=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of DGER
!
!
      IF(INCY > 0) THEN
         JY = 1
      ELSE
         JY = 1 - ( N - 1 )*INCY
      END IF
!
      IF(INCX == 1) THEN

         FORALL (j = 1:N, Y(JY + (j-1)*INCY) /= 0.0d0)
            A(1:M,j) = A(1:M,j) + ALPHA*X(1:M)*Y(JY + (j-1)*INCY)
         END FORALL
!
      ELSE
!
         IF(INCX > 0) THEN
            KX = 1
         ELSE
            KX = 1 - ( M - 1 )*INCX
         END IF
!
         FORALL (j = 1:N, Y(JY + (j-1)*INCY) /= 0.0d0)
            A(1:M,j) = A(1:M,j) + ALPHA*X(KX:KX+(M-1)*INCX:INCX)*Y(JY + (j-1)*INCY)
         END FORALL
!
      END IF
!
!
!=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of DGER
!
!
      RETURN
      END SUBROUTINE DGER
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE DTBSV(N,K,A,LDA,X)
!
! ----------
! ... Statement of implicitness
! ----------
!
      IMPLICIT NONE
!
! ----------
! ... Integer variables
! ----------
!
      INTEGER :: N, K, LDA
      INTEGER :: j, L, i_b, i_e, Kp1
!
! ----------
! ... Double precision variables
! ----------
!
      REAL(KIND = 8) :: TEMP
!
! ----------
! ... Double precision arrays
! ----------
!
      REAL(KIND = 8), DIMENSION(N)     :: X
      REAL(KIND = 8), DIMENSION(LDA,N) :: A
!
!
!=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of DTBSV
!
!
      Kp1  = K + 1
!
      DO j = N, 1, -1
!
         IF_Pos: IF(X(j) /= 0.0d0) THEN
!
            L    = Kp1 - j
            X(j) = X(j)/A(Kp1,j)
            TEMP = X(j)
!
            i_b = j-1
            i_e = MAX(1,j-K)
!
            X(i_b:i_e:-1) = X(i_b:i_e:-1) - TEMP*A(L+i_b:L+i_e:-1,j)
!
         END IF IF_Pos
!
      END DO
!
!
!=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of DTBSV
!
!
      RETURN
      END SUBROUTINE DTBSV
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE DGEMM(M,N,K,ALPHA,A,LDA,B,LDB,C,LDC)
!
! ----------
! ... Statement of implicitness
! ----------
!
      IMPLICIT NONE
!
! ----------
! ... Integer variables
! ----------
!
      INTEGER :: M, N, K, LDA, LDB, LDC
      INTEGER :: i, j, l
!
! ----------
! ... Double precision variables
! ----------
!
      REAL(KIND = 8) :: ALPHA
!
! ----------
! ... Double precision arrays
! ----------
!
      REAL(KIND = 8), DIMENSION (LDA,K) :: A
      REAL(KIND = 8), DIMENSION (LDB,N) :: B
      REAL(KIND = 8), DIMENSION (LDC,N) :: C
!
!
!=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of DGEMM
!
!
      DO l=1,K
        FORALL(j=1:N, i=1:M, B(l,j) /= 0.0d0)  C(i,j) = C(i,j) - ALPHA*B(l,j)*A(i,l)
      END DO
!
!
!=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of DGEMM
!
!
      RETURN
      END SUBROUTINE DGEMM
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE DTRSM(M,N,A,LDA,B,LDB)
!
! ----------
! ... Statement of implicitness
! ----------
!
      IMPLICIT NONE
!
! ----------
! ... Integer variables
! ----------
!
      INTEGER :: M, N, LDA, LDB
      INTEGER :: j, k
!
! ----------
! ... Double precision arrays
! ----------
!
      REAL(KIND = 8), DIMENSION (LDA,M) :: A
      REAL(KIND = 8), DIMENSION (LDB,N) :: B
!
!
!=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of DTRSM
!
!
      FORALL(j=1:N, k=1:M, B(k,j) /= 0.0d0)
         B(k+1:M,j) = B(k+1:M,j) - B(k,j)*A(k+1:M,k)
      END FORALL
!
!
!=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of DTRSM
!
!
      RETURN
      END SUBROUTINE DTRSM
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE DSWAP(n,dx,incx,dy,incy)
!
! ----------
! ... Statement of implicitness
! ----------
!
      IMPLICIT NONE
!
! ----------
! ... Integer variables
! ----------
!
      INTEGER :: n, incx, incy
!
! ----------
! ... Double precision arrays
! ----------
!
      REAL(KIND = 8), DIMENSION(N) :: dtemp
!
      REAL(KIND = 8), DIMENSION(1+(N-1)*incx) :: dx
      REAL(KIND = 8), DIMENSION(1+(N-1)*incy) :: dy
!
!
!=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of DSWAP
!
!
      IF(n <= 0) RETURN
!
      IF(incx == 1 .AND. incy == 1) THEN
         dtemp = dx
         dx    = dy
         dy    = dtemp
         RETURN
      END IF
!
!     Code for unequal increments or equal increments not equal to 1
!
      dtemp(1:N)              = dx(1:1+(N-1)*incx:incx)
      dx(1:1+(N-1)*incx:incx) = dy(1:1+(N-1)*incy:incy)
      dy(1:1+(N-1)*incy:incy) = dtemp(1:N)
!
!
!=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of DSWAP
!
!
      RETURN
      END SUBROUTINE DSWAP
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      END MODULE Matrix_Solvers
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!      I T E R A T I V E    S O L V E R S
!
!
!
!***************************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE DBCG(N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MTTVEC, MSOLVE, MTSOLV, ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT,  &
     &                R, Z, P, RR, ZZ, PP, DZ, RWORK, IWORK)
!
!
!
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX
      INTEGER ITER, IERR, IUNIT, IWORK(*)
      DOUBLE PRECISION B(N), X(N), A(NELT), TOL, ERR, R(N), Z(N), P(N)
      DOUBLE PRECISION RR(N), ZZ(N), PP(N), DZ(N), RWORK(*)

      DOUBLE PRECISION DMACH(5)
      DATA DMACH(3) / 1.1101827117665D-16 /

      EXTERNAL MATVEC, MTTVEC, MSOLVE, MTSOLV
!
      ITER = 0
      IERR = 0
      IF( N.LT.1 ) THEN
         IERR = 3
         RETURN
      ENDIF
      FUZZ = DMACH(3)
      TOLMIN = 500.0*FUZZ
      FUZZ = FUZZ*FUZZ
      IF( TOL.LT.TOLMIN ) THEN
         TOL = TOLMIN
         IERR = 4
      ENDIF
!
      CALL MATVEC(N, X, R, NELT, IA, JA, A, ISYM)
      DO 10 I = 1, N
         R(I)  = B(I) - R(I)
         RR(I) = R(I)
 10   CONTINUE

      CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)
      CALL MTSOLV(N, RR, ZZ, NELT, IA, JA, A, ISYM, RWORK, IWORK)
!
      IF( ISDBCG(N, B, X, NELT, IA, JA, A, ISYM, MSOLVE, ITOL, TOL,  &
     &     ITMAX, ITER, ERR, IERR, IUNIT, R, Z, P, RR, ZZ, PP,       &
     &     DZ, RWORK, IWORK, AK, BK, BNRM, SOLNRM) .NE. 0 )    GO TO 200

      IF( IERR.NE.0 ) RETURN
!
!
      DO 100 K=1,ITMAX
         ITER = K
!
         DDOT = 0.D0
         DO 15 I = 1,N
           DDOT = DDOT + Z(I)*RR(I)
   15    CONTINUE
         BKNUM = DDOT

         IF( ABS(BKNUM).LE.FUZZ ) THEN
            IERR = 6
            RETURN
         ENDIF
         IF(ITER .EQ. 1) THEN
            DO 17 I = 1,N
              P(I) = Z(I)
   17       CONTINUE

            DO 18 I = 1,N
              PP(I) = ZZ(I)
   18       CONTINUE
         ELSE
            BK = BKNUM/BKDEN
            DO 20 I = 1, N
               P(I) = Z(I) + BK*P(I)
               PP(I) = ZZ(I) + BK*PP(I)
 20         CONTINUE
         ENDIF
         BKDEN = BKNUM
!
         CALL MATVEC(N, P, Z, NELT, IA, JA, A, ISYM)

         DDOT = 0.D0
         DO 25 I = 1,N
           DDOT = DDOT + PP(I)*Z(I)
   25    CONTINUE
         AKDEN = DDOT

         AK = BKNUM/AKDEN
         IF( ABS(AKDEN).LE.FUZZ ) THEN
            IERR = 6
            RETURN
         ENDIF

          DO 26 I = 1,N
             X(I) = X(I) + AK*P(I)
   26     CONTINUE

          DO 27 I = 1,N
             R(I) = R(I) - AK*Z(I)
   27     CONTINUE

         CALL MTTVEC(N, PP, ZZ, NELT, IA, JA, A, ISYM)

          DO 28 I = 1,N
             RR(I) = RR(I) - AK*ZZ(I)
   28     CONTINUE

         CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)
         CALL MTSOLV(N, RR, ZZ, NELT, IA, JA, A, ISYM, RWORK, IWORK)
!
         IF( ISDBCG(N, B, X, NELT, IA, JA, A, ISYM, MSOLVE, ITOL, TOL, &
     &        ITMAX, ITER, ERR, IERR, IUNIT, R, Z, P, RR, ZZ,          &
     &        PP, DZ, RWORK, IWORK, AK, BK, BNRM, SOLNRM) .NE. 0 ) GO TO 200
!
 100  CONTINUE
!
!
      ITER = ITMAX + 1
      IERR = 2
!
 200  RETURN
      END
!
!
!
!***************************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***************************************************************************
!
!
      SUBROUTINE DSLUBC(N, B, X, NELT, IA, JA, A, ISYM, ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, RWORK, LENW, IWORK, LENIW )
!
!
!
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX, ITER
      INTEGER IERR, IUNIT, LENW, IWORK(LENIW), LENIW
      DOUBLE PRECISION B(N), X(N), A(NELT), TOL, ERR, RWORK(LENW)
      EXTERNAL DSMV, DSMTV, DSLUI, DSLUTI
      PARAMETER (LOCRB=1, LOCIB=11)
!
      IERR = 0
      IF( N.LT.1 .OR. NELT.LT.1 ) THEN
         IERR = 3
         RETURN
      ENDIF
      CALL DS2Y( N, NELT, IA, JA, A, ISYM, IUNIT )
!
      NL = 0
      NU = 0
      DO 20 ICOL = 1, N
!         Don't count diagonal.
         JBGN = JA(ICOL)+1
         JEND = JA(ICOL+1)-1
         IF( JBGN.LE.JEND ) THEN
            DO 10 J = JBGN, JEND
               IF( IA(J).GT.ICOL ) THEN
                  NL = NL + 1
                  IF( ISYM.NE.0 ) NU = NU + 1
               ELSE
                  NU = NU + 1
               ENDIF
 10         CONTINUE
         ENDIF
 20   CONTINUE
!
      LOCIL = LOCIB
      LOCJL = LOCIL + N+1
      LOCIU = LOCJL + NL
      LOCJU = LOCIU + NU
      LOCNR = LOCJU + N+1
      LOCNC = LOCNR + N
      LOCIW = LOCNC + N
!
      LOCL = LOCRB
      LOCDIN = LOCL + NL
      LOCU = LOCDIN + N
      LOCR = LOCU + NU
      LOCZ = LOCR + N
      LOCP = LOCZ + N
      LOCRR = LOCP + N
      LOCZZ = LOCRR + N
      LOCPP = LOCZZ + N
      LOCDZ = LOCPP + N
      LOCW = LOCDZ + N
!
      CALL DCHKW( 'DSLUBC', LOCIW, LENIW, LOCW, LENW, IERR, ITER, ERR )
      IF( IERR.NE.0 ) RETURN
!
      IWORK(1) = LOCIL
      IWORK(2) = LOCJL
      IWORK(3) = LOCIU
      IWORK(4) = LOCJU
      IWORK(5) = LOCL
      IWORK(6) = LOCDIN
      IWORK(7) = LOCU
      IWORK(9) = LOCIW
      IWORK(10) = LOCW
!
      CALL DSILUS( N, NELT, IA, JA, A, ISYM, NL, IWORK(LOCIL),         &
     &     IWORK(LOCJL), RWORK(LOCL), RWORK(LOCDIN), NU, IWORK(LOCIU), &
     &     IWORK(LOCJU), RWORK(LOCU), IWORK(LOCNR), IWORK(LOCNC) )
!
      CALL DBCG(N, B, X, NELT, IA, JA, A, ISYM, DSMV, DSMTV,           &
     &     DSLUI, DSLUTI, ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT,    &
     &     RWORK(LOCR), RWORK(LOCZ), RWORK(LOCP), RWORK(LOCRR), RWORK(LOCZZ), &
     &     RWORK(LOCPP), RWORK(LOCDZ), RWORK, IWORK )
      RETURN
      END
!
!
!
!***************************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***************************************************************************
!
!
      FUNCTION ISDBCG(N, B, X, NELT, IA, JA, A, ISYM, MSOLVE, ITOL, &
     &     TOL, ITMAX, ITER, ERR, IERR, IUNIT, R, Z, P, RR, ZZ, PP, DZ, &
     &     RWORK, IWORK, AK, BK, BNRM, SOLNRM)
!
!
!
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX
      INTEGER ITER, IERR, IUNIT, IWORK(*)
      DOUBLE PRECISION B(N), X(N), A(NELT), TOL, ERR, R(N), Z(N), P(N)
      DOUBLE PRECISION RR(N), ZZ(N), PP(N), DZ(N), RWORK(*)
      DOUBLE PRECISION AK, BK, BNRM, SOLNRM
      COMMON /SOLBLK/ SOLN(1)
      EXTERNAL MSOLVE
!
      ISDBCG = 0
!
      IF( ITOL.EQ.1 ) THEN
         IF(ITER .EQ. 0) BNRM = DNRM2(N, B, 1)
         ERR = DNRM2(N, R, 1)/BNRM
      ELSE IF( ITOL.EQ.2 ) THEN
         IF(ITER .EQ. 0) THEN
            CALL MSOLVE(N, B, DZ, NELT, IA, JA, A, ISYM, RWORK, IWORK)
            BNRM = DNRM2(N, DZ, 1)
         ENDIF
         ERR = DNRM2(N, Z, 1)/BNRM
      ELSE IF( ITOL.EQ.11 ) THEN
         IF(ITER .EQ. 0) SOLNRM = DNRM2(N, SOLN, 1)
         DO 10 I = 1, N
            DZ(I) = X(I) - SOLN(I)
 10      CONTINUE
         ERR = DNRM2(N, DZ, 1)/SOLNRM
      ELSE
         ERR = 1.0E10
         IERR = 3
      ENDIF
!
      IF(IUNIT .NE. 0) THEN
         IF( ITER.EQ.0 ) THEN
            WRITE(IUNIT,1000) N, ITOL
         ENDIF
         WRITE(IUNIT,1010) ITER, ERR, AK, BK
      ENDIF
      IF(ERR .LE. TOL) ISDBCG = 1
!
      RETURN
 1000 FORMAT(' Preconditioned BiConjugate Gradient for N, ITOL = ', I5,I5,/, &
     &       ' ITER','   Error Estimate','            Alpha             Beta')
 1010 FORMAT(1X,I4,1X,E16.7,1X,E16.7,1X,E16.7)
      END
!
!
!
!***************************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***************************************************************************
!
!
      SUBROUTINE DCGS(N, B, X, NELT, IA, JA, A, ISYM, MATVEC,  &
     &     MSOLVE, ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT,   &
     &     R, R0, P, Q, U, V1, V2, RWORK, IWORK)
!
!
!
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX
      INTEGER ITER, IERR, IUNIT, IWORK(*)
      DOUBLE PRECISION B(N), X(N), A(NELT), TOL, ERR, R(N), R0(N), P(N)
      DOUBLE PRECISION Q(N), U(N), V1(N), V2(N), RWORK(*)

      DOUBLE PRECISION DMACH(5)
      DATA DMACH(3) / 1.1101827117665D-16 /

      EXTERNAL MATVEC, MSOLVE
!
      ITER = 0
      IERR = 0
      IF( N.LT.1 ) THEN
         IERR = 3
         RETURN
      ENDIF
      TOLMIN = 500.0*DMACH(3)
      IF( TOL.LT.TOLMIN ) THEN
         TOL = TOLMIN
         IERR = 4
      ENDIF
!
      CALL MATVEC(N, X, R, NELT, IA, JA, A, ISYM)
      DO 10 I = 1, N
         V1(I)  = R(I) - B(I)
 10   CONTINUE
      CALL MSOLVE(N, V1, R, NELT, IA, JA, A, ISYM, RWORK, IWORK)
!
      IF( ISDCGS(N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MSOLVE,   &
     &     ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, R, R0, P, Q,  &
     &     U, V1, V2, RWORK, IWORK, AK, BK, BNRM, SOLNRM) .NE. 0 ) GO TO 200
      IF( IERR.NE.0 ) RETURN
!
!
      FUZZ = DMACH(3)**2
      DO 20 I = 1, N
         R0(I) = R(I)
 20   CONTINUE
      RHONM1 = 1.0
!
!
      DO 100 K=1,ITMAX
         ITER = K
!
         DDOT = 0.D0
         DO 15 I = 1,N
           DDOT = DDOT + R0(I)*R(I)
   15    CONTINUE
         RHON = DDOT

         IF( ABS(RHONM1).LT.FUZZ ) GOTO 998
         BK = RHON/RHONM1
         IF( ITER.EQ.1 ) THEN
            DO 30 I = 1, N
               U(I) = R(I)
               P(I) = R(I)
 30         CONTINUE
         ELSE
            DO 40 I = 1, N
               U(I) = R(I) + BK*Q(I)
               V1(I) = Q(I) + BK*P(I)
 40         CONTINUE
            DO 50 I = 1, N
               P(I) = U(I) + BK*V1(I)
 50         CONTINUE
         ENDIF
!
         CALL MATVEC(N, P, V2, NELT, IA, JA, A, ISYM)
         CALL MSOLVE(N, V2, V1, NELT, IA, JA, A, ISYM, RWORK, IWORK)

         DDOT = 0.D0
         DO 25 I = 1,N
           DDOT = DDOT + R0(I)*V1(I)
   25    CONTINUE
         SIGMA = DDOT

         IF( ABS(SIGMA).LT.FUZZ ) GOTO 999
         AK = RHON/SIGMA
         AKM = -AK
         DO 60 I = 1, N
            Q(I) = U(I) + AKM*V1(I)
 60      CONTINUE

         DO 70 I = 1, N
            V1(I) = U(I) + Q(I)
 70      CONTINUE

          DO 72 I = 1,N
             X(I) = X(I) + AKM*V1(I)
   72     CONTINUE
!                     -1
         CALL MATVEC(N, V1, V2, NELT, IA, JA, A, ISYM)
         CALL MSOLVE(N, V2, V1, NELT, IA, JA, A, ISYM, RWORK, IWORK)

          DO 73 I = 1,N
             R(I) = R(I) + AKM*V1(I)
   73     CONTINUE
!
         IF( ISDCGS(N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MSOLVE,   &
     &        ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, R, R0, P, Q,   &
     &        U, V1, V2, RWORK, IWORK, AK, BK, BNRM, SOLNRM) .NE. 0 ) GO TO 200
!
         RHONM1 = RHON
 100  CONTINUE
!
      ITER = ITMAX + 1
      IERR = 2
 200  RETURN
!
 998  IERR = 5
      RETURN
!
 999  IERR = 6
      RETURN
      END
!
!
!
!***************************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***************************************************************************
!
!
      SUBROUTINE DSLUCS(N, B, X, NELT, IA, JA, A, ISYM, ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, RWORK, LENW, IWORK, LENIW )
!
!
!
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX, ITER
      INTEGER IERR, IUNIT, LENW, IWORK(LENIW), LENIW
      DOUBLE PRECISION B(N), X(N), A(NELT), TOL, ERR, RWORK(LENW)
      EXTERNAL DSMV, DSLUI
      PARAMETER (LOCRB=1, LOCIB=11)
!
      IERR = 0
      IF( N.LT.1 .OR. NELT.LT.1 ) THEN
         IERR = 3
         RETURN
      ENDIF
      CALL DS2Y( N, NELT, IA, JA, A, ISYM, IUNIT )
!
      NL = 0
      NU = 0
      DO 20 ICOL = 1, N
         JBGN = JA(ICOL)+1
         JEND = JA(ICOL+1)-1
         IF( JBGN.LE.JEND ) THEN
            DO 10 J = JBGN, JEND
               IF( IA(J).GT.ICOL ) THEN
                  NL = NL + 1
                  IF( ISYM.NE.0 ) NU = NU + 1
               ELSE
                  NU = NU + 1
               ENDIF
 10         CONTINUE
         ENDIF
 20   CONTINUE
!
      LOCIL = LOCIB
      LOCJL = LOCIL + N+1
      LOCIU = LOCJL + NL
      LOCJU = LOCIU + NU
      LOCNR = LOCJU + N+1
      LOCNC = LOCNR + N
      LOCIW = LOCNC + N
!
      LOCL   = LOCRB
      LOCDIN = LOCL + NL
      LOCUU  = LOCDIN + N
      LOCR   = LOCUU + NU
      LOCR0  = LOCR + N
      LOCP   = LOCR0 + N
      LOCQ   = LOCP + N
      LOCU   = LOCQ + N
      LOCV1  = LOCU + N
      LOCV2  = LOCV1 + N
      LOCW   = LOCV2 + N
!
      CALL DCHKW( 'DSLUCS', LOCIW, LENIW, LOCW, LENW, IERR, ITER, ERR )
      IF( IERR.NE.0 ) RETURN
!
      IWORK(1) = LOCIL
      IWORK(2) = LOCJL
      IWORK(3) = LOCIU
      IWORK(4) = LOCJU
      IWORK(5) = LOCL
      IWORK(6) = LOCDIN
      IWORK(7) = LOCUU
      IWORK(9) = LOCIW
      IWORK(10) = LOCW
!
      CALL DSILUS( N, NELT, IA, JA, A, ISYM, NL, IWORK(LOCIL),  &
     &     IWORK(LOCJL), RWORK(LOCL), RWORK(LOCDIN), NU, IWORK(LOCIU),  &
     &     IWORK(LOCJU), RWORK(LOCUU), IWORK(LOCNR), IWORK(LOCNC) )
!
      CALL DCGS(N, B, X, NELT, IA, JA, A, ISYM, DSMV,  &
     &     DSLUI, ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT,  &
     &     RWORK(LOCR), RWORK(LOCR0), RWORK(LOCP),  &
     &     RWORK(LOCQ), RWORK(LOCU), RWORK(LOCV1), RWORK(LOCV2), RWORK, IWORK )
      RETURN
      END
!
!
!
!***************************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***************************************************************************
!
!
      FUNCTION ISDCGS(N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MSOLVE,   &
     &     ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, R, R0, P, Q, U,   &
     &     V1, V2, RWORK, IWORK, AK, BK, BNRM, SOLNRM)
!
!
!
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX
      INTEGER ITER, IERR, IUNIT, IWORK(*)
      DOUBLE PRECISION B(N), X(N), A(NELT), TOL, ERR, R(N), R0(N), P(N)
      DOUBLE PRECISION Q(N), U(N), V1(N), V2(N), RWORK(*)
      DOUBLE PRECISION AK, BK, BNRM, SOLNRM
      COMMON /SOLBLK/ SOLN(1)
      EXTERNAL MATVEC, MSOLVE
!
      ISDCGS = 0
!
      IF( ITOL.EQ.1 ) THEN
         IF(ITER .EQ. 0) BNRM = DNRM2(N, B, 1)
         CALL MATVEC(N, X, V2, NELT, IA, JA, A, ISYM )
         DO 5 I = 1, N
            V2(I) = V2(I) - B(I)
 5       CONTINUE
         ERR = DNRM2(N, V2, 1)/BNRM
      ELSE IF( ITOL.EQ.2 ) THEN
         IF(ITER .EQ. 0) THEN
            CALL MSOLVE(N, B, V2, NELT, IA, JA, A, ISYM, RWORK, IWORK)
            BNRM = DNRM2(N, V2, 1)
         ENDIF
         ERR = DNRM2(N, R, 1)/BNRM
      ELSE IF( ITOL.EQ.11 ) THEN
         IF(ITER .EQ. 0) SOLNRM = DNRM2(N, SOLN, 1)
         DO 10 I = 1, N
            V2(I) = X(I) - SOLN(I)
 10      CONTINUE
         ERR = DNRM2(N, V2, 1)/SOLNRM
      ELSE
         ERR = 1.0E10
         IERR = 3
      ENDIF
!
      IF(IUNIT .NE. 0) THEN
         IF( ITER.EQ.0 ) THEN
            WRITE(IUNIT,1000) N, ITOL
         ENDIF
         WRITE(IUNIT,1010) ITER, ERR, AK, BK
      ENDIF
      IF(ERR .LE. TOL) ISDCGS = 1
!
      RETURN
 1000 FORMAT(' Preconditioned BiConjugate Gradient Squared for N, ITOL = ',I5, I5,/, &
     &' ITER','   Error Estimate','            Alpha             Beta')
 1010 FORMAT(1X,I4,1X,E16.7,1X,E16.7,1X,E16.7)
      END
!
!
!
!***************************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***************************************************************************
!
!
      SUBROUTINE DSMV( N, X, Y, NELT, IA, JA, A, ISYM )
!
!
!
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM
      DOUBLE PRECISION A(NELT), X(N), Y(N)
!
      DO 10 I = 1, N
         Y(I) = 0.0D0
 10   CONTINUE
!
      DO 30 ICOL = 1, N
         IBGN = JA(ICOL)
         IEND = JA(ICOL+1)-1
         DO 20 I = IBGN, IEND
            Y(IA(I)) = Y(IA(I)) + A(I)*X(ICOL)
 20      CONTINUE
 30   CONTINUE
!
      IF( ISYM.EQ.1 ) THEN
!
         DO 50 IROW = 1, N
            JBGN = JA(IROW)+1
            JEND = JA(IROW+1)-1
            IF( JBGN.GT.JEND ) GOTO 50
            DO 40 J = JBGN, JEND
               Y(IROW) = Y(IROW) + A(J)*X(IA(J))
 40         CONTINUE
 50      CONTINUE
      ENDIF
      RETURN
      END
!
!
!
!***************************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***************************************************************************
!
!
      SUBROUTINE DSMTV( N, X, Y, NELT, IA, JA, A, ISYM )
!
!
!
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM
      DOUBLE PRECISION X(N), Y(N), A(NELT)
!
      DO 10 I = 1, N
         Y(I) = 0.0D0
 10   CONTINUE
!
      DO 30 IROW = 1, N
         IBGN = JA(IROW)
         IEND = JA(IROW+1)-1
         DO 20 I = IBGN, IEND
            Y(IROW) = Y(IROW) + A(I)*X(IA(I))
 20      CONTINUE
 30   CONTINUE
!
      IF( ISYM.EQ.1 ) THEN
         DO 50 ICOL = 1, N
            JBGN = JA(ICOL)+1
            JEND = JA(ICOL+1)-1
            IF( JBGN.GT.JEND ) GOTO 50
            DO 40 J = JBGN, JEND
               Y(IA(J)) = Y(IA(J)) + A(J)*X(ICOL)
 40         CONTINUE
 50      CONTINUE
      ENDIF
      RETURN
      END
!
!
!
!***************************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***************************************************************************
!
!
      SUBROUTINE DSLUI(N, B, X, NELT, IA, JA, A, ISYM, RWORK, IWORK)
!
!
!
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, IWORK(*)
      DOUBLE PRECISION B(N), X(N), A(NELT), RWORK(*)
!
      LOCIL = IWORK(1)
      LOCJL = IWORK(2)
      LOCIU = IWORK(3)
      LOCJU = IWORK(4)
      LOCL = IWORK(5)
      LOCDIN = IWORK(6)
      LOCU = IWORK(7)
!
      CALL DSLUI2(N, B, X, IWORK(LOCIL), IWORK(LOCJL), RWORK(LOCL), RWORK(LOCDIN), IWORK(LOCIU), IWORK(LOCJU), RWORK(LOCU) )
!
      RETURN
      END
!
!
!
!***************************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***************************************************************************
!
!
      SUBROUTINE DSLUI2(N, B, X, IL, JL, L, DINV, IU, JU, U )
!
!
!
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER N, IL(*), JL(*), IU(*), JU(*)
      DOUBLE PRECISION B(N), X(N), L(*), DINV(N), U(*)
!
      DO 10 I = 1, N
         X(I) = B(I)
 10   CONTINUE
      DO 30 IROW = 2, N
         JBGN = IL(IROW)
         JEND = IL(IROW+1)-1
         IF( JBGN.LE.JEND ) THEN
            DO 20 J = JBGN, JEND
               X(IROW) = X(IROW) - L(J)*X(JL(J))
 20         CONTINUE
         ENDIF
 30   CONTINUE
!
      DO 40 I=1,N
         X(I) = X(I)*DINV(I)
 40   CONTINUE
!
      DO 60 ICOL = N, 2, -1
         JBGN = JU(ICOL)
         JEND = JU(ICOL+1)-1
         IF( JBGN.LE.JEND ) THEN
            DO 50 J = JBGN, JEND
               X(IU(J)) = X(IU(J)) - U(J)*X(ICOL)
 50         CONTINUE
         ENDIF
 60   CONTINUE
!
      RETURN
      END
!
!
!
!***************************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***************************************************************************
!
!
      SUBROUTINE DSLUTI(N, B, X, NELT, IA, JA, A, ISYM, RWORK, IWORK)
!
!
!
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, IWORK(*)
      DOUBLE PRECISION B(N), X(N), A(N), RWORK(*)
!
      LOCIL = IWORK(1)
      LOCJL = IWORK(2)
      LOCIU = IWORK(3)
      LOCJU = IWORK(4)
      LOCL = IWORK(5)
      LOCDIN = IWORK(6)
      LOCU = IWORK(7)
!
      CALL DSLUI4(N, B, X, IWORK(LOCIL), IWORK(LOCJL), RWORK(LOCL), RWORK(LOCDIN), IWORK(LOCIU), IWORK(LOCJU), RWORK(LOCU))
!
      RETURN
      END
!
!
!
!***************************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***************************************************************************
!
!
      SUBROUTINE DSLUI4(N, B, X, IL, JL, L, DINV, IU, JU, U )

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER N, IL(*), JL(*), IU(*), JU(*)
      DOUBLE PRECISION B(N), X(N), L(*), DINV(N), U(*)
!
      DO 10 I=1,N
         X(I) = B(I)
 10   CONTINUE
!
      DO 80 IROW = 2, N
         JBGN = JU(IROW)
         JEND = JU(IROW+1) - 1
         IF( JBGN.LE.JEND ) THEN
            DO 70 J = JBGN, JEND
               X(IROW) = X(IROW) - U(J)*X(IU(J))
 70         CONTINUE
         ENDIF
 80   CONTINUE
!
      DO 90 I = 1, N
         X(I) = X(I)*DINV(I)
 90   CONTINUE
!
      DO 110 ICOL = N, 2, -1
         JBGN = IL(ICOL)
         JEND = IL(ICOL+1) - 1
         IF( JBGN.LE.JEND ) THEN
            DO 100 J = JBGN, JEND
               X(JL(J)) = X(JL(J)) - L(J)*X(ICOL)
 100        CONTINUE
         ENDIF
 110  CONTINUE
      RETURN
      END
!
!
!
!***************************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***************************************************************************
!
!
      SUBROUTINE DSILUS(N, NELT, IA, JA, A, ISYM, NL, IL, JL, L, DINV, NU, IU, JU, U, NROW, NCOL)
!
!
!
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!
      INTEGER N, NELT,IA(NELT),JA(NELT),ISYM,NL,IL(NL+2),JL(NL+2)
      INTEGER NU, IU(NU+2), JU(NU+2), NROW(N), NCOL(N)
!
      DOUBLE PRECISION A(NELT), L(NL), DINV(N), U(NU)
!
      DO 10 I=1,N
         NROW(I) = 0
         NCOL(I) = 0
 10   CONTINUE
      DO 30 ICOL = 1, N
         JBGN = JA(ICOL)+1
         JEND = JA(ICOL+1)-1
         IF( JBGN.LE.JEND ) THEN
            DO 20 J = JBGN, JEND
               IF( IA(J).LT.ICOL ) THEN
                  NCOL(ICOL) = NCOL(ICOL) + 1
               ELSE
                  NROW(IA(J)) = NROW(IA(J)) + 1
                  IF( ISYM.NE.0 ) NCOL(IA(J)) = NCOL(IA(J)) + 1
               ENDIF
 20         CONTINUE
         ENDIF
 30   CONTINUE
      JU(1) = 1
      IL(1) = 1
!
      DO 40 ICOL = 1, N
         IL(ICOL+1) = IL(ICOL) + NROW(ICOL)
         JU(ICOL+1) = JU(ICOL) + NCOL(ICOL)
         NROW(ICOL) = IL(ICOL)
         NCOL(ICOL) = JU(ICOL)
 40   CONTINUE
!
      DO 60 ICOL = 1, N
         DINV(ICOL) = A(JA(ICOL))
         JBGN = JA(ICOL)+1
         JEND = JA(ICOL+1)-1
         IF( JBGN.LE.JEND ) THEN
            DO 50 J = JBGN, JEND
               IROW = IA(J)
               IF( IROW.LT.ICOL ) THEN
                  IU(NCOL(ICOL)) = IROW
                  U(NCOL(ICOL)) = A(J)
                  NCOL(ICOL) = NCOL(ICOL) + 1
               ELSE
                  JL(NROW(IROW)) = ICOL
                  L(NROW(IROW)) = A(J)
                  NROW(IROW) = NROW(IROW) + 1
                  IF( ISYM.NE.0 ) THEN
                     IU(NCOL(IROW)) = ICOL
                     U(NCOL(IROW)) = A(J)
                     NCOL(IROW) = NCOL(IROW) + 1
                  ENDIF
               ENDIF
 50         CONTINUE
         ENDIF
 60   CONTINUE
!
      DO 110 K = 2, N
         JBGN = JU(K)
         JEND = JU(K+1)-1
         IF( JBGN.LT.JEND ) THEN
            DO 80 J = JBGN, JEND-1
               DO 70 I = J+1, JEND
                  IF( IU(J).GT.IU(I) ) THEN
                     ITEMP = IU(J)
                     IU(J) = IU(I)
                     IU(I) = ITEMP
                     TEMP = U(J)
                     U(J) = U(I)
                     U(I) = TEMP
                  ENDIF
 70            CONTINUE
 80         CONTINUE
         ENDIF
         IBGN = IL(K)
         IEND = IL(K+1)-1
         IF( IBGN.LT.IEND ) THEN
            DO 100 I = IBGN, IEND-1
               DO 90 J = I+1, IEND
                  IF( JL(I).GT.JL(J) ) THEN
                     JTEMP = JU(I)
                     JU(I) = JU(J)
                     JU(J) = JTEMP
                     TEMP = L(I)
                     L(I) = L(J)
                     L(J) = TEMP
                  ENDIF
 90            CONTINUE
 100        CONTINUE
         ENDIF
 110  CONTINUE
!
      DO 300 I=2,N
!
         INDX1 = IL(I)
         INDX2 = IL(I+1) - 1
         IF(INDX1 .GT. INDX2) GO TO 200
         DO 190 INDX=INDX1,INDX2
            IF(INDX .EQ. INDX1) GO TO 180
            INDXR1 = INDX1
            INDXR2 = INDX - 1
            INDXC1 = JU(JL(INDX))
            INDXC2 = JU(JL(INDX)+1) - 1
            IF(INDXC1 .GT. INDXC2) GO TO 180
 160        KR = JL(INDXR1)
 170        KC = IU(INDXC1)
            IF(KR .GT. KC) THEN
               INDXC1 = INDXC1 + 1
               IF(INDXC1 .LE. INDXC2) GO TO 170
            ELSEIF(KR .LT. KC) THEN
               INDXR1 = INDXR1 + 1
               IF(INDXR1 .LE. INDXR2) GO TO 160
            ELSEIF(KR .EQ. KC) THEN
               L(INDX) = L(INDX) - L(INDXR1)*DINV(KC)*U(INDXC1)
               INDXR1 = INDXR1 + 1
               INDXC1 = INDXC1 + 1
               IF(INDXR1 .LE. INDXR2 .AND. INDXC1 .LE. INDXC2) GO TO 160
            ENDIF
 180        L(INDX) = L(INDX)/DINV(JL(INDX))
 190     CONTINUE
!
 200     INDX1 = JU(I)
         INDX2 = JU(I+1) - 1
         IF(INDX1 .GT. INDX2) GO TO 260
         DO 250 INDX=INDX1,INDX2
            IF(INDX .EQ. INDX1) GO TO 240
            INDXC1 = INDX1
            INDXC2 = INDX - 1
            INDXR1 = IL(IU(INDX))
            INDXR2 = IL(IU(INDX)+1) - 1
            IF(INDXR1 .GT. INDXR2) GO TO 240
 210        KR = JL(INDXR1)
 220        KC = IU(INDXC1)
            IF(KR .GT. KC) THEN
               INDXC1 = INDXC1 + 1
               IF(INDXC1 .LE. INDXC2) GO TO 220
            ELSEIF(KR .LT. KC) THEN
               INDXR1 = INDXR1 + 1
               IF(INDXR1 .LE. INDXR2) GO TO 210
            ELSEIF(KR .EQ. KC) THEN
               U(INDX) = U(INDX) - L(INDXR1)*DINV(KC)*U(INDXC1)
               INDXR1 = INDXR1 + 1
               INDXC1 = INDXC1 + 1
               IF(INDXR1 .LE. INDXR2 .AND. INDXC1 .LE. INDXC2) GO TO 210
            ENDIF
 240        U(INDX) = U(INDX)/DINV(IU(INDX))
 250     CONTINUE
!
 260     INDXR1 = IL(I)
         INDXR2 = IL(I+1) - 1
         IF(INDXR1 .GT. INDXR2) GO TO 300
         INDXC1 = JU(I)
         INDXC2 = JU(I+1) - 1
         IF(INDXC1 .GT. INDXC2) GO TO 300
 270     KR = JL(INDXR1)
 280     KC = IU(INDXC1)
         IF(KR .GT. KC) THEN
            INDXC1 = INDXC1 + 1
            IF(INDXC1 .LE. INDXC2) GO TO 280
         ELSEIF(KR .LT. KC) THEN
            INDXR1 = INDXR1 + 1
            IF(INDXR1 .LE. INDXR2) GO TO 270
         ELSEIF(KR .EQ. KC) THEN
            DINV(I) = DINV(I) - L(INDXR1)*DINV(KC)*U(INDXC1)
            INDXR1 = INDXR1 + 1
            INDXC1 = INDXC1 + 1
            IF(INDXR1 .LE. INDXR2 .AND. INDXC1 .LE. INDXC2) GO TO 270
         ENDIF
!
 300  CONTINUE
!
      DO 430 I=1,N
         DINV(I) = 1./DINV(I)
 430  CONTINUE
!
      RETURN
      END
!
!
!
!***************************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***************************************************************************
!
!
      SUBROUTINE DCHKW( NAME, LOCIW, LENIW, LOCW, LENW, IERR, ITER, ERR )
!
!
!
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*(*) NAME
      CHARACTER*72 MESG
      INTEGER LOCIW, LENIW, LOCW, LENW, IERR, ITER
      DOUBLE PRECISION ERR

      DOUBLE PRECISION DMACH(5)
      DATA DMACH(2) / 1.79769313486231D+308 /
!
      IERR = 0
      IF( LOCIW.GT.LENIW ) THEN
         IERR = 1
         ITER = 0
         ERR = DMACH(2)
         MESG = NAME // ': INTEGER work array too short. '//' IWORK needs i1: have allocated i2.'
      ENDIF
!
      IF( LOCW.GT.LENW ) THEN
         IERR = 1
         ITER = 0
         ERR = DMACH(2)
         MESG = NAME // ': DOUBLE PRECISION work array too short. '//' RWORK needs i1: have allocated i2.'
      ENDIF

      RETURN
      END
!
!
!
!***************************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***************************************************************************
!
!
      DOUBLE PRECISION FUNCTION DNRM2 ( N, DX, INCX)
      INTEGER          NEXT
      DOUBLE PRECISION   DX(*), CUTLO, CUTHI, HITEST, SUM, XMAX,ZERO,ONE
      DATA   ZERO, ONE /0.0D0, 1.0D0/
      DATA   CUTLO, CUTHI / 4.441d-16,  1.304d19 /
!
!
!
      IF(N .GT. 0) GO TO 10
         DNRM2  = ZERO
         GO TO 300
!
   10 ASSIGN 30 TO NEXT
      SUM = ZERO
      NN = N * INCX
      I = 1
   20    GO TO NEXT,(30, 50, 70, 110)
   30 IF( DABS(DX(I)) .GT. CUTLO) GO TO 85
      ASSIGN 50 TO NEXT
      XMAX = ZERO
!
!
   50 IF( DX(I) .EQ. ZERO) GO TO 200
      IF( DABS(DX(I)) .GT. CUTLO) GO TO 85
!
      ASSIGN 70 TO NEXT
      GO TO 105
!
!
  100 I = J
      ASSIGN 110 TO NEXT
      SUM = (SUM / DX(I)) / DX(I)
  105 XMAX = DABS(DX(I))
      GO TO 115
!
   70 IF( DABS(DX(I)) .GT. CUTLO ) GO TO 75
!
  110 IF( DABS(DX(I)) .LE. XMAX ) GO TO 115
         SUM = ONE + SUM * (XMAX / DX(I))**2
         XMAX = DABS(DX(I))
         GO TO 200
!
  115 SUM = SUM + (DX(I)/XMAX)**2
      GO TO 200
!
   75 SUM = (SUM * XMAX) * XMAX
!
!
   85 HITEST = CUTHI/FLOAT( N )
!
      DO 95 J =I,NN,INCX
      IF(DABS(DX(J)) .GE. HITEST) GO TO 100
   95    SUM = SUM + DX(J)**2
      DNRM2 = DSQRT( SUM )
      GO TO 300
!
  200 CONTINUE
      I = I + INCX
      IF ( I .LE. NN ) GO TO 20
!
      DNRM2 = XMAX * DSQRT(SUM)
  300 CONTINUE
      RETURN
      END
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
      SUBROUTINE QS2I1D( IA, JA, A, N, KFLAG, IUNIT )
!
!
!
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION IL(21),IU(21)
      INTEGER   IA(N),JA(N),IT,IIT,JT,JJT
      DOUBLE PRECISION A(N), TA, TTA
!
      NN = N
      IF (NN.LT.1) THEN
         WRITE(IUNIT,6100)
 6100 FORMAT(/,'QS2I1D- the number of values to be sorted was NOT POSITIVE.')
         RETURN
      ENDIF
      IF( N.EQ.1 ) RETURN
      KK = IABS(KFLAG)
      IF ( KK.NE.1 ) THEN
         WRITE(IUNIT,6101)
 6101 FORMAT(/,'QS2I1D- the sort control parameter, k, was not 1 OR -1.')
         RETURN
      ENDIF
!
      IF( KFLAG.LT.1 ) THEN
         DO 20 I=1,NN
            IA(I) = -IA(I)
 20      CONTINUE
      ENDIF
!
      M = 1
      I = 1
      J = NN
      R = .375
 210  IF( R.LE.0.5898437 ) THEN
         R = R + 3.90625E-2
      ELSE
         R = R-.21875
      ENDIF
 225  K = I
!
!
      IJ = I + IDINT( DBLE(J-I)*R )
      IT = IA(IJ)
      JT = JA(IJ)
      TA = A(IJ)
!
!
      IF( IA(I).GT.IT ) THEN
         IA(IJ) = IA(I)
         IA(I)  = IT
         IT     = IA(IJ)
         JA(IJ) = JA(I)
         JA(I)  = JT
         JT     = JA(IJ)
         A(IJ)  = A(I)
         A(I)   = TA
         TA     = A(IJ)
      ENDIF
      L=J
!
!
      IF( IA(J).LT.IT ) THEN
         IA(IJ) = IA(J)
         IA(J)  = IT
         IT     = IA(IJ)
         JA(IJ) = JA(J)
         JA(J)  = JT
         JT     = JA(IJ)
         A(IJ)  = A(J)
         A(J)   = TA
         TA     = A(IJ)
!
!
         IF ( IA(I).GT.IT ) THEN
            IA(IJ) = IA(I)
            IA(I)  = IT
            IT     = IA(IJ)
            JA(IJ) = JA(I)
            JA(I)  = JT
            JT     = JA(IJ)
            A(IJ)  = A(I)
            A(I)   = TA
            TA     = A(IJ)
         ENDIF
      ENDIF
!
!
  240 L=L-1
      IF( IA(L).GT.IT ) GO TO 240
!
!
  245 K=K+1
      IF( IA(K).LT.IT ) GO TO 245
!
!
      IF( K.LE.L ) THEN
         IIT   = IA(L)
         IA(L) = IA(K)
         IA(K) = IIT
         JJT   = JA(L)
         JA(L) = JA(K)
         JA(K) = JJT
         TTA   = A(L)
         A(L)  = A(K)
         A(K)  = TTA
         GOTO 240
      ENDIF
!
!
      IF( L-I.GT.J-K ) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 260
!
!
  255 M = M-1
      IF( M.EQ.0 ) GO TO 300
      I = IL(M)
      J = IU(M)
  260 IF( J-I.GE.1 ) GO TO 225
      IF( I.EQ.J ) GO TO 255
      IF( I.EQ.1 ) GO TO 210
      I = I-1
  265 I = I+1
      IF( I.EQ.J ) GO TO 255
      IT = IA(I+1)
      JT = JA(I+1)
      TA =  A(I+1)
      IF( IA(I).LE.IT ) GO TO 265
      K=I
  270 IA(K+1) = IA(K)
      JA(K+1) = JA(K)
      A(K+1)  =  A(K)
      K = K-1
      IF( IT.LT.IA(K) ) GO TO 270
      IA(K+1) = IT
      JA(K+1) = JT
      A(K+1)  = TA
      GO TO 265
!
!
  300 IF( KFLAG.LT.1 ) THEN
         DO 310 I=1,NN
            IA(I) = -IA(I)
 310     CONTINUE
      ENDIF
      RETURN
      END
!
!
!
!***************************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***************************************************************************
!
!
      SUBROUTINE DS2Y(N, NELT, IA, JA, A, ISYM,IUNIT )
!
!
!
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM
      DOUBLE PRECISION A(NELT)
!
      IF( JA(N+1).EQ.NELT+1 ) RETURN
      CALL QS2I1D( JA, IA, A, NELT, 1,IUNIT )
      JA(1) = 1
      DO 20 ICOL = 1, N-1
         DO 10 J = JA(ICOL)+1, NELT
            IF( JA(J).NE.ICOL ) THEN
               JA(ICOL+1) = J
               GOTO 20
            ENDIF
 10      CONTINUE
 20   CONTINUE
      JA(N+1) = NELT+1
!
      JA(N+2) = 0
!
      DO 70 ICOL = 1, N
         IBGN = JA(ICOL)
         IEND = JA(ICOL+1)-1
         DO 30 I = IBGN, IEND
            IF( IA(I).EQ.ICOL ) THEN
               ITEMP = IA(I)
               IA(I) = IA(IBGN)
               IA(IBGN) = ITEMP
               TEMP = A(I)
               A(I) = A(IBGN)
               A(IBGN) = TEMP
               GOTO 40
            ENDIF
 30      CONTINUE
 40      IBGN = IBGN + 1
         IF( IBGN.LT.IEND ) THEN
            DO 60 I = IBGN, IEND
               DO 50 J = I+1, IEND
                  IF( IA(I).GT.IA(J) ) THEN
                     ITEMP = IA(I)
                     IA(I) = IA(J)
                     IA(J) = ITEMP
                     TEMP = A(I)
                     A(I) = A(J)
                     A(J) = TEMP
                  ENDIF
 50            CONTINUE
 60         CONTINUE
         ENDIF
 70   CONTINUE
      RETURN
      END
!
!
!***************************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***************************************************************************
!
!
        SUBROUTINE LIS(N, R, work, Non0, RowNum, ColNum, CO)
!
!
!
!#include "omp.h"
        implicit none
              INTEGER N, Non0, RowNum(Non0), ColNum(Non0), Row_New(Non0), Col_New(Non0)!!!!,index(N)-----save memory
              DOUBLE PRECISION R(N), work(N),CO(Non0),start_time,end_time
!
#include "lisf.h"
              LIS_MATRIX A
              LIS_VECTOR b,x
              LIS_SOLVER solver
              LIS_INTEGER i,gn,ln,ierr
              LIS_INTEGER matrix_type,comm_world
              LIS_INTEGER omp_get_num_procs,omp_get_max_threads
!
        Row_New = RowNum-1
        Col_New = ColNum-1
!
!
!
              call lis_initialize(ierr)
!
              matrix_type = LIS_MATRIX_COO
              comm_world = LIS_COMM_WORLD
              ln = 0
!
#ifdef _OPENMP
!         write(*,*) 'max number of threads = ',omp_get_num_procs()
!         write(*,*) 'number of threads = ', omp_get_max_threads()
#endif
!
              call lis_matrix_create(comm_world,A,ierr)
              call lis_matrix_set_size(A,ln,N,ierr)
              call lis_matrix_set_type(A,matrix_type,ierr)
              call lis_matrix_set_coo(Non0,Row_New,Col_New,CO,A,ierr)
              call lis_matrix_assemble(A,ierr)
              call lis_vector_duplicate(A,b,ierr)
              call lis_vector_duplicate(A,x,ierr)
!
!-------------------alternative section for assigning vector values---------------
!        do i=1,N
!                index(i)=i
!        enddo
!
!        call lis_vector_set_values(LIS_INS_VALUE,N,index,R,b,ierr)
!        call lis_vector_set_values(LIS_INS_VALUE,N,index,work,x,ierr)
!
!-------------------alternative section for assigning vector values---------------
        do i=1,N
            call lis_vector_set_value(LIS_INS_VALUE,i,R(i),b,ierr)
        enddo
!
        call lis_vector_set_all(0.0d0,x,ierr)
        call lis_solver_create(solver,ierr)
        call lis_solver_set_option('-i 5 -ell 2 -p 9',solver,ierr)
        call lis_solver_set_option('-tol 1.0e-8',solver,ierr)
!

        call lis_solve(A,b,x,solver,ierr)
!        call lis_solver_get_precon(solver,gn,ierr)
!        PRINT *, 'preconditioner used is  ',gn
!
        call lis_vector_get_values(x,1,N,work,ierr)
              call lis_solver_destroy(solver,ierr)
        call lis_matrix_unset(A,ierr)
        call lis_matrix_destroy(A,ierr)
        call lis_vector_destroy(x,ierr)
        call lis_vector_destroy(b,ierr)
!
        call lis_finalize(ierr)
              RETURN
              END

!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!
!
!
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
