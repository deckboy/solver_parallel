!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>                                                                     >
!>                                                                     >
!>      T_UtilityFunctions.f95: Code unit including the routines       >
!>                describing general utility functions                 >
!>                                                                     >
!>                                                                     >
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
      MODULE Utility_Functions
!
         PRIVATE
!
         PUBLIC :: Polynomial, ExponSeries, SineSeries,                  &
      &            Integral_poly, Integral_ExponSrs,  Integral_SineSrs,  &
      &            N_Character, Cubic_Equation_Roots
!
!
!
         CONTAINS
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
         REAL(KIND = 8) FUNCTION Polynomial(argument,n_poly,A)
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*                 Routine for computing a polynomial                  *
!*                                                                     *
!*                     Version 1.0 - June 6, 2007                      *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
            IMPLICIT NONE
!
! -------------
! ......... Integer input variables
! -------------
!
            INTEGER, INTENT(IN) :: n_poly
!
! -------------
! ......... Integer variables
! -------------
!
            INTEGER :: i
!
! -------------
! ......... Real input variables
! -------------
!
            REAL(KIND = 8), INTENT(IN)                      :: argument
            REAL(KIND = 8), INTENT(IN), DIMENSION(0:n_poly) :: A
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=> Main body of Q_Polynomial
!
!
            Polynomial = A(n_poly)
!
            DO i = n_poly-1,0,-1
!
               Polynomial = Polynomial*argument + A(i)
!
            END DO
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=> End of Polynomial
!
!
            RETURN
!
         END FUNCTION Polynomial
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
         REAL(KIND = 8) FUNCTION Integral_poly(argument, A)
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*         Routine for computing the integral of a polynomial          *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
            IMPLICIT NONE
!
! -------------
! ......... Integer input variables
! -------------
!
            INTEGER :: n_poly
!
! -------------
! ......... Integer variables
! -------------
!
            INTEGER :: i
!
! -------------
! ......... Real input variables
! -------------
!
            REAL(KIND = 8), INTENT(IN)                :: argument
            REAL(KIND = 8), INTENT(IN), DIMENSION(0:) :: A
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of Integral_poly
!
!
            n_poly = SIZE(A) - 1
            Integral_poly = A(n_poly)/(DBLE(n_poly) + 1.0d0)
!
            DO i = n_poly-1,0,-1
!
               Integral_poly = Integral_poly*argument + A(i)/(DBLE(i) + 1.0d0)
!
            END DO
!
            Integral_poly = Integral_poly*argument
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End of Integral_poly
!
!
            RETURN
!
         END FUNCTION Integral_poly
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
         REAL(KIND = 8) FUNCTION ExponSeries(argument,n_expon,A,B,C,acceptable)
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*             Routine for computing an exponential series             *
!*                                                                     *
!*                     Version 1.0 - June 6, 2007                      *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
            IMPLICIT NONE
!
! -------------
! ......... Integer input variables
! -------------
!
            INTEGER, INTENT(IN) :: n_expon
!
! -------------
! ......... Real input variables
! -------------
!
            REAL(KIND = 8), INTENT(IN)                       :: argument
            REAL(KIND = 8), INTENT(IN), DIMENSION(0:n_expon) :: A,B,C
!
! -------------
! ......... Logical variable
! -------------
!
            LOGICAL, INTENT(OUT) :: acceptable
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of ExponSeries
!
!
            IF(ANY((B + C*argument) > HUGE(1.0d0))) THEN
               acceptable = .FALSE.
            ELSE
               acceptable = .TRUE.
               ExponSeries = SUM(A*(EXP(B + C*argument)))
            END IF
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End of ExponSeries
!
!
            RETURN
!
         END FUNCTION ExponSeries
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
         REAL(KIND = 8) FUNCTION Integral_ExponSrs(argument,n_expon,A,B,C,acceptable)
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*     Routine for computing the integral of an exponential series     *
!*                                                                     *
!*                     Version 1.0 - June 6, 2007                      *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
            IMPLICIT NONE
!
! -------------
! ......... Integer input variables
! -------------
!
            INTEGER, INTENT(IN) :: n_expon
!
! -------------
! ......... Real input variables
! -------------
!
            REAL(KIND = 8), INTENT(IN)                       :: argument
            REAL(KIND = 8), INTENT(IN), DIMENSION(0:n_expon) :: A,B,C
!
! -------------
! ......... Logical variable
! -------------
!
            LOGICAL, INTENT(OUT) :: acceptable
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of Integral_ExponSrs
!
!
            IF(ANY((B + C*argument) > HUGE(1.0d0))) THEN
               acceptable = .FALSE.
            ELSE
               acceptable = .TRUE.
               Integral_ExponSrs = SUM((A/C)*(EXP(B + C*argument)))
            END IF
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End of Integral_ExponSrs
!
!
            RETURN
!
         END FUNCTION Integral_ExponSrs
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
         REAL(KIND = 8) FUNCTION SineSeries(argument,n_sine,A,B,C)
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*                Routine for computing a sine series                  *
!*                                                                     *
!*                     Version 1.0 - June 6, 2007                      *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
            IMPLICIT NONE
!
! -------------
! ......... Integer input variables
! -------------
!
            INTEGER, INTENT(IN) :: n_sine
!
! -------------
! ......... Real input variables
! -------------
!
            REAL(KIND = 8), INTENT(IN)                      :: argument
            REAL(KIND = 8), INTENT(IN), DIMENSION(0:n_sine) :: A,B,C
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=> Main body of SineSeries
!
!
            SineSeries = SUM(A*(SIN(B + C*argument)))
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=> End of SineSeries
!
!
            RETURN
!
         END FUNCTION SineSeries
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
         REAL(KIND = 8) FUNCTION Integral_SineSrs(argument,n_sine,A,B,C)
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*         Routine for computing the integral of a sine series         *
!*                                                                     *
!*                     Version 1.0 - June 6, 2007                      *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
            IMPLICIT NONE
!
! -------------
! ......... Integer input variables
! -------------
!
            INTEGER, INTENT(IN) :: n_sine
!
! -------------
! ......... Real input variables
! -------------
!
            REAL(KIND = 8), INTENT(IN)                      :: argument
            REAL(KIND = 8), INTENT(IN), DIMENSION(0:n_sine) :: A,B,C
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=> Main body of Integral_SineSrs
!
!
            Integral_SineSrs = SUM((-A/C)*(COS(B + C*argument)))
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=> End of Integral_SineSrs
!
!
            RETURN
!
         END FUNCTION Integral_SineSrs
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
      CHARACTER(LEN=1) FUNCTION N_Character(N)
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*         ROUTINE FOR ASSIGNING CCHARACTER VALUES TO NUMBERS          *
!*                          IN ELEMENT NAMES                           *
!*                                                                     *
!*                     Version 1.00, July 14, 2007                     *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
      IMPLICIT NONE
!
! ----------
! ... Integer variables
! ----------
!
      INTEGER, INTENT(IN) :: N
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  BEGIN N_Character
!
!
      SELECT CASE (N)
      CASE (0)
         N_Character = '0'
      CASE (1)
         N_Character = '1'
      CASE (2)
         N_Character = '2'
      CASE (3)
         N_Character = '3'
      CASE (4)
         N_Character = '4'
      CASE (5)
         N_Character = '5'
      CASE (6)
         N_Character = '6'
      CASE (7)
         N_Character = '7'
      CASE (8)
         N_Character = '8'
      CASE (9)
         N_Character = '9'
      END SELECT
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  END N_Character
!
!
      RETURN
!
      END FUNCTION N_Character
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
      SUBROUTINE Cubic_Equation_Roots(A2,A1,A0,XMX,XMN,ixx)
!
!*************************************************************************
!*************************************************************************
!*                                                                       *
!*                     ROOTS OF THE CUBIC EQUATION                       *
!*                        X^3+A2*X^2+A1*X+A0=0                           *
!*                                                                       *
!*************************************************************************
!*************************************************************************
!
      IMPLICIT NONE
!
! -------
! ... Double precision parameters
! -------
!
      REAL(KIND = 8), PARAMETER :: PI=3.14159265358979324D0
!
! -------
! ... Double precision variables
! -------
!
      REAL(KIND = 8), INTENT(IN)  :: A2,A1,A0
      REAL(KIND = 8), INTENT(OUT) :: XMX,XMN
!
      REAL(KIND = 8) :: G,H,GH,S,TT,XA,X1,X2,X3,THETA
      REAL(KIND = 8) :: X1mn,X2mn,X3mn
!
! -------
! ... Integer variables
! -------
!
      INTEGER, INTENT(OUT) :: ixx
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of Cubic_Equation_Roots
!
! -------
! ... Intermediate variables
! -------
!
      G  =  A1-(A2*A2/3.0d0)
      H  = (2.0d0*A2*A2*A2-9.0d0*A2*A1+2.7d1*A0)/2.7d1
      GH = (G*G*G)/2.7d1+(H*H)/4.0d0
!
!
!
      IF_NoRoots: IF (GH > 0.0d0) THEN
!
! -------------
! ...... Single Real Root
! -------------
!
         S  = -5.0d-1*H+SQRT(GH)
         TT = -5.0d-1*H-SQRT(GH)
!
         IF (S >= 0.0d0) THEN
            S = S**(1.0d0/3.0d0)
         ELSE
            S = -((-S)**(1.0d0/3.0d0))
         END IF
!
         IF (TT >= 0.0d0) THEN
            TT = TT**(1.0d0/3.0d0)
         ELSE
            TT = -((-TT)**(1.0d0/3.0d0))
         END IF
!
         XMX = S+TT-A2/3.0d0
         XMN = XMX
         ixx = 1
!
      ELSE
!
! -------------
! ...... Three Real Roots
! -------------
!
         XA    = -5.0d-1*H/SQRT(-(G*G*G/2.7d1))
         THETA = (PI/2.0d0-ATAN(XA/(SQRT(1.0d0-XA*XA))))/3.0d0
!
         X1 = 2.0d0*SQRT(-G/3.0d0)*COS(THETA)
         X2 = 2.0d0*SQRT(-G/3.0d0)*COS(THETA+2.0d0*PI/3.0d0)
         X3 = 2.0d0*SQRT(-G/3.0d0)*COS(THETA+4.0d0*PI/3.0d0)
!
         XMX  = MAX(X1,X2,X3)-A2/3.0d0
!
         X1mn = X1-A2/3.0d0
         X2mn = X2-A2/3.0d0
         X3mn = X3-A2/3.0d0
!
         IF(X1mn < 0.0d0) X1mn = 1.0d10
         IF(X2mn < 0.0d0) X2mn = 1.0d10
         IF(X3mn < 0.0d0) X3mn = 1.0d10
!
         XMN  = MIN(X1mn,X2mn,X3mn)
         ixx  = 3
!
      END IF IF_NoRoots
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of Cubic_Equation_Roots
!
!
      RETURN
!
      END SUBROUTINE Cubic_Equation_Roots
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
      END MODULE Utility_Functions
!
!
!
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
