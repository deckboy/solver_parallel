!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>                                                                     >
!>                                                                     >
!>            GAS_EOS.f95: Code unit of GAS EOS,                       >
!>       the corresponding Jacobian, and associated procedures         >
!>                                                                     >
!>                                                                     >
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  Segment that describes the equation of state of the gas fluid,
!  assigns initial conditions, computes the thermophysical properties of
!  the gas-bearing medium, and determines phase changes and the state of
!  the system from the various possible options.  This segment also
!  includes the procedures that computes the elements of the Jacobian
!  matrix for the Newton-Raphson iteration.
!
!
      SUBROUTINE Equation_of_State
!
! ...... Modules to be used
!
         USE Universal_Parameters
         USE Basic_Parameters
         USE EOS_Default_Parameters
!
         USE General_External_File_Units
         USE General_Control_Parameters
!
         USE Grid_Geometry
         USE Element_Attributes
!
         USE Geologic_Media_Parameters
!
         USE Solution_Matrix_Arrays, ONLY: X, DX
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*    DETERMINES PHASE TRANSITIONS AND COMPUTES ALL THERMOPHYSICAL     *
!*           PROPERTIES OF TWO COMPONENTS IN 1,2, OR 3 PHASES          *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
      IMPLICIT NONE
!
! -------
! ... Double precision variables
! -------
!
      REAL(KIND = 8) :: XINCR
!
! -------
! ... Integer variables
! -------
!
      INTEGER(KIND = 2) :: nmat
      INTEGER(KIND = 4) :: n, nloc, nx
!
      INTEGER(KIND = 1) :: State_X
!
! -------
! ... Logical variables
! -------
!
      LOGICAL :: First_call = .TRUE., array_size_uncertain = .TRUE.
!
! -------
! ... Saving variables
! -------
!
      SAVE First_call
!
!
!  =>=>=>=>=>=>=>=>   Main body of Equation_of_State
!
!
      IF(First_call) THEN
         First_call = .FALSE.
         WRITE(UNIT = VERS_Unit, FMT = 6000)
      END IF
!
      IF (Option_Print_EOSInfo >= 7) WRITE(*,6001) Tot_NumTimeSteps, Tot_NumNRIterations
!
!***********************************************************************
!*                                                                     *
!*                      SET UP INITIAL CONDITIONS                      *
!*                                                                     *
!***********************************************************************
!
      IF (NumTimeSteps == 0) CALL EOS_Initial_Assignment
!
      IF (Option_Print_EOSInfo >= 7) WRITE(*,6101)
!
! >>>>>>>>>>>>>>>>>>
! ...
! ... Determining the number of elements to include in the computations
! ...
! >>>>>>>>>>>>>>>>>>
!
! Check that properties are assigned to all gridblocks (including boundaries)
      IF_ArraySize: IF(array_size_uncertain) THEN
!
         IF_Check: IF(NumTimeSteps > 0) THEN
!
            ElemArraySize = NumElem
!
            array_size_uncertain = .FALSE.
!
         ELSE
!
            ElemArraySize = NumElemTot
!
         END IF IF_Check
!
      END IF IF_ArraySize
!
!
!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>                                                                     >
!>      Main loop determining conditions and assigning properties      >
!>                                                                     >
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
      DO_NumEle: DO n=1,ElemArraySize
!
!***********************************************************************
!*                                                                     *
!*        Set updated primary variables after last NR-iteration        *
!*                                                                     *
!***********************************************************************
!
         nmat = elem(n)%MatNum
         nloc = (n-1)*NumComPlus1
!
! ----------
! ...... Regular variable update in active gridblocks
! ----------
!
         DO_UpdatePV: DO nx=1,NumComPlus1
!
            XINCR = 0.0d0
            IF(NumNRIterations /= 0) XINCR = DX(nloc+nx)
            XX(nx) = X(nloc+nx) + XINCR
!
         END DO DO_UpdatePV
!
! ----------
! ...... Dealing with the possibility of non-physical or "NaN" primary variables
! ----------
!
         IF ( ( ABS(XX(1)) > 5.0d8 ) .OR. ( XX(1) < 0.0d0 ) .OR. ( XX(1) /= XX(1) ) ) unrealistic_conditions = .TRUE.
!
         IF (Option_Print_EOSInfo >= 7) WRITE(*,6102) elem(n)%name, XX(1:NumComPlus1)
!
         IF_NotGood: IF(unrealistic_conditions) THEN
            CALL Warning(elem(n)%name, XX, NumComPlus1, NumTimeSteps)
            RETURN
         END IF IF_NotGood
!
!
! >>>>>>>>>>>>>>>>>>
! ......
! ...... Choice of phase conditions after time step reduction
! ......
! >>>>>>>>>>>>>>>>>>
!
!
         IF_Ihaf_NE_0: IF(NumTimeStepCuts /= 0 .AND. Tot_NumNRIterations == 0) THEN
            State_X = ElemState%index(n,previous)
         ELSE
            State_X = ElemState%index(n,current)
         END IF IF_Ihaf_NE_0
!
!
! >>>>>>>>>>>>>>>>>>
! ......
! ...... Choice of phase conditions in a "normal" case
! ......
! >>>>>>>>>>>>>>>>>>
!
         Case_State: SELECT CASE(State_X)
!
! ----------
! ...... Single phase: Aqueous
! ----------
!
         CASE(1_1)
!
            CALL Single_Phase_Gas(n,nloc,nmat)
!
            IF(unrealistic_conditions) THEN
              RETURN
            ELSE
               CYCLE DO_NumEle
            END IF
!
         END SELECT Case_State
!
!
      END DO DO_NumEle
!
!
!***********************************************************************
!*                                                                     *
!*    PRINT SECONDARY PARAMETERS (When Option_Print_EOSInfo >= 8)      *
!*                                                                     *
!***********************************************************************
!
!
      IF(Option_Print_EOSInfo >= 8) CALL Print_Secondary_Variables
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 6000 FORMAT(/,'Equation_of_State',T50,'[v1.0,   01 June      2008]',/,  &
     &         ':::::   Selection of the appropriate state in the Gas system')
!
 6001 FORMAT(/,' !!!!!!!!!!!   EOS (Gas)  :::::::::::::  [Timestep, NR-iterations] = [',I6,',',I7,']'/)
!
 6101 FORMAT(/,' PRIMARY VARIABLES',/)
 6102 FORMAT(' AT ELEMENT "',A8,'" ::: ',(4(2X,1pE20.13)))
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of Equation_of_State
!
!
      RETURN
!
      END SUBROUTINE Equation_of_State
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE EOS_Initial_Assignment
!
! ... Modules to be used
!
         USE Universal_Parameters
         USE Basic_Parameters
         USE EOS_Default_Parameters
!
         USE General_External_File_Units
         USE General_Control_Parameters
!
         USE Gas_Parameters
!
         USE Grid_Geometry, ONLY: elem
         USE Element_Attributes
!
         USE Geologic_Media_Parameters
!
         USE Gas_Thermophysical_Properties
!
         USE Solution_Matrix_Arrays, ONLY: X, DX, DELX
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*            Assigns initial thermophysical properties                *
!*                     in a Gas system                                 *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
      IMPLICIT NONE
!
! -------
! ... Double precision variables
! -------
!
      REAL(KIND = 8) :: PX, TX
!
! -------
! ... Integer variables
! -------
!
      INTEGER :: n, nloc
!
! -------
! ... Logical variables
! -------
!
      LOGICAL :: First_call = .TRUE.
!
! -------
! ... Saving variables
! -------
!
      SAVE First_call
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of EOS_Initial_Assignment
!
!
      IF(First_call) THEN
         WRITE(VERS_Unit,6000)
         First_call = .FALSE.
      END IF
!
      IF (Option_Print_EOSInfo >= 7) WRITE(*,6001) Tot_NumTimeSteps, Tot_NumNRIterations
!
!***********************************************************************
!*                                                                     *
!*                      SET UP INITIAL CONDITIONS                      *
!*                                                                     *
!***********************************************************************
!
      DO_InConSt: DO n = 1,NumElemTot
!
         nloc = (n-1) * NumComPlus1
!
!***********************************************************************
!*                                                                     *
!*   Check if the element is under single-phase conditions: Gas        *
!*                                                                     *
!***********************************************************************
!
         Case_State: SELECT CASE(ElemState%index(n,current))
!
         CASE(1_1)
!
            PX = X(nloc+1)
            TX = X(nloc+2)
!
! -------------
! ......... Erroneous 1-phase conditions are specified : Temperature < freezing
! -------------
!
            IF (TX < 0.0d0) THEN
               WRITE(*,6501) elem(n)%name, TX         ! Print a clarifying message about the error
               STOP                                   ! Stop the simulation
            END IF
!
! -------------
! ......... Erroneous 1-phase conditions are specified : Pressure too high or too low
! -------------
!
            IF ( PX < 0.0d0 .OR. PX > 5.0d8 ) THEN
               WRITE(*,6502) elem(n)%name, PX         ! Print a clarifying message about the error
               STOP                                   ! Stop the simulation
            END IF
!
! -------------
! ......... Initializing the pressure, temperature, saturations
! -------------
!
            ElemState%pres(n,initial:previous) = PX
            ElemState%temp(n,initial:previous) = TX
!
            ElemProp(n,0)%satur(GasPhase) = 1.0d0
!
! ----------
! ..... Unrealistic/incorrect initial conditions
! ----------
!
        CASE DEFAULT
!
           CALL Warning(elem(n)%name, XX, NumComPlus1, NumTimeSteps)
!
        END SELECT Case_State
!
!
      END DO DO_InConSt
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 6000 FORMAT(/,'EOS_Initial_Assignment',T50,'[v1.0,   01 June      2008]',/,  &
     &         ':::::   Initial assignment of thermophysical properties in a Gas system')
!
 6001 FORMAT(/,' !!!!!!!!!!!   EOS (Gas)  :::::::::::::  [Timestep, NR-iterations] = [',I6,',',I7,']'/)
!
!
!
 6501 FORMAT(//,20('ERROR-'),//,  &
     &         T20,'R U N   A B O R T E D',/,  &
     &         T20,'The temperature specified in element "',a8,'" is erroneous !!!',/,  &
     &         T20,'The temperature T =',1pe11.4,' oC is outside the allowed range ( T > 0 oC)',//,  &
     &         T20,'CORRECT AND TRY AGAIN',//,  &
     &          20('ERROR-'))
!
 6502 FORMAT(//,20('ERROR-'),//,  &
     &         T20,'R U N   A B O R T E D',/,  &
     &         T20,'The pressure specified in element "',a8,'" is erroneous !!!',/,  &
     &         T20,'The pressure P =',1pe11.4,' Pa is outside the allowed range ( 0 < P < 1.0E8 Pa)',//,  &
     &         T20,'CORRECT AND TRY AGAIN',//,  &
     &          20('ERROR-'))
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End of EOS_Initial_Assignment
!
!
      RETURN
!
      END SUBROUTINE EOS_Initial_Assignment
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE Single_Phase_Gas(N,nloc,nmat)
!
! ...... Modules to be used
!
         USE Universal_Parameters
         USE Basic_Parameters
         USE EOS_Default_Parameters
         USE Gas_Parameters
!
         USE General_External_File_Units
         USE General_Control_Parameters
!
         USE Grid_Geometry, ONLY : elem
         USE Element_Attributes
!
         USE Geologic_Media_Parameters
!
         USE Gas_Thermophysical_Properties
!
         USE Solution_Matrix_Arrays, ONLY: X, DX, DELX
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*       ROUTINE FOR COMPUTING AND ASSIGNING ALL THE SECONDARY         *
!*         VARIABLES UNDER SINGLE-PHASE (Gas) CONDITIONS               *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
      IMPLICIT NONE
!
! -------
! ... Double precision variables
! -------
!
      REAL(KIND = 8) :: PX, TX, T_k
!
      REAL(KIND = 8) :: gas_phase_density, gas_phase_viscosity
      REAL(KIND = 8) :: gas_phase_IntEnergy, gas_phase_enthalpy
      REAL(KIND = 8) :: gas_compressibility, TR, RK, alpha, AxAl, AUP, BUP
      REAL(KIND = 8) :: gas_phase_IdealEnthalpy, gas_phase_IdealIntEnergy
      REAL(KIND = 8) :: gas_phase_DepEnthalpy, gas_phase_DepIntEnergy
!
! -------
! ... Integer variables
! -------
!
      INTEGER(KIND = 2), INTENT(IN) :: nmat
      INTEGER(KIND = 4), INTENT(IN) :: N, nloc
!
      INTEGER :: k
!
! -------
! ... Logical variables
! -------
!
      LOGICAL :: First_call = .TRUE.
!
! -------
! ... Saving variables
! -------
!
      SAVE First_call
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of Single_Phase_gas
!
!
      IF(First_call) THEN
         First_call = .FALSE.
         WRITE(UNIT = VERS_Unit, FMT = 6000)
      END IF
!
!
!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>                                                                     >
!>      Main loop determining conditions and assigning properties      >
!>                                                                     >
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
      DO_NumEqs: DO k = 1,NumEquPlus1
!
! ----------
! ...... Initial assignments of primary variable values
! ...... Primary variable: P
! ----------
!
         PX = XX(1)
         TX = XX(2)
!
! ...... Ensure a realistic pressure
!
         IF(PX > 5.0d8 .OR. PX < 0.0d0) THEN
            unrealistic_conditions = .TRUE.
            RETURN
         END IF
!
! ...... Ensure a realistic temperature
!
         IF(TX > 1.0d5 .OR. TX < -25.0d0) THEN
            unrealistic_conditions = .TRUE.
            RETURN
         END IF
!
!
! >>>>>>>>>>>>>>>>>>
! ......
! ...... NOTE: In the first iteration      (k = 1): Unincremented variables
! ......       In the remaining iterations (k > 1): Incremented variables
! ......
! >>>>>>>>>>>>>>>>>>
!
!
         CASE_Incr: SELECT CASE(k)
!
! ...... Incrementing pressure
!
         CASE(2)
!
            DELX(nloc+1) = derivative_increment*XX(1)
            PX           = XX(1) + DELX(nloc+1)
!
         CASE(3)
!
            DELX(nloc+2) = derivative_increment*XX(2)
            TX           = XX(2) + DELX(nloc+2)
!
         END SELECT CASE_Incr
!
!***********************************************************************
!*                                                                     *
!*             P R O P E R T I E S  +  P A R A M E T E R S             *
!*                                                                     *
!***********************************************************************
!
  500    T_k = TX + 273.15d0
!
! -------------
! ...... Density of the gas phase
! -------------
!
         IF( Compute_Gas_Density(PX, T_k, gas_compressibility, gas_phase_density, TR, RK, alpha,  &
        &                                   AxAl, AUP, BUP) .EQV. .FALSE.) THEN
            unrealistic_conditions = .TRUE.
            RETURN
         END IF
!
! ----------
! ...... Viscosity of the gas phase
! ----------
!
         IF( Compute_Gas_Viscosity(TX, gas_phase_density, gas_phase_viscosity) .EQV. .FALSE.) THEN
            unrealistic_conditions = .TRUE.
            RETURN
         END IF
!
! ----------
! ...... Enthalpy of the gas phase
! ----------
!
         IF(nonisothermal_conditions) THEN
            IF(     (Compute_Ideal_Gas_Enthalpy(273.15d0, T_k, gas_phase_IdealEnthalpy, gas_phase_IdealIntEnergy) .EQV. .FALSE.)                      &
          &    .OR. (Compute_Departure_Enthalpy( T_k, gas_compressibility, gas_phase_DepEnthalpy, gas_phase_DepIntEnergy, TR,   &
          &          RK, alpha, AxAL, AUP, BUP) .EQV. .FALSE.)) THEN
!
               unrealistic_conditions = .TRUE.
               RETURN
         END IF
!
         gas_phase_enthalpy   = 1.0d-3*(gas_phase_IdealEnthalpy  - gas_phase_DepEnthalpy  + gas%reference_departure_enthalpy)/gas%MolWeight
         gas_phase_IntEnergy  = 1.0d-3*(gas_phase_IdealIntEnergy - gas_phase_DepIntEnergy + gas%reference_departure_IntEnergy)/gas%MolWeight
         END IF
!
! ----------
! ...... Assign value to state index
! ----------
!
         IF(k == 1) ElemState%index(n,current) = 1_1
!
! ----------
! ...... Assign values of secondary parameters
! ----------
!
         ElemProp(n,k-1)%satur(GasPhase)        = 1.0d0
         ElemProp(n,k-1)%MassFrac(gas_2,GasPhase) = 1.0d0
!
         ElemProp(n,k-1)%viscosity(GasPhase) = gas_phase_viscosity
         ElemProp(n,k-1)%density(GasPhase)   = gas_phase_density
!
         ElemProp(n,k-1)%RelPerm(GasPhase) = 1.0d0

         IF(nonisothermal_conditions) THEN
            ElemProp(n,k-1)%enthalpy(GasPhase)  = gas_phase_enthalpy
            ElemProp(n,k-1)%IntEnergy(GasPhase) = gas_phase_IntEnergy
         ELSE
            ElemProp(n,k-1)%enthalpy(GasPhase)  = 0.0d0
            ElemProp(n,k-1)%IntEnergy(GasPhase) = 0.0d0
         END IF
!
! ----------
! ...... Store the incremented pressure and temperature
! ----------
!
         ElemState%pres(n,k-1) = PX
         ElemState%temp(n,k-1) = TX
!
!
!
      END DO DO_NumEqs
!
!
!
      RETURN     ! Task completed - Exit the routine
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 6000 FORMAT(/,'Single_Phase_Gas',T50,'[v1.0,   01 June      2008]',/,  &
     &         ':::::   Computes the secondary variables under 1-phase (gas) conditions in a 100% gas system')
!
 6001 FORMAT(' X_aA in <Phase1_Aqueous> : Element "',A8,'" ==> Px =',1pE13.6,   &
     &       ' X_aA =',1pE13.6,' Tx =',1pE13.6,'oC  PSat_H2O =',1pE13.6)
 6004 FORMAT(' !!!!!! <Phase1_Aqueous>: Gas phase evolves at element "',A8,'" ==> Px  = ',1pE15.8,'  P_gas = ',1pE15.8)
 6005 FORMAT(' !!!!!! <Phase1_Aqueous>: Ice evolves at element "',A8,'" ==> TX = ',1pE15.8,' oC   T_fusion = ',1pE15.8,' oC')
!
 6101 FORMAT(' !!!!!!!!!!   Routine <Phase1_Aqueous>, Element "',A8,&
     &       '": ERRONEOUS DATA INITIALIZATION  ==>  STOP EXECUTION <<<<<<<<<<')
!
 6102 FORMAT(' !!!!!!!!!!   Routine <Phase1_Aqueous>: ERRONEOUS PHASE CHANGE   ==>  STOP EXECUTION <<<<<<<<<<')
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End of Single_Phase_Gas
!
!
      RETURN
!
      END SUBROUTINE Single_Phase_Gas
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE JACOBIAN_SetUp
!
! ... Modules to be used
!
         USE EOS_Routine_Selector
!
         USE Universal_Parameters
         USE Basic_Parameters
         USE EOS_Default_Parameters
!
         USE General_External_File_Units
         USE General_Control_Parameters
!
         USE Gas_Parameters
         USE Gas_Thermophysical_Properties
!
         USE Grid_Geometry
         USE Element_Attributes
         USE Connection_Attributes
!
         USE Geologic_Media_Parameters
         USE Geological_Media_Properties
!
         USE Solution_Matrix_Arrays
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*             ROUTINE FOR SETTING UP THE JACOBIAN MATRIX              *
!*                                                                     *
!*                          GAS System                                 *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
      IMPLICIT NONE
!
! -------
! ... Double precision arrays
! -------
!
      REAL(KIND = 8), DIMENSION(NumPerturb,NumPerturb) :: Flow
!
      REAL(KIND = 8), DIMENSION(NumPhases,NumEquPlus1) :: phase_mass
!
      REAL(KIND = 8), DIMENSION(NumEqu,NumEquPlus1)    :: accum
!
      REAL(KIND = 8), DIMENSION(NumPhases)             :: phase_therm_conductivities1, phase_therm_conductivities2
!
      REAL(KIND = 8), ALLOCATABLE, DIMENSION(:) :: DistInv, AdjGrav
!
! -------
! ... Double precision variables
! -------
!
      REAL(KIND = 8) :: D1,D2,FlowArea,FAC1,FAC2,dPdL
      REAL(KIND = 8) :: perm1,perm2,phi, rock_heat_capacity
!
      REAL(KIND = 8) :: REL1,S1,VIS1,RHO1,RHO10,W1,Wup_1
      REAL(KIND = 8) :: REL2,S2,VIS2,RHO2,RHO20,W2,Wup_2
!
      REAL(KIND = 8) :: AvgDensity, P_gradient, XM1, XM2, mobility
      REAL(KIND = 8) :: PhaseFlow, CompInPhase, enthalpy, thermal_conductivity, gas_therm_conductivity
!
      REAL(KIND = 8) :: LogRatio1, LogRatio2
! -------
! ... Integer variables
! -------
!
      INTEGER(KIND = 2) :: nmat1, nmat2
!
      INTEGER :: n,nloc,m,np,k, mn
      INTEGER :: n1,n2,N1LOC,N2LOC,N1LOCP,N2LOCP,ISO
      INTEGER :: col, row
!
      INTEGER :: NN1, KK1, N1KL, N2KL, m1, mm1, m2, mm2, nm1, nm2
!
! -------
! ... Logical variables
! -------
!
      LOGICAL :: First_call = .TRUE.
!
! -------
! ... Saving variables
! -------
!
      SAVE First_call, DistInv, AdjGrav
!
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of JACOBIAN_SetUp
!
!
! ... If first time through subroutine, ALLOCATE array memory, etc.
      IF(First_call) THEN
!
         First_call = .FALSE.
         WRITE(VERS_Unit,6000)
!
! ...... Allocate local working arrays
!
         ALLOCATE( wt1(NumConx), wt2(NumConx), DistInv(NumConx), AdjGrav(NumConx) )

! ...... Define elements of local arrays - ARRAY operations
!
         DistInv = 1.0d0/(conx(1:NumConx)%d1 + conx(1:NumConx)%d2)
         DO n=1,NumConx
            wt1(n)     = conx(n)%d2*DistInv(n)
            wt2(n)     = 1.0d0 - wt1(n)
         END DO
         AdjGrav = conx(1:NumConx)%beta*gravity
!
         IF(coordinate_system .EQ. 'CYL') THEN
!
            DO n=1,NumConx
!
               n1 = conx(n)%n1
               n2 = conx(n)%n2
!
               IF(n1 == 0 .OR. n2 == 0)            CYCLE
               IF(n1 > NumElem .AND. n2 > NumElem) CYCLE
!
               IF(ABS(elem(n1)%coord(1) - elem(n2)%coord(1)) >= 1.0d-4 ) THEN
!
                  wt1(n) = LOG( elem(n2)%coord(1) /(elem(n1)%coord(1) + conx(n)%d1)) / &
                &          LOG( elem(n1)%coord(1) / (elem(n1)%coord(1) + noise) )
                  wt2(n) = 1.0d0 - wt1(n)
!
               END IF
!
            END DO

         END IF
      END IF
!
! -------
! ... Some print-out for large Option_Print_JacobianInfo values
! -------
!
      IF_PrOut: IF(Option_Print_JacobianInfo >= 1) THEN
         WRITE(*,6001) Tot_NumTimeSteps, NumNRIterations
         WRITE(*,6002)
      END IF IF_PrOut
!
!
!***********************************************************************
!*                                                                     *
!*                           INITIALIZATIONS                           *
!*                                                                     *
!***********************************************************************
!
!
      Non0 = 0
!
! ... Whole array operations
!
      phase_mass = 0.0d0
      accum      = 0.0d0
!
!
!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>                                                                     >
!>                          ELEMENT LOOPS                              >
!>  Compute accumulation-related properties, parameters and variables  >
!>                                                                     >
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
      DO_NumEle: DO n = 1,NumElem
!
         nloc = Loc(n)
!
         IF(nonisothermal_conditions) THEN
            mn   = elem(n)%MatNum
            rock_heat_capacity = (1.0d0-ElemMedia(n,original)%porosity)*media(mn)%DensG*media(mn)%SpcHt
         END IF
!
!
! >>>>>>>>>>>>
! >>>>>>>>>>>>
! >>>>>>>>>>>>
! ......
! ...... PERTURBATION LOOP: Computation of the accumulation terms at state points (m = 1) and using incremented variables (m > 1)
! ......
! >>>>>>>>>>>>
! >>>>>>>>>>>>
! >>>>>>>>>>>>
!
!
! +1 for state equation (initial condition), NumEqu for # increments, NumEqu = 2*NumEqu when non-isothermal
         DO_NumEq: DO m = 1,NumEquPlus1
!
!
! ......... Initializations
!
            mm1 = m-1
!
            phi = ElemMedia(n,mm1)%porosity
!
!
! >>>>>>>>>>>>>>>
! .........
! ......... PHASE LOOP: Computation of the mass of each phase
! .........
! >>>>>>>>>>>>>>>
!
!
            DO_NumPhase: DO np=1,NumPhases
!
! ............ Computation of < Phi * S * Rho >
!
! Calculate phase mass
               phase_mass(np,m) = phi * ElemProp(n,mm1)%satur(np) * ElemProp(n,mm1)%density(np)
!
            END DO DO_NumPhase
!
!
! >>>>>>>>>>>>>>>
! .........
! ......... COMPONENT LOOP: Computation of the mass of each component in the various phases
! .........
! >>>>>>>>>>>>>>>
!
!
            DO_NumComp: DO k=1,NumCom
!
! ......... Computation of < Phi * S * Rho > * < X >
!
! Compute mass of component, m is perturbation
               accum(k,m) = SUM( phase_mass(1:NumPhases,m) * ElemProp(n,mm1)%MassFrac(k,1:NumPhases) )
!
            END DO DO_NumComp
!
! >>>>>>>>>>>>>>>
! .........
! ......... HEAT Accumulation
! .........
! >>>>>>>>>>>>>>>
!
            IF(nonisothermal_conditions) THEN
               accum(heat,m) =  rock_heat_capacity * ElemState%temp(n,mm1) &
             &                + SUM( phase_mass(1:NumPhases,m) * ElemProp(n,mm1)%IntEnergy(1:NumPhases) )
            END IF

!***********************************************************************
!*                                                                     *
!*    Assign accumulation terms at the beginning of each time step     *
!*                                                                     *
!***********************************************************************
!
            IF_BeginDt: IF(NumNRIterations == 1 .AND. m == 1) THEN
!
! ........... <old_accum(Loc(n)+K)> holds the quantity of component k in element N at the end of the last completed time step
!
               old_accum(nloc+1 : nloc+NumEqu) = accum(1:NumEqu,m)              ! Whole array operations
!
            END IF IF_BeginDt
!
            IF(m == 1 .AND. Option_Print_JacobianInfo >= 3) THEN
               WRITE(*,6005) elem(n)%name, (accum(k,m), k=1,NumEqu)
            END IF
!
! <<<<<<
! <<<<<< End of the EQUATIONS LOOP
! <<<<<<
!
         END DO DO_NumEq
!
!
!***********************************************************************
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*             ASSIGNING ALL SINGLE-ELEMENT-RELATED TERMS              *
!*                       IN THE JACOBIAN MATRIX                        *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!***********************************************************************
!
!
! ...... <row> is the row index in matrix block (N,N)
!
         DO_NumEq1: DO row = 1,NumEqu
!
! ......... Initialize residuals
!
! Accumulation residual
            R(nloc+row) = accum(row,1) - old_accum(nloc + row)
!
! -------------
! ......... <col> is the column index in matrix block (N,N)
! ......... Compute matrix element arising from dependence of component "row" upon primary variable "col"
! -------------
!
            DO_NumEq2: DO col = 1,NumEqu
!
               RowNum(Non0+1) = nloc + row
               ColNum(Non0+1) = nloc + col
!
               IF(ElemMedia(n,0)%porosity == 0.0d0 .AND. row == col .AND. row /= NumComPlus1) THEN
                  CO(Non0+1) = 1.0d0
               ELSE
                  CO(Non0+1) = -(Accum(row,col+1) - Accum(row,1))/DELX(Locp(n)+col)
               END IF
!
               Non0 = Non0+1
!
            END DO DO_NumEq2
!
! <<<<<<
! <<<<<< End of the single-element assignments
! <<<<<<
!
         END DO DO_NumEq1
!
! <<<
! <<< End of the MAIN ELEMENT LOOP
! <<<
!
      END DO DO_NumEle
!
!
!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>                                                                     >
!>                         Sinks & Sources                             >
!>                                                                     >
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
      IF(NumSS /= 0) CALL SourceSink_Equation_Terms
!
!
!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>                                                                     >
!>                      MAIN CONNECTION LOOP                           >
!>                                                                     >
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
     FORALL (n=1:NumConx) ConxFlow(n)%CompInPhase(1:NumCom,1:NumMobPhases) = 0.0d0   ! ... Initialization - Whole array operation
!
!
! >>>>>>>>>>>>>>>>>>
! ...
! ... At each interface, we have one flux for each component, which is a sum of contributions from all phases.
! ... For computing derivatives, we need also the fluxes obtained by incrementing
! ... any of the <NumEqu> variables for the "1st" element, and any of the <NumEqu> variables for the "2nd" element
! ...
! >>>>>>>>>>>>>>>>>>
!
!
      DO_NumConX: DO n=1,NumConx
!
!
!
         n1 = conx(n)%n1
         n2 = conx(n)%n2
!
         IF(n1 == 0 .OR. n2 == 0)            CYCLE DO_NumConX
         IF(n1 > NumElem .AND. n2 > NumElem) CYCLE DO_NumConX
!
! ...... Identify location of matrix blocks in element arrays
!
         N1LOC  = Loc(n1)
         N2LOC  = Loc(n2)
!
         N1LOCP = Locp(n1)
         N2LOCP = Locp(n2)
!
! ...... Rock type/material and associated properties
!
         nmat1 = elem(n1)%MatNum
         nmat2 = elem(n2)%MatNum
!
! ----------
! ...... Basic connection parameters & quantities
! ----------
!
         D1 = conx(n)%d1
         D2 = conx(n)%d2
!
         ISO = conx(n)%ki
         FlowArea = conx(n)%area
!
! ...... Assign factors for flux terms
!
         IF(n1 <= NumElem) THEN
            FAC1 = TimeIncrement/elem(n1)%vol
         ELSE
            FAC1 = 0.0d0
         END IF
!
         IF(n2 <= NumElem) THEN
            FAC2 = TimeIncrement/elem(n2)%vol
         ELSE
            FAC2 = 0.0d0
         END IF
!
         IF(Option_Print_JacobianInfo >= 4) WRITE(*,6008) n, elem(n1)%name, elem(n2)%name, FAC1, FAC2
!
! ----------
! ...... Flow initialization - Whole array operation
! ----------
!
         Flow = 0.0d0
!
!
! >>>>>>>>>>>>>>>>>>
! ......
! ...... PERTURBATION LOOP
! ......
! ...... M=1 corresponds to the parameter values at the state point = latest updated variables
! ...... M=2, .......... ,  NumEquPlus1 are the NumEqu perturbations corresponding to incrementing the variables for element N1
! ...... M=NumEqu+2, ... ,2*NumEquPlus1 are the NEQ perturbations corresponding to incrementing the variables for element N2
! ......
! ......
! >>>>>>>>>>>>>>>>>>
!
         DO_NumFlo: DO m=1,NumPerturb
!
! -------------
! ......... Determine indices of secondary variables
! ......... m1 denotes fluxes incremented at N1, m2 denotes fluxes incremented at N2
! -------------
!
            IF_Perturations: IF(m >= 2 .AND. m <= NumEquPlus1) THEN
                         m1 = m
                         m2 = 1
                      ELSE IF(m == 1) THEN
                         m1 = 1
                         m2 = 1
                      ELSE
                         m1 = 1
                         m2 = m - NumEqu
                      END IF IF_Perturations
!
                      mm1 = m1 - 1                                      !Perturation for grid block 1
                      mm2 = m2 - 1                                      !Perturation for grid block 2
!
! >>>>>>>>>>>>>
! .........
! ......... Compute pressure gradients
! .........
! >>>>>>>>>>>>>
!
            dPdL = ( ElemState%pres(n2,mm2) - ElemState%pres(n1,mm1) )*DistInv(n)
!
! >>>>>>>>>>>>>
! .........
! ......... Compute conductive heat transport
! .........
! >>>>>>>>>>>>>
!
        IF(nonisothermal_conditions) THEN
!
! ... Calculate gas conductivity if the linear moddel is to be used
!
           IF_Linear:  IF(Option_ThermalConductivity > 1) THEN
!
              IF((Compute_Gas_Therm_Conductivity( (ElemState%temp(n1,mm1) + 273.15d0),                &
       &                                        phase_therm_conductivities1(1) ) .EQV. .FALSE.) .OR.   &
       &         (Compute_Gas_Therm_Conductivity( (ElemState%temp(n2,mm2) + 273.15d0),                &
       &                                        phase_therm_conductivities2(1) ) .EQV. .FALSE.))       &
              THEN
                 unrealistic_conditions = .TRUE.
                 RETURN
              END IF
!
              thermal_conductivity = Composite_Thermal_Conductivity ( n1 = n1, n2 = n2, m1 = mm1, m2 = mm2,                 &
       &                                                             nm1 = nmat1, nm2 = nmat2, wt1in = wt1(n), wt2in = wt2(n),  &
       &                                                             gas_thermal_cond1 = phase_therm_conductivities1,       &
       &                                                             gas_thermal_cond2 = phase_therm_conductivities2 )
           ELSE
              thermal_conductivity = Composite_Thermal_Conductivity ( n1 = n1, n2 = n2, m1 = mm1, m2 = mm2,                 &
       &                                                             nm1 = nmat1, nm2 = nmat2, wt1in = wt1(n), wt2in = wt2(n),  &
       &                                                             liquid_phases = (/ 1 /) )
!
           END IF IF_Linear
!
           flow(heat, m) = FlowArea * thermal_conductivity * (ElemState%temp(n2,mm2) - Elemstate%temp(n1,mm1) ) * DistInv(n)
!
        END IF
!
! >>>>>>>>>>>>>
! .........
! ......... Compute the interface intrinsic permeability (harmonic average)
! .........
! >>>>>>>>>>>>>
!
! ......... Permeability adjustment - Klinkenberg factors
!
            perm1 = ElemMedia(n1,mm1)%perm(iso)*(1.0d0 + media(nmat1)%Klink/ElemState%pres(n1,mm1))
!
            perm2 = ElemMedia(n2,mm2)%perm(iso)*(1.0d0 + media(nmat2)%Klink/ElemState%pres(n2,mm2))
!
! -------------
! ......... Print-outs for large Option_Print_JacobianInfo values
! -------------
!
            IF((Option_Print_JacobianInfo >= 6 .AND. m == 1) .OR. Option_Print_JacobianInfo >= 8) THEN
               WRITE(*, FMT = 6010) m, ElemState%pres(n1,mm1), ElemState%pres(n2,mm2), dPdL
            END IF
!
!
! >>>>>>>>>>>>>>>>>>
! .........
! ......... PHASE LOOP (Obtain Flow for Each Phase)
! .........
! >>>>>>>>>>>>>>>>>>
!
!
            DO_NumPhFlo: DO np=1,NumMobPhases
!
! ............ Initializations
!
               IF(m == 1) ConxFlow(n)%PoreVel(np) = 0.0d0
               PhaseFlow = 0.0d0
!
! ............ Relative Permeabilities
!
               REL1 = ElemProp(n1,mm1)%RelPerm(np)
               REL2 = ElemProp(n2,mm2)%RelPerm(np)
!
! ............ Check phase mobility & realistic perm. regime
!
               IF(REL1 == 0.0d0 .AND. REL2 == 0.0d0) GO TO 500
!
! ----------------
! ............ Retrieve properties and parameters of the two connection elements
! ----------------
!
               S1     = ElemProp(n1,mm1)%satur(np)
               VIS1   = ElemProp(n1,mm1)%viscosity(np)
               RHO1   = ElemProp(n1,mm1)%density(np)
               RHO10  = ElemProp(n1,0)%density(np)
!
               S2     = ElemProp(n2,mm2)%satur(np)
               VIS2   = ElemProp(n2,mm2)%viscosity(np)
               RHO2   = ElemProp(n2,mm2)%density(np)
               RHO20  = ElemProp(n2,0)%density(np)
!
! ............ Compute weighted interface density
!
               IF_W1Val: IF(RHO1 == 0.0d0) THEN
                            W1 = 0.0d0
                         ELSE IF(RHO2 == 0.0d0) THEN
                            W1 = 1.0d0
                         ELSE
                            W1 = 5.0d-1
                         END IF IF_W1Val
!
               W2 = 1.0d0 - W1
!
               AvgDensity = W1*RHO1 + W2*RHO2
!
! ............ Effective pressure gradient
!
               P_gradient  = dPdL -  AvgDensity  * AdjGrav(n)
!
! ----------------
! ............ Determine upstream weighting factors
! ----------------
!
               IF(m == 1) THEN
!
                  IF_Upstream: IF(P_gradient > 0.0d0) THEN
                     Wup_1 = 1.0d0 - W_upstream
                  ELSE
                     Wup_1 = W_upstream
                  END IF IF_Upstream
!
                  Wup_2 = 1.0d0 - Wup_1
!
               END IF
!
! ----------------
! ............ Compute interface mobilities - harmonic mean
! ----------------
!
               XM1 = perm1*REL1/(VIS1 + 1.0d-20)
               XM2 = perm2*REL2/(VIS2 + 1.0d-20)
!
               mobility = XM1*XM2/( WT1(n)*XM1 + WT2(n)*XM2 + 1.0d-20 )
!
! ----------------
! ............ Compute pore velocities for fluid phases
! ----------------
!
               IF_PoreVel: IF(m == 1) THEN
!
                  ConxFlow(n)%PoreVel(np) = mobility * P_gradient/(  Wup_1*ElemMedia(n1,0)%porosity*S1           &
     &                                                             + Wup_2*ElemMedia(n2,0)%porosity*S2 + 1.0d-20 )
!
               END IF IF_PoreVel
!
! ----------------
! ............ Compute flow of phase NP in perturmation m
! ----------------
!
               PhaseFlow = mobility * AvgDensity * P_gradient * FlowArea
!
!
! >>>>>>>>>>>>>>>>>>>>>
! ............
! ............ COMPONENT LOOP (Obtain Flow for Each Component in Each Phase)
! ............
! >>>>>>>>>>>>>>>>>>>>>
!
!
               DO_NumCompn: DO k=1,NumCom
!
                  CompInPhase = Wup_1 * ElemProp(n1,mm1)%MassFrac(k,np) + Wup_2 * ElemProp(n2,mm2)%MassFrac(k,np)
!
                  IF(m == 1) ConxFlow(n)%CompInPhase(k,np) = CompInPhase * PhaseFlow
!
! ............... Compute fluid component flows
!
                  Flow(k,m) = Flow(k,m) + CompInPhase * PhaseFlow
!
!
! <<<<<<<<<<<<
! <<<<<<<<<<<< End of the COMPONENT LOOP
! <<<<<<<<<<<<
!
               END DO DO_NumCompn
!
! ----------------
! ............ Compute advective heat transport
! ----------------
!
               IF(nonisothermal_conditions) THEN
                  Enthalpy = Wup_1*ElemProp(n1,mm1)%enthalpy(np) + Wup_2*ElemProp(n2,mm1)%enthalpy(np)
                  flow(heat,m) = flow(heat,m) + PhaseFlow * Enthalpy
               END IF
!
!
  500          CONTINUE
!
! ----------------
! ............ Print-outs for Option_Print_JacobianInfo >= 7
! ----------------
!
               IF_PrintOut: IF((Option_Print_JacobianInfo == 7 .AND. m == 1) .OR.   &
     &                         (Option_Print_JacobianInfo == 8 .AND. m == 1) .OR.   &
     &                          Option_Print_JacobianInfo == 9                    ) &
     &         THEN
!
                  PRINT 6012, NP,perm1,perm2
                  PRINT 6014, S1,REL1,VIS1,RHO1
                  PRINT 6016, S2,REL2,VIS2,RHO2
                  PRINT 6018, mobility, AvgDensity, P_gradient, PhaseFlow
                  PRINT 6020, (Flow(K,M), k=1,heat)
!
               END IF IF_PrintOut
!
! ----------------
! ............ Store the fluxes in each phase (at the state point)
! ----------------
!
               IF(m == 1) ConxFlow(n)%rate(np) = PhaseFlow
!
! <<<<<<<<<
! <<<<<<<<< End of the PHASE LOOP
! <<<<<<<<<
!
            END DO DO_NumPhFlo
!
! <<<<<<
! <<<<<< End of the FLUX LOOP
! <<<<<<
!!
         END DO DO_NumFlo
!
!
!***********************************************************************
!*                                                                     *
!*                     ASSIGN ALL INTERFACE TERMS                      *
!*                                                                     *
!***********************************************************************
!
!

!
!
! >>>>>>>>>>>>>>>>>>>>>
! ......
! ...... EQUATION-1 LOOP (Filling-up the Jacobian matrix)
! ...... K is the row index within a block pertaining to element N1 or N2
! ......
! >>>>>>>>>>>>>>>>>>>>>
!
!
         DO_NumEqs1: DO row = 1,NumEqu
!
! -------------
! ......... Compute flux contributions to residuals
! -------------
!
!
            IF(n1 <= NumElem) R(N1LOC+row) = R(N1LOC+row) - FAC1*Flow(row,1)
!
            IF(n2 <= NumElem) R(N2LOC+row) = R(N2LOC+row) + FAC2*Flow(row,1)
!
!
! >>>>>>>>>>>>>>>>>>>>>
! .........
! ......... EQUATION-2 LOOP (Filling-up the Jacobian matrix)
! ......... L is the column index within a block (element N1 or N2)
! .........
! >>>>>>>>>>>>>>>>>>>>>
!
!
            DO_NumEqs2: DO col = 1,NumEqu
!
               IF_N1N2Act: IF(n1 <= NumElem .AND. n2 <= NumElem) THEN
!
! ............... Off-diagonal term in equation for element N1, arising from dependence
! ............... of the flux of component "row" upon variable "col" in element N2
!
                  RowNum(Non0+1) = N1LOC + row
                  ColNum(Non0+1) = N2LOC + col
!
                  CO(Non0+1) = FAC1*(Flow(row,col+1+NumEqu) - Flow(row,1))/DELX(N2LOCP+col)
                  Non0       = Non0+1
!
! ............... Off-diagonal term in equation for element N2, arising from dependence
! ............... of the flux of component "row" upon variable "col" in element N1
!
                  RowNum(Non0+1) = N2LOC + row
                  ColNum(Non0+1) = N1LOC + col
!
                  CO(Non0+1) = -FAC2*(Flow(row,col+1) - Flow(row,1))/DELX(N1LOCP+col)
                  Non0       = Non0+1
!
               END IF IF_N1N2Act
!
!
! ............ Diagonal term in equation for element N1, arising from dependence
! ............ of the flux of component "row" upon variable "col" in element N1
!
!
               IF_N1Act: IF(n1 <= NumElem) THEN
!
                  N1KL     = (n1-1)*NumEquSquare + (row-1)*NumEqu + col
!
                  CO(N1KL) = CO(N1KL) + FAC1*(Flow(row,col+1) - Flow(row,1))/DELX(N1LOCP+col)
!
               END IF IF_N1Act
!
!
! ............ Diagonal term in equation for element N2, arising from dependence
! ............ of the flux of component "row" upon variable "col" in element N2
!
!
               IF_N2Act: IF(n2 <= NumElem) THEN
!
                  N2KL     = (n2-1)*NumEquSquare + (row-1)*NumEqu + col
!
                  CO(N2KL) = CO(N2KL) - FAC2*(Flow(row,col+1+NumEqu) - Flow(row,1))/DELX(N2LOCP+col)
!
               END IF IF_N2Act
!
! <<<<<<<<<
! <<<<<<<<< End of the EQUATION-2 LOOP
! <<<<<<<<<
!
            END DO DO_NumEqs2
!
! <<<<<<
! <<<<<< End of the EQUATION-1 LOOP
! <<<<<<
!
            END DO DO_NumEqs1
!
! <<<
! <<< End of the CONNECTIONS LOOP
! <<<
!
      END DO DO_NumConX
!
!
!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>                                                                     >
!>                        CONVERGENCE TEST                             >
 !>                                                                     >
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
      IF(Option_Print_JacobianInfo >= 2) THEN
         PRINT 6022
         DO n = 1,NumElem
            PRINT 6024, elem(n)%name,(R((N-1)*NumEqu+K),K=1,NumEqu)
         END DO
         DO n = 1,Non0
            !PRINT 6824, n,RowNum(n),ColNum(n),CO(n)
 !6824 FORMAT(T5,'n=',i5,1x,'Row=',i5,1x,'Col=',i5,1x,'CO=',1pe15.8)
         END DO
      END IF
!
! -------
! ... Compute residuals - CAREFUL: Whole array/subarray operations
! -------
!
      WHERE(ABS(old_accum(1:NumUnknowns)) > abs_convergence_crit)
         rwork(1:NumUnknowns) = ABS(R(1:NumUnknowns)/old_accum(1:NumUnknowns))
      ELSEWHERE
         rwork(1:NumUnknowns) = ABS(R(1:NumUnknowns)/abs_convergence_crit)
      END WHERE
!
! -------
! ... Determine maximum residual and its location in the 1-D array - Whole array/subarray operations
! ... NOTE: No need to reinitialize <rwork> and <iwork> because they are reset in routine <Solve_Jacobian_Matrix_Equation>
! -------
!
      iwork(1:1) = MAXLOC(rwork(1:NumUnknowns))
!
      MaxResidual = rwork(iwork(1))
!
! -------
! ... Determine the location of maximum residual (Cell # and Equation #)
! -------
!
      NN1 = iwork(1)/NumEqu
      KK1 = MOD(iwork(1),NumEqu)
!
      IF(KK1 == 0) THEN
         MaxResElemNum = NN1
         MaxResEquNum  = NumEqu
      ELSE
         MaxResElemNum = NN1+1
         MaxResEquNum  = KK1
      END IF
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 6000 FORMAT(/,'JACOBIAN_SetUp',T50,'[v1.0,   01 June      2008]',6X,/,   &
     &         ':::::   Assemble all accumulation and flow terms; Populate the Jacobian matrix - Gas system')
!
 6001 FORMAT(/,' !!!!!!!!!!   ==> SUBROUTINE "JACOBIAN_SetUp":   [Timestep, NR-iteration] = [',I5,',',I3,']',/)
 6002 FORMAT(' <=><=><=><=><=>  ACCUMULATION TERMS:  Mass balances first, Energy balances last',/)
!
 6005 FORMAT('    At element "',A8,'":  ',5(1X,1pe17.10))
 6008 FORMAT(/,' ::::: Connection #',I5,'   Elements: (',A8,',',A8,')','  =>  (Dt/V1) = ',1pE12.5,'   (Dt/V2) = ',1pE12.5)
!
 6010 FORMAT(/,'   *** FLUX NO.',I3,15X,'      P1 = ',1pE12.5,'  P2 = ',1pE12.5,'    dPdL = ',1pE12.5)
 6011 FORMAT(21X,'       CONI = ',1pE12.5,'    T1 = ',1pE12.5,'    T2 = ',1pE12.5,' DTX = ',1pE12.5,'   FNK1 = ',1pE12.5)
!
 6012 FORMAT(9X,'NP =',I2,'    PER1 = ',1pE12.5,' PER2= ',1pE12.5)
 6014 FORMAT('       S1   = ',1pE12.5,'  REL1 = ',1pE12.5,'  VIS1 = ',1pE12.5,'  RHO1 = ',1pE12.5)
 6016 FORMAT('       S2   = ',1pE12.5,'  REL2 = ',1pE12.5,'  VIS2 = ',1pE12.5,'  RHO2 = ',1pE12.5)
 6018 FORMAT('       DMOBI= ',1pE12.5,'  RHOX = ',1pE12.5,'    DR = ',1pE12.5,'  PhaseFlow = ',1pE12.5)
 6020 FORMAT('     :: FLOW TERMS',8(1X,1pE12.5))
!
 6022 FORMAT(/,' => => => => =>  RESIDUALS:  Mass balances first, Energy balance last',/)
 6024 FORMAT('    AT ELEMENT "',A8,'":   ',4(1X,1pe17.10))
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End of JACOBIAN_SetUp
!
!
      RETURN
!
      END SUBROUTINE JACOBIAN_SetUp
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
