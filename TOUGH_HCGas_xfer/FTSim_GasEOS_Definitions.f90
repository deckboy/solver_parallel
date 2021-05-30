!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>                                                                     >
!>                                                                     >
!>   FTSim_GasEOS_Definitions.f95: Code unit including definitions:    >
!>      assignment of default parameters associated exclusively        >
!>                   with the gas equation of state                    >
!>                                                                     >
!>                                                                     >
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  Segment providing default parameter values describing the basic
!  attributes of the equation of state (i.e., number of components,
!  number of phases, etc.)
!
!
      MODULE EOS_Default_Parameters
!
         SAVE
!
! ----------
! ...... Integer parameters
! ----------
!
         INTEGER, PARAMETER :: Number_of_states = 1
!
         INTEGER, PARAMETER :: GasPhase = 1
         INTEGER, PARAMETER :: gas_2 = 1, heat = 2
!
! ----------
! ...... Character parameter arrays
! ----------
!
         CHARACTER(LEN=7), PARAMETER :: EOS_ID = "GAS_EOS"
!
         CHARACTER(LEN=3), PARAMETER, DIMENSION(Number_of_states) :: State_name = (/ "Gas" /)
!
! ----------
! ...... Character arrays
! ----------
!
         CHARACTER(LEN = 35), DIMENSION(Number_of_states) :: EOS_Variables
!
         CHARACTER(LEN = 400), DIMENSION(Number_of_states) :: Format_Generic_InitialCond
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
!
         SUBROUTINE EOS_Defaults( EOS_name, NumCom, NumEqu, binary_diffusion,                                         &
     &                                      Max_NumMassComp, Min_NumMassComp, Max_NumEquations, Min_NumEquations,     &
     &                                      Max_NumPhases, Max_NumMobPhases, NumStateParam, NumAuxParam,              &
     &                                      VERS_Unit, MemoryAlloc_Unit )
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*        Routine for default (maximum) numbers of components,         *
!*          equations, and phases of the EOS used by FTSim             *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
         IMPLICIT NONE
!
! -------------
! ...... Integer variables
! -------------
!
         INTEGER, INTENT(IN)  :: NumCom, NumEqu
         INTEGER, INTENT(IN)  :: VERS_Unit, MemoryAlloc_Unit
!
         INTEGER, INTENT(OUT) :: Max_NumMassComp, Min_NumMassComp, Max_NumEquations, Min_NumEquations
         INTEGER, INTENT(OUT) :: Max_NumPhases, Max_NumMobPhases, NumStateParam, NumAuxParam
!
! -------------
! ...... Character variables
! -------------
!
         CHARACTER(LEN = 15) :: EOS_name
!
! -------------
! ...... Logical variables
! -------------
!
         LOGICAL, INTENT(IN) :: binary_diffusion
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of EOS_Defaults
!
!
         WRITE(UNIT = VERS_Unit, FMT = 6000)
!
!
!***********************************************************************
!*                                                                     *
!*            DEFAULT NUMBERS FOR THE  EQUATION OF STATE               *
!*                                                                     *
!***********************************************************************
!
!
         SELECT CASE(EOS_name(1:7))
!
! ...... Parameters for the GAS Equation of State
!
         CASE('Gas_Eos', 'GAS_EOS', 'gas_eos')
!
            Max_NumMassComp  = 1   ! Maximum number of mass components
            Min_NumMassComp  = 1   ! Minimum number of mass components
            Max_NumEquations = 2   ! Maximum number of equations
            Min_NumEquations = 1   ! Minimum number of equations
            Max_NumPhases    = 1   ! Maximum number of phases
! CMF - Max_NumMobPhases = 1 causes problems with gfortran
            Max_NumMobPhases = 1   ! Maximum number of mobile phases
!
            NumAuxParam      = 0      ! The size of the auxilliary parameter array
!
            NumStateParam    = 2      ! The number of parameters stored at the state point
!
! ...... ERROR: Unavailable Equation of State
!
         CASE DEFAULT
!
            WRITE(*, 6001)               EOS_name(1:7)
            WRITE(MemoryAlloc_Unit,6001) EOS_name(1:7)
            STOP
!
         END SELECT
!
! ----------
! ...... Ensure that NumCom, NumEqu, NumPhases do not exceed the maximum numbers
! ----------
!
         IF( (NumCom > Max_NumMassComp) .OR. (NumEqu > Max_NumEquations) ) THEN
!
            WRITE(*,6010) NumCom, NumEqu, Max_NumMassComp, Max_NumEquations
            STOP
!
         END IF
!
! ----------
! ...... Ensure that NumCom, NumEqu, NumPhases are not smaller that the minimum numbers
! ----------
!
         IF( (NumCom < Min_NumMassComp) .OR. (NumEqu < Min_NumEquations) ) THEN
!
            WRITE(*,6011) NumCom, NumEqu, Min_NumMassComp, Min_NumEquations
            STOP
!
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
 6000 FORMAT(/,'EOS_Defaults',T50,'[v1.0,   7 April     2006]',/,   &
     &         ':::::   Default parameters of the equation of state for a pure gas system')
!
 6001 FORMAT(//,20('ERROR-'),//,T33,  &
     &             '       S I M U L A T I O N   A B O R T E D',/,  &
     &         T20,'The equation of state specified by EOS_Name = ',A7,' is unavailable'  &
     &       /,T32,'               CORRECT AND TRY AGAIN',  &
     &       //,20('ERROR-'))
!
 6010 FORMAT(//,20('ERROR-'),//,T33,                                                   &
     &             '       S I M U L A T I O N   A B O R T E D',//,                    &
     &         T6,'One or more of the parameters NumCom, NumEqu = (',i2,',',i2,')',/,  &
     &         T6,'exceeds the maximum values ',                                       &
     &            '<Max_NumMassComp>, <Max_NumEquations>, = (',i2,',',i2,')',/,        &
     &       /,T32,'               CORRECT AND TRY AGAIN',                             &
     &       //,20('ERROR-'))
!
 6011 FORMAT(//,20('ERROR-'),//,T33,                                                   &
     &             '       S I M U L A T I O N   A B O R T E D',//,                    &
     &         T6,'One or more of the parameters NumCom, NumEqu = (',i2,',',i2,')',/,  &
     &         T6,'exceeds the minimum values ',                                       &
     &            '<Min_NumMassComp>, <Min_NumEquations>, = (',i2,',',i2,')',/,        &
     &       /,T32,'               CORRECT AND TRY AGAIN',                             &
     &       //,20('ERROR-'))
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of EOS_Defaults
!
!
            RETURN
!
         END SUBROUTINE EOS_Defaults
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
         SUBROUTINE EOS_Primary_Variable_List(EOS_name,NumCom,NumEqu,N_EOS,EOS_variables,VERS_Unit,MemoryAlloc_Unit)
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*         Routine for identifying the primary variable list           *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
         IMPLICIT NONE
!
! -------------
! ...... Integer variables
! -------------
!
         INTEGER, INTENT(IN) :: NumCom,NumEqu,N_EOS
         INTEGER, INTENT(IN) :: VERS_Unit,MemoryAlloc_Unit
!
! -------------
! ...... Character variables
! -------------
!
         CHARACTER(LEN = 15) :: EOS_name
!
! -------------
! ...... Character variables
! -------------
!
         CHARACTER(LEN = 35), DIMENSION(N_EOS) :: EOS_variables
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of EOS_Primary_Variable_List
!
!
         WRITE(VERS_Unit,6000)
!
!
!***********************************************************************
!*                                                                     *
!*             DEFAULT NUMBERS FOR THE EQUATION OF STATE               *
!*                                                                     *
!***********************************************************************
!
!
         SELECT CASE(EOS_name(1:7))
!
! ...... Parameters for the Gas Equation of State
!
         CASE('Gas_Eos', 'GAS_EOS', 'gas_eos')
!
            EOS_variables(1) = ': P, T '
!
        Format_Generic_InitialCond(1) = '(">>>INITIAL_CONDITIONS",/,                            &
     &      "&Uniform_Initial_Conditions     state                   = ",a1,a3,a1,",",/,        &
     &      "                                primary_variable_names  = ",a1,a5,a1,",",/,        &
     &      "                                primary_variable_values = ",2(1pe20.13,1x,","),/,  &
     &      "                                /")'
!
! ...... ERROR: Unavailable Equation of State
!
         CASE DEFAULT
!
            WRITE(*, 6001)               EOS_name(1:7)
            WRITE(MemoryAlloc_Unit,6001) EOS_name(1:7)
            STOP
!
         END SELECT
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 6000 FORMAT(/,'EOS_Primary_Variable_List',T50,'[v1.0,  01 June 2008]',/,   &
     &         ':::::   Determines the list of primary variables for the GAS equation of state')
!
 6001 FORMAT(//,20('ERROR-'),//,T33,  &
     &             '       S I M U L A T I O N   A B O R T E D',/,  &
     &         T20,'The equation of state specified by EOS_Name = ',A7 ,' is unavailable'  &
     &       /,T32,'               CORRECT AND TRY AGAIN',  &
     &       //,20('ERROR-'))
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End of EOS_Primary_Variable_List
!
!
            RETURN
!
         END SUBROUTINE EOS_Primary_Variable_List
!
!
!
      END MODULE EOS_Default_Parameters
!
!
!
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
