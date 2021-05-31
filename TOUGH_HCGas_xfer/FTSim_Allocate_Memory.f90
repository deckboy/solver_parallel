!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>                                                                     >
!>                                                                     >
!>       T_AllocateMem.f95: Coode unit defining main variables,        >
!>    main objects and main arrays, and allocating memory to arrays    >
!>                                                                     >
!>                                                                     >
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  Segment responsible for the dynamic memory allocation (following
!  input describing the size of the problem) and dimensioning of most
!  arrays needed by the code, in addition to memory deallocation of
!  unnecessary arrays.
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      MODULE EOS_Parameters
!
         SAVE
!
! ----------
! ...... Integer parameters
! ----------
!
         INTEGER :: Max_NumMassComp   ! Maximum number of mass components
         INTEGER :: Min_NumMassComp   ! Minimum number of mass components
         INTEGER :: Max_NumEquations  ! Maximum number of equations
         INTEGER :: Min_NumEquations  ! Minimum number of equations
         INTEGER :: Max_NumPhases     ! Maximum number of phases
         INTEGER :: Max_NumMobPhases  ! Maximum number of mobile phases
!
!
      END MODULE EOS_Parameters
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      MODULE Basic_Parameters
!
         SAVE
!
! ----------
! ...... Double precision arrays
! ----------
!
         REAL(KIND = 8), DIMENSION(10) :: default_initial_cond
         REAL(KIND = 8), DIMENSION(10) :: XX
!
! ----------
! ...... Double precision parameters
! ----------
!
         REAL(KIND = 8), PARAMETER :: zero = 1.0d-6, noise = 1.0d-20
         REAL(KIND = 8), PARAMETER :: ZA   = 1.0d00
!
! ----------
! ...... Integer arrays
! ----------
!
         INTEGER, ALLOCATABLE, DIMENSION(:) :: Loc, Locp
!
! ----------
! ...... Integer parameters
! ----------
!
         INTEGER :: current = 0, initial = -2, previous = -1, original = -2
!
! ----------
! ...... Integer variables
! ----------
!
         INTEGER :: NumEqu, NumCom, NumPhases, NumMobPhases, NumStateParam, NumAuxParam
!
         INTEGER :: ElemNameLength
!
         INTEGER :: NumEquSquare, NumEquPlus1, NumComPlus1, NumPerturb
!
         INTEGER :: Max_NumElem, NumElemTot, NumElem, ElemArraySize
!
         INTEGER :: Max_NumConx, NumConx
!
         INTEGER :: Max_NumPrimaryVar, NumSecondaryVar
!
         INTEGER :: MatrixOrder, Max_NumNonZeros, Non0, NumAboveMain, NumBelowMain, NumUnknowns
!
         INTEGER :: real_work_array_size, intg_work_array_size
!
         INTEGER :: Max_NumSS, NumSS
!
         INTEGER :: Max_NumMedia, NumMedia
!
         INTEGER :: NumActiveDimensions, Active2ndDimension, longest(1)
!
! ----------
! ...... Character parameters
! ----------
!
         CHARACTER(LEN = 3), DIMENSION(6) :: Available_Equation = (/"Pol","Exp","Pow","Log","Sin","Cos"/)
!
! ----------
! ...... Character variables
! ----------
!
         CHARACTER(LEN =   3) :: coordinate_system, DataBlockEndSpecifier
         CHARACTER(LEN =  15) :: EOS_Name
         CHARACTER(LEN = 120) :: title
!
! ----------
! ...... Logical variables
! ----------
!
         LOGICAL :: nonisothermal_conditions
         LOGICAL :: binary_diffusion, file_VERS_exists
         LOGICAL :: variable_porosity, porosity_perm_dependence
         LOGICAL :: element_by_element_properties
!
!
      END MODULE Basic_Parameters
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      MODULE General_External_File_Units
!
         SAVE
!
! ----------
! ...... Integer parameters
! ----------
!
!
         INTEGER, PARAMETER :: MemoryAlloc_Unit = 10
         INTEGER, PARAMETER :: Active_Conx_Unit = 11
!
         INTEGER, PARAMETER :: MESH_Unit = 12
         INTEGER, PARAMETER :: MINC_Unit = 14
!
         INTEGER, PARAMETER :: INCON_Unit = 15
         INTEGER, PARAMETER :: GENER_Unit = 16
!
         INTEGER, PARAMETER :: VERS_Unit  = 17
         INTEGER, PARAMETER :: LINEQ_Unit = 18
!
         INTEGER, PARAMETER :: SAVE_Unit = 21
!
         INTEGER, PARAMETER :: Elem_TimeSeries_Unit = 22
         INTEGER, PARAMETER :: Conx_TimeSeries_Unit = 23
         INTEGER, PARAMETER :: SS_TimeSeries_Unit   = 24
!
         INTEGER, PARAMETER :: SinkSource_Unit = 25
!
         INTEGER, PARAMETER :: Plot_coord_Unit = 26
!
         INTEGER, PARAMETER :: Plot_Elem_Unit  = 27
         INTEGER, PARAMETER :: Plot_Conx_Unit  = 28
!
         INTEGER, PARAMETER :: Parameter_Update_Unit = 30
!
!
!
      END MODULE General_External_File_Units
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      MODULE Universal_Parameters
!
! ----------
! ...... Double precision parameters
! ----------
!
         REAL(KIND = 8), PARAMETER :: pi    = 3.14159265358979324D0
         REAL(KIND = 8), PARAMETER :: R_gas = 8.31456d3
!
         REAL(KIND = 8), PARAMETER :: StefanBoltzmanConst = 5.6687d-8
!
!
      END MODULE Universal_Parameters
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      MODULE Time_Series_Parameters
!
         SAVE
!
! ----------
! ...... Derived types
! ----------
!
         TYPE Observation_Elem
!
            INTEGER, DIMENSION(100)            :: num      ! Element # of observation elements
            CHARACTER(LEN = 5), DIMENSION(100) :: name     ! Names of observation elements
!
         END TYPE Observation_Elem
!
!
         TYPE Observation_Conx
!
            INTEGER, DIMENSION(100)             :: num     ! Connection # of observation elements
            CHARACTER(LEN = 16), DIMENSION(100) :: name    ! Names of observation connections
!
         END TYPE Observation_Conx
!
!
         TYPE Observation_ss
!
            INTEGER, DIMENSION(100)             :: num     ! Source/sink # of observation elements
            CHARACTER(LEN = 16), DIMENSION(100) :: name    ! Names of observation sources/sinks
!
         END TYPE Observation_ss
!
! ----------
! ...... Derived type variables
! ----------
!
         TYPE(Observation_Elem) :: ObsElem  ! The observation element list
         TYPE(Observation_Conx) :: ObsConx  !  the observation connection list
         TYPE(Observation_SS)   :: ObsSS    !  the observation SS list
!
! ----------
! ...... Integer variables
! ----------
!
         INTEGER :: NumobsElem,  &          ! The # of elements in the observation element loop
        &           NumObsConx,  &          ! The # of connections in the observation connections loop
        &           NumObsSS                ! The # of sources/sinks in the observation sources/sinks loop
!
! ..................MORE DECLARATIONS NECESSARY ???????????????
!
      END MODULE Time_Series_Parameters
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      MODULE General_Control_Parameters
!
         SAVE
!
! ----------
! ...... Double precision arrays
! ----------
!
         REAL(KIND = 8), DIMENSION(100) :: UserTimeStep, PrintTime
!
! ----------
! ...... Integer arrays
! ----------
!
         INTEGER, DIMENSION(24) :: MOP
!
! ----------
! ...... Double precision variables
! ----------
!
         REAL(KIND = 8) :: CPU_InitialTime, CPU_MatrixSolTime
         REAL(KIND = 8) :: time, TimeOrigin, TimeShift, SimulationTimeEnd
         REAL(KIND = 8) :: NewTimeStep, TimeStep, TimeIncrement, InitialTimeStep, MaxTimeStep, TimeStepMax_PostPrint, Dt_reducer
         REAL(KIND = 8) :: rel_convergence_crit, abs_convergence_crit, MaxResidual, derivative_increment
         REAL(KIND = 8) :: gravity, W_upstream, W_implicitness, W_NRIteration
!
! ----------
! ...... Integer variables
! ----------
!
! ...... MOP replacement
!
         INTEGER ::  Option_OutputAmount,            &  ! Describes the amount of output to be printed
        &            Option_Print_NRIterationInfo,   &  ! Was MOP(1)
        &            Option_Print_CyclingInfo,       &  ! Was MOP(2)
        &            Option_Print_JacobianInfo,      &  ! Was MOP(3)
        &            Option_Print_SourceSinkInfo,    &  ! Was MOP(4)
        &            Option_Print_EOSInfo,           &  ! Was MOP(5)
        &            Option_Print_SolverInfo,        &  ! Was MOP(6)
        &            Option_ThermalConductivity,     &
        &            Option_UpstreamWeighting,       &
        &            Option_InterfaceDensity,        &
        &            Option_RelPerm_CapPress,        &  ! Was MOP(8)
        &            Option_FluidComposition,        &  ! Was MOP(9)
        &            Option_AbsolutePerm,            &  ! Was MOP(11)
        &            Option_GenerationRates,         &  ! Was MOP(12)
        &            Option_PermAdjustment,          &  ! Was MOP(13)
        &            Option_GasSolubility,           &  ! Was MOP(14)
        &            Option_HeatExchange,            &  ! Was MOP(15)
        &            DoublingDtCriterion,            &  ! Was MOP(16)
        &            Option_BinGasDiffus,            &  ! Was MOP(17)
        &            Option_InConCheck,              &  ! Was MOP(20)
        &            Option_LineEqSolver,            &  ! Was MOP(21)
        &            Option_AddOutputFile,           &  ! Was MOP(22)
        &            Option_SelectDiffOption,        &  ! Was MOP(23)
        &            Option_MultiDiffuseFlux            ! Was MOP(24)

!
!
         INTEGER :: NumUserTimeSteps, NumPrintTimes
         INTEGER :: NumNRIterations, Max_NumNRIterations, Tot_NumNRIterations
         INTEGER :: NumTimeSteps, Tot_NumTimeSteps, Max_NumTimeSteps, PRINT_frequency
         INTEGER :: DefaultStateIndex
!
         INTEGER :: MaxResElemNum, MaxResEquNum
         INTEGER :: IGOOD, TrackElemNum
         INTEGER :: NumTimeStepCuts
!
         INTEGER :: SAVE_frequency, NumParamUpdates
!
! ----------
! ...... Character variables
! ----------
!
         CHARACTER(LEN = 4) :: Option_OutputFormat
         CHARACTER(LEN = 5) :: TrackElemName
!
! ----------
! ...... Logical variables
! ----------
!
         LOGICAL :: CoordNeed, NoFlowSimulation, NoVersion
         LOGICAL :: Flag_CheckInitialConditions, Flag_PrintInputs
         LOGICAL :: PrintOutNow, convergence
         LOGICAL :: RestartSimulation, RadiativeHeat
!
         LOGICAL :: unrealistic_conditions, realistic_result, impossible_solution
!
!
      END MODULE General_Control_Parameters
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      MODULE Solver_Parameters
!
         SAVE
!
! ----------
! ...... Double precision parameters
! ----------
!
         REAL(KIND = 8), PARAMETER :: seed = 1.0d-25
!
! ----------
! ...... Double precision variables
! ----------
!
         REAL(KIND = 8) :: Max_CGIterationRatio, CG_convergence_crit
!
! ----------
! ...... Integer parameters
! ----------
!
         INTEGER, PARAMETER :: iunit = 0
!
! ----------
! ...... Common integer variables
! ----------
!
         INTEGER :: MatrixSolver, Max_NumCGIterations
!
! ----------
! ...... Common logical variables
! ----------
!
         LOGICAL :: SolverBlok
!
!
      END MODULE Solver_Parameters
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      MODULE Input_Comments
!
         SAVE
!
! ----------
! ...... Double precision arrays
! ----------
!
         CHARACTER(LEN = 80), ALLOCATABLE, DIMENSION(:) :: comm
!
!
      END MODULE Input_Comments
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      MODULE Grid_Geometry
!
         SAVE
!
! ----------
! ...... Derived-type variable: Element geometry and property assignment
! ----------
!
         TYPE Element_Geometry
!
            CHARACTER(LEN = 5), POINTER, DIMENSION(:) :: name                       ! Element name
!
            CHARACTER(LEN = 1), POINTER, DIMENSION(:) :: activity                   ! Activity flag
!
            INTEGER(KIND = 2), POINTER, DIMENSION(:) :: MatNum                      ! Element material #
!
            REAL(KIND = 8), POINTER, DIMENSION(:) :: vol                            ! Element volume
!
            REAL(KIND = 4), POINTER, DIMENSION(:,:) :: coord            ! Coordinates of element center
!
         END TYPE Element_Geometry
!
! ----------
! ...... Derived-type variables: Connection geometry and property assignment
! ----------
!
         TYPE Connection_Geometry
!
            CHARACTER(LEN = 5), POINTER, DIMENSION(:) :: name1               ! Element names in the connection
            CHARACTER(LEN = 5), POINTER, DIMENSION(:) :: name2               ! Element names in the connection
!
            INTEGER, POINTER, DIMENSION(:) :: n1                                 ! Element numbers in the connection
            INTEGER, POINTER, DIMENSION(:) :: n2                                 ! Element numbers in the connection
            INTEGER, POINTER, DIMENSION(:) :: ki                                    ! Permeability index (direction on which the connection permeability is based)
!
            REAL(KIND = 8), POINTER, DIMENSION(:) :: d1                          ! Distance from interface (m)
            REAL(KIND = 8), POINTER, DIMENSION(:) :: d2                          ! Distance from interface (m)
            REAL(KIND = 8), POINTER, DIMENSION(:) :: area                           ! Interface area (m^2)
            REAL(KIND = 8), POINTER, DIMENSION(:) :: beta                           ! Cos(angle) with vertical
!
         END TYPE Connection_Geometry
!
! ----------
! ...... Derived-type arrays
! ----------
!
         TYPE(Element_Geometry) :: elem
!
         TYPE(Connection_Geometry) :: conx
!
! ----------
! ...... Double precision allocatable arrays
! ----------
!
         REAL(KIND = 4), ALLOCATABLE, DIMENSION(:) :: wt1, wt2
!
!
!***********************************************************************
!*                                                                     *
!*                         OBJECT PROCEDURES                           *
!*                                                                     *
!***********************************************************************
!
!
         PUBLIC :: Renumber_Elements, Determine_Connection_Elements
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
         SUBROUTINE Renumber_Elements( num_inact_elem, NumElem, NumElemTot, Max_NumElem, MemoryAlloc_Unit )
!
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*        Routine that renumbers elements by placing inactive          *
!*           ones after the end of the active element list             *
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
            INTEGER, INTENT(IN) :: num_inact_elem, NumElem, NumElemTot, Max_NumElem, MemoryAlloc_Unit
!
            TYPE(Element_Geometry) :: temp_elem
!
! -------------
! ......... Integer variables
! -------------
!
            INTEGER :: ier
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of <Renumber_Elements>
!

            ALLOCATE(temp_elem%name(num_inact_elem), temp_elem%activity(num_inact_elem), temp_elem%MatNum(num_inact_elem), temp_elem%vol(num_inact_elem), temp_elem%coord(num_inact_elem, 3), STAT=ier)          ! First, allocate some memory
!
! ......... Print explanatory comments
!
            IF(ier == 0) THEN
               WRITE(MemoryAlloc_Unit, FMT = 6101)
            ELSE
               WRITE(*, FMT = 6102)
               WRITE(MemoryAlloc_Unit, FMT = 6102)
               STOP
            END IF
!
!
! ......... Use the temporary array for storage of information
!
!
            temp_elem%name(1 : num_inact_elem : 1)     = elem%name(Max_NumElem : Max_NumElem-num_inact_elem+1 : -1)
            temp_elem%activity(1 : num_inact_elem : 1) = elem%activity(Max_NumElem : Max_NumElem-num_inact_elem+1 : -1)
            temp_elem%MatNum(1 : num_inact_elem : 1)   = elem%MatNum(Max_NumElem : Max_NumElem-num_inact_elem+1 : -1)
            temp_elem%vol(1 : num_inact_elem : 1)      = elem%vol(Max_NumElem : Max_NumElem-num_inact_elem+1 : -1)
            temp_elem%coord(1 : num_inact_elem : 1, :) = elem%coord(Max_NumElem : Max_NumElem-num_inact_elem+1 : -1, :)
!
!
! ......... Renumber inactive element starting on element NumElem+1 and assign properties
!
!
            !elem(NumElem+1 : NumElemTot) = temp_elem(1:num_inact_elem)%elem
            elem%name(NumElem+1 : NumElemTot)     = temp_elem%name(1 : num_inact_elem)
            elem%activity(NumElem+1 : NumElemTot) = temp_elem%activity(1 : num_inact_elem)
            elem%MatNum(NumElem+1 : NumElemTot)   = temp_elem%MatNum(1 : num_inact_elem)
            elem%vol(NumElem+1 : NumElemTot)      = temp_elem%vol(1 : num_inact_elem)
            elem%coord(NumElem+1 : NumElemTot, :) = temp_elem%coord(1 : num_inact_elem, :)
!
!
! ......... Finally, deallocate memory
!
            DEALLOCATE(temp_elem%name, temp_elem%activity, temp_elem%MatNum, temp_elem%vol, temp_elem%coord, STAT=ier)
!
! ......... Print explanatory comments
!
            IF(ier == 0) THEN
               WRITE(MemoryAlloc_Unit, FMT = 6103)
            ELSE
               WRITE(*, 6104)
               WRITE(MemoryAlloc_Unit, FMT = 6104)
               STOP
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
 6101 FORMAT(T2,'Memory allocation to array <temp_elem> in subroutine <Renumber_Elements> was successful')
!
 6102 FORMAT(//,20('ERROR-'),//,   &
     &       T2,'Memory allocation to arrays <temp_elem> in subroutine <Renumber_Elements> was unsuccessful',//,   &
     &       T2,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6103 FORMAT(T2,'Memory deallocation to array <temp_elem> in subroutine <Renumber_Elements> was successful')
!
 6104 FORMAT(//,20('ERROR-'),//,   &
     &       T2,'Memory deallocation to arrays <temp_elem> in subroutine <Renumber_Elements> was unsuccessful',//,   &
     &       T2,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End of <Renumber_Elements>
!
!
            RETURN
!
         END SUBROUTINE Renumber_Elements
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
         SUBROUTINE Determine_Connection_Elements( NumElemTot, NumConx )
!
#ifdef USE_CPP
            use, intrinsic :: iso_c_binding
#endif
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*        Routine for determining the neigbors of each element         *
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
            INTEGER, INTENT(IN) :: NumElemTot, NumConx
!
! -------------
! ......... Integer variables
! -------------
!
            INTEGER :: N_invalid_elem, n
!
#ifdef USE_TIMER
            real(KIND = 8) :: start, finish
#endif

#ifdef USE_CPP
            INTEGER ( c_int ) :: n1 ( NumConx )
            INTEGER ( c_int ) :: n2 ( NumConx )

            CHARACTER( LEN=5 ), target :: name1( NumConx )
            CHARACTER( LEN=5 ), target :: name2( NumConx )
            CHARACTER( LEN=5 ), target :: elems( NumElemTot )

            TYPE( c_ptr ), DIMENSION( NumConx ) :: name1_p
            TYPE( c_ptr ), DIMENSION( NumConx ) :: name2_p
            TYPE( c_ptr ), DIMENSION( NumElemTot ) :: elems_p

            interface
              subroutine check_edges ( num_elem, num_conx, n1, n2, name1_p, name2_p, elems_p) bind ( c )
                use iso_c_binding
                integer ( c_int ), VALUE :: num_elem
                integer ( c_int ), VALUE :: num_conx
                integer ( c_int ), dimension ( NumConx ) :: n1(*)
                integer ( c_int ), dimension ( NumConx ) :: n2(*)
                type ( c_ptr ), dimension ( NumConx ) :: name1_p(*)
                type ( c_ptr ), dimension ( NumConx ) :: name2_p(*)
                type ( c_ptr ), dimension ( NumElemTot ) :: elems_p(*)
              end subroutine check_edges
            end interface
#endif
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of <Determine_Connection_Elements>
!
!
! ......... Determine the element numbers of the connections
!
#ifdef USE_TIMER
            call CPU_Timing_Routine(start)
#endif

#ifdef USE_CPP
#ifdef USE_OMP
!$OMP PARALLEL
!$OMP DO schedule(auto)
#endif
            DO n = 1,NumConx
               name1_p(n) = C_LOC( conx%name1(n) )
               name2_p(n) = C_LOC( conx%name2(n) )
            END DO
#ifdef USE_OMP
!$OMP END DO NOWAIT
!$OMP DO schedule(auto)
#endif
            DO n = 1,NumElemTot
               elems_p(n) = C_LOC( elem%name(n) )
            END DO
#ifdef USE_OMP
!$OMP END DO
!$OMP END PARALLEL
#endif

            call check_edges(NumElemTot, NumConx, conx%n1, conx%n2, name1_p, name2_p, elems_p)
#else
            DO n = 1,NumElemTot
               WHERE(conx%name1(1:NumConx) == elem%name(n)) conx%n1(1:NumConx) = n   ! Identify the 1st elements in the connections
               WHERE(conx%name2(1:NumConx) == elem%name(n)) conx%n2(1:NumConx) = n   ! Identify the 2nd elements in the connections
            END DO

            N_invalid_elem = COUNT(conx%n1(1:NumConx) == 0) + COUNT(conx%n2(1:NumConx) == 0)  ! Total # of invalid elements in connections

            IF(N_invalid_elem /= 0) THEN
               DO n = 1,NumConx
                  IF(conx%n1(n) == 0) WRITE(*, FMT = 6008) conx%name1(n), n
                  IF(conx%n2(n) == 0) WRITE(*, FMT = 6008) conx%name2(n), n
               END DO
            END IF
#endif

#ifdef USE_TIMER
            call CPU_Timing_Routine(finish)

            write (*,*) __FILE__, ":", __LINE__, "edgelist check time: ", finish-start
#endif
!
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 6008 FORMAT(' Reference to unknown element "',A8,'" at connection number ',I6,' ==> Will ignore connection')
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End of <Determine_Connection_Elements>
!
!
            RETURN
!
         END SUBROUTINE Determine_Connection_Elements
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
      END MODULE Grid_Geometry
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      MODULE Element_Attributes
!
         SAVE
!
! ----------
! ...... Derived-type variable: Pressure and temperature in the elements
! ----------
!
         TYPE State_Conditions
!
            INTEGER(KIND = 1), DIMENSION(:,:), ALLOCATABLE :: index       ! Element phase state index
!
            REAL(KIND = 8), DIMENSION(:,:), ALLOCATABLE :: pres     ! Pressure
!
            REAL(KIND = 8), DIMENSION(:,:), ALLOCATABLE :: temp     ! Temperature
!
            REAL(KIND = 8), DIMENSION(:,:), ALLOCATABLE :: param    ! Auxilliary parameters stored at state point
!
         END TYPE State_Conditions
!
! ----------
! ...... Derived-type variable: Hydraulic properties
! ----------
!
         TYPE Hydraulic_properties
!
            REAL(KIND = 8)               :: porosity    ! Element porosity
            REAL(KIND = 8), DIMENSION(3) :: perm        ! Element permeability tensor components
!
         END TYPE Hydraulic_properties
!
! ----------
! ...... Derived-type variables: Secondary variables
! ----------
!
         TYPE Secondary_Variable
!
! ......... Accumulation-related
!
            REAL(KIND = 8), DIMENSION(:),   POINTER :: satur     ! Phase saturation
            REAL(KIND = 8), DIMENSION(:,:), POINTER :: MassFrac  ! Mass fraction of components in the phases
!
! ......... Thermophysical
!
            REAL(KIND = 8), POINTER, DIMENSION(:) :: density     ! Phase density (kg/m^3)
            REAL(KIND = 8), POINTER, DIMENSION(:) :: enthalpy    ! Phase specific enthalpy (J/Kg)
            REAL(KIND = 8), POINTER, DIMENSION(:) :: IntEnergy   ! Phase specific internal energy (J/Kg)
!
! ......... Mobility-related
!
            REAL(KIND = 8), POINTER, DIMENSION(:) :: RelPerm     ! Phase relative permeability
            REAL(KIND = 8), POINTER, DIMENSION(:) :: viscosity   ! Phase viscosity (Pa.s)
!
         END TYPE Secondary_Variable
!
! ----------
! ...... Derived-type arrays
! ----------
!
         TYPE(State_Conditions)                                  :: ElemState
!
         TYPE(Hydraulic_properties), ALLOCATABLE, DIMENSION(:,:) :: ElemMedia
!
         TYPE(Secondary_Variable), ALLOCATABLE, DIMENSION(:,:)   :: ElemProp
!
!
      END MODULE Element_Attributes
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      MODULE Connection_Attributes
!
         SAVE
!
! ----------
! ...... Derived-type variables: Connection
! ----------
!
         TYPE Connection_Flow
!
            REAL(KIND = 8) :: heat                                   ! Heat flow (W)
!
            REAL(KIND = 8), POINTER, DIMENSION(:) :: rate            ! Phase flow rate (kg/s)
            REAL(KIND = 8), POINTER, DIMENSION(:) :: DarcyVel        ! Phase Darcy velocities (m/s)
            REAL(KIND = 8), POINTER, DIMENSION(:) :: PoreVel         ! Phase pore velocities (m/s)
!
            REAL(KIND = 8), POINTER, DIMENSION(:,:) :: CompInPhase   ! Component flow in the phases (kg/s)
!
         END TYPE Connection_Flow
!
! ----------
! ...... Derived-type arrays
! ----------
!
         TYPE(Connection_Flow), ALLOCATABLE, DIMENSION(:) :: ConxFlow
!
!
      END MODULE Connection_Attributes
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!

      MODULE Initial_Conditions
!
         SAVE
!
! ----------
! ...... Double precision arrays
! ----------
!
         REAL(KIND = 8), DIMENSION(10) :: default_initial_cond
         REAL(KIND = 8), DIMENSION(10) :: XX
!
!
!
      END MODULE Initial_Conditions
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      MODULE Solution_Matrix_Arrays
!
         SAVE
!
! ----------
! ...... Double precision allocatable arrays
! ----------
!
!
         REAL(KIND = 8), ALLOCATABLE, DIMENSION(:) :: X, R, DX, DELX       !X=change in pressure, R=residual, DX=dist. b/w nodes,  DELX = derivate increment
!
         REAL(KIND = 8), ALLOCATABLE, DIMENSION(:) :: old_accum            !Time-step accumulation term
!
         REAL(KIND = 8), ALLOCATABLE, DIMENSION(:) :: CO
         REAL(KIND = 8), ALLOCATABLE, DIMENSION(:) :: rwork
!
         REAL(KIND = 8), ALLOCATABLE, DIMENSION(:,:) :: AB
!
! ----------
! ...... Integer allocatable arrays
! ----------
!
         INTEGER, ALLOCATABLE, DIMENSION(:) :: RowNum,ColNum
         INTEGER, ALLOCATABLE, DIMENSION(:) :: iwork
!
!
      END MODULE Solution_Matrix_Arrays
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      MODULE Sources_and_Sinks
!
         SAVE
!
! ----------
! ...... Derived-type variables
! ----------
!
         TYPE SourceSink
!
! ... Source/sink general attributes
!
            CHARACTER(LEN = 5) :: name          !i.e., Name of well
            CHARACTER(LEN = 5) :: ElemName      !i.e., Element name where well is
            CHARACTER(LEN = 5) :: type          !i.e., Type of well
!
            CHARACTER(LEN = 1) :: Rate_vs_time  !i.e., Production rate
!
            CHARACTER(LEN = 4) :: RateType      !i.e., Table of equation
!
            INTEGER :: ElemNum                  !Corresponds w/ ElemName
            INTEGER :: TypeIndx
!
! ......... Source/sink general parameters
!
            REAL(KIND = 8) :: rate, enth
!
! ......... Phase fractional flows at the well
!
            REAL(KIND = 8), DIMENSION(:), POINTER :: PhaseFracFlow
            REAL(KIND = 8), DIMENSION(:), POINTER :: PhaseVolRate
!
         END TYPE SourceSink
!
! ----------
! ...... Derived-type arrays
! ----------
!
         TYPE(SourceSink), ALLOCATABLE, DIMENSION(:) :: SS
!
!
      END MODULE Sources_and_Sinks
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      MODULE Geologic_Media_Parameters
!
         SAVE
!
! ----------
! ...... Derived-type variables
! ----------
!
         TYPE PFMedium                             ! Porous-fractured media
!
            REAL(KIND = 8) :: DensG,  &            ! The PF medium grain density (kg/m3)
     &                        Poros,  &            ! The PF medium porosity
     &                        SpcHt,  &            ! The PF medium specific heat (J/kg/C)
     &                        KThrW,  &            ! The thermal conductivity of the saturated PF medium (W/m/C)
     &                        KThrD                ! The thermal conductivity of the dry PF medium (W/m/C)
!
            REAL(KIND = 8) :: Compr,  &            ! The pore compressibility (1/Pa)
     &                        Expan,  &            ! The pore expansivity (1/C)
     &                        Tortu,  &            ! The PF medium tortuosity factor - diffusion
     &                        Klink                ! The Kinkenberg parameter (1/Pa) - low-flow rates
!
            REAL(KIND = 8), DIMENSION(3) :: Perm   ! Intrinsic permeabilities along principal axes (m^2)
!
            LOGICAL :: CompAsSolidSatFunction      ! Compressibility function
!
         END TYPE PFMedium
!
! ----------
! ...... Derived-type arrays
! ----------
!
         TYPE(PFMedium),  ALLOCATABLE, DIMENSION(:) :: media     !i.e., media(m)%Poros
!
! ----------
! ...... Integer allocatable arrays
! ----------
!
         INTEGER(KIND = 1), ALLOCATABLE, DIMENSION(:) :: StateI
         INTEGER(KIND = 1), ALLOCATABLE, DIMENSION(:) :: StateIndex
!
! ----------
! ...... Character allocatable arrays
! ----------
!
         CHARACTER(LEN = 5), ALLOCATABLE, DIMENSION(:) :: MediumName, Rk_name
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
      END MODULE Geologic_Media_Parameters
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE Core_Memory_Allocation
!
         USE EOS_Parameters
         USE EOS_Default_Parameters
!
         USE Basic_Parameters
         USE General_External_File_Units
!
         USE General_Control_Parameters
!
         USE Grid_Geometry
         USE Element_Attributes
         USE Connection_Attributes
!
         USE Solution_Matrix_Arrays
         USE Sources_and_Sinks
!
         USE Geologic_Media_Parameters
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*      MAIN ROUTINE FOR ALLOCATING MEMORY TO MOST FTSim ARRAYS        *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
         IMPLICIT NONE
!
! ----------
! ...... Integer Arrays
! ----------
!
         INTEGER, DIMENSION(0:30) :: ierror = 99999
!
! ----------
! ...... Integer variables
! ----------
!
         INTEGER :: ii, i, j
!
         INTEGER :: number_of_components, number_of_equations
         INTEGER :: Max_number_of_elements, Max_number_of_connections
         INTEGER :: Max_number_of_sources_and_sinks, Max_number_of_geologic_media
!
! ----------
! ...... Character variables
! ----------
!
         CHARACTER(LEN = 15) :: B_name = '               '
         CHARACTER(LEN = 26) :: HEADR
         CHARACTER(LEN =  9) :: Memory_header
!
! ----------
! ...... Logical variables
! ----------
!
         LOGICAL exists, NewInputFormat
!
         LOGICAL Flag_binary_diffusion
#ifdef USE_TIMER
         real(KIND = 8) :: start, finish
#endif
!
! ----------
! ...... Namelists
! ----------
!
         NAMELIST/Basic_Parameter_Definitions/ EOS_Name,                      &  ! The name of the EOS segment compiled with FTSim
     &                                         number_of_components,          &  ! Number of mass components
     &                                         number_of_equations,           &  ! Number of equations per gridblock
     &                                         Flag_binary_diffusion             ! Logical variable determining whether diffusion is to be considered
!
         NAMELIST/System_Specifications/ coordinate_system,                   &  ! Coordinate system
     &                                   Max_number_of_elements,              &  ! Maximum number of elements (cells)
     &                                   Max_number_of_connections,           &  ! Maximum number of connections (cells)
     &                                   Max_number_of_sources_and_sinks,     &  ! Maximum number of sources and sinks
     &                                   Max_number_of_geologic_media            ! Maximum number of geological media
!
         SAVE
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of Core_Memory_Allocation
!
!
         WRITE(VERS_Unit,6000)
!
!
!***********************************************************************
!*                                                                     *
!*              OPEN FILE TO DOCUMENT MEMORY ALLOCATION                *
!*                                                                     *
!***********************************************************************
!
!
         INQUIRE(FILE='ALLOC',EXIST=exists)
         IF(exists) THEN
            OPEN(UNIT = MemoryAlloc_Unit, FILE = 'ALLOC', STATUS = 'OLD')
            REWIND (UNIT = MemoryAlloc_Unit)
         ELSE
            OPEN(UNIT = MemoryAlloc_Unit, FILE = 'ALLOC', STATUS = 'NEW')
         END IF
!
!
!***********************************************************************
!*                                                                     *
!*              READ NUMBERS OF ELEMENTS AND CONNECTIONS               *
!*                                                                     *
!***********************************************************************
!
!
         READ(*,*) HEADR
         B_name        = TRIM(ADJUSTL(HEADR))
         Memory_header = B_name(1:9)
!
!
         IF( Memory_header /= '>>>memory'   .AND. Memory_header /= '>>>Memory'   .AND. Memory_header /= '>>>MEMORY' .AND.   &
    &        Memory_header(1:6) /= 'memory' .AND. Memory_header(1:6) /= 'Memory' .AND. Memory_header(1:6) /= 'MEMORY'   )   &
    &    THEN
            WRITE(*, FMT = 6002)
            WRITE(UNIT = MemoryAlloc_Unit, FMT = 6002)
            STOP
         END IF
!
! ...... Determine the input format style
!
         IF( Memory_header == '>>>memory' .OR. Memory_header == '>>>Memory' .OR. Memory_header == '>>>MEMORY' ) THEN
            NewInputFormat = .TRUE.
         ELSE
            NewInputFormat = .FALSE.
         END IF
!
!
!***********************************************************************
!*                                                                     *
!*  READ THE NUMBER OF COMPONENTS, EQUATIONS, PHASES, SEC. VARIABLES   *
!*                                                                     *
!***********************************************************************
!
!
         IF_InputFormat1: IF( NewInputFormat .EQV. .FALSE. ) THEN
!
! ......... Reading data using the old input format
!
            READ(*,*) B_name
            EOS_Name = TRIM(ADJUSTL(B_name))
!
            READ(*,*) NumCom, NumEqu, binary_diffusion
!
         ELSE
!
! ......... Reading data using the new input format (NAMELIST-based)
!
            READ(*, NML = Basic_Parameter_Definitions, IOSTAT = ierror(0))
!
            NumCom = number_of_components                                     !If isothermal, NumEqu = NumCom
            NumEqu = number_of_equations                                      !If nonisothermal, NumEqu = NumCom+1 (per gridblock)
!
            binary_diffusion = Flag_binary_diffusion
!
! ......... Stop if there is a problem reading the fundamental parameters
!
            IF (ierror(0) /= 0) THEN
               WRITE (*, FMT = 6100)
               STOP
            END IF
!
         END IF IF_InputFormat1
!
!
!***********************************************************************
!*                                                                     *
!*   DETERMINE DEFAULT PARAMETERS FOR THE VARIOUS EQUATIONS OF STATE   *
!*                                                                     *
!***********************************************************************
!
!
         CALL EOS_Defaults(ADJUSTL(EOS_Name), NumCom, NumEqu, binary_diffusion,                                        &
     &                                        Max_NumMassComp, Min_NumMassComp, Max_NumEquations, Min_NumEquations,    &
     &                                        Max_NumPhases, Max_NumMobPhases, NumStateParam, NumAuxParam,             &
     &                                        VERS_Unit, MemoryAlloc_Unit )
!
! ...... Determine the number of mobile phases
!
         NumPhases    = Max_NumPhases
         NumMobPhases = Max_NumMobPhases
!
! ...... Define some additional parameters
!
         NumEquPlus1  = NumEqu + 1
         NumEquSquare = NumEqu*NumEqu
         NumPerturb   = 2*NumEqu + 1
         NumComPlus1  = NumCom + 1
!
! ...... Define whether the simulation is conducted isothermally
!
         IF(NumEqu == NumCom) THEN
            nonisothermal_conditions = .FALSE.
         ELSE
            nonisothermal_conditions = .TRUE.
         END IF
!
! ----------
! ...... Determine the primary variable list (to be included in the INCON and SAVE files)
! ----------
!
         CALL EOS_Primary_Variable_List( ADJUSTL(EOS_Name), NumCom, NumEqu, Number_of_states, &
     &                                   EOS_variables, VERS_Unit, MemoryAlloc_Unit )
!
! ----------
! ...... Number of secondary parameters (to be stored in the DP array)
! ----------
!
         IF(binary_diffusion) THEN
            NumSecondaryVar = 10*Max_NumPhases + 1 + NumAuxParam
         ELSE
            NumSecondaryVar = 8*Max_NumPhases + 1 + NumAuxParam
         END IF
!
!
!***********************************************************************
!*                                                                     *
!*       SYSTEM DESCRIPTION: GRID, SOURCES & SINKS, MEDIA, ETC.        *
!*                                                                     *
!***********************************************************************
!
!
         IF_InputFormat2: IF( NewInputFormat .EQV. .FALSE. ) THEN
!
! -------------
! ......... Old Input Style: Read the number of elements, connections and number of characters in the element name
! -------------
!
            READ(*,*) coordinate_system, Max_NumElem, Max_NumConx
            READ(*,*) Max_NumSS
            READ(*,*) Max_NumMedia
!
         ELSE
!
! -------------
! ......... New Input Style: Read the number of elements, connections and number of characters in the element name
! -------------
!
            READ(*, NML = System_Specifications, IOSTAT = ierror(0))
!
! ......... Stop if there is a problem reading the fundamental parameters
!
            IF (ierror(0) /= 0) THEN
               WRITE (*, FMT = 6105)
               STOP
            END IF
!
            Max_NumElem  = Max_number_of_elements
            Max_NumConx  = Max_number_of_connections
            Max_NumSS    = Max_number_of_sources_and_sinks
            Max_NumMedia = Max_number_of_geologic_media
!
         END IF IF_InputFormat2
!
         SELECT CASE(coordinate_system(1:3))
         CASE('Car', 'car', 'CAR')
            coordinate_system(1:3) = 'Car'
         CASE( 'Cyl', 'cyl', 'CYL' )
            coordinate_system(1:3) = 'Cyl'
         CASE DEFAULT
            WRITE(*,6010) coordinate_system
         END SELECT
!
! ...... Add a couple of extra elements and connections (to avoid array boundary errors)
!
         Max_NumElem = Max_NumElem + 2
         Max_NumConx = Max_NumConx + 2
!
!
!***********************************************************************
!*                                                                     *
!*            MEMORY ALLOCATION TO ELEMENT-RELATED ARRAYS              *
!*                                                                     *
!***********************************************************************
!
!
         ElemNameLength = 5
!
! ----------
! ...... Allocate memory to derived-type arrays
! ----------
!
!
         ALLOCATE(elem%name(Max_NumElem), elem%activity(Max_NumElem), elem%MatNum(Max_NumElem), elem%vol(Max_NumElem), elem%coord(Max_NumElem, 3), ElemMedia(Max_NumElem, -2:NumEqu), STAT=ierror(1))    !Define arrays for elements, element state, element media (n,k) at original and successive timesteps
!
! ----------
! ...... Allocate memory for arrays within the variable <ElemState> (derived type)
! ----------
!
         ALLOCATE(ElemState%index(Max_NumElem, -1:0), ElemState%pres(Max_NumElem,-2:NumEqu), ElemState%temp(Max_NumElem,-2:NumEqu), ElemState%param(Max_NumElem,NumStateParam), STAT = ierror(2))
! ----------
! ...... Allocate memory to derived-type arrays - Secondary variables
! ----------
!
         ALLOCATE(ElemProp(Max_NumElem,0:NumEqu), STAT=ierror(3))                   !Record only at current state and next increment
!
! ----------
! ...... Allocate memory for arrays within the variable <ElemProp> (derived type)
! ----------
!
         DO_j1 : DO j = 0,NumEqu
            DO_i1 : DO i = 1,Max_NumElem
!
               ALLOCATE(ElemProp(i,j)%satur(Max_NumPhases),                      &
     &                  ElemProp(i,j)%MassFrac(Max_NumMassComp, Max_NumPhases),  &
     &                  ElemProp(i,j)%density(Max_NumPhases),                    &
     &                  ElemProp(i,j)%enthalpy(Max_NumPhases),                   &
     &                  ElemProp(i,j)%IntEnergy(Max_NumPhases),                  &
     &                  ElemProp(i,j)%RelPerm(Max_NumMobPhases),                 &
     &                  ElemProp(i,j)%viscosity(Max_NumMobPhases),               &
     &                  STAT = ierror(4))
!
            END DO DO_i1
         END DO DO_j1
!
! ----------
! ...... Allocate memory to the integer arrays that provide the element starting locations in the arrays of the unknowns
! ----------
!
         ALLOCATE(Loc(Max_NumElem), Locp(Max_NumElem), STAT=ierror(8))
!
! ...... Define the locator arrays
!
         FORALL (i = 1 : Max_NumElem)
            Loc(i)  = (i-1)*NumEqu
            Locp(i) = (i-1)*NumComPlus1
         END FORALL
!
!
!***********************************************************************
!*                                                                     *
!*          MEMORY ALLOCATION TO CONNECTION-RELATED ARRAYS             *
!*                                                                     *
!***********************************************************************
!
! ----------
! ...... Allocate memory to derived-type arrays
! ----------
!
         ALLOCATE( conx%name1(Max_NumConx), conx%name2(Max_NumConx), conx%n1(Max_NumConx), conx%n2(Max_NumConx), conx%ki(Max_NumConx), conx%d1(Max_NumConx), conx%d2(Max_NumConx), conx%area(Max_NumConx), conx%beta(Max_NumConx), ConxFlow(Max_NumConx), STAT=ierror(9) )
!
! ----------
! ...... Allocate memory to arrays within the derived type: Conx
! ----------
!
         DO_Conx1 : DO i = 1,Max_NumConx
!
            ALLOCATE(ConxFlow(i)%rate(Max_NumMobPhases),     STAT = ierror(10))
            ALLOCATE(ConxFlow(i)%DarcyVel(Max_NumMobPhases), ConxFlow(i)%PoreVel(Max_NumMobPhases),  STAT = ierror(11))
!
            ALLOCATE(ConxFlow(i)%CompInPhase(Max_NumMassComp, Max_NumMobPhases),  STAT = ierror(12))
!
         END DO DO_Conx1
!
!
!***********************************************************************
!*                                                                     *
!*            MEMORY ALLOCATION TO JACOBIAN MATRIX ARRAYS              *
!*                                                                     *
!***********************************************************************
!
!
         Max_NumPrimaryVar = (Max_NumMassComp + 1)*Max_NumElem
         MatrixOrder       =  NumEqu*Max_NumElem
!
! ----------
! ...... Allocate memory to double precision allocatable arrays
! ----------
!
         ALLOCATE(DX(Max_NumPrimaryVar),    STAT=ierror(15))
         ALLOCATE(DELX(Max_NumPrimaryVar),  STAT=ierror(16))
         ALLOCATE(old_accum(MatrixOrder),   STAT=ierror(17))
!
!
!
         ALLOCATE(X(Max_NumPrimaryVar), R(MatrixOrder + 1), STAT = ierror(18))         !X= Primary variables (pres,temp,sat)
!
!
!***********************************************************************
!*                                                                     *
!*          MEMORY ALLOCATION TO SOURCE/SINK-RELATED ARRAYS            *
!*                                                                     *
!***********************************************************************
!
!
         Max_NumSS = Max_NumSS + 2                                          ! Add a couple of extra sources/sinks (to avoid array boundary errors)
!
! ----------
! ...... Allocate memory to derived-type
! ----------
!
         ALLOCATE(SS(Max_NumSS), STAT = ierror(20))
!
         IF(Max_NumPhases > 1) THEN
            DO_SS : DO i = 1,Max_NumSS
               ALLOCATE(SS(i)%PhaseFracFlow(Max_NumMobPhases), STAT = ierror(21))
               ALLOCATE(SS(i)%PhaseVolRate(Max_NumMobPhases),  STAT = ierror(22))
            END DO DO_SS
         END IF
!
!
!***********************************************************************
!*                                                                     *
!*           MEMORY ALLOCATION TO ARRAYS OF ROCK PROPERTIES            *
!*                                                                     *
!***********************************************************************
!
!
         Max_NumMedia = Max_NumMedia + 2                                  !  Add a couple of extra media (to avoid array boundary errors)
!
! ----------
! ...... Allocate memory to derived-type allocatable arrays
! ----------
!
         ALLOCATE(media(Max_NumMedia), STAT=ierror(25))
!
!
!
!***********************************************************************
!*                                                                     *
!*         ENSURE PROPER MEMORY ALLOCATION - STOP IF PROBEMS           *
!*                                                                     *
!***********************************************************************
!
! Output memory allocation problem if exists
! Cannot put write statement in array process, so use do loop
!
         DO_InitCheck: DO ii=0,27
!
            IF(ierror(ii) == 0) THEN
               WRITE(MemoryAlloc_Unit,6015) ii
            ELSE IF(ierror(ii) == 99999) THEN
               CONTINUE
            ELSE
               WRITE(*,6020)  ii
               WRITE(MemoryAlloc_Unit,6020) ii
               STOP
            END IF
!
         END DO DO_InitCheck
!
!
!***********************************************************************
!*                                                                     *
!*                      END OF THE DATA BLOCK                          *
!*                                                                     *
!***********************************************************************
!
!
      IF_FormatStyle: IF( NewInputFormat ) THEN
!
         READ( *, FMT = '(A3)', IOSTAT = ierror(0) ) DataBlockEndSpecifier
!
!-----------
! ...... Stop if there is a problem reading the data block end specifies
!-----------
!
         IF( (ierror(0) == 0) .AND. (DataBlockEndSpecifier == '<<<') ) THEN
            RETURN
         ELSE
            WRITE (*, FMT = 6500)
            STOP
         END IF
!
      END IF IF_FormatStyle
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 6000 FORMAT(/,'Core_Memory_Allocation',T50,'[v1.0,   29 May      2008]',/,        &
     &         ':::::   Allocate memory to most arrays based on size provided by input parameters')
!
 6002 FORMAT(//,20('ERROR-'),//,T33,                                              &
     &             '       S I M U L A T I O N   A B O R T E D',/,                &
     &         T30,'The header <MEMORY> is missing at the top of the input file', &
     &       /,T32,'               CORRECT AND TRY AGAIN',                        &
     &       //,20('ERROR-'))
!
 6010 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'The coordinate system (variable <coordinate_system> = ',a11, &
     &          ' in subroutine <Core_Memory_Allocation>) is unavailable',//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6015 FORMAT(T2,' Memory allocation at point ',i3,' in subroutine "Core_Memory_Allocation" was successful')
!
 6020 FORMAT(//,20('ERROR-'),//,                                          &
     &       T2,' Memory allocation at point ',i3,' in subroutine "Core_Memory_Allocation" was unsuccessful',//,   &
     &       T2,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',        &
     &       //,20('ERROR-'))
!
 6100 FORMAT(//,22('ERROR-'),//,T33,  &
     &             '       S I M U L A T I O N   A B O R T E D',/,                              &
     &         T12,'There is a problem reading the namelist <Basic_Parameter_Definitions> in the <Memory_Allocation> data block', &
     &       /,T32,'               CORRECT AND TRY AGAIN',                                      &
     &       //,22('ERROR-'))
!
 6105 FORMAT(//,22('ERROR-'),//,T33,  &
     &             '       S I M U L A T I O N   A B O R T E D',/,                              &
     &         T18,'There is a problem reading the namelist <System_Specifications> in the <Memory_Allocation> data block', &
     &       /,T32,'               CORRECT AND TRY AGAIN',                                      &
     &       //,22('ERROR-'))
!
 6500 FORMAT(//,22('ERROR-'),//,T33,  &
     &             '       S I M U L A T I O N   A B O R T E D',/,                                     &
     &         T10,'There is a problem reading the data block end specifier <DataBlockEndSpecifier> ', &
     &             'in the <Memory_Allocation> block',                                                 &
     &       /,T32,'               CORRECT AND TRY AGAIN',                                             &
     &       //,22('ERROR-'))
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End of Core_Memory_Allocation
!
!
      END SUBROUTINE Core_Memory_Allocation
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE Solver_Memory_Allocation
!
! ... Modules to be used
!
         USE Basic_Parameters
         USE General_External_File_Units

         USE Solution_Matrix_Arrays
         USE General_Control_Parameters
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*  ROUTINE FOR MEMORY ALLOCATION TO SOLUTION MATRIX AND WORK ARRAYS   *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
         IMPLICIT NONE
!
! ----------
! ...... Integer Arrays
! ----------
!
         INTEGER, DIMENSION(26:30) :: ierror = 99999
!
! ----------
! ...... Integer Variables
! ----------
!
!
         INTEGER :: ii
!
         SAVE
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of Solver_Memory_Allocation
!
!
         WRITE(VERS_Unit,6000)
         WRITE(MemoryAlloc_Unit,6005)
!
! ...... Account for the solver memory needs
!
         Max_NumNonZeros = (Max_NumElem+2*Max_NumConx)*NumEquSquare
!
! ...... Account for the solver needs
!
         real_work_array_size = MAX( 1000 + Max_NumNonZeros + 10*MatrixOrder, 9*Max_NumConx+1000 )
         intg_work_array_size = 100 + Max_NumNonZeros + 5*MatrixOrder
!
!***********************************************************************
!*                                                                     *
!*        MEMORY ALLOCATION TO SOLUTION MATRIX AND WORK ARRAYS         *
!*                                                                     *
!***********************************************************************
!
! ----------
! ...... Allocate memory to double precision allocatable arrays
! ----------
!
! CO = Jacobian matrix coeffs., i.e., CO is Jacobian matrix in 1-D array
! 1-D most efficient for banded matrix
!
         ALLOCATE(CO(Max_NumNonZeros+1),       STAT=ierror(26))
         ALLOCATE(RWORK(real_work_array_size), STAT=ierror(27))
!
! ----------
! ...... Allocate memory to integer allocatable arrays
! ----------
!
         ALLOCATE(RowNum(Max_NumNonZeros+1), STAT=ierror(28))
         ALLOCATE(ColNum(Max_NumNonZeros+1), STAT=ierror(29))
!
         ALLOCATE(iwork(intg_work_array_size), STAT=ierror(30))
!
!
!***********************************************************************
!*                                                                     *
!*         ENSURE PROPER MEMORY ALLOCATION - STOP IF PROBEMS           *
!*                                                                     *
!***********************************************************************
!
!
         DO_InitCheck: DO ii=26,30
!
            IF(ierror(ii) == 0) THEN
               WRITE(MemoryAlloc_Unit,6001) ii
            ELSE IF(ierror(ii) == 99999) THEN
               CONTINUE
            ELSE
               WRITE(*,6002)  ii
               WRITE(MemoryAlloc_Unit,6002) ii
               STOP
            END IF
!
         END DO DO_InitCheck
!
! ----------
! ...... Print information on the size of the arrays
! ----------
!
         WRITE(*,6010) Max_NumNonZeros, Max_NumNonZeros+1, real_work_array_size, intg_work_array_size
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 6000 FORMAT(/,'Solver_Memory_Allocation',T50,'[v1.0,   29 May      2008]',/,   &
     &         ':::::   Allocate memory to the Jacobian  and solution matrix arrays')
!
 6001 FORMAT(T2,' Memory allocation at point ',i3,' in subroutine "Allocate_Solver_Memory" was successful')
!
 6002 FORMAT(//,20('ERROR-'),//,   &
     &       T2,' Memory allocation at point ',i3,' in subroutine "Allocate_Solver_Memory" was unsuccessful',//,   &
     &       T2,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
 6005 FORMAT(/)
!
 6010 FORMAT(' ',131('='),//,   &
     &       ' Maximum number of Jacobian matrix elements:',16X,' Max_NumNonZeros  =',I9,//,   &
     &       ' Large linear equation arrays - Length of "Row", "Col" and "CO" is             "Max_NumNonZeros+1"      = ',I10,//, &
     &       ' Maximum length of the Conjugate Gradient REAL*8 work array "real_work_array" is "real_work_array_size" = ',i10,/, &
     &       ' Maximum length of the Conjugate Gradient integer array "intg_work_array is" "intg_work_array_size"     = ',i10,//,&
     &       ' ',131('='))
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End of Solver_Memory_Allocation
!
!
      END SUBROUTINE Solver_Memory_Allocation
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE Allocate_Temp_Input_Memory
!
! ... Modules to be used
!
         USE Input_Comments
         USE Basic_Parameters
         USE General_External_File_Units

!
         USE General_Control_Parameters
!
         USE Geologic_Media_Parameters
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*       ROUTINE FOR MEMORY ALLOCATION TO TEMPORARY INPUT ARRAYS       *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
         IMPLICIT NONE
!
! ----------
! ...... Integer Arrays
! ----------
!
         INTEGER, DIMENSION(41:45) :: ierror = 99999
!
! ----------
! ...... Integer Variables
! ----------
!
!
         INTEGER :: ii
!
         SAVE
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of Allocate_Temp_Input_Memory
!
!
         WRITE(VERS_Unit,6000)
         WRITE(MemoryAlloc_Unit,6005)
!
         ALLOCATE(Rk_name(Max_NumMedia), MediumName(Max_NumMedia), STAT=ierror(41))
!
         ALLOCATE(StateI(Max_NumMedia),  StateIndex(Max_NumMedia), STAT=ierror(42))
!
         ALLOCATE(comm(50), STAT=ierror(43))
!
!
!***********************************************************************
!*                                                                     *
!*         ENSURE PROPER MEMORY ALLOCATION - STOP IF PROBEMS           *
!*                                                                     *
!***********************************************************************
!
!
         DO_InitCheck: DO ii=41,45
!
            IF(ierror(ii) == 0) THEN
               WRITE(MemoryAlloc_Unit,6001) ii
            ELSE IF(ierror(ii) == 99999) THEN
               CONTINUE
            ELSE
               WRITE(*,6002)  ii
               WRITE(MemoryAlloc_Unit,6002) ii
               STOP
            END IF
!
         END DO DO_InitCheck
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 6000 FORMAT(/,'Allocate_Temp_Input_Memory',T50,'[v1.0,   29 May      2008]',/,   &
     &         ':::::   Allocate memory to to temporary input arrays')
!
 6001 FORMAT(T2,' Memory allocation at point ',i3,' in subroutine "Allocate_Temp_Input_Memory" was successful')
 6002 FORMAT(//,20('ERROR-'),//,   &
     &       T2,' Memory allocation at point ',i3,   &
     &          ' in subroutine "Allocate_Temp_Input_Memory" was unsuccessful',   &
     &       //,T2,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
 6005 FORMAT(/)
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of Allocate_Temp_Input_Memory
!
!
      END SUBROUTINE Allocate_Temp_Input_Memory
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE Deallocate_Unneeded_Arrays
!
! ... Modules to be used
!
         USE Input_Comments
!
         USE Basic_Parameters
         USE General_External_File_Units

         USE General_Control_Parameters
!
         USE Grid_Geometry
!
         USE Geologic_Media_Parameters
         USE Sources_and_Sinks
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*     ROUTINE FOR MEMORY DEALLOCATION FROM TEMPORARY INPUT ARRAYS     *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
         IMPLICIT NONE
!
! ----------
! ...... Integer Arrays
! ----------
!
         INTEGER, DIMENSION(41:45) :: ierror = 99999
!
! ----------
! ...... Integer Variables
! ----------
!
!
         INTEGER :: ii,ierG = 0
!
         SAVE
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of Deallocate_Unneeded_Arrays
!
!
         WRITE(VERS_Unit,6000)
         WRITE(MemoryAlloc_Unit,6005)
!
         DEALLOCATE(Rk_name,MediumName, STAT=ierror(41))  ! PF medium name
         DEALLOCATE(StateI,StateIndex,  STAT=ierror(42))  ! PF state index
!
         DEALLOCATE(comm,    STAT=ierror(43))             ! Comments in input file
!
!
!***********************************************************************
!*                                                                     *
!*        ENSURE PROPER MEMORY DEALLOCATION - STOP IF PROBEMS          *
!*                                                                     *
!***********************************************************************
!
!
         DO_InitCheck: DO ii=41,45
!
            IF(ierror(ii) == 0) THEN
               WRITE(MemoryAlloc_Unit,6001) ii
            ELSE IF(ierror(ii) == 99999) THEN
               CONTINUE
            ELSE
               WRITE(*,6002)  ii                         ! Unsuccesful deallocation
               WRITE(MemoryAlloc_Unit,6002) ii
               STOP
            END IF
!
         END DO DO_InitCheck
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 6000 FORMAT(/,'Deallocate_Unneeded_Arrays',T50,'[v1.0,   29 May      2008]',/,   &
     &         ':::::   Deallocation of memory from all temporary and unnecessary arrays')
!
 6001 FORMAT(T2,' Memory deallocation at point ',i3,' in subroutine "Deallocate_Unneeded_Arrays" was successful')
 6002 FORMAT(//,20('ERROR-'),//,   &
     &       T2,' Memory deallocation at point ',i3,' in subroutine "Deallocate_Unneeded_Arrays" was unsuccessful',   &
     &       //,T2,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
 6005 FORMAT(/)
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of Deallocate_Unneeded_Arrays
!
!
      END SUBROUTINE Deallocate_Unneeded_Arrays
!
!
!
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
