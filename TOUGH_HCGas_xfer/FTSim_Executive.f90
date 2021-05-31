!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>                                                                     >
!>                                                                     >
!>         T_Executive.f95: Executive code unit including the          >
!>          executive routines common to all TFx simulations           >
!>                                                                     >
!>                                                                     >
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  The executive segment of FTSim.  It includes the procedures that
!  advance the time in the simulation process, estimate the time-step
!  for optimum performance, populate the matrix arrays and invoke the
!  solvers of the Jacobian, invoke special linear algebra for matrix
!  pre-processing in cases of very demanding linear algebra problems,
!  compute mass and energy balances, compute rates in sources and sinks,
!  compute binary diffusion coefficients, write special output files,
!  and conduct other miscellaneous operations.
!
!
      SUBROUTINE Simulation_Cycle
!
! ...... Modules to be used
!
         USE EOS_Routine_Selector
!
         USE Basic_Parameters
         USE General_Control_Parameters
         USE General_External_File_Units
!
         USE Grid_Geometry
         USE Element_Attributes
!
         USE Solution_Matrix_Arrays
!
         USE Geological_Media_Properties
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*          EXECUTIVE ROUTINE THAT PERFORMS THE SIMULATION             *
!*                     WHILE ADVANCING IN TIME                         *
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
      REAL(KIND = 8) :: CPU_time, Dt_solv, T_dd
!
! -------
! ... Integer variables
! -------
!
      INTEGER :: NSTL, KIT
      INTEGER :: n_t, n_t0, i_n, NIT, NIT1
!
! -------
! ... Logical variables
! -------
!
      LOGICAL :: ConvergenceFailure, DtCut_Flag, PrintNow_Flag
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of <Simulation_Cycle>
!
!
      WRITE(VERS_Unit,6000)
!
!
!***********************************************************************
!*                                                                     *
!*                           INITIALIZATIONS                           *
!*                                                                     *
!***********************************************************************
!
!
      NumTimeSteps = 0
!
      IGOOD = 0
!
      time = TimeShift
!
      NewTimeStep = InitialTimeStep
!
      NSTL = (TrackElemNum-1)*NumComPlus1
!
      NumNRIterations = 0
!
      KIT  = 0
!
      convergence             = .FALSE.       ! Flag denoting convergence of the Newton-Raphson iteration - Simulation continuing
!
      ConvergenceFailure      = .FALSE.       ! Flag denoting convergence failure leading to abrupt simulation end
!
      unrealistic_conditions  = .FALSE.       ! Flag denoting that the conditions in any of elements are unrealistic
!
      impossible_solution     = .FALSE.       ! Flag denoting that the solution of the Jacobian matrix is impossible
!
! ----------
! ... Initialize thermophysical properties
! ----------
!
      CALL Equation_of_State                  ! Calling the Equation-of-State routines
!
! ... Computing changes to porosity, permeability and capillary pressure
!
      CALL Hydraulic_Property_Changes
!
      CALL Mass_and_Energy_Balance(0)         ! Compute mass and heat balance
!
!
!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>                                                                     >
!>                       MAIN TIME STEP LOOP                           >
!>                                                                     >
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
      n_t0 = Tot_NumTimeSteps + 1                                            ! Beginning value of the cumulative number of time steps (including restarts)
!
      DO_TimeSteps: DO n_t = n_t0, Max_NumTimeSteps
!
         Tot_NumTimeSteps = n_t                                              ! Update <Tot_NumTimeSteps>
!
         NumTimeSteps = NumTimeSteps + 1                                     ! Update the number of time steps (new time step)
!
         IF(Option_Print_NRIterationInfo /= 0) THEN
            WRITE(*,6001) NumTimeSteps, Tot_NumTimeSteps, NewTimeStep, time  ! Print explanatory messsage
         END IF
!
         NumTimeStepCuts = 0                                                 ! Initializing the timestep reduction counter
!
!
!***********************************************************************
!*                                                                     *
!*                 INITIALIZATIONS FOR NEW TIME STEP                   *
!*                                                                     *
!***********************************************************************
!
!
 2000    NumNRIterations = 0                                       ! The counter of the Newtonian iterations in this time step
!
         PrintOutNow   = .FALSE.                                   ! Flag indicating print-out at a requested time
!
         PrintNow_Flag = .FALSE.                                   ! Flag indicating immediate print-out
!
! -------------
! ...... Set the time-step (Dt) size & adjust related variables
! -------------
!
         IF(NumPrintTimes /= 0) CALL Adjust_Timestep_for_Printout  ! Adjust Dt to coincide with desired print-out times
!
         TimeStep = NewTimeStep                                    ! Set Dt for the ensueing computations

         convergence = .FALSE.                                     ! Set the convergence flag
!
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*                   THE NEWTONIAN ITERATION LOOP                      *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
!
         DO_NRIterations: DO i_n = 1,Max_NumNRIterations+1
!
! ----------------
! ......... Initialize changes in the primary variables for this iteration
! ----------------
!
            IF(NumNRIterations == 0) DX = 0.0d0  ! Whole array operation
!
! ----------------
! ......... Update counters
! ----------------
!
            NumNRIterations     = NumNRIterations + 1      ! Update the Newtonian iteration counter in this time step
            Tot_NumNRIterations = Tot_NumNRIterations + 1  ! Update the cumulative Newtonian iteration counter since t = 0
!
            TimeIncrement = W_implicitness*TimeStep        ! Compute the implicitness-weighted time increment
!
!
!***********************************************************************
!*                                                                     *
!*        Compute the accumulation, source/sink and flow terms         *
!*                    SET UP THE JACOBIAN EQUATIONS                    *
!*                                                                     *
!***********************************************************************
!
!
            CALL JACOBIAN_SetUp
!
            IF(impossible_solution) THEN
!
               CALL Timestep_Reduction(KIT,DtCut_Flag,ConvergenceFailure)  ! Cut Dt if necessary
               IF(DtCut_Flag .EQV. .TRUE.) GO TO 2000                      ! If so, restart the Newtonian iteration loop
!
            END IF
!
!
!***********************************************************************
!*                                                                     *
!*                        UPON CONVERGENCE ...                         *
!*                                                                     *
!***********************************************************************
!
!
            IF_Convrg: IF(MaxResidual <= rel_convergence_crit) THEN
!
               CALL Upon_Convergence(NSTL,KIT,ConvergenceFailure,PrintNow_Flag)
!
               IF(PrintNow_Flag .EQV. .TRUE.) THEN  ! If the print-out flag is set, ...
                  CALL Print_Standard_Output        ! ...  produce a print-out and ...
                  CALL Mass_and_Energy_Balance(0)   ! ...  compute mass and heat balances
               END IF
!
               EXIT DO_NRIterations                 ! End NR iterations, and continue timestepping
!
            END IF IF_Convrg
!
!
!***********************************************************************
!*                                                                     *
!*               IF CONVERGENCE IS NOT YET ATTAINED ...                *
!*                                                                     *
!***********************************************************************
!
!
            IF(Option_Print_NRIterationInfo /= 0) THEN
!
! ......... Print basic information for Option_Print_NRIterationInfo /= 0
!
               WRITE(*,6006) Tot_NumTimeSteps, NumNRIterations, TimeStep, MaxResidual, elem%name(MaxResElemNum), MaxResEquNum
            END IF
!
! ----------------
! ......... Print extra information for Option_Print_NRIterationInfo > 2
! ----------------
!
            IF_Option_Print_NRIterationInfo: IF(Option_Print_NRIterationInfo >= 2) THEN
!
               IF(TrackElemNum /= 0) THEN
                  NIT = TrackElemNum
               ELSE
                  NIT = MaxResElemNum
               END IF
!
               NIT1 = (NIT-1)*NumComPlus1
!
               WRITE(*,6008) elem%name(nit), DX(NIT1+1), DX(NIT1+2), ElemState%temp(nit,current), &
     &                       ElemState%pres(nit,current), ElemProp(nit,0)%satur(1)
!
            END IF IF_Option_Print_NRIterationInfo
!
! ----------------
! ......... Check if the maximum number of Newtonian iterations has been reached
! ----------------
!
            IF(NumNRIterations > Max_NumNRIterations) THEN
!
               CALL Timestep_Reduction(KIT,DtCut_Flag,ConvergenceFailure)  ! Cut Dt if necessary
               IF(DtCut_Flag .EQV. .TRUE.) GO TO 2000                      ! If so, restart the Newtonian iteration loop
!
            END IF
!
!
!***********************************************************************
!*                                                                     *
!*                         SOLVE THE JACOBIAN                          *
!*                                                                     *
!***********************************************************************
!
!
            CALL CPU_Elapsed_Time(0,Dt_solv,T_dd)   ! Set timer
!
            CALL Solve_Jacobian_Matrix_Equation
!
            CALL CPU_Elapsed_Time(1,Dt_solv,T_dd)   ! Determine execution time
!
            CPU_MatrixSolTime = CPU_MatrixSolTime + Dt_solv
!
            IF(impossible_solution) THEN
!
               CALL Timestep_Reduction(KIT,DtCut_Flag,ConvergenceFailure)  ! Cut Dt if necessary
               IF(DtCut_Flag .EQV. .TRUE.) GO TO 2000                      ! If so, restart the Newtonian iteration loop
!
            END IF
!
! ----------------
! ......... Update the thermophysical properties
! ----------------
!
            CALL Equation_of_State
!
! ----------------
! ......... If unrealistic conditions are encountered, reduce time step
! ----------------
!
            IF(unrealistic_conditions) THEN
!
               CALL Timestep_Reduction(KIT,DtCut_Flag,ConvergenceFailure)  ! Cut Dt if necessary
               IF(DtCut_Flag .EQV. .TRUE.) GO TO 2000                      ! If so, restart the Newtonian iteration loop
!
            END IF
!
! ----------------
! ......... Computing changes to porosity, permeability and capillary pressure
! ----------------
!
            CALL Hydraulic_Property_Changes
!
! ----------------
! ......... Check if print-out is needed
! ----------------
!
            IF(Option_OutputAmount >= 10) CALL Print_Standard_Output
!
! <<<<<<
! <<<<<< End of the NEWTONIAN ITERATION LOOP
! <<<<<<
!
         END DO DO_NRIterations
!
!
!***********************************************************************
!*                                                                     *
!*          CHECK IF ANY OF THE VARIOUS TERMINATION CRITERIA           *
!*                         HAVE BEEN FULFILLED                         *
!*                                                                     *
!***********************************************************************
!
!
         CALL CPU_Timing_Routine(CPU_time)                                                   ! Check the CPU time up to this point
!
         IF(  (ConvergenceFailure) .OR.                                                   &  ! The termination flag is invoked
     &        (SimulationTimeEnd /= 0.0d0 .AND. time >= SimulationTimeEnd) .OR.           &  ! The simulation time is exceeded
     &        (Tot_NumTimeSteps >= Max_NumTimeSteps)                                      &  ! The number of time steps is exceeded
     &     ) THEN
!
            CALL Write_File_SAVE                                                             ! Write results into file SAVE
!
            EXIT                                                                             ! Done!  Exit the subroutine
!
         END IF
!
! <<<
! <<< End of the TIME STEP LOOP
! <<<
!
      END DO DO_TimeSteps
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 6000 FORMAT(/,'Simulation_Cycle',T50,'[v1.0,  18 August    2007]',/,   &
     &         ':::::   Executive routine that performs the simulation while marching in time')
!
!
 6001 FORMAT(' <NumTimeSteps> = ',I4,' <Tot_NumTimeSteps> =',I4, &
     &         '  <TimeStep> = ',1pE13.6,'  <time> = ',1pE13.6)
!
 6002 FORMAT(' Write coefficients for semi-analytical heat loss calculation onto file "TABLE"')
!
 6005 FORMAT(4E20.13)
 6006 FORMAT(' .... ITERATING:  At [',I6,',',I2,']: TimeStep = ',1pE13.6,  &
     &       '   MAX{Residual} = ',1pE13.6,'  at element "',A8,'", equation ',I2)
!
 6008 FORMAT(31X,'     AT "',A8,'" ...   DX1= ',1pE12.5,' DX2= ',1pE12.5,' T = ',F7.3,' P = ',F9.0,' S = ',1pE12.5)
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End of Simulation_Cycle
!
!
      RETURN
!
      END SUBROUTINE Simulation_Cycle
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE Timestep_Reduction(k_it,DtCut_Flag,ConvergenceFailure)
!
! ...... Modules to be used
!
         USE Basic_Parameters
         USE General_External_File_Units
         USE General_Control_Parameters
!
         USE Geological_Media_Properties
!
         USE Element_Attributes
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*               ROUTINE FOR REDUCING THE TIME-STEP SIZE               *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
      IMPLICIT NONE
!
! -------
! ... Integer variables
! -------
!
      INTEGER, INTENT(IN) :: k_it
!
! -------
! ... Logical variables
! -------
!
      LOGICAL, INTENT(OUT) :: ConvergenceFailure, DtCut_Flag
!
      LOGICAL  :: First_call = .TRUE.
!
! -------
! ... Saving variables
! -------
!
      SAVE First_call
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of Timestep_Reduction
!
!
      IF(First_call) THEN
         First_call = .FALSE.
         WRITE(VERS_Unit,6000)
      END IF
!
! -------
! ... Reduce time-step size
! -------
!
      NewTimeStep = TimeStep/Dt_reducer
!
      unrealistic_conditions = .FALSE.                   ! Reset the <unrealistic_conditions> flag
!
      DtCut_Flag = .FALSE.                               ! Initialize the flag for transfer of control on exit
!
!***********************************************************************
!*                                                                     *
!*  Check on number of previous time steps that converged on ITER = 1  *
!*                                                                     *
!***********************************************************************
!
      IF(k_it > 2 .AND. DoublingDtCriterion /= 0) THEN
         WRITE(*,6001) Tot_NumTimeSteps, TimeStep        ! Print explanatory comment
         ConvergenceFailure = .TRUE.                     ! Set flag to stop execution after following Dt
      END IF
!
      WRITE(*,6002) Tot_NumTimeSteps, NumNRIterations, NewTimeStep   ! Print explanatory comment
!
! -------
! ... Determine the number of Dt cutbacks
! -------
!
      NumTimeStepCuts = NumTimeStepCuts + 1
!
!***********************************************************************
!*                                                                     *
!*             STOP AFTER 25 ATTEMPTS TO REDUCE TIME STEP              *
!*                                                                     *
!***********************************************************************
!
      IF_IhafMax: IF(NumTimeStepCuts > 25) THEN
!
         WRITE(*,6003) TimeStep            ! Print warning message
         CALL Print_Standard_Output        ! Write output before stopping
         CALL Write_File_SAVE              ! Write the SAVE file before stopping
         STOP
!
      ELSE
!
         NumNRIterations = 0               ! Reset the counter of Newtonian iterations for new Dt
!
! ----------
! ... Update element pressures and temperatures
! ----------
!
         CALL Equation_of_State            ! Re-estimate thermophysical properties
!
         CALL Hydraulic_Property_Changes   ! Computing changes in hydraulic properties
!
         DtCut_Flag = .TRUE.               ! Set the flag for a new NI iteration with a shorter Dt on exit
!
      END IF IF_IhafMax
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 6000 FORMAT(/,'Timestep_Reduction',T50,'[v1.0,  18 December  2007]',/, &
     &         ':::::   Routine to reduce the timestep size')
!
 6001 FORMAT(/,' >>>>>>>>  CONVERGENCE FAILURE on time step #',I6,   &
     &         ' with Dt = ',1pE13.6,' seconds, following two steps that converged on <NumNRIterations> = 1',/,   &
     &         '      ==>  STOP EXECUTION AFTER NEXT TIME STEP')
 6002 FORMAT(' >>>>>  REDUCE TIME STEP AT (',I6,',',I2,') ==>  New timestep = ',1pE13.6)
!
 6003 FORMAT(//,' NO CONVERGENCE AFTER 25 REDUCTIONS OF TIME STEP SIZE',/,   &
     &          ' LAST TIME STEP WAS ',1pE13.6,' SECONDS ',//,   &
     &          ' !!!!!!!!!!  S T O P   E X E C U T I O N  !!!!!!!!!!')
 6004 FORMAT(' Write coefficients for semi-analytical heat loss calculation onto file "TABLE"')
!
 6005 FORMAT(4E20.13)
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>    End of Timestep_Reduction
!
!
      RETURN
!
      END SUBROUTINE Timestep_Reduction
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE Upon_Convergence(NSTL,KIT,ConvergenceFailure,PrintNow_Flag)
!
! ...... Modules to be used
!
         USE Basic_Parameters
         USE EOS_Default_Parameters
         USE General_External_File_Units
!
         USE General_Control_Parameters
!
         USE Grid_Geometry, ONLY: elem
         USE Element_Attributes
         USE Sources_and_Sinks
!
         USE Solution_Matrix_Arrays
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*        ROUTINE FOR UPDATING PRIMARY VARIABLES AND SELECTING         *
!*          AFTER SUCCESSFUL CONVERGENCE OF THE NI ITERATIONS          *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
      IMPLICIT NONE
!
! -------
! ... Integer variables
! -------
!
      INTEGER, INTENT(IN)  :: NSTL
      INTEGER, INTENT(OUT) :: KIT
!
      INTEGER :: N_up, n
!
! -------
! ... Logical variables
! -------
!
      LOGICAL, INTENT(INOUT) :: ConvergenceFailure
      LOGICAL, INTENT(OUT)   :: PrintNow_Flag
!
      LOGICAL :: First_call = .TRUE.

#ifdef USE_TIMER
         real(KIND = 8) :: start, finish
#endif
!
! -------
! ... Saving variables
! -------
!
      SAVE First_call, N_up
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of Upon_Convergence
!
      IF(First_call) THEN
         First_call = .FALSE.
         N_up       = (NumElem-1)*NumComPlus1 + 1
         WRITE(VERS_Unit,6000)
      END IF
!
! -------
! ... Reset the convergence flag to indicate convergence
! -------
!
      convergence = .TRUE.
!
!
!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>                                                                     >
!>                       ACTIVE ELEMENT LOOP                           >
!>       Compute Changes In Properties, Parameters And Variables       >
!>                        After Convergence                            >
!>                                                                     >
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
! ----------
! ... Compute changes in porosity: Array operations
! ----------
!
      ElemMedia(1:NumElem,previous)%porosity = ElemMedia(1:NumElem,current)%porosity
!
! ----------
! ... Update element pressures and temperatures
! ----------
!
      DO n = 1,NumElem
         ElemState%pres(n, previous) = ElemState%pres(n, current)
         ElemState%temp(n, previous) = ElemState%temp(n, current)
      END DO
!
! ----------
! ... Update primary variables: CAREFUL! Array operations!
! ----------
!
! This is the solution
!
 1000 X = X + DX                   !!! <<<<<< Careful !
!
!
!***********************************************************************
!*                                                                     *
!*             Store the state index in case of Dt cutback             *
!*                                                                     *
!***********************************************************************
!
      ElemState%index(n, previous) = ElemState%index(n, current)   ! CAREFUL! Whole array operation
!
! -------
! ... Udpate current time
! -------
!
      time = time + TimeStep                                 ! The new current time
!
!
!***********************************************************************
!*                                                                     *
!*             SET NEW TIMESTEP; UPDATE OTHER PARAMETERS               *
!*                                                                     *
!***********************************************************************
!
!
      IF(SimulationTimeEnd /= 0.0d0 .AND. time == SimulationTimeEnd) THEN
         PrintOutNow = .TRUE.
      END IF
!
! ----------
! ... For internally-determined time step sizes
! ----------
!
      IF((NumUserTimeSteps == 0) .OR. (NumTimeSteps+1 > 8*NumUserTimeSteps .OR. UserTimeStep(NumTimeSteps+1) == 0.0d0)) THEN
!
         IF_DtDouble: IF(NumNRIterations <= DoublingDtCriterion) THEN     ! If convergence within <DoublingDtCriterion> iterations, ...
            NewTimeStep = 2.0d0*TimeStep                                  !    ... double time step, ...
         ELSE                                                             !    ... otherwise, ...
            NewTimeStep = TimeStep                                        !    ... use old time step
         END IF IF_DtDouble
!
! ----------
! ... For internally-determined time step sizes after exhaustion of user-supplied Dt list
! ----------
!
      ELSE IF(NumTimeSteps+1 < NumUserTimeSteps) THEN
!
         IF(UserTimeStep(NumTimeSteps + 1) == 0.0d0) THEN
         IF_DtDouble2: IF(NumNRIterations <= DoublingDtCriterion) THEN    ! If convergence within <DoublingDtCriterion> iterations, ...
            NewTimeStep = 2.0d0*TimeStep                                  !    ... double time step, ...
         ELSE                                                             !    ... otherwise, ...
            NewTimeStep = TimeStep                                        !    ... use old time step
         END IF IF_DtDouble2
!
      END IF
!
! ----------
! ... For given (predetermined) time step sizes
! ----------
!
      ELSE
         NewTimeStep = UserTimeStep(NumTimeSteps+1)                       ! Use specified time step
      END IF
!
! ----------
! ... Adjust Dt to ensure conformance to the maximum simulation period <SimulationTimeEnd>
! ----------
!
      IF(SimulationTimeEnd > 0.0d0) NewTimeStep = MIN( NewTimeStep, SimulationTimeEnd-time )
!
! ----------
! ... Adjust Dt to ensure conformance with max{Dt} = MaxTimeStep
! ----------
!
      IF(MaxTimeStep > 0.0d0) NewTimeStep = MIN(NewTimeStep,MaxTimeStep)
!
!***********************************************************************
!*                                                                     *
!*            Determine flags and counters for print-outs              *
!*                                                                     *
!***********************************************************************
!
      IF(NumNRIterations == 1) THEN    ! Count number of consecutive time steps that converge on num_NR_iter = 1
         KIT = KIT+1
      ELSE
         KIT = 0                       ! For NumNRIterations > 1, do not consider
      END IF
!
! -------
! ... Print messages if KIT >= 10
! -------
!
      IF(KIT >= 10) THEN
         WRITE(*,6001)
         ConvergenceFailure = .TRUE.   ! Set the flag indicating convergence failure
      END IF
!
!***********************************************************************
!*                                                                     *
!*             Print out important data upon convergence               *
!*                                                                     *
!***********************************************************************
!
      IF_Conver: IF(convergence .EQV. .TRUE.) THEN
!
         IF_NST: IF(TrackElemNum /= 0) THEN
!
            WRITE(*,6002) elem%name(TrackElemNum), Tot_NumTimeSteps, NumNRIterations, time,    &
     &                    TimeStep, DX(NSTL+1), DX(NSTL+2),                                    &
     &                    ElemState%temp(TrackElemNum,0), ElemState%pres(TrackElemNum,0),    &
     &                    ElemProp(TrackElemNum,0)%satur(GasPhase)
!
         ELSE
!
            WRITE(*,6002) elem%name(MaxResElemNum), Tot_NumTimeSteps, NumNRIterations, time,                      &
     &                    TimeStep,DX( (MaxResElemNum-1)*NumComPlus1+1 ),DX( (MaxResElemNum-1)*NumComPlus1+2 ),   &
     &                    ElemState%temp(MaxResElemNum,0), ElemState%pres(MaxResElemNum,0),                     &
     &                    ElemProp(MaxResElemNum,0)%satur(GasPhase)
!
         END IF IF_NST
!
      END IF IF_Conver
!
!***********************************************************************
!*                                                                     *
!*              Determine whether printout is required                 *
!*                                                                     *
!***********************************************************************
!
      IF( (ConvergenceFailure) .OR.                                                  &
     &    ( (PrintOutNow) .AND. (convergence) ) .OR.                                 &
     &    ( (convergence) .AND. (MOD(Tot_NumTimeSteps,PRINT_frequency) == 0) ) .OR.  &
     &    (Tot_NumTimeSteps >= Max_NumTimeSteps) .OR.                                &
     &    (Option_OutputAmount >= 10)                                                &
     &  ) PrintNow_Flag = .TRUE.
!
!***********************************************************************
!*                                                                     *
!*             Determine whether to update the SAVE file,              *
!*           as well as important computational parameters             *
!*                                                                     *
!***********************************************************************
!
      IF(SAVE_frequency > 0) THEN
!
         IF_CheckPrint: IF(MOD(Tot_NumTimeSteps, SAVE_frequency) == 0) THEN
!
! ......... First, print an interim SAVE file
!
            CALL Write_File_SAVE
!
! -------------
! ......... Modify simulation control parameters without stopping execution
! -------------
!
            CALL Update_Simulation_Parameters
!
            IF( (Tot_NumTimeSteps >= Max_NumTimeSteps) .OR. (time >= SimulationTimeEnd) ) PrintNow_Flag = .TRUE.
!
         END IF IF_CheckPrint
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
 6000 FORMAT(/,'Upon_Convergence',T50,'[v1.0,  29 September 2006]',/,   &
     &         ':::::   Update primary variables and select timestep after convergence of the NI iterations (generic)')
!
 6001 FORMAT(' !!!!!!!!!  FOR 10 CONSECUTIVE TIME STEPS HAVE CONVERGENCE ON <NumNRIterations> = 1',/,   &
     &       '       ==>  WRITE OUT CURRENT DATA, THEN STOP EXECUTION')
!
 6002 FORMAT(1X,A8,'(',I6,',',I2,')',' ST =',1pE12.5,' DT =',1pE12.5,' DX1= ',1pE12.5,' DX2= ',1pE12.5,' T = ',1pE11.4, &
     &             ' P = ',1pE13.6,' S = ',1pE12.5)
!
 6005 FORMAT(/,' ! >>>>>>>>>>  ',/, &
     &         ' ! >>>>>>>>>>  ',/, &
     &         ' ! ..........  The pressure in cell <',a8,'> exceeds the maximum level (=',1pE13.6,  &
     &         ' Pa) defined for source/sink <',a5,'>',/, &
     &         ' ! ..........  The rate factor is reset to zero, and the simulation continues (SS_response = <ZERO>)',/, &
     &         ' ! >>>>>>>>>>  ',/, &
     &         ' ! >>>>>>>>>>  ',/)
!
 6006 FORMAT(/,' ! >>>>>>>>>>  ',/, &
     &         ' ! >>>>>>>>>>  ',/, &
     &         ' ! ..........  The pressure in cell <',a8,'> is below the minimum level (=',1pE13.6,  &
     &         ' Pa) defined for source/sink <',a5,'>',/, &
     &         ' ! ..........  The rate factor is reset to zero, and the simulation continues (SS_response = <ZERO>)',/, &
     &         ' ! >>>>>>>>>>  ',/, &
     &         ' ! >>>>>>>>>>  ',/)
!
 6015 FORMAT(/,' ! >>>>>>>>>>  ',/, &
     &         ' ! >>>>>>>>>>  ',/, &
     &         ' ! ..........  The pressure in cell <',a8,'> exceeds the maximum level (=',1pE13.6,  &
     &         ' Pa) defined for source/sink <',a5,'>',/, &
     &         ' ! ..........  The rate factor is adjusted from',1pE13.6,' to',1pE13.6, &
     &         ', and the simulation continues (SS_response = <ADJU>)',/, &
     &         ' ! >>>>>>>>>>  ',/, &
     &         ' ! >>>>>>>>>>  ',/)
!
 6016 FORMAT(/,' ! >>>>>>>>>>  ',/, &
     &         ' ! >>>>>>>>>>  ',/, &
     &         ' ! ..........  The pressure in cell <',a8,'> is below the minimum level (=',1pE13.6,  &
     &         ' Pa) defined for source/sink <',a5,'>',/, &
     &         ' ! ..........  The rate factor is adjusted from',1pE13.6,' to',1pE13.6, &
     &         ', and the simulation continues (SS_response = <ADJU>)',/, &
     &         ' ! >>>>>>>>>>  ',/, &
     &         ' ! >>>>>>>>>>  ',/)
!
 6025 FORMAT(/,' ! >>>>>>>>>>  ',/, &
     &         ' ! >>>>>>>>>>  ',/, &
     &         ' ! ..........  The pressure in cell <',a8,'> exceeds the maximum level (=',1pE13.6,  &
     &         ' Pa) defined for source/sink <',a5,'>',/, &
     &         ' ! ..........  The results at the end of this time step are printed,', &
     &         ' and the simulation is halted (SS_response = <STOP>)',/, &
     &         ' ! >>>>>>>>>>  ',/, &
     &         ' ! >>>>>>>>>>  ',/)
!
 6026 FORMAT(/,' ! >>>>>>>>>>  ',/, &
     &         ' ! >>>>>>>>>>  ',/, &
     &         ' ! ..........  The pressure in cell <',a8,'> is below the minimum level (=',1pE13.6,  &
     &         ' Pa) defined for source/sink <',a5,'>',/, &
     &         ' ! ..........  The results at the end of this time step are printed,', &
     &         ' and the simulation is halted (SS_response = <STOP>)',/, &
     &         ' ! >>>>>>>>>>  ',/, &
     &         ' ! >>>>>>>>>>  ',/)
!
 6200 FORMAT(8(1pe15.8,1x))
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End of Upon_Convergence
!
!
      RETURN
!
      END SUBROUTINE Upon_Convergence
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE Update_Simulation_Parameters
!
! ...... Modules to be used
!
         USE Basic_Parameters
         USE General_External_File_Units
         USE General_Control_Parameters
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*       ROUTINE TO ALLOW MODIFICATION OF SIMULATION PARAMETERS        *
!*                 WITHOUT INTERRUPTION OF EXECUTION                   *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
      IMPLICIT NONE
!
! -------
! ... Integer variables
! -------
!
      INTEGER :: n, time_step_selection
!
! -------
! ... Character variables
! -------
!
      CHARACTER (LEN = 28) :: UpdateHeader
      CHARACTER (LEN = 28) :: FirstLineUnit, FirstLine
!
! -------
! ... Logical variables
! -------
!
      LOGICAL :: First_call = .TRUE., exists
!
! -------
! ... Namelists
! -------
!
      NAMELIST/Real_Parameters_To_Update/ SimulationTimeEnd, MaxTimeStep, rel_convergence_crit, abs_convergence_crit
      NAMELIST/Integer_Parameters_To_Update/ Max_NumTimeSteps, Max_NumNRIterations, time_step_selection, SAVE_frequency
!
! -------
! ... Saving variables
! -------
!
      SAVE First_call
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of Update_Simulation_Parameters
!
!
      IF(First_call) THEN
!
         First_call = .FALSE.
         NumParamUpdates = 0
!
         WRITE(VERS_Unit,6000)
!
      END IF
!
! -------
! ... Determine the status of the file <Parameter_Update_File>
! -------
!
      INQUIRE(FILE='Parameter_Update_File', EXIST = exists)
!
      IF(exists) THEN                                                   ! If it already exists, ...
         OPEN(Parameter_Update_Unit, FILE='Parameter_Update_File')      !    open it as an old file and ...
      ELSE                                                              ! Otherwise, ...
         RETURN                                                         !    exit the routine
      END IF
!
! -------
! ... Check if important simulation parameters are to be updated
! -------
!
      READ(UNIT = Parameter_Update_Unit, FMT = *, IOSTAT = n) UpdateHeader
!
      IF(n < 0) THEN
         CLOSE (UNIT = Parameter_Update_Unit)
         RETURN
      END IF
!
! ...... Check heading first
!
      IF_Update: IF(UpdateHeader == 'Update_Simulation_Parameters' .OR. UpdateHeader == 'update_simulation_parameters' .OR. &
     &              UpdateHeader == 'UPDATE_SIMULATION_PARAMETERS') THEN
!
! ...... Read and update real parameters
!
         READ(UNIT = Parameter_Update_Unit, NML = Real_Parameters_To_Update, IOSTAT = n)
!
         IF(n < 0) THEN
            CLOSE (UNIT = Parameter_Update_Unit)
            RETURN
         END IF
!
! ...... Read and update integer parameters
!
         READ(UNIT = Parameter_Update_Unit, NML = Integer_Parameters_To_Update, IOSTAT = n)
!
         CLOSE (UNIT = Parameter_Update_Unit)
!
         IF(n < 0) RETURN
!
! ...... Update number of updates, some adjustments, printouts
!
         NumParamUpdates = NumParamUpdates + 1
!
         WRITE(*, FMT = 6010) NumParamUpdates,                                                            &
     &                        SimulationTimeEnd, MaxTimeStep, rel_convergence_crit, abs_convergence_crit, &
     &                        Max_NumTimeSteps, Max_NumNRIterations, DoublingDtCriterion, SAVE_frequency
!
! ...... Modify the update file as a DIRECT access file (allows over-writing info in the file)
!
         OPEN(Parameter_Update_Unit, FILE='Parameter_Update_File', ACCESS = 'DIRECT', RECL = 28)
!
         WRITE(UNIT = FirstLineUnit, FMT = 6005)  NumParamUpdates  !    Use internal file to store the modified 1st line as Character + Integer ...
!
         READ (UNIT = FirstLineUnit, FMT = 6006)  FirstLine        !    Retrieve the 1st line as a character variable ...
!
! ...... Modify the top line of the update file to avoid accidentally affecting future or continuation simulations
!
         WRITE (UNIT = Parameter_Update_Unit, REC = 1) FirstLine
!
      END IF IF_Update
!
      CLOSE (UNIT = Parameter_Update_Unit)
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 6000 FORMAT(/,'Update_Simulation_Parameters',T50,'[v1.0,  15 February  2007]',/,   &
     &         ':::::   Update simulation contol parameters without interrupting execution')
!
 6005 FORMAT('==> Completed Update #',i3,'   ')
 6006 FORMAT(a28)
!
 6010 FORMAT(//,' ! >>>>>>>>>>  ',/, &
     &          ' ! >>>>>>>>>>  ',/, &
     &          ' ! ..........  UPDATING SIMULATION CONTROL PARAMETERS - Update #',i3, /, &
     &          ' ! ..........  ',/,   &
     &          ' ! ..........       Length of simulation period <SimulationTimeEnd>               =',1pe12.5,/,   &
     &          ' ! ..........       Maximum allowable timstep size <MaxTimeStep>                  =',1pe12.5,/,   &
     &          ' ! ..........       Relative convergence criterion <rel_convergence_crit>         =',1pe12.5,/,   &
     &          ' ! ..........       Absolute convergence criterion <rel_convergence_crit>         =',1pe12.5,/,   &
     &          ' ! ..........  ',/,   &
     &          ' ! ..........       Maximum number of timesteps <Max_NumTimeSteps>                     =',i7,/,   &
     &          ' ! ..........       Maximum number of NR iterations per timestep <Max_NumNRIterations> =',i7,/,   &
     &          ' ! ..........       NR iteration-to-convergence criterion for doubling Dt              =',i7,/,   &
     &          ' ! ..........       Printing frequency of SAVE file <SAVE_frequency>                   =',i7,/,   &
     &          ' ! ..........  ',/,   &
     &          ' ! >>>>>>>>>>  ',/, &
     &          ' ! >>>>>>>>>>  ',//)
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Update_Simulation_Parameters
!
!
      RETURN
!
      END SUBROUTINE Update_Simulation_Parameters
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE Adjust_Timestep_for_Printout
!
! ...... Modules to be used
!
         USE Basic_Parameters
         USE General_External_File_Units
         USE General_Control_Parameters
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*        ROUTINE MODIFYING THE TIME-STEP SIZE TO CONFORM WITH         *
!*     USER-SUPPLIED INPUT TIMES AT WHICH PRINT-OUTS ARE DESIRED       *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
      IMPLICIT NONE
!
! -------
! ... Integer variables
! -------
!
      INTEGER :: i
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
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of Adjust_Timestep_for_Printout
!
!
      IF(First_call) THEN
         First_call = .FALSE.
         WRITE(VERS_Unit,6000)
      END IF
!
! -------
! ... Check if the current time exceeds the maximum print-out time
! -------
!
      IF(PrintTime(NumPrintTimes) <= time)  RETURN                     ! Done ! Exit the routine
!
!
!***********************************************************************
!*                                                                     *
!*                        PRINT-OUT TIME LOOP                          *
!*                                                                     *
!***********************************************************************
!
!
      DO_NPTimes: DO  i=1,NumPrintTimes
!
! ----------
! ...... If the desired print-out time is reached
! ----------
!
         IF_ChkTIS: IF(PrintTime(i) == time) THEN
!
            IF(TimeStepMax_PostPrint > 0.0d0) NewTimeStep = MIN(NewTimeStep,TimeStepMax_PostPrint)  ! Obtain an estimate for the new Dt
!
! ......... Check the Dt estimate against the next print-out time
!
            IF(PrintTime(I+1) > PrintTime(i)+NewTimeStep) THEN              ! If Dt < PrintTime(I+1)-PrintTime(I)
               RETURN                                                       !    Done ! Exit the routine
            ELSE                                                            ! Otherwise
               NewTimeStep = MIN(NewTimeStep,PrintTime(i+1)-PrintTime(i))   !    Readjust Dt
            END IF
!
            PrintOutNow = .TRUE.                                            ! Reset the print-out flag
!
            RETURN                                                          ! Done ! Exit the routine
!
! ----------
! ...... If the desired print-out time has already been passsed
! ----------
!
         ELSE IF(PrintTime(i) < time) THEN
!
            CYCLE                                                          ! Continue the loop
!
! ----------
! ...... If the desired print-out time has not been passsed
! ----------
!
         ELSE
!
            IF(time + NewTimeStep >= PrintTime(i)) THEN                    ! If Dt results in current t > PrintTime(I), ...
               NewTimeStep   = PrintTime(i)-time                           ! ... Dt is adjusted and
               PrintOutNow = .TRUE.                                        ! ... the print-out flag is reset
            END IF
!
            RETURN                                                         ! Done ! Exit the routine
!
         END IF IF_ChkTIS
!
! <<<
! <<< End of the PRINT-OUT TIME LOOP
! <<<
!
      END DO DO_NPTimes
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 6000 FORMAT(/,'Adjust_Timestep_for_Printout',T50,'[v1.0,   4 November  2007]',/,   &
     &         ':::::   Adjust time steps to conform with user-defined times for a print-out')
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Adjust_Timestep_for_Printout
!
!
      RETURN
!
      END SUBROUTINE Adjust_Timestep_for_Printout
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE Solve_Jacobian_Matrix_Equation
!
! ...... Modules to be used
!
         USE Basic_Parameters
         USE General_External_File_Units
!
         USE General_Control_Parameters
         USE Solver_Parameters
!
         USE Grid_Geometry
         USE Element_Attributes
         USE Solution_Matrix_Arrays
!
         USE Matrix_Solvers
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*     ROUTINE CALLING THE SOLVERS FOR THE JACOBIAN MATRIX EQUATION    *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
      IMPLICIT NONE
!
! -------
! ... double-precision work array (automatic creation & destruction)
! -------
!
      REAL(KIND = 8), DIMENSION(SIZE(R)) :: work
      REAL(KIND = 8) :: ERR
!
! -------
! ... Integer variables
! -------
!
      INTEGER :: ICALL = 0, iteruc = 0, iprpro =0, izero0 = 0
!
      INTEGER :: N, nn, k, N_k, N_e, MainDiagonal, ITERU, IERR
      INTEGER :: izerod, ierG, info
      INTEGER :: iddfup, iddfdn, matord, nsupdg, nsubdg, ntotd
!
      SAVE ICALL, N, iteruc, iprpro, izero0
      SAVE matord, nsupdg, nsubdg, ntotd, N_k, N_e, MainDiagonal
!
!
#ifdef USE_TIMER
         real(KIND = 8) :: start, finish
#endif
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of Solve_Jacobian_Matrix_Equation
!
!
      ICALL = ICALL+1
!
      IF(ICALL == 1) THEN
!
         WRITE(VERS_Unit,6000)
!
         N   =  NumEqu*NumElem
         N_k = (NumElem-1)*NumComPlus1
         N_e = (NumElem-1)*NumEqu
!
         MainDiagonal = NumEquSquare*NumElem
!
      END IF
!
!
!***********************************************************************
!*                                                                     *
!*              ACCOUNTING FOR ZEROs ON THE MAIN DIAGONAL              *
!*                                                                     *
!***********************************************************************
!
!
!
! ... Initializations
!
      izerod = 0
      iprpro = 0
!
! ... DIRECT SOLVERS: Print solver data for Option_Print_SolverInfo > 0
!
      IF(Option_Print_SolverInfo /= 0) WRITE(*,6006) Tot_NumTimeSteps,NumNRIterations,icall,N
!
! -------
! ... Print matrix when Option_Print_SolverInfo >= 7
! ... Whole array operations - CAREFUL!!!
! -------
!
      IF(Option_Print_SolverInfo >= 7) THEN
         WRITE(*,6015)
         WRITE(*,6010) (RowNum(NN),ColNum(NN),CO(NN),NN=1,Non0)
      END IF
!
!
!***********************************************************************
!*                                                                     *
!*       DETERMINATION OF THE BANDWIDTHS OF THE U AND L MATRICES       *
!*                IN THE DIRECT LU DECOMPOSITION SOLVER                *
!*                                                                     *
!***********************************************************************
!
      IF(icall == 1) THEN
         CO(Max_NumNonZeros+1) = 0.0d0
      END IF
!
! -------
! ... Determine bandwidths and memory requirements for direct solver
! -------
!
      IF_Direct: IF(MatrixSolver == 1 .AND. icall == 1) THEN
!
         iddfup = MAXVAL( ColNum(1:Non0) - RowNum(1:Non0) )
         iddfdn = MAXVAL( RowNum(1:Non0) - ColNum(1:Non0) )
!
         matord = NumEqu*NumElem
         nsupdg = iddfup
         nsubdg = iddfdn
         ntotd  = 2*nsubdg+nsupdg+1
!
! ...... Deallocate memory for iterative solver, then re-allocate it with far less memory
!
         DEALLOCATE(rwork, iwork, STAT = ierG)
!
! ...... Allocate memory for direct solver
!
         ALLOCATE(rwork(MAX(NumComPlus1*NumPhases*Max_NumElem,19*Max_NumElem+Max_NumConx*NumPerturb)), &
     &            iwork(2*MatrixOrder), STAT = ierG)
!
! ...... Check if properly deallocated
!
         IF(ierG == 0) THEN
            WRITE(MemoryAlloc_Unit,6101) 102   ! Succesful deallocation
         ELSE
            WRITE(*,6102) 102                  ! Unsuccesful deallocation
            WRITE(MemoryAlloc_Unit,6102) 102
            STOP
         END IF
!
! ...... Allocate memory for direct solver
!
         ALLOCATE(AB(ntotd,matord), STAT = ierG)
!
! ...... Check if properly allocated
!
         IF(ierG == 0) THEN
            WRITE(MemoryAlloc_Unit,6103) 104   ! Succesful deallocation
         ELSE
            WRITE(*,6104)  104                 ! Unsuccesful deallocation
            WRITE(MemoryAlloc_Unit,6104) 104
            STOP
         END IF
!
      END IF IF_Direct
!
!
!***********************************************************************
!*                                                                     *
!*                          MATRIX SOLUTION                            *
!*                                                                     *
!***********************************************************************
!
      IF_MatrixSolution: IF(MatrixSolver == 1) THEN
!
! ----------
! ...... Direct solver
! ----------
!
         CALL LUD_Direct_Solver(matord,Non0,nsubdg,nsupdg,ntotd,matord,info)
!
      ELSE
!
         work(1:N) = 0.0d0 ! whole array operation
!
!
         Max_NumCGIterations = MAX(20, INT(N * Max_CGIterationRatio))
!
         IF(MatrixSolver == 2) THEN
            CALL DSLUBC(N, R, work, Non0, RowNum, ColNum, CO, 0, 2,                         &
          &             CG_convergence_crit, Max_NumCGIterations, ITERU, ERR, IERR,         &
          &             LINEQ_Unit, RWORK, real_work_array_size, IWORK, intg_work_array_size )
!
         ELSE IF (MatrixSolver == 3) THEN
!
            CALL DSLUCS(N, R, work, Non0, RowNum, ColNum, CO, 0, 2,                         &
          &             CG_convergence_crit, Max_NumCGIterations, ITERU, ERR, IERR,         &
          &             LINEQ_UNIT, RWORK, real_work_array_size, IWORK, intg_work_array_size )
!
         ELSE IF (MatrixSolver == 4) THEN
!
            CALL LIS(N, R, work, Non0, RowNum, ColNum, CO)
!
         END IF
!
! ...... Retrieve information from work
         R(1:N) = work(1:N)
!
      END IF IF_MatrixSolution
!
! -------
! ... Print-out for Option_Print_SolverInfo >= 5
! -------
!
      IF_Option_Print_SolverInfo: IF(Option_Print_SolverInfo >= 5) THEN
         WRITE(*,6050)
!
         DO_NumEle: DO nn=1,NumElem
            WRITE(*,6055) elem%name(nn),(R((nn-1)*NumEqu+k),k=1,NumEqu)
         END DO DO_NumEle
!
      END IF IF_Option_Print_SolverInfo
!
!
!*********************************************************************
!*                                                                   *
!*            UPDATE CHANGES - Whole array operations                *
!*                                                                   *
!*********************************************************************
!
!
#ifdef USE_TIMER
         call CPU_Timing_Routine(start)
#endif

#ifdef USE_OMP
!$OMP PARALLEL
!$OMP DO schedule(auto)
#endif
      DO k=1,NumEqu
         DX(k:N_k+k:NumComPlus1) = DX(k:N_k+k:NumComPlus1) + W_NRIteration * R(k:N_e+k:NumEqu)
      END DO
#ifdef USE_OMP
!$OMP END DO
!$OMP END PARALLEL
#endif

#ifdef USE_TIMER
         call CPU_Timing_Routine(finish)

         write (*,*) __FILE__, ":", __LINE__, " time: ", finish-start
#endif
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 6000 FORMAT(/,'Solve_Jacobian_Matrix_Equation',T50,'[v1.0   24 August    2007]',/,   &
     &         ':::::   Solves the Jacobian matrix equation through an interface with linear equation solvers ',/,   &
     &         '        Can call a direct solver or a package of conjugate gradient solvers')
!
 6003 FORMAT(1X,'At <Tot_NumTimeSteps> =',i6,' and <NumNRIterations> =',i7,', IZEROD=',I6,    &
     &              ', <Z_preprocessing> = ',a2,' and <O_preprocessing> = ',a2)
 6004 FORMAT(//,' ',15('WARNING-'),//T14,F5.2,'% of the elements on ',   &
     &          'the main diagonal of the Jacobian matrix are zeros',/,   &
     &       T14,'The matrix preprocessor <O_preprocessing> = ',A2,' cannot be used',/,   &
     &       T14,'Action taken: reset <O_preprocessing> to "O0" (no O-preprocessing): Continue execution')
 6005 FORMAT(' <<< Solve_Jacobian_Matrix_Equation >>>  at [Tot_NumTimeSteps, NR_iter] = [',I6,',',I3,   &
     &       ']','   ICALL =',I5,'   n =',I7,'   Non0 =',I8,'   IZEROD =',I7,'   or ',F6.2,' % zeros')
 6006 FORMAT(' <<< Solve_Jacobian_Matrix_Equation >>>  at [Tot_NumTimeSteps, NR_iter] = [',I6,',',I3,']',  &
     &       '   ICALL =',I5,'   N =',I7,'   (Direct solver)')
!
 6010 FORMAT(5(1X,I5,1x,I5,1x,1pE14.7))
 6015 FORMAT(/,' MATRIX OF COEFFICIENTS',/)
!
 6020 FORMAT(//,20('ERROR-'),//,T33,   &
     &             '       S I M U L A T I O N   A B O R T E D',   &
     &       /, T40,'DECLARED BANDWIDTH SMALLER THAN NEEDED',      &
     &       /, T33,'        PLEASE CORRECT AND TRY AGAIN',        &
     &       //,T20,'THE NUMBER OF NEEDED DIAGONALS IS', I4,' WHILE THE AVAILABLE NUMBER IS ', I4,   &
     &       //,T20,'THE PARAMETER real_work_array_size MUST BE INCREASED TO AT LEAST ',I10,   &
     &       //,20('ERROR-'))
 6030 FORMAT('  EPS = ',1pE11.4,  '  RMIN = ',1pE11.4,   &
     &       '  RESID = ',1pE11.4,'  IRNCP = ',I6,       &
     &       '  MINIRN = ',I6,  '  MINICN = ',I6,'IRANK = ',I4)
!
 6040 FORMAT('      SOLUTION TIME = ',1PE12.5,'  SECONDS')
 6045 FORMAT(T5,' At [',I6,',',I3,']',' <NewTimeStep> =',1pE12.5,' IERR=',I1,'& ERR=',1PE12.6,' IT=',I5,' ITC=',I10)
!
 6050 FORMAT(/,' ===== INCREMENTS ====   in order of primary variables',/)
 6055 FORMAT('       AT ELEMENT "',A5,'"   ',8(1X,1pE15.8))
!
 6101 FORMAT(/,T3,'Memory allocation at point ',i3,' in <Solve_Jacobian_Matrix_Equation> was successful')
 6102 FORMAT(//,20('ERROR-'),//,   &
     &       T3,'Memory allocation at point ',i3,' in <Solve_Jacobian_Matrix_Equation> was unsuccessful',//,   &
     &       T2,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6103 FORMAT(T3,'Memory allocation at point ',i3,' in <Solve_Jacobian_Matrix_Equation> was successful')
 6104 FORMAT(//,20('ERROR-'),//,   &
     &       T3,'Memory allocation at point ',i3,' in <Solve_Jacobian_Matrix_Equation> was unsuccessful',//,   &
     &       T2,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End of Solve_Jacobian_Matrix_Equation
!
!
      RETURN
!
      END SUBROUTINE Solve_Jacobian_Matrix_Equation
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE LUD_Direct_Solver(N,NZx,KL,KU,LDAB,LDB,INFO)
!
! ...... Modules to be used
!
         USE Basic_Parameters
         USE General_External_File_Units
!
         USE General_Control_Parameters
         USE Solver_Parameters
!
         USE Solution_Matrix_Arrays
!
         USE Matrix_Solvers
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*         ROUTINE FOR POPULATING THE 2-D ARRAY FOR THE DIRECT         *
!*             BANDED MATRIX SOLVER USING LU DECOMPOSITION             *
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
      REAL(KIND = 8) :: TS,TT,TD,TSS
!
! -------
! ... Integer variables
! -------
!
      INTEGER, INTENT(IN)  :: N,NZx,KL,KU,LDAB,LDB
      INTEGER, INTENT(OUT) :: INFO
!
      INTEGER :: i
!
! -------
! ... Logical variables
! -------
!
      LOGICAL :: First_call = .TRUE.
!
#ifdef USE_TIMER
      real(KIND = 8) :: start, finish
#endif
! -------
! ... Saving variables
! -------
!
      SAVE First_call
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of LUD_Direct_Solver
!
!
      IF(First_call) THEN
         First_call = .FALSE.
         WRITE(VERS_Unit,6000)
      END IF
!
!
!***********************************************************************
!*                                                                     *
!*                    I N I T I A L I Z A T I O N                      *
!*                                                                     *
!***********************************************************************
!
      info           = 0
#ifdef USE_TIMER
            call CPU_Timing_Routine(start)
#endif

#ifdef USE_OMP
!$OMP SIMD
#endif
      DO i=1,N
         iwork(i)     = 0                                    !  CAREFUL! Whole array operation
         AB(1:LDAB,i) = 0.0d0                                !  CAREFUL! Whole array operation
      END DO

#ifdef USE_TIMER
            call CPU_Timing_Routine(finish)

            write (*,*) __FILE__, ":", __LINE__, " time: ", finish-start
#endif
!
!***********************************************************************
!*                                                                     *
!*            PLACEMENT OF MATRIX ELEMENTS INTO THE AB ARRAY           *
!*                                                                     *
!***********************************************************************
!
      FORALL (i=1:NZx)
         AB(KL+KU+1+RowNum(i)-ColNum(i),ColNum(i)) = co(i)  ! CAREFUL! Whole array operation
      END FORALL
!
      IF(Option_Print_SolverInfo /= 0) CALL CPU_Elapsed_Time(0,TS,TT)
!
!***********************************************************************
!*                                                                     *
!*                          LU DECOMPOSITION                           *
!*                                                                     *
!***********************************************************************
!
      CALL DGBTRF(N,N,KL,KU,AB,LDAB,iwork,INFO)
!
      IF(Option_Print_SolverInfo /= 0) THEN
         CALL CPU_Elapsed_Time(1,TD,TT)
         WRITE(*,6001) TD,1
      END IF
!
!***********************************************************************
!*                                                                     *
!*       SOLUTION USING THE LU FACTORIZATION COMPUTED BY DGBTRF        *
!*                                                                     *
!***********************************************************************
!
      IF(info == 0) THEN
!
         IF(Option_Print_SolverInfo /= 0) CALL CPU_Elapsed_Time(0,TS,TT)
!
         CALL DGBTRS(N,KL,KU,AB,LDAB,iwork,R,LDB)
!
         IF(Option_Print_SolverInfo /= 0) THEN
            CALL CPU_Elapsed_Time(1,TSS,TT)
            WRITE(*,6002) TSS
         END IF
!
      ELSE
         WRITE(*,6003) info
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
 6000 FORMAT(/,'LUD_Direct_Solver',T50,'[v1.0   24 February  2007]',/,   &
     &         ':::::   Direct banded matrix solver using LU decomposition')
!
 6001 FORMAT(' LU decomposition time = ',1PE12.4,'  seconds','   ==>  INFO =',I4)
!
 6002 FORMAT('      SOLUTION TIME = ',1PE12.4,'  SECONDS')
 6003 FORMAT(//,20('ERROR-'),//,T33,   &
     &             '       S I M U L A T I O N   A B O R T E D',   &
     &       /,T24,'THE UPPER TRIANGULAR MATRIX ELEMENT U(I,I) ',   &
     &             'IS EXACTLY ZERO AT I = ',I6,   &
     &       /,T38,'THE FACTORIZATION IS COMPLETED BUT DIVISION',   &
     &       /,T38,'BY ZERO WILL OCCUR IF SOLUTION IS ATTEMPTED',   &
     &       //,20('ERROR-'))
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End of LUD_Direct_Solver
!
!
      RETURN
!
      END SUBROUTINE LUD_Direct_Solver
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************!
!
!
      SUBROUTINE Print_Secondary_Variables
!
! ...... Modules to be used
!
         USE Basic_Parameters
         USE General_External_File_Units
         USE General_Control_Parameters
!
         USE Grid_Geometry, ONLY: elem
         USE Element_Attributes
         USE Solution_Matrix_Arrays, ONLY: DELX
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*          ROUTINE FOR PRINTING ALL THE SECONDARY PARAMETERS          *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
      IMPLICIT NONE
!
! -------
! ... Double precision local arrays
! -------
!
      REAL(KIND = 8), DIMENSION(NumSecondaryVar) :: DP
!
! -------
! ... Integer variables
! -------
!
      INTEGER :: NLOC,n,k,m,i,kk
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
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of Print_Secondary_Variables
!
!
      IF(First_call) THEN
         First_call = .FALSE.
         WRITE(VERS_Unit,6000)
      END IF
!
! ... Write a heading
!
      WRITE(*,6001)
!
!
!***********************************************************************
!*                                                                     *
!*                          MAIN ELEMENT LOOP                          *
!*                                                                     *
!***********************************************************************
!
!
      DO_NoElem: DO n = 1,NumElemTot
!
! -------------
! ......... Print parameters at state point
! -------------
!
         WRITE(*,6002)  elem%name(n),   &
     &                 (ElemProp(n,0)%satur(i),                                   &
     &                  ElemProp(n,0)%density(i),                                 &
     &                  ElemProp(n,0)%enthalpy(i),                                &
     &                 (ElemProp(n,0)%MassFrac(m,i),m=1,NumCom),i=1,NumPhases),   &
     &                 (ElemProp(n,0)%RelPerm(i),                                 &
     &                  ElemProp(n,0)%viscosity(i),i=1,NumMobPhases),             &
     &                  ElemState%temp(n,0)
!
! -------------
!  ........ Print incremented parameters
! -------------
!
         IF_Optionx: IF(Flag_PrintInputs .EQV. .TRUE.) THEN               ! 'Flag_PrintInputs' might need to be 'Flag_print_input_data'
!
            DO_NEq_1: DO k=2,NumEquPlus1
               WRITE(*,6003)  (ElemProp(n,k-1)%satur(i),                                   &
     &                         ElemProp(n,k-1)%RelPerm(i),                                 &
     &                         ElemProp(n,k-1)%enthalpy(i),                                &
     &                        (ElemProp(n,k-1)%MassFrac(m,i),m=1,NumCom),i=1,NumPhases),   &
     &                        (ElemProp(n,k-1)%RelPerm(i),                                 &
     &                         ElemProp(n,k-1)%viscosity(i),i=1,NumMobPhases),             &
     &                         ElemState%temp(n,k-1)
            END DO DO_NEq_1
!
! -------------
!  ........ Print derivatives
! -------------
!
         ELSE
!
            NLOC = Locp(n)                   ! = (n-1)*NumComPlus1
!
            DO_1: DO k=2,NumEquPlus1
!
               kk = 0
!
! ............ Basic parameters
!
               DO_2: DO i=1,NumPhases
!
                  kk = kk+1
                  DP(kk) = ( ElemProp(n,k-1)%satur(i) - ElemProp(n,0)%satur(i) )/DELX(NLOC+K-1)
!
                  kk = kk+1
                  DP(kk) = ( ElemProp(n,k-1)%density(i) - ElemProp(n,0)%density(i) )/DELX(NLOC+K-1)
!
! ............... Parameters related to mass fractions of components in phases
!
                  DO_2a: DO m=1,NumCom
                     kk = kk+1
                     DP(kk) = ( ElemProp(n,k-1)%MassFrac(m,i) - ElemProp(n,0)%MassFrac(m,i) )/DELX(NLOC+K-1)
                  END DO DO_2a
!
               END DO DO_2
!
! ............ Parameters related to phase mobility
!
               DO_3: DO i=1,NumMobPhases
!
                  kk = kk+1
                  DP(kk) = ( ElemProp(n,k-1)%RelPerm(i) - ElemProp(n,0)%RelPerm(i) )/DELX(NLOC+K-1)
!
                  kk = kk+1
                  DP(kk) = ( ElemProp(n,k-1)%viscosity(i) - ElemProp(n,0)%viscosity(i) )/DELX(NLOC+K-1)
!
               END DO DO_3
!
! ............ Temperature
!
               kk = kk+1
               DP(kk) = ( ElemState%temp(n,k-1) - ElemState%temp(n,0) )/DELX(NLOC+K-1)
!
               WRITE(*,6003) (DP(m),m=1,NumSecondaryVar)
!
            END DO DO_1
!
         END IF IF_Optionx
!
      END DO DO_NoElem
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
 6000 FORMAT(/,'Print_Secondary_Variables',T50,'[v1.0,  23 August    2007]',/,   &
     &         ':::::   Routine for printing all secondary parameters (stored in the derived-type array "cell")')
!
 6001 FORMAT(/,' Secondary parameters')
 6002 FORMAT(/,' Element ',A8,/,(10(1X,1pE12.5)))
 6003 FORMAT(/,(10(1X,1pE12.5)))
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of Print_Secondary_Variables
!
!
      RETURN
!
      END SUBROUTINE Print_Secondary_Variables
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE Write_File_SAVE
!
! ...... Modules to be used
!
         USE EOS_Default_Parameters
         USE Basic_Parameters
         USE General_External_File_Units

         USE General_Control_Parameters
!
         USE Grid_Geometry
         USE Element_Attributes
!
         USE Sources_and_Sinks
!
         USE Solution_Matrix_Arrays, ONLY: X
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*     ROUTINE GENERATING A FILE <SAVE> TO BE USED FOR RESTARTING      *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
      IMPLICIT NONE
!
! -------
! ... Integer variables
! -------
!
      INTEGER :: n, i, n2, ier
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
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of Write_File_SAVE
!
!
      IF(First_call) THEN
         First_call = .FALSE.
         WRITE(VERS_Unit,6000)
      END IF
!
      WRITE(*,6001) Tot_NumTimeSteps,time                          ! Write time and timestep into main output file
!
! -------
! ... Associate and print data into file SAVE
! -------
!
      REWIND (UNIT = SAVE_Unit)                                    ! Open file SAVE
      WRITE  (UNIT = SAVE_Unit, FMT = 6002) NumElemTot,time        ! Print headings
!
! -------
! ... Print primary variables into file SAVE
! -------
!
         DO_NumEle2: DO n=1,NumElemTot
!
            n2 = ElemState%Index(n,current)
!
! ...Original code written to SAVE file - delete after new code works!
    !           WRITE(SAVE_Unit,6004)  elem%name(n), ElemMedia(n,current)%porosity, State_name(n2), EOS_variables(n2), &
    ! &                               (X((n-1)*NumComPlus1+i),i=1,NumComPlus1)
!
! ...Write namelist data to SAVE file
!
         WRITE(SAVE_Unit, 6005)


!
         END DO DO_NumEle2
!
! -------
! ... Print continuation data into file SAVE using a namelist
! -------
!
      WRITE(UNIT = SAVE_Unit, FMT = 6350)                          ! Print continuation keyword ':::'
!
      WRITE(UNIT = SAVE_Unit, FMT = 6400, IOSTAT = ier) Tot_NumTimeSteps, Tot_NumNRIterations, TimeOrigin, time
!
! ... Stop if there is a problem writing the namelist
!
      IF( ier /= 0 ) THEN
         WRITE (*, FMT = 6500)
         STOP
      END IF
!
      WRITE(UNIT = SAVE_Unit, FMT = 6450)                          ! Print end-of-datablock indicator '<<<'
!
!
!
      ENDFILE (UNIT = SAVE_Unit)                                   ! Close the SAVE file
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 6000 FORMAT(/,'Write_File_SAVE',T50,'[v1.0,  14 January   2007]',/,   &
     &         ':::::   At the completion of a FTSim run, write primary variables into file "SAVE"')
!
 6001 FORMAT(/,' >>>>> Write file "SAVE" after ',i7,' time steps  ==>  The time is ',1pe13.6,' seconds'/)
 6002 FORMAT('INCON: Initial conditions for ',i7,' elements at time ',1pE14.7)
!
 6003 FORMAT(A8, 7X,1pE15.8,2x,a3,a35,/,6(1pE20.13))
!!!!!!!!!! 6004 FORMAT(A5,10X,1pE15.8,2x,a3,a35,/,6(1pE20.13)) !!!! old format for SAVE format - delete after new code works
!

 6005 FORMAT('"',a5,'"',/, &
     &       '&InCond  state   = "',a3,'",',/,              &
     &       '         phi     = ',1pe15.9,',',/            &
     &       '         perm    = ',3(1pe15.9,',',1x),/,     &
     &       '         PVList  = "',a5,'",',/,              &
     &       '         PV      = ',2(1pe20.13,1x),/,        &
     &       '         /')



!
 6013 FORMAT(A8, 7X,1pE15.8,2x,a3,a35,1x,3(e15.8),/,6(1pE20.13))
 6014 FORMAT(A5,10X,1pE15.8,2x,a3,a35,1x,3(e15.8),/,6(1pE20.13))
!
 6350 FORMAT(':::  ')
!
 6400 FORMAT('&Data_For_Continuation_Run',/,                          &
     &       '       timesteps_to_this_point     = ',i10,',',/,       &
     &       '       NR_iterations_to_this_point = ',i10,',',/,       &
     &       '       origin_of_time              = ',1pe21.14,',',/,  &
     &       '       time_to_this_point          = ',1pe21.14,',',/,  &
     &       '                                   /' )
!
 6450 FORMAT('<<<  ')
!
 6500 FORMAT(//,22('ERROR-'),//,T33,  &
     &             '       S I M U L A T I O N   A B O R T E D',/,                              &
     &         T23,'There is a problem writing the namelist <Data_For_Continuation_Run> in the <SAVE> data file', &
     &       /,T32,'               CORRECT AND TRY AGAIN',                                      &
     &       //,22('ERROR-'))
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End of Write_File_SAVE
!
!
      RETURN
!
      END SUBROUTINE Write_File_SAVE
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE SourceSink_Equation_Terms
!
! ...... Modules to be used
!
         USE Basic_Parameters
         USE General_External_File_Units

         USE General_Control_Parameters
         USE Solver_Parameters
!
         USE Grid_Geometry, ONLY: elem
         USE Sources_and_Sinks
!
         USE Solution_Matrix_Arrays
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*     ROUTINE COMPUTING ALL TERMS ARISING FROM SINKS AND SOURCES      *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
      IMPLICIT NONE
!
! -------
! ... Double precision array
! -------
!
      REAL(KIND = 8), DIMENSION(NumComPlus1 , NumEquPlus1) :: Q_contr
!
! -------
! ... Double precision variables
! -------
!
      REAL(KIND = 8) :: FAC,Gn
!
! -------
! ... Integer variables
! -------
!
      INTEGER :: n,j,k,m,jkm,np
      INTEGER :: JLOC,JLOCP,MN,JMN,JNK1,TableLength
!
! -------
! ... Character variables
! -------
!
      CHARACTER(LEN = 1) :: ITABA
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
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of SourceSink_Equation_Terms
!
!
      IF(First_call) THEN
         First_call = .FALSE.
         WRITE(VERS_Unit,6000)
      END IF
!
! ... Print additional info
!
      IF(Option_Print_SourceSinkInfo >= 1) WRITE(*,6002) Tot_NumTimeSteps, Tot_NumNRIterations
!
!
!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>                                                                     >
!>                         SOURCE/SINK LOOP                            >
!>                                                                     >
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
      DO_NumSS: DO n=1,NumSS
!
         j = SS(n)%ElemNum
         IF(j == 0 .OR. j > NumElem) GO TO 1100
         FAC = TimeIncrement/elem%vol(j)
!
         JLOC  = (j-1)*NumEqu
         JLOCP = (j-1)*NumComPlus1
!
         mn   = SS(n)%TypeIndx
         JMN  = JLOC+MN
         JNK1 = JLOC+NumComPlus1
!
!
!***********************************************************************
!*                                                                     *
!*             PRODUCTION OR INJECTION OF A MASS COMPONENT             *
!*  THE EFFECTS OF HEAT OF INJECTED/PRODUCED COMPONENT ARE CONSIDERED  *
!*                                                                     *
!***********************************************************************
!
!
         IF_MassC: IF(mn <= NumCom) THEN     ! Mass components only
!
!
! >>>>>>>>>>>>>>>>>>
! .........
! ......... Injection/production at a prescribed fixed rate (TableLength <= 1)
! .........
! >>>>>>>>>>>>>>>>>>
!
!
!            IF_QFixed: IF(TableLength <= 1) THEN
!
               Gn = SS(n)%rate
!
! -------------------
! ............ Mass injection (GN > 0), or
! ............     production (GN < 0) with prescribed enthalpy
! -------------------
!
!               IF_Enth1: IF( Gn >= 0.0d0 .OR. (GN <  0.0d0 .AND. itaba /= ' ') ) THEN
               IF_Enth1: IF( Gn >= 0.0d0) THEN
                  R(jmn) = R(jmn) - FAC*Gn
!
                  IF(NumEqu == NumComPlus1) R(JNK1) = R(JNK1) - FAC * Gn * SS(n)%enth
!
                  GO TO 1100
!
! -------------------
! ............ Mass production (GN < 0), with enthalpy to be determined
! ............    from conditions in producing block
! -------------------
!
               ELSE
!
                  CALL Source_Phase_EnthComp(n,j,Gn,FAC,Q_contr) ! Calculate phase composition and flowing enthalpy of multiphase source fluid
!
                  GO TO 1000                                     ! Go to modify the Jacobian
!
               END IF IF_Enth1
!
! <<<<<<<<<
! <<<<<<<<< End of the Fixed Rate IF
! <<<<<<<<<
!
!            END IF IF_QFixed
!
! <<<<<<
! <<<<<< End of the Mass Component IF
! <<<<<<
!
         END IF IF_MassC
!
!
!***********************************************************************
!*                                                                     *
!*                 DIRECT HEAT INJECTION/PRODUCTION                    *
!*                                                                     *
!***********************************************************************
!
!
         IF_Heat: IF(mn == NumComPlus1) THEN
!
!
! >>>>>>>>>>>>>>>>>>
! .........
! ......... Heat injection/production at a prescribed fixed rate
! .........
! >>>>>>>>>>>>>>>>>>
!
!
            IF(NumEqu == NumComPlus1) THEN            ! Ignore heat sinks/sources when not solving the energy equation
!
               r(jmn) = r(jmn) - fac*SS(n)%rate       ! Adjusting the RHS
!
            END IF
!
            GO TO 1100
!
! <<<<<<
! <<<<<< End of the Mass Component IF
! <<<<<<
!
         END IF IF_Heat
!
!
!***********************************************************************
!*                                                                     *
!*                MODIFY THE JACOBIAN MATRIX EQUATION                  *
!*                                                                     *
!***********************************************************************
!
!
 1000    DO_NeqOuter: DO k=1,NumEqu
!
            R(JLOC+k) = R(JLOC+k) + Q_contr(k,1)                                         ! Adjusting the RHS
!
            DO_NeqInner: DO m=1,NumEqu
               jkm     = (j-1)*NumEquSquare + (k-1)*NumEqu + m                           ! Locating the CO index
               CO(jkm) = CO(jkm) - ( Q_contr(k,m+1) - Q_contr(k,1) )/DELX(JLOCP+m)       ! Adjusting the CO coefficient
            END DO DO_NeqInner
!
         END DO DO_NeqOuter
!
! -------
! ... Print supporting information (when Option_Print_SourceSinkInfo >=2 only)
! -------
!
 1100    IF(Option_Print_SourceSinkInfo >= 2) WRITE(*,6004) SS(n)%ElemName, SS(n)%name, SS(n)%rate, SS(n)%enth, &
     &                                (SS(n)%PhaseFracFlow(np), np=1,NumMobPhases)
!
! <<<
! <<<
! <<< End of the SOURCE/SINK LOOP
! <<<
! <<<
!
      END DO DO_NumSS
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 6000 FORMAT(/,'SourceSink_Equation_Terms',T50,'[v1.0,  18 November  2006]',/,                             &
     &         ':::::   Assembling all source and sink terms in the equations (generic version)',/,        &
     &         '        Rigorous  step rate capability for Option_GenerationRates = 2, and capability for', &
     &         '        flowing wellbore pressure corrections')
!
!
 6002 FORMAT(/,' ==============> SUBROUTINE QU <==============   [Timestep, NR-iterations] = [',I6,',',I7,']')
!
 6004 FORMAT(' ELEMENT ',A8,'   SOURCE ',A5,'   :::   FLOW RATE = ',1pE13.6,'   SPECIFIC ENTHALPY = ',1pE13.6,/,   &
     &       ' FNP1 = ',1pE13.6,'   FNP2 = ',1pE13.6,'   FNP3 = ',1pE13.6,'   FNP4 = ',1pE13.6)
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End of SourceSink_Equation_Terms
!
!
      RETURN
!
!
      END SUBROUTINE SourceSink_Equation_Terms
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
      SUBROUTINE Source_Phase_EnthComp(n,j,Gn,fac,Q_contr)
!
! ...... Modules to be used
!
         USE Basic_Parameters
         USE General_External_File_Units
!
         USE General_Control_Parameters
!
         USE Grid_Geometry
         USE Element_Attributes
!
         USE Sources_and_Sinks
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*  ROUTINE DETERMINING PHASE COMPOSITION & ENTHALPY OF SOURCE FLUID   *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
      IMPLICIT NONE
!
! -------
! ... Quad precision array
! -------
!
      REAL(KIND = 8), INTENT(OUT), DIMENSION(NumComPlus1, NumEquPlus1) :: Q_contr
!
! -------
! ... Double precision variables
! -------
!
      REAL(KIND = 8), INTENT(IN) :: Gn,fac
!
      REAL(KIND = 8) :: FFS
!
! -------
! ... Integer variables
! -------
!
      INTEGER, INTENT(IN) :: n,j
!
      INTEGER :: ic,m,np,i
!
! -------
! ... Logical variables
! -------
!
      LOGICAL :: First_call = .TRUE.
!
#ifdef USE_TIMER
      real(KIND = 8) :: start, finish
#endif
! -------
! ... Saving variables
! -------
!
      SAVE First_call
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of Source_Phase_EnthComp
!
!
      IF(First_call) THEN
         First_call = .FALSE.
         WRITE(VERS_Unit,6000)
      END IF
!
! -------
! ... Assignment of parameter values
! -------
!
      !ic = Option_FluidComposition + 1
      ic = 1
!
!
!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>                                                                     >
!>              LOOP OVER NumEqu+1 SETS OF PARAMETERS                  >
!>                                                                     >
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
#ifdef USE_TIMER
         call CPU_Timing_Routine(start)
#endif

#ifdef USE_OMP
!$OMP PARALLEL PRIVATE (FFS)
!$OMP DO schedule(auto)
#endif
      DO_Neq1: DO m=1,NumEquPlus1
! ----------
! ...... Initialize the block Q_contr(j,j)
! ----------
#ifdef USE_OMP
!$OMP SIMD
#endif
         DO i=1,NumComPlus1
            Q_contr(i,m) = 0.0d0   ! CAREFUL - Whole array operation
         END DO
! ----------
! ...... Initialization
! ----------
         FFS        = 0.0d0                 ! = SUM(mobility*density) in all phases
         SS(n)%enth = 0.0d0                 ! Flow enthalpy
!
!
!***********************************************************************
!*                                                                     *
!*           COMPUTE FRACTIONAL FLOWS IN THE VARIOUS PHASES            *
!*                                                                     *
!***********************************************************************
!
         DO_NumPhase1: DO np=1,NumMobPhases
!
            SS(n)%PhaseFracFlow(np) = 0.0d0   ! Initialization of fractional flows
!
            IF(ic == 1 .AND. ElemProp(j,m-1)%viscosity(np) /= 0.0d0) THEN
               SS(n)%PhaseFracFlow(np) = ElemProp(j,m-1)%RelPerm(np) * ElemProp(j,m-1)%density(np)   &
     &                                  /ElemProp(j,m-1)%viscosity(np)
            END IF
!
            IF(ic == 2) THEN
               SS(n)%PhaseFracFlow(np) = ElemProp(j,m-1)%satur(np) * ElemProp(j,m-1)%density(np)
            END IF
!
! -------------
! ......... Adjsust the sum in FFS
! -------------
            FFS = FFS + SS(n)%PhaseFracFlow(np)
!
! -------------
! ......... For Option_Print_SourceSinkInfo >= 5, print info
! -------------
            IF(Option_Print_SourceSinkInfo >= 5) WRITE(*,6002) m, np, ElemProp(j,m-1)%RelPerm(np), ElemProp(j,m-1)%viscosity(np),   &
     &                                           ElemProp(j,m-1)%density(np), SS(n)%PhaseFracFlow(np), FFS
         END DO DO_NumPhase1
!
!***********************************************************************
!*                                                                     *
!*           RENORMALIZE FRACTIONAL FLOWS TO 1, AND COMPUTE            *
!*               CONTRIBUTIONS OF COMPONENTS IN PHASES                 *
!*                                                                     *
!***********************************************************************
!
         DO_NumPhase2: DO np=1,NumMobPhases
!
! ......... Normalization of fractional flows
            IF(FFS /= 0.0d0) SS(n)%PhaseFracFlow(np) = SS(n)%PhaseFracFlow(np)/FFS
!
! ......... Contribution of components to phase enthalpy
            SS(n)%enth =  SS(n)%enth + SS(n)%PhaseFracFlow(np) * ElemProp(j,m-1)%enthalpy(np)
!
! ......... Adjustment of the Jacobian submatrices:  CAREFUL! Whole array operation
#ifdef USE_OMP
! Probably won't do anything unless advanced gather instructions available
!$OMP SIMD
#endif
            DO i=1,NumCom
               Q_contr(i, m) = Q_contr(i, m) - FAC*GN*SS(n)%PhaseFracFlow(np) * ElemProp(j,m-1)%MassFrac(i,np)
            END DO
!
         END DO DO_NumPhase2
!
! ----------
! ...... Adjust the heat-related Jacobian submatrices
! ----------
         Q_contr(NumComPlus1,m) = Q_contr(NumComPlus1,m) - FAC*GN*SS(n)%enth
!
! <<<
! <<< End of the NumEquPlus1 LOOP
! <<<
!
      END DO DO_Neq1
#ifdef USE_OMP
!$OMP END DO
!$OMP END PARALLEL
#endif

#ifdef USE_TIMER
         call CPU_Timing_Routine(finish)

         write (*,*) __FILE__, ":", __LINE__, " time: ", finish-start
#endif
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 6000 FORMAT(/,'Source_Phase_EnthComp',T50,'[v1.0,  27 April     2007]',/,   &
     &         ':::::   Determine phase composition and enthalpy of source fluid')
!
 6002 FORMAT(' M=',I2,' NP=',I1,' K=',1pE13.6,' VIS=',1pE13.6,   &
     &       ' D=',1pE13.6,' FF=',1pE13.6,' FFS=',1pE13.6)
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End of Source_Phase_EnthComp
!
!
      RETURN
!
!
      END SUBROUTINE Source_Phase_EnthComp
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE CPU_Elapsed_Time(IFLAG,DT,TOTAL)
!
! ...... Modules to be used
!
         USE General_External_File_Units
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*                 ROUTINE COMPUTING ELAPSED CPU TIME                  *
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
      REAL(KIND = 8), INTENT(OUT) :: DT,TOTAL
!
      REAL(KIND = 8) :: T1,T2
!
! -------
! ... Integer variables
! -------
!
      INTEGER, INTENT(IN) :: IFLAG
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
      SAVE First_call,T1
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of CPU_Elapsed_Time
!
!
      IF(First_call) THEN
         First_call = .FALSE.
         WRITE(VERS_Unit,6000)
      END IF
!
! -------
! ... Call the CPU timing routine
! -------
!
      CALL CPU_Timing_Routine(T2)
!
! -------
! ... Compute elapsed time
! -------
!
      IF(IFLAG /= 0) THEN
         DT    = T2-T1
         TOTAL = T2
      ELSE
         T1 = T2
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
 6000 FORMAT(/,'CPU_Elapsed_Time',T50,'[v1.0,  15 February  2007]',6X,/,   &
     &         ':::::   Compute elapsed CPU time')
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of CPU_Elapsed_Time
!
!
      RETURN
!
!
      END SUBROUTINE CPU_Elapsed_Time
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE Warning(ELEM,XX,NK1,KC)
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*       ROUTINE FOR PRINTING WARNINGS AND RELATED INFORMATION         *
!*   WHEN THE RANGE IS EXCEEDED AND/OR INITILIALIZATION IS INCORRECT   *
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
      INTEGER, INTENT(IN) :: NK1,KC
!
      INTEGER :: m
!
! ----------
! ... Double precision arrays
! ----------
!
      REAL(KIND=8), DIMENSION(NK1), INTENT(IN) :: XX
!
! ----------
! ... Character variables
! ----------
!
      CHARACTER(LEN=5), INTENT(IN) :: ELEM
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  BEGIN Warning
!
!
      WRITE(*,6100) ELEM,(XX(m), m=1,NK1)
!
      IF(KC == 0) THEN
         WRITE(*,6101)
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
 6100 FORMAT(' !!!!!!!!! Cannot find parameters at element "',A8,'" with primary variables  XX(M) =',10(1X,1pE12.5))
!
 6101 FORMAT(' !!!!!!!!!!!  Erroneous Data Initialization  !!!!!!!!!!',11X,'STOP EXECUTION  <=========')
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  END Warning
!
!
      RETURN
!
      END SUBROUTINE Warning
!
!
!
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
