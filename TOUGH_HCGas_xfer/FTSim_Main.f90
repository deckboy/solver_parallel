!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>                                                                     >
!>                                                                     >
!>        FTSim_Main.f95: Code unit including the main program         >
!>                      and all related routines                       >
!>                                                                     >
!>                                                                     >
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  Main program that organizes the calling sequence of the high-level
!  events in the simulation process, and includes the writing of
!  important general comments in the standard output files, timing
!  procedures, and handling of files needed by the code and/or created
!  during the code execution
!
!
      PROGRAM FT_Simulator
!
! ...... Modules to be used
!
         USE Basic_Parameters
         USE General_External_File_Units
!
         USE General_Control_Parameters
!
!***********************************************************************
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*       FTSim Main Program: Calls several lower-level routines        *
!*                    which execute computations                       *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!***********************************************************************
!
      IMPLICIT NONE
!
! -------
! ... Character allocatable arrays
! -------
!
      CHARACTER(LEN = 5), ALLOCATABLE, DIMENSION(:) :: VV
!
! -------
! ... Double Precision Variables
! -------
!
      REAL(KIND = 8) :: ELT, ELT1, ELTC
!
! -------
! ... Integer Variables
! -------
!
      INTEGER :: i, ierG
!
! -------
! ... Logical Variables
! -------
!
      LOGICAL :: Flag_MINC
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of <FTSim>
!
!
!
! -------
! ... Open file VERS to record subroutine use
! -------
!
      INQUIRE(FILE = 'VERS',EXIST = file_VERS_exists)
!
      IF_VERS: IF(file_VERS_exists) THEN
                  OPEN(VERS_Unit, FILE = 'VERS', STATUS = 'OLD')
               ELSE
                  OPEN(VERS_Unit, FILE = 'VERS', STATUS = 'NEW')
                  ENDFILE (UNIT = VERS_Unit)
               END IF IF_VERS
!
      REWIND (UNIT = VERS_Unit)
!
      WRITE(VERS_Unit,6000)
!
! -------
! ... Reading the heading/title of the input data file
! -------
!
      READ (*, FMT = 6002) title
      WRITE(*, FMT = 6003) title
!
!
!***********************************************************************
!*                                                                     *
!*                Allocate memory for most arrays                      *
!*                                                                     *
!***********************************************************************
!
!
      CALL Core_Memory_Allocation
!
!
!***********************************************************************
!*                                                                     *
!*          Print headings and important general information           *
!*                                                                     *
!***********************************************************************
!
!
      CALL Write_Headings
!
!
!
!***********************************************************************
!*                                                                     *
!*               Check the directory for various files                 *
!*                                                                     *
!***********************************************************************
!
!
      CALL Open_External_IO_Files
!
!
!***********************************************************************
!*                                                                     *
!*                       Initialize the clock                          *
!*                                                                     *
!***********************************************************************
!
!
      CALL CPU_Timing_Routine(CPU_InitialTime)
!
!
!***********************************************************************
!*                                                                     *
!*             Allocate memory to temporary input arrays               *
!*                                                                     *
!***********************************************************************
!
!
      CALL Allocate_Temp_Input_Memory
!
!
!
!***********************************************************************
!*                                                                     *
!*             Read data from the main input data file                 *
!*                                                                     *
!***********************************************************************
!
!
      CALL Read_Main_Input_File(Flag_MINC)
!
      IF(NoFlowSimulation .EQV. .TRUE.) THEN
         ELT1 = 0.0d0
         GO TO 1000
      END IF
!
!
! -------
! ... Determine the floating point accuracy of the processor
! -------
!
      CALL Significant_Digits
!
!
!***********************************************************************
!*                                                                     *
!*      Read input data from the various files, both pre-existing      *
!*                  and created by the INPUT routine                   *
!*                                                                     *
!***********************************************************************
!
!
      CALL Read_External_Input_Files(Flag_MINC)
!
!
!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!*                                                                     *
!*                      CONDUCT THE SIMULATION                         *
!*                                                                     *
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
      IF_NoFloSimul:  IF(NoFlowSimulation .EQV. .FALSE.) THEN   ! Check if the simulation ...
!                                                               ! ... is to be bypassed
! ----------
! ... Print additional information
! ----------
!
         WRITE(*,6004) NumElemTot, NumElem, NumConx, NumSS
!
!
!***********************************************************************
!*                                                                     *
!*             Allocate memory to matrix solver arrays                 *
!*                                                                     *
!***********************************************************************
!
!
         CALL Solver_Memory_Allocation
!
!
! ----------
! ...... Print input data (when Flag_print_input_data /= 0)
! ----------
!
         IF(Flag_PrintInputs .NEQV. .TRUE.) CALL Print_Input_Data                       !'Flag_PrintInputs' might be 'Flag_print_input_data'
!
! ----------
! ...... Deterine elapsed time for data input
! ----------
!
         CALL CPU_Timing_Routine(ELT1)
         ELT1 = ELT1 - CPU_InitialTime
!
         WRITE(*,6006) ELT1
!
! ----------
! ...... Deallocate unnecessary memory (Temporary input arrays)
! ----------
!
         CALL Deallocate_Unneeded_Arrays
!
!
! ----------
! ...... Print additional information
! ----------
!
         CALL Write_More_Generic_Info
!
         CALL Write_EOS_Specific_Info
!
!
!
!***********************************************************************
!*                                                                     *
!*     Perform simulation by solving the mass and heat equations       *
!*                                                                     *
!***********************************************************************
!
!
         CALL Simulation_Cycle
!
!
! <<<
! <<< End of the "IF_NoFloSimul" construct
! <<<
!
      END IF IF_NoFloSimul
!
!
!***********************************************************************
!*                                                                     *
!*                   Print-out of version information                  *
!*                                                                     *
!***********************************************************************
!
!
 1000 IF_NoVersion: IF(NoVersion .EQV. .FALSE.) THEN  ! Check whether to bypass
!
! ----------
! ...... Print information
! ----------
!
         WRITE(*, FMT = 6008)
!
! ...... Close the VERS file
!
         ENDFILE (UNIT = VERS_Unit)
!
! ...... Go to the top of the VERS file
!
         REWIND (UNIT = VERS_Unit)
!
! ----------
! ...... Allocate memory to the VV local array
! ----------
!
         ALLOCATE(VV(26), STAT=ierG)
!
! ...... Check if properly allocated
!
         IF(ierG /= 0) THEN   ! Unsuccesful allocation
            WRITE(*, FMT = 6101)
            STOP
         END IF
!
! ----------
! ...... Initialization of the VV local array
! ----------
!
         DO_Version: DO
!
! ......... Initialization of the VV local array
!
            VV = '     '     !  CAREFUL - Whole array operation
!
! ......... Read version information
!
            READ(UNIT = VERS_Unit, FMT = 6010, IOSTAT = ierG) (VV(i),i=1,26)
!
! ......... For EOF/EOR, exit the loop
!
            IF(ierG /= 0) EXIT DO_Version
!
! ......... Write version information
!
            WRITE(*, FMT = 6010) (VV(i),I=1,26)
!
         END DO DO_Version
!
! ----------
! ...... Allocate memory from the VV local array
! ----------
!
         DEALLOCATE (VV, STAT=ierG)
!
! ...... Check if properly deallocated
!
         IF(ierG /= 0) THEN   ! Unsuccesful deallocation
            WRITE(*, FMT = 6102)
            STOP
         END IF
!
!
      END IF IF_NoVersion
!
!
!***********************************************************************
!*                                                                     *
!*               Compute execution timing information                  *
!*                                                                     *
!***********************************************************************
!
!
      WRITE(*, FMT = 6015)
!
      CALL CPU_Timing_Routine(ELT)
!
      ELT  = ELT-CPU_InitialTime
      ELTC = ELT-ELT1
!
      WRITE(*, FMT = 6020) ELT, ELT1, ELTC, CPU_MatrixSolTime   !  Print execution time data
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 6000 FORMAT(/,'FT_Simulator',T50,'[v1.0,  28 October   2007]',/,   &
     &         ':::::   The <FT_Simulator> main program (calls all the high-level routines)')
!
 6002 FORMAT(A120)
 6003 FORMAT(/1X,131('=')//5X,'PROBLEM TITLE:  ',A120//)


 6004 FORMAT(/' MESH HAS',I7,' ELEMENTS (',I7,' ACTIVE) AND',I7,' CONNECTIONS (INTERFACES)',   &
     &        ' BETWEEN THEM'/' GENER HAS',I6,' SINKS/SOURCES')
!
 6006 FORMAT(//,' END OF <FT_Simulator> INPUT JOB: Elapsed Time = ',1pe15.8,' seconds'/)
!
 6008 FORMAT(/,   &
     &       132('*'),/,'*',130X,'*',/,                                   &
     &       '*',50X,'SUMMARY OF PROGRAM UNITS USED',51X,'*',/,           &
     &       '*',130X,'*',/,132('*'),/)
!
 6010 FORMAT(26A5)
 6015 FORMAT(//,132('*'))
 6020 FORMAT(//, ' END OF <FTSim> SIMULATION RUN:',   &
     &       T35,' Elapsed Time             = ',1pe15.8,' sec',/,   &
     &       T35,' Data Input Time          = ',1pe15.8,' sec',/,   &
     &       T35,' Calculation Time         = ',1pe15.8,' sec',/,   &
     &       T40,' Matrix Solving Time = ',1pe15.8,' sec')
!
 6101 FORMAT(//,20('ERROR-'),//,   &
     &       T2,'Memory allocation of array VV in the main program <FTSim> was unsuccessful',//,   &
     &       T2,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
 6102 FORMAT(//,20('ERROR-'),//,   &
     &       T2,'Memory deallocation of array VV in <FTSim> was unsuccessful',//,   &
     &       T2,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 8000 FORMAT(///,T5,'!!!!! >>>>>  ',6('QUE PASSA HOMBRES?  '),///, &
     &           T5,'!!!!! >>>>>   I   W A N T   T O   G O   O N   B U T   Y O U   D O   N O T   L E T   M E')
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End of <FTSim>
!
!
      STOP
      END PROGRAM FT_Simulator
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE Write_Headings
!
! ...... Modules to be used
!
         USE EOS_Parameters
!
         USE Basic_Parameters
         USE General_External_File_Units
         USE General_Control_Parameters
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*          ROUTINE FOR PRINTING HEADINGS IN THE OUTPUT FILE           *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
      IMPLICIT NONE
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of Write_Headings
!
!
      WRITE(VERS_Unit,6000)
!
! ... Format
!
 6000 FORMAT(/,'Write_Headings',T50,'[v1.0,  23 April     2007]',/,  &
     &         ':::::   Routine for printing the headings in the standard <FTSim> output file')
!
! -------
! ... Print headers
! -------
!
      WRITE(*, 6001) ' '
!
! ... Format
!
 6001 FORMAT(A1,/,   &
     &       T21,'@@@@@  @@@@@          @@@  @  @   @  @  @  @      @@   @@@@@   @@   @@@      @@@   @  @  @   @',/, &
     &       T21,'@        @           @     @  @@ @@  @  @  @     @  @    @    @  @  @  @     @  @  @  @  @@  @',/, &
     &       T21,'@@@      @    @@@@@   @@   @  @ @ @  @  @  @     @@@@    @    @  @  @@@      @@@   @  @  @ @ @',/, &
     &       T21,'@        @              @  @  @   @  @  @  @     @  @    @    @  @  @ @      @ @   @  @  @  @@',/, &
     &       T21,'@        @           @@@   @  @   @   @@   @@@@  @  @    @     @@   @  @     @  @   @@   @   @',//)
!
!
!
      WRITE(*, 6002)
!
! ... Format
!
 6002 FORMAT(20X,'FTSim IS A PROGRAM FOR MULTIPHASE MULTICOMPONENT FLOW IN PERMEABLE MEDIA, INCLUDING HEAT FLOW',/,   &
     &       28X,'IT WAS DEVELOPED BY THE AWESOME DUDES OF THE PETROLEUM ENGINEERING DEPARTMENT ',/,    &
     &       53X,'AT TEXAS A&M UNIVERSITY')
!
      WRITE(*, 6003)
 6003 FORMAT(/,132('*'),/)
!
!
!
      WRITE(*, 6005)
!
! ... Format
!
 6005 FORMAT(/26X,80('*'),/,26X,20('*'),40X,20('*'),/,                   &
     &        26X,20('*'),7X,'FTSim V1.0 (MONTH___ 20__)',7X,20('*')/,   &
     &        26X,20('*'),40X,20('*'),/,26X,80('*'),/)
!
!
!
      WRITE(*, 6008) Max_NumElem, Max_NumConx, Max_NumEquations, Max_NumMassComp, Max_NumPhases, &
     &               Max_NumSS
!
! ... Format
!
 6008 FORMAT(//,' PARAMETERS FOR DYNAMIC DIMENSIONING OF MAJOR ARRAYS (MAIN PROGRAM) ARE AS FOLLOWS',//, &
     &          T10,'Max_NumElem      =',I8,/,T10,'Max_NumConx      =',I8,/,    &
     &          T10,'Max_NumEquations =',I8,/,T10,'Max_NumMassComp  =',I8,/,    &
     &          T10,'Max_NumPhases    =',I8,/,T10,'Max_NumSS        =',I8,/,    &
     &      ' ',131('='))
!
!
!
      WRITE(*, 6009) Max_NumElem, Max_NumConx, Max_NumPrimaryVar, Max_NumSS, Max_NumMedia
!
! ... Format
!
 6009 FORMAT(/,' MAXIMUM NUMBER OF VOLUME ELEMENTS (GRID BLOCKS)   :',9X,'Max_NumElem        =',I9,/,   &
     &         ' MAXIMUM NUMBER OF CONNECTIONS (INTERFACES)        :',9X,'Max_NumConx        =',I9,/,   &
     &         ' MAXIMUM LENGTH OF PRIMARY VARIABLE ARRAYS         :',9X,'Max_NumPrimaryVar  =',I9,/,   &
     &         ' MAXIMUM NUMBER OF GENERATION ITEMS (SINKS/SOURCES):',9X,'Max_NumSS          =',I9,/,   &
     &         ' MAXIMUM NUMBER OF POROUS/FRACTURED MEDIA:          ',9X,'Max_NumMedia       =',I9,//,' ',131('='))
!
!
!
      WRITE(*, 6012)
!
! ... Format
!
 6012 format(/,' Array dimensioning is dynamic.  The use of direct solvers will reduce the size of the maximum solvable problem')
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End of Write_Headings
!
!
      RETURN
!
      END SUBROUTINE Write_Headings
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE Write_More_Generic_Info
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
!*       ROUTINE FOR PRINTING ADDITIONAL INFO IN THE OUTPUT FILE       *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
      IMPLICIT NONE
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of Write_More_Generic_Info
!
!
      WRITE(VERS_Unit,6000)
!
! ... Format
!
 6000 FORMAT(/,'Write_More_Generic_Info',T50,'[v1.0,  23 April     2007]',/,  &
     &         ':::::   Routine for printing additional information in the output file')
!
! -------
! ... Print headers
! -------
!
      WRITE(*,6002)
!
! ... Format
!
 6002 FORMAT(' ',/' ',131('*'))
!
!
!
      WRITE(*,6003)
!
! ... Format
!
 6003 FORMAT(' *          ARRAY "MOP"',   &
     &       ' ALLOWS TO GENERATE MORE PRINTOUT IN VARIOUS SUBROUTINES, AND TO MAKE SOME CALCULATIONAL CHOICES            *')
!
!
!
      WRITE(*,6004)
!
! ... Format
!
 6004 FORMAT(' ',131('*'),/)
!
!
!
      WRITE(*,6005) Option_Print_NRIterationInfo
!
! ... Format
!
 6005 FORMAT(' ','  Option_Print_NRIterationInfo  =',I2,                                       &
     &       ' ==> ALLOWS TO GENERATE A SHORT PRINTOUT FOR EACH NEWTON-RAPHSON ITERATION',/,   &
     &       '           = 0, 1, OR 2: Generate 0, 1, or 2 lines of printout',//,12X)
!
      WRITE(*,6006) Option_Print_CyclingInfo, Option_Print_JacobianInfo, Option_Print_SourceSinkInfo,   &
      &             Option_Print_EOSInfo, Option_Print_SolverInfo
!
! ... Format
!
 6006 FORMAT('   Option_Print_CyclingInfo  =',I2,' ==> CYCIT;  Option_Print_JacobianInfo =',I2,' ==> MULTI;  ',                  &
     &       'Option_Print_SourceSinkInfo =',I2,' ==> QU;  Option_Print_EOSInfo =',I2,' ==> EOS;  Option_Print_SolverInfo =',I2, &
     &       ' ==> LINEQ   ',/)
!
!
!
      WRITE(*,6008) Flag_PrintInputs               !'Flag_PrintInputs' might need to be 'Flag_print_input_data'
!
! ... Format
!
 6008 FORMAT('   Flag_print_input_data  =',L2,' ==> IF UNEQUAL ZERO, THE RUN WILL GENERATE A PRINTOUT OF INPUT DATA',/)
!
!
!
      WRITE(*,6009) Option_RelPerm_CapPress,Option_FluidComposition,Option_ThermalConductivity
!
! ... Format
!
 6009 FORMAT(' ',8X,'   CALCULATIONAL CHOICES OFFERED BY Option_ ARE AS FOLLOWS:',//,                          &
     & '   Option_RelPerm_CapPress  =',I2,' ==> METHODS FOR RELATIVE PERMEABILITY',                            &
     &                       ' AND CAPILLARY PRESSURE ESTIMATION IN THE PRESSENCE OF SOLID PHASES',/,          &
     & '           = 0: Relative permeability model based on the Original Porous Medium (OPM) concept, ',      &
     &                 'scaling of capillary pressures (EPM #1)',/,                                            &
     & '           = 1: Relative permeability model based on the Evolving Porous Medium (EPM) #1 concept, ',   &
     &                 'scaling of capillary pressures (EPM #1)',/,                                            &
     & '           = 2: Relative permeability model based on the Evolving Porous Medium (EPM) #2 concept, ',   &
     &                 'scaling of capillary pressures (EPM #2)',/,                                            &
     & '           = 3: Relative permeability model based on the Evolving Porous Medium (EPM) #1 concept, ',   &
     &                 'no scaling of capillary pressures ',/,                                                 &
     & '           = 4: Relative permeability model based on the Evolving Porous Medium (EPM) #2 concept, ',   &
     &                 'no scaling of capillary pressures ',/,                                                 &
     & '           = 9: Relative permeability model based on the Original Porous Medium (OPM) concept, ',      &
     &                 'no scaling of capillary pressures',//,                                                 &
     & '   Option_FluidComposition  =',I2,' ==> CHOOSES FLUID COMPOSITION ON WITHDRAWAL (PRODUCTION)',/,       &
     & '           = 0: According to relative mobilities',/,                                                   &
     & '           = 1: According to composition in producing element',//,                                     &
     & '   Option_ThermalConductivity =',I2,' ==> CHOOSES INTERPOLATION FORMULA FOR DEPENDENCE OF THERMAL CONDUCTIVITY ON PHASE SATURATIONS ',/,  &
     & '           = 0: k = k_dry + [sqrt(S_aqu) + sqrt(S_hyd)]*(k_wet - k_dry) + phi*S_ice*k_ice',/,          &
     & '           = 1: k = k_dry + (S_aqu + S_hyd)*(k_wet - k_dry) + phi*S_ice*k_ice',/,                      &
     & '           = 2: k computed from component thermal conductivities (gas contribution is neglected)',/,   &
     & '           = 3: k computed from all component thermal conductivities (gas contribution is included)',/)
!
!
      WRITE(*,6010) Option_AbsolutePerm
!
! ... Format
!
 6010 FORMAT(  &
   & '   Option_AbsolutePerm =',I2,' ==> CHOOSES EVALUATION OF MOBILITY AND ABSOLUTE PERMEABILITY AT INTERFACES ' ,/,            &
   & '           = 0: Mobilities are upstream weighted with <W_upstream> (default is 1.0). Permeability is upstream weighted',/, &
   & '           = 1: Mobilities are averaged between adjacent elements. Permeability is upstream weighted',/,                   &
   & '           = 2: Mobilities are upstream weighted with <W_upstream> (default is 1.0). Permeability is harmonic weighted,',/,&
   & '           = 3: Mobilities are averaged between adjacent elements. Permeability is harmonic weighted',/,                   &
   & '           = 4: Mobility * permeability product is harmonic weighted',/)
!
!
!
      WRITE(*,6011) Option_GenerationRates,Option_PermAdjustment
!
! ... Format
!
 6011 FORMAT('   Option_GenerationRates =',I2,' ==> CHOOSES PROCEDURE FOR INTERPOLATING GENERATION RATES FROM A TIME TABLE.',/,  &
     &'           = 0: Triple linear interpolation',/,                                                                           &
     &'           = 1: "Step function" option ',//,                                                                              &
     &       '   Option_PermAdjustment =',I2,' ==> CHOOSES METHOD OF PERMEABILITY ADJUSTMENT IN CASES OF VARIABLE POROSITY ',    &
     &       'AND/OR GEOMECHANICAL CHANGES',/,                                                                                   &
     &       '           = 0: No permeability adjustment'/                                                                       &
     &       '           = 1: Permeability adjusted according to an empirical exponential model',/,                              &
     &       '           = 2: Permeability adjusted according to a full geomechanical model',/)
!
!
!
      WRITE(*,6012) Option_GasSolubility, Option_HeatExchange
!
! ... Format
!
 6012 FORMAT("   Option_GasSolubility =",I2," ==> SPECIFIES THE HANDLING OF GAS SOLUBILITY",/,                                  &
     &       "           = 0: Use Henry's constants (T-indepenent)",/,                                                          &
     &       "           = 1: Use appropriate equations of Henry's constants (T-depenent)",/,                                   &
     &       "           > 1: Use fugacities (T- and P-depenent)",//,                                                           &
     &       "   HeatExchange =",I2," ==> ALLOWS TO SELECT A SEMI-ANALYTICAL HEAT EXCHANGE CALCULATION WITH CONFINING BEDS",/,  &
     &       "           = 0: No semi-analytical heat exchange",/,                                                              &
     &       "           > 0: Semi-analytical heat exchange engaged (when a special sub-routine 'qloss' is present)",/)
!
!
!
      WRITE(*,6014) DoublingDtCriterion, Option_BinGasDiffus, Option_InterfaceDensity
!
! ... Format
!
 6014 FORMAT(' ',   &
     &       '  DoublingDtCriterion =',I2,' ==> PERMITS TO CHOOSE TIME STEP SELECTION OPTION',/,              &
     &       '           = 0: Use time steps explicitly provided as input',/,                                   &
     &       '           > 0: Increase time step by at least a factor 2, if convergence occurs in',             &
     &       '                < Option_TimeStepSelect iterations',//,                                           &
     &       '  Option_BinGasDiffus =',I2,' ==> SPECIFIES THE HANDLING OF BINARY GAS DIFFUSIVITIES',/,          &
     &       '           = 0: Compute using the method of Fuller et al. [1969]',/,                              &
     &       '           > 7: Adjust for high pressures using the method of Riazi and Whitson [1993]',//,       &
     &       '  Option_InterfaceDensity =',I2,' ==> ALLOWS TO SELECT HANDLING OF INTERFACE DENSITY',/,          &
     &       '           = 0: Perform upstream weighting for interface density',/,                              &
     &       '           > 0: Compute interface density as average of the two grid block densities',/,          &
     &       '                however, when one of the two phase saturations is zero, do upstream weighting',/)
!
!
!
      WRITE(*,6016) Option_OutputAmount,Option_InConCheck,Option_LineEqSolver
!
! ... Format
!
 6016 FORMAT(                                                                                                                  &
     & '   Option_OutputAmount =',I2,' ==> CONTROLS THE TYPE(S) OF OUTPUT',/,                                                  &
     & '          < 8: The run produces a standard <FTSim> output',/,                                                          &
     & '          = 8: In addition to the standard <FTSim> output, the run produces the files',                                &
     &               ' "Plot_Data_Elem" and "Plot_Data_Conx" ',/,                                                              &
     &               ' (of the most important element and connection variables) in plotter-ready format',/,                    &
     & '          = 9: No standard <FTSim> output - The run produces only the files "Plot_Data_Elem" and "Plot_Data_Conx" ',/, &
     &               ' (of the most important element and connection variables) in plotter-ready format',//,                   &
     & '   Option_InConCheck =',I2,' ==> CONTROLS THE CHECKING OF INITIAL CONDITIONS ',//,                                     &
     &           ' < 9: The initial conditions are checked to ensure physically meaningful primary variable values',/,         &
     &           ' = 9: No checking of initial conditions is performed ',//,                                                   &
     & '   Option_LineEqSolver =',I2,' ==> PERMITS TO SELECT LINEAR EQUATION SOLVER',/,                                        &
     & '          = 0: DEFAULTS TO Option_LineEqSolver = 3',/,                                                                 &
     & '          = 1: SUBROUTINE LUBAND: Direct solver using LU decomposition',/,                                             &
     & '          = 2: SUBROUTINE DSLUBC: Bi-conjugate gradient solver; Preconditioner: Incomplete LU factorization'/,         &
     & '          = 3: SUBROUTINE DSLUCS: Bi-conjugate gradient solver - Lanczos type; Preconditioner:',                       &
     &                                  ' Incomplete LU factorization',/,                                                      &
     & '          = 4: SUBROUTINE DSLUGM: Generalized minimum residual conjugate gradients; Preconditioner:',                  &
     &                                  ' Incomplete LU factorization',/,                                                      &
     & '          = 5: SUBROUTINE DLUSTB: Stabilized bi-conjugate gradient solver; Preconditioner: Incomplete LU factorization',//)
!
!
      WRITE(*,6018) Option_AddOutputFile,Option_SelectDiffOption
!
! ... Format
!
 6018 FORMAT('   Option_AddOutputFile =',I2,' ==> DETERMINES IF AN ADDITIONAL OUTPUT FILE OF max{P} WILL BE PRINTED',/,   &
     & '           = 9: A time series of the maximum pressure in the domain and its ',                                    &
     &                 'location is printed in the file "MaxP_Time_Series"',/,                                            &
     & '           < 9: No max{P} output is provided ',//,                                                                &
     &       '   Option_SelectDiffOption =',I2,' ==> Used in some EOS-modules for selecting different options',//)
!
!
!
      WRITE(*,6020) Option_MultiDiffuseFlux
!
! ... Format
!
 6020 format('   Option_MultiDiffuseFlux  =',I2,' ==> PERMITS TO SELECT HANDLING OF MULTIPHASE DIFFUSIVE FLUXES AT INTERFACES',/, &
     &       '           = 0: Harmonic weighting of fully-coupled effective multiphase diffusivity',/,                            &
     &       '           = 1: Separate harmonic weighting of gas and liquid phase diffusivities'//                                &
     &       ' ',131('*'))
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=> End of Write_More_Generic_Info
!
!
      RETURN
!
      END SUBROUTINE Write_More_Generic_Info
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE Open_External_IO_Files
!
! ...... Modules to be used
!
         USE Basic_Parameters, ONLY: file_VERS_exists
         USE General_External_File_Units
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*      ROUTINE FOR CHECKING/OPENING THE FILES NEEDED BY <FTSim>       *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
      IMPLICIT NONE
!
! ----------
! ... Logical variables
! ----------
!
      LOGICAL :: exists
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of Open_External_IO_Files
!
!
      WRITE(*,6001)            ! Write comment in the output file
!
      WRITE(VERS_Unit,6000)    ! Write comment in file VERS
!
      IF_VERS: IF(file_VERS_exists) THEN
                  WRITE(*,6002)
               ELSE
                  WRITE(*,6003)
               END IF IF_VERS
!
!
! -------
! ... File MESH
! -------
!
      INQUIRE(FILE='MESH',EXIST = exists)
!
      IF_MESH: IF(exists) THEN
                  WRITE(*,6004)
                  OPEN(MESH_Unit, FILE = 'MESH', STATUS = 'OLD')
               ELSE
                  WRITE(*,6005)
                  OPEN(MESH_Unit, FILE = 'MESH', STATUS = 'NEW')
               END IF IF_MESH
!
! -------
! ... File INCON
! -------
!
      INQUIRE(FILE='INCON',EXIST = exists)
!
      IF_INCON: IF(exists) THEN
                   WRITE(*,6006)
                   OPEN(INCON_Unit, FILE = 'INCON', STATUS = 'OLD')
                ELSE
                   WRITE(*,6007)
                   OPEN(INCON_Unit, FILE = 'INCON', STATUS = 'NEW')
                   ENDFILE (UNIT = INCON_Unit)
                END IF IF_INCON
!
! -------
! ... File GENER
! -------
!
      INQUIRE(FILE='SINKS&SOURCES', EXIST = exists)
!
      IF_SSExists: IF(exists .EQV. .FALSE.) THEN
!
         INQUIRE(FILE='GENER', EXIST = exists)
         IF_GENER: IF(exists) THEN
            WRITE(*,6008)
            OPEN(UNIT = GENER_Unit, FILE = 'GENER', STATUS = 'OLD')
         ELSE
            WRITE(*,6009)
            OPEN(GENER_Unit, FILE = 'GENER', STATUS = 'NEW')
            ENDFILE (UNIT = GENER_Unit)
         END IF IF_GENER
!
      END IF IF_SSExists
!
! -------
! ... File SAVE
! -------
!
      INQUIRE(FILE='SAVE',EXIST = exists)
!
      IF_SAVE: IF(exists) THEN
                  WRITE(*,6010)
                  OPEN(SAVE_Unit, FILE = 'SAVE', STATUS = 'OLD')
               ELSE
                  WRITE(*,6011)
                  OPEN(SAVE_Unit, FILE = 'SAVE', STATUS = 'NEW')
               END IF IF_SAVE
!
! -------
! ... File LINEQ
! -------
!
      INQUIRE(FILE='LINEQ',EXIST = exists)
!
      IF_LINEQ: IF(exists) THEN
                   WRITE(*,6012)
                   OPEN(LINEQ_Unit, FILE = 'LINEQ', STATUS = 'OLD')
                ELSE
                   WRITE(*,6013)
                   OPEN(LINEQ_Unit, FILE = 'LINEQ', STATUS = 'NEW')
                END IF IF_LINEQ
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
 6000 FORMAT(/,'Open_External_IO_Files',T50,'[v1.0,  14 January   2006]',/,  &
     &         ':::::   Open files <VERS>, <MESH>, <INCON>, <GENER>, <SAVE>, <LINEQ>, and <TABLE>')
!
 6001 FORMAT(//' SUMMARY OF DISK FILES'/)
!
 6002 FORMAT(' File <VERS> exists ==> Open as an old file')
 6003 FORMAT(' File <VERS> does not exist ==> Open as a new file')
!
 6004 FORMAT(' File <MESH> exists ==> Open as an old file')
 6005 FORMAT(' File <MESH> does not exist ==> Open as a new file')
!
 6006 FORMAT(' File <INCON> exists ==> Open as an old file')
 6007 FORMAT(' File <INCON> does not exist ==> Open as a new file')
!
 6008 FORMAT(' File <GENER> exists ==> Open as an old file')
 6009 FORMAT(' File <GENER> does not exist ==> Open as a new file')
!
 6010 FORMAT(' File <SAVE> exists ==> Open as an old file')
 6011 FORMAT(' File <SAVE> does not exist ==> Open as a new file')
!
 6012 FORMAT(' File <LINEQ> exists ==> Open as an old file')
 6013 FORMAT(' File <LINEQ> does not exist ==> Open as a new file')
!
 6014 FORMAT(' File <TABLE> exists ==> Open as an old file')
 6015 FORMAT(' File <TABLE> does not exist ==> Open as a new file')
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End of Open_External_IO_Files
!
!
      RETURN
!
      END SUBROUTINE Open_External_IO_Files
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE Significant_Digits
!
! ...... Modules to be used
!
         USE General_Control_Parameters
!
         USE General_External_File_Units
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*    ROUTINE FOR CALCULATING THE NUMBER OF SIGNIFICANT DIGITS FOR     *
!*           FLOATING POINT PROCESSING; ASSIGN DEFAULT FOR             *
!*     <derivative_increment>, AND PRINT APPROPRIATE WARNING WHEN      *
!*                 MACHINE ACCURACY IS INSUFFICIENT                    *
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
      REAL(KIND = 8) :: DF
!
! -------
! ... Integer variables
! -------
!
      INTEGER :: N10
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=> Main body of Significant_Digits
!
!
      WRITE(VERS_Unit,6000)
!
! -------
! ... Machine accuracy from FORTRAN95 intrinsic functions
! -------
!
      N10 = PRECISION(1.0d-12)
      DF  = SQRT(EPSILON(1.0d-12))
!
      WRITE(*,6003) N10,DF
!
! -------
! ... Determine the derivative increment
! -------
!
      IF_DFAC: IF(derivative_increment == 0.0d0) THEN
                  derivative_increment = DF
                  WRITE(*,6010)
               ELSE
                  WRITE(*,6011) derivative_increment
               END IF IF_DFAC
!
      WRITE(*,6016)
!
! -------
! ... Additional print-outs
! -------
!
      IF(N10 <= 12 .AND. N10 > 8) WRITE(*,6024)
      IF(N10 <= 8) WRITE(*,6025)
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 6000 FORMAT(/,'Significant_Digits',T50,'[v1.0,  11 March     2006]',/,  &
     &         ':::::   Calculate number of significant digits for floating point arithmetic')
!
 6003 FORMAT(/,24X,84('*'),/,24X,'*',23X,'EVALUATE FLOATING POINT ARITHMETIC',25X,'*',/,24X,84('*'),/,24X,'*',82X,'*',/,   &
     &       24X,'*         FLOATING POINT PROCESSOR HAS APPROXIMATELY',I3,' SIGNIFICANT DIGITS',T108,'*',/, &
     &       24X,'*',82X,'*',/,      &
     &       24X,'*               The default value of increment factor for numerical',T108,'*',   /, &
     &       24X,'*                   derivatives is <derivative_increment> =',1pE11.4,T108,'*',/, &
     &       24X,'*',82X,'*')
 6010 FORMAT(24X,'*               Default value for <derivative_increment> will be used',T108,'*')
!
 6011 FORMAT(24X,'*          User-specified value <derivative_increment> = ',1pe11.4,' will be used',T108,'*')
 6016 FORMAT(24X,'*',82X,'*'/24X,84('*')/)
!
 6024 FORMAT(' !!!!!!!!!!  WARNING  >>>>>>>>>> NUMBER OF SIGNIFICANT DIGITS IS MARGINAL;',   &
     &       ' EXPECT DETERIORATED CONVERGENCE BEHAVIOR',/)
 6025 FORMAT(' !!!!!!!!!!  WARNING  >>>>>>>>>> NUMBER OF SIGNIFICANT DIGITS IS INSUFFICIENT;',   &
     &       ' CONVERGENCE WILL BE POOR OR FAIL',/,   &
     &       ' WWWWWWWWWWWWWWWWWWWWWWWWWWWWW: CODE SHOULD BE RUN IN DOUBLE PRECISION!',/)
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of Significant_Digits
!
!
      RETURN
!
      END SUBROUTINE Significant_Digits
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE CPU_Timing_Routine(time)
!
         USE General_External_File_Units
!
      IMPLICIT NONE
!
      REAL(KIND = 8), INTENT(OUT) :: time
!
      INTEGER :: itime,irate
!
      LOGICAL :: First_call = .TRUE.
!
      SAVE First_call
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Begin CPU_Timing_Routine
!
!
      IF(First_call) THEN
         WRITE(VERS_Unit,6000)
         First_call = .FALSE.
      END IF
!
 6000 FORMAT(/,'CPU_Timing_Routine',T50,'[v1.0,  10 March     2006]',/,  &
     &         ':::::   CPU time computation: Uses FORTRAN95 intrinsic timing functions')
!
!
      call SYSTEM_CLOCK(COUNT=itime, COUNT_RATE=irate)
      TIME = dble(itime)/dble(irate)
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End CPU_Timing_Routine
!
!
      RETURN
!
      END SUBROUTINE CPU_Timing_Routine
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE Print_Input_Data
!
! ...... Modules to be used
!
         USE Basic_Parameters
         USE General_External_File_Units
         USE General_Control_Parameters
!
         USE Grid_Geometry
         USE Element_Attributes, ONLY: ElemMedia
         USE Solution_Matrix_Arrays, ONLY: X
!
         USE Geologic_Media_Parameters
         USE Sources_and_Sinks
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
! ... Integer variables
! -------
!
      INTEGER :: i,j,n
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of Print_Input_Data
!
!
      WRITE(VERS_Unit,6000)
!
! -------
! ... Print headings
! -------
!
      WRITE(*,100) ' '
      WRITE(*,150) title
      WRITE(*,1000)
      WRITE(*,1000)
!
! -------
! ... Print basic computational parameters
! -------
!
      WRITE(*,6100)
      WRITE(*,6102) Max_NumNRIterations,Max_NumTimeSteps,PRINT_frequency,Option_OutputAmount,SAVE_frequency,                       &
     &              derivative_increment,rel_convergence_crit,abs_convergence_crit,W_upstream,W_NRIteration,W_implicitness, &
     &              TimeOrigin,SimulationTimeEnd,InitialTimeStep,MaxTimeStep,Dt_reducer,gravity
!
!
 6100 FORMAT(/,T5,'                         P R O B L E M    S P E C I F I C A T I O N S',//)
 6102 FORMAT(T5,'                          COMPUTATIONAL/SIMULATION CONTROL PARAMETERS',/, &
     &       T5,'            Maximun number of Newton-Raphson iterations per time step,   Max_NumNRIterations  = ',i6,/,  &
     &       T5,'                                         Maximun number of time steps,   Max_NumTimeSteps     = ',i6,/,  &
     &       T5,'                           Output printing frequency (# of timesteps),   PRINT_frequency      = ',i6,/,  &
     &       T5,'                                              Output printing options,   Option_OutputAmount  = ',i6,/,  &
     &       T5,'                 Frequency of updating the SAVE file (# of timesteps),   SAVE_frequency       = ',i6,//, &
     &       T5,'       Perturbation increment for computing the numerical derivatives,   derivative_increment = ',1pe15.8,/,  &
     &       T5,'  Convergence criterion #1 for relative error of Newtonian iterations,   rel_convergence_crit = ',1pe15.8,/,  &
     &       T5,'  Convergence criterion #2 for absolute error of Newtonian iterations,   abs_convergence_crit = ',1pe15.8,/,  &
     &       T5,'                                             Upstream weighing factor,   W_upstream           = ',1pe15.8,/,  &
     &       T5,'                             Newton-Raphson iteration weighing factor,   W_NRIteration        = ',1pe15.8,/,  &
     &       T5,'                Weighing factor determining the level of implicitness,   W_implicitness       = ',1pe15.8,//, &
     &       T5,'                                        TIME AND TIME-STEP PARAMETERS',/, &
     &       T5,'                         Time at beginning of simulation period (sec),   TimeOrigin           = ',1pe15.8,/,  &
     &       T5,'                               Time at end of simulation period (sec),   SimulationTimeEnd    = ',1pe15.8,/,  &
     &       T5,'                                              Initial time step (sec),   InitialTimeStep      = ',1pe15.8,/,  &
     &       T5,'                                              Maximum time step (sec),   MaxTimeStep          = ',1pe15.8,/,  &
     &       T5,'                               Time step reduction factor in cutbacks,   Dt_reducer           = ',1pe15.8,/,  &
     &       T5,'                                      Acceleration of gravity (m/s^2),   gravity              = ',1pe15.8,//, &
     &       T6,64('='),/)
!
! -------
! ... Print user-specified Dt info
! -------
!
      IF(NumUserTimeSteps /= 0) THEN
         WRITE(*,450)
         WRITE(*,451) ((UserTimeStep(J+8*(N-1)),J=1,8), n=1,NumUserTimeSteps)
      END IF
!
! -------
! ... Print the initial values of the primary variables (to be applied uniformly)
! -------
!
      WRITE(*,6104)
      WRITE(*,6105) (default_initial_cond(n), n = 1,NumEqu)
!
!
 6104 FORMAT(T6,'U N I F O R M    P R I M A R Y    V A R I A B L E    V A L U E S',/, &
     &       T6,64('_'),//, &
     &       T6,'Primary variable # 1',4X,'Primary variable # 2',4X,'Primary variable # 3',  &
     &       4X,'Primary variable # 4',4X,'Primary variable # 5',/)
 6105 FORMAT(T6,6(1pe20.13,4x))
!
!
      WRITE(*,6500)
 6500 FORMAT(/)
!
!
      WRITE(*,1000)
      WRITE(*,1000)
!
      WRITE(*,700)
!
! -------
! ... Print the properties of the porous/fractured medium
! -------
!
      DO_NumRok: DO i=1,NumMedia
         WRITE(*,750)
         WRITE(*,800) i,MediumName(i),    &
     &                  media(i)%DensG,   &
     &                  media(i)%Poros,   &
     &                  media(i)%KThrW,   &
     &                  media(i)%SpcHt,   &
     &                  media(i)%Compr,   &
     &                  media(i)%Expan
         WRITE(*,810)
         WRITE(*,820) media(i)%Perm(1),media(i)%Perm(2),media(i)%Perm(3)
      END DO DO_NumRok
!
      WRITE(*,1000)
      WRITE(*,830)
      WRITE(*,840)
!
! -------
! ... Print the element-specific data
! -------
!
      DO_NumEle: DO i=1,NumElemTot
         WRITE(*,850) elem(i)%name,elem(i)%MatNum,elem(i)%vol
      END DO DO_NumEle
!
      WRITE(*,1000)
      WRITE(*,860)
      WRITE(*,870)
!
! -------
! ... Print the connection-specific data
! -------
!
      DO_NumConx: DO i=1,NumConx
         WRITE(*,880) conx(i)%name1, conx(i)%name2, conx(i)%ki, conx(i)%d1, conx(i)%d2, conx(i)%area, conx(i)%beta
      END DO DO_NumConx
!
! -------
! ... Print the source/sink data
! -------
!
      IF(NumSS <= 0 ) GO TO 9000
!
      WRITE(*,1000)
      WRITE(*,890)
!
!
!
      DO_NumSS: DO i=1,NumSS
!
! ...... Print flow rates of source/sink
!
         WRITE(*,895)
         WRITE(*,900) SS(i)%ElemName, SS(i)%name, SS(i)%rate, SS(i)%enth
!
! ...... Print comment for fixed rates or tabular data
!
         WRITE(*,896)
!
      END DO DO_NumSS
!
! -------
! ... Print initial conditions by element
! -------
!
 9000 WRITE(*,1000)
      WRITE(*,910)
      WRITE(*,920)
!
      DO_NumElem: DO i=1,NumElemTot
         WRITE(*,930) elem(i)%name, ElemMedia(i,current)%porosity, (X((i-1)*NumComPlus1+j),j=1,NumComPlus1)
      END DO DO_NumElem
!
! -------
! ... Print initial closing lines
! -------
!
      WRITE(*,1000)
      WRITE(*,940)
      WRITE(*,1000)
      WRITE(*,1000)
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 6000 FORMAT(/,'Print_Input_Data',T50,'[v1.0,  22 November  2006]',6X,/,   &
     &         ':::::   Provide printout of most data provided through the <FTSim> main input file')
!
  100 FORMAT(A1,/,'     <FTSim> INPUT DATA'/)
  150 FORMAT(5X,'PROBLEM TITLE: ',A120,/)
!
  400 FORMAT(5X,8(5X,1pE10.3))
  401 FORMAT(1X,4(4X,1pE10.4),9X,A5,3(4X,1pE10.4),4x,i5)
  450 FORMAT(/,5X,' VARIABLE TIME STEPS ARE PRESCRIBED')
  451 FORMAT(5X,8(5X,E10.4))
!
  700 FORMAT(/,5X,'ROCK PROPERTIES',/)
  750 FORMAT(5X,'DOMAIN     MAT        DENSITY        POROSITY     CONDUCTIVITY     HEAT CAP       COMPR          EXPAN')
!
  800 FORMAT(5X,I4,6X,A5,6(5X,1pE10.3),/)
  810 FORMAT(2X,'          PERM1          PERM2          PERM3')
  820 FORMAT(5X,3(5X,1pE10.3),/)
  830 FORMAT(/,5X,'ELEMENTS',/)
  840 FORMAT(5X,'        ELEMENT       MATERIAL         VOLUME')
  850 FORMAT(14X,A5,8X,I5,10X,1pE10.3)
!
  885 FORMAT(5X,'WELL ON DELIVERABILITY   ::::::::   OPEN IN',I3,' LAYERS   :::::::',///,   &
     &       5X,'   ELEMENT         SOURCE           PI             PWB            DEL(Z)',/)
  860 FORMAT(/,5X,'CONNECTIONS',/)
  870 FORMAT(5X,'          ELEM1          ELEM2           ISOT           DEL1           DEL2           AREA           BETA',/)
  880 FORMAT(5X,2(10X,A5),9X,I5,6X,4(5X,1pE10.3))
  890 FORMAT(/,5X,'GENERATION DATA')
  895 FORMAT('        ELEMENT         SOURCE           RATE           ENTHALPY'/)
  896 FORMAT(5X,'CONSTANT RATE IS SPECIFIED',/)
!
  900 FORMAT(5X,5X,A5,10X,A5,4X,4(5X,1pE10.3))
  910 FORMAT(/,5X,'INITIAL CONDITIONS',/)
  920 FORMAT(5X,'        ELEMENT       POROSITY         X1             X2             X3             X4             X5',/)
  930 FORMAT(5X,10X,A5,6(5X,1pE10.3))
  940 FORMAT(/,5X,'END OF INPUT DATA',/)
!
 1000 FORMAT(5X,117('*'))
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End of Print_Input_Data
!
!
      RETURN
!
!
      END SUBROUTINE Print_Input_Data
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
