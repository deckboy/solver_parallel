!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>                                                                     >
!>                                                                     >
!>        T_Inputs.f95: Code unit including the routines that          >
!>           read the inputs from the various input files              >
!>                                                                     >
!>                                                                     >
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! This segment includes the procedures involved in the reading of the
! general input files needed for FTSim simulations.  It does not include
! any procedure reading gas-related data (this is accomplished in the
! FTSim_Gas_Specifics.f95 segment).
!
!
      SUBROUTINE Read_Main_Input_File(Flag_MINC)
!
! ...... Modules to be used
!
         USE EOS_Routine_Selector
!
         USE EOS_Default_Parameters
         USE EOS_Parameters
!
         USE Input_Comments
!
         USE Basic_Parameters
         USE General_External_File_Units
         USE General_Control_Parameters
         USE Solver_Parameters
!
         USE Grid_Geometry
         USE Connection_Attributes
!
         USE Geologic_Media_Parameters
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*         ROUTINE FOR READING THE DATA IN THE MAIN INPUT FILE         *
!*             DIRECTLY FROM THE <FTsim> INPUT DATA BLOCK              *
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
      INTEGER :: icom = 0, ier = 0
!
      INTEGER :: N_DomInit
!
      INTEGER :: ico,i
!
! -------
! ... Character variables
! -------
!
      CHARACTER(LEN = 5)  :: KeyWord
      CHARACTER(LEN = 8)  :: LongKeyWord
      CHARACTER(LEN = 72) :: w72
      CHARACTER(LEN = 75) :: w75
!
! -------
! ... Logical variables
! -------
!
      LOGICAL, INTENT(OUT) :: Flag_MINC
!
! -------
! ... Namelists
! -------
!
      NAMELIST/Solver_Specifications/  MatrixSolver,          &
     &                                 Max_CGIterationRatio,  &
     &                                 CG_convergence_crit
!
! -------
! ... Saving variables
! -------
!
      SAVE icom
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of Read_Main_Input_File
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
      CALL Initializations_and_Defaults
!

      N_DomInit = 0       ! Number of domains for domain-wise initialization
!
      Flag_MINC = .FALSE. ! <Flag_MINC> = .FALSE. : A regular grid is involved; otherwise, a fractured medium grid is involved
!
!
!***********************************************************************
!*                                                                     *
!*                 READING AND SEARCHING FOR KEYWORDS                  *
!*                                                                     *
!***********************************************************************
!
!
 1000 READ(*, FMT = '(A8,A72)') LongKeyword,w72
!
!--------
! ... Check if the input involves FTSim input formats
!--------
!
      IF (LongKeyword(1:3) /= '>>>') THEN
         KeyWord = LongKeyword(1:5)
         w75     = LongKeyword(6:8)//w72
      ELSE
         KeyWord = LongKeyword(4:8)
         w75     = w72//'   '
      END IF
!
!
!
      CASE_KeyWord: SELECT CASE(KeyWord)
!
!
!***********************************************************************
!*                                                                     *
!*               KeyWord = 'ROCKS': READ ROCK PROPERTIES               *
!*                                                                     *
!***********************************************************************
!
!
      CASE ('MEDIA','ROCKS')
!
         CALL Read_Geologic_Media_Properties
!
         GO TO 1000
!
!
!***********************************************************************
!*                                                                     *
!*          KeyWord = 'DTIME': READ TIME DISCRETIZATION DATA           *
!*                                                                     *
!***********************************************************************
!
!
      CASE ('TIME_')
!
         CALL Read_Time_Discretization_Data
!
         GO TO 1000
!
!
!***********************************************************************
!*                                                                     *
!*          KeyWord = 'NRI_D': READ Newton-Raphson Iteration DATA      *
!*                                                                     *
!***********************************************************************
!
!
      CASE ('NRI_D')
!
         CALL Read_NR_Iteration_Data
!
         GO TO 1000
!
!
!***********************************************************************
!*                                                                     *
!*          KeyWord = 'OUTPU': READ Output Options                     *
!*                                                                     *
!***********************************************************************
!
!
      CASE ('OUTPU')
!
         CALL Read_Output_Options
!
         GO TO 1000
!
!
!***********************************************************************
!*                                                                     *
!*          KeyWord = 'COMPU': READ COMPUTATIONAL PARAMETERS           *
!*                                                                     *
!***********************************************************************
!
!
      CASE ('COMPU')
!
         CALL Read_Computational_Parameters
!
         GO TO 1000
!
!
!***********************************************************************
!*                                                                     *
!*          KeyWord = 'INITI': READ INITIAL CONDITIONS                 *
!*                                                                     *
!***********************************************************************
!
!
      CASE ('INITI')
!
         CALL Read_Initial_Conditions
!
         GO TO 1000
!
!
!***********************************************************************
!*                                                                     *
!*            KeyWord = 'ELEME': READ ELEMENT-RELATED DATA             *
!*                                                                     *
!***********************************************************************
!
!
      CASE ('ELEME')
!
         CALL Read_Element_Data
!
         GO TO 1000
!
!
!***********************************************************************
!*                                                                     *
!*          KeyWord = 'CONNE': READ CONNECTION-RELATED DATA            *
!*                                                                     *
!***********************************************************************
!
!
      CASE ('CONNE')
!
         CALL Read_Connection_Data
!
         GO TO 1000
!
!
!***********************************************************************
!*                                                                     *
!*         KeyWord = 'GENER': READ SINK/SOURCE-RELATED DATA            *
!*                                                                     *
!***********************************************************************
!
!
      CASE ('GENER')
!
         CALL Read_Generation_Data
!
         GO TO 1000
!
!
!***********************************************************************
!*                                                                     *
!*  KeyWord = 'SOLVR': READ SOLVER TYPE AND CORRESPONDING PARAMETERS   *
!*                                                                     *
!***********************************************************************
!
!
      CASE ('SOLVE')
!
         SolverBlok = .TRUE.
!
         READ(*, NML = Solver_Specifications, IOSTAT = ier)
!
! ...Stop if there is a problem reading the namelist
!
         IF(ier /= 0) THEN
            WRITE (*, FMT = 6800)
            STOP
         END IF
!
         GO TO 1000
!
!
!***********************************************************************
!*                                                                     *
!*            KeyWord = READ EOS-SPECIFIC DATA/PARAMETERS              *
!*                                                                     *
!***********************************************************************
!
!
      CASE (EOS_ID(1:5))
!
         CALL READ_EOS_Specific_Data
!
         GO TO 1000
!
!
!***********************************************************************
!*                                                                     *
!*            KeyWord = 'TIMES': READ TIMES FOR PRINT-OUTS             *
!*                                                                     *
!***********************************************************************
!
!
      CASE ('TIMES')
!
         CALL Read_Printout_Time_Data
!
         GO TO 1000
!
!
!***********************************************************************
!*                                                                     *
!*        KeyWord = 'ENDCY': READING OF INPUT FILE IS COMPLETE         *
!*                                                                     *
!***********************************************************************
!
!
      CASE ('ENDCY')
!
! ...... Initialization - CAREFUL! Array operations
!
         StateIndex = 0_1
!
! ----------
! ...... Print accumulated comments
! ----------
!
         IF_ComPrint: IF(icom /= 0) THEN
!
            ico = MIN(50,icom)
            WRITE(*,6010) ' ',icom
!
            DO i=1,ico
               WRITE(*,6011) comm(i)
            END DO
!
            WRITE(*,6012)
!
         END IF IF_ComPrint
!
!
!***********************************************************************
!*                                                                     *
!*      Default case: The lines that start with unknown keywords       *
!*                    are stored as comments (<= 50 lines)             *
!*                                                                     *
!***********************************************************************
!
!
      CASE DEFAULT
!
! ...... Storing lines with unknown leading keywords
!
         icom = icom+1
         IF(icom <= 50) comm(icom) = KeyWord//w75
!
         GO TO 1000
!
!
!
      END SELECT CASE_KeyWord
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 5004 FORMAT(I5,5X,7E10.4)
!
 5006 FORMAT(16I5)
 5008 FORMAT(i1,2x,a2,3x,a2,2(e10.4))
!
 5010 FORMAT(8E10.4)
!
 5012 FORMAT(A5,5x,a3)
 5022 FORMAT(A8,2x,a3)
!
 5014 FORMAT(A10,5x,a3)
 5024 FORMAT(A16,4x,a3)
!
!
!
 6000 FORMAT(/,'Read_Main_Input_File',T50,'[v1.0,  29 September 2007]', /, &
               ':::::   Read all data provided in the main FTSim input file')
!
 6004 FORMAT(' !!!!! WARNING !!!!! ==>  IE(1) =',I5,' in block <SELEC> is too large - Can read at most 64 records',/)
!
 6006 FORMAT(/,' Have read keyword "ENDFI" in file <INPUT> ==> Bypass flow simulation   ',62('*'),/,   &
     &         ' Optional  printout of input data (for Flag_print_input_data /= 0) is available')
!
 6010 FORMAT(A1,/,23X,86('*'),/,23X,'*',84X,'*',/,   &
     &       23X,'*',5X,'Comments in Input Data (have encountered',I4,' records with unknown keywords)',   &
     &       4X,'*',/,   &
     &       23X,'*',30X,'(Print up to 50 of them)',30X,'*',/,   &
     &       23X,'*',84X,'*'/23X,86('*'),/,23X,'*',84X,'*')
 6011 FORMAT(23X,'*  ',A80,'  *')
 6012 FORMAT(23X,'*',84X,'*'/23X,86('*'))
!
 6100 FORMAT(//,20('ERROR-'),//,T33,  &
     &             '       S I M U L A T I O N   A B O R T E D',/,                              &
     &         T27,'There is a problem reading one of the data block keywords <KeyWord> in the input file', &
     &       /,T32,'               CORRECT AND TRY AGAIN',                                      &
     &       //,20('ERROR-'))
!
 6800 FORMAT(//,20('ERROR-'),//,T33,  &
     &             '       S I M U L A T I O N   A B O R T E D',/,                                   &
     &         T27,'There is a problem reading namelist <Solver_Specifications> in the input file',  &
     &       /,T32,'               CORRECT AND TRY AGAIN',                                           &
     &       //,20('ERROR-'))
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of Read_Main_Input_File
!
!
      RETURN
!
      END SUBROUTINE Read_Main_Input_File
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE Initializations_and_Defaults
!
! ...... Modules to be used
!
         USE Basic_Parameters
         USE General_External_File_Units
         USE General_Control_Parameters
!
         USE Solver_Parameters, ONLY: SolverBlok
!
         USE Connection_Attributes
!
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*      ROUTINE FOR INITIALIZING IMPORTANT WIDELY USED PARAMETERS      *
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
      INTEGER :: n
!
#ifdef USE_TIMER
         real :: start, finish
#endif
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of Read_Main_Input_File
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
      MaxTimeStep   = 0.0d0            ! Maximum timestep size
!
      NumPrintTimes = 0                ! NumPrintTimes : # of print-out time data specified
!
! ----------
! ... Initialize computational parameters
! ----------
!
      MOP = 0                          ! WHOLE ARRAY operation
!
! ----------
! ... Initialize logical variables
! ----------
!
      NoFlowSimulation = .FALSE.       ! <NoFlowSimulation> = .FALSE. : Flow simulation will be conducted.
                                       !    Reset to .TRUE. if flow simulation is bypassed.
!
      NoVersion     = .FALSE.          ! <NoVersion> = .FALSE. : Version info is to be printed.
                                       !    Reset to .TRUE. if reset by the NOVER keyword.
!
      RestartSimulation = .FALSE.      ! <RestartSimulation> = .FALSE. : This is not a restart simulation.
                                       !    The initial conditions are not read from the renamed SAVE file
                                       !    from an earlier run.  Will be internally reset to .TRUE. if otherwise
!
      variable_porosity = .FALSE.      !  The porosity and permeability are assumed to remain constant
                                       !   despite changes in pressure and temperature.  Will be reset to true
                                       !   if geomechanics/geochemistry are considered, or if pore compressibilities
                                       !   and/or pore expansivities are nonzero
!
! ----------
! ... Initialize primary variable values - Array operations
! ----------
!
      default_initial_cond(1:NumComPlus1) = 0.0d0
      xx                                  = 0.0d0
!
! ----------
! ... Initialize flows: Whole array operation
! ----------
!
#ifdef USE_TIMER
         call cpu_time(start)
#endif

#ifdef USE_OMP
!$OMP WORKSHARE
#endif
      FORALL (n=1:Max_NumConx)
!
         ConxFlow(n)%rate(1:NumMobPhases)     = 0.0d0
         ConxFlow(n)%DarcyVel(1:NumMobPhases) = 0.0d0
         ConxFlow(n)%PoreVel(1:NumMobPhases)  = 0.0d0
!
         ConxFlow(n)%CompInPhase(1:NumCom,1:NumMobPhases) = 0.0d0
!
      END FORALL
!
#ifdef USE_OMP
!$OMP END WORKSHARE
#endif

#ifdef USE_TIMER
         call cpu_time(finish)

         write (*,*) __FILE__, ":", __LINE__, " time: ", finish-start
#endif
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 6000 FORMAT(/,'Initializations_and_Defaults',T50,'[v1.0,  30 September 2007]', /, &
               ':::::   Initialize important parameters and variables ')
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of <Initializations_and_Defaults>
!
!
      RETURN
!
      END SUBROUTINE Initializations_and_Defaults
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE Read_Geologic_Media_Data
!
! ...... Modules to be used
!
         USE Universal_Parameters
         USE Basic_Parameters
         USE General_External_File_Units
!
         USE Geologic_Media_Parameters
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*   ROUTINE FOR READING THE DATA DESCRIBING THE DOMAIN CONNECTIONS    *
!*             DIRECTLY FROM THE FTSim INPUT DATA BLOCK                *
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
      INTEGER :: NAD,mn,i,n_char
!
! -------
! ... Character variables
! -------
!
      CHARACTER(LEN = 5) :: MAT_name, w_temp
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of Read_Geologic_Media_Data
!
!
      WRITE(VERS_Unit,6000)
!
!
!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>                                                                     >
!>                         ROCK TYPE LOOP                              >
!>                                                                     >
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
      DO_NumRock: DO mn = 1,Max_NumMedia+1
!
!***********************************************************************
!*                                                                     *
!*             Checking the adequacy of the dimensioning               *
!*                                                                     *
!***********************************************************************
!
         IF_Limit: IF(mn == Max_NumMedia+1) THEN
!
             READ(*,5002) MAT_name      ! Reading the Max_NumMedia+1 record
!
! .......... If 'MAT_name' is not blank, print an error message & stop
!
             IF_Last: IF(MAT_name /= '     ') THEN
                WRITE(*,6001) Max_NumMedia
                STOP
!
! ......... Otherwise, the rock data read-in is completed
!
             ELSE
                NumMedia = Max_NumMedia  ! The number of media is set
                EXIT DO_NumRock
             END IF IF_Last
!
         END IF IF_Limit
!
!***********************************************************************
!*                                                                     *
!*              READING THE ROCK TYPE DATA: 1st Record                 *
!*                                                                     *
!***********************************************************************
!
         READ(*,5002) MediumName(mn),           &  ! Rock type name
     &                NAD,                      &  ! Number of additional lines with media properties
     &                media(mn)%DensG,          &  ! Rock bulk density (kg/m^3)
     &                media(mn)%Poros,          &  ! Porosity
     &               (media(mn)%Perm(i),i=1,3), &  ! Intrinsic permeability along principal axes (m^2)
     &                media(mn)%KThrW,          &  ! Thermal conductivity of water-saturated rock (W/m^2)
     &                media(mn)%SpcHt              ! Rock specific heat (J/kg/C)
!
! --------
! ... End of the rock records
! --------
!
         IF_RokEnd: IF(MediumName(mn) == '     ') THEN
!
! ......... Determine the number of media in the system
!
            NumMedia = mn-1
!
! ......... Otherwise, the rock data read-in is completed
!
            EXIT DO_NumRock
!
         END IF IF_RokEnd
!
!***********************************************************************
!*                                                                     *
!*             INITIALIZATION OF DATA IN SECOND RECORD                 *
!*                                                                     *
!***********************************************************************
!
!
         media(mn)%Compr     = 0.0D0
         media(mn)%Expan     = 0.0D0
         media(mn)%KThrD     = 0.0D0
         media(mn)%Tortu     = 0.0D0
         media(mn)%Klink     = 0.0D0
!
! ...... Ensuring that the rock name is left justified
!
         w_temp = '     '
         n_char = LEN_TRIM(ADJUSTL(MediumName(mn)))
!
         w_temp(1:n_char) = TRIM(ADJUSTL(MediumName(mn)))
         MediumName(mn)   = w_temp
!
! -----------
! ...... Account for the number of additional records for the cases of
! ...... (a) reference conditions or (b) seeds for permeability multipliers
! -----------
!
         IF(MediumName(mn) == 'REFCO' .OR. MediumName(mn) == 'SEED') THEN
            IF(NAD /= 0) NAD = 0
         END IF
!
!***********************************************************************
!*                                                                     *
!*              READING THE ROCK TYPE DATA: 2nd Record                 *
!*                                                                     *
!***********************************************************************
!
         IF_RokInLine2: IF(NAD > 0) THEN
!
            READ(*,5004) media(mn)%Compr,       &  ! Rock compressibility (1/Pa)
     &                   media(mn)%Expan,       &  ! Rock expansivity     (1/Pa)
     &                   media(mn)%KThrD,       &  ! Dry rock thermal conductivity (W/m^2)
     &                   media(mn)%Tortu,       &  ! Tortuosity
     &                   media(mn)%Klink           ! Klinkenberg parameter
         END IF IF_RokInLine2
!
! -----------
! ...... Assign default value to dry thermal conductivity if needed
! -----------
!
         IF(ABS(media(mn)%KThrD) < 1.0D-6) media(mn)%KThrD = media(mn)%KThrW
!
! ...... Printing the rock type number and name
!
         WRITE(*, FMT = 6002) mn, MediumName(mn)
!
! <<<
! <<< End of the ROCK TYPE LOOP
! <<<
!
      END DO DO_NumRock
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 5002 FORMAT(A5,I5,8E10.4)
 5004 FORMAT(10E10.4)
 5008 FORMAT(I5,5X,7E10.4)
 5018 FORMAT(I5,5X,7E20.13)
!
 6000 FORMAT(/,'Read_Geologic_Media_Data',T50,'[v1.0,   6 January   2007]',/,   &
     &         ':::::   Read the rock type data from the <ROCK> block of the FTSim input data file')
!
 6001 FORMAT(//,20('ERROR-'),//,T43,   &
     &             'S I M U L A T I O N   A B O R T E D',/,   &
     &         T14,'The number of input rock types = ',i4,' exceeds the maximum declared number <Max_NumMedia> = ',i4, &
     &       /,T50,'CORRECT AND TRY AGAIN',   &
     &       //,20('ERROR-'))
!
 6002 FORMAT(/,' Domain No:',I3,'     Material Name : ',A5)
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of Read_Geologic_Media_Data
!
!
      RETURN
!
      END SUBROUTINE Read_Geologic_Media_Data
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE Read_Geologic_Media_Properties
!
! ...... Modules to be used
!
         USE Universal_Parameters
         USE Basic_Parameters
         USE General_External_File_Units
!
         USE Geologic_Media_Parameters
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*   ROUTINE FOR READING THE DATA DESCRIBING THE DOMAIN CONNECTIONS    *
!*             DIRECTLY FROM THE FTSim INPUT DATA BLOCK                *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
      IMPLICIT NONE
!
! -------
! ... Integer Arrays
! -------
!
      INTEGER, DIMENSION(0:30) :: ierror = 99999
!
! -------
! ... Integer variables
! -------
!
      INTEGER :: mn, i, n_char, ier = 0
!
! -------
! ... Character variables
! -------
!
      CHARACTER(LEN = 5) :: MAT_name, w_temp
!
! ------
! ... Floating Point
! ------
      REAL(kind = 8) :: grain_density, porosity, permeability(3), sat_therm_conductivity, unsat_therm_conductivity,  &
      &                 specific_heat, compressibility, pore_expansivity, tortuosity, klinkenberg
!
! ------
! ... Namelist
! ------
      NAMELIST/Media_Properties/ grain_density,            &	    ![kg/m^3]
        &                        porosity,                 &
        &                        permeability,	           &        !Absolute perms along three principal axes
        &                        sat_therm_conductivity,   &	    ![W/m/degC], Formation heat conductivity under fully liquid-saturated conditions
        &                        unsat_therm_conductivity, &        ![W/m/degC], Formation heat conductivity under desaturated conditions
        &                        specific_heat,            &	    ![J/kg/degC]
        &                        compressibility,          &
        &                        pore_expansivity,         &
        &                        tortuosity,               &        !Tortuosity factor for binary diffusion
        &                        klinkenberg                        !Klinkenberg parameter for enhancing gas phase permeability
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of Read_Geologic_Media_Properties
!
!
      WRITE(VERS_Unit,6000)
!

      DO_NumRock: DO mn = 1, Max_NumMedia+1

         READ(*,'(A5)',IOSTAT = ierror(0)) MAT_name

         IF(mn == Max_NumMedia+1) THEN
         IF(MAT_name(1:3) == '<<<') THEN
            WRITE(*,6001) Max_NumMedia
            STOP
         ELSE
           NumMedia=Max_NumMedia
           EXIT DO_NumRock
         END IF
         END IF

         IF(MAT_name(1:3) == '<<<') THEN
            NumMedia = mn - 1
            RETURN
         ELSE
           WRITE(*,6002) mn, MAT_name
         END IF
!
!
!***********************************************************************
!*                                                                     *
!*                          INITIALIZATIONS                            *
!*                                                                     *
!***********************************************************************
!
	  grain_density 	        = 0.0d0
	  porosity		            = 0.0d0
	  permeability(1:3)	        = 0.0d0
	  sat_therm_conductivity	= 0.0d0
	  unsat_therm_conductivity  = 0.0d0
	  specific_heat		        = 0.0d0
	  compressibility           = 0.0d0
	  pore_expansivity          = 0.0d0
	  tortuosity                = 0.0d0
	  klinkenberg               = 0.0d0
!
      media(mn)%DensG           = 0.0d0
      media(mn)%Poros           = 0.0d0
      media(mn)%Perm(1:3)       = 0.0d0
      media(mn)%KThrW           = 0.0d0
      media(mn)%KThrD           = 0.0d0
      media(mn)%SpcHt           = 0.0d0
      media(mn)%Compr           = 0.0d0
      media(mn)%Expan           = 0.0d0
      media(mn)%Tortu           = 0.0d0
      media(mn)%Klink           = 0.0d0
!
!***********************************************************************
!*                                                                     *
!*        READING THE BASIC DATA DESCRIBING THE GAS PROPERTIES         *
!*                                                                     *
!***********************************************************************
!
      READ(*, NML = Media_Properties, IOSTAT = ier)    ! Reading the properties of the porous medium
!
!--------
! ... Stop if there is a problem reading the medium properties
!--------
!
      IF (ier /= 0) THEN
         WRITE (*, FMT = 6100)
         STOP
      END IF
!
!
!***********************************************************************
!*                                                                     *
!*      Assigning properties and parameters to the object <gas>        *
!*                                                                     *
!***********************************************************************
!
!
     MediumName(mn) = MAT_name
!
!
     IF(grain_density /= 0.0d0) THEN
        media(mn)%DensG = grain_density
        END IF
     IF(porosity /= 0.0d0) THEN
        media(mn)%Poros = porosity
        END IF
     IF(permeability(1) /= 0.0d0) THEN
        media(mn)%Perm(1) = permeability(1)
        END IF
     IF(permeability(2) /= 0.0d0) THEN
        media(mn)%Perm(2) = permeability(2)
        END IF
     IF(permeability(3) /= 0.0d0) THEN
        media(mn)%Perm(3) = permeability(3)
        END IF
     IF(sat_therm_conductivity /= 0.0d0) THEN
        media(mn)%KThrW = sat_therm_conductivity
        END IF
     IF(unsat_therm_conductivity /= 0.0d0) THEN
        media(mn)%KThrD = unsat_therm_conductivity
        END IF
     IF(specific_heat /= 0.0d0) THEN
        media(mn)%SpcHt = specific_heat
        END IF
     IF(compressibility /= 0.0d0) THEN
        media(mn)%Compr = compressibility
        END IF
     IF(pore_expansivity /= 0.0d0) THEN
        media(mn)%Expan = pore_expansivity
        END IF
     IF(tortuosity /= 0.0d0) THEN
        media(mn)%Tortu = tortuosity
        END IF
     IF(klinkenberg /= 0.0d0) THEN
        media(mn)%Klink = klinkenberg
        END IF
!
!
   END DO DO_NumRock
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
 5002 FORMAT(A5,I5,8E10.4)
 5004 FORMAT(10E10.4)
 5008 FORMAT(I5,5X,7E10.4)
 5018 FORMAT(I5,5X,7E20.13)
!
 6000 FORMAT(/,'Read_Geologic_Media_Properties',T50,'[v1.0,   1 June   2008]',/,   &
     &         ':::::   Read the rock type data from the <ROCK> block of the FTSim input data file')
!
 6001 FORMAT(//,20('ERROR-'),//,T43,   &
     &             'S I M U L A T I O N   A B O R T E D',/,   &
     &         T14,'The number of input rock types = ',i4,' exceeds the maximum declared number <Max_NumMedia> = ',i4, &
     &       /,T50,'CORRECT AND TRY AGAIN',   &
     &       //,20('ERROR-'))
!
 6002 FORMAT(/,' Domain No:',I3,'     Material Name : ',A5)
!
 6100 FORMAT(//,22('ERROR-'),//,T33,  &
     &             '       S I M U L A T I O N   A B O R T E D',/,                                             &
     &         T22,'There is a problem reading the namelist <Media_Properties> in the <DO_NumRock> data block',  &
     &       /,T32,'               CORRECT AND TRY AGAIN',                                                     &
     &       //,22('ERROR-'))
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of Read_Geologic_Media_Properties
!
!
      RETURN
!
      END SUBROUTINE Read_Geologic_Media_Properties
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE Read_Element_Data
!
! ...... Modules to be used
!
         USE Basic_Parameters
         USE General_External_File_Units
!
         USE Geologic_Media_Parameters, ONLY: MediumName
!
         USE Utility_Functions, ONLY: N_Character
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*     ROUTINE FOR READING THE DATA DESCRIBING THE DOMAIN ELEMENTS     *
!*             DIRECTLY FROM THE FTSim INPUT DATA BLOCK                *
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
      REAL(KIND = 8) :: VOLX, X, Y, Z
!
! -------
! ... Integer variables
! -------
!
      INTEGER :: n,i,m
!
! -------
! ... Character variables
! -------
!
      CHARACTER(LEN = 1) :: activity
      CHARACTER(LEN = 5) :: MA12, ElName5C
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of Read_Element_Data
!
!
      WRITE(VERS_Unit,6000)
!
! --------
! ... Printing headings
! --------
!
      WRITE(*,5002)
      REWIND (UNIT = MESH_Unit)
!
! --------
! ... Printing according to the number of characters in the element names
! --------
!
      WRITE(*, 6004)
      WRITE(MESH_Unit,6002)
!
!***********************************************************************
!*                                                                     *
!*                       READING THE ELEMENT DATA                      *
!*                                                                     *
!***********************************************************************
!
 1000 CONTINUE
!
! ...... Read data for 5-Character elements
!
         READ(*,5010) ElName5C,   &           ! 5-Character element name
     &                MA12,       &           ! Name of the corresponding rock type
     &                VOLX,       &           ! Element volume (m^3)
     &                X,          &           ! X-coordinate of element center (m)
     &                Y,          &           ! Y-coordinate of element center (m)
     &                Z,          &           ! Z-coordinate of element center (m)
     &                activity                ! Activity flag
!
! ... Reset the activity coefficient
!
      IF(VOLX <= 1.0d-10 .OR. VOLX > 1.0d10) activity = 'I'
!
! --------
! ... End of the element records (Check it and then write)
! --------
!
      IF_ElEnd: IF(ElName5C == '     ') THEN
                   WRITE(MESH_Unit,5002)
                   RETURN
                END IF IF_ElEnd
!
       WRITE(MESH_Unit,5010) ElName5C,M,VOLX,X,Y,Z,activity
!
! ......... Continue reading data
!
       GO TO 1000
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 5002 FORMAT('     ')
!
 5010 FORMAT(A5,10x,A5,E10.4,20x,3E10.4,1x,a1)
!
!
 6000 FORMAT(/,'Read_Element_Data',T50,'[v1.0,   6 January   2007]',/,   &
     &         ':::::   Read the element data from the <ELEME> block of the FTSim input data file')
!
 6002 FORMAT('ELEME')
 6004 FORMAT(' Write file <MESH> from input data')
 6006 FORMAT('ELEME',a4)
 6008 FORMAT(' Write file <MESH> from input data - Mesh data follow the ',a4,' option')
!
 6401 FORMAT(' Reference to unknown material <',A3,A2,'> at element <',A5,'> ==> Will ignore element')
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of Read_Element_Data
!
!
      RETURN
!
      END SUBROUTINE Read_Element_Data
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE Read_Connection_Data
!
! ...... Modules to be used
!
         USE Basic_Parameters
         USE General_External_File_Units
!
         USE Grid_Geometry, ONLY: conx
!
         USE Utility_Functions, ONLY: N_Character
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*   ROUTINE FOR READING THE DATA DESCRIBING THE DOMAIN CONNECTIONS    *
!*              DIRECTLY FROM THE FTsim INPUT DATA BLOCK               *
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
      REAL(KIND = 8) :: D1,D2,AREAX,BETAX
!
! -------
! ... Integer variables
! -------
!
      INTEGER :: NE1,NE2,ISOT
      INTEGER :: N,N1,N2,N_1,N_2
!
! -------
! ... Character variables
! -------
!
      CHARACTER(LEN = 3) :: EL1,EL2
!
      CHARACTER(LEN = 1) :: F1_Format,F2_Format
      CHARACTER(LEN = 5) :: Elm1_5,Elm2_5,Elem1_W5,Elem2_W5
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of Read_Connection_Data
!
!
      WRITE(VERS_Unit,6000)
!
! --------
! ... Printing headings and initializing
! --------
!
      WRITE(MESH_Unit,6002)
!
! ... Initializing the connection counter to zero
!
      NumConx  = 0
!
!***********************************************************************
!*                                                                     *
!*                     READING THE CONNECTION DATA                     *
!*                                                                     *
!***********************************************************************
!
 1000 CONTINUE
!
! ............... Read connection data for 5-Character elements
!
         READ(*,5004) Elm1_5, &  ! 5-Character name of the 1st element in the connection
     &                ELm2_5, &  ! 5-Character name of the 2nd element in the connection
     &                ISOT,   &  ! =1,2,3: Permeability between the elements given by the ISOT perm. component in ROCKS
     &                D1,     &  ! Distance of interface from the center of 1st element (m)
     &                D2,     &  ! Distance of interface from the center of 2nd element (m)
     &                AREAX,  &  ! Interface area (m^2)
     &                BETAX      ! = Cos(b), b: angle between vertical and line connecting the element centers
!
         EL1       = Elm1_5(1:3)
         F1_Format = Elm1_5(4:4)
!
         EL2       = Elm2_5(1:3)
         F2_Format = Elm2_5(4:4)
!
         READ(Elm1_5(4:5),'(I2)') NE1        ! Read NE1(target) from Elm1_5(4:5)(string)
         READ(Elm2_5(4:5),'(I2)') NE2        ! Read NE2(target) from Elm2_5(4:5)(string)
!
! --------
! ... End of the connection records - No info on connection elements
! --------
!
      IF_El1End: IF(EL1 == '   ' .AND. NE1 == 0) THEN
!
                    WRITE(MESH_Unit,5001)
!
! ......The connection data read-in is completed
!
                    RETURN
!
                 END IF IF_El1End
!
! --------
! ... End of the connection records - Info on connection elements
! --------
!
      IF_El2End: IF(EL1 == '+++') THEN
!
                    WRITE(MESH_Unit,5002)
!
! ......Read the connection-associated element numbers
!
                    READ(*,5010)          (conx(N)%n1,conx(N)%n2, N=1,NumConx)
                    WRITE(MESH_Unit,5010) (conx(N)%n1,conx(N)%n2, N=1,NumConx)
!
! ......The connection data read-in is completed
!
                    WRITE(MESH_Unit,5001)
                    RETURN
!
                 END IF IF_El2End
!
!***********************************************************************
!*                                                                     *
!*          WRITING THE CONNECTION DATA INTO THE ÒMESHÓ FILE           *
!*                                                                     *
!***********************************************************************
!
         NumConx = NumConx + 1    ! ...... Count connections
!
! ...... Determine element numbers in connection NumConx
!
         N1 = NE1
!
         N2 = NE2
!
! ...... Writing the element info into MESH
!
         IF(N1 < 10 .AND. F1_Format == ' ') THEN
           Elem1_W5 = EL1//' '//N_Character(N1)
         ELSE IF(N1 < 10 .AND. F1_Format == '0') THEN
            Elem1_W5 = EL1//'0'//N_Character(N1)
         ELSE
            N_1     = N1/10
            N_2     = MOD(N1,10)
            Elem1_W5 = EL1//N_Character(N_1)//N_Character(N_2)
         END IF
!
         IF(N2 < 10 .AND. F2_Format == ' ') THEN
            Elem2_W5 = EL2//' '//N_Character(N2)
         ELSE IF(N2 < 10 .AND. F2_Format == '0') THEN
            Elem2_W5 = EL2//'0'//N_Character(N2)
         ELSE
            N_1     = N2/10
            N_2     = MOD(N2,10)
            Elem2_W5 = EL2//N_Character(N_1)//N_Character(N_2)
         END IF
!
         WRITE(MESH_Unit,6004) Elem1_W5,Elem2_W5,ISOT,D1,D2,AREAX,BETAX
!
! --------
! ... Continue reading connection data
! --------
!
      GO TO 1000
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 5001 FORMAT('     ')
 5002 FORMAT('+++  ')
!
 5004 FORMAT(A5,A5,15X,I5,5E10.4)
!
 5010 FORMAT(16I6)
!
! 5010 FORMAT(16I5)
!
 6000 FORMAT(/,'Read_Connection_Data',T50,'[v1.0,   6 January   2007]',/, &
     &         ':::::   Read the connection data from the <CONNE> block of the FTSim input data file')
!
 6002 FORMAT('CONNE')
!
 6004 FORMAT(A5,A5,15X,I5,5E10.4)
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of Read_Connection_Data
!
!
      RETURN
!
      END SUBROUTINE Read_Connection_Data
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE Read_Time_Discretization_Data
!
! ...... Modules to be used
!
         USE EOS_Default_Parameters
!
         USE Basic_Parameters
         USE General_External_File_Units
         USE General_Control_Parameters
!
         USE Solver_Parameters, ONLY: MatrixSolver
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*       ROUTINE FOR READING THE TIME DISCRETIZATION DATA              *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
      IMPLICIT NONE
!
! -------
! ... Real variables
! -------
!
      REAL(KIND = 8) :: time_at_simulation_beginning, time_at_simulation_end, initial_timestep
      REAL(KIND = 8) :: maximum_timestep, Dt_reduction_factor
!
! -------
! ... Integer variables
! -------
!
      INTEGER :: i, ier = 0
      INTEGER :: Max_number_of_timesteps, num_of_user_defined_timesteps
!
! -------
! ... Character variables
! -------
!
      CHARACTER(LEN = 3) :: units_of_time, StateIndx
!
! -------
! ... Character variables
! -------
!
      REAL(KIND = 8), DIMENSION(100) :: user_defined_timestep_array
!
! -------
! ... Namelists
! -------
!
      NAMELIST/Time_Discretization_Data/  Max_number_of_timesteps,      &  !Maximum number of time steps
     &                                    time_at_simulation_beginning, &  !Time at the beginning of the simulation
     &                                    time_at_simulation_end,       &  !Time at the end of the simulation
     &                                    initial_timestep,             &  !Initial time step size (sec)
     &                                    maximum_timestep,             &  !Maximum time step size (sec)
     &                                    Dt_reduction_factor,          &  !Reduction factor for cutting Dt
     &                                    units_of_time,                &  ! = <min>, <hrs>, or <days>; Default is <sec>
     &                                    num_of_user_defined_timesteps    !Number of the user-specified timesteps
!
      NAMELIST/User_Specified_Timesteps/  UserTimeStep                     !Listing of the user-specified timesteps
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of Read_Computational_Parameters
!
      WRITE(VERS_Unit,6000)
!
!
!*************************
!...Initializations
!*************************
!
      max_number_of_timesteps        = 0
      time_at_simulation_beginning   = 0.0d0
      time_at_simulation_end         = 0.0d0
      initial_timestep               = 0.0d0
      maximum_timestep               = 0.0d0
      Dt_reduction_factor            = 4.0d0
      units_of_time                  = 'sec'
      num_of_user_defined_timesteps  = 0
      user_defined_timestep_array    = 0.0d0
!
      Max_NumTimeSteps     = 0
      TimeOrigin           = 0.0d0
      SimulationTimeEnd    = 0.0d0
      InitialTimeStep      = 0
      MaxTimeStep          = 0.0d0
      Dt_reducer           = 0.0d0
      NumUserTimeSteps     = 0
      UserTimeStep         = 0.0d0
!
! ------
! ...READ THE TIME DISCRETIZATION DATA
! ------
!
     READ(*, NML = Time_Discretization_Data, IOSTAT = ier)
!
! ......Stop if there is a problem reading the medium properties
!
     IF(ier /= 0) THEN
        WRITE(*, FMT = 6100)
        STOP
     END IF
!
! ......Assign properties
!
      Max_NumTimeSteps   = max_number_of_timesteps
      NumUserTimeSteps   = num_of_user_defined_timesteps
      TimeOrigin         = time_at_simulation_beginning
      SimulationTimeEnd  = time_at_simulation_end
      Dt_reducer         = Dt_reduction_factor
      InitialTimeStep    = initial_timestep
      MaxTimeStep        = maximum_timestep
!
! ......Assign default values
!
      IF(MaxTimeStep == 0.0d0)   MaxTimeStep      = SimulationTimeEnd
      IF(PRINT_frequency  == 0)  PRINT_frequency  = Max_NumTimeSteps
      IF(Dt_reducer  <= 0.0d0)   Dt_reducer       = 4.0d0
!
! ......Convert all time-units to 'seconds'
!
       IF(units_of_time == 'min') THEN
         TimeOrigin         = TimeOrigin * 60.0d0
         SimulationTimeEnd  = SimulationTimeEnd * 60.0d0
         InitialTimeStep    = InitialTimeStep * 60.0d0
         MaxTimeStep        = MaxTimeStep * 60.0d0
       END IF
!
       IF(units_of_time == 'hrs') THEN
         TimeOrigin         = TimeOrigin * 3600.0d0
         SimulationTimeEnd  = SimulationTimeEnd * 3600.0d0
         InitialTimeStep    = InitialTimeStep * 3600.0d0
         MaxTimeStep        = MaxTimeStep * 3600.0d0
       END IF
!
       IF(units_of_time == 'day') THEN
         TimeOrigin         = TimeOrigin * 86400.0d0
         SimulationTimeEnd  = SimulationTimeEnd * 86400.0d0
         InitialTimeStep    = InitialTimeStep * 86400.0d0
         MaxTimeStep        = MaxTimeStep * 86400.0d0
       END IF
!
       Tot_NumTimeSteps     = 0
       NumNRIterations      = 0
       Tot_NumNRIterations  = 0
       TimeShift            = TimeOrigin
!
! ......If there are no user-defined time-steps, we are done
!
       IF(num_of_user_defined_timesteps <= 0) GO TO 1000
!
!
! ...READ USER SPECIFIED TIMESTEPS
!
          READ(*, NML = User_Specified_Timesteps, IOSTAT = ier)
!
! ...Stop if there is a problem reading the time metrics
!
          IF(ier /= 0) THEN
             WRITE (*, FMT = 6200)
             STOP
          END IF
!
! ... Check if the user-proved time-steps are valid
!
       IF( ANY(UserTimeStep(1:NumUserTimeSteps) < 1.0d-16) ) THEN
          WRITE(*, FMT = 6150)
          STOP
       END IF
!
! ... Convert all times to seconds
!
       IF(units_of_time == 'min') THEN
          UserTimeStep = UserTimeStep * 60.0d0
       ELSE IF (units_of_time == 'hrs') THEN
          UserTimeStep = UserTimeStep * 3600.0d0
       ELSE IF( units_of_time == 'DAY') THEN
          UserTimeStep = UserTimeStep * 86400.0d0
       END IF
!
! ... Initialization
!
      InitialTimeStep = UserTimeStep(1)
!
!***********************************************************************
!*                                                                     *
!*                      END OF THE DATA BLOCK                          *
!*                                                                     *
!***********************************************************************
!
 1000 READ(*, FMT = '(A3)', IOSTAT = ier) DataBlockEndSpecifier
!
!-----------
! ... Stop if there is a problem reading the data block end specifies
!-----------
!
      IF( (ier == 0) .AND. (DataBlockEndSpecifier == '<<<') ) THEN
         RETURN
      ELSE
         WRITE (*, FMT = 6500)
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
 5002 FORMAT(2I2,3I4,24I1,35x,i5)
 5004 FORMAT(4E10.4,A5,5X,3E10.4)
!
 5006 FORMAT(8E10.4)
 5010 FORMAT(7E10.4,2x,a3,5x,3E10.4)
 5020 FORMAT(6E20.13)
!
 6000 FORMAT(/,'Read_Time_Discretization_Data',T50,'[v1.0,  01 June 2008]',/,   &
     &         ':::::   Read the time discretization data from the TIME_DISCRETIZATION_DATA block of the FTSim input data file')
!
 6100 FORMAT(//,22('ERROR-'),//,T33,  &
     &             '       S I M U L A T I O N   A B O R T E D',/,                         &
     &         T15,'There is a problem reading the namelist <Time_Discretization_Data> ',  &
     &             'in the <Read_Time_Discretization_Data> data block',                    &
     &       /,T32,'               CORRECT AND TRY AGAIN',                                 &
     &       //,22('ERROR-'))
!
 6150 FORMAT(//,22('ERROR-'),//,T33,  &
     &             '       S I M U L A T I O N   A B O R T E D',/,                                             &
     &         T22,'There is a problem reading the user-specified time-steps <UserTimeStep> <1.0E-16 sec',/,   &
     &       /,T32,'               CORRECT AND TRY AGAIN',                                                     &
     &       //,22('ERROR-'))
!
 6200 FORMAT(//,22('ERROR-'),//,T33,  &
     &             '       S I M U L A T I O N   A B O R T E D',/,                                     &
     &         T22,'There is a problem reading the namelist <User_Specified_Timesteps> in subroutine   &
          &              Read_Time_Discretization_Data ',                                              &
     &       /,T32,'               CORRECT AND TRY AGAIN',                                             &
     &       //,22('ERROR-'))
!
 6500 FORMAT(//,22('ERROR-'),//,T33,  &
     &             '       S I M U L A T I O N   A B O R T E D',/,                                     &
     &         T04,'There is a problem reading the data block end specifier <DataBlockEndSpecifier> ', &
     &             'in the <Read_Time_Discretization_Data> block',                                     &
     &       /,T32,'               CORRECT AND TRY AGAIN',                                             &
     &       //,22('ERROR-'))
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End of Read_Time_Discretization_Data
!
!
      RETURN
!
      END SUBROUTINE Read_Time_Discretization_Data
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE Read_NR_Iteration_Data
!
! ...... Modules to be used
!
         USE EOS_Default_Parameters
!
         USE Basic_Parameters
         USE General_External_File_Units
         USE General_Control_Parameters
!
         USE Solver_Parameters, ONLY: MatrixSolver
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*       ROUTINE FOR READING THE MAIN NR Iteration Data                *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
      IMPLICIT NONE
!
! -------
! ... Real variables
! -------
!
      REAL(KIND = 8) :: relative_convergence_criterion, absolute_convergence_criterion
      REAL(KIND = 8) :: NR_solution_weighing_factor, criterion_for_Dt_doubling
      REAL(KIND = 8) :: W_NRSolution
!
! -------
! ... Integer variables
! -------
!
      INTEGER :: ier = 0
      INTEGER :: Max_number_of_NR_iterations
!
! -------
! ... Character variables
! -------
!
      CHARACTER(LEN = 3) :: StateIndx
!
! -------
! ... Namelists
! -------
!
      NAMELIST/NR_Iteration_Data/  Max_number_of_NR_iterations,     &  !Maximum number of Newtonian iterations
      &                            criterion_for_Dt_doubling,       &  !If NR convergence attained in fewer iterations than this number, then Dt is doubled on the next timestep
      &                            relative_convergence_criterion,  &  !Convergence criterion #1 for relative error of newtonian iterations
      &                            absolute_convergence_criterion,  &  !Convergence criterion #2 for absolute error of Newtonian iterations
      &                            NR_solution_weighing_factor         !Weighing factor of most recent Newton-Raphson solution in updating
                                                                       !    primary variables
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of Read_Computational_Parameters
!
      WRITE(VERS_Unit,6000)
!
!***********************************************************************
!*                                                                     *
!*       READING THE Newton-Raphson Iteration Data                     *
!*                                                                     *
!***********************************************************************
!
! ...Initializations
!
      Max_number_of_NR_iterations     = 0
      criterion_for_Dt_doubling       = 0
      relative_convergence_criterion  = 1.0d-5
      absolute_convergence_criterion  = 1.0d0
      NR_solution_weighing_factor     = 0.0d0
!
      Max_NumNRIterations             = 0.0d0
      DoublingDtCriterion             = 0.0d0
      rel_convergence_crit            = 0.0d0
      abs_convergence_crit            = 0.0d0
      W_NRSolution                    = 0.0d0
!
!***********************************************************************
!*                                                                     *
!*       READING THE NR Iteration Data                                 *
!*
!***********************************************************************
!
     READ(*, NML = NR_Iteration_Data, IOSTAT = ier)
!
     IF(ier /= 0) THEN
        WRITE(*, FMT = 6100)
        STOP
     END IF
!
! ...Relationships to code-used variables:
!
      Max_NumNRIterations   = Max_number_of_NR_iterations
      rel_convergence_crit  = relative_convergence_criterion
      abs_convergence_crit  = absolute_convergence_criterion
      W_NRSolution          = NR_solution_weighing_factor
!
! --------
! ... Assign default values
! --------
!
      IF(Max_NumNRIterations  == 0)      Max_NumNRIterations  = 8
      IF(rel_convergence_crit <= 0.0d0)  rel_convergence_crit = 1.0d-5
      IF(abs_convergence_crit <= 0.0d0)  abs_convergence_crit = 1.0d0
      IF(W_NRIteration        == 0.0d0)  W_NRIteration        = 1.0d0
      IF(DoublingDtCriterion  <= 0)      DoublingDtCriterion  = 4
!
!***********************************************************************
!*                                                                     *
!*                      END OF THE DATA BLOCK                          *
!*                                                                     *
!***********************************************************************
!
 1000 READ(*, FMT = '(A3)', IOSTAT = ier) DataBlockEndSpecifier
!
!-----------
! ... Stop if there is a problem reading the data block end specifies
!-----------
!
      IF( (ier == 0) .AND. (DataBlockEndSpecifier == '<<<') ) THEN
         RETURN
      ELSE
         WRITE (*, FMT = 6500)
         STOP
      END IF
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 5002 FORMAT(2I2,3I4,24I1,35x,i5)
 5004 FORMAT(4E10.4,A5,5X,3E10.4)
!
 5006 FORMAT(8E10.4)
 5010 FORMAT(7E10.4,2x,a3,5x,3E10.4)
 5020 FORMAT(6E20.13)
!
 6000 FORMAT(/,'Read_Computational_Parameters',T50,'[v1.0,  21 September 2007]',/,   &
     &         ':::::   Read the computational parameter data from the NRI block of the FTSim input data file')
!
 6100 FORMAT(//,20('ERROR-'),//,T43,   &
     &             'S I M U L A T I O N   A B O R T E D',/,   &
     &         T25,'There is a problem reading the namelist <NR_Iteration_Data> in subroutine <Read_NR_Iteration_Data>',/,   &
     &       //,T50,'CORRECT AND TRY AGAIN',   &
     &       //,20('ERROR-'))
!
 6101 FORMAT(/,' FILE <MaxP_Time_Series>  EXISTS ==> OPEN AS AN OLD FILE')
 6102 FORMAT(/,' FILE <MaxP_Time_Series> DOES NOT EXIST ==> OPEN AS A NEW FILE')
!
 6500 FORMAT(//,22('ERROR-'),//,T33,  &
     &             '       S I M U L A T I O N   A B O R T E D',/,                                     &
     &         T04,'There is a problem reading the data block end specifier <DataBlockEndSpecifier> ', &
     &             'in the <Read_NR_Iteration_Data> block',                                            &
     &       /,T32,'               CORRECT AND TRY AGAIN',                                             &
     &       //,22('ERROR-'))
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End of Read_NR_Iteration_Data
!
!
      RETURN
!
      END SUBROUTINE Read_NR_Iteration_Data
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
      SUBROUTINE Read_Output_Options
!
! ...... Modules to be used
!
         USE EOS_Default_Parameters
!
         USE Basic_Parameters
         USE General_External_File_Units
         USE General_Control_Parameters
!
         USE Solver_Parameters, ONLY: MatrixSolver
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*       ROUTINE FOR READING THE MAIN OUTPUT OPTIONS                   *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
      IMPLICIT NONE
!
!
! -------
! ... Integer variables
! -------
!
      INTEGER :: ier
      INTEGER :: Option_amount_of_output
      INTEGER :: output_frequency
      INTEGER :: frequency_of_storing_SAVE_file
      INTEGER :: number_specified_output_times
      INTEGER :: number_output_times
!
! -------
! ... CHARACTER variables
! -------
!
      CHARACTER(LEN = 3) :: StateIndx
      CHARACTER(LEN = 3) :: units_of_time
      CHARACTER(LEN = 4) :: Option_format_of_output
      CHARACTER(LEN = 5) :: name_of_tracked_element
!
! -------
! ... Logical variables
! -------
!
      LOGICAL :: Flag_print_input_data
!
! -------
! ... Namelists
! -------
!
      NAMELIST/Output_Options/  Option_amount_of_output,        &  !Option determining the amount of output
      &                         Option_format_of_output,        &  != Determining the format of the output (standard, plotting-package, or both)
      &                         Option_Print_NRIterationInfo,   &  != MOP(1): Determing the printout amount after each NR iteration
      &                         Option_Print_CyclingInfo,       &  != MOP(2): Determing the printout amount in the simulation cycling subroutine
      &                         Option_Print_JacobianInfo,      &  != MOP(3): Determing the printout amount in the subroutine setting up the Jacobian
      &                         Option_Print_SourceSinkInfo,    &  != MOP(4): Determing the printout amount in the source/sink subroutine
      &                         Option_Print_EOSInfo,           &  != MOP(5): Determing the printout amount in the equation-of-state subroutine
      &                         Option_Print_SolverInfo,        &  != MOP(6): Determing the printout amount in the linear equation solvers
      &                         Flag_print_input_data,          &  != MOP(7): Flag determining whether the inputs will be printed in the output file
      &                         PRINT_frequency,                &  !Print-outs every <output_frequency> time steps
      &                         frequency_of_storing_SAVE_file, &  !Frequency of storing the SAVE file
      &                         name_of_tracked_element,        &  !Name of tracked element
      &                         number_specified_output_times      !Number of provided rpint-out times (NumPrintTimes <= 100)
!
      NAMELIST/Specified_Output_Times/ PrintTime,               &  !Listing of the user-specified times...
      &                                units_of_time               ! = <min>, <hrs>, or <days>; Default is <sec>
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of Read_Computational_Parameters
!
      WRITE(VERS_Unit,6000)
!
!***********************************************************************
!*                                                                     *
!*              READING THE OPTIONAL OUTPUT                            *
!*                                                                     *
!***********************************************************************
!
!
!...Initializations
!
      Option_amount_of_output        = 0
      Option_format_of_output        = 'STAN'
      Option_Print_NRIterationInfo   = 0
      Option_Print_CyclingInfo       = 0
      Option_Print_JacobianInfo      = 0
      Option_Print_SourceSinkInfo    = 0
      Option_Print_EOSInfo           = 0
      Option_Print_SolverInfo        = 0
      Flag_print_input_data          = .FALSE.
      PRINT_frequency                = 0
      frequency_of_storing_SAVE_file = 0
      name_of_tracked_element        = '     '
      number_specified_output_times  = 0
!
! ....Initialization of variables declared in <FTSim_Allocate_Memory.f90>
!
      Option_OutputAmount            = 0         !Option_OutputAmount
      Option_OutputFormat            = 'STAN'    !Currently parameter Plot_Elem_Unit
      Flag_PrintInputs               = .FALSE.
      SAVE_frequency                 = 0
      NumPrintTimes                  = 0
      TrackElemName                  = '     '
!
      READ(*, NML = Output_Options, IOSTAT = ier)
!
! ... Stop if there is a problem reading the namelist
!
      IF(ier /= 0) THEN
         WRITE(*, FMT = 6105)
         STOP
      END IF
!
!...Relationships to code-used variables:
!
      Option_OutputAmount   = Option_amount_of_output
      Option_OutputFormat   = Option_format_of_output
!
      Flag_PrintInputs      = Flag_print_input_data
      SAVE_frequency        = frequency_of_storing_SAVE_file
      NumPrintTimes         = number_specified_output_times
      TrackElemName         = name_of_tracked_element
!
! --------
! ... Assign default values (for some)
! --------
!
      IF(Option_amount_of_output <= 0 .OR. Option_amount_of_output > 3) Option_OutputAmount = 1
!
      IF(Option_OutputFormat /= 'PLOT' .AND. Option_OutputFormat /= 'BOTH'       &
    &    .AND. Option_OutputFormat /= 'Plot' .AND. Option_OutputFormat /= 'Both' &
    &    .AND. Option_OutputFormat /= 'plot' .AND. Option_OutputFormat /= 'both') Option_OutputFormat = 'STAN'
!
      IF(PRINT_frequency <= 0) PRINT_frequency = 10000000
!
!... Read user-specified output times
!
      IF(NumPrintTimes <= 0) GO TO 1000
!
!***********************************************************************
!*                                                                     *
!*                 Read user-specified output times                    *
!*                                                                     *
!***********************************************************************
!
      READ(*, NML = Specified_Output_Times, IOSTAT = ier)
!
!--------
! ...... Stop if there is a problem reading the sim time metrics
!--------
!
      IF(ier /= 0) THEN
         WRITE (*, FMT = 6200)
         STOP
      END IF
!
!***********************************************************************
!*                                                                     *
!*                      END OF THE DATA BLOCK                          *
!*                                                                     *
!***********************************************************************
!
 1000 READ(*, FMT = '(A3)', IOSTAT = ier) DataBlockEndSpecifier
!
!-----------
! ... Stop if there is a problem reading the data block end specifies
!-----------
!
      IF( (ier == 0) .AND. (DataBlockEndSpecifier == '<<<') ) THEN
         RETURN
      ELSE
         WRITE (*, FMT = 6500)
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
 5002 FORMAT(2I2,3I4,24I1,35x,i5)
 5004 FORMAT(4E10.4,A5,5X,3E10.4)
!
 5006 FORMAT(8E10.4)
 5010 FORMAT(7E10.4,2x,a3,5x,3E10.4)
 5020 FORMAT(6E20.13)
!
 6000 FORMAT(/,'Read_Output_Options',T50,'[v1.0,  21 September 2007]',/,   &
     &         ':::::   Read the data from the OUTPUT_Options block of the FTSim input data file')
!
 6101 FORMAT(/,' FILE <MaxP_Time_Series>  EXISTS ==> OPEN AS AN OLD FILE')
 6102 FORMAT(/,' FILE <MaxP_Time_Series> DOES NOT EXIST ==> OPEN AS A NEW FILE')
!
 6105 FORMAT(//,22('ERROR-'),//,T33,  &
     &             '       S I M U L A T I O N   A B O R T E D',/,               &
     &         T15,'There is a problem reading the namelist <Output_Options> ',  &
     &             'in the <Read_Output_Options> data block',                    &
     &       /,T32,'               CORRECT AND TRY AGAIN',                       &
     &       //,22('ERROR-'))
!
 6200 FORMAT(//,22('ERROR-'),//,T33,  &
     &             '       S I M U L A T I O N   A B O R T E D',/,                                            &
     &         T22,'There is a problem reading the namelist <Output_Options> in subroutine Read_Output_Options ',   &
     &       /,T32,'               CORRECT AND TRY AGAIN',                                                    &
     &       //,22('ERROR-'))
!
 6500 FORMAT(//,22('ERROR-'),//,T33,  &
     &             '       S I M U L A T I O N   A B O R T E D',/,                                     &
     &         T04,'There is a problem reading the data block end specifier <DataBlockEndSpecifier> ', &
     &             'in the <Read_Output_Options> block',                                            &
     &       /,T32,'               CORRECT AND TRY AGAIN',                                             &
     &       //,22('ERROR-'))
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End of Read_Output_Options
!
!
      RETURN
!
      END SUBROUTINE Read_Output_Options
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE Read_Computational_Parameters
!
! ...... Modules to be used
!
         USE EOS_Default_Parameters
!
         USE Basic_Parameters
         USE General_External_File_Units
         USE General_Control_Parameters
!
         USE Solver_Parameters, ONLY: MatrixSolver
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*       ROUTINE FOR READING THE MAIN COMPUTATIONAL PARAMETERS         *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
      IMPLICIT NONE
!
! -------
! ... Real variables
! -------
!
      REAL(KIND = 8) :: acceleration_of_gravity, upstream_weighting_factor, implicitness_weighting_factor
!
! -------
! ... Integer variables
! -------
!
      INTEGER :: ier
      INTEGER :: Option_thermal_conductivity, Option_upstream_weighting, Option_interface_density
!
! -------
! ... CHARACTER variables
! -------
!
      CHARACTER(LEN = 3) :: line_feed
!
! -------
! ... LOGICAL variables
! -------
!
      LOGICAL :: Flag_check_initial_conditions, Flag_analytical_heat_exchange
!
! -------
! ... Namelists
! -------
!
      NAMELIST/Computational_Parameters/  upstream_weighting_factor,       &  !Upstream weighing factor
      &                                   implicitness_weighting_factor,   &  !Weighing factor determining the level of implicitness
      &                                   derivative_increment,            &  !Increment for estimation of numerical derivative
      &                                   acceleration_of_gravity,         &  !Acceleration of gravity
      &                                   Option_thermal_conductivity,     &  !Option describing the method for computing the thermal cond. of composite systems
      &                                   Option_upstream_weighting,       &  !Option describing the method of upstream weighting of mobilities and permeabilities
      &                                   Option_interface_density,        &  !Option describing the method of estimation of the interface density
      &                                   Flag_check_initial_conditions       !Flag for checking the accuracy/correctness/feasibility of
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of Read_Computational_Parameters
!
      WRITE(VERS_Unit,6000)
!
!
!...Initializations
!
      upstream_weighting_factor       = 0.0d0
      implicitness_weighting_factor   = 0.0d0
      derivative_increment            = 0.0d0
      acceleration_of_gravity         = 0.0d0
      Option_thermal_conductivity     = 0
      Option_upstream_weighting       = 0
      Option_interface_density        = 0
      Flag_check_initial_conditions   = .TRUE.
!
! ....Initialization of variables declared in <FTSim_Allocate_Memory.f90>
!
      W_upstream                     = 1.0d0
      W_implicitness                 = 1.0d0
      gravity                        = 0.0d0
      Option_ThermalConductivity     = 0
      Option_UpstreamWeighting       = 0
      Option_InterfaceDensity        = 0
      Flag_CheckInitialConditions    = .TRUE.
!
!--------
! ...... Read simulation runtime metrics
!--------
!
      READ(*, NML = Computational_Parameters, IOSTAT = ier)
!
      IF(ier /= 0) THEN
         WRITE(*, FMT = 6100)
         STOP
      END IF
!
!--------
! ...... Assigning parameter values
!--------
!
      W_upstream                     = upstream_weighting_factor
      W_implicitness                 = implicitness_weighting_factor
      gravity                        = acceleration_of_gravity
      Option_ThermalConductivity     = Option_thermal_conductivity
      Option_UpstreamWeighting       = Option_upstream_weighting
      Option_InterfaceDensity        = Option_InterfaceDensity
      Flag_CheckInitialConditions    = Flag_check_initial_conditions
!
! --------
! ... Assign default values (for some)
! --------
!
      IF(gravity < 0.0d0)                                     gravity        = 0.0d0
      IF(W_implicitness < 0.0d0 .OR. W_implicitness > 1.0d0)  W_implicitness = 1.0d0
      IF(W_upstream < 0.0d0 .OR. W_upstream > 1.0d0)          W_upstream     = 1.0d0
!
!***********************************************************************
!*                                                                     *
!*                      END OF THE DATA BLOCK                          *
!*                                                                     *
!***********************************************************************
!
 1000 READ(*, FMT = '(A3)', IOSTAT = ier) DataBlockEndSpecifier
!
!-----------
! ... Stop if there is a problem reading the data block end specifies
!-----------
!
      IF( (ier == 0) .AND. (DataBlockEndSpecifier == '<<<') ) THEN
         RETURN

      ELSE
         WRITE (*, FMT = 6500)
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
 5002 FORMAT(2I2,3I4,24I1,35x,i5)
 5004 FORMAT(4E10.4,A5,5X,3E10.4)
!
 5006 FORMAT(8E10.4)
 5010 FORMAT(7E10.4,2x,a3,5x,3E10.4)
 5020 FORMAT(6E20.13)
!
 6000 FORMAT(/,'Read_Computational_Parameters',T50,'[v1.0,  21 September 2007]',/,   &
     &         ':::::   Read the computational parameter data from the PARAM block of the FTSim input data file')
!
 6100 FORMAT(//,22('ERROR-'),//,T33,  &
     &             '       S I M U L A T I O N   A B O R T E D',/,                                            &
     &         T22,'There is a problem reading the namelist <Computational_Parameters> ',  &
     &             'in subroutine <Read_Computational_Parameters>', &
     &       /,T32,'               CORRECT AND TRY AGAIN',                                                    &
     &       //,22('ERROR-'))
!
 6500 FORMAT(//,22('ERROR-'),//,T33,  &
     &             '       S I M U L A T I O N   A B O R T E D',/,                                     &
     &         T04,'There is a problem reading the data block end specifier <DataBlockEndSpecifier> ', &
     &             'in the <Read_Computational_Parameters> block',                                            &
     &       /,T32,'               CORRECT AND TRY AGAIN',                                             &
     &       //,22('ERROR-'))
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End of Read_Computational_Parameters
!
!
      RETURN
!
      END SUBROUTINE Read_Computational_Parameters
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
      SUBROUTINE Read_Generation_Data
!
! ...... Modules to be used
!
         USE Basic_Parameters
         USE General_External_File_Units
!
         USE Sources_and_Sinks
!
         USE Utility_Functions, ONLY: N_Character
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*     ROUTINE FOR READING THE GENERATION (SOURCE AND SINK) DATA       *
!*             DIRECTLY FROM THE FTSim INPUT DATA BLOCK                *
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
      REAL(KIND = 8) :: GX, EX
!
! -------
! ... Integer variables
! -------
!
      INTEGER :: NE, NS, n, N_1, N_2
!
! -------
! ... Character variables
! -------
!
      CHARACTER(LEN = 3) :: EL,SL
      CHARACTER(LEN = 4) :: SS_Type
!
      CHARACTER(LEN = 1) :: F_Format,S_Format
      CHARACTER(LEN = 5) :: ElName5C, Elem_W5, SS_name
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of Read_Generation_Data
!
!
      WRITE(VERS_Unit,6000)
!
! --------
! ... Printing headings and attaching GENER file
! --------
!
      WRITE(*,5001)
      REWIND (UNIT = GENER_Unit)
!
      WRITE(*, 6002)
      WRITE(GENER_Unit,6004)
!
! ... Seeting the source/sink counter to zero
!
      NumSS = 0
!
!***********************************************************************
!*                                                                     *
!*                    READING THE SOURCE/SINK DATA                     *
!*                                                                     *
!***********************************************************************
!
 1000 CONTINUE
!
! ...... Read data for 5-Character elements
!
         READ(*,5004) ElName5C,       &  ! 5-Character name of the source/sink element
     &                SS_name,        &  ! 5-Character name of the source/sink
     &                SS_Type,        &  ! Type of source/sink
     &                GX,             &  ! Generation rate (kg/s)
     &                EX                 ! Corresponding specific enthalpy (J/kg, Injection ONLY !!!)
!
         EL       = ElName5C(1:3)
         F_Format = ElName5C(4:4)
!
         SL       = SS_name(1:3)
         S_Format = SS_name(4:4)
!
         READ(ElName5C(4:5),'(I2)') NE   ! Read NE(target) from ElName5C(4:5)(string); NE is the 2-Digit number component of the element
!
         READ(SS_name(4:5), '(I2)') NS   ! Read NS(target) from SS_name(4:5)(string)
!
! --------
! ... End of the source/sink records - No info on correlation to elements
! --------
!
      IF_El1End: IF( (EL == '   ' .AND. NE == 0) .OR. (EL(1:3) == '<<<') ) THEN
!
                    WRITE(UNIT = GENER_Unit, FMT = 5001)
!
! ......The source/sink data read-in is completed
!
                    RETURN
!
                 END IF IF_El1End
!
! --------
! ... End of the source/sink records - Info on correlation to elements
! --------
!
      IF_El2End: IF(EL == '+++') THEN
!
                    WRITE(UNIT = GENER_Unit, FMT = 5002)
!
! ......Read the source/sink-associated element numbers
!
                    READ(*,5010)                         (SS(n)%ElemNum, n=1,NumSS)
                    WRITE(UNIT = GENER_Unit, FMT = 5010) (SS(n)%ElemNum, n=1,NumSS)
!
! ......The source/sink data read-in is completed
!
                    WRITE(UNIT = GENER_Unit, FMT = 5001)
                    RETURN
!
                 END IF IF_El2End
!
!
!***********************************************************************
!*                                                                     *
!*         WRITING THE SOURCE/SINK DATA INTO THE 'GENER' FILE          *
!*                                                                     *
!***********************************************************************
!
      NumSS = NumSS+1               ! ...... Count sources/sinks
!
! ... Determine the appropriate format for writing the source/sink name in the GENER file
!
      IF(NS < 10 .AND. S_Format == ' ') THEN
         SS_name = SL//' '//N_Character(NS)
      ELSE IF(NS < 10 .AND. S_Format == '0') THEN
         SS_name = SL//'0'//N_Character(NS)
      ELSE
         N_1     = NS/10
         N_2     = MOD(NS,10)
         SS_name = SL//N_Character(N_1)//N_Character(N_2)
      END IF
!
! ... Writing the element info into file GENER
!
      IF(NE < 10 .AND. F_Format == ' ') THEN
         Elem_W5 = EL//' '//N_Character(NE)
      ELSE IF(NE < 10 .AND. F_Format == '0') THEN
         Elem_W5 = EL//'0'//N_Character(NE)
      ELSE
         N_1     = NE/10
         N_2     = MOD(NE,10)
         Elem_W5 = EL//N_Character(N_1)//N_Character(N_2)
      END IF
!
      WRITE(GENER_Unit,5004) Elem_W5,SS_name,SS_Type,GX,EX
!
! -----------
! ...... Continue reading source/sink data
! -----------
!
      GO TO 1000
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 5001 FORMAT('     ')
 5002 FORMAT('+++  ')
 5003 FORMAT('<<<  ')
!
 5004 FORMAT(A5,A5,25X,A4,1x,3E10.4)
!
 5010 FORMAT(16I5)
!
 5020 FORMAT(4E14.7)
!
 6000 FORMAT(/,'Read_Generation_Data',T50,'[v1.0,  11 September 2007]',/,   &
     &         ':::::   Read the source/sink data from the <GENER> block of the FTSim input data file')
!
 6002 FORMAT(' Write file <GENER> from input data')
 6004 FORMAT('GENER')
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of Read_Generation_Data
!
!
      RETURN
!
      END SUBROUTINE Read_Generation_Data
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE Read_Initial_Conditions
!
! ...... Modules to be used
!
         USE EOS_Default_Parameters
!
         USE Basic_Parameters
         USE General_External_File_Units
         USE General_Control_Parameters
         USE Grid_Geometry
         USE Element_Attributes
!
         USE Utility_Functions, ONLY: N_Character
         USE Sources_and_Sinks
         USE Solution_Matrix_Arrays, ONLY: X
         USE Solver_Parameters, ONLY: MatrixSolver
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*          ROUTINE FOR READING THE INITIAL CONDITION DATA             *
!*             DIRECTLY FROM THE FTSim INPUT DATA BLOCK                *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
      IMPLICIT NONE
!
      REAL(KIND = 8) :: phi
      REAL(KIND = 8) :: origin_of_time, time_to_this_point
!
      REAL(KIND = 8), DIMENSION(3)  :: perm
      REAL(KIND = 8), DIMENSION(10) :: PV, primary_variable_values
!
! -------
! ... Integer variables
! -------
!
      INTEGER :: ier = 0, i, num
      INTEGER :: timesteps_to_this_point, NR_iterations_to_this_point
!
      INTEGER, DIMENSION(0:30) :: ierror = 99999
!
! -------
! ... CHARACTER variables
! -------
!
      CHARACTER(LEN = 3) :: StateIndx, state
      CHARACTER(LEN = 5) :: name
!
      CHARACTER(LEN = 35) :: PVList, primary_variable_names
!
! -------
! ... Namelists
! -------
!
      NAMELIST/Uniform_Initial_Conditions/  state,                  &
     &                                      primary_variable_names, &
     &                                      primary_variable_values
!
      NAMELIST/InCond/  name,    &
     &                  state,   &
     &                  PVList,  &
     &                  PV,      &
     &                  phi,     &
     &                  perm
!
      NAMELIST/ Data_For_Continuation_Run / timesteps_to_this_point,      &    ! Number of total timesteps up to the beginning of the continuation run
     &                                      NR_iterations_to_this_point,  &    ! Number of total NR iterations up to the beginning of the continuation run
     &                                      origin_of_time,               &    ! Time at beginning of simulation
     &                                      time_to_this_point                 ! Time at beginning of continuation run
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of Read_Initial_Conditions
!
      WRITE(VERS_Unit, FMT = 6000)
!
! ... Print headings and write to the INCON file
!
      WRITE(*,5001)
      WRITE(*,6002)
!
      REWIND(UNIT = INCON_Unit)
      WRITE(INCON_Unit, 6004)
!
!***********************************************************************
!*                                                                     *
!*                       Initializations Values                        *
!*                                                                     *
!***********************************************************************
!
      state                          = '   '   ! Uniform initial thermodynamic state
      primary_variable_names         = '   '   ! Definition of properties/conditions used as primary variables
      primary_variable_values        = 0.0d0   ! Values of the initial uniform primary variable values
!
!--------
! ... Read namelist uniform initial conditions
!--------
!
      READ(*, NML = Uniform_Initial_Conditions, IOSTAT = ier)
!
!--------
! ... Stop if there is a problem reading the namelist
!--------
!
      IF (ier /= 0) THEN
         WRITE (*, FMT = 6100)
         STOP
      END IF
!
!--------
! ... Checking validity of state
!--------
!
      IF( ALL(state_name(1:Number_of_states) /= state) ) THEN
         WRITE (*, FMT = 6105) state
         STOP
      END IF
!
      DO i=1, Number_of_states
         IF(state == state_name(i)) THEN
            DefaultStateIndex = i
            EXIT
         END IF
      END DO
!
!--------
! ... Write to file <INCON>
!--------
!
      REWIND (UNIT = INCON_Unit)
!
      WRITE (UNIT = INCON_Unit, FMT = Format_Generic_InitialCond(DefaultStateIndex)), char(34), state, char(34),         &
      &                      char(34), primary_variable_names, char(34), (primary_variable_values(i), i=1,NumComPlus1)
!
!--------
! ... Check for block-end identifier or additional initial conditions
!--------
!
 1000 READ(*,*) name
!
      IF(name(1:3) == '<<<') THEN
!
        WRITE (UNIT = INCON_Unit, FMT = 6106)
!
        RETURN
!
!--------
! ... Check for sub-block indicating element-by-element initialization
!--------
!
      ELSE IF(name == ':::>>') THEN
!
         WRITE(UNIT = INCON_Unit, FMT = 6008)
!
         state  = '   '   ! Element initial thermodynamic state (= phase coexistence)
         PV     = 0.0d0   ! Values of the initial primary variable values in this element
         phi    = 0.0d0   ! Element initial uniform porosity
         perm   = 0.0d0   ! Element initial uniform permeabilities
!
!--------
! ...... Read element name and corresponding namelist
!--------
!
 2000    READ(*,*) name
!
         IF(name == ':::<<') THEN
            WRITE (UNIT = INCON_Unit, FMT = 6108 )
            GO TO 1000
         END IF
!
! ...... Blank lines are not tolerated
!
         IF(name == '     ') THEN
            WRITE(*, FMT = 6109)
            STOP
         END IF
!
!-----------
! ...... Read element specific namelist
!-----------
!
         READ(*, NML = InCond, IOSTAT = ier)
!
!-----------
! ...... Stop if there is a problem reading element specific namelist
!-----------
!
         IF (ier /= 0) THEN
            WRITE (*, FMT = 6110)
            STOP
         END IF
!
!-----------
! ...... Checking validity of state
!-----------
!
         IF( ALL(state_name(1:Number_of_states) /= state) ) THEN
            WRITE (*, FMT = 6115) name, state
            STOP
         END IF
!
!-----------
! ...... Writing into file <INCON>
!-----------
!
         WRITE (UNIT = INCON_Unit, FMT = 6010) name, state, phi, perm, PVList, (PV(i), i=1,NumComPlus1)
!
! ...... Read more data
!
         GO TO 2000
!
      ELSE
!
         WRITE (*, FMT = 6300)
         STOP
!
      END IF
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
 5991 FORMAT(/,'Read_Initial_Conditions',T50,'[v1.0,  26 September 2007]',/,   &
     &         ':::::   Read the initial condition data from the <INITIAL_CONDITIONS> block of the FTSim input data file')
 5001 FORMAT('     ')
!
 5995 FORMAT('<<<  ')
!
 6000 FORMAT(/,'Read_Initial_Conditions',T50,'[v1.0,   04 June      2008]',/,  &
     &         ':::::   Read data specific to the <INITIAL CONDITIONS> data block of the input data file')
!
 6002 FORMAT(' Write file <INCON> from input data')
!
 6004 FORMAT('INCON')
!
 6005 FORMAT('>>>INITIAL_CONDITIONS',/, &
     &       '&Uniform_Initial_Conditions   state                   = "',a3,'",',/,           &
     &       '                              primary_variable_names  = "',a5,'",',/,           &
     &       '                              primary_variable_values = ',2(1pe20.13,1x),/,     &
     &       '                              /')
!
 6008 FORMAT(':::>>>Element-by-element initial conditions')
!
 6010 FORMAT('"',a5,'"',/, &
     &       '&InCond  state   = "',a3,'",',/,              &
     &       '         phi     = ',1pe15.9,',',/            &
     &       '         perm    = ',3(1pe15.9,',',1x),/,     &
     &       '         PVList  = "',a5,'",',/,              &
     &       '         PV      = ',2(1pe20.13,1x),/,        &
     &       '         /')
!
 6100 FORMAT(//,22('ERROR-'),//,T33,  &
     &             '       S I M U L A T I O N   A B O R T E D',/,                                                                 &
     &         T22,'There is a problem reading the namelist <Uniform_Initial_Conditions> in subroutine <Read_Initial_Conditions> ',&
     &       /,T32,'               CORRECT AND TRY AGAIN',                                                                         &
     &       //,22('ERROR-'))
!
 6105 FORMAT(//,22('ERROR-'),//,T43,   &
     &             'S I M U L A T I O N   A B O R T E D',/,   &
     &         T25,'The state identifier of uniform initial conditions <state> = "',a3,'"',/,   &
     &         T31,'This is not an available option for this Equation of State',   &
     &       /,T50,'CORRECT AND TRY AGAIN',   &
     &       //,22('ERROR-'))
!
 6106 FORMAT('<<<END of INITIAL_CONDITIONS')
 6108 FORMAT(':::<<<END of Element-by-element initial conditions')
 !
 6109 FORMAT(//,22('ERROR-'),//,T33,  &
     &             '       S I M U L A T I O N   A B O R T E D',/,                                                                 &
     &         T22,'There is a blank line within the element-specific namelist <InCond> in subroutine <Read_Initial_Conditions> ', &
     &       /,T32,'               CORRECT AND TRY AGAIN',                                                                         &
     &       //,22('ERROR-'))
!
 6110 FORMAT(//,22('ERROR-'),//,T33,  &
     &             '       S I M U L A T I O N   A B O R T E D',/,                                                                 &
     &         T22,'There is a problem reading the element-specific namelist <InCond> in subroutine <Read_Initial_Conditions> ',&
     &       /,T32,'               CORRECT AND TRY AGAIN',                                                                         &
     &       //,22('ERROR-'))
!
 6115 FORMAT(//,22('ERROR-'),//,T43,   &
     &             'S I M U L A T I O N   A B O R T E D',/,   &
     &         T25,'In element "',a5,'", the state identifier of initial conditions <state> = "',a3,'"',/,   &
     &         T31,'This is not an available option for this Equation of State',   &
     &       /,T50,'CORRECT AND TRY AGAIN',   &
     &       //,22('ERROR-'))
!
 6300 FORMAT(//,22('ERROR-'),//,T33,  &
     &             '       S I M U L A T I O N   A B O R T E D',/,                                                         &
     &         T22,'There is a problem exiting the <INITIAL CONDITIONS> input block (check for ":::" and "<<<" formats) ', &
     &       /,T32,'               CORRECT AND TRY AGAIN',                                                                 &
     &       //,22('ERROR-'))
!
 6400 FORMAT(F10.2,6(5X,ES10.2E2))
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of Read_Initial_Conditions
!
      RETURN
!
      END SUBROUTINE Read_Initial_Conditions
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE Read_Printout_Time_Data
!
! ...... Modules to be used
!
         USE General_Control_Parameters, ONLY: NumPrintTimes, TimeStepMax_PostPrint, PrintTime
!
         USE General_External_File_Units
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*             ROUTINE FOR READING THE PRINT-OUT TIME DATA             *
!*              DIRECTLY FROM THE FTSim INPUT DATA BLOCK               *
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
      REAL(KIND = 8) :: PrintTimeIncrement
!
! -------
! ... Integer variables
! -------
!
      INTEGER :: Max_NumPrintTimes,i
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of Read_Printout_Time_Data
!
!
      WRITE(VERS_Unit,6000)
!
!***********************************************************************
!*                                                                     *
!*      READING THE TIME DATA: 1st record - Basic time parameters      *
!*                                                                     *
!***********************************************************************
!
      READ(*,5001) NumPrintTimes,         &  ! Number of provided print-out times (NumPrintTimes<=100)
     &             Max_NumPrintTimes,     &  ! Number of desired print-out times  (NumPrintTimes<=Max_NumPrintTimes<=100)
     &             TimeStepMax_PostPrint, &  ! Maximum Dt size after any of the proscribed times have been reached (sec)
     &             PrintTimeIncrement        ! Time increment for times NumPrintTimes,NumPrintTimes+1,NumPrintTimes+2,...
!
!***********************************************************************
!*                                                                     *
!*         READING THE TIME DATA: 2nd record - Print-out times         *
!*                                                                     *
!***********************************************************************
!
      READ(*,5002) (PrintTime(i), i=1,NumPrintTimes)
!
! ... Reading of print-out times is completed
!
      IF(Max_NumPrintTimes <= NumPrintTimes .OR. PrintTimeIncrement <= 0.0d0) RETURN
!
!***********************************************************************
!*                                                                     *
!*    ASSIGNING INCREMENTED PRINT-OUT TIMES FOR ITE > NumPrintTimes    *
!*                                                                     *
!***********************************************************************
!
      DO i = NumPrintTimes+1, Max_NumPrintTimes
         PrintTime(i) = PrintTime(i-1) + PrintTimeIncrement
      END DO
!
! ... Reset NumPrintTimes
!
      NumPrintTimes = Max_NumPrintTimes
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 5001 FORMAT(2I5,3E10.4)
 5002 FORMAT(8E10.4)
!
 5007 FORMAT(A6)
!
 6000 FORMAT(/,'Read_Printout_Time_Data',T50,'[v1.0,   7 January   2007]',/,   &
     &         ':::::   Read the print-out times from the <TIMES> block of the FTSim input data file')
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of Read_Printout_Time_Data
!
!
      RETURN
!
      END SUBROUTINE Read_Printout_Time_Data
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE Read_External_Input_Files(Flag_MINC)
!
! ...... Modules to be used
!
         USE EOS_Routine_Selector
!
         USE Basic_Parameters, ONLY: EOS_Name
         USE General_External_File_Units
!
         USE General_Control_Parameters, ONLY: NoFlowSimulation
!
         USE Grid_Geometry, ONLY: elem
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*  ROUTINE FOR READING INPUT DATA FROM DISK FILES, WHICH ARE EITHER   *
!*        PROVIDED BY THE USER, OR ARE INTERNALLY GENERATED IN         *
!*        SUBROUTINE INPUT FROM DATA GIVEN IN THE JOB DECK             *
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
      INTEGER :: fileno,n2x,IM
!
! -------
! ... Character variables
! -------
!
      CHARACTER(LEN = 9)  :: suffix
      CHARACTER(LEN = 14) :: filenm
!
! -------
! ... Logical variables
! -------
!
      LOGICAL :: exists,file_GENER_open
!
      LOGICAL, INTENT(IN) :: Flag_MINC
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of Read_External_Input_Files
!
!
      WRITE(VERS_Unit,6000)
!
! ----------
! ... Opening files for writing outputs for plotting
! ----------
!
      DO_PlotF: DO n2x=1,3
!
         IF(n2x == 1) THEN
            suffix = 'coords_EC'               ! Suffix for the coordinate file
         ELSE IF(n2x == 2) THEN
            suffix = 'Data_Elem'               ! Suffix for the element-related output data file
         ELSE
            suffix = 'Data_Conx'               ! Suffix for the connection-related output data file
         END IF
!
         filenm = 'Plot_'//suffix              ! The new file name
         fileno = Plot_coord_Unit - 1 + n2x    ! The new file number
!
! -------------
! ...... Determine if the files exist
! -------------
!
         INQUIRE(FILE=filenm, EXIST=exists)
!
         IF(exists) THEN                                 ! If the file exists ...
            WRITE(*,6001) filenm                         !    print info
            OPEN(fileno, FILE = filenm, STATUS = 'OLD')  !    open as an old file
            REWIND (UNIT = fileno)                       !    go to the top of the file
         ELSE                                            ! If the file does not exist ...
            WRITE(*,6002) filenm                         !    print info
            OPEN(fileno, FILE = filenm, STATUS='NEW')    !    open as a new file
         END IF
!
      END DO DO_PlotF
!
! ... Printing for outout file formatting purposes (cosmetic)
!
      WRITE(*,6005)
!
!
!***********************************************************************
!*                                                                     *
!*          Read data from the ELEME block of the file MESH            *
!*                                                                     *
!***********************************************************************
!
!
      IF_FileNo: IF(Flag_MINC) THEN
                    IM = MINC_Unit        ! MINC-type mesh file (fractured medium)
                 ELSE
                    IM = MESH_Unit        ! Regular MESH file (non-fractured medium)
                 END IF IF_FileNo
!
!
!
      CALL Read_Elements_From_MESH(IM)
      IF(NoFlowSimulation .EQV. .TRUE.) RETURN
!
!
!***********************************************************************
!*                                                                     *
!*          Read data from the CONNE block of the file MESH            *
!*                                                                     *
!***********************************************************************
!
!
      CALL Read_Connections_From_MESH(IM)
      IF(NoFlowSimulation .EQV. .TRUE.) RETURN
!
!
!***********************************************************************
!*                                                                     *
!*               Read sink/source data from file GENER                 *
!*                                                                     *
!***********************************************************************
!
!
      INQUIRE(FILE = 'GENER', OPENED = file_GENER_open)
!
!     IF(file_GENER_open) THEN
         CALL Read_File_GENER
!     ELSE
!        CALL READ_File_SinkSource(EOS_Specific_Flag)
!     END IF
!
      IF(NoFlowSimulation .EQV. .TRUE.) RETURN
!
!
!***********************************************************************
!*                                                                     *
!*            Read initial condition data from file INCON              *
!*                                                                     *
!***********************************************************************
!
!
      CALL Read_Initial_Conditions_File
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
 5007 FORMAT(A6)
!
 6000 FORMAT(/,'Read_External_Input_Files',T50,'[v1.0,  12 January   2007]',/,   &
     &         ':::::   Initialize data from files <MESH> or <MINC>, <GENER> and <INCON>',/,   &
     &         '        Also initializes permeability modifiers and coordinate arrays',    &
     &         ' and optionally reads tables with flowing wellbore pressures')
!
 6001 FORMAT(' File <',a10,'> does not exist ==> Open as a new file')
 6002 FORMAT(' File <',a10,'> exists ==> Open as an old file')
 6005 FORMAT(' ')
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End of Read_External_Input_Files
!
!
      RETURN
!
      END SUBROUTINE Read_External_Input_Files
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE Read_Elements_From_MESH(IM)
!
! ...... Modules to be used
!
         USE Basic_Parameters
         USE General_External_File_Units
         USE General_Control_Parameters
!
         USE Grid_Geometry
         USE Element_Attributes
!
         USE Geologic_Media_Parameters
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*     ROUTINE FOR READING THE DATA DESCRIBING THE DOMAIN ELEMENTS     *
!*                   DIRECTLY FROM THE FILE 'MESH'                     *
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
      REAL(KIND = 8) :: elem_vol, elem_vol_prev, X_co, Y_co, Z_co
!
! -------
! ... Integer variables
! -------
!
      INTEGER, INTENT(IN) :: IM
!
      INTEGER :: ier = 0
!
      INTEGER :: N,M,matmat,i,i_ASCII
      INTEGER :: elem_MatNum,elem_number
      INTEGER :: num_act_elem = 0, num_inact_elem = 0
!
! -------
! ... Character variables
! -------
!
      CHARACTER(LEN = 1) :: elem_activity
      CHARACTER(LEN = 5) :: MA12,DENT
      CHARACTER(LEN = 8) :: elem_name
!
! -------
! ... Logical variables
! -------
!
      LOGICAL :: Medium_ByName
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of Read_Elements_From_MESH
!
!
      WRITE(VERS_Unit,6000)
!
! ----------
! ... Reading the headings of the file ÒMESHÓ
! ----------
!
      REWIND (UNIT = IM)
      READ(IM,5001) DENT
!
! ----------
! ... Initialization
! ----------
!
      NumElem = 0             ! The number of active gridblocks
!
!
!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>                                                                     >
!>                       MAIN ELEMENT LOOPS                            >
!>                                                                     >
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
      DO_NumEle1: DO n=1,10000000
!
! -------------
! ...... Check if "Max_NumElem" is exceeded
! -------------
!
         IF(n > Max_NumElem) THEN
            WRITE(*,6001) Max_NumElem
            NoFlowSimulation = .TRUE.        ! Bypass simulation
            RETURN
         END IF
!
! ...... Initialize the element name
!
         elem_name =    '        '
!
!
! >>>>>>>>>>>>>>>>>>
! ......
! ...... Read the element data from the ELEME block of MESH
! ......
! >>>>>>>>>>>>>>>>>>
!
         READ(IM,5001) elem_name(1:5), &  ! 5-Character element name
     &                 MA12,           &  ! 5-Character name of corresponding rock type
     &                 elem_vol,       &  ! Element volume (m^3)
     &                 X_co,           &  ! X-coordinate element center (m)
     &                 Y_co,           &  ! Y-coordinate of element center (m)
     &                 Z_co,           &  ! Z-coordinate of element center (m)
     &                 elem_activity      ! Activity flag of the element
!
!
!***********************************************************************
!*                                                                     *
!*                Processing the element-specific data                 *
!*                                                                     *
!***********************************************************************
!
! -------------
! ...... If the element name is blank, the element list is exhausted
! -------------
!
         IF_ElEnd: IF(elem_name(1:5) == '     ') THEN
!
            NumElemTot = n-1
            NumElem    = num_act_elem
!
            NumUnknowns = NumElem*NumEqu
!
            EXIT DO_NumEle1
!
         END IF IF_ElEnd
!
! -------------
! ...... Determine the number of active elements
! -------------
!
         IF_Inactive: IF( elem_activity == 'I' ) THEN
!
            num_inact_elem    = num_inact_elem + 1
            elem_number       = Max_NumElem - num_inact_elem + 1
!
         ELSE
!
            IF( (elem_vol <= 1.0d-10) .OR. (elem_vol >= 1.0d+10) ) THEN
               elem_activity  = 'I'
               num_inact_elem = num_inact_elem + 1
               elem_number    = Max_NumElem - num_inact_elem + 1
            ELSE
               elem_activity = ' '
               num_act_elem  = num_act_elem + 1
               elem_number   = num_act_elem
            END IF
!
         END IF IF_Inactive
!
! -------------
! ...... Storing the volume of current element
! -------------
!
         elem_vol_prev = elem_vol
!
! -------------
! ...... Determine if the medium type is by name or by number
! -------------
!
         Medium_ByName = .FALSE.                     ! Default: Medium by number
!
         DO_MedName: DO i=1,5
!
            i_ASCII = ICHAR(MA12(i:1))               ! Determine the ASCII value of each character of the medium name
!
            IF(i_ASCII == 32) CYCLE DO_MedName       ! Continue iteration for blanks
!
            IF(i_ASCII > 57 .OR. i_ASCII <48) THEN   ! If any character is a letter, ...
               Medium_ByName = .TRUE.                ! ... the flag is reset
               EXIT DO_MedName
            END IF
!
         END DO DO_MedName
!
! -------------
! ...... Assign rock type to elements
! -------------
!
         IF_Mat1: IF(Medium_ByName .EQV. .FALSE.) THEN
!
            READ(MA12,'(I5)') elem_MatNum     ! Elements with a domain number
                                              ! Read MATX(n)(target) from MA12(string)
!
         ELSE
!
! ......... Elements with a domain name
!
            matmat = 0
            DO_NumRok1: DO m=1,NumMedia
!
               IF(MA12 == MediumName(m)) THEN
                  elem_MatNum = m             ! Find material index
                  matmat      = 1
                  EXIT DO_NumRok1
               END IF
!
            END DO DO_NumRok1
!
! ......... On failure to match, print warning
!
            IF_NoMatch1: IF(matmat == 0) THEN
               WRITE(*,6004) MA12,elem_name
               elem_MatNum = 0
            END IF IF_NoMatch1
!
         END IF IF_Mat1
!
! -------------
! ...... Elements without domain assignment are assigned to Domain #1 by default
! -------------
!
         IF_DomDef1: IF(elem_MatNum == 0) THEN
            elem_MatNum = 1
            WRITE(*,6005) elem_name, MediumName(1)
         END IF IF_DomDef1
!
!***********************************************************************
!*                                                                     *
!*   Renumbering based on activity status & assignment of properties   *
!*        NOTE: Inactive element numbering sequence is reverse         *
!*                       to order of appearance                        *
!*                                                                     *
!***********************************************************************
!
         elem(elem_number)%name     = elem_name
         elem(elem_number)%MatNum   = elem_MatNum
         elem(elem_number)%vol      = elem_vol
         elem(elem_number)%activity = elem_activity
!
         elem(elem_number)%coord(1) = X_co
         elem(elem_number)%coord(2) = Y_co
         elem(elem_number)%coord(3) = Z_co
!
! <<<
! <<< End of the ELEMENT LOOP
! <<<
!
      END DO DO_NumEle1
!
!***********************************************************************
!*                                                                     *
!*                    RESET THE NUMEBERING SEQUENCE:                   *
!* First the active, then the inactive elements in order of appearance *
!*                                                                     *
!***********************************************************************
!
      CALL Renumber_Elements( num_inact_elem, NumElem, NumElemTot, Max_NumElem, MemoryAlloc_Unit )
!
! ----------
! ... Print the coordinates of active elements in file "Plot_Coord"
! ----------
!
      WRITE(Plot_coord_Unit,6050)
      DO n = 1,NumElemTot
         WRITE(Plot_coord_Unit,6002) n,elem(n)%coord(1),elem(n)%coord(2),elem(n)%coord(3)
      END DO
!
!***********************************************************************
!*                                                                     *
!*           DETERMINE NUMBER OF FLAG ELEMENT (FOR PRINT-OUT)          *
!*                                                                     *
!***********************************************************************
!
      TrackElemNum = 0
!
      IF_ActEl2: IF(TrackElemName(1:5) /= '     ') THEN
!
         DO_NumEleB: DO n = 1,NumElemTot
!
            IF(elem(n)%name == TrackElemName) THEN
               TrackElemNum = n
               EXIT
            END IF
!
         END DO DO_NumEleB
!
      END IF IF_ActEl2
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 5001 FORMAT(A5,10X,A5,E10.4,20x,3E10.4,1x,a1)
!
 5007 FORMAT(A6)
!
 5101 FORMAT(A5,10X,i5,E10.4,20x,3E10.4,1x,a1)
!
!
 6000 FORMAT(/,'Read_Elements_From_MESH',T50,'[v1.0,  19 August    2007]',/,   &
     &         ':::::   Read the element data from the <ELEME> block of the file <MESH>')
!
 6001 FORMAT(' NUMBER OF ELEMENTS SPECIFIED IN DATA BLOCK "ELEME" EXCEEDS ALLOWABLE MAXIMUM OF ',I7,/, &
     &       ' INCREASE PARAMETER <Max_NumElem> IN MAIN PROGRAM,',//,   &
     &       ' >>>>>>>>>>   F L O W   S I M U L A T I O N   A B O R T E D   <<<<<<<<<<')
 6002 FORMAT(i6,2x,3(1pE15.8,1x))
!
 6004 FORMAT(' REFERENCE TO UNKNOWN MATERIAL "',A5,'" AT ELEMENT "',A8,'" ==> WILL IGNORE ELEMENT')
 6005 FORMAT(' !!!!! ==> WARNING <== !!!!! :',   &
     &       ' ELEMENT "',A8,'" HAS NO DOMAIN ASSIGNMENT; ASSIGN TO DOMAIN # 1, "',A5,'"')
!
 6050 FORMAT('Gridblock center coordinates: Elem # - X[m] - Y[m] - Z [m] at cell center')
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of Read_Elements_From_MESH
!
!
      RETURN
!
      END SUBROUTINE Read_Elements_From_MESH
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE Read_Connections_From_MESH(IM)
!
! ...... Modules to be used
!
         USE Basic_Parameters
         USE General_External_File_Units
         USE General_Control_Parameters
!
         USE Grid_Geometry
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*   ROUTINE FOR READING THE DATA DESCRIBING THE DOMAIN CONNECTIONS    *
!*                   DIRECTLY FROM THE FILE ÒMESHÓ                     *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
      IMPLICIT NONE
!
! -------
! ... Integer arrays
! -------
!
      INTEGER, DIMENSION(3) :: isumd = 0
!
! -------
! ... Integer variables
! -------
!
      INTEGER, INTENT(IN) :: IM
!
      INTEGER :: i, N, ipluss
!
! -------
! ... Character variables
! -------
!
      CHARACTER(LEN = 5) :: EL5_CR,DENT
!
! -------
! ... Logical variables
! -------
!
      LOGICAL :: known_ElemNum_in_connections = .FALSE.
!
#ifdef USE_TIMER
         real :: start, finish
#endif
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of Read_Connections_From_MESH
!
!
      WRITE(UNIT = VERS_Unit, FMT = 6000)
!
! ----------
! ... Initializations
! ----------
!
#ifdef USE_TIMER
         call cpu_time(start)
#endif

#ifdef USE_OMP
!$OMP WORKSHARE
#endif
      NumConx = 0
      ipluss  = 0
!
      isumd  = 0  ! Whole array operation
#ifdef USE_OMP
!$OMP END WORKSHARE
#endif

#ifdef USE_TIMER
         call cpu_time(finish)

         write (*,*) __FILE__, ":", __LINE__, " time: ", finish-start
#endif
!
! ----------
! ... Reading the headings of the block ÒCONNEÓ
! ----------
!
      i = 0
 1000 READ(UNIT = IM, FMT = 5001) DENT
      IF(DENT /= 'CONNE') THEN
         i = i + 1
         IF( i > 5 ) THEN
            WRITE (*, FMT = 6001)
            STOP
         END IF
         GO TO 1000
      END IF
!
! ... Initialization - CAREFUL! Whole array operations
!
#ifdef USE_TIMER
         call cpu_time(start)
#endif

#ifdef USE_OMP
!$OMP WORKSHARE
#endif
         conx%n1 = 0
         conx%n2 = 0
#ifdef USE_OMP
!$OMP END WORKSHARE
#endif

#ifdef USE_TIMER
         call cpu_time(finish)

         write (*,*) __FILE__, ":", __LINE__, " time: ", finish-start
#endif
!
         conx%name1 = '        '
         conx%name2 = '        '
!
!
!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>                                                                     >
!>                       MAIN CONNECTION LOOP                          >
!>                                                                     >
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
      DO_NumCon: DO n=1,10000000
!
! -------------
! ...... Check if "Max_NumConx" is exceeded
! -------------
!
         IF(n > Max_NumConx) THEN
            WRITE(*,6002) Max_NumConx
            NoFlowSimulation = .TRUE.        ! Bypass simulation
            RETURN
         END IF
!
! ----------------
! ......... Read the element data from the CONNE block of MESH - 8-character elements
! ----------------
!
         READ(IM,5001) conx(n)%name1,  &  ! 5-Character name of connection element #1
     &                 conx(n)%name2,  &  ! 5-Character name of connection element #2
     &                 conx(n)%ki,     &  ! Directional indicator of the permeability tensor
     &                 conx(n)%d1,     &  ! Distance of the center of element #1 from interface (m)
     &                 conx(n)%d2,     &  ! Distance of the center of element #2 from interface (m)
     &                 conx(n)%area,   &  ! Interface area (m^2)
     &                 conx(n)%beta       ! = COS(b), b = angle between the vertical and the connection line
!
         IF(conx(n)%name1 == conx(n)%name2 .AND. conx(n)%name1(1:5) /= '     ') THEN
            WRITE(*,6075) n
            STOP
         END IF
!
         EL5_CR = conx(n)%name1(1:5)
!
! -------------
! ...... If EL5_CR = Ô+++  Ô, read element numbers in the connection list
! -------------
!
         IF_ConEnd1: IF(EL5_CR == '+++  ') THEN
!
            NumConx = n-1                                          ! Define actual NumConx
!
            READ(IM,5003) (conx(i)%n1,conx(i)%n2, i=1,NumConx)     ! Read connection-element # data, 5-ch. names
!
            known_ElemNum_in_connections = .TRUE.
!
            EXIT                                                   ! Done !
!
         END IF IF_ConEnd1
!
! -------------
! ...... If EL5_CR is blank, the element list is exhausted
! -------------
!
         IF_ConEnd2: IF(EL5_CR == '     ') THEN
!
            NumConx = n-1                             ! Define actual NumConx
            BACKSPACE IM                              ! Move back one line in MESH
            EXIT                                      ! Done !
!
         END IF IF_ConEnd2
!
! -------------
! ...... Determine the number of connections along the principal directions
! -------------
!
         IF(conx(n)%ki == 1) isumd(1) = isumd(1)+1
         IF(conx(n)%ki == 2) isumd(2) = isumd(2)+1
         IF(conx(n)%ki == 3) isumd(3) = isumd(3)+1
!
         IF_3to2D: IF(conx(n)%ki < 0) THEN
                      IF(conx(n)%ki == -2) THEN
                         isumd(3) = isumd(3)+1
                      ELSE
                         isumd(1) = isumd(1)+1
                      END IF
                   END IF IF_3to2D
!
! <<<
! <<< End of the CONNECTION LOOP
! <<<
!
      END DO DO_NumCon
!
!***********************************************************************
!*                                                                     *
!*                     DETERMINE ACTIVE DIMENSIONS                     *
!*                                                                     *
!***********************************************************************
!
      NumActiveDimensions = COUNT(isumd > 0)         ! Number of active dimensions, array operation
!
      IF(NumActiveDimensions /= 1 .AND. NumActiveDimensions /= 2 .AND. NumActiveDimensions /= 3) THEN
         WRITE(*,6078) NumActiveDimensions
         STOP
      END IF
!
! -------------
! ...... Determine the 2nd active direction in a 2D system (X,R is always the 1st)
! -------------
!
      IF_2D: IF(NumActiveDimensions == 2) THEN
!
         IF(isumd(2) == 0) THEN
            Active2ndDimension = 3                   ! The 3rd direction is active
         ELSE IF(isumd(3) == 0) THEN
            Active2ndDimension = 2                   ! The 2nd direction is active
         ELSE
            WRITE(*,6079) Active2ndDimension
         END IF
!
      END IF IF_2D
!
      WRITE(*,6080) NumActiveDimensions
      WRITE(*,6081) isumd(1),isumd(2),isumd(3)
!
      longest(1:1)   = MAXLOC(isumd)
!
!***********************************************************************
!*                                                                     *
!*    DETERMINE IF THE ELEMENT NUMBERS OF THE CONNECTIONS ARE KNOWN    *
!*                                                                     *
!***********************************************************************
!
      IF(known_ElemNum_in_connections) RETURN
!
!***********************************************************************
!*                                                                     *
!*          DETERMINE THE ELEMENT NUMBERS OF THE CONNECTIONS           *
!*                                                                     *
!***********************************************************************
!
      CALL Determine_Connection_Elements( NumElemTot, NumConx )
!
!***********************************************************************
!*                                                                     *
!*           WRITE THE ELEMENT NUMBERS OF THE CONNECTIONS              *
!*                  AT THE END OF THE <MESH> FILE                      *
!*                                                                     *
!***********************************************************************
!
      WRITE(IM,6003)      ! Write '+++  ' at the end of MESH
!
      WRITE(IM,5003) (conx(i)%n1,conx(i)%n2, i=1,NumConx)
!
      ENDFILE IM
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 5001 FORMAT(2A5,15X,I5,5(E10.4))
!
 5003 FORMAT(16I6)
!
 5007 FORMAT(A6)
!
 6000 FORMAT(/,'Read_Connections_From_MESH',T50,'[v1.0,  20 August    2007]',/,   &
     &         ':::::   Read the connection data from the <CONNE> block of the file <MESH>')
!
 6001 FORMAT(//, ' !!!!! ==> WARNING <== !!!!! :',  &
     &        /, ' There are least 5 lines between the bottom of the ELEME block',   &
     &           ' and the beginning of the CONNE block in the MESH file ',/,        &
     &           ' >>> The MESH file may be corrupted - The simulation is ABORTED ', &
     &       //, ' CHECK, CORRECT AND TRY AGAIN !!! ')
!
 6002 FORMAT('NUMBER OF CONNECTIONS SPECIFIED IN DATA BLOCK <CONNE> EXCEEDS ALLOWABLE MAXIMUM OF ',I8,/, &
     &       ' INCREASE PARAMETER <Max_NumConx> IN MAIN PROGRAM',//,   &
     &       ' !!!!! ==> ==> ==> FLOW SIMULATION IS ABORTED')
 6003 FORMAT('+++  ')
!
 6004 FORMAT(' REFERENCE TO UNKNOWN MATERIAL ',A5,' AT ELEMENT "',A8,'" ==> WILL IGNORE ELEMENT')
 6005 FORMAT(' !!!!! ==> WARNING <== !!!!! :',   &
     &       ' ELEMENT "',A8,'" HAS NO DOMAIN ASSIGNMENT; ASSIGN TO DOMAIN # 1, "',A5,'"')
 6008 FORMAT(' REFERENCE TO UNKNOWN ELEMENT ',A8,' AT CONNECTION ',I6,' ==> WILL IGNORE CONNECTION')
!
 6075 FORMAT(//,20('ERROR-'),//,T33,   &
     &             '       S I M U L A T I O N   A B O R T E D',/,   &
     &         T31,'At connection number ',i7,', the two connection elements are the same',   &
     &       /,T33,'               CORRECT AND TRY AGAIN',   &
     &       //,20('ERROR-'))
 6078 FORMAT(//,20('ERROR-'),//,T33,   &
     &             '       S I M U L A T I O N   A B O R T E D',/,   &
     &         T5,'The number of active dimensions <NumActiveDimensions> = ',i1,   &
     &            ' is not an acceptable number (must be 1, 2 or 3)',   &
     &       /,T33,'               CORRECT AND TRY AGAIN',   &
     &       //,20('ERROR-'))
 6079 FORMAT(//,20('ERROR-'),//,T33,   &
     &             '       S I M U L A T I O N   A B O R T E D',/,   &
     &         T3, 'The index of the second active dimension <Active2ndDimension> = ',i2,   &
     &             ' is not an acceptable number (must be 2 or 3)',   &
     &       /,T33,'               CORRECT AND TRY AGAIN',   &
     &       //,20('ERROR-'))
 6080 FORMAT(//,130('*'),/,'*',T130,'*',/,'*',   &
     &       T10,'The number of active dimensions <NumActiveDimensions> =',i2,   &
     &       T130,'*',/,'*',T130,'*')
 6081 FORMAT('*',T10,'The numbers of connections in the X-, Y- and Z-directions are ',i6,', ',i6,' and ',i6,' respectively', &
     &           T130,'*',/,'*',T130,'*',/,130('*'),//)
!
 6090 FORMAT(//,'Connection coordinates: Conx # - elem_1 # - elem_2 # - X[m] - Y[m] - Z [m] at interface')
 6092 FORMAT(i7,2x,i6,2x,i6,2x,3(1pE15.8,1x))
!
 6101 FORMAT(T2,' Memory deallocation of array <emissivity> at point ',i3, &
     &          ' in subroutine <Read_Connections_From_MESH> was successful')
 6102 FORMAT(//,20('ERROR-'),//,   &
     &       T2,' Memory deallocation of array <emissivity> at point ',i3, &
     &          ' in subroutine <Read_Connections_From_MESH> was unsuccessful', &
     &       //,T2,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6801 FORMAT(//,20('ERROR-'),//,T33,   &
     &             '       S I M U L A T I O N   A B O R T E D',/,   &
     &         T2,'The number of characters of the element names <ElemNameLength> = ',i3,   &
     &            ' is not an acceptable number (must be either 5 or 8)',   &
     &       /,T33,'               CORRECT AND TRY AGAIN',   &
     &       //,20('ERROR-'))
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of Read_Connections_From_MESH
!
!
      RETURN
!
      END SUBROUTINE Read_Connections_From_MESH
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE Read_File_GENER
!
! ...... Modules to be used
!
         USE Basic_Parameters
         USE General_External_File_Units
         USE General_Control_Parameters
!
         USE Grid_Geometry, ONLY: elem
!
         USE Sources_and_Sinks
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*        ROUTINE FOR READING THE SOURCE/SINK GENERATION DATA          *
!*                   DIRECTLY FROM THE FILE ÒGENERÓ                    *
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
      INTEGER :: NGL  = 0
!
      INTEGER :: i, j, n, kf, non_common_SS
!
! -------
! ... Character variables
! -------
!
      CHARACTER(LEN = 4) :: SS_Typ
      CHARACTER(LEN = 5) :: EL5_CR, DENT
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of Read_File_GENER
!
!
      WRITE(VERS_Unit,6000)
!
! ----------
! ... Initializations
! ----------
!
      NumSS = 0
      kf    = 0
!
! ... CAREFUL! Whole array operations
!
      SS%ElemName = '     '
      SS%ElemNum  = 0
!
! ----------
! ... Reading the headings of the file ÒGENERÓ
! ----------
!
      REWIND (UNIT = GENER_Unit)
      READ(GENER_Unit,5001) DENT
!
!
!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>                                                                     >
!>                         SINK/SOURCE LOOP                            >
!>                                                                     >
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
      DO_NumSS: DO n=1,Max_NumSS+1
!
! -------------
! ...... Check if Max_NumSS is exceeded
! -------------
!
         IF_Limit: IF(n == Max_NumSS+1) THEN
!
            READ(UNIT = GENER_Unit, FMT = 5001) EL5_CR
!
            IF(EL5_CR /= '+++  ' .AND. EL5_CR /= '     ') THEN
               WRITE(*,6001) Max_NumSS
               NoFlowSimulation = .TRUE.
               RETURN
            ELSE IF (EL5_CR == '+++  ') THEN
               GO TO 1000
            ELSE IF (EL5_CR == '     ') THEN
               GO TO 2000
            END IF
!
         END IF IF_Limit
!
! -------------
! ...... Reading Source/Sink data
! -------------
!
         READ(GENER_Unit,5001) SS(n)%ElemName,   &  ! 5-Character name of cell containing
     &                         SS(n)%name,       &  ! 5-Character name of source/sink (well name)
     &                         SS(n)%type,       &  ! Type of source/sink
     &                         SS(n)%rate,       &  ! Injection/production rate (kg/s) or PI for SS(n)%type = 'DELV'
     &                         SS(n)%enth           ! Specific enthalpy of injected fluid (J/kg) or bottomhole pressure for SS(n)%type = ÔDELVÕ
!
         EL5_CR = SS(n)%ElemName(1:5)
!
!
!
         SS_Typ = SS(n)%type(1:4)
!
!
!***********************************************************************
!*                                                                     *
!*  If EL5_CR = Ô+++  Ô, read element numbers in the source/sink list  *
!*                                                                     *
!***********************************************************************
!
 1000    IF_SSEnd1: IF(EL5_CR == '+++  ') THEN
!
            NumSS = n-1                                       ! Define actual NumSS
!
            READ(GENER_Unit,5003) (SS(i)%ElemNum, i=1,NumSS)  ! Read connection-element # data, 5-ch. names
!
            RETURN                                            ! End of the GENER file reached. Done !
!
         END IF IF_SSEnd1
!
!***********************************************************************
!*                                                                     *
!*         If EL5_CR is blank, the element list is exhausted           *
!*                                                                     *
!***********************************************************************
!
 2000    IF_SSEnd2: IF( (EL5_CR == '     ') .OR. (EL5_CR(1:3) == '<<<') ) THEN
!
            NumSS = n-1                           ! Define actual NumSS
!
            BACKSPACE (UNIT = GENER_Unit)         ! Move back one line in GENER
            WRITE(GENER_Unit,6002)                ! Write Ô+++  Ô
!
            IF(NumSS <= 0) RETURN
!
! ----------------
! ......... Determine the element number in the cell that contains the source/sink
! ----------------
!
            DO j=1,NumElemTot
               FORALL (i=1:NumSS, elem(j)%name == SS(i)%ElemName) SS(i)%ElemNum = j
            END DO
!
! ----------------
! ......... Determining the sources not corresponding to actual elements
! ----------------
!
            non_common_SS = COUNT(SS(1:NumSS)%ElemNum == 0)   ! CAREFUL! Whole array operations
!
! ----------------
! ......... Printing a warning, and continuing
! ----------------
!
            IF(non_common_SS > 0) THEN
               DO_ObsSS: DO i=1,NumSS
                  IF(SS(i)%ElemNum == 0) THEN
                     WRITE(*,6004) ADJUSTL(SS(i)%ElemName),i
                  END IF
               END DO DO_ObsSS
            END IF
!
! ----------------
! ......... Print the source/sink element numbers at the end of the GENER file
! ----------------
!
            WRITE(UNIT = GENER_Unit, FMT = 5003) (SS(i)%ElemNum, i=1,NumSS)
!
            ENDFILE(UNIT = GENER_Unit)             ! Close the GENER file
!
            RETURN                                 ! The processing of data in GENER is complete
!
         END IF IF_SSEnd2
!
!***********************************************************************
!*                                                                     *
!*       READING AND PROCESSING SINK/SOURCE DATA IS IN PROGRESS        *
!*                                                                     *
!***********************************************************************
!

!
! -------------
! ...... Specifying source/sink option
! -------------
!
         SELECT_SS: SELECT CASE(SS_Typ)
!
! -------------------------
! ................. Source/sink fluid component is MASS
! -------------------------
!
                    CASE('MASS')
                       SS(n)%TypeIndx = 1
                       SS(n)%RateType = 'MASS'
!
! -------------------------
! ................. Source/sink fluid component is COM1 = H2O
! -------------------------
!
                    CASE('COM1')
                       SS(n)%TypeIndx = 1
                       SS(n)%RateType = 'COM1'
!
! -------------------------
! ................. Source/sink fluid component is COM2 = Air
! -------------------------
!
                    CASE('COM2')
                       SS(n)%TypeIndx = 2
                       SS(n)%RateType = 'COM2'
!
! -------------------------
! ................. Source/sink fluid component is COM3 = ...
! -------------------------
!
                    CASE('COM3')
                       SS(n)%TypeIndx = 3
                       SS(n)%RateType = 'COM3'
!
! -------------------------
! ................. Source/sink fluid component is COM4 = ...
! -------------------------
!
                    CASE('COM4')
                       SS(n)%TypeIndx = 4
                       SS(n)%RateType = 'COM4'
!
! -------------------------
! ................. Source/sink fluid component is COM5 = ...
! -------------------------
!
                    CASE('COM5')
                       SS(n)%TypeIndx = 5
                       SS(n)%RateType = 'COM5'
!
! -------------------------
! ................. Source/sink fluid component is COM6 = ...
! -------------------------
!
                    CASE('COM6')
                       SS(n)%TypeIndx = 6
                       SS(n)%RateType = 'COM6'
!
! -------------------------
! ................. Source/sink fluid component is HEAT = ...
! -------------------------
!
                    CASE('HEAT')
                       SS(n)%TypeIndx = 2
                       SS(n)%RateType = 'HEAT'
!
! -------------------------
! ................. Default case: No available option
! -------------------------
!
                    CASE DEFAULT
!
                       SS(n)%TypeIndx = 0
!
                       WRITE(*,6005) SS_Typ,SS(n)%ElemName(1:5),SS(n)%name  !       Print a warning
                       SS(n)%TypeIndx = 1                                   !       Assign the COM1 option
                       SS(n)%rate     = 0.0d0                               !       Set the rate at zero
                       CYCLE                                                !       Continue the source/sink loop
!
                    END SELECT SELECT_SS
!
! <<<
! <<< End of the SOURCE/SINK LOOP
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
 5001 FORMAT(2A5,25X,A5,3G10.4)
!
 5003 FORMAT(16I5)
!
 5006 FORMAT(4E14.7)
!
 5007 FORMAT(A6)
!
 6000 FORMAT(/,'Read_Generic_File_GENER',T50,'[v1.0,   8 March      2007]',/, &
     &         ':::::   Read the sink/source data from file GENER (generic)')
!
!
 6001 FORMAT(' Number of sinks/sources specified in data block <GENER> exceeds allowable maximum of ',I4,/,  &
     &       ' Increase parameter <Max_NumSS> in the main program and recompile',//,   &
     &       ' >>>>>>>>>>  S K I P   F L O W   S I M U L A T I O N   <<<<<<<<<<')
 6002 FORMAT('+++  ')
!
 6004 FORMAT(' Reference to unknown element <',A8,'> at source ',I2,' ==> Will ignore source')
 6005 FORMAT(' Ignore unknown generation option <',A4,'> at element <',A8,'>, source <',A5,'>')
!
 6103 FORMAT(T2,'Memory allocation of <SS(n)%TimeColumn>, <SS(n)%RateColumn> in ', &
     &          'subroutine <Read_Generic_File_GENER> was successful')
 6104 FORMAT(//,22('ERROR-'),//,   &
     &       T2,'Memory allocation of <SS(n)%TimeColumn>, <SS(n)%RateColumn> in ', &
     &          'subroutine <Read_Generic_File_GENER> was unsuccessful',//,  &
     &       T2,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6105 FORMAT(T2,'Memory allocation of <SS(n)%EnthColumn> in subroutine <Read_Generic_File_GENER> was successful')
 6106 FORMAT(//,22('ERROR-'),//,   &
     &       T2,'Memory allocation of <SS(n)%EnthColumn> in subroutine <Read_Generic_File_GENER> was unsuccessful',//,   &
     &       T2,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End of Read_File_GENER
!
!
      RETURN
!
      END SUBROUTINE Read_File_GENER
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE Read_Initial_Conditions_File
!
! ...... Modules to be used
!
         USE EOS_Default_Parameters
!
         USE Basic_Parameters
         USE General_External_File_Units
         USE General_Control_Parameters
!
         USE Grid_Geometry, ONLY: elem
         USE Element_Attributes
         USE Solution_Matrix_Arrays, ONLY: X,R
!
         USE Geologic_Media_Parameters
!
         USE Sources_and_Sinks
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*     ROUTINE FOR READING THE DATA DESCRIBING THE DOMAIN ELEMENTS     *
!*                   DIRECTLY FROM THE FILE "INCON"                    *
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
      REAL(KIND = 8) :: origin_of_time, time_to_this_point
      REAL(KIND = 8) :: phi
!
      REAL(KIND = 8), DIMENSION(10) :: PV, primary_variable_values
!
! -------
! ... Double precision arrays
! -------
!
      REAL(KIND = 8), DIMENSION(3) :: perm
!
! -------
! ... Integer variables
! -------
!
      INTEGER           :: i, j, n, NMX, nloc, JLOC, No_UnkN, ier
      INTEGER(KIND = 1) :: StIndx
!
      INTEGER :: timesteps_to_this_point, NR_iterations_to_this_point
!
! -------
! ... Character variables
! -------
!
      CHARACTER(LEN = 3) :: state
      CHARACTER(LEN = 5) :: DENT, name
!
      CHARACTER(LEN = 35) :: PVList, primary_variable_names
!
! -------
! ... Logical variables
! -------
!
      LOGICAL :: end_of_subblock_reached = .FALSE.
!
! -------
! ... Namelist 1: Continuation data
! -------
!
      NAMELIST/ Data_For_Continuation_Run / timesteps_to_this_point,      &    ! Number of total timesteps up to the beginning of the continuation run
     &                                      NR_iterations_to_this_point,  &    ! Number of total NR iterations up to the beginning of the continuation run
     &                                      origin_of_time,               &    ! Time at beginning of simulation
     &                                      time_to_this_point                 ! Time at beginning of continuation run
!
! -------
! ... Namelist 2: Uniform initial conditions
! -------
!
      NAMELIST/Uniform_Initial_Conditions/ state,                   &
     &                                     primary_variable_names,  &
     &                                     primary_variable_values
!
! -------
! ... Namelist 3: Element-specific initial conditions
! -------
!
      NAMELIST/InCond/ name,    &
     &                 state,   &
     &                 PVList,  &
     &                 PV,      &
     &                 phi,     &
     &                 perm
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of Read_Initial_Conditions_File
!
!
      WRITE(VERS_Unit,6000)
!
!
      No_UnkN = 0
!
! ----------
! ... Reading the headings of the file <INCON>
! ----------
!
      REWIND (UNIT = INCON_Unit)
      READ(INCON_Unit,5001) DENT
!
!***********************************************************************
!*                                                                     *
!*     Initialize generic hydraulic properties & related options       *
!*                                                                     *
!***********************************************************************
!
!
      ElemMedia(1:NumElemTot,current)%porosity = media(elem(1:NumElemTot)%MatNum)%Poros      ! Assign porosities
!
!**********************************************************************
!*****Parallelize Initial generic hydraulic properties and related
!*****options                                                       ***
!**********************************************************************
      DO i = 1,3
         ElemMedia(1:NumElemTot,current)%perm(i) = media(elem(1:NumElemTot)%MatNum)%Perm(i)  ! Assign permeabilities
      END DO
!
!
! ... Determine whether porosity and permeability will be varying as functions of P, T
!
      IF(ANY(media(1:NumMedia)%Expan /= 0.0d0) .OR. ANY(media(1:NumMedia)%Compr /= 0.0d0)) THEN
         variable_porosity = .TRUE.
      END IF
!
!***********************************************************************
!*                                                                     *
!*         Read generic initial conditions from file <INCON>           *
!*                                                                     *
!***********************************************************************
!
!
      READ(UNIT = INCON_Unit, NML = Uniform_Initial_Conditions, IOSTAT = ier)
!
!--------
! ... Stop if there is a problem reading the namelist
!--------
!
      IF (ier /= 0) THEN
         WRITE (*, FMT = 6100)
         STOP
      END IF
!
!--------
! ... Checking validity of state
!--------
!
      IF( ALL(state_name(1:Number_of_states) /= state) ) THEN
         WRITE (*, FMT = 6105) state
         STOP
      END IF
!
      DO i=1,Number_of_states
         IF(state == state_name(i)) THEN
            DefaultStateIndex = i
            EXIT
         END IF
      END DO
!
!--------
! ... Assigning initial conditions
!--------
!
      default_initial_cond(1:NumComPlus1) = primary_variable_values(1:NumComPlus1)
!
!***********************************************************************
!*                                                                     *
!*                 Assign generic initial conditions                   *
!*                                                                     *
!***********************************************************************
!
!Cannot be parallelized
      DO_NumEle: DO n=1,NumElemTot
!
         nloc = Locp(n)      ! Determine the location in the arrays
!
! -----------
! ...... Assigning uniform initial conditions
! ----------
!
         ElemState%index(n,current)  = DefaultStateIndex                    ! Assign initial state index
         X(nloc+1 : nloc+NumComPlus1) = default_initial_cond(1:NumComPlus1)  ! Assign initial conditions: Whole array operation !!!
!
      END DO DO_NumEle
!
!
!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>                                                                     >
!>           ELEMENT LOOP FOR ASSIGNING INITIAL CONDITIONS             >
!>                                                                     >
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
      DO_NumEle1: DO n=1,100000000
!
         READ(UNIT = INCON_Unit, FMT = *) name
!
! -------------
! ...... If EL is blank, the element list is exhausted
! -------------
!
         IF(name(1:3) == '<<<') THEN
!
            IF(No_UnkN > 0) WRITE(*,6001) No_UnkN
!
            EXIT DO_NumEle1
!
         END IF
!
! -------------
! ...... If EL = ':::  ', read restart data using a namelist
! -------------
!
         IF(name(1:3) == ':::' .AND. name(4:5) /= '>>' .AND. name(4:5) /= '<<') THEN
!
            READ(UNIT = INCON_Unit, NML = Data_For_Continuation_Run, IOSTAT = ier)
!
! ......... Stop if there is a problem reading the data block end specifies
!
            IF( ier /= 0 ) THEN
               WRITE (*, FMT = 6500)
               STOP
            END IF
!
! ......... Assign parameter values
!
            Tot_NumTimeSteps    = timesteps_to_this_point
            Tot_NumNRIterations = NR_iterations_to_this_point
            TimeOrigin          = origin_of_time
            TimeShift           = time_to_this_point
!
            IF(SimulationTimeEnd /= 0.0D0 .AND. TimeShift >= SimulationTimeEnd) TimeShift = 0.0D0
!
! ......... Reset the restart file indicator flag
!
            RestartSimulation = .TRUE.
!
            EXIT DO_NumEle1                          ! Done; exit the loop
!
         END IF
!
! -------------
! ...... If <name> = ':::>>', skip to the beginning
! -------------
!
         IF(name == ':::>>') THEN
            IF(end_of_subblock_reached) THEN
               WRITE(*, FMT = 6205)
               STOP
            ELSE
               CYCLE DO_NumEle1
            END IF
         END IF
!
! -------------
! ...... If <name> = ':::<<', continue reading to reach the identifiers of end or continuation data
! -------------
!
         IF(name == ':::<<') THEN
            end_of_subblock_reached = .TRUE.
            CYCLE DO_NumEle1
         END IF
!
! -------------
! ...... If a name is read after having reached the end of the subblock, then there is a mistake
! -------------
!
         IF(end_of_subblock_reached) THEN
            WRITE(*, FMT = 6210)
            STOP
         END IF
!
! -------------
! ...... For a regular element name, begin reading element-by-element init
! -------------
!
         state  = '   '   ! Element initial thermodynamic state (= phase coexistence)
         PV     = 0.0d0   ! Values of the initial primary variable values in this element
         phi    = 0.0d0   ! Element initial uniform porosity
         perm   = 0.0d0   ! Element initial uniform permeabilities
!
! -------------
! ...... Read the element data from the INCON file
! -------------
!
         READ(UNIT = INCON_Unit, NML = InCond, IOSTAT = ier)
!
!-----------
! ...... Stop if there is a problem reading the namelist
!-----------
!
         IF (ier /= 0) THEN
            WRITE (*, FMT = 6110)
            STOP
         END IF
!
!-----------
! ...... Checking validity of state
!-----------
!
         IF( ALL(state_name(1:Number_of_states) /= state) ) THEN
            WRITE (*, FMT = 6115) name, state
            STOP
         END IF
!
! -------------
! ...... Determine the state index
! -------------
!
         DO_State: DO j = 1,Number_of_states
!
! ......... If match, assign the state index number
!
            IF(state == State_name(j)) THEN
               StIndx = j
               EXIT DO_State
            END IF
!
         END DO DO_State
!
! -------------
! ...... Assign the initial conditions - primary variables
! -------------
!
         DO_NumEle2: DO j=1,NumElemTot
!
            IF(elem(j)%name(1:5) /= name) CYCLE DO_NumEle2
!
            jloc = Locp(j)                                   ! Determine pointer
!
            IF(phi /= 0.0d0) THEN
               ElemMedia(j,current)%porosity = phi           ! Override porosity
            END IF
!
            IF(ANY(perm(1:3) /= 0.0d0)) THEN
               ElemMedia(j,current)%perm(1:3) = perm(1:3)    ! Override permeability
            END IF
!
            ElemState%index(j,current)  = StIndx            ! Assign state index
!
            X(jloc+1 : jloc+NumComPlus1) = PV(1:NumComPlus1) ! Assign/override initial conditions - Whole array operation
!
            GO TO 1000
!
         END DO DO_NumEle2
!
! -------------
! ...... If the element number is unknown, print message & continue
! -------------
!
         No_UnkN = No_UnkN+1                    ! Count the unknown elements
         IF(No_UnkN <= 50) WRITE(*,6008) name   ! Print warning message
!
! <<<
! <<< End of the ELEMENT LOOP
! <<<
!
 1000 END DO DO_NumEle1
!
!
!***********************************************************************
!*                                                                     *
!*           Store the initial/original state of the system            *
!*                                                                     *
!***********************************************************************
!
!
      ElemState%index(:,previous) = ElemState%index(:,current)                                            ! CAREFUL! Whole array operation
!
      ElemMedia(1:NumElemTot,original)%porosity = ElemMedia(1:NumElemTot,current)%porosity
      ElemMedia(1:NumElemTot,previous)%porosity = ElemMedia(1:NumElemTot,current)%porosity
!
      FORALL (i=1:3)
         ElemMedia(1:NumElemTot,original)%perm(i) = ElemMedia(1:NumElemTot,current)%perm(i)
         ElemMedia(1:NumElemTot,previous)%perm(i) = ElemMedia(1:NumElemTot,current)%perm(i)
      END FORALL
!
      DO j = 1,NumEqu
         ElemMedia(1:NumElemTot,j)%porosity               = ElemMedia(1:NumElemTot,current)%porosity  ! Initialize the incremented-state porosity
         FORALL (i=1:3) ElemMedia(1:NumElemTot,j)%perm(i) = ElemMedia(1:NumElemTot,current)%perm(i)   ! Initialize the incremented-state permeability
      END DO
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 5001 FORMAT(A5)
!
 5002 FORMAT(A5,10X,E15.8,2x,a3,36x,3(E15.8))
 5005 FORMAT(6E20.13)
 5008 FORMAT(2I5,8E15.8)
!
 6000 FORMAT(/,'Read_Initial_Conditions_File',T50,'[v1.0,   1 October   2007]',/,        &
     &         ':::::   Read the initial condition data directly from the file <INCON>')
!
 6001 FORMAT(/,' !!!!! WARNING !!!!!! ==>  INCON data at ',I5,' unknown elements have been ignored')
!
 6008 FORMAT(' WILL IGNORE INITIAL CONDITION AT UNKNOWN ELEMENT ',A8)
!
 6100 FORMAT(//,22('ERROR-'),//,T33,                                                                                               &
     &             '       S I M U L A T I O N   A B O R T E D',/,                                                                 &
     &         T22,'There is a problem reading the namelist <Uniform_Initial_Conditions> in subroutine <Read_Initial_Conditions> ',&
     &       /,T32,'               CORRECT AND TRY AGAIN',                                                                         &
     &       //,22('ERROR-'))
!
 6105 FORMAT(//,22('ERROR-'),//,T43,                                                            &
     &             'S I M U L A T I O N   A B O R T E D',/,                                     &
     &         T25,'The state identifier of uniform initial conditions <state> = "',a3,'"',/,   &
     &         T31,'This is not an available option for this Equation of State',                &
     &       /,T50,'CORRECT AND TRY AGAIN',                                                     &
     &       //,22('ERROR-'))
!
 6106 FORMAT('<<<END of INITIAL_CONDITIONS')
 6108 FORMAT(':::<<<END of Element-by-element initial conditions')
 !
 6109 FORMAT(//,22('ERROR-'),//,T33,                                                                                               &
     &             '       S I M U L A T I O N   A B O R T E D',/,                                                                 &
     &         T22,'There is a blank line within the element-specific namelist <InCond> in subroutine <Read_Initial_Conditions> ', &
     &       /,T32,'               CORRECT AND TRY AGAIN',                                                                         &
     &       //,22('ERROR-'))
!
 6110 FORMAT(//,22('ERROR-'),//,T33,                                                                                               &
     &             '       S I M U L A T I O N   A B O R T E D',/,                                                                 &
     &         T22,'There is a problem reading the element-specific namelist <InCond> in subroutine <Read_Initial_Conditions> ',   &
     &       /,T32,'               CORRECT AND TRY AGAIN',                                                                         &
     &       //,22('ERROR-'))
!
 6115 FORMAT(//,22('ERROR-'),//,T43,                                                                         &
     &             'S I M U L A T I O N   A B O R T E D',/,                                                  &
     &         T25,'In element "',a5,'", the state identifier of initial conditions <state> = "',a3,'"',/,   &
     &         T31,'This is not an available option for this Equation of State',   &
     &       /,T50,'CORRECT AND TRY AGAIN',   &
     &       //,22('ERROR-'))
!
 6205 FORMAT(//,22('ERROR-'),//,T33,                                                                          &
     &             '       S I M U L A T I O N   A B O R T E D',/,                                            &
     &         T22,'There is a sub-block initialization identifier after an end-block identifier ":::<<" ',   &
     &       /,T32,'               CORRECT AND TRY AGAIN',                                                    &
     &       //,22('ERROR-'))
!
 6210 FORMAT(//,22('ERROR-'),//,T33,                                                                                 &
     &             '       S I M U L A T I O N   A B O R T E D',/,                                                   &
     &         T22,'There is an attempt to provide initial conditions to an after an end-block identifier ":::<<"',  &
     &       /,T32,'               CORRECT AND TRY AGAIN',                                                           &
     &       //,22('ERROR-'))
!
 6500 FORMAT(//,22('ERROR-'),//,T33,                                                                               &
     &             '       S I M U L A T I O N   A B O R T E D',/,                                                 &
     &         T22,'There is a problem reading the namelist <Data_For_Continuation_Run> in the <INCON> data file', &
     &       /,T32,'               CORRECT AND TRY AGAIN',                                                         &
     &       //,22('ERROR-'))
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End of Read_Initial_Conditions_File
!
!
      RETURN
!
      END SUBROUTINE Read_Initial_Conditions_File
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
