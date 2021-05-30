!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>                                                                     >
!>                                                                     >
!>   GAS_Specifics.f95: Code unit of gas-specific routines,            >
!>           routines overriding FTSim generic routines,               >
!>                     and associated interfaces                       >
!>                                                                     >
!>                                                                     >
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  Segment that includes procedures specific to the gas EOS simulation,
!  such as the reading of gas-specific inputs and the preparation of
!  gas-specific output files.  Generic procedures and operator
!  extension - which override (overload) the standard procedures used by
!  FTSim for the simulation of the most-common problems - are defined in
!  this segment, which does not include any procedures describing the EOS
!
!
      MODULE EOS_Routine_Selector
!
         SAVE
!
!***********************************************************************
!*                                                                     *
!*       Interface: selection of appropriate routine for reading       *
!*                  the EOS-specific data                              *
!*                                                                     *
!***********************************************************************
!
         INTERFACE READ_EOS_Specific_Data
!
! -------------
! ......... The GAS EOS-specific routine
! -------------
!
            SUBROUTINE Read_Gas_Specific_Data
!
               USE General_Control_Parameters
               USE Gas_Parameters
               USE General_External_File_Units
!
               IMPLICIT NONE
!
               REAL(KIND = 8) :: gas_reference_density, gas_compressibility
               REAL(KIND = 8) :: P_at_reference_density, T_at_reference_density
               REAL(KIND = 8) :: gas_reference_specific_heat, gas_reference_thermal_cond
!
               REAL(KIND = 8), DIMENSION(10)   :: Cp_coefficients, Cv_coefficients, viscosity_coefficients, ThermCond_coefficients
!
               INTEGER :: ier, VisCoeffNumber, ThermCondCoeffNumber
!
               CHARACTER (LEN = 20) gas_name
!
            END SUBROUTINE Read_Gas_Specific_Data
!
         END INTERFACE
!
!
      END MODULE EOS_Routine_Selector
!
!
!
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!
!
!
      SUBROUTINE Read_Gas_Specific_Data
!
! ...... Modules to be used
!
         USE General_Control_Parameters
!
         USE Gas_Parameters
!
         USE General_External_File_Units
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*           ROUTINE FOR READING THE GAS PROPERTY DATA                 *
!*           DIRECTLY FROM THE FTSim INPUT DATA BLOCK                  *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
      IMPLICIT NONE
!
      REAL(KIND = 8) :: gas_reference_density
      REAL(KIND = 8) :: gas_compressibility
      REAL(KIND = 8) :: gas_reference_viscosity
      REAL(KIND = 8) :: P_at_reference_density
      REAL(KIND = 8) :: T_at_reference_density
      REAL(KIND = 8) :: gas_reference_specific_heat
      REAL(KIND = 8) :: gas_reference_thermal_cond
!
      REAL(KIND = 8), DIMENSION(10)  :: viscosity_coefficients, ThermCond_coefficients
      REAL(KIND = 8), DIMENSION(0:4) :: Cp_coefficients
      REAL(KIND = 8) :: TCrit, PCrit, VCrit, ZCrit, omega, MolWeight
!
      INTEGER :: ier, VisCoeffNumber, ThermCondCoeffNumber
!
      CHARACTER(LEN = 20) gas_name
      CHARACTER(LEN = 3)  DataBlockEndSpecifier
!
! -------
! ... Namelists
! -------
!
      NAMELIST/Gas_EOS_Specific_Data/ gas_name, gas_reference_density, gas_compressibility,    &
     &                                P_at_reference_density, T_at_reference_density,          &
     &                                gas_reference_viscosity, viscosity_coefficients,         &
     &                                ThermCond_coefficients, VisCoeffNumber,                  &
     &                                ThermCondCoeffNumber, Cp_coefficients,                   &
     &                                gas_reference_specific_heat, gas_reference_thermal_cond, &
     &                                TCrit, PCrit, VCrit, ZCrit, Omega, MolWeight
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of Read_Gas_Specific_Data
!
!
      WRITE(UNIT = VERS_Unit, FMT = 6000)
!
!
!***********************************************************************
!*                                                                     *
!*                          Initializations                            *
!*                                                                     *
!***********************************************************************
!
      gas_name = 'GAS'
!
      MolWeight = 0.0d0
      omega      = 0.0d0
!
! ... Density-related properties and parameters
!
      gas_reference_density  = 0.0d0
      P_at_reference_density = 0.0d0
      T_at_reference_density = 0.0d0
      gas_compressibility    = 0.0d0
!
! ... Viscosity-related properties and parameters
!
      gas_reference_viscosity = 0.0d0
!
      constant_gas_viscosity  = .TRUE.
      viscosity_coefficients  = 0.0d0
      VisCoeffNumber          = 0
!
! ... Thermal properties and parameters
!
      gas_reference_specific_heat        = 0.0d0
      gas_reference_thermal_cond         = 0.0d0
      ThermCondCoeffNumber               = 0
      Cp_coefficients                    = 0.0d0
      ThermCond_coefficients             = 0.0d0
!
! ... Critical properties and parameters
!
      PCrit                       = 0.0d0
      TCrit                       = 0.0d0
      VCrit                       = 0.0d0
      ZCrit                       = 0.0d0
!
!***********************************************************************
!*                                                                     *
!*        READING THE BASIC DATA DESCRIBING THE GAS PROPERTIES         *
!*                                                                     *
!***********************************************************************
!
      READ(*, NML = Gas_EOS_Specific_Data, IOSTAT = ier)        ! Reading the properties of the porous medium
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
!--------
! ... check for unrealistic critical properties
!--------
!
      IF((PCrit == 0.0d0) .OR. (TCrit == 0.0d0) .OR. (VCrit == 0.0d0) .OR. (ZCrit == 0.0d0)) THEN
         WRITE (*, FMT = 6110)
         STOP
      END IF
!
!--------
! ... check for unrealistic critical properties
!--------
!
      IF(omega <= 0.0d0 .OR. MolWeight <= 0.0d0) THEN
         WRITE (*, FMT = 6120)
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
      gas%name             = gas_name
!
      gas%PresAtRefDensity = P_at_reference_density
      gas%TempAtRefDensity = T_at_reference_density
      gas%RefDensity       = gas_reference_density
      gas%RefViscosity     = gas_reference_viscosity
      gas%compressibility  = gas_compressibility
      gas%RefThermalCond   = gas_reference_thermal_cond
      gas%PCrit            = PCrit
      gas%TCrit            = TCrit
      gas%VCrit            = VCrit
      gas%ZCrit            = ZCrit
      gas%omega            = omega
      gas%MolWeight        = MolWeight
!
      gas%Cp_coefficients  = Cp_coefficients
!
!--------
! ... determine viscosity status
!--------
!
      gas%RefViscosity  = gas_reference_viscosity
!
      IF(VisCoeffNumber == 0 .OR. gas_reference_viscosity == 0.0d0) THEN
         WRITE (*, FMT = 6130)
         STOP
      END IF
!
      IF(ANY(viscosity_coefficients /= 0.0d0)) THEN
         constant_gas_viscosity     = .FALSE.
         ALLOCATE(gas%viscosity_coefficients(1:VisCoeffNumber))
         gas%viscosity_coefficients = viscosity_coefficients       ! whole-array operation
      END IF
!
!--------
! ... determine gas thermal conductivity status
!--------
!
      gas%RefThermalCond = gas_reference_thermal_cond
!
      IF(ThermCondCoeffNumber <= 0 .AND. gas_reference_thermal_cond <= 0.0d0) THEN
         WRITE (*, FMT = 6140)
         STOP
      END IF
!
      IF(ANY(ThermCond_coefficients /= 0.0d0)) THEN
         constant_gas_thermal_cond = .FALSE.
         ALLOCATE(gas%ThermCond_coefficients(1:ThermCondCoeffNumber))
         gas%ThermCond_coefficients = ThermCond_coefficients ! whole-array operation
      END IF
!
!--------
! ... Write gas properties to standard output
!--------
!
      WRITE (*, FMT = 6390) gas%name
      WRITE (*, FMT = 6400) gas%RefDensity, gas%PresAtRefDensity, gas%TempAtRefDensity,      &
     &                      gas%compressibility, gas%RefViscosity, gas%RefThermalCond
!
!--------
! ... Read past the end-of-the-data block
!--------
!
      READ( *, FMT = '(A3)', IOSTAT = ier ) DataBlockEndSpecifier
!
!-----------
! ...... Stop if there is a problem reading the data block end specifies
!-----------
!
      IF_FormatStyle: IF( (ier == 0) .AND. (DataBlockEndSpecifier == '<<<') ) THEN
         RETURN
      ELSE
         WRITE (*, FMT = 6500)
         STOP
      END IF IF_FormatStyle
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 6000 FORMAT(/,'Read_Gas_Specific_Data',T50,'[v1.0,   29 May      2008]',/,  &
     &         ':::::   Read data specific to the Air-H2O system from the <Air+H2O> data block of the input data file')
!
 6100 FORMAT(//,22('ERROR-'),//,T33,  &
     &             '       S I M U L A T I O N   A B O R T E D',/,                                             &
     &         T22,'There is a problem reading the namelist <Gas_EOS_Specific_Data> in the <Gas> data block',  &
     &       /,T32,'               CORRECT AND TRY AGAIN',                                                     &
     &       //,22('ERROR-'))
!
 6110 FORMAT(//,22('ERROR-'),//,T33,  &
     &             '       S I M U L A T I O N   A B O R T E D',/,                                             &
     &         T22,'One or more of the GAS critical properties in namelist <Gas_EOS_Specific_Data> is zero',   &
     &       /,T32,'               CORRECT AND TRY AGAIN',                                                     &
     &       //,22('ERROR-'))
!
 6120 FORMAT(//,22('ERROR-'),//,T33,  &
     &             '       S I M U L A T I O N   A B O R T E D',/,                                                &
     &         T22,'Either the gas molecular weight or omega in namelist <Gas_EOS_Specific_Data> is unrealistic', &
     &       /,T32,'               CORRECT AND TRY AGAIN',                                                        &
     &       //,22('ERROR-'))
!
 6130 FORMAT(//,22('ERROR-'),//,T33,  &
     &             '       S I M U L A T I O N   A B O R T E D',/,                                                    &
     &         T22,'Improper definition of gas viscosity and coefficients in <Gas_EOS_Specific_Data> is unrealistic', &
     &       /,T32,'               CORRECT AND TRY AGAIN',                                                            &
     &       //,22('ERROR-'))
!
 6140 FORMAT(//,22('ERROR-'),//,T33,  &
     &             '       S I M U L A T I O N   A B O R T E D',/,                                                       &
     &         T22,'Improper definition of gas therm. cond. and coefficients in <Gas_EOS_Specific_Data> is unrealistic', &
     &       /,T32,'               CORRECT AND TRY AGAIN',                                                               &
     &       //,22('ERROR-'))
!
 6390 FORMAT(/,/,                                                                                                                                    &
     &       A15,' Reference state and properties:',/,/,                                                                                                       &
     &       T5,'Density           Pres                Temp                Compressibility     Viscosity           Specific_heat       Thermal_cond',&
     &       /,T5,'(kg/m^3)          (Pa)                (C)                 (1/Pa)             (Pa.s)               (J/kg/C)            (W/m/C)')
!
 6400 FORMAT(F10.2,6(10X,ES10.2E2))
!
 6500 FORMAT(//,22('ERROR-'),//,T33,  &
     &             '       S I M U L A T I O N   A B O R T E D',/,                                     &
     &         T10,'There is a problem reading the data block end specifier <DataBlockEndSpecifier> ', &
     &             'in the <GAS_EOS> block',                                                           &
     &       /,T32,'               CORRECT AND TRY AGAIN',                                             &
     &       //,22('ERROR-'))
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End of Read_Gas_Specific_Data
!
!
      RETURN
!
      END SUBROUTINE Read_Gas_Specific_Data
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE Write_EOS_Specific_Info
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
!*      PRINTING ADDITIONAL EOS-SPECIFIC INFO IN THE OUTPUT FILE       *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
      IMPLICIT NONE
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of Write_EOS_Specific_Info
!
!
      WRITE(VERS_Unit,6000)
!
! ... Format
!
 6000 FORMAT(/,'Write_EOS_Specific_Info',T50,'[v1.0,   29 May      2008]',/,   &
     &         ':::::   Routine for printing additional EOS-specific information in the standard <FTSim> output file')
!
!
!
      WRITE(*,6022)
!
! ... Format
!
 6022 FORMAT(' ',/,          &
     &       ' ',139('*'),/, &
     &       '     * FTSim+Gas: EQUATION OF STATE FOR PURE GAS IN A SINGLE PHASE (Gas)',T140,'*',/,   &
     &            ' ',139('*'))
!
!
!
      WRITE(*,6026) NumCom, NumEqu, NumPhases, NumMobPhases
!
! ... Format
!
 6026 FORMAT(/,' OPTIONS SELECTED ARE: (NumCom, NumEqu, NumPhases, NumMobPhases) = (',3(I1,','),I1,')',/)
!
!
!
      WRITE(*,6028) NumCom, NumEqu, NumPhases, NumMobPhases
!
! ... Format
!
 6028 FORMAT(25X,'NumCom           = ',I2,' - NUMBER OF FLUID COMPONENTS',/,   &
     &       25X,'NumEqu           = ',I2,' - NUMBER OF EQUATIONS PER GRID BLOCK',/,   &
     &       25X,'NumPhases        = ',I2,' - NUMBER OF POSSIBLE PHASES ',/,   &
     &       25X,'NumMobPhases     = ',I2,' - NUMBER OF MOBILE PHASES')
!
!
!
      WRITE(*,6030)
!
! ... Format
!
 6030 FORMAT(/,' AVAILABLE OPTIONS ARE: (NumCom, NumEqu, NumPhases, NumMobPhases) = (1,1,1,1)- GAS, ISOTHERMAL')
!
!
!
      WRITE(*,6032)
!
! ... Format
!
 6032 FORMAT(/,' THE POSSIBLE PRIMARY VARIABLES ARE:',/,/,   &
     &         23X,'     P:   PRESSURE (Pa)',/   &
     &         23X,'     T:   TEMPERATURE (C)',/)
!
!
!
!      WRITE(*,6033)
!
! ... Format
!
 6033 FORMAT(/,132('='))
!
!
!
!      WRITE(*,6034)
!
! ... Format
!
 6034 FORMAT(/,' FOR (NumCom,NumEqu,NumPhases) = (2,3,4) - EQUILIBRIUM REACTION',//,   &
     &   ' ******************************  ******************************  ',  &
     &   ' *********************************************************',/,       &
     &   ' *         COMPONENTS         *  *           PHASES           *  ',  &
     &   ' *   FLUID PHASE CONDITIONS       PRIMARY VARIABLES      *',/,       &
     &   ' ******************************  ******************************  ',  &
     &   ' *********************************************************',/,       &
     &   ' *                            *  *                            *  ',  &
     &   ' *                                 X1    X2    Index     *',/,       &
     &   ' *       # 1  -  GAS          *  *       # 1  -  Gas      *  ',  &
     &   ' *                                --------------------   *',/,       &
     &   ' *                            *  *                            *  ',  &
     &   ' *                                                       *',/,       &
     &   ' ******************************  ******************************  ',  &
     &   ' *   SINGLE-PHASE Gas          P,    T,     GAS      *',/,       &
     &   '                                                                 ',  &
     &   ' *                                                       *',/,       &
     &   '                                                                 ',  &
     &   ' *********************************************************')
!
!
!
       WRITE(*,6040)
!
! ... Format
!
 6040 FORMAT(123('*'))
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End of Write_EOS_Specific_Info
!
!
      RETURN
!
      END SUBROUTINE Write_EOS_Specific_Info
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE Print_Standard_Output
!
! ...... Modules to be used
!
         USE Universal_Parameters
         USE Basic_Parameters
!
         USE EOS_Default_Parameters
         USE General_External_File_Units
         USE General_Control_Parameters
!
         USE Grid_Geometry
         USE Element_Attributes
         USE Connection_Attributes
         USE Sources_and_Sinks
!
         USE Solution_Matrix_Arrays, ONLY: X, DX
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*      ROUTINE FOR GENERATING PRINT-OUT FOR THE Gas SYSTEM            *
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
      REAL(KIND = 8), DIMENSION(4) :: DXM
!
! -------
! ... Double precision variables
! -------
!
      REAL(KIND = 8) :: DAY, GC, HGC
!
! -------
! ... Integer variables
! -------
!
      INTEGER :: n,NLOC,i,j
!
! -------
! ... Character variables
! -------
!
      CHARACTER(LEN = 1) :: HB = ' ', H0 = ' '
!
      CHARACTER(LEN = 1) :: ConvInd
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
      SAVE First_call,HB,H0
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of Print_Standard_Output
!
!
      IF_1stCall: IF(First_call) THEN
!
         First_call = .FALSE.
         WRITE(VERS_Unit,6000)
!
         IF_Option_OutputAmount: IF(Option_OutputFormat == 'BOTH') THEN
!
            IF(coordinate_system == 'cyl' .OR. coordinate_system == 'CYL' .OR. coordinate_system == 'Cyl') THEN
               WRITE(Plot_Elem_Unit,6101)
            ELSE
               WRITE(Plot_Elem_Unit,6102)
            END IF
!
         ELSE IF(Option_OutputFormat == 'PLOT') THEN
!
            IF(coordinate_system == 'cyl' .OR. coordinate_system == 'CYL' .OR. coordinate_system == 'Cyl') THEN
               WRITE(Plot_Elem_Unit,6101)
            ELSE
               WRITE(Plot_Elem_Unit,6102)
            END IF
!
         END IF IF_Option_OutputAmount
!
      END IF IF_1stCall
!
!
      DAY = time/8.64d4
!
! -------
! ... Option for printing short <FTSim> outputs
! -------
!
      IF_output: IF(Option_OutputFormat == 'BOTH') THEN
!
         WRITE(Plot_Elem_Unit,6105) time, NumElemTot
         WRITE(Plot_Conx_Unit,6115) time
!
      ELSE IF(Option_OutputFormat == 'PLOT') THEN
!
         WRITE(Plot_Elem_Unit,6105) time, NumElemTot
         WRITE(Plot_Conx_Unit,6115) time
         GO TO 101
!
      END IF IF_output
!
!***********************************************************************
!*                                                                     *
!*                       COMPUTE MAXIMUM CHANGES                       *
!*                                                                     *
!***********************************************************************
!
      DXM(1) = MAXVAL(ABS(DX(1 : (NumElemTot-1)*NumComPlus1+1 : NumComPlus1)))
      DXM(2) = MAXVAL(ABS(DX(2 : (NumElemTot-1)*NumComPlus1+2 : NumComPlus1)))
      DXM(3) = MAXVAL(ABS(DX(3 : (NumElemTot-1)*NumComPlus1+3 : NumComPlus1)))
!
!***********************************************************************
!*                                                                     *
!*                    PRINT HEADER INFORMATION - 1                     *
!*                                                                     *
!***********************************************************************
!
      WRITE(*, 6002) H0, title, Tot_NumTimeSteps, Tot_NumNRIterations, DAY
      WRITE(*, 6004) HB
!
      IF(convergence .EQV. .TRUE.) THEN
         ConvInd = 'Y'
      ELSE
         ConvInd = 'N'
      END IF
!
      WRITE(*, 6005)
      WRITE(*, 6006) time, Tot_NumTimeSteps, Tot_NumNRIterations, ConvInd, (DXM(I),I=1,3), &
     &               MaxResidual, MaxResElemNum, MaxResEquNum, TimeStep
!
      WRITE(*, 6004) HB
      WRITE(*, 6008) H0
      WRITE(*, 6010)
!
!
!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>                                                                     >
!>                         ELEMENT LOOP - 1                            >
!>             Print primary variables & related information           >
!>                                                                     >
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
 101 DO_NumEle1: DO n=1,NumElemTot
!
! ----------
! ...... Print standard output
! ----------
!
      IF_PrintOptions: IF(Option_OutputFormat == 'STAN') THEN
!
         IF_act1: IF(elem(n)%activity /= 'I' .AND. elem(n)%activity /= 'V') THEN     ! Active elements
!
!
            WRITE(*, 6012) ADJUSTR(elem(n)%name),N,               &
     &                     ElemState%pres(n,0),                  &
     &                     ElemState%temp(n,0),                  &
     &                     elem(n)%coord(3),                      &
     &                     ElemProp(n,0)%density(GasPhase),       &
     &                     ElemProp(n,0)%viscosity(GasPhase),     &
     &                     ElemProp(n,0)%MassFrac(gas_2,GasPhase),  &
     &                     elem(n)%coord(1),                      &
     &                     ElemProp(n,0)%enthalpy(GasPhase),      &
     &                     ElemProp(n,0)%IntEnergy(GasPhase)
!
         ELSE
!
            WRITE(*, 6011) ADJUSTR(elem(n)%name),N,               &
     &                     ElemState%pres(n,0),                  &
     &                     ElemState%temp(n,0),                  &
     &                     elem(n)%coord(3),                      &
     &                     ElemProp(n,0)%density(GasPhase),       &
     &                     ElemProp(n,0)%viscosity(GasPhase),     &
     &                     ElemProp(n,0)%MassFrac(gas_2,GasPhase),  &
     &                     elem(n)%coord(1),                      &
     &                     ElemProp(n,0)%enthalpy(GasPhase),      &
     &                     ElemProp(n,0)%IntEnergy(GasPhase)
!
         END IF IF_act1
!
! ----------
! ...... Print standard output & plotting files
! -------
!
      ELSEIF(Option_OutputFormat == 'BOTH') THEN

         IF_act2: IF(elem(n)%activity /= 'I' .AND. elem(n)%activity /= 'V') THEN     ! Active elements
!
!
            WRITE(*, 6012) ADJUSTR(elem(n)%name),N,               &
     &                     ElemState%pres(n,0),                  &
     &                     ElemState%temp(n,0),                  &
     &                     ElemProp(n,0)%satur(GasPhase),         &
     &                     ElemProp(n,0)%density(GasPhase),       &
     &                     ElemProp(n,0)%viscosity(GasPhase),     &
     &                     ElemProp(n,0)%MassFrac(gas_2,GasPhase),  &
     &                     ElemProp(n,0)%enthalpy(GasPhase),      &
     &                     ElemProp(n,0)%IntEnergy(GasPhase)
!
         ELSE
!
            WRITE(*, 6011) ADJUSTR(elem(n)%name),N,               &
     &                     ElemState%pres(n,0),                  &
     &                     ElemState%temp(n,0),                  &
     &                     ElemProp(n,0)%satur(GasPhase),         &
     &                     ElemProp(n,0)%density(GasPhase),       &
     &                     ElemProp(n,0)%viscosity(GasPhase),     &
     &                     ElemProp(n,0)%MassFrac(gas_2,GasPhase),  &
     &                     ElemProp(n,0)%enthalpy(GasPhase),      &
     &                     ElemProp(n,0)%IntEnergy(GasPhase)
!
         END IF IF_act2
!
!
         IF_Coord1: IF(coordinate_system == 'cyl' .OR. coordinate_system == 'CYL' .OR. coordinate_system == 'Cyl') THEN
!
               WRITE(Plot_Elem_Unit,6013) ADJUSTR(elem(n)%name),N,                 &
     &                                  elem(n)%coord(1), elem(n)%coord(3),      &
     &                     ElemState%pres(n,0),                  &
     &                     ElemState%temp(n,0),                  &
     &                                  ElemProp(n,0)%satur(GasPhase),           &
     &                                  ElemProp(n,0)%density(GasPhase),         &
     &                                  ElemProp(n,0)%viscosity(GasPhase),       &
     &                                  ElemProp(n,0)%MassFrac(gas_2,GasPhase),    &
     &                                  ElemProp(n,0)%enthalpy(GasPhase),        &
     &                                  ElemProp(n,0)%IntEnergy(GasPhase)
            ELSE
!
               WRITE(Plot_Elem_Unit,6013) ADJUSTR(elem(n)%name),N,                               &
     &                                  elem(n)%coord(1), elem(n)%coord(2), elem(n)%coord(3),  &
     &                     ElemState%pres(n,0),                  &
     &                     ElemState%temp(n,0),                  &
     &                                  ElemProp(n,0)%satur(GasPhase),                         &
     &                                  ElemProp(n,0)%density(GasPhase),                       &
     &                                  ElemProp(n,0)%viscosity(GasPhase),                     &
     &                                  ElemProp(n,0)%MassFrac(gas_2,GasPhase),                  &
     &                                  ElemProp(n,0)%enthalpy(GasPhase),                      &
     &                                  ElemProp(n,0)%IntEnergy(GasPhase)
!
            END IF IF_Coord1
!
! ----------
! ...... Print plotting files only
! -------
!
      ELSEIF(Option_OutputFormat == 'PLOT') THEN
!
            IF_Coord2: IF(coordinate_system == 'cyl' .OR. coordinate_system == 'CYL' .OR. coordinate_system == 'Cyl') THEN
!
               WRITE(Plot_Elem_Unit,6013) ADJUSTR(elem(n)%name),N,  &
     &                     elem(n)%coord(1), elem(n)%coord(3),    &
     &                     ElemState%pres(n,0),                  &
     &                     ElemState%temp(n,0),                  &
     &                     ElemProp(n,0)%satur(GasPhase),         &
     &                     ElemProp(n,0)%density(GasPhase),       &
     &                     ElemProp(n,0)%viscosity(GasPhase),     &
     &                     ElemProp(n,0)%MassFrac(gas_2,GasPhase),  &
     &                     ElemProp(n,0)%enthalpy(GasPhase),      &
     &                     ElemProp(n,0)%IntEnergy(GasPhase)
!
            ELSE
!
               WRITE(Plot_Elem_Unit,6014) ADJUSTR(elem(n)%name),N,                &
     &                     elem(n)%coord(1), elem(n)%coord(2), elem(n)%coord(3),  &
     &                     ElemState%pres(n,0),                  &
     &                     ElemState%temp(n,0),                  &
     &                     ElemProp(n,0)%satur(GasPhase),                         &
     &                     ElemProp(n,0)%density(GasPhase),                       &
     &                     ElemProp(n,0)%viscosity(GasPhase),                     &
     &                     ElemProp(n,0)%MassFrac(gas_2,GasPhase),                  &
     &                     ElemProp(n,0)%enthalpy(GasPhase),                      &
     &                     ElemProp(n,0)%IntEnergy(GasPhase)
!
            END IF IF_Coord2
!
      END IF IF_PrintOptions
!
      END DO DO_NumEle1
!
! -------
! ... Option for printing short FTSim outputs
! -------
!
      IF(Option_OutputFormat == 'BOTH') THEN
         WRITE(Plot_Elem_Unit,6106)
      ELSE IF(Option_OutputFormat == 'PLOT') THEN
         WRITE(Plot_Elem_Unit,6106)
         GO TO 1001
      END IF
!
!
!***********************************************************************
!*                                                                     *
!*                    PRINT HEADER INFORMATION - 2                     *
!*                                                                     *
!***********************************************************************
!
      WRITE(*, 6015)
      IF (MOD(Option_OutputAmount,10) < 2 .OR. NumConx == 0) GO TO 1001
!
      WRITE(*, 6016) ' ', title, Tot_NumTimeSteps, Tot_NumNRIterations, time
      WRITE(*, 6024) H0
      WRITE(*, 6026)
!
!
!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>                                                                     >
!>                      MAIN CONNECTION LOOP                           >
!>                        Print flow terms                             >
!>                                                                     >
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
      IF_ConxPrint: IF(Option_OutputFormat == 'STAN') THEN              ! Standard output
!
      DO_NumCon1: DO n=1,NumConx
!
         IF(conx(n)%n1 == 0 .OR. conx(n)%n2 == 0) CYCLE
!
         WRITE(*,6028)  ADJUSTR(conx(n)%name1), ADJUSTR(conx(n)%name2), N,  &
     &                  ConxFlow(n)%rate(GasPhase),                         &
     &                  ConxFlow(n)%CompInPhase(gas_2,GasPhase),              &
     &                  ConxFlow(n)%PoreVel(GasPhase)
!
      END DO DO_NumCon1
!
      ELSEIF(Option_OutputFormat == 'BOTH') THEN                        ! Standard output & plotting
!
      DO_NumCon2: DO n=1,NumConx
!
         IF(conx(n)%n1 == 0 .OR. conx(n)%n2 == 0) CYCLE
!
         WRITE(*,6028)  ADJUSTR(conx(n)%name1), ADJUSTR(conx(n)%name2), N,  &
     &                  ConxFlow(n)%rate(GasPhase),                         &
     &                  ConxFlow(n)%CompInPhase(gas_2,GasPhase),              &
     &                  ConxFlow(n)%PoreVel(GasPhase)
!
            WRITE(Plot_Conx_Unit,6014) ConxFlow(n)%heat,                                                             &
     &                                 ConxFlow(n)%rate(GasPhase),            ConxFlow(n)%rate(GasPhase),            &
     &                                 ConxFlow(n)%CompinPhase(gas_2,GasPhase), ConxFlow(n)%CompinPhase(gas_2,GasPhase), &
     &                                 ConxFlow(n)%PoreVel(GasPhase),         ConxFlow(n)%PoreVel(GasPhase)
!
      END DO DO_NumCon2
!
      ELSEIF(Option_OutputFormat == 'PLOT') THEN                        ! Standard plotting only
!
         DO_NumCon3: DO n=1,NumConx
!
            IF(conx(n)%n1 == 0 .OR. conx(n)%n2 == 0) CYCLE DO_NumCon3
!
            WRITE(Plot_Conx_Unit,6014) ConxFlow(n)%heat,                                                             &
     &                                 ConxFlow(n)%rate(GasPhase),            ConxFlow(n)%rate(GasPhase),            &
     &                                 ConxFlow(n)%CompinPhase(gas_2,GasPhase), ConxFlow(n)%CompinPhase(gas_2,GasPhase), &
     &                                 ConxFlow(n)%PoreVel(GasPhase),         ConxFlow(n)%PoreVel(GasPhase)
!
         END DO DO_NumCon3
!
! <<<
! <<< End of the CONNECTIONS LOOP
! <<<
!
      END IF IF_ConxPrint
!
! -------
! ... Option for printing short outputs
! -------
!
      IF(Option_OutputFormat == 'PLOT' .OR. Option_OutputFormat == 'BOTH') WRITE(Plot_Conx_Unit,6106)
!
!***********************************************************************
!*                                                                     *
!*                    PRINT HEADER INFORMATION - 2                     *
!*                                                                     *
!***********************************************************************
!
      WRITE(*, 6015)
!
      IF(MOD(Option_OutputAmount,10) < 3) GO TO 1001
!
      WRITE(*, 6016) ' ', title, Tot_NumTimeSteps, Tot_NumNRIterations, time
!
      WRITE(*, 6034) H0
!
!
!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>                                                                     >
!>                         ELEMENT LOOP - 3                            >
!>                Primary variables and their changes                  >
!>                                                                     >
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
      DO_NumEle3: DO n = 1,NumElemTot
!
! ...... Locations in the X array
!
         NLOC = Locp(n)
!
         IF_act4: IF(elem(n)%activity /= 'I') THEN    ! Active elements
!
            WRITE(*, FMT = 6040) ADJUSTR(elem(n)%name), N, ( X(NLOC+i), i = 1,NumComPlus1 ), &
     &                                                     ( DX(NLOC+i), i = 1,NumComPlus1 ), ElemMedia(n,current)%porosity
!
         ELSE
!
            WRITE(*, FMT = 6041) ADJUSTR(elem(n)%name), N, ( X(NLOC+i), i = 1,NumComPlus1 ), &
     &                                                     ( DX(NLOC+i), i = 1,NumComPlus1 ), ElemMedia(n,current)%porosity
!
         END IF IF_act4
!
! <<<
! <<< End of the ELEMENT LOOP
! <<<
!
      END DO DO_NumEle3
!
      WRITE(*, 6015)
!
!
!
 1001 IF(NumSS == 0) RETURN
!
! ... Initializations - only if active wells exist
!
      GC       = 0.0d0
      HGC      = 0.0d0
!
!***********************************************************************
!*                                                                     *
!*                    PRINT HEADER INFORMATION - 4                     *
!*                                                                     *
!***********************************************************************
!
      WRITE(*, 6016) ' ', title, Tot_NumTimeSteps, Tot_NumNRIterations, time
      WRITE(*, 6042) H0
      WRITE(*, 6044)
!
!
!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>                                                                     >
!>                         SINKS & SOURCES                             >
!>                                                                     >
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
      DO_NumSoSi: DO n=1,NumSS
!
! ...... Locations of the S/S element
!
         j = SS(n)%ElemNum
!
         IF (j == 0) CYCLE DO_NumSoSi
!
! ----------
! ...... Print-out options
! ----------
!
!
! ...... For injection
!
         IF_GPO: IF(SS(n)%rate > 0.0d0) THEN
!
            WRITE(*,6045) ADJUSTR(SS(n)%ElemName), SS(n)%name, N, SS(n)%rate, SS(n)%enth, SS(n)%type(1:4)
!
! ...... For production
!
         ELSE
!
            WRITE(*,6046)  ADJUSTR(SS(n)%ElemName), SS(n)%name, n, SS(n)%rate
!
         END IF IF_GPO
!
! <<<
! <<< End of the SINKS & SOURCES LOOP
! <<<
!
      END DO DO_NumSoSi
!
!
      WRITE(*, 6015)
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 6000 FORMAT(/,'Print_Standard_Output',T50,'[v1.0,   29 May      2008]',/,   &
     &         ':::::   Print results for elements, connections, and sinks & sources in the standard <FTSim> output file')
!
 6002 FORMAT(A1,/,' ',A120,//,   &
     &       10X,'Output data after (',I5,',',I6,') time steps ',53('.'),' The time is ',1pE12.5,' days')
 6004 FORMAT(/151('='),/,A1)
 6005 FORMAT(T3,  '  TOTAL    ',T19,'          ',T31,'             ',   &
     &       T46, ' Solution  ',T63,'      ',T77,'      ',T91,'  ',         &
     &       T103,'            ',T117,'. . . at',T129,'. . . at',T141,'        ')
 6006 FORMAT(T3,  ' TIME (s)',  T19,'Timestep #',T31,'NR-iterations',   &
     &       T46, 'Convergence',T63,'Max DP',T77,'Max DT',T91,'Max **',         &
     &       T103,'Max Residual',T117,'Elem. Num.',T129,'Eqn. Num.',T141,'TimeStep',/,   &
     &       1X,1pE12.5,T20,i7,T33,i7,T51,a1,T58,1pE13.6,   &
     &       T72,1pE13.6,T86,1pE13.6,T102,1pE13.6,T120,i4,T127,i6,T139,1pE13.6)
!
 6008 FORMAT(A1,/,T4,'ELEM',T11,'INDEX',T20,'Pressure',T33,'Temperature',T49,'S_Gas',T61,'Density',T74,   &
     &       'Viscosity',T88,'X_o_O',T100,'k_rel_o',T114,'*',T126,'*',T139,'*')
 6010 FORMAT(T20,'  (Pa)  ',T34,'(Deg. C)',T61,'(kg/m^3)',T76,'(Pa.s)',T114,'',/,151('_'),/)
!
 6011 FORMAT(A8,'*',I6,1X,1pe14.7,10(1X,1pE12.5))
 6012 FORMAT(A8,1X,I6,1X,1pe14.7,10(1X,1pE12.5))
 6013 FORMAT(A8,I4,12(1pE15.8,1x))
 6014 FORMAT(10(1pE15.8,1x))
 6015 FORMAT(/,151('@'))
 6016 FORMAT(A1,/,10X,A120,//,   &
     &       74X,'Timesteps = ',I6,' :|: NR iterations = ',I7,' :|: TIME = ',1pE13.6)
!
 6018 FORMAT(A1,/,145('_'),//,   &
     &       T4, 'ELEM',T11,'INDEX',T18,'C-AIRinGAS',T30,'C-AIRinAQU',   &
     &       T42,'Gas Density Aqu Density Ice Density  Gas Visco   Aqu Visco   Gas_Enth    Aqu_Enth    Ice_Enth')
!
 6020 FORMAT(T19,'(Kg/m^3)    (Kg/m^3)    (Kg/m^3)    (Kg/m^3)    (Kg/m^3)      ', &
     &           '(Pa*s)      (P*s)      (J/kg)      (J/kg)      (J/kg)',/,145('_'),/)
!
 6021 FORMAT(A8,'*',I6,1x,10(1pE12.5))
 6022 FORMAT(A8,1x,I6,1x,10(1pE12.5))
!
 6024 FORMAT(A1,/,124('_'),//,   &
     &       T4, 'ELEM1    ELEM2   INDEX',T31,'GasRate',T43,'gasInGasRate',T58,'GasPoreVel',T75,'*',   &
     &       T89,'*',T103,'*',T117,'*')
!
 6026 FORMAT(T33,'(W)',T46,'(Kg/s)',T60,'(Kg/s)',T74,'',T88,'',T102,'',T117,'',/,124('_'),/)
 6028 FORMAT(A8,1x,A8,1x,I7,1x,(7(1X,1pE13.6)))
!
 6031 FORMAT(A1,/,151('_'),//,   &
     &       T4,'ELEM   INDEX    Enthl-Gas    Enthl-Aqu     Enthl-Hydr',4x,'Enthl-Ice',4x,'Reactn_Rate',3x,'Reactn_Heat',/,   &
     &       T13,(4(8X,'(J/Kg)')),7x,'(Kg/s)',8x,'(J/kg)',/,151('_'),/)
!
 6034 FORMAT(A1,/,151('_'),//,   &
     &       T4,'ELEM   INDEX    Pres(Pa)       Temp(C)       deltaP       deltaT       Porosity                               ',   &
     &       /,151('_'),/)
!
!
 6040 FORMAT(A8,1x,I6,8(1x,1pE13.6))
 6041 FORMAT(A8,'*',I6,8(1x,1pE13.6))
!
 6042 FORMAT(A1,/,151('_'),//,   &
     &       ' ELEMENT SOURCE INDEX   Generation Rate    Enthalpy     GasFraction   AquFraction   ',   &
     &       'Air-MasRate  Air_VRateTot  Air_VRateGas  Air_VRateAqu   Inhib-MRate')
!
 6044 FORMAT(T26, '(Kg/s) or (W)      (J/Kg)',T87,'(kg/s)',   &
     &       T100, 'St m^3/s',T114,'St m^3/s',T128,'St m^3/s',T143,'(kg/s)',/,151('_'),/)
!
 6045 FORMAT(A8,2X,A5,2X,I4,4X,1pE12.5,4X,1pE12.5,'  (Source/Sink Type:',a4,')')
 6046 FORMAT(A8,2X,A5,2X,I4,4X,1pE12.5,4X,1pE12.5,7(2X,1pE12.5))
!
 6050 FORMAT(' >>>>> SOURCE  "',A5,'":     Rate =',1pE12.5,   &
     &       ' kg/s     Flowing Enthalpy =',1pE12.5,' J/kg     CH4 Mass Fraction =',1pE12.5/)
!
 6101 FORMAT(T1,'Variables = Name  N  r  z  P  T  S_gas  rho_gas  mu_gas  X_gas Enth.  IntEn.')
 6102 FORMAT(T1,'Variables = Name  N  x  y  z  P  T  S_gas  rho_gas  mu_gas  X_gas Enth.  IntEn.')
!
 6105 FORMAT(T1,'ZONE T(ime)= "',1pe15.8,'", F=POINT I(elem. total)=',i6)
!
!6105 FORMAT(T1,'TIME = ',1pe15.8,' : P[Pa] - T[C] - S_gas - S_aqu - S_ice - X_wG - X_wA - P_cap - k_rg - k_rw')
!
 6106 FORMAT(//)
!
 6115 FORMAT(T1,'TIME = ',1pe15.8,' : Heat_Flo [W] - G_Flo [kg/s] - A_Flo [kg/s] - H2OinG_Flo [kg/s] - ', &
     &                            'H2OinA_Flo [kg/s] - G_vel [m/s] - A_vel [m/s]')
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End of Print_Standard_Output
!
!
      RETURN
!
      END SUBROUTINE Print_Standard_Output
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE Mass_and_Energy_Balance(print_head_flag)
!
! ... Modules to be used
!
         USE Universal_Parameters
         USE General_External_File_Units
         USE Basic_Parameters
         USE General_Control_Parameters
!
         USE Grid_Geometry, ONLY: elem
         USE Element_Attributes
!
         USE Solution_Matrix_Arrays, ONLY: CO
!
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*          ROUTINE FOR COMPUTING VOLUME AND MASS BALANCES             *
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
      REAL(KIND = 8), DIMENSION(4)   :: vol_phase,mass_phase,mass_comp
      REAL(KIND = 8), DIMENSION(4,4) :: mass_CinPh
!
! -------
! ... Double precision variables
! -------
!
      REAL(KIND = 8) :: DAY,Sumx
!
! -------
! ... Integer variables
! -------
!
      INTEGER :: print_head_flag,n,i,k
!
! -------
! ... Logical variables
! -------
!
      LOGICAL :: First_call = .TRUE.

#ifdef USE_TIMER
         real :: start, finish
#endif
!
! -------
! ... Saving variables
! -------
!
      SAVE First_call
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of Mass_and_Energy_Balance
!
!
      IF(First_call) THEN
         First_call = .FALSE.
         WRITE(VERS_Unit,6000)
      END IF
!
!***********************************************************************
!*                                                                     *
!*                        BASIC INITIALIZATIONS                        *
!*                                                                     *
!***********************************************************************
!
      DAY = time/8.64E4
!
! ----------
! ... Initialize arrays - CAREFUL: Whole array operation
! ----------
!
#ifdef USE_TIMER
         call cpu_time(start)
#endif

#ifdef USE_OMP
!$OMP WORKSHARE
#endif
      vol_phase  = 0.0d0
      mass_phase = 0.0d0
      mass_comp  = 0.0d0
      mass_CinPh = 0.0d0
#ifdef USE_OMP
!$OMP END WORKSHARE
#endif

#ifdef USE_TIMER
         call cpu_time(finish)

         write (*,*) __FILE__, ":", __LINE__, " time: ", finish-start
#endif
!
! ... Printing out headings
!
      IF(print_head_flag == 0) THEN
         WRITE(*,6001)
         WRITE(*,6002) Tot_NumTimeSteps, Tot_NumNRIterations, time, DAY
      END IF
!
!
!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>                                                                     >
!>                PHASE LOOP - ! Whole Array Operations                >
!>                                                                     >
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
      DO_NumPhase: DO i=1,NumPhases
!
! ...... Determine the product (Vol*Phi*Saturation)
!
         FORALL (n=1:NumElem) CO(n) = ElemProp(n,0)%satur(i)
!
!array operation
         CO(NumElem+1 : 2*NumElem) = ElemMedia(1:NumElem,current)%porosity * elem(1:NumElem)%vol * CO(1:NumElem)
!CO is a 1-D array w/: sat in 1st block, phi*V*sat in 2nd block
!
! ...... Determine the product (Vol*Phi*Saturation*Dens)
!
!replace sat in 1st block w/ density
         FORALL (n=1:NumElem) CO(n) = ElemProp(n,0)%density(i)
!
!add a 3rd array block w/ Vol*Phi*Saturation*Dens product
         CO(2*NumElem+1 : 3*NumElem) = CO(NumElem+1:2*NumElem)*CO(1:NumElem)
!
! ...... Compute phase volumes and masses
!
         vol_phase(i)  = SUM(CO(  NumElem+1 : 2*NumElem))
         mass_phase(i) = SUM(CO(2*NumElem+1 : 3*NumElem))
!
         IF(print_head_flag /= 0) CYCLE DO_NumPhase
!
!
! >>>>>>>>>>>>>>>>>>
! ......
! ...... COMPONENT LOOP - Summation of component masses in each phase
! ......
! >>>>>>>>>>>>>>>>>>
!
!
         DO_NumComp: DO k=1,NumCom
!
! ......... Determine the product (Vol*Phi*Saturation*Dens*Mass fraction)
!
            FORALL (n=1:NumElem) CO(n) = ElemProp(n,0)%MassFrac(k,i)
            Sumx = SUM(CO(2*NumElem+1 : 3*NumElem)*CO(1:NumElem))
!
! ......... Compute mass of component k
!
            mass_comp(k)    = mass_comp(k)    + Sumx    ! Total mass in all phases
            mass_CinPh(k,i) = mass_CinPh(k,i) + Sumx    ! Mass in phase i
!
! <<<<<<
! <<<<<< End of the COMPONENT LOOP
! <<<<<<
!
         END DO DO_NumComp
!
! <<<
! <<< End of the PHASE LOOP
! <<<
!
      END DO DO_NumPhase
!
!
!***********************************************************************
!*                                                                     *
!*                        P R I N T - O U T S                          *
!*                                                                     *
!***********************************************************************
!
      IF(print_head_flag /= 0) THEN
         write(Elem_TimeSeries_Unit,6010)  time,(mass_phase(i),i=1,NumPhases)
         RETURN
      END IF
!
! -------
! ... Phase-based mass and volume balances
! -------
!
      WRITE(*, FMT = 6003) (vol_phase(i), i=1,NumPhases), (mass_phase(i),i=1,NumPhases)
!
! -------
! ... Component and component-in-phase mass balances
! -------
!
      WRITE(*, FMT = 6004) ((mass_CinPh(k,i),k=1,NumCom), i=1,NumPhases),(mass_comp(k),k=1,NumCom)
!
      WRITE(*, FMT = 6005)
!
! -------
! ... Reset the CO array - CAREFUL! Array operation
! -------
!
      CO(1:3*NumElem) = 0.0d0

!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 6000 FORMAT(/,'Mass_and_Energy_Balance',T50,'[v1.0,   29 May      2008]', /, &
               ':::::   Perform summary balances for volume, mass, and energy')
!
 6001 FORMAT(//,130('*'),/,'*',T130,'*',/,'*',   &
     &       T43,'V O L U M E   A N D   M A S S   B A L A N C E S',   &
     &       T130,'*',/,'*',T130,'*',/,130('*'),/)
 6002 FORMAT(' =========  [Timestep, NR-iterations] = [',I6,',',I7,']  =========',13X,   &
     &       'The time is ',1pE12.5,' sec, or ',1pe12.5,' days'/)
!
 6003 FORMAT(/,36X,'PHASES PRESENT',/,   &
     &       ' ',67('='),/,   &
     &       '  PHASES       |                       Gas',/,   &
     &       ' ',67('='),/,15x,'|',/,   &
     &       '  VOLUME (m^3) |',18X,1pE15.8,/,   &
     &       '  MASS   (Kg)  |',18X,1pE15.8,/,   &
     &       15x,'|',/,' ',67('='))
!
 6004 FORMAT(//,95X,   &
     &       'COMPONENT MASS IN PLACE (Kg)',/,   &
     &       73x,' ',53('='),/,   &
     &       73x,'   COMPONENTS    |                gas ',/,   &
     &       73x,' ',53('='),/,   &
     &       76x,'PHASES',8X,'|',/,74x,16('-'),'|',/,   &
     &       76X,'Gas phase |',10X,1pE15.8,/,   &
     &       73X,' ',53('-'),/,   &
     &       76X,'TOTAL         |',10X,1pE15.8,/,   &
     &       90x,'|',/,73x,' ',53('='))
!
 6005 FORMAT(/,1X,129('*'),/,1X,129('*'),//)
!
 6010 FORMAT(1pe12.5,1x,11(1pe13.6,1x))
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of Mass_and_Energy_Balance
!
!
      RETURN
!
      END SUBROUTINE Mass_and_Energy_Balance
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
