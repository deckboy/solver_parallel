!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>                                                                     >
!>                                                                     >
!>            GAS_Properties.f95: Module including all                 >
!>               gas-related properties and processes                  >
!>                                                                     >
!>                                                                     >
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  Segment that includes (a) all the gas-related constants (parameter),
!  and (b) procedures describing the water behavior and thermophysical
!  properties/processes in its entire thermodynamic phase diagram.
!
!
!
      MODULE Gas_Parameters
!
         SAVE
!
! ----------
! ...... Derived types
! ----------
!
         TYPE Gas_Basic_Parameters
!
            CHARACTER(LEN = 20) :: name                                         ! ...... Gas name
!
! ......... Reference gas properties
!
            REAL(KIND = 8) :: MolWeight                    ! ...... Molecular weight
!
            REAL(KIND = 8) :: RefThermalCond               ! ...... Thermal conductivity
            REAL(KIND = 8) :: StandardDensity              ! ...... Reference density
            REAL(KIND = 8) :: RefViscosity                 ! ...... Gas reference viscosity
            REAL(KIND = 8) :: RefDensity
            REAL(KIND = 8) :: compressibility
            REAL(KIND = 8) :: TempAtRefDensity
            REAL(KIND = 8) :: PresAtRefDensity
!
            REAL(KIND = 8) :: reference_departure_IntEnergy, reference_departure_enthalpy
!
            REAL(KIND = 8) :: PCrit, TCrit, VCrit, ZCrit, omega
!
            INTEGER :: VisCoeffNumber, ThermCondCoeffNumber
!
            REAL(KIND = 8), DIMENSION(0:4) :: Cp_coefficients
!
            REAL(KIND = 8), ALLOCATABLE, DIMENSION(:)   :: ThermCond_coefficients, viscosity_coefficients
!
         END TYPE Gas_Basic_Parameters
!
! ----------
! ...... Derived-type variables
! ----------
!
         TYPE(Gas_Basic_Parameters) :: gas         !The gas object add ,DIMENSION(x) to have a library of x # of gases (i.e. methane, ethane, etc.)
!
! -------
! ... LOGICAL variables
! -------
!
! ... If any of the coefficients in the above-defined coefficient array are nonzero, then LOGICAL = TRUE
         LOGICAL :: constant_gas_Cp, constant_gas_Cv
         LOGICAL :: constant_gas_thermal_cond
         LOGICAL :: constant_gas_viscosity
!
      END MODULE Gas_Parameters
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      MODULE Gas_Thermophysical_Properties
!
         IMPLICIT NONE
!
         PRIVATE
!
         PUBLIC :: Compute_Gas_Density, Compute_Ideal_Gas_Enthalpy, Compute_Departure_Enthalpy,  &
        &          Compute_Gas_Viscosity, Compute_Gas_Therm_Conductivity
!
         CONTAINS
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!    &        "CH4   ",                                 &                ! Name, AtomDV                    xxxxxPROPERTIES TO BE ADDED INTO GAS NAMELIST
!    &         1.9056d+02,  4.600155d+06,  9.9000d-02,  2.880d-01,   &   ! TCrit,PCrit,VCrit,ZCrit
!    &         1.1000d-02,  1.604300d-02,                            &   ! Omega,MolWt
!    &         4.568d+00,   -8.9750d-03,  3.631000d-05, -3.4070d-08,  1.091d-11,   &   ! A0,A1,A2,A3,A4
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      LOGICAL FUNCTION Compute_Gas_Density( pressure, temperature, compressibility, gas_density, &
     &                                      TR, RK, alpha, AxAL, AUP, BUP )
!
! ...... Modules to be used
!
         USE Universal_Parameters
!
         USE Gas_Parameters
         USE Utility_Functions
         USE General_External_File_Units
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*                  Computation of compressibility Z                   *
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
      REAL(KIND = 8), INTENT(IN)  :: pressure, temperature
!
      REAL(KIND = 8), INTENT(OUT) :: compressibility, gas_density
      REAL(KIND = 8), INTENT(OUT) :: TR, RK, alpha, AxAL, AUP, BUP
!
      REAL(KIND = 8) :: BB, AA0,AA1,AA2
      REAL(KIND = 8) :: MaxSolution, MinSolution
!
! -------
! ... Integer variables
! -------
!
      INTEGER :: ixx
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
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of <Compute_Gas_Density>
!
!
      IF(First_call) THEN
         First_call = .FALSE.
         WRITE(UNIT = VERS_Unit, FMT = 6000)
      END IF
!
! -------
! ... Initialization
! -------
!
      Compute_Gas_Density = .TRUE.
!
!***********************************************************************
!*                                                                     *
!*                     For the Peng-Robinson EOS                       *
!*                                                                     *
!***********************************************************************
!
      TR = temperature/gas%TCrit
!
      RK = 3.7464d-1 + 1.54226d0*gas%Omega - 2.6992d-1 * gas%Omega * gas%Omega
!
      alpha  = (1.0d0+RK*(1.0d0-SQRT(TR)))*(1.0d0 + RK*(1.0d0 - SQRT(TR)))
!
      AxAL = 4.5724d-1 * alpha * R_gas * R_gas * gas%TCrit * gas%TCrit/gas%PCrit
!
      BB   = 7.78d-2 * R_gas * gas%TCrit/gas%PCrit
!
      AUP = AxAL * pressure/(R_gas * R_gas * temperature * temperature)
      BUP = BB * pressure/(R_gas * temperature)
!
! -------
! ... The parameters of the cubic equation
! -------
!
      AA2 = -(1.0d0 - BUP)
      AA1 =  AUP - 3.0d0*BUP*BUP - 2.0d0*BUP
      AA0 = -BUP*( AUP - BUP*(1.0d0 + BUP) )
!
!***********************************************************************
!*                                                                     *
!*                 SOLVE THE CUBIC EQUATION OF STATE                   *
!*                                                                     *
!***********************************************************************
!
      compressibility = 0.0d0  ! Initialization
!
      CALL Cubic_Equation_Roots(AA2, AA1, AA0, MaxSolution, MinSolution, ixx)
!
! -------
! ... Determine the compressibility factor
! -------
!
      compressibility = MaxSolution
!
! -------
! ... Compute the densities
! -------
!
      gas_density = 1.0d3*pressure*gas%MolWeight/(compressibility*R_gas*temperature)
!
      IF(compressibility <= 0.0d0 ) THEN
         Compute_Gas_Density = .FALSE.
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
 6000 FORMAT(/,'Compute_Gas_Density',T50,'[v1.0,  20 June      2008]',/,      &
     &         ':::::   Computation of the compressibility factor Z and density of a real gas')
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of <Compute_Gas_Density>
!
!
      RETURN
!
      END FUNCTION Compute_Gas_Density
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      LOGICAL FUNCTION Compute_Ideal_Gas_Enthalpy(T0, T1, enthalpy, internal_energy)
!
! ...... Modules to be used
!
         USE Universal_Parameters
!
         USE Gas_Parameters
!
         USE Utility_Functions
         USE General_External_File_Units
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*           Computation of ideal enthalpy and entropy changes         *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
!
      IMPLICIT NONE
!
! -------
! ... Double precision variables
! -------
!
      REAL(KIND = 8), INTENT(IN)  :: T0, T1
      REAL(KIND = 8), INTENT(OUT) :: enthalpy, internal_energy
!
      REAL(KIND = 8) :: Hcom
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
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of <Compute_Ideal_Gas_Enthalpy>
!
!
      IF(First_call) THEN
         First_call = .FALSE.
         WRITE(UNIT = VERS_Unit, FMT = 6000)
      END IF
!
! -------
! ... Initialization
! -------
!
      Compute_Ideal_Gas_Enthalpy = .TRUE.
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*    COMPUTE ENTHALPY AND ENTROPY CHANGES BASED ON IDEAL BEHAVIOR     *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
      Hcom = Integral_poly(T1,gas%Cp_coefficients(0:4)) - Integral_poly(T0,gas%Cp_coefficients(0:4))
!
      enthalpy = Hcom*R_gas
!
      internal_energy = enthalpy - (T1-T0)*R_gas
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 6000 FORMAT(/,'Compute_Ideal_Gas_Enthalpy',T50,'[v1.0,  11 April     2004]',/,  &
     &         ':::::   Computation of the ideal enthalpy and internal energy of ideal gases')
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of <Compute_Ideal_Gas_Enthalpy>
!
!
      RETURN
!
      END FUNCTION Compute_Ideal_Gas_Enthalpy
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      LOGICAL FUNCTION Compute_Departure_Enthalpy( temperature, compressibility, departure_enthalpy, departure_internal_energy, &
     &                                             TR, RK, alpha, AxAL, AUP, BUP )
!
! ...... Modules to be used
!
         USE Universal_Parameters
!
         USE Gas_Parameters
         USE General_External_File_Units
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*            Computation of departure enthalpy and entrop             *
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
      REAL(KIND = 8), INTENT(IN)  :: temperature, compressibility
      REAL(KIND = 8), INTENT(IN)  :: TR, RK, alpha, AxAL, AUP, BUP
!
      REAL(KIND = 8), INTENT(OUT) :: departure_enthalpy, departure_internal_energy
!
      REAL(KIND = 8) :: Di0, Z0, D0, TR0, RK0, alpha0, AxAL0, AUP0, BUP0, Di
!
! -------
! ... Logical variables
! -------
!
      LOGICAL :: status, First_call = .TRUE.
!
! -------
! ... Saving variables
! -------
!
      SAVE First_call
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of <Compute_Departure_Enthalpy>
!
!
      IF(First_call) THEN
!
         First_call = .FALSE.
         WRITE(UNIT = VERS_Unit, FMT = 6000)
!
         status = Compute_Gas_Density( 1.0d4, 273.15d0, Z0, D0, TR0, RK0, alpha0, AxAL0, AUP0, BUP0 )
!
         Di0 = RK0 * AxAL0 * SQRT(273.15d0/gas%Tcrit/alpha0)
!
         gas%reference_departure_IntEnergy = R_gas*273.15d0*(AUP0/BUP0/2.828d0)*(1.0d0 + Di0/AxAL0)   &
     &                                        *LOG((Z0 + 2.414d0*BUP0)/(Z0 - 0.414d0*BUP0))
!
         gas%reference_departure_enthalpy  =  gas%reference_departure_IntEnergy + R_gas * 273.15d0 * (1.0d0 - Z0)
!
      END IF
!
! -------
! ... Initialization
! -------
!
      Compute_Departure_Enthalpy = .TRUE.
!
!***********************************************************************
!*                                                                     *
!*                     For the Peng-Robinson EOS                       *
!*                                                                     *
!***********************************************************************
!
      Di = RK * AxAL * SQRT(TR/alpha)
!
      departure_internal_energy = R_gas*temperature*(AUP/BUP/2.828d0)*(1.0d0 + Di/AxAL)   &
     &                                 *LOG((compressibility + 2.414d0*BUP)/(compressibility - 0.414d0*BUP))
!
      departure_enthalpy =  departure_internal_energy + R_gas * temperature * (1.0d0 - compressibility)
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 6000 FORMAT(/,'Compute_Departure_Enthalpy',T50,'[v1.0,  11 April     2004]',/,  &
     &         ':::::   Computation of the enthalpy, internal energy and entropy departures of real gases')
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of <Compute_Departure_Enthalpy>
!
!
      RETURN
!
      END FUNCTION Compute_Departure_Enthalpy
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
         LOGICAL FUNCTION Compute_Gas_Viscosity(temperature, gas_density, gas_viscosity)
!
! ......... Modules to be used
!
            USE General_External_File_Units
!
            USE Gas_Parameters
!
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*                Compute a single gas viscosity                       *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
!
         IMPLICIT NONE
!
! ----------
! ...... Double precision variables
! ----------
!
         REAL(KIND = 8), INTENT(IN)  :: temperature, gas_density      !Temp in C
!
         REAL(KIND = 8), INTENT(OUT) :: gas_viscosity                 ! Density in Pa.s
!
         REAL(KIND = 8)              :: T_k
!
! ----------
! ...... LOGICAL variables
! ----------
!
         LOGICAL :: First_call = .TRUE.
!
! ----------
! ...... Saving variables
! ----------
!
         SAVE First_call
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of <Compute_Gas_Viscosity>
!
!
         IF(First_call) THEN
            First_call = .FALSE.
            WRITE(UNIT = VERS_Unit, FMT = 6000)
         END IF
!
! -------------

! ...... Initialization
! -------------
!
         Compute_Gas_Viscosity = .TRUE.
!
! -------------
! ...... Gas viscosity computation
! -------------
!
      T_k = temperature + 273.15d0
!
      IF_TDepend: IF(constant_gas_viscosity) THEN
         gas_viscosity = gas%RefViscosity
      ELSE
!
! ... check equation applicability...eqn. here
         gas_viscosity = 0.001*(gas%viscosity_coefficients(1)          &
        &              + gas%viscosity_coefficients(2)*T_k             &
        &              + gas%viscosity_coefficients(3)*T_k**2          &
        &              + gas%viscosity_coefficients(4)*T_k**3          &
        &              + gas%viscosity_coefficients(5)*gas_density     &
        &              + gas%viscosity_coefficients(6)*gas_density**2  &
        &              + gas%viscosity_coefficients(7)*gas_density**3  &
        &              + gas%viscosity_coefficients(8)*gas_density**4)
!
      END IF IF_TDepend
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 6000 FORMAT(/,'Compute_Gas_Viscosity',T50,'[v1.0,  12 June 2008]',6X,/,   &
     &         ':::::   Calculate the gas viscosity (scalar version)')
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End of <Compute_Gas_Viscosity>
!
!
         RETURN
!
         END FUNCTION Compute_Gas_Viscosity
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
         LOGICAL FUNCTION Compute_Gas_Therm_Conductivity(temperature, gas_therm_conductivity)
!
! ......... Modules to be used
!
            USE General_External_File_Units
!
            USE Gas_Parameters
!
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*            Compute a single gas thermal conductivity                *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
!
         IMPLICIT NONE
!
! ----------
! ...... Double precision variables
! ----------
!
         REAL(KIND = 8), INTENT(IN)  :: temperature               ! Temp in C
         REAL(KIND = 8), INTENT(OUT) :: gas_therm_conductivity    ! Thermal Conductivity in XX(units)
!
         REAL(KIND = 8)              :: exponent, T_k
!
! ----------
! ...... LOGICAL variables
! ----------
!
         LOGICAL :: First_call = .TRUE.
!
! ----------
! ...... Saving variables
! ----------
!
         SAVE First_call
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of <Compute_Gas_Thermal_Conductivity>
!
!
         IF(First_call) THEN
            First_call = .FALSE.
            WRITE(UNIT = VERS_Unit, FMT = 6000)
         END IF
!
! -------------
! ...... Initialization
! -------------
!
         Compute_Gas_Therm_Conductivity = .TRUE.
!
! -------------
! ...... Gas thermal conductivity computation
! -------------
!
      IF(constant_gas_thermal_cond) THEN
         gas_therm_conductivity = gas%RefThermalCond
      ELSE
!
! ...Check equation applicability
!
         gas_therm_conductivity =   gas%ThermCond_coefficients(1) + gas%ThermCond_coefficients(2) * temperature  &
        &                         + gas%ThermCond_coefficients(3) * temperature * temperature
!
      END IF
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 6000 FORMAT(/,'Compute_Gas_Thermal_Conductivity',T50,'[v1.0,   01 June      2008]',6X,/,   &
     &         ':::::   Calculate the gas thermal conductivity (scalar version)')
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End of <Compute_gas_Thermal_Conductivity>
!
!
         RETURN
!
         END FUNCTION Compute_Gas_Therm_Conductivity
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
      END MODULE Gas_Thermophysical_Properties
!
!
!
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

