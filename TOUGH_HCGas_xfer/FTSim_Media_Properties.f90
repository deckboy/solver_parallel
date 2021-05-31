!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>                                                                     >
!>                                                                     >
!>   FTSim_Media_Properties.f95: Code unit including the routines      >
!>       describing the porous media properties and processes          >
!>                                                                     >
!>                                                                     >
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  Segment that describes the hydraulic and thermal behavior of the
!  geologic medium (porous or fractured), i.e., capillary pressure and
!  relative permeability under multiphase conditions, interface
!  permeability and mobility, and interface thermal conductivity
!
!
      MODULE Geological_Media_Properties
!
! ...... Modules to be used
!
         IMPLICIT NONE
!
         PRIVATE
!
         PUBLIC :: Hydraulic_Property_Changes, Composite_Thermal_Conductivity
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
!
      SUBROUTINE Hydraulic_Property_Changes
!
! ...... Modules to be used
!
         USE General_External_File_Units
!
         USE Basic_Parameters, N1 => ElemArraySize
         USE EOS_Default_Parameters
!
         USE General_External_File_Units
         USE General_Control_Parameters
!
         USE Grid_Geometry
!
         USE Element_Attributes
!
         USE Geologic_Media_Parameters
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*  Computes the wettability properties of the gas-bearing medium      *
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
      INTEGER :: n, k
!
      INTEGER(KIND = 2) :: mn
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
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of Hydraulic_Property_Changes
!
!
      IF(First_call) THEN
!
         First_call = .FALSE.
         WRITE(VERS_Unit,6000)
!
      END IF
!
!
!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>                                                                     >
!>             Compute the incremented state of porosity               >
!>                                                                     >
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
      IF_VarPorosity: IF(variable_porosity) THEN
!
#ifdef USE_TIMER
            call CPU_Timing_Routine(start)
#endif

#ifdef USE_OMP
!$OMP PARALLEL PRIVATE(k, n, mn)
!$OMP DO schedule(auto)
#endif
         DO_NumEquA: DO k = 0,NumEqu
!
!***********************************************************************
!************ Parallelize the computation of the state of Porosity *****
!***********************************************************************
            DO n = 1,N1
!
               mn = elem%MatNum(n)
!
! ............ Standard treatment: Porosity as an exponential function of P, T, compressibility and expansivity
!
               ElemMedia(n,k)%porosity =  ElemMedia(n,original)%porosity                                                 &
     &                                  * EXP(  media(mn)%Compr * ( ElemState%pres(n,k) - ElemState%pres(n,initial) )  &
     &                                        + media(mn)%Expan * ( ElemState%temp(n,k) - ElemState%temp(n,initial) )  )
!
            END DO
!
         END DO DO_NumEquA
#ifdef USE_OMP
!$OMP END DO
!$OMP END PARALLEL
#endif

#ifdef USE_TIMER
            call CPU_Timing_Routine(finish)

            write (*,*) __FILE__, ":", __LINE__, " time: ", finish-start
#endif
!
      END IF IF_VarPorosity
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 6000 FORMAT(/,'Hydraulic_Property_Changes',T50,'[v1.0,   01 June      2008]',/,  &
     &         ':::::     Transient porosity and permeability as a function of P, T, stress and strain')
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of Hydraulic_Property_Changes
!
!
      RETURN
!
      END SUBROUTINE Hydraulic_Property_Changes
!
!
!
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!
    REAL(KIND = 8) FUNCTION Composite_Thermal_Conductivity(n1, n2, m1, m2, nm1, nm2, wt1in, wt2in, gas_thermal_cond1, &
   &                                                       gas_thermal_cond2, liquid_phases)
!
! ...... Modules to be used
!
         USE General_External_File_Units
!
         USE Basic_Parameters, ONLY: NumPhases, noise
         USE EOS_Default_Parameters
!
!         USE General_External_File_Units
         USE General_Control_Parameters
!
         USE Grid_Geometry
!
         USE Element_Attributes
!
         USE Geologic_Media_Parameters
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*  Computes the wettability properties of the gas-bearing medium      *
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
      REAL(KIND = 8), INTENT(IN), OPTIONAL, DIMENSION(NumPhases) :: gas_thermal_cond1, gas_thermal_cond2
!
! -------
! ... Double precision variables
! -------
!
      REAL(KIND = 4), INTENT(IN) :: wt1in, wt2in
      REAL(KIND = 8) :: liquid_saturation, Thermal_Conductivity1, Thermal_Conductivity2
!
! -------
! ... Integer variables
! -------
!
      INTEGER, INTENT(IN) :: n1, n2, m1, m2
      INTEGER(KIND = 2), INTENT(IN) :: nm1, nm2
      INTEGER :: number_liquid_phases
!
! -------
! ... Integer arrays
! -------
!
       INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:) :: liquid_phases
!
! -------
! ... Logical variables
! -------
!
      LOGICAL :: First_call = .TRUE.
      LOGICAL :: dry_system = .FALSE.
!
! -------
! ... Saving variables
! -------
!
      SAVE First_call, number_liquid_phases, dry_system
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of Composite_thermal_conductivity
!
!
      IF(First_call) THEN
!
         IF(Option_ThermalConductivity < 2)  THEN
            number_liquid_phases = SIZE(liquid_phases)
!
            IF(number_liquid_phases == 1 .AND. liquid_phases(1) == 0) THEN
               dry_system = .TRUE.
            ENDIF
         END IF
!
         First_call = .FALSE.
         WRITE(UNIT = VERS_Unit, FMT = 6000)
!
      END IF
!
!
!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>                                                                     >
!>             Compute the incremented state of porosity               >
!>                                                                     >
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
      Case_Option: SELECT CASE (Option_ThermalConductivity)
!
         CASE (2)   ! Linear model
!
            Thermal_Conductivity1 =   media(nm1)%KThrD &
           &                        + ElemMedia(n1,m1)%porosity * SUM(ElemProp(n1,m1)%satur(1:NumPhases) * gas_thermal_cond1(1:NumPhases))
            Thermal_Conductivity2 =   media(nm2)%KThrD &
           &                        + ElemMedia(n2,m2)%porosity * SUM(ElemProp(n2,m2)%satur(1:NumPhases) * gas_thermal_cond2(1:NumPhases))
!
         CASE (1)   ! Sommerton #1 Model
!
            IF(dry_system) THEN
               Thermal_Conductivity1 = media(nm1)%KThrD
               Thermal_Conductivity2 = media(nm2)%KThrD
            ELSE
!
               liquid_saturation     = SUM(ElemProp(n1,m1)%satur(liquid_phases(1:number_liquid_phases)))
               Thermal_Conductivity1 = media(nm1)%KThrD + SQRT(liquid_saturation) * (media(nm1)%KThrW - media(nm1)%KThrD)
!
               liquid_saturation     = SUM(ElemProp(n2,m2)%satur(liquid_phases(1:number_liquid_phases)))
               Thermal_Conductivity2 = media(nm2)%KThrD + SQRT(liquid_saturation) * (media(nm2)%KThrW - media(nm2)%KThrD)
!
            END IF
!
         CASE (0)   ! Sommerton #2 Model
!
            IF(dry_system) THEN
               Thermal_Conductivity1 = media(nm1)%KThrD
               Thermal_Conductivity2 = media(nm2)%KThrD
            ELSE
!
               liquid_saturation     = SUM(ElemProp(n1,m1)%satur(liquid_phases(1:number_liquid_phases)))
               Thermal_Conductivity1 = media(nm1)%KThrD + liquid_saturation * (media(nm1)%KThrW - media(nm1)%KThrD)
!
               liquid_saturation     = SUM(ElemProp(n2,m2)%satur(liquid_phases(1:number_liquid_phases)))
               Thermal_Conductivity2 = media(nm2)%KThrD + liquid_saturation * (media(nm2)%KThrW - media(nm2)%KThrD)
!
            END IF
!
         CASE DEFAULT   !Unavailable Model - print an error message
            WRITE(*,6100) Option_ThermalConductivity
            STOP
!
      END SELECT Case_Option
!
! ......Compoupte the interblock thermal conductivity
!
      Composite_Thermal_Conductivity =  Thermal_Conductivity1 * Thermal_Conductivity2 /  &
     &                                 (Thermal_Conductivity1 * wt1in + Thermal_Conductivity2 * wt2in + noise)
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 6000 FORMAT(/,'Composite_Thermal_Conductivity',T50,'[v1.0,   01 June      2008]',/,                        &
     &         ':::::     Transient porosity and permeability as a function of P, T, stress and strain')
!
 6100 FORMAT(//,22('ERROR-'),//,T33,                                                                                          &
     &             '       S I M U L A T I O N   A B O R T E D',/,                                                            &
     &         T22,'There is a problem reading the case <Option_ThermalConductivity> in the <Composite_Thermal_Conductivity>  &
     &         data block', /,T32,'               CORRECT AND TRY AGAIN', //,22('ERROR-'))
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of Composite_thermal_conductivity
!
!
      RETURN
!
      END FUNCTION Composite_thermal_conductivity
!
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!
      END MODULE Geological_Media_Properties
!
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
