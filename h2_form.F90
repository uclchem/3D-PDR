!=======================================================================
!
!  Calculate the rate of molecular hydrogen (H2) formation on grains
!  using the treatment of Cazaux & Tielens (2002, ApJ, 575, L29) and
!  Cazaux & Tielens (2004, ApJ, 604, 222).
!
!-----------------------------------------------------------------------
FUNCTION H2_FORMATION_RATE(GAS_TEMPERATURE,GRAIN_TEMPERATURE) RESULT(RATE)

   USE DEFINITIONS
   USE HEALPIX_TYPES
   USE GLOBAL_MODULE, ONLY: metallicity, g2d
   IMPLICIT NONE

   REAL(KIND=DP) :: RATE
   REAL(KIND=DP), INTENT(IN) :: GAS_TEMPERATURE,GRAIN_TEMPERATURE

   REAL(KIND=DP) :: THERMAL_VELOCITY,STICKING_COEFFICIENT,TOTAL_CROSS_SECTION
   REAL(KIND=DP) :: FLUX,FACTOR1,FACTOR2,EPSILON
   REAL(KIND=DP) :: SILICATE_FORMATION_EFFICIENCY,GRAPHITE_FORMATION_EFFICIENCY
   REAL(KIND=DP) :: SILICATE_CROSS_SECTION,SILICATE_MU,SILICATE_E_S,SILICATE_E_H2
   REAL(KIND=DP) :: SILICATE_E_HP,SILICATE_E_HC,SILICATE_NU_H2,SILICATE_NU_HC
   REAL(KIND=DP) :: GRAPHITE_CROSS_SECTION,GRAPHITE_MU,GRAPHITE_E_S,GRAPHITE_E_H2
   REAL(KIND=DP) :: GRAPHITE_E_HP,GRAPHITE_E_HC,GRAPHITE_NU_H2,GRAPHITE_NU_HC

!  Mean thermal velocity of hydrogen atoms (cm s^-1)
   THERMAL_VELOCITY=1.45D5*SQRT(GAS_TEMPERATURE/1.0D2)

!  Calculate the thermally averaged sticking coefficient of hydrogen atoms on grains,
!  as given by Hollenbach & McKee (1979, ApJS, 41, 555, eqn 3.7)
   STICKING_COEFFICIENT=1.0D0/(1.0D0+0.04D0*SQRT(GAS_TEMPERATURE+GRAIN_TEMPERATURE) &
                    & + 0.2D0*(GAS_TEMPERATURE/1.0D2)+0.08D0*(GAS_TEMPERATURE/1.0D2)**2)

   FLUX=1.0D-10 ! Flux of H atoms in monolayers per second (mLy s^-1)

   TOTAL_CROSS_SECTION=6.273D-22 ! Total mixed grain cross section per H nucleus (cm^-2/nucleus)
   SILICATE_CROSS_SECTION=8.473D-22 ! Silicate grain cross section per H nucleus (cm^-2/nucleus)
   GRAPHITE_CROSS_SECTION=7.908D-22 ! Graphite grain cross section per H nucleus (cm^-2/nucleus)

!  Silicate grain properties
   SILICATE_MU=0.005D0   ! Fraction of newly formed H2 that stays on the grain surface
   SILICATE_E_S=110.0D0  ! Energy of the saddle point between a physisorbed and a chemisorbed site (K)
   SILICATE_E_H2=320.0D0 ! Desorption energy of H2 molecules (K)
   SILICATE_E_HP=450.0D0 ! Desorption energy of physisorbed H atoms (K)
   SILICATE_E_HC=3.0D4   ! Desorption energy of chemisorbed H atoms (K)
   SILICATE_NU_H2=3.0D12 ! Vibrational frequency of H2 molecules in surface sites (s^-1)
   SILICATE_NU_HC=1.3D13 ! Vibrational frequency of H atoms in their surface sites (s^-1)

   FACTOR1=SILICATE_MU*FLUX/(2*SILICATE_NU_H2*EXP(-SILICATE_E_H2/GRAIN_TEMPERATURE))

   FACTOR2=1.0D0*(1.0D0+SQRT((SILICATE_E_HC-SILICATE_E_S)/(SILICATE_E_HP-SILICATE_E_S)))**2 &
        & /4.0D0*EXP(-SILICATE_E_S/GRAIN_TEMPERATURE)

   EPSILON=1.0D0/(1.0D0+SILICATE_NU_HC/(2*FLUX)*EXP(-1.5*SILICATE_E_HC/GRAIN_TEMPERATURE) &
              & *(1.0D0+SQRT((SILICATE_E_HC-SILICATE_E_S)/(SILICATE_E_HP-SILICATE_E_S)))**2)

   SILICATE_FORMATION_EFFICIENCY=1.0D0/(1.0D0+FACTOR1+FACTOR2)*EPSILON

!  Graphite grain properties
   GRAPHITE_MU=0.005D0   ! Fraction of newly formed H2 that stays on the grain surface
   GRAPHITE_E_S=260.0D0  ! Energy of the saddle point between a physisorbed and a chemisorbed site (K)
   GRAPHITE_E_H2=520.0D0 ! Desorption energy of H2 molecules (K)
   GRAPHITE_E_HP=800.0D0 ! Desorption energy of physisorbed H atoms (K)
   GRAPHITE_E_HC=3.0D4   ! Desorption energy of chemisorbed H atoms (K)
   GRAPHITE_NU_H2=3.0D12 ! Vibrational frequency of H2 molecules in surface sites (s^-1)
   GRAPHITE_NU_HC=1.3D13 ! Vibrational frequency of H atoms in their surface sites (s^-1)

   FACTOR1=GRAPHITE_MU*FLUX/(2*GRAPHITE_NU_H2*EXP(-GRAPHITE_E_H2/GRAIN_TEMPERATURE))

   FACTOR2=1.0D0*(1.0D0+SQRT((GRAPHITE_E_HC-GRAPHITE_E_S)/(GRAPHITE_E_HP-GRAPHITE_E_S)))**2 &
        & /4.0D0*EXP(-GRAPHITE_E_S/GRAIN_TEMPERATURE)

   EPSILON=1.0D0/(1.0D0+GRAPHITE_NU_HC/(2*FLUX)*EXP(-1.5*GRAPHITE_E_HC/GRAIN_TEMPERATURE) &
              & *(1.0D0+SQRT((GRAPHITE_E_HC-GRAPHITE_E_S)/(GRAPHITE_E_HP-GRAPHITE_E_S)))**2)

   GRAPHITE_FORMATION_EFFICIENCY=1.0D0/(1.0D0+FACTOR1+FACTOR2)*EPSILON

!!$!  Use the tradional rate, with a simple temperature dependence based on the
!!$!  thermal velocity of the H atoms in the gas and neglecting any temperature
!!$!  dependency of the formation and sticking efficiencies
!!$   RATE=3.0D-18*SQRT(GAS_TEMPERATURE)

!!$!  Use the treatment of de Jong (1977, A&A, 55, 137, p140 right column).
!!$!  The second exponential dependence on the gas temperature reduces the
!!$!  efficiency at high temperatures and so prevents runaway H2 formation
!!$!  heating at high temperatures:
!!$!
!!$!  k_H2 = 3E-18 * T^0.5 * exp(-T/1000)   [cm3/s]
!!$!
!!$   RATE=3.0D-18*SQRT(GAS_TEMPERATURE)*EXP(-(GAS_TEMPERATURE/1.0D3))

!!$!  Use the treatment of Tielens & Hollenbach (1985, ApJ, 291, 722, eqn 4)
!!$   RATE=0.5D0*THERMAL_VELOCITY*TOTAL_CROSS_SECTION*STICKING_COEFFICIENT

!  Use the treatment of Cazaux & Tielens (2002, ApJ, 575, L29) and
!  Cazaux & Tielens (2004, ApJ, 604, 222)
   RATE=0.5D0*THERMAL_VELOCITY*(SILICATE_CROSS_SECTION*SILICATE_FORMATION_EFFICIENCY &
    & + GRAPHITE_CROSS_SECTION*GRAPHITE_FORMATION_EFFICIENCY)*STICKING_COEFFICIENT*METALLICITY*100./g2d

   !RATE = RATE / 3.0

!!$!  Use the expression given by Markus Rollig during the February 2012 Leiden workshop
!!$   RATE=0.5D0*THERMAL_VELOCITY &
!!$      & *(SILICATE_CROSS_SECTION/((1.0D0 + 6.998D24/EXP(1.5*SILICATE_E_HC/GRAIN_TEMPERATURE)) &
!!$      & *(1.0D0 + 1.0D0/(EXP(SILICATE_E_HP/GRAIN_TEMPERATURE) &
!!$      & *(0.427D0/EXP((SILICATE_E_HP-SILICATE_E_S)/GRAIN_TEMPERATURE) + 2.5336D-14*SQRT(GRAIN_TEMPERATURE))))) &
!!$      & + GRAPHITE_CROSS_SECTION/((1.0D0 + 4.610D24/EXP(1.5*GRAPHITE_E_HC/GRAIN_TEMPERATURE)) &
!!$      & *(1.0D0 + 1.0D0/(EXP(GRAPHITE_E_HP/GRAIN_TEMPERATURE) &
!!$      & *(0.539D0/EXP((GRAPHITE_E_HP-GRAPHITE_E_S)/GRAIN_TEMPERATURE) + 5.6334D-14*SQRT(GRAIN_TEMPERATURE)))))) &
!!$      & *STICKING_COEFFICIENT

   RETURN
END FUNCTION H2_FORMATION_RATE
!=======================================================================
