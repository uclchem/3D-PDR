
!=======================================================================
!
!     Calculate the total heating rate at the current grid point.
!
!-----------------------------------------------------------------------
      SUBROUTINE CALC_HEATING(DENSITY,GAS_TEMPERATURE,DUST_TEMPERATURE, &
                 & UV_FIELD,V_TURB,NSPEC,ABUNDANCE,NREAC,RATE,HEATING_RATE,&
                 & NRGR,NRH2,NRHD,NRCO,NRCI,NRSI) 

      USE DEFINITIONS
      USE HEALPIX_TYPES
      USE GLOBAL_MODULE
      USE MAINCODE_MODULE, ONLY : UV_FAC, ZETA

      IMPLICIT NONE

      INTEGER(KIND=I4B), INTENT(IN) :: NSPEC,NREAC
      REAL(KIND=DP), INTENT(IN)     :: DENSITY,GAS_TEMPERATURE,DUST_TEMPERATURE,UV_FIELD,V_TURB
      REAL(KIND=DP), INTENT(IN)     :: ABUNDANCE(1:NSPEC),RATE(1:NREAC)
      INTEGER(KIND=I4B), INTENT(IN) :: NRGR,NRH2,NRHD,NRCO,NRCI,NRSI
      REAL(KIND=DP), INTENT(OUT)    :: HEATING_RATE(1:12)

      REAL(KIND=DP) :: TOTAL_HEATING,PHOTOELECTRIC_HEATING,PAHPHOTOELEC_HEATING,WEINGARTNER_HEATING, &
                     & CIONIZATION_HEATING,H2FORMATION_HEATING,H2PHOTODISS_HEATING,FUVPUMPING_HEATING, &
                     & COSMICRAY_HEATING,TURBULENT_HEATING,CHEMICAL_HEATING,GASGRAIN_HEATING,SOFTXRAY_HEATING

!     Dust PE heating
      INTEGER(KIND=I4B) :: ITERATION
      REAL(KIND=DP) :: HABING_FIELD
      REAL(KIND=DP) :: X,XX,XK,XD,GAMMA,DELTA
      REAL(KIND=DP) :: DELTAD,DELTAUV,Y,HNUD,HNUH

!     PAH heating/cooling
      REAL(KIND=DP) :: EPSILON,ALPHA,BETA,PAH_HEATING,PAH_COOLING
      REAL(KIND=DP) :: PHI_PAH

!     Weingartner & Draine treatment of photoelectric heating (grains+PAHs)
      REAL(KIND=DP) :: C0,C1,C2,C3,C4,C5,C6

!     H2* FUV pumping heating
      REAL(KIND=DP) :: NCR_H2

!     Turbulent heating
      REAL(KIND=DP) :: L_TURB

!     Gas-grain coupling heating/cooling
      REAL(KIND=DP) :: ACCOMMODATION,NGRAIN,CGRAIN

!!     Soft X-ray heating
!      REAL(KIND=DP) :: PP1,PP2,F6
!     Convert the FUV field (in Draine units) to the Habing equivalent
      HABING_FIELD=1.68D0*UV_FIELD

!-----------------------------------------------------------------------
!     Dust photoelectric heating
!
!     Use the treatment of Tielens & Hollenbach, 1985, ApJ, 291, 722,
!     which follows de Jong (1977,1980)
!
!     The charge of a dust grain can be found by equating the rate of
!     photo-ejection of electrons from the dust grain to the rate of
!     recombination of electrons with the dust grain (Spitzer)
!
!     The various parameter values are taken from Table 2 of the paper
!-----------------------------------------------------------------------

      DELTAD=1.0D0
      DELTAUV=1.8D0
      Y=0.1D0
      HNUD=6.0D0
      HNUH=13.6D0

      XK=KB*GAS_TEMPERATURE/(HNUH*EV)
      XD=HNUD/HNUH
      GAMMA=2.9D-4*Y*SQRT(GAS_TEMPERATURE)*HABING_FIELD/(ABUNDANCE(NELECT)*DENSITY)
      DELTA=XK-XD+GAMMA

!     Iterate to determine X by finding the zero of the function F
      X=0.5D0
      DO ITERATION=1,100
         XX=X-(F(X,DELTA,GAMMA)/FF(X,DELTA))
         IF(ABS(XX-X).LT.1.0D-2) EXIT
         X=XX
      ENDDO
      X=XX

     IF(ITERATION.GE.100) THEN
        WRITE(10,*)'WARNING! Grain parameter X not found in PE heating'
        WRITE(10,*)'Using final value from interation loop: X =',X
     ENDIF

!     Assume the dust PE heating scales linearly with metallicity
      PHOTOELECTRIC_HEATING=2.7D-25*DELTAUV*DELTAD*DENSITY*Y*HABING_FIELD &
             & *(((1.0D0-X)**2)/X + XK*((X**2)-1.0D0)/(X**2))*METALLICITY


!==================PAH PHOTOELECTRIC HEATING===============================
!-----------------------------------------------------------------------
!  Grain + PAH photoelectric heating (MRN size distribution; r = 3-100
!  Å)
!
!  Use the treatment of Bakes & Tielens (1994, ApJ, 427, 822) with the
!  modifications suggested by Wolfire et al. (2003, ApJ, 587, 278) to
!  account for the revised PAH abundance estimate from Spitzer data.
!
!  See also:
!  Wolfire et al. (1995, ApJ, 443, 152)
!  Le Page, Snow & Bierbaum (2001, ApJS, 132, 233)
!-----------------------------------------------------------------------

!  Adopt the PAH rate scaling factor of Wolfire et al. (2008, ApJ, 680,
!  384)
!  Setting this factor to 1.0 gives the standard Bakes & Tielens
!  expression
   PHI_PAH=1.0D0!0.4D0

   ALPHA=0.944D0
   BETA=0.735D0/GAS_TEMPERATURE**0.068
   DELTA=HABING_FIELD*SQRT(GAS_TEMPERATURE)/(ABUNDANCE(NELECT)*DENSITY*PHI_PAH)
   EPSILON=4.87D-2/(1.0D0+4.0D-3*DELTA**0.73) + 3.65D-2*(GAS_TEMPERATURE/1.0D4)**0.7/(1.0D0+2.0D-4*DELTA)

   PAH_HEATING=1.30D-24*EPSILON*HABING_FIELD*DENSITY
   PAH_COOLING=4.65D-30*GAS_TEMPERATURE**ALPHA*(DELTA**BETA)*ABUNDANCE(NELECT)*DENSITY*PHI_PAH*DENSITY

!  Assume the PE heating rate scales linearly with metallicity
   PAHPHOTOELEC_HEATING=(PAH_HEATING-PAH_COOLING)*METALLICITY
!
!!==============================================================================
!

!-----------------------------------------------------------------------
!     Weingartner & Draine, 2001, ApJS, 134, 263
!
!     Includes photoelectric heating due to PAHs, VSGs and larger grains
!     Assumes a gas-to-dust mass ratio of 100:1
!-----------------------------------------------------------------------

      C0=5.72D+0
      C1=3.45D-2
      C2=7.08D-3
      C3=1.98D-2
      C4=4.95D-1
      C5=6.92D-1
      C6=5.20D-1

      WEINGARTNER_HEATING=METALLICITY*1.0D-26*(HABING_FIELD*DENSITY)*(C0+C1*GAS_TEMPERATURE**C4) &
             & /(1.0D0+C2*(HABING_FIELD*SQRT(GAS_TEMPERATURE)/(ABUNDANCE(NELECT)*DENSITY))**C5  &
             & *(1.0D0+C3*(HABING_FIELD*SQRT(GAS_TEMPERATURE)/(ABUNDANCE(NELECT)*DENSITY))**C6))

!-----------------------------------------------------------------------
!     Carbon photoionization heating
!
!     1 eV on average per carbon ionization
!     Use the C photoionization rate determined by the CALCULATE_REACTION_RATES subroutine: RATE(NRCI) [units: s^-1]
!-----------------------------------------------------------------------
      CIONIZATION_HEATING=(1.0*EV)*RATE(NRCI)*ABUNDANCE(NC)*DENSITY
!-----------------------------------------------------------------------
!     H2 formation heating
!
!     Assume 1.5 eV liberated as heat during H2 formation
!     See: Hollenbach & Tielens, Review of Modern Physics, 1999, 71, 173
!     Use the rate determined by the CALCULATE_REACTION_RATES subroutine: RATE(NRGR) [units: cm^3.s^-1]
!-----------------------------------------------------------------------

      H2FORMATION_HEATING=(1.5*EV)*RATE(NRGR)*DENSITY*ABUNDANCE(NH)*DENSITY

!-----------------------------------------------------------------------
!     H2 photodissociation heating
!
!     0.4 eV on average per photodissociated molecule
!     Use the rate determined by the CALCULATE_REACTION_RATES subroutine: RATE(NRH2) [units: s^-1]
!-----------------------------------------------------------------------

      H2PHOTODISS_HEATING=(0.4*EV)*RATE(NRH2)*ABUNDANCE(NH2)*DENSITY

!-----------------------------------------------------------------------
!     H2 FUV pumping heating
!
!     2.2 eV on average per vibrationally excited H2* molecule
!     See: Hollenbach & McKee (1979)
!     Use the H2 photodissociation rate determined by CALCULATE_REACTION_RATES: RATE(NRH2) [units: s^-1]
!     Use the H2 critical density calculation from Hollenbach & McKee (1979)
!-----------------------------------------------------------------------

      NCR_H2=1.0D6/SQRT(GAS_TEMPERATURE)/(1.6D0*ABUNDANCE(NH)*EXP(-((400.0D0/GAS_TEMPERATURE)**2)) &
                                       & + 1.4D0*ABUNDANCE(NH2)*EXP(-(18100.0D0/(GAS_TEMPERATURE+1200.0D0))))

      FUVPUMPING_HEATING=(2.2*EV)*9.0D0*RATE(NRH2)*ABUNDANCE(NH2)*DENSITY/(1.0D0+NCR_H2/DENSITY)

!-----------------------------------------------------------------------
!     Cosmic-ray ionization heating
!
!     8.0 eV of heat deposited per primary ionization (plus some from He ionization)
!     Use the treatment of Tielens & Hollenbach, 1985, ApJ, 291, 772
!     See also: Shull & Van Steenberg, 1985, ApJ, 298, 268
!               Clavel et al. (1978), Kamp & van Zadelhoff (2001)
!-----------------------------------------------------------------------

!      COSMICRAY_HEATING=(20.0*EV)*(1.3D-17*ZETA)*DENSITY*ABUNDANCE(NH2) !20.0 -> 9.4 eV
      COSMICRAY_HEATING=(9.4*EV)*(1.3D-17*ZETA)*DENSITY*ABUNDANCE(NH2) !20.0 -> 9.4 eV
!      COSMICRAY_HEATING=(9.4*EV)*(1.3D-17*ZETA)*DENSITY

!-----------------------------------------------------------------------
!     Supersonic turbulent decay heating
!
!     Most relevant for the inner parsecs of galaxies (Black)
!     Black, in Interstellar Processes, 1987, p731
!     See also: Rodriguez-Fernandez et al., 2001, A&A, 365, 174
!
!     V_TURB = turbulent velocity (km/s); Galactic center ~ 15 km/s
!     L_TURB = turbulent scale length (pc); typically 5 pc
!-----------------------------------------------------------------------

      L_TURB=5.0D0
      TURBULENT_HEATING=3.5D-28*((V_TURB/1.0D5)**3)*(1.0D0/L_TURB)*DENSITY

!-----------------------------------------------------------------------
!     Exothermic chemical reaction heating
!
!     See: Clavel et al., 1978, A&A, 65, 435
!     Recombination reactions: HCO+ (7.51 eV); H3+ (4.76+9.23 eV); H3O+ (1.16+5.63+6.27 eV)
!     Ion-neutral reactions  : He+ + H2 (6.51 eV); He+ + CO (2.22 eV)
!     For each reaction, the heating rate should be: n(1) * n(2) * K * E
!     with n(1) and n(2) the densities, K the rate coefficient [cm^3.s^-1], and E the energy [erg]
!-----------------------------------------------------------------------

#ifdef REDUCED
      CHEMICAL_HEATING=ABUNDANCE(NH2x)*DENSITY*ABUNDANCE(NELECT)*RATE(216)*10.9*EV& !H2+ + e-
                   & + ABUNDANCE(NH2x)*DENSITY*ABUNDANCE(NH)*RATE(155)*0.94*EV& !H2+ + H
                   & + ABUNDANCE(NHCOx)*DENSITY*ABUNDANCE(NELECT)*DENSITY*(RATE(240)*(7.51*EV)) &                                         ! HCO+ + e-
                   & + ABUNDANCE(NH3x)*DENSITY*ABUNDANCE(NELECT)*DENSITY*(RATE(217)*(4.76*EV)+RATE(218)*(9.23*EV)) &                      ! H3+  + e-
                   & + ABUNDANCE(NH3Ox)*DENSITY*ABUNDANCE(NELECT)*DENSITY*(RATE(236)*(1.16*EV)+&
                   &RATE(237)*(5.63*EV)+RATE(238)*(6.27*EV)) &                                                                            ! H3O+ + e-
                    & + ABUNDANCE(NHEx)*DENSITY*ABUNDANCE(NH2)*DENSITY*(RATE(50)*(6.51*EV)+RATE(170)*(6.51*EV)) &                          ! He+  + H2
                   & + ABUNDANCE(NHEx)*DENSITY*ABUNDANCE(NCO)*DENSITY*(RATE(89)*(2.22*EV)+RATE(90)*(2.22*EV)+&
                   &RATE(91)*(2.22*EV))          ! He+  + CO
#endif

#ifdef FULL
      CHEMICAL_HEATING=ABUNDANCE(NHCOx)*DENSITY*ABUNDANCE(NELECT)*DENSITY*(RATE(730)*(7.51*EV)) &                                         ! HCO+ + e-
                   & + ABUNDANCE(NH3x)*DENSITY*ABUNDANCE(NELECT)*DENSITY*(RATE(706)*(4.76*EV)+RATE(705)*(9.23*EV)) &                      ! H3+  + e-
                   & + ABUNDANCE(NH3Ox)*DENSITY*ABUNDANCE(NELECT)*DENSITY*(RATE(717)*(1.16*EV)+&
                   &RATE(716)*(5.63*EV)+RATE(714)*(6.27*EV)) &                                                                            ! H3O+ + e-
                   & + ABUNDANCE(NHEx)*DENSITY*ABUNDANCE(NH2)*DENSITY*(RATE(1227)*(6.51*EV)+RATE(265)*(6.51*EV)) &                          ! He+  + H2
                   & + ABUNDANCE(NHEx)*DENSITY*ABUNDANCE(NCO)*DENSITY*(RATE(1541)*(2.22*EV))
#endif

#ifdef MYNETWORK
      STOP "CHEMICAL_HEATING function has to be declared at &
             & [sub_calculate_heating.F90] &
             If you are using the pre-set 'mynetwork' network comment &
             & out this STOP [sub_calculate_heating.F90]"
     CHEMICAL_HEATING=ABUNDANCE(NHCOx)*DENSITY*ABUNDANCE(NELECT)*DENSITY*(RATE(185)*(7.51*EV)) &                                         ! HCO+ + e-
                   & + ABUNDANCE(NH3x)*DENSITY*ABUNDANCE(NELECT)*DENSITY*(RATE(173)*(4.76*EV)+RATE(172)*(9.23*EV)) &                      ! H3+  + e-
                   & + ABUNDANCE(NH3Ox)*DENSITY*ABUNDANCE(NELECT)*DENSITY*(RATE(179)*(1.16*EV)+&
                   &RATE(178)*(5.63*EV)+RATE(176)*(6.27*EV)) &                                                                            ! H3O+ + e-
                   & + ABUNDANCE(NHEx)*DENSITY*ABUNDANCE(NH2)*DENSITY*(RATE(297)*(6.51*EV)+RATE(67)*(6.51*EV)) &                          ! He+  + H2
                   & + ABUNDANCE(NHEx)*DENSITY*ABUNDANCE(NCO)*DENSITY*(RATE(364)*(2.22*EV))
#endif

!-----------------------------------------------------------------------
!     Gas-grain collisional heating
!
!     Use the treatment of Burke & Hollenbach, 1983, ApJ, 265, 223, and
!     accommodation fitting formula of Groenewegen, 1994, A&A, 290, 531
!
!     Other relevant references:
!     Hollenbach & McKee, 1979, ApJS, 41,555
!     Tielens & Hollenbach, 1985, ApJ, 291,722
!     Goldsmith, 2001, ApJ, 557, 736
!
!     This process is insignificant for the energy balance of the dust
!     but can influence the gas temperature. If the dust temperature is
!     lower than the gas temperature, this becomes a cooling mechanism
!
!     In Burke & Hollenbach (1983) the factor:
!
!     (8*kb/(pi*mass_proton))**0.5*2*kb = 4.003D-12
!
!     This value has been used in the expression below
!-----------------------------------------------------------------------
!
      ACCOMMODATION=0.35D0*EXP(-SQRT((DUST_TEMPERATURE+GAS_TEMPERATURE)/5.0D2))+0.1D0
      NGRAIN=1.998D-12*DENSITY*METALLICITY*100./g2d
      CGRAIN=PI*GRAIN_RADIUS**2

      GASGRAIN_HEATING=4.003D-12*DENSITY*NGRAIN*CGRAIN*ACCOMMODATION*SQRT(GAS_TEMPERATURE) &
                    & *(DUST_TEMPERATURE-GAS_TEMPERATURE)

! Gas-grain collisional heating (Tielens 2005 ISM book)
!      GASGRAIN_HEATING=-1d-33*DENSITY**2*sqrt(GAS_TEMPERATURE)*(GAS_TEMPERATURE-DUST_TEMPERATURE)


!-----------------------------------------------------------------------
!     Soft X-ray heating
!
!     Use the treatment of from Wolfire et al., 1995, ApJ, 443, 152
!-----------------------------------------------------------------------

!      IF(DEPTH.EQ.0) THEN
!         PP1=0.0D0
!      ELSE
!         PP1=DLOG10(DENSITY*(DIST(DEPTH)-DIST(DEPTH-1))/1.0D18)
!      ENDIF
!
!      PP2=DLOG10(ABUNDANCE(NELECT))
!      F6=0.990D0-2.74D-3*PP2+1.13D-3*PP2**2
!
!      SOFTXRAY_HEATING=10.0D0**(F6*(-26.5D0-0.920D0*PP1+5.89D-2*PP1**2)
!     *                +F6*0.96D0*EXP(-(((PP1-0.38D0)/0.87D0)**2)))



!-----------------------------------------------------------------------
!     Total heating rate (sum of all contributions)
!-----------------------------------------------------------------------

      TOTAL_HEATING=&
!           & + PHOTOELECTRIC_HEATING &
           & + PAHPHOTOELEC_HEATING &
!           & + WEINGARTNER_HEATING &
           & + CIONIZATION_HEATING &
           & + H2FORMATION_HEATING &
           & + H2PHOTODISS_HEATING &
           & + FUVPUMPING_HEATING &
           & + COSMICRAY_HEATING &
           & + TURBULENT_HEATING &
           & + CHEMICAL_HEATING &
!           & + SOFTXRAY_HEATING &
           & + GASGRAIN_HEATING

      HEATING_RATE(1)=PHOTOELECTRIC_HEATING
      HEATING_RATE(2)=PAHPHOTOELEC_HEATING
      HEATING_RATE(3)=WEINGARTNER_HEATING
      HEATING_RATE(4)=CIONIZATION_HEATING
      HEATING_RATE(5)=H2FORMATION_HEATING
      HEATING_RATE(6)=H2PHOTODISS_HEATING
      HEATING_RATE(7)=FUVPUMPING_HEATING
      HEATING_RATE(8)=COSMICRAY_HEATING
      HEATING_RATE(9)=TURBULENT_HEATING
      HEATING_RATE(10)=CHEMICAL_HEATING
      HEATING_RATE(11)=GASGRAIN_HEATING
      HEATING_RATE(12)=TOTAL_HEATING

!-----------------------------------------------------------------------

CONTAINS ! Dust photoelectric heating functions...

!=======================================================================
!     X is the grain charge parameter and is the solution to F(X)=0
!-----------------------------------------------------------------------
      FUNCTION F(X,DELTA,GAMMA)

      USE DEFINITIONS
      USE HEALPIX_TYPES

      IMPLICIT NONE

      REAL(KIND=DP) :: F
      REAL(KIND=DP), INTENT(IN) :: X,DELTA,GAMMA

      F=(X**3)+DELTA*(X**2)-GAMMA

      END FUNCTION F
!-----------------------------------------------------------------------

!=======================================================================
!     FF(X) is the derivative of F(X) with respect to X
!-----------------------------------------------------------------------
      FUNCTION FF(X,DELTA)

      USE DEFINITIONS
      USE HEALPIX_TYPES

      IMPLICIT NONE

      REAL(KIND=DP) :: FF
      REAL(KIND=DP), INTENT(IN) :: X,DELTA

      FF=3*(X**2)+DELTA*(2*X)

      END FUNCTION FF
!-----------------------------------------------------------------------

      END SUBROUTINE CALC_HEATING
!=======================================================================



