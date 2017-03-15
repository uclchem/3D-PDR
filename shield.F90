!C=======================================================================
!C
!C     H2 photodissociation rate taking into account
!C     self-shielding and grain extinction
!C
!C-----------------------------------------------------------------------
!C
!C     Input parameters:
!C     K0  = Unattenuated photodissociation rate (in cm^3/s)
!C     G0  = Incident FUV field (in Draine units)
!C     AV  = visual extinction (in magnitudes)
!C     NH2 = H2 column density (in cm^-2)
!C
!C     Program variables:
!C     H2PDRATE = H2 photodissociation rate taking into
!C                account self-shielding and grain extinction
!C     DOPW     = Doppler linewidth (in Hz) of a typical transition
!C                (assuming turbulent broadening with b=3 km/s)
!C     RADW     = radiative linewidth (in Hz) of a typical transition
!C     LAMBDA   = wavelength (in Å) of a typical transition
!C
!C     Functions called:
!C     H2SHIELD = H2 self-shielding function
!C     S!CATTER  = attenuation due to scattering by dust
!C
!C-----------------------------------------------------------------------
      FUNCTION H2PDRATE(K0,G0,AV,NH2)
!      IMPLICIT NONE
!      real(kind=dp) :: K0,G0,AV,NH2
!      real(kind=dp) :: LAMBDA,SCATTER
!!C      DOUBLE PRECISION DOPW,RADW
!!C      DOUBLE PRECISION H2SHIELD1
!      real(kind=dp) :: H2SHIELD2
     use definitions
     use healpix_types
     use maincode_module, only : v_turb
!     use global_module, only : nh2
    implicit none
     real(kind=dp) :: H2PDRATE
     real(kind=dp), intent(in) :: k0, g0, av
!     integer(kind=i4b), intent(in) :: nh2
     real(kind=dp), intent(in) :: nh2
     real(kind=dp) :: lambda, scatter!, h2shield2
     real(kind=dp) :: dopw, radw, h2shield1

      LAMBDA=1000.0D0
      DOPW=V_TURB/(LAMBDA*1.0D-8)
      RADW=8.0D7

!     Calculate the H2 photodissociation rate (H2PDRATE)
      H2PDRATE=K0*G0*H2SHIELD1(NH2,DOPW,RADW)*SCATTER(AV,LAMBDA)/2.0
!      H2PDRATE=K0*G0*H2SHIELD2(NH2)*SCATTER(AV,LAMBDA)/2.0

      RETURN
      END
!C=======================================================================

!C=======================================================================
!C
!C     !CO photodissociation rate taking into
!C     account shielding and grain extinction
!C
!C-----------------------------------------------------------------------
!C
!C     Input parameters:
!C     K0  = Unattenuated photodissociation rate (in cm^3/s)
!C     G0  = Incident FUV field (in Draine units)
!C     AV  = visual extinction (in magnitudes)
!C     N!CO = !CO column density (in cm^-2)
!C     NH2 = H2 column density (in cm^-2)
!C
!C     Program variables:
!C     !COPDRATE = !CO photodissociation rate taking into
!C                account self-shielding and grain extinction
!C     LAMBDA   = wavelength (in Å) of a typical transition
!C
!C     Functions called:
!C     LBAR     = function to determine the wavelength
!C     !COSHIELD = !CO shielding function
!C     S!CATTER  = attenuation due to scattering by dust
!C
!C-----------------------------------------------------------------------
      FUNCTION COPDRATE(K0,G0,AV,NCO,NH2)

!      IMPLICIT NONE
!      real(kind=dp) :: K0,G0,AV,NCO,NH2
!      real(kind=dp) :: LAMBDA,LBAR
!      real(kind=dp) :: COSHIELD,SCATTER

     use definitions
     use healpix_types
!     use global_module, only : nh2
     implicit none
     real(kind=dp) :: copdrate
     real(kind=dp), intent(in) :: k0, g0, av, nco
!     integer(kind=i4b), intent(in) :: nh2
     real(kind=dp), intent(in) :: nh2
     real(kind=dp) :: lambda, lbar, coshield, scatter


      LAMBDA=LBAR(NCO,NH2)

!C     Calculate the CO photodissociation rate (COPDRATE)
      COPDRATE=K0*G0*COSHIELD(NCO,NH2)*SCATTER(AV,LAMBDA)/2.0
      RETURN
      END
!C=======================================================================

!C=======================================================================
!C
!C     !CI photoionization rate taking into account grain extinction
!C     and shielding by !CI and H2 lines, adopting the treatment of
!C     Kamp & Bertoldi (2000, A&A, 353, 276, Equation 8)
!C
!C-----------------------------------------------------------------------
!C
!C     Input parameters:
!C     K0   = Unattenuated photoionization rate (in cm^3/s)
!C     G0   = Incident FUV field (in Draine units)
!C     AV   = visual extinction (in magnitudes)
!C     KAV  = tau(λ)/tau(V) correction factor
!C     N!CI  = !CI column density (in cm^-2)
!C     NH2  = H2 column density (in cm^-2)
!C     TGAS = gas temperature (in K)
!C
!C     Program variables:
!C     !CIPDRATE = !CI photoionization rate taking into
!C                account shielding and grain extinction
!C     TAU!C     = optical depth in the !CI absorption band
!C
!C-----------------------------------------------------------------------
       FUNCTION CIPDRATE(K0,G0,AV,KAV,NCI,NH2,TGAS)

!      IMPLICIT NONE
!      real(kind=dp) :: K0,G0,AV,KAV,NCI,NH2,TGAS
!      real(kind=dp) :: TAUC

     use definitions
     use healpix_types
!     use global_module, only : nh2
     implicit none
     real(kind=dp) :: cipdrate
     real(kind=dp), intent(in) :: K0,G0,AV,KAV,NCI,TGAS
!     integer(kind=i4b), intent(in) :: nh2
     real(kind=dp), intent(in) :: nh2
     real(kind=dp) :: tauc


!C     !Calculate the optical depth in the !CI absorption band, accounting
!C     for grain extinction and shielding by !CI and overlapping H2 lines
      TAUC=KAV*AV+1.1D-17*NCI+(0.9D0*TGAS**0.27D0*(NH2/1.59D21)**0.45D0)

!C     Calculate the CI photoionization rate (CIPDRATE)
      CIPDRATE=K0*G0*EXP(-TAUC)/2.0

      RETURN
      END
!C=======================================================================

!C=======================================================================
!C
!C     SI photoionization rate -- needs to be implemented!
!C     For now, use the standard expression for photorates
!C
!C-----------------------------------------------------------------------
!C
!C     Input parameters:
!C     K0   = Unattenuated photoionization rate (in cm^3/s)
!C     G0   = Incident FUV field (in Draine units)
!C     AV   = visual extinction (in magnitudes)
!C     KAV  = tau(λ)/tau(V) correction factor
!C     NSI  = SI column density (in cm^-2)
!C
!C     Program variables:
!C     SIPDRATE = SI photoionization rate taking into
!C                account shielding and grain extinction
!C     TAUS     = optical depth in the SI absorption band
!C
!C-----------------------------------------------------------------------
      FUNCTION SIPDRATE(K0,G0,AV,KAV)!,NSI)

!      IMPLICIT NONE
!      real(kind=dp) K0,G0,AV,KAV,NSI
!      real(kind=dp) TAUS

     use definitions
     use healpix_types
     implicit none
     real(kind=dp) :: sipdrate
     real(kind=dp), intent(in) :: K0,G0,AV,KAV!,NSI
     real(kind=dp) :: taus
!C     Calculate the optical depth in the SI absorption band, accounting
!C     for grain extinction and shielding by ???
      TAUS=KAV*AV

!C     Calculate the SI photoionization rate (SIPDRATE)
      SIPDRATE=K0*G0*EXP(-TAUS)/2.0

      RETURN
      END
!C=======================================================================

!C=======================================================================
!C
!C     H2 line self-shielding, adopting the treatment of
!C     Federman, Glassgold & Kwan (1979, ApJ, 227, 466)
!C
!C-----------------------------------------------------------------------
!C
!C     Input parameters:
!C     NH2  = H2 column density (in cm^-2)
!C     DOPW = Doppler linewidth (in Hz)
!C     RADW = radiative linewidth (in Hz)
!C
!C     Program variables:
!C     H2SHIELD1 = total self-shielding function containing
!C                 both Doppler and radiative contributions
!C     FPARA     = fraction of H2 in para state: 1/(1+o/p ratio)
!C     FOS!C      = oscillator strength of a typical transition
!C     TAUD      = parameter tauD (eq. A7) in Federman's paper
!C                 (optical depth at line centre)
!C     R         = parameter r  (eq. A2) in Federman's paper
!C     T         = parameter t1 (eq. A6) in Federman's paper
!C     U         = parameter u1 (eq. A6) in Federman's paper
!C     JD        = parameter JD (eq. A8) in Federman's paper
!C                 (Doppler contribution to self-shielding)
!C     JR        = parameter JR (eq. A9) in Federman's paper
!C                (radiative contribution to self-shielding)
!C
!C-----------------------------------------------------------------------
      FUNCTION H2SHIELD1(NH2,DOPW,RADW)

!      IMPLICIT NONE
!      real(kind=dp) ::  NH2,DOPW,RADW
!      real(kind=dp) ::  FPARA,FOSC,TAUD
!      real(kind=dp) ::  R,T,U,JD,JR

     use definitions
     use healpix_types
!     use global_module, only : nh2
     implicit none
     real(kind=dp) :: h2shield1
!     integer(kind=i4b), intent(in) :: nh2
     real(kind=dp), intent(in) :: nh2
     real(kind=dp), intent(in) ::DOPW,RADW
     real(kind=dp) :: FPARA, FOSC, TAUD, R, T, U, JD, JR


!C     Calculate the optical depth at line centre = N(H2)*f_para*(πe^2/mc)*f/(√πß) ≈ N(H2)*f_para*(1.5E-2)*f/ß
      FPARA=0.5D0 ! (assume o/p ratio=1)
      FOSC=1.0D-2
      TAUD=NH2*FPARA*(1.497358985D-2)*FOSC/DOPW

!C     Calculate the Doppler core contribution to the self-shielding (JD)
      IF(TAUD.EQ.0.0D0) THEN
         JD=1.0D0
      ELSE IF(TAUD.LT.2.0D0) THEN
         JD=EXP(-(0.666666667D0*TAUD))
      ELSE IF(TAUD.LT.10.0D0) THEN
         JD=0.638D0*TAUD**(-1.25D0)
      ELSE IF(TAUD.LT.100.0D0) THEN
         JD=0.505D0*TAUD**(-1.15D0)
      ELSE
         JD=0.344D0*TAUD**(-1.0667D0)
      ENDIF

!C     Calculate the radiative wing contribution to self-shielding (JR)
      IF(RADW.EQ.0.0D0) THEN
         JR=0.0D0
      ELSE
         R=RADW/(1.772453851D0*DOPW)
         T=3.02D0*((R*1.0D3)**(-0.064D0))
         U=SQRT(TAUD*R)/T
         JR=R/(T*SQRT(0.785398163D0+U**2))
      ENDIF

!C     Calculate the total self-shielding function (H2SHIELD1)
      H2SHIELD1=JD+JR

      RETURN
      END
!C=======================================================================

!C=======================================================================
!C
!C     H2 line shielding, using the computed values listed in
!C     Lee et al. (1996, A&A, 311, 690, Table 10)
!C
!C-----------------------------------------------------------------------
!C
!C     Input parameters:
!C     NH2 = H2 column density (in cm^-2)
!C
!C     Program variables:
!C     H2SHIELD2 = total H2 shielding factor containing
!C                 contributions from both H2 and H lines
!C                 from spline interpolation over the grid
!C     SH2_GRID  = H2 shielding factors from Lee et al. (1996)
!C                 as a function of H2 column density
!C     SH2_DERIV = 2nd derivative of SH2_GRID values from SPLINE
!C     !COL_GRID  = H2 column densities (in cm^-2)
!C     NUMH2     = number of entries in the table
!C     START     = .TRUE. when H2SHIELD2 is first called
!C
!C     Functions called:
!C     SPLINE =
!C     SPLINT =
!C
!C-----------------------------------------------------------------------
      FUNCTION H2SHIELD2(NH2)

!      IMPLICIT NONE
!      LOGICAL :: START
!      INTEGER(kind=i4b) :: NUMH2
!      real(kind=dp) ::  NH2
!      real(kind=dp) ::  COL_GRID(105),SH2_GRID(105),SH2_DERIV(105)
!      COMMON /STATUS/START
!      COMMON /H2GRID/SH2_GRID,SH2_DERIV,COL_GRID,NUMH2
     use definitions
     use healpix_types
     use uclpdr_module, only : start, numh2, COL_GRID, SH2_GRID, SH2_DERIV
!     use global_module, only : nh2
     implicit none
     real(kind=dp) :: h2shield2
!     integer(kind=i4b), intent(in) :: nh2
     real(kind=dp), intent(inout) :: nh2

      IF(START) CALL SPLINE(COL_GRID,SH2_GRID,NUMH2, &
     &                      1.0D30,1.0D30,SH2_DERIV)
      IF(NH2.LT.COL_GRID(1))     NH2=COL_GRID(1)
      IF(NH2.GT.COL_GRID(NUMH2)) NH2=COL_GRID(NUMH2)
      CALL SPLINT(COL_GRID,SH2_GRID,SH2_DERIV,NUMH2,NH2,H2SHIELD2)
      IF(H2SHIELD2.LT.0.0D0) H2SHIELD2=0.0D0

      RETURN
      END
!C=======================================================================

!C=======================================================================
!C
!C     12!CO line shielding, using the computed values listed in
!C     van Dishoeck & Black (1988, ApJ, 334, 771, Table 5)
!C
!C     Appropriate shielding factors are determined by performing a
!C     2-dimensional spline interpolation over the values listed in
!C     Table 5 of van Dishoeck & Black, which include contributions
!C     from self-shielding and H2 screening
!C
!C-----------------------------------------------------------------------
!C
!C     Input parameters:
!C     N!CO = !CO column density (in cm^-2)
!C     NH2 = H2 column density (in cm^-2)
!C
!C     Program variables:
!C     !COSHIELD  = total 12!CO shielding factor containing
!C                 contributions from both H2 and !CO lines
!C                 from 2D spline interpolation over the grid
!C     S!CO_GRID  = log10 values of the 12!CO shielding factors
!C                 from van Dishoeck & Black (1988) as a function
!C                 of !CO column density (1st index) and H2 column
!C                 density (2nd index)
!C     S!CO_DERIV = 2nd derivative of S!CO_GRID values from SPLIE2
!C     N!CO_GRID  = log10 values of !CO column densities (in cm^-2)
!C     NH2_GRID  = log10 values of H2 column densities (in cm^-2)
!C     DIM!CO     = number of !CO column densities
!C     DIMH2     = number of H2 column densities
!C     START     = .TRUE. when !COSHIELD is first called
!C
!C     Functions called:
!C     SPLIE2 =
!C     SPLIN2 =
!C
!C-----------------------------------------------------------------------
       FUNCTION COSHIELD(NCO,NH2)

!      IMPLICIT NONE
!      LOGICAL :: START
!      INTEGER(kind=i4b) :: DIMCO,DIMH2
!      real(kind=dp) ::  NCO,NH2
!      real(kind=dp) ::  LOGNCO,LOGNH2
!      real(kind=dp) ::  NCO_GRID(8),NH2_GRID(6)
!      real(kind=dp) ::  SCO_GRID(8,6),SCO_DERIV(8,6)
!      COMMON /STATUS/START
!      COMMON /COGRID/SCO_GRID,SCO_DERIV,NCO_GRID,NH2_GRID,DIMCO,DIMH2

     use definitions
     use healpix_types
!     use global_module, only : NCO, NH2
     use uclpdr_module, only : start, NCO_GRID, NH2_GRID, SCO_GRID, &
            & DIMCO, SCO_DERIV, DIMH2, SCO_DERIV
     implicit none
     real(kind=dp) :: COSHIELD
     real(kind=dp) :: LOGNCO, LOGNH2
!     integer(kind=i4b), intent(in) :: NCO, NH2
     real(kind=dp), intent(in) :: nh2, nco

      IF(START) THEN
         CALL SPLIE2(NCO_GRID,NH2_GRID,SCO_GRID,DIMCO,DIMH2,SCO_DERIV)
         START=.FALSE.
      ENDIF

      LOGNCO=DLOG10(NCO+1.0D0)
      LOGNH2=DLOG10(NH2+1.0D0)

      IF(LOGNCO.LT.NCO_GRID(1)) LOGNCO=NCO_GRID(1)
      IF(LOGNH2.LT.NH2_GRID(1)) LOGNH2=NH2_GRID(1)
      IF(LOGNCO.GT.NCO_GRID(DIMCO)) LOGNCO=NCO_GRID(DIMCO)
      IF(LOGNH2.GT.NH2_GRID(DIMH2)) LOGNH2=NH2_GRID(DIMH2)

      CALL SPLIN2(NCO_GRID,NH2_GRID,SCO_GRID,SCO_DERIV,&
     &            DIMCO,DIMH2,LOGNCO,LOGNH2,COSHIELD)
      COSHIELD=10.0D0**COSHIELD

      RETURN
      END
!C=======================================================================

!C=======================================================================
!C
!C     Scattering by dust grains, adopting the treatment of
!C     Wagenblast & Hartquist (1989, MNRAS, 237, 1019) and
!C     Flannery, Roberge & Rybicki (1980, ApJ, 236, 598)
!C
!C-----------------------------------------------------------------------
!C
!C     Input parameters:
!C     AV     = visual extinction (in magnitudes)
!C     LAMBDA = wavelength (in Å) of incident radiation
!C
!C     Program variables:
!C     S!CATTER = attenuation factor describing the influence of
!C               grain scattering on the FUV flux, dependening
!C               on the total column density and wavelength of
!C               light (assuming albedo=0.3 gscat=0.8)
!C     TAUV    = optical depth at visual wavelength (λ=5500Å)
!C     TAUL    = optical depth at wavelength LAMBDA
!C     A(0)    = a(0)*exp(-k(0)*tau)
!C             = relative intensity decrease for 0 < tau < 1
!C     A(I)    = ∑ a(i)*exp(-k(i)*tau) for i=1,5
!C               relative intensity decrease for tau ≥ 1
!C     K(0)    = see A0
!C     K(I)    = see A(I)
!C
!C     Functions called:
!C     XLAMBDA = function to determine tau(λ)/tau(V)
!C
!C-----------------------------------------------------------------------
      FUNCTION SCATTER(AV,LAMBDA)

!      IMPLICIT NONE
!      INTEGER(kind=i4b) :: I
!      real(kind=dp) ::  AV,LAMBDA
!      real(kind=dp) ::  TAUV,TAUL
!      real(kind=dp) ::  A(0:5),K(0:5),EXPONENT
!      real(kind=dp) ::  XLAMBDA
!      DATA A/1.000D0,2.006D0,-1.438D0,0.7364D0,-0.5076D0,-0.0592D0/
!      DATA K/0.7514D0,0.8490D0,1.013D0,1.282D0,2.005D0,5.832D0/

     use definitions
     use healpix_types
     implicit none
     real(kind=dp) :: scatter, LAMBDAVAR
     real(kind=dp), intent(in) :: AV, LAMBDA
     real(kind=dp), dimension(0:5), save :: A = (/&
           &1.000D0,2.006D0,-1.438D0,0.7364D0,-0.5076D0,-0.0592D0/)
     real(kind=dp), dimension(0:5), save :: K = (/&
           &0.7514D0,0.8490D0,1.013D0,1.282D0,2.005D0,5.832D0/)
     real(kind=dp) :: EXPONENT, XLAMBDA
     integer(kind=i4b) :: i
     real(kind=dp) :: TAUL, TAUV

!C     Calculate the optical depth at visual wavelength
      TAUV=AV/1.086D0

!C     Convert the optical depth to that at the desired wavelength
      LAMBDAVAR=LAMBDA
      TAUL=TAUV*XLAMBDA(LAMBDAVAR)

!C     Calculate the attenuation due to scattering by dust (SCATTER)
      SCATTER=0.0D0
      IF(TAUL.LT.1.0D0) THEN
         EXPONENT=K(0)*TAUL
         IF(EXPONENT.LT.100.0D0) THEN
            SCATTER=A(0)*EXP(-EXPONENT)
         ENDIF
      ELSE
         DO I=1,5
            EXPONENT=K(I)*TAUL
            IF(EXPONENT.LT.100.0D0) THEN
               SCATTER=SCATTER+A(I)*EXP(-EXPONENT)
            ENDIF
         ENDDO
      ENDIF

      RETURN
      END
!C=======================================================================

!C=======================================================================
!C
!C     Determine the ratio of the optical depth at a given wavelength to
!C     that at visual wavelength (λ=5500Å) using the extinction curve of
!C     Savage & Mathis (1979, ARA&A, 17, 73, Table 2)
!C
!C-----------------------------------------------------------------------
!C
!C     Input parameters:
!C     LAMBDA  = wavelength (in Å)
!C
!C     Program variables:
!C     XLAMBDA = value of tau(λ)/tau(V) at the desired wavelength
!C               (by spline interpolation over a table of values)
!C     L_GRID  = wavelengths listed in Table 2 of Savage & Mathis
!C     X_GRID  = tau(λ)/tau(V) values, determined by dividing the
!C               Aλ/E(B-V) values in Table 2 by R=AV/E(B-V)=3.1
!C     X_DERIV = 2nd derivative of X_GRID values from SPLINE
!C     N_GRID  = number of wavelengths
!C
!C     Functions called:
!C     SPLIE =
!C     SPLIN =
!C
!C-----------------------------------------------------------------------
      FUNCTION XLAMBDA(LAMBDA)

!      IMPLICIT NONE
!      LOGICAL :: START
!      INTEGER(kind=i4b) ::  N_GRID
!      real(kind=dp) ::  LAMBDA
!      real(kind=dp) ::  L_GRID(30),X_GRID(30),X_DERIV(30)
!      COMMON /STATUS/START
!      COMMON /TAUGRID/L_GRID,X_GRID,X_DERIV,N_GRID

     use definitions
     use healpix_types
     use uclpdr_module, only : start, N_GRID, L_GRID, X_GRID, X_DERIV
     implicit none
     real(kind=dp) :: xlambda
     real(kind=dp), intent(inout) :: lambda

!C     Find the appropriate value for XLAMBDA using spline interpolation
      IF(START) CALL SPLINE(L_GRID,X_GRID,N_GRID,1.0D30,1.0D30,X_DERIV)
      IF(LAMBDA.LT.L_GRID(1))      LAMBDA=L_GRID(1)
      IF(LAMBDA.GT.L_GRID(N_GRID)) LAMBDA=L_GRID(N_GRID)
      CALL SPLINT(L_GRID,X_GRID,X_DERIV,N_GRID,LAMBDA,XLAMBDA)

      RETURN
      END
!C=======================================================================

!C=======================================================================
!C
!C     !Calculate the mean wavelength (in Å) of the 33 dissociating bands,
!C     weighted by their fractional contribution to the total shielding
!C     van Dishoeck & Black (1988, ApJ, 334, 771, Equation 4)
!C
!C-----------------------------------------------------------------------
!C
!C     Input parameters:
!C     N!CO = !CO column density (in cm^-2)
!C     NH2 = H2 column density (in cm^-2)
!C
!C     Program variables:
!C     LBAR = mean wavelength (in Å)
!C     U    = log10(N!CO)
!C     W    = log10(NH2)
!C
!C-----------------------------------------------------------------------
      FUNCTION LBAR(NCO,NH2)

!      IMPLICIT NONE
!      real(kind=dp) ::  NCO,NH2
!      real(kind=dp) ::  U,W


     use definitions
     use healpix_types
!     use global_module, only : NCO, NH2
     implicit none
     real(kind=dp) :: lbar
     real(kind=dp) :: U,W
!     integer(kind=i4b), intent(in) :: NCO, NH2
     real(kind=dp), intent(in) :: nh2,nco

      U=DLOG10(NCO+1.0D0)
      W=DLOG10(NH2+1.0D0)

      LBAR=(5675.0D0 - 200.6D0*W) &
     &    - (571.6D0 - 24.09D0*W)*U &
     &   + (18.22D0 - 0.7664D0*W)*U**2

!C     LBAR cannot be larger than the wavelength of band 33 (1076.1Å)
!C     and cannot be smaller than the wavelength of band 1 (913.6Å)
      IF(LBAR.LT.913.6D0)  LBAR=913.6D0
      IF(LBAR.GT.1076.1D0) LBAR=1076.1D0

      RETURN
      END
!C=======================================================================
