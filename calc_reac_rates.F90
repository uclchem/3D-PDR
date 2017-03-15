!=======================================================================
!
!  Calculate the rate coefficients for all reactions at the specified
!  temperature and visual extinction A_V. The photodissociation of H2
!  and CO and the photoionization of CI and SI are treated separately
!  in detail (see the routines in photorates.f90). Multiple rates for
!  the same reaction (duplicates) are allowed in the ratefile and are
!  activated based on their minimum and maximum temperature specified
!  in that file. Negative gamma factors are ignored below the minimum
!  temperature at which the reaction rate is valid.
!
!  X-ray induced reaction rates are calculated following the detailed
!  treatment of Meijerink & Spaans (2005, A&A, 436, 397).
!
!-----------------------------------------------------------------------
      SUBROUTINE CALCULATE_REACTION_RATES(TEMPERATURE,DUST_TEMPERATURE,NRAYS,RAD_SURFACE, &
               & AV,COLUMN,NREAC,REACTANT,PRODUCT,ALPHA,BETA,GAMMA,RATE,RTMIN,RTMAX,DUPLICATE,NSPEC, &
               & NRGR,NRH2,NRHD,NRCO,NRCI,NRSI)

 use definitions
 use healpix_types
 use global_module
 use functions_module
 use maincode_module , only : zeta, AV_fac

      IMPLICIT NONE

      INTEGER(kind=i4b), intent(in) :: NRAYS, NSPEC
      real(kind=dp),intent(in) :: TEMPERATURE, DUST_TEMPERATURE
      real(kind=dp),intent(in) :: RAD_SURFACE(0:nrays-1),AV(0:nrays-1),COLUMN(0:nrays-1,1:nspec)
      INTEGER(kind=i4b),intent(in) :: NREAC,DUPLICATE(1:nreac)
      real(kind=dp),intent(in) :: ALPHA(1:nreac),BETA(1:nreac),GAMMA(1:nreac),RTMIN(1:nreac),RTMAX(1:nreac)
      real(kind=dp), intent(out) :: rate(1:nreac)
      INTEGER(kind=i4b),intent(out):: NRGR,NRH2,NRHD,NRCO,NRCI,NRSI
      CHARACTER(len=10),intent(in) ::  REACTANT(1:nreac,1:3),PRODUCT(1:nreac,1:4)

      real(kind=dp) :: PHI_PAH,CION,STICKING,FLUX,YIELD
      INTEGER(kind=i4b) :: I,J,K

!     Initialize the rate coefficients.
      RATE=0.0D0

!     Initialize the stored reaction numbers. If they are not assigned
!     subsequently, any attempt to access that reaction will generate an
!     error and the code will crash. This is a useful bug catch.
      NRGR=0
      NRH2=0
      NRHD=0
      NRCO=0
      NRCI=0
      NRSI=0

      DO I=1,NREAC
!        Determine the type of reaction
         IF(REACTANT(I,2).EQ."PHOTON") GOTO 1
         IF(REACTANT(I,2).EQ."CRP   ") GOTO 2
         IF(REACTANT(I,2).EQ."CRPHOT") GOTO 3
         IF(REACTANT(I,2).EQ."FREEZE") GOTO 4
         IF(REACTANT(I,2).EQ."ELFRZE") GOTO 5
         IF(REACTANT(I,2).EQ."CRH   ") GOTO 6
         IF(REACTANT(I,2).EQ."PHOTD ") GOTO 7
         IF(REACTANT(I,2).EQ."THERM ") GOTO 8
         IF(REACTANT(I,2)(1:1).EQ."#") GOTO 9
         IF(REACTANT(I,2).EQ."XRAY  ") GOTO 110
         IF(REACTANT(I,2).EQ."XRSEC ") GOTO 120
         IF(REACTANT(I,2).EQ."XRLYA ") GOTO 130
         IF(REACTANT(I,2).EQ."XRPHOT") GOTO 140

!-----------------------------------------------------------------------

!     Thermal reactions:

!     The rate of H2 formation on grains is calculated separately
!     by the function H2_FORMATION_RATE (see function for details)
         IF((REACTANT(I,1).EQ."H  " .AND. REACTANT(I,2).EQ."H  "   .AND. &
          & (REACTANT(I,3).EQ."   " .OR.  REACTANT(I,3).EQ."#  ")) .AND. &
          &  (PRODUCT(I,1).EQ."H2 " .AND. &
          &  (PRODUCT(I,2).EQ."   " .OR.  PRODUCT(I,2).EQ."#  "))) THEN
#ifdef H2FORM
            RATE(I)=H2_FORMATION_RATE(TEMPERATURE,DUST_TEMPERATURE)
#else
!            RATE(I)=3.0D-18*SQRT(TEMPERATURE)*EXP(-(TEMPERATURE/1.0D3))
            RATE(I)=3.0D-18*SQRT(TEMPERATURE)
#endif
            NRGR=I
            GOTO 10
         ENDIF

!     Rates for reactions involving PAHs are calculated according to the
!     treatment of Wolfire et al. (2003, ApJ, 587, 278; 2008, ApJ, 680, 384)
         IF(ANY(REACTANT(I,:).EQ."PAH  ") .OR. ANY(REACTANT(I,:).EQ."PAH0 ") .OR. &
          & ANY(REACTANT(I,:).EQ."PAH+ ") .OR. ANY(REACTANT(I,:).EQ."PAH- ")) THEN
            PHI_PAH=0.4D0
            RATE(I)=ALPHA(I)*(TEMPERATURE/100.0D0)**BETA(I)*PHI_PAH
            GOTO 10
         END IF

!C     Check for large negative gamma values that might cause discrepant
!C     rates at low temperatures. Set these rates to zero when T < RTMIN.

!CODE RESPONSINBLE FOR O + H+ --> O+ + H
         IF(DUPLICATE(I).EQ.0) THEN
            IF(GAMMA(I).LT.-200.0D0 .AND. TEMPERATURE.LT.RTMIN(I)) THEN
               RATE(I)=0.0D0
            ELSE
!if (i==76) then
!rate(i)=alpha(i)
!write(6,*) '1'
!else
               RATE(I)=ALPHA(I)*(TEMPERATURE/300.0D0)**BETA(I)*EXP(-(GAMMA(I)/TEMPERATURE))
!endif
            ENDIF
         ELSE IF(DUPLICATE(I).EQ.1) THEN
            J=I
            DO
               IF(TEMPERATURE.LE.RTMAX(J)) THEN
                  IF(GAMMA(J).LT.-200.0D0 .AND. TEMPERATURE.LT.RTMIN(J)) THEN
                  RATE(J)=0.0D0
                  ELSE
!if (i==76) then
!rate(j)=alpha(j)
!write(6,*) '2'
!else
                  RATE(J)=ALPHA(J)*(TEMPERATURE/300.0D0)**BETA(J)*EXP(-(GAMMA(J)/TEMPERATURE))
!endif
                  ENDIF
                  EXIT
               ELSE IF(DUPLICATE(J+1).LT.DUPLICATE(J)) THEN
                  IF(GAMMA(J).LT.-200.0D0 .AND. TEMPERATURE.LT.RTMIN(J)) THEN
                  RATE(J)=0.0D0
                  ELSE
!if (i==76) then
!rate(j)=alpha(j)
!write(6,*) '3'
!else
                  RATE(J)=ALPHA(J)*(TEMPERATURE/300.0D0)**BETA(J)*EXP(-(GAMMA(J)/TEMPERATURE))
!endif
                  ENDIF
                  EXIT
               ELSE
                  RATE(J)=0.0D0
                  J=J+1
               ENDIF
            ENDDO
         ENDIF
         GOTO 10

!C-----------------------------------------------------------------------

!C     Photoreactions:

!C     Store the reaction number for H2 photodissociation. The rate itself
!C     is calculated separately by the function H2PDRATE (within shield.f)
 1       IF(REACTANT(I,1).EQ."H2 " .AND. REACTANT(I,3).EQ."   ") THEN
!C           Loop over all rays
            DO K=0,NRAYS-1
               RATE(I)=RATE(I) + H2PDRATE(ALPHA(I),RAD_SURFACE(K),AV(K),COLUMN(K,NH2))
            ENDDO
            IF(PRODUCT(I,1).EQ."H " .AND. PRODUCT(I,2).EQ."H ") NRH2=I
            GOTO 10
         ENDIF

!C     Store the reaction number for HD photodissociation. The rate itself
!C     is calculated separately by the function H2PDRATE (within shield.f)
         IF(REACTANT(I,1).EQ."HD " .AND. REACTANT(I,3).EQ."   ") THEN
!C           Loop over all rays
            DO K=0,NRAYS-1
               RATE(I)=RATE(I) + H2PDRATE(ALPHA(I),RAD_SURFACE(K),AV(K),COLUMN(K,NHD))
            ENDDO
            IF(ANY(PRODUCT(I,:).EQ."H ") .AND. ANY(PRODUCT(I,:).EQ."D ")) NRHD=I
            GOTO 10
         ENDIF

!C     Store the reaction number for !CO photodissociation. The rate itself
!C     is calculated separately by the function !COPDRATE (within shield.f)
         IF(REACTANT(I,1).EQ."CO " .AND. REACTANT(I,3).EQ."   " .AND. &
       & ANY(PRODUCT(I,:).EQ."C ") .AND. ANY(PRODUCT(I,:).EQ."O ")) THEN
!C           Loop over all rays
            DO K=0,NRAYS-1
               RATE(I)=RATE(I) + COPDRATE(ALPHA(I),RAD_SURFACE(K),AV(K),COLUMN(K,NCO),COLUMN(K,NH2))
            ENDDO
            NRCO=I
            GOTO 10

         ENDIF

!C     Store the reaction number for !CI photoionization. The rate itself
!C     is calculated separately by the function CIPDRATE (within shield.f)
         IF(REACTANT(I,1).EQ."C  " .AND. REACTANT(I,3).EQ."   " .AND.&
      &    ((PRODUCT(I,1).EQ."C+ " .AND. PRODUCT(I,2).EQ."e- ") .OR.&
      &     (PRODUCT(I,1).EQ."e- " .AND. PRODUCT(I,2).EQ."C+ "))) THEN
!C           Loop over all rays
            DO K=0,NRAYS-1
               RATE(I)=RATE(I) + CIPDRATE(ALPHA(I),RAD_SURFACE(K),AV(K),GAMMA(I),COLUMN(K,NC),COLUMN(K,NH2),TEMPERATURE)
            ENDDO
            NRCI=I
            GOTO 10

         ENDIF

!C     Store the reaction number for SI photoionization. The rate itself
!C     is calculated separately by the function SIPDRATE (within shield.f)
         IF(REACTANT(I,1).EQ."S  " .AND. REACTANT(I,3).EQ."   " .AND.&
      &    ((PRODUCT(I,1).EQ."S+ " .AND. PRODUCT(I,2).EQ."e- ") .OR.&
      &     (PRODUCT(I,1).EQ."e- " .AND. PRODUCT(I,2).EQ."S+ "))) THEN
!C           Loop over all rays
            DO K=0,NRAYS-1
               RATE(I)=RATE(I) + SIPDRATE(ALPHA(I),RAD_SURFACE(K),AV(K),GAMMA(I),COLUMN(K,NS))
            ENDDO
            NRSI=I
            GOTO 10
         ENDIF

         IF(DUPLICATE(I).EQ.0) THEN
!C           Loop over all rays
            DO K=0,NRAYS-1
                RATE(I)=RATE(I) + ALPHA(I)*RAD_SURFACE(K)*EXP(-(GAMMA(I)*AV(K)))/2.0
            ENDDO
         ELSE IF(DUPLICATE(I).EQ.1) THEN
            J=I
            DO
               IF(TEMPERATURE.LE.RTMAX(J)) THEN
!C                 Loop over all rays
                  DO K=0,NRAYS-1
                     RATE(J)=RATE(J) + ALPHA(J)*RAD_SURFACE(K)*EXP(-(GAMMA(J)*AV(K)))/2.0
                  ENDDO
                  EXIT
               ELSE IF(DUPLICATE(J+1).LT.DUPLICATE(J)) THEN
!C                 Loop over all rays
                  DO K=0,NRAYS-1
                     RATE(J)=RATE(J) + ALPHA(J)*RAD_SURFACE(K)*EXP(-(GAMMA(J)*AV(K)))/2.0
                  ENDDO
                  EXIT
               ELSE
                  RATE(J)=0.0D0
                  J=J+1
               ENDIF
            ENDDO
         ENDIF
         GOTO 10

!C-----------------------------------------------------------------------

!C     Cosmic ray-induced ionization:

 2       IF(DUPLICATE(I).EQ.0) THEN
            RATE(I)=ALPHA(I)*ZETA
         ELSE IF(DUPLICATE(I).EQ.1) THEN
            J=I
            DO
               IF(TEMPERATURE.LE.RTMAX(J)) THEN
                  RATE(J)=ALPHA(J)*ZETA
                  EXIT
               ELSE IF(DUPLICATE(J+1).LT.DUPLICATE(J)) THEN
                  RATE(J)=ALPHA(J)*ZETA
                  EXIT
               ELSE
                  RATE(J)=0.0D0
                  J=J+1
               ENDIF
            ENDDO
         ENDIF
         GOTO 10

110     RATE(I)=0.0D0; GOTO 10
120     RATE(I)=0.0D0; GOTO 10
130     RATE(I)=0.0D0; GOTO 10
140     RATE(I)=0.0D0; GOTO 10

!C-----------------------------------------------------------------------

!C     Photoreactions due to cosmic ray-induced secondary photons:

 3       IF(DUPLICATE(I).EQ.0) THEN
            RATE(I)=ALPHA(I)*ZETA*(TEMPERATURE/300.0D0)**BETA(I)&
     &           *GAMMA(I)/(1.0D0-OMEGA)
         ELSE IF(DUPLICATE(I).EQ.1) THEN
            J=I
            DO
               IF(TEMPERATURE.LE.RTMAX(J)) THEN
                  RATE(J)=ALPHA(J)*ZETA*(TEMPERATURE/300.0D0)**BETA(J)&
     &                 *GAMMA(J)/(1.0D0-OMEGA)
                  EXIT
               ELSE IF(DUPLICATE(J+1).LT.DUPLICATE(J)) THEN
                  RATE(J)=ALPHA(J)*ZETA*(TEMPERATURE/300.0D0)**BETA(J)&
     &                 *GAMMA(J)/(1.0D0-OMEGA)
                  EXIT
               ELSE
                  RATE(J)=0.0D0
                  J=J+1
               ENDIF
            ENDDO
         ENDIF
         GOTO 10

!C-----------------------------------------------------------------------

!C     Freeze-out of neutral species:

 4       IF(BETA(I).EQ.0.0D0) THEN
            CION=1.0D0
         ELSE IF(BETA(I).EQ.1.0D0) THEN
            CION=1.0D0+16.71D-4/(GRAIN_RADIUS*TEMPERATURE)
         ELSE
            CION=0.0D0
         ENDIF
         STICKING=0.3D0
         RATE(I)=ALPHA(I)*4.57D4*2.4D-22*SQRT(TEMPERATURE/GAMMA(I))*CION*STICKING
         GOTO 10

!C-----------------------------------------------------------------------

!C     Freeze-out of singly charged positive ions:

 5       IF(BETA(I).EQ.0.0D0) THEN
            CION=1.0D0
         ELSE IF(BETA(I).EQ.1.0D0) THEN
            CION=1.0D0+16.71D-4/(GRAIN_RADIUS*TEMPERATURE)
         ELSE
            CION=0.0D0
         ENDIF
         STICKING=0.3D0
         RATE(I)=ALPHA(I)*4.57D4*2.4D-22*SQRT(TEMPERATURE/GAMMA(I))*CION*STICKING
         GOTO 10

!C-----------------------------------------------------------------------

!C     Desorption due to cosmic ray heating:

!CC     Treatment of Hasegawa & Herbst (1993, MNRAS, 261, 83, Equation 15)
!C 6       RATE(I)=ALPHA(I)*ZETA

!C     Treatment of Roberts et al. (2007, MNRAS, 382, 773, Equation 3)
 6       IF(GAMMA(I).LE.1210.0D0) THEN
            YIELD=1.0D5 ! Number of adsorbed molecules released per cosmic ray impact
         ELSE
            YIELD=0.0D0
         ENDIF
         FLUX=2.06D-3 ! Flux of iron nuclei cosmic rays (in cm^-2 s^-1)
         RATE(I)=FLUX*ZETA*2.4D-22*YIELD
         GOTO 10

!C-----------------------------------------------------------------------

!C     Photodesorption:

 7       IF(TEMPERATURE.LT.50.0D0) THEN
            YIELD=3.5D-3
         ELSE IF(TEMPERATURE.LT.85.0D0) THEN
            YIELD=4.0D-3
         ELSE IF(TEMPERATURE.LT.100.0D0) THEN
            YIELD=5.5D-3
         ELSE
            YIELD=7.5D-3
         ENDIF
!C         FLUX=1.0D8 ! Flux of FUV photons in the unattenuated Habing field (in photons cm^-2 s^-1)
         FLUX=1.7D8 ! Flux of FUV photons in the unattenuated Draine field (in photons cm^-2 s^-1)
!C        Loop over all rays
         DO K=0,NRAYS-1
            RATE(I)=RATE(I) + FLUX*RAD_SURFACE(K)*EXP(-(1.8D0*AV(K)))*2.4D-22*YIELD
         ENDDO
         GOTO 10

!C-----------------------------------------------------------------------

!C     Thermal desorption:

!C     Treatment of Hasegawa, Herbst & Leung (1992, ApJS, 82, 167, Equations 2 & 3)
 8       RATE(I)=SQRT(2.0D0*1.5D15*KB/(PI**2*AU)*ALPHA(I)/GAMMA(I))&
     &        *EXP(-(ALPHA(I)/DUST_TEMPERATURE))
         GOTO 10

!C-----------------------------------------------------------------------

!C     Grain mantle reactions:

 9       RATE(I)=ALPHA(I)
         GOTO 10

!C-----------------------------------------------------------------------

!C     Check that the rate is physical (0<RATE(I)<1) and produce an error
!C     message if not. Impose a lower cut-off on all rate coefficients to
!C     prevent the problem becoming too stiff. Rates less than 1E-99 are
!C     set to zero. Grain-surface reactions and desorption mechanisms are
!C     allowed rates greater than 1.
 10      IF(RATE(I).LT.0.0D0) THEN
           PRINT *,'ERROR! Negative rate for reaction',I
           STOP
         ENDIF
         IF(RATE(I).GT.1.0D0 .AND. REACTANT(I,1)(1:1).NE."G") THEN
           WRITE(10,*)'WARNING! Rate is too large for reaction',I
           WRITE(10,*)'RATE =',RATE(I)
           RATE(I)=1.0D0
         ENDIF
         IF(RATE(I).LT.1.0D-99) RATE(I)=0.0D0
!C     End of loop over rates
      ENDDO

      RETURN
      END SUBROUTINE CALCULATE_REACTION_RATES
