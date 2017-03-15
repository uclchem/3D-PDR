SUBROUTINE find_Ccoeff(NTEMP,NLEV,TEMPERATURE,TEMPERATURES,H_COL,HP_COL,EL_COL,HE_COL, &
                     & H2_COL,PH2_COL,OH2_COL,C_COEFFS,H_abd,Hp_abd,elec_abd,He_abd,H2_abd,&
                     & coolant)!, step_check, coef)
!T.Bell

 use definitions
 use healpix_types
 use global_module
 implicit none

 integer(kind=i4b), intent(in) :: NTEMP, NLEV, coolant
 real(kind=dp), intent(in)::TEMPERATURE
 real(kind=dp), intent(in)::TEMPERATURES(1:7,1:NTEMP)
 real(kind=dp), intent(in)::H_COL(1:NLEV,1:NLEV,1:NTEMP)
 real(kind=dp), intent(in)::HP_COL(1:NLEV,1:NLEV,1:NTEMP)
 real(kind=dp), intent(in)::EL_COL(1:NLEV,1:NLEV,1:NTEMP)
 real(kind=dp), intent(in)::HE_COL(1:NLEV,1:NLEV,1:NTEMP)
 real(kind=dp), intent(in)::H2_COL(1:NLEV,1:NLEV,1:NTEMP)
 real(kind=dp), intent(in)::PH2_COL(1:NLEV,1:NLEV,1:NTEMP)
 real(kind=dp), intent(in)::OH2_COL(1:NLEV,1:NLEV,1:NTEMP)

 real(kind=dp), intent(in)::H_abd, Hp_abd, elec_abd, He_abd, H2_abd

 real(kind=dp), intent(out)::C_COEFFS(1:NLEV,1:NLEV)
! real(kind=dp), intent(out)::coef(1:nlev,1:nlev,1:7), step_check(1:7)

 integer(kind=i4b) :: i,j,k,klo,khi,partner_id
 real(kind=dp) :: step, tmp
 real(kind=dp) :: FPARA, FORTHO

!  Initialize the collisional rates
   C_COEFFS=0.0D0

!  Calculate the H2 ortho/para ratio at equilibrium for the specified
!  temperature and the resulting fractions of H2 in para & ortho form
   FPARA=0.0D0 ; FORTHO=0.0D0
   IF(H2_abd.GT.0.0D0) THEN
      FPARA=1.0D0/(1.0D0+9.0D0*EXP(-170.5D0/TEMPERATURE))
      FORTHO=1.0D0-FPARA
   ENDIF

   DO PARTNER_ID=1,7 ! Loop over collision partners

!     Skip the collision partner if no rates are available
      IF(TEMPERATURES(PARTNER_ID,1).EQ.0.0D0) CYCLE

!     Determine the two nearest temperature values
!     present within the list of collisional rates
      KLO=0; KHI=0
      DO K=1,NTEMP ! Loop over temperatures
         IF(TEMPERATURES(PARTNER_ID,K).GT.TEMPERATURE) THEN
            KLO=K-1
            KHI=K
            EXIT
         ELSE IF(TEMPERATURES(PARTNER_ID,K).EQ.0.0D0) THEN
            KLO=K-1
            KHI=K-1
            EXIT
         END IF
      END DO

!     If the required temperature is above or below the range of available
!     temperature values then use the highest or lowest value in the range
      IF(KHI.EQ.0) THEN
         KLO=NTEMP
         KHI=NTEMP
      ELSE IF(KHI.EQ.1) THEN
         KLO=1
         KHI=1
      END IF

!     Calculate the "distance" between the two temperature
!     values, to be used in the linear interpolation below
      IF(KLO.EQ.KHI) THEN
         STEP=0.0D0
      ELSE
         STEP=(TEMPERATURE-TEMPERATURES(PARTNER_ID,KLO)) &
           & /(TEMPERATURES(PARTNER_ID,KHI)-TEMPERATURES(PARTNER_ID,KLO))
      END IF
!if (coolant.eq.1) step_check(partner_id)=step

!     Linearly interpolate the collisional rate coefficients
!     for each collision partner at the required temperature
      IF(PARTNER_ID.EQ.1) THEN ! Collisions with H2
         DO I=1,NLEV
            DO J=1,NLEV
               TMP=H2_COL(I,J,KLO)+(H2_COL(I,J,KHI)-H2_COL(I,J,KLO))*STEP
               C_COEFFS(I,J)=C_COEFFS(I,J)+TMP*H2_abd!*DENSITY*ABUNDANCE(nH2)
!if (coolant.eq.1) coef(i,j,partner_id) = TMP*H2_abd
            END DO
         END DO
      END IF
      IF(PARTNER_ID.EQ.2) THEN ! Collisions with para-H2
         DO I=1,NLEV
            DO J=1,NLEV
               TMP=PH2_COL(I,J,KLO)+(PH2_COL(I,J,KHI)-PH2_COL(I,J,KLO))*STEP
               C_COEFFS(I,J)=C_COEFFS(I,J)+TMP*H2_abd*FPARA!*DENSITY*ABUNDANCE(nH2)*PARA_FRACTION
!if (coolant.eq.1) coef(i,j,partner_id) = TMP*H2_abd*FPARA
            END DO
         END DO
      END IF
      IF(PARTNER_ID.EQ.3) THEN ! Collisions with ortho-H2
         DO I=1,NLEV
            DO J=1,NLEV
               TMP=OH2_COL(I,J,KLO)+(OH2_COL(I,J,KHI)-OH2_COL(I,J,KLO))*STEP
               C_COEFFS(I,J)=C_COEFFS(I,J)+TMP*H2_abd*FORTHO!*DENSITY*ABUNDANCE(nH2)*ORTHO_FRACTION
!if (coolant.eq.1) coef(i,j,partner_id) = TMP*H2_abd*FORTHO
            END DO
         END DO
      END IF
      IF(PARTNER_ID.EQ.4) THEN ! Collisions with electrons
         DO I=1,NLEV
            DO J=1,NLEV
               TMP=EL_COL(I,J,KLO)+(EL_COL(I,J,KHI)-EL_COL(I,J,KLO))*STEP
               C_COEFFS(I,J)=C_COEFFS(I,J)+TMP*elec_abd!*DENSITY*ABUNDANCE(nelect)
!if (coolant.eq.1) coef(i,j,partner_id) = TMP*elec_abd
            END DO
         END DO
      END IF
      IF(PARTNER_ID.EQ.5) THEN ! Collisions with H
         DO I=1,NLEV
            DO J=1,NLEV
               TMP=H_COL(I,J,KLO)+(H_COL(I,J,KHI)-H_COL(I,J,KLO))*STEP
               C_COEFFS(I,J)=C_COEFFS(I,J)+TMP*H_abd!*DENSITY*ABUNDANCE(nH)
!if (coolant.eq.1) coef(i,j,partner_id) = TMP*H_abd
            END DO
         END DO
      END IF
      IF(PARTNER_ID.EQ.6) THEN ! Collisions with He
         DO I=1,NLEV
            DO J=1,NLEV
               TMP=HE_COL(I,J,KLO)+(HE_COL(I,J,KHI)-HE_COL(I,J,KLO))*STEP
               C_COEFFS(I,J)=C_COEFFS(I,J)+TMP*He_abd!*DENSITY*ABUNDANCE(nHe)
!if (coolant.eq.1) coef(i,j,partner_id) = TMP*He_abd
            END DO
         END DO
      END IF
      IF(PARTNER_ID.EQ.7) THEN ! Collisions with protons
         DO I=1,NLEV
            DO J=1,NLEV
               TMP=HP_COL(I,J,KLO)+(HP_COL(I,J,KHI)-HP_COL(I,J,KLO))*STEP
               C_COEFFS(I,J)=C_COEFFS(I,J)+TMP*Hp_abd!*DENSITY*ABUNDANCE(nHx)
!if (coolant.eq.1) coef(i,j,partner_id) = TMP*Hp_abd
            END DO
         END DO
      END IF

   END DO ! End of loop over collision partners
   RETURN
END SUBROUTINE
