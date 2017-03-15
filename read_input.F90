!-----------------------------------------------------------------------
!     Read atomic/molecular datafile (in LAMDA/RADEX format)
!     The standard LAMDA/RADEX datafile format is used:
!     http://www.strw.leidenuniv.nl/~moldata/molformat.html
!-----------------------------------------------------------------------
      SUBROUTINE READINPUT(FILENAME,NLEV,NTEMP,ENERGIES,WEIGHTS,&
     &                     A_COEFFS,B_COEFFS,FREQUENCIES,TEMPERATURES,&
     &                     H_COL,HP_COL,EL_COL,HE_COL,H2_COL,PH2_COL,OH2_COL)

!T.Bell

  use healpix_types

      implicit none
      INTEGER(kind=I4B), INTENT(IN)::NLEV
      INTEGER(kind=I4B), INTENT(IN)::NTEMP
      CHARACTER(len=*),INTENT(IN)::FILENAME

      real(kind=dp), intent(out)::ENERGIES(1:NLEV), WEIGHTS(1:NLEV)
      real(kind=dp), intent(out)::A_COEFFS(1:NLEV,1:NLEV), B_COEFFS(1:NLEV,1:NLEV)
      real(kind=dp), intent(out)::FREQUENCIES(1:NLEV,1:NLEV), TEMPERATURES(1:7,1:NTEMP)
      real(kind=dp), intent(out)::H_COL(1:NLEV,1:NLEV,1:NTEMP)
      real(kind=dp), intent(out)::HP_COL(1:NLEV,1:NLEV,1:NTEMP)
      real(kind=dp), intent(out)::EL_COL(1:NLEV,1:NLEV,1:NTEMP)
      real(kind=dp), intent(out)::HE_COL(1:NLEV,1:NLEV,1:NTEMP)
      real(kind=dp), intent(out)::H2_COL(1:NLEV,1:NLEV,1:NTEMP)
      real(kind=dp), intent(out)::PH2_COL(1:NLEV,1:NLEV,1:NTEMP)
      real(kind=dp), intent(out)::OH2_COL(1:NLEV,1:NLEV,1:NTEMP)

      INTEGER(kind=I4B)::NLIN,NPART,NCOL
      INTEGER(kind=I4B)::I,J,K,L,M,N,P
      real(kind=dp):: ENERGY,WEIGHT,EINSTEINA,FREQUENCY
      real(kind=dp)::COEFF(1:NTEMP)


!     Initialize all variables to zero before reading in the data
      DO I=1,NLEV
         ENERGIES(I)=0.0D0
         WEIGHTS(I)=0.0D0
         DO J=1,NLEV
            A_COEFFS(I,J)=0.0D0
            B_COEFFS(I,J)=0.0D0
            FREQUENCIES(I,J)=0.0D0
            DO K=1,NTEMP
               TEMPERATURES(:,K)=0.0D0
               H_COL(I,J,K)=0.0D0
               EL_COL(I,J,K)=0.0D0
               HE_COL(I,J,K)=0.0D0
               H2_COL(I,J,K)=0.0D0
               PH2_COL(I,J,K)=0.0D0
               OH2_COL(I,J,K)=0.0D0
            ENDDO
         ENDDO
      ENDDO

      OPEN(8,FILE=FILENAME,STATUS='OLD')
      READ(8,'(////)') !empty line
      READ(8,*) N !number of levels of the file
      IF(N.NE.NLEV) STOP "ERROR! Incorrect number of energy levels, N>NLEV"
      READ(8,*)
      N=0
      DO WHILE(N.LT.NLEV)
         READ(8,*) N,ENERGY,WEIGHT
         ENERGIES(N)=ENERGY*C*HP ! Convert from cm^-1 to erg
         WEIGHTS(N)=WEIGHT
      ENDDO
      READ(8,*)
      READ(8,*) NLIN
      READ(8,*)
      N=0
      DO WHILE(N.LT.NLIN)
         READ(8,*) N,I,J,EINSTEINA,FREQUENCY
         FREQUENCIES(I,J)=FREQUENCY*1.0D9 ! Convert from GHz to Hz
         FREQUENCIES(J,I)=FREQUENCIES(I,J)
         A_COEFFS(I,J)=EINSTEINA
!        Calculate the Einstein B coefficients using Bij = Aij/(2.h.nu^3/c^2)
         B_COEFFS(I,J)=A_COEFFS(I,J)&
     &                 /(2.0D0*HP*(FREQUENCIES(I,J)**3)/(C**2))
         B_COEFFS(J,I)=B_COEFFS(I,J)*(WEIGHTS(I)/WEIGHTS(J))
      ENDDO
!     Calculate the transition frequencies between all levels (even if forbidden)
      DO I=1,NLEV
         DO J=1,NLEV
            FREQUENCY=ABS(ENERGIES(I)-ENERGIES(J))/HP
            IF(FREQUENCIES(I,J).NE.0.0D0) THEN
!              Check if the calculated and measured frequencies differ by >1%
               IF(ABS(FREQUENCY-FREQUENCIES(I,J))&
     &                /FREQUENCIES(I,J).GT.1.0D-2) THEN
                  WRITE(6,*) 'ERROR! Calculated frequency differs by >1%:'
                  WRITE(6,*) FREQUENCY,' Hz vs',FREQUENCIES(I,J),' Hz'
                  STOP
               ENDIF
            ELSE
               FREQUENCIES(I,J)=FREQUENCY
            ENDIF
         ENDDO
      ENDDO
      READ(8,*)
      READ(8,*) NPART
!     Read the collisional rate coefficients (cm^3 s^-1) for each collision partner
      DO L=1,NPART
         READ(8,*)
         READ(8,*) P
         READ(8,*)
         READ(8,*) NCOL
         READ(8,*)
         READ(8,*) M
         IF(M.GT.NTEMP) THEN
            WRITE(6,*) 'ERROR! Too many temperature values (>NTEMP):',M
            STOP
         ENDIF
         IF(P.EQ.1) THEN
            READ(8,*)
            READ(8,*) (TEMPERATURES(P,K),K=1,M)
            READ(8,*)
            N=0
            DO WHILE(N.LT.NCOL)
               READ(8,*) N,I,J,(COEFF(K),K=1,M)
               DO K=1,M
                  H2_COL(I,J,K)=COEFF(K)
!                 Calculate the reverse (excitation) rate coefficient
!                 from detailed balance: Cji = Cij*gi/gj*exp(-(Ei-Ej)/kT)
                  IF(H2_COL(I,J,K).NE.0.0D0 .AND. H2_COL(J,I,K).EQ.0.0D0) THEN
                     H2_COL(J,I,K)=H2_COL(I,J,K)*(WEIGHTS(I)/WEIGHTS(J)) &
     &              *EXP(-(ENERGIES(I)-ENERGIES(J))/(KB*TEMPERATURES(P,K)))
                  ENDIF
               ENDDO
            ENDDO
         ELSE IF(P.EQ.2) THEN
            READ(8,*)
            READ(8,*) (TEMPERATURES(P,K),K=1,M)
            READ(8,*)
            N=0
            DO WHILE(N.LT.NCOL)
               READ(8,*) N,I,J,(COEFF(K),K=1,M)
               DO K=1,M
                  PH2_COL(I,J,K)=COEFF(K)
!                 Calculate the reverse (excitation) rate coefficient
!                 from detailed balance: Cji = Cij*gi/gj*exp(-(Ei-Ej)/kT)
                  IF(PH2_COL(I,J,K).NE.0.0D0 .AND. PH2_COL(J,I,K).EQ.0.0D0) THEN
                     PH2_COL(J,I,K)=PH2_COL(I,J,K)*(WEIGHTS(I)/WEIGHTS(J)) &
     &              *EXP(-(ENERGIES(I)-ENERGIES(J))/(KB*TEMPERATURES(P,K)))
                  ENDIF
               ENDDO
            ENDDO
!         ELSE IF(P.EQ.3) THEN
!            READ(8,*)
!            READ(8,*) (TEMPERATURES(P,K),K=1,M)
!            READ(8,*)
!            N=0
!            DO WHILE(N.LT.NCOL)
!               READ(8,*) N,I,J,(COEFF(K),K=1,M)
!               DO K=1,M
!                  OH2_COL(I,J,K)=COEFF(K)   
!!                 Calculate the reverse (excitation) rate coefficient
!!                 from detailed balance: Cji = Cij*gi/gj*exp(-(Ei-Ej)/kT)
!                  IF(PH2_COL(I,J,K).NE.0.0D0 .AND. PH2_COL(J,I,K).EQ.0.0D0) THEN
!                     PH2_COL(J,I,K)=PH2_COL(I,J,K)*(WEIGHTS(I)/WEIGHTS(J)) &
!     &              *EXP(-(ENERGIES(I)-ENERGIES(J))/(KB*TEMPERATURES(P,K)))
!                  ENDIF
!               ENDDO
!            ENDDO
         ELSE IF(P.EQ.3) THEN
            READ(8,*)
            READ(8,*) (TEMPERATURES(P,K),K=1,M)
            READ(8,*)
            N=0
            DO WHILE(N.LT.NCOL)
               READ(8,*) N,I,J,(COEFF(K),K=1,M)
               DO K=1,M
                  OH2_COL(I,J,K)=COEFF(K)
!                 Calculate the reverse (excitation) rate coefficient
!                 from detailed balance: Cji = Cij*gi/gj*exp(-(Ei-Ej)/kT)
                  IF(OH2_COL(I,J,K).NE.0.0D0 .AND. OH2_COL(J,I,K).EQ.0.0D0) THEN
                     OH2_COL(J,I,K)=OH2_COL(I,J,K)*(WEIGHTS(I)/WEIGHTS(J)) &
     &              *EXP(-(ENERGIES(I)-ENERGIES(J))/(KB*TEMPERATURES(P,K)))
                  ENDIF
               ENDDO
            ENDDO
         ELSE IF(P.EQ.4) THEN
            READ(8,*)
            READ(8,*) (TEMPERATURES(P,K),K=1,M)
            READ(8,*)
            N=0
            DO WHILE(N.LT.NCOL)
               READ(8,*) N,I,J,(COEFF(K),K=1,M)
               DO K=1,M
                  EL_COL(I,J,K)=COEFF(K)
!                 Calculate the reverse (excitation) rate coefficient
!                 from detailed balance: Cji = Cij*gi/gj*exp(-(Ei-Ej)/kT)
                  IF(EL_COL(I,J,K).NE.0.0D0 .AND. EL_COL(J,I,K).EQ.0.0D0) THEN
                     EL_COL(J,I,K)=EL_COL(I,J,K)*(WEIGHTS(I)/WEIGHTS(J)) &
     &              *EXP(-(ENERGIES(I)-ENERGIES(J))/(KB*TEMPERATURES(P,K)))
                  ENDIF
               ENDDO
            ENDDO
         ELSE IF(P.EQ.5) THEN
            READ(8,*)
            READ(8,*) (TEMPERATURES(P,K),K=1,M)
            READ(8,*)
            N=0
            DO WHILE(N.LT.NCOL)
               READ(8,*) N,I,J,(COEFF(K),K=1,M)
               DO K=1,M
                  H_COL(I,J,K)=COEFF(K)
!                 Calculate the reverse (excitation) rate coefficient
!                 from detailed balance: Cji = Cij*gi/gj*exp(-(Ei-Ej)/kT)
                  IF(H_COL(I,J,K).NE.0.0D0 .AND. H_COL(J,I,K).EQ.0.0D0) THEN
                     H_COL(J,I,K)=H_COL(I,J,K)*(WEIGHTS(I)/WEIGHTS(J)) &
     &              *EXP(-(ENERGIES(I)-ENERGIES(J))/(KB*TEMPERATURES(P,K)))
                  ENDIF
               ENDDO
            ENDDO
         ELSE IF(P.EQ.6) THEN
            READ(8,*)
            READ(8,*) (TEMPERATURES(P,K),K=1,M)
            READ(8,*)
            N=0
            DO WHILE(N.LT.NCOL)
               READ(8,*) N,I,J,(COEFF(K),K=1,M)
               DO K=1,M
                  HE_COL(I,J,K)=COEFF(K)
!                 Calculate the reverse (excitation) rate coefficient
!                 from detailed balance: Cji = Cij*gi/gj*exp(-(Ei-Ej)/kT)
                  IF(HE_COL(I,J,K).NE.0.0D0 .AND. HE_COL(J,I,K).EQ.0.0D0) THEN
                     HE_COL(J,I,K)=HE_COL(I,J,K)*(WEIGHTS(I)/WEIGHTS(J)) &
     &              *EXP(-(ENERGIES(I)-ENERGIES(J))/(KB*TEMPERATURES(P,K)))
                  ENDIF
               ENDDO
            ENDDO
         ELSE IF(P.EQ.7) THEN
            READ(8,*)
            READ(8,*) (TEMPERATURES(P,K),K=1,M)
            READ(8,*)
            N=0
            DO WHILE(N.LT.NCOL)
               READ(8,*) N,I,J,(COEFF(K),K=1,M)
               DO K=1,M
                  HP_COL(I,J,K)=COEFF(K)
!                 Calculate the reverse (excitation) rate coefficient
!                 from detailed balance: Cji = Cij*gi/gj*exp(-(Ei-Ej)/kT)
                  IF(HP_COL(I,J,K).NE.0.0D0 .AND. HP_COL(J,I,K).EQ.0.0D0) THEN
                     HP_COL(J,I,K)=HP_COL(I,J,K)*(WEIGHTS(I)/WEIGHTS(J)) &
     &              *EXP(-(ENERGIES(I)-ENERGIES(J))/(KB*TEMPERATURES(P,K)))
                  ENDIF
               ENDDO
            ENDDO
         ELSE
            WRITE(6,*) 'ERROR! Unrecognized collision partner ID:',P
            STOP
         ENDIF
      ENDDO

      CLOSE(8)
      WRITE(6,*) 'Cooling datafile: ',FILENAME,' read successfully'
      RETURN

      END subroutine

