!=======================================================================
      SUBROUTINE ANALYSE_CHEMISTRY(GRIDPOINT,TIME,DENSITY,TEMPERATURE, &
     &                             NRAYS,AV,NSPEC,SPECIES,ABUNDANCE, &
     &                             NREAC,REACTANT,PRODUCT,RATE)

!T.Bell

      USE DEFINITIONS
      USE HEALPIX_TYPES
      USE GLOBAL_MODULE
      use maincode_module, only : alpha, beta

      IMPLICIT NONE

      INTEGER(KIND=I4B), INTENT(IN) :: GRIDPOINT,NRAYS,NSPEC,NREAC
      REAL(KIND=DP),     INTENT(IN) :: TIME,DENSITY,TEMPERATURE
      !REAL(KIND=DP),     INTENT(IN) :: AV(0:NRAYS-1),ABUNDANCE(1:NSPEC),RATE(1:NREAC)
      REAL(KIND=DP),     INTENT(IN) :: AV,ABUNDANCE(1:NSPEC),RATE(1:NREAC)
      CHARACTER(LEN=10), INTENT(IN) :: SPECIES(1:NSPEC),REACTANT(1:NREAC,1:3),PRODUCT(1:NREAC,1:4)

      INTEGER(KIND=I4B) :: I,J,L,M,N,IP,ID
      INTEGER(KIND=I4B) :: NPR(1:NREAC),NDR(1:NREAC)
      INTEGER(KIND=I4B) :: NMAX(1),TOTAL
      REAL(KIND=DP)     :: X1,X2,RMULT,PERCENT
      REAL(KIND=DP)     :: P(1:NREAC),PTOT
      REAL(KIND=DP)     :: D(1:NREAC),DTOT

!     Specify the number of species to examine and their names
      INTEGER(KIND=I4B), SAVE :: NLIST=7
      CHARACTER(LEN=10), SAVE :: SPECLIST(1:7)=&
     &      (/"H2        ","CO      ","C         ",&
     & "C+        ","CH        ","OH        ","e-        "/)
      integer::pl

      WRITE(98,6)
      DO L=1,NLIST
!        Find index I of the species that corresponds
!        to the current entry in the SPECLIST array
         DO I=1,NSPEC
            IF(SPECIES(I).EQ.SPECLIST(L)) EXIT
         ENDDO

!        Check that species I corresponds to the desired species
!        Continue to the next entry in SPECLIST if it does not
         IF(SPECIES(I).NE.SPECLIST(L)) CYCLE

         !WRITE(98,4) GRIDPOINT,TIME,MINVAL(AV),DENSITY,TEMPERATURE,SPECIES(I)
         WRITE(98,4) GRIDPOINT,TIME,AV,DENSITY,TEMPERATURE,SPECIES(I)

!        Reset the formation/destruction rate counters
         IP=0
         ID=0

!        Reset the total formation/destruction rates
         PTOT=0.0D0
         DTOT=0.0D0

!        Check all reactions to find relevant ones
         DO J=1,NREAC
!           Reset the reactant abundances
            X1=0.0D0
            X2=0.0D0

!-----------------------------------------------------------------------
!           Formation routes for species I
!-----------------------------------------------------------------------
            IF(PRODUCT(J,1).EQ.SPECIES(I) .OR. PRODUCT(J,2).EQ.SPECIES(I) .OR. &
     &         PRODUCT(J,3).EQ.SPECIES(I) .OR. PRODUCT(J,4).EQ.SPECIES(I)) THEN

               IP=IP+1

!              Store the reactant abundances
               DO N=1,NSPEC
                  IF(SPECIES(N).EQ.REACTANT(J,1)) X1=ABUNDANCE(N)
                  IF(SPECIES(N).EQ.REACTANT(J,2)) X2=ABUNDANCE(N)
               ENDDO

               RMULT=0.0D0
               IF(PRODUCT(J,1).EQ.SPECIES(I)) RMULT=RMULT+1.0D0
               IF(PRODUCT(J,2).EQ.SPECIES(I)) RMULT=RMULT+1.0D0
               IF(PRODUCT(J,3).EQ.SPECIES(I)) RMULT=RMULT+1.0D0
               IF(PRODUCT(J,4).EQ.SPECIES(I)) RMULT=RMULT+1.0D0

!              Calculate the reaction rate
               IF(REACTANT(J,2).EQ."PHOTON" .OR. &
     &            REACTANT(J,2).EQ."CRP   " .OR. &
     &            REACTANT(J,2).EQ."CRPHOT" .OR. &
     &            REACTANT(J,2).EQ."FREEZE" .OR. &
     &            REACTANT(J,2).EQ."ELFRZE" .OR. &
     &            REACTANT(J,2).EQ."CRH   " .OR. &
     &            REACTANT(J,2).EQ."PHOTD " .OR. &
     &            REACTANT(J,2).EQ."THERM ") THEN
                  P(IP)=RATE(J)*X1*RMULT
               ELSE
                  P(IP)=RATE(J)*X1*X2*DENSITY*RMULT
               ENDIF

               PTOT=PTOT+P(IP)
               NPR(IP)=J
            ENDIF

!-----------------------------------------------------------------------
!           Destruction routes for species I
!-----------------------------------------------------------------------

            IF(REACTANT(J,1).EQ.SPECIES(I) .OR. REACTANT(J,2).EQ.SPECIES(I)) THEN

               ID=ID+1

!              Store the reactant abundances
               DO N=1,NSPEC
                  IF(SPECIES(N).EQ.REACTANT(J,1)) X1=ABUNDANCE(N)
                  IF(SPECIES(N).EQ.REACTANT(J,2)) X2=ABUNDANCE(N)
               ENDDO

               RMULT=0.0D0
               IF(REACTANT(J,1).EQ.SPECIES(I)) RMULT=RMULT+1.0D0
               IF(REACTANT(J,2).EQ.SPECIES(I)) RMULT=RMULT+1.0D0

!              Calculate the reaction rate
               IF(REACTANT(J,2).EQ."PHOTON" .OR. &
     &            REACTANT(J,2).EQ."CRP   " .OR. &
     &            REACTANT(J,2).EQ."CRPHOT" .OR. &
     &            REACTANT(J,2).EQ."FREEZE" .OR. &
     &            REACTANT(J,2).EQ."ELFRZE" .OR. &
     &            REACTANT(J,2).EQ."CRH   " .OR. &
     &            REACTANT(J,2).EQ."PHOTD " .OR. &
     &            REACTANT(J,2).EQ."THERM ") THEN
                  D(ID)=RATE(J)*X1*RMULT
               ELSE
                  D(ID)=RATE(J)*X1*X2*DENSITY*RMULT
               ENDIF

               DTOT=DTOT+D(ID)
               NDR(ID)=J
            ENDIF

!        End of loop over all reactions
         ENDDO

!        Prevent divide-by-zero errors by setting minimum finite
!        values for the total formation and destruction rates
         IF(PTOT.LT.1.0D-99) PTOT=1.0D-99
         IF(DTOT.LT.1.0D-99) DTOT=1.0D-99

!-----------------------------------------------------------------------
!        Output the formation reactions and their rates
!-----------------------------------------------------------------------

         NMAX=0
         TOTAL=0

!        List the formation reactions in order of decreasing importance
         DO M=1,IP
!           Find the location of the maximum value
            NMAX=MAXLOC(P(1:IP))
            N=NMAX(1)
!           Exit the loop once the reaction rates reach zero
            IF(P(N).EQ.0.0D0) EXIT

!           Calculate the percentage of the total formation
!           rate that is contributed by the current reaction
            PERCENT=1.0D2*(P(N)/PTOT)
            IF(PERCENT.LT.0.5D0) PERCENT=1.0D0
            TOTAL=TOTAL+NINT(PERCENT)


            IF(PRODUCT(NPR(N),3).EQ." ") THEN
               WRITE(98,1) (REACTANT(NPR(N),J),J=1,2),(PRODUCT(NPR(N),J),J=1,2),NINT(PERCENT)
            ELSE IF(PRODUCT(NPR(N),4).EQ." ") THEN
               WRITE(98,2) (REACTANT(NPR(N),J),J=1,2),(PRODUCT(NPR(N),J),J=1,3),NINT(PERCENT)
            ELSE
               WRITE(98,3) (REACTANT(NPR(N),J),J=1,2),(PRODUCT(NPR(N),J),J=1,4),NINT(PERCENT)
            ENDIF

!           Exit the loop once the sum of the reaction rates reaches 100%
            IF(TOTAL.GE.100) EXIT

!           Set the formation rate of this reaction to zero
!           to prevent it from being included more than once
            P(N)=0.0D0

!        End of loop over formation reactions
         ENDDO
         WRITE(98,*)

!-----------------------------------------------------------------------
!        Output the destruction reactions and their rates
!-----------------------------------------------------------------------

         NMAX=0
         TOTAL=0

!        List the destruction reactions in order of decreasing importance
         DO M=1,ID
!           Find the location of the maximum value
            NMAX=MAXLOC(D(1:ID))
            N=NMAX(1)

!           Exit the loop once the reaction rates reach zero
            IF(D(N).EQ.0.0D0) EXIT

!           Calculate the percentage of the total destruction
!           rate that is contributed by the current reaction
            PERCENT=1.0D2*(D(N)/DTOT)
            IF(PERCENT.LT.0.5D0) PERCENT=1.0D0
            TOTAL=TOTAL+NINT(PERCENT)

            IF(PRODUCT(NDR(N),3).EQ." ") THEN
               WRITE(98,1) (REACTANT(NDR(N),J),J=1,2),(PRODUCT(NDR(N),J),J=1,2),-NINT(PERCENT)
            ELSE IF(PRODUCT(NDR(N),4).EQ." ") THEN
               WRITE(98,2) (REACTANT(NDR(N),J),J=1,2),(PRODUCT(NDR(N),J),J=1,3),-NINT(PERCENT)
            ELSE
               WRITE(98,3) (REACTANT(NDR(N),J),J=1,2),(PRODUCT(NDR(N),J),J=1,4),-NINT(PERCENT)
            ENDIF

!           Exit the loop once the sum of the reaction rates reaches 100%
            IF(TOTAL.GE.100) EXIT

!           Set the destruction rate of this reaction to zero
!           to prevent it from being included more than once
            D(N)=0.0D0

!        End of loop over destruction reactions
         ENDDO
         WRITE(98,*)

!-----------------------------------------------------------------------

!        Output the abundance and total formation/destruction rates
         WRITE(98,5) ABUNDANCE(I),PTOT,DTOT
         WRITE(98,6)

!        End of loop over SPECLIST array
      ENDDO
      WRITE(98,7)

 1    FORMAT(3X,A10,'+',3X,A10,'-->',3X,A10,'+',3X,A10,12X,'Rate:',I4,'%')
 2    FORMAT(3X,A10,'+',3X,A10,'-->',3X,A10,'+',3X,A10,'+',3X,A5,3X,'Rate:',I4,'%')
 3    FORMAT(3X,A10,'+',3X,A10,'-->',3X,A10,'+',3X,A10,'+',3X,A5,3X,'+',3X,A5,3X,'Rate:',I4,'%')
 4    FORMAT('Gridpoint =',I7,', t =',1PE7.1E1,' yr, Av =',1PE7.1E1,' mag',/, &
     &       'n_H =',1PE7.1E1,' cm-3, T_gas =',1PE7.1E1,' K',/, &
     &       'Species = ',A10,/)
 5    FORMAT('Abundance        =',1PE10.3,/, &
     &       'Formation Rate   =',1PE10.3,' s-1',/, &
     &       'Destruction Rate =',1PE10.3,' s-1',/)
 6    FORMAT(80('-'),/)
 7    FORMAT(80('='),/)

      RETURN
      END SUBROUTINE ANALYSE_CHEMISTRY
!=======================================================================
