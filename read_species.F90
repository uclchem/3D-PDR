!C***********************************************************************
!C     Read in the chemical reaction rates and the species, masses and
!C     initial abundances (if specified). The  rates and species files
!C     are assumed to have comma separated values (CSV) format. This is
!C     in line with the Rate05 formatting, removing the need for file-
!C     dependent FORMAT statements.
!C***********************************************************************
!
!C-----------------------------------------------------------------------
!C     Read in the species data, including initial fractional abundances
!C     (3rd column) and their masses (4th column). Check that the value
!C     of NSPEC agrees with the number of species in the file and produce
!C     an error message if not.
!C-----------------------------------------------------------------------
      SUBROUTINE READ_SPECIES(NSPEC,SPECIES,ABUNDANCE,MASS)

!T.Bell
use definitions
use healpix_types
use global_module

      IMPLICIT NONE
      INTEGER(kind=i4b),intent(in) :: NSPEC
      real(kind=dp), intent(out) :: ABUNDANCE(1:nspec),MASS(1:nspec)
      CHARACTER(len=10), intent(out) :: SPECIES(1:nspec)

      INTEGER(kind=i4b) :: I,INDEX,SPECIESFILE

      SPECIESFILE = 3

!C     Initialize the variables and read in the species data. Check that
!C     the value of NSPEC agrees with the number of species in the file
!C     and produce an error message if not.
      SPECIES="          "
      ABUNDANCE=0.0D0
      MASS=0.0D0

!C     Initialize all the species index labels. If they are not assigned
!C     subsequently, any attempt to access that species will generate an
!C     error and the code will crash. This is a useful bug catch.
      NH=0
      ND=0
      NH2=0
      NHD=0
      NH2x=0
      NPROTON=0
      NC=0
      NCx=0
      NO=0
      NOx=0
      NN=0
      NNx=0
      NS=0
      NSx=0
      NHE=0
      NHEx=0
      NNA=0
      NNAx=0
      NMG=0
      NMGx=0
      NSI=0
      NSIx=0
      NFE=0
      NFEx=0
      NCL=0
      NCLx=0
      NCA=0
      NCAx=0
      NCAxx=0
      NCO=0
      NCH=0
      NCH2=0
      NOH=0
      NO2=0
      NCS=0
      NH2O=0
      NELECT=0
      NH3x=0
      NH3Ox=0
      NHCOx=0
#ifdef REDUCED
      OPEN(SPECIESFILE,FILE="species_reduced.d",STATUS="OLD")
#endif
#ifdef FULL
      OPEN(SPECIESFILE,FILE="species_full.d",STATUS="OLD")
#endif
#ifdef MYNETWORK
      OPEN(SPECIESFILE,FILE="species_mynetwork.d",STATUS="OLD")
#endif
      REWIND(SPECIESFILE)
      DO I=1,NSPEC
         READ(SPECIESFILE,*,END=1) INDEX,SPECIES(I),ABUNDANCE(I),MASS(I)

!C        Assign the various index labels to their correct species.
         IF(SPECIES(I).EQ."H         ") NH      = I
         IF(SPECIES(I).EQ."D         ") ND      = I
         IF(SPECIES(I).EQ."H2        ") NH2     = I
         IF(SPECIES(I).EQ."HD        ") NHD     = I
         IF(SPECIES(I).EQ."H2+       ") NH2x    = I
         IF(SPECIES(I).EQ."H3+       ") NH3x    = I
         IF(SPECIES(I).EQ."H+        ") NPROTON = I
         IF(SPECIES(I).EQ."C         ") NC      = I
         IF(SPECIES(I).EQ."C+        ") NCx     = I
         IF(SPECIES(I).EQ."O         ") NO      = I
         IF(SPECIES(I).EQ."O+        ") NOx     = I
         IF(SPECIES(I).EQ."N         ") NN      = I
         IF(SPECIES(I).EQ."N+        ") NNx     = I
         IF(SPECIES(I).EQ."S         ") NS      = I
         IF(SPECIES(I).EQ."S+        ") NSx     = I
         IF(SPECIES(I).EQ."He        ") NHE     = I
         IF(SPECIES(I).EQ."HE        ") NHE     = I
         IF(SPECIES(I).EQ."He+       ") NHEx    = I
         IF(SPECIES(I).EQ."HE+       ") NHEx    = I
         IF(SPECIES(I).EQ."Na        ") NNA     = I
         IF(SPECIES(I).EQ."NA        ") NNA     = I
         IF(SPECIES(I).EQ."Na+       ") NNAx    = I
         IF(SPECIES(I).EQ."NA+       ") NNAx    = I
         IF(SPECIES(I).EQ."Mg        ") NMG     = I
         IF(SPECIES(I).EQ."MG        ") NMG     = I
         IF(SPECIES(I).EQ."Mg+       ") NMGx    = I
         IF(SPECIES(I).EQ."MG+       ") NMGx    = I
         IF(SPECIES(I).EQ."Si        ") NSI     = I
         IF(SPECIES(I).EQ."SI        ") NSI     = I
         IF(SPECIES(I).EQ."Si+       ") NSIx    = I
         IF(SPECIES(I).EQ."SI+       ") NSIx    = I
         IF(SPECIES(I).EQ."Fe        ") NFE     = I
         IF(SPECIES(I).EQ."FE        ") NFE     = I
         IF(SPECIES(I).EQ."Fe+       ") NFEx    = I
         IF(SPECIES(I).EQ."FE+       ") NFEx    = I
         IF(SPECIES(I).EQ."Cl        ") NCL     = I
         IF(SPECIES(I).EQ."CL        ") NCL     = I
         IF(SPECIES(I).EQ."Cl+       ") NCLx    = I
         IF(SPECIES(I).EQ."CL+       ") NCLx    = I
         IF(SPECIES(I).EQ."Ca        ") NCA     = I
         IF(SPECIES(I).EQ."CA        ") NCA     = I
         IF(SPECIES(I).EQ."Ca+       ") NCAx    = I
         IF(SPECIES(I).EQ."CA+       ") NCAx    = I
         IF(SPECIES(I).EQ."Ca++      ") NCAxx   = I
         IF(SPECIES(I).EQ."CA++      ") NCAxx   = I
         IF(SPECIES(I).EQ."CO        ") NCO     = I
         IF(SPECIES(I).EQ."CH        ") NCH     = I
         IF(SPECIES(I).EQ."CH2       ") NCH2    = I
         IF(SPECIES(I).EQ."OH        ") NOH     = I
         IF(SPECIES(I).EQ."O2        ") NO2     = I
         IF(SPECIES(I).EQ."CS        ") NCS     = I
         IF(SPECIES(I).EQ."H2O       ") NH2O    = I
         IF(SPECIES(I).EQ."H3O+      ") NH3Ox   = I
         IF(SPECIES(I).EQ."HCO+      ") NHCOx   = I
         IF(SPECIES(I).EQ."e-        ") NELECT  = I
         IF(SPECIES(I).EQ."ELECTR    ") NELECT  = I
      ENDDO

      I=I-1
      READ(SPECIESFILE,*,END=1)
      I=I+1
 1    IF(I.NE.NSPEC) THEN
         write(6,*) 'ERROR! Number of species (NSPEC) does not match ',&
     &           'the number of entries in the species file'
         STOP
      ENDIF

!C     Check that the final species in the file is e-. Print a warning
!C     message to screen and logfile if not.
      IF(SPECIES(NSPEC).NE."e-") THEN
         write(6,*) 'WARNING! Last entry in species file is not e-'
         WRITE(10,*)'WARNING! Last entry in species file is not e-'
      ENDIF


!C     Check that the total hydrogen nuclei abundance adds up to 1.
!C     If not, modify the abundance of H2 (only consider H, H+ & H2)
      IF((ABUNDANCE(NH)+ABUNDANCE(NPROTON)+2.0D0*ABUNDANCE(NH2)).NE.1.0D0) THEN
         ABUNDANCE(NH2)=0.5D0*(1.0D0-ABUNDANCE(NH)-ABUNDANCE(NPROTON))
      ENDIF

!C     Calculate the intial electron abundance, if not 
!C     specified, as the sum of the metal ion abundances
      IF(ABUNDANCE(NELECT).LE.0.0D0) THEN
       ABUNDANCE(NELECT)=0.0D0
       IF(NCx.NE.0) ABUNDANCE(NELECT)=ABUNDANCE(NELECT)+ABUNDANCE(NCx)
       IF(NSx.NE.0) ABUNDANCE(NELECT)=ABUNDANCE(NELECT)+ABUNDANCE(NSx)
       IF(NNAx.NE.0) ABUNDANCE(NELECT)=ABUNDANCE(NELECT)+ABUNDANCE(NNAx)
       IF(NMGx.NE.0) ABUNDANCE(NELECT)=ABUNDANCE(NELECT)+ABUNDANCE(NMGx)
       IF(NSIx.NE.0) ABUNDANCE(NELECT)=ABUNDANCE(NELECT)+ABUNDANCE(NSIx)
       IF(NFEx.NE.0) ABUNDANCE(NELECT)=ABUNDANCE(NELECT)+ABUNDANCE(NFEx)
       IF(NCLx.NE.0) ABUNDANCE(NELECT)=ABUNDANCE(NELECT)+ABUNDANCE(NCLx)
       IF(NCAx.NE.0) ABUNDANCE(NELECT)=ABUNDANCE(NELECT)+ABUNDANCE(NCAx)
       IF(NCAxx.NE.0) ABUNDANCE(NELECT)=ABUNDANCE(NELECT)+2.0D0*ABUNDANCE(NCAxx)
       IF(NPROTON.NE.0) ABUNDANCE(NELECT)=ABUNDANCE(NELECT)+ABUNDANCE(NPROTON)
      ENDIF

      CLOSE(SPECIESFILE)
      RETURN
      END SUBROUTINE
