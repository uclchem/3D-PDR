!C***********************************************************************
!C     Read in the chemical reaction rates and the species, masses and
!C     initial abundances (if specified). The  rates and species files
!C     are assumed to have comma separated values (CSV) format. This is
!C     in line with the Rate05 formatting, removing the need for file-
!C     dependent FORMAT statements.
!C***********************************************************************
!
      SUBROUTINE READ_RATES(NREAC,REAC,PROD,ALPHA,BETA,GAMMA,RATE,&
                           &DUPLICATE,RTMIN,RTMAX)
!T.Bell
use definitions
use healpix_types
use global_module

      IMPLICIT NONE
      INTEGER(kind=i4b), intent(in) :: NREAC
      integer(kind=i4b), intent(out) :: DUPLICATE(1:nreac)
      real(kind=dp), intent(out) :: ALPHA(1:nreac),BETA(1:nreac),&
              &GAMMA(1:nreac),RATE(1:nreac),RTMIN(1:nreac),RTMAX(1:nreac)
      CHARACTER(len=10), intent(out) :: REAC(1:nreac,1:3),PROD(1:nreac,1:4)

      INTEGER(kind=i4b) :: I,J,N,RATEFILE
      CHARACTER(len=1) :: CLEM
      RATEFILE = 2

!C     Initialize the variables and read in the ratefile data. Check that
!C     the value of NREAC agrees with the number of reactions in the file
!C     and produce an error message if not.
      REAC="          "
      PROD="          "
      ALPHA=0.0D0
      BETA=0.0D0
      GAMMA=0.0D0
      RTMIN=0.0D0
      RTMAX=0.0D0
      DUPLICATE=0

      RATE=0.0D0

       
#ifdef REDUCED
      OPEN(RATEFILE,FILE="rates_reduced.d",STATUS="OLD")
#endif
#ifdef FULL
      OPEN(RATEFILE,FILE="rates_full.d",STATUS="OLD")
#endif
#ifdef MYNETWORK
      OPEN(RATEFILE,FILE="rates_mynetwork.d",STATUS="OLD")
#endif
      REWIND(RATEFILE)
      DO I=1,NREAC
         READ(RATEFILE,*,END=1) N,(REAC(I,J),J=1,3),(PROD(I,J),J=1,4),&
     &                          ALPHA(I),BETA(I),GAMMA(I), &
     &                          CLEM,RTMIN(I),RTMAX(I)
         IF(CLEM.NE."") CLEM=""

!C     Check for duplicate reactions and set the DUPLICATE counter to the
!C     appropriate value. Adjust their minimum temperatures so that the
!C     temperature ranges are adjacent.
         IF(I.GT.1) THEN
           IF(REAC(I,1).EQ.REAC(I-1,1) .AND. &
     &     REAC(I,2).EQ.REAC(I-1,2) .AND. REAC(I,3).EQ.REAC(I-1,3) .AND. &
     &     PROD(I,1).EQ.PROD(I-1,1) .AND. PROD(I,2).EQ.PROD(I-1,2) .AND. &
     &     PROD(I,3).EQ.PROD(I-1,3) .AND. PROD(I,4).EQ.PROD(I-1,4)) THEN 
            IF(DUPLICATE(I-1).EQ.0) DUPLICATE(I-1)=1
            DUPLICATE(I)=DUPLICATE(I-1)+1
            RTMIN(I)=RTMAX(I-1)
           ELSE
            DUPLICATE(I)=0
           ENDIF
         ELSE
            DUPLICATE(I)=0
         ENDIF

!C     Check for negative gamma values as they could cause problems when
!C     calculating abundances. Produce a warning message if they occur.
         IF(GAMMA(I).LT.0.0D0) THEN
          write(6,*) 'Negative gamma factor in rate',N
          WRITE(10,"('Negative gamma factor in rate',I5,' (',F8.1,')')")&
     &         N,GAMMA(I)
         ENDIF
      ENDDO
      I=I-1
      READ(RATEFILE,*,END=1)
      I=I+1
 1    IF(I.NE.NREAC) THEN
         write(6,*) 'ERROR! Number of reactions (NREAC) does not match ', &
     &           'the number of entries in the ratefile'
         STOP
      ENDIF

      CLOSE(RATEFILE)
      RETURN
      END SUBROUTINE
!C-----------------------------------------------------------------------
