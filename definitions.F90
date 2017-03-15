! Definition of precision kind variables
! ============================================================================

!#include "macros.h"

! ============================================================================
!T.Bisbas -- taken from SEREN/HEALPix

MODULE definitions

  integer, parameter :: DP = selected_real_kind(p=15) ! double precision
  integer, parameter :: SP = selected_real_kind(p=6)  ! single precision

#ifdef DOUBLE_PRECISION
  integer, parameter :: PR = DP                       ! particle precision
#else
  integer, parameter :: PR = SP                       ! default = single
#endif

  integer, parameter :: ILP = 4                       ! Integer long precision

END MODULE definitions


! ============================================================================
