!-----------------------------------------------------------------------------
!
!  Copyright (C) 1997-2005 Krzysztof M. Gorski, Eric Hivon, 
!                          Benjamin D. Wandelt, Anthony J. Banday, 
!                          Matthias Bartelmann, Hans K. Eriksen, 
!                          Frode K. Hansen, Martin Reinecke
!
!
!  This file is part of HEALPix.
!
!  HEALPix is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!
!  HEALPix is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with HEALPix; if not, write to the Free Software
!  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
!
!  For more information about HEALPix see http://healpix.jpl.nasa.gov
!
!-----------------------------------------------------------------------------
MODULE healpix_types
  ! This module sets the types used in the Fortran 90 modules
  ! of the HEALPIX distribution and follows the example of Numerical Recipes
  !
  ! Benjamin D. Wandelt October 1997
  ! Eric Hivon June 1998
  ! Eric Hivon Oct  2001, edited to be compatible with 'F' compiler
  ! Eric Hivon July 2002, addition of i8b, i2b, i1b
  !                       addition of max_i8b, max_i2b and max_i1b
  !            Jan 2005, explicit form of max_i1b because of ifc 8.1.021
  !            June 2005, redefine i8b as 16 digit integer because of Nec f90 compiler

! Include definitions module for universal definition of DP, PR, SP
  use definitions

  INTEGER, PARAMETER, public :: i8b = SELECTED_INT_KIND(16)
  INTEGER, PARAMETER, public :: i4b = SELECTED_INT_KIND(9)
  INTEGER, PARAMETER, public :: i2b = SELECTED_INT_KIND(4)
  INTEGER, PARAMETER, public :: i1b = SELECTED_INT_KIND(2)
  INTEGER, PARAMETER, public :: lgt = KIND(.TRUE.)
  INTEGER, PARAMETER, public :: spc = KIND((1.0_sp, 1.0_sp))
  INTEGER, PARAMETER, public :: dpc = KIND((1.0_dp, 1.0_dp))
  !
  INTEGER(I8B),  PARAMETER, public :: max_i8b = HUGE(1_i8b)
  INTEGER,       PARAMETER, public :: max_i4b = HUGE(1_i4b)
  INTEGER,       PARAMETER, public :: max_i2b = HUGE(1_i2b)
  INTEGER,       PARAMETER, public :: max_i1b = 127
  REAL(kind=sp), PARAMETER, public :: max_sp  = HUGE(1.0_sp)
  REAL(kind=dp), PARAMETER, public :: max_dp  = HUGE(1.0_dp)

  ! Numerical Constant (Double precision)
  REAL(kind=dp), PARAMETER, public :: QUARTPI=0.785398163397448309615660845819875721049_dp
  REAL(kind=dp), PARAMETER, public :: HALFPI= 1.570796326794896619231321691639751442099_dp
  REAL(kind=dp), PARAMETER, public :: PI    = 3.141592653589793238462643383279502884197_dp
  REAL(kind=dp), PARAMETER, public :: TWOPI = 6.283185307179586476925286766559005768394_dp
  REAL(kind=dp), PARAMETER, public :: FOURPI=12.56637061435917295385057353311801153679_dp
  REAL(kind=dp), PARAMETER, public :: SQRT2 = 1.41421356237309504880168872420969807856967_dp
  REAL(kind=dp), PARAMETER, public :: SQ4PI_INV = 0.2820947917738781434740397257803862929220_dp
  REAL(kind=dp), PARAMETER, public :: TWOTHIRD = 0.6666666666666666666666666666666666666666_dp

  real(kind=DP), parameter, public :: RAD2DEG = 180.0_DP / PI
  real(kind=DP), parameter, public :: DEG2RAD = PI / 180.0_DP
  real(kind=SP), parameter, public :: hpx_sbadval = -1.6375e30_sp
  real(kind=DP), parameter, public :: hpx_dbadval = -1.6375e30_dp

  ! Maximum length of filenames
  integer, parameter :: filenamelen = 1024


  ! ---- Normalisation and convention ----
  ! normalisation of spin weighted functions
  real(kind=dp), parameter, public ::  KvS = 1.0_dp ! 1.0 : CMBFAST (Healpix 1.2)
  ! sign of Q
  real(kind=dp), parameter, public :: sgQ = -1.0_dp ! -1 : CMBFAST (Healpix 1.2)
  ! sign of spin weighted function !
  real(kind=dp), parameter, public :: SW1 = -1.0_dp ! -1 : Healpix 1.2, bug correction
  real(kind=dp), parameter, public :: iKvS = 1.0_dp / KvS  ! inverse of KvS

!parameters for 3DPDR
  real(kind=dp), parameter, public ::  KB = 1.38065040D-16 !Boltzmann constant cgs
  real(kind=dp), parameter, public ::  C  = 2.99792458D+10 !speed of light cgs
  real(kind=dp), parameter, public ::  MP = 1.67262164D-24 !proton mass cgs
  real(kind=dp), parameter, public ::  HP = 6.62606896D-27 !Planck's constant cgs
  real(kind=dp), parameter, public ::  HB = 1.05457163D-27 !Planck's constant / 2pi
  real(kind=dp), parameter, public ::  HK = 4.79923734D-11 !Planck's constant / Boltzmann constant
  real(kind=dp), parameter, public ::  NA = 6.02214179D+23 !Avogadro's number
  real(kind=dp), parameter, public ::  AU = 1.66053878D-24 !atomic mass unit
  real(kind=dp), parameter, public ::  MH = 1.67372346D-24 !hydrogen mass cgs
  real(kind=dp), parameter, public ::  ME = 9.10938215D-28 !electron mass cgs
  real(kind=dp), parameter, public ::  EC = 4.80320427D-10 !elementary charge in esu
  real(kind=dp), parameter, public ::  PC = 3.08568025D+18 !pc in cm
  real(kind=dp), parameter, public ::  EV = 1.60217646D-12 !electron volt in erg


END MODULE healpix_types
