MODULE healpix_module

!T.Bisbas
  use definitions
  use healpix_types
  integer(kind=i4b):: ns_max                      ! ..
  integer(kind=i4b), dimension(0:1023) :: pix2x   ! ..
  integer(kind=i4b), dimension(0:1023) :: pix2y   ! ..
  real(kind=DP) :: x2pix(0:1023)          ! ..
  real(kind=DP) :: y2pix(0:1023)          ! ..
  real(kind=DP), allocatable :: vertex(:,:)      ! ..
  real(kind=DP), allocatable :: vector(:)      ! ..
END MODULE healpix_module

MODULE maincode_module
  USE ISO_C_BINDING
  use definitions
  use healpix_types

  integer(kind=I4B) :: i            ! counter
  integer :: ii           ! counter
  integer(kind=I4B) :: k            ! counter
  integer(kind=I4B) :: kk           ! counter
  integer(kind=I4B) :: p            ! counter
  integer(kind=I4B) :: ipix         ! pix id
  integer(kind=I4B) :: level        ! current level
  integer(kind=I4B) :: nrays        ! no. of rays on current level
  integer(kind=I4B) :: nside        ! refer to healpix manual
  integer(kind=I4B) :: itot         ! total number of grid points 
  integer(kind=I4B) :: ktot         ! total number of grid points 
  integer(kind=I4B) :: tot_eval     ! total number of evaluation points
  INTEGER(KIND=I4B) :: NLEV,NTEMP
  integer(kind=i4b) :: iteration, ITERTOT
  integer(kind=i4b) :: iterstep     ! output interval (per how many iterations)
  integer(kind=i4b) :: NSPEC, NREAC
  integer(kind=i4b) :: CII_NLEV, CII_NTEMP   !CII cooling variables
  integer(kind=i4b) :: CI_NLEV, CI_NTEMP     !CI cooling variables
  integer(kind=i4b) :: OI_NLEV, OI_NTEMP     !OI cooling variables
  integer(kind=i4b) :: C12O_NLEV, C12O_NTEMP !12CO cooling variables
  integer(kind=i4b) :: maxpoints
  integer(kind=i4b) :: suma
  integer :: CHEMITERATIONS
  integer :: pdr_tot
  
  logical::writeiterations
  logical::SPH

  integer(kind=I4B), allocatable :: rb(:)       ! ray with the IDs of the input file
  integer(kind=I4B), allocatable :: rrb(:)       ! ray with the IDs of the input file

  integer(kind=I4B) :: iwork(1:1000)
  real(kind=dp) :: rwork(1:1000000)

  real(kind=dp) :: ZETA
  real(kind=DP) :: angle_los   ! line-of-sight angle
  real(kind=dp) :: radius      ! maximum distance from the origin(1:3)
  real(kind=dp) :: theta_crit  ! critical theta angle to produce evaluation point
  real(kind=DP) :: rvec(1:3)   ! local 3 column array for the x,y,z coordinates of the grid point
  real(kind=DP) :: theta, phi  ! used from healpix to convert cartesian -> spherical
  real(kind=dp) :: adaptive_step !adaptive step for calculations of UVfield and column density.
  real(kind=dp) :: rhs1        ! right-hand-side 1 (to make integration calculation simpler)
  real(kind=DP) :: rhs2        ! -ditto-
  real(kind=DP) :: points      ! dummy variable used to find itot
  real(kind=DP) :: origin(1:3) ! grid point at which we perform calculation
  real(kind=DP) :: healpixvector(1:3) ! Healpix rays
  real(kind=dp) :: ENERGY,WEIGHT,EINSTEINA,FREQUENCY
  real(kind=dp) :: Tguess
  real(kind=dp) :: n_H
  real(kind=dp) :: CII_Z_function
  real(kind=dp) :: CI_Z_function
  real(kind=dp) :: OI_Z_function
  real(kind=dp) :: C12O_Z_function
  real(kind=dp) :: Z_increment
  real(kind=dp) :: frac1, frac2, frac3, tau_increment
  real(kind=dp) :: beta_ij_ray, beta_ij_sum, beta_ij
  real(kind=dp) :: TPOP, TMP2, BB_ij
  real(kind=dp) :: relch, v_turb, v_turb_inp
  real(kind=dp) :: Tlow0
  real(kind=dp) :: Thigh0
  real(kind=dp) :: Tmin
  real(kind=dp) :: Tmax
  real(kind=dp) :: Fcrit
  real(kind=dp) :: Tdiff
  real(kind=dp) :: dust_temperature
!  real(kind=dp) :: h
  real(kind=dp) :: avmax

  character(len=3) :: fieldchoice
  real(kind=dp) :: Gext(1)
  real(kind=dp) :: AV_fac, UV_fac

  real(kind=DP), allocatable :: vectors(:,:) ! Healpix rays
  real(kind=DP), allocatable :: ep(:,:)            ! evaluation point along each ray (local)
  real(kind=dp), allocatable :: ra(:)              ! distance    
  real(kind=dp), allocatable :: density(:)         ! density of each grid point - from input
  real(kind=DP), allocatable :: c_dens(:)          ! column density
  real(kind=dp), allocatable :: rra(:)             ! distance    
  real(kind=dp), allocatable :: COEFF(:)
  real(kind=dp), allocatable :: ENERGIES(:), WEIGHTS(:)
  real(kind=dp), allocatable :: A_COEFFS(:,:), B_COEFFS(:,:), C_COEFFS(:,:)
  real(kind=dp), allocatable :: FREQUENCIES(:,:), TEMPERATURES(:)
  real(kind=dp), allocatable :: H_COL(:,:,:)
  real(kind=dp), allocatable :: EL_COL(:,:,:)
  real(kind=dp), allocatable :: HE_COL(:,:,:)
  real(kind=dp), allocatable :: H2_COL(:,:,:)
  real(kind=dp), allocatable :: PH2_COL(:,:,:)
  real(kind=dp), allocatable :: OH2_COL(:,:,:)
  real(kind=dp), allocatable :: tau_ij(:)
  real(kind=dp), allocatable :: field(:,:)
  real(kind=dp), allocatable :: transition_CII(:,:),transition_CI(:,:)
  real(kind=dp), allocatable :: transition_OI(:,:),transition_C12O(:,:)
!CII cooling variables
  real(kind=dp) :: CII_COOLING
  real(kind=dp), allocatable :: CII_ENERGIES(:),CII_WEIGHTS(:)
  real(kind=dp), allocatable :: CII_A_COEFFS(:,:),CII_B_COEFFS(:,:),CII_C_COEFFS(:,:)
  real(kind=dp), allocatable :: CII_FREQUENCIES(:,:)
  real(kind=dp), allocatable :: CII_TEMPERATURES(:,:), CII_HP(:,:,:)
  real(kind=dp), allocatable :: CII_H(:,:,:),CII_EL(:,:,:)
  real(kind=dp), allocatable :: CII_HE(:,:,:),CII_H2(:,:,:)
  real(kind=dp), allocatable :: CII_PH2(:,:,:),CII_OH2(:,:,:)
!CI cooling variables
  real(kind=dp) :: CI_COOLING
  real(kind=dp), allocatable :: CI_ENERGIES(:),CI_WEIGHTS(:)
  real(kind=dp), allocatable :: CI_A_COEFFS(:,:),CI_B_COEFFS(:,:),CI_C_COEFFS(:,:)
  real(kind=dp), allocatable :: CI_FREQUENCIES(:,:)
  real(kind=dp), allocatable :: CI_TEMPERATURES(:,:),CI_HP(:,:,:)
  real(kind=dp), allocatable :: CI_H(:,:,:),CI_EL(:,:,:)
  real(kind=dp), allocatable :: CI_HE(:,:,:),CI_H2(:,:,:)
  real(kind=dp), allocatable :: CI_PH2(:,:,:),CI_OH2(:,:,:)
!OI cooling variables
  real(kind=dp) :: OI_COOLING
  real(kind=dp), allocatable :: OI_ENERGIES(:),OI_WEIGHTS(:)
  real(kind=dp), allocatable :: OI_A_COEFFS(:,:),OI_B_COEFFS(:,:),OI_C_COEFFS(:,:)
  real(kind=dp), allocatable :: OI_FREQUENCIES(:,:)
  real(kind=dp), allocatable :: OI_TEMPERATURES(:,:),OI_HP(:,:,:)
  real(kind=dp), allocatable :: OI_H(:,:,:),OI_EL(:,:,:)
  real(kind=dp), allocatable :: OI_HE(:,:,:),OI_H2(:,:,:)
  real(kind=dp), allocatable :: OI_PH2(:,:,:),OI_OH2(:,:,:)
!12CO cooling variables
  real(kind=dp) :: C12O_COOLING
  real(kind=dp), allocatable :: C12O_ENERGIES(:),C12O_WEIGHTS(:)
  real(kind=dp), allocatable :: C12O_A_COEFFS(:,:),C12O_B_COEFFS(:,:),C12O_C_COEFFS(:,:)
  real(kind=dp), allocatable :: C12O_FREQUENCIES(:,:)
  real(kind=dp), allocatable :: C12O_TEMPERATURES(:,:),C12O_HP(:,:,:)
  real(kind=dp), allocatable :: C12O_H(:,:,:),C12O_EL(:,:,:)
  real(kind=dp), allocatable :: C12O_HE(:,:,:),C12O_H2(:,:,:)
  real(kind=dp), allocatable :: C12O_PH2(:,:,:),C12O_OH2(:,:,:)

  character(len=20) :: input
  character(len=20) :: inputchem
  character(len=20) :: C12Oinput, CIIinput, CIinput, OIinput
  character(len=20)  :: output

  type columndens_node
     real(kind=dp), pointer :: columndens_point(:,:)
  end type columndens_node
  type(columndens_node), allocatable :: column(:)
  integer(kind=i4b), allocatable :: DUPLICATE(:)
  
  real(kind=dp),allocatable :: ALPHA(:),BETA(:),GAMMA(:),RATE(:),RTMIN(:),RTMAX(:)
  CHARACTER(len=10), allocatable :: REACTANT(:,:),PRODUCT(:,:)
  real(kind=dp), allocatable :: MASS(:),dummyabundance(:)!,abundances(:,:)
  CHARACTER(len=10), allocatable :: SPECIES(:)

  real(kind=dp),bind(c,name='maincode_module_mp_start_time_'):: start_time
  real(kind=dp),bind(c,name='maincode_module_mp_end_time_'):: end_time
  integer(kind=i4b) :: status

  real(kind=dp), allocatable :: dusttemperature(:)
  real(kind=dp), allocatable :: gastemperature(:)
  real(kind=dp), allocatable :: previousgastemperature(:)

#ifdef THERMALBALANCE
  real(kind=dp), allocatable :: Fmean(:)
  real(kind=dp), allocatable :: Fratio(:)
  real(kind=dp) :: dummytemperature
  real(kind=dp), allocatable :: Tlow(:)
  real(kind=dp), allocatable :: Thigh(:)
#endif

  type pdr_node
     integer(kind=i4b), pointer :: epray(:)        !population of evaluation points per ray
     integer(kind=i4b), pointer :: projected(:,:)  !ID of projected grid points in the line of sight
     integer(kind=i4b), pointer :: raytype(:)      !raytype
     real(kind=dp), pointer :: epoint(:,:,:)       !co-ordinates of each evaluation point
     real(kind=dp), pointer :: columndensity(:)    !column density
     real(kind=dp), pointer :: AV(:)               !AV
     real(kind=dp), pointer :: rad_surface(:)      !rad_surface
     real(kind=dp), pointer :: abundance(:)        !abundance of species
     real(kind=dp), pointer :: CII_pop(:)          !CII population density
     real(kind=dp), pointer :: CI_pop(:)           !CI population density
     real(kind=dp), pointer :: OI_pop(:)           !OI population density
     real(kind=dp), pointer :: C12O_pop(:)         !CO population density
     real(kind=dp), pointer :: CII_line(:,:)       !CII line cooling
     real(kind=dp), pointer :: CI_line(:,:)        !CI line cooling
     real(kind=dp), pointer :: OI_line(:,:)        !OI line cooling
     real(kind=dp), pointer :: C12O_line(:,:)      !C12O line cooling
!===============
     real(kind=dp), pointer :: CII_optdepth(:,:,:)       !CII line cooling
     real(kind=dp), pointer :: CI_optdepth(:,:,:)        !CI line cooling
     real(kind=dp), pointer :: OI_optdepth(:,:,:)        !OI line cooling
     real(kind=dp), pointer :: C12O_optdepth(:,:,:)      !C12O line cooling
     real(kind=dp), pointer :: C12O_beta(:,:,:)      !C12O line escape probability
     real(kind=dp), pointer :: CII_beta(:,:,:)      !C12O line escape probability
     real(kind=dp), pointer :: CI_beta(:,:,:)      !C12O line escape probability
     real(kind=dp), pointer :: OI_beta(:,:,:)      !C12O line escape probability
     real(kind=dp), pointer :: CII_CCOEFFS(:,:)
     real(kind=dp), pointer :: check_step(:)
     real(kind=dp), pointer :: check_coef(:,:,:)
!===============
     real(kind=dp) :: UVfield                      !attenuated UV field of element
     real(kind=dp) :: rho                          !density of element
     real(kind=dp) :: smoo                         !smoothing length (if SPH = .TRUE. in params.dat)
     real(kind=dp) :: x,y,z                        !position of element
     integer(kind=i4b) :: etype                    !element type (1 = PDR, 2 = ION, 3 = DARK)
#ifdef DUST2
     real(kind=dp) :: dust_t
#endif
  end type pdr_node 
  type (pdr_node), allocatable :: pdr(:)      !main 3DPDR array for each grid point p
  real(kind=dp), allocatable :: CII_solution(:,:), CIIsolution(:)
  real(kind=dp), allocatable :: CI_solution(:,:), CIsolution(:)
  real(kind=dp), allocatable :: OI_solution(:,:), OIsolution(:)
  real(kind=dp), allocatable :: C12O_solution(:,:), C12Osolution(:)
  real(kind=dp), allocatable :: CIIevalpop(:,:,:)
  real(kind=dp), allocatable :: CIevalpop(:,:,:)
  real(kind=dp), allocatable :: OIevalpop(:,:,:)
  real(kind=dp), allocatable :: C12Oevalpop(:,:,:)

  integer(kind=i4b)::levpop_iteration
  real(kind=dp) :: rho_min
  real(kind=dp) :: rho_max
#ifdef GUESS_TEMP
  real(kind=dp),allocatable :: Tmin_array(:)
  real(kind=dp),allocatable :: Tmax_array(:)
#endif

!================================
!================================
!================================
real::time_evalpoints,t2,t3,t3b,t4 !times for the cpu_time subroutine
real(kind=dp)::percentage,thermal_percentage
real(kind=dp)::OI_percentage,C12O_percentage
real(kind=dp)::CII_percentage,CI_percentage
#ifdef CO_FIX
real(kind=dp)::levpop_percentage
#endif
real(kind=dp) :: rad_tot
character(len=32) :: out_file, out_file2
character(len=7) :: file_ext
character(len=6) :: file_numb
real(kind=dp),allocatable :: CII_cool(:)
real(kind=dp),allocatable :: CI_cool(:)
real(kind=dp),allocatable :: OI_cool(:)
real(kind=dp),allocatable :: C12O_cool(:)
real(kind=dp),allocatable :: total_cooling_rate(:)
real(kind=dp),allocatable :: dummyarray_CII(:,:), dummyarray_CI(:,:)
real(kind=dp),allocatable :: dummyarray_OI(:,:), dummyarray_C12O(:,:)
!=================
real(kind=dp),allocatable :: dummyarray_CII_tau(:,:,:), dummyarray_CI_tau(:,:,:)
real(kind=dp),allocatable :: dummyarray_OI_tau(:,:,:), dummyarray_C12O_tau(:,:,:)
real(kind=dp),allocatable :: dummyarray_CII_beta(:,:,:),dummyarray_CI_beta(:,:,:)
real(kind=dp),allocatable :: dummyarray_OI_beta(:,:,:),dummyarray_C12O_beta(:,:,:)
real(kind=dp),allocatable :: dummycoef(:,:,:), dummystep(:)
!=================
real(kind=dp),allocatable :: total_heating(:)
integer(kind=i4b) :: referee, id
#ifdef THERMALBALANCE
logical,allocatable :: converged(:),doleveltmin(:)
!logical :: show_balance
#endif
logical,allocatable :: level_converged(:)
logical :: level_conv,first_time
logical,allocatable :: CII_conv(:),CI_conv(:),OI_conv(:),C12O_conv(:)
real(kind=dp),allocatable :: CII_RELCH(:,:),CI_RELCH(:,:),OI_RELCH(:,:),C12O_RELCH(:,:)
logical :: relch_conv
integer(kind=i4b) :: CII_i,CI_i,OI_i,C12O_i
#ifdef OPENMP
integer::CPUs
#endif
logical,allocatable::expanded(:)
integer(kind=i4b) :: grand_ptot
integer(kind=i4b) :: pdr_ptot,ion_ptot,dark_ptot
integer(kind=i4b),allocatable :: IDlist_pdr(:),IDlist_ion(:),IDlist_dark(:)
real(kind=DP), allocatable :: pdrpoint(:,:)      ! coordinates of pdr element 
integer(kind=i4b) :: pp
real(kind=dp) :: xpos,ypos,zpos,denst
!================================
!================================
!================================

END MODULE maincode_module


module uclpdr_module

  use definitions
  use healpix_types

  integer(kind=i4b), save :: NUMH2=105
  real(kind=dp), dimension(105), save :: COL_GRID=(/&
     &              0.000D+00,3.690D+11,3.715D+12,3.948D+13,1.233D+14, &
     &              2.536D+14,4.342D+14,6.653D+14,6.689D+14,9.075D+14, &
     &              1.234D+15,1.631D+15,2.105D+15,2.363D+15,2.899D+15, &
     &              3.207D+15,3.848D+15,4.636D+15,5.547D+15,6.604D+15, &
     &              7.855D+15,9.368D+15,1.122D+16,1.352D+16,1.643D+16, &
     &              2.017D+16,2.515D+16,3.190D+16,4.128D+16,5.439D+16, &
     &              7.315D+16,1.009D+17,1.432D+17,2.092D+17,3.123D+17, &
     &              4.738D+17,5.388D+17,8.935D+17,1.381D+18,2.164D+18, &
     &              3.330D+18,5.024D+18,7.404D+18,9.029D+18,1.316D+19, &
     &              1.813D+19,2.453D+19,3.248D+19,4.216D+19,5.370D+19, &
     &              6.722D+19,8.277D+19,9.894D+19,1.186D+20,1.404D+20, &
     &              1.644D+20,1.908D+20,2.197D+20,2.510D+20,2.849D+20, &
     &              3.214D+20,3.604D+20,4.019D+20,4.456D+20,4.915D+20, &
     &              5.393D+20,5.886D+20,6.392D+20,6.909D+20,7.433D+20, &
     &              7.965D+20,8.505D+20,9.056D+20,9.627D+20,1.011D+21, &
     &              1.068D+21,1.125D+21,1.185D+21,1.250D+21,1.327D+21, &
     &              1.428D+21,1.578D+21,1.851D+21,2.128D+21,2.298D+21, &
     &              2.389D+21,2.459D+21,2.519D+21,2.571D+21,2.618D+21, &
     &              2.707D+21,2.790D+21,2.887D+21,3.001D+21,3.139D+21, &
     &              3.303D+21,3.497D+21,3.722D+21,3.983D+21,4.283D+21, &
     &              4.644D+21,5.127D+21,5.945D+21,8.205D+21,1.015D+22/)
  real(kind=dp), dimension(105), save :: SH2_GRID=(/&
     &              1.000D+00,9.983D-01,9.853D-01,8.761D-01,7.199D-01, &
     &              5.728D-01,4.455D-01,3.431D-01,3.418D-01,2.732D-01, &
     &              2.110D-01,1.619D-01,1.236D-01,1.084D-01,8.447D-02, &
     &              7.410D-02,5.774D-02,4.416D-02,3.390D-02,2.625D-02, &
     &              2.048D-02,1.606D-02,1.264D-02,9.987D-03,7.937D-03, &
     &              6.343D-03,5.088D-03,4.089D-03,3.283D-03,2.640D-03, &
     &              2.130D-03,1.725D-03,1.397D-03,1.129D-03,9.097D-04, &
     &              7.340D-04,6.883D-04,5.377D-04,4.352D-04,3.475D-04, &
     &              2.771D-04,2.205D-04,1.753D-04,1.549D-04,1.210D-04, &
     &              9.666D-05,7.705D-05,6.148D-05,4.904D-05,3.909D-05, &
     &              3.112D-05,2.473D-05,1.997D-05,1.578D-05,1.244D-05, &
     &              9.769D-06,7.634D-06,5.932D-06,4.581D-06,3.515D-06, &
     &              2.679D-06,2.029D-06,1.527D-06,1.144D-06,8.523D-07, &
     &              6.332D-07,4.693D-07,3.475D-07,2.574D-07,1.907D-07, &
     &              1.413D-07,1.047D-07,7.739D-08,5.677D-08,4.386D-08, &
     &              3.227D-08,2.385D-08,1.750D-08,1.248D-08,8.389D-09, &
     &              5.026D-09,2.382D-09,6.259D-10,1.653D-10,7.399D-11, &
     &              4.824D-11,3.474D-11,2.633D-11,2.069D-11,1.663D-11, &
     &              1.099D-11,7.506D-12,4.825D-12,2.864D-12,1.534D-12, &
     &              7.324D-13,3.087D-13,1.135D-13,3.591D-14,9.689D-15, &
     &              2.045D-15,2.618D-16,8.918D-18,3.041D-21,1.739D-23/)
  logical :: start
  integer(kind=i4b), save :: DIMH2=6
  integer(kind=i4b), save :: DIMCO=8
  real(kind=dp), dimension(8), save :: NCO_GRID=(/&
    &12.0D0,13.0D0,14.0D0,15.0D0, 16.0D0,17.0D0,18.0D0,19.0D0/)
  real(kind=dp), dimension(6), save :: NH2_GRID=(/&
    &18.0D0,19.0D0,20.0D0,21.0D0,22.0D0,23.0D0/)
  real(kind=dp) :: SCO_GRID(1:8,1:6)
  integer(kind=i4b), save :: N_GRID=30
  real(kind=dp), dimension(30), save :: L_GRID=(/&
     &             910.0D0, 950.0D0,1000.0D0,1050.0D0,1110.0D0, &
     &            1180.0D0,1250.0D0,1390.0D0,1490.0D0,1600.0D0, &
     &            1700.0D0,1800.0D0,1900.0D0,2000.0D0,2100.0D0, &
     &            2190.0D0,2300.0D0,2400.0D0,2500.0D0,2740.0D0, &
     &            3440.0D0,4000.0D0,4400.0D0,5500.0D0,7000.0D0, &
     &            9000.0D0,12500.0D0,22000.0D0,34000.0D0,1.0D9/)
  real(kind=dp), dimension(30), save :: X_GRID=(/&
     &            5.76D0,5.18D0,4.65D0,4.16D0,3.73D0, &
     &            3.40D0,3.11D0,2.74D0,2.63D0,2.62D0, &
     &            2.54D0,2.50D0,2.58D0,2.78D0,3.01D0, &
     &            3.12D0,2.86D0,2.58D0,2.35D0,2.00D0, &
     &            1.58D0,1.42D0,1.32D0,1.00D0,0.75D0, &
     &            0.48D0,0.28D0,0.12D0,0.05D0,0.00D0/)

  real(kind=dp) :: SH2_DERIV(1:105), SCO_DERIV(1:8,1:6), X_DERIV(1:30)

end module uclpdr_module

module global_module

  USE ISO_C_BINDING
  use definitions
  use healpix_types

  INTEGER(kind=i4b) :: NH,ND,NH2,NHD,NC,NCx,NCO,NO,NPROTON,NH2O,NHe, &
    &        NMG,NMGx,NN,NFE,NFEx,NSI,NSIx,NCA,NCAx,NCAxx,NS,NSx,NCS, &
    &        NOSH,NCL,NCLx,NH2x,NHEx,NOx,NNx,NNA,NNAx,NCH,NCH2,NOH,NO2, &
    &        NH3x, NH3Ox, NHCOx
  integer(kind=i4b),bind(c,name='global_module_mp_nelect_')::NELECT

  !INTEGER(kind=i4b) :: NRGR,NRH2,NRHD,NRCO,NRCI,NRSI

! REAL(kind=dp), save :: ZETA=3.85D0,OMEGA=0.42D0,GRAIN_RADIUS=1.0D-5,METALLICITY=1.0D0
 REAL(kind=dp) :: g2d
 REAL(kind=dp) :: metallicity
 REAL(kind=dp) :: omega
 REAL(kind=dp) :: grain_radius
! REAL(kind=dp), save :: OMEGA=0.42D0,GRAIN_RADIUS=1.0D-7,METALLICITY=1.0D0
! REAL(kind=dp), save :: ZETA=1.0D0,OMEGA=0.42D0,GRAIN_RADIUS=1.0D-5,METALLICITY=1.0D0

 real(kind=dp),allocatable :: all_heating(:,:)
 real(kind=dp),allocatable :: allheating(:) 
 
end module global_module

module functions_module

 use definitions
 use healpix_types
 implicit none
  
 interface
   function H2PDRATE(K0,G0,AV,NH2)
     use definitions
     use healpix_types
     real(kind=dp) :: H2PDRATE
     real(kind=dp), intent(in) :: k0, g0, av
     real(kind=dp), intent(in) :: nh2
     real(kind=dp) :: lambda, scatter, h2shield2
   end function H2PDRATE

   function COPDRATE(K0,G0,AV,NCO,NH2)
     use definitions
     use healpix_types
     real(kind=dp) :: copdrate
     real(kind=dp), intent(in) :: k0, g0, av, nco
     real(kind=dp), intent(in) :: nh2
     real(kind=dp) :: lambda, lbar, coshield, scatter
   end function copdrate

   function CIPDRATE(K0,G0,AV,KAV,NCI,NH2,TGAS)
     use definitions
     use healpix_types
     real(kind=dp) :: cipdrate
     real(kind=dp), intent(in) :: K0,G0,AV,KAV,NCI,TGAS
     real(kind=dp), intent(in) :: nh2
     real(kind=dp) :: tauc
   end function cipdrate

   function SIPDRATE(K0,G0,AV,KAV,NSI)
     use definitions
     use healpix_types
     real(kind=dp) :: sipdrate
     real(kind=dp), intent(in) :: K0,G0,AV,KAV,NSI
     real(kind=dp) :: taus
   end function sipdrate

   function H2SHIELD1(NH2,DOPW,RADW)
     use definitions
     use healpix_types
     real(kind=dp) :: h2shield1
     real(kind=dp), intent(in) :: nh2
     real(kind=dp), intent(in) ::DOPW,RADW
     real(kind=dp) :: FPARA, FOSC, TAUD, R, T, U, JD, JR
   end function h2shield1

   function h2shield2(nh2)
     use definitions
     use healpix_types
     use uclpdr_module, only : start, numh2, COL_GRID, SH2_GRID, SH2_DERIV
     real(kind=dp) :: h2shield2
     real(kind=dp), intent(in) :: nh2
   end function h2shield2

   function COSHIELD(NCO,NH2)  
     use definitions
     use healpix_types
     use uclpdr_module, only : start, NCO_GRID, NH2_GRID, SCO_GRID, SCO_DERIV
     real(kind=dp) :: COSHIELD
     real(kind=dp) :: LOGNCO, LOGNH2
     real(kind=dp), intent(in) :: NCO, NH2
   end function COSHIELD

   function SCATTER(AV,LAMBDA)
     use definitions
     use healpix_types
     real(kind=dp) :: scatter
     real(kind=dp), intent(in) :: AV, LAMBDA
     real(kind=dp), dimension(0:5), save :: A = (/&
           &1.000D0,2.006D0,-1.438D0,0.7364D0,-0.5076D0,-0.0592D0/)
     real(kind=dp), dimension(0:5), save :: K = (/&
           &0.7514D0,0.8490D0,1.013D0,1.282D0,2.005D0,5.832D0/)
     real(kind=dp) :: EXPONENT, XLAMBDA
   end function scatter

   function XLAMBDA(LAMBDA)
     use definitions
     use healpix_types
     use uclpdr_module, only : start, N_GRID, L_GRID, X_GRID, X_DERIV
     real(kind=dp) :: xlambda
     real(kind=dp), intent(in) :: lambda
   end function xlambda

   function LBAR(NCO,NH2)
     use definitions
     use healpix_types
     real(kind=dp) :: lbar
     real(kind=dp) :: U,W 
     real(kind=dp) :: NCO, NH2
   end function LBAR

    function calculate_heating(density, gas_temperature, dust_temperature, UV_field, &
             & v_turb, nspec, dummyabundance, nreac, rate)
      use definitions
      use healpix_types
      real(kind=dp) :: calculate_heating
      integer(kind=i4b) :: nspec, nreac
      real(kind=dp) :: density, gas_temperature, dust_temperature, UV_field, v_turb
      real(kind=dp) :: dummyabundance(1:nspec), rate(1:nreac)
    end function calculate_heating

#ifdef H2FORM

      FUNCTION H2_FORMATION_RATE(GAS_TEMPERATURE,GRAIN_TEMPERATURE) RESULT(RATE)
         USE DEFINITIONS
         USE HEALPIX_TYPES
         IMPLICIT NONE
         REAL(KIND=DP) :: RATE
         REAL(KIND=DP), INTENT(IN) :: GAS_TEMPERATURE,GRAIN_TEMPERATURE
      END FUNCTION H2_FORMATION_RATE

#endif




 end interface

end module functions_module

module chemistry_module
  USE ISO_C_BINDING
  use healpix_types

  real(kind=dp),bind(c,name='chemistry_module_mp_relative_abundance_tolerance_') :: relative_abundance_tolerance
  real(kind=dp),bind(c,name='chemistry_module_mp_absolute_abundance_tolerance_') :: absolute_abundance_tolerance
end module chemistry_module


