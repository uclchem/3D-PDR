program threedpdr

use definitions
use healpix_types
use healpix_module
use maincode_module
use uclpdr_module, only : start, SCO_GRID
use global_module
#ifdef OPENMP
use omp_lib
#endif
use chemistry_module
  
logical::write_output
logical,allocatable::dobinarychop(:)
character(len=1),allocatable::previouschange(:)
!real(kind=dp),allocatable :: Rrate(:)
real(kind=dp)::neutralfraction,maximum_density,minimum_density
integer(kind=i4b)::mmpdr,j,newmaxpoints,ilevel,jlevel
real(kind=dp),allocatable::dummy_abundance(:,:), dummy_density(:)
real(kind=dp),allocatable::dummy_rate(:,:), dummy_temperature(:)
real(kind=dp)::tt1,tt2
integer(kind=i4b)::NRGR,NRH2,NRHD,NRCO,NRCI,NRSI
real(kind=dp)::uvfieldaux
real(kind=dp),allocatable::prev_cooling(:)
!reversing rows declaration
real(kind=dp),allocatable::x_rev(:),y_rev(:),z_rev(:),n_rev(:)

write(6,*) '=============================================================================='
write(6,*) '*********     *********             *********     *********      *********'
write(6,*) '        **    **       **           **       **   **       **    **       **'
write(6,*) '         **   **        **          **        **  **        **   **        **'
write(6,*) '        **    **         **         **       **   **         **  **       **'
write(6,*) '  *******     **         **  *****  *********     **         **  *********'
write(6,*) '        **    **         **         **            **         **  **     **'
write(6,*) '         **   **        **          **            **        **   **      **'
write(6,*) '        **    **       **           **            **       **    **       **'
write(6,*) '*********     *********             **            *********      **        **'
write(6,*) '=============================================================================='
write(6,*) '********************   Coders:   T.G.Bisbas, T.A.Bell   **********************'
write(6,*) '*************   Collaborators:   S.Viti, J.Yates, M.Barlow   *****************'
write(6,*) '*************************        Version 1.0          ************************'
write(6,*) '=============================================================================='
write(6,*) ''
write(6,*) 'Reading input [params.dat]'
call readparams
v_turb=v_turb_inp*1.0D5
write(6,*) 'Input file:               ',input
write(6,*) 'HEALPix level:            ',level
#ifdef PSEUDO_1D
if (level.gt.0) STOP "HEALPix level must be set to 0 in PSEUDO_1D mode"
#endif
#ifdef PSEUDO_2D
if (level.gt.0) STOP "HEALPix level must be set to 0 in PSEUDO_1D mode"
#endif
write(6,*) 'Theta critical:           ',theta_crit
write(6,*) 'Angle between rays:       ',sqrt(pi/3.0D0/4.0D0**(real(level)))
write(6,*) 'Maxpoints                 ',maxpoints
if (theta_crit.ge.pi/2.0D0) stop 'theta_crit must be less than pi/2'
write(6,*) 'Guess Temperature (K):    ',Tguess
write(6,*) 'Dust  Temperature (K):    ',dust_temperature
write(6,*) 'Turbulent velocity (cm/s):',v_turb
write(6,*) 'minimum density (cm^-3):  ',rho_min
write(6,*) 'maximum density (cm^-3):  ',rho_max
#ifdef THERMALBALANCE
write(6,*) 'Tlow:                     ',Tlow0
write(6,*) 'Thigh:                    ',Thigh0
#endif
write(6,*) 'Tmin:                     ',Tmin
write(6,*) 'Tmax:                     ',Tmax
write(6,*) 'Fcrit:                    ',Fcrit
write(6,*) 'Tdiff:                    ',Tdiff
start = .true.
write(6,*) 'Form of field:            ',fieldchoice
write(6,*) 'Gext:                     ',Gext(1)
write(6,*) 'AV factor:                ',AV_fac
write(6,*) 'UV factor:                ',UV_fac
#ifdef REDUCED
write(6,*) 'Chemical network:           REDUCED'
#endif
#ifdef FULL
write(6,*) 'Chemical network:           FULL'
#endif
#ifdef MYNETWORK
write(6,*) 'Chemical network:           MYNETWORK'
#endif
write(6,*) 'Number of species:        ',nspec
write(6,*) 'Number of reactions:      ',nreac
write(6,*) 'Total iterations:         ',itertot
write(6,*) 'Output interval / iter.:  ',iterstep
write(6,*) 'Chemiterations:           ',chemiterations
write(6,*) 'Zeta:                     ',zeta*1.3d-17
write(6,*) 'Gas-to-dust               ',g2d!*100.0
write(6,*) 'Metallicity               ',metallicity
write(6,*) 'Omega                     ',omega
write(6,*) 'Grain radius              ',grain_radius
close(12)
write(6,*) '============================================='
write(6,*) ''
write(6,*) '------FLAGS-----'
#ifdef THERMALBALANCE
write(6,*) 'THERMALBALANCE'
#endif
#ifdef OPENMP
write(6,*) 'OPENMP'
#endif
#ifdef PSEUDO_1D
write(6,*) 'PSEUDO_1D'
#elif PSEUDO_2D
write(6,*) 'PSEUDO_2D'
#else
write(6,*) 'FULL_3D'
#endif
#ifdef DUST
write(6,*) 'DUST'
#endif
#ifdef DUST2
write(6,*) 'DUST2'
#endif
#ifdef CO_FIX
write(6,*) 'CO_FIX'
#endif
#ifdef H2FORM
write(6,*) 'H2FORM'
#endif
#ifdef TEMP_FIX
write(6,*) 'TEMP_FIX'
#endif
#ifdef GUESS_TEMP
write(6,*) 'GUESS_TEMP'
#endif
write(6,*) ''
write(6,*) 'Reading initial conditions file'



open(unit=2,file=input,status='old')


!finds how many grid points are in the file
i=0
do 
  read(2,*,end=100) points
  i=i+1
enddo
100 continue

grand_ptot=i
write(6,*) 'Total elements: ',grand_ptot

allocate(prev_cooling(1:grand_ptot))
prev_cooling=-1.0D10

close(2)
write(6,*) ''

allocate(pdr(1:grand_ptot))
allocate(IDlist_pdr(1:grand_ptot))
allocate(IDlist_ion(1:grand_ptot))
allocate(IDlist_dark(1:grand_ptot))
pdr_ptot=0
ion_ptot=0
dark_ptot=0
maximum_density=0.0D0
minimum_density=1.0D10
open(unit=2,file=input,status='old')

do p=1,grand_ptot
    read(2,*) xpos,ypos,zpos,denst
    if (denst.le.rho_min) then
      ion_ptot = ion_ptot + 1
      pdr(p)%etype = 2 !IONIZED
      pdr(p)%x=xpos
      pdr(p)%y=ypos
      pdr(p)%z=zpos
      pdr(p)%rho=denst
      IDlist_ion(ion_ptot)=p
    endif
    if ((denst.gt.rho_min).AND.(denst.le.rho_max)) then
      pdr_ptot = pdr_ptot + 1
      pdr(p)%etype = 1 !PDR
      pdr(p)%x=xpos
      pdr(p)%y=ypos
      pdr(p)%z=zpos
      pdr(p)%rho=denst
      if (denst.gt.maximum_density) maximum_density=denst
      if (denst.lt.minimum_density) minimum_density=denst
      IDlist_pdr(pdr_ptot)=p
    endif
    if (denst.gt.rho_max) then
      dark_ptot = dark_ptot + 1
      pdr(p)%etype = 3 !DARK MOLECULAR
      pdr(p)%x=xpos
      pdr(p)%y=ypos
      pdr(p)%z=zpos
      pdr(p)%rho=denst
      IDlist_dark(dark_ptot)=p
    endif
enddo
write(6,*) 'PDR elements       = ',pdr_ptot
write(6,*) 'IONIZED elements   = ',ion_ptot
write(6,*) 'MOLECULAR elements = ',dark_ptot
write(6,*) 'Maximum PDR density = ',maximum_density
write(6,*) 'Minimum PDR density = ',minimum_density
write(6,*) 'Density used in DMR = ',2.0D0*rho_max


!reversing the last two rows
allocate(x_rev(1:pdr_ptot))
allocate(y_rev(1:pdr_ptot))
allocate(z_rev(1:pdr_ptot))
allocate(n_rev(1:pdr_ptot))
do pp=1,pdr_ptot-2
  p=IDlist_pdr(pp)
  x_rev(pp)=pdr(p)%x
  y_rev(pp)=pdr(p)%y
  z_rev(pp)=pdr(p)%z
  n_rev(pp)=pdr(p)%rho
enddo
x_rev(pdr_ptot-1)=pdr(IDlist_pdr(pdr_ptot))%x
y_rev(pdr_ptot-1)=pdr(IDlist_pdr(pdr_ptot))%y
z_rev(pdr_ptot-1)=pdr(IDlist_pdr(pdr_ptot))%z
n_rev(pdr_ptot-1)=pdr(IDlist_pdr(pdr_ptot))%rho

x_rev(pdr_ptot)=pdr(IDlist_pdr(pdr_ptot-1))%x
y_rev(pdr_ptot)=pdr(IDlist_pdr(pdr_ptot-1))%y
z_rev(pdr_ptot)=pdr(IDlist_pdr(pdr_ptot-1))%z
n_rev(pdr_ptot)=pdr(IDlist_pdr(pdr_ptot-1))%rho

!updates..
do pp=1,pdr_ptot
  p=IDlist_pdr(pp)
  pdr(p)%x=x_rev(pp)
  pdr(p)%y=y_rev(pp)
  pdr(p)%z=z_rev(pp)
  pdr(p)%rho=n_rev(pp)
enddo
!end reversing...


allocate(rra(0:pdr_ptot))
allocate(rrb(1:pdr_ptot))
allocate(pdrpoint(1:3,0:pdr_ptot)) !0 is for the ONE molecular element
close(2)
do pp=1,pdr_ptot
  p=IDlist_pdr(pp)
  rra(pp) = sqrt(pdr(p)%x**2+pdr(p)%y**2+pdr(p)%z**2)
  rrb(pp) = p
  pdrpoint(1,pp) = pdr(p)%x
  pdrpoint(2,pp) = pdr(p)%y
  pdrpoint(3,pp) = pdr(p)%z
enddo
!assigning x,y,z co-ordinates for the ONE molecular element
if (dark_ptot.gt.0) then
  pdrpoint(1,0) = pdr(IDlist_dark(1))%x
  pdrpoint(2,0) = pdr(IDlist_dark(1))%y
  pdrpoint(3,0) = pdr(IDlist_dark(1))%z
endif

call allocations

call readinput(C12Oinput,C12O_NLEV,C12O_NTEMP,C12O_ENERGIES,C12O_WEIGHTS,&
      &     C12O_A_COEFFS,C12O_B_COEFFS,C12O_FREQUENCIES,C12O_TEMPERATURES,&
      &     C12O_H,C12O_HP,C12O_EL,C12O_HE,C12O_H2,C12O_PH2,C12O_OH2)
call readinput(CIIinput,CII_NLEV,CII_NTEMP,CII_ENERGIES,CII_WEIGHTS,&
      &     CII_A_COEFFS,CII_B_COEFFS,CII_FREQUENCIES,CII_TEMPERATURES,&
      &     CII_H,CII_HP,CII_EL,CII_HE,CII_H2,CII_PH2,CII_OH2)
call readinput(CIinput,CI_NLEV,CI_NTEMP,CI_ENERGIES,CI_WEIGHTS,&
      &     CI_A_COEFFS,CI_B_COEFFS,CI_FREQUENCIES,CI_TEMPERATURES,&
      &     CI_H,CI_HP,CI_EL,CI_HE,CI_H2,CI_PH2,CI_OH2)
call readinput(OIinput,OI_NLEV,OI_NTEMP,OI_ENERGIES,OI_WEIGHTS,&
      &     OI_A_COEFFS,OI_B_COEFFS,OI_FREQUENCIES,OI_TEMPERATURES,&
      &     OI_H,OI_HP,OI_EL,OI_HE,OI_H2,OI_PH2,OI_OH2)

write(6,*) ''

call read_species(nspec, species, dummyabundance, mass)

call READ_RATES(NREAC,REACTANT,PRODUCT,ALPHA,BETA,GAMMA,rate,DUPLICATE,RTMIN,RTMAX)

dusttemperature = dust_temperature

#ifndef GUESS_TEMP
gastemperature = Tguess
previousgastemperature = Tguess
#ifdef THERMALBALANCE
Tlow = Tlow0
Thigh = Thigh0
#endif
#endif

do p=1,grand_ptot
  allocate(pdr(p)%abundance(1:nspec))
  do j=1,nspec
    if (species(j).eq.'H2'.OR.species(j).eq.'H'.OR.species(j).eq.'He') then
        pdr(p)%abundance(j) = dummyabundance(j)
    else
        pdr(p)%abundance(j) = dummyabundance(j)*metallicity
    endif
  enddo
enddo

!Sorting with increasing the distance from the (0,0,0)
!Direction of the PDR zone.
call heapsort(pdr_ptot,rrb,rra) !returns sorted the ID (rb) and the distance (ra)

!guess population density
allocate(total_heating(0:pdr_ptot))

!HEALPix stuff
nside=2**level
nrays=12*nside**2
ns_max=8192


allocate(vector(1:3))
allocate(vertex(1:3,1:4))
allocate(vectors(1:3,0:nrays-1))
!allocate(expanded(0:pdr_ptot));expanded=.false.
allocate(dobinarychop(0:pdr_ptot));dobinarychop=.false.
allocate(previouschange(0:pdr_ptot))

!HEALPix routine. 
call mk_xy2pix 


!creating HEALPix vectors.
write(6,*) 'Building HEALPix vectors...'
open(unit=77,file='HEALPix_vectors.dat',status='replace')
do i=1,nrays
  ipix=i-1 !ipix is the ID of a HEALPix ray. Runs with values 0:nrays-1
  call pix2vec_nest(nside,ipix,pix2x,pix2y,vector,vertex)
  vectors(1:3,ipix)=vector(1:3) !Store in memory
  write(77,'(3ES11.3,I7)') vectors(1:3,ipix),ipix
enddo

write(6,*) 'Done!';write(6,*) ''

write(6,*) 'allocating memory...'
do p=1,grand_ptot
    allocate(pdr(p)%CII_pop(1:CII_nlev))
    allocate(pdr(p)%CI_pop(1:CI_nlev))
    allocate(pdr(p)%OI_pop(1:OI_nlev))
    allocate(pdr(p)%C12O_pop(1:C12O_nlev))
enddo
do pp=1,pdr_ptot
    p=IDlist_pdr(pp)
    allocate(pdr(p)%epray(0:nrays-1))
    allocate(pdr(p)%epoint(1:3,0:nrays-1,0:maxpoints))
    allocate(pdr(p)%projected(0:nrays-1,0:maxpoints))
    allocate(pdr(p)%columndensity(0:nrays-1))
    allocate(pdr(p)%AV(0:nrays-1))
    allocate(pdr(p)%rad_surface(0:nrays-1))
    allocate(pdr(p)%CII_line(1:CII_nlev,1:CII_nlev))
    allocate(pdr(p)%CI_line(1:CI_nlev,1:CI_nlev))
    allocate(pdr(p)%OI_line(1:OI_nlev,1:OI_nlev))
    allocate(pdr(p)%C12O_line(1:C12O_nlev,1:C12O_nlev))
    allocate(pdr(p)%raytype(0:nrays-1))
!==========
    allocate(pdr(p)%CII_optdepth(1:CII_nlev,1:CII_nlev,0:nrays-1))
    allocate(pdr(p)%CI_optdepth(1:CI_nlev,1:CI_nlev,0:nrays-1))
    allocate(pdr(p)%OI_optdepth(1:OI_nlev,1:OI_nlev,0:nrays-1))
    allocate(pdr(p)%C12O_optdepth(1:C12O_nlev,1:C12O_nlev,0:nrays-1))
!==========
enddo
!Allocating for the ONE molecular element------
if (dark_ptot.gt.0) then
  allocate(pdr(IDlist_dark(1))%epray(0:nrays-1))
  allocate(pdr(IDlist_dark(1))%epoint(1:3,0:nrays-1,0:maxpoints))
  allocate(pdr(IDlist_dark(1))%projected(0:nrays-1,0:maxpoints))
  allocate(pdr(IDlist_dark(1))%columndensity(0:nrays-1))
  allocate(pdr(IDlist_dark(1))%AV(0:nrays-1))
  allocate(pdr(IDlist_dark(1))%rad_surface(0:nrays-1))
  allocate(pdr(IDlist_dark(1))%CII_line(1:CII_nlev,1:CII_nlev))
  allocate(pdr(IDlist_dark(1))%CI_line(1:CI_nlev,1:CI_nlev))
  allocate(pdr(IDlist_dark(1))%OI_line(1:OI_nlev,1:OI_nlev))
  allocate(pdr(IDlist_dark(1))%C12O_line(1:C12O_nlev,1:C12O_nlev))
  allocate(pdr(p)%raytype(0:nrays-1))
endif
!----------------------------------------------
allocate(CII_solution(0:pdr_ptot,1:CII_nlev))
allocate(CI_solution(0:pdr_ptot,1:CI_nlev))
allocate(OI_solution(0:pdr_ptot,1:OI_nlev))
allocate(C12O_solution(0:pdr_ptot,1:C12O_nlev))
write(6,*) 'Memory OK';write(6,*) ''


write(6,*) 'Building evaluation points...'
call evaluation_points

#ifndef PSEUDO_1D
write(6,*) 'subtracting 1 evaluation point in raytype(j)=-2'
do pp=1,pdr_ptot
  p=IDlist_pdr(pp)
  do j=0,nrays-1
    if (pdr(p)%epray(j).eq.0) cycle
    if (pdr(p)%raytype(j).eq.-2) pdr(p)%epray(j)=pdr(p)%epray(j)-1
  enddo
enddo
#endif

call cpu_time(time_evalpoints)
#ifdef OPENMP
write(6,*) 'Time [PARALLEL] = ',time_evalpoints/real(CPUs),' seconds'
#else
write(6,*) 'Time = ',time_evalpoints,' seconds'
#endif

maxpoints = 0
do pp=1,pdr_ptot
   p=IDlist_pdr(pp)
   newmaxpoints = maxval(pdr(p)%epray)
   if (newmaxpoints.gt.maxpoints) maxpoints = newmaxpoints
enddo

!Dark Molecular element-------
if (dark_ptot.gt.0) then
p=IDlist_dark(1)
newmaxpoints = maxval(pdr(p)%epray)
if (newmaxpoints.gt.maxpoints) maxpoints = newmaxpoints
endif
!-----------------------------
write(6,*) '';write(6,*) 'new maxpoints = ',maxpoints

write(6,*) 'Calculating UV field ...'

call calc_UVfield

#ifdef DUST2
call calculate_dust_temperatures
!updates the dust temperature for the PDR particles
do pp=1,pdr_ptot
   p=IDlist_pdr(pp)
   dusttemperature(pp)=pdr(p)%dust_t
enddo
#endif

#ifdef THERMALBALANCE
allocate(converged(0:pdr_ptot))
allocate(doleveltmin(0:pdr_ptot))
doleveltmin=.false.
converged=.false.
level_conv=.false.
first_time=.true.
#endif
allocate(column(0:pdr_ptot))
write(6,*) 'Calculating column densities...'
referee=0
call calc_columndens
referee=1

start_time = 0.0D0

if (dark_ptot.gt.0) call DARK_MOLECULAR_REGION

ITERATION = 0
!======== LTE LEVEL POPULATIONS ============
  write(6,*) ''; write(6,*) 'Calculating LTE level populations...' 
  call cpu_time(t3b)


relative_abundance_tolerance = 1.0D-8
absolute_abundance_tolerance = 1.0D-30
#ifdef THERMALBALANCE
dobinarychop=.false.
#endif

 allocate(dummy_rate(1:nreac,1:pdr_ptot))
 allocate(dummy_abundance(1:nspec,1:pdr_ptot))
 allocate(dummy_density(1:pdr_ptot))
 allocate(dummy_temperature(1:pdr_ptot))



call cpu_time(tt1)
DO II=1,CHEMITERATIONS
 write(6,*) 'Chemical iteration ',ii

#ifdef OPENMP
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE (pp,p,rate)&
!$OMP PRIVATE(NRGR,NRH2,NRHD,NRCO,NRCI,NRSI)
#endif
  do pp=1,pdr_ptot
    p=IDlist_pdr(pp)
!    if (allocated(rate)) deallocate(rate); allocate(rate(1:nreac))
    CALL CALCULATE_REACTION_RATES(gastemperature(pp),dusttemperature(pp),nrays,pdr(p)%rad_surface(0:nrays-1),&
          &pdr(p)%AV(0:nrays-1),column(pp)%columndens_point(0:nrays-1,1:nspec),&
          &nreac, reactant, product, alpha, beta, gamma, rate, rtmin, rtmax, duplicate, nspec,&
          &NRGR,NRH2,NRHD,NRCO,NRCI,NRSI)
      DUMMY_RATE(:,pp) = rate
      DUMMY_ABUNDANCE(:,pp) = pdr(p)%abundance
      DUMMY_DENSITY(pp) = pdr(p)%rho
      DUMMY_TEMPERATURE(pp) = gastemperature(pp)
   enddo
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif
  CALL CALCULATE_ABUNDANCES(DUMMY_ABUNDANCE,DUMMY_RATE,DUMMY_DENSITY,DUMMY_TEMPERATURE,pdr_ptot,NSPEC,NREAC)

#ifdef OPENMP
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE (pp,p)
#endif
  DO pp=1,pdr_ptot
   p=IDlist_pdr(pp)
   pdr(p)%abundance = dummy_abundance(:,pp)
  enddo
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

   call calc_columndens

ENDDO
 deallocate(dummy_abundance)
 deallocate(dummy_density)
 deallocate(dummy_temperature)



#ifdef OPENMP
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(pp,p,CII_Z_function)&
!$OMP PRIVATE(CI_Z_function,OI_Z_function,C12O_Z_function)
#endif
 do pp=1,pdr_ptot
    p=IDlist_pdr(pp)
    !Calculate the partition functions
    CALL CALCULATE_PARTITION_FUNCTION(CII_Z_FUNCTION,CII_NLEV,CII_ENERGIES,CII_WEIGHTS,GASTEMPERATURE(PP))
    CALL CALCULATE_PARTITION_FUNCTION(CI_Z_FUNCTION,CI_NLEV,CI_ENERGIES,CI_WEIGHTS,GASTEMPERATURE(PP))
    CALL CALCULATE_PARTITION_FUNCTION(OI_Z_FUNCTION,OI_NLEV,OI_ENERGIES,OI_WEIGHTS,GASTEMPERATURE(PP))
    CALL CALCULATE_PARTITION_FUNCTION(C12O_Z_FUNCTION,C12O_NLEV,C12O_ENERGIES,C12O_WEIGHTS,GASTEMPERATURE(PP))
#ifndef GUESS_TEMP
    IF (pp.EQ.1) THEN
       WRITE(6,*) ''
       WRITE(6,*) 'Z(CII)  = ',CII_Z_FUNCTION
       WRITE(6,*) 'Z(CI)   = ',CI_Z_FUNCTION
       WRITE(6,*) 'Z(OI)   = ',OI_Z_FUNCTION
       WRITE(6,*) 'Z(C12O) = ',C12O_Z_FUNCTION
       WRITE(6,*) ''
    ENDIF
#endif
    !
    ! Calculate the LTE level populations
    CALL CALCULATE_LTE_POPULATIONS(CII_NLEV,PDR(P)%CII_POP,CII_ENERGIES,&
        &CII_WEIGHTS,CII_Z_FUNCTION,PDR(P)%ABUNDANCE(NCx)*PDR(P)%RHO,GASTEMPERATURE(PP))
    CALL CALCULATE_LTE_POPULATIONS(CI_NLEV, PDR(P)%CI_POP, CI_ENERGIES, &
        &CI_WEIGHTS, CI_Z_FUNCTION, PDR(P)%ABUNDANCE(NC)*PDR(P)%RHO, GASTEMPERATURE(PP))
    CALL CALCULATE_LTE_POPULATIONS(OI_NLEV, PDR(P)%OI_POP, OI_ENERGIES, &
        &OI_WEIGHTS, OI_Z_FUNCTION, PDR(P)%ABUNDANCE(NO)*PDR(P)%RHO, GASTEMPERATURE(PP))
    CALL CALCULATE_LTE_POPULATIONS(C12O_NLEV,PDR(P)%C12O_POP,C12O_ENERGIES,&
        &C12O_WEIGHTS,C12O_Z_FUNCTION,PDR(P)%ABUNDANCE(NCO)*PDR(P)%RHO,GASTEMPERATURE(PP))
 enddo !ii=1,itot
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

call cpu_time(tt2)
#ifdef OPENMP
write(6,*) 'No. CPUs = ',CPUs
write(6,*) 'Time for LTE level populations = ',(tt2-tt1)/real(CPUs,kind=dp),' seconds'
#else
write(6,*) 'Time for LTE level populations (SERIAL)= ',(tt2-tt1),' seconds'
#endif

call cpu_time(t3)
write(6,*) 'Total time = ',t3,' seconds'

write(6,*) ''

allocate(level_converged(0:pdr_ptot))
level_converged=.false.
allocate(CII_conv(0:pdr_ptot))
allocate(CI_conv(0:pdr_ptot))
allocate(OI_conv(0:pdr_ptot))
allocate(C12O_conv(0:pdr_ptot))
CII_conv=.false.
CI_conv=.false.
OI_conv=.false.
C12O_conv=.false.
allocate(CII_RELCH(0:pdr_ptot,1:CII_NLEV))
allocate(CI_RELCH(0:pdr_ptot,1:CI_NLEV))
allocate(OI_RELCH(0:pdr_ptot,1:OI_NLEV))
allocate(C12O_RELCH(0:pdr_ptot,1:C12O_NLEV))

write(6,*) 'Begin iterations...'

levpop_iteration=0
cii_percentage=0
ci_percentage=0
oi_percentage=0
DO ITERATION=1,ITERTOT
  write_output=.false.
  write(6,*) ''
  write(6,*) 'Iteration ',iteration
  call cpu_time(t3b)
!  show_balance=.true.
  levpop_iteration=levpop_iteration+1
  write(6,*) 'Level population iteration ',levpop_iteration
    IF (iteration.gt.1.and.levpop_iteration.eq.1) THEN
 allocate(dummy_abundance(1:nspec,1:pdr_ptot))
 allocate(dummy_density(1:pdr_ptot))
 allocate(dummy_temperature(1:pdr_ptot))
      DO II=1,3
         write(6,*) 'Chemical iteration ',ii
#ifdef OPENMP
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE (pp,p,rate)&
!$OMP PRIVATE(NRGR,NRH2,NRHD,NRCO,NRCI,NRSI)
#endif
        do pp=1,pdr_ptot
          p=IDlist_pdr(pp)
!          if (allocated(rate)) deallocate(rate); allocate(rate(1:nreac))
          CALL CALCULATE_REACTION_RATES(gastemperature(pp),dusttemperature(pp),nrays,pdr(p)%rad_surface(0:nrays-1),&
             &pdr(p)%AV(0:nrays-1),column(pp)%columndens_point(0:nrays-1,1:nspec),&
             &nreac, reactant, product, alpha, beta, gamma, rate, rtmin, rtmax, duplicate, nspec,&
             &NRGR,NRH2,NRHD,NRCO,NRCI,NRSI)
          DUMMY_RATE(:,pp) = rate
          DUMMY_ABUNDANCE(:,pp) = pdr(p)%abundance
          DUMMY_DENSITY(pp) = pdr(p)%rho
          DUMMY_TEMPERATURE(pp) = gastemperature(pp)
       enddo
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif
       CALL CALCULATE_ABUNDANCES(DUMMY_ABUNDANCE,DUMMY_RATE,DUMMY_DENSITY,DUMMY_TEMPERATURE,pdr_ptot,NSPEC,NREAC)
#ifdef OPENMP
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE (pp,p)
#endif
       DO pp=1,pdr_ptot
       	  if (converged(pp)) cycle
          p=IDlist_pdr(pp)
          pdr(p)%abundance = dummy_abundance(:,pp)
       enddo
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif
       call calc_columndens
     ENDDO !CHEMICAL ITERATION
 deallocate(dummy_abundance)
 deallocate(dummy_density)
 deallocate(dummy_temperature)


!Perform LTE...
#ifdef OPENMP
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(p,pp) &
!$OMP PRIVATE(CII_Z_FUNCTION,CI_Z_FUNCTION) &
!$OMP PRIVATE(OI_Z_FUNCTION,C12O_Z_FUNCTION)
#endif
    do pp=1,pdr_ptot
       p=IDlist_pdr(pp)
#ifdef THERMALBALANCE
        if (level_converged(pp).or.converged(pp)) cycle
#else
        if (level_converged(pp)) cycle
#endif
        !Calculate the partition functions
        CALL CALCULATE_PARTITION_FUNCTION(CII_Z_FUNCTION,CII_NLEV,CII_ENERGIES,CII_WEIGHTS,GASTEMPERATURE(PP))
        CALL CALCULATE_PARTITION_FUNCTION(CI_Z_FUNCTION,CI_NLEV,CI_ENERGIES,CI_WEIGHTS,GASTEMPERATURE(PP))
        CALL CALCULATE_PARTITION_FUNCTION(OI_Z_FUNCTION,OI_NLEV,OI_ENERGIES,OI_WEIGHTS,GASTEMPERATURE(PP))
        CALL CALCULATE_PARTITION_FUNCTION(C12O_Z_FUNCTION,C12O_NLEV,C12O_ENERGIES,C12O_WEIGHTS,GASTEMPERATURE(PP))
        !Calculate the LTE level populations
        CALL CALCULATE_LTE_POPULATIONS(CII_NLEV,PDR(P)%CII_POP,CII_ENERGIES,CII_WEIGHTS,&
             &CII_Z_FUNCTION,PDR(P)%ABUNDANCE(NCx)*PDR(P)%RHO,GASTEMPERATURE(PP))
        CALL CALCULATE_LTE_POPULATIONS(CI_NLEV, PDR(P)%CI_POP, CI_ENERGIES, CI_WEIGHTS, &
             &CI_Z_FUNCTION, PDR(P)%ABUNDANCE(NC)*PDR(P)%RHO, GASTEMPERATURE(PP))
        CALL CALCULATE_LTE_POPULATIONS(OI_NLEV, PDR(P)%OI_POP, OI_ENERGIES, OI_WEIGHTS, &
             &OI_Z_FUNCTION, PDR(P)%ABUNDANCE(NO)*PDR(P)%RHO, GASTEMPERATURE(PP))
        CALL CALCULATE_LTE_POPULATIONS(C12O_NLEV,PDR(P)%C12O_POP,C12O_ENERGIES,C12O_WEIGHTS,&
             &C12O_Z_FUNCTION,PDR(P)%ABUNDANCE(NCO)*PDR(P)%RHO,GASTEMPERATURE(PP))

     enddo !second particle loop
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

    ELSE
#ifdef OPENMP
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(pp,p) &
!$OMP PRIVATE(j,i,ilevel) &
!$OMP PRIVATE(dummystep,dummycoef) &
!$OMP PRIVATE(CII_C_COEFFS,CI_C_COEFFS,OI_C_COEFFS,C12O_C_COEFFS) &
!$OMP PRIVATE(transition_CII,transition_CI,transition_OI,transition_C12O) &
!$OMP PRIVATE(dummyarray_CII,dummyarray_CI,dummyarray_OI,dummyarray_C12O) &
!$OMP PRIVATE(dummyarray_CII_tau,dummyarray_CI_tau,dummyarray_OI_tau,dummyarray_C12O_tau)&
!$OMP PRIVATE(dummyarray_CII_beta,dummyarray_CI_beta,dummyarray_OI_beta,dummyarray_C12O_beta) &
!$OMP PRIVATE(CIIsolution,CIsolution,OIsolution,C12Osolution) &
!$OMP PRIVATE(CIIevalpop,CIevalpop,OIevalpop,C12Oevalpop)
#endif
    do pp=1,pdr_ptot
       p=IDlist_pdr(pp)
!-----------------------------------------------
#ifdef THERMALBALANCE
    if (level_converged(pp).or.converged(pp)) cycle
#else
    if (level_converged(pp)) cycle
#endif
!-----------------------------------------------
    if (level_converged(pp)) cycle
          pdr(p)%projected(:,0) = p

ALLOCATE(CII_C_COEFFS(1:CII_NLEV,1:CII_NLEV))
ALLOCATE(CI_C_COEFFS(1:CI_NLEV,1:CI_NLEV))
ALLOCATE(OI_C_COEFFS(1:OI_NLEV,1:OI_NLEV))
ALLOCATE(C12O_C_COEFFS(1:C12O_NLEV,1:C12O_NLEV))
allocate(transition_CII(1:CII_nlev,1:CII_nlev))
allocate(transition_CI(1:CI_nlev,1:CI_nlev))
allocate(transition_OI(1:OI_nlev,1:OI_nlev))
allocate(transition_C12O(1:C12O_nlev,1:C12O_nlev))
allocate(dummyarray_CII(1:CII_nlev,1:CII_nlev))
allocate(dummyarray_CI(1:CI_nlev,1:CI_nlev))
allocate(dummyarray_OI(1:OI_nlev,1:OI_nlev))
allocate(dummyarray_C12O(1:C12O_nlev,1:C12O_nlev))

!========
allocate(dummyarray_CII_tau(1:CII_nlev,1:CII_nlev,0:nrays-1))
allocate(dummyarray_CI_tau(1:CI_nlev,1:CI_nlev,0:nrays-1))
allocate(dummyarray_OI_tau(1:OI_nlev,1:OI_nlev,0:nrays-1))
allocate(dummyarray_C12O_tau(1:C12O_nlev,1:C12O_nlev,0:nrays-1))
allocate(dummyarray_C12O_beta(1:C12O_nlev,1:C12O_nlev,0:nrays-1))
allocate(dummyarray_CII_beta(1:CII_nlev,1:CII_nlev,0:nrays-1))
allocate(dummyarray_CI_beta(1:CI_nlev,1:CI_nlev,0:nrays-1))
allocate(dummyarray_OI_beta(1:OI_nlev,1:OI_nlev,0:nrays-1))
allocate(dummystep(1:7))
allocate(dummycoef(1:CII_nlev,1:CII_Nlev,1:7))
!=======
allocate(CIIsolution(1:CII_nlev))
allocate(CIsolution(1:CI_nlev))
allocate(OIsolution(1:OI_nlev))
allocate(C12Osolution(1:C12O_nlev))
allocate(CIIevalpop(0:nrays-1,0:maxpoints,1:CII_nlev))
allocate(CIevalpop(0:nrays-1,0:maxpoints,1:CI_nlev))
allocate(OIevalpop(0:nrays-1,0:maxpoints,1:OI_nlev))
allocate(C12Oevalpop(0:nrays-1,0:maxpoints,1:C12O_nlev))
CIIevalpop=0.0D0; CIevalpop=0.0D0; OIevalpop=0.0D0; C12Oevalpop=0.0D0
   
       ! Specify the evaluation points along each ray from the current pdrpoint
       do j=0,nrays-1
          do i=0,pdr(p)%epray(j)
             do ilevel=1,CII_nlev 
                CIIevalpop(j,i,ilevel)=pdr(int(pdr(p)%projected(j,i)))%CII_pop(ilevel)
             enddo
             do ilevel=1,CI_nlev 
                CIevalpop(j,i,ilevel)=pdr(int(pdr(p)%projected(j,i)))%CI_pop(ilevel)
             enddo
             do ilevel=1,OI_nlev 
                OIevalpop(j,i,ilevel)=pdr(int(pdr(p)%projected(j,i)))%OI_pop(ilevel)
             enddo
             do ilevel=1,C12O_nlev 
                C12Oevalpop(j,i,ilevel)=pdr(int(pdr(p)%projected(j,i)))%C12O_pop(ilevel)
             enddo
          enddo !pdr(p)%epray(j)
       enddo
       !
       ! Use the LVG (escape probability) method to determine the transition matrices and solve for the level populations
! CII calculations -------------------------------------------------
       ! Calculate the collisional rate coefficients
       CALL FIND_CCOEFF(CII_NTEMP,CII_NLEV,gastemperature(pp),CII_TEMPERATURES,&
          & CII_H,CII_HP,CII_EL,CII_HE,CII_H2,CII_PH2,CII_OH2,&
          & CII_C_COEFFS,pdr(p)%abundance(NH)*pdr(p)%rho,pdr(p)%abundance(NPROTON)*pdr(p)%rho, &
          & pdr(p)%abundance(NELECT)*pdr(p)%rho, pdr(p)%abundance(NHE)*pdr(p)%rho,pdr(p)%abundance(NH2)*pdr(p)%rho,&
          & 1)
       call escape_probability(transition_CII, dusttemperature(pp), nrays, CII_nlev, &
              &CII_A_COEFFS, CII_B_COEFFS, CII_C_COEFFS, &
              &CII_frequencies, CIIevalpop, maxpoints, &
              &gastemperature(pp), v_turb, pdr(p)%epray, pdr(p)%CII_pop, &
              &pdr(p)%epoint, CII_weights,CII_cool(pp),dummyarray_CII,&
              &dummyarray_CII_tau,1,pdr(p)%rho,metallicity,dummyarray_CII_beta)
       pdr(p)%CII_line=dummyarray_CII
       pdr(p)%CII_optdepth=dummyarray_CII_tau
       call solvlevpop(CII_nlev,transition_CII,pdr(p)%abundance(NCx)*pdr(p)%rho,CIIsolution)!,1)
       CII_solution(pp,:)=CIIsolution
#ifdef CO_FIX
       if (levpop_iteration.ge.120) then
          CII_solution(pp,:)=pdr(p)%CII_pop
       else if (levpop_iteration.ge.75) then
          CII_solution(pp,:)=0.5*(CII_solution(pp,:) + pdr(p)%CII_pop)
       endif
#endif
!-------------------------------------------------------------------

! CI calculations --------------------------------------------------
       ! Calculate the collisional rate coefficients
       CALL FIND_CCOEFF(CI_NTEMP,CI_NLEV,gastemperature(pp),CI_TEMPERATURES,&
          & CI_H,CI_HP,CI_EL,CI_HE,CI_H2,CI_PH2,CI_OH2,&
          & CI_C_COEFFS,pdr(p)%abundance(NH)*pdr(p)%rho,pdr(p)%abundance(NPROTOn)*pdr(p)%rho, &
          & pdr(p)%abundance(NELECT)*pdr(p)%rho, pdr(p)%abundance(NHE)*pdr(p)%rho,pdr(p)%abundance(NH2)*pdr(p)%rho,2)
       call escape_probability(transition_CI, dusttemperature(pp), nrays, CI_nlev, &
              &CI_A_COEFFS, CI_B_COEFFS, CI_C_COEFFS, &
              &CI_frequencies, CIevalpop, maxpoints, &
              &gastemperature(pp), v_turb, pdr(p)%epray, pdr(p)%CI_pop, &
              &pdr(p)%epoint,CI_weights,CI_cool(pp),dummyarray_CI,dummyarray_CI_tau,2,pdr(p)%rho,metallicity,dummyarray_CI_beta)
       pdr(p)%CI_line=dummyarray_CI
       pdr(p)%CI_optdepth=dummyarray_CI_tau
       call solvlevpop(CI_nlev,transition_CI,pdr(p)%abundance(NC)*pdr(p)%rho,CIsolution)!,2)
       CI_solution(pp,:)=CIsolution
#ifdef CO_FIX
       if (levpop_iteration.ge.120) then
          CI_solution(pp,:)=pdr(p)%CI_pop
       else if (levpop_iteration.ge.75) then
          CI_solution(pp,:)=0.5*(CI_solution(pp,:) + pdr(p)%CI_pop)
       endif
#endif
!-------------------------------------------------------------------

! OI calculations --------------------------------------------------
       ! Calculate the collisional rate coefficients
       CALL FIND_CCOEFF(OI_NTEMP,OI_NLEV,gastemperature(pp),OI_TEMPERATURES,&
          & OI_H,OI_HP,OI_EL,OI_HE,OI_H2,OI_PH2,OI_OH2,&
          & OI_C_COEFFS,pdr(p)%abundance(NH)*pdr(p)%rho,pdr(p)%abundance(NPROTON)*pdr(p)%rho, &
          & pdr(p)%abundance(NELECT)*pdr(p)%rho, pdr(p)%abundance(NHE)*pdr(p)%rho,pdr(p)%abundance(NH2)*pdr(p)%rho,3)
       call escape_probability(transition_OI, dusttemperature(pp), nrays, OI_nlev, &
              &OI_A_COEFFS, OI_B_COEFFS, OI_C_COEFFS, &
              &OI_frequencies, OIevalpop, maxpoints, &
              &gastemperature(pp), v_turb, pdr(p)%epray, pdr(p)%OI_pop, &
              &pdr(p)%epoint,OI_weights,OI_cool(pp),dummyarray_OI,dummyarray_OI_tau,3,pdr(p)%rho,metallicity,dummyarray_OI_beta)
       pdr(p)%OI_line=dummyarray_OI
       pdr(p)%OI_optdepth=dummyarray_OI_tau
       call solvlevpop(OI_nlev,transition_OI,pdr(p)%abundance(NO)*pdr(p)%rho,OIsolution)!,3)
       OI_solution(pp,:)=OIsolution
#ifdef CO_FIX
       if (levpop_iteration.ge.120) then
          OI_solution(pp,:)=pdr(p)%OI_pop
       else if (levpop_iteration.ge.75) then
          OI_solution(pp,:)=0.5*(OI_solution(pp,:) + pdr(p)%OI_pop)
       endif
#endif
!-------------------------------------------------------------------

! C12O calculations ------------------------------------------------
       ! Calculate the collisional rate coefficients
       CALL FIND_CCOEFF(C12O_NTEMP,C12O_NLEV,gastemperature(pp),C12O_TEMPERATURES,&
          & C12O_H,C12O_HP,C12O_EL,C12O_HE,C12O_H2,C12O_PH2,C12O_OH2,&
          & C12O_C_COEFFS,pdr(p)%abundance(NH)*pdr(p)%rho,pdr(p)%abundance(NPROTON)*pdr(p)%rho, &
          & pdr(p)%abundance(NELECT)*pdr(p)%rho, pdr(p)%abundance(NHE)*pdr(p)%rho,pdr(p)%abundance(NH2)*pdr(p)%rho,4)
       call escape_probability(transition_C12O, dusttemperature(pp), nrays, C12O_nlev, &
              &C12O_A_COEFFS, C12O_B_COEFFS, C12O_C_COEFFS, &
              &C12O_frequencies, C12Oevalpop, maxpoints, &
              &gastemperature(pp), v_turb, pdr(p)%epray, pdr(p)%C12O_pop, &
              &pdr(p)%epoint,C12O_weights,C12O_cool(pp),dummyarray_C12O,&
              &dummyarray_C12O_tau,4,pdr(p)%rho,metallicity,dummyarray_C12O_beta)
       pdr(p)%C12O_line=dummyarray_C12O
       pdr(p)%C12O_optdepth=dummyarray_C12O_tau
       call solvlevpop(C12O_nlev,transition_C12O,pdr(p)%abundance(NCO)*pdr(p)%rho,C12Osolution)!,4)
       C12O_solution(pp,:)=C12Osolution
!-------------------------------------------------------------------

#ifdef CO_FIX
       ! If the level populations are oscillating and not converging, try to suppress the oscillations
       ! by taking the average of the level populations from this iteration and the previous iteration
       if (CII_percentage.eq.100.and.CI_percentage.eq.100.and.&
           &OI_percentage.eq.100) then
          if (levpop_iteration.ge.120) then
             C12O_solution(pp,:)=pdr(p)%C12O_pop
          else if (levpop_iteration.ge.75) then
             C12O_solution(pp,:)=0.5*(C12O_solution(pp,:) + pdr(p)%C12O_pop)
          endif
       endif
#endif

deALLOCATE(CII_C_COEFFS)
deALLOCATE(CI_C_COEFFS)
deALLOCATE(OI_C_COEFFS)
deALLOCATE(C12O_C_COEFFS)
deallocate(transition_CII)
deallocate(transition_CI)
deallocate(transition_OI)
deallocate(transition_C12O)
deallocate(dummyarray_CII)
deallocate(dummyarray_CI)
deallocate(dummyarray_OI)
deallocate(dummyarray_C12O)
!============
deallocate(dummyarray_CII_tau)
deallocate(dummyarray_CI_tau)
deallocate(dummyarray_OI_tau)
deallocate(dummyarray_C12O_tau)
deallocate(dummyarray_C12O_beta)
deallocate(dummyarray_CII_beta)
deallocate(dummyarray_CI_beta)
deallocate(dummyarray_OI_beta)
deallocate(dummystep)
deallocate(dummycoef)
!============
deallocate(CIIsolution)
deallocate(CIsolution)
deallocate(OIsolution)
deallocate(C12Osolution)
deallocate(CIIevalpop)
deallocate(CIevalpop)
deallocate(OIevalpop)
deallocate(C12Oevalpop)


enddo !particles
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif
total_cooling_rate=CII_cool+CI_cool+OI_cool+C12O_cool

!================
   endif ! LTE or LVG calculation of the level populations
!================
#ifdef OPENMP
#ifdef THERMALBALANCE
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(p,pp,allheating,dummytemperature,rate) &
#else
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(p,pp,allheating,rate)&
#endif
!$OMP PRIVATE(NRGR,NRH2,NRHD,NRCO,NRCI,NRSI)
#endif
do pp=1,pdr_ptot
   p=IDlist_pdr(pp)
#ifdef THERMALBALANCE
! Skip this pdrpoint if the temperature has already converged
   if (converged(pp)) cycle
#endif
!   if (allocated(rate)) deallocate(rate); allocate(rate(1:nreac))
   if (allocated(allheating)) deallocate(allheating); allocate(allheating(1:12))
    CALL CALCULATE_REACTION_RATES(gastemperature(pp),dusttemperature(pp),nrays,pdr(p)%rad_surface(0:nrays-1),&
          &pdr(p)%AV(0:nrays-1),column(pp)%columndens_point(0:nrays-1,1:nspec),&
          &nreac, reactant, product, alpha, beta, gamma, rate, rtmin, rtmax, duplicate, nspec,&
          &NRGR,NRH2,NRHD,NRCO,NRCI,NRSI)
    call calc_heating(pdr(p)%rho,gastemperature(pp),dusttemperature(pp),pdr(p)%UVfield, &
          &v_turb,nspec,pdr(p)%abundance(:),nreac,rate,allheating,&
          &NRGR,NRH2,NRHD,NRCO,NRCI,NRSI)

   all_heating(pp,:)=allheating

#ifdef THERMALBALANCE
! Calculate the difference between the total heating and total cooling rates (Fmean)
! and the absolute value of the relative difference between the two rates (Fratio)
   Fmean(pp) = all_heating(pp,12) - total_cooling_rate(pp)
   Fratio(pp) = 2.0D0*abs(Fmean(pp))/abs(all_heating(pp,12) + total_cooling_rate(pp))

! Store the current temperature in a dummy variable
! Do not start testing the thermal balance until enough iterations have passed for the level populations to begin to converge...
   if (level_conv.and.first_time) then
   dummytemperature = gastemperature(pp)
! Determine the temperature bracket to begin searching within...
      if (Fmean(pp).eq.0) then ! Handle the (very rare) case when the initial guess temperature is the correct value
         Tlow(pp) = gastemperature(pp)  ! Update the value of Tlow
         Thigh(pp) = gastemperature(pp) ! Update the value of Thigh
      else if (Fmean(pp).gt.0) then !---> HEATING
         Tlow(pp) = gastemperature(pp)  ! Update the value of Tlow
         gastemperature(pp) = 1.3D0*gastemperature(pp) !increase 30%
         previouschange(pp) = "H" !we increase
      else if (Fmean(pp).lt.0) then !---> COOLING
         Thigh(pp) = gastemperature(pp) ! Update the value of Thigh
         gastemperature(pp) = 0.7D0*gastemperature(pp) !decrease 30%
         if (gastemperature(pp).lt.Tmin) gastemperature(pp)=Tmin
         previouschange(pp) = "C" !we decrease
      endif
   previousgastemperature(pp) = dummytemperature
   if (previousgastemperature(pp).lt.Tmin) previousgastemperature(pp)=Tmin

   else if (level_conv.and..not.first_time) then
   dummytemperature = gastemperature(pp)
! Check for convergence in both the heating-cooling imbalance and the temperature difference between iterations
      if (Fratio(pp).le.Fcrit) converged(pp) = .true.

if (.not.dobinarychop(pp)) then
   !if we *still* need to heat, increase by 30% 
   if (Fmean(pp).gt.0.and.previouschange(pp).eq."H") then
       Tlow(pp) = gastemperature(pp)
       gastemperature(pp) = 1.3D0*gastemperature(pp)
       Thigh(pp) = gastemperature(pp)
       previouschange(pp) = "H"

   endif
   !if we *still* need to cool, decrease by 30%
   if (Fmean(pp).lt.0.and.previouschange(pp).eq."C") then
       Thigh(pp) = gastemperature(pp)
       gastemperature(pp) = 0.7D0*gastemperature(pp)
       Tlow(pp) = gastemperature(pp)
       previouschange(pp) = "C"
       if (gastemperature(pp).lt.Tmin) then 
           gastemperature(pp)=Tmin
           Tlow(pp)=Tmin
           Thigh(pp)=Tmin
       endif

   endif
   !For all other cases do binary chop and flag the process as .true.
   !Needs heating but previously it was decreased by 30%. 
   if (Fmean(pp).gt.0.and.previouschange(pp).eq."C") then
     gastemperature(pp) = (Thigh(pp) + Tlow(pp))/2.0D0
     dobinarychop(pp)=.true.  !from now on

   endif
   !Needs cooling but previously it was increased by 30%. Now do binary chop
   if (Fmean(pp).lt.0.and.previouschange(pp).eq."H") then 
           gastemperature(pp) = (Thigh(pp) + Tlow(pp))/2.0D0
           dobinarychop(pp)=.true.  !from now on
   endif


else !from now on only binary chop (we found the min-max by the 30% increase/decrease)
   if (Fmean(pp).gt.0) then 
         Tlow(pp) = gastemperature(pp)
         gastemperature(pp) = (gastemperature(pp) + Thigh(pp)) / 2.0D0
   endif
   if (Fmean(pp).lt.0) then
         Thigh(pp) = gastemperature(pp)
         gastemperature(pp) = (gastemperature(pp) + Tlow(pp)) / 2.0D0
   endif
endif !dobinarychop

#ifdef TEMP_FIX
      ! If the temperatures are converging in a value that thermal balance can't be reached
      ! double the high temperature and half the low one. If the temperatures are still not converging, force convergence.
      if ((abs(gastemperature(pp)-previousgastemperature(pp)).le.Tdiff).and.(Fratio(pp).gt.Fcrit)) then
         converged(pp)=.true.
   !      if (expanded(pp)) converged(pp) = .true.
   !      if (Fmean(pp).gt.0.and..not.converged(pp)) then 
   !        Tlow(pp) = gastemperature(pp)
   !        gastemperature(pp) = 4.0D0*gastemperature(pp)
   !      endif
   !      if (Fmean(pp).lt.0.and..not.converged(pp)) then 
   !        Thigh(pp) = gastemperature(pp)
   !        gastemperature(pp) = 0.25D0*gastemperature(pp)
   !      endif
   !      expanded(pp) = .true.
      endif
#endif

! Replace the previous temperature with the current value
      previousgastemperature(pp) = dummytemperature

      if ((dummytemperature.lt.Tmin).and.(Fmean(pp).lt.0)) converged(pp)=.true.
      if ((dummytemperature.gt.Tmax).and.(Fmean(pp).gt.0)) converged(pp)=.true.
 
      if (converged(pp)) then
         if (dummytemperature.lt.Tmin) then 
             previousgastemperature(pp) = Tmin
             gastemperature(pp) = Tmin
             if (doleveltmin(pp)) then 
                converged(pp)=.true.
             else
                converged(pp)=.false.
                level_converged(pp)=.false.
                dummytemperature=Tmin
             endif
             doleveltmin(pp)=.true.
         endif
         if (dummytemperature.gt.Tmax) then
              previousgastemperature(pp) = Tmax
              gastemperature(pp) = Tmax
         endif
      endif

   endif !level_conv first_time
!endif THERMALBALANCE
#endif 

enddo !ii=1,itot
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif


do pp=1,pdr_ptot
  prev_cooling(pp)=total_cooling_rate(pp)
enddo

 write(6,*) 'Checking for convergence...'
#ifdef THERMALBALANCE
   if (level_conv.and.first_time) first_time=.false.
   i=0
   do p=1,pdr_ptot
      if (converged(p)) i=i+1
   enddo
   if (level_conv) then
      write(6,*) 'Resetting [level_conv=.false.]'
      level_conv=.false.
      write_output=.true.
   endif
   thermal_percentage = 100.D0*real(i,kind=dp)/real(pdr_ptot,kind=dp)
   write(*,'(" Thermal balance is ",F5.1,"% converged.")') thermal_percentage
   write(*,'(" [",I6,"/",I6,"]")') i,pdr_ptot
   if (i.eq.pdr_ptot) then
     write(6,*) '#### Converged through thermal balance ####'
     goto 2
   endif
#endif
   RELCH_conv=.true.
   level_converged=.false.
#ifdef OPENMP
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(pp,p,ilevel,RELCH)
#endif 
   do pp=1,pdr_ptot
      p=IDlist_pdr(pp)
      CII_conv(pp)=.true.;CI_conv(pp)=.true.;OI_conv(pp)=.true.;C12O_conv(pp)=.true.
      DO ilevel=1,CII_NLEV
       IF(CII_SOLUTION(pp,ilevel).GE.pdr(p)%abundance(NCx)*1.0D-10) THEN
         IF(CII_SOLUTION(pp,ilevel).EQ.0.0D0 .AND. pdr(p)%CII_pop(ilevel).EQ.0.0D0) THEN
            RELCH=0.0D0
         ELSE
            RELCH=2.0D0*ABS((CII_SOLUTION(pp,ilevel)-pdr(p)%CII_pop(ilevel))&
                &/(CII_SOLUTION(pp,ilevel)+pdr(p)%CII_pop(ilevel)))
         ENDIF
         CII_RELCH(pp,ilevel)=RELCH
         IF(RELCH.GT.1.0D-2) then 
            RELCH_conv=.false.
            CII_conv(pp) = .false.
         ENDIF        
       ENDIF
      ENDDO
      DO ilevel=1,CI_NLEV
       IF(CI_SOLUTION(pp,ilevel).GE.pdr(p)%abundance(NC)*1.0D-10) THEN
         IF(CI_SOLUTION(pp,ilevel).EQ.0.0D0 .AND. pdr(p)%CI_pop(ilevel).EQ.0.0D0) THEN
            RELCH=0.0D0
         ELSE
            RELCH=2.0D0*ABS((CI_SOLUTION(pp,ilevel)-pdr(p)%CI_pop(ilevel))&
                &/(CI_SOLUTION(pp,ilevel)+pdr(p)%CI_pop(ilevel)))
         ENDIF
         CI_RELCH(pp,ilevel)=RELCH
         IF(RELCH.GT.1.0D-2) then
             RELCH_conv=.false.
             CI_conv(pp)=.false.
         endif
       ENDIF
      ENDDO
      DO ilevel=1,OI_NLEV
       IF(OI_SOLUTION(pp,ilevel).GE.pdr(p)%abundance(NO)*1.0D-10) THEN
         IF(OI_SOLUTION(pp,ilevel).EQ.0.0D0 .AND. pdr(p)%OI_pop(ilevel).EQ.0.0D0) THEN
            RELCH=0.0D0
         ELSE
            RELCH=2.0D0*ABS((OI_SOLUTION(pp,ilevel)-pdr(p)%OI_pop(ilevel))&
                &/(OI_SOLUTION(pp,ilevel)+pdr(p)%OI_pop(ilevel)))
         ENDIF
         OI_RELCH(pp,ilevel)=RELCH
         IF(RELCH.GT.1.0D-2) then 
            RELCH_conv=.false.
            OI_conv(pp) = .false.
         endif
       ENDIF
      ENDDO
      DO ilevel=1,C12O_NLEV
       IF(C12O_SOLUTION(pp,ilevel).GE.pdr(p)%abundance(NCO)*1.0D-10) THEN
         IF(C12O_SOLUTION(pp,ilevel).EQ.0.0D0 .AND. pdr(p)%C12O_pop(ilevel).EQ.0.0D0) THEN
            RELCH=0.0D0
         ELSE
            RELCH=2.0D0*ABS((C12O_SOLUTION(pp,ilevel)-pdr(p)%C12O_pop(ilevel))&
                &/(C12O_SOLUTION(pp,ilevel)+pdr(p)%C12O_pop(ilevel)))
         ENDIF
         C12O_RELCH(pp,ilevel)=RELCH
         IF(RELCH.GT.1.0D-2) then 
            RELCH_conv=.false.
            C12O_conv(pp) = .false.
         endif
       ENDIF
      ENDDO

      if (CII_conv(pp).and.CI_conv(pp).and.OI_conv(pp).and.C12O_conv(pp))  then
        level_converged(pp) = .true.
      endif
1001 continue
   enddo
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif
if (.not.relch_conv) then 
   goto 1
else
   write(6,*) '#### Converged through level populations ####'
   levpop_iteration=0

#ifdef THERMALBALANCE
   write(6,*) 'Enabling thermal balance routine in next iteration'
   level_conv=.true.
   goto 1
#else
   goto 2
#endif
endif

1 continue 

   i=0;CII_i=0;CI_i=0;OI_i=0;C12O_i=0
   do p=1,pdr_ptot
      if (level_converged(p)) i=i+1
      if (CII_conv(p)) CII_i=CII_i+1
      if (CI_conv(p)) CI_i=CI_i+1
      if (OI_conv(p)) OI_i=OI_i+1
      if (C12O_conv(p)) C12O_i=C12O_i+1
   enddo
   levpop_percentage  = 100.D0*real(i,kind=dp)/real(pdr_ptot,kind=dp)
   CII_percentage = int(100.D0*real(CII_i,kind=dp)/real(pdr_ptot,kind=dp),kind=i4b)
   write(*,'(" CII is ",F5.1,"% converged.")') CII_percentage
   write(*,'(" [",I6,"/",I6,"]")') CII_i,pdr_ptot
   CI_percentage = int(100.D0*real(CI_i,kind=dp)/real(pdr_ptot,kind=dp),kind=i4b)
   write(*,'(" CI is ",F5.1,"% converged.")') CI_percentage
   write(*,'(" [",I6,"/",I6,"]")') CI_i,pdr_ptot
   OI_percentage = int(100.D0*real(OI_i,kind=dp)/real(pdr_ptot,kind=dp),kind=i4b)
   write(*,'(" OI is ",F5.1,"% converged.")') OI_percentage
   write(*,'(" [",I6,"/",I6,"]")') OI_i,pdr_ptot
   C12O_percentage = int(100.D0*real(C12O_i,kind=dp)/real(pdr_ptot,kind=dp),kind=i4b)
   write(*,'(" CO is ",F5.1,"% converged.")') C12O_percentage
   write(*,'(" [",I6,"/",I6,"]")') C12O_i,pdr_ptot
   write(*,'(" Level populations are ",F5.1,"% converged.")') levpop_percentage
   write(*,'(" [",I6,"/",I6,"]")') i,pdr_ptot

#ifdef THERMALBALANCE
   if (int(levpop_percentage,kind=i4b).ge.100) then
     level_converged = .false.
     write(6,*) 'Resetting [level_converged=.false.] array'
   endif
#endif
   write(6,*) 'Updating population densities...'
#ifdef OPENMP
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(pp,p,ilevel)
#endif
   do pp=1,pdr_ptot
     p=IDlist_pdr(pp)
     DO ilevel=1,CII_NLEV
       pdr(p)%CII_pop(ilevel) = CII_solution(pp,ilevel)
     ENDDO
     DO ilevel=1,CI_NLEV
       pdr(p)%CI_pop(ilevel) = CI_solution(pp,ilevel)
     ENDDO
     DO ilevel=1,OI_NLEV
       pdr(p)%OI_pop(ilevel) = OI_solution(pp,ilevel)
     ENDDO
     DO ilevel=1,C12O_NLEV
       pdr(p)%C12O_pop(ilevel) = C12O_solution(pp,ilevel)
     ENDDO

   enddo
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif
   call cpu_time(t3)
#ifdef OPENMP
   write(6,*) "Number of CPUs: ",CPUs
   write(6,*) 'Iteration time [PARALLEL] = ',(t3-t3b)/real(CPUs),' seconds.'
   write(6,*) 'Total time [PARALLEL] = ',t3/real(CPUs),' seconds.'
#else
   write(6,*) 'Iteration time = ',t3-t3b,' seconds.'
   write(6,*) 'Total time = ',t3,' seconds.'
#endif

END DO !ITERATIONS

2  continue
if (iteration.ge.1) then

   if (iteration.lt.itertot) then
     WRITE(6,*) '3DPDR converged after ',ITERATION-1,' iterations'
   else
     write(6,*) 'Reached maximum number of iterations without convergence.'
     write(6,*) 'To reach convergence, increase the relative number in [params.dat]'
   endif
   write(6,*) 'Writing final outputs'


!-------------------------------------
!OUTPUT FOR CHEMICAL ANALYSIS
!-------------------------------------
 
!   do pp=1,pdr_ptot
!      p=IDlist_pdr(pp)
!      call analyse_chemistry(p, end_time, pdr(p)%rho, previousgastemperature(pp), &
!        &12, pdr(p)%AV(6), nspec, species,pdr(p)%abundance(1:nspec),nreac, reactant, &
!        & product, dummy_rate(:,pp))
!!      write(6,*)
!   enddo

!-------------------------------------
!OUTPUT FOR ABUNDANCES AND TEMPERATURE
!-------------------------------------
   out_file = trim(adjustl(output))//".pdr.fin"
   out_file2 = trim(adjustl(out_file))//"]"
   write(6,'(" Writing file [",A)') out_file2
   open(unit=21,file=out_file,status='replace')

  do pp=1,pdr_ptot-2
     p=IDlist_pdr(pp)
#ifdef PSEUDO_1D
     write(21,'(I7,4ES11.3,I5,300ES11.3)') p,pdr(p)%x, pdr(p)%AV(6), previousgastemperature(pp),dusttemperature(pp),pdr(p)%etype,&
     &pdr(p)%rho,pdr(p)%UVfield,pdr(p)%abundance
#else
     write(21,'(I7,5ES11.3,I5,300ES11.3)') p,pdr(p)%x, pdr(p)%y, pdr(p)%z, previousgastemperature(pp),dusttemperature(pp),&
     &pdr(p)%etype,pdr(p)%rho,pdr(p)%UVfield,pdr(p)%abundance,pdr(p)%AV
#endif
  enddo
  pp=pdr_ptot
  p=IDlist_pdr(pp)
#ifdef PSEUDO_1D
  write(21,'(I7,4ES11.3,I5,300ES11.3)') p,pdr(p)%x, pdr(p)%AV(6), previousgastemperature(pp),dusttemperature(pp),pdr(p)%etype,&
  &pdr(p)%rho,pdr(p)%UVfield,pdr(p)%abundance
#else
  write(21,'(I7,5ES11.3,I5,300ES11.3)') p,pdr(p)%x, pdr(p)%y, pdr(p)%z, previousgastemperature(pp),dusttemperature(pp),&
  &pdr(p)%etype,pdr(p)%rho,pdr(p)%UVfield,pdr(p)%abundance,pdr(p)%AV
#endif
  pp=pdr_ptot-1
  p=IDlist_pdr(pp)
#ifdef PSEUDO_1D
  write(21,'(I7,4ES11.3,I5,300ES11.3)') p,pdr(p)%x, pdr(p)%AV(6), previousgastemperature(pp),dusttemperature(pp),pdr(p)%etype,&
  &pdr(p)%rho,pdr(p)%UVfield,pdr(p)%abundance
#else
  write(21,'(I7,5ES11.3,I5,300ES11.3)') p,pdr(p)%x, pdr(p)%y, pdr(p)%z, previousgastemperature(pp),dusttemperature(pp),&
  &pdr(p)%etype,pdr(p)%rho,pdr(p)%UVfield,pdr(p)%abundance,pdr(p)%AV
#endif


if (ion_ptot.gt.0) then

   out_file = trim(adjustl(output))//".ion.fin"
   out_file2 = trim(adjustl(out_file))//"]"
   write(6,'(" Writing file [",A)') out_file2
   close(21);open(unit=21,file=out_file,status='replace')

  do pp=1,ion_ptot
     p=IDlist_ion(pp)
#ifdef PSEUDO_1D
     write(21,'(I7,3ES11.3,I5,300ES11.3)') p,pdr(p)%x, pdr(p)%AV(6), previousgastemperature(pp),pdr(p)%etype,&
     &pdr(p)%rho,pdr(p)%UVfield,pdr(p)%abundance
#else
     write(21,'(I7,4ES11.3,I5,400ES11.3)') p,pdr(p)%x, pdr(p)%y, pdr(p)%z, previousgastemperature(pp),pdr(p)%etype,&
     &pdr(p)%rho,pdr(p)%UVfield,pdr(p)%abundance,pdr(p)%AV(:)
#endif
  enddo
endif

if (dark_ptot.gt.0) then

   out_file = trim(adjustl(output))//".mol.fin"
   out_file2 = trim(adjustl(out_file))//"]"
   write(6,'(" Writing file [",A)') out_file2
   close(21);open(unit=21,file=out_file,status='replace')

  do pp=1,dark_ptot
     p=IDlist_dark(pp)
#ifdef PSEUDO_1D
     write(21,'(I7,3ES11.3,I5,300ES11.3)') p,pdr(p)%x, pdr(p)%AV(6), previousgastemperature(pp),pdr(p)%etype,&
     &pdr(p)%rho,pdr(p)%UVfield,pdr(p)%abundance
#else
     write(21,'(I7,4ES11.3,I5,400ES11.3)') p,pdr(p)%x, pdr(p)%y, pdr(p)%z, previousgastemperature(pp),pdr(p)%etype,&
     &pdr(p)%rho,pdr(p)%UVfield,pdr(p)%abundance,pdr(p)%AV(:)
#endif
  enddo
endif

close(21)
!-----------------------------------------
!END OUTPUT FOR ABUNDANCES AND TEMPERATURE
!-----------------------------------------


!---------------------------
!OUTPUT FOR COOLING FUNCTION
!---------------------------
   out_file = trim(adjustl(output))//trim(adjustl(".cool"))//".fin"
   out_file2 = trim(adjustl(out_file))//"]"
   write(6,'(" Writing file [",A)') out_file2
   open(unit=13,file=out_file,status='replace')

   do pp=1,pdr_ptot
      p=IDlist_pdr(pp)
#ifdef PSEUDO_1D
      write(13,'(I7,200ES11.3)') p, pdr(p)%x, pdr(p)%AV(6), CII_cool(pp),CI_cool(pp), & 
                               &OI_cool(pp),C12O_cool(pp), total_cooling_rate(pp)
#else
      write(13,'(I7,200ES11.3)') p, pdr(p)%x, pdr(p)%y, pdr(p)%z, CII_cool(pp), CI_cool(pp),&
                     &OI_cool(pp), C12O_cool(pp), total_cooling_rate(pp), pdr(p)%AV(:)
#endif
   enddo
!-------------------------------
!END OUTPUT FOR COOLING FUNCTION
!-------------------------------


!---------------------------
!OUTPUT FOR HEATING FUNCTION
!---------------------------
   out_file = trim(adjustl(output))//trim(adjustl(".heat"))//".fin"
   out_file2 = trim(adjustl(out_file))//"]"
   write(6,'(" Writing file [",A)') out_file2
   open(unit=14,file=out_file,status='replace')


   do pp=1,pdr_ptot
      p=IDlist_pdr(pp)
#ifdef PSEUDO_1D
      write(14,'(I7,200ES11.3)') p, pdr(p)%x, pdr(p)%AV(6), all_heating(pp,:)
#else
      write(14,'(I7,200ES11.3)') p, pdr(p)%x, pdr(p)%y, pdr(p)%z, all_heating(pp,:), pdr(p)%AV(:)
#endif
   enddo
   close(14)
!-------------------------------
!END OUTPUT FOR HEATING FUNCTION
!-------------------------------

!---------------------------
!OUTPUT FOR TRANSITION LINES
!---------------------------
   out_file = trim(adjustl(output))//trim(adjustl(".line"))//".fin"
   out_file2 = trim(adjustl(out_file))//"]"
   write(6,'(" Writing file [",A)') out_file2
   open(unit=16,file=out_file,status='replace')

   do pp=1,pdr_ptot-2
      p=IDlist_pdr(pp)
#ifdef PSEUDO_1D
      write(16,'(I9,200ES11.3)') pp, pdr(p)%x, pdr(p)%AV(6), &    !ID,x,AV(6)
      &pdr(p)%CII_line(2,1),&                                          !CII line
      &pdr(p)%CI_line(2,1), pdr(p)%CI_line(3,1), pdr(p)%CI_line(3,2),& !CI line
      &pdr(p)%OI_line(2,1), pdr(p)%OI_line(3,1), pdr(p)%OI_line(3,2),& !OI line
      &pdr(p)%C12O_line(2,1), pdr(p)%C12O_line(3,2), pdr(p)%C12O_line(4,3), pdr(p)%C12O_line(5,4), &
      &pdr(p)%C12O_line(6,5), pdr(p)%C12O_line(7,6), pdr(p)%C12O_line(8,7), pdr(p)%C12O_line(9,8), &
      &pdr(p)%C12O_line(10,9), pdr(p)%C12O_line(11,10)     !CO line
#else
      write(16,'(I9,200ES11.3)') pp, pdr(p)%x, pdr(p)%y, pdr(p)%z, &    !ID,x,y,z
      &pdr(p)%CII_line(2,1),&                                          !CII line
      &pdr(p)%CI_line(2,1), pdr(p)%CI_line(3,1), pdr(p)%CI_line(3,2),& !CI line
      &pdr(p)%OI_line(2,1), pdr(p)%OI_line(3,1), pdr(p)%OI_line(3,2),& !OI line
      &pdr(p)%C12O_line(2,1), pdr(p)%C12O_line(3,2), pdr(p)%C12O_line(4,3), pdr(p)%C12O_line(5,4), pdr(p)%C12O_line(6,5), &
      &pdr(p)%C12O_line(7,6), pdr(p)%C12O_line(8,7), pdr(p)%C12O_line(9,8), pdr(p)%C12O_line(10,9), pdr(p)%C12O_line(11,10),&     !CO line
      &pdr(p)%AV(:)
#endif
   enddo
      pp=pdr_ptot
      p=IDlist_pdr(pp)
#ifdef PSEUDO_1D
      write(16,'(I9,200ES11.3)') pp, pdr(p)%x, pdr(p)%AV(6), &    !ID,x,AV(6)
      &pdr(p)%CII_line(2,1),&                                          !CII line
      &pdr(p)%CI_line(2,1), pdr(p)%CI_line(3,1), pdr(p)%CI_line(3,2),& !CI line
      &pdr(p)%OI_line(2,1), pdr(p)%OI_line(3,1), pdr(p)%OI_line(3,2),& !OI line
      &pdr(p)%C12O_line(2,1), pdr(p)%C12O_line(3,2), pdr(p)%C12O_line(4,3), pdr(p)%C12O_line(5,4), &
      &pdr(p)%C12O_line(6,5), pdr(p)%C12O_line(7,6), pdr(p)%C12O_line(8,7), pdr(p)%C12O_line(9,8), &
      &pdr(p)%C12O_line(10,9), pdr(p)%C12O_line(11,10)     !CO line
#else
      write(16,'(I9,200ES11.3)') pp, pdr(p)%x, pdr(p)%y, pdr(p)%z, &    !ID,x,y,z
      &pdr(p)%CII_line(2,1),&                                          !CII line
      &pdr(p)%CI_line(2,1), pdr(p)%CI_line(3,1), pdr(p)%CI_line(3,2),& !CI line
      &pdr(p)%OI_line(2,1), pdr(p)%OI_line(3,1), pdr(p)%OI_line(3,2),& !OI line
      &pdr(p)%C12O_line(2,1), pdr(p)%C12O_line(3,2), pdr(p)%C12O_line(4,3), pdr(p)%C12O_line(5,4), pdr(p)%C12O_line(6,5), &
      &pdr(p)%C12O_line(7,6), pdr(p)%C12O_line(8,7), pdr(p)%C12O_line(9,8), pdr(p)%C12O_line(10,9), pdr(p)%C12O_line(11,10),&     !CO line
      &pdr(p)%AV(:)
#endif
      pp=pdr_ptot-1
      p=IDlist_pdr(pp)
#ifdef PSEUDO_1D
      write(16,'(I9,200ES11.3)') pp, pdr(p)%x, pdr(p)%AV(6), &    !ID,x,AV(6)
      &pdr(p)%CII_line(2,1),&                                          !CII line
      &pdr(p)%CI_line(2,1), pdr(p)%CI_line(3,1), pdr(p)%CI_line(3,2),& !CI line
      &pdr(p)%OI_line(2,1), pdr(p)%OI_line(3,1), pdr(p)%OI_line(3,2),& !OI line
      &pdr(p)%C12O_line(2,1), pdr(p)%C12O_line(3,2), pdr(p)%C12O_line(4,3), pdr(p)%C12O_line(5,4), &
      &pdr(p)%C12O_line(6,5), pdr(p)%C12O_line(7,6), pdr(p)%C12O_line(8,7), pdr(p)%C12O_line(9,8), &
      &pdr(p)%C12O_line(10,9), pdr(p)%C12O_line(11,10)     !CO line
#else
      write(16,'(I9,200ES11.3)') pp, pdr(p)%x, pdr(p)%y, pdr(p)%z, &    !ID,x,y,z
      &pdr(p)%CII_line(2,1),&                                          !CII line
      &pdr(p)%CI_line(2,1), pdr(p)%CI_line(3,1), pdr(p)%CI_line(3,2),& !CI line
      &pdr(p)%OI_line(2,1), pdr(p)%OI_line(3,1), pdr(p)%OI_line(3,2),& !OI line
      &pdr(p)%C12O_line(2,1), pdr(p)%C12O_line(3,2), pdr(p)%C12O_line(4,3), pdr(p)%C12O_line(5,4), pdr(p)%C12O_line(6,5), &
      &pdr(p)%C12O_line(7,6), pdr(p)%C12O_line(8,7), pdr(p)%C12O_line(9,8), pdr(p)%C12O_line(10,9), pdr(p)%C12O_line(11,10),&     !CO line
      &pdr(p)%AV(:)
#endif
   close(16)
!-------------------------------
!END OUTPUT FOR TRANSITION LINES
!-------------------------------

!---------------------------
!OUTPUT FOR OPTICAL DEPTHS
!---------------------------
   out_file = trim(adjustl(output))//trim(adjustl(".spop"))//".fin"
   out_file2 = trim(adjustl(out_file))//"]"
   write(6,'(" Writing file [",A)') out_file2
   open(unit=16,file=out_file,status='replace')

   do pp=1,pdr_ptot
      p=IDlist_pdr(pp)
      write(16,'(I9,200ES11.3)') pp, pdr(p)%Av(6), pdr(p)%CII_pop, pdr(p)%CI_pop, pdr(p)%OI_pop, pdr(p)%C12O_pop
   enddo
   close(16)
!-------------------------------
!END OUTPUT FOR OPTICAL DEPTHS
!-------------------------------



#ifdef PSEUDO_1D
!---------------------------
!OUTPUT FOR OPTICAL DEPTHS
!---------------------------
   out_file = trim(adjustl(output))//trim(adjustl(".opdp"))//".fin"
   out_file2 = trim(adjustl(out_file))//"]"
   write(6,'(" Writing file [",A)') out_file2
   open(unit=16,file=out_file,status='replace')

   do pp=1,pdr_ptot
      p=IDlist_pdr(pp)
      write(16,'(I9,200ES11.3)') pp, pdr(p)%x,&     !ID,x,AV(6)
      &pdr(p)%CII_optdepth(2,1,6),&                                          !CII line
      &pdr(p)%CI_optdepth(2,1,6), pdr(p)%CI_optdepth(3,1,6), pdr(p)%CI_optdepth(3,2,6),& !CI line
      &pdr(p)%OI_optdepth(2,1,6), pdr(p)%OI_optdepth(3,1,6), pdr(p)%OI_optdepth(3,2,6),& !OI line
      &pdr(p)%C12O_optdepth(2,1,6), pdr(p)%C12O_optdepth(3,2,6), pdr(p)%C12O_optdepth(4,3,6), pdr(p)%C12O_optdepth(5,4,6), &
      &pdr(p)%C12O_optdepth(6,5,6), pdr(p)%C12O_optdepth(7,6,6), pdr(p)%C12O_optdepth(8,7,6), pdr(p)%C12O_optdepth(9,8,6), &
      &pdr(p)%C12O_optdepth(10,9,6), pdr(p)%C12O_optdepth(11,10,6)     !CO line
   enddo
   close(16)
!-------------------------------
!END OUTPUT FOR OPTICAL DEPTHS
!-------------------------------
#endif


!==============
endif

write(6,*) ''
call cpu_time(t4)
#ifdef OPENMP
write(6,*) "Number of CPUs: ",CPUs
write(6,*) 'Simulation time [PARALLEL] = ',(t4-t2)/real(CPUs),' seconds.'
#else
write(6,*) 'Simulation time = ',t4-t2,' seconds.'
#endif
write(6,*) 'Finished !'

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~END OF MAIN PROGRAM~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
contains
 

subroutine allocations

!load SCO_GRID data [UCL_PDR]
SCO_GRID(1:8,1) = (/0.000D+00,-1.408D-02,-1.099D-01,-4.400D-01,-1.154D+00,-1.888D+00,-2.760D+00,-4.001D+00/)
SCO_GRID(1:8,2) = (/-8.539D-02,-1.015D-01,-2.104D-01,-5.608D-01,-1.272D+00,-1.973D+00,-2.818D+00,-4.055D+00/)
SCO_GRID(1:8,3) = (/-1.451D-01,-1.612D-01,-2.708D-01,-6.273D-01,-1.355D+00,-2.057D+00,-2.902D+00,-4.122D+00/)
SCO_GRID(1:8,4) = (/-4.559D-01,-4.666D-01,-5.432D-01,-8.665D-01,-1.602D+00,-2.303D+00,-3.146D+00,-4.421D+00/)
SCO_GRID(1:8,5) = (/-1.303D+00,-1.312D+00,-1.367D+00,-1.676D+00,-2.305D+00,-3.034D+00,-3.758D+00,-5.077D+00/)
SCO_GRID(1:8,6) = (/-3.883D+00,-3.888D+00,-3.936D+00,-4.197D+00,-4.739D+00,-5.165D+00,-5.441D+00,-6.446D+00/)

!CII
ALLOCATE(CII_ENERGIES(1:CII_NLEV))
ALLOCATE(CII_WEIGHTS(1:CII_NLEV))
ALLOCATE(CII_A_COEFFS(1:CII_NLEV,1:CII_NLEV))
ALLOCATE(CII_B_COEFFS(1:CII_NLEV,1:CII_NLEV))
ALLOCATE(CII_FREQUENCIES(1:CII_NLEV,1:CII_NLEV))
ALLOCATE(CII_TEMPERATURES(1:7,1:CII_NTEMP))
ALLOCATE(CII_HP(1:CII_NLEV,1:CII_NLEV,1:CII_NTEMP))
ALLOCATE(CII_H(1:CII_NLEV,1:CII_NLEV,1:CII_NTEMP))
ALLOCATE(CII_EL(1:CII_NLEV,1:CII_NLEV,1:CII_NTEMP))
ALLOCATE(CII_HE(1:CII_NLEV,1:CII_NLEV,1:CII_NTEMP))
ALLOCATE(CII_H2(1:CII_NLEV,1:CII_NLEV,1:CII_NTEMP))
ALLOCATE(CII_PH2(1:CII_NLEV,1:CII_NLEV,1:CII_NTEMP))
ALLOCATE(CII_OH2(1:CII_NLEV,1:CII_NLEV,1:CII_NTEMP))
!CI
ALLOCATE(CI_ENERGIES(1:CI_NLEV))
ALLOCATE(CI_WEIGHTS(1:CI_NLEV))
ALLOCATE(CI_A_COEFFS(1:CI_NLEV,1:CI_NLEV))
ALLOCATE(CI_B_COEFFS(1:CI_NLEV,1:CI_NLEV))
ALLOCATE(CI_FREQUENCIES(1:CI_NLEV,1:CI_NLEV))
ALLOCATE(CI_TEMPERATURES(1:7,1:CI_NTEMP))
ALLOCATE(CI_HP(1:CI_NLEV,1:CI_NLEV,1:CI_NTEMP))
ALLOCATE(CI_H(1:CI_NLEV,1:CI_NLEV,1:CI_NTEMP))
ALLOCATE(CI_EL(1:CI_NLEV,1:CI_NLEV,1:CI_NTEMP))
ALLOCATE(CI_HE(1:CI_NLEV,1:CI_NLEV,1:CI_NTEMP))
ALLOCATE(CI_H2(1:CI_NLEV,1:CI_NLEV,1:CI_NTEMP))
ALLOCATE(CI_PH2(1:CI_NLEV,1:CI_NLEV,1:CI_NTEMP))
ALLOCATE(CI_OH2(1:CI_NLEV,1:CI_NLEV,1:CI_NTEMP))
!OI
ALLOCATE(OI_ENERGIES(1:OI_NLEV))
ALLOCATE(OI_WEIGHTS(1:OI_NLEV))
ALLOCATE(OI_A_COEFFS(1:OI_NLEV,1:OI_NLEV))
ALLOCATE(OI_B_COEFFS(1:OI_NLEV,1:OI_NLEV))
ALLOCATE(OI_FREQUENCIES(1:OI_NLEV,1:OI_NLEV))
ALLOCATE(OI_TEMPERATURES(1:7,1:OI_NTEMP))
ALLOCATE(OI_HP(1:OI_NLEV,1:OI_NLEV,1:OI_NTEMP))
ALLOCATE(OI_H(1:OI_NLEV,1:OI_NLEV,1:OI_NTEMP))
ALLOCATE(OI_EL(1:OI_NLEV,1:OI_NLEV,1:OI_NTEMP))
ALLOCATE(OI_HE(1:OI_NLEV,1:OI_NLEV,1:OI_NTEMP))
ALLOCATE(OI_H2(1:OI_NLEV,1:OI_NLEV,1:OI_NTEMP))
ALLOCATE(OI_PH2(1:OI_NLEV,1:OI_NLEV,1:OI_NTEMP))
ALLOCATE(OI_OH2(1:OI_NLEV,1:OI_NLEV,1:OI_NTEMP))
!12CO
ALLOCATE(C12O_ENERGIES(1:C12O_NLEV))
ALLOCATE(C12O_WEIGHTS(1:C12O_NLEV))
ALLOCATE(C12O_A_COEFFS(1:C12O_NLEV,1:C12O_NLEV))
ALLOCATE(C12O_B_COEFFS(1:C12O_NLEV,1:C12O_NLEV))
ALLOCATE(C12O_FREQUENCIES(1:C12O_NLEV,1:C12O_NLEV))
ALLOCATE(C12O_TEMPERATURES(1:7,1:C12O_NTEMP))
ALLOCATE(C12O_HP(1:C12O_NLEV,1:C12O_NLEV,1:C12O_NTEMP))
ALLOCATE(C12O_H(1:C12O_NLEV,1:C12O_NLEV,1:C12O_NTEMP))
ALLOCATE(C12O_EL(1:C12O_NLEV,1:C12O_NLEV,1:C12O_NTEMP))
ALLOCATE(C12O_HE(1:C12O_NLEV,1:C12O_NLEV,1:C12O_NTEMP))
ALLOCATE(C12O_H2(1:C12O_NLEV,1:C12O_NLEV,1:C12O_NTEMP))
ALLOCATE(C12O_PH2(1:C12O_NLEV,1:C12O_NLEV,1:C12O_NTEMP))
ALLOCATE(C12O_OH2(1:C12O_NLEV,1:C12O_NLEV,1:C12O_NTEMP))

allocate(C_COEFFS(1:NLEV,1:NLEV))

allocate(species(1:nspec))
allocate(dummyabundance(1:nspec))
allocate(mass(1:nspec))
allocate(reactant(1:nreac,1:3))
allocate(product(1:nreac,1:4))
allocate(rate(1:nreac))
allocate(alpha(1:nreac))
allocate(beta(1:nreac))
allocate(gamma(1:nreac))
allocate(rtmin(1:nreac))
allocate(rtmax(1:nreac))
allocate(duplicate(1:nreac))

!allocations start from 0 to cope with the ONE dark molecular element
allocate(total_cooling_rate(0:pdr_ptot))
allocate(CII_cool(0:pdr_ptot))
allocate(CI_cool(0:pdr_ptot))
allocate(OI_cool(0:pdr_ptot))
allocate(C12O_cool(0:pdr_ptot))

allocate(all_heating(0:pdr_ptot,1:12))

allocate(dusttemperature(0:pdr_ptot))
allocate(gastemperature(0:pdr_ptot))
allocate(previousgastemperature(0:pdr_ptot))

#ifdef THERMALBALANCE
allocate(Fratio(0:pdr_ptot))
allocate(Fmean(0:pdr_ptot));Fmean=0
allocate(Tlow(0:pdr_ptot))
allocate(Thigh(0:pdr_ptot))
#endif

return
end subroutine allocations



SUBROUTINE CALC_UVFIELD
! Calculation of the UVfield
! Updated to allow multiple sources of UV radiation
! Added option for point source - assumes source is at coords (0,0,0)
! T.Bisbas, T.Bell

      REAL(KIND=DP) :: Gpoint(1:3)
      real(kind=dp)::maxdotproduct, auxmaxdotproduct
      integer(kind=i4b) :: jj_min,jj
      integer::epray_new


do p=1,grand_ptot
      pdr(p)%UVfield = 0.0D0
enddo

DO pp=1,pdr_ptot
  p=IDlist_pdr(pp)
  pdr(p)%rad_surface = 0.0D0
  pdr(p)%columndensity = 0.0D0
  pdr(p)%AV = 0.0D0
  pdr(p)%projected(:,0)=p        
  DO j=0,nrays-1
     IF (pdr(p)%epray(j).GT.0) THEN
        DO i=1,pdr(p)%epray(j)
           adaptive_step = SQRT((pdr(p)%epoint(1,j,i-1)-pdr(p)%epoint(1,j,i))**2 + &
                              & (pdr(p)%epoint(2,j,i-1)-pdr(p)%epoint(2,j,i))**2 + &
                              & (pdr(p)%epoint(3,j,i-1)-pdr(p)%epoint(3,j,i))**2)
           pdr(p)%columndensity(j) = pdr(p)%columndensity(j) + &
              &((pdr(INT(pdr(p)%projected(j,i-1)))%rho + &
              &pdr(INT(pdr(p)%projected(j,i)))%rho)/2.)*adaptive_step*pc
        ENDDO
     ENDIF
     pdr(p)%AV(j) = pdr(p)%columndensity(j)*AV_fac
  ENDDO

  DO j=0,nrays-1
    ! Isotropic illumination, with strength given by Gx
    IF (fieldchoice.EQ."ISO") pdr(p)%rad_surface(j) = Gext(1)/real(nrays,kind=dp)

    ! Uni-directional illumination, with vector defined by (Gx,Gy,Gz)
    IF (fieldchoice.EQ."UNI") THEN
           if (j==0) then 
           !calculate which is the minimum negative dot_product of Gext and
           !vectors. The minimum dot_product corresponds to the direction which
           !the plane parallel radiation is impinging. Store the ID of the ray
           !which has this smallest dot_product (jj_min).
           maxdotproduct=0
            do jj=0,nrays-1
             auxmaxdotproduct=-dot_product(Gext(:),vectors(:,jj))
             if (auxmaxdotproduct.ge.maxdotproduct) then 
                maxdotproduct=auxmaxdotproduct
                jj_min=jj
             endif
            enddo
          endif 
        if (j.eq.jj_min)  pdr(p)%rad_surface(j) = -DOT_PRODUCT(Gext(:),vectors(:,j))
    ENDIF

   ENDDO ! End of loop over rays
   DO j=0,nrays-1
            ! Impose a lower limit on the value of UVfield to prevent numerical issues
            pdr(p)%UVfield = pdr(p)%UVfield + pdr(p)%rad_surface(j)*EXP(-pdr(p)%AV(j)*UV_fac)
            IF (pdr(p)%UVfield.LT.1.0D-50) pdr(p)%UVfield = 0.0D0
   ENDDO ! End of loop over rays
ENDDO ! End of loop over pdrpoints

#ifdef THERMALBALANCE
#ifdef GUESS_TEMP
ALLOCATE(Tmin_array(0:pdr_ptot))
ALLOCATE(Tmax_array(0:pdr_ptot))
DO p=1,pdr_ptot
   Tguess = 10.0D0*(1.0D0+(2.*pdr(IDlist_pdr(p))%UVfield)**(1.0D0/3.0D0))
   gastemperature(p) = Tguess
   previousgastemperature(p) = Tguess

   Tlow(p) = Tguess/2.0D0
   Thigh(p) = Tguess*1.5D0

   IF (Tlow(p).LT.Tmin)  Tlow(p)  = Tmin
   IF (Thigh(p).GT.Tmax) Thigh(p) = Tmax

   Tmin_array(p) =  Tlow(p)/3.0D0 ! Bound minimum
   Tmax_array(p) = Thigh(p)*2.0D0 ! Bound maximum
   IF (Tmin_array(p).LT.Tmin) Tmin_array(p) = Tmin
   IF (Tmax_array(p).GT.Tmax) Tmax_array(p) = Tmax
ENDDO
#endif
#endif
if (dark_ptot.gt.0) then
!Dark Molecular element-------------------------------------------------------------------------------------------------
p=IDlist_dark(1)
pdr(p)%rad_surface = 0.0D0
pdr(p)%columndensity = 0.0D0
pdr(p)%AV = 0.0D0
pdr(p)%projected(:,0)=p        
DO j=0,nrays-1
  IF (pdr(p)%epray(j).GT.0) THEN
     DO i=1,pdr(p)%epray(j)
        adaptive_step = SQRT((pdr(p)%epoint(1,j,i-1)-pdr(p)%epoint(1,j,i))**2 + &
                           & (pdr(p)%epoint(2,j,i-1)-pdr(p)%epoint(2,j,i))**2 + &
                           & (pdr(p)%epoint(3,j,i-1)-pdr(p)%epoint(3,j,i))**2)
        pdr(p)%columndensity(j) = pdr(p)%columndensity(j) + &
                 &((pdr(INT(pdr(p)%projected(j,i-1)))%rho + &
                 &pdr(INT(pdr(p)%projected(j,i)))%rho)/2.)*adaptive_step*pc
     ENDDO
  ENDIF
  pdr(p)%AV(j) = pdr(p)%columndensity(j)*AV_fac
ENDDO

DO j=0,nrays-1
   ! Isotropic illumination, with strength given by Gx
    IF (fieldchoice.EQ."ISO") pdr(p)%rad_surface(j) = Gext(1)/real(nrays,kind=dp)
    ! Uni-directional illumination, with vector defined by (Gx,Gy,Gz)
    IF (fieldchoice.EQ."UNI") THEN
        if (j==0) then 
        !calculate which is the minimum negative dot_product of Gext and
        !vectors. The minimum dot_product corresponds to the direction which
        !the plane parallel radiation is impinging. Store the ID of the ray
        !which has this smallest dot_product (jj_min).
        maxdotproduct=0
        do jj=0,nrays-1
          auxmaxdotproduct=-dot_product(Gext(:),vectors(:,jj))
          if (auxmaxdotproduct.ge.maxdotproduct) then 
             maxdotproduct=auxmaxdotproduct
             jj_min=jj
          endif
        enddo
     endif
     if (j.eq.jj_min)  pdr(p)%rad_surface(j) = pdr(p)%rad_surface(j)-DOT_PRODUCT(Gext(:),vectors(:,j))
  ENDIF !over fieldchoices
ENDDO ! End of loop over rays

DO j=0,nrays-1
    ! Impose a lower limit on the value of UVfield to prevent numerical issues
    pdr(p)%UVfield = pdr(p)%UVfield + pdr(p)%rad_surface(j)*EXP(-pdr(p)%AV(j)*UV_fac) 
   IF (pdr(p)%UVfield.LT.1.0D-50) pdr(p)%UVfield = 0.0D0
ENDDO

gastemperature(0)=Tmin
previousgastemperature(0)=Tmin
#ifdef THERMALBALANCE
Tlow(0)=Tmin
Thigh(0)=Tmin
#endif
endif
!-------------------------------------
RETURN

END SUBROUTINE CALC_UVFIELD


subroutine calc_columndens
!calculation of column density
!T.Bisbas


!DARK MOLECULAR ELEMENT -------------------------------------------
if (dark_ptot.gt.0) then
  p=IDlist_dark(1)
  if (referee.eq.0) allocate(column(0)%columndens_point(0:nrays-1,1:nspec))
  column(0)%columndens_point = 0.0D0
  do j=0,nrays-1
     if (pdr(p)%epray(j).eq.0) cycle
     do i=1,pdr(p)%epray(j)!nb_projected_points
       adaptive_step = sqrt((pdr(p)%epoint(1,j,i-1)-pdr(p)%epoint(1,j,i))**2+&
                  &(pdr(p)%epoint(2,j,i-1)-pdr(p)%epoint(2,j,i))**2 + &
                  &(pdr(p)%epoint(3,j,i-1)-pdr(p)%epoint(3,j,i))**2)
       do k=1,nspec
          column(0)%columndens_point(j,k) = column(0)%columndens_point(j,k) + adaptive_step*PC*&
             & (pdr(int(pdr(p)%projected(j,i-1)))%rho*pdr(int(pdr(p)%projected(j,i-1)))%abundance(k) +&
             & pdr(int(pdr(p)%projected(j,i)))%rho*pdr(int(pdr(p)%projected(j,i)))%abundance(k))/2.
       enddo ! End of loop over species (n)
     enddo ! End of loop over evaluation points along ray (i)
  enddo ! End of j loop over rays (j)

endif 
!-------------------------------------------------------------------------

#ifdef OPENMP
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(pp,p,adaptive_step,j,i,k)
#endif
do pp=1,pdr_ptot
   p=IDlist_pdr(pp)
   pdr(p)%projected(:,0)=p
#ifdef THERMALBALANCE
  if (converged(pp)) cycle
#endif
  if (referee.eq.0) allocate(column(pp)%columndens_point(0:nrays-1,1:nspec))
  column(pp)%columndens_point = 0.0D0
  do j=0,nrays-1
     if (pdr(p)%epray(j).eq.0) cycle
     do i=1,pdr(p)%epray(j)!nb_projected_points
       adaptive_step = sqrt((pdr(p)%epoint(1,j,i-1)-pdr(p)%epoint(1,j,i))**2+&
                  &(pdr(p)%epoint(2,j,i-1)-pdr(p)%epoint(2,j,i))**2 + &
                  &(pdr(p)%epoint(3,j,i-1)-pdr(p)%epoint(3,j,i))**2)
       do k=1,nspec
          column(pp)%columndens_point(j,k) = column(pp)%columndens_point(j,k) + adaptive_step*PC*&
             & (pdr(int(pdr(p)%projected(j,i-1)))%rho*pdr(int(pdr(p)%projected(j,i-1)))%abundance(k) +&
             & pdr(int(pdr(p)%projected(j,i)))%rho*pdr(int(pdr(p)%projected(j,i)))%abundance(k))/2.
       enddo ! End of loop over species (n)
     enddo ! End of loop over evaluation points along ray (i)
  enddo ! End of j loop over rays (j)
enddo ! End of ii loop over pdrpoints (ii)
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif
return
end subroutine calc_columndens


! Calculate the partition function for the given species
SUBROUTINE CALCULATE_PARTITION_FUNCTION(PARTITION_FUNCTION,NLEV,ENERGIES,WEIGHTS,TEMPERATURE)

      USE HEALPIX_TYPES
      IMPLICIT NONE

      INTEGER(KIND=I4B), INTENT(IN) :: NLEV
      REAL(KIND=DP), INTENT(IN)     :: ENERGIES(1:NLEV),WEIGHTS(1:NLEV)
      REAL(KIND=DP), INTENT(IN)     :: TEMPERATURE
      REAL(KIND=DP), INTENT(OUT)    :: PARTITION_FUNCTION

      INTEGER(KIND=I4B) :: ILEVEL

      PARTITION_FUNCTION=0.0D0
      DO ILEVEL=1,NLEV
         PARTITION_FUNCTION=PARTITION_FUNCTION + WEIGHTS(ILEVEL)*EXP(-ENERGIES(ILEVEL)/KB/TEMPERATURE)
      ENDDO

RETURN
END SUBROUTINE CALCULATE_PARTITION_FUNCTION

! Calculate the LTE level populations for the given species
SUBROUTINE CALCULATE_LTE_POPULATIONS(NLEV,LEVEL_POP,ENERGIES,WEIGHTS,PARTITION_FUNCTION,DENSITY,TEMPERATURE)

      USE HEALPIX_TYPES
      IMPLICIT NONE

      INTEGER(KIND=I4B), INTENT(IN) :: NLEV
      REAL(KIND=DP), INTENT(IN)     :: ENERGIES(1:NLEV),WEIGHTS(1:NLEV)
      REAL(KIND=DP), INTENT(IN)     :: PARTITION_FUNCTION
      REAL(KIND=DP), INTENT(IN)     :: DENSITY,TEMPERATURE
      REAL(KIND=DP), INTENT(OUT)    :: LEVEL_POP(1:NLEV)

      INTEGER(KIND=I4B) :: ILEVEL
      REAL(KIND=DP) :: TOTAL_POP

  
      TOTAL_POP=0.0D0
      DO ILEVEL=1,NLEV
         LEVEL_POP(ILEVEL)=DENSITY*WEIGHTS(ILEVEL)*EXP(-ENERGIES(ILEVEL)/KB/TEMPERATURE)/PARTITION_FUNCTION
         TOTAL_POP=TOTAL_POP + LEVEL_POP(ILEVEL)
      ENDDO

      ! Check that the sum of the level populations adds up to the total density
      IF(ABS(TOTAL_POP-DENSITY)/DENSITY.GT.1.0D-3) THEN
         WRITE(6,*)'ERROR! Sum of LTE level populations differs from the total density by ', &
		 & INT(1.0D2*ABS(TOTAL_POP-DENSITY)/DENSITY),'%'
         STOP
      ENDIF

RETURN
END SUBROUTINE CALCULATE_LTE_POPULATIONS

subroutine dark_molecular_region

write(6,*) '' 
write(6,*) '*** Dark Molecular Region ***'
write(6,*) 'Calculating LTE level populations...'
n_H = 2.0D0*rho_max!maxval(pdr(:)%rho) !needed for ODEs
write(6,*) 'Density = ',n_H
DO II=1,CHEMITERATIONS
   write(6,*) 'Chemical iteration ',ii

 allocate(dummy_rate(1:nreac,1))

   p=IDlist_dark(1)
    CALL CALCULATE_REACTION_RATES(gastemperature(0),dusttemperature(0),nrays,pdr(p)%rad_surface(0:nrays-1),&
          &pdr(p)%AV(0:nrays-1),column(0)%columndens_point(0:nrays-1,1:nspec),&
          &nreac, reactant, product, alpha, beta, gamma, rate, rtmin, rtmax, duplicate, nspec,&
          &NRGR,NRH2,NRHD,NRCO,NRCI,NRSI)
      DUMMY_RATE(:,1) = rate

   allocate(dummy_abundance(1:nspec,1))
   allocate(dummy_density(1))
   allocate(dummy_temperature(1))
    p=IDlist_dark(1)

!     Set the abundances to their appropriate starting values
      DUMMY_ABUNDANCE(:,1) = pdr(p)%abundance
      DUMMY_DENSITY(1) = pdr(P)%rho
      DUMMY_TEMPERATURE(1) = gastemperature(0)

  CALL CALCULATE_ABUNDANCES(DUMMY_ABUNDANCE,DUMMY_RATE,DUMMY_DENSITY,DUMMY_TEMPERATURE,1,NSPEC,NREAC)
   p=IDlist_dark(1)
   pdr(p)%abundance = dummy_abundance(:,1)
  deallocate(dummy_abundance)
  deallocate(dummy_density)
  deallocate(dummy_temperature)
  deallocate(dummy_rate)

   call calc_columndens
ENDDO

p=IDlist_dark(1)
!Calculate the partition functions
CALL CALCULATE_PARTITION_FUNCTION(CII_Z_FUNCTION,CII_NLEV,CII_ENERGIES,CII_WEIGHTS,GASTEMPERATURE(0))
CALL CALCULATE_PARTITION_FUNCTION(CI_Z_FUNCTION,CI_NLEV,CI_ENERGIES,CI_WEIGHTS,GASTEMPERATURE(0))
CALL CALCULATE_PARTITION_FUNCTION(OI_Z_FUNCTION,OI_NLEV,OI_ENERGIES,OI_WEIGHTS,GASTEMPERATURE(0))
CALL CALCULATE_PARTITION_FUNCTION(C12O_Z_FUNCTION,C12O_NLEV,C12O_ENERGIES,C12O_WEIGHTS,GASTEMPERATURE(0))
! Calculate the LTE level populations
CALL CALCULATE_LTE_POPULATIONS(CII_NLEV,PDR(P)%CII_POP,CII_ENERGIES,&
        &CII_WEIGHTS,CII_Z_FUNCTION,PDR(P)%ABUNDANCE(NCx)*N_H,GASTEMPERATURE(0))
CALL CALCULATE_LTE_POPULATIONS(CI_NLEV, PDR(P)%CI_POP, CI_ENERGIES, &
        &CI_WEIGHTS, CI_Z_FUNCTION, PDR(P)%ABUNDANCE(NC)*N_H,GASTEMPERATURE(0))
CALL CALCULATE_LTE_POPULATIONS(OI_NLEV, PDR(P)%OI_POP, OI_ENERGIES, &
        &OI_WEIGHTS, OI_Z_FUNCTION, PDR(P)%ABUNDANCE(NO)*N_H,GASTEMPERATURE(0))
CALL CALCULATE_LTE_POPULATIONS(C12O_NLEV,PDR(P)%C12O_POP,C12O_ENERGIES,&
        &C12O_WEIGHTS,C12O_Z_FUNCTION,PDR(P)%ABUNDANCE(NCO)*N_H,GASTEMPERATURE(0))

write(6,*) 'Assigning properties to all dark particles'
do pp=2,dark_ptot
   p=IDlist_dark(pp)
   pdr(p)%cii_pop = pdr(IDlist_dark(1))%cii_pop
   pdr(p)%ci_pop = pdr(IDlist_dark(1))%ci_pop
   pdr(p)%oi_pop = pdr(IDlist_dark(1))%oi_pop
   pdr(p)%c12o_pop = pdr(IDlist_dark(1))%c12o_pop
   pdr(p)%abundance = pdr(IDlist_dark(1))%abundance
enddo
write(6,*) 'Done! Proceeding with PDR...'

return
end subroutine dark_molecular_region


end Program
