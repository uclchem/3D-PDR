!T.Bisbas, T.Bell

subroutine escape_probability(transition, dust_temperature, nrays, nlev, &
                   &A_COEFFS, B_COEFFS, C_COEFFS, &
                   &frequencies,s_evalpop, maxpoints, Tguess, v_turb,&
                   &s_jjr, s_pop, s_evalpoint, weights,cooling_rate,line,tau,coolant,density,metallicity,bbeta)


use definitions
use maincode_module, only : p,pdr,vectors
use healpix_types
use healpix_module
use global_module, only: g2d

implicit none

integer(kind=i4b), intent(in) :: nrays
integer(kind=i4b), intent(in) :: nlev
integer(kind=i4b), intent(in) :: maxpoints
integer(kind=i4b), intent(in) :: s_jjr(0:nrays-1)
integer(kind=i4b), intent(in) :: coolant
real(kind=dp), intent(inout) :: transition(1:nlev,1:nlev)
real(kind=dp), intent(in) :: A_COEFFS(1:nlev, 1:nlev)
real(kind=dp), intent(in) :: B_COEFFS(1:nlev, 1:nlev)
real(kind=dp), intent(in) :: C_COEFFS(1:nlev, 1:nlev)
real(kind=dp), intent(in) :: frequencies(1:nlev, 1:nlev)
real(kind=dp), intent(in) :: s_evalpop(0:nrays-1,0:maxpoints,1:nlev)
real(kind=dp), intent(in) :: s_evalpoint(1:3,0:nrays-1,0:maxpoints)
real(kind=dp), intent(in) :: Tguess, v_turb
real(kind=dp), intent(in) :: weights(1:nlev)
real(kind=dp), intent(in) :: s_pop(1:nlev)
real(kind=dp), intent(in) :: dust_temperature,density,metallicity

integer(kind=i4b) :: i, j
integer(kind=i4b) :: ilevel, jlevel
real(kind=dp) :: beta_ij, beta_ij_sum
real(kind=dp) :: frac1, frac2, frac3, rhs2
real(kind=dp) :: tpop, tmp2
real(kind=dp) :: S_ij, BB_ij
real(kind=dp) :: tau_increment
real(kind=dp), allocatable :: tau_ij(:)
real(kind=dp), allocatable :: field(:,:)
real(kind=dp) :: beta_ij_ray(0:nrays-1)
real(kind=dp), intent(out) :: line(1:nlev,1:nlev)
real(kind=dp), intent(out) :: cooling_rate
real(kind=dp), intent(out) :: tau(1:nlev,1:nlev,0:nrays-1)
real(kind=dp),intent(out) :: bbeta(1:nlev,1:nlev,0:nrays-1)
real(kind=dp) :: emissivity, bb_ij_dust, ngrain, rho_grain


line=0.0D0
cooling_rate = 0.0D0
    allocate(tau_ij(0:nrays-1))
    allocate(field(1:nlev,1:nlev))
    field=0.0D0
    frac2=1.0D0/sqrt(8.0*KB*Tguess/PI/MP + v_turb**2)
    !frac2=1.0D0/sqrt(KB*Tguess/PI/MP+v_turb**2/2.)/sqrt(2.*PI)
    do ilevel=1,nlev
       do jlevel=1,nlev !i>j
         if (jlevel.ge.ilevel) exit
         tau_ij=0.0D0; beta_ij=0.0D0; beta_ij_ray=0.0D0; beta_ij_sum=0.0D0
         frac1=(A_COEFFS(ilevel,jlevel)*(C**3))/(8.0*pi*(frequencies(ilevel,jlevel)**3))
         TMP2=2.0D0*HP*(FREQUENCIES(ilevel,jlevel)**3)/(C**2)

         !Planck function !2.7D0 is the CMBR temperature
         BB_ij = TMP2*(1.0D0/(EXP(HP*frequencies(ilevel,jlevel)/KB/2.7D0)-1.0D0))
         NGRAIN=2.0D-12*density*metallicity*100./g2d!densityofgas depth depented
         rho_grain=2.0D0
         EMISSIVITY=(RHO_GRAIN*NGRAIN)*(0.01*(1.3*FREQUENCIES(ilevel,jlevel)/3.0D11))
         BB_ij_dust = TMP2*(1.0D0/(EXP(HP*frequencies(ilevel,jlevel)/KB/DUST_TEMPERATURE)-1.D0)*EMISSIVITY)
         BB_ij = BB_ij + BB_ij_dust
         if (s_pop(ilevel).eq.0) then
            S_ij=0.0D0
            beta_ij=1.0D0
            goto 2
         endif
         TPOP=(s_pop(jlevel)*WEIGHTS(ilevel))/(s_pop(ilevel)*WEIGHTS(jlevel))-1.0D0
         IF(abs(TPOP).lt.1.0D-50) then
              S_ij=HP*FREQUENCIES(ilevel,jlevel)*s_pop(ilevel)*A_COEFFS(ilevel,jlevel)/4./pi
              beta_ij=1.0D0
              goto 1
         else
         !calculation of source function (taken from UCL_PDR)
              S_ij=TMP2/TPOP
         endif
         do j=0,nrays-1
#ifdef PSEUDO_1D
         if (j.ne.6) then
           tau_ij(j) = 1.0D50
         else
#endif
#ifdef PSEUDO_2D
         if (abs(vectors(3,j).gt.1d-10) tau_ij(j) = 1.0D50 !Not in Equator
#endif


 do i=1,s_jjr(j)
             !calculations of tau_ij
     frac3=((s_evalpop(j,i-1,jlevel)*weights(ilevel)-s_evalpop(j,i-1,ilevel)*weights(jlevel))+&
      &(s_evalpop(j,i,jlevel)*weights(ilevel)-s_evalpop(j,i,ilevel)*weights(jlevel)))/2./weights(jlevel)
     rhs2=sqrt((s_evalpoint(1,j,i-1)-s_evalpoint(1,j,i))**2+&
              &(s_evalpoint(2,j,i-1)-s_evalpoint(2,j,i))**2+&
              &(s_evalpoint(3,j,i-1)-s_evalpoint(3,j,i))**2) !adaptive step
     tau_increment=frac1*frac2*frac3*rhs2*PC
     tau_ij(j)=tau_ij(j)+tau_increment !optical depth
 enddo !i=1,jr(j)
#ifdef PSEUDO_1D
         endif
#endif
#ifdef PSEUDO_2D
         endif
#endif


           ! Prevent exploding beta values caused by strong masing (tau < -10)
           ! Assume tau = -10 and calculate the escape probability accordingly
           if (tau_ij(j).lt.-5.0D0) then
              beta_ij_ray(j)=(1.0D0-EXP(5.0D0))/(-5.0D0)
!           ! Treat weak masing using the standard escape probability formalism
!           else if (tau_ij(j).lt.0.0D0) then
!              beta_ij_ray(j)=(1.0D0-EXP(-tau_ij(j)))/tau_ij(j)
           ! Prevent floating point overflow caused by very low opacity (tau < 1e-6)
           else if (abs(tau_ij(j)).lt.1.0D-8) then !was D-6
              beta_ij_ray(j)=1.0D0
           ! For all other cases use the standard escape probability formalism
           else
              beta_ij_ray(j)=(1.0D0-EXP(-tau_ij(j)))/tau_ij(j)
           endif

        !=============
        tau(ilevel,jlevel,j)=tau_ij(j)
        bbeta(ilevel,jlevel,j)=beta_ij_ray(j)
        !=============

         enddo !j=0,nrays-1
         beta_ij_sum=sum(beta_ij_ray)
         !calculation of average beta_ij in the origin grid point
#ifdef PSEUDO_1D
         beta_ij = beta_ij_sum
#elif PSEUDO_2D
         beta_ij = beta_ij_sum / 4.
#else
         beta_ij = beta_ij_sum / real(nrays,kind=DP) 
#endif

1 continue
         line(ilevel,jlevel) = A_COEFFS(ilevel,jlevel)*HP*frequencies(ilevel,jlevel) * &
                             & s_pop(ilevel)*beta_ij*(S_ij-BB_ij)/S_ij
         cooling_rate = cooling_rate + line(ilevel,jlevel)
2 continue
         !<J_ij>
         field(ilevel,jlevel) = (1.0D0-beta_ij)*S_ij + beta_ij*BB_ij
         field(jlevel,ilevel) = field(ilevel,jlevel)
       enddo !jlevel=1,nlev
     enddo !ilevel=1,nlev
 

    !R_IJ CALCULATIONS
    !Update the transition matrix: Rij = Aij + Bij.<J> + Cij				    	 
    DO ilevel=1,NLEV
      DO jlevel=1,NLEV
        TRANSITION(ilevel,jlevel)=A_COEFFS(ilevel,jlevel)&
        & +B_COEFFS(ilevel,jlevel)*FIELD(ilevel,jlevel)&
        & +C_COEFFS(ilevel,jlevel)
        IF(ABS(TRANSITION(ilevel,jlevel)).LT.1.0D-50) TRANSITION(ilevel,jlevel)=0.0D0
      ENDDO !jlevel=1,nlev
    ENDDO !ilevel=1,nlev

    deallocate(tau_ij)
    deallocate(field)

  return

end subroutine escape_probability
