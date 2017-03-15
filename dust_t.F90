!=======================================================================
!
!  Calculate the dust temperature for each particle using the treatment
!  of Hollenbach, Takahashi & Tielens (1991, ApJ, 377, 192, eqns 5 & 6)
!  for the heating due to the incident FUV photons and the treatment of
!  Meijerink & Spaans (2005, A&A, 436, 397, eqn B.6) for heating due to
!  the incident flux of X-ray photons.
!
!  Among other things, the dust temperature can influence:
!
!     1) Cooling budget by emitting FIR photons that
!        interact with the line radiative transfer;
!     2) Gas-grain collisional heating or cooling rate;
!     3) H2 formation by changing the sticking probability;
!     4) Evaporation and condensation of molecules on grains.
!
!  The formula derived by Hollenbach, Takahashi & Tielens (1991) has
!  been modified to include the attenuation of the IR radiation. The
!  incident FUV radiation is absorbed and re-emitted in the infrared
!  by dust at the surface of the cloud (up to Av ~ 1mag). In the HTT
!  derivation, this IR radiation then serves as a second heat source
!  for dust deeper into the cloud. However, in their treatment, this
!  second re-radiated component is not attenuated with distance into
!  the cloud so it is *undiluted* with depth, leading to higher dust
!  temperatures deep within the cloud which in turn heat the gas via
!  collisions to unrealistically high temperatures. Models with high
!  gas densities and high incident FUV fluxes (e.g. n_H = 10^5 cm-3,
!  X_0 = 10^8 Draine) can produce T_gas ~ 100 K at Av ~ 50 mag!
!
!  Attenuation of the FIR radiation has therefore been introduced by
!  using an approximation for the infrared-only dust temperature from
!  Rowan-Robinson (1980, eqn 30b):
!
!  T_dust = T_0*(r/r_0)^(-0.4)
!
!  where r_0 is the cloud depth at which T_dust = T_0, corresponding
!  to an A_V of ~ 1 mag, the assumed size of the outer region of the
!  cloud that processes the incident FUV radiation and then re-emits
!  it in the FIR (see the original HTT 1991 paper for details). This
!  should prevent the dust temperature from dropping off too rapidly
!  with distance and maintain a larger warm dust region (~50-100 K).
!
!-----------------------------------------------------------------------
SUBROUTINE CALCULATE_DUST_TEMPERATURES

   USE HEALPIX_TYPES
   USE MAINCODE_MODULE

   IMPLICIT NONE

   INTEGER(KIND=I4B) :: J

   REAL(KIND=DP) :: NU_0,R_0,T_0,TAU_100
   REAL(KIND=DP) :: T_CMB

!  Parameters used in the HHT equations (see their paper for details)
   NU_0=2.65D15
   TAU_100=1.0D-3
   R_0=1.0D0/AV_FAC
   T_CMB=2.73D0

   DO PP=1,pdr_ptot ! Loop over particles
      P=IDlist_pdr(pp)
!     Calculate the contribution to the dust temperature from the local FUV flux and the CMB background
      PDR(P)%DUST_T=8.9D-11*NU_0*(1.71D0*PDR(P)%UVfield)+T_CMB**5

      DO J=0,NRAYS-1 ! Loop over rays

!        The minimum dust temperature is related to the incident FUV flux along each ray
!        Convert the incident FUV flux from Draine to Habing units by multiplying by 1.7
         T_0=12.2*(1.71D0*PDR(P)%RAD_SURFACE(J))**0.2

!        Add the contribution to the dust temperature from the FUV flux incident along this ray
         IF(T_0.GT.0) PDR(P)%DUST_T=PDR(P)%DUST_T &
                      & + (0.42-LOG(3.45D-2*TAU_100*T_0))*(3.45D-2*TAU_100*T_0)*T_0**5

      END DO ! End of loop over rays

!     Convert from total dust emission intensity to dust temperature
      PDR(P)%DUST_T=PDR(P)%DUST_T**0.2

!     Impose a lower limit on the dust temperature, since values below 10 K can dramatically
!     limit the rate of H2 formation on grains (the molecule cannot desorb from the surface)
      IF(PDR(P)%DUST_T.LT.10.0D0) THEN
         PDR(P)%DUST_T=10.0D0
      END IF

!     Check that the dust temperature is physical
      IF(PDR(P)%DUST_T.GT.1000) THEN
         WRITE(6,*) 'ERROR! Calculated dust temperature exceeds 1000 K'
         STOP
      END IF

   END DO ! End of loop over particles

   RETURN
END SUBROUTINE CALCULATE_DUST_TEMPERATURES
!=======================================================================
