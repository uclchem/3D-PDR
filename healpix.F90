
!#include "macros.h" 
 subroutine mk_xy2pix
    !=======================================================================
    !     sets the array giving the number of the pixel lying in (x,y)
    !     x and y are in {1,128}
    !     the pixel number is in {0,128**2-1}
    !
    !     if  i-1 = sum_p=0  b_p * 2^p
    !     then ix = sum_p=0  b_p * 4^p
    !          iy = 2*ix
    !     ix + iy in {0, 128**2 -1}
    !=======================================================================
    use definitions
    use healpix_types
    use healpix_module, only : x2pix, y2pix
    implicit none

!    REAL(KIND=DP), INTENT(INOUT) :: x2pix(0:1023),y2pix(0:1023)
    INTEGER(KIND=I4B):: kk,ipp,ii,jj,idd
    !=======================================================================

    do ii = 1,128           !for converting x,y into
       jj  = ii-1            !pixel numbers
       kk  = 0
       ipp = 1

       do
          if (jj==0) then
             x2pix(ii) = kk
             y2pix(ii) = 2*kk
             exit
          else
             idd = MODULO(JJ,2)
             jj  = jj/2
             kk  = ipp*idd+kk
             ipp = ipp*4
          endif
       enddo

    enddo

    RETURN
  END subroutine mk_xy2pix

! ============================================================================
  subroutine pix2vec_nest(nside, ipix, pix2x, pix2y, vector, vertex)
    use definitions
    use healpix_types
    use healpix_module, only : ns_max
    implicit none
    !=======================================================================
    !     renders vector (x,y,z) coordinates of the nominal pixel center
    !     for the pixel number ipix (NESTED scheme)
    !     given the map resolution parameter nside
    !     also returns the (x,y,z) position of the 4 pixel vertices (=corners)
    !     in the order N,W,S,E
    !=======================================================================
    INTEGER(KIND=I4B), INTENT(IN) :: nside, ipix
    REAL(KIND=DP), INTENT(OUT) :: vector(1:3)
    REAL(KIND=DP),     INTENT(OUT) :: vertex(1:3,1:4)

    integer(kind=i4b),intent(in)::pix2x(0:1023),pix2y(0:1023)

    INTEGER(KIND=I4B) :: npix, npface, &
         &     ipf, ip_low, ip_trunc, ip_med, ip_hi, &
         &     jrt, jr, nr, jpt, jp, kshift, nl4
    REAL(KIND=DP) :: z, fn, fact1, fact2, sth, phi

    INTEGER(KIND=I4B) ::  ix, iy, face_num
!     common /xy_nest/ ix, iy, face_num ! can be useful to calling routine

    ! coordinate of the lowest corner of each face
    INTEGER(KIND=I4B), dimension(1:12) :: jrll = (/ 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4 /) ! in unit of nside
    INTEGER(KIND=I4B), dimension(1:12) :: jpll = (/ 1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7 /) ! in unit of nside/2

    real(kind=DP) :: phi_nv, phi_wv, phi_sv, phi_ev, phi_up, phi_dn
    real(kind=DP) :: z_nv, z_sv, sth_nv, sth_sv
    real(kind=DP) :: hdelta_phi
    integer(kind=I4B) :: iphi_mod, iphi_rat
    logical(kind=LGT) :: do_vertex
    !-----------------------------------------------------------------------

ns_max=8192
      
    if (nside<1 .or. nside>ns_max) stop 'a'
    npix = 12 * nside**2
    if (ipix <0 .or. ipix>npix-1) stop 'b'

    !     initiates the array for the pixel number -> (x,y) mapping
    if (pix2x(1023) <= 0) call mk_pix2xy()

    fn = real(nside,kind=dp)
    fact1 = 1.0_dp/(3.0_dp*fn*fn)
    fact2 = 2.0_dp/(3.0_dp*fn)
    nl4   = 4*nside

    !     finds the face, and the number in the face
    npface = nside**2

    face_num = ipix/npface  ! face number in {0,11}
    ipf = MODULO(ipix,npface)  ! pixel number in the face {0,npface-1}

    !     finds the x,y on the face (starting from the lowest corner)
    !     from the pixel number
    ip_low = MODULO(ipf,1024)       ! content of the last 10 bits
    ip_trunc =   ipf/1024        ! truncation of the last 10 bits
    ip_med = MODULO(ip_trunc,1024)  ! content of the next 10 bits
    ip_hi  =     ip_trunc/1024   ! content of the high weight 10 bits

    ix = 1024*pix2x(ip_hi) + 32*pix2x(ip_med) + pix2x(ip_low)
    iy = 1024*pix2y(ip_hi) + 32*pix2y(ip_med) + pix2y(ip_low)

    !     transforms this in (horizontal, vertical) coordinates
    jrt = ix + iy  ! 'vertical' in {0,2*(nside-1)}
    jpt = ix - iy  ! 'horizontal' in {-nside+1,nside-1}

    !     computes the z coordinate on the sphere
    jr =  jrll(face_num+1)*nside - jrt - 1   ! ring number in    !deallocate(ra)
    !deallocate(rb)
    !deallocate(dis)

    !write(7,*) ifront(1:3),evalpoint(1:3,id)
    !write(4,*) ifront(1:3) {1,4*nside-1}

    nr = nside                  ! equatorial region (the most frequent)
    z  = (2*nside-jr)*fact2
    kshift = MODULO(jr - nside, 2)

    do_vertex=.false.

    if (do_vertex) then
       z_nv = (2*nside-jr+1)*fact2
       z_sv = (2*nside-jr-1)*fact2
       if (jr == nside) then ! northern transition
          z_nv =  1.0_dp - (nside-1)**2 * fact1
       elseif (jr == 3*nside) then  ! southern transition
          z_sv = -1.0_dp + (nside-1)**2 * fact1!pix2ang
       endif
    endif
    if (jr < nside) then     ! north pole region
       nr = jr
       z = 1.0_dp - nr*nr*fact1
       kshift = 0
       if (do_vertex) then
          z_nv = 1.0_dp - (nr-1)**2*fact1
          z_sv = 1.0_dp - (nr+1)**2*fact1
       endif
    else if (jr > 3*nside) then ! south pole region
       nr = nl4 - jr
       z = - 1.0_dp + nr*nr*fact1
       kshift = 0
       if (do_vertex) then
          z_nv = - 1.0_dp + (nr+1)**2*fact1
          z_sv = - 1.0_dp + (nr-1)**2*fact1
       endif
    endif
    !     computes the phi coordinate on the sphere, in [0,2Pi]
    jp = (jpll(face_num+1)*nr + jpt + 1 + kshift)/2  ! 'phi' number in the ring in {1,4*nr}
    if (jp > nl4) jp = jp - nl4
    if (jp < 1)   jp = jp + nl4

    phi = (jp - (kshift+1)*0.5_dp) * (halfpi / nr)

    sth = SQRT((1.0_dp-z)*(1.0_dp+z))
    vector(1) = sth * COS(phi)
    vector(2) = sth * SIN(phi)
    vector(3) = z

!two lines added by T.Bisbas due to a problem with precision.
if (abs(vector(2)).lt.1d-10) vector(2)=0.0_dp
if (abs(vector(1)).lt.1d-10) vector(1)=0.0_dp

    if (do_vertex) then
       phi_nv = phi
       phi_sv = phi

       phi_up = 0.0_dp
       iphi_mod = MODULO(jp-1, nr) ! in {0,1,... nr-1}
       iphi_rat = (jp-1) / nr      ! in {0,1,2,3}
       if (nr > 1) phi_up = HALFPI * (iphi_rat +  iphi_mod   /real(nr-1,kind=dp))
       phi_dn             = HALFPI * (iphi_rat + (iphi_mod+1)/real(nr+1,kind=dp))
       if (jr < nside) then            ! North polar cap
          phi_nv = phi_up
          phi_sv = phi_dn
       else if (jr > 3*nside) then     ! South polar cap
          phi_nv = phi_dn
          phi_sv = phi_up
       else if (jr == nside) then      ! North transition
          phi_nv = phi_up
       else if (jr == 3*nside) then    ! South transition
          phi_sv = phi_up
       endif

       hdelta_phi = PI / (4.0_dp*nr)

       ! west vertex
       phi_wv      = phi - hdelta_phi
       vertex(1,2) = sth * COS(phi_wv)
       vertex(2,2) = sth * SIN(phi_wv)
       vertex(3,2) = z

       ! east vertex
       phi_ev      = phi + hdelta_phi
       vertex(1,4) = sth * COS(phi_ev)
       vertex(2,4) = sth * SIN(phi_ev)
       vertex(3,4) = z

       ! north vertex
       sth_nv = SQRT((1.0_dp-z_nv)*(1.0_dp+z_nv))
       vertex(1,1) = sth_nv * COS(phi_nv)
       vertex(2,1) = sth_nv * SIN(phi_nv)
       vertex(3,1) = z_nv

       ! south vertex
       sth_sv = SQRT((1.0_dp-z_sv)*(1.0_dp+z_sv))
       vertex(1,3) = sth_sv * COS(phi_sv)
       vertex(2,3) = sth_sv * SIN(phi_sv)
       vertex(3,3) = z_sv
    endif

    return
  end subroutine pix2vec_nest
 
subroutine mk_pix2xy()
   use definitions
   use healpix_types
   use healpix_module
   implicit none
    !=======================================================================
    !     constructs the array giving x and y in the face from pixel number
    !     for the nested (quad-cube like) ordering of pixels    !
    !     the bits corresponding to x and y are interleaved in the pixel number
    !     one breaks up the pixel number by even and odd bits
    !=======================================================================
    INTEGER ::  kpix, jpix, ix, iy, ip, id

    !cc cf block data      data      pix2x(1023) /0/
    !-----------------------------------------------------------------------
    !      print *, 'initiate pix2xy'
    do kpix=0,1023          ! pixel number
       jpix = kpix
       IX = 0
       IY = 0
       IP = 1               ! bit position (in x and y)
!        do while (jpix/=0) ! go through all the bits
       do
          if (jpix == 0) exit ! go through all the bits
          ID = MODULO(jpix,2)  ! bit value (in kpix), goes in ix
          jpix = jpix/2
          IX = ID*IP+IX

          ID = MODULO(jpix,2)  ! bit value (in kpix), goes in iy
          jpix = jpix/2
          IY = ID*IP+IY

          IP = 2*IP         ! next bit (in x and y)
       enddo
       pix2x(kpix) = IX     ! in 0,31
       pix2y(kpix) = IY     ! in 0,31
    enddo

    return
  end subroutine mk_pix2xy

! ============================================================================
  subroutine vec2ang(rvec, theta, phi)
    !=======================================================================
    !     renders the angles theta, phi corresponding to vector (x,y,z) (rvec)
    !     theta (co-latitude measured from North pole, in [0,Pi] radians)
    !     and phi (longitude measured eastward, in [0,2Pi[ radians)
    !     North pole is (x,y,z)=(0,0,1)
    !     added by EH, Feb 2000
    !=======================================================================
    use definitions
    use healpix_types
    implicit none

    REAL(KIND=DP), INTENT(IN) :: rvec(1:3)
    REAL(KIND=DP), INTENT(OUT) :: theta, phi

    REAL(KIND=DP) :: dnorm, z
    !=======================================================================

    dnorm = SQRT(rvec(1)**2+rvec(2)**2+rvec(3)**2)
    if (dnorm.eq.0) then
        phi=0
        theta=0
        z=0
!       write(6,*) 'found dnorm = 0'
       return
    end if
    z = rvec(3) / dnorm
    theta = ACOS(z)

    phi = 0.0_dp
    if (rvec(1) /= 0.0_dp .or. rvec(2) /= 0.0_dp) &
         &     phi = ATAN2(rvec(2),rvec(1)) ! phi in ]-pi,pi]
    if (phi < 0.0)     phi = phi + twopi ! phi in [0,2pi[

    return
  end subroutine vec2ang

! ============================================================================
 subroutine ang2pix_nest_id(nside, theta, phi, ipix)
    !=======================================================================
    !     renders the pixel number ipix (NESTED scheme) for a pixel which contains
    !     a point on a sphere at coordinates theta and phi, given the map
    !     resolution parametr nside
    !
    !     the computation is made to the highest resolution available (nside=8192)
    !     and then degraded to that required (by integer division)
    !     this doesn't cost more, and it makes sure
    !     that the treatement of round-off will be consistent
    !     for every resolution
    !=======================================================================
    use definitions
    use healpix_types
    use healpix_module, only : x2pix, y2pix, ns_max
    implicit none

    INTEGER(KIND=I4B), INTENT(IN) :: nside
    INTEGER(KIND=I4B), INTENT(OUT) :: ipix
    REAL(KIND=DP), INTENT(IN) ::  theta, phi

    REAL(KIND=DP) ::  z, za, tt, tp, tmp
    INTEGER(KIND=I4B) :: jp, jm, ifp, ifm, face_num, &
         &     ix, iy, ix_low, ix_hi, iy_low, iy_hi, ipf, ntt

    !-----------------------------------------------------------------------
    if (x2pix(128) <= 0) call mk_xy2pix()

    z  = COS(theta)
    za = ABS(z)
    tt = MODULO(phi, twopi) / halfpi  ! in [0,4[

    if (za <= twothird) then ! equatorial region

       !        (the index of edge lines increase when the longitude=phi goes up)
       jp = INT(ns_max*(0.5_dp + tt - z*0.75_dp)) !  ascending edge line index
       jm = INT(ns_max*(0.5_dp + tt + z*0.75_dp)) ! descending edge line index

       !        finds the face
       if (ns_max.eq.0) stop 'insane ns_max @ _IF'
       ifp = jp / ns_max  ! in {0,4}
       ifm = jm / ns_max
       if (ifp == ifm) then          ! faces 4 to 7
          face_num = MODULO(ifp,4) + 4
       else if (ifp < ifm) then     ! (half-)faces 0 to 3
          face_num = MODULO(ifp,4)
       else                            ! (half-)faces 8 to 11
          face_num = MODULO(ifm,4) + 8
       endif

       ix = MODULO(jm, ns_max)
       iy = ns_max - MODULO(jp, ns_max) - 1

    else ! polar region, za > 2/3
        ntt = INT(tt)
       if (ntt >= 4) ntt = 3
       tp = tt - ntt
       tmp = SQRT( 3.0_dp*(1.0_dp - za) )  ! in ]0,1]

       !        (the index of edge lines increase when distance from the closest pole goes up)
       jp = INT( ns_max * tp * tmp ) ! line going toward the pole as phi increases
       jm = INT( ns_max * (1.0_dp - tp) * tmp ) ! that one goes away of the closest pole
       jp = MIN(ns_max-1, jp) ! for points too close to the boundary
       jm = MIN(ns_max-1, jm)

       !        finds the face and pixel's (x,y)
       if (z >= 0) then
          face_num = ntt  ! in {0,3}
          ix = ns_max - jm - 1
          iy = ns_max - jp - 1
       else
          face_num = ntt + 8 ! in {8,11}
          ix =  jp
          iy =  jm
       endif

       !         print*,z,face_num,ix,iy
    endif

    ix_low = MODULO(ix,128)
    ix_hi  =     ix/128
    iy_low = MODULO(iy,128)
    iy_hi  =     iy/128

    ipf =  (x2pix(ix_hi +1)+y2pix(iy_hi +1)) * (128 * 128) &
         &     + (x2pix(ix_low+1)+y2pix(iy_low+1))

    if (ns_max.eq.0.or.nside.eq.0) stop 'insane ns_max and/or nside @ _IF'
    !ipf = (ipf / ( ns_max/nside ))/(ns_max/nside)  ! in {0, nside**2 - 1}
    ipf = ipf / ( ns_max/nside ) **2
    ipix = ipf + face_num* nside **2    ! in {0, 12*nside**2 - 1}

    return
end subroutine ang2pix_nest_id
