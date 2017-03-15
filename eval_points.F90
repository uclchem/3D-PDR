subroutine evaluation_points
!calculation of evaluation points
!T.Bisbas
use definitions
use healpix_types
use healpix_module
use maincode_module
#ifdef OPENMP
use omp_lib
#endif

double precision::adaptivemin
logical::killray(0:nrays-1)
integer::j

if (dark_ptot.gt.0) then

write(6,*) 'Creating evaluation points for the Dark Molecular element'
!BUILDING HEALPIX VECTORS FOR THE ONE DARK MOLECULAR ELEMENT ======================================
!SERIAL PROCESS -----------------
!defines the origin to transfer all the domain in the original co-ordinate system
killray=.false.
origin(1:3) = pdrpoint(1:3,0)
allocate(ra(0:grand_ptot-1)) !needs one extra place for sorting in heapsort
allocate(rb(1:grand_ptot-1)) !-1 to avoid overlapping origin & pdrpoint
allocate(ep(1:3,0:nrays-1))
!calculating distances from the origin(1:3)
kk=0
pdr(IDlist_dark(1))%epray = 0
do i=1,grand_ptot
   if (i.eq.IDlist_dark(1)) cycle
   kk=kk+1
   !locates the grid point in the new computational domain
   rvec(1)=pdr(i)%x-origin(1)
   rvec(2)=pdr(i)%y-origin(2)
   rvec(3)=pdr(i)%z-origin(3)
   !next two lines return the ipix ray that the rvec(1:3) point belongs to.
   call vec2ang(rvec,theta,phi)
   call ang2pix_nest_id(nside,theta,phi,ipix)
   ra(kk)=sqrt(rvec(1)**2+rvec(2)**2+rvec(3)**2)
   rb(kk)=i !stores the identifier of each grid point
enddo
ktot=kk !ktot should be grand_ptot-1
if (ktot.ne.(grand_ptot-1)) then 
   write(6,*) 'ktot = ',ktot,' grand_ptot-1 = ',grand_ptot-1
   stop 'ktot is not equal to grand_ptot-1 !!'
endif
!calling heapsort and sorting with increasing the distance from the origin(1:3)
call heapsort(ktot,rb,ra)
!maximum distance from origin(1:3). This is the radius at 
!which the HEALPix vectors should expand.
radius=ra(ktot)
!gives values for the first evaluation point which is the origin (1:3)
ep=0.
do j=0,nrays-1
   pdr(IDlist_dark(1))%epoint(1:3,j,0) = origin(1:3)
end do
!loops over all the domain and finds evaluation points. [straight N loop]
do k=1,ktot
   !locates the grid point in the new computational domain
   rvec(1)=pdr(rb(k))%x-origin(1)
   rvec(2)=pdr(rb(k))%y-origin(2)
   rvec(3)=pdr(rb(k))%z-origin(3)
   !next two lines return the ipix ray that the rvec(1:3) point belongs to.
   call vec2ang(rvec,theta,phi)
   call ang2pix_nest_id(nside,theta,phi,ipix)
   if (killray(ipix)) cycle
   healpixvector(1:3) = 1.1_DP*radius*vectors(1:3,ipix) !expand unit healpix vectors
   !calculates the angle along the line of sight of EVALUATION POINT -- HEALPIXVECTOR
   angle_los=acos(dot_product(rvec(1:3)-ep(1:3,ipix),healpixvector(1:3)-ep(1:3,ipix))/&
       &(sqrt((rvec(1)-ep(1,ipix))**2+(rvec(2)-ep(2,ipix))**2+(rvec(3)-ep(3,ipix))**2) * &
       &sqrt((healpixvector(1)-ep(1,ipix))**2+(healpixvector(2)-ep(2,ipix))**2+&
       &(healpixvector(3)-ep(3,ipix))**2)))
   !if the angle is less than the critical theta (user defined), then we have a new
   !evaluation point. This point is the projection of the grid point in the above line of sight.
   !if ((angle_los.le.theta_crit).and.(angle_los.ge.0D0)) then 
   if (angle_los.le.theta_crit) then
      !All next if-statements are conditions to avoid division by zero (i.e. x-plane, y-plane, z-plane)
      if (healpixvector(3).ne.0.0_dp) then
         ep(3,ipix) = (healpixvector(1)*healpixvector(3)*rvec(1) + healpixvector(2)*healpixvector(3)*&
                  &rvec(2) + (healpixvector(3)**2)*rvec(3))/(healpixvector(1)**2+healpixvector(2)**2+&
                  &healpixvector(3)**2)
         ep(2,ipix) = ep(3,ipix)*healpixvector(2)/healpixvector(3)
         ep(1,ipix) = ep(3,ipix)*healpixvector(1)/healpixvector(3)
      else
         if (healpixvector(1).eq.0.0_dp) then
            ep(1,ipix) = 0.0_dp
            ep(2,ipix) = rvec(2)
            ep(3,ipix) = 0.0_dp
         else if (healpixvector(2).eq.0.0_dp) then
            ep(1,ipix) = rvec(1)
            ep(2,ipix) = 0.0_dp
            ep(3,ipix) = 0.0_dp
         else
            ep(3,ipix) = 0.0_dp
            ep(1,ipix) = ((healpixvector(1)**2)*rvec(1) + healpixvector(1)*healpixvector(2)*rvec(2))/&
                  &(healpixvector(1)**2 + healpixvector(2)**2)
            ep(2,ipix) = ep(1,ipix) * healpixvector(2)/healpixvector(1)
         endif
      endif !healpixvector(3)
    !updates memory and stores the evaluation point in the original computational domain (so evaluation point + origin)
    pdr(IDlist_dark(1))%epray(ipix) = pdr(IDlist_dark(1))%epray(ipix)+1
    id = pdr(IDlist_dark(1))%epray(ipix)
    if (pdr(IDlist_dark(1))%epray(ipix).gt.maxpoints) STOP 'Increase maxpoints!'
    pdr(IDlist_dark(1))%epoint(1:3,ipix,id)=ep(1:3,ipix)+origin(1:3)
    pdr(IDlist_dark(1))%projected(ipix,id)=rb(k)
    if (pdr(rb(k))%etype.eq.2) killray(ipix)=.true. !if the projected is ionized, stop propagating the ray (it has hit the HII region)
    endif !angle_los
enddo !k=1,ktot
deallocate(ra)
deallocate(rb)
deallocate(ep)
!==========================================================================================
endif
#ifdef OPENMP
!$OMP PARALLEL
!$OMP MASTER
CPUs = OMP_GET_NUM_THREADS()
write(6,*) "Proceeding for the PDR (PARALLEL)..."
write(6,*) "Number of CPUs: ", OMP_GET_NUM_THREADS()
!$OMP END MASTER
!$OMP END PARALLEL
#else
write(6,*) 'Proceeding for the PDR (SERIAL)...'
#endif

!BUILDING HEALPIX VECTORS FOR ALL PDR ELEMENTS. 
!PARALLEL PROCESS-----------------------------
#ifdef OPENMP
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(origin, p, ra, rb, ep) &
!$OMP PRIVATE(i, rvec, theta, phi, ipix, ktot, radius) &
!$OMP PRIVATE(j, healpixvector, angle_los, id, killray) REDUCTION (+ : kk)
#endif 
do p=1,pdr_ptot
    !defines the origin to transfer all the domain in the original co-ordinate system
    killray=.false.
    origin(1:3) = pdrpoint(1:3,p)
    allocate(ra(0:grand_ptot-1)) !needs one extra place for sorting in heapsort
    allocate(rb(1:grand_ptot-1)) !-1 to avoid overlapping origin & pdrpoint
    allocate(ep(1:3,0:nrays-1))

    !calculating distances from the origin(1:3)
    kk=0
    pdr(IDlist_pdr(p))%epray = 0

    do i=1,grand_ptot
      if (i.eq.IDlist_pdr(p)) cycle
      kk=kk+1
      !locates the grid point in the new computational domain
      rvec(1)=pdr(i)%x-origin(1)
      rvec(2)=pdr(i)%y-origin(2)
      rvec(3)=pdr(i)%z-origin(3)
      !next two lines return the ipix ray that the rvec(1:3) point belongs to.
      call vec2ang(rvec,theta,phi)
      call ang2pix_nest_id(nside,theta,phi,ipix)
      ra(kk)=sqrt(rvec(1)**2+rvec(2)**2+rvec(3)**2)
      rb(kk)=i !stores the identifier of each grid point
    enddo
    ktot=kk !ktot should be grand_ptot-1
    if (ktot.ne.(grand_ptot-1)) then 
      write(6,*) 'ktot = ',ktot,' grand_ptot-1 = ',grand_ptot-1
      stop 'ktot is not equal to grand_ptot-1 !!'
    endif

    !calling heapsort and sorting with increasing the distance from the origin(1:3)
    call heapsort(ktot,rb,ra)
    !maximum distance from origin(1:3). This is the radius at 
    !which the HEALPix vectors should expand.
    radius=ra(ktot)

    !gives values for the first evaluation point which is the origin (1:3)
    ep=0.
    do j=0,nrays-1
      pdr(IDlist_pdr(p))%epoint(1:3,j,0) = origin(1:3)
    end do
    !loops over all the domain and finds evaluation points. [straight N loop]
    do k=1,ktot
      !locates the grid point in the new computational domain
      rvec(1)=pdr(rb(k))%x-origin(1)
      rvec(2)=pdr(rb(k))%y-origin(2)
      rvec(3)=pdr(rb(k))%z-origin(3)
      !next two lines return the ipix ray that the rvec(1:3) point belongs to.
      call vec2ang(rvec,theta,phi)
      call ang2pix_nest_id(nside,theta,phi,ipix)
      if (killray(ipix)) cycle
      healpixvector(1:3) = 1.1_DP*radius*vectors(1:3,ipix) !expand unit healpix vectors
      !calculates the angle along the line of sight of EVALUATION POINT -- HEALPIXVECTOR
      angle_los=acos(dot_product(rvec(1:3)-ep(1:3,ipix),healpixvector(1:3)-ep(1:3,ipix))/&
          &(sqrt((rvec(1)-ep(1,ipix))**2+(rvec(2)-ep(2,ipix))**2+(rvec(3)-ep(3,ipix))**2) * &
          &sqrt((healpixvector(1)-ep(1,ipix))**2+(healpixvector(2)-ep(2,ipix))**2+&
          &(healpixvector(3)-ep(3,ipix))**2)))
      !if the angle is less than the critical theta (user defined), then we have a new
      !evaluation point. This point is the projection of the grid point in the above line of sight.
      !if ((angle_los.le.theta_crit).and.(angle_los.ge.0D0)) then 
      if (angle_los.le.theta_crit) then
         !All next if-statements are conditions to avoid division by zero (i.e. x-plane, y-plane, z-plane)
         if (healpixvector(3).ne.0.0_dp) then
            ep(3,ipix) = (healpixvector(1)*healpixvector(3)*rvec(1) + healpixvector(2)*healpixvector(3)*&
                         &rvec(2) + (healpixvector(3)**2)*rvec(3))/(healpixvector(1)**2+healpixvector(2)**2+&
                         &healpixvector(3)**2)
            ep(2,ipix) = ep(3,ipix)*healpixvector(2)/healpixvector(3)
            ep(1,ipix) = ep(3,ipix)*healpixvector(1)/healpixvector(3)
         else
            if (healpixvector(1).eq.0.0_dp) then
               ep(1,ipix) = 0.0_dp
               ep(2,ipix) = rvec(2)
               ep(3,ipix) = 0.0_dp
            else if (healpixvector(2).eq.0.0_dp) then
               ep(1,ipix) = rvec(1)
               ep(2,ipix) = 0.0_dp
               ep(3,ipix) = 0.0_dp
            else
               ep(3,ipix) = 0.0_dp
               ep(1,ipix) = ((healpixvector(1)**2)*rvec(1) + healpixvector(1)*healpixvector(2)*rvec(2))/&
                            &(healpixvector(1)**2 + healpixvector(2)**2)
               ep(2,ipix) = ep(1,ipix) * healpixvector(2)/healpixvector(1)
            endif
         endif !healpixvector(3)
!       !updates memory and stores the evaluation point in the original computational domain (so evaluation point + origin)
       pdr(IDlist_pdr(p))%epray(ipix) = pdr(IDlist_pdr(p))%epray(ipix)+1
       id = pdr(IDlist_pdr(p))%epray(ipix)
       if (pdr(IDlist_pdr(p))%epray(ipix).gt.maxpoints) STOP 'Increase maxpoints!'
       pdr(IDlist_pdr(p))%epoint(1:3,ipix,id)=ep(1:3,ipix)+origin(1:3)
       pdr(IDlist_pdr(p))%projected(ipix,id)=rb(k)
       if (pdr(rb(k))%etype.eq.2) killray(ipix)=.true.
       endif !angle_los
     enddo !k=1,ktot
    deallocate(ra)
    deallocate(rb)
    deallocate(ep)
enddo !pp/p=1,pdr_ptot
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif


suma=0
do pp=1,pdr_ptot
  p=IDlist_pdr(pp)
  suma = suma + sum(pdr(p)%epray(:))
enddo
if (dark_ptot.gt.0) then
!Include the Dark Molecular element
suma = suma + sum(pdr(IDlist_dark(1))%epray(:))
endif
write(6,*) 'No. evaluation points:',suma
write(6,*) 'Done!';write(6,*) ''

write(6,*) 'Checking for negative steps...'
adaptivemin=100.0D0
do pp=1,pdr_ptot
   p=IDlist_pdr(pp)
  do j=0,nrays-1
    if (pdr(p)%epray(j).gt.0) then
    do i=1,pdr(p)%epray(j)
       adaptive_step = sqrt((pdr(p)%epoint(1,j,0)-pdr(p)%epoint(1,j,i))**2+&
                           &(pdr(p)%epoint(2,j,0)-pdr(p)%epoint(2,j,i))**2+&
                           &(pdr(p)%epoint(3,j,0)-pdr(p)%epoint(3,j,i))**2)-&
                      &sqrt((pdr(p)%epoint(1,j,0)-pdr(p)%epoint(1,j,i-1))**2+&
                           &(pdr(p)%epoint(2,j,0)-pdr(p)%epoint(2,j,i-1))**2+&
                           &(pdr(p)%epoint(3,j,0)-pdr(p)%epoint(3,j,i-1))**2)
       if (adaptive_step.lt.0) stop 'found negative adaptive step!'
       if (adaptive_step.lt.adaptivemin) adaptivemin = adaptive_step
    enddo
    endif
   enddo
enddo
if (dark_ptot.gt.0) then
!Checking for the Dark Molecular element
p=IDlist_dark(1)
do j=0,nrays-1
  if (pdr(p)%epray(j).gt.0) then
     do i=1,pdr(p)%epray(j)
        adaptive_step = sqrt((pdr(p)%epoint(1,j,0)-pdr(p)%epoint(1,j,i))**2+&
                         &(pdr(p)%epoint(2,j,0)-pdr(p)%epoint(2,j,i))**2+&
                         &(pdr(p)%epoint(3,j,0)-pdr(p)%epoint(3,j,i))**2)-&
                    &sqrt((pdr(p)%epoint(1,j,0)-pdr(p)%epoint(1,j,i-1))**2+&
                         &(pdr(p)%epoint(2,j,0)-pdr(p)%epoint(2,j,i-1))**2+&
                         &(pdr(p)%epoint(3,j,0)-pdr(p)%epoint(3,j,i-1))**2)
        if (adaptive_step.lt.0) stop 'found negative adaptive step!'
        if (adaptive_step.lt.adaptivemin) adaptivemin = adaptive_step
     enddo
  endif
enddo
endif

write(6,*) 'No negative steps found'
write(6,*) 'Minimum adaptive step = ',adaptivemin


write(6,*) 'Assigning raytypes'
do pp=1,pdr_ptot
   p=IDlist_pdr(pp)
   allocate(pdr(p)%raytype(0:nrays-1))
   do j=0,nrays-1
     if (pdr(p)%epray(j).gt.0) then
       pdr(p)%raytype(j) = -pdr(pdr(p)%projected(j,pdr(p)%epray(j)))%etype
    else
       pdr(p)%raytype(j) = -pdr(p)%etype
     endif
  enddo
enddo

!Assigning raytype for the Dark Molecular element
if (dark_ptot.gt.0) then
  p=IDlist_dark(1)
  allocate(pdr(p)%raytype(0:nrays-1))
  do j=0,nrays-1
    if (pdr(p)%epray(j).gt.0) then
       pdr(p)%raytype(j) = -pdr(pdr(p)%projected(j,pdr(p)%epray(j)))%etype
    else
       pdr(p)%raytype(j) = -pdr(p)%etype
    endif
  enddo
endif

return
end subroutine evaluation_points
