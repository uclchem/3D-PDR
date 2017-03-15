
!#include "macros.h" 
! ============================================================================

!
! HEAPSORT.F90
! TAKEN FROM SEREN SPH CODE
! C. P. Batty - 16/5/2006
! Sorts a list of values (and its identifiers) using heapsort (bottom-up)
! ----------------------------------------------------------------------------

SUBROUTINE heapsort(pp_pot, plist, rlist)
  use definitions
  use healpix_types
  implicit none

  integer(kind=i4b), intent(in) :: pp_pot       ! number of potential neighbours found
  integer(kind=i4b), dimension(1:pp_pot+1), intent(inout) :: plist
                                      ! list of potential neighbours
  real(kind=dp), dimension(0:pp_pot+1), intent(inout) :: rlist
                                      ! list of drsqd values for reference
  integer :: start
  integer :: end

  integer :: sift_start
  integer :: sift_end  ! limit of how far down the heap to sift

  integer :: root
  integer :: child

  integer       :: swapi  ! storage for array element swap
  real(kind=PR) :: swapr  ! storage for array element swap

!  integer::pp_limit


  ! Place array(s) in max-heap order
  start = pp_pot / 2  ! assign index in array of last parent node
  do while (start > 0)
     ! Sift down the node at index start to place such that all nodes below
     ! start index are in heap order
     sift_start = start
     sift_end = pp_pot-1

     root = sift_start
     do
        child = root * 2  ! ...point to the left child
        if (child > sift_end) exit  ! if root has no children then exit
        ! If child has greater sibling then point to right child instead
        if ((child < sift_end) .and. (rlist(child) < rlist(child+1))) child = child + 1
        if (rlist(root) < rlist(child)) then  ! out of max-heap order
           swapr = rlist(root);  rlist(root) = rlist(child);  rlist(child) = swapr
           swapi = plist(root);  plist(root) = plist(child);  plist(child) = swapi
           root = child  ! repeat to continue sifting down child
        else
           exit
        end if
     end do

     start = start - 1
     ! After sifting down root, all elements/nodes are in heap order
  end do

  end = pp_pot
  do while (end > 1)
     ! Swap root of heap (maximum value) with last element of heap
     swapr = rlist(end);  rlist(end) = rlist(1);  rlist(1) = swapr
     swapi = plist(end);  plist(end) = plist(1);  plist(1) = swapi
     ! Decrease size of heap by one (previous max value stays in proper place)
     end = end - 1

     ! Put the heap back in max-heap order
     sift_start = 1
     sift_end = end

     root = sift_start
     do
        child = root * 2  ! ...point to the left child
        if (child > sift_end) exit  ! if root has no children then exit
        ! If child has greater sibling then point to right child instead
        if ((child < sift_end) .and. (rlist(child) < rlist(child+1))) child = child + 1
        if (rlist(root) < rlist(child)) then  ! out of max-heap order
           swapr = rlist(root);  rlist(root) = rlist(child);  rlist(child) = swapr
           swapi = plist(root);  plist(root) = plist(child);  plist(child) = swapi
           root = child  ! repeat to continue sifting down child
        else
           exit
        end if
     end do
  end do

  return
END SUBROUTINE heapsort

