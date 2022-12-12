module cbl_masks_mod
  public fveg_mask
  public fsunlit_mask
  public fsunlit_veg_mask

contains

!=============================================================================

subroutine fveg_mask( veg_mask, mp, lai_thresh, reducedLAIdue2snow )

implicit none
integer ::  mp
LOGICAL :: veg_mask(mp)
real :: lai_thresh
real :: reducedLAIdue2snow(mp)
integer :: i
 
veg_mask = .FALSE.
! Define vegetation mask:
do i=1, mp  
  if( reducedLAIdue2snow(i) > lai_thresh ) veg_mask(i) = .true.
end do

End subroutine fveg_mask

!=============================================================================

subroutine fsunlit_mask( sunlit_mask, mp, coszen_tols, coszen )

implicit none
integer ::  mp
LOGICAL :: sunlit_mask(mp)
real :: coszen_tols
real :: coszen(mp)   
integer :: i

sunlit_mask = .FALSE.
! Define sunlit mask:
do i=1, mp  
  if( coszen(i) > coszen_tols ) sunlit_mask(i) = .true.  
end do

End subroutine fsunlit_mask

!=============================================================================

subroutine fsunlit_veg_mask( sunlit_veg_mask, veg_mask, sunlit_mask, mp )

implicit none
integer ::  mp
LOGICAL :: veg_mask(mp)
LOGICAL :: sunlit_mask(mp)
LOGICAL :: sunlit_veg_mask(mp)
integer :: i

sunlit_veg_mask = .FALSE.
! Define sunlit AND vegetation mask:
do i=1, mp  
  if( veg_mask(i) .AND.  sunlit_mask(i) ) sunlit_veg_mask(i) = .true.
end do

End subroutine fsunlit_veg_mask

!=============================================================================


End module cbl_masks_mod
