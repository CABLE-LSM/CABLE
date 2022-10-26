!jhan:Althoughthis isonly calling subrs and not using data - still needs revision to only call once and use L_tile_pts from here
module cbl_masks_mod
  public L_tile_pts
  public fveg_mask
  public fsunlit_mask
  public fsunlit_veg_mask
  public veg_mask
  public sunlit_mask
  public sunlit_veg_mask
!H! remove SAVE attr later 
  !mask TRUE where tile fraction is greater than zero
  LOGICAL, allocatable,SAVE :: L_tile_pts(:,:)
  LOGICAL, allocatable, SAVE :: veg_mask(:), sunlit_mask(:), sunlit_veg_mask(:)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fveg_mask( veg_mask, mp, lai_thresh, reducedLAIdue2snow )

implicit none
LOGICAL, allocatable :: veg_mask(:)
integer ::  mp
real :: lai_thresh
real :: reducedLAIdue2snow(mp)
integer :: i
 
IF ( .NOT. ALLOCATED(veg_mask)) ALLOCATE( veg_mask(mp) )

veg_mask = .FALSE.
! Define vegetation mask:
do i=1, mp  
  if( reducedLAIdue2snow(i) > lai_thresh ) veg_mask(i) = .true.
end do

End subroutine fveg_mask

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fsunlit_mask( sunlit_mask, mp, coszen_tols, coszen )

implicit none
LOGICAL, allocatable :: sunlit_mask(:)
integer ::  mp
real :: coszen_tols
real :: coszen(mp)   
integer :: i

IF ( .NOT. ALLOCATED(sunlit_mask)) ALLOCATE( sunlit_mask(mp) )

sunlit_mask = .FALSE.
! Define sunlit mask:
do i=1, mp  
  if( coszen(i) > coszen_tols ) sunlit_mask(i) = .true.  
end do

End subroutine fsunlit_mask

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fsunlit_veg_mask( sunlit_veg_mask, mp )

implicit none
LOGICAL, allocatable :: sunlit_veg_mask(:)
integer ::  mp
integer :: i

IF ( .NOT. ALLOCATED(sunlit_veg_mask) ) ALLOCATE( sunlit_veg_mask(mp) )

sunlit_veg_mask = .FALSE.
! Define sunlit AND vegetation mask:
do i=1, mp  
  if( veg_mask(i) .AND.  sunlit_mask(i) ) sunlit_veg_mask(i) = .true.
end do

End subroutine fsunlit_veg_mask

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


End module cbl_masks_mod
