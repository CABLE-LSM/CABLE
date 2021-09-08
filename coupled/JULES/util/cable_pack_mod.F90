MODULE cable_pack_mod

CONTAINS

SUBROUTINE cable_pack_rr(jules_var, cable_var)

USE cable_types_mod,          ONLY: mp, l_tile_pts
USE ancil_info,               ONLY: row_length, rows, land_index, land_pts,   &
                                    surft_pts, surft_index
USE jules_surface_types_mod,  ONLY: ntype

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   JULES met forcing vars needed by CABLE commonly have JULES dimensions
!   (row_length,rows), which are no good for CABLE. These have to be
!   re-packed in a single vector of active tiles. Hence we use
!   conditional "mask" l_tile_pts(land_pts,ntiles) which is .true.
!   if the land point is/has an active tile
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!
!-----------------------------------------------------------------------------


REAL, INTENT(IN) :: jules_var(row_length, rows)
REAL, INTENT(INOUT) :: cable_var(mp)
REAL :: fvar(land_pts, ntype)
INTEGER :: n, k, l, j, i

fvar(:, :) = 0.0
DO n = 1, ntype
  ! loop over number of points per tile
  DO k = 1, surft_pts(n)
    l = surft_index(k, n)
    j = (land_index(l) - 1) / row_length + 1
    i = land_index(l) - (j-1) * row_length
    fvar(l, n) = jules_var(i, j)
  END DO
END DO
cable_var =  PACK(fvar, l_tile_pts)

END SUBROUTINE cable_pack_rr


SUBROUTINE cable_pack_lp(umvar, defaultin, cablevar, soiltype, skip )

USE cable_types_mod,          ONLY: mp, l_tile_pts
USE jules_surface_types_mod,  ONLY: ntype
USE cable_surface_types_mod,  ONLY: mstype, ICE_soiltype_cbl
                              ! number of non-veg (i.e. soil) parameters

USE ancil_info,               ONLY: land_pts, surft_pts, surft_index

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   UM met forcing vars needed by CABLE which have UM dimensions
!   (land_points)[_lp], which is no good to cable. These have to be
!   re-packed in a single vector of active tiles. Hence we use
!   conditional "mask" l_tile_pts(land_pts,ntiles) which is .true.
!   if the land point is/has an active tile
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!
!-----------------------------------------------------------------------------

REAL, INTENT(IN) :: umvar(land_pts)
REAL, INTENT(IN) :: defaultin(mstype)
REAL, INTENT(INOUT) :: cablevar(mp)
INTEGER, INTENT(INOUT) :: soiltype(mp)
REAL, ALLOCATABLE:: fvar(:,:)
LOGICAL, OPTIONAL :: skip
LOGICAL :: to_skip
INTEGER :: n,k,l

IF ( PRESENT(skip) ) THEN
  to_skip = skip
ELSE
  to_skip = .FALSE.
END IF


ALLOCATE( fvar(land_pts, ntype) )

! loop over ntype
fvar = 0.0
!skip_present = PRESENT(skip)
DO n = 1, ntype
  DO k = 1, surft_pts(n)
    ! index of each point per tile in an array of dim=(land_pts,ntype)
    l = surft_index(k,n)
     ! at this point fvar=umvar, ELSE=0.0
    fvar(l,n) = umvar(l)
    ! unless explicitly skipped by including argument in subroutine call
    IF ( to_skip ) THEN
       ! on perma frost tile, set fvar=defaultin
      IF ( n == ntype ) THEN
        fvar(l,n) =  defaultin(ICE_soiltype_cbl)
      END IF
    END IF
  END DO
END DO

! unless explicitly skipped by including argument in subroutine call
IF ( to_skip ) THEN
  cablevar = PACK(fvar, l_tile_pts)
  DO n = 1, mp
    IF (soiltype(n) == ICE_soiltype_cbl)                                     &
      cablevar(n) =  defaultin(ICE_soiltype_cbl)
  END DO
END IF

DEALLOCATE(fvar)

END SUBROUTINE cable_pack_lp

END MODULE cable_pack_mod
