!******************************COPYRIGHT********************************************
! (c) CSIRO 2022.
! All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms and
! conditions set out therein.
!
! [Met Office Ref SC0237]
!******************************COPYRIGHT********************************************
MODULE cable_rad_unpack_mod

!-----------------------------------------------------------------------------
! Description:
!  Unpack JULES variables into CABLE variables as needed for the
!  radiation calculations.
!
! This MODULE is USEd in:
!     cable_land_albedo_mod_cbl.F90 (JULES)
!
! This MODULE contains 1 public Subroutine:
!     cable_rad_unpack
!
! Code owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

IMPLICIT NONE
PUBLIC :: cable_rad_unpack

CONTAINS

SUBROUTINE cable_rad_unpack( land_albedo, alb_surft,                           &
                             mp, nrs, row_length, rows, land_pts,              &
                             nsurft, tile_pts, tile_index,                     &
                             land_index, tile_frac, L_tile_pts,                &
                             EffSurfRefl_dif, EffSurfRefl_beam )

! Description:
!   Nothing further to add to module description.

IMPLICIT NONE

! Model(field) dimensions
INTEGER, INTENT(IN) :: mp                       !total number of "tiles"
INTEGER, INTENT(IN) :: nrs                      !# rad bands VIS,NIR. 3rd WAS LW
INTEGER, INTENT(IN) :: row_length               !grid cell x
INTEGER, INTENT(IN) :: rows                     !grid cell y
INTEGER, INTENT(IN) :: land_pts                 !grid cell land points -x,y grid
INTEGER, INTENT(IN) :: nsurft                   !grid cell # surface types

! Return Albedos CABLE to fulfill contract with JULES
REAL, INTENT(OUT) :: land_albedo(row_length,rows,nrs)
REAL, INTENT(OUT) :: alb_surft(land_pts,nsurft,nrs)

! Inherited dimensions from JULES
INTEGER, INTENT(IN) :: tile_pts(nsurft)         !Number of land points per PFT
INTEGER, INTENT(IN) :: land_index(land_pts)     !land point Index in (x,y) array
INTEGER, INTENT(IN) :: tile_index(land_pts,nsurft) !Index of land point in (land_pts) array

!recieved as spatial maps from the UM.
REAL,    INTENT(IN) :: tile_frac(land_pts,nsurft)     !fraction of each surface type per land point
LOGICAL, INTENT(IN) :: L_tile_pts( land_pts, nsurft ) !mask:=TRUE where tile_frac>0, else FALSE. pack mp according to this mask

! Albedos
REAL, INTENT(IN) :: EffSurfRefl_dif(mp,nrs)     !Effective Surface Relectance as seen by atmosphere [Diffuse SW]  (rad%reffdf)
REAL, INTENT(IN) :: EffSurfRefl_beam(mp,nrs)    !Effective Surface Relectance as seen by atmosphere [Direct Beam SW] (rad%reffbm)

!___ local vars
INTEGER :: i,j,k,l,n
REAL :: miss = 0.0

! std template args
CHARACTER(LEN=*), PARAMETER :: subr_name = "cable_rad_unpack"
REAL :: Sumreffbm(mp)
REAL :: Sumreffdf(mp)

! UNPACK Albedo (per rad stream) per surface tile
alb_surft(:,:,:) = 0.0        ! guarantee flushed
! Direct beam, visible / near-IR
alb_surft(:,:,1) = UNPACK(EffSurfRefl_beam(:,1),l_tile_pts, miss)
alb_surft(:,:,3) = UNPACK(EffSurfRefl_beam(:,2),l_tile_pts, miss)
! Diffuse, visible / near-IR
alb_surft(:,:,2) = UNPACK(EffSurfRefl_dif(:,1),l_tile_pts, miss)
alb_surft(:,:,4) = UNPACK(EffSurfRefl_dif(:,2),l_tile_pts, miss)

! ERROR trap: Model stopped as albedo is unphysical
DO i = 1,land_pts
  DO j = 1,nsurft

    IF ( alb_surft(i,j,1) > 1.0 .OR. alb_surft(i,j,1) < 0.0) THEN
      WRITE(6,*) 'albedo(i,j,1) is unphysical ',alb_surft(i,j,1)
      STOP 'CABLE ERROR'
    ELSE IF ( alb_surft(i,j,2) > 1.0 .OR. alb_surft(i,j,2) < 0.0) THEN
      WRITE(6,*) 'albedo(i,j,2) is unphysical ',alb_surft(i,j,2)
      STOP 'CABLE ERROR'
    ELSE IF ( alb_surft(i,j,3) > 1.0 .OR. alb_surft(i,j,3) < 0.0) THEN
      WRITE(6,*) 'albedo(i,j,3) is unphysical ',alb_surft(i,j,3)
      STOP 'CABLE ERROR'
    ELSE IF ( alb_surft(i,j,4) > 1.0 .OR. alb_surft(i,j,4) < 0.0) THEN
      WRITE(6,*) 'albedo(i,j,4) is unphysical ',alb_surft(i,j,4)
      STOP 'CABLE ERROR'
    END IF

  END DO
END DO

! Aggregate albedo (per rad stream) OVER surface tiles to get per cell value
land_albedo(:,:,:) = 0.0        ! guarantee flushed
DO n = 1,nsurft
  DO k = 1,tile_pts(n)

    l = tile_index(k,n)
    j=(land_index(l) - 1) / row_length + 1
    i = land_index(l) - (j-1) * row_length

    ! Direct beam, visible
    land_albedo(i,j,1) = land_albedo(i,j,1) + tile_frac(l,n) * ALB_surft(l,n,1)
    ! Diffuse, visible
    land_albedo(i,j,2) = land_albedo(i,j,2) + tile_frac(l,n) * ALB_surft(l,n,2)
    ! Direct beam, nearinfrared
    land_albedo(i,j,3) = land_albedo(i,j,3) + tile_frac(l,n) * ALB_surft(l,n,3)
    ! Diffuse, nearinfrared
    land_albedo(i,j,4) = land_albedo(i,j,4) + tile_frac(l,n) * ALB_surft(l,n,4)

  END DO
END DO

RETURN

END SUBROUTINE cable_rad_unpack

END MODULE cable_rad_unpack_mod
