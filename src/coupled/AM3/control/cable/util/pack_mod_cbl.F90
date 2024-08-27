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

MODULE cable_pack_mod

!-----------------------------------------------------------------------------
! Description:
!   JULES met forcing vars needed by CABLE commonly have JULES dimensions
!   (row_length,rows), which are no good for CABLE. These have to be
!   re-packed in a single vector of active tiles. Hence we use
!   conditional "mask" l_tile_pts(land_pts,ntiles) which is .true.
!   if the land point is/has an active tile. A packing routine for
!   land_point dimensioned JULES fields is to be included as required on the
!   subsequent (e.g. explicit) CALLs to CABLE
!
! This MODULE is USEd by:
!      cable_land_albedo_mod.F90
!
! This MODULE contains 1 public Subroutine:
!      cable_pack_rr,
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

IMPLICIT NONE
PUBLIC :: cable_pack_rr, pack_landpts2mp_ICE, pack_landpts2mp
PRIVATE

CONTAINS

SUBROUTINE cable_pack_rr( cable_var, jules_var, mp, l_tile_pts, row_length,    &
                          rows, ntype, land_pts, land_index, surft_pts,        &
                          surft_index )

! Description:
!   Nothing further to add to module description.

IMPLICIT NONE

INTEGER, INTENT(IN) :: mp, row_length, rows, ntype, land_pts
REAL, INTENT(OUT)   :: cable_var(mp)
LOGICAL, INTENT(IN) :: l_tile_pts(land_pts, ntype)
INTEGER, INTENT(IN) :: land_index(land_pts)           !Index in (x,y) array
INTEGER, INTENT(IN) :: surft_pts(ntype)              !# land points per PFT
INTEGER, INTENT(IN) :: surft_index(land_pts,ntype)   !Index in land_pts array
REAL, INTENT(IN)    :: jules_var(row_length, rows)
!local vars
REAL    :: fvar(land_pts, ntype)
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

RETURN
END SUBROUTINE cable_pack_rr

!--- UM met forcing vars needed by CABLE which have UM dimensions
!---(land_points)[_lp], which is no good to cable. These have to be 
!--- re-packed in a single vector of active tiles. Hence we use 
!--- conditional "mask" l_tile_pts(land_pts,ntiles) which is .true.
!--- if the land point is/has an active tile
SUBROUTINE pack_landpts2mp_ICE( nsurft, land_pts, mp, nsoil_max, ICE_soiltype, &
                                tile_pts, tile_index, L_tile_pts, umvar,       &
                                soiltype, ICE_value, cablevar )
  
IMPLICIT NONE
 
INTEGER, INTENT(IN) :: nsurft 
INTEGER, INTENT(IN) :: land_pts
INTEGER, INTENT(IN) :: mp 
INTEGER, INTENT(IN) :: nsoil_max 
INTEGER, INTENT(IN) :: ICE_soiltype
INTEGER, INTENT(IN) :: soiltype(mp)
INTEGER, INTENT(IN) :: tile_pts(nsurft)   
INTEGER, INTENT(IN) :: tile_index(land_pts, nsurft)   
LOGICAL, INTENT(IN) :: L_tile_pts(land_pts, nsurft)
REAL, INTENT(IN)    :: umvar(land_pts)
REAL, INTENT(IN)    :: ICE_value(nsoil_max)
REAL, INTENT(INOUT) :: cablevar(mp)
!local vars
REAL :: fvar(land_pts, nsurft)   
INTEGER :: n,k,l,i

fvar(:,:) = 0.0
DO n=1, nsurft

  ! loop over number of points per nth tile
  DO k=1,tile_pts(n)
    ! land pt index of point 
    L = tile_index(k,n)
    ! at this point fvar=umvar, ELSE=0.0 
    fvar(l,n) = umvar(l)
  ENDDO

ENDDO

cablevar = PACK(fvar,L_tile_pts)

DO i=1,mp
  
  IF(soiltype(i)==ICE_soiltype) THEN
    cablevar(i) =  ICE_value(ICE_soiltype)
  END IF
    
ENDDO        

END SUBROUTINE pack_landpts2mp_ICE

SUBROUTINE pack_landpts2mp( nsurft, land_pts, mp, tile_pts, tile_index,        &
                            L_tile_pts, umvar, cablevar )
  
IMPLICIT NONE
 
INTEGER, INTENT(IN) :: nsurft 
INTEGER, INTENT(IN) :: land_pts
INTEGER, INTENT(IN) :: mp 
LOGICAL, INTENT(IN) :: L_tile_pts(land_pts, nsurft)   
INTEGER, INTENT(IN) :: tile_pts(nsurft)   
INTEGER, INTENT(IN) :: tile_index(land_pts, nsurft)   
REAL, INTENT(IN)    :: umvar(land_pts) 
REAL, INTENT(INOUT) :: cablevar(mp)
!local vars
REAL :: fvar(land_pts, nsurft)   
INTEGER :: n,k,l
 
fvar(:,:) = 0.0
DO n=1, nsurft
  ! loop over number of points per nth tile
  DO k=1,tile_pts(n)
    ! land pt index of point 
    l = tile_index(k,n)
    ! at this point fvar=umvar, ELSE=0.0 
    fvar(l,n) = umvar(l)
  ENDDO
ENDDO

cablevar = PACK(fvar,L_tile_pts)
 
END SUBROUTINE pack_landpts2mp


END MODULE cable_pack_mod
