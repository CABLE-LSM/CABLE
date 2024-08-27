MODULE cable_um_init_respiration_mod
   
IMPLICIT NONE
PUBLIC init_respiration

CONTAINS
 
SUBROUTINE init_respiration( land_pts, ntiles, npft, l_tilepts,                &
                             npp_pft_acc, resp_w_pft_acc, canopy )
USE cable_def_types_mod, ONLY: canopy_type

INTEGER, INTENT(IN) :: land_pts          ! # land points being processed
INTEGER, INTENT(IN) :: ntiles            ! # tiles 
INTEGER, INTENT(IN) :: npft              ! # plant functional types
LOGICAL, INTENT(IN) :: L_tilepts(land_pts, ntiles)  ! TRUE if active tile
TYPE(canopy_type), INTENT(OUT) :: canopy
REAL, INTENT(IN)               :: npp_pft_acc(land_pts, npft)
REAL, INTENT(IN)               :: resp_w_pft_acc(land_pts, npft)

REAL  :: fnpp_pft_acc(land_pts, ntiles)
REAL  :: fresp_w_pft_acc(land_pts, ntiles)

! make ntile versions of npp_pft_acc/resp_w_pft_acc
fnpp_pft_acc(:,1:ntiles)    = 0.0
fresp_w_pft_acc(:,1:ntiles) = 0.0
fnpp_pft_acc(:,1:npft)      = npp_pft_acc(:,1:npft)
fresp_w_pft_acc(:,1:npft)   = resp_w_pft_acc(:,1:npft)  

!---set soil & plant respiration (now in dim(land_pts,ntiles))
canopy%frs = PACK(fnpp_pft_acc   , l_tilepts)
canopy%frp = PACK(fresp_w_pft_acc, l_tilepts)

!---convert units to g C m-2 s-1
canopy%frs = canopy%frs * 1000.
canopy%frp = canopy%frp * 1000.

RETURN
END SUBROUTINE init_respiration

END MODULE cable_um_init_respiration_mod




