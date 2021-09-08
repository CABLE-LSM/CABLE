MODULE cable_rad_unpack_mod
  
CONTAINS
 
SUBROUTINE cable_rad_unpack(                                                  &
land_albedo,                                                                  &
alb_surft,                                                                    &
mp,                                                                           &
nrb,                                                                          &
row_length,                                                                   &
rows,                                                                         &
land_pts,                                                                     &
nsurft,                                                                       &
sm_levels,                                                                    &
tile_pts,                                                                     &
tile_index,                                                                   &
land_index,                                                                   &
tile_frac,                                                                    &
L_tile_pts,                                                                   &
EffSurfRefl_dif,                                                              &
EffSurfRefl_beam,                                                             &
snow_flag_cable,                                                              &
SnowFlag_3L                                                                   &


                           )

IMPLICIT NONE

!___ re-decl input args

!model dimensions
!-------------------------------------------------------------------------------
!JaC:todo:ultimatelty get this from JaC~
INTEGER :: mp                       !total number of "tiles"  
INTEGER :: nrb                      !number of radiation bands [per legacy=3, but really=2 VIS,NIR. 3rd dim was for LW]
INTEGER :: row_length                       !grid cell x
INTEGER :: rows                             !grid cell y
INTEGER :: land_pts                         !grid cell land points on the x,y grid
INTEGER :: nsurft                           !grid cell number of surface types 
INTEGER :: sm_levels                        !grid cell number of soil levels 
INTEGER :: tile_pts(nsurft)                 !Number of land points per PFT 
INTEGER :: tile_index(land_pts,nsurft)      !Index of land point in (land_pts) array
INTEGER :: land_index(land_pts)             !Index of land points in (x,y) array - see below
!-------------------------------------------------------------------------------

!Return from CABLE to vulvill contract with JULES
!-------------------------------------------------------------------------------
REAL :: land_albedo(row_length,rows,4)      
REAL :: alb_surft(Land_pts,nsurft,4)        
!-------------------------------------------------------------------------------

!recieved as spatial maps from the UM. 
!-------------------------------------------------------------------------------
REAL :: tile_frac(land_pts,nsurft)          !fraction of each surface type per land point 
LOGICAL :: L_tile_pts( land_pts, nsurft )   !mask:=TRUE where tile_frac>0, else FALSE. pack mp according to this mask
!-------------------------------------------------------------------------------

!recieved as spatial maps from the UM. remapped to "mp"
!-------------------------------------------------------------------------------
INTEGER:: surface_type(mp)          ! Integer index of Surface type (veg%iveg)
REAL :: LAI_pft_cbl(mp)             !LAI -  "limited" and remapped
REAL :: HGT_pft_cbl(mp)             !canopy height -  "limited" and remapped
!-------------------------------------------------------------------------------

!Prognostics
!-------------------------------------------------------------------------------
INTEGER:: SnowFlag_3L(mp)           !Flag to treat snow as 3 layer  - if enough present. Updated depending on total depth (ssnow%isflag)
REAL :: snow_flag_cable(land_pts,nsurft) !Flag to treat snow as 3 layer  - REALized and in JULES dimensioonal format
!-------------------------------------------------------------------------------
                                                                                           
! Albedos
!-------------------------------------------------------------------------------
REAL :: EffSurfRefl_dif(mp,nrb)     !Effective Surface Relectance as seen by atmosphere [Diffuse SW]  (rad%reffdf)
REAL :: EffSurfRefl_beam(mp,nrb)    !Effective Surface Relectance as seen by atmosphere [Direct Beam SW] (rad%reffbm)
!-------------------------------------------------------------------------------


!___ local vars
INTEGER :: i,j,k,l,n
REAL :: miss = 0.0
 
! std template args 
CHARACTER(LEN=*), PARAMETER :: subr_name = "cable_rad_unpack"
REAL :: Sumreffbm(mp) 
REAL :: Sumreffdf(mp) 
  
! only for land points, at present do not have a method for treating 
! mixed land/sea or land/seaice points as yet.
!   Albedo for surface tiles
!     (:,:,1) direct beam visible
!     (:,:,2) diffuse visible
!     (:,:,3) direct beam near-IR
!     (:,:,4) diffuse near-IR
alb_surft(:,:,:) = 0.0
alb_surft(:,:,1) = UNPACK(EffSurfRefl_beam(:,1),l_tile_pts, miss)
alb_surft(:,:,3) = UNPACK(EffSurfRefl_beam(:,2),l_tile_pts, miss)
alb_surft(:,:,2) = UNPACK(EffSurfRefl_dif(:,1),l_tile_pts, miss)
alb_surft(:,:,4) = UNPACK(EffSurfRefl_dif(:,2),l_tile_pts, miss)

DO i = 1,land_pts
  DO j = 1,nsurft
              
    IF ( alb_surft(i,j,1)> 1.0) THEN
      PRINT  * , 'alb > 1',alb_surft(i,j,1)
      STOP
    END IF
 
  END DO
END DO

land_albedo = 0

DO n = 1,nsurft
  DO k = 1,tile_pts(n)
    l = tile_index(k,n)
    j=(land_index(l) - 1) / row_length + 1
    i = land_index(l) - (j-1) * row_length
    
    ! direct beam visible
    land_albedo(i,j,1) = land_albedo(i,j,1) +                                 &
                               tile_frac(l,n) * ALB_surft(l,n,1)

    ! diffuse beam visible
    land_albedo(i,j,2) = land_albedo(i,j,2) +                                 &
                               tile_frac(l,n) * ALB_surft(l,n,2)

    ! direct beam nearinfrared 
    land_albedo(i,j,3) = land_albedo(i,j,3) +                                 &
                               tile_frac(l,n) * ALB_surft(l,n,3)

    ! diffuse beam nearinfrared
    land_albedo(i,j,4) = land_albedo(i,j,4) +                                 &
                               tile_frac(l,n) * ALB_surft(l,n,4)
  END DO
END DO

RETURN
     
END SUBROUTINE cable_rad_unpack
END MODULE cable_rad_unpack_mod
