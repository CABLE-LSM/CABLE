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

MODULE init_cable_parms_mod

!------------------------------------------------------------------------------
! Description:
!   Initialises CABLE parameter variables from values read in namelist
!
! This MODULE is USEd in:
!      cable_land_albedo_mod_cbl.F90
!
! This MODULE contains 1 public Subroutine:
!      init_cable_parms_rad
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!------------------------------------------------------------------------------

IMPLICIT NONE

PUBLIC init_cable_parms_rad

CONTAINS

SUBROUTINE init_cable_parms_rad( VegXfang, VegTaul, VegRefl, SurfaceType,      &
                                 SoilType, mp, nrb, l_tile_pts,                &
                                 ICE_SurfaceType, ICE_SoilType, VeginXfang,    &
                                 VeginTaul, VeginRefl, land_pts, nsurft,       &
                                 tile_frac )
! Description:
!   Nothing further to add to the module description.

IMPLICIT NONE

INTEGER, INTENT(OUT), ALLOCATABLE :: SurfaceType(:) !CABLE surface tile PFT/nveg
INTEGER, INTENT(OUT), ALLOCATABLE :: SoilType(:)    !CABLE soil type per tile
REAL,    INTENT(OUT), ALLOCATABLE :: VegXfang(:)    !Leaf Angle (CABLE)
REAL,    INTENT(OUT), ALLOCATABLE :: VegTaul(:,:)   !Leaf Transmisivity (CABLE)
REAL,    INTENT(OUT), ALLOCATABLE :: VegRefl(:,:)   !Leaf Reflectivity (CABLE)

INTEGER, INTENT(IN) :: mp
INTEGER, INTENT(IN) :: nrb
INTEGER, INTENT(IN) :: nsurft
INTEGER, INTENT(IN) :: land_pts
INTEGER, INTENT(IN) :: ICE_SurfaceType              !index ICE surface type
INTEGER, INTENT(IN) :: ICE_SoilType                 !index soil type
LOGICAL, INTENT(IN) :: l_tile_pts(land_pts,nsurft)  ! TRUE if active tile
REAL,    INTENT(IN) :: tile_frac(land_pts,nsurft)

!---IN: CABLE Vegetation parameters. decl in params_io_cbl.F90
REAL, INTENT(IN) :: VeginXfang(nsurft)              !Leaf Angle
REAL, INTENT(IN) :: VeginTaul(nrb, nsurft )         !Leaf Transmisivity
REAL, INTENT(IN) :: VeginRefl(nrb, nsurft )         !Leaf Reflectivity

!local vars
INTEGER :: JSurfaceTypeID(land_pts,nsurft)
INTEGER :: i, j, h

IF (.NOT. ALLOCATED(SurfaceType) ) ALLOCATE (SurfaceType(mp) )
IF (.NOT. ALLOCATED(VegXfang) )    ALLOCATE (VegXfang(mp) )
IF (.NOT. ALLOCATED(VegRefl) )     ALLOCATE (VegRefl(mp, nrb) )
IF (.NOT. ALLOCATED(VegTaul) )     ALLOCATE (VegTaul(mp, nrb) )
IF (.NOT. ALLOCATED(SoilType))     ALLOCATE( SoilType(mp) )

!local var to pack surface type:
JSurfaceTypeID = 0
DO j = 1, land_pts
  DO i = 1, nsurft
    IF ( tile_frac(j,i) > 0 ) JSurfaceTypeID(j,i) = i
  END DO
END DO

SurfaceType = PACK( JSurfaceTypeID, L_tile_pts)

SoilType(:)=2
DO j=1,mp
  IF (SurfaceType(j) == ICE_SurfaceType) THEN
    SoilType(j)= ICE_SoilType
  END IF
END DO

! Prescribe parameters for current gridcell based on veg/soil type
! (which may have loaded from default value file or met file):
DO h = 1, mp          ! over each patch in current grid
  VegTaul(h,1)   = VeginTaul(1,SurfaceType(h))
  VegTaul(h,2)   = VeginTaul(2,SurfaceType(h))
  VegRefl(h,1)   = VeginRefl(1,SurfaceType(h))
  VegRefl(h,2)   = VeginRefl(2,SurfaceType(h))
  VegXfang(h)    = VeginXfang( SurfaceType(h))
END DO ! over each veg patch in land point

RETURN
END SUBROUTINE init_cable_parms_rad

END MODULE init_cable_parms_mod
