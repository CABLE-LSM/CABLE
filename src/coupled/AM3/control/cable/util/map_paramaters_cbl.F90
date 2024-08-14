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

MODULE map_cable_parms_mod

!------------------------------------------------------------------------------
! Description:
!   Initialises CABLE parameter variables from values read in namelist
!
! This MODULE is USEd in:
!      cable_land_albedo_mod_cbl.F90
!
! This MODULE contains 1 public Subroutine:
!      map_cable_parms_rad
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!------------------------------------------------------------------------------

IMPLICIT NONE

PUBLIC map_cable_parms

CONTAINS

SUBROUTINE map_cable_parms( mp, ms, nrb, land_pts, nsurft, l_tile_pts,         &
                             ICE_SurfaceType, ICE_SoilType, soil_zse, veg,     &
                             soil, pars, tile_frac )
                                 
! Description:
!   Nothing further to add to the module description.

USE params_io_mod_cbl,   ONLY: params_io_data_type
USE cable_def_types_mod, ONLY: veg_parameter_type
USE cable_def_types_mod, ONLY: soil_parameter_type

IMPLICIT NONE

INTEGER, INTENT(IN) :: mp
INTEGER, INTENT(IN) :: ms
INTEGER, INTENT(IN) :: nrb
INTEGER, INTENT(IN) :: nsurft
INTEGER, INTENT(IN) :: land_pts
INTEGER, INTENT(IN) :: ICE_SurfaceType              !index ICE surface type
INTEGER, INTENT(IN) :: ICE_SoilType                 !index soil type
LOGICAL, INTENT(IN) :: l_tile_pts(land_pts,nsurft)  ! TRUE if active tile
REAL,    INTENT(IN) :: soil_zse(ms)                 ! soil depth per layer 
REAL,    INTENT(IN) :: tile_frac(land_pts,nsurft)

TYPE(veg_parameter_type),  INTENT(OUT) :: veg        ! vegetation parameters
TYPE(soil_parameter_type), INTENT(OUT) :: soil       ! soil parameters
TYPE(params_io_data_type), INTENT(IN)  :: pars


!local vars
INTEGER :: JSurfaceTypeID(land_pts,nsurft)
INTEGER :: i, j, h

CALL DefinePatchType( veg%iveg, soil%isoilm, mp, land_pts, nsurft,             &
                      ICE_SurfaceType, ICE_SoilType, tile_frac, L_tile_pts )

CALL init_ALLveg( veg%iveg, mp, ms, nrb, veg, pars )

CALL init_ALLsoil( soil%isoilm, mp, ms, soil_zse, soil, pars )

RETURN
END SUBROUTINE map_cable_parms

SUBROUTINE DefinePatchType( SurfaceType, SoilType, mp, land_pts, nsurft,       &
                           ICE_SurfaceType, ICE_SoilType, tile_frac, L_tile_pts)
                            
IMPLICIT NONE

INTEGER, INTENT(IN) :: mp
INTEGER, INTENT(OUT) :: SurfaceType(mp) !CABLE surface tile PFT/nveg
INTEGER, INTENT(OUT) :: SoilType(mp)    !CABLE soil type per tile

INTEGER, INTENT(IN) :: nsurft
INTEGER, INTENT(IN) :: land_pts
INTEGER, INTENT(IN) :: ICE_SurfaceType              ! index ICE surface type
INTEGER, INTENT(IN) :: ICE_SoilType                 ! index soil type
LOGICAL, INTENT(IN) :: l_tile_pts(land_pts,nsurft)  ! TRUE if active tile
REAL,    INTENT(IN) :: tile_frac(land_pts,nsurft)

!local vars
INTEGER :: JSurfaceTypeID(land_pts,nsurft)
INTEGER :: i, j

!local var to pack surface type:
JSurfaceTypeID = 0
DO j = 1, land_pts
  DO i = 1, nsurft
    IF ( tile_frac(j,i) > 0 ) JSurfaceTypeID(j,i) = i
  END DO
END DO

SurfaceType = 0
SurfaceType = PACK( JSurfaceTypeID, L_tile_pts)

SoilType(:)=2
DO j=1,mp
  IF (SurfaceType(j) == ICE_SurfaceType) THEN
    SoilType(j)= ICE_SoilType
  END IF
END DO

END SUBROUTINE DefinePatchType

SUBROUTINE init_ALLveg( SurfaceType, mp, ms, nrb, veg, pars )
                                 
! Description: Map per PFT parameters onto CABLE vectors of length mp
!              Note: Canopy height (%hc) & LAI (%vlai) are clobbered by
!              range limited UM/JULES spatial fields in limit_HGT_LAI()

USE params_io_mod_cbl,   ONLY: params_io_data_type
USE cable_def_types_mod, ONLY: veg_parameter_type, r_2

IMPLICIT NONE

INTEGER, INTENT(IN) :: mp
INTEGER, INTENT(IN) :: SurfaceType(mp) ! CABLE surface tile PFT/nveg
INTEGER, INTENT(IN) :: ms
INTEGER, INTENT(IN) :: nrb

TYPE(veg_parameter_type),  INTENT(OUT) :: veg        ! vegetation parameters
TYPE(params_io_data_type), INTENT(IN)  :: pars

!local vars
INTEGER :: i, h

DO h = 1, mp          ! over each patch in current grid

  veg%frac4(h)    = pars%vegin_frac4(SurfaceType(h))
  veg%taul(h,1)   = pars%vegin_taul(1,SurfaceType(h))
  veg%taul(h,2)   = pars%vegin_taul(2,SurfaceType(h))
  veg%refl(h,1)   = pars%vegin_refl(1,SurfaceType(h))
  veg%refl(h,2)   = pars%vegin_refl(2,SurfaceType(h))
  veg%canst1(h)   = pars%vegin_canst1(SurfaceType(h))
  veg%dleaf(h)    = pars%vegin_dleaf(SurfaceType(h))
  veg%vcmax(h)    = pars%vegin_vcmax(SurfaceType(h))
  veg%ejmax(h)    = pars%vegin_ejmax(SurfaceType(h))
  veg%hc(h)       = pars%vegin_hc(SurfaceType(h))
  veg%xfang(h)    = pars%vegin_xfang(SurfaceType(h))
  veg%vbeta(h)    = pars%vegin_vbeta(SurfaceType(h))
  veg%xalbnir(h)  = pars%vegin_xalbnir(SurfaceType(h))
  veg%rp20(h)     = pars%vegin_rp20(SurfaceType(h))
  veg%rpcoef(h)   = pars%vegin_rpcoef(SurfaceType(h))
  veg%rs20(h)     = pars%vegin_rs20(SurfaceType(h))
  veg%shelrb(h)   = pars%vegin_shelrb(SurfaceType(h))
  veg%wai(h)      = pars%vegin_wai(SurfaceType(h))
  veg%a1gs(h)     = pars%vegin_a1gs(SurfaceType(h))
  veg%d0gs(h)     = pars%vegin_d0gs(SurfaceType(h))
  veg%vegcf(h)    = pars%vegin_vegcf(SurfaceType(h))
  veg%extkn(h)    = pars%vegin_extkn(SurfaceType(h))
  veg%tminvj(h)   = pars%vegin_tminvj(SurfaceType(h))
  veg%tmaxvj(h)   = pars%vegin_tmaxvj(SurfaceType(h))
  veg%g0(h)       = pars%vegin_g0(SurfaceType(h)) ! Ticket #56
  veg%g1(h)       = pars%vegin_g1(SurfaceType(h)) ! Ticket #56
  veg%a1gs(h)     = pars%vegin_a1gs(SurfaceType(h))
  veg%d0gs(h)     = pars%vegin_d0gs(SurfaceType(h))
  veg%alpha(h)    = pars%vegin_alpha(SurfaceType(h))
  veg%convex(h)   = pars%vegin_convex(SurfaceType(h))
  veg%cfrd(h)     = pars%vegin_cfrd(SurfaceType(h))
  veg%gswmin(h)   = pars%vegin_gswmin(SurfaceType(h))
  veg%conkc0(h)   = pars%vegin_conkc0(SurfaceType(h))
  veg%conko0(h)   = pars%vegin_conko0(SurfaceType(h))
  veg%ekc(h)      = pars%vegin_ekc(SurfaceType(h))
  veg%eko(h)      = pars%vegin_eko(SurfaceType(h))
  veg%rootbeta(h) = pars%vegin_rootbeta(SurfaceType(h))
  veg%zr(h)       = pars%vegin_zr(SurfaceType(h))
  veg%clitt(h)    = pars%vegin_clitt(SurfaceType(h))

END DO ! over each veg patch in land point
 
RETURN
END SUBROUTINE init_ALLveg

SUBROUTINE init_ALLsoil( SoilType, mp, ms, soil_zse, soil, pars )

! Description: UM/JULES spatial fields are used later in initialize_soil() to
!              define soil properties. However note these are per landpoint as 
!              JULES doesn't have soil per tile. We use PFT=ICE to define 
!              Soil=ICE and use values from here. Also, there is no soil density
!              field in the UM, we effectively use density @ soil%isoilm=2 

USE params_io_mod_cbl,   ONLY: params_io_data_type
USE cable_def_types_mod, ONLY: soil_parameter_type

IMPLICIT NONE

INTEGER, INTENT(IN) :: mp
INTEGER, INTENT(IN) :: ms
INTEGER, INTENT(IN) :: SoilType(mp)    !CABLE soil type per tile
REAL,    INTENT(IN) :: soil_zse(ms)                 ! soil depth per layer 

TYPE(soil_parameter_type), INTENT(OUT) :: soil       ! vegetation parameters
TYPE(params_io_data_type), INTENT(IN)  :: pars

!local vars
INTEGER :: h

soil%zse(:) = soil_zse(:)

DO h = 1, mp          ! over each patch in current grid

  soil%swilt(h)   =  pars%soilin_swilt(SoilType(h))
  soil%sfc(h)     =  pars%soilin_sfc(SoilType(h))
  soil%ssat(h)    =  pars%soilin_ssat(SoilType(h))
  soil%bch(h)     =  pars%soilin_bch(SoilType(h))
  soil%hyds(h)    =  pars%soilin_hyds(SoilType(h))
  soil%sucs(h)    =  pars%soilin_sucs(SoilType(h))
  soil%rhosoil(h) =  pars%soilin_rhosoil(SoilType(h))
  soil%css(h)     =  pars%soilin_css(SoilType(h))
  soil%silt(h)    =  pars%soilin_silt(SoilType(h))
  soil%clay(h)    =  pars%soilin_clay(SoilType(h))
  soil%sand(h)    =  pars%soilin_sand(SoilType(h))
 
END DO

END SUBROUTINE init_ALLsoil

END MODULE map_cable_parms_mod
