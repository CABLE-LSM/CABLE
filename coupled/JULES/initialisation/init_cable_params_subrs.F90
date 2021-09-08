#if !defined(UM_JULES) 

MODULE init_cable_pftparms_mod

IMPLICIT NONE

CONTAINS
 
SUBROUTINE init_cable_veg()

USE cable_types_mod,   ONLY: mp, l_tile_pts
USE cable_params_mod,  ONLY: veg => veg_cbl, vegin
USE ancil_info,        ONLY: nsurft, land_pts, frac_surft

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   init_cables veg parameters using values read from namelist
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!-----------------------------------------------------------------------------

INTEGER :: JSurfaceTypeID(land_pts,nsurft)  
INTEGER :: i
INTEGER :: h

!local var to pack surface type:
JSurfaceTypeID = 0
DO i = 1,nsurft
  IF ( frac_surft(1,i) > 0 ) JSurfaceTypeID(:,i) = i
END DO

veg%iveg = PACK( JSurfaceTypeID, L_tile_pts)

! Prescribe parameters for current gridcell based on veg/soil type
! (which may have loaded from default value file or met file):
DO h = 1, mp          ! over each patch in current grid
  veg%taul(h,1)   = vegin%taul(1,veg%iveg(h))
  veg%taul(h,2)   = vegin%taul(2,veg%iveg(h))
  veg%refl(h,1)   = vegin%refl(1,veg%iveg(h))
  veg%refl(h,2)   = vegin%refl(2,veg%iveg(h))
  veg%cplant(h,1)   = vegin%cplant(1,veg%iveg(h))
  veg%cplant(h,2)   = vegin%cplant(2,veg%iveg(h))
  veg%cplant(h,3)   = vegin%cplant(3,veg%iveg(h))
  veg%csoil(h,1)   = vegin%csoil(1,veg%iveg(h))
  veg%csoil(h,2)   = vegin%csoil(2,veg%iveg(h))
  veg%ratecp(h,1)   = vegin%ratecp(1,veg%iveg(h))
  veg%ratecp(h,2)   = vegin%ratecp(2,veg%iveg(h))
  veg%ratecp(h,3)   = vegin%ratecp(3,veg%iveg(h))
  veg%ratecs(h,1)   = vegin%ratecs(1,veg%iveg(h))
  veg%ratecs(h,2)   = vegin%ratecs(2,veg%iveg(h))
  veg%hc(h)       = vegin%hc(veg%iveg(h))
  veg%xfang(h)    = vegin%xfang(veg%iveg(h))
  veg%frac4(h)    = vegin%frac4(veg%iveg(h))
  veg%canst1(h)   = vegin%canst1(veg%iveg(h))
  veg%dleaf(h)    = vegin%dleaf(veg%iveg(h))
  veg%vcmax(h)    = vegin%vcmax(veg%iveg(h))
  veg%ejmax(h)    = vegin%ejmax(veg%iveg(h))
  veg%vbeta(h)    = vegin%vbeta(veg%iveg(h))
  veg%xalbnir(h)  = vegin%xalbnir(veg%iveg(h))
  veg%rp20(h)     = vegin%rp20(veg%iveg(h))
  veg%rpcoef(h)   = vegin%rpcoef(veg%iveg(h))
  veg%rs20(h)     = vegin%rs20(veg%iveg(h))
  veg%shelrb(h)   = vegin%shelrb(veg%iveg(h))
  veg%wai(h)      = vegin%wai(veg%iveg(h))
  veg%a1gs(h)     = vegin%a1gs(veg%iveg(h))
  veg%d0gs(h)     = vegin%d0gs(veg%iveg(h))
  veg%vegcf(h)    = vegin%vegcf(veg%iveg(h))
  veg%extkn(h)    = vegin%extkn(veg%iveg(h))
  veg%tminvj(h)   = vegin%tminvj(veg%iveg(h))
  veg%tmaxvj(h)   = vegin%tmaxvj(veg%iveg(h))
  veg%g0(h)       = vegin%g0(veg%iveg(h)) ! Ticket #56
  veg%g1(h)       = vegin%g1(veg%iveg(h)) ! Ticket #56
  veg%a1gs(h)   = vegin%a1gs(veg%iveg(h))
  veg%d0gs(h)   = vegin%d0gs(veg%iveg(h))
  veg%alpha(h)  = vegin%alpha(veg%iveg(h))
  veg%convex(h) = vegin%convex(veg%iveg(h))
  veg%cfrd(h)   = vegin%cfrd(veg%iveg(h))
  veg%gswmin(h) = vegin%gswmin(veg%iveg(h))
  veg%conkc0(h) = vegin%conkc0(veg%iveg(h))
  veg%conko0(h) = vegin%conko0(veg%iveg(h))
  veg%ekc(h)    = vegin%ekc(veg%iveg(h))
  veg%eko(h)    = vegin%eko(veg%iveg(h))
  veg%rootbeta(h)  = vegin%rootbeta(veg%iveg(h))
  veg%zr(h)       = vegin%zr(veg%iveg(h))
  veg%clitt(h)    = vegin%clitt(veg%iveg(h))
END DO ! over each veg patch in land point

END SUBROUTINE init_cable_veg

END MODULE init_cable_pftparms_mod

#endif
#if !defined(UM_JULES) 

MODULE init_cable_soilparms_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE init_cable_soil()

!H!USE cable_types_mod,    ONLY: mp, l_tile_pts
USE cable_params_mod,  ONLY: soil => soil_cbl, soilin, veg => veg_cbl
!H!USE ancil_info,         ONLY : nsurft, land_pts, frac_surft

IMPLICIT NONE

soil%isoilm  =  2
   
WHERE ( veg%iveg == 17 ) soil%isoilm = 9
           
END SUBROUTINE init_cable_soil

END MODULE init_cable_soilparms_mod

#endif
