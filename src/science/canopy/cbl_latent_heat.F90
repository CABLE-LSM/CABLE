MODULE cable_latent_heat_module
  !* This MODULE contains one SUBROUTINE which evaluates the latent heat from
  ! ground given the previoulsy evaluated rate of potential evaporation.
  ! This is a component of the calcualtion of the energy balance and is carried out
  ! once per cycle in the MO iteration in [[define_canopy]].

  IMPLICIT NONE

  PUBLIC latent_heat_flux 
  PRIVATE

CONTAINS

SUBROUTINE Latent_heat_flux( mp, CTFRZ, dels, soil_zse, soil_swilt,           &
                             cable_user_l_new_reduce_soilevp, pwet, air_rlam,  &
                             ssnow_snowd, ssnow_wb, ssnow_wbice,             &
                             ssnow_pudsto, ssnow_pudsmx, ssnow_potev,          &
                             ssnow_wetfac, ssnow_evapfbl, ssnow_cls,          & 
                             ssnow_tss, canopy_fes, canopy_fess, canopy_fesp  )
  !*## Purpose
  !
  ! This SUBROUTINE converts the previously evaluated rate of potential
  ! evaporation `ssnow_potev` from the soil to the associated latent heat flux
  ! `canopy_fes`. The evaluation accounts for whether the water flux is
  ! from/to a wet or ice surface,
  ! the direction of the flux (to or from the soil), any limits due to a lack of
  ! water the surface layer of the soil, and
  ! the partition of the latent heat flux between soil `canopy_fess` and
  ! puddle `canopy_fesp` contributions.
  !
  !## Warning
  !
  ! This current implementation is possibly confusing.  The variable `ssnow_potev`
  ! is not potential evaporation (in kg water m\(^{-2}\)s\(^{-1}\)) but is instead
  ! the latent heat flux (in Wm\(^{-2}\)) associated with the potential evaporation.
  ! `ssnow_potev` assumes that the water flux described by the potential
  ! evaporation is of liquid water regardless of temperature
  ! and/or whether the surface is covered by snow. Additionally one option for the
  ! evaluation of `ssnow_potev` (the [[Penman_Monteith]] option) does not function
  ! correctly over frozen surfaces. A bug fix will require changes to
  ! [[define_canopy]], [[Penman_Monteith]] and this SUBROUTINE; markers have been
  ! inserted into the code accordingly.  

USE cable_def_types_mod, ONLY : r_2
IMPLICIT NONE

INTEGER :: mp
REAL(KIND=r_2), INTENT(OUT) :: canopy_fes(mp) !! latent heat flux from the ground (Wm\(^{-2}\)) 
REAL(KIND=r_2), INTENT(OUT) :: canopy_fess(mp)!! latent heat flux from the soil (Wm\(^{-2}\))
REAL(KIND=r_2), INTENT(OUT) :: canopy_fesp(mp)!! latent heat flux from any puddles (Wm\(^{-2}\))
REAL, INTENT(OUT) :: ssnow_cls(mp)
!! factor denoting whether latent flux is to/from water or ice (=1 if water, =1.1335 if ice)
REAL, INTENT(IN OUT) :: pwet(mp)
!! factor to reduce soil evaporation due to presence of a puddle (-)
REAL, INTENT(IN OUT) :: ssnow_wetfac(mp)      !! wetness factor for soil (between 0 and 1)


REAL, INTENT(IN) :: CTFRZ                     !! temperature at freezing point (K) 
REAL, INTENT(IN) :: dels                      !! time step of CABLE (s)
REAL, INTENT(IN) :: soil_zse                  !! thickness of topmost soil layer (m)
REAL, INTENT(IN) :: soil_swilt(mp)
!! soil moisture content at wilting point (m\(^3\) water m\(^{-3}\) volume of soil)
LOGICAL , INTENT(IN) :: cable_user_l_new_reduce_soilevp
!! NAMELIST switch to use alternate soil evaporation scheme

REAL, INTENT(IN) :: air_rlam(mp)              !! density of air (kgm\(^{-3}\))
REAL, INTENT(IN) :: ssnow_snowd(mp)
!! depth of snow in liquid water equivalent (mm m\(^{-2}\))
REAL, INTENT(IN) :: ssnow_pudsto(mp)          !! amount of water in puddles (kgm\(^{-2}\))
REAL, INTENT(IN) :: ssnow_pudsmx(mp)
!! maximum amount of water possible in puddles (kgm\(^{-2}\))
REAL, INTENT(IN) :: ssnow_potev(mp)
!! latent heat flux associated potential evaporation (Wm\(^{-2}\))
REAL, INTENT(IN) :: ssnow_evapfbl(mp)         !! flux of water from soil surface (kg m\(^{-2}\))
REAL(KIND=r_2), INTENT(IN) :: ssnow_wb(mp)
!! water content in surface soil layer (m\(^{3}\) liquid water m\(^{-3}\) volume of soil)
REAL(KIND=r_2), INTENT(IN) :: ssnow_wbice(mp)
!! ice content in surface soil layer (m\(^{3}\) frozen water m\(^{-3}\) volume of soil)
REAL, INTENT(IN) :: ssnow_tss(mp)             !! temperature of surface soil/snow layer (K)

REAL, DIMENSION(mp) ::                                                      &
  frescale,  flower_limit, fupper_limit

INTEGER :: j

!>## Method
!
! 'Potential evaporation' quantifies the theoretical value of evaporation from
! a fully saturated surface (with no limits on water availability) at a given
! temperature and humidity, in response to the available energy.
! This SUBROUTINE evaluates the value of the latent heat flux from the
! ground/soil/snow packgiven the potential evaporation according to one of four cases
! (see Ticket 137)
!
!  1. evaporation from/dew onto surfaces with no snow when soil temperature > freezing point
!  2. sublimation from/deposition of frost onto surfaces with snow cover
!  3. evaporation of liquid water from a frozen soil column if there is no snow present
!     (sublimation of frozen water from within the soil column is not possible)
!  4. deposition of frost onto frozen soils if there is no snow cover
!
! while satisfying any constraints due to limits on soil moisture.
!
!<br><\br>
! There are six outputs from this SUBROUTINE
!
! * three latent heat fluxes (soil, puddle and total) noting that these are referenced
! to the area of the land point (not to the area of the puddle/non-puddle split).
! * `pwet` quantifying the impact of puddles in reducing the latent heat flux from the soil.
! * `ssnow_wetfac` quantifying control of soil moisture on the rate of evaporation.
! * `ssnow_cls` quantifying whether the latent heat flux represents a change between
! liquid/vapour of ice/vapour phases. 
!
! **IMPORTANT** the value of `ssnow_cls` set in this SUBROUTINE is used to control whether
! the water fluxes from the snow pack/soil column in [[soil_snow]] are of water or ice.
!
!<br><\br>
!## Workflow
! The underpinning equation linking the \(\beta\)-formulation for the soil latent heat flux
! to the potential evaporation is
!
! \(\lambda E_s = \beta c_{ls} \lambda E_{pot} \)
!
! where \(\lambda\) =`air_rlam` is the latent heat of evaporation, \(c_{ls}\)=`ssnow_cls`
! is a factor set by whether the flux of water is to/from water or to/from ice
! (\(c_{ls}=1.1335\) if ice - 1.1335 is the ratio of the latent heat of sublimation of water
! to the latent heat of evaporation of water), \(\beta\)=`ssnow_wetfac` is a coefficient controlling
! the availability of water in the soil surface layer and \(\lamba E_{pot}\)=`ssnow_potev`
! is the latent heat flux associated with potential evaporation.
!
! When the latent heat flux is negative, i.e. the flux of water is from the atmosphere
! onto the soil/snow there is no limit due to availability of water in the soil column
! and `ssnow_wetfac` is overwritten to a value of 1.
! Otherwise a first estimate for the soil latent heat flux is calcualted as above without
! further limits.
WHERE (ssnow_potev < 0. ) ssnow_wetfac(:) = 1.0
canopy_fess= ssnow_wetfac * ssnow_potev

!> When there is puddle a fraction of the soil column is covered by water, `pwet`.
! As a result the evaporation from the soil is reduced by the change in area fraction.
pwet = MAX(0.,MIN(0.2,ssnow_pudsto/MAX(1.,ssnow_pudsmx)))
canopy_fess = canopy_fess * (1.-pwet)

! frescale is a factor used to convert an amount of water (in m3/m3)
! in the surface layer of the soil into a corresponding latent heat flux.
! 1000 is the density of water in kg/m3
frescale = soil_zse * 1000. * air_rlam / dels

!> The main workflow loops over the land points in the CABLE array and determines
! additional limits and the value for \(c_{ls}\) depending on which
! of the four cases that land point falls in (Ticket 137)
! Cases 1 and 3 represent a conversion between liquid and vapour phases of water.
DO j=1,mp

  IF(ssnow_snowd(j) < 0.1 .AND. canopy_fess(j) .GT. 0. ) THEN

     !>For cases 1 and 3 there must be sufficient water in the surface soil layer
     ! to provide the water for evaporation, this sets `fupper_limit`
     ! the maximum latent heat flux possible.  Two options
     ! are provided, set by `cable_user_l_new_reduce_soilevap`, which
     ! differ in the amount of water that must remain at the end of the time step.
     IF (.NOT.cable_user_l_new_reduce_soilevp) THEN
        flower_limit(j) = REAL(ssnow_wb(j))-soil_swilt(j)/2.0
     ELSE
        ! E.Kowalczyk 2014 - reduces the soil evaporation
        flower_limit(j) = REAL(ssnow_wb(j))-soil_swilt(j)
     ENDIF
     fupper_limit(j) = MAX( 0.,                                        &
          flower_limit(j) * frescale(j)                       &
          - ssnow_evapfbl(j)*air_rlam(j)/dels)

     canopy_fess(j) = MIN(canopy_fess(j), REAL(fupper_limit(j),r_2))

     !> Additionally (for case 3) the evaporation of liquid water from within the
     ! frozen soil column must not reduce the liquid water fraction to the
     ! point that that ice fraction of soil moisture exceeds the
     ! frozen_limit=0.85.  This provides a second upper limit on the evaporation and
     ! latent flux. **WARNING** frozen_limit=0.85 is hard coded - if changed
     ! then the corresponding limit in [[cbl_soilsnow]] must also be changed.
     fupper_limit(j) = REAL(ssnow_wb(j)-ssnow_wbice(j)/0.85)*frescale(j)
     fupper_limit(j) = MAX(REAL(fupper_limit(j),r_2),0.)

     canopy_fess(j) = MIN(canopy_fess(j), REAL(fupper_limit(j),r_2))

  END IF   

  ssnow_cls(j)=1.

  !> If the surface has snow cover then the latent heat flux represents a
  ! conversion between solid and vapour phases of water and (\(c_{ls}\)=1.1335.
  ! The first estimate is updated by the change in value of (\(c_{ls}\).
  IF (ssnow_snowd(j) >=0.1 ) THEN
     ssnow_cls(j) = 1.1335
     canopy_fess(j) = ssnow_cls(j)*ssnow_potev(j)
  ENDIF

  !> For case 4 deposition of frost onto frozen soil (temperature below freezing)
  ! occurs even if there is no snow (so \(c_{ls}\)=1.1335) and there is no limit.  
  IF (ssnow_snowd(j) < 0.1 .AND. ssnow_potev(j) < 0. .AND. &
       ssnow_tss(j)<CTFRZ) THEN
     ssnow_cls(j)=1.1335
     canopy_fess(j) = ssnow_cls(j)*ssnow_potev(j)
  ENDIF

  !> For case 2, if `ssnow_potev`>0 then there needs to be sufficient snow to
  ! sublimate over the time step.  This places an upper limit on the latent
  ! heat flux possible.
  IF (ssnow_snowd(j) >= 0.1 .AND. ssnow_potev(j) > 0.) THEN

     ssnow_cls(j) = 1.1335

     !INH - if changes to PM routine then matching changes here
     canopy_fess(j) = MIN( (ssnow_wetfac(j)*ssnow_potev(j))*ssnow_cls(j), &
          ssnow_snowd(j)/dels*air_rlam(j)*ssnow_cls(j))

  ENDIF

ENDDO
!*  Note that the case of dew fall onto a surface with little/no snow while
! the soil surface temperature is above freezing is unrestricted and
! the latent heat flux defaults to the first estimate.

!> The latent heat flux associated with evaporation from puddles is set
! to the area fraction of `ssnow_potev` but then limited by the amount of water
! available in the puddles.
! Finally the total latant heat lux is obtained by summing the soil and
! puddle contributions.
canopy_fesp = MIN(ssnow_pudsto/dels*air_rlam,MAX(pwet*ssnow_potev,0.))
canopy_fes = canopy_fess + canopy_fesp

END SUBROUTINE latent_heat_flux

END MODULE cable_latent_heat_module
