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
  ! the partition of the latent heat flux between soil and puddle contributions.
  !
  !## Warning
  !
  ! This current implementation is possibly confusing.  The variable `ssnow_potev`
  ! is not potential evaporation (in kg water m\(^{-2}\)s\(^{-1}\)) but the
  ! latent heat flux (in Wm\(^{-2}\) associated with the potential evaporation.
  ! `ssnow_potev` also assumes that the water flux described by potential
  ! evaporation is of liquid water regardless of temperature
  ! and whether the surface is covered by snow. Additionally one option for the
  ! evlaution of `ssnow_potev` (the [[Penman_Monteith]] option) does not function
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
EAL, INTENT(IN) :: soil_zse                   !! thickness of topmost soil layer (m)
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

!*## Method
!
! 'Potential evaporation' quantifies the theoretical value of evaporation from
! a fully saturated surface (with no limits of water) at a given temperature
! and avaialble energy.
! This SUBROUTINE evaluates a true value of the ground latent heat flux in one of four cases
! (Ticket 137)
!
!  1. evaporation from/dew onto surfaces with no snow, where soil temperature > freezing point
!  2. sublimation from/frost onto surfaces with snow cover
!  3. evaporation of liquid water from frozen soil column if there is no snow present
!     (sublimation of frozen water from within the soil column is not possible)
!  4. deposition of frost onto a frozen soils if no snow cover
!
!IMPORTANTLY the value of _cls set here is used to control whether
!water fluxes are from the snow pack or soil column in _soilsnow

! Soil latent heat:
WHERE (ssnow_potev < 0. ) ssnow_wetfac(:) = 1.0
canopy_fess= ssnow_wetfac * ssnow_potev

! Reduce soil evap due to presence of puddle
pwet = MAX(0.,MIN(0.2,ssnow_pudsto/MAX(1.,ssnow_pudsmx)))
canopy_fess = canopy_fess * (1.-pwet)

frescale = soil_zse * 1000. * air_rlam / dels

DO j=1,mp

  IF(ssnow_snowd(j) < 0.1 .AND. canopy_fess(j) .GT. 0. ) THEN

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

     !Ticket 137 - case iii)
     !evaporation from frozen soils needs to respect the assumption that
     !ice fraction of soil moisture cannot exceed frozen_limit=0.85
     !see soilsnow: if frozen_limit changes need to be consistent
     fupper_limit(j) = REAL(ssnow_wb(j)-ssnow_wbice(j)/0.85)*frescale(j)
     fupper_limit(j) = MAX(REAL(fupper_limit(j),r_2),0.)

     canopy_fess(j) = MIN(canopy_fess(j), REAL(fupper_limit(j),r_2))

  END IF

  ssnow_cls(j)=1.

  !Ticket 137 - case ii) deposition of frost onto snow
  ! case of sublimation of snow overwrites later
  IF (ssnow_snowd(j) >=0.1 ) THEN
     ssnow_cls(j) = 1.1335
     canopy_fess(j) = ssnow_cls(j)*ssnow_potev(j)
  ENDIF

  !Ticket 137 - case iv) deposition of frost onto frozen soil, no snow
  IF (ssnow_snowd(j) < 0.1 .AND. ssnow_potev(j) < 0. .AND. &
       ssnow_tss(j)<CTFRZ) THEN
     ssnow_cls(j)=1.1335
     canopy_fess(j) = ssnow_cls(j)*ssnow_potev(j)
  ENDIF

  !Ticket 137 - case ii) sublimation of snow
  IF (ssnow_snowd(j) >= 0.1 .AND. ssnow_potev(j) > 0.) THEN

     ssnow_cls(j) = 1.1335

     !INH - if changes to PM routine then matching changes here
     canopy_fess(j) = MIN( (ssnow_wetfac(j)*ssnow_potev(j))*ssnow_cls(j), &
          ssnow_snowd(j)/dels*air_rlam(j)*ssnow_cls(j))

  ENDIF

ENDDO

! Evaporation from soil puddle
canopy_fesp = MIN(ssnow_pudsto/dels*air_rlam,MAX(pwet*ssnow_potev,0.))
canopy_fes = canopy_fess + canopy_fesp

END SUBROUTINE latent_heat_flux

END MODULE cable_latent_heat_module
