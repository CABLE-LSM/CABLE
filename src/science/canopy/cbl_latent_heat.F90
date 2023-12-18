MODULE cable_latent_heat_module
  !* This MODULE contains one SUBROUTINE which evaluates the latent heat from
  ! ground/soil/snow pack given the previously evaluated rate of potential
  ! evaporation. This is a component of the calculation of the surface energy
  ! balance and is called twice per cycle in the MO iteration section of
  ! [[define_canopy]].

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
  ! evaporation `ssnow_potev` from the soil into the associated latent heat flux
  ! `canopy_fes`. The evaluation accounts for whether the water flux is
  ! in liquid or solid phase (evaporation or sublimation, dew or frost),
  ! the direction of the flux (to or from the soil), and any limits due to a lack of
  ! water in the surface layer of the soil column.  The SUBROUTINE also evaluates
  ! the partition of the latent heat flux between soil `canopy_fess` and
  ! puddle `canopy_fesp` contributions.
  !
  !## Warning
  !
  ! This current CABLE implementation is possibly confusing.
  ! The variable `ssnow_potev` is not potential evaporation (in kg water
  ! m\(^{-2}\)s\(^{-1}\)) but is instead the latent heat flux (in Wm\(^{-2}\))
  ! associated with the potential evaporation.
  ! `ssnow_potev` assumes that the water flux described by the potential
  ! evaporation is of liquid water regardless of temperature
  ! and/or whether the surface is covered by snow. Additionally one option for the
  ! evaluation of `ssnow_potev` (the [[Penman_Monteith]] option) does not function
  ! correctly over frozen surfaces. A bug fix will require changes to
  ! [[define_canopy]], [[Penman_Monteith]] and this SUBROUTINE; comments have been
  ! inserted into the code accordingly.  
  !
  ! The subroutine implicitly sets `ssnow_wetfac=1` over snow but does not 
  ! actually change the value of `ssnow_wetfac` for these points.
  ! Later on, other calculations use `ssnow_wetfac` and may use a different
  ! `ssnow_wetfac` value over snow. It is important to note,
  ! [[Surf_wetness_fact]] calculates `ssnow_wetfac` and sets it to 0.9
  ! over snow. And the value for points with new snow could be smaller.


USE cable_def_types_mod, ONLY : r_2
IMPLICIT NONE

INTEGER :: mp
REAL(KIND=r_2), INTENT(OUT) :: canopy_fes(mp)
!! latent heat flux from the ground (Wm\(^{-2}\)) 
REAL(KIND=r_2), INTENT(OUT) :: canopy_fess(mp)
!! latent heat flux from the soil (Wm\(^{-2}\))
REAL(KIND=r_2), INTENT(OUT) :: canopy_fesp(mp)
!! latent heat flux from any puddles (Wm\(^{-2}\))
REAL, INTENT(OUT) :: ssnow_cls(mp)
!! factor denoting phase of water flux (=1 if liquid, =1.1335 if ice)
REAL, INTENT(IN OUT) :: pwet(mp)
!! factor to reduce soil evaporation due to presence of a puddle (-)
REAL, INTENT(IN OUT) :: ssnow_wetfac(mp)
!! wetness factor for soil (between 0 and 1)


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

!|## Method
!
! 'Potential evaporation' quantifies the theoretical value of evaporation from
! a fully saturated surface (i.e. with no limits on water availability) at a given
! temperature and humidity, in response to the available energy.
! This SUBROUTINE evaluates the value of the latent heat flux from the
! ground/soil/snow pack given the potential evaporation according to one of four cases
! (see Ticket 137)
!
!  1. evaporation from/dew onto surfaces where there is no snow and
!     when soil temperature > freezing point
!  2. sublimation from/deposition of frost onto surfaces where there is snow cover
!  3. evaporation of liquid water from a frozen soil column if there is no snow present
!     (sublimation of frozen water from within the soil column is not possible per CABLE)
!  4. deposition of frost onto frozen soils if there is no snow cover
!
! while satisfying any constraints due to limits from soil moisture availability.
!<br></br>
! There are six outputs from this SUBROUTINE
!
! * three latent heat fluxes (soil, puddle and total). Note that these fluxes
! are referenced to the area of the land point and not to the area of the
! puddle/non-puddle split.
! * `pwet` quantifying the impact of puddles in reducing the latent heat flux from the soil.
! * `ssnow_wetfac` quantifying the control of soil moisture on the rate of evaporation.
! * `ssnow_cls` quantifying whether the latent heat flux represents a change between
! liquid/vapour of ice/vapour phases.
!<br></br>
! `ssnow_cls` takes one of two values - 1.0 if the flux of water is of liquid, and
! 1.1335 if the flux of water is of ice.  1.1335 = latent heat of sublimation of ice /
! latent heat of evaporation of liquid water.  
!
! **IMPORTANT** the value of `ssnow_cls` set in this SUBROUTINE is used to control whether
! the water fluxes from the snow pack/soil column are of liquid or solid water
! in [[soil_snow]] (and its SUBROUTINES).
! Inconsistencies between the two sets of SUBROUTINES
! would lead to loss of energy closure and/or loss of mass.
!
!<br></br>
!## Workflow
! The underpinning equation linking the \(\beta\)-formulation for the soil latent heat flux
! to the potential evaporation is
!
! \(\lambda E_s = \beta c_{ls} \lambda E_{pot} \)
!
! where \(\lambda\) =`air_rlam` is the latent heat of evaporation
! (which can theoretically vary so is an mp-array), \(c_{ls}\)=`ssnow_cls`
! is a factor set by the phase of the water being fluxed,
! \(\beta\)=`ssnow_wetfac` is a coefficient quantifying the control that soil surface
! layer moisture has on evaporation, and \(\lambda E_{pot}\)=`ssnow_potev`
! is the latent heat flux associated with the potential evaporation.
!
! <br></br>
! The workflow is
!
! 1. A first estimate for the soil latent heat flux is calculated
! assuming that the flux is of liquid water. If `ssnow_potev` is negative
! any limitation due to soil moisture availability is negected.
! <br></br>
WHERE (ssnow_potev < 0. ) ssnow_wetfac(:) = 1.0
canopy_fess= ssnow_wetfac * ssnow_potev

!| 2. If there is water in the puddle store `ssnow_pudsto`, the fraction
! of surface area covered by the puddle is quantified, `pwet` (between 0-0.2).
! The soil latent heat flux is reduced by the change in area fraction.
! <br></br>
pwet = MAX(0.,MIN(0.2,ssnow_pudsto/MAX(1.,ssnow_pudsmx)))
canopy_fess = canopy_fess * (1.-pwet)

! frescale is a factor used to convert an amount of water (in m3/m3)
! in the surface layer of the soil into a limit on the soil latent heat flux.
! 1000 is the density of water in kg/m3
frescale = soil_zse * 1000. * air_rlam / dels

!| 3. (the main loop) The value for \(c_{ls}\) and additional limits
! on the latent heat flux(es) are applied, according to which of the four cases
! is occurring.  Inside the loop the workflow is as follows:
!
DO j=1,mp

  IF(ssnow_snowd(j) < 0.1 .AND. canopy_fess(j) .GT. 0. ) THEN

!|    - For cases 1 and 3 the flux is of liquid water, \(c_{ls} =1\).
!     
!     - Limit 1: For cases 1 and 3 there must be sufficient liquid water in the surface
!       soil layer to provide the water for evaporation. This sets
!       the maximum soil latent heat flux possible, `fupper_limit`.  Two options
!       are provided for `fupper_limit`, set by the `cable_user_l_new_reduce_soilevap`
!       switch.
!       The options differ in the amount of water that remains at the end of the time step.
!
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

!|     - Limit 2: Additionally for case 3, the evaporation of liquid water
!        from within the frozen soil column must not reduce the liquid water fraction
!        to the point that that ice fraction of soil moisture exceeds an upper limit
!        frozen_limit=0.85.  This provides a second upper limit on the evaporation and
!        soil latent flux. **WARNING** frozen_limit=0.85 is hard coded - if it is changed
!        then the corresponding limit in [[cbl_soilsnow]] must also be changed.
!
     fupper_limit(j) = REAL(ssnow_wb(j)-ssnow_wbice(j)/0.85)*frescale(j)
     fupper_limit(j) = MAX(REAL(fupper_limit(j),r_2),0.)

     canopy_fess(j) = MIN(canopy_fess(j), REAL(fupper_limit(j),r_2))

  END IF   
!*     - The case of dew fall onto a surface with little/no snow while
!        the soil surface temperature is above freezing is unrestricted and
!        the soil latent heat flux defaults to the first estimate.
  
  ssnow_cls(j)=1.

!|     - If the surface has snow cover, or `ssnow_potev`<0 **and** the soil
!        temperature is below freezing (frost), then the latent heat flux represents a
!        conversion between solid and vapour phases of water, \(c_{ls}\)=1.1335.
!        The first estimate for the soil latent heat flux is updated by the change
!        in value of \(c_{ls}\).
!
  IF (ssnow_snowd(j) >=0.1 ) THEN
     ssnow_cls(j) = 1.1335
     canopy_fess(j) = ssnow_cls(j)*ssnow_potev(j)
  ENDIF

!|     - For case 4, deposition of frost onto frozen surface (temperature below
!        freezing) occurs even if there is no snow - there is no limit
!        on the magnitude of the soil latent heat flux.
!
  IF (ssnow_snowd(j) < 0.1 .AND. ssnow_potev(j) < 0. .AND. &
       ssnow_tss(j)<CTFRZ) THEN
     ssnow_cls(j)=1.1335
     canopy_fess(j) = ssnow_cls(j)*ssnow_potev(j)
  ENDIF

!|     - Limit 3: For case 2, if `ssnow_potev`>0 then there needs to be sufficient
!        snow to sublimate over the time step.  This places an upper limit on
!        the soil latent heat flux.
!
  IF (ssnow_snowd(j) >= 0.1 .AND. ssnow_potev(j) > 0.) THEN

     ssnow_cls(j) = 1.1335

     !INH - if changes to PM routine then matching changes here
     canopy_fess(j) = MIN( (ssnow_wetfac(j)*ssnow_potev(j))*ssnow_cls(j), &
          ssnow_snowd(j)/dels*air_rlam(j)*ssnow_cls(j))

  ENDIF

ENDDO

!|  4. The latent heat flux associated with evaporation from puddles is set
! to the area fraction of the potential evaporation (`pwet` * `ssnow_potev`).
! If there is insufficient water in the puddle to support this flux then an upper
! limit is applied.
! <br></br>
canopy_fesp = MIN(ssnow_pudsto/dels*air_rlam,MAX(pwet*ssnow_potev,0.))

!| 5. The total latent heat flux is obtained by summing the soil and
! puddle contributions. 
canopy_fes = canopy_fess + canopy_fesp

END SUBROUTINE latent_heat_flux

END MODULE cable_latent_heat_module
