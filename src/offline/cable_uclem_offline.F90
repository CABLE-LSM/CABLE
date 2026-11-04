!==============================================================================
! UCLEM urban canopy model for the offline drivers (serial and MPI worker).
!
! UCLEM runs on the urban patches (veg%iveg == urban_cable) after cbm has
! advanced the whole mp vector.  CABLE still integrates its own soil and
! vegetation on those patches, as on the ACCESS3 urban tile, but the fluxes,
! radiative temperature, albedo and runoff that reach the output are replaced
! by the UCLEM values.
!
! UCLEM is switched on and parameterised by the namelist group uclem_nml in the
! CABLE namelist file.  When the group is absent, or uclem_on is false, nothing
! in this module changes the CABLE state.
!
! The coupling routines in uclem_coupling_mod read the CABLE forcing through one
! contiguous slice, which matches the ACCESS3 tile-major packing.  Offline, mp
! is ordered by land point, so the urban patches of different land points are
! separated by the other patches.  This module therefore gathers the urban
! patches, through the index array upts, into nu-length copies of the met,
! radiation, roughness and air types, calls the coupling routines on those
! copies as one contiguous slice, and scatters the results back.
!==============================================================================
MODULE cable_uclem_offline_mod

USE cable_def_types_mod, ONLY: met_type, air_type, radiation_type,             &
                               roughness_type, canopy_type, soil_snow_type,    &
                               veg_parameter_type, alloc_cbm_var
USE uclem_state_mod,     ONLY: uclem_type, uclem_data_type,                    &
                               alloc_uclem_type, assoc_uclem_type,             &
                               update_uclem_geometry, update_uclem_optical,    &
                               normalise_uclem_road_volcp,                     &
                               rebuild_uclem_derived_state
USE uclem_science,       ONLY: urbtemp

IMPLICIT NONE

PRIVATE
PUBLIC :: uclem_offline_init, uclem_offline_step, uclem_on

! Namelist uclem_nml.  A value <= 0 keeps the UCLEM default for urban type 1.
LOGICAL :: uclem_on  = .FALSE.
REAL    :: bldheight = -1.0   ! building height (m)
REAL    :: hwratio   = -1.0   ! building height to canyon width ratio
REAL    :: sigmabld  = -1.0   ! plan area fraction of buildings
REAL    :: sigvegc   = -1.0   ! plan area fraction of canyon vegetation
REAL    :: alb_roof  = -1.0
REAL    :: alb_wall  = -1.0
REAL    :: alb_road  = -1.0
REAL    :: alb_vegc  = -1.0
REAL    :: emis_roof = -1.0
REAL    :: emis_wall = -1.0
REAL    :: emis_road = -1.0
REAL    :: emis_vegc = -1.0

NAMELIST /uclem_nml/ uclem_on, bldheight, hwratio, sigmabld, sigvegc,          &
                     alb_roof, alb_wall, alb_road, alb_vegc,                   &
                     emis_roof, emis_wall, emis_road, emis_vegc

TYPE(uclem_data_type), TARGET, SAVE :: uclem_data
TYPE(uclem_type),              SAVE :: uclem
INTEGER, ALLOCATABLE,          SAVE :: upts(:)   ! mp index of each urban patch
INTEGER,                       SAVE :: nu = 0
LOGICAL,                       SAVE :: initialised = .FALSE.
LOGICAL,                       SAVE :: first_step  = .TRUE.

! Packed nu-length copies of the CABLE types read by the coupling routines.
TYPE(met_type),       SAVE :: umet
TYPE(radiation_type), SAVE :: urad
TYPE(roughness_type), SAVE :: urough
TYPE(air_type),       SAVE :: uair

CONTAINS

!==============================================================================
! Read uclem_nml, locate the urban patches and set up the UCLEM state.
! Called once, after veg%iveg is known on this process.
!==============================================================================
SUBROUTINE uclem_offline_init(nml_file, mp, veg, logn)

USE cable_surface_types_mod, ONLY: urban_cable

CHARACTER(LEN=*),           INTENT(IN) :: nml_file
INTEGER,                    INTENT(IN) :: mp, logn
TYPE(veg_parameter_type),   INTENT(IN) :: veg

INTEGER :: unit, ios, i
REAL, ALLOCATABLE :: wrr(:)

IF (initialised) RETURN
initialised = .TRUE.

OPEN(NEWUNIT=unit, FILE=nml_file, STATUS='OLD', ACTION='READ', IOSTAT=ios)
IF (ios /= 0) RETURN
READ(unit, NML=uclem_nml, IOSTAT=ios)
CLOSE(unit)
! A missing group leaves UCLEM off.  A malformed group must not pass silently.
IF (ios > 0) STOP 'UCLEM: error reading namelist uclem_nml'
IF (ios < 0) uclem_on = .FALSE.
IF (.NOT. uclem_on) RETURN

nu = COUNT(veg%iveg(1:mp) == urban_cable)
IF (nu == 0) THEN
  ! An MPI worker can hold no urban patches while others do.
  WRITE(logn,*) 'UCLEM: no urban patches on this process'
  uclem_on = .FALSE.
  RETURN
END IF
upts = PACK([(i, i = 1, mp)], veg%iveg(1:mp) == urban_cable)

CALL alloc_cbm_var(umet,   nu)
CALL alloc_cbm_var(urad,   nu)
CALL alloc_cbm_var(urough, nu)
CALL alloc_cbm_var(uair,   nu)

CALL alloc_uclem_type(uclem_data, nu)
CALL assoc_uclem_type(uclem, uclem_data)

! update_uclem_geometry takes the road width fraction, 1 - sigmabld.
ALLOCATE(wrr(nu))
wrr = -1.0
IF (sigmabld > 0.0) wrr = 1.0 - sigmabld
CALL update_uclem_geometry(uclem_data, SPREAD(bldheight, 1, nu),               &
                           SPREAD(hwratio, 1, nu), wrr, nu)
DEALLOCATE(wrr)

! update_uclem_optical applies sigvegc only where alb_roof > 0.
CALL update_uclem_optical(uclem_data,                                          &
                          SPREAD(alb_roof, 1, nu),  SPREAD(alb_wall, 1, nu),   &
                          SPREAD(alb_road, 1, nu),  SPREAD(sigvegc, 1, nu),    &
                          SPREAD(alb_vegc, 1, nu),  SPREAD(emis_roof, 1, nu),  &
                          SPREAD(emis_wall, 1, nu), SPREAD(emis_road, 1, nu),  &
                          SPREAD(emis_vegc, 1, nu), nu)
CALL normalise_uclem_road_volcp(uclem_data, nu)
CALL rebuild_uclem_derived_state(uclem_data, nu)

! The patch fraction already carries the urban area of the grid cell, so the
! whole urban patch is urban.
uclem%fp%sigmau = 1.0

WRITE(logn,*) 'UCLEM: on for ', nu, ' urban patches, first and last mp index ', &
              upts(1), upts(nu)
WRITE(logn,*) 'UCLEM: bldheight ', uclem%fp%bldheight(1),                     &
              ' hwratio ', uclem%fp%hwratio(1),                               &
              ' sigmabld ', uclem%fp%sigmabld(1)

END SUBROUTINE uclem_offline_init

!==============================================================================
! Advance UCLEM one step and replace the urban patch output fields.
! Called right after cbm, while the runoff fields are still rates (kg/m2/s).
!==============================================================================
SUBROUTINE uclem_offline_step(mp, dels, met, air, rad, rough, canopy, ssnow)

USE uclem_coupling_mod,       ONLY: run_uclem_solar_geometry,                  &
                                    run_uclem_beam_fraction,                   &
                                    run_uclem_access_coupling
USE cable_phys_constants_mod, ONLY: CSBOLTZ => SBOLTZ

INTEGER,               INTENT(IN)    :: mp
REAL,                  INTENT(IN)    :: dels
TYPE(met_type),        INTENT(IN)    :: met
TYPE(air_type),        INTENT(IN)    :: air
TYPE(radiation_type),  INTENT(INOUT) :: rad
TYPE(roughness_type),  INTENT(IN)    :: rough
TYPE(canopy_type),     INTENT(INOUT) :: canopy
TYPE(soil_snow_type),  INTENT(INOUT) :: ssnow

INTEGER :: j
REAL    :: utc_hour
REAL    :: wind_v(nu)

IF (.NOT. uclem_on) RETURN

! Gather the forcing that the coupling routines read.
umet%tk        = met%tk(upts)
umet%qv        = met%qv(upts)
umet%ua        = met%ua(upts)
umet%pmb       = met%pmb(upts)
umet%fld       = met%fld(upts)
umet%precip    = met%precip(upts)
umet%precip_sn = met%precip_sn(upts)
umet%coszen    = met%coszen(upts)
umet%fsd       = met%fsd(upts,:)
urad%fbeam     = rad%fbeam(upts,:)
urad%latitude  = rad%latitude(upts)
urad%longitude = rad%longitude(upts)
urough%za_uv   = rough%za_uv(upts)
uair%rho       = air%rho(upts)

! Start the outdoor facets at the air temperature of the first step.  The
! indoor room, slab and internal mass keep the UCLEM default.
IF (first_step) THEN
  first_step = .FALSE.
  DO j = 0, UBOUND(uclem%roof%nodetemp, 2)
    uclem%roof%nodetemp(:,j)  = REAL(umet%tk - urbtemp, 8)
    uclem%road%nodetemp(:,j)  = REAL(umet%tk - urbtemp, 8)
    uclem%walle%nodetemp(:,j) = REAL(umet%tk - urbtemp, 8)
    uclem%wallw%nodetemp(:,j) = REAL(umet%tk - urbtemp, 8)
  END DO
END IF

! The offline met%hod is local solar time, set from the longitude when the
! forcing is in GMT.  UCLEM takes UTC and adds the longitude back itself.  UTC
! is the same at every patch, so it is taken from the first urban patch, and so
! is the day of year, which only sets the solar declination.
utc_hour = MODULO(met%hod(upts(1)) - rad%longitude(upts(1))/15.0, 24.0)

! Offline radiation is evaluated every step, so the sun position is updated
! every step, unlike the ACCESS3 path where it follows the radiation calls.
CALL run_uclem_solar_geometry(nu, umet%coszen, INT(met%doy(upts(1))), utc_hour, &
                              1, nu, urad, uclem)
CALL run_uclem_beam_fraction(nu, 1, nu, umet, urad, uclem)

! The offline forcing keeps only the wind speed, so the wind is taken to blow
! along the x axis.  run_uclem_access_coupling also forms umet%ofsd.
wind_v = 0.0
CALL run_uclem_access_coupling(nu, dels, 1, nu, umet%ua, wind_v,               &
                               umet, uair, urough, uclem, 0)

! Replace the urban patch fields read by the output.  Outgoing longwave is
! sigma*u_ts**4 because UCLEM defines u_ts from the outgoing longwave.
canopy%fh(upts)   = uclem%u_fg
canopy%fe(upts)   = uclem%u_eg
rad%swnet(upts)   = umet%ofsd * (1.0 - uclem%pd%u_alb)
rad%lwnet(upts)   = umet%fld - CSBOLTZ * uclem%u_ts**4
rad%rnet(upts)    = rad%swnet(upts) + rad%lwnet(upts)
canopy%rnet(upts) = rad%rnet(upts)
! Storage heat flux, net of the anthropogenic heat, as the residual, which keeps
! Qg = Rnet - Qh - Qle in the output.  The CABLE energy balance check (bal%ebal)
! is built from the vegetation and soil terms and does not close on the urban
! patch.
canopy%ga(upts)   = rad%rnet(upts) - canopy%fh(upts) - canopy%fe(upts)
rad%trad(upts)    = uclem%u_ts
rad%albedo(upts,1) = uclem%pd%u_alb
rad%albedo(upts,2) = uclem%pd%u_alb
rad%albedo_T(upts) = uclem%pd%u_alb
! u_rn is the runoff over the step (kg/m2).  UCLEM has no deep drainage.
ssnow%rnof1(upts)  = uclem%u_rn / dels
ssnow%rnof2(upts)  = 0.0
ssnow%runoff(upts) = ssnow%rnof1(upts)

END SUBROUTINE uclem_offline_step

END MODULE cable_uclem_offline_mod
