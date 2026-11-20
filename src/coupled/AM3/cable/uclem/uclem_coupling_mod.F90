MODULE uclem_coupling_mod

IMPLICIT NONE

PRIVATE
PUBLIC :: run_uclem_beam_fraction, run_uclem_access_coupling
PUBLIC :: run_uclem_solar_geometry, run_uclem_albedo_coupling

CONTAINS

!===============================================================================
! Set the UCLEM solar geometry (zenith, hour angle, local time).
!
! Called only from run_uclem_albedo_coupling, i.e. only when the host radiation
! scheme runs.  uclem_eval then absorbs shortwave with the same sun position that
! produced the albedo the host was given.  CCAM calls uclem_ccangle only from
! seaesfrad for the same reason.
!===============================================================================
SUBROUTINE run_uclem_solar_geometry(mp, coszen, day_of_year, utc_hour, iu0, nu, rad, uclem)

USE uclem_state_mod, ONLY: calculate_solar_angles, uclem_newangle, uclem_type
USE cable_def_types_mod, ONLY : radiation_type

IMPLICIT NONE

INTEGER, INTENT(IN) :: mp
REAL,    INTENT(IN) :: coszen(mp)
INTEGER, INTENT(IN) :: day_of_year, iu0, nu
REAL,    INTENT(IN) :: utc_hour
TYPE (radiation_type), INTENT(INOUT) :: rad
TYPE (uclem_type),     INTENT(INOUT) :: uclem

INTEGER :: j, iu1
REAL, ALLOCATABLE :: azimuth(:)
REAL, ALLOCATABLE :: uclem_ctime(:)

IF (nu == 0) RETURN

! UCLEM arrays are allocated on the urban points only.  The mp-length CABLE
! forcing is read through the contiguous slice iu0:iu1 located by
! set_uclem_packing and verified by the probe in init_uclem_restart.
iu1 = iu0 + nu - 1

allocate(azimuth(nu))
allocate(uclem_ctime(nu))

! Local time as a fraction of the day, which is what uclem_newangle expects
! (its ctimein argument is documented as 0 <= ctime <= 1).  Held in a local
! rather than in met%hod: that field is named for an hour of day, and writing
! a day fraction into it would leave any future reader of CABLE's met type
! silently a factor of 24 out.  UCLEM is currently its only user.
uclem_ctime = MOD(utc_hour + rad%longitude(iu0:iu1)/15.0, 24.0) / 24.0

DO j = 1, nu
  CALL calculate_solar_angles(rad%latitude(iu0+j-1), rad%longitude(iu0+j-1), &
                              REAL(day_of_year), utc_hour, azimuth(j))
END DO

call uclem_newangle(uclem,nu,coszen(iu0:iu1),uclem_ctime,azimuth)

deallocate(azimuth)
deallocate(uclem_ctime)

RETURN
END SUBROUTINE run_uclem_solar_geometry


!===============================================================================
! Return the UCLEM canyon albedo to CABLE's effective surface reflectance, i.e.
! the quantity that cable_rad_unpack passes to the UM radiation scheme as
! alb_surft/land_albedo.  Called from cable_land_albedo on the radiation
! pathway, after cable_rad_driver has filled EffSurfRefl_* with the
! vegetation/soil/snow values and before they are unpacked.
!
! Notes
!  - split=1 and split=2 return the direct-beam and diffuse canyon albedo.
!    uclem_alb overrides fp%fbeam internally for these two modes, so the beam
!    fraction does not need to be set here.
!  - UCLEM (and U-Surf) carry a single broadband facet albedo, so the same
!    value is applied to the VIS and NIR bands.  Band 3 is the legacy LW slot
!    and is left untouched.
!  - sigmau is used only as an urban point mask, not as an area weight; the
!    tile fraction weighting is applied later by the host.
!===============================================================================
SUBROUTINE run_uclem_albedo_coupling(mp, nrb, coszen, day_of_year, utc_hour,  &
                                     iu0, nu, rad, uclem,                    &
                                     EffSurfRefl_beam, EffSurfRefl_dif)

USE uclem_state_mod, ONLY: uclem_alb, uclem_type
USE cable_def_types_mod, ONLY : radiation_type

IMPLICIT NONE

INTEGER, INTENT(IN) :: mp, nrb
REAL,    INTENT(IN) :: coszen(mp)
INTEGER, INTENT(IN) :: day_of_year, iu0, nu
REAL,    INTENT(IN) :: utc_hour
TYPE (radiation_type), INTENT(INOUT) :: rad
TYPE (uclem_type),     INTENT(INOUT) :: uclem
REAL, INTENT(INOUT) :: EffSurfRefl_beam(mp,nrb)
REAL, INTENT(INOUT) :: EffSurfRefl_dif(mp,nrb)

INTEGER :: iu1
REAL, ALLOCATABLE :: urban_alb_direct(:)
REAL, ALLOCATABLE :: urban_alb_diffuse(:)

CALL run_uclem_solar_geometry(mp, coszen, day_of_year, utc_hour, iu0, nu, rad, uclem)

IF (nu == 0) RETURN

iu1 = iu0 + nu - 1

! No first-call guard is needed here any more.  cable_land_albedo packs
! rad%latitude and rad%longitude before this call, and run_uclem_solar_geometry
! takes the hour of day from the model clock, so the sun position is defined on
! the first radiation step of every run and every restart.
allocate(urban_alb_direct(nu))
allocate(urban_alb_diffuse(nu))

call uclem_alb(uclem, nu, urban_alb_direct,  split=1)
call uclem_alb(uclem, nu, urban_alb_diffuse, split=2)

urban_alb_direct  = min(max(urban_alb_direct,  0.0), 1.0)
urban_alb_diffuse = min(max(urban_alb_diffuse, 0.0), 1.0)

where (uclem%fp%sigmau > 0.0)
  EffSurfRefl_beam(iu0:iu1,1) = urban_alb_direct
  EffSurfRefl_beam(iu0:iu1,2) = urban_alb_direct
  EffSurfRefl_dif(iu0:iu1,1)  = urban_alb_diffuse
  EffSurfRefl_dif(iu0:iu1,2)  = urban_alb_diffuse
end where

deallocate(urban_alb_direct)
deallocate(urban_alb_diffuse)

RETURN
END SUBROUTINE run_uclem_albedo_coupling


!===============================================================================
! Hand UCLEM the beam fraction of the downward shortwave for the current step.
!
! This routine deliberately does not touch the solar geometry.  The sun position
! belongs to the radiation pathway: run_uclem_albedo_coupling sets vangle, hangle
! and ctime when the host radiation scheme is called, and uclem_eval then absorbs
! shortwave with the same sun position that produced the albedo the host was
! given.  Re-setting the angles here, once per model timestep, made the two
! disagree between radiation calls and broke the shortwave budget by the
! difference between the two albedos.  CCAM, the reference host for UCLEM, calls
! uclem_ccangle only from seaesfrad for the same reason.
!
! The beam fraction is still set every model timestep because rad%fbeam and
! met%fsd are packed from the host radiation fields on the explicit pathway and
! are current there, whereas work%sw_down_ij is filled after surf_couple_radiation
! has already run.  The ratio is independent of the cosine-of-zenith rescaling,
! so this is the current radiation step's value either way.
!===============================================================================
SUBROUTINE run_uclem_beam_fraction(mp, iu0, nu, met, rad, uclem)

USE uclem_state_mod, ONLY: uclem_fbeam, uclem_type
USE cable_def_types_mod, ONLY : met_type, radiation_type

IMPLICIT NONE

INTEGER, INTENT(IN) :: mp
INTEGER, INTENT(IN) :: iu0, nu
TYPE (met_type),       INTENT(IN)    :: met
TYPE (radiation_type), INTENT(IN)    :: rad
TYPE (uclem_type),     INTENT(INOUT) :: uclem

INTEGER :: iu1
REAL, ALLOCATABLE :: ufbeam(:)

IF (nu == 0) RETURN

iu1 = iu0 + nu - 1

allocate(ufbeam(nu))

! Weight the two shortwave bands by their downward fluxes rather than equally.
! The host radiation scheme applies the direct and diffuse canyon albedos to the
! actual direct and diffuse fluxes, so a flux-weighted beam fraction is what
! makes the single blended albedo used inside uclem_eval consistent with it.
! CCAM forms the same quantity in seaesfrad as
! fbeamvis*swrsave + fbeamnir*(1-swrsave).
ufbeam = ( rad%fbeam(iu0:iu1,1)*met%fsd(iu0:iu1,1)                              &
         + rad%fbeam(iu0:iu1,2)*met%fsd(iu0:iu1,2) )                           &
       / MAX( met%fsd(iu0:iu1,1) + met%fsd(iu0:iu1,2), TINY(1.0) )
call uclem_fbeam(uclem,nu,ufbeam)

deallocate(ufbeam)

RETURN
END SUBROUTINE run_uclem_beam_fraction


! Note: canopy is no longer an argument.  It was only ever written by the roughness sync
! that this routine used to perform, which had no effect on the atmosphere.  The urban
! aerodynamic override now happens in cable_expl_unpack.
SUBROUTINE run_uclem_access_coupling(mp, dels, iu0, nu, wind_u, wind_v,       &
                                     met, air, rough, uclem, progcalc)

USE uclem_science, ONLY: uclem_eval, uclem_prepinterior,                     &
                         refheight, urbtemp
USE uclem_state_mod, ONLY: uclem_type
USE cable_def_types_mod, ONLY : met_type, roughness_type, air_type

IMPLICIT NONE

INTEGER, INTENT(IN) :: mp, progcalc, iu0, nu
REAL,    INTENT(IN) :: dels
REAL,    INTENT(IN) :: wind_u(mp), wind_v(mp)

TYPE (met_type),       INTENT(INOUT) :: met
TYPE (roughness_type), INTENT(INOUT) :: rough
TYPE (uclem_type),     INTENT(INOUT) :: uclem
TYPE (air_type),       INTENT(IN)    :: air

INTEGER :: iu1

REAL, ALLOCATABLE :: a_udir(:)
REAL, ALLOCATABLE :: a_umag(:)
REAL, ALLOCATABLE :: a_umag_src(:)
REAL, ALLOCATABLE :: a_umagroof(:)
REAL, ALLOCATABLE :: a_temp_uclem(:)
REAL, ALLOCATABLE :: a_temproof_uclem(:)
REAL, ALLOCATABLE :: a_mixr_uclem(:)
REAL, ALLOCATABLE :: a_mixrroof_uclem(:)
REAL, ALLOCATABLE :: a_rnd_uclem(:)
REAL, ALLOCATABLE :: a_snd_uclem(:)
REAL, ALLOCATABLE :: a_zmin_uclem(:)
REAL, ALLOCATABLE :: a_zroof_uclem(:)

! met%ofsd is CABLE working state, so it is formed over the whole mp vector even
! though only the urban slice is passed on to uclem_eval.
met%ofsd = met%fsd(:,1) + met%fsd(:,2)

IF (nu == 0) RETURN

! UCLEM arrays are allocated on the urban points only.  The mp-length CABLE
! forcing is read through the contiguous slice iu0:iu1 located by
! set_uclem_packing and verified by the probe in init_uclem_restart.
iu1 = iu0 + nu - 1

allocate(a_udir(nu))
allocate(a_umag(nu))
allocate(a_umag_src(nu))
allocate(a_umagroof(nu))
allocate(a_temp_uclem(nu))
allocate(a_temproof_uclem(nu))
allocate(a_mixr_uclem(nu))
allocate(a_mixrroof_uclem(nu))
allocate(a_rnd_uclem(nu))
allocate(a_snd_uclem(nu))
allocate(a_zmin_uclem(nu))
allocate(a_zroof_uclem(nu))

a_umag_src = max(met%ua(iu0:iu1), 0.01)
a_udir = atan2(wind_v(iu0:iu1), wind_u(iu0:iu1))
a_rnd_uclem = MAX(met%precip(iu0:iu1) - met%precip_sn(iu0:iu1), 0.0) / MAX(dels, TINY(dels))
a_snd_uclem = MAX(met%precip_sn(iu0:iu1), 0.0) / MAX(dels, TINY(dels))

! UCLEM takes a reference height and a set of reference air properties for each of
! its two exchange paths.  The two heights are measured from different origins,
! refheight*bldheight for the canyon and bldheight for the roof, and inside
! uclem_eval they set the two path lengths pd%cndzmin and d_rfdzmin.  Those two
! path lengths, together with the separate canyon and roof roughness lengths and
! the separate reference pressures d_sigd and d_sigr, are what keeps the two
! facets aerodynamically distinct, so both heights are still required here.
! UCLEM expects absolute heights above ground and subtracts the displacement
! internally.  The MAX against bldheight is the single-level substitute for the
! level search CCAM performs in sflux; without it d_rfdzmin collapses to
! zoroof+0.2 whenever the buildings are taller than the lowest model level.
a_zmin_uclem  = MAX(rough%za_uv(iu0:iu1), refheight*uclem%fp%bldheight + 2.0)
a_zroof_uclem = MAX(rough%za_uv(iu0:iu1),           uclem%fp%bldheight + 2.0)

! Both paths are forced with the lowest UM level.  The two sets of air properties
! can only differ if the host resolves the layer between the canyon and roof
! origins, and the UM surface exchange hands down a single level.  CCAM behaves
! the same way whenever its own lowest level already clears bldheight + 2, in
! which case its level search returns one level for both paths.
!
! An earlier version interpolated this level to the two reference heights with
! uclem_most_profile, anchored on ssnow%tss and ssnow%qstss.  On the urban tile
! ssnow%tss is CABLE's surface temperature, which the override design leaves
! without physical meaning, and ssnow%qstss is the saturation humidity at that
! temperature rather than the actual surface humidity.  The interpolation placed
! roughly a quarter of its weight on that anchor for typical building heights,
! which inflated the canyon forcing humidity, and it added no information that
! was not already in the single level it started from.
a_temp_uclem     = met%tk(iu0:iu1)
a_temproof_uclem = met%tk(iu0:iu1)
a_mixr_uclem     = met%qv(iu0:iu1)
a_mixrroof_uclem = met%qv(iu0:iu1)
a_umag           = a_umag_src
a_umagroof       = a_umag_src

call uclem_prepinterior(1,uclem%fp,uclem%pd,nu)

call uclem_eval(uclem%u_fg, uclem%u_eg, uclem%u_ts, uclem%u_wf, &
                uclem%u_rn, uclem%u_evspsbl, uclem%u_sbl, dels, &
                met%ofsd(iu0:iu1), met%fld(iu0:iu1), air%rho(iu0:iu1), &
                a_temp_uclem, a_temproof_uclem, &
                a_mixr_uclem, a_mixrroof_uclem, met%pmb(iu0:iu1)*100, &
                a_umag, a_umagroof, a_udir, &
                a_rnd_uclem, a_snd_uclem, a_zmin_uclem, a_zroof_uclem, &
                uclem%fp, uclem%fp_intm,uclem%fp_road, uclem%fp_roof, uclem%fp_slab,uclem%fp_wall, &
                uclem%intm, uclem%pd, uclem%rdhyd, uclem%rfhyd, uclem%rfveg, &
                uclem%road, uclem%roof, uclem%room, uclem%slab, uclem%walle, &
                uclem%wallw, uclem%cnveg, uclem%intl, nu, progcalc, 0, &
                uclem%qf,uclem%u_temp9,uclem%u_shf9)

uclem%u_ts = uclem%u_ts + urbtemp

! The urban aerodynamic properties are no longer written back into CABLE's working state.
! Writing them here had no effect on the atmosphere: cable_expl_unpack had already handed
! the UM its values by the time the prognostic call ran, and ruff_resist overwrote
! rough%z0m again at the start of the next timestep.  The urban tile's roughness, drag and
! transfer coefficients are now applied directly to the UM tile arrays in cable_expl_unpack,
! which leaves CABLE's working state, and therefore the natural tiles and the soil column,
! untouched.

deallocate(a_udir)
deallocate(a_umag)
deallocate(a_umag_src)
deallocate(a_umagroof)
deallocate(a_temp_uclem)
deallocate(a_temproof_uclem)
deallocate(a_mixr_uclem)
deallocate(a_mixrroof_uclem)
deallocate(a_rnd_uclem)
deallocate(a_snd_uclem)
deallocate(a_zmin_uclem)
deallocate(a_zroof_uclem)

RETURN
END SUBROUTINE run_uclem_access_coupling

END MODULE uclem_coupling_mod
