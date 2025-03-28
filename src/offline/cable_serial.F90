!==============================================================================
! This source code is part of the
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CSIRO Open Source Software License
! Agreement (variation of the BSD / MIT License).
!
! You may not use this file except in compliance with this License.
! A copy of the License (CSIRO_BSD_MIT_License_v2.0_CABLE.txt) is located
! in each directory containing CABLE code.
!
! ==============================================================================
! Purpose: Offline driver for CABLE
! Contact: Bernard.Pak@csiro.au
!
! History: Since 1.4b, capability to run global offline (ncciy = YEAR),
!          inclusion of call to CASA-CNP (icycle>0)
!          exclusion of call to cbm (icycle>10)
!          soil_snow_type now ssnow (instead of ssoil)
!
!
! ========================================untitled======================================
! Uses:           cable_def_types_mod
!                 cable_IO_vars_module
!                 cable_common_module
!                 cable_input_module
!                 cable_output_module
!                 cable_cbm_module
!                 casadimension
!                 casavariable
!
! CALLs:       open_met_file
!              load_parameters
!              open_output_file
!              get_met_data
!              casa_feedback
!              cbm
!              bgcdriver
!              sumcflux
!              write_output
!              casa_poolout
!              casa_fluxout
!              create_restart
!              close_met_file
!              close_output_file
!              prepareFiles
!
!
! input  file: [SiteName].nc
!              poolcnpIn[SiteName].csv -- for CASA-CNP only
!              gridinfo_CSIRO_1x1.nc
!              def_veg_params.txt
!              def_soil_params.txt -- nearly redundant, can be switched on
!              restart_in.nc -- not strictly required
!
! output file: log_cable.txt
!              out_cable.nc
!              restart_out.nc
!              poolcnpOut.csv -- from CASA-CNP
!==============================================================================

MODULE cable_serial
  !! Offline serial driver for CABLE.
  USE cable_driver_common_mod, ONLY : &
    cable_driver_init,                &
    vegparmnew,                       &
    spinup,                           &
    spincasa,                         &
    CASAONLY,                         &
    l_casacnp,                        &
    l_landuse,                        &
    l_laiFeedbk,                      &
    l_vcmaxFeedbk,                    &
    delsoilM,                         &
    delsoilT,                         &
    delgwM,                           &
    LALLOC,                           &
    prepareFiles,                     &
    prepareFiles_princeton,           &
    LUCdriver,                        &
    compare_consistency_check_values
  USE cable_def_types_mod
  USE cable_IO_vars_module, ONLY: logn,gswpfile,ncciy,leaps,                  &
       fixedCO2,output,check,&
       patch_type,landpt,&
       defaultLAI, sdoy, smoy, syear, timeunits, calendar, &
       NO_CHECK
  USE casa_ncdf_module, ONLY: is_casa_time
  USE cable_common_module,  ONLY: ktau_gl, kend_gl, knode_gl, cable_user,     &
       filename, myhome,            &
       CurYear,    &
       IS_LEAPYEAR, &
       kwidth_gl

! physical constants
USE cable_phys_constants_mod, ONLY : CTFRZ   => TFRZ
USE cable_phys_constants_mod, ONLY : CEMLEAF => EMLEAF
USE cable_phys_constants_mod, ONLY : CEMSOIL => EMSOIL
USE cable_phys_constants_mod, ONLY : CSBOLTZ => SBOLTZ

  USE cable_input_module,   ONLY: open_met_file,load_parameters,              &
       get_met_data,close_met_file,                &
       ncid_rain,       &
       ncid_snow,       &
       ncid_lw,         &
       ncid_sw,         &
       ncid_ps,         &
       ncid_qa,         &
       ncid_ta,         &
       ncid_wd,ncid_mask
  USE cable_output_module,  ONLY: create_restart,open_output_file,            &
       write_output,close_output_file
   USE cable_checks_module, ONLY: constant_check_range
  USE cable_write_module,   ONLY: nullify_write
  USE cable_IO_vars_module, ONLY: timeunits,calendar
   USE cable_cbm_module, ONLY : cbm
  !mpidiff
  USE cable_climate_mod
    
  ! modules related to CASA-CNP
  USE casadimension,        ONLY: icycle
  USE casavariable,         ONLY: casafile, casa_biome, casa_pool, casa_flux,  &
                                !mpidiff
       casa_met, casa_balance, zero_sum_casa, update_sum_casa
  USE phenvariable,         ONLY: phen_variable

  ! vh_js !
  ! modules related to POP
  USE POP_Types,            ONLY: POP_TYPE
  USE POPLUC_Types, ONLY : POPLUC_Type
  USE POPLUC_Module, ONLY:  WRITE_LUC_OUTPUT_NC, WRITE_LUC_OUTPUT_GRID_NC, &
       POP_LUC_CASA_transfer,  WRITE_LUC_RESTART_NC, POPLUC_set_patchfrac
  USE POP_Constants,        ONLY: HEIGHT_BINS, NCOHORT_MAX

  ! PLUME-MIP only
  USE CABLE_PLUME_MIP,      ONLY: PLUME_MIP_TYPE, PLUME_MIP_GET_MET,&
       PLUME_MIP_INIT

  USE CABLE_CRU,            ONLY: CRU_TYPE, CRU_GET_SUBDIURNAL_MET
  USE CABLE_site,           ONLY: site_TYPE, site_GET_CO2_Ndep

  ! LUC_EXPT only
  USE CABLE_LUC_EXPT, ONLY: LUC_EXPT_TYPE, LUC_EXPT_INIT
#ifdef NAG
  USE F90_UNIX
#endif
  USE casa_inout_module
  USE casa_cable
USE cbl_soil_snow_init_special_module
USE landuse_constant, ONLY: mstate,mvmax,mharvw
USE landuse_variable
USE bgcdriver_mod, ONLY : bgcdriver
USE casa_offline_inout_module, ONLY : WRITE_CASA_RESTART_NC, WRITE_CASA_OUTPUT_NC 
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: serialdrv

CONTAINS

SUBROUTINE serialdrv(NRRRR, dels, koffset, kend, GSWP_MID, PLUME, CRU, site)
  !! Offline serial driver.
  INTEGER, INTENT(IN) :: NRRRR !! Number of repeated spin-up cycles
  REAL, INTENT(INOUT) :: dels !! Time step size in seconds
  INTEGER, INTENT(INOUT) :: koffset !! Timestep to start at
  INTEGER, INTENT(INOUT) :: kend !! No. of time steps in run
  INTEGER, ALLOCATABLE, INTENT(INOUT) :: GSWP_MID(:,:) !! NetCDF file IDs for GSWP met forcing
  TYPE(PLUME_MIP_TYPE), INTENT(IN) :: PLUME
  TYPE(CRU_TYPE), INTENT(IN) :: CRU
  TYPE (site_TYPE), INTENT(IN) :: site

  ! timing variables
  INTEGER, PARAMETER ::  kstart = 1   ! start of simulation
  INTEGER, PARAMETER ::  mloop  = 30   ! CASA-CNP PreSpinup loops

  INTEGER        ::                                                           &
       ktau,       &  ! increment equates to timestep, resets if spinning up
       ktau_tot = 0, &  ! NO reset when spinning up, total timesteps by model
       koffset_met = 0, &  !offfset for site met data ('site' only)
       ktauday,    &  ! day counter for CASA-CNP
       idoy,       &  ! day of year (1:365) counter for CASA-CNP
       nyear,      &  ! year counter for CASA-CNP
       casa_it,    &  ! number of calls to CASA-CNP
       YYYY,       &  !
       RYEAR = 0,  &  ! current repeat year
       RRRR,       &  !
       ctime = 0,  &  ! day count for casacnp
       LOY, &         ! days in year
       count_sum_casa ! number of time steps over which casa pools &
                                !and fluxes are aggregated (for output)

  CHARACTER     :: dum*9

  ! CABLE variables
  TYPE (met_type)       :: met     ! met input variables
  TYPE (air_type)       :: air     ! air property variables
  TYPE (canopy_type)    :: canopy  ! vegetation variables
  TYPE (radiation_type) :: rad     ! radiation variables
  TYPE (roughness_type) :: rough   ! roughness varibles
  TYPE (balances_type)  :: bal     ! energy and water balance variables
  TYPE (soil_snow_type) :: ssnow   ! soil and snow variables
  !mpidiff
  TYPE (climate_type)   :: climate     ! climate variables

  ! CABLE parameters
  TYPE (soil_parameter_type) :: soil ! soil parameters
  TYPE (veg_parameter_type)  :: veg  ! vegetation parameters

  TYPE (sum_flux_type)  :: sum_flux ! cumulative flux variables
  TYPE (bgc_pool_type)  :: bgc  ! carbon pool variables

  ! CASA-CNP variables
  TYPE (casa_biome)     :: casabiome
  TYPE (casa_pool)      :: casapool
  TYPE (casa_flux)      :: casaflux
  TYPE (casa_pool)      :: sum_casapool
  TYPE (casa_flux)      :: sum_casaflux
  TYPE (casa_met)       :: casamet
  TYPE (casa_balance)   :: casabal
  TYPE (phen_variable)  :: phen
  ! vh_js !
  TYPE (POP_TYPE)       :: POP
  TYPE(POPLUC_TYPE) :: POPLUC
  TYPE (LUC_EXPT_TYPE) :: LUC_EXPT
  TYPE (landuse_mp)     :: lucmp
  CHARACTER             :: cyear*4
  CHARACTER             :: ncfile*99

  LOGICAL, SAVE           :: &
       spinConv = .FALSE.,         & ! has spinup converged?
       CALL1 = .TRUE.,             &
       SPINon= .TRUE.

  INTEGER :: Metyear, Y, LOYtmp

  ! temporary storage for soil moisture/temp. in spin up mode
  REAL, ALLOCATABLE, DIMENSION(:,:)  :: &
       soilMtemp,                         &
       soilTtemp

  REAL, ALLOCATABLE, DIMENSION(:) :: GWtemp

  ! timing
  REAL:: etime ! Declare the type of etime(), For receiving user and system time, total time
  REAL, allocatable  :: heat_cap_lower_limit(:,:)

  !mpidiff
  INTEGER :: i,x,kk,m,np,ivt

  DOUBLE PRECISION ::                                                                     &
       new_sumbal = 0.0, &
       new_sumfpn = 0.0, &
       new_sumfe = 0.0
  !For consistency w JAC
  REAL,ALLOCATABLE, SAVE :: c1(:,:)
  REAL,ALLOCATABLE, SAVE :: rhoch(:,:)
  REAL,ALLOCATABLE, SAVE :: xk(:,:)

  INTEGER :: nkend=0
  INTEGER :: count_bal = 0

  ! for landuse
  integer     mlon,mlat, mpx
  real(r_2), dimension(:,:,:),   allocatable,  save  :: luc_atransit
  real(r_2), dimension(:,:),     allocatable,  save  :: luc_fharvw
  real(r_2), dimension(:,:,:),   allocatable,  save  :: luc_xluh2cable
  real(r_2), dimension(:),       allocatable,  save  :: arealand
  integer,   dimension(:,:),     allocatable,  save  :: landmask
  integer,   dimension(:),       allocatable,  save  :: cstart,cend,nap
  real(r_2), dimension(:,:,:),   allocatable,  save  :: patchfrac_new

! END header

! INISTUFF


  ! outer loop - spinup loop no. ktau_tot :
  ktau     = 0
  SPINLOOP:DO WHILE ( SPINon )

    NREP: DO RRRR = 1, NRRRR

      YEAR: DO YYYY= CABLE_USER%YearStart,  CABLE_USER%YearEnd

        CurYear = YYYY

        !ccc Set "calendar" for netcdf time attribute and
        ! number of days in the year
        calendar = "noleap"
        LOY = 365
        IF ( leaps ) THEN
          calendar = "standard"
        END IF
        IF ( leaps .AND. IS_LEAPYEAR( YYYY ) ) THEN
          LOY = 366
        END IF

        ! Check for gswp run
        SELECT CASE (TRIM(cable_user%MetType))
        CASE ('gswp', 'prin')
          ncciy = CurYear

          IF (cable_user%MetType == 'gswp') THEN
            CALL prepareFiles(ncciy)
          ELSE ! cable_user%MetType == 'prin'
            CALL prepareFiles_princeton(ncciy) ! MMY
          END IF

          IF ( RRRR .EQ. 1 ) THEN
            CALL open_met_file( dels, koffset, kend, spinup, CTFRZ )
            IF (leaps.AND.is_leapyear(YYYY).AND.kend.EQ.2920) THEN
              STOP 'LEAP YEAR INCOMPATIBILITY WITH INPUT MET !'
            ENDIF
            IF ( NRRRR .GT. 1 ) THEN
              GSWP_MID(1,YYYY) = ncid_rain
              GSWP_MID(2,YYYY) = ncid_snow
              GSWP_MID(3,YYYY) = ncid_lw
              GSWP_MID(4,YYYY) = ncid_sw
              GSWP_MID(5,YYYY) = ncid_ps
              GSWP_MID(6,YYYY) = ncid_qa
              GSWP_MID(7,YYYY) = ncid_ta
              GSWP_MID(8,YYYY) = ncid_wd
            ENDIF
          ELSE
            ncid_rain = GSWP_MID(1,YYYY)
            ncid_snow = GSWP_MID(2,YYYY)
            ncid_lw   = GSWP_MID(3,YYYY)
            ncid_sw   = GSWP_MID(4,YYYY)
            ncid_ps   = GSWP_MID(5,YYYY)
            ncid_qa   = GSWP_MID(6,YYYY)
            ncid_ta   = GSWP_MID(7,YYYY)
            ncid_wd   = GSWP_MID(8,YYYY)
            kend      = ktauday * LOY
          ENDIF

        CASE ('plum')
          ! PLUME experiment setup using WATCH
          IF ( .NOT. PLUME%LeapYears ) LOY = 365
          kend = NINT(24.0*3600.0/dels) * LOY

        CASE ('cru')
          ! TRENDY experiment using CRU-NCEP
          LOY = 365
          kend = NINT(24.0*3600.0/dels) * LOY

        CASE ('gswp3')
          ncciy = CurYear
          CALL open_met_file( dels, koffset, kend, spinup, CTFRZ )

        CASE ('site')
          ! site experiment eg AmazonFace (spinup or  transient run type)
          kend = NINT(24.0*3600.0/dels) * LOY
          ! get koffset to add to time-step of sitemet
          SELECT CASE (TRIM(site%RunType))
          CASE ('historical')
            MetYear = CurYear
          CASE ('spinup', 'transient')
            ! setting met year so we end the spin-up at the end of the site data-years.
            MetYear = site%spinstartyear + &
                MOD(CurYear- &
                (site%spinstartyear-(site%spinendyear-site%spinstartyear +1)*100), &
                (site%spinendyear-site%spinstartyear +1))
          END SELECT
          WRITE(*,*) 'MetYear: ', MetYear
          WRITE(*,*) 'Simulation Year: ', CurYear
          koffset_met = 0
          IF (MetYear .GT. site%spinstartyear) THEN
            DO Y = site%spinstartyear, MetYear-1
              LOYtmp = 365
              IF (IS_LEAPYEAR(Y)) LOYtmp = 366
              koffset_met = koffset_met + INT( REAL(LOYtmp) * 86400./REAL(dels) )
            ENDDO
          ENDIF

        END SELECT

        ! somethings (e.g. CASA-CNP) only need to be done once per day
        ktauday=INT(24.0*3600.0/dels)

        ! Checks where parameters and initialisations should be loaded from.
        ! If they can be found in either the met file or restart file, they will
        ! load from there, with the met file taking precedence. Otherwise, they'll
        ! be chosen from a coarse global grid of veg and soil types, based on
        ! the lat/lon coordinates. Allocation of CABLE's main variables also here.
        IF ( CALL1 ) THEN

          IF (cable_user%POPLUC) THEN
            CALL LUC_EXPT_INIT (LUC_EXPT)
          ENDIF
          ! vh_js !
          CALL load_parameters( met, air, ssnow, veg,climate,bgc,           &
               soil, canopy, rough, rad, sum_flux,                   &
               bal, logn, vegparmnew, casabiome, casapool,           &
               casaflux, sum_casapool, sum_casaflux, &
               casamet, casabal, phen, POP, spinup,        &
               CEMSOIL, CTFRZ, LUC_EXPT, POPLUC )

          IF (check%ranges /= NO_CHECK) THEN
            WRITE (*, *) "Checking parameter ranges"
            CALL constant_check_range(soil, veg, 0, met)
          END IF

          IF ( CABLE_USER%POPLUC .AND. TRIM(CABLE_USER%POPLUC_RunType) .EQ. 'static') &
               CABLE_USER%POPLUC= .FALSE.

          ! Open output file:
          IF (.NOT.CASAONLY) THEN
            IF ( TRIM(filename%out) .EQ. '' ) THEN
              IF ( CABLE_USER%YEARSTART .GT. 0 ) THEN
                WRITE( dum, FMT="(I4,'_',I4)")CABLE_USER%YEARSTART, &
                     CABLE_USER%YEAREND
                filename%out = TRIM(filename%path)//'/'//&
                    TRIM(cable_user%RunIden)//'_'//&
                    TRIM(dum)//'_cable_out.nc'
              ELSE
                filename%out = TRIM(filename%path)//'/'//&
                    TRIM(cable_user%RunIden)//'_cable_out.nc'
              ENDIF
            ENDIF
            CALL nullify_write() ! nullify pointers
            CALL open_output_file( dels, soil, veg, bgc, rough, met)
          ENDIF

          ssnow%otss_0 = ssnow%tgg(:,1)
          ssnow%otss = ssnow%tgg(:,1)
          ssnow%tss = ssnow%tgg(:,1)
          canopy%fes_cor = 0.
          canopy%fhs_cor = 0.
          met%ofsd = 0.1


          CALL zero_sum_casa(sum_casapool, sum_casaflux)
          count_sum_casa = 0

          IF (cable_user%call_climate) CALL climate_init ( climate, mp, ktauday )
          IF (cable_user%call_climate .AND.(.NOT.cable_user%climate_fromzero)) &
              CALL READ_CLIMATE_RESTART_NC (climate, ktauday)

          spinConv = .FALSE. ! initialise spinup convergence variable
          IF (.NOT.spinup)  spinConv=.TRUE.
          IF( icycle>0 .AND. spincasa) THEN
            PRINT *, 'EXT spincasacnp enabled with mloop= ', mloop
            CALL spincasacnp(dels,kstart,kend,mloop,veg,soil,casabiome,casapool, &
                  casaflux,casamet,casabal,phen,POP,climate,LALLOC)
            SPINon = .FALSE.
            SPINconv = .FALSE.

          ELSEIF ( casaonly .AND. (.NOT. spincasa) ) THEN !.AND. cable_user%popluc) THEN

            CALL CASAONLY_LUC(dels,kstart,kend,veg,soil,casabiome,casapool, &
                 casaflux,casamet,casabal,phen,POP,climate,LALLOC, LUC_EXPT, POPLUC, &
                 sum_casapool, sum_casaflux)
            SPINon = .FALSE.
            SPINconv = .FALSE.
            ktau = kend

          ENDIF



        ENDIF ! CALL 1

        ! globally (WRT code) accessible kend through USE cable_common_module
        kwidth_gl = INT(dels)
        kend_gl  = kend
        knode_gl = 0

        IF (casaonly) THEN
          EXIT
        ENDIF

        if( .NOT. allocated(heat_cap_lower_limit) ) then
          allocate( heat_cap_lower_limit(mp,ms) )
          heat_cap_lower_limit = 0.01
        end if

        call spec_init_soil_snow(dels, soil, ssnow, canopy, met, bal, veg, heat_cap_lower_limit)

        ! time step loop over ktau
        DO ktau=kstart, kend

          WRITE(logn,*) 'Progress -',REAL(ktau)/REAL(kend)*100.0

          ! increment total timstep counter
          ktau_tot = ktau_tot + 1

          ! globally (WRT code) accessible kend through USE cable_common_module
          ktau_gl = ktau_tot

          idoy =INT( MOD(REAL(CEILING(REAL((ktau+koffset)/ktauday))),REAL(LOY)))
          IF ( idoy .EQ. 0 ) idoy = LOY

          ! needed for CASA-CNP
          nyear     =INT((kend+koffset)/(LOY*ktauday))

          ! Get met data and LAI, set time variables.
          ! Rainfall input may be augmented for spinup purposes:
          IF (CALL1 .AND. CASAONLY) THEN
            SELECT CASE (TRIM(cable_user%MetType))
            CASE ('plum')
              CALL PLUME_MIP_GET_MET(PLUME, MET, YYYY, ktau, kend, &
                   (YYYY.EQ.CABLE_USER%YearEnd .AND. ktau.EQ.kend))

            CASE ('cru')
              CALL CRU_GET_SUBDIURNAL_MET(CRU, met, &
                   YYYY, ktau, kend, &
                   YYYY.EQ.CABLE_USER%YearEnd)
            END SELECT
          END IF

          SELECT CASE (TRIM(cable_user%MetType))
          CASE ('plum')
            IF (.NOT. CASAONLY) THEN
              CALL PLUME_MIP_GET_MET(PLUME, MET, YYYY, ktau, kend, &
                   (YYYY.EQ.CABLE_USER%YearEnd .AND. ktau.EQ.kend))
            END IF

          CASE ('cru')
            IF (.NOT. CASAONLY) THEN
              CALL CRU_GET_SUBDIURNAL_MET(CRU, met, &
                   YYYY, ktau, kend, &
                   YYYY.EQ.CABLE_USER%YearEnd)
            END IF

          CASE ('site')
            CALL get_met_data( spinup, spinConv, met, soil,            &
                 rad, veg, kend, dels, CTFRZ, ktau+koffset_met,             &
                 kstart+koffset_met )

            CALL site_get_CO2_Ndep(site)

            ! Two options: (i) if we have sub-annual varying CO2, then
            ! these data should be put into the met file and in site.nml
            ! CO2 should be set to -9999; or (ii) if we only have annual
            ! CO2 numbers then these should be read from the site csv file
            WHERE (met%ca .EQ. fixedCO2/1000000.0)
              met%ca = site%CO2 / 1.e+6
            END WHERE

            met%Ndep = site%Ndep  *1000./10000./365. ! kg ha-1 y-1 > g m-2 d-1
            met%Pdep = site%Pdep  *1000./10000./365. ! kg ha-1 y-1 > g m-2 d-1
            met%fsd = MAX(met%fsd,0.0)

          CASE DEFAULT
            CALL get_met_data( spinup, spinConv, met, soil,            &
              rad, veg, kend, dels, CTFRZ, ktau+koffset,                 &
              kstart+koffset )

          END SELECT

          IF (TRIM(cable_user%MetType).EQ.'' ) THEN
            CurYear = met%year(1)
            IF ( leaps .AND. IS_LEAPYEAR( CurYear ) ) THEN
              LOY = 366
            ELSE
              LOY = 365
            ENDIF
          ENDIF
          met%ofsd = met%fsd(:,1) + met%fsd(:,2)
          canopy%oldcansto=canopy%cansto
          ! Zero out lai where there is no vegetation acc. to veg. index
          WHERE ( veg%iveg(:) .GE. 14 ) veg%vlai = 0.

          ! At first time step of year, set tile area according to updated LU areas
          ! and zero casa fluxes
          IF (ktau == 1) THEN
            IF (icycle>1) CALL casa_cnpflux(casaflux,casapool,casabal,.TRUE.)
            IF ( CABLE_USER%POPLUC) CALL POPLUC_set_patchfrac(POPLUC,LUC_EXPT)
          ENDIF

          IF ( .NOT. CASAONLY ) THEN

            ! Feedback prognostic vcmax and daily LAI from casaCNP to CABLE
            IF (l_vcmaxFeedbk) CALL casa_feedback( ktau, veg, casabiome,    &
                casapool, casamet )

            IF (l_laiFeedbk.AND.icycle>0) veg%vlai(:) = casamet%glai(:)

            IF (.NOT. allocated(c1)) ALLOCATE( c1(mp,nrb), rhoch(mp,nrb), xk(mp,nrb) )
            ! Call land surface scheme for this timestep, all grid points:
            CALL cbm( ktau, dels, air, bgc, canopy, met, bal,                             &
                 rad, rough, soil, ssnow, sum_flux, veg, climate, xk, c1, rhoch )

            IF (cable_user%CALL_climate) &
              CALL cable_climate(ktau_tot,kstart,kend,ktauday,idoy,LOY,met, &
                   climate, canopy, air, rad, dels, mp)


            ssnow%smelt = ssnow%smelt*dels
            ssnow%rnof1 = ssnow%rnof1*dels
            ssnow%rnof2 = ssnow%rnof2*dels
            ssnow%runoff = ssnow%runoff*dels





          ELSE IF ( IS_CASA_TIME("dread", yyyy, ktau, kstart, &
                koffset, kend, ktauday, logn) ) THEN                ! CLN READ FROM FILE INSTEAD !
            WRITE(CYEAR,FMT="(I4)")CurYear + INT((ktau-kstart+koffset)/(LOY*ktauday))
            ncfile  = TRIM(casafile%c2cdumppath)//'c2c_'//CYEAR//'_dump.nc'
            casa_it = NINT( REAL(ktau / ktauday) )

            CALL read_casa_dump( ncfile, casamet, casaflux,phen, climate, casa_it, kend, .FALSE. )
          ENDIF

          !jhan this is insufficient testing. condition for
          !spinup=.false. & we want CASA_dump.nc (spinConv=.true.)
          IF(icycle >0 .OR.  CABLE_USER%CASA_DUMP_WRITE ) THEN
            ! vh_js !

            CALL bgcdriver( ktau, kstart, kend, dels, met,                &
                 ssnow, canopy, veg, soil, climate, casabiome,                     &
                 casapool, casaflux, casamet, casabal,                    &
                 phen, pop, spinConv, spinup, ktauday, idoy, loy,         &
                 CABLE_USER%CASA_DUMP_READ, CABLE_USER%CASA_DUMP_WRITE,   &
                 LALLOC )

            IF(MOD((ktau-kstart+1),ktauday)==0) THEN

              !mpidiff
              ! update time-aggregates of casa pools and fluxes
              CALL update_sum_casa(sum_casapool, sum_casaflux, casapool, casaflux, &
                   & .TRUE. , .FALSE., 1)
              count_sum_casa = count_sum_casa + 1
            ENDIF


            IF( ((MOD((ktau-kstart+1),ktauday)==0) .AND.  &
                  MOD((ktau-kstart+1)/ktauday,LOY)==0) )THEN ! end of year
              IF (CABLE_USER%POPLUC) THEN
                ! Dynamic LUC
                CALL LUCdriver( casabiome,casapool,casaflux,POP,     &
                     LUC_EXPT, POPLUC, veg )
              ENDIF

              ! one annual time-step of POP
              CALL POPdriver(casaflux,casabal,veg, POP)

              IF (CABLE_USER%POPLUC) THEN
                ! Dynamic LUC: update casa pools according to LUC transitions
                CALL POP_LUC_CASA_transfer(POPLUC,POP,LUC_EXPT,casapool,casabal,casaflux,ktauday)
                ! Dynamic LUC: write output
                CALL WRITE_LUC_OUTPUT_NC( POPLUC, YYYY, ( YYYY.EQ.cable_user%YearEnd ))

              ENDIF
            ENDIF

          ENDIF
          ! WRITE CASA OUTPUT
          IF(icycle >0) THEN


            IF ( IS_CASA_TIME("write", yyyy, ktau, kstart, &
                 koffset, kend, ktauday, logn) ) THEN
              ctime = ctime +1
              !mpidiff
              CALL update_sum_casa(sum_casapool, sum_casaflux, casapool, casaflux, &
                   .FALSE. , .TRUE. , count_sum_casa)
              CALL WRITE_CASA_OUTPUT_NC (veg, casamet, sum_casapool, casabal, sum_casaflux, &
                   CASAONLY, ctime, ( ktau.EQ.kend .AND. YYYY .EQ.               &
                   cable_user%YearEnd.AND. RRRR .EQ.NRRRR ) )
              !mpidiff
              count_sum_casa = 0
              CALL zero_sum_casa(sum_casapool, sum_casaflux)
            ENDIF


            IF (((.NOT.spinup).OR.(spinup.AND.spinConv)).AND. &
                 MOD((ktau-kstart+1),ktauday)==0) THEN
              IF ( CABLE_USER%CASA_DUMP_WRITE )  THEN

                SELECT CASE (TRIM(cable_user%MetType))
                CASE ('', 'site')
                  WRITE(CYEAR,FMT="(I4)") CurYear
                CASE DEFAULT
                  !CLN CHECK FOR LEAP YEAR
                  WRITE(CYEAR,FMT="(I4)") CurYear + INT((ktau-kstart)/(LOY*ktauday))
                END SELECT
                ncfile = TRIM(casafile%c2cdumppath)//'c2c_'//CYEAR//'_dump.nc'

                IF (TRIM(cable_user%MetType).EQ.'' ) THEN
                  !jhan:assuming doy for mp=1 is same as ....
                  CALL write_casa_dump( ncfile, casamet , casaflux, phen, climate,&
                       INT(met%doy(1)), LOY )
                ELSE
                  CALL write_casa_dump( ncfile, casamet , casaflux, &
                       phen, climate, idoy, kend/ktauday )
                ENDIF

              ENDIF
            ENDIF

          ENDIF


          IF ( .NOT. CASAONLY ) THEN
            ! sumcflux is pulled out of subroutine cbm
            ! so that casaCNP can be called before adding the fluxes
            ! (Feb 2008, YP)
            CALL sumcflux( ktau, kstart, kend, dels, bgc,              &
                 canopy, soil, ssnow, sum_flux, veg,                   &
                 met, casaflux, l_vcmaxFeedbk )


          ENDIF

          ! Write timestep's output to file if either: we're not spinning up
          ! or we're spinning up and the spinup has converged:

          IF ( (.NOT. CASAONLY) .AND. spinConv ) THEN
            !mpidiff
            SELECT CASE (TRIM(cable_user%MetType))
            CASE ('plum', 'cru', 'bios', 'gswp', 'gswp3', 'site')
              CALL write_output( dels, ktau_tot, met, canopy, casaflux, casapool, casamet, &
                   ssnow, rad, bal, air, soil, veg, CSBOLTZ, CEMLEAF, CEMSOIL )
            CASE DEFAULT
              CALL write_output( dels, ktau, met, canopy, casaflux, casapool, casamet, &
                   ssnow, rad, bal, air, soil, veg, CSBOLTZ, CEMLEAF, CEMSOIL )
            END SELECT
          ENDIF


          ! Check triggered by cable_user%consistency_check = .TRUE. in cable.nml
          IF(cable_user%consistency_check) THEN

            count_bal = count_bal +1;
            new_sumbal = new_sumbal + SUM(bal%wbal)/mp +  SUM(bal%ebal)/mp
            new_sumfpn = new_sumfpn + SUM(canopy%fpn)/mp
            new_sumfe = new_sumfe + SUM(canopy%fe)/mp
            IF (ktau == kend) PRINT*
            IF (ktau == kend) PRINT*, "time-space-averaged energy & water balances"
            IF (ktau == kend) PRINT*,"Ebal_tot[Wm-2], Wbal_tot[mm per timestep]", &
                 SUM(bal%ebal_tot)/mp/count_bal, SUM(bal%wbal_tot)/mp/count_bal
            IF (ktau == kend) PRINT*, "time-space-averaged latent heat and net photosynthesis"
            IF (ktau == kend) PRINT*, "sum_fe[Wm-2], sum_fpn[umol/m2/s]",  &
                 new_sumfe/count_bal, new_sumfpn/count_bal
            IF (ktau == kend) WRITE(logn,*)
            IF (ktau == kend) WRITE(logn,*), "time-space-averaged energy & water balances"
            IF (ktau == kend) WRITE(logn,*),"Ebal_tot[Wm-2], Wbal_tot[mm per timestep]", &
                 SUM(bal%ebal_tot)/mp/count_bal, SUM(bal%wbal_tot)/mp/count_bal
            IF (ktau == kend) WRITE(logn,*), "time-space-averaged latent heat and net photosynthesis"
            IF (ktau == kend) WRITE(logn,*), "sum_fe[Wm-2], sum_fpn[umol/m2/s]",  &
                 new_sumfe/count_bal, new_sumfpn/count_bal


            ! vh ! commented code below detects Nans in evaporation flux and stops if there are any.
!$          do kk=1,mp
!$            if( canopy%fe(kk).NE.( canopy%fe(kk))) THEN
!$              write(*,*) 'fe nan', kk, ktau,met%qv(kk), met%precip(kk),met%precip_sn(kk), &
!$                   met%fld(kk), met%fsd(kk,:), met%tk(kk), met%ua(kk), ssnow%potev(kk), met%pmb(kk), &
!$                   canopy%ga(kk), ssnow%tgg(kk,:), canopy%fwsoil(kk)
!$
!$
!$              stop
!$            endif
!$            if ( casaflux%cnpp(kk).NE. casaflux%cnpp(kk)) then
!$              write(*,*) 'npp nan', kk, ktau,  casaflux%cnpp(kk)
!$              stop
!$
!$            endif
!$
!$
!$            !if (canopy%fwsoil(kk).eq.0.0) then
!$            !   write(*,*) 'zero fwsoil', ktau, canopy%fpn(kk)
!$            !endif
!$
!$
!$          enddo

            IF( ktau == kend ) THEN
              nkend = nkend+1
              PRINT *, "NB. Offline-serial runs spinup cycles:", nkend
              CALL compare_consistency_check_values(new_sumbal)
            ENDIF

          ENDIF

          CALL1 = .FALSE.

        END DO ! END Do loop over timestep ktau

        CALL1 = .FALSE.

        !jhan this is insufficient testing. condition for
        !spinup=.false. & we want CASA_dump.nc (spinConv=.true.)
        ! see if spinup (if conducting one) has converged:
        !IF(spinup.AND..NOT.spinConv) THEN
        IF(spinup.AND.(.NOT.spinConv).AND.(.NOT.CASAONLY) ) THEN

          ! Write to screen and log file:
          WRITE(*,'(A18,I3,A24)') ' Spinning up: run ',INT(ktau_tot/kend),&
               ' of data set complete...'
          WRITE(logn,'(A18,I3,A24)') ' Spinning up: run ',                &
               INT(ktau_tot/kend), ' of data set complete...'


          ! IF not 1st run through whole dataset:
!$        IF( MOD( ktau_tot, kend ) .EQ. 0 .AND. ktau_Tot .GT. kend .AND. &
!$             YYYY.EQ. CABLE_USER%YearEnd .OR. ( NRRRR .GT. 1 .AND. &
!$             RRRR.EQ. NRRRR) ) THEN

          IF( MOD( ktau_tot, kend ) .EQ. 0 .AND. ktau_Tot .GT. kend .AND. &
               YYYY.EQ. CABLE_USER%YearEnd ) THEN

            ! evaluate spinup
            IF( ANY( ABS(ssnow%wb-soilMtemp)>delsoilM).OR.               &
                 ANY( ABS(ssnow%tgg-soilTtemp)>delsoilT) .OR. &
                 MAXVAL(ABS(ssnow%GWwb-GWtemp),dim=1) > delgwM) THEN

              ! No complete convergence yet
              PRINT *, 'ssnow%wb : ', ssnow%wb
              PRINT *, 'soilMtemp: ', soilMtemp
              PRINT *, 'ssnow%tgg: ', ssnow%tgg
              PRINT *, 'soilTtemp: ', soilTtemp
              IF (cable_user%gw_model) THEN
                PRINT *, 'ssnow%GWwb : ', ssnow%GWwb
                PRINT *, 'GWtemp: ', GWtemp
              ENDIF

            ELSE ! spinup has converged

              spinConv = .TRUE.
              ! Write to screen and log file:
              WRITE(*,'(A33)') ' Spinup has converged - final run'
              WRITE(logn,'(A52)')                                       &
                   ' Spinup has converged - final run - writing all data'
              WRITE(logn,'(A37,F8.5,A28)')                              &
                   ' Criteria: Change in soil moisture < ',             &
                   delsoilM, ' in any layer over whole run'
              WRITE(logn,'(A40,F8.5,A28)' )                             &
                   '           Change in soil temperature < ',           &
                   delsoilT, ' in any layer over whole run'
            END IF

          ELSE ! allocate variables for storage

            IF (.NOT.ALLOCATED(soilMtemp)) ALLOCATE(  soilMtemp(mp,ms) )
            IF (.NOT.ALLOCATED(soilTtemp)) ALLOCATE(  soilTtemp(mp,ms) )
            IF (.NOT.ALLOCATED(GWtemp)   ) ALLOCATE(  GWtemp(mp)  )

          END IF

          IF (cable_user%max_spins .GT. 0) THEN
            IF ((ktau_tot/kend .GE. cable_user%max_spins)) THEN
              spinConv = .TRUE.
              WRITE(logn,*) 'Past ',cable_user%max_spins,' spin cycles, running anyways'
            END IF
          END IF

          ! store soil moisture and temperature
          IF ( YYYY.EQ. CABLE_USER%YearEnd ) THEN
            soilTtemp = ssnow%tgg
            soilMtemp = REAL(ssnow%wb)
            GWtemp = ssnow%GWwb
          ENDIF

        ELSE

          ! if not spinning up, or spin up has converged, exit:
          IF ( SpinOn ) THEN
            PRINT*,"setting SPINON -> FALSE", YYYY, RRRR
            SPINon = .FALSE.
          END IF

        END IF


        IF((.NOT.(spinup.OR.casaonly)).OR.(spinup.AND.spinConv)) THEN
          IF (icycle > 0) THEN
            CALL casa_fluxout( nyear, veg, soil, casabal, casamet)
          END IF

        ENDIF

        IF ( .NOT. (spinup.OR.casaonly) .OR. spinconv ) THEN
          IF ( NRRRR .GT. 1 ) THEN
            RYEAR = YYYY + ( CABLE_USER%YearEnd - CABLE_USER%YearStart + 1 ) &
                 * ( RRRR - 1 )
          ELSE
            RYEAR = YYYY
          END IF
          IF ( cable_user%CALL_POP.AND.POP%np.GT.0 ) THEN

            IF (TRIM(cable_user%POP_out).EQ.'epi') THEN
              CALL POP_IO( pop, casamet, RYEAR, 'WRITE_EPI', &
                   (YYYY.EQ.CABLE_USER%YearEnd .AND. RRRR.EQ.NRRRR) )

            ENDIF
          ENDIF

        ENDIF
        ! Close met data input file:

        IF ( TRIM(cable_user%MetType) .EQ. "gswp" .AND. &
             RRRR .EQ. NRRRR ) THEN
          CALL close_met_file
          IF ( YYYY .EQ. CABLE_USER%YearEnd .AND. &
               NRRRR .GT. 1 ) DEALLOCATE ( GSWP_MID )
        ENDIF
        IF (TRIM(cable_user%MetType) == "gswp3") CALL close_met_file

        IF ((icycle.GT.0).AND.(.NOT.casaonly)) THEN
          ! re-initalise annual flux sums
          casabal%FCgppyear=0.0
          casabal%FCrpyear=0.0
          casabal%FCnppyear=0
          casabal%FCrsyear=0.0
          casabal%FCneeyear=0.0
        ENDIF
        CALL CPU_TIME(etime)
        PRINT *, 'Finished. ', etime, ' seconds needed for year'

      END DO YEAR


    END DO NREP

  END DO SPINLOOP

  l_landuse=.false.

  IF ( SpinConv .AND. .NOT. CASAONLY ) THEN
    ! Close output file and deallocate main variables:
    CALL close_output_file( bal, air, bgc, canopy, met,                      &
         rad, rough, soil, ssnow,                                            &
         sum_flux, veg )
  ENDIF

  IF ( cable_user%CALL_POP.AND.POP%np.GT.0 ) THEN
    !mpidiff
    IF ( CASAONLY .OR. cable_user%pop_fromzero &
         .OR.TRIM(cable_user%POP_out).EQ.'ini' ) THEN
      CALL POP_IO( pop, casamet, RYEAR+1, 'WRITE_INI', .TRUE.)
    ELSE
      CALL POP_IO( pop, casamet, RYEAR+1, 'WRITE_RST', .TRUE.)
    ENDIF
  ENDIF



  IF (icycle > 0.and. .not.l_landuse) THEN

    !CALL casa_poolout( ktau, veg, soil, casabiome,              &
    ! casapool, casaflux, casamet, casabal, phen )
    CALL write_casa_restart_nc ( casamet, casapool,casaflux,phen, CASAONLY )

  END IF

  IF (l_landuse .AND. .NOT. CASAONLY ) THEN
    mlon = maxval(landpt(1:mp)%ilon)
    mlat = maxval(landpt(1:mp)%ilat)
    print *, 'before landuse: mlon mlat ', mlon,mlat
    allocate(luc_atransit(mland,mvmax,mvmax))
    allocate(luc_fharvw(mland,mharvw))
    allocate(luc_xluh2cable(mland,mvmax,mstate))
    allocate(landmask(mlon,mlat))
    allocate(arealand(mland))
    allocate(patchfrac_new(mlon,mlat,mvmax))
    allocate(cstart(mland),cend(mland),nap(mland))

    do m=1,mland
      cstart(m) = landpt(m)%cstart
      cend(m)   = landpt(m)%cend
      nap(m)    = landpt(m)%nap
    enddo

    call landuse_data(mlon,mlat,landmask,arealand,luc_atransit,luc_fharvw,luc_xluh2cable)
    call  landuse_driver(mlon,mlat,landmask,arealand,ssnow,soil,veg,bal,canopy,  &
                         phen,casapool,casabal,casamet,casabiome,casaflux,bgc,rad, &
                         cstart,cend,nap,lucmp)

    do m=1,mland
      do np=cstart(m),cend(m)
        ivt=lucmp%iveg(np)
        if(ivt<1.or.ivt>mvmax) then
          print *, 'landuse: error in vegtype',m,np,ivt
          stop
        endif
        patchfrac_new(landpt(m)%ilon,landpt(m)%ilat,ivt) = lucmp%patchfrac(np)
      enddo
    enddo

    call create_new_gridinfo(filename%type,filename%gridnew,mlon,mlat,landmask,patchfrac_new)

    print *, 'writing casapools: land use'
    call WRITE_LANDUSE_CASA_RESTART_NC(cend(mland), lucmp, CASAONLY )

    print *, 'writing cable restart: land use'
    call create_landuse_cable_restart(logn, dels, ktau, soil, cend(mland),lucmp,cstart,cend,nap, met)

    print *, 'deallocating'
    call landuse_deallocate_mp(cend(mland),ms,msn,nrb,mplant,mlitter,msoil,mwood,lucmp)

  ENDIF

  IF (cable_user%POPLUC .AND. .NOT. CASAONLY ) THEN
    CALL WRITE_LUC_RESTART_NC ( POPLUC, YYYY )
  ENDIF

  IF ( .NOT. CASAONLY.and. .not. l_landuse ) THEN
    ! Write restart file if requested:
    IF(output%restart)                                           &
      CALL create_restart( logn, dels, ktau, soil, veg, ssnow,  &
           canopy, rough, rad, bgc, bal, met )
    !mpidiff
    IF (cable_user%CALL_climate) &
      CALL WRITE_CLIMATE_RESTART_NC ( climate, ktauday )

    !--- LN ------------------------------------------[
  ENDIF



  IF ( TRIM(cable_user%MetType) .NE. "gswp" .AND. &
       TRIM(cable_user%MetType) .NE. "gswp3" .AND. &
       TRIM(cable_user%MetType) .NE. "plum" .AND. &
       TRIM(cable_user%MetType) .NE. "cru" ) CALL close_met_file

  !WRITE(logn,*) bal%wbal_tot, bal%ebal_tot, bal%ebal_tot_cncheck
  CALL CPU_TIME(etime)
  WRITE(logn,*) 'Finished. ', etime, ' seconds needed for ', kend,' hours'
  ! Close log file
  CLOSE(logn)
  CALL CPU_TIME(etime)
  PRINT *, 'Finished. ', etime, ' seconds needed for ', kend,' hours'

END SUBROUTINE serialdrv

END MODULE cable_serial
