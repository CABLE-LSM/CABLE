MODULE crop_def
  !! Definitions of constants and parameters and initialisation routines for
  !! crops.

  use cable_def_types_mod, only: dp=>r_2
  use cable_common_module, only: get_unit, filenames_type
  use casadimension,       only: mplant

  implicit none

  save

  public


  ! variables/parameters
  !! Number of crop types in grid cell.
  integer :: nc = 1
  !! Total number of crop functional types (CFTs). Read from crop parameter
  !! file (filename%crop).
  integer :: ncmax

  !! Crop stage
  integer, parameter :: baresoil = 0
  !! Crop stage
  integer, parameter :: planted = 1
  !! Crop stage
  integer, parameter :: emergent = 2
  !! Crop stage
  integer, parameter :: growing = 3


  ! Parameters related to planting and germination.
  !! Critical soil moisture relative to field capacity for germination.
  real, parameter :: fcrit_germination = 0.8
  !! Number of days checked before possible planting date.
  integer, parameter :: nrdays_before=7
  !! Number of days checked after planting date.
  integer, parameter :: nrdays_after=3

  ! parameters related to carbon allocation etc.
  !! Maximum days since sowing above which germination is assumed to have
  !! failed.
  integer, parameter :: maxdays_ger = 40
  !! conversion dry matter to carbon
  real(dp), parameter :: DMtoC = 0.475_dp
  !! fPHU at flowering.
  real(dp), parameter :: fPHU_flowering = 0.5_dp
  !! Ratio of carbon content in senescent and green leaves.
  real(dp), parameter :: rCsen_Cgr = 0.666_dp
  !! Growth respiration coefficient (see Lokupitiya et al. 2009).
  real(dp), parameter :: Rgcoeff = 0.22_dp
  !! Fraction of leaf C that is mobilised to products at senesc.
  real(dp), parameter :: fCleaf_mobile = 0.5_dp

  ! parameters related to management
  ! irrigation
  !! Soil depth used for calculation of irrigation demand (m).
  real, parameter :: irrig_depth = 0.5
  !! Fraction of available soil moisture below which irrigation is triggered.
  real, parameter :: irrig_trigger = 0.9
  !! Irrigation loss factor (0-1), integrative parameter. Maybe rename to
  !! refill/irrigation efficiency or similar.
  real, parameter :: Firrig_loss = 0.5

  ! types
  type crop_type
    !! Derrived type representing crops.
    !!
    !! Crop stages:
    !! - 0: bare soil, nothing planted/sown
    !! - 1: planted/sown (but not yet germinated)
    !! - 2: emergent
    !! - 3: growing

    !! Crop state.
    integer, dimension (:), pointer :: state
    !! Number of days elapsed in each state.
    integer, dimension (:,:), pointer :: state_nrdays

    ! Switches related to physiology/development
    !! Switch for whether vernalisation occurs
    logical, dimension (:), pointer :: vernalisation
    !! Switch for whether crop responds to photoperiod
    logical, dimension (:), pointer :: photoperiodism

    ! Crop temperature requirements
    !! Crop-specific base temperature.
    real, dimension (:), pointer :: Tbase
    !! Crop-specific max. temperature (no further PHU accumulation above Tmax)
    real, dimension (:), pointer :: Tmax
    !! Planting temperature
    real, dimension (:), pointer :: Tplant

    ! Phenological heat units (PHU)
    !! Phenological heat units (soil) until germination in the absence of
    !! water stress.
    real (dp), dimension (:), pointer :: PHU_germination
    !! fPHU (air) until emergence (use of seed storage) is completed.
    real (dp), dimension (:), pointer :: fPHU_emergence
    !! Phenological heat units (air) until maturity.
    real (dp), dimension (:), pointer :: PHU_maturity
    !! Actual phenological heat units relative to maturity (growth)
    real (dp), dimension (:), pointer :: fPHU

    ! vernalisation requirements
    !! Accrued vernalisation units.
    real (dp), dimension (:), pointer :: VU
    !! Fraction of required vernalisation units.
    real (dp), dimension (:), pointer :: fVU
    !! True if vernalisation effects have been taken into account
    logical, dimension (:), pointer :: vacc

    ! Carbon allocation
    !! Dynamic C allocation? (i.e. allocation depends on soil moisture and
    !! nutrients).
    logical, dimension (:), pointer :: dynamic_allocation
    !! Initial C allocation coefficient to roots at emergence.
    real (dp), dimension (:), pointer :: fCalloc_root_init
    !! Initial C allocation coefficient to leaves at emergence.
    real (dp), dimension (:), pointer :: fCalloc_leaf_init
    !! C allocation coefficient to leaves at first breakpoint.
    real (dp), dimension (:), pointer :: fCalloc_leaf_bpt
    !! fPHU at which C allocation to roots stops.
    real (dp), dimension (:), pointer :: Calloc_root_end
    !! fPHU at first breakpoint (leaves).
    real (dp), dimension (:), pointer :: Calloc_leaf_bpt
    !! fPHU at which C allocation to leaves stops.
    real (dp), dimension (:), pointer :: Calloc_leaf_end
    !! fPHU at which C allocation coefficient to products is 1.
    real (dp), dimension (:), pointer :: Calloc_prod_max
    !! Carbon pool in the stem that is used as mobile reserves.
    real (dp), dimension (:), pointer :: Cstem_mobile

    ! Senescence and remobilisation
    !! Fraction of C in the stem that goes to mobile reserves.
    real (dp), dimension(:), pointer :: fCstem_mobile
    !! Maximum C stored in plant pools.
    real (dp), dimension(:,:), pointer :: Cmax
    !! Maximum C stored in mobile reserves of stem.
    real (dp), dimension(:), pointer :: Cmaxstemmob
    !! Fraction of C mobilised from stem to products.
    real (dp), dimension(:), pointer :: fmobilise
    !! Fraction of the crop canopy that is senescent (0-1).
    real (dp), dimension(:), pointer :: fsenesc
    !! Total (green + brown) LAI.
    real (dp), dimension(:), pointer :: LAItot
    !! Senescent (brown) LAI.
    real (dp), dimension(:), pointer :: LAIsen
    !! Exponent in leaf senescence function.
    real (dp), dimension(:), pointer :: drsen

    ! Planting/Sowing dates (DOY)
    !! Earliest sowing date.
    integer, dimension (:), pointer :: sowing_doymin
    !! Latest sowing date.
    integer, dimension (:), pointer :: sowing_doymax
    !! Actual sowing date.
    integer, dimension (:), pointer :: sowing_doy

    ! Germination requirements
    !! Total germination requirements (0-1)
    real (dp), dimension (:), pointer :: fgermination

    ! Management options
    !! Carbon content in seeds (gC m-2) corresponds to planting density.
    real (dp), dimension (:), pointer :: Cseed
    !! Sowing depth (m).
    real (dp), dimension (:), pointer :: sowing_depth
    !! Fraction of aboveground biomass (stem and leaves) removed at harvest.
    !! Remainder goes to litter pool.
    real (dp), dimension (:), pointer :: Cplant_remove

    !! Soil layer corresponding to snow depth.
    integer, dimension (:), pointer :: sl

    ! Harvest variables
    !! Yield at harvest (gC m-2).
    real (dp), dimension (:), pointer :: yield
    !! Harvest index (Cproduct/total Cplant).
    real (dp), dimension (:), pointer :: harvest_index

    ! Leaf structural variables
    !! Specific leaf area (m2/gDM) at plant maturity.
    real (dp), dimension (:), pointer :: sla_maturity
    !! Extinction cofficient in specific leaf area formula.
    real (dp), dimension (:), pointer :: sla_beta

    ! Management variables
    ! file names of spatial data sets
    !! File containing 'area equipped for irrigation'.
    character (len=200) :: AEI_file
    !! File containing fertilisation data; not yet implemented.
    character (len=200) :: Fertilisation_file


    ! Crop-related variables in other types:
    ! Water added as irrigation directly to soil surface (mm)
    ! canopy%irrig_surface
    ! Water added as irrigation above canopy (mm)
    ! canopy%irrig_sprinkler

  end type crop_type


contains


  subroutine allocate_init_cropvars(crop, filename)
    !! Allocates and initialises crop variables in one routine.

    type (crop_type), intent(inout) :: crop
    type (filenames_type), intent(in) :: filename

    !local
    integer :: j, jcrop   ! loop counter
    integer :: ioerror    ! input error integer
    integer :: vernalisation, photoperiodism ! to be converted to logical
    character (len=25) :: cropnametmp  ! not sure if needed/useful

    open (40, FILE=filename%crop, STATUS='old', ACTION='READ',IOSTAT=ioerror)

    if (ioerror/=0) then
      STOP 'CABLE_log: Cannot open crop parameter file!'
    end if

    read (40,*)
    read (40,*) ncmax ! Number of crop functional types (CFTs).

    ! allocate crop structure
    allocate (crop%state(ncmax),       &
        crop%vernalisation(ncmax),     &
        crop%photoperiodism(ncmax),    &
        crop%Tbase(ncmax),             &
        crop%Tmax(ncmax),              &
        crop%Tplant(ncmax),            &
        crop%PHU_germination(ncmax),   &
        crop%fPHU_emergence(ncmax),    &
        crop%PHU_maturity(ncmax),      &
        crop%fPHU(ncmax),              &
        crop%VU(ncmax),                &
        crop%fVU(ncmax),               &
        crop%vacc(ncmax),              &
        crop%dynamic_allocation(ncmax),&
        crop%fCalloc_root_init(ncmax), &
        crop%fCalloc_leaf_init(ncmax), &
        crop%fCalloc_leaf_bpt(ncmax),  &
        crop%Calloc_root_end(ncmax),   &
        crop%Calloc_leaf_bpt(ncmax),   &
        crop%Calloc_leaf_end(ncmax),   &
        crop%Calloc_prod_max(ncmax),   &
        crop%Cstem_mobile(ncmax),      &
        crop%fCstem_mobile(ncmax),     &
        crop%Cmax(ncmax,mplant),       &
        crop%Cmaxstemmob(ncmax),       &
        crop%fmobilise(ncmax),         &
        crop%fsenesc(ncmax),           &
        crop%LAItot(ncmax),            &
        crop%LAIsen(ncmax),            &
        crop%drsen(ncmax),             &
        crop%sowing_doymin(ncmax),     &
        crop%sowing_doymax(ncmax),     &
        crop%sowing_doy(ncmax),        &
        crop%fgermination(ncmax),      &
        crop%Cseed(ncmax),             &
        crop%sowing_depth(ncmax),      &
        crop%Cplant_remove(ncmax),     &
        crop%sl(ncmax),                &
        crop%yield(ncmax),             &
        crop%harvest_index(ncmax),     &
        crop%sla_maturity(ncmax),      &
        crop%sla_beta(ncmax)           &
        )

    allocate (crop%state_nrdays(ncmax,size(crop%state)))

    ! initialise parameter values given in crop parameter file
    do j=1,ncmax
      read (40,*) jcrop, cropnametmp

      if (jcrop > ncmax) STOP 'jcrop out of range in parameter file'

      crop%vernalisation(jcrop) = .FALSE.
      crop%vernalisation(jcrop) = .FALSE.

      ! Read actual parameter values
      read (40,*) vernalisation, photoperiodism
      read (40,*) crop%Tbase(jcrop), crop%Tmax(jcrop), crop%Tplant(jcrop)
      read (40,*) crop%PHU_germination(jcrop), crop%fPHU_emergence(jcrop)
      read (40,*) crop%PHU_maturity(jcrop)
      read (40,*) crop%fCalloc_root_init(jcrop), &
          crop%fCalloc_leaf_init(jcrop),&
          crop%fCalloc_leaf_bpt(jcrop)
      read (40,*) crop%Calloc_root_end(jcrop), &
          crop%Calloc_leaf_bpt(jcrop), &
          crop%Calloc_leaf_end(jcrop), &
          crop%Calloc_prod_max(jcrop)
      read (40,*) crop%fCstem_mobile(jcrop)
      read (40,*) crop%Cseed(jcrop), crop%sowing_depth(jcrop), &
          crop%Cplant_remove(jcrop)
      read (40,*) crop%sla_maturity(jcrop), crop%sla_beta(jcrop), &
          crop%drsen(jcrop)

      if (vernalisation > 0) then
         crop%vernalisation(jcrop) = .TRUE.
      endif
      if (photoperiodism > 0) then
         crop%photoperiodism(jcrop) = .TRUE.
      endif

    end do ! loop over CFTs

    close (40)

    ! initialise all other variables
    crop%state = 0
    crop%state_nrdays = 0
    crop%fPHU = 0.0_dp
    crop%VU = 0.0_dp
    crop%fVU = 0.0_dp
    crop%vacc = .FALSE.
    crop%dynamic_allocation = .FALSE. ! to be moved to cable.nml, if possible!
    crop%Cstem_mobile = 0.0_dp
    crop%Cmax = 0.0_dp
    crop%Cmaxstemmob = 0.0_dp
    crop%fmobilise = 0.0_dp
    crop%fsenesc = 0.0_dp
    crop%LAItot = 0.0_dp
    crop%LAIsen = 0.0_dp
    crop%sowing_doymin = 272
    crop%sowing_doymax = 335
    crop%sowing_doy = 0
    crop%fgermination = 0.0_dp
    crop%sl = 0
    crop%yield = 0.0_dp
    crop%harvest_index = 0.0_dp
  end subroutine allocate_init_cropvars


  subroutine init_crop_data(crop)
    !! Read spatial datasets needed for crop modeling: irrigation,
    !! fertilisation, crop distribution etc.

    type (crop_type), intent (inout) :: crop

    ! local
    ! Error status returned by nc routines zero=ok, non-zero=error
    integer :: ErrStatus 
    integer :: nmlunit ! Unit number for reading namelist file
    character (len=200) :: AEI_file
    character (len=200) :: Fertilisation_file

    ! define CROP namelist
    NAMELIST /CROP_NML/ AEI_file, Fertilisation_file

    ! Read CROP namelist settings
    CALL get_unit(nmlunit) ! CABLE routine finds spare unit number
    OPEN (nmlunit, FILE="crop.nml", STATUS='OLD', ACTION='READ')
    READ (nmlunit, NML=CROP_NML)
    CLOSE (nmlunit)

    ! Assign namelist settings to corresponding elements in the crop structure
    crop%AEI_file = AEI_file
    crop%Fertilisation_file = Fertilisation_file

    ! Open AEI file
    !ErrStatus = NF90_OPEN(TRIM(crop%AEI_file), NF90_NOWRITE, FID)
  end subroutine init_crop_data


END MODULE crop_def
