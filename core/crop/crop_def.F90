MODULE crop_def

  use cable_def_types_mod, only: dp => r_2
  use cable_common_module, only: filenames_type
  use casadimension,       only: mplant
  
  implicit none

  save

  public

  
  ! variables/parameters
  integer :: nc=1     ! number of crop types in grid cell
  integer :: ncmax    ! total number of crop functional types (CFTs)
                      ! read from crop parameter file (filename%crop)

  ! crop stages
  integer, parameter :: baresoil=0
  integer, parameter :: sown=1
  integer, parameter :: emergent=2
  integer, parameter :: growing=3

  
  integer, parameter :: maxdays_ger=40 ! maximum days since sowing above which germination
                                       ! is assumed to have failed
  real(dp), parameter :: DMtoC=0.475_dp         ! conversion dry matter to carbon
  real(dp), parameter :: fPHU_flowering=0.5_dp  ! fPHU at flowering
  real(dp), parameter :: rCsen_Cgr=0.666_dp     ! ratio of carbon content in senescent and green leaves 
  real(dp), parameter :: Rgcoeff=0.22_dp        ! growth respiration coefficient (see Lokupitiya et al. 2009)
  real(dp), parameter :: fCleaf_mobile=0.5_dp   ! fraction of leaf C that is mobilised to products at senesc.
  
  ! types 
  type crop_type

    ! crop stages:
    ! 0: bare soil, not sown
    ! 1: sown (but not yet germinated)
    ! 2: emergent
    ! 3: growing
    integer, dimension(:),  pointer :: state

    ! Switches related to physiology/development
    logical, dimension(:), pointer :: vernalisation   ! does vernalisation occur?
    logical, dimension(:), pointer :: photoperiodism  ! crop responds to photoperiod?
    
    ! Crop temperature requirements
    real, dimension(:), pointer :: Tbase     ! crop-specific base temperature
    real, dimension(:), pointer :: Tmax      ! crop-specific max. temperature
                                             ! (no further PHU accumulation above Tmax)

    ! Phenological heat units (PHU)
    real(dp), dimension(:), pointer :: PHU_germination  ! PHU (soil) until germination
                                                        ! in the absence of water stress
    real(dp), dimension(:), pointer :: fPHU_emergence   ! fPHU (air) until emergence (use of seed storage)
                                                        ! is completed
    real(dp), dimension(:), pointer :: PHU_maturity     ! PHU (air) until maturity
    real(dp), dimension(:), pointer :: fPHU             ! actual PHU relative to maturity (growth)

    ! vernalisation requirements
    real(dp), dimension(:), pointer :: VU     ! accrued vernalisation units
    real(dp), dimension(:), pointer :: fVU    ! fraction of required vernalisation units
    logical,  dimension(:), pointer :: vacc   ! true if vernalisation effects have been taken into account

    ! Carbon allocation
    real(dp), dimension(:), pointer :: fCalloc_root_init  ! initial C allocation coefficient to roots at emergence
    real(dp), dimension(:), pointer :: fCalloc_leaf_init  ! initial C allocation coefficient to leaves at emergence
    real(dp), dimension(:), pointer :: fCalloc_leaf_bpt   ! C allocation coefficient to leaves at first breakpoint
    real(dp), dimension(:), pointer :: Calloc_root_end    ! fPHU at which C allocation to roots stops
    real(dp), dimension(:), pointer :: Calloc_leaf_bpt    ! fPHU at first breakpoint (leaves)
    real(dp), dimension(:), pointer :: Calloc_leaf_end    ! fPHU at which C allocation to leaves stops
    real(dp), dimension(:), pointer :: Calloc_prod_max    ! fPHU at which C allocation coefficient to products is 1
    real(dp), dimension(:), pointer :: Cstem_mobile       ! Carbon pool in the stem that is used as mobile reserves
    

    ! Senescence and remobilisation
    real(dp), dimension(:),   pointer :: fCstem_mobile   ! fraction of C in the stem that goes to mobile reserves 
    real(dp), dimension(:,:), pointer :: Cmax        ! max. C stored in plant pools
    real(dp), dimension(:),   pointer :: Cmaxstemmob ! max. C stored in mobile reserves of stem
    real(dp), dimension(:),   pointer :: fmobilise ! fraction of C mobilised from stem to products
    real(dp), dimension(:),   pointer :: fsenesc  ! fraction of the crop canopy that is senescent (0-1)
    real(dp), dimension(:),   pointer :: LAItot   ! total (green + brown) LAI
    real(dp), dimension(:),   pointer :: LAIsen   ! senescent (brown) LAI
    real(dp), dimension(:),   pointer :: drsen    ! exponent in leaf senescence function
    
    ! Sowing dates (DOY)
    integer, dimension(:), pointer  :: sowing_doymin ! earliest sowing date
    integer, dimension(:), pointer  :: sowing_doymax ! latest sowing date
    integer, dimension(:), pointer  :: sowing_doy    ! actual sowing date

    ! Germination requirements
    real(dp), dimension(:), pointer :: fgermination  ! total germination requirements (0-1)
    
    ! Management options
    real(dp), dimension(:), pointer :: Cseed         ! carbon content in seeds (gC m-2)
                                                     ! --> corresponds to planting density
    real(dp), dimension(:), pointer :: sowing_depth  ! sowing depth (m)
    real(dp), dimension(:), pointer :: Cplant_remove ! fraction of aboveground biomass (stem and leaves)
                                                     ! removed at harvest. Remainder goes to litter pool
    
    !
    integer, dimension(:), pointer :: sl  ! soil layer corresponding to sowing depth

    ! Harvest variables
    real(dp), dimension(:), pointer :: yield          ! yield at harvest (gC m-2)
    real(dp), dimension(:), pointer :: harvest_index  ! harvest index (Cproduct/total Cplant)

    
    ! Leaf structural variables
    real(dp), dimension(:), pointer :: sla_maturity     ! specific leaf area (m2/gDM) at plant maturity
    real(dp), dimension(:), pointer :: sla_beta         ! extinction cofficient in SLA formula


    
  end type crop_type


  
Contains


  
  subroutine allocate_init_cropvars(crop,filename)
  ! allocates and initialises crop variables in one routine
    
    type(crop_type),      intent(inout)   :: crop
    type(filenames_type), intent(in)      :: filename

    !local
    integer :: j, jcrop   ! loop counter
    integer :: ioerror    ! input error integer
    integer :: vernalisation, photoperiodism ! to be converted to logical
    character(len=25) :: cropnametmp  ! not sure if needed/useful


    open(40,FILE=filename%crop,STATUS='old',ACTION='READ',IOSTAT=ioerror)

    if (ioerror/=0) then
      STOP 'CABLE_log: Cannot open crop parameter file!'
    end if 

    read(40,*)
    read(40,*) ncmax   ! number of crop functional types (CFTs)

    
    ! allocate crop structure  
    allocate(crop%state(ncmax),             &
             crop%vernalisation(ncmax),     &
             crop%photoperiodism(ncmax),    &
             crop%Tbase(ncmax),             &
             crop%Tmax(ncmax),              &
             crop%PHU_germination(ncmax),   &
             crop%fPHU_emergence(ncmax),    &
             crop%PHU_maturity(ncmax),      &
             crop%fPHU(ncmax),              &
             crop%VU(ncmax),                &
             crop%fVU(ncmax),               &
             crop%vacc(ncmax),              &
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


      
    ! initialise parameter values given in crop parameter file
    do j=1,ncmax

      read(40,*) jcrop, cropnametmp

      if (jcrop > ncmax) STOP 'jcrop out of range in parameter file'


      crop%vernalisation(jcrop) = .FALSE.
      crop%vernalisation(jcrop) = .FALSE.

      ! Read actual parameter values
      read(40,*) vernalisation, photoperiodism
      read(40,*) crop%Tbase(jcrop), crop%Tmax(jcrop)
      read(40,*) crop%PHU_germination(jcrop), crop%fPHU_emergence(jcrop)
      read(40,*) crop%PHU_maturity(jcrop)
      read(40,*) crop%fCalloc_root_init(jcrop), crop%fCalloc_leaf_init(jcrop), crop%fCalloc_leaf_bpt(jcrop)
      read(40,*) crop%Calloc_root_end(jcrop), crop%Calloc_leaf_bpt(jcrop), crop%Calloc_leaf_end(jcrop), crop%Calloc_prod_max(jcrop)
      read(40,*) crop%fCstem_mobile(jcrop)
      read(40,*) crop%Cseed(jcrop), crop%sowing_depth(jcrop), crop%Cplant_remove(jcrop)
      read(40,*) crop%sla_maturity(jcrop), crop%sla_beta(jcrop), crop%drsen(jcrop)

      if (vernalisation > 0) then
         crop%vernalisation(jcrop) = .TRUE.
      endif
      if (photoperiodism > 0) then
         crop%photoperiodism(jcrop) = .TRUE.
      endif
      
    end do ! loop over CFTs

    close(40)

    ! initialise all other variables
    crop%state         = 0
    crop%fPHU          = 0.0_dp
    crop%VU            = 0.0_dp
    crop%fVU           = 0.0_dp
    crop%vacc          = .FALSE.
    crop%Cstem_mobile  = 0.0_dp
    crop%Cmax          = 0.0_dp
    crop%Cmaxstemmob   = 0.0_dp
    crop%fmobilise     = 0.0_dp
    crop%fsenesc       = 0.0_dp
    crop%LAItot        = 0.0_dp
    crop%LAIsen        = 0.0_dp
    crop%sowing_doymin = 270
    crop%sowing_doymax = 290
    crop%sowing_doy    = 0
    crop%fgermination  = 0.0_dp
    crop%sl            = 0
    crop%yield         = 0.0_dp
    crop%harvest_index = 0.0_dp
    

  end subroutine allocate_init_cropvars

  
END MODULE crop_def
