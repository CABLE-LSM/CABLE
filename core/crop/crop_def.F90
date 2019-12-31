MODULE crop_def

  use cable_def_types_mod, only: dp => r_2
  use cable_common_module, only: filenames_type
  
  
  implicit none

  save

  public

  
  ! variables/parameters
  integer :: nc=1     ! number of crop types in grid cell
  integer :: ncmax    ! total number of crop functional types (CFTs)
                      ! read from crop parameter file (filename%crop)

  integer, parameter :: maxdays_ger=40 ! maximum days since sowing above which germination
                                       ! is assumed to have failed
  real(dp), parameter :: Cinit_root=0.4_dp  ! initial C allocation to roots at emergence
  real(dp), parameter :: Cinit_stem=0.0_dp  ! initial C allocation to stems at emergence
  real(dp), parameter :: Cinit_leaf=0.6_dp  ! initial C allocation to leaves at emergence
  real(dp), parameter :: DMtoC=0.475_dp     ! conversion dry matter to carbon

  
  ! types 
  type crop_type

    ! crop stages:
    ! 0: bare soil, not sown
    ! 1: sown (but not yet germinated)
    ! 2: emergent
    ! 3: growing
    integer, dimension(:),  pointer :: state
    ! Crop temperature requirements
    real, dimension(:), pointer :: Tbase     ! crop-specific base temperature

    ! Phenological heat units (PHU)
    real(dp), dimension(:), pointer :: PHU_germination  ! PHU (soil) until germination
                                                        ! in the absence of water stress
    real(dp), dimension(:), pointer :: fPHU_emergence   ! fPHU (air) until emergence (use of seed storage)
                                                        ! is completed
    real(dp), dimension(:), pointer :: PHU_maturity     ! PHU (air) until maturity
    real(dp), dimension(:), pointer :: fPHU             ! actual PHU relative to maturity (growth)

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
    character(len=25) :: cropnametmp  ! not sure if needed/useful


    open(40,FILE=filename%crop,STATUS='old',ACTION='READ',IOSTAT=ioerror)

    if (ioerror/=0) then
      STOP 'CABLE_log: Cannot open crop parameter file!'
    end if 

    read(40,*)
    read(40,*) ncmax   ! number of crop functional types (CFTs)

    
    ! allocate crop structure  
    allocate(crop%state(ncmax),           &
             crop%Tbase(ncmax),           &
             crop%PHU_germination(ncmax), &
             crop%fPHU_emergence(ncmax),  &
             crop%PHU_maturity(ncmax),    &
             crop%fPHU(ncmax),            &
             crop%sowing_doymin(ncmax),   &
             crop%sowing_doymax(ncmax),   &
             crop%sowing_doy(ncmax),      &
             crop%fgermination(ncmax),    &
             crop%Cseed(ncmax),           &
             crop%sowing_depth(ncmax),    &
             crop%Cplant_remove(ncmax),   &
             crop%sl(ncmax),              &
             crop%yield(ncmax),           &
             crop%harvest_index(ncmax),   &
             crop%sla_maturity(ncmax),    &
             crop%sla_beta(ncmax)         &
            )

      
    ! initialise parameter values given in crop parameter file
    do j=1,ncmax

      read(40,*) jcrop, cropnametmp

      if (jcrop > ncmax) STOP 'jcrop out of range in parameter file'

      ! Read actual parameter values
      read(40,*) crop%Tbase(jcrop)
      read(40,*) crop%PHU_germination(jcrop), crop%fPHU_emergence(jcrop)
      read(40,*) crop%PHU_maturity(jcrop)
      read(40,*) crop%Cseed(jcrop), crop%sowing_depth(jcrop), crop%Cplant_remove(jcrop)
      read(40,*) crop%sla_maturity(jcrop), crop%sla_beta(jcrop)

    end do ! loop over CFTs

    close(40)

      
    ! initialise all other variables
    crop%state         = 0
    crop%fPHU          = 0.0_dp
    crop%sowing_doymin = 270
    crop%sowing_doymax = 290
    crop%sowing_doy    = 0
    crop%fgermination  = 0.0_dp
    crop%sl            = 0
    crop%yield         = 0.0_dp
    crop%harvest_index = 0.0_dp
    

  end subroutine allocate_init_cropvars

  
END MODULE crop_def
