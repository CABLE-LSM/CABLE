! ==============================================================================
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
! Purpose: defines/allocates variables for CASA-CNP
!
! Contact: Yingping.Wang@csiro.au
!
! History: Developed for offline CASA-CNP, code revision likely to better
!          suit ACCESS and to merge more consistently with CABLE code
!
!
! ==============================================================================
! casa_variable.f90
!
! the following modules are used when "casacnp" is coupled to "cable"
!   casadimension
!   casaparm
!   casavariable with subroutine alloc_casa_var
!   phenvariable with subroutine alloc_phenvariable

module casadimension

  use cable_def_types_mod, only: r_2

  implicit none

  public

  integer, parameter :: mdyear = 365  ! days per year
  integer, parameter :: mdmonth = 30  ! days per month
  integer, parameter :: mdweek = 7    ! days per week
  integer, parameter :: mmyear = 12   ! month per year
  integer, parameter :: mt = 36500    ! integration time step
  integer, parameter :: mpftmax = 2   ! max. PFT/cell
  integer, parameter :: mplant   = 3  ! plant pools
  integer, parameter :: mlitter  = 3  ! litter pools
  integer, parameter :: msoil    = 3  ! soil pools
  integer, parameter :: mso      = 12 ! soil order number
  integer, parameter :: mhwp     = 1  ! harvested wood pools
  integer, parameter :: mclear   = 1  ! forest clearing pools
  integer, parameter :: mheights = 10 ! height clas
  ! BP put icycle into namelist file
  integer            :: icycle
  ! INTEGER, PARAMETER :: icycle=3    ! =1 for C, =2 for C+N; =3 for C+N+P
  integer, parameter :: mstart = 1    ! starting time step
  integer, parameter :: mphase = 4    ! phen. phases
  real(r_2), parameter :: deltcasa = 1.0_r_2/365.0_r_2 ! fraction 1 day of year
  real(r_2), parameter :: deltpool = 1.0_r_2           ! pool delt(1day)

end module casadimension


! ------------------------------------------------------------------


module casaparm

  use cable_def_types_mod, only: r_2
  ! use casadimension, only: deltcasa

  implicit none

  public

  integer, parameter :: initcasa = 1  ! =0 spin; 1 restart file
  integer, parameter :: iceland  = 17 !=13 for casa vegtype =15 for IGBP vegtype
  integer, parameter :: cropland = 9  ! 12 and 14 for IGBP vegtype
  integer, parameter :: croplnd2 = 10 ! ditto
  integer, parameter :: forest   = 3
  integer, parameter :: shrub    = 2
  integer, parameter :: grass    = 1
  integer, parameter :: icewater = 0
  integer, parameter :: LEAF     = 1
  integer, parameter :: WOOD     = 2
  integer, parameter :: FROOT    = 3
  ! integer, parameter :: LABILE = 4
  integer, parameter :: METB    = 1
  integer, parameter :: STR     = 2
  integer, parameter :: CWD     = 3
  integer, parameter :: MIC     = 1
  integer, parameter :: SLOW    = 2
  integer, parameter :: pass    = 3
  integer, parameter :: PLAB    = 1
  integer, parameter :: PSORB   = 2
  integer, parameter :: POCC    = 3
  !! vh_js !! LALLOC moved to bgcdriver to allow for value to be switchable
  ! integer, parameter :: LALLOC = 0  ! 0 constant; 1 variable
  real(r_2), parameter :: z30 = 0.3_r_2
  real(r_2), parameter :: R0 = 0.3_r_2
  real(r_2), parameter :: S0 = 0.3_r_2
  real(r_2), parameter :: fixed_stem = 1.0_r_2 / 3.0_r_2
  real(r_2), parameter :: Q10alloc = 2.0_r_2
  real(r_2), parameter :: ratioNCstrfix = 1.0_r_2 / 150.0_r_2
  real(r_2), parameter :: ratioNPstrfix = 25.0_r_2
  real(r_2), parameter :: fracCbiomass = 0.50_r_2
  real(r_2), parameter :: tsoilrefc = 25.0_r_2
  real(r_2), parameter :: tkzeroc = 273.15_r_2
  real(r_2), parameter :: frootparma = 0.3192_r_2
  real(r_2), parameter :: frootparmb =-0.0485_r_2
  real(r_2), parameter :: frootparmc = 0.1755_r_2
  real(r_2), parameter :: xweightalloc = 0.2_r_2
  !  real(r_2), parameter :: xkplab  = 0.5_r_2 * deltcasa
  !  real(r_2), parameter :: xkpsorb = 0.01_r_2 * deltcasa
  !  real(r_2), parameter :: xkpocc  = 0.01_r_2 * deltcasa

end module casaparm


! ------------------------------------------------------------------


module casavariable

  use cable_def_types_mod, only: r_2

  implicit none

  private

  ! types
  public :: casa_balance
  public :: casa_biome
  public :: casa_flux
  public :: casa_met
  public :: casa_pool
  public :: casafiles_type

  ! routines on types
  public :: alloc_casa_var
  public :: print_casa_var
  public :: read_netcdf_casa_var
  public :: write_netcdf_casa_var
  public :: zero_casa_var

  ! routines on sum variables
  public :: alloc_sum_casa
  public :: update_sum_casa
  public :: zero_sum_casa

  ! public variables
  public :: casa_timeunits
  public :: casafile

  ! number of variables in type definitions
  ! used in write_netcdf and in MPI code
  integer, parameter, public :: ncasa_biome = 53
  integer, parameter, public :: ncasa_pool = 42
  integer, parameter, public :: ncasa_flux = 91
  integer, parameter, public :: ncasa_met = 47
  integer, parameter, public :: ncasa_bal = 47

  ! private section

  character(len=200) :: casa_timeunits

  type casa_biome
     integer, dimension(:), pointer :: &
          ivt2 => null()
     real(r_2), dimension(:), pointer :: &
          xkleafcoldmax => null(), &
          xkleafcoldexp => null(), &
          xkleafdrymax => null(), &
          xkleafdryexp => null(), &
          glaimax => null(), &
          glaimin => null(), &
          sla => null(), &
          ratiofrootleaf => null(), &
          kroot => null(), &
          krootlen => null(), &
          rootdepth => null(), &
          kuptake => null(), &
          kminN => null(), &
          kuplabP => null(), &
          kclabrate => null(), &
          xnpmax => null(), &
          q10soil => null(), &
          xkoptlitter => null(), &
          xkoptsoil => null(), &
          xkplab => null(), &
          xkpsorb => null(), &
          xkpocc => null(), &
          prodptase => null(), &
          costnpup => null(), &
          maxfinelitter => null(), &
          maxcwd => null(), &
          nintercept => null(), &
          nslope => null(), &
          la_to_sa => null(), &
          vcmax_scalar => null(), &
          disturbance_interval => null(), &
          DAMM_EnzPool => null(), &
          DAMM_KMO2 => null(), &
          DAMM_KMcp => null(), &
          DAMM_Ea => null(), &
          DAMM_alpha => null()
     real(r_2), dimension(:,:), pointer :: &
          plantrate => null(), &
          rmplant => null(), &
          fracnpptoP => null(), &
          fraclignin => null(), &
          fraclabile => null(), &
          ratioNCplantmin => null(), &
          ratioNCplantmax => null(), &
          ratioNPplantmin => null(), &
          ratioNPplantmax => null(), &
          fracLigninplant => null(), &
          ftransNPtoL => null(), &
          ftransPPtoL => null(), &
          litterrate => null(), &
          ratioPcplantmin => null(), &
          ratioPcplantmax => null()
     real(r_2), dimension(:,:), pointer :: &
          soilrate => null()
  end type casa_biome


  type casa_pool
     real(r_2), dimension(:), pointer :: &
          Clabile => null(), &
          dClabiledt => null(), &
          Ctot => null(), &          !! vh_js !!
          Ctot_0 => null()
     real(r_2), dimension(:,:), pointer :: &
          Cplant => null(), &
          Nplant => null(), &
          Pplant => null(), &
          dCplantdt => null(), &
          dNplantdt => null(), &
          dPplantdt => null(), &
          ratioNCplant => null(), &
          ratioNPplant => null()
     real(r_2), dimension(:), pointer :: &
          Nsoilmin => null(), &
          Psoillab => null(), &
          Psoilsorb => null(), &
          Psoilocc => null(), &
          dNsoilmindt => null(), &
          dPsoillabdt => null(), &
          dPsoilsorbdt => null(), &
          dPsoiloccdt => null()
     real(r_2), dimension(:,:), pointer :: &
          Clitter => null(), &
          Nlitter => null(), &
          Plitter => null(), &
          dClitterdt => null(), &
          dNlitterdt => null(), &
          dPlitterdt => null(), &
          ratioNClitter => null(), &
          ratioNPlitter => null()
     real(r_2), dimension(:,:), pointer :: &
          Csoil => null(), &
          Nsoil => null(), &
          Psoil => null(), &
          dCsoildt => null(), &
          dNsoildt => null(), &
          dPsoildt => null(), &
          ratioNCsoil => null(), &
          ratioNCsoilnew => null(), &
          ratioNPsoil => null(), &
          ratioNCsoilmin => null(), &
          ratioNCsoilmax => null(), &
          ratioPCsoil => null(), &
          ratioPCplant => null(), &
          ratioPClitter => null()
  end type casa_pool


  type casa_flux
     real(r_2), dimension(:), pointer :: &
          Cgpp => null(), &
          Cnpp => null(), &
          Crp => null(), &
          Crgplant => null(), &
          Nminfix => null(), &
          Nminuptake => null(), &
          Plabuptake => null(), &
          Clabloss => null(), &
          fracClabile => null(), &
          !! vh_js !! the 3 variables below are needed for POP coupling to CASA
          stemnpp => null(), &
          frac_sapwood => null(), &
          sapwood_area => null(), &
          Charvest => null(), &  ! leaf biomass removed due to crop or pasture management
          Nharvest => null(), & ! leaf N removed due to crop or pasture management
          Pharvest => null(), & ! leaf P removed due to crop or pasture management
          fHarvest => null(), &  ! fraction leaf biomass removed due to crop or pasture management
          fcrop => null()        ! fraction of 'grass' that is crop
     real(r_2), dimension(:,:), pointer :: &
          fracCalloc => null(), &
          fracNalloc => null(), &
          fracPalloc => null(), &
          Crmplant => null(), &
          kplant => null(), &
          !! vh_js !! additional diagnostic
          Cplant_turnover => null()
     real(r_2), dimension(:,:,:), pointer :: &
          fromPtoL => null()
     real(r_2), dimension(:),pointer :: &
          Cnep => null(), &
          Crsoil => null(), &
          Nmindep => null(), &
          Nminloss => null(), &
          Nminleach => null(), &
          Nupland => null(), &
          Nlittermin => null(), &
          Nsmin => null(), &
          Nsimm => null(), &
          Nsnet => null(), &
          fNminloss => null(), &
          fNminleach => null(), &
          Pdep => null(), &
          Pwea => null(), &
          Pleach => null(), &
          Ploss => null(), &
          Pupland => null(), &
          Plittermin => null(), &
          Psmin => null(), &
          Psimm => null(), &
          Psnet => null(), &
          fPleach => null(), &
          kplab => null(), &
          kpsorb => null(), &
          kpocc => null(), &
          kmlabp => null(), &
          Psorbmax => null(), &
          !! additional diagnostics for partitioning biomass turnover
          Cplant_turnover_disturbance => null(), &
          Cplant_turnover_crowding  => null(), &
          Cplant_turnover_resource_limitation => null()
     real(r_2), dimension(:,:),   pointer :: klitter => null()
     real(r_2), dimension(:,:),   pointer :: ksoil => null()
     real(r_2), dimension(:,:,:), pointer :: fromLtoS => null()
     real(r_2), dimension(:,:,:), pointer :: fromStoS => null()
     real(r_2), dimension(:,:),   pointer :: fromLtoCO2 => null()
     real(r_2), dimension(:,:),   pointer :: fromStoCO2 => null()
     real(r_2), dimension(:,:),   pointer :: FluxCtolitter => null()
     real(r_2), dimension(:,:),   pointer :: FluxNtolitter => null()
     real(r_2), dimension(:,:),   pointer :: FluxPtolitter => null()
     real(r_2), dimension(:,:),   pointer :: FluxCtosoil => null()
     real(r_2), dimension(:,:),   pointer :: FluxNtosoil => null()
     real(r_2), dimension(:,:),   pointer :: FluxPtosoil => null()
     real(r_2), dimension(:),     pointer :: FluxCtoCO2 => null()
     real(r_2), dimension(:),     pointer :: FluxCtohwp => null()
     real(r_2), dimension(:),     pointer :: FluxNtohwp => null()
     real(r_2), dimension(:),     pointer :: FluxPtohwp => null()
     real(r_2), dimension(:),     pointer :: FluxCtoclear => null()
     real(r_2), dimension(:),     pointer :: FluxNtoclear => null()
     real(r_2), dimension(:),     pointer :: FluxPtoclear => null()
     real(r_2), dimension(:),     pointer :: CtransferLUC => null()
     !CVH variables inherited from BLAZE
     real(r_2), dimension(:,:,:), pointer :: fromPtoL_fire => null()
     real(r_2), dimension(:,:), pointer   :: klitter_fire => null()
     ! sum of fire turnover and non-fire turnover (litter)
     real(r_2), dimension(:,:), pointer   :: klitter_tot => null()
     real(r_2), dimension(:,:), pointer   :: kplant_fire => null()
     ! sum of fire turnover and non-fire turnover (plants)
     real(r_2), dimension(:,:), pointer   :: kplant_tot => null()
     !CVH diagnostic: CO2 emissions from fire
     real(r_2), dimension(:), pointer     :: fluxCtoCO2_plant_fire => null()
     real(r_2), dimension(:), pointer     :: fluxCtoCO2_litter_fire => null()
     ! contribution to fire emissions from individual plant pools
     real(r_2), dimension(:,:), pointer   :: fluxfromPtoCO2_fire => null()
     ! contribution to fire emissions from individual litter pools
     real(r_2), dimension(:,:), pointer   :: fluxfromLtoCO2_fire => null()
     real(r_2), dimension(:), pointer     :: fluxNtoAtm_fire => null()
     ! real(r_2), dimension(:,:,:),pointer  :: fire_mortality_vs_height => null()
     ! Diagnostic fluxes for use in 13C
     real(r_2), dimension(:,:,:), pointer :: FluxFromPtoL => null()
     real(r_2), dimension(:,:,:), pointer :: FluxFromLtoS => null()
     real(r_2), dimension(:,:,:), pointer :: FluxFromStoS => null()
     real(r_2), dimension(:,:),   pointer :: FluxFromPtoCO2 => null()
     real(r_2), dimension(:,:),   pointer :: FluxFromLtoCO2 => null()
     real(r_2), dimension(:,:),   pointer :: FluxFromStoCO2 => null()
     real(r_2), dimension(:),     pointer :: FluxFromPtoHarvest => null()
  end type casa_flux


  type casa_met
     real(r_2), dimension(:), pointer :: &
          glai => null(), &
          Tairk => null(), &
          precip => null(), &
          tsoilavg => null(), &
          moistavg => null(), &
          btran => null()
     integer, dimension(:), pointer :: &
          lnonwood => null()
     real(r_2), dimension(:,:), pointer :: &
          Tsoil => null(), &
          moist => null()
     integer, dimension(:), pointer :: &
          iveg2 => null(), &
          ijgcm => null(), &
          isorder => null()
     real(r_2), dimension(:), pointer :: &
          lat => null(), &
          lon => null(), &
          areacell => null()
     ! added yp wang 5/nov/2012
     real(r_2), dimension(:,:), pointer :: &
          Tairkspin => null(), &
          cgppspin => null(), &
          crmplantspin_1 => null(), &
          crmplantspin_2 => null(), &
          crmplantspin_3 => null(), &
          Tsoilspin_1 => null(), &
          Tsoilspin_2 => null(), &
          Tsoilspin_3 => null(), &
          Tsoilspin_4 => null(), &
          Tsoilspin_5 => null(), &
          Tsoilspin_6 => null(), &
          moistspin_1 => null(), &
          moistspin_2 => null(), &
          moistspin_3 => null(), &
          moistspin_4 => null(), &
          moistspin_5 => null(), &
          moistspin_6 => null(), &
          mtempspin => null(), &
          frecspin => null()
     ! 13C
     real(r_2), dimension(:,:), pointer :: &
          ! daily cumulated total 12CO2 net assimilation in [g(C)/m2]
          cAn12spin => null()
     real(r_2), dimension(:,:), pointer :: &
          ! daily cumulated total 13CO2 net assimilation in [g(13C)/m2]
          cAn13spin => null()
     ! BLAZE
     real(r_2), dimension(:,:), pointer :: &
          dprecip_spin => null(), &
          aprecip_av20_spin => null(), &
          du10_max_spin     => null(), &
          drhum_spin        => null(), &
          dtemp_max_spin    => null(), &
          dtemp_min_spin    => null(), &
          KBDI_spin         => null(), &
          D_MacArthur_spin  => null(), &
          FFDI_spin         => null(), &
          last_precip_spin  => null()
     integer, dimension(:,:), pointer :: &
          DSLR_spin => null()
  end type casa_met


  type casa_balance
     real(r_2), dimension(:), pointer :: &
          FCgppyear => null(), &
          FCnppyear => null(), &
          FCrmleafyear => null(), &
          FCrmwoodyear => null(), &
          FCrmrootyear => null(), &
          FCrgrowyear => null(), &
          FCrpyear => null(), &
          FCrsyear => null(), &
          FCneeyear => null(), &
          dCdtyear => null(), &
          LAImax => null(), &
          Cleafmean => null(), &
          Crootmean => null(), &
          FNdepyear => null(), &
          FNfixyear => null(), &
          FNsnetyear => null(), &
          FNupyear => null(), &
          FNleachyear => null(), &
          FNlossyear => null(), &
          FPweayear => null(), &
          FPdustyear => null(), &
          FPsnetyear => null(), &
          FPupyear => null(), &
          FPleachyear => null(), &
          FPlossyear => null()
     real(r_2), dimension(:,:),pointer :: &
          glaimon => null(), &
          glaimonx => null()
     real(r_2), dimension(:,:), pointer :: &
          cplantlast => null(), &
          nplantlast => null(), &
          pplantlast => null()
     real(r_2), dimension(:,:), pointer :: &
          clitterlast => null(), &
          nlitterlast => null(), &
          plitterlast => null()
     real(r_2), dimension(:,:), pointer :: &
          csoillast => null(), &
          nsoillast => null(), &
          psoillast => null()
     real(r_2), dimension(:), pointer :: &
          nsoilminlast => null(), &
          psoillablast => null(), &
          psoilsorblast => null(), &
          psoilocclast => null(), &
          cbalance => null(), &
          nbalance => null(), &
          pbalance => null(), &
          sumcbal => null(), &
          sumnbal => null(), &
          sumpbal => null()
     real(r_2), dimension(:), pointer :: &
          clabilelast => null()
  end type casa_balance


  ! The following declarations are removed and have to be passed using
  ! parameter list for each subroutine (BP apr2010)
  !  TYPE (casa_biome)              :: casabiome
  !  TYPE (casa_pool)               :: casapool
  !  TYPE (casa_flux)               :: casaflux
  !  TYPE (casa_met)                :: casamet
  !  TYPE (casa_balance)            :: casabal


  ! Added filename type for casaCNP (BP apr2010)
  type casafiles_type
     character(LEN=200) :: cnpbiome    ! file for biome-specific BGC parameters
     character(LEN=200) :: cnppoint    ! file for point-specific BGC inputs
     character(LEN=200) :: cnpepool    ! file for end-of-run pool sizes
     character(LEN=200) :: cnpipool=''    ! file for inital pool sizes
     character(LEN=200) :: cnpmetin      ! met file for spin up
     character(LEN=200) :: cnpmetout     ! met file for spin up
     character(LEN=200) :: ndep          ! N deposition input file
     ! added yp wang
     character(LEN=200) :: cnpspin       ! input file for spin up
     character(LEN=200) :: dump_cnpspin  ! name of dump file for spinning casa-cnp

     character(LEN=200) :: phen        ! leaf phenology datafile
     character(LEN=200) :: cnpflux     ! modelled mean yearly CNP fluxes
     logical            :: l_ndep
     ! added vh
     character(LEN=200) :: c2cdumppath='' ! cable2casa dump for casa spinup
     character(LEN=200) :: out=''    ! casa output file
  end type casafiles_type
  type(casafiles_type) :: casafile

  interface alloc_casa_var
     module procedure &
          alloc_casabiome, &
          alloc_casapool, &
          alloc_casaflux, &
          alloc_casamet, &
          alloc_casabal
  end interface alloc_casa_var

  interface zero_casa_var
     module procedure &
          zero_casabiome, &
          zero_casapool, &
          zero_casaflux, &
          zero_casamet, &
          zero_casabal
  end interface zero_casa_var

  interface print_casa_var
     module procedure &
          print_casabiome, &
          print_casapool, &
          print_casaflux, &
          print_casamet, &
          print_casabal
  end interface print_casa_var

  interface read_netcdf_casa_var
     module procedure &
          read_netcdf_casabiome, &
          read_netcdf_casapool, &
          read_netcdf_casaflux, &
          read_netcdf_casamet, &
          read_netcdf_casabal
  end interface read_netcdf_casa_var

  interface write_netcdf_casa_var
     module procedure &
          write_netcdf_casabiome, &
          write_netcdf_casapool, &
          write_netcdf_casaflux, &
          write_netcdf_casamet, &
          write_netcdf_casabal
  end interface write_netcdf_casa_var

contains

  ! ------------------------------------------------------------------

  subroutine alloc_casabiome(casabiome)

    use cable_def_types_mod, only: mvtype
    use casadimension, only: mplant, mlitter, msoil, mso

    implicit none

    type(casa_biome), intent(inout) :: casabiome

    allocate( &
         casabiome%ivt2(mvtype), &
         casabiome%xkleafcoldmax(mvtype), &
         casabiome%xkleafcoldexp(mvtype), &
         casabiome%xkleafdrymax(mvtype), &
         casabiome%xkleafdryexp(mvtype), &
         casabiome%glaimax(mvtype), &
         casabiome%glaimin(mvtype), &
         casabiome%sla(mvtype), &
         casabiome%ratiofrootleaf(mvtype), &
         casabiome%kroot(mvtype), &
         casabiome%krootlen(mvtype), &
         casabiome%rootdepth(mvtype), &
         casabiome%kuptake(mvtype), &
         casabiome%kminN(mvtype), &
         casabiome%KuplabP(mvtype), &
         casabiome%kclabrate(mvtype), &
         casabiome%xnpmax(mvtype), &
         casabiome%q10soil(mvtype), &
         casabiome%xkoptlitter(mvtype), &
         casabiome%xkoptsoil(mvtype), &
         casabiome%xkplab(mso), &
         casabiome%xkpsorb(mso), &
         casabiome%xkpocc(mso), &
         casabiome%prodptase(mvtype), &
         casabiome%costnpup(mvtype), &
         casabiome%maxfinelitter(mvtype), &
         casabiome%maxcwd(mvtype), &
         casabiome%nintercept(mvtype), &
         casabiome%nslope(mvtype), &
         casabiome%la_to_sa(mvtype), &
         casabiome%vcmax_scalar(mvtype), &
         casabiome%disturbance_interval(mvtype), &
         casabiome%DAMM_EnzPool(mvtype), &
         casabiome%DAMM_KMO2(mvtype), &
         casabiome%DAMM_KMcp(mvtype), &
         casabiome%DAMM_Ea(mvtype), &
         casabiome%DAMM_alpha(mvtype), &
         casabiome%plantrate(mvtype,mplant), &
         casabiome%rmplant(mvtype,mplant), &
         casabiome%fracnpptoP(mvtype,mplant), &
         casabiome%fraclignin(mvtype,mplant), &
         casabiome%fraclabile(mvtype,mplant), &
         casabiome%ratioNCplantmin(mvtype,mplant), &
         casabiome%ratioNCplantmax(mvtype,mplant), &
         casabiome%ratioNPplantmin(mvtype,mplant), &
         casabiome%ratioNPplantmax(mvtype,mplant), &
         casabiome%fracLigninplant(mvtype,mplant), &
         casabiome%ftransNPtoL(mvtype,mplant), &
         casabiome%ftransPPtoL(mvtype,mplant), &
         casabiome%litterrate(mvtype,mlitter), &
         !  casabiome%ratioPcplantmax(mvtype,leaf), &
         !  casabiome%ratioPcplantmin(mvtype,leaf) &
         !! vh_js !!
         casabiome%ratioPcplantmax(mvtype,mplant), &
         casabiome%ratioPcplantmin(mvtype,mplant), &
         casabiome%soilrate(mvtype,msoil) &
         )

  end subroutine alloc_casabiome


  subroutine alloc_casapool(casapool, arraysize)

    use casadimension, only: mplant, mlitter, msoil

    implicit none

    type(casa_pool), intent(inout) :: casapool
    integer,         intent(in)    :: arraysize

    allocate( &
         casapool%Clabile(arraysize), &
         casapool%dClabiledt(arraysize), &
         casapool%Ctot(arraysize), &
         casapool%Ctot_0(arraysize), &
         casapool%Cplant(arraysize,mplant), &
         casapool%Nplant(arraysize,mplant), &
         casapool%Pplant(arraysize,mplant), &
         casapool%dCplantdt(arraysize,mplant), &
         casapool%dNplantdt(arraysize,mplant), &
         casapool%dPplantdt(arraysize,mplant), &
         casapool%ratioNCplant(arraysize,mplant), &
         casapool%ratioNPplant(arraysize,mplant), &
         casapool%Nsoilmin(arraysize), &
         casapool%Psoillab(arraysize), &
         casapool%Psoilsorb(arraysize), &
         casapool%Psoilocc(arraysize), &
         casapool%dNsoilmindt(arraysize), &
         casapool%dPsoillabdt(arraysize), &
         casapool%dPsoilsorbdt(arraysize), &
         casapool%dPsoiloccdt(arraysize), &
         casapool%Clitter(arraysize,mlitter), &
         casapool%Nlitter(arraysize,mlitter), &
         casapool%Plitter(arraysize,mlitter), &
         casapool%dClitterdt(arraysize,mlitter), &
         casapool%dNlitterdt(arraysize,mlitter), &
         casapool%dPlitterdt(arraysize,mlitter), &
         casapool%ratioNClitter(arraysize,mlitter), &
         casapool%ratioNPlitter(arraysize,mlitter), &
         casapool%Csoil(arraysize,msoil), &
         casapool%Nsoil(arraysize,msoil), &
         casapool%Psoil(arraysize,msoil), &
         casapool%dCsoildt(arraysize,msoil), &
         casapool%dNsoildt(arraysize,msoil), &
         casapool%dPsoildt(arraysize,msoil), &
         casapool%ratioNCsoil(arraysize,msoil), &
         casapool%ratioNPsoil(arraysize,msoil), &
         casapool%ratioNCsoilnew(arraysize,msoil), &
         casapool%ratioNCsoilmin(arraysize,msoil), &
         casapool%ratioNCsoilmax(arraysize,msoil), &
         casapool%ratioPCsoil(arraysize,msoil), &
         casapool%ratioPCplant(arraysize,mplant), &
         casapool%ratioPClitter(arraysize,mlitter) &
         )

  end subroutine alloc_casapool


  subroutine alloc_casaflux(casaflux, arraysize)

    use casadimension, only: mplant, mlitter, msoil

    implicit none

    type(casa_flux), intent(inout) :: casaflux
    integer,         intent(in)    :: arraysize

    allocate( &
         casaflux%Cgpp(arraysize), &
         casaflux%Cnpp(arraysize), &
         casaflux%Crp(arraysize), &
         casaflux%Crgplant(arraysize), &
         casaflux%Nminfix(arraysize), &
         casaflux%Nminuptake(arraysize), &
         casaflux%Plabuptake(arraysize), &
         casaflux%Clabloss(arraysize), &
         casaflux%fracClabile(arraysize), &
         casaflux%stemnpp(arraysize), &
         casaflux%frac_sapwood(arraysize), &
         casaflux%sapwood_area(arraysize), &
         casaflux%Charvest(arraysize), &
         casaflux%Nharvest(arraysize), &
         casaflux%Pharvest(arraysize), &
         casaflux%fharvest(arraysize), &
         casaflux%fcrop(arraysize), &
         casaflux%fracCalloc(arraysize,mplant), &
         casaflux%fracNalloc(arraysize,mplant), &
         casaflux%fracPalloc(arraysize,mplant), &
         casaflux%Crmplant(arraysize,mplant), &
         casaflux%kplant(arraysize,mplant), &
         casaflux%Cplant_turnover(arraysize,mplant), &
         casaflux%fromPtoL(arraysize,mlitter,mplant), &
         casaflux%Cnep(arraysize), &
         casaflux%Crsoil(arraysize), &
         casaflux%Nmindep(arraysize), &
         casaflux%Nminloss(arraysize), &
         casaflux%Nminleach(arraysize), &
         casaflux%Nupland(arraysize), &
         casaflux%Nlittermin(arraysize), &
         casaflux%Nsmin(arraysize), &
         casaflux%Nsimm(arraysize), &
         casaflux%Nsnet(arraysize), &
         casaflux%fNminloss(arraysize), &
         casaflux%fNminleach(arraysize), &
         casaflux%Pdep(arraysize), &
         casaflux%Pwea(arraysize), &
         casaflux%Pleach(arraysize), &
         casaflux%Ploss(arraysize), &
         casaflux%Pupland(arraysize), &
         casaflux%Plittermin(arraysize), &
         casaflux%Psmin(arraysize), &
         casaflux%Psimm(arraysize), &
         casaflux%Psnet(arraysize), &
         casaflux%fPleach(arraysize), &
         casaflux%kplab(arraysize), &
         casaflux%kpsorb(arraysize), &
         casaflux%kpocc(arraysize), &
         casaflux%kmlabP(arraysize), &
         casaflux%Psorbmax(arraysize), &
         casaflux%Cplant_turnover_disturbance(arraysize), &
         casaflux%Cplant_turnover_crowding(arraysize), &
         casaflux%Cplant_turnover_resource_limitation(arraysize), &
         casaflux%klitter(arraysize,mlitter), &
         casaflux%ksoil(arraysize,msoil), &
         casaflux%fromLtoS(arraysize,msoil,mlitter), &
         casaflux%fromStoS(arraysize,msoil,msoil), &
         casaflux%fromLtoCO2(arraysize,mlitter), &
         casaflux%fromStoCO2(arraysize,msoil), &
         casaflux%FluxCtolitter(arraysize,mlitter), &
         casaflux%FluxNtolitter(arraysize,mlitter), &
         casaflux%FluxPtolitter(arraysize,mlitter), &
         casaflux%FluxCtosoil(arraysize,msoil), &
         casaflux%FluxNtosoil(arraysize,msoil), &
         casaflux%FluxPtosoil(arraysize,msoil), &
         casaflux%FluxCtoCO2(arraysize), &
         casaflux%FluxCtohwp(arraysize), &
         casaflux%FluxNtohwp(arraysize), &
         casaflux%FluxPtohwp(arraysize), &
         casaflux%FluxCtoclear(arraysize), &
         casaflux%FluxNtoclear(arraysize), &
         casaflux%FluxPtoclear(arraysize), &
         casaflux%CtransferLUC(arraysize), &
         casaflux%fromPtoL_fire(arraysize,mlitter,mplant), &
         casaflux%klitter_fire(arraysize,mlitter), &
         casaflux%klitter_tot(arraysize,mlitter), &
         casaflux%kplant_fire(arraysize,mplant), &
         casaflux%kplant_tot(arraysize,mplant), &
         casaflux%FluxCtoCO2_plant_fire(arraysize), &
         casaflux%FluxCtoCO2_litter_fire(arraysize), &
         casaflux%fluxfromPtoCO2_fire(arraysize,mplant), &
         casaflux%fluxfromLtoCO2_fire(arraysize,mlitter), &
         casaflux%FluxNtoAtm_fire(arraysize), &
         ! casabal%fire_mortality_vs_height(arraysize,mheights,2) )
         casaflux%FluxFromPtoL(arraysize,mplant,mlitter), &
         casaflux%FluxFromLtoS(arraysize,mlitter,msoil), &
         casaflux%FluxFromStoS(arraysize,msoil,msoil), &
         casaflux%FluxFromPtoCO2(arraysize,mplant), &
         casaflux%FluxFromLtoCO2(arraysize,mlitter), &
         casaflux%FluxFromStoCO2(arraysize,msoil), &
         casaflux%FluxFromPtoHarvest(arraysize))

  end subroutine alloc_casaflux


  subroutine alloc_casamet(casamet, arraysize)

    use cable_def_types_mod, only: ms
    use casadimension, only: mdyear

    implicit none

    type(casa_met), intent(inout) :: casamet
    integer,        intent(in)    :: arraysize

    allocate( &
         casamet%glai(arraysize), &
         casamet%Tairk(arraysize), &
         casamet%precip(arraysize), &
         casamet%tsoilavg(arraysize), &
         casamet%moistavg(arraysize), &
         casamet%btran(arraysize), &
         casamet%lnonwood(arraysize), &
         casamet%Tsoil(arraysize,ms), &
         casamet%moist(arraysize,ms), &
         casamet%iveg2(arraysize), &
         casamet%ijgcm(arraysize), &
         casamet%isorder(arraysize), &
         casamet%lat(arraysize), &
         casamet%lon(arraysize), &
         casamet%areacell(arraysize), &
         casamet%Tairkspin(arraysize,mdyear), &
         casamet%cgppspin(arraysize,mdyear), &
         casamet%crmplantspin_1(arraysize,mdyear), &
         casamet%crmplantspin_2(arraysize,mdyear), &
         casamet%crmplantspin_3(arraysize,mdyear), &
         casamet%Tsoilspin_1(arraysize,mdyear), &
         casamet%Tsoilspin_2(arraysize,mdyear), &
         casamet%Tsoilspin_3(arraysize,mdyear), &
         casamet%Tsoilspin_4(arraysize,mdyear), &
         casamet%Tsoilspin_5(arraysize,mdyear), &
         casamet%Tsoilspin_6(arraysize,mdyear), &
         casamet%moistspin_1(arraysize,mdyear), &
         casamet%moistspin_2(arraysize,mdyear), &
         casamet%moistspin_3(arraysize,mdyear), &
         casamet%moistspin_4(arraysize,mdyear), &
         casamet%moistspin_5(arraysize,mdyear), &
         casamet%moistspin_6(arraysize,mdyear), &
         casamet%mtempspin(arraysize,mdyear), &
         casamet%frecspin(arraysize,mdyear), &
         casamet%cAn12spin(arraysize,mdyear), &
         casamet%cAn13spin(arraysize,mdyear), &
         casamet%dprecip_spin(arraysize,mdyear), &
         casamet%aprecip_av20_spin(arraysize,mdyear), &
         casamet%du10_max_spin(arraysize,mdyear), &
         casamet%drhum_spin(arraysize,mdyear), &
         casamet%dtemp_max_spin(arraysize,mdyear), &
         casamet%dtemp_min_spin(arraysize,mdyear), &
         casamet%KBDI_spin(arraysize,mdyear), &
         casamet%D_MacArthur_spin(arraysize,mdyear), &
         casamet%FFDI_spin(arraysize,mdyear), &
         casamet%last_precip_spin(arraysize,mdyear), &
         casamet%DSLR_spin(arraysize,mdyear))

  end subroutine alloc_casamet


  subroutine alloc_casabal(casabal, arraysize)

    use casadimension, only: mplant, mlitter, msoil

    implicit none

    type(casa_balance), intent(inout) :: casabal
    integer,            intent(in)    :: arraysize

    allocate( &
         casabal%FCgppyear(arraysize), &
         casabal%FCnppyear(arraysize), &
         casabal%FCrmleafyear(arraysize), &
         casabal%FCrmwoodyear(arraysize), &
         casabal%FCrmrootyear(arraysize), &
         casabal%FCrgrowyear(arraysize), &
         casabal%FCrpyear(arraysize), &
         casabal%FCrsyear(arraysize), &
         casabal%FCneeyear(arraysize), &
         casabal%dCdtyear(arraysize), &
         casabal%LAImax(arraysize), &
         casabal%Cleafmean(arraysize), &
         casabal%Crootmean(arraysize), &
         casabal%FNdepyear(arraysize), &
         casabal%FNfixyear(arraysize), &
         casabal%FNsnetyear(arraysize), &
         casabal%FNupyear(arraysize), &
         casabal%FNleachyear(arraysize), &
         casabal%FNlossyear(arraysize), &
         casabal%FPweayear(arraysize), &
         casabal%FPdustyear(arraysize), &
         casabal%FPsnetyear(arraysize), &
         casabal%FPupyear(arraysize), &
         casabal%FPleachyear(arraysize), &
         casabal%FPlossyear(arraysize), &
         casabal%glaimon(arraysize,12), &
         casabal%glaimonx(arraysize,12), &
         casabal%cplantlast(arraysize,mplant), &
         casabal%nplantlast(arraysize,mplant), &
         casabal%pplantlast(arraysize,mplant), &
         casabal%clitterlast(arraysize,mlitter), &
         casabal%nlitterlast(arraysize,mlitter), &
         casabal%plitterlast(arraysize,mlitter), &
         casabal%csoillast(arraysize,msoil), &
         casabal%nsoillast(arraysize,msoil), &
         casabal%psoillast(arraysize,msoil), &
         casabal%nsoilminlast(arraysize), &
         casabal%psoillablast(arraysize), &
         casabal%psoilsorblast(arraysize), &
         casabal%psoilocclast(arraysize), &
         casabal%cbalance(arraysize), &
         casabal%nbalance(arraysize), &
         casabal%pbalance(arraysize), &
         casabal%sumcbal(arraysize), &
         casabal%sumnbal(arraysize), &
         casabal%sumpbal(arraysize), &
         casabal%clabilelast(arraysize))

  end subroutine alloc_casabal

  ! ------------------------------------------------------------------

  subroutine alloc_sum_casa(sum_casapool, sum_casaflux, arraysize)

    implicit none

    type(casa_pool), intent(inout) :: sum_casapool
    type(casa_flux), intent(inout) :: sum_casaflux
    integer,         intent(in)    :: arraysize

    call alloc_casapool(sum_casapool, arraysize)
    call alloc_casaflux(sum_casaflux, arraysize)

  end subroutine alloc_sum_casa

  ! ------------------------------------------------------------------

  subroutine zero_casabiome(casabiome)

    use cable_def_types_mod, only: r_2

    implicit none

    type(casa_biome), intent(inout) :: casabiome

    casabiome%ivt2 = 0
    casabiome%xkleafcoldmax = 0.0_r_2
    casabiome%xkleafcoldexp = 0.0_r_2
    casabiome%xkleafdrymax = 0.0_r_2
    casabiome%xkleafdryexp = 0.0_r_2
    casabiome%glaimax = 0.0_r_2
    casabiome%glaimin = 0.0_r_2
    casabiome%sla = 0.0_r_2
    casabiome%ratiofrootleaf = 0.0_r_2
    casabiome%kroot = 0.0_r_2
    casabiome%krootlen = 0.0_r_2
    casabiome%rootdepth = 0.0_r_2
    casabiome%kuptake = 0.0_r_2
    casabiome%kminN = 0.0_r_2
    casabiome%kuplabP = 0.0_r_2
    casabiome%kclabrate = 0.0_r_2
    casabiome%xnpmax = 0.0_r_2
    casabiome%q10soil = 0.0_r_2
    casabiome%xkoptlitter = 0.0_r_2
    casabiome%xkoptsoil = 0.0_r_2
    casabiome%xkplab = 0.0_r_2
    casabiome%xkpsorb = 0.0_r_2
    casabiome%xkpocc = 0.0_r_2
    casabiome%prodptase = 0.0_r_2
    casabiome%costnpup = 0.0_r_2
    casabiome%maxfinelitter = 0.0_r_2
    casabiome%maxcwd = 0.0_r_2
    casabiome%nintercept = 0.0_r_2
    casabiome%nslope = 0.0_r_2
    casabiome%la_to_sa = 0.0_r_2
    casabiome%vcmax_scalar = 0.0_r_2
    casabiome%disturbance_interval = 0.0_r_2
    casabiome%DAMM_EnzPool = 0.0_r_2
    casabiome%DAMM_KMO2 = 0.0_r_2
    casabiome%DAMM_KMcp = 0.0_r_2
    casabiome%DAMM_Ea = 0.0_r_2
    casabiome%DAMM_alpha = 0.0_r_2
    casabiome%plantrate = 0.0_r_2
    casabiome%rmplant = 0.0_r_2
    casabiome%fracnpptoP = 0.0_r_2
    casabiome%fraclignin = 0.0_r_2
    casabiome%fraclabile = 0.0_r_2
    casabiome%ratioNCplantmin = 0.0_r_2
    casabiome%ratioNCplantmax = 0.0_r_2
    casabiome%ratioNPplantmin = 0.0_r_2
    casabiome%ratioNPplantmax = 0.0_r_2
    casabiome%fracLigninplant = 0.0_r_2
    casabiome%ftransNPtoL = 0.0_r_2
    casabiome%ftransPPtoL = 0.0_r_2
    casabiome%litterrate = 0.0_r_2
    casabiome%ratioPcplantmin = 0.0_r_2
    casabiome%ratioPcplantmax = 0.0_r_2
    casabiome%soilrate = 0.0_r_2

  end subroutine zero_casabiome


  subroutine zero_casapool(casapool)

    use cable_def_types_mod, only: r_2

    implicit none

    type(casa_pool), intent(inout) :: casapool

    casapool%Clabile = 0.0_r_2
    casapool%dClabiledt = 0.0_r_2
    casapool%Ctot = 0.0_r_2
    casapool%Ctot_0 = 0.0_r_2
    casapool%Cplant = 0.0_r_2
    casapool%Nplant = 0.0_r_2
    casapool%Pplant = 0.0_r_2
    casapool%dCplantdt = 0.0_r_2
    casapool%dNplantdt = 0.0_r_2
    casapool%dPplantdt = 0.0_r_2
    casapool%ratioNCplant = 0.0_r_2
    casapool%ratioNPplant = 0.0_r_2
    casapool%Nsoilmin = 0.0_r_2
    casapool%Psoillab = 0.0_r_2
    casapool%Psoilsorb = 0.0_r_2
    casapool%Psoilocc = 0.0_r_2
    casapool%dNsoilmindt = 0.0_r_2
    casapool%dPsoillabdt = 0.0_r_2
    casapool%dPsoilsorbdt = 0.0_r_2
    casapool%dPsoiloccdt = 0.0_r_2
    casapool%Clitter = 0.0_r_2
    casapool%Nlitter = 0.0_r_2
    casapool%Plitter = 0.0_r_2
    casapool%dClitterdt = 0.0_r_2
    casapool%dNlitterdt = 0.0_r_2
    casapool%dPlitterdt = 0.0_r_2
    casapool%ratioNClitter = 0.0_r_2
    casapool%ratioNPlitter = 0.0_r_2
    casapool%Csoil = 0.0_r_2
    casapool%Nsoil = 0.0_r_2
    casapool%Psoil = 0.0_r_2
    casapool%dCsoildt = 0.0_r_2
    casapool%dNsoildt = 0.0_r_2
    casapool%dPsoildt = 0.0_r_2
    casapool%ratioNCsoil = 0.0_r_2
    casapool%ratioNCsoilnew = 0.0_r_2
    casapool%ratioNPsoil = 0.0_r_2
    casapool%ratioNCsoilmin = 0.0_r_2
    casapool%ratioNCsoilmax = 0.0_r_2
    casapool%ratioPCsoil = 0.0_r_2
    casapool%ratioPCplant = 0.0_r_2
    casapool%ratioPClitter = 0.0_r_2

  end subroutine zero_casapool


  subroutine zero_casaflux(casaflux)

    use cable_def_types_mod, only: r_2

    implicit none

    type(casa_flux), intent(inout) :: casaflux

    casaflux%Cgpp = 0.0_r_2
    casaflux%Cnpp = 0.0_r_2
    casaflux%Crp = 0.0_r_2
    casaflux%Crgplant = 0.0_r_2
    casaflux%Nminfix = 0.0_r_2
    casaflux%Nminuptake = 0.0_r_2
    casaflux%Plabuptake = 0.0_r_2
    casaflux%Clabloss = 0.0_r_2
    casaflux%fracClabile = 0.0_r_2
    casaflux%stemnpp = 0.0_r_2
    casaflux%frac_sapwood = 0.0_r_2
    casaflux%sapwood_area = 0.0_r_2
    casaflux%Charvest = 0.0_r_2
    casaflux%Nharvest = 0.0_r_2
    casaflux%Pharvest = 0.0_r_2
    casaflux%fHarvest = 0.0_r_2
    casaflux%fcrop = 0.0_r_2
    casaflux%fracCalloc = 0.0_r_2
    casaflux%fracNalloc = 0.0_r_2
    casaflux%fracPalloc = 0.0_r_2
    casaflux%Crmplant = 0.0_r_2
    casaflux%kplant = 0.0_r_2
    casaflux%Cplant_turnover = 0.0_r_2
    casaflux%fromPtoL = 0.0_r_2
    casaflux%Cnep = 0.0_r_2
    casaflux%Crsoil = 0.0_r_2
    casaflux%Nmindep = 0.0_r_2
    casaflux%Nminloss = 0.0_r_2
    casaflux%Nminleach = 0.0_r_2
    casaflux%Nupland = 0.0_r_2
    casaflux%Nlittermin = 0.0_r_2
    casaflux%Nsmin = 0.0_r_2
    casaflux%Nsimm = 0.0_r_2
    casaflux%Nsnet = 0.0_r_2
    casaflux%fNminloss = 0.0_r_2
    casaflux%fNminleach = 0.0_r_2
    casaflux%Pdep = 0.0_r_2
    casaflux%Pwea = 0.0_r_2
    casaflux%Pleach = 0.0_r_2
    casaflux%Ploss = 0.0_r_2
    casaflux%Pupland = 0.0_r_2
    casaflux%Plittermin = 0.0_r_2
    casaflux%Psmin = 0.0_r_2
    casaflux%Psimm = 0.0_r_2
    casaflux%Psnet = 0.0_r_2
    casaflux%fPleach = 0.0_r_2
    casaflux%kplab = 0.0_r_2
    casaflux%kpsorb = 0.0_r_2
    casaflux%kpocc = 0.0_r_2
    casaflux%kmlabp = 0.0_r_2
    casaflux%Psorbmax = 0.0_r_2
    casaflux%Cplant_turnover_disturbance = 0.0_r_2
    casaflux%Cplant_turnover_crowding  = 0.0_r_2
    casaflux%Cplant_turnover_resource_limitation = 0.0_r_2
    casaflux%klitter = 0.0_r_2
    casaflux%ksoil = 0.0_r_2
    casaflux%fromLtoS = 0.0_r_2
    casaflux%fromStoS = 0.0_r_2
    casaflux%fromLtoCO2 = 0.0_r_2
    casaflux%fromStoCO2 = 0.0_r_2
    casaflux%FluxCtolitter = 0.0_r_2
    casaflux%FluxNtolitter = 0.0_r_2
    casaflux%FluxPtolitter = 0.0_r_2
    casaflux%FluxCtosoil = 0.0_r_2
    casaflux%FluxNtosoil = 0.0_r_2
    casaflux%FluxPtosoil = 0.0_r_2
    casaflux%FluxCtoCO2 = 0.0_r_2
    casaflux%FluxCtohwp = 0.0_r_2
    casaflux%FluxNtohwp = 0.0_r_2
    casaflux%FluxPtohwp = 0.0_r_2
    casaflux%FluxCtoclear = 0.0_r_2
    casaflux%FluxNtoclear = 0.0_r_2
    casaflux%FluxPtoclear = 0.0_r_2
    casaflux%CtransferLUC = 0.0_r_2
    casaflux%fromPtoL_fire = 0.0_r_2
    casaflux%klitter_fire = 0.0_r_2
    casaflux%klitter_tot = 0.0_r_2
    casaflux%kplant_fire = 0.0_r_2
    casaflux%kplant_tot = 0.0_r_2
    casaflux%fluxCtoCO2_plant_fire = 0.0_r_2
    casaflux%fluxCtoCO2_litter_fire = 0.0_r_2
    casaflux%fluxfromPtoCO2_fire = 0.0_r_2
    casaflux%fluxfromLtoCO2_fire = 0.0_r_2
    casaflux%fluxNtoAtm_fire = 0.0_r_2
    ! casaflux%fire_mortality_vs_height = 0.0_r_2
    casaflux%FluxFromPtoL = 0.0_r_2
    casaflux%FluxFromLtoS = 0.0_r_2
    casaflux%FluxFromStoS = 0.0_r_2
    casaflux%FluxFromPtoCO2 = 0.0_r_2
    casaflux%FluxFromLtoCO2 = 0.0_r_2
    casaflux%FluxFromStoCO2 = 0.0_r_2
    casaflux%FluxFromPtoHarvest = 0.0_r_2

  end subroutine zero_casaflux


  subroutine zero_casamet(casamet)

    use cable_def_types_mod, only: r_2

    implicit none

    type(casa_met), intent(inout) :: casamet

    casamet%glai = 0.0_r_2
    casamet%Tairk = 0.0_r_2
    casamet%precip = 0.0_r_2
    casamet%tsoilavg = 0.0_r_2
    casamet%moistavg = 0.0_r_2
    casamet%btran = 0.0_r_2
    casamet%lnonwood = 0
    casamet%Tsoil = 0.0_r_2
    casamet%moist = 0.0_r_2
    casamet%iveg2 = 0
    casamet%ijgcm = 0
    casamet%isorder = 0
    casamet%lat = 0.0_r_2
    casamet%lon = 0.0_r_2
    casamet%areacell = 0.0_r_2
    casamet%Tairkspin = 0.0_r_2
    casamet%cgppspin = 0.0_r_2
    casamet%crmplantspin_1 = 0.0_r_2
    casamet%crmplantspin_2 = 0.0_r_2
    casamet%crmplantspin_3 = 0.0_r_2
    casamet%Tsoilspin_1 = 0.0_r_2
    casamet%Tsoilspin_2 = 0.0_r_2
    casamet%Tsoilspin_3 = 0.0_r_2
    casamet%Tsoilspin_4 = 0.0_r_2
    casamet%Tsoilspin_5 = 0.0_r_2
    casamet%Tsoilspin_6 = 0.0_r_2
    casamet%moistspin_1 = 0.0_r_2
    casamet%moistspin_2 = 0.0_r_2
    casamet%moistspin_3 = 0.0_r_2
    casamet%moistspin_4 = 0.0_r_2
    casamet%moistspin_5 = 0.0_r_2
    casamet%moistspin_6 = 0.0_r_2
    casamet%mtempspin = 0.0_r_2
    casamet%frecspin = 0.0_r_2
    casamet%cAn12spin = 0.0_r_2
    casamet%cAn13spin = 0.0_r_2
    casamet%dprecip_spin = 0.0_r_2
    casamet%aprecip_av20_spin = 0.0_r_2
    casamet%du10_max_spin     = 0.0_r_2
    casamet%drhum_spin        = 0.0_r_2
    casamet%dtemp_max_spin    = 0.0_r_2
    casamet%dtemp_min_spin    = 0.0_r_2
    casamet%KBDI_spin         = 0.0_r_2
    casamet%D_MacArthur_spin  = 0.0_r_2
    casamet%FFDI_spin         = 0.0_r_2
    casamet%last_precip_spin  = 0.0_r_2
    casamet%DSLR_spin = 0

  end subroutine zero_casamet


  subroutine zero_casabal(casabal)

    use cable_def_types_mod, only: r_2

    implicit none

    type(casa_balance), intent(inout) :: casabal

    casabal%FCgppyear = 0.0_r_2
    casabal%FCnppyear = 0.0_r_2
    casabal%FCrmleafyear = 0.0_r_2
    casabal%FCrmwoodyear = 0.0_r_2
    casabal%FCrmrootyear = 0.0_r_2
    casabal%FCrgrowyear = 0.0_r_2
    casabal%FCrpyear = 0.0_r_2
    casabal%FCrsyear = 0.0_r_2
    casabal%FCneeyear = 0.0_r_2
    casabal%dCdtyear = 0.0_r_2
    casabal%LAImax = 0.0_r_2
    casabal%Cleafmean = 0.0_r_2
    casabal%Crootmean = 0.0_r_2
    casabal%FNdepyear = 0.0_r_2
    casabal%FNfixyear = 0.0_r_2
    casabal%FNsnetyear = 0.0_r_2
    casabal%FNupyear = 0.0_r_2
    casabal%FNleachyear = 0.0_r_2
    casabal%FNlossyear = 0.0_r_2
    casabal%FPweayear = 0.0_r_2
    casabal%FPdustyear = 0.0_r_2
    casabal%FPsnetyear = 0.0_r_2
    casabal%FPupyear = 0.0_r_2
    casabal%FPleachyear = 0.0_r_2
    casabal%FPlossyear = 0.0_r_2
    casabal%glaimon = 0.0_r_2
    casabal%glaimonx = 0.0_r_2
    casabal%cplantlast = 0.0_r_2
    casabal%nplantlast = 0.0_r_2
    casabal%pplantlast = 0.0_r_2
    casabal%clitterlast = 0.0_r_2
    casabal%nlitterlast = 0.0_r_2
    casabal%plitterlast = 0.0_r_2
    casabal%csoillast = 0.0_r_2
    casabal%nsoillast = 0.0_r_2
    casabal%psoillast = 0.0_r_2
    casabal%nsoilminlast = 0.0_r_2
    casabal%psoillablast = 0.0_r_2
    casabal%psoilsorblast = 0.0_r_2
    casabal%psoilocclast = 0.0_r_2
    casabal%cbalance = 0.0_r_2
    casabal%nbalance = 0.0_r_2
    casabal%pbalance = 0.0_r_2
    casabal%sumcbal = 0.0_r_2
    casabal%sumnbal = 0.0_r_2
    casabal%sumpbal = 0.0_r_2
    casabal%clabilelast = 0.0_r_2

  end subroutine zero_casabal

  ! ------------------------------------------------------------------

  subroutine print_casabiome(casabiome)

    implicit none

    type(casa_biome), intent(in) :: casabiome

    write(*,*) 'ivt2 ', casabiome%ivt2
    write(*,*) 'xkleafcoldmax ', casabiome%xkleafcoldmax
    write(*,*) 'xkleafcoldexp ', casabiome%xkleafcoldexp
    write(*,*) 'xkleafdrymax ', casabiome%xkleafdrymax
    write(*,*) 'xkleafdryexp ', casabiome%xkleafdryexp
    write(*,*) 'glaimax ', casabiome%glaimax
    write(*,*) 'glaimin ', casabiome%glaimin
    write(*,*) 'sla ', casabiome%sla
    write(*,*) 'ratiofrootleaf ', casabiome%ratiofrootleaf
    write(*,*) 'kroot ', casabiome%kroot
    write(*,*) 'krootlen ', casabiome%krootlen
    write(*,*) 'rootdepth ', casabiome%rootdepth
    write(*,*) 'kuptake ', casabiome%kuptake
    write(*,*) 'kminN ', casabiome%kminN
    write(*,*) 'KuplabP ', casabiome%KuplabP
    write(*,*) 'kclabrate ', casabiome%kclabrate
    write(*,*) 'xnpmax ', casabiome%xnpmax
    write(*,*) 'q10soil ', casabiome%q10soil
    write(*,*) 'xkoptlitter ', casabiome%xkoptlitter
    write(*,*) 'xkoptsoil ', casabiome%xkoptsoil
    write(*,*) 'xkplab ', casabiome%xkplab
    write(*,*) 'xkpsorb ', casabiome%xkpsorb
    write(*,*) 'xkpocc ', casabiome%xkpocc
    write(*,*) 'prodptase ', casabiome%prodptase
    write(*,*) 'costnpup ', casabiome%costnpup
    write(*,*) 'maxfinelitter ', casabiome%maxfinelitter
    write(*,*) 'maxcwd ', casabiome%maxcwd
    write(*,*) 'nintercept ', casabiome%nintercept
    write(*,*) 'nslope ', casabiome%nslope
    write(*,*) 'plantrate ', casabiome%plantrate
    write(*,*) 'rmplant ', casabiome%rmplant
    write(*,*) 'fracnpptoP ', casabiome%fracnpptoP
    write(*,*) 'fraclignin ', casabiome%fraclignin
    write(*,*) 'fraclabile ', casabiome%fraclabile
    write(*,*) 'ratioNCplantmin ', casabiome%ratioNCplantmin
    write(*,*) 'ratioNCplantmax ', casabiome%ratioNCplantmax
    write(*,*) 'ratioNPplantmin ', casabiome%ratioNPplantmin
    write(*,*) 'ratioNPplantmax ', casabiome%ratioNPplantmax
    write(*,*) 'fracLigninplant ', casabiome%fracLigninplant
    write(*,*) 'ftransNPtoL ', casabiome%ftransNPtoL
    write(*,*) 'ftransPPtoL ', casabiome%ftransPPtoL
    write(*,*) 'litterrate ', casabiome%litterrate
    write(*,*) 'soilrate ', casabiome%soilrate
    write(*,*) 'ratioPcplantmax ', casabiome%ratioPcplantmax
    write(*,*) 'ratioPcplantmin ', casabiome%ratioPcplantmin
    write(*,*) 'la_to_sa ', casabiome%la_to_sa
    write(*,*) 'vcmax_scalar ', casabiome%vcmax_scalar
    write(*,*) 'disturbance_interval ', casabiome%disturbance_interval
    write(*,*) 'DAMM_EnzPool ', casabiome%DAMM_EnzPool
    write(*,*) 'DAMM_KMO2 ', casabiome%DAMM_KMO2
    write(*,*) 'DAMM_KMcp ', casabiome%DAMM_KMcp
    write(*,*) 'DAMM_Ea ', casabiome%DAMM_Ea
    write(*,*) 'DAMM_alpha ', casabiome%DAMM_alpha

  end subroutine print_casabiome


  subroutine print_casapool(casapool)

    implicit none

    type(casa_pool), intent(in) :: casapool

    write(*,*) 'Clabile ', casapool%Clabile
    write(*,*) 'dClabiledt ', casapool%dClabiledt
    write(*,*) 'Cplant ', casapool%Cplant
    write(*,*) 'Nplant ', casapool%Nplant
    write(*,*) 'Pplant ', casapool%Pplant
    write(*,*) 'dCplantdt ', casapool%dCplantdt
    write(*,*) 'dNplantdt ', casapool%dNplantdt
    write(*,*) 'dPplantdt ', casapool%dPplantdt
    write(*,*) 'ratioNCplant ', casapool%ratioNCplant
    write(*,*) 'ratioNPplant ', casapool%ratioNPplant
    write(*,*) 'Nsoilmin ', casapool%Nsoilmin
    write(*,*) 'Psoillab ', casapool%Psoillab
    write(*,*) 'Psoilsorb ', casapool%Psoilsorb
    write(*,*) 'Psoilocc ', casapool%Psoilocc
    write(*,*) 'dNsoilmindt ', casapool%dNsoilmindt
    write(*,*) 'dPsoillabdt ', casapool%dPsoillabdt
    write(*,*) 'dPsoilsorbdt ', casapool%dPsoilsorbdt
    write(*,*) 'dPsoiloccdt ', casapool%dPsoiloccdt
    write(*,*) 'Clitter ', casapool%Clitter
    write(*,*) 'Nlitter ', casapool%Nlitter
    write(*,*) 'Plitter ', casapool%Plitter
    write(*,*) 'dClitterdt ', casapool%dClitterdt
    write(*,*) 'dNlitterdt ', casapool%dNlitterdt
    write(*,*) 'dPlitterdt ', casapool%dPlitterdt
    write(*,*) 'ratioNClitter ', casapool%ratioNClitter
    write(*,*) 'ratioNPlitter ', casapool%ratioNPlitter
    write(*,*) 'Csoil ', casapool%Csoil
    write(*,*) 'Nsoil ', casapool%Nsoil
    write(*,*) 'Psoil ', casapool%Psoil
    write(*,*) 'dCsoildt ', casapool%dCsoildt
    write(*,*) 'dNsoildt ', casapool%dNsoildt
    write(*,*) 'dPsoildt ', casapool%dPsoildt
    write(*,*) 'ratioNCsoil ', casapool%ratioNCsoil
    write(*,*) 'ratioNPsoil ', casapool%ratioNPsoil
    write(*,*) 'ratioNCsoilnew ', casapool%ratioNCsoilnew
    write(*,*) 'ratioNCsoilmin ', casapool%ratioNCsoilmin
    write(*,*) 'ratioNCsoilmax ', casapool%ratioNCsoilmax
    write(*,*) 'ratioPCsoil ', casapool%ratioPCsoil
    write(*,*) 'ratioPCplant ', casapool%ratioPCplant
    write(*,*) 'ratioPClitter ', casapool%ratioPClitter
    write(*,*) 'Ctot_0 ', casapool%Ctot_0
    write(*,*) 'Ctot ', casapool%Ctot

  end subroutine print_casapool


  subroutine print_casaflux(casaflux)

    implicit none

    type(casa_flux), intent(in) :: casaflux

    write(*,*) 'Cgpp ', casaflux%Cgpp
    write(*,*) 'Cnpp ', casaflux%Cnpp
    write(*,*) 'Crp ', casaflux%Crp
    write(*,*) 'Crgplant ', casaflux%Crgplant
    write(*,*) 'Nminfix ', casaflux%Nminfix
    write(*,*) 'Nminuptake ', casaflux%Nminuptake
    write(*,*) 'Plabuptake ', casaflux%Plabuptake
    write(*,*) 'Clabloss ', casaflux%Clabloss
    write(*,*) 'fracClabile ', casaflux%fracClabile
    write(*,*) 'fracCalloc ', casaflux%fracCalloc
    write(*,*) 'fracNalloc ', casaflux%fracNalloc
    write(*,*) 'fracPalloc ', casaflux%fracPalloc
    write(*,*) 'kplant ', casaflux%kplant
    write(*,*) 'Crmplant ', casaflux%Crmplant
    write(*,*) 'fromPtoL ', casaflux%fromPtoL
    write(*,*) 'Cnep ', casaflux%Cnep
    write(*,*) 'Crsoil ', casaflux%Crsoil
    write(*,*) 'Nmindep ', casaflux%Nmindep
    write(*,*) 'Nminloss ', casaflux%Nminloss
    write(*,*) 'Nminleach ', casaflux%Nminleach
    write(*,*) 'Nupland ', casaflux%Nupland
    write(*,*) 'Nlittermin ', casaflux%Nlittermin
    write(*,*) 'Nsmin ', casaflux%Nsmin
    write(*,*) 'Nsimm ', casaflux%Nsimm
    write(*,*) 'Nsnet ', casaflux%Nsnet
    write(*,*) 'fNminloss ', casaflux%fNminloss
    write(*,*) 'fNminleach ', casaflux%fNminleach
    write(*,*) 'Pdep ', casaflux%Pdep
    write(*,*) 'Pwea ', casaflux%Pwea
    write(*,*) 'Pleach ', casaflux%Pleach
    write(*,*) 'Ploss ', casaflux%Ploss
    write(*,*) 'Pupland ', casaflux%Pupland
    write(*,*) 'Plittermin ', casaflux%Plittermin
    write(*,*) 'Psmin ', casaflux%Psmin
    write(*,*) 'Psimm ', casaflux%Psimm
    write(*,*) 'Psnet ', casaflux%Psnet
    write(*,*) 'fPleach ', casaflux%fPleach
    write(*,*) 'kplab ', casaflux%kplab
    write(*,*) 'kpsorb ', casaflux%kpsorb
    write(*,*) 'kpocc ', casaflux%kpocc
    write(*,*) 'kmlabP ', casaflux%kmlabP
    write(*,*) 'Psorbmax ', casaflux%Psorbmax
    write(*,*) 'klitter ', casaflux%klitter
    write(*,*) 'ksoil ', casaflux%ksoil
    write(*,*) 'fromLtoS ', casaflux%fromLtoS
    write(*,*) 'fromStoS ', casaflux%fromStoS
    write(*,*) 'fromLtoCO2 ', casaflux%fromLtoCO2
    write(*,*) 'fromStoCO2 ', casaflux%fromStoCO2
    write(*,*) 'stemnpp ', casaflux%stemnpp
    write(*,*) 'frac_sapwood ', casaflux%frac_sapwood
    write(*,*) 'sapwood_area ', casaflux%sapwood_area
    write(*,*) 'fharvest ', casaflux%fharvest
    write(*,*) 'Charvest ', casaflux%Charvest
    write(*,*) 'Nharvest ', casaflux%Nharvest
    write(*,*) 'Pharvest ', casaflux%Pharvest
    write(*,*) 'fcrop ', casaflux%fcrop
    write(*,*) 'Cplant_turnover ', casaflux%Cplant_turnover
    write(*,*) 'Cplant_turnover_disturbance ', casaflux%Cplant_turnover_disturbance
    write(*,*) 'Cplant_turnover_crowding ', casaflux%Cplant_turnover_crowding
    write(*,*) 'Cplant_turnover_resource_limitation ', casaflux%Cplant_turnover_resource_limitation

    write(*,*) 'fromPtoL_fire ', casaflux%fromPtoL_fire
    write(*,*) 'kplant_fire ', casaflux%kplant_fire
    write(*,*) 'klitter_fire ', casaflux%klitter_fire
    write(*,*) 'kplant_tot ', casaflux%kplant_tot
    write(*,*) 'klitter_tot ', casaflux%klitter_tot
    write(*,*) 'FluxCtoCO2_plant_fire ', casaflux%FluxCtoCO2_plant_fire
    write(*,*) 'FluxCtoCO2_litter_fire ', casaflux%FluxCtoCO2_litter_fire
    write(*,*) 'fluxfromPtoCO2_fire ', casaflux%fluxfromPtoCO2_fire
    write(*,*) 'fluxfromLtoCO2_fire ', casaflux%fluxfromLtoCO2_fire
    write(*,*) 'FluxNtoAtm_fire ', casaflux%FluxNtoAtm_fire

    write(*,*) 'FluxCtolitter ', casaflux%FluxCtolitter
    write(*,*) 'FluxNtolitter ', casaflux%FluxNtolitter
    write(*,*) 'FluxPtolitter ', casaflux%FluxPtolitter

    write(*,*) 'FluxCtosoil ', casaflux%FluxCtosoil
    write(*,*) 'FluxNtosoil ', casaflux%FluxNtosoil
    write(*,*) 'FluxPtosoil ', casaflux%FluxPtosoil

    write(*,*) 'FluxCtohwp ', casaflux%FluxCtohwp
    write(*,*) 'FluxNtohwp ', casaflux%FluxNtohwp
    write(*,*) 'FluxPtohwp ', casaflux%FluxPtohwp

    write(*,*) 'FluxCtoclear ', casaflux%FluxCtoclear
    write(*,*) 'FluxNtoclear ', casaflux%FluxNtoclear
    write(*,*) 'FluxPtoclear ', casaflux%FluxPtoclear

    write(*,*) 'CtransferLUC ', casaflux%CtransferLUC

    write(*,*) 'FluxCtoco2 ', casaflux%FluxCtoco2

    write(*,*) 'FluxFromPtoL ', casaflux%FluxFromPtoL
    write(*,*) 'FluxFromLtoS ', casaflux%FluxFromLtoS
    write(*,*) 'FluxFromStoS ', casaflux%FluxFromStoS
    write(*,*) 'FluxFromPtoCO2 ', casaflux%FluxFromPtoCO2
    write(*,*) 'FluxFromLtoCO2 ', casaflux%FluxFromLtoCO2
    write(*,*) 'FluxFromStoCO2 ', casaflux%FluxFromStoCO2
    write(*,*) 'FluxFromPtoHarvest ', casaflux%FluxFromPtoHarvest

  end subroutine print_casaflux


  subroutine print_casamet(casamet)

    implicit none

    type(casa_met), intent(in) :: casamet

    write(*,*) 'glai ', casamet%glai
    write(*,*) 'lnonwood ', casamet%lnonwood
    write(*,*) 'Tairk ', casamet%Tairk
    write(*,*) 'precip ', casamet%precip
    write(*,*) 'tsoilavg ', casamet%tsoilavg
    write(*,*) 'moistavg ', casamet%moistavg
    write(*,*) 'btran ', casamet%btran
    write(*,*) 'Tsoil ', casamet%Tsoil
    write(*,*) 'moist ', casamet%moist
    write(*,*) 'iveg2 ', casamet%iveg2
    write(*,*) 'ijgcm ', casamet%ijgcm
    write(*,*) 'isorder ', casamet%isorder
    write(*,*) 'lat ', casamet%lat
    write(*,*) 'lon ', casamet%lon
    write(*,*) 'areacell ', casamet%areacell
    write(*,*) 'Tairkspin ', casamet%Tairkspin
    write(*,*) 'cgppspin ', casamet%cgppspin
    write(*,*) 'crmplantspin_1 ', casamet%crmplantspin_1
    write(*,*) 'crmplantspin_2 ', casamet%crmplantspin_2
    write(*,*) 'crmplantspin_3 ', casamet%crmplantspin_3
    write(*,*) 'Tsoilspin_1 ', casamet%Tsoilspin_1
    write(*,*) 'Tsoilspin_2 ', casamet%Tsoilspin_2
    write(*,*) 'Tsoilspin_3 ', casamet%Tsoilspin_3
    write(*,*) 'Tsoilspin_4 ', casamet%Tsoilspin_4
    write(*,*) 'Tsoilspin_5 ', casamet%Tsoilspin_5
    write(*,*) 'Tsoilspin_6 ', casamet%Tsoilspin_6
    write(*,*) 'moistspin_1 ', casamet%moistspin_1
    write(*,*) 'moistspin_2 ', casamet%moistspin_2
    write(*,*) 'moistspin_3 ', casamet%moistspin_3
    write(*,*) 'moistspin_4 ', casamet%moistspin_4
    write(*,*) 'moistspin_5 ', casamet%moistspin_5
    write(*,*) 'moistspin_6 ', casamet%moistspin_6
    write(*,*) 'mtempspin ', casamet%mtempspin
    write(*,*) 'cAn12spin ', casamet%cAn12spin
    write(*,*) 'cAn13spin ', casamet%cAn13spin

  end subroutine print_casamet


  subroutine print_casabal(casabal)

    implicit none

    type(casa_balance), intent(in) :: casabal

    write(*,*) 'FCgppyear ', casabal%FCgppyear
    write(*,*) 'FCnppyear ', casabal%FCnppyear
    write(*,*) 'FCrpyear ', casabal%FCrpyear
    write(*,*) 'FCrmleafyear ', casabal%FCrmleafyear
    write(*,*) 'FCrmwoodyear ', casabal%FCrmwoodyear
    write(*,*) 'FCrmrootyear ', casabal%FCrmrootyear
    write(*,*) 'FCrgrowyear ', casabal%FCrgrowyear
    write(*,*) 'FCrsyear ', casabal%FCrsyear
    write(*,*) 'FCneeyear ', casabal%FCneeyear
    write(*,*) 'FNdepyear ', casabal%FNdepyear
    write(*,*) 'FNfixyear ', casabal%FNfixyear
    write(*,*) 'FNsnetyear ', casabal%FNsnetyear
    write(*,*) 'FNupyear ', casabal%FNupyear
    write(*,*) 'FNleachyear ', casabal%FNleachyear
    write(*,*) 'FNlossyear ', casabal%FNlossyear
    write(*,*) 'FPweayear ', casabal%FPweayear
    write(*,*) 'FPdustyear ', casabal%FPdustyear
    write(*,*) 'FPsnetyear ', casabal%FPsnetyear
    write(*,*) 'FPupyear ', casabal%FPupyear
    write(*,*) 'FPleachyear ', casabal%FPleachyear
    write(*,*) 'FPlossyear ', casabal%FPlossyear
    write(*,*) 'dCdtyear ', casabal%dCdtyear
    write(*,*) 'LAImax ', casabal%LAImax
    write(*,*) 'Cleafmean ', casabal%Cleafmean
    write(*,*) 'Crootmean ', casabal%Crootmean

    write(*,*) 'glaimon ', casabal%glaimon
    write(*,*) 'glaimonx ', casabal%glaimonx

    write(*,*) 'cplantlast ', casabal%cplantlast
    write(*,*) 'nplantlast ', casabal%nplantlast
    write(*,*) 'pplantlast ', casabal%pplantlast

    write(*,*) 'clitterlast ', casabal%clitterlast
    write(*,*) 'nlitterlast ', casabal%nlitterlast
    write(*,*) 'plitterlast ', casabal%plitterlast

    write(*,*) 'csoillast ', casabal%csoillast
    write(*,*) 'nsoillast ', casabal%nsoillast
    write(*,*) 'psoillast ', casabal%psoillast

    write(*,*) 'nsoilminlast ', casabal%nsoilminlast
    write(*,*) 'psoillablast ', casabal%psoillablast
    write(*,*) 'psoilsorblast ', casabal%psoilsorblast
    write(*,*) 'psoilocclast ', casabal%psoilocclast
    write(*,*) 'cbalance ', casabal%cbalance
    write(*,*) 'nbalance ', casabal%nbalance
    write(*,*) 'pbalance ', casabal%pbalance
    write(*,*) 'sumcbal ', casabal%sumcbal
    write(*,*) 'sumnbal ', casabal%sumnbal
    write(*,*) 'sumpbal ', casabal%sumpbal
    write(*,*) 'clabilelast ', casabal%clabilelast

  end subroutine print_casabal

  ! ------------------------------------------------------------------

  subroutine zero_sum_casa(sum_casapool, sum_casaflux)

    implicit none

    type(casa_pool), intent(inout) :: sum_casapool
    type(casa_flux), intent(inout) :: sum_casaflux

    call zero_casapool(sum_casapool)
    call zero_casaflux(sum_casaflux)

  end subroutine zero_sum_casa


  subroutine update_sum_casa(sum_casapool, sum_casaflux, casapool, casaflux, sum_now, average_now, nsteps)

    use cable_def_types_mod, only: r_2

    implicit none

    type(casa_pool), intent(inout) :: sum_casapool
    type(casa_flux), intent(inout) :: sum_casaflux
    type(casa_pool), intent(in)    :: casapool
    type(casa_flux), intent(in)    :: casaflux
    logical,         intent(in)    :: sum_now, average_now
    integer,         intent(in)    :: nsteps

    real(r_2) :: rnsteps

    rnsteps = 1.0 / real(nsteps, r_2)

    if (sum_now) then
       sum_casapool%Clabile        = sum_casapool%Clabile        + casapool%Clabile
       sum_casapool%dClabiledt     = sum_casapool%Clabile        + casapool%Clabile
       sum_casapool%Cplant         = sum_casapool%Cplant         + casapool%Cplant
       sum_casapool%Nplant         = sum_casapool%Nplant         + casapool%Nplant
       sum_casapool%Pplant         = sum_casapool%Pplant         + casapool%Pplant
       sum_casapool%dCplantdt      = sum_casapool%dCplantdt      + casapool%dCplantdt
       sum_casapool%dNplantdt      = sum_casapool%dNplantdt      + casapool%dNplantdt
       sum_casapool%dPplantdt      = sum_casapool%dPplantdt      + casapool%dPplantdt
       sum_casapool%ratioNCplant   = sum_casapool%ratioNCplant   + casapool%ratioNCplant
       sum_casapool%ratioNPplant   = sum_casapool%ratioNPplant   + casapool%ratioNPplant
       sum_casapool%Nsoilmin       = sum_casapool%Nsoilmin       + casapool%Nsoilmin
       sum_casapool%Psoillab       = sum_casapool%Psoillab       + casapool%Psoillab
       sum_casapool%Psoilsorb      = sum_casapool%Psoilsorb      + casapool%Psoilsorb
       sum_casapool%Psoilocc       = sum_casapool%Psoilocc       + casapool%Psoilocc
       sum_casapool%dNsoilmindt    = sum_casapool%dNsoilmindt    + casapool%dNsoilmindt
       sum_casapool%dPsoillabdt    = sum_casapool%dPsoillabdt    + casapool%dPsoillabdt
       sum_casapool%dPsoilsorbdt   = sum_casapool%dPsoilsorbdt   + casapool%dPsoilsorbdt
       sum_casapool%dPsoiloccdt    = sum_casapool%dPsoiloccdt    + casapool%dPsoiloccdt
       sum_casapool%Clitter        = sum_casapool%Clitter        + casapool%Clitter
       sum_casapool%Nlitter        = sum_casapool%Nlitter        + casapool%Nlitter
       sum_casapool%Plitter        = sum_casapool%Plitter        + casapool%Plitter
       sum_casapool%dClitterdt     = sum_casapool%dClitterdt     + casapool%dClitterdt
       sum_casapool%dNlitterdt     = sum_casapool%dNlitterdt     + casapool%dNlitterdt
       sum_casapool%dPlitterdt     = sum_casapool%dPlitterdt     + casapool%dPlitterdt
       sum_casapool%ratioNClitter  = sum_casapool%ratioNClitter  + casapool%ratioNClitter
       sum_casapool%ratioNPlitter  = sum_casapool%ratioNPlitter  + casapool%ratioNPlitter
       sum_casapool%Csoil          = sum_casapool%Csoil          + casapool%Csoil
       sum_casapool%Nsoil          = sum_casapool%Nsoil          + casapool%Nsoil
       sum_casapool%Psoil          = sum_casapool%Psoil          + casapool%Psoil
       sum_casapool%dCsoildt       = sum_casapool%dCsoildt       + casapool%dCsoildt
       sum_casapool%dNsoildt       = sum_casapool%dNsoildt       + casapool%dNsoildt
       sum_casapool%dPsoildt       = sum_casapool%dPsoildt       + casapool%dPsoildt
       sum_casapool%ratioNCsoil    = sum_casapool%ratioNCsoil    + casapool%ratioNCsoil
       sum_casapool%ratioNPsoil    = sum_casapool%ratioNPsoil    + casapool%ratioNPsoil
       sum_casapool%ratioNCsoilnew = sum_casapool%ratioNCsoilnew + casapool%ratioNCsoilnew
       sum_casapool%ratioNCsoilmin = sum_casapool%ratioNCsoilmin + casapool%ratioNCsoilmin
       sum_casapool%ratioNCsoilmax = sum_casapool%ratioNCsoilmax + casapool%ratioNCsoilmax
       sum_casapool%ratioPcsoil    = sum_casapool%ratioPcsoil    + casapool%ratioPcsoil
       sum_casapool%ratioPcplant   = sum_casapool%ratioPcplant   + casapool%ratioPcplant
       sum_casapool%ratioPclitter  = sum_casapool%ratioPclitter  + casapool%ratioPclitter

       sum_casaflux%Cgpp            = sum_casaflux%Cgpp            + casaflux%Cgpp
       sum_casaflux%Cnpp            = sum_casaflux%Cnpp            + casaflux%Cnpp
       sum_casaflux%Crp             = sum_casaflux%Crp             + casaflux%Crp
       sum_casaflux%Crgplant        = sum_casaflux%Crgplant        + casaflux%Crgplant
       sum_casaflux%Nminfix         = sum_casaflux%Nminfix         + casaflux%Nminfix
       sum_casaflux%Nminuptake      = sum_casaflux%Nminuptake      + casaflux%Nminuptake
       sum_casaflux%Plabuptake      = sum_casaflux%Plabuptake      + casaflux%Plabuptake
       sum_casaflux%Clabloss        = sum_casaflux%Clabloss        + casaflux%Clabloss
       sum_casaflux%fracClabile     = sum_casaflux%fracClabile     + casaflux%fracClabile
       ! sum_casaflux%fracCalloc(:,1) = sum_casaflux%fracCalloc(:,1) + casaflux%fracCalloc(:,1)
       ! sum_casaflux%fracCalloc(:,2) = sum_casaflux%fracCalloc(:,2) + casaflux%fracCalloc(:,2)
       ! sum_casaflux%fracCalloc(:,3) = sum_casaflux%fracCalloc(:,3) + casaflux%fracCalloc(:,3)
       sum_casaflux%fracCalloc(:,1) = sum_casaflux%fracCalloc(:,1) + casaflux%fracCalloc(:,1) * casaflux%Cnpp
       sum_casaflux%fracCalloc(:,2) = sum_casaflux%fracCalloc(:,2) + casaflux%fracCalloc(:,2) * casaflux%Cnpp
       sum_casaflux%fracCalloc(:,3) = sum_casaflux%fracCalloc(:,3) + casaflux%fracCalloc(:,3) * casaflux%Cnpp
       sum_casaflux%fracNalloc      = sum_casaflux%fracNalloc      + casaflux%fracNalloc
       sum_casaflux%fracPalloc      = sum_casaflux%fracPalloc      + casaflux%fracPalloc
       sum_casaflux%Crmplant        = sum_casaflux%Crmplant        + casaflux%Crmplant
       ! sum_casaflux%kplant          = sum_casaflux%kplant          + casaflux%kplant
       sum_casaflux%kplant          = sum_casaflux%kplant          + casaflux%kplant * casapool%Cplant
       sum_casaflux%fromPtoL        = sum_casaflux%fromPtoL        + casaflux%fromPtoL
       sum_casaflux%Cnep            = sum_casaflux%Cnep            + casaflux%Cnep
       sum_casaflux%Crsoil          = sum_casaflux%Crsoil          + casaflux%Crsoil
       sum_casaflux%Nmindep         = sum_casaflux%Nmindep         + casaflux%Nmindep
       sum_casaflux%Nminloss        = sum_casaflux%Nminloss        + casaflux%Nminloss
       sum_casaflux%Nminleach       = sum_casaflux%Nminleach       + casaflux%Nminleach
       sum_casaflux%Nupland         = sum_casaflux%Nupland         + casaflux%Nupland
       sum_casaflux%Nlittermin      = sum_casaflux%Nlittermin      + casaflux%Nlittermin
       sum_casaflux%Nsmin           = sum_casaflux%Nsmin           + casaflux%Nsmin
       sum_casaflux%Nsimm           = sum_casaflux%Nsimm           + casaflux%Nsimm
       sum_casaflux%Nsnet           = sum_casaflux%Nsnet           + casaflux%Nsnet
       sum_casaflux%fNminloss       = sum_casaflux%fNminloss       + casaflux%fNminloss
       sum_casaflux%fNminleach      = sum_casaflux%fNminleach      + casaflux%fNminleach
       sum_casaflux%Pdep            = sum_casaflux%Pdep            + casaflux%Pdep
       sum_casaflux%Pwea            = sum_casaflux%Pwea            + casaflux%Pwea
       sum_casaflux%Pleach          = sum_casaflux%Pleach          + casaflux%Pleach
       sum_casaflux%Ploss           = sum_casaflux%Ploss           + casaflux%Ploss
       sum_casaflux%Pupland         = sum_casaflux%Pupland         + casaflux%Pupland
       sum_casaflux%Plittermin      = sum_casaflux%Plittermin      + casaflux%Plittermin
       sum_casaflux%Psmin           = sum_casaflux%Psmin           + casaflux%Psmin
       sum_casaflux%Psimm           = sum_casaflux%Psimm           + casaflux%Psimm
       sum_casaflux%Psnet           = sum_casaflux%Psnet           + casaflux%Psnet
       sum_casaflux%fPleach         = sum_casaflux%fPleach         + casaflux%fPleach
       sum_casaflux%kplab           = sum_casaflux%kplab           + casaflux%kplab
       sum_casaflux%kpsorb          = sum_casaflux%kpsorb          + casaflux%kpsorb
       sum_casaflux%kpocc           = sum_casaflux%kpocc           + casaflux%kpocc
       sum_casaflux%kmlabP          = sum_casaflux%kmlabP          + casaflux%kmlabP
       sum_casaflux%Psorbmax        = sum_casaflux%Psorbmax        + casaflux%Psorbmax
       sum_casaflux%klitter         = sum_casaflux%klitter         + casaflux%klitter
       sum_casaflux%ksoil           = sum_casaflux%ksoil           + casaflux%ksoil
       sum_casaflux%fromLtoS        = sum_casaflux%fromLtoS        + casaflux%fromLtoS
       sum_casaflux%fromStoS        = sum_casaflux%fromStoS        + casaflux%fromStoS
       sum_casaflux%fromLtoCO2      = sum_casaflux%fromLtoCO2      + casaflux%fromLtoCO2
       sum_casaflux%fromStoCO2      = sum_casaflux%fromStoCO2      + casaflux%fromStoCO2
       sum_casaflux%stemnpp         = sum_casaflux%stemnpp         + casaflux%stemnpp
       sum_casaflux%frac_sapwood    = sum_casaflux%frac_sapwood    + casaflux%frac_sapwood
       sum_casaflux%sapwood_area    = sum_casaflux%sapwood_area    + casaflux%sapwood_area
       sum_casaflux%Cplant_turnover = &
            sum_casaflux%Cplant_turnover + casaflux%Cplant_turnover
       sum_casaflux%Cplant_turnover_disturbance = &
            sum_casaflux%Cplant_turnover_disturbance + casaflux%Cplant_turnover_disturbance
       sum_casaflux%Cplant_turnover_crowding = &
            sum_casaflux%Cplant_turnover_crowding + casaflux%Cplant_turnover_crowding
       sum_casaflux%Cplant_turnover_resource_limitation = &
            sum_casaflux%Cplant_turnover_resource_limitation + casaflux%Cplant_turnover_resource_limitation

       sum_casaflux%FluxCtolitter = sum_casaflux%FluxCtolitter + casaflux%FluxCtolitter
       sum_casaflux%FluxNtolitter = sum_casaflux%FluxNtolitter + casaflux%FluxNtolitter
       sum_casaflux%FluxPtolitter = sum_casaflux%FluxPtolitter + casaflux%FluxPtolitter

       sum_casaflux%FluxCtosoil = sum_casaflux%FluxCtosoil + casaflux%FluxCtosoil
       sum_casaflux%FluxNtosoil = sum_casaflux%FluxNtosoil + casaflux%FluxNtosoil
       sum_casaflux%FluxPtosoil = sum_casaflux%FluxPtosoil + casaflux%FluxPtosoil

       sum_casaflux%FluxCtoco2 = sum_casaflux%FluxCtoco2 + casaflux%FluxCtoco2

       sum_casaflux%fromPtoL_fire = sum_casaflux%fromPtoL_fire + casaflux%fromPtoL_fire
       sum_casaflux%kplant_fire   = sum_casaflux%kplant_fire   + casaflux%kplant_fire
       sum_casaflux%klitter_fire  = sum_casaflux%klitter_fire  + casaflux%klitter_fire
       sum_casaflux%kplant_tot    = sum_casaflux%kplant_tot    + casaflux%kplant_tot
       sum_casaflux%klitter_tot   = sum_casaflux%klitter_tot   + casaflux%klitter_tot
       sum_casaflux%FluxCtoCO2_plant_fire  = sum_casaflux%FluxCtoCO2_plant_fire + &
            casaflux%FluxCtoCO2_plant_fire
       sum_casaflux%FluxCtoCO2_litter_fire = sum_casaflux%FluxCtoCO2_litter_fire + &
            casaflux%FluxCtoCO2_litter_fire
       sum_casaflux%fluxfromPtoCO2_fire    = sum_casaflux%fluxfromPtoCO2_fire + &
            casaflux%fluxfromPtoCO2_fire
       sum_casaflux%fluxfromLtoCO2_fire    = sum_casaflux%fluxfromLtoCO2_fire + &
            casaflux%fluxfromLtoCO2_fire

       sum_casaflux%FluxNtoAtm_fire = sum_casaflux%FluxNtoAtm_fire + &
            casaflux%FluxNtoAtm_fire

       sum_casaflux%FluxFromPtoL   = sum_casaflux%FluxFromPtoL   + casaflux%FluxFromPtoL
       sum_casaflux%FluxFromLtoS   = sum_casaflux%FluxFromLtoS   + casaflux%FluxFromLtoS
       sum_casaflux%FluxFromStoS   = sum_casaflux%FluxFromStoS   + casaflux%FluxFromStoS
       sum_casaflux%FluxFromPtoCO2 = sum_casaflux%FluxFromPtoCO2 + casaflux%FluxFromPtoCO2
       sum_casaflux%FluxFromLtoCO2 = sum_casaflux%FluxFromLtoCO2 + casaflux%FluxFromLtoCO2
       sum_casaflux%FluxFromStoCO2 = sum_casaflux%FluxFromStoCO2 + casaflux%FluxFromStoCO2
       sum_casaflux%FluxFromPtoHarvest = sum_casaflux%FluxFromPtoHarvest + &
            casaflux%FluxFromPtoHarvest
    endif ! sum_now

    if (average_now) then
       ! sum_casaflux%fracCalloc = sum_casaflux%fracCalloc * rnsteps
       where (sum_casaflux%Cnpp > 1.e-12_r_2)
          sum_casaflux%fracCalloc(:,1) = sum_casaflux%fracCalloc(:,1) / sum_casaflux%Cnpp
          sum_casaflux%fracCalloc(:,2) = sum_casaflux%fracCalloc(:,2) / sum_casaflux%Cnpp
          sum_casaflux%fracCalloc(:,3) = sum_casaflux%fracCalloc(:,3) / sum_casaflux%Cnpp
       elsewhere
          sum_casaflux%fracCalloc(:,1) = 0.0_r_2
          sum_casaflux%fracCalloc(:,2) = 0.0_r_2
          sum_casaflux%fracCalloc(:,3) = 0.0_r_2
       endwhere
       ! sum_casaflux%kplant = sum_casaflux%kplant * rnsteps
       where (sum_casapool%Cplant > 1.e-12_r_2)
          sum_casaflux%kplant = sum_casaflux%kplant / sum_casapool%Cplant
       elsewhere
          sum_casaflux%kplant = 0.0_r_2
       endwhere
       sum_casapool%Clabile    = sum_casapool%Clabile * rnsteps
       sum_casapool%dClabiledt = sum_casapool%Clabile * rnsteps
       sum_casapool%Cplant         = sum_casapool%Cplant         * rnsteps
       sum_casapool%Nplant         = sum_casapool%Nplant         * rnsteps
       sum_casapool%Pplant         = sum_casapool%Pplant         * rnsteps
       sum_casapool%dCplantdt      = sum_casapool%dCplantdt      * rnsteps
       sum_casapool%dNplantdt      = sum_casapool%dNplantdt      * rnsteps
       sum_casapool%dPplantdt      = sum_casapool%dPplantdt      * rnsteps
       sum_casapool%ratioNCplant   = sum_casapool%ratioNCplant   * rnsteps
       sum_casapool%ratioNPplant   = sum_casapool%ratioNPplant   * rnsteps
       sum_casapool%Nsoilmin       = sum_casapool%Nsoilmin       * rnsteps
       sum_casapool%Psoillab       = sum_casapool%Psoillab       * rnsteps
       sum_casapool%Psoilsorb      = sum_casapool%Psoilsorb      * rnsteps
       sum_casapool%Psoilocc       = sum_casapool%Psoilocc       * rnsteps
       sum_casapool%dNsoilmindt    = sum_casapool%dNsoilmindt    * rnsteps
       sum_casapool%dPsoillabdt    = sum_casapool%dPsoillabdt    * rnsteps
       sum_casapool%dPsoilsorbdt   = sum_casapool%dPsoilsorbdt   * rnsteps
       sum_casapool%dPsoiloccdt    = sum_casapool%dPsoiloccdt    * rnsteps
       sum_casapool%Clitter        = sum_casapool%Clitter        * rnsteps
       sum_casapool%Nlitter        = sum_casapool%Nlitter        * rnsteps
       sum_casapool%Plitter        = sum_casapool%Plitter        * rnsteps
       sum_casapool%dClitterdt     = sum_casapool%dClitterdt     * rnsteps
       sum_casapool%dNlitterdt     = sum_casapool%dNlitterdt     * rnsteps
       sum_casapool%dPlitterdt     = sum_casapool%dPlitterdt     * rnsteps
       sum_casapool%ratioNClitter  = sum_casapool%ratioNClitter  * rnsteps
       sum_casapool%ratioNPlitter  = sum_casapool%ratioNPlitter  * rnsteps
       sum_casapool%Csoil          = sum_casapool%Csoil          * rnsteps
       sum_casapool%Nsoil          = sum_casapool%Nsoil          * rnsteps
       sum_casapool%Psoil          = sum_casapool%Psoil          * rnsteps
       sum_casapool%dCsoildt       = sum_casapool%dCsoildt       * rnsteps
       sum_casapool%dNsoildt       = sum_casapool%dNsoildt       * rnsteps
       sum_casapool%dPsoildt       = sum_casapool%dPsoildt       * rnsteps
       sum_casapool%ratioNCsoil    = sum_casapool%ratioNCsoil    * rnsteps
       sum_casapool%ratioNPsoil    = sum_casapool%ratioNPsoil    * rnsteps
       sum_casapool%ratioNCsoilnew = sum_casapool%ratioNCsoilnew * rnsteps
       sum_casapool%ratioNCsoilmin = sum_casapool%ratioNCsoilmin * rnsteps
       sum_casapool%ratioNCsoilmax = sum_casapool%ratioNCsoilmax * rnsteps
       sum_casapool%ratioPCsoil    = sum_casapool%ratioPCsoil    * rnsteps
       sum_casapool%ratioPCplant   = sum_casapool%ratioPCplant   * rnsteps
       sum_casapool%ratioPClitter  = sum_casapool%ratioPClitter  * rnsteps

       sum_casaflux%Cgpp        = sum_casaflux%Cgpp        * rnsteps
       sum_casaflux%Cnpp        = sum_casaflux%Cnpp        * rnsteps
       sum_casaflux%Crp         = sum_casaflux%Crp         * rnsteps
       sum_casaflux%Crgplant    = sum_casaflux%Crgplant    * rnsteps
       sum_casaflux%Nminfix     = sum_casaflux%Nminfix     * rnsteps
       sum_casaflux%Nminuptake  = sum_casaflux%Nminuptake  * rnsteps
       sum_casaflux%Plabuptake  = sum_casaflux%Plabuptake  * rnsteps
       sum_casaflux%Clabloss    = sum_casaflux%Clabloss    * rnsteps
       sum_casaflux%fracClabile = sum_casaflux%fracClabile * rnsteps
       sum_casaflux%fracNalloc = sum_casaflux%fracNalloc * rnsteps
       sum_casaflux%fracPalloc = sum_casaflux%fracPalloc * rnsteps

       sum_casaflux%Crmplant     = sum_casaflux%Crmplant     * rnsteps
       sum_casaflux%fromPtoL     = sum_casaflux%fromPtoL     * rnsteps
       sum_casaflux%Cnep         = sum_casaflux%Cnep         * rnsteps
       sum_casaflux%Crsoil       = sum_casaflux%Crsoil       * rnsteps
       sum_casaflux%Nmindep      = sum_casaflux%Nmindep      * rnsteps
       sum_casaflux%Nminloss     = sum_casaflux%Nminloss     * rnsteps
       sum_casaflux%Nminleach    = sum_casaflux%Nminleach    * rnsteps
       sum_casaflux%Nupland      = sum_casaflux%Nupland      * rnsteps
       sum_casaflux%Nlittermin   = sum_casaflux%Nlittermin   * rnsteps
       sum_casaflux%Nsmin        = sum_casaflux%Nsmin        * rnsteps
       sum_casaflux%Nsimm        = sum_casaflux%Nsimm        * rnsteps
       sum_casaflux%Nsnet        = sum_casaflux%Nsnet        * rnsteps
       sum_casaflux%fNminloss    = sum_casaflux%fNminloss    * rnsteps
       sum_casaflux%fNminleach   = sum_casaflux%fNminleach   * rnsteps
       sum_casaflux%Pdep         = sum_casaflux%Pdep         * rnsteps
       sum_casaflux%Pwea         = sum_casaflux%Pwea         * rnsteps
       sum_casaflux%Pleach       = sum_casaflux%Pleach       * rnsteps
       sum_casaflux%Ploss        = sum_casaflux%Ploss        * rnsteps
       sum_casaflux%Pupland      = sum_casaflux%Pupland      * rnsteps
       sum_casaflux%Plittermin   = sum_casaflux%Plittermin   * rnsteps
       sum_casaflux%Psmin        = sum_casaflux%Psmin        * rnsteps
       sum_casaflux%Psimm        = sum_casaflux%Psimm        * rnsteps
       sum_casaflux%Psnet        = sum_casaflux%Psnet        * rnsteps
       sum_casaflux%fPleach      = sum_casaflux%fPleach      * rnsteps
       sum_casaflux%kplab        = sum_casaflux%kplab        * rnsteps
       sum_casaflux%kpsorb       = sum_casaflux%kpsorb       * rnsteps
       sum_casaflux%kpocc        = sum_casaflux%kpocc        * rnsteps
       sum_casaflux%kmlabP       = sum_casaflux%kmlabP       * rnsteps
       sum_casaflux%Psorbmax     = sum_casaflux%Psorbmax     * rnsteps
       sum_casaflux%klitter      = sum_casaflux%klitter      * rnsteps
       sum_casaflux%ksoil        = sum_casaflux%ksoil        * rnsteps
       sum_casaflux%fromLtoS     = sum_casaflux%fromLtoS     * rnsteps
       sum_casaflux%fromStoS     = sum_casaflux%fromStoS     * rnsteps
       sum_casaflux%fromLtoCO2   = sum_casaflux%fromLtoCO2   * rnsteps
       sum_casaflux%fromStoCO2   = sum_casaflux%fromStoCO2   * rnsteps
       sum_casaflux%stemnpp      = sum_casaflux%stemnpp      * rnsteps
       sum_casaflux%frac_sapwood = sum_casaflux%frac_sapwood * rnsteps
       sum_casaflux%sapwood_area = sum_casaflux%sapwood_area * rnsteps
       sum_casaflux%Cplant_turnover = &
            sum_casaflux%Cplant_turnover                     * rnsteps
       sum_casaflux%Cplant_turnover_disturbance = &
            casaflux%Cplant_turnover_disturbance             * rnsteps
       sum_casaflux%Cplant_turnover_crowding = &
            sum_casaflux%Cplant_turnover_crowding            * rnsteps
       sum_casaflux%Cplant_turnover_resource_limitation = &
            sum_casaflux%Cplant_turnover_resource_limitation * rnsteps
       sum_casaflux%FluxCtolitter = sum_casaflux%FluxCtolitter * rnsteps
       sum_casaflux%FluxNtolitter = sum_casaflux%FluxNtolitter * rnsteps
       sum_casaflux%FluxPtolitter = sum_casaflux%FluxPtolitter * rnsteps

       sum_casaflux%FluxCtosoil = sum_casaflux%FluxCtosoil * rnsteps
       sum_casaflux%FluxNtosoil = sum_casaflux%FluxNtosoil * rnsteps
       sum_casaflux%FluxPtosoil = sum_casaflux%FluxPtosoil * rnsteps

       sum_casaflux%FluxCtoco2 =  sum_casaflux%FluxCtoco2 * rnsteps

       sum_casaflux%fromPtoL_fire          = sum_casaflux%fromPtoL_fire          * rnsteps
       sum_casaflux%kplant_fire            = sum_casaflux%kplant_fire            * rnsteps
       sum_casaflux%klitter_fire           = sum_casaflux%klitter_fire           * rnsteps
       sum_casaflux%kplant_tot             = sum_casaflux%kplant_tot             * rnsteps
       sum_casaflux%klitter_tot            = sum_casaflux%klitter_tot            * rnsteps
       sum_casaflux%FluxCtoCO2_plant_fire  = sum_casaflux%FluxCtoCO2_plant_fire  * rnsteps
       sum_casaflux%FluxCtoCO2_litter_fire = sum_casaflux%FluxCtoCO2_litter_fire * rnsteps
       sum_casaflux%fluxfromPtoCO2_fire    = sum_casaflux%fluxfromPtoCO2_fire    * rnsteps
       sum_casaflux%fluxfromLtoCO2_fire    = sum_casaflux%fluxfromLtoCO2_fire    * rnsteps

       sum_casaflux%FluxNtoAtm_fire = sum_casaflux%FluxNtoAtm_fire * rnsteps

       sum_casaflux%FluxFromPtoL       = sum_casaflux%FluxFromPtoL       * rnsteps
       sum_casaflux%FluxFromLtoS       = sum_casaflux%FluxFromLtoS       * rnsteps
       sum_casaflux%FluxFromStoS       = sum_casaflux%FluxFromStoS       * rnsteps
       sum_casaflux%FluxFromPtoCO2     = sum_casaflux%FluxFromPtoCO2     * rnsteps
       sum_casaflux%FluxFromLtoCO2     = sum_casaflux%FluxFromLtoCO2     * rnsteps
       sum_casaflux%FluxFromStoCO2     = sum_casaflux%FluxFromStoCO2     * rnsteps
       sum_casaflux%FluxFromPtoHarvest = sum_casaflux%FluxFromPtoHarvest * rnsteps
    endif ! average_now

  end subroutine update_sum_casa


  ! ------------------------------------------------------------------


  subroutine read_netcdf_casabiome(filename, casabiome)

    use netcdf, only: nf90_open, nf90_nowrite, &
         nf90_inq_varid, nf90_get_var, nf90_close
#ifdef __MPI__
    use mpi, only: MPI_Abort
#endif
    use cable_def_types_mod, only: nc_err

    implicit none

    character(len=*), intent(in)    :: filename
    type(casa_biome), intent(inout) :: casabiome

    logical :: existfile
    integer :: fid, vid
#ifdef __MPI__
    integer :: ierr
#endif

    ! open netCDF file
    inquire(file=trim(filename), exist=existfile)
    if (.not. existfile) then
       write(*,*) filename, ' does not exist!'
#ifdef __MPI__
       call MPI_Abort(0, 173, ierr)
#else
       stop 173
#endif
    endif

    ! open netCDF file
    call nc_err(nf90_open(trim(filename), nf90_nowrite, fid))

    ! read variables
    ! integer vectors
    call nc_err(nf90_inq_varid(fid, 'ivt2', vid))
    call nc_err(nf90_get_var(fid, vid, casabiome%ivt2))
    ! double vectors
    call nc_err(nf90_inq_varid(fid, 'xkleafcoldmax', vid))
    call nc_err(nf90_get_var(fid, vid, casabiome%xkleafcoldmax))
    call nc_err(nf90_inq_varid(fid, 'xkleafcoldexp', vid))
    call nc_err(nf90_get_var(fid, vid, casabiome%xkleafcoldexp))
    call nc_err(nf90_inq_varid(fid, 'xkleafdrymax', vid))
    call nc_err(nf90_get_var(fid, vid, casabiome%xkleafdrymax))
    call nc_err(nf90_inq_varid(fid, 'xkleafdryexp', vid))
    call nc_err(nf90_get_var(fid, vid, casabiome%xkleafdryexp))
    call nc_err(nf90_inq_varid(fid, 'glaimax', vid))
    call nc_err(nf90_get_var(fid, vid, casabiome%glaimax))
    call nc_err(nf90_inq_varid(fid, 'glaimin', vid))
    call nc_err(nf90_get_var(fid, vid, casabiome%glaimin))
    call nc_err(nf90_inq_varid(fid, 'sla', vid))
    call nc_err(nf90_get_var(fid, vid, casabiome%sla))
    call nc_err(nf90_inq_varid(fid, 'ratiofrootleaf', vid))
    call nc_err(nf90_get_var(fid, vid, casabiome%ratiofrootleaf))
    call nc_err(nf90_inq_varid(fid, 'kroot', vid))
    call nc_err(nf90_get_var(fid, vid, casabiome%kroot))
    call nc_err(nf90_inq_varid(fid, 'krootlen', vid))
    call nc_err(nf90_get_var(fid, vid, casabiome%krootlen))
    call nc_err(nf90_inq_varid(fid, 'rootdepth', vid))
    call nc_err(nf90_get_var(fid, vid, casabiome%rootdepth))
    call nc_err(nf90_inq_varid(fid, 'kuptake', vid))
    call nc_err(nf90_get_var(fid, vid, casabiome%kuptake))
    call nc_err(nf90_inq_varid(fid, 'kminn', vid))
    call nc_err(nf90_get_var(fid, vid, casabiome%kminn))
    call nc_err(nf90_inq_varid(fid, 'kuplabp', vid))
    call nc_err(nf90_get_var(fid, vid, casabiome%kuplabp))
    call nc_err(nf90_inq_varid(fid, 'kclabrate', vid))
    call nc_err(nf90_get_var(fid, vid, casabiome%kclabrate))
    call nc_err(nf90_inq_varid(fid, 'xnpmax', vid))
    call nc_err(nf90_get_var(fid, vid, casabiome%xnpmax))
    call nc_err(nf90_inq_varid(fid, 'q10soil', vid))
    call nc_err(nf90_get_var(fid, vid, casabiome%q10soil))
    call nc_err(nf90_inq_varid(fid, 'xkoptlitter', vid))
    call nc_err(nf90_get_var(fid, vid, casabiome%xkoptlitter))
    call nc_err(nf90_inq_varid(fid, 'xkoptsoil', vid))
    call nc_err(nf90_get_var(fid, vid, casabiome%xkoptsoil))
    call nc_err(nf90_inq_varid(fid, 'xkplab', vid))
    call nc_err(nf90_get_var(fid, vid, casabiome%xkplab))
    call nc_err(nf90_inq_varid(fid, 'xkpsorb', vid))
    call nc_err(nf90_get_var(fid, vid, casabiome%xkpsorb))
    call nc_err(nf90_inq_varid(fid, 'xkpocc', vid))
    call nc_err(nf90_get_var(fid, vid, casabiome%xkpocc))
    call nc_err(nf90_inq_varid(fid, 'prodptase', vid))
    call nc_err(nf90_get_var(fid, vid, casabiome%prodptase))
    call nc_err(nf90_inq_varid(fid, 'costnpup', vid))
    call nc_err(nf90_get_var(fid, vid, casabiome%costnpup))
    call nc_err(nf90_inq_varid(fid, 'maxfinelitter', vid))
    call nc_err(nf90_get_var(fid, vid, casabiome%maxfinelitter))
    call nc_err(nf90_inq_varid(fid, 'maxcwd', vid))
    call nc_err(nf90_get_var(fid, vid, casabiome%maxcwd))
    call nc_err(nf90_inq_varid(fid, 'nintercept', vid))
    call nc_err(nf90_get_var(fid, vid, casabiome%nintercept))
    call nc_err(nf90_inq_varid(fid, 'nslope', vid))
    call nc_err(nf90_get_var(fid, vid, casabiome%nslope))
    call nc_err(nf90_inq_varid(fid, 'la_to_sa', vid))
    call nc_err(nf90_get_var(fid, vid, casabiome%la_to_sa))
    call nc_err(nf90_inq_varid(fid, 'vcmax_scalar', vid))
    call nc_err(nf90_get_var(fid, vid, casabiome%vcmax_scalar))
    call nc_err(nf90_inq_varid(fid, 'disturbance_interval', vid))
    call nc_err(nf90_get_var(fid, vid, casabiome%disturbance_interval))
    call nc_err(nf90_inq_varid(fid, 'damm_enzpool', vid))
    call nc_err(nf90_get_var(fid, vid, casabiome%damm_enzpool))
    call nc_err(nf90_inq_varid(fid, 'damm_kmo2', vid))
    call nc_err(nf90_get_var(fid, vid, casabiome%damm_kmo2))
    call nc_err(nf90_inq_varid(fid, 'damm_kmcp', vid))
    call nc_err(nf90_get_var(fid, vid, casabiome%damm_kmcp))
    call nc_err(nf90_inq_varid(fid, 'damm_ea', vid))
    call nc_err(nf90_get_var(fid, vid, casabiome%damm_ea))
    call nc_err(nf90_inq_varid(fid, 'damm_alpha', vid))
    call nc_err(nf90_get_var(fid, vid, casabiome%damm_alpha))
    ! double arrays
    call nc_err(nf90_inq_varid(fid, 'plantrate', vid))
    call nc_err(nf90_get_var(fid, vid, casabiome%plantrate))
    call nc_err(nf90_inq_varid(fid, 'rmplant', vid))
    call nc_err(nf90_get_var(fid, vid, casabiome%rmplant))
    call nc_err(nf90_inq_varid(fid, 'fracnpptop', vid))
    call nc_err(nf90_get_var(fid, vid, casabiome%fracnpptop))
    call nc_err(nf90_inq_varid(fid, 'fraclignin', vid))
    call nc_err(nf90_get_var(fid, vid, casabiome%fraclignin))
    call nc_err(nf90_inq_varid(fid, 'fraclabile', vid))
    call nc_err(nf90_get_var(fid, vid, casabiome%fraclabile))
    call nc_err(nf90_inq_varid(fid, 'rationcplantmin', vid))
    call nc_err(nf90_get_var(fid, vid, casabiome%rationcplantmin))
    call nc_err(nf90_inq_varid(fid, 'rationcplantmax', vid))
    call nc_err(nf90_get_var(fid, vid, casabiome%rationcplantmax))
    call nc_err(nf90_inq_varid(fid, 'rationpplantmin', vid))
    call nc_err(nf90_get_var(fid, vid, casabiome%rationpplantmin))
    call nc_err(nf90_inq_varid(fid, 'rationpplantmax', vid))
    call nc_err(nf90_get_var(fid, vid, casabiome%rationpplantmax))
    call nc_err(nf90_inq_varid(fid, 'fracligninplant', vid))
    call nc_err(nf90_get_var(fid, vid, casabiome%fracligninplant))
    call nc_err(nf90_inq_varid(fid, 'ftransnptol', vid))
    call nc_err(nf90_get_var(fid, vid, casabiome%ftransnptol))
    call nc_err(nf90_inq_varid(fid, 'ftranspptol', vid))
    call nc_err(nf90_get_var(fid, vid, casabiome%ftranspptol))
    call nc_err(nf90_inq_varid(fid, 'litterrate', vid))
    call nc_err(nf90_get_var(fid, vid, casabiome%litterrate))
    call nc_err(nf90_inq_varid(fid, 'ratiopcplantmin', vid))
    call nc_err(nf90_get_var(fid, vid, casabiome%ratiopcplantmin))
    call nc_err(nf90_inq_varid(fid, 'ratiopcplantmax', vid))
    call nc_err(nf90_get_var(fid, vid, casabiome%ratiopcplantmax))
    call nc_err(nf90_inq_varid(fid, 'soilrate', vid))
    call nc_err(nf90_get_var(fid, vid, casabiome%soilrate))

    ! close NetCDF file
    call nc_err(nf90_close(fid))

  end subroutine read_netcdf_casabiome


  subroutine read_netcdf_casapool(filename, casapool)

    use netcdf, only: nf90_open, nf90_nowrite, &
         nf90_inq_varid, nf90_get_var, nf90_close
#ifdef __MPI__
    use mpi, only: MPI_Abort
#endif
    use cable_def_types_mod, only: nc_err

    implicit none

    character(len=*), intent(in)    :: filename
    type(casa_pool),  intent(inout) :: casapool

    logical :: existfile
    integer :: fid, vid
#ifdef __MPI__
    integer :: ierr
#endif

    ! open netCDF file
    inquire(file=trim(filename), exist=existfile)
    if (.not. existfile) then
       write(*,*) filename, ' does not exist!'
#ifdef __MPI__
       call MPI_Abort(0, 174, ierr)
#else
       stop 174
#endif
    endif

    ! open netCDF file
    call nc_err(nf90_open(trim(filename), nf90_nowrite, fid))

    ! read variables
    call nc_err(nf90_inq_varid(fid, 'clabile', vid))
    call nc_err(nf90_get_var(fid, vid, casapool%clabile))
    call nc_err(nf90_inq_varid(fid, 'dclabiledt', vid))
    call nc_err(nf90_get_var(fid, vid, casapool%dclabiledt))
    call nc_err(nf90_inq_varid(fid, 'ctot', vid))
    call nc_err(nf90_get_var(fid, vid, casapool%ctot))
    call nc_err(nf90_inq_varid(fid, 'ctot_0', vid))
    call nc_err(nf90_get_var(fid, vid, casapool%ctot_0))
    call nc_err(nf90_inq_varid(fid, 'cplant', vid))
    call nc_err(nf90_get_var(fid, vid, casapool%cplant))
    call nc_err(nf90_inq_varid(fid, 'nplant', vid))
    call nc_err(nf90_get_var(fid, vid, casapool%nplant))
    call nc_err(nf90_inq_varid(fid, 'pplant', vid))
    call nc_err(nf90_get_var(fid, vid, casapool%pplant))
    call nc_err(nf90_inq_varid(fid, 'dcplantdt', vid))
    call nc_err(nf90_get_var(fid, vid, casapool%dcplantdt))
    call nc_err(nf90_inq_varid(fid, 'dnplantdt', vid))
    call nc_err(nf90_get_var(fid, vid, casapool%dnplantdt))
    call nc_err(nf90_inq_varid(fid, 'dpplantdt', vid))
    call nc_err(nf90_get_var(fid, vid, casapool%dpplantdt))
    call nc_err(nf90_inq_varid(fid, 'rationcplant', vid))
    call nc_err(nf90_get_var(fid, vid, casapool%rationcplant))
    call nc_err(nf90_inq_varid(fid, 'rationpplant', vid))
    call nc_err(nf90_get_var(fid, vid, casapool%rationpplant))
    call nc_err(nf90_inq_varid(fid, 'nsoilmin', vid))
    call nc_err(nf90_get_var(fid, vid, casapool%nsoilmin))
    call nc_err(nf90_inq_varid(fid, 'psoillab', vid))
    call nc_err(nf90_get_var(fid, vid, casapool%psoillab))
    call nc_err(nf90_inq_varid(fid, 'psoilsorb', vid))
    call nc_err(nf90_get_var(fid, vid, casapool%psoilsorb))
    call nc_err(nf90_inq_varid(fid, 'psoilocc', vid))
    call nc_err(nf90_get_var(fid, vid, casapool%psoilocc))
    call nc_err(nf90_inq_varid(fid, 'dnsoilmindt', vid))
    call nc_err(nf90_get_var(fid, vid, casapool%dnsoilmindt))
    call nc_err(nf90_inq_varid(fid, 'dpsoillabdt', vid))
    call nc_err(nf90_get_var(fid, vid, casapool%dpsoillabdt))
    call nc_err(nf90_inq_varid(fid, 'dpsoilsorbdt', vid))
    call nc_err(nf90_get_var(fid, vid, casapool%dpsoilsorbdt))
    call nc_err(nf90_inq_varid(fid, 'dpsoiloccdt', vid))
    call nc_err(nf90_get_var(fid, vid, casapool%dpsoiloccdt))
    call nc_err(nf90_inq_varid(fid, 'clitter', vid))
    call nc_err(nf90_get_var(fid, vid, casapool%clitter))
    call nc_err(nf90_inq_varid(fid, 'nlitter', vid))
    call nc_err(nf90_get_var(fid, vid, casapool%nlitter))
    call nc_err(nf90_inq_varid(fid, 'plitter', vid))
    call nc_err(nf90_get_var(fid, vid, casapool%plitter))
    call nc_err(nf90_inq_varid(fid, 'dclitterdt', vid))
    call nc_err(nf90_get_var(fid, vid, casapool%dclitterdt))
    call nc_err(nf90_inq_varid(fid, 'dnlitterdt', vid))
    call nc_err(nf90_get_var(fid, vid, casapool%dnlitterdt))
    call nc_err(nf90_inq_varid(fid, 'dplitterdt', vid))
    call nc_err(nf90_get_var(fid, vid, casapool%dplitterdt))
    call nc_err(nf90_inq_varid(fid, 'rationclitter', vid))
    call nc_err(nf90_get_var(fid, vid, casapool%rationclitter))
    call nc_err(nf90_inq_varid(fid, 'rationplitter', vid))
    call nc_err(nf90_get_var(fid, vid, casapool%rationplitter))
    call nc_err(nf90_inq_varid(fid, 'csoil', vid))
    call nc_err(nf90_get_var(fid, vid, casapool%csoil))
    call nc_err(nf90_inq_varid(fid, 'nsoil', vid))
    call nc_err(nf90_get_var(fid, vid, casapool%nsoil))
    call nc_err(nf90_inq_varid(fid, 'psoil', vid))
    call nc_err(nf90_get_var(fid, vid, casapool%psoil))
    call nc_err(nf90_inq_varid(fid, 'dcsoildt', vid))
    call nc_err(nf90_get_var(fid, vid, casapool%dcsoildt))
    call nc_err(nf90_inq_varid(fid, 'dnsoildt', vid))
    call nc_err(nf90_get_var(fid, vid, casapool%dnsoildt))
    call nc_err(nf90_inq_varid(fid, 'dpsoildt', vid))
    call nc_err(nf90_get_var(fid, vid, casapool%dpsoildt))
    call nc_err(nf90_inq_varid(fid, 'rationcsoil', vid))
    call nc_err(nf90_get_var(fid, vid, casapool%rationcsoil))
    call nc_err(nf90_inq_varid(fid, 'rationcsoilnew', vid))
    call nc_err(nf90_get_var(fid, vid, casapool%rationcsoilnew))
    call nc_err(nf90_inq_varid(fid, 'rationpsoil', vid))
    call nc_err(nf90_get_var(fid, vid, casapool%rationpsoil))
    call nc_err(nf90_inq_varid(fid, 'rationcsoilmin', vid))
    call nc_err(nf90_get_var(fid, vid, casapool%rationcsoilmin))
    call nc_err(nf90_inq_varid(fid, 'rationcsoilmax', vid))
    call nc_err(nf90_get_var(fid, vid, casapool%rationcsoilmax))
    call nc_err(nf90_inq_varid(fid, 'ratiopcsoil', vid))
    call nc_err(nf90_get_var(fid, vid, casapool%ratiopcsoil))
    call nc_err(nf90_inq_varid(fid, 'ratiopcplant', vid))
    call nc_err(nf90_get_var(fid, vid, casapool%ratiopcplant))
    call nc_err(nf90_inq_varid(fid, 'ratiopclitter', vid))
    call nc_err(nf90_get_var(fid, vid, casapool%ratiopclitter))

    ! close NetCDF file
    call nc_err(nf90_close(fid))

  end subroutine read_netcdf_casapool


  subroutine read_netcdf_casaflux(filename, casaflux)

    use netcdf, only: nf90_open, nf90_nowrite, &
         nf90_inq_varid, nf90_get_var, nf90_close
#ifdef __MPI__
    use mpi, only: MPI_Abort
#endif
    use cable_def_types_mod, only: nc_err

    implicit none

    character(len=*), intent(in)    :: filename
    type(casa_flux),  intent(inout) :: casaflux

    logical :: existfile
    integer :: fid, vid
#ifdef __MPI__
    integer :: ierr
#endif

    ! open netCDF file
    inquire(file=trim(filename), exist=existfile)
    if (.not. existfile) then
       write(*,*) filename, ' does not exist!'
#ifdef __MPI__
       call MPI_Abort(0, 175, ierr)
#else
       stop 175
#endif
    endif

    ! open netCDF file
    call nc_err(nf90_open(trim(filename), nf90_nowrite, fid))

    ! read variables
    call nc_err(nf90_inq_varid(fid, 'cgpp', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%cgpp))
    call nc_err(nf90_inq_varid(fid, 'cnpp', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%cnpp))
    call nc_err(nf90_inq_varid(fid, 'crp', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%crp))
    call nc_err(nf90_inq_varid(fid, 'crgplant', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%crgplant))
    call nc_err(nf90_inq_varid(fid, 'nminfix', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%nminfix))
    call nc_err(nf90_inq_varid(fid, 'nminuptake', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%nminuptake))
    call nc_err(nf90_inq_varid(fid, 'plabuptake', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%plabuptake))
    call nc_err(nf90_inq_varid(fid, 'clabloss', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%clabloss))
    call nc_err(nf90_inq_varid(fid, 'fracclabile', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%fracclabile))
    call nc_err(nf90_inq_varid(fid, 'stemnpp', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%stemnpp))
    call nc_err(nf90_inq_varid(fid, 'frac_sapwood', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%frac_sapwood))
    call nc_err(nf90_inq_varid(fid, 'sapwood_area', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%sapwood_area))
    call nc_err(nf90_inq_varid(fid, 'charvest', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%charvest))
    call nc_err(nf90_inq_varid(fid, 'nharvest', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%nharvest))
    call nc_err(nf90_inq_varid(fid, 'pharvest', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%pharvest))
    call nc_err(nf90_inq_varid(fid, 'fharvest', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%fharvest))
    call nc_err(nf90_inq_varid(fid, 'fcrop', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%fcrop))
    call nc_err(nf90_inq_varid(fid, 'fraccalloc', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%fraccalloc))
    call nc_err(nf90_inq_varid(fid, 'fracnalloc', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%fracnalloc))
    call nc_err(nf90_inq_varid(fid, 'fracpalloc', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%fracpalloc))
    call nc_err(nf90_inq_varid(fid, 'crmplant', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%crmplant))
    call nc_err(nf90_inq_varid(fid, 'kplant', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%kplant))
    call nc_err(nf90_inq_varid(fid, 'cplant_turnover', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%cplant_turnover))
    call nc_err(nf90_inq_varid(fid, 'fromptol', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%fromptol))
    call nc_err(nf90_inq_varid(fid, 'cnep', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%cnep))
    call nc_err(nf90_inq_varid(fid, 'crsoil', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%crsoil))
    call nc_err(nf90_inq_varid(fid, 'nmindep', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%nmindep))
    call nc_err(nf90_inq_varid(fid, 'nminloss', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%nminloss))
    call nc_err(nf90_inq_varid(fid, 'nminleach', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%nminleach))
    call nc_err(nf90_inq_varid(fid, 'nupland', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%nupland))
    call nc_err(nf90_inq_varid(fid, 'nlittermin', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%nlittermin))
    call nc_err(nf90_inq_varid(fid, 'nsmin', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%nsmin))
    call nc_err(nf90_inq_varid(fid, 'nsimm', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%nsimm))
    call nc_err(nf90_inq_varid(fid, 'nsnet', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%nsnet))
    call nc_err(nf90_inq_varid(fid, 'fnminloss', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%fnminloss))
    call nc_err(nf90_inq_varid(fid, 'fnminleach', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%fnminleach))
    call nc_err(nf90_inq_varid(fid, 'pdep', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%pdep))
    call nc_err(nf90_inq_varid(fid, 'pwea', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%pwea))
    call nc_err(nf90_inq_varid(fid, 'pleach', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%pleach))
    call nc_err(nf90_inq_varid(fid, 'ploss', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%ploss))
    call nc_err(nf90_inq_varid(fid, 'pupland', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%pupland))
    call nc_err(nf90_inq_varid(fid, 'plittermin', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%plittermin))
    call nc_err(nf90_inq_varid(fid, 'psmin', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%psmin))
    call nc_err(nf90_inq_varid(fid, 'psimm', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%psimm))
    call nc_err(nf90_inq_varid(fid, 'psnet', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%psnet))
    call nc_err(nf90_inq_varid(fid, 'fpleach', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%fpleach))
    call nc_err(nf90_inq_varid(fid, 'kplab', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%kplab))
    call nc_err(nf90_inq_varid(fid, 'kpsorb', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%kpsorb))
    call nc_err(nf90_inq_varid(fid, 'kpocc', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%kpocc))
    call nc_err(nf90_inq_varid(fid, 'kmlabp', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%kmlabp))
    call nc_err(nf90_inq_varid(fid, 'psorbmax', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%psorbmax))
    call nc_err(nf90_inq_varid(fid, 'cplant_turnover_disturbance', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%cplant_turnover_disturbance))
    call nc_err(nf90_inq_varid(fid, 'cplant_turnover_crowding', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%cplant_turnover_crowding))
    call nc_err(nf90_inq_varid(fid, 'cplant_turnover_resource_limitation', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%cplant_turnover_resource_limitation))
    call nc_err(nf90_inq_varid(fid, 'klitter', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%klitter))
    call nc_err(nf90_inq_varid(fid, 'ksoil', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%ksoil))
    call nc_err(nf90_inq_varid(fid, 'fromltos', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%fromltos))
    call nc_err(nf90_inq_varid(fid, 'fromstos', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%fromstos))
    call nc_err(nf90_inq_varid(fid, 'fromltoco2', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%fromltoco2))
    call nc_err(nf90_inq_varid(fid, 'fromstoco2', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%fromstoco2))
    call nc_err(nf90_inq_varid(fid, 'fluxctolitter', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%fluxctolitter))
    call nc_err(nf90_inq_varid(fid, 'fluxntolitter', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%fluxntolitter))
    call nc_err(nf90_inq_varid(fid, 'fluxptolitter', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%fluxptolitter))
    call nc_err(nf90_inq_varid(fid, 'fluxctosoil', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%fluxctosoil))
    call nc_err(nf90_inq_varid(fid, 'fluxntosoil', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%fluxntosoil))
    call nc_err(nf90_inq_varid(fid, 'fluxptosoil', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%fluxptosoil))
    call nc_err(nf90_inq_varid(fid, 'fluxctoco2', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%fluxctoco2))
    call nc_err(nf90_inq_varid(fid, 'fluxctohwp', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%fluxctohwp))
    call nc_err(nf90_inq_varid(fid, 'fluxntohwp', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%fluxntohwp))
    call nc_err(nf90_inq_varid(fid, 'fluxptohwp', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%fluxptohwp))
    call nc_err(nf90_inq_varid(fid, 'fluxctoclear', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%fluxctoclear))
    call nc_err(nf90_inq_varid(fid, 'fluxntoclear', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%fluxntoclear))
    call nc_err(nf90_inq_varid(fid, 'fluxptoclear', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%fluxptoclear))
    call nc_err(nf90_inq_varid(fid, 'ctransferluc', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%ctransferluc))
    call nc_err(nf90_inq_varid(fid, 'fromptol_fire', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%fromptol_fire))
    call nc_err(nf90_inq_varid(fid, 'klitter_fire', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%klitter_fire))
    call nc_err(nf90_inq_varid(fid, 'klitter_tot', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%klitter_tot))
    call nc_err(nf90_inq_varid(fid, 'kplant_fire', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%kplant_fire))
    call nc_err(nf90_inq_varid(fid, 'kplant_tot', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%kplant_tot))
    call nc_err(nf90_inq_varid(fid, 'fluxctoco2_plant_fire', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%fluxctoco2_plant_fire))
    call nc_err(nf90_inq_varid(fid, 'fluxctoco2_litter_fire', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%fluxctoco2_litter_fire))
    call nc_err(nf90_inq_varid(fid, 'fluxfromptoco2_fire', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%fluxfromptoco2_fire))
    call nc_err(nf90_inq_varid(fid, 'fluxfromltoco2_fire', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%fluxfromltoco2_fire))
    call nc_err(nf90_inq_varid(fid, 'fluxntoatm_fire', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%fluxntoatm_fire))
    call nc_err(nf90_inq_varid(fid, 'fluxfromptol', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%fluxfromptol))
    call nc_err(nf90_inq_varid(fid, 'fluxfromltos', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%fluxfromltos))
    call nc_err(nf90_inq_varid(fid, 'fluxfromstos', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%fluxfromstos))
    call nc_err(nf90_inq_varid(fid, 'fluxfromptoco2', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%fluxfromptoco2))
    call nc_err(nf90_inq_varid(fid, 'fluxfromltoco2', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%fluxfromltoco2))
    call nc_err(nf90_inq_varid(fid, 'fluxfromstoco2', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%fluxfromstoco2))
    call nc_err(nf90_inq_varid(fid, 'fluxfromptoharvest', vid))
    call nc_err(nf90_get_var(fid, vid, casaflux%fluxfromptoharvest))

    ! close NetCDF file
    call nc_err(nf90_close(fid))

  end subroutine read_netcdf_casaflux


  subroutine read_netcdf_casamet(filename, casamet)

    use netcdf, only: nf90_open, nf90_nowrite, &
         nf90_inq_varid, nf90_get_var, nf90_close
#ifdef __MPI__
    use mpi, only: MPI_Abort
#endif
    use cable_def_types_mod, only: nc_err

    implicit none

    character(len=*), intent(in)    :: filename
    type(casa_met),   intent(inout) :: casamet

    logical :: existfile
    integer :: fid, vid
#ifdef __MPI__
    integer :: ierr
#endif

    ! open netCDF file
    inquire(file=trim(filename), exist=existfile)
    if (.not. existfile) then
       write(*,*) filename, ' does not exist!'
#ifdef __MPI__
       call MPI_Abort(0, 176, ierr)
#else
       stop 176
#endif
    endif

    ! open netCDF file
    call nc_err(nf90_open(trim(filename), nf90_nowrite, fid))

    ! read variables
    call nc_err(nf90_inq_varid(fid, 'glai', vid))
    call nc_err(nf90_get_var(fid, vid, casamet%glai))
    call nc_err(nf90_inq_varid(fid, 'tairk', vid))
    call nc_err(nf90_get_var(fid, vid, casamet%tairk))
    call nc_err(nf90_inq_varid(fid, 'precip', vid))
    call nc_err(nf90_get_var(fid, vid, casamet%precip))
    call nc_err(nf90_inq_varid(fid, 'tsoilavg', vid))
    call nc_err(nf90_get_var(fid, vid, casamet%tsoilavg))
    call nc_err(nf90_inq_varid(fid, 'moistavg', vid))
    call nc_err(nf90_get_var(fid, vid, casamet%moistavg))
    call nc_err(nf90_inq_varid(fid, 'btran', vid))
    call nc_err(nf90_get_var(fid, vid, casamet%btran))
    call nc_err(nf90_inq_varid(fid, 'lnonwood', vid))
    call nc_err(nf90_get_var(fid, vid, casamet%lnonwood))
    call nc_err(nf90_inq_varid(fid, 'tsoil', vid))
    call nc_err(nf90_get_var(fid, vid, casamet%tsoil))
    call nc_err(nf90_inq_varid(fid, 'moist', vid))
    call nc_err(nf90_get_var(fid, vid, casamet%moist))
    call nc_err(nf90_inq_varid(fid, 'iveg2', vid))
    call nc_err(nf90_get_var(fid, vid, casamet%iveg2))
    call nc_err(nf90_inq_varid(fid, 'ijgcm', vid))
    call nc_err(nf90_get_var(fid, vid, casamet%ijgcm))
    call nc_err(nf90_inq_varid(fid, 'isorder', vid))
    call nc_err(nf90_get_var(fid, vid, casamet%isorder))
    call nc_err(nf90_inq_varid(fid, 'lat', vid))
    call nc_err(nf90_get_var(fid, vid, casamet%lat))
    call nc_err(nf90_inq_varid(fid, 'lon', vid))
    call nc_err(nf90_get_var(fid, vid, casamet%lon))
    call nc_err(nf90_inq_varid(fid, 'areacell', vid))
    call nc_err(nf90_get_var(fid, vid, casamet%areacell))
    call nc_err(nf90_inq_varid(fid, 'tairkspin', vid))
    call nc_err(nf90_get_var(fid, vid, casamet%tairkspin))
    call nc_err(nf90_inq_varid(fid, 'cgppspin', vid))
    call nc_err(nf90_get_var(fid, vid, casamet%cgppspin))
    call nc_err(nf90_inq_varid(fid, 'crmplantspin_1', vid))
    call nc_err(nf90_get_var(fid, vid, casamet%crmplantspin_1))
    call nc_err(nf90_inq_varid(fid, 'crmplantspin_2', vid))
    call nc_err(nf90_get_var(fid, vid, casamet%crmplantspin_2))
    call nc_err(nf90_inq_varid(fid, 'crmplantspin_3', vid))
    call nc_err(nf90_get_var(fid, vid, casamet%crmplantspin_3))
    call nc_err(nf90_inq_varid(fid, 'tsoilspin_1', vid))
    call nc_err(nf90_get_var(fid, vid, casamet%tsoilspin_1))
    call nc_err(nf90_inq_varid(fid, 'tsoilspin_2', vid))
    call nc_err(nf90_get_var(fid, vid, casamet%tsoilspin_2))
    call nc_err(nf90_inq_varid(fid, 'tsoilspin_3', vid))
    call nc_err(nf90_get_var(fid, vid, casamet%tsoilspin_3))
    call nc_err(nf90_inq_varid(fid, 'tsoilspin_4', vid))
    call nc_err(nf90_get_var(fid, vid, casamet%tsoilspin_4))
    call nc_err(nf90_inq_varid(fid, 'tsoilspin_5', vid))
    call nc_err(nf90_get_var(fid, vid, casamet%tsoilspin_5))
    call nc_err(nf90_inq_varid(fid, 'tsoilspin_6', vid))
    call nc_err(nf90_get_var(fid, vid, casamet%tsoilspin_6))
    call nc_err(nf90_inq_varid(fid, 'moistspin_1', vid))
    call nc_err(nf90_get_var(fid, vid, casamet%moistspin_1))
    call nc_err(nf90_inq_varid(fid, 'moistspin_2', vid))
    call nc_err(nf90_get_var(fid, vid, casamet%moistspin_2))
    call nc_err(nf90_inq_varid(fid, 'moistspin_3', vid))
    call nc_err(nf90_get_var(fid, vid, casamet%moistspin_3))
    call nc_err(nf90_inq_varid(fid, 'moistspin_4', vid))
    call nc_err(nf90_get_var(fid, vid, casamet%moistspin_4))
    call nc_err(nf90_inq_varid(fid, 'moistspin_5', vid))
    call nc_err(nf90_get_var(fid, vid, casamet%moistspin_5))
    call nc_err(nf90_inq_varid(fid, 'moistspin_6', vid))
    call nc_err(nf90_get_var(fid, vid, casamet%moistspin_6))
    call nc_err(nf90_inq_varid(fid, 'mtempspin', vid))
    call nc_err(nf90_get_var(fid, vid, casamet%mtempspin))
    call nc_err(nf90_inq_varid(fid, 'frecspin', vid))
    call nc_err(nf90_get_var(fid, vid, casamet%frecspin))
    call nc_err(nf90_inq_varid(fid, 'can12spin', vid))
    call nc_err(nf90_get_var(fid, vid, casamet%can12spin))
    call nc_err(nf90_inq_varid(fid, 'can13spin', vid))
    call nc_err(nf90_get_var(fid, vid, casamet%can13spin))
    call nc_err(nf90_inq_varid(fid, 'dprecip_spin', vid))
    call nc_err(nf90_get_var(fid, vid, casamet%dprecip_spin))
    call nc_err(nf90_inq_varid(fid, 'aprecip_av20_spin', vid))
    call nc_err(nf90_get_var(fid, vid, casamet%aprecip_av20_spin))
    call nc_err(nf90_inq_varid(fid, 'du10_max_spin', vid))
    call nc_err(nf90_get_var(fid, vid, casamet%du10_max_spin))
    call nc_err(nf90_inq_varid(fid, 'drhum_spin', vid))
    call nc_err(nf90_get_var(fid, vid, casamet%drhum_spin))
    call nc_err(nf90_inq_varid(fid, 'dtemp_max_spin', vid))
    call nc_err(nf90_get_var(fid, vid, casamet%dtemp_max_spin))
    call nc_err(nf90_inq_varid(fid, 'dtemp_min_spin', vid))
    call nc_err(nf90_get_var(fid, vid, casamet%dtemp_min_spin))
    call nc_err(nf90_inq_varid(fid, 'kbdi_spin', vid))
    call nc_err(nf90_get_var(fid, vid, casamet%kbdi_spin))
    call nc_err(nf90_inq_varid(fid, 'd_macarthur_spin', vid))
    call nc_err(nf90_get_var(fid, vid, casamet%d_macarthur_spin))
    call nc_err(nf90_inq_varid(fid, 'ffdi_spin', vid))
    call nc_err(nf90_get_var(fid, vid, casamet%ffdi_spin))
    call nc_err(nf90_inq_varid(fid, 'last_precip_spin', vid))
    call nc_err(nf90_get_var(fid, vid, casamet%last_precip_spin))
    call nc_err(nf90_inq_varid(fid, 'dslr_spin', vid))
    call nc_err(nf90_get_var(fid, vid, casamet%dslr_spin))

    ! close NetCDF file
    call nc_err(nf90_close(fid))

  end subroutine read_netcdf_casamet


  subroutine read_netcdf_casabal(filename, casabal)

    use netcdf, only: nf90_open, nf90_nowrite, &
         nf90_inq_varid, nf90_get_var, nf90_close
#ifdef __MPI__
    use mpi, only: MPI_Abort
#endif
    use cable_def_types_mod, only: nc_err

    implicit none

    character(len=*),   intent(in)    :: filename
    type(casa_balance), intent(inout) :: casabal

    logical :: existfile
    integer :: fid, vid
#ifdef __MPI__
    integer :: ierr
#endif

    ! open netCDF file
    inquire(file=trim(filename), exist=existfile)
    if (.not. existfile) then
       write(*,*) filename, ' does not exist!'
#ifdef __MPI__
       call MPI_Abort(0, 177, ierr)
#else
       stop 177
#endif
    endif

    ! open netCDF file
    call nc_err(nf90_open(trim(filename), nf90_nowrite, fid))

    ! read variables
    call nc_err(nf90_inq_varid(fid, 'fcgppyear', vid))
    call nc_err(nf90_get_var(fid, vid, casabal%fcgppyear))
    call nc_err(nf90_inq_varid(fid, 'fcnppyear', vid))
    call nc_err(nf90_get_var(fid, vid, casabal%fcnppyear))
    call nc_err(nf90_inq_varid(fid, 'fcrmleafyear', vid))
    call nc_err(nf90_get_var(fid, vid, casabal%fcrmleafyear))
    call nc_err(nf90_inq_varid(fid, 'fcrmwoodyear', vid))
    call nc_err(nf90_get_var(fid, vid, casabal%fcrmwoodyear))
    call nc_err(nf90_inq_varid(fid, 'fcrmrootyear', vid))
    call nc_err(nf90_get_var(fid, vid, casabal%fcrmrootyear))
    call nc_err(nf90_inq_varid(fid, 'fcrgrowyear', vid))
    call nc_err(nf90_get_var(fid, vid, casabal%fcrgrowyear))
    call nc_err(nf90_inq_varid(fid, 'fcrpyear', vid))
    call nc_err(nf90_get_var(fid, vid, casabal%fcrpyear))
    call nc_err(nf90_inq_varid(fid, 'fcrsyear', vid))
    call nc_err(nf90_get_var(fid, vid, casabal%fcrsyear))
    call nc_err(nf90_inq_varid(fid, 'fcneeyear', vid))
    call nc_err(nf90_get_var(fid, vid, casabal%fcneeyear))
    call nc_err(nf90_inq_varid(fid, 'dcdtyear', vid))
    call nc_err(nf90_get_var(fid, vid, casabal%dcdtyear))
    call nc_err(nf90_inq_varid(fid, 'laimax', vid))
    call nc_err(nf90_get_var(fid, vid, casabal%laimax))
    call nc_err(nf90_inq_varid(fid, 'cleafmean', vid))
    call nc_err(nf90_get_var(fid, vid, casabal%cleafmean))
    call nc_err(nf90_inq_varid(fid, 'crootmean', vid))
    call nc_err(nf90_get_var(fid, vid, casabal%crootmean))
    call nc_err(nf90_inq_varid(fid, 'fndepyear', vid))
    call nc_err(nf90_get_var(fid, vid, casabal%fndepyear))
    call nc_err(nf90_inq_varid(fid, 'fnfixyear', vid))
    call nc_err(nf90_get_var(fid, vid, casabal%fnfixyear))
    call nc_err(nf90_inq_varid(fid, 'fnsnetyear', vid))
    call nc_err(nf90_get_var(fid, vid, casabal%fnsnetyear))
    call nc_err(nf90_inq_varid(fid, 'fnupyear', vid))
    call nc_err(nf90_get_var(fid, vid, casabal%fnupyear))
    call nc_err(nf90_inq_varid(fid, 'fnleachyear', vid))
    call nc_err(nf90_get_var(fid, vid, casabal%fnleachyear))
    call nc_err(nf90_inq_varid(fid, 'fnlossyear', vid))
    call nc_err(nf90_get_var(fid, vid, casabal%fnlossyear))
    call nc_err(nf90_inq_varid(fid, 'fpweayear', vid))
    call nc_err(nf90_get_var(fid, vid, casabal%fpweayear))
    call nc_err(nf90_inq_varid(fid, 'fpdustyear', vid))
    call nc_err(nf90_get_var(fid, vid, casabal%fpdustyear))
    call nc_err(nf90_inq_varid(fid, 'fpsnetyear', vid))
    call nc_err(nf90_get_var(fid, vid, casabal%fpsnetyear))
    call nc_err(nf90_inq_varid(fid, 'fpupyear', vid))
    call nc_err(nf90_get_var(fid, vid, casabal%fpupyear))
    call nc_err(nf90_inq_varid(fid, 'fpleachyear', vid))
    call nc_err(nf90_get_var(fid, vid, casabal%fpleachyear))
    call nc_err(nf90_inq_varid(fid, 'fplossyear', vid))
    call nc_err(nf90_get_var(fid, vid, casabal%fplossyear))
    call nc_err(nf90_inq_varid(fid, 'glaimon', vid))
    call nc_err(nf90_get_var(fid, vid, casabal%glaimon))
    call nc_err(nf90_inq_varid(fid, 'glaimonx', vid))
    call nc_err(nf90_get_var(fid, vid, casabal%glaimonx))
    call nc_err(nf90_inq_varid(fid, 'cplantlast', vid))
    call nc_err(nf90_get_var(fid, vid, casabal%cplantlast))
    call nc_err(nf90_inq_varid(fid, 'nplantlast', vid))
    call nc_err(nf90_get_var(fid, vid, casabal%nplantlast))
    call nc_err(nf90_inq_varid(fid, 'pplantlast', vid))
    call nc_err(nf90_get_var(fid, vid, casabal%pplantlast))
    call nc_err(nf90_inq_varid(fid, 'clitterlast', vid))
    call nc_err(nf90_get_var(fid, vid, casabal%clitterlast))
    call nc_err(nf90_inq_varid(fid, 'nlitterlast', vid))
    call nc_err(nf90_get_var(fid, vid, casabal%nlitterlast))
    call nc_err(nf90_inq_varid(fid, 'plitterlast', vid))
    call nc_err(nf90_get_var(fid, vid, casabal%plitterlast))
    call nc_err(nf90_inq_varid(fid, 'csoillast', vid))
    call nc_err(nf90_get_var(fid, vid, casabal%csoillast))
    call nc_err(nf90_inq_varid(fid, 'nsoillast', vid))
    call nc_err(nf90_get_var(fid, vid, casabal%nsoillast))
    call nc_err(nf90_inq_varid(fid, 'psoillast', vid))
    call nc_err(nf90_get_var(fid, vid, casabal%psoillast))
    call nc_err(nf90_inq_varid(fid, 'nsoilminlast', vid))
    call nc_err(nf90_get_var(fid, vid, casabal%nsoilminlast))
    call nc_err(nf90_inq_varid(fid, 'psoillablast', vid))
    call nc_err(nf90_get_var(fid, vid, casabal%psoillablast))
    call nc_err(nf90_inq_varid(fid, 'psoilsorblast', vid))
    call nc_err(nf90_get_var(fid, vid, casabal%psoilsorblast))
    call nc_err(nf90_inq_varid(fid, 'psoilocclast', vid))
    call nc_err(nf90_get_var(fid, vid, casabal%psoilocclast))
    call nc_err(nf90_inq_varid(fid, 'cbalance', vid))
    call nc_err(nf90_get_var(fid, vid, casabal%cbalance))
    call nc_err(nf90_inq_varid(fid, 'nbalance', vid))
    call nc_err(nf90_get_var(fid, vid, casabal%nbalance))
    call nc_err(nf90_inq_varid(fid, 'pbalance', vid))
    call nc_err(nf90_get_var(fid, vid, casabal%pbalance))
    call nc_err(nf90_inq_varid(fid, 'sumcbal', vid))
    call nc_err(nf90_get_var(fid, vid, casabal%sumcbal))
    call nc_err(nf90_inq_varid(fid, 'sumnbal', vid))
    call nc_err(nf90_get_var(fid, vid, casabal%sumnbal))
    call nc_err(nf90_inq_varid(fid, 'sumpbal', vid))
    call nc_err(nf90_get_var(fid, vid, casabal%sumpbal))
    call nc_err(nf90_inq_varid(fid, 'clabilelast', vid))
    call nc_err(nf90_get_var(fid, vid, casabal%clabilelast))

    ! close NetCDF file
    call nc_err(nf90_close(fid))

  end subroutine read_netcdf_casabal


  ! ------------------------------------------------------------------


  subroutine write_netcdf_casabiome(filename, casabiome)

    use netcdf, only: nf90_create, nf90_clobber, nf90_64bit_offset, &
         nf90_def_dim, nf90_def_var, nf90_int, nf90_double, &
         nf90_enddef, nf90_put_var, nf90_close
    use cable_def_types_mod, only: nc_err

    implicit none

    character(len=*), intent(in) :: filename
    type(casa_biome), intent(in) :: casabiome

    integer :: fid
    integer :: dimid1, dimid2, dimid3, dimid4
    integer :: i
    integer, dimension(ncasa_biome) :: vid

    ! create netCDF file
    call nc_err(nf90_create(trim(filename), ior(nf90_clobber, nf90_64bit_offset), fid))

    ! define dimensions
    ! vegetation types
    call nc_err(nf90_def_dim(fid, 'dim1', size(casabiome%plantrate, 1), dimid1))
    ! mplant
    call nc_err(nf90_def_dim(fid, 'dim2', size(casabiome%plantrate, 2), dimid2))
    ! mlitter
    call nc_err(nf90_def_dim(fid, 'dim3', size(casabiome%litterrate, 2), dimid3))
    ! msoil
    call nc_err(nf90_def_dim(fid, 'dim4', size(casabiome%soilrate, 2), dimid4))

    ! define variables
    i = 1
    ! define integer vector variables [dim1]
    call nc_err(nf90_def_var(fid, 'ivt2', nf90_int, &
         [dimid1], vid(i)), i)

    ! define double vector variables [dim1]
    call nc_err(nf90_def_var(fid, 'xkleafcoldmax', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'xkleafcoldexp', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'xkleafdrymax', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'xkleafdryexp', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'glaimax', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'glaimin', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'sla', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'ratiofrootleaf', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'kroot', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'krootlen', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'rootdepth', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'kuptake', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'kminn', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'kuplabp', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'kclabrate', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'xnpmax', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'q10soil', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'xkoptlitter', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'xkoptsoil', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'xkplab', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'xkpsorb', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'xkpocc', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'prodptase', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'costnpup', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'maxfinelitter', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'maxcwd', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'nintercept', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'nslope', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'la_to_sa', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'vcmax_scalar', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'disturbance_interval', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'damm_enzpool', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'damm_kmo2', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'damm_kmcp', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'damm_ea', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'damm_alpha', nf90_double, &
         [dimid1], vid(i)), i)

    ! define double array variables [dim1, dim2]
    call nc_err(nf90_def_var(fid, 'plantrate', nf90_double, &
         [dimid1, dimid2], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'rmplant', nf90_double, &
         [dimid1, dimid2], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fracnpptop', nf90_double, &
         [dimid1, dimid2], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fraclignin', nf90_double, &
         [dimid1, dimid2], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fraclabile', nf90_double, &
         [dimid1, dimid2], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'rationcplantmin', nf90_double, &
         [dimid1, dimid2], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'rationcplantmax', nf90_double, &
         [dimid1, dimid2], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'rationpplantmin', nf90_double, &
         [dimid1, dimid2], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'rationpplantmax', nf90_double, &
         [dimid1, dimid2], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fracligninplant', nf90_double, &
         [dimid1, dimid2], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'ftransnptol', nf90_double, &
         [dimid1, dimid2], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'ftranspptol', nf90_double, &
         [dimid1, dimid2], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'ratiopcplantmax', nf90_double, &
         [dimid1, dimid2], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'ratiopcplantmin', nf90_double, &
         [dimid1, dimid2], vid(i)), i)

    ! define double array variables [dim1, dim3]
    call nc_err(nf90_def_var(fid, 'litterrate', nf90_double, &
         [dimid1, dimid3], vid(i)), i)

    ! define double array variables [dim1, dim4]
    call nc_err(nf90_def_var(fid, 'soilrate', nf90_double, &
         [dimid1, dimid4], vid(i)), i)

    ! end define mode
    call nc_err(nf90_enddef(fid))

    ! put variables
    i = 1
    call nc_err(nf90_put_var(fid, vid(i), casabiome%ivt2), i)
    call nc_err(nf90_put_var(fid, vid(i), casabiome%xkleafcoldmax), i)
    call nc_err(nf90_put_var(fid, vid(i), casabiome%xkleafcoldexp), i)
    call nc_err(nf90_put_var(fid, vid(i), casabiome%xkleafdrymax), i)
    call nc_err(nf90_put_var(fid, vid(i), casabiome%xkleafdryexp), i)
    call nc_err(nf90_put_var(fid, vid(i), casabiome%glaimax), i)
    call nc_err(nf90_put_var(fid, vid(i), casabiome%glaimin), i)
    call nc_err(nf90_put_var(fid, vid(i), casabiome%sla), i)
    call nc_err(nf90_put_var(fid, vid(i), casabiome%ratiofrootleaf), i)
    call nc_err(nf90_put_var(fid, vid(i), casabiome%kroot), i)
    call nc_err(nf90_put_var(fid, vid(i), casabiome%krootlen), i)
    call nc_err(nf90_put_var(fid, vid(i), casabiome%rootdepth), i)
    call nc_err(nf90_put_var(fid, vid(i), casabiome%kuptake), i)
    call nc_err(nf90_put_var(fid, vid(i), casabiome%kminN), i)
    call nc_err(nf90_put_var(fid, vid(i), casabiome%KuplabP), i)
    call nc_err(nf90_put_var(fid, vid(i), casabiome%kclabrate), i)
    call nc_err(nf90_put_var(fid, vid(i), casabiome%xnpmax), i)
    call nc_err(nf90_put_var(fid, vid(i), casabiome%q10soil), i)
    call nc_err(nf90_put_var(fid, vid(i), casabiome%xkoptlitter), i)
    call nc_err(nf90_put_var(fid, vid(i), casabiome%xkoptsoil), i)
    call nc_err(nf90_put_var(fid, vid(i), casabiome%xkplab), i)
    call nc_err(nf90_put_var(fid, vid(i), casabiome%xkpsorb), i)
    call nc_err(nf90_put_var(fid, vid(i), casabiome%xkpocc), i)
    call nc_err(nf90_put_var(fid, vid(i), casabiome%prodptase), i)
    call nc_err(nf90_put_var(fid, vid(i), casabiome%costnpup), i)
    call nc_err(nf90_put_var(fid, vid(i), casabiome%maxfinelitter), i)
    call nc_err(nf90_put_var(fid, vid(i), casabiome%maxcwd), i)
    call nc_err(nf90_put_var(fid, vid(i), casabiome%nintercept), i)
    call nc_err(nf90_put_var(fid, vid(i), casabiome%nslope), i)
    call nc_err(nf90_put_var(fid, vid(i), casabiome%la_to_sa), i)
    call nc_err(nf90_put_var(fid, vid(i), casabiome%vcmax_scalar), i)
    call nc_err(nf90_put_var(fid, vid(i), casabiome%disturbance_interval), i)
    call nc_err(nf90_put_var(fid, vid(i), casabiome%DAMM_EnzPool), i)
    call nc_err(nf90_put_var(fid, vid(i), casabiome%DAMM_KMO2), i)
    call nc_err(nf90_put_var(fid, vid(i), casabiome%DAMM_KMcp), i)
    call nc_err(nf90_put_var(fid, vid(i), casabiome%DAMM_Ea), i)
    call nc_err(nf90_put_var(fid, vid(i), casabiome%DAMM_alpha), i)
    call nc_err(nf90_put_var(fid, vid(i), casabiome%plantrate), i)
    call nc_err(nf90_put_var(fid, vid(i), casabiome%rmplant), i)
    call nc_err(nf90_put_var(fid, vid(i), casabiome%fracnpptoP), i)
    call nc_err(nf90_put_var(fid, vid(i), casabiome%fraclignin), i)
    call nc_err(nf90_put_var(fid, vid(i), casabiome%fraclabile), i)
    call nc_err(nf90_put_var(fid, vid(i), casabiome%ratioNCplantmin), i)
    call nc_err(nf90_put_var(fid, vid(i), casabiome%ratioNCplantmax), i)
    call nc_err(nf90_put_var(fid, vid(i), casabiome%ratioNPplantmin), i)
    call nc_err(nf90_put_var(fid, vid(i), casabiome%ratioNPplantmax), i)
    call nc_err(nf90_put_var(fid, vid(i), casabiome%fracLigninplant), i)
    call nc_err(nf90_put_var(fid, vid(i), casabiome%ftransNPtoL), i)
    call nc_err(nf90_put_var(fid, vid(i), casabiome%ftransPPtoL), i)
    call nc_err(nf90_put_var(fid, vid(i), casabiome%ratioPcplantmax), i)
    call nc_err(nf90_put_var(fid, vid(i), casabiome%ratioPcplantmin), i)
    call nc_err(nf90_put_var(fid, vid(i), casabiome%litterrate), i)
    call nc_err(nf90_put_var(fid, vid(i), casabiome%soilrate), i)

    ! close NetCDF file
    call nc_err(nf90_close(fid))

  end subroutine write_netcdf_casabiome


  subroutine write_netcdf_casapool(filename, casapool)

    use netcdf, only: nf90_create, nf90_clobber, nf90_64bit_offset, &
         nf90_def_dim, nf90_def_var, nf90_double, &
         nf90_enddef, nf90_put_var, nf90_close
    use cable_def_types_mod, only: nc_err

    implicit none

    character(len=*), intent(in) :: filename
    type(casa_pool),  intent(in) :: casapool

    integer :: fid
    integer :: dimid1, dimid2, dimid3, dimid4
    integer :: i
    integer, dimension(ncasa_pool) :: vid

    ! create netCDF file
    call nc_err(nf90_create(trim(filename), ior(nf90_clobber, nf90_64bit_offset), fid))

    ! define dimensions
    ! land
    call nc_err(nf90_def_dim(fid, 'dim1', size(casapool%Cplant, 1), dimid1))
    ! mplant
    call nc_err(nf90_def_dim(fid, 'dim2', size(casapool%Cplant, 2), dimid2))
    ! mlitter
    call nc_err(nf90_def_dim(fid, 'dim3', size(casapool%Clitter, 2), dimid3))
    ! msoil
    call nc_err(nf90_def_dim(fid, 'dim4', size(casapool%Csoil, 2), dimid4))

    ! define variables
    i = 1
    ! define double vector variables [dim1]
    call nc_err(nf90_def_var(fid, 'clabile', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'dclabiledt', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'nsoilmin', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'psoillab', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'psoilsorb', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'psoilocc', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'dnsoilmindt', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'dpsoillabdt', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'dpsoilsorbdt', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'dpsoiloccdt', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'ctot_0', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'ctot', nf90_double, &
         [dimid1], vid(i)), i)

    ! define double array variables [dim1, dim2]
    call nc_err(nf90_def_var(fid, 'cplant', nf90_double, &
         [dimid1, dimid2], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'nplant', nf90_double, &
         [dimid1, dimid2], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'pplant', nf90_double, &
         [dimid1, dimid2], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'dcplantdt', nf90_double, &
         [dimid1, dimid2], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'dnplantdt', nf90_double, &
         [dimid1, dimid2], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'dpplantdt', nf90_double, &
         [dimid1, dimid2], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'rationcplant', nf90_double, &
         [dimid1, dimid2], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'rationpplant', nf90_double, &
         [dimid1, dimid2], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'ratiopcplant', nf90_double, &
         [dimid1, dimid2], vid(i)), i)

    ! define double array variables [dim1, dim3]
    call nc_err(nf90_def_var(fid, 'clitter', nf90_double, &
         [dimid1, dimid3], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'nlitter', nf90_double, &
         [dimid1, dimid3], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'plitter', nf90_double, &
         [dimid1, dimid3], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'dclitterdt', nf90_double, &
         [dimid1, dimid3], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'dnlitterdt', nf90_double, &
         [dimid1, dimid3], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'dplitterdt', nf90_double, &
         [dimid1, dimid3], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'rationclitter', nf90_double, &
         [dimid1, dimid3], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'rationplitter', nf90_double, &
         [dimid1, dimid3], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'ratiopclitter', nf90_double, &
         [dimid1, dimid3], vid(i)), i)

    ! define double array variables [dim1, dim4]
    call nc_err(nf90_def_var(fid, 'csoil', nf90_double, &
         [dimid1, dimid4], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'nsoil', nf90_double, &
         [dimid1, dimid4], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'psoil', nf90_double, &
         [dimid1, dimid4], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'dcsoildt', nf90_double, &
         [dimid1, dimid4], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'dnsoildt', nf90_double, &
         [dimid1, dimid4], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'dpsoildt', nf90_double, &
         [dimid1, dimid4], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'rationcsoil', nf90_double, &
         [dimid1, dimid4], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'rationpsoil', nf90_double, &
         [dimid1, dimid4], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'rationcsoilnew', nf90_double, &
         [dimid1, dimid4], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'rationcsoilmin', nf90_double, &
         [dimid1, dimid4], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'rationcsoilmax', nf90_double, &
         [dimid1, dimid4], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'ratiopcsoil', nf90_double, &
         [dimid1, dimid4], vid(i)), i)

    ! end define mode
    call nc_err(nf90_enddef(fid))

    ! put variables
    i = 1
    call nc_err(nf90_put_var(fid, vid(i), casapool%Clabile), i)
    call nc_err(nf90_put_var(fid, vid(i), casapool%dClabiledt), i)
    call nc_err(nf90_put_var(fid, vid(i), casapool%Nsoilmin), i)
    call nc_err(nf90_put_var(fid, vid(i), casapool%Psoillab), i)
    call nc_err(nf90_put_var(fid, vid(i), casapool%Psoilsorb), i)
    call nc_err(nf90_put_var(fid, vid(i), casapool%Psoilocc), i)
    call nc_err(nf90_put_var(fid, vid(i), casapool%dNsoilmindt), i)
    call nc_err(nf90_put_var(fid, vid(i), casapool%dPsoillabdt), i)
    call nc_err(nf90_put_var(fid, vid(i), casapool%dPsoilsorbdt), i)
    call nc_err(nf90_put_var(fid, vid(i), casapool%dPsoiloccdt), i)
    call nc_err(nf90_put_var(fid, vid(i), casapool%Ctot_0), i)
    call nc_err(nf90_put_var(fid, vid(i), casapool%Ctot), i)
    call nc_err(nf90_put_var(fid, vid(i), casapool%Cplant), i)
    call nc_err(nf90_put_var(fid, vid(i), casapool%Nplant), i)
    call nc_err(nf90_put_var(fid, vid(i), casapool%Pplant), i)
    call nc_err(nf90_put_var(fid, vid(i), casapool%dCplantdt), i)
    call nc_err(nf90_put_var(fid, vid(i), casapool%dNplantdt), i)
    call nc_err(nf90_put_var(fid, vid(i), casapool%dPplantdt), i)
    call nc_err(nf90_put_var(fid, vid(i), casapool%ratioNCplant), i)
    call nc_err(nf90_put_var(fid, vid(i), casapool%ratioNPplant), i)
    call nc_err(nf90_put_var(fid, vid(i), casapool%ratioPCplant), i)
    call nc_err(nf90_put_var(fid, vid(i), casapool%Clitter), i)
    call nc_err(nf90_put_var(fid, vid(i), casapool%Nlitter), i)
    call nc_err(nf90_put_var(fid, vid(i), casapool%Plitter), i)
    call nc_err(nf90_put_var(fid, vid(i), casapool%dClitterdt), i)
    call nc_err(nf90_put_var(fid, vid(i), casapool%dNlitterdt), i)
    call nc_err(nf90_put_var(fid, vid(i), casapool%dPlitterdt), i)
    call nc_err(nf90_put_var(fid, vid(i), casapool%ratioNClitter), i)
    call nc_err(nf90_put_var(fid, vid(i), casapool%ratioNPlitter), i)
    call nc_err(nf90_put_var(fid, vid(i), casapool%ratioPClitter), i)
    call nc_err(nf90_put_var(fid, vid(i), casapool%Csoil), i)
    call nc_err(nf90_put_var(fid, vid(i), casapool%Nsoil), i)
    call nc_err(nf90_put_var(fid, vid(i), casapool%Psoil), i)
    call nc_err(nf90_put_var(fid, vid(i), casapool%dCsoildt), i)
    call nc_err(nf90_put_var(fid, vid(i), casapool%dNsoildt), i)
    call nc_err(nf90_put_var(fid, vid(i), casapool%dPsoildt), i)
    call nc_err(nf90_put_var(fid, vid(i), casapool%ratioNCsoil), i)
    call nc_err(nf90_put_var(fid, vid(i), casapool%ratioNPsoil), i)
    call nc_err(nf90_put_var(fid, vid(i), casapool%ratioNCsoilnew), i)
    call nc_err(nf90_put_var(fid, vid(i), casapool%ratioNCsoilmin), i)
    call nc_err(nf90_put_var(fid, vid(i), casapool%ratioNCsoilmax), i)
    call nc_err(nf90_put_var(fid, vid(i), casapool%ratioPCsoil), i)

    ! close NetCDF file
    call nc_err(nf90_close(fid))

  end subroutine write_netcdf_casapool


  subroutine write_netcdf_casaflux(filename, casaflux)

    use netcdf, only: nf90_create, nf90_clobber, nf90_64bit_offset, &
         nf90_def_dim, nf90_def_var, nf90_double, &
         nf90_enddef, nf90_put_var, nf90_close
    use cable_def_types_mod, only: nc_err

    implicit none

    character(len=*), intent(in) :: filename
    type(casa_flux),  intent(in) :: casaflux

    integer :: fid
    integer :: dimid1, dimid2, dimid3, dimid4
    integer :: i
    integer, dimension(ncasa_flux) :: vid

    ! create netCDF file
    call nc_err(nf90_create(trim(filename), ior(nf90_clobber, nf90_64bit_offset), fid))

    ! define dimensions
    ! land
    call nc_err(nf90_def_dim(fid, 'dim1', size(casaflux%kplant, 1), dimid1))
    ! mplant
    call nc_err(nf90_def_dim(fid, 'dim2', size(casaflux%kplant, 2), dimid2))
    ! mlitter
    call nc_err(nf90_def_dim(fid, 'dim3', size(casaflux%klitter, 2), dimid3))
    ! msoil
    call nc_err(nf90_def_dim(fid, 'dim4', size(casaflux%ksoil, 2), dimid4))

    ! define variables
    i = 1
    ! define double vector variables [dim1]
    call nc_err(nf90_def_var(fid, 'cgpp', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'cnpp', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'crp', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'crgplant', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'nminfix', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'nminuptake', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'plabuptake', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'clabloss', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fracclabile', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'cnep', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'crsoil', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'nmindep', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'nminloss', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'nminleach', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'nupland', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'nlittermin', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'nsmin', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'nsimm', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'nsnet', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fnminloss', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fnminleach', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'pdep', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'pwea', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'pleach', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'ploss', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'pupland', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'plittermin', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'psmin', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'psimm', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'psnet', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fpleach', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'kplab', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'kpsorb', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'kpocc', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'kmlabp', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'psorbmax', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'stemnpp', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'frac_sapwood', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'sapwood_area', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fharvest', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'charvest', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'nharvest', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'pharvest', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fcrop', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'cplant_turnover_disturbance', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'cplant_turnover_crowding', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'cplant_turnover_resource_limitation', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fluxctoco2_plant_fire', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fluxctoco2_litter_fire', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fluxntoatm_fire', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fluxctohwp', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fluxntohwp', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fluxptohwp', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fluxctoclear', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fluxntoclear', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fluxptoclear', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'ctransferluc', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fluxctoco2', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fluxfromptoharvest', nf90_double, &
         [dimid1], vid(i)), i)

    ! define double array variables [dim1, dim2]
    call nc_err(nf90_def_var(fid, 'fraccalloc', nf90_double, &
         [dimid1, dimid2], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fracnalloc', nf90_double, &
         [dimid1, dimid2], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fracpalloc', nf90_double, &
         [dimid1, dimid2], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'kplant', nf90_double, &
         [dimid1, dimid2], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'crmplant', nf90_double, &
         [dimid1, dimid2], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'cplant_turnover', nf90_double, &
         [dimid1, dimid2], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'kplant_fire', nf90_double, &
         [dimid1, dimid2], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'kplant_tot', nf90_double, &
         [dimid1, dimid2], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fluxfromptoco2_fire', nf90_double, &
         [dimid1, dimid2], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fluxfromptoco2', nf90_double, &
         [dimid1, dimid2], vid(i)), i)

    ! define double array variables [dim1, dim3]
    call nc_err(nf90_def_var(fid, 'klitter', nf90_double, &
         [dimid1, dimid3], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fromltoco2', nf90_double, &
         [dimid1, dimid3], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'klitter_fire', nf90_double, &
         [dimid1, dimid3], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'klitter_tot', nf90_double, &
         [dimid1, dimid3], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fluxctolitter', nf90_double, &
         [dimid1, dimid3], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fluxntolitter', nf90_double, &
         [dimid1, dimid3], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fluxptolitter', nf90_double, &
         [dimid1, dimid3], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fluxfromltoco2_fire', nf90_double, &
         [dimid1, dimid3], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fluxfromltoco2', nf90_double, &
         [dimid1, dimid3], vid(i)), i)

    ! define double array variables [dim1, dim4]
    call nc_err(nf90_def_var(fid, 'ksoil', nf90_double, &
         [dimid1, dimid4], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fromstoco2', nf90_double, &
         [dimid1, dimid4], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fluxctosoil', nf90_double, &
         [dimid1, dimid4], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fluxntosoil', nf90_double, &
         [dimid1, dimid4], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fluxptosoil', nf90_double, &
         [dimid1, dimid4], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fluxfromstoco2', nf90_double, &
         [dimid1, dimid4], vid(i)), i)

    ! define double array variables [dim1, dim2, dim3]
    call nc_err(nf90_def_var(fid, 'fluxfromptol', nf90_double, &
         [dimid1, dimid2, dimid3], vid(i)), i)

    ! define double array variables [dim1, dim3, dim2]
    call nc_err(nf90_def_var(fid, 'fromptol', nf90_double, &
         [dimid1, dimid3, dimid2], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fromptol_fire', nf90_double, &
         [dimid1, dimid3, dimid2], vid(i)), i)

    ! define double array variables [dim1, dim3, dim4]
    call nc_err(nf90_def_var(fid, 'fluxfromltos', nf90_double, &
         [dimid1, dimid3, dimid4], vid(i)), i)

    ! define double array variables [dim1, dim4, dim3]
    call nc_err(nf90_def_var(fid, 'fromltos', nf90_double, &
         [dimid1, dimid4, dimid3], vid(i)), i)

    ! define double array variables [dim1, dim4, dim4]
    call nc_err(nf90_def_var(fid, 'fromstos', nf90_double, &
         [dimid1, dimid4, dimid4], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fluxfromstos', nf90_double, &
         [dimid1, dimid4, dimid4], vid(i)), i)

    ! end define mode
    call nc_err(nf90_enddef(fid))

    ! put variables
    i = 1
    call nc_err(nf90_put_var(fid, vid(i), casaflux%Cgpp), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%Cnpp), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%Crp), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%Crgplant), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%Nminfix), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%Nminuptake), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%Plabuptake), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%Clabloss), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%fracClabile), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%Cnep), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%Crsoil), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%Nmindep), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%Nminloss), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%Nminleach), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%Nupland), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%Nlittermin), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%Nsmin), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%Nsimm), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%Nsnet), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%fNminloss), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%fNminleach), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%Pdep), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%Pwea), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%Pleach), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%Ploss), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%Pupland), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%Plittermin), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%Psmin), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%Psimm), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%Psnet), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%fPleach), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%kplab), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%kpsorb), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%kpocc), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%kmlabP), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%Psorbmax), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%stemnpp), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%frac_sapwood), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%sapwood_area), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%fharvest), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%Charvest), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%Nharvest), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%Pharvest), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%fcrop), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%Cplant_turnover_disturbance), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%Cplant_turnover_crowding), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%Cplant_turnover_resource_limitation), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%FluxCtoCO2_plant_fire), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%FluxCtoCO2_litter_fire), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%FluxNtoAtm_fire), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%FluxCtohwp), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%FluxNtohwp), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%FluxPtohwp), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%FluxCtoclear), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%FluxNtoclear), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%FluxPtoclear), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%CtransferLUC), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%FluxCtoco2), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%FluxFromPtoHarvest), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%fracCalloc), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%fracNalloc), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%fracPalloc), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%kplant), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%Crmplant), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%Cplant_turnover), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%kplant_fire), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%kplant_tot), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%fluxfromPtoCO2_fire), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%FluxFromPtoCO2), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%klitter), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%fromLtoCO2), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%klitter_fire), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%klitter_tot), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%FluxCtolitter), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%FluxNtolitter), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%FluxPtolitter), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%fluxfromLtoCO2_fire), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%FluxFromLtoCO2), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%ksoil), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%fromStoCO2), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%FluxCtosoil), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%FluxNtosoil), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%FluxPtosoil), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%FluxFromStoCO2), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%FluxFromPtoL), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%fromPtoL), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%fromPtoL_fire), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%FluxFromLtoS), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%fromLtoS), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%fromStoS), i)
    call nc_err(nf90_put_var(fid, vid(i), casaflux%FluxFromStoS), i)

    ! close NetCDF file
    call nc_err(nf90_close(fid))

  end subroutine write_netcdf_casaflux


  subroutine write_netcdf_casamet(filename, casamet)

    use netcdf, only: nf90_create, nf90_clobber, nf90_64bit_offset, &
         nf90_def_dim, nf90_def_var, nf90_int, nf90_double, &
         nf90_enddef, nf90_put_var, nf90_close
    use cable_def_types_mod, only: nc_err

    implicit none

    character(len=*), intent(in) :: filename
    type(casa_met),   intent(in) :: casamet

    integer :: fid
    integer :: dimid1, dimid2, dimid3
    integer :: i
    integer, dimension(ncasa_met) :: vid

    ! create netCDF file
    call nc_err(nf90_create(trim(filename), ior(nf90_clobber, nf90_64bit_offset), fid))

    ! define dimensions
    ! land
    call nc_err(nf90_def_dim(fid, 'dim1', size(casamet%Tsoil, 1), dimid1))
    ! soil layer
    call nc_err(nf90_def_dim(fid, 'dim2', size(casamet%Tsoil, 2), dimid2))
    ! days of year = 365
    call nc_err(nf90_def_dim(fid, 'dim3', size(casamet%Tairkspin, 2), dimid3))

    ! define variables
    i = 1
    ! define integer vector variables [dim1]
    call nc_err(nf90_def_var(fid, 'lnonwood', nf90_int, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'iveg2', nf90_int, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'ijgcm', nf90_int, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'isorder', nf90_int, &
         [dimid1], vid(i)), i)

    ! define double vector variables [dim1]
    call nc_err(nf90_def_var(fid, 'glai', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'tairk', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'precip', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'tsoilavg', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'moistavg', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'btran', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'lat', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'lon', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'areacell', nf90_double, &
         [dimid1], vid(i)), i)

    ! define double array variables [dim1, dim2]
    call nc_err(nf90_def_var(fid, 'tsoil', nf90_double, &
         [dimid1, dimid2], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'moist', nf90_double, &
         [dimid1, dimid2], vid(i)), i)

    ! define integer array variables [dim1, dim3]
    call nc_err(nf90_def_var(fid, 'dslr_spin', nf90_int, &
         [dimid1, dimid3], vid(i)), i)

    ! define double array variables [dim1, dim3]
    call nc_err(nf90_def_var(fid, 'tairkspin', nf90_double, &
         [dimid1, dimid3], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'cgppspin', nf90_double, &
         [dimid1, dimid3], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'crmplantspin_1', nf90_double, &
         [dimid1, dimid3], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'crmplantspin_2', nf90_double, &
         [dimid1, dimid3], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'crmplantspin_3', nf90_double, &
         [dimid1, dimid3], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'tsoilspin_1', nf90_double, &
         [dimid1, dimid3], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'tsoilspin_2', nf90_double, &
         [dimid1, dimid3], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'tsoilspin_3', nf90_double, &
         [dimid1, dimid3], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'tsoilspin_4', nf90_double, &
         [dimid1, dimid3], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'tsoilspin_5', nf90_double, &
         [dimid1, dimid3], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'tsoilspin_6', nf90_double, &
         [dimid1, dimid3], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'moistspin_1', nf90_double, &
         [dimid1, dimid3], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'moistspin_2', nf90_double, &
         [dimid1, dimid3], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'moistspin_3', nf90_double, &
         [dimid1, dimid3], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'moistspin_4', nf90_double, &
         [dimid1, dimid3], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'moistspin_5', nf90_double, &
         [dimid1, dimid3], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'moistspin_6', nf90_double, &
         [dimid1, dimid3], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'mtempspin', nf90_double, &
         [dimid1, dimid3], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'frecspin', nf90_double, &
         [dimid1, dimid3], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'can12spin', nf90_double, &
         [dimid1, dimid3], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'can13spin', nf90_double, &
         [dimid1, dimid3], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'dprecip_spin', nf90_double, &
         [dimid1, dimid3], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'aprecip_av20_spin', nf90_double, &
         [dimid1, dimid3], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'du10_max_spin', nf90_double, &
         [dimid1, dimid3], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'drhum_spin', nf90_double, &
         [dimid1, dimid3], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'dtemp_max_spin', nf90_double, &
         [dimid1, dimid3], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'dtemp_min_spin', nf90_double, &
         [dimid1, dimid3], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'kbdi_spin', nf90_double, &
         [dimid1, dimid3], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'd_macarthur_spin', nf90_double, &
         [dimid1, dimid3], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'ffdi_spin', nf90_double, &
         [dimid1, dimid3], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'last_precip_spin', nf90_double, &
         [dimid1, dimid3], vid(i)), i)

    ! end define mode
    call nc_err(nf90_enddef(fid))

    ! put variables
    i = 1
    call nc_err(nf90_put_var(fid, vid(i), casamet%lnonwood), i)
    call nc_err(nf90_put_var(fid, vid(i), casamet%iveg2), i)
    call nc_err(nf90_put_var(fid, vid(i), casamet%ijgcm), i)
    call nc_err(nf90_put_var(fid, vid(i), casamet%isorder), i)
    call nc_err(nf90_put_var(fid, vid(i), casamet%glai), i)
    call nc_err(nf90_put_var(fid, vid(i), casamet%Tairk), i)
    call nc_err(nf90_put_var(fid, vid(i), casamet%precip), i)
    call nc_err(nf90_put_var(fid, vid(i), casamet%tsoilavg), i)
    call nc_err(nf90_put_var(fid, vid(i), casamet%moistavg), i)
    call nc_err(nf90_put_var(fid, vid(i), casamet%btran), i)
    call nc_err(nf90_put_var(fid, vid(i), casamet%lat), i)
    call nc_err(nf90_put_var(fid, vid(i), casamet%lon), i)
    call nc_err(nf90_put_var(fid, vid(i), casamet%areacell), i)
    call nc_err(nf90_put_var(fid, vid(i), casamet%Tsoil), i)
    call nc_err(nf90_put_var(fid, vid(i), casamet%moist), i)
    call nc_err(nf90_put_var(fid, vid(i), casamet%DSLR_spin), i)
    call nc_err(nf90_put_var(fid, vid(i), casamet%Tairkspin), i)
    call nc_err(nf90_put_var(fid, vid(i), casamet%cgppspin), i)
    call nc_err(nf90_put_var(fid, vid(i), casamet%crmplantspin_1), i)
    call nc_err(nf90_put_var(fid, vid(i), casamet%crmplantspin_2), i)
    call nc_err(nf90_put_var(fid, vid(i), casamet%crmplantspin_3), i)
    call nc_err(nf90_put_var(fid, vid(i), casamet%Tsoilspin_1), i)
    call nc_err(nf90_put_var(fid, vid(i), casamet%Tsoilspin_2), i)
    call nc_err(nf90_put_var(fid, vid(i), casamet%Tsoilspin_3), i)
    call nc_err(nf90_put_var(fid, vid(i), casamet%Tsoilspin_4), i)
    call nc_err(nf90_put_var(fid, vid(i), casamet%Tsoilspin_5), i)
    call nc_err(nf90_put_var(fid, vid(i), casamet%Tsoilspin_6), i)
    call nc_err(nf90_put_var(fid, vid(i), casamet%moistspin_1), i)
    call nc_err(nf90_put_var(fid, vid(i), casamet%moistspin_2), i)
    call nc_err(nf90_put_var(fid, vid(i), casamet%moistspin_3), i)
    call nc_err(nf90_put_var(fid, vid(i), casamet%moistspin_4), i)
    call nc_err(nf90_put_var(fid, vid(i), casamet%moistspin_5), i)
    call nc_err(nf90_put_var(fid, vid(i), casamet%moistspin_6), i)
    call nc_err(nf90_put_var(fid, vid(i), casamet%mtempspin), i)
    call nc_err(nf90_put_var(fid, vid(i), casamet%frecspin), i)
    call nc_err(nf90_put_var(fid, vid(i), casamet%cAn12spin), i)
    call nc_err(nf90_put_var(fid, vid(i), casamet%cAn13spin), i)
    call nc_err(nf90_put_var(fid, vid(i), casamet%dprecip_spin), i)
    call nc_err(nf90_put_var(fid, vid(i), casamet%aprecip_av20_spin), i)
    call nc_err(nf90_put_var(fid, vid(i), casamet%du10_max_spin), i)
    call nc_err(nf90_put_var(fid, vid(i), casamet%drhum_spin), i)
    call nc_err(nf90_put_var(fid, vid(i), casamet%dtemp_max_spin), i)
    call nc_err(nf90_put_var(fid, vid(i), casamet%dtemp_min_spin), i)
    call nc_err(nf90_put_var(fid, vid(i), casamet%KBDI_spin), i)
    call nc_err(nf90_put_var(fid, vid(i), casamet%D_MacArthur_spin), i)
    call nc_err(nf90_put_var(fid, vid(i), casamet%FFDI_spin), i)
    call nc_err(nf90_put_var(fid, vid(i), casamet%last_precip_spin), i)

    ! close NetCDF file
    call nc_err(nf90_close(fid))

  end subroutine write_netcdf_casamet


  subroutine write_netcdf_casabal(filename, casabal)

    use netcdf, only: nf90_create, nf90_clobber, nf90_64bit_offset, &
         nf90_def_dim, nf90_def_var, nf90_float, nf90_double, &
         nf90_enddef, nf90_put_var, nf90_close
    use cable_def_types_mod, only: nc_err

    implicit none

    character(len=*),   intent(in) :: filename
    type(casa_balance), intent(in) :: casabal

    integer :: fid
    integer :: dimid1, dimid2, dimid3, dimid4, dimid5
    integer :: i
    integer, dimension(ncasa_bal) :: vid

    ! create netCDF file
    call nc_err(nf90_create(trim(filename), ior(nf90_clobber, nf90_64bit_offset), fid))

    ! define dimensions
    ! land
    call nc_err(nf90_def_dim(fid, 'dim1', size(casabal%cplantlast, 1), dimid1))
    ! mplant
    call nc_err(nf90_def_dim(fid, 'dim2', size(casabal%cplantlast, 2), dimid2))
    ! mlitter
    call nc_err(nf90_def_dim(fid, 'dim3', size(casabal%clitterlast, 2), dimid3))
    ! msoil
    call nc_err(nf90_def_dim(fid, 'dim4', size(casabal%csoillast, 2), dimid4))
    ! number of months = 12
    call nc_err(nf90_def_dim(fid, 'dim5', size(casabal%glaimon, 2), dimid5))

    ! define variables
    i = 1
    ! define double vector variables [dim1]
    call nc_err(nf90_def_var(fid, 'fcgppyear', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fcnppyear', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fcrpyear', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fcrmleafyear', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fcrmwoodyear', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fcrmrootyear', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fcrgrowyear', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fcrsyear', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fcneeyear', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fndepyear', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fnfixyear', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fnsnetyear', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fnupyear', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fnleachyear', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fnlossyear', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fpweayear', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fpdustyear', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fpsnetyear', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fpupyear', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fpleachyear', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fplossyear', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'dcdtyear', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'laimax', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'cleafmean', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'crootmean', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'nsoilminlast', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'psoillablast', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'psoilsorblast', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'psoilocclast', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'cbalance', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'nbalance', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'pbalance', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'sumcbal', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'sumnbal', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'sumpbal', nf90_double, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'clabilelast', nf90_double, &
         [dimid1], vid(i)), i)

    ! define double array variables [dim1, dim2]
    call nc_err(nf90_def_var(fid, 'cplantlast', nf90_double, &
         [dimid1, dimid2], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'nplantlast', nf90_double, &
         [dimid1, dimid2], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'pplantlast', nf90_double, &
         [dimid1, dimid2], vid(i)), i)

    ! define double array variables [dim1, dim3]
    call nc_err(nf90_def_var(fid, 'clitterlast', nf90_double, &
         [dimid1, dimid3], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'nlitterlast', nf90_double, &
         [dimid1, dimid3], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'plitterlast', nf90_double, &
         [dimid1, dimid3], vid(i)), i)

    ! define double array variables [dim1, dim4]
    call nc_err(nf90_def_var(fid, 'csoillast', nf90_double, &
         [dimid1, dimid4], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'nsoillast', nf90_double, &
         [dimid1, dimid4], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'psoillast', nf90_double, &
         [dimid1, dimid4], vid(i)), i)

    ! define double array variables [dim1, dim5]
    call nc_err(nf90_def_var(fid, 'glaimon', nf90_double, &
         [dimid1, dimid5], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'glaimonx', nf90_double, &
         [dimid1, dimid5], vid(i)), i)

    ! end define mode
    call nc_err(nf90_enddef(fid))

    ! put variables
    i = 1
    call nc_err(nf90_put_var(fid, vid(i), casabal%FCgppyear), i)
    call nc_err(nf90_put_var(fid, vid(i), casabal%FCnppyear), i)
    call nc_err(nf90_put_var(fid, vid(i), casabal%FCrpyear), i)
    call nc_err(nf90_put_var(fid, vid(i), casabal%FCrmleafyear), i)
    call nc_err(nf90_put_var(fid, vid(i), casabal%FCrmwoodyear), i)
    call nc_err(nf90_put_var(fid, vid(i), casabal%FCrmrootyear), i)
    call nc_err(nf90_put_var(fid, vid(i), casabal%FCrgrowyear), i)
    call nc_err(nf90_put_var(fid, vid(i), casabal%FCrsyear), i)
    call nc_err(nf90_put_var(fid, vid(i), casabal%FCneeyear), i)
    call nc_err(nf90_put_var(fid, vid(i), casabal%FNdepyear), i)
    call nc_err(nf90_put_var(fid, vid(i), casabal%FNfixyear), i)
    call nc_err(nf90_put_var(fid, vid(i), casabal%FNsnetyear), i)
    call nc_err(nf90_put_var(fid, vid(i), casabal%FNupyear), i)
    call nc_err(nf90_put_var(fid, vid(i), casabal%FNleachyear), i)
    call nc_err(nf90_put_var(fid, vid(i), casabal%FNlossyear), i)
    call nc_err(nf90_put_var(fid, vid(i), casabal%FPweayear), i)
    call nc_err(nf90_put_var(fid, vid(i), casabal%FPdustyear), i)
    call nc_err(nf90_put_var(fid, vid(i), casabal%FPsnetyear), i)
    call nc_err(nf90_put_var(fid, vid(i), casabal%FPupyear), i)
    call nc_err(nf90_put_var(fid, vid(i), casabal%FPleachyear), i)
    call nc_err(nf90_put_var(fid, vid(i), casabal%FPlossyear), i)
    call nc_err(nf90_put_var(fid, vid(i), casabal%dCdtyear), i)
    call nc_err(nf90_put_var(fid, vid(i), casabal%LAImax), i)
    call nc_err(nf90_put_var(fid, vid(i), casabal%Cleafmean), i)
    call nc_err(nf90_put_var(fid, vid(i), casabal%Crootmean), i)
    call nc_err(nf90_put_var(fid, vid(i), casabal%nsoilminlast), i)
    call nc_err(nf90_put_var(fid, vid(i), casabal%psoillablast), i)
    call nc_err(nf90_put_var(fid, vid(i), casabal%psoilsorblast), i)
    call nc_err(nf90_put_var(fid, vid(i), casabal%psoilocclast), i)
    call nc_err(nf90_put_var(fid, vid(i), casabal%cbalance), i)
    call nc_err(nf90_put_var(fid, vid(i), casabal%nbalance), i)
    call nc_err(nf90_put_var(fid, vid(i), casabal%pbalance), i)
    call nc_err(nf90_put_var(fid, vid(i), casabal%sumcbal), i)
    call nc_err(nf90_put_var(fid, vid(i), casabal%sumnbal), i)
    call nc_err(nf90_put_var(fid, vid(i), casabal%sumpbal), i)
    call nc_err(nf90_put_var(fid, vid(i), casabal%clabilelast), i)
    call nc_err(nf90_put_var(fid, vid(i), casabal%cplantlast), i)
    call nc_err(nf90_put_var(fid, vid(i), casabal%nplantlast), i)
    call nc_err(nf90_put_var(fid, vid(i), casabal%pplantlast), i)
    call nc_err(nf90_put_var(fid, vid(i), casabal%clitterlast), i)
    call nc_err(nf90_put_var(fid, vid(i), casabal%nlitterlast), i)
    call nc_err(nf90_put_var(fid, vid(i), casabal%plitterlast), i)
    call nc_err(nf90_put_var(fid, vid(i), casabal%csoillast), i)
    call nc_err(nf90_put_var(fid, vid(i), casabal%nsoillast), i)
    call nc_err(nf90_put_var(fid, vid(i), casabal%psoillast), i)
    call nc_err(nf90_put_var(fid, vid(i), casabal%glaimon), i)
    call nc_err(nf90_put_var(fid, vid(i), casabal%glaimonx), i)

    ! close NetCDF file
    call nc_err(nf90_close(fid))

  end subroutine write_netcdf_casabal


end module casavariable


! ------------------------------------------------------------------


module phenvariable

  use cable_def_types_mod, only: mvtype, r_2
  use casadimension, only: mdyear, mphase

  implicit none

  private

  ! type
  public :: phen_variable

  ! routines on type
  public :: alloc_phenvariable
  public :: dealloc_phenvariable
  public :: print_phenvariable
  public :: read_netcdf_phen_var
  public :: write_netcdf_phen_var
  public :: zero_phenvariable

  ! number of variables in type definitions
  ! used in write_netcdf and in MPI code
  integer, parameter, public :: ncasa_phen = 10

  type phen_variable
     integer,   dimension(:),   pointer :: phase => null()
     real(r_2), dimension(:),   pointer :: TKshed => null()
     integer,   dimension(:,:), pointer :: doyphase => null()
     ! fraction of max LAI
     real,      dimension(:),   pointer :: phen => null()
     ! annual leaf on sum
     real,      dimension(:),   pointer :: aphen => null()
     integer,   dimension(:,:), pointer :: phasespin => null()
     integer,   dimension(:,:), pointer :: doyphasespin_1 => null()
     integer,   dimension(:,:), pointer :: doyphasespin_2 => null()
     integer,   dimension(:,:), pointer :: doyphasespin_3 => null()
     integer,   dimension(:,:), pointer :: doyphasespin_4 => null()
  end type phen_variable

contains

  subroutine alloc_phenvariable(phen,arraysize)

    use cable_def_types_mod, only: mvtype
    use casadimension, only: mdyear, mphase

    implicit none

    type(phen_variable), intent(inout) :: phen
    integer,             intent(in) :: arraysize

    allocate(phen%Tkshed(mvtype))
    allocate( &
         phen%phase(arraysize), &
         phen%doyphase(arraysize,mphase), &
         phen%phen(arraysize), &
         phen%aphen(arraysize), &
         phen%phasespin(arraysize,mdyear), &
         phen%doyphasespin_1(arraysize,mdyear), &
         phen%doyphasespin_2(arraysize,mdyear), &
         phen%doyphasespin_3(arraysize,mdyear), &
         phen%doyphasespin_4(arraysize,mdyear))

  end subroutine alloc_phenvariable


  subroutine dealloc_phenvariable(phen)

    implicit none

    type(phen_variable), intent(inout) :: phen

    deallocate(phen%Tkshed)
    deallocate(phen%phase)
    deallocate(phen%doyphase)
    deallocate(phen%phen)
    deallocate(phen%aphen)
    deallocate(phen%phasespin)
    deallocate(phen%doyphasespin_1)
    deallocate(phen%doyphasespin_2)
    deallocate(phen%doyphasespin_3)
    deallocate(phen%doyphasespin_4)

  end subroutine dealloc_phenvariable


  subroutine print_phenvariable(phen)

    implicit none

    type(phen_variable), intent(in) :: phen

    write(*,*) 'phen%Tkshed ', phen%Tkshed
    write(*,*) 'phen%phase ', phen%phase
    write(*,*) 'phen%doyphase ', phen%doyphase
    write(*,*) 'phen%phen ', phen%phen
    write(*,*) 'phen%aphen ', phen%aphen
    write(*,*) 'phen%phasespin ', phen%phasespin
    write(*,*) 'phen%doyphasespin_1 ', phen%doyphasespin_1
    write(*,*) 'phen%doyphasespin_2 ', phen%doyphasespin_2
    write(*,*) 'phen%doyphasespin_3 ', phen%doyphasespin_3
    write(*,*) 'phen%doyphasespin_4 ', phen%doyphasespin_4

  end subroutine print_phenvariable


  subroutine zero_phenvariable(phen)

    use cable_def_types_mod, only: r_2

    implicit none

    type(phen_variable), intent(inout) :: phen

    phen%phase          = 0
    phen%Tkshed         = 0.0_r_2
    phen%doyphase       = 0
    phen%phen           = 0.0
    phen%aphen          = 0.0
    phen%phasespin      = 0
    phen%doyphasespin_1 = 0
    phen%doyphasespin_2 = 0
    phen%doyphasespin_3 = 0
    phen%doyphasespin_4 = 0

  end subroutine zero_phenvariable


  subroutine read_netcdf_phen_var(filename, phen)

    use netcdf, only: nf90_open, nf90_nowrite, &
         nf90_inq_varid, nf90_get_var, nf90_close
#ifdef __MPI__
    use mpi, only: MPI_Abort
#endif
    use cable_def_types_mod, only: nc_err

    implicit none

    character(len=*),    intent(in)    :: filename
    type(phen_variable), intent(inout) :: phen

    logical :: existfile
    integer :: fid, vid
#ifdef __MPI__
    integer :: ierr
#endif

    ! open netCDF file
    inquire(file=trim(filename), exist=existfile)
    if (.not. existfile) then
       write(*,*) filename, ' does not exist!'
#ifdef __MPI__
       call MPI_Abort(0, 178, ierr)
#else
       stop 178
#endif
    endif

    ! open netCDF file
    call nc_err(nf90_open(trim(filename), nf90_nowrite, fid))

    ! read variables
    ! integer vectors
    call nc_err(nf90_inq_varid(fid, 'phase', vid))
    call nc_err(nf90_get_var(fid, vid, phen%phase))
    ! integer arrays
    call nc_err(nf90_inq_varid(fid, 'doyphase', vid))
    call nc_err(nf90_get_var(fid, vid, phen%doyphase))
    call nc_err(nf90_inq_varid(fid, 'phasespin', vid))
    call nc_err(nf90_get_var(fid, vid, phen%phasespin))
    call nc_err(nf90_inq_varid(fid, 'doyphasespin_1', vid))
    call nc_err(nf90_get_var(fid, vid, phen%doyphasespin_1))
    call nc_err(nf90_inq_varid(fid, 'doyphasespin_2', vid))
    call nc_err(nf90_get_var(fid, vid, phen%doyphasespin_2))
    call nc_err(nf90_inq_varid(fid, 'doyphasespin_3', vid))
    call nc_err(nf90_get_var(fid, vid, phen%doyphasespin_3))
    call nc_err(nf90_inq_varid(fid, 'doyphasespin_4', vid))
    call nc_err(nf90_get_var(fid, vid, phen%doyphasespin_4))
    ! real vectors
    call nc_err(nf90_inq_varid(fid, 'phen', vid))
    call nc_err(nf90_get_var(fid, vid, phen%phen))
    call nc_err(nf90_inq_varid(fid, 'aphen', vid))
    call nc_err(nf90_get_var(fid, vid, phen%aphen))
    ! double vectors
    call nc_err(nf90_inq_varid(fid, 'tkshed', vid))
    call nc_err(nf90_get_var(fid, vid, phen%tkshed))

    ! close NetCDF file
    call nc_err(nf90_close(fid))

  end subroutine read_netcdf_phen_var


  subroutine write_netcdf_phen_var(filename, phen)

    use netcdf, only: nf90_create, nf90_clobber, nf90_64bit_offset, &
         nf90_def_dim, nf90_def_var, nf90_int, nf90_float, nf90_double, &
         nf90_enddef, nf90_put_var, nf90_close
    use cable_def_types_mod, only: nc_err

    implicit none

    character(len=*),    intent(in) :: filename
    type(phen_variable), intent(in) :: phen

    integer :: fid
    integer :: dimid1, dimid2, dimid3, dimid4
    integer :: i
    integer, dimension(ncasa_phen) :: vid

    ! create netCDF file
    call nc_err(nf90_create(trim(filename), ior(nf90_clobber, nf90_64bit_offset), fid))

    ! define dimensions
    ! land
    call nc_err(nf90_def_dim(fid, 'dim1', size(phen%phase, 1), dimid1))
    ! phenology phases
    call nc_err(nf90_def_dim(fid, 'dim2', size(phen%doyphase, 2), dimid2))
    ! days of the year = 365
    call nc_err(nf90_def_dim(fid, 'dim3', size(phen%phasespin, 2), dimid3))
    ! vegetation types
    call nc_err(nf90_def_dim(fid, 'dim4', size(phen%TKshed, 1), dimid4))

    ! define variables
    i = 1
    ! define integer vector variables [dim1]
    call nc_err(nf90_def_var(fid, 'phase', nf90_int, &
         [dimid1], vid(i)), i)

    ! define integer array variables [dim1, dim2]
    call nc_err(nf90_def_var(fid, 'doyphase', nf90_int, &
         [dimid1, dimid2], vid(i)), i)

    ! define integer array variables [dim1, dim3]
    call nc_err(nf90_def_var(fid, 'phasespin', nf90_int, &
         [dimid1, dimid3], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'doyphasespin_1', nf90_int, &
         [dimid1, dimid3], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'doyphasespin_2', nf90_int, &
         [dimid1, dimid3], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'doyphasespin_3', nf90_int, &
         [dimid1, dimid3], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'doyphasespin_4', nf90_int, &
         [dimid1, dimid3], vid(i)), i)

    ! define real vector variables [dim1]
    call nc_err(nf90_def_var(fid, 'phen', nf90_float, &
         [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'aphen', nf90_float, &
         [dimid1], vid(i)), i)

    ! define double vector variables [dim4]
    call nc_err(nf90_def_var(fid, 'tkshed', nf90_double, &
         [dimid4], vid(i)), i)

    ! end define mode
    call nc_err(nf90_enddef(fid))

    ! put variables
    i = 1
    call nc_err(nf90_put_var(fid, vid(i), phen%phase), i)
    call nc_err(nf90_put_var(fid, vid(i), phen%doyphase), i)
    call nc_err(nf90_put_var(fid, vid(i), phen%phasespin), i)
    call nc_err(nf90_put_var(fid, vid(i), phen%doyphasespin_1), i)
    call nc_err(nf90_put_var(fid, vid(i), phen%doyphasespin_2), i)
    call nc_err(nf90_put_var(fid, vid(i), phen%doyphasespin_3), i)
    call nc_err(nf90_put_var(fid, vid(i), phen%doyphasespin_4), i)
    call nc_err(nf90_put_var(fid, vid(i), phen%phen), i)
    call nc_err(nf90_put_var(fid, vid(i), phen%aphen), i)
    call nc_err(nf90_put_var(fid, vid(i), phen%TKshed), i)

    ! close NetCDF file
    call nc_err(nf90_close(fid))

  end subroutine write_netcdf_phen_var

end module phenvariable
