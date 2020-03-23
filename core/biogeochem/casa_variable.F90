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
!   casavariable with subroutine alloc_casavariable
!   phenvariable with subroutine alloc_phenvariable

MODULE casadimension

  USE cable_def_types_mod, ONLY : mp, r_2, mvtype, ms

  IMPLICIT NONE

  INTEGER, PARAMETER :: mdyear=365         ! days per year
  INTEGER, PARAMETER :: mdmonth=30         ! days per month
  INTEGER, PARAMETER :: mdweek=7           ! days per week
  INTEGER, PARAMETER :: mmyear=12          ! month per year
  INTEGER, PARAMETER :: mt=36500           ! integration time step
  INTEGER, PARAMETER :: mpftmax=2          ! max. PFT/cell
  INTEGER, PARAMETER :: mplant   = 3       ! plant pools
  INTEGER, PARAMETER :: mlitter  = 3       ! litter pools
  INTEGER, PARAMETER :: msoil    = 3       ! soil pools
  INTEGER, PARAMETER :: mso      = 12      ! soil order number
  INTEGER, PARAMETER :: mhwp     = 1       ! harvested wood pools
  INTEGER, PARAMETER :: mclear   = 1       ! forest clearing pools
  INTEGER, PARAMETER :: mheights = 10      ! height clas
  ! BP put icycle into namelist file
  INTEGER            :: icycle
  ! INTEGER, PARAMETER :: icycle=3           ! =1 for C, =2 for C+N; =3 for C+N+P
  INTEGER, PARAMETER :: mstart=1           ! starting time step
  INTEGER, PARAMETER :: mphase=4           ! phen. phases
  REAL(r_2), PARAMETER :: deltcasa = 1.0_r_2/365.0_r_2 ! fraction 1 day of year
  REAL(r_2), PARAMETER :: deltpool = 1.0_r_2           ! pool delt(1day)

END MODULE casadimension


MODULE casaparm

  USE casadimension

  IMPLICIT NONE
  
  INTEGER, PARAMETER :: initcasa = 1  ! =0 spin; 1 restart file
  INTEGER, PARAMETER :: iceland  = 17 !=13 for casa vegtype =15 for IGBP vegtype
  INTEGER, PARAMETER :: cropland = 9  ! 12 and 14 for IGBP vegtype
  INTEGER, PARAMETER :: croplnd2 = 10 ! ditto
  INTEGER, PARAMETER :: forest   = 3
  INTEGER, PARAMETER :: shrub    = 2
  INTEGER, PARAMETER :: grass    = 1
  INTEGER, PARAMETER :: icewater = 0
  INTEGER, PARAMETER :: LEAF     = 1
  INTEGER, PARAMETER :: WOOD     = 2
  INTEGER, PARAMETER :: FROOT    = 3
  !  INTEGER, PARAMETER :: LABILE = 4
  INTEGER, PARAMETER :: METB    = 1
  INTEGER, PARAMETER :: STR     = 2
  INTEGER, PARAMETER :: CWD     = 3
  INTEGER, PARAMETER :: MIC     = 1
  INTEGER, PARAMETER :: SLOW    = 2
  INTEGER, PARAMETER :: PASS    = 3
  INTEGER, PARAMETER :: PLAB    = 1
  INTEGER, PARAMETER :: PSORB   = 2
  INTEGER, PARAMETER :: POCC    = 3
  !! vh_js !! LALLOC moved to bgcdriver to allow for value to be switchable
  ! INTEGER, PARAMETER :: LALLOC  = 0      !=0 constant; 1 variable
  REAL(r_2), PARAMETER :: z30=0.3_r_2
  REAL(r_2), PARAMETER :: R0=0.3_r_2
  REAL(r_2), PARAMETER :: S0=0.3_r_2
  REAL(r_2), PARAMETER :: fixed_stem=1.0_r_2/3.0_r_2
  REAL(r_2), PARAMETER :: Q10alloc=2.0_r_2
  REAL(r_2), PARAMETER :: ratioNCstrfix = 1.0_r_2/150.0_r_2
  REAL(r_2), PARAMETER :: ratioNPstrfix = 25.0_r_2
  REAL(r_2), PARAMETER :: fracCbiomass = 0.50_r_2
  REAL(r_2), PARAMETER :: tsoilrefc=25.0_r_2
  REAL(r_2), PARAMETER :: tkzeroc=273.15_r_2
  REAL(r_2), PARAMETER :: frootparma = 0.3192_r_2
  REAL(r_2), PARAMETER :: frootparmb =-0.0485_r_2
  REAL(r_2), PARAMETER :: frootparmc = 0.1755_r_2
  REAL(r_2), PARAMETER :: xweightalloc = 0.2_r_2
  !  REAL(r_2), PARAMETER :: xkplab=0.5_r_2*deltcasa
  !  REAL(r_2), PARAMETER :: xkpsorb=0.01_r_2*deltcasa
  !  REAL(r_2), PARAMETER :: xkpocc =0.01_r_2*deltcasa
  
END MODULE casaparm


MODULE casavariable

  USE casadimension

  IMPLICIT NONE

  SAVE

  character(len=200) :: casa_timeunits

  TYPE casa_biome
     INTEGER,   DIMENSION(:),POINTER :: ivt2 => null()
     REAL(r_2), DIMENSION(:),POINTER :: xkleafcoldmax => null(),  &
          xkleafcoldexp => null(),  &
          xkleafdrymax => null(),   &
          xkleafdryexp => null(),   &
          glaimax => null(),        &
          glaimin => null(),        &
          sla => null(),            &
          ratiofrootleaf => null(), &
          kroot => null(),          &
          krootlen => null(),       &
          rootdepth => null(),      &
          kuptake => null(),        &
          kminN => null(),          &
          kuplabP => null(),        &
          kclabrate => null(),      &
          xnpmax => null(),         &
          q10soil => null(),        &
          xkoptlitter => null(),    &
          xkoptsoil => null(),      &
          xkplab => null(),         &
          xkpsorb => null(),        &
          xkpocc => null(),         &
          prodptase => null(),      &
          costnpup => null(),       &
          maxfinelitter => null(),  &
          maxcwd => null(),         &
          nintercept => null(),     &
          nslope => null(),         &
          la_to_sa => null(),       &
          vcmax_scalar => null(),   &
          disturbance_interval => null(), &
          DAMM_EnzPool => null(), &
          DAMM_KMO2 => null(), &
          DAMM_KMcp => null(), &
          DAMM_Ea => null(), &
          DAMM_alpha => null()

     REAL(r_2), DIMENSION(:,:),POINTER :: plantrate => null(),     &
          rmplant => null(),         &
          fracnpptoP => null(),      &
          fraclignin => null(),      &
          fraclabile => null(),      &
          ratioNCplantmin => null(), &
          ratioNCplantmax => null(), &
          ratioNPplantmin => null(), &
          ratioNPplantmax => null(), &
          fracLigninplant => null(), &
          ftransNPtoL => null(),     &
          ftransPPtoL => null(),     &
          litterrate => null(),      &
          ratioPcplantmin => null(), &
          ratioPcplantmax => null()
     REAL(r_2), DIMENSION(:,:),POINTER :: soilrate => null()
  END TYPE casa_biome

  
  TYPE casa_pool
     REAL(r_2), DIMENSION(:),POINTER :: Clabile => null(),       &
          dClabiledt => null(),    &
          Ctot => null(),         &          !! vh_js !!
          Ctot_0 => null()
     REAL(r_2), DIMENSION(:,:),POINTER :: Cplant => null(),      &
          Nplant => null(),        &
          Pplant => null(),        &
          dCplantdt => null(),     &
          dNplantdt => null(),     &
          dPplantdt => null(),     &
          ratioNCplant => null(),  &
          ratioNPplant => null()
     REAL(r_2), DIMENSION(:),POINTER :: Nsoilmin => null(),      &
          Psoillab => null(),      &
          Psoilsorb => null(),     &
          Psoilocc => null(),      &
          dNsoilmindt => null(),   &
          dPsoillabdt => null(),   &
          dPsoilsorbdt => null(),  &
          dPsoiloccdt => null()
     REAL(r_2), DIMENSION(:,:), POINTER :: Clitter => null(),    &
          Nlitter => null(),       &
          Plitter => null(),       &
          dClitterdt => null(),    &
          dNlitterdt => null(),    &
          dPlitterdt => null(),    &
          ratioNClitter => null(), &
          ratioNPlitter => null()
     REAL(r_2), DIMENSION(:,:),POINTER :: Csoil => null(),       &
          Nsoil => null(),         &
          Psoil => null(),         &
          dCsoildt => null(),      &
          dNsoildt => null(),      &
          dPsoildt => null(),      &
          ratioNCsoil => null(),   &
          ratioNCsoilnew => null(),&
          ratioNPsoil => null(),   &
          ratioNCsoilmin => null(),&
          ratioNCsoilmax => null(),&
          ratioPCsoil => null(),   &
          ratioPCplant => null(),  &
          ratioPClitter => null()
  END TYPE casa_pool

  
  TYPE casa_flux
     REAL(r_2), DIMENSION(:),POINTER :: Cgpp => null(),          &
          Cnpp => null(),          &
          Crp => null(),           &
          Crgplant => null(),      &
          Nminfix => null(),       &
          Nminuptake => null(),    &
          Plabuptake => null(),    &
          Clabloss => null(),      &
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
     REAL(r_2), DIMENSION(:,:),POINTER :: fracCalloc => null(),  &
          fracNalloc => null(),    &
          fracPalloc => null(),    &
          Crmplant => null(),      &
          kplant => null(),        &
                                !! vh_js !! additional diagnostic
          Cplant_turnover => null()
     REAL(r_2), DIMENSION(:,:,:),POINTER :: fromPtoL => null()
     REAL(r_2), DIMENSION(:),POINTER :: Cnep => null(),        &
          Crsoil => null(),      &
          Nmindep => null(),     &
          Nminloss => null(),    &
          Nminleach => null(),   &
          Nupland => null(),     &
          Nlittermin => null(),  &
          Nsmin => null(),       &
          Nsimm => null(),       &
          Nsnet => null(),       &
          fNminloss => null(),   &
          fNminleach => null(),  &
          Pdep => null(),        &
          Pwea => null(),        &
          Pleach => null(),      &
          Ploss => null(),       &
          Pupland => null(),     &
          Plittermin => null(),  &
          Psmin => null(),       &
          Psimm => null(),       &
          Psnet => null(),       &
          fPleach => null(),     &
          kplab => null(),       &
          kpsorb => null(),      &
          kpocc => null(),       &
          kmlabp => null(),      &
          Psorbmax => null(),    &
          !! additional diagnostics for partitioning biomass turnover
          Cplant_turnover_disturbance => null(), &
          Cplant_turnover_crowding  => null(), &
          Cplant_turnover_resource_limitation => null()

     REAL(r_2), DIMENSION(:,:),POINTER    :: klitter => null()
     REAL(r_2), DIMENSION(:,:),POINTER    :: ksoil => null()
     REAL(r_2), DIMENSION(:,:,:),POINTER  :: fromLtoS => null()
     REAL(r_2), DIMENSION(:,:,:),POINTER  :: fromStoS => null()
     REAL(r_2), DIMENSION(:,:),POINTER    :: fromLtoCO2 => null()
     REAL(r_2), DIMENSION(:,:),POINTER    :: fromStoCO2 => null()
     REAL(r_2), DIMENSION(:,:),POINTER    :: FluxCtolitter => null()
     REAL(r_2), DIMENSION(:,:),POINTER    :: FluxNtolitter => null()
     REAL(r_2), DIMENSION(:,:),POINTER    :: FluxPtolitter => null()
     REAL(r_2), DIMENSION(:,:),POINTER    :: FluxCtosoil => null()
     REAL(r_2), DIMENSION(:,:),POINTER    :: FluxNtosoil => null()
     REAL(r_2), DIMENSION(:,:),POINTER    :: FluxPtosoil => null()
     REAL(r_2), DIMENSION(:),POINTER      :: FluxCtoCO2 => null()
     REAL(r_2), DIMENSION(:),POINTER    :: FluxCtohwp => null()
     REAL(r_2), DIMENSION(:),POINTER    :: FluxNtohwp => null()
     REAL(r_2), DIMENSION(:),POINTER    :: FluxPtohwp => null()
     REAL(r_2), DIMENSION(:),POINTER    :: FluxCtoclear => null()
     REAL(r_2), DIMENSION(:),POINTER    :: FluxNtoclear => null()
     REAL(r_2), DIMENSION(:),POINTER    :: FluxPtoclear => null()
     REAL(r_2), DIMENSION(:),POINTER    :: CtransferLUC => null()

     !CVH variables inherited from BLAZE
     REAL(r_2), DIMENSION(:,:,:),POINTER  :: fromPtoL_fire => null()
     REAL(r_2), DIMENSION(:,:),POINTER    :: klitter_fire => null()
     REAL(r_2), DIMENSION(:,:),POINTER    :: klitter_tot => null()  ! sum of fire turnover and non-fire turnover (litter)
     REAL(r_2), DIMENSION(:,:),POINTER    :: kplant_fire => null()
     REAL(r_2), DIMENSION(:,:),POINTER    :: kplant_tot => null()  ! sum of fire turnover and non-fire turnover (plants)

     !CVH diagnostic: CO2 emissions from fire
     REAL(r_2), DIMENSION(:),POINTER      :: fluxCtoCO2_plant_fire => null()
     REAL(r_2), DIMENSION(:),POINTER      :: fluxCtoCO2_litter_fire => null()
     ! contribution to fire emissions from individual plant pools
     REAL(r_2), DIMENSION(:,:),POINTER      :: fluxfromPtoCO2_fire => null()
     ! contribution to fire emissions from individual litter pools
     REAL(r_2), DIMENSION(:,:),POINTER      :: fluxfromLtoCO2_fire => null()
     REAL(r_2), DIMENSION(:),POINTER      :: fluxNtoAtm_fire => null()
     !REAL(r_2), DIMENSION(:,:,:),POINTER  :: fire_mortality_vs_height => null()

     ! Diagnostic fluxes for use in 13C
     REAL(r_2), DIMENSION(:,:,:), POINTER :: FluxFromPtoL => null()
     REAL(r_2), DIMENSION(:,:,:), POINTER :: FluxFromLtoS => null()
     REAL(r_2), DIMENSION(:,:,:), POINTER :: FluxFromStoS => null()
     REAL(r_2), DIMENSION(:,:),   POINTER :: FluxFromPtoCO2 => null()
     REAL(r_2), DIMENSION(:,:),   POINTER :: FluxFromLtoCO2 => null()
     REAL(r_2), DIMENSION(:,:),   POINTER :: FluxFromStoCO2 => null()
     REAL(r_2), DIMENSION(:),     POINTER :: FluxFromPtoHarvest => null()
  END TYPE casa_flux

  
  TYPE casa_met
     REAL(r_2), DIMENSION(:),POINTER    :: glai => null(),     &
          Tairk => null(),    &
          precip => null(),   &
          tsoilavg => null(), &
          moistavg => null(), &
          btran => null()
     INTEGER, DIMENSION(:), POINTER     :: lnonwood => null()
     REAL(r_2), DIMENSION(:,:), POINTER :: Tsoil => null(),    &
          moist => null()
     INTEGER, DIMENSION(:), POINTER     :: iveg2 => null(),    &
          ijgcm => null(),    &
          isorder => null()
     REAL(r_2), DIMENSION(:), POINTER   :: lat => null(),      &
          lon => null(),      &
          areacell => null()
     ! added yp wang 5/nov/2012
     REAL(r_2), DIMENSION(:,:), POINTER :: Tairkspin => null(),&
          cgppspin => null(),&
          crmplantspin_1 => null(),&
          crmplantspin_2 => null(),&
          crmplantspin_3 => null(),&
          Tsoilspin_1 => null(),&
          Tsoilspin_2 => null(),&
          Tsoilspin_3 => null(),&
          Tsoilspin_4 => null(),&
          Tsoilspin_5 => null(),&
          Tsoilspin_6 => null(),&
          moistspin_1 => null(),&
          moistspin_2 => null(),&
          moistspin_3 => null(),&
          moistspin_4 => null(),&
          moistspin_5 => null(),&
          moistspin_6 => null(), &
          mtempspin => null()
     ! 13C
     real(r_2), dimension(:,:), pointer :: cAn12spin => null() ! daily cumulated total 12CO2 net assimilation in [g(C)/m2]
     real(r_2), dimension(:,:), pointer :: cAn13spin => null() ! daily cumulated total 13CO2 net assimilation in [g(13C)/m2]
  END TYPE casa_met

  
  TYPE casa_balance
     REAL(r_2), DIMENSION(:),POINTER   :: FCgppyear => null(), FCnppyear => null(), &
          FCrmleafyear => null(), FCrmwoodyear => null(), FCrmrootyear => null(), FCrgrowyear => null(), &
          FCrpyear => null(), FCrsyear => null(),FCneeyear => null(),  dCdtyear => null(), &
          LAImax => null(), Cleafmean => null(), Crootmean => null(), &
          FNdepyear => null(), FNfixyear => null(), FNsnetyear => null(), FNupyear => null(), &
          FNleachyear => null(), FNlossyear => null(), &
          FPweayear => null(), FPdustyear => null(), FPsnetyear => null(), &
          FPupyear => null(), FPleachyear => null(), FPlossyear => null()

     REAL(r_2), DIMENSION(:,:),POINTER :: glaimon => null(), glaimonx => null()
     REAL(r_2), DIMENSION(:,:),POINTER :: cplantlast => null(), nplantlast => null(), pplantlast => null()
     REAL(r_2), DIMENSION(:,:),POINTER :: clitterlast => null(), nlitterlast => null(), plitterlast => null()
     REAL(r_2), DIMENSION(:,:),POINTER :: csoillast => null(), nsoillast => null(), psoillast => null()
     REAL(r_2), DIMENSION(:),  POINTER :: nsoilminlast => null(), psoillablast => null(),  &
          psoilsorblast => null(), psoilocclast => null(), &
          cbalance => null(), nbalance => null(), pbalance => null(), &
          sumcbal => null(), sumnbal => null(), sumpbal => null()
     REAL(r_2), DIMENSION(:),POINTER   :: clabilelast => null()
  END TYPE casa_balance

  
  ! The following declarations are removed and have to be passed using
  ! parameter list for each subroutine (BP apr2010)
  !  TYPE (casa_biome)              :: casabiome
  !  TYPE (casa_pool)               :: casapool
  !  TYPE (casa_flux)               :: casaflux
  !  TYPE (casa_met)                :: casamet
  !  TYPE (casa_balance)            :: casabal

  
  ! Added filename type for casaCNP (BP apr2010)
  TYPE casafiles_type
     CHARACTER(LEN=200) :: cnpbiome    ! file for biome-specific BGC parameters
     CHARACTER(LEN=200) :: cnppoint    ! file for point-specific BGC inputs
     CHARACTER(LEN=200) :: cnpepool    ! file for end-of-run pool sizes
     CHARACTER(LEN=200) :: cnpipool=''    ! file for inital pool sizes
     CHARACTER(LEN=200) :: cnpmetin      ! met file for spin up
     CHARACTER(LEN=200) :: cnpmetout     ! met file for spin up
     CHARACTER(LEN=200) :: ndep          ! N deposition input file
     ! added yp wang
     CHARACTER(LEN=200) :: cnpspin       ! input file for spin up
     CHARACTER(LEN=200) :: dump_cnpspin  ! name of dump file for spinning casa-cnp

     CHARACTER(LEN=200) :: phen        ! leaf phenology datafile
     CHARACTER(LEN=200) :: cnpflux     ! modelled mean yearly CNP fluxes
     LOGICAL            :: l_ndep
     ! added vh
     CHARACTER(LEN=200) :: c2cdumppath='' ! cable2casa dump for casa spinup
     CHARACTER(LEN=200) :: out=''    ! casa output file
  END TYPE casafiles_type
  TYPE(casafiles_type) :: casafile

  
Contains

  
  SUBROUTINE alloc_casavariable(casabiome, casapool, casaflux, casamet, casabal, arraysize)
    
    IMPLICIT NONE
    
    TYPE (casa_biome)  , INTENT(INOUT) :: casabiome
    TYPE (casa_pool)   , INTENT(INOUT) :: casapool
    TYPE (casa_flux)   , INTENT(INOUT) :: casaflux
    TYPE (casa_met)    , INTENT(INOUT) :: casamet
    TYPE (casa_balance), INTENT(INOUT) :: casabal
    INTEGER,             INTENT(IN) :: arraysize

    ALLOCATE(casabiome%ivt2(mvtype),               &
         casabiome%xkleafcoldmax(mvtype),          &
         casabiome%xkleafcoldexp(mvtype),          &
         casabiome%xkleafdrymax(mvtype),           &
         casabiome%xkleafdryexp(mvtype),           &
         casabiome%glaimax(mvtype),                &
         casabiome%glaimin(mvtype),                &
         casabiome%sla(mvtype),                    &
         casabiome%ratiofrootleaf(mvtype),         &
         casabiome%kroot(mvtype),                  &
         casabiome%krootlen(mvtype),               &
         casabiome%rootdepth(mvtype),              &
         casabiome%kuptake(mvtype),                &
         casabiome%kminN(mvtype),                  &
         casabiome%KuplabP(mvtype),                &
         casabiome%kclabrate(mvtype),              &
         casabiome%xnpmax(mvtype),                 &
         casabiome%q10soil(mvtype),                &
         casabiome%xkoptlitter(mvtype),            &
         casabiome%xkoptsoil(mvtype),              &
         casabiome%xkplab(mso),                    &
         casabiome%xkpsorb(mso),                   &
         casabiome%xkpocc(mso),                    &
         casabiome%prodptase(mvtype),              &
         casabiome%costnpup(mvtype),               &
         casabiome%maxfinelitter(mvtype),          &
         casabiome%maxcwd(mvtype),                 &
         casabiome%nintercept(mvtype),             &
         casabiome%nslope(mvtype),                 &
         casabiome%plantrate(mvtype,mplant),       &
         casabiome%rmplant(mvtype,mplant),         &
         casabiome%fracnpptoP(mvtype,mplant),      &
         casabiome%fraclignin(mvtype,mplant),      &
         casabiome%fraclabile(mvtype,mplant),      &
         casabiome%ratioNCplantmin(mvtype,mplant), &
         casabiome%ratioNCplantmax(mvtype,mplant), &
         casabiome%ratioNPplantmin(mvtype,mplant), &
         casabiome%ratioNPplantmax(mvtype,mplant), &
         casabiome%fracLigninplant(mvtype,mplant), &
         casabiome%ftransNPtoL(mvtype,mplant),     &
         casabiome%ftransPPtoL(mvtype,mplant),     &
         casabiome%litterrate(mvtype,mlitter),     &
         casabiome%soilrate(mvtype,msoil),         &
         !  casabiome%ratioPcplantmax(mvtype,leaf),   &
         !  casabiome%ratioPcplantmin(mvtype,leaf)    &
         !! vh_js !!
         casabiome%ratioPcplantmax(mvtype,mplant),   &
         casabiome%ratioPcplantmin(mvtype,mplant),   &
         casabiome%la_to_sa(mvtype),                 &
         casabiome%vcmax_scalar(mvtype),             &
         casabiome%disturbance_interval(mvtype),     &
         casabiome%DAMM_EnzPool(mvtype),     &
         casabiome%DAMM_KMO2(mvtype),     &
         casabiome%DAMM_KMcp(mvtype),     &
         casabiome%DAMM_Ea(mvtype),       &
         casabiome%DAMM_alpha(mvtype)     &
         )

    ALLOCATE(casapool%Clabile(arraysize),           &
         casapool%dClabiledt(arraysize),            &
         casapool%Cplant(arraysize,mplant),         &
         casapool%Nplant(arraysize,mplant),         &
         casapool%Pplant(arraysize,mplant),         &
         casapool%dCplantdt(arraysize,mplant),      &
         casapool%dNplantdt(arraysize,mplant),      &
         casapool%dPplantdt(arraysize,mplant),      &
         casapool%ratioNCplant(arraysize,mplant),   &
         casapool%ratioNPplant(arraysize,mplant),   &
         casapool%Nsoilmin(arraysize),              &
         casapool%Psoillab(arraysize),              &
         casapool%Psoilsorb(arraysize),             &
         casapool%Psoilocc(arraysize),              &
         casapool%dNsoilmindt(arraysize),           &
         casapool%dPsoillabdt(arraysize),           &
         casapool%dPsoilsorbdt(arraysize),          &
         casapool%dPsoiloccdt(arraysize),           &
         casapool%Clitter(arraysize,mlitter),       &
         casapool%Nlitter(arraysize,mlitter),       &
         casapool%Plitter(arraysize,mlitter),       &
         casapool%dClitterdt(arraysize,mlitter),    &
         casapool%dNlitterdt(arraysize,mlitter),    &
         casapool%dPlitterdt(arraysize,mlitter),    &
         casapool%ratioNClitter(arraysize,mlitter), &
         casapool%ratioNPlitter(arraysize,mlitter), &
         casapool%Csoil(arraysize,msoil),           &
         casapool%Nsoil(arraysize,msoil),           &
         casapool%Psoil(arraysize,msoil),           &
         casapool%dCsoildt(arraysize,msoil),        &
         casapool%dNsoildt(arraysize,msoil),        &
         casapool%dPsoildt(arraysize,msoil),        &
         casapool%ratioNCsoil(arraysize,msoil),     &
         casapool%ratioNPsoil(arraysize,msoil),     &
         casapool%ratioNCsoilnew(arraysize,msoil),  &
         casapool%ratioNCsoilmin(arraysize,msoil),  &
         casapool%ratioNCsoilmax(arraysize,msoil),  &
         casapool%ratioPCsoil(arraysize,msoil),     &
         casapool%ratioPCplant(arraysize,mplant),   &
         casapool%ratioPClitter(arraysize,mlitter), &
         casapool%Ctot_0(arraysize),                &
         casapool%Ctot(arraysize)   )               

    ALLOCATE(casaflux%Cgpp(arraysize),                 &
         casaflux%Cnpp(arraysize),                     &
         casaflux%Crp(arraysize),                      &
         casaflux%Crgplant(arraysize),                 &
         casaflux%Nminfix(arraysize),                  &
         casaflux%Nminuptake(arraysize),               &
         casaflux%Plabuptake(arraysize),               &
         casaflux%Clabloss(arraysize),                 &
         casaflux%fracClabile(arraysize),              &
         casaflux%fracCalloc(arraysize,mplant),        &
         casaflux%fracNalloc(arraysize,mplant),        &
         casaflux%fracPalloc(arraysize,mplant),        &
         casaflux%kplant(arraysize,mplant),            &
         casaflux%Crmplant(arraysize,mplant),          &
         casaflux%fromPtoL(arraysize,mlitter,mplant),  &
         casaflux%Cnep(arraysize),                     &
         casaflux%Crsoil(arraysize),                   &
         casaflux%Nmindep(arraysize),                  &
         casaflux%Nminloss(arraysize),                 &
         casaflux%Nminleach(arraysize),                &
         casaflux%Nupland(arraysize),                  &
         casaflux%Nlittermin(arraysize),               &
         casaflux%Nsmin(arraysize),                    &
         casaflux%Nsimm(arraysize),                    &
         casaflux%Nsnet(arraysize),                    &
         casaflux%fNminloss(arraysize),                &
         casaflux%fNminleach(arraysize),               &
         casaflux%Pdep(arraysize),                     &
         casaflux%Pwea(arraysize),                     &
         casaflux%Pleach(arraysize),                   &
         casaflux%Ploss(arraysize),                    &
         casaflux%Pupland(arraysize),                  &
         casaflux%Plittermin(arraysize),               &
         casaflux%Psmin(arraysize),                    &
         casaflux%Psimm(arraysize),                    &
         casaflux%Psnet(arraysize),                    &
         casaflux%fPleach(arraysize),                  &
         casaflux%kplab(arraysize),                    &
         casaflux%kpsorb(arraysize),                   &
         casaflux%kpocc(arraysize),                    &
         casaflux%kmlabP(arraysize),                   &
         casaflux%Psorbmax(arraysize),                 &
         casaflux%klitter(arraysize,mlitter),          &
         casaflux%ksoil(arraysize,msoil),              &
         casaflux%fromLtoS(arraysize,msoil,mlitter),   &
         casaflux%fromStoS(arraysize,msoil,msoil),     &
         casaflux%fromLtoCO2(arraysize,mlitter),       &
         casaflux%fromStoCO2(arraysize,msoil),         &
         casaflux%stemnpp(arraysize),                  &
         casaflux%frac_sapwood(arraysize),             &
         casaflux%sapwood_area(arraysize), &
         casaflux%fharvest(arraysize), &
         casaflux%Charvest(arraysize), &
         casaflux%Nharvest(arraysize), &
         casaflux%Pharvest(arraysize), &
         casaflux%fcrop(arraysize), &
         casaflux%Cplant_turnover(arraysize,mplant) , &
         casaflux%Cplant_turnover_disturbance(arraysize) , &
         casaflux%Cplant_turnover_crowding(arraysize) , &
         casaflux%Cplant_turnover_resource_limitation(arraysize))

    !CVH Alllocate fire turnover rates and plant-to-litter partitioning coefficients
    ALLOCATE(casaflux%fromPtoL_fire(arraysize,mlitter,mplant), &
         casaflux%kplant_fire(arraysize,mplant),        &
         casaflux%klitter_fire(arraysize,mlitter),           &
         casaflux%kplant_tot(arraysize,mplant),        &
         casaflux%klitter_tot(arraysize,mlitter),           &
         casaflux%FluxCtoCO2_plant_fire(arraysize),             &
         casaflux%FluxCtoCO2_litter_fire(arraysize),             &
         casaflux%fluxfromPtoCO2_fire(arraysize,mplant),         &
         casaflux%fluxfromLtoCO2_fire(arraysize,mlitter),         &
         casaflux%FluxNtoAtm_fire(arraysize))
    !casaflux%fire_mortality_vs_height(arraysize,mheights,2))

    ALLOCATE(casaflux%FluxCtolitter(arraysize,mlitter),    &
         casaflux%FluxNtolitter(arraysize,mlitter),    &
         casaflux%FluxPtolitter(arraysize,mlitter))

    ALLOCATE(casaflux%FluxCtosoil(arraysize,msoil),        &
         casaflux%FluxNtosoil(arraysize,msoil),        &
         casaflux%FluxPtosoil(arraysize,msoil))

    ALLOCATE(casaflux%FluxCtohwp(arraysize),    &
         casaflux%FluxNtohwp(arraysize),    &
         casaflux%FluxPtohwp(arraysize))

    ALLOCATE(casaflux%FluxCtoclear(arraysize),    &
         casaflux%FluxNtoclear(arraysize),    &
         casaflux%FluxPtoclear(arraysize))

    ALLOCATE(casaflux%CtransferLUC(arraysize))

    ALLOCATE(casaflux%FluxCtoco2(arraysize))

    ALLOCATE(casaflux%FluxFromPtoL(arraysize,mplant,mlitter), &
         casaflux%FluxFromLtoS(arraysize,mlitter,msoil), &
         casaflux%FluxFromStoS(arraysize,msoil,msoil), &
         casaflux%FluxFromPtoCO2(arraysize,mplant), &
         casaflux%FluxFromLtoCO2(arraysize,mlitter), &
         casaflux%FluxFromStoCO2(arraysize,msoil), &
         casaflux%FluxFromPtoHarvest(arraysize))

    ALLOCATE(casamet%glai(arraysize),                 &
         casamet%lnonwood(arraysize),             &
         casamet%Tairk(arraysize),                &
         casamet%precip(arraysize),               &
         casamet%tsoilavg(arraysize),             &
         casamet%moistavg(arraysize),             &
         casamet%btran(arraysize),                &
         casamet%Tsoil(arraysize,ms),             &
         casamet%moist(arraysize,ms),             &
         casamet%iveg2(arraysize),                &
         casamet%ijgcm(arraysize),                &
         casamet%isorder(arraysize),              &
         casamet%lat(arraysize),                  &
         casamet%lon(arraysize),                  &
         casamet%areacell(arraysize),             &
         casamet%Tairkspin(arraysize,mdyear),     &
         casamet%cgppspin(arraysize,mdyear),      &
         casamet%crmplantspin_1(arraysize,mdyear),&
         casamet%crmplantspin_2(arraysize,mdyear),&
         casamet%crmplantspin_3(arraysize,mdyear),&
         casamet%Tsoilspin_1(arraysize,mdyear),   &
         casamet%Tsoilspin_2(arraysize,mdyear),   &
         casamet%Tsoilspin_3(arraysize,mdyear),   &
         casamet%Tsoilspin_4(arraysize,mdyear),   &
         casamet%Tsoilspin_5(arraysize,mdyear),   &
         casamet%Tsoilspin_6(arraysize,mdyear),   &
         casamet%moistspin_1(arraysize,mdyear),   &
         casamet%moistspin_2(arraysize,mdyear),   &
         casamet%moistspin_3(arraysize,mdyear),   &
         casamet%moistspin_4(arraysize,mdyear),   &
         casamet%moistspin_5(arraysize,mdyear),   &
         casamet%moistspin_6(arraysize,mdyear),  &
         casamet%mtempspin(arraysize,mdyear))
    allocate(casamet%cAn12spin(arraysize,mdyear), &
         casamet%cAn13spin(arraysize,mdyear))

    ALLOCATE(casabal%FCgppyear(arraysize),           &
         casabal%FCnppyear(arraysize),           &
         casabal%FCrpyear(arraysize),            &
         casabal%FCrmleafyear(arraysize),        &
         casabal%FCrmwoodyear(arraysize),        &
         casabal%FCrmrootyear(arraysize),        &
         casabal%FCrgrowyear(arraysize),         &
         casabal%FCrsyear(arraysize),            &
         casabal%FCneeyear(arraysize),           &
         casabal%FNdepyear(arraysize),           &
         casabal%FNfixyear(arraysize),           &
         casabal%FNsnetyear(arraysize),          &
         casabal%FNupyear(arraysize),            &
         casabal%FNleachyear(arraysize),         &
         casabal%FNlossyear(arraysize),          &
         casabal%FPweayear(arraysize),           &
         casabal%FPdustyear(arraysize),          &
         casabal%FPsnetyear(arraysize),          &
         casabal%FPupyear(arraysize),            &
         casabal%FPleachyear(arraysize),         &
         casabal%FPlossyear(arraysize),          &
         casabal%dCdtyear(arraysize),            & 
         casabal%LAImax(arraysize),              &  
         casabal%Cleafmean(arraysize),           & 
         casabal%Crootmean(arraysize)) !,           &
    !casabal%fire_mortality_vs_height(arraysize,mheights,2) )


    ALLOCATE(casabal%glaimon(arraysize,12),          &
         casabal%glaimonx(arraysize,12))

    ALLOCATE(casabal%cplantlast(arraysize,mplant),   &
         casabal%nplantlast(arraysize,mplant),   &
         casabal%pplantlast(arraysize,mplant))

    ALLOCATE(casabal%clitterlast(arraysize,mlitter), &
         casabal%nlitterlast(arraysize,mlitter), &
         casabal%plitterlast(arraysize,mlitter))

    ALLOCATE(casabal%csoillast(arraysize,msoil),     &
         casabal%nsoillast(arraysize,msoil),     &
         casabal%psoillast(arraysize,msoil))

    ALLOCATE(casabal%nsoilminlast(arraysize),        &
         casabal%psoillablast(arraysize),        &
         casabal%psoilsorblast(arraysize),       &
         casabal%psoilocclast(arraysize),        &
         casabal%cbalance(arraysize),            &
         casabal%nbalance(arraysize),            &
         casabal%pbalance(arraysize),            &
         casabal%sumcbal(arraysize),             &
         casabal%sumnbal(arraysize),             &
         casabal%sumpbal(arraysize),             &
         casabal%clabilelast(arraysize))       

  END SUBROUTINE alloc_casavariable


  SUBROUTINE alloc_sum_casavariable(sum_casapool, sum_casaflux, arraysize)

    IMPLICIT NONE

    TYPE(casa_pool), INTENT(INOUT) :: sum_casapool
    TYPE(casa_flux), INTENT(INOUT) :: sum_casaflux
    INTEGER,         INTENT(IN)    :: arraysize

    ALLOCATE(sum_casapool%Clabile(arraysize),           &
         sum_casapool%dClabiledt(arraysize),            &
         sum_casapool%Cplant(arraysize,mplant),         &
         sum_casapool%Nplant(arraysize,mplant),         &
         sum_casapool%Pplant(arraysize,mplant),         &
         sum_casapool%dCplantdt(arraysize,mplant),      &
         sum_casapool%dNplantdt(arraysize,mplant),      &
         sum_casapool%dPplantdt(arraysize,mplant),      &
         sum_casapool%ratioNCplant(arraysize,mplant),   &
         sum_casapool%ratioNPplant(arraysize,mplant),   &
         sum_casapool%Nsoilmin(arraysize),              &
         sum_casapool%Psoillab(arraysize),              &
         sum_casapool%Psoilsorb(arraysize),             &
         sum_casapool%Psoilocc(arraysize),              &
         sum_casapool%dNsoilmindt(arraysize),           &
         sum_casapool%dPsoillabdt(arraysize),           &
         sum_casapool%dPsoilsorbdt(arraysize),          &
         sum_casapool%dPsoiloccdt(arraysize),           &
         sum_casapool%Clitter(arraysize,mlitter),       &
         sum_casapool%Nlitter(arraysize,mlitter),       &
         sum_casapool%Plitter(arraysize,mlitter),       &
         sum_casapool%dClitterdt(arraysize,mlitter),    &
         sum_casapool%dNlitterdt(arraysize,mlitter),    &
         sum_casapool%dPlitterdt(arraysize,mlitter),    &
         sum_casapool%ratioNClitter(arraysize,mlitter), &
         sum_casapool%ratioNPlitter(arraysize,mlitter), &
         sum_casapool%Csoil(arraysize,msoil),           &
         sum_casapool%Nsoil(arraysize,msoil),           &
         sum_casapool%Psoil(arraysize,msoil),           &
         sum_casapool%dCsoildt(arraysize,msoil),        &
         sum_casapool%dNsoildt(arraysize,msoil),        &
         sum_casapool%dPsoildt(arraysize,msoil),        &
         sum_casapool%ratioNCsoil(arraysize,msoil),     &
         sum_casapool%ratioNPsoil(arraysize,msoil),     &
         sum_casapool%ratioNCsoilnew(arraysize,msoil),  &
         sum_casapool%ratioNCsoilmin(arraysize,msoil),  &
         sum_casapool%ratioNCsoilmax(arraysize,msoil),  &
         sum_casapool%ratioPCsoil(arraysize,msoil),     &
         sum_casapool%ratioPCplant(arraysize,mplant),   &
         sum_casapool%ratioPClitter(arraysize,mlitter)  &
         )

    ALLOCATE(sum_casaflux%Cgpp(arraysize),                 &
         sum_casaflux%Cnpp(arraysize),                     &
         sum_casaflux%Crp(arraysize),                      &
         sum_casaflux%Crgplant(arraysize),                 &
         sum_casaflux%Nminfix(arraysize),                  &
         sum_casaflux%Nminuptake(arraysize),               &
         sum_casaflux%Plabuptake(arraysize),               &
         sum_casaflux%Clabloss(arraysize),                 &
         sum_casaflux%fracClabile(arraysize),              &
         sum_casaflux%fracCalloc(arraysize,mplant),        &
         sum_casaflux%fracNalloc(arraysize,mplant),        &
         sum_casaflux%fracPalloc(arraysize,mplant),        &
         sum_casaflux%kplant(arraysize,mplant),            &
         sum_casaflux%Crmplant(arraysize,mplant),          &
         sum_casaflux%fromPtoL(arraysize,mlitter,mplant),  &
         sum_casaflux%Cnep(arraysize),                     &
         sum_casaflux%Crsoil(arraysize),                   &
         sum_casaflux%Nmindep(arraysize),                  &
         sum_casaflux%Nminloss(arraysize),                 &
         sum_casaflux%Nminleach(arraysize),                &
         sum_casaflux%Nupland(arraysize),                  &
         sum_casaflux%Nlittermin(arraysize),               &
         sum_casaflux%Nsmin(arraysize),                    &
         sum_casaflux%Nsimm(arraysize),                    &
         sum_casaflux%Nsnet(arraysize),                    &
         sum_casaflux%fNminloss(arraysize),                &
         sum_casaflux%fNminleach(arraysize),               &
         sum_casaflux%Pdep(arraysize),                     &
         sum_casaflux%Pwea(arraysize),                     &
         sum_casaflux%Pleach(arraysize),                   &
         sum_casaflux%Ploss(arraysize),                    &
         sum_casaflux%Pupland(arraysize),                  &
         sum_casaflux%Plittermin(arraysize),               &
         sum_casaflux%Psmin(arraysize),                    &
         sum_casaflux%Psimm(arraysize),                    &
         sum_casaflux%Psnet(arraysize),                    &
         sum_casaflux%fPleach(arraysize),                  &
         sum_casaflux%kplab(arraysize),                    &
         sum_casaflux%kpsorb(arraysize),                   &
         sum_casaflux%kpocc(arraysize),                    &
         sum_casaflux%kmlabP(arraysize),                   &
         sum_casaflux%Psorbmax(arraysize),                 &
         sum_casaflux%klitter(arraysize,mlitter),          &
         sum_casaflux%ksoil(arraysize,msoil),              &
         sum_casaflux%fromLtoS(arraysize,msoil,mlitter),   &
         sum_casaflux%fromStoS(arraysize,msoil,msoil),     &
         sum_casaflux%fromLtoCO2(arraysize,mlitter),       &
         sum_casaflux%fromStoCO2(arraysize,msoil),         &
         sum_casaflux%stemnpp(arraysize),                  &
         sum_casaflux%frac_sapwood(arraysize),             &
         sum_casaflux%sapwood_area(arraysize), &
         sum_casaflux%Cplant_turnover(arraysize,mplant) , &
         sum_casaflux%Cplant_turnover_disturbance(arraysize) , &
         sum_casaflux%Cplant_turnover_crowding(arraysize) , &
         sum_casaflux%Cplant_turnover_resource_limitation(arraysize))

    ALLOCATE(sum_casaflux%FluxCtolitter(arraysize,mlitter),    &
         sum_casaflux%FluxNtolitter(arraysize,mlitter),    &
         sum_casaflux%FluxPtolitter(arraysize,mlitter))

    ALLOCATE(sum_casaflux%FluxCtosoil(arraysize,msoil),        &
         sum_casaflux%FluxNtosoil(arraysize,msoil),        &
         sum_casaflux%FluxPtosoil(arraysize,msoil))

    ALLOCATE(sum_casaflux%fromPtoL_fire(arraysize,mlitter,mplant), &
         sum_casaflux%kplant_fire(arraysize,mplant),        &
         sum_casaflux%klitter_fire(arraysize,mlitter),           &
         sum_casaflux%kplant_tot(arraysize,mplant),        &
         sum_casaflux%klitter_tot(arraysize,mlitter),           &
         sum_casaflux%FluxCtoCO2_plant_fire(arraysize),             &
         sum_casaflux%FluxCtoCO2_litter_fire(arraysize),             &
         sum_casaflux%fluxfromPtoCO2_fire(arraysize,mplant),         &
         sum_casaflux%fluxfromLtoCO2_fire(arraysize,mlitter),        &
         sum_casaflux%FluxNtoAtm_fire(arraysize))

    ALLOCATE(sum_casaflux%FluxFromPtoL(arraysize,mplant,mlitter), &
         sum_casaflux%FluxFromLtoS(arraysize,mlitter,msoil), &
         sum_casaflux%FluxFromStoS(arraysize,msoil,msoil), &
         sum_casaflux%FluxFromPtoCO2(arraysize,mplant), &
         sum_casaflux%FluxFromLtoCO2(arraysize,mlitter), &
         sum_casaflux%FluxFromStoCO2(arraysize,msoil), &
         sum_casaflux%FluxFromPtoHarvest(arraysize))

    ALLOCATE(sum_casaflux%FluxCtoco2(arraysize))

  END SUBROUTINE alloc_sum_casavariable


  subroutine zero_casavariable(casabiome, casapool, casaflux, casamet, casabal)
    
    implicit none
    
    type(casa_biome)  , intent(inout) :: casabiome
    type(casa_pool)   , intent(inout) :: casapool
    type(casa_flux)   , intent(inout) :: casaflux
    type(casa_met)    , intent(inout) :: casamet
    type(casa_balance), intent(inout) :: casabal

    casabiome%ivt2                 = 0.0_r_2
    casabiome%xkleafcoldmax        = 0.0_r_2
    casabiome%xkleafcoldexp        = 0.0_r_2
    casabiome%xkleafdrymax         = 0.0_r_2
    casabiome%xkleafdryexp         = 0.0_r_2
    casabiome%glaimax              = 0.0_r_2
    casabiome%glaimin              = 0.0_r_2
    casabiome%sla                  = 0.0_r_2
    casabiome%ratiofrootleaf       = 0.0_r_2
    casabiome%kroot                = 0.0_r_2
    casabiome%krootlen             = 0.0_r_2
    casabiome%rootdepth            = 0.0_r_2
    casabiome%kuptake              = 0.0_r_2
    casabiome%kminN                = 0.0_r_2
    casabiome%KuplabP              = 0.0_r_2
    casabiome%kclabrate            = 0.0_r_2
    casabiome%xnpmax               = 0.0_r_2
    casabiome%q10soil              = 0.0_r_2
    casabiome%xkoptlitter          = 0.0_r_2
    casabiome%xkoptsoil            = 0.0_r_2
    casabiome%xkplab               = 0.0_r_2
    casabiome%xkpsorb              = 0.0_r_2
    casabiome%xkpocc               = 0.0_r_2
    casabiome%prodptase            = 0.0_r_2
    casabiome%costnpup             = 0.0_r_2
    casabiome%maxfinelitter        = 0.0_r_2
    casabiome%maxcwd               = 0.0_r_2
    casabiome%nintercept           = 0.0_r_2
    casabiome%nslope               = 0.0_r_2
    casabiome%plantrate            = 0.0_r_2
    casabiome%rmplant              = 0.0_r_2
    casabiome%fracnpptoP           = 0.0_r_2
    casabiome%fraclignin           = 0.0_r_2
    casabiome%fraclabile           = 0.0_r_2
    casabiome%ratioNCplantmin      = 0.0_r_2
    casabiome%ratioNCplantmax      = 0.0_r_2
    casabiome%ratioNPplantmin      = 0.0_r_2
    casabiome%ratioNPplantmax      = 0.0_r_2
    casabiome%fracLigninplant      = 0.0_r_2
    casabiome%ftransNPtoL          = 0.0_r_2
    casabiome%ftransPPtoL          = 0.0_r_2
    casabiome%litterrate           = 0.0_r_2
    casabiome%soilrate             = 0.0_r_2
    casabiome%ratioPcplantmax      = 0.0_r_2
    casabiome%ratioPcplantmin      = 0.0_r_2
    casabiome%la_to_sa             = 0.0_r_2
    casabiome%vcmax_scalar         = 0.0_r_2
    casabiome%disturbance_interval = 0.0_r_2
    casabiome%DAMM_EnzPool         = 0.0_r_2
    casabiome%DAMM_KMO2            = 0.0_r_2
    casabiome%DAMM_KMcp            = 0.0_r_2
    casabiome%DAMM_Ea              = 0.0_r_2
    casabiome%DAMM_alpha           = 0.0_r_2

    casapool%Clabile        = 0.0_r_2
    casapool%dClabiledt     = 0.0_r_2
    casapool%Cplant         = 0.0_r_2
    casapool%Nplant         = 0.0_r_2
    casapool%Pplant         = 0.0_r_2
    casapool%dCplantdt      = 0.0_r_2
    casapool%dNplantdt      = 0.0_r_2
    casapool%dPplantdt      = 0.0_r_2
    casapool%ratioNCplant   = 0.0_r_2
    casapool%ratioNPplant   = 0.0_r_2
    casapool%Nsoilmin       = 0.0_r_2
    casapool%Psoillab       = 0.0_r_2
    casapool%Psoilsorb      = 0.0_r_2
    casapool%Psoilocc       = 0.0_r_2
    casapool%dNsoilmindt    = 0.0_r_2
    casapool%dPsoillabdt    = 0.0_r_2
    casapool%dPsoilsorbdt   = 0.0_r_2
    casapool%dPsoiloccdt    = 0.0_r_2
    casapool%Clitter        = 0.0_r_2
    casapool%Nlitter        = 0.0_r_2
    casapool%Plitter        = 0.0_r_2
    casapool%dClitterdt     = 0.0_r_2
    casapool%dNlitterdt     = 0.0_r_2
    casapool%dPlitterdt     = 0.0_r_2
    casapool%ratioNClitter  = 0.0_r_2
    casapool%ratioNPlitter  = 0.0_r_2
    casapool%Csoil          = 0.0_r_2
    casapool%Nsoil          = 0.0_r_2
    casapool%Psoil          = 0.0_r_2
    casapool%dCsoildt       = 0.0_r_2
    casapool%dNsoildt       = 0.0_r_2
    casapool%dPsoildt       = 0.0_r_2
    casapool%ratioNCsoil    = 0.0_r_2
    casapool%ratioNPsoil    = 0.0_r_2
    casapool%ratioNCsoilnew = 0.0_r_2
    casapool%ratioNCsoilmin = 0.0_r_2
    casapool%ratioNCsoilmax = 0.0_r_2
    casapool%ratioPCsoil    = 0.0_r_2
    casapool%ratioPCplant   = 0.0_r_2
    casapool%ratioPClitter  = 0.0_r_2
    casapool%Ctot_0         = 0.0_r_2
    casapool%Ctot           = 0.0_r_2

    casaflux%Cgpp         = 0.0_r_2
    casaflux%Cnpp         = 0.0_r_2
    casaflux%Crp          = 0.0_r_2
    casaflux%Crgplant     = 0.0_r_2
    casaflux%Nminfix      = 0.0_r_2
    casaflux%Nminuptake   = 0.0_r_2
    casaflux%Plabuptake   = 0.0_r_2
    casaflux%Clabloss     = 0.0_r_2
    casaflux%fracClabile  = 0.0_r_2
    casaflux%fracCalloc   = 0.0_r_2
    casaflux%fracNalloc   = 0.0_r_2
    casaflux%fracPalloc   = 0.0_r_2
    casaflux%kplant       = 0.0_r_2
    casaflux%Crmplant     = 0.0_r_2
    casaflux%fromPtoL     = 0.0_r_2
    casaflux%Cnep         = 0.0_r_2
    casaflux%Crsoil       = 0.0_r_2
    casaflux%Nmindep      = 0.0_r_2
    casaflux%Nminloss     = 0.0_r_2
    casaflux%Nminleach    = 0.0_r_2
    casaflux%Nupland      = 0.0_r_2
    casaflux%Nlittermin   = 0.0_r_2
    casaflux%Nsmin        = 0.0_r_2
    casaflux%Nsimm        = 0.0_r_2
    casaflux%Nsnet        = 0.0_r_2
    casaflux%fNminloss    = 0.0_r_2
    casaflux%fNminleach   = 0.0_r_2
    casaflux%Pdep         = 0.0_r_2
    casaflux%Pwea         = 0.0_r_2
    casaflux%Pleach       = 0.0_r_2
    casaflux%Ploss        = 0.0_r_2
    casaflux%Pupland      = 0.0_r_2
    casaflux%Plittermin   = 0.0_r_2
    casaflux%Psmin        = 0.0_r_2
    casaflux%Psimm        = 0.0_r_2
    casaflux%Psnet        = 0.0_r_2
    casaflux%fPleach      = 0.0_r_2
    casaflux%kplab        = 0.0_r_2
    casaflux%kpsorb       = 0.0_r_2
    casaflux%kpocc        = 0.0_r_2
    casaflux%kmlabP       = 0.0_r_2
    casaflux%Psorbmax     = 0.0_r_2
    casaflux%klitter      = 0.0_r_2
    casaflux%ksoil        = 0.0_r_2
    casaflux%fromLtoS     = 0.0_r_2
    casaflux%fromStoS     = 0.0_r_2
    casaflux%fromLtoCO2   = 0.0_r_2
    casaflux%fromStoCO2   = 0.0_r_2
    casaflux%stemnpp      = 0.0_r_2
    casaflux%frac_sapwood = 0.0_r_2
    casaflux%sapwood_area = 0.0_r_2
    casaflux%fharvest     = 0.0_r_2
    casaflux%Charvest     = 0.0_r_2
    casaflux%Nharvest     = 0.0_r_2
    casaflux%Pharvest     = 0.0_r_2
    casaflux%fcrop        = 0.0_r_2
    casaflux%Cplant_turnover                     = 0.0_r_2
    casaflux%Cplant_turnover_disturbance         = 0.0_r_2
    casaflux%Cplant_turnover_crowding            = 0.0_r_2
    casaflux%Cplant_turnover_resource_limitation = 0.0_r_2

    casaflux%fromPtoL_fire          = 0.0_r_2
    casaflux%kplant_fire            = 0.0_r_2
    casaflux%klitter_fire           = 0.0_r_2
    casaflux%kplant_tot             = 0.0_r_2
    casaflux%klitter_tot            = 0.0_r_2
    casaflux%FluxCtoCO2_plant_fire  = 0.0_r_2
    casaflux%FluxCtoCO2_litter_fire = 0.0_r_2
    casaflux%fluxfromPtoCO2_fire    = 0.0_r_2
    casaflux%fluxfromLtoCO2_fire    = 0.0_r_2
    casaflux%FluxNtoAtm_fire        = 0.0_r_2

    casaflux%FluxCtolitter = 0.0_r_2
    casaflux%FluxNtolitter = 0.0_r_2
    casaflux%FluxPtolitter = 0.0_r_2

    casaflux%FluxCtosoil = 0.0_r_2
    casaflux%FluxNtosoil = 0.0_r_2
    casaflux%FluxPtosoil = 0.0_r_2

    casaflux%FluxCtohwp = 0.0_r_2
    casaflux%FluxNtohwp = 0.0_r_2
    casaflux%FluxPtohwp = 0.0_r_2

    casaflux%FluxCtoclear = 0.0_r_2
    casaflux%FluxNtoclear = 0.0_r_2
    casaflux%FluxPtoclear = 0.0_r_2

    casaflux%CtransferLUC = 0.0_r_2

    casaflux%FluxCtoco2 = 0.0_r_2

    casaflux%FluxFromPtoL = 0.0_r_2
    casaflux%FluxFromLtoS       = 0.0_r_2
    casaflux%FluxFromStoS       = 0.0_r_2
    casaflux%FluxFromPtoCO2     = 0.0_r_2
    casaflux%FluxFromLtoCO2     = 0.0_r_2
    casaflux%FluxFromStoCO2     = 0.0_r_2
    casaflux%FluxFromPtoHarvest = 0.0_r_2

    casamet%glai           = 0.0_r_2
    casamet%lnonwood       = 0.0_r_2
    casamet%Tairk          = 0.0_r_2
    casamet%precip         = 0.0_r_2
    casamet%tsoilavg       = 0.0_r_2
    casamet%moistavg       = 0.0_r_2
    casamet%btran          = 0.0_r_2
    casamet%Tsoil          = 0.0_r_2
    casamet%moist          = 0.0_r_2
    casamet%iveg2          = 0.0_r_2
    casamet%ijgcm          = 0.0_r_2
    casamet%isorder        = 0.0_r_2
    casamet%lat            = 0.0_r_2
    casamet%lon            = 0.0_r_2
    casamet%areacell       = 0.0_r_2
    casamet%Tairkspin      = 0.0_r_2
    casamet%cgppspin       = 0.0_r_2
    casamet%crmplantspin_1 = 0.0_r_2
    casamet%crmplantspin_2 = 0.0_r_2
    casamet%crmplantspin_3 = 0.0_r_2
    casamet%Tsoilspin_1    = 0.0_r_2
    casamet%Tsoilspin_2    = 0.0_r_2
    casamet%Tsoilspin_3    = 0.0_r_2
    casamet%Tsoilspin_4    = 0.0_r_2
    casamet%Tsoilspin_5    = 0.0_r_2
    casamet%Tsoilspin_6    = 0.0_r_2
    casamet%moistspin_1    = 0.0_r_2
    casamet%moistspin_2    = 0.0_r_2
    casamet%moistspin_3    = 0.0_r_2
    casamet%moistspin_4    = 0.0_r_2
    casamet%moistspin_5    = 0.0_r_2
    casamet%moistspin_6    = 0.0_r_2
    casamet%mtempspin      = 0.0_r_2
    casamet%cAn12spin      = 0.0_r_2
    casamet%cAn13spin      = 0.0_r_2

    casabal%FCgppyear    = 0.0_r_2
    casabal%FCnppyear    = 0.0_r_2
    casabal%FCrpyear     = 0.0_r_2
    casabal%FCrmleafyear = 0.0_r_2
    casabal%FCrmwoodyear = 0.0_r_2
    casabal%FCrmrootyear = 0.0_r_2
    casabal%FCrgrowyear  = 0.0_r_2
    casabal%FCrsyear     = 0.0_r_2
    casabal%FCneeyear    = 0.0_r_2
    casabal%FNdepyear    = 0.0_r_2
    casabal%FNfixyear    = 0.0_r_2
    casabal%FNsnetyear   = 0.0_r_2
    casabal%FNupyear     = 0.0_r_2
    casabal%FNleachyear  = 0.0_r_2
    casabal%FNlossyear   = 0.0_r_2
    casabal%FPweayear    = 0.0_r_2
    casabal%FPdustyear   = 0.0_r_2
    casabal%FPsnetyear   = 0.0_r_2
    casabal%FPupyear     = 0.0_r_2
    casabal%FPleachyear  = 0.0_r_2
    casabal%FPlossyear   = 0.0_r_2
    casabal%dCdtyear     = 0.0_r_2
    casabal%LAImax       = 0.0_r_2
    casabal%Cleafmean    = 0.0_r_2
    casabal%Crootmean    = 0.0_r_2

    casabal%glaimon  = 0.0_r_2
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

    casabal%nsoilminlast  = 0.0_r_2
    casabal%psoillablast  = 0.0_r_2
    casabal%psoilsorblast = 0.0_r_2
    casabal%psoilocclast  = 0.0_r_2
    casabal%cbalance      = 0.0_r_2
    casabal%nbalance      = 0.0_r_2
    casabal%pbalance      = 0.0_r_2
    casabal%sumcbal       = 0.0_r_2
    casabal%sumnbal       = 0.0_r_2
    casabal%sumpbal       = 0.0_r_2
    casabal%clabilelast   = 0.0_r_2

  END SUBROUTINE zero_casavariable


  SUBROUTINE zero_sum_casa(sum_casapool, sum_casaflux)

    IMPLICIT NONE

    TYPE(casa_pool), INTENT(INOUT) :: sum_casapool
    TYPE(casa_flux), INTENT(INOUT) :: sum_casaflux

    sum_casapool%Clabile        = 0.0_r_2
    sum_casapool%dClabiledt     = 0.0_r_2
    sum_casapool%Cplant         = 0.0_r_2
    sum_casapool%Nplant         = 0.0_r_2
    sum_casapool%Pplant         = 0.0_r_2
    sum_casapool%dCplantdt      = 0.0_r_2
    sum_casapool%dNplantdt      = 0.0_r_2
    sum_casapool%dPplantdt      = 0.0_r_2
    sum_casapool%ratioNCplant   = 0.0_r_2
    sum_casapool%ratioNPplant   = 0.0_r_2
    sum_casapool%Nsoilmin       = 0.0_r_2
    sum_casapool%Psoillab       = 0.0_r_2
    sum_casapool%Psoilsorb      = 0.0_r_2
    sum_casapool%Psoilocc       = 0.0_r_2
    sum_casapool%dNsoilmindt    = 0.0_r_2
    sum_casapool%dPsoillabdt    = 0.0_r_2
    sum_casapool%dPsoilsorbdt   = 0.0_r_2
    sum_casapool%dPsoiloccdt    = 0.0_r_2
    sum_casapool%Clitter        = 0.0_r_2
    sum_casapool%Nlitter        = 0.0_r_2
    sum_casapool%Plitter        = 0.0_r_2
    sum_casapool%dClitterdt     = 0.0_r_2
    sum_casapool%dNlitterdt     = 0.0_r_2
    sum_casapool%dPlitterdt     = 0.0_r_2
    sum_casapool%ratioNClitter  = 0.0_r_2
    sum_casapool%ratioNPlitter  = 0.0_r_2
    sum_casapool%Csoil          = 0.0_r_2
    sum_casapool%Nsoil          = 0.0_r_2
    sum_casapool%Psoil          = 0.0_r_2
    sum_casapool%dCsoildt       = 0.0_r_2
    sum_casapool%dNsoildt       = 0.0_r_2
    sum_casapool%dPsoildt       = 0.0_r_2
    sum_casapool%ratioNCsoil    = 0.0_r_2
    sum_casapool%ratioNPsoil    = 0.0_r_2
    sum_casapool%ratioNCsoilnew = 0.0_r_2
    sum_casapool%ratioNCsoilmin = 0.0_r_2
    sum_casapool%ratioNCsoilmax = 0.0_r_2
    sum_casapool%ratioPCsoil    = 0.0_r_2
    sum_casapool%ratioPCplant   = 0.0_r_2
    sum_casapool%ratioPClitter  = 0.0_r_2

    sum_casaflux%Cgpp         = 0.0_r_2
    sum_casaflux%Cnpp         = 0.0_r_2
    sum_casaflux%Crp          = 0.0_r_2
    sum_casaflux%Crgplant     = 0.0_r_2
    sum_casaflux%Nminfix      = 0.0_r_2
    sum_casaflux%Nminuptake   = 0.0_r_2
    sum_casaflux%Plabuptake   = 0.0_r_2
    sum_casaflux%Clabloss     = 0.0_r_2
    sum_casaflux%fracClabile  = 0.0_r_2
    sum_casaflux%fracCalloc   = 0.0_r_2
    sum_casaflux%fracNalloc   = 0.0_r_2
    sum_casaflux%fracPalloc   = 0.0_r_2
    sum_casaflux%kplant       = 0.0_r_2
    sum_casaflux%Crmplant     = 0.0_r_2
    sum_casaflux%fromPtoL     = 0.0_r_2
    sum_casaflux%Cnep         = 0.0_r_2
    sum_casaflux%Crsoil       = 0.0_r_2
    sum_casaflux%Nmindep      = 0.0_r_2
    sum_casaflux%Nminloss     = 0.0_r_2
    sum_casaflux%Nminleach    = 0.0_r_2
    sum_casaflux%Nupland      = 0.0_r_2
    sum_casaflux%Nlittermin   = 0.0_r_2
    sum_casaflux%Nsmin        = 0.0_r_2
    sum_casaflux%Nsimm        = 0.0_r_2
    sum_casaflux%Nsnet        = 0.0_r_2
    sum_casaflux%fNminloss    = 0.0_r_2
    sum_casaflux%fNminleach   = 0.0_r_2
    sum_casaflux%Pdep         = 0.0_r_2
    sum_casaflux%Pwea         = 0.0_r_2
    sum_casaflux%Pleach       = 0.0_r_2
    sum_casaflux%Ploss        = 0.0_r_2
    sum_casaflux%Pupland      = 0.0_r_2
    sum_casaflux%Plittermin   = 0.0_r_2
    sum_casaflux%Psmin        = 0.0_r_2
    sum_casaflux%Psimm        = 0.0_r_2
    sum_casaflux%Psnet        = 0.0_r_2
    sum_casaflux%fPleach      = 0.0_r_2
    sum_casaflux%kplab        = 0.0_r_2
    sum_casaflux%kpsorb       = 0.0_r_2
    sum_casaflux%kpocc        = 0.0_r_2
    sum_casaflux%kmlabP       = 0.0_r_2
    sum_casaflux%Psorbmax     = 0.0_r_2
    sum_casaflux%klitter      = 0.0_r_2
    sum_casaflux%ksoil        = 0.0_r_2
    sum_casaflux%fromLtoS     = 0.0_r_2
    sum_casaflux%fromStoS     = 0.0_r_2
    sum_casaflux%fromLtoCO2   = 0.0_r_2
    sum_casaflux%fromStoCO2   = 0.0_r_2
    sum_casaflux%stemnpp      = 0.0_r_2
    sum_casaflux%frac_sapwood = 0.0_r_2
    sum_casaflux%sapwood_area = 0.0_r_2
    sum_casaflux%Cplant_turnover                     = 0.0_r_2
    sum_casaflux%Cplant_turnover_disturbance         = 0.0_r_2
    sum_casaflux% Cplant_turnover_crowding           = 0.0_r_2
    sum_casaflux%Cplant_turnover_resource_limitation = 0.0_r_2

    sum_casaflux%FluxCtolitter = 0.0_r_2
    sum_casaflux%FluxNtolitter = 0.0_r_2
    sum_casaflux%FluxPtolitter = 0.0_r_2

    sum_casaflux%FluxCtosoil = 0.0_r_2
    sum_casaflux%FluxNtosoil = 0.0_r_2
    sum_casaflux%FluxPtosoil = 0.0_r_2

    sum_casaflux%FluxCtoco2 = 0.0_r_2

    sum_casaflux%fromPtoL_fire          = 0.0_r_2
    sum_casaflux%kplant_fire            = 0.0_r_2
    sum_casaflux%klitter_fire           = 0.0_r_2
    sum_casaflux%kplant_tot             = 0.0_r_2
    sum_casaflux%klitter_tot            = 0.0_r_2
    sum_casaflux%FluxCtoCO2_plant_fire  = 0.0_r_2
    sum_casaflux%FluxCtoCO2_litter_fire = 0.0_r_2
    sum_casaflux%fluxfromPtoCO2_fire    = 0.0_r_2
    sum_casaflux%fluxfromLtoCO2_fire    = 0.0_r_2
    sum_casaflux%FluxNtoAtm_fire        = 0.0_r_2

    sum_casaflux%FluxFromPtoL       = 0.0_r_2
    sum_casaflux%FluxFromLtoS       = 0.0_r_2
    sum_casaflux%FluxFromStoS       = 0.0_r_2
    sum_casaflux%FluxFromPtoCO2     = 0.0_r_2
    sum_casaflux%FluxFromLtoCO2     = 0.0_r_2
    sum_casaflux%FluxFromStoCO2     = 0.0_r_2
    sum_casaflux%FluxFromPtoHarvest = 0.0_r_2

  END SUBROUTINE zero_sum_casa


  SUBROUTINE update_sum_casa(sum_casapool, sum_casaflux, casapool, casaflux, sum_now, average_now, nsteps)

    IMPLICIT NONE
    
    TYPE(casa_pool), INTENT(INOUT) :: sum_casapool
    TYPE(casa_flux), INTENT(INOUT) :: sum_casaflux
    TYPE(casa_pool), INTENT(IN)    :: casapool
    TYPE(casa_flux), INTENT(IN)    :: casaflux
    LOGICAL,         INTENT(IN)    :: sum_now, average_now
    INTEGER,         INTENT(IN)    :: nsteps

    real :: rnsteps

    rnsteps = 1.0 / real(nsteps)

    IF (sum_now) then
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
       sum_casaflux%Cplant_turnover                     = sum_casaflux%Cplant_turnover + &
            casaflux%Cplant_turnover
       sum_casaflux%Cplant_turnover_disturbance         = sum_casaflux%Cplant_turnover_disturbance + &
            casaflux%Cplant_turnover_disturbance
       sum_casaflux%Cplant_turnover_crowding            = sum_casaflux%Cplant_turnover_crowding + &
            casaflux%Cplant_turnover_crowding
       sum_casaflux%Cplant_turnover_resource_limitation = sum_casaflux%Cplant_turnover_resource_limitation +  &
            casaflux%Cplant_turnover_resource_limitation

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
       where (sum_casaflux%Cnpp .gt. 1.e-12_r_2)
          sum_casaflux%fracCalloc(:,1) = sum_casaflux%fracCalloc(:,1) / sum_casaflux%Cnpp
          sum_casaflux%fracCalloc(:,2) = sum_casaflux%fracCalloc(:,2) / sum_casaflux%Cnpp
          sum_casaflux%fracCalloc(:,3) = sum_casaflux%fracCalloc(:,3) / sum_casaflux%Cnpp
       elsewhere
          sum_casaflux%fracCalloc(:,1) = 0.0_r_2
          sum_casaflux%fracCalloc(:,2) = 0.0_r_2
          sum_casaflux%fracCalloc(:,3) = 0.0_r_2
       endwhere
       ! sum_casaflux%kplant = sum_casaflux%kplant * rnsteps
       where (sum_casapool%Cplant .gt. 1.e-12_r_2)
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
       sum_casaflux%Cplant_turnover             = sum_casaflux%Cplant_turnover          * rnsteps
       sum_casaflux%Cplant_turnover_disturbance = casaflux%Cplant_turnover_disturbance  * rnsteps
       sum_casaflux%Cplant_turnover_crowding    = sum_casaflux%Cplant_turnover_crowding * rnsteps
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

  END SUBROUTINE update_sum_casa


END MODULE casavariable


MODULE phenvariable
  
  USE casadimension
  
  IMPLICIT NONE
  
  TYPE phen_variable
     INTEGER,   DIMENSION(:),  POINTER :: phase => null()
     REAL(r_2), DIMENSION(:),  POINTER :: TKshed => null()
     INTEGER,   DIMENSION(:,:),POINTER :: doyphase => null()
     REAL, DIMENSION(:),  POINTER :: phen => null()   ! fraction of max LAI
     REAL, DIMENSION(:),  POINTER :: aphen => null()  ! annual leaf on sum
     INTEGER,   DIMENSION(:,:),POINTER :: phasespin => null()
     INTEGER,   DIMENSION(:,:),POINTER :: doyphasespin_1 => null()
     INTEGER,   DIMENSION(:,:),POINTER :: doyphasespin_2 => null()
     INTEGER,   DIMENSION(:,:),POINTER :: doyphasespin_3 => null()
     INTEGER,   DIMENSION(:,:),POINTER :: doyphasespin_4 => null()
  END type phen_variable

  
CONTAINS

  
  SUBROUTINE alloc_phenvariable(phen,arraysize)

    IMPLICIT NONE
    
    TYPE(phen_variable), INTENT(INOUT) :: phen
    INTEGER,             INTENT(IN) :: arraysize

    ALLOCATE(phen%Tkshed(mvtype))
    ALLOCATE(phen%phase(arraysize),         &
         phen%doyphase(arraysize,mphase))
    ALLOCATE(phen%phen(arraysize), &
         phen%aphen(arraysize), &
         phen%phasespin(arraysize,mdyear), &
         phen%doyphasespin_1(arraysize,mdyear), &
         phen%doyphasespin_2(arraysize,mdyear), &
         phen%doyphasespin_3(arraysize,mdyear), &
         phen%doyphasespin_4(arraysize,mdyear))
    
  END SUBROUTINE alloc_phenvariable

  
  subroutine zero_phenvariable(phen)

    implicit none
    
    type(phen_variable), intent(inout) :: phen

    phen%Tkshed         = 0
    phen%phase          = 0.0_r_2
    phen%doyphase       = 0
    phen%phen           = 0.0
    phen%aphen          = 0.0
    phen%phasespin      = 0
    phen%doyphasespin_1 = 0
    phen%doyphasespin_2 = 0
    phen%doyphasespin_3 = 0
    phen%doyphasespin_4 = 0
    
  end subroutine zero_phenvariable
  
End MODULE phenvariable
