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
  INTEGER, PARAMETER :: mplant = 3         ! plant pools
  INTEGER, PARAMETER :: mlitter= 3         ! litter pools
  INTEGER, PARAMETER :: msoil  = 3         ! soil pools
  INTEGER, PARAMETER :: mso    = 12        ! soil order number
! BP put icycle into namelist file
  INTEGER            :: icycle
!  INTEGER, PARAMETER :: icycle=3           ! =1 for C, =2 for C+N; =3 for C+N+P
  INTEGER, PARAMETER :: mstart=1           ! starting time step
  INTEGER, PARAMETER :: mphase=4           ! phen. phases
  REAL(r_2),    PARAMETER :: deltcasa=1.0/365.0 ! year
  REAL(r_2),    PARAMETER :: deltpool=1.0       ! pool delt(1day)

END MODULE casadimension

MODULE casaparm
  USE casadimension

  IMPLICIT NONE
  INTEGER, PARAMETER :: initcasa= 1   ! =0 spin; 1 restart file
  INTEGER, PARAMETER :: iceland  = 17 !=13 for casa vegtype =15 for IGBP vegtype
  INTEGER, PARAMETER :: cropland = 9  ! 12 and 14 for IGBP vegtype 
  INTEGER, PARAMETER :: croplnd2 =10  ! ditto
  INTEGER, PARAMETER :: forest  = 3
  INTEGER, PARAMETER :: shrub   = 2
  INTEGER, PARAMETER :: grass   = 1
  INTEGER, PARAMETER :: icewater= 0
  INTEGER, PARAMETER :: LEAF    = 1
  INTEGER, PARAMETER :: WOOD    = 2
  INTEGER, PARAMETER :: FROOT   = 3
!  INTEGER, PARAMETER :: LABILE  = 4
  INTEGER, PARAMETER :: METB    = 1
  INTEGER, PARAMETER :: STR     = 2
  INTEGER, PARAMETER :: CWD     = 3
  INTEGER, PARAMETER :: MIC     = 1
  INTEGER, PARAMETER :: SLOW    = 2
  INTEGER, PARAMETER :: PASS    = 3
  INTEGER, PARAMETER :: PLAB    = 1
  INTEGER, PARAMETER :: PSORB   = 2
  INTEGER, PARAMETER :: POCC    = 3
  INTEGER, PARAMETER :: LALLOC  = 0      !=0 constant; 1 variable
  REAL(r_2), PARAMETER :: z30=0.3
  REAL(r_2), PARAMETER :: R0=0.3
  REAL(r_2), PARAMETER :: S0=0.3
  REAL(r_2), PARAMETER :: fixed_stem=1.0/3.0
  REAL(r_2), PARAMETER :: Q10alloc=2.0
  REAL(r_2), PARAMETER :: ratioNCstrfix = 1.0/150.0
  REAL(r_2), PARAMETER :: ratioPCstrfix = ratioNCstrfix/25.0
  REAL(r_2), PARAMETER :: fracCbiomass = 0.50
  REAL(r_2), PARAMETER :: tsoilrefc=25.0
  REAL(r_2), PARAMETER :: tkzeroc=273.15
  REAL(r_2), PARAMETER :: frootparma = 0.3192
  REAL(r_2), PARAMETER :: frootparmb =-0.0485
  REAL(r_2), PARAMETER :: frootparmc = 0.1755
  REAL(r_2), PARAMETER :: xweightalloc = 0.2
!  REAL(r_2), PARAMETER :: xkplab=0.5*deltcasa
!  REAL(r_2), PARAMETER :: xkpsorb=0.01*deltcasa
!  REAL(r_2), PARAMETER :: xkpocc =0.01*deltcasa
END MODULE casaparm

MODULE casavariable
  USE casadimension
  IMPLICIT NONE

  TYPE casa_biome
    INTEGER,   DIMENSION(:),POINTER :: ivt2
    REAL(r_2), DIMENSION(:),POINTER :: xkleafcoldmax,  &
                                       xkleafcoldexp,  &
                                       xkleafdrymax,   &
                                       xkleafdryexp,   &
                                       glaimax,        &
                                       glaimin,        &
                                       sla,            &
                                       ratiofrootleaf, &
                                       kroot,          &
                                       krootlen,       &
                                       rootdepth,      &
                                       kuptake,        &
                                       kminN,          &
                                       kuplabP,        &
                                       kclabrate,      &
                                       xnpmax,         &
                                       q10soil,        &
                                       xkoptlitter,    &
                                       xkoptsoil,      &
                                       xkplab,         &
                                       xkpsorb,        &
                                       xkpocc,         &
                                       prodptase,      &
                                       costnpup,       &
                                       maxfinelitter,  &
                                       maxcwd,         &             
                                       nintercept,     &  
                                       nslope             

    REAL(r_2), DIMENSION(:,:),POINTER :: plantrate,     &
                                       rmplant,         &
                                       fracnpptoP,      &
                                       fraclignin,      &
                                       fraclabile,      &
                                       ratioNCplantmin, &
                                       ratioNCplantmax, &
                                       ratioPCplantmin, &
                                       ratioPCplantmax, &
                                       fracLigninplant, &
                                       ftransNPtoL,     &
                                       ftransPPtoL,     &
                                       litterrate
    REAL(r_2), DIMENSION(:,:),POINTER :: soilrate
  END TYPE casa_biome

  TYPE casa_pool
    REAL(r_2), DIMENSION(:),POINTER :: Clabile,       &
                                       dClabiledt               
    REAL(r_2), DIMENSION(:,:),POINTER :: Cplant,      &
                                       Nplant,        &
                                       Pplant,        &
                                       dCplantdt,     &
                                       dNplantdt,     &
                                       dPplantdt,     &
                                       ratioNCplant,  &
                                       ratioPCplant
    REAL(r_2), DIMENSION(:),POINTER :: Nsoilmin,      &
                                       Psoillab,      &
                                       Psoilsorb,     &
                                       Psoilocc,      & 
                                       dNsoilmindt,   &
                                       dPsoillabdt,   &
                                       dPsoilsorbdt,  &
                                       dPsoiloccdt          
    REAL(r_2), DIMENSION(:,:), POINTER :: Clitter,    &
                                       Nlitter,       &
                                       Plitter,       &
                                       dClitterdt,    &
                                       dNlitterdt,    &
                                       dPlitterdt,    &
                                       ratioNClitter, &
                                       ratioPClitter
    REAL(r_2), DIMENSION(:,:),POINTER :: Csoil,       &
                                       Nsoil,         &
                                       Psoil,         &
                                       dCsoildt,      &
                                       dNsoildt,      &
                                       dPsoildt,      &
                                       ratioNCsoil,   &
                                       ratioNCsoilnew,&
                                       ratioPCsoil,   &
                                       ratioNCsoilmin,&
                                       ratioNCsoilmax
  END TYPE casa_pool

  TYPE casa_flux
    REAL(r_2), DIMENSION(:),POINTER :: Cgpp,          &
                                       Cnpp,          &
                                       Crp,           &
                                       Crgplant,      &
                                       Nminfix,       &
                                       Nminuptake,    &
                                       Plabuptake,    &
                                       Clabloss,      &
                                       fracClabile
    REAL(r_2), DIMENSION(:,:),POINTER :: fracCalloc,  &
                                       fracNalloc,    &
                                       fracPalloc,    &
                                       Crmplant,      &
                                       kplant
    REAL(r_2), DIMENSION(:,:,:),POINTER :: fromPtoL
    REAL(r_2), DIMENSION(:),POINTER :: Cnep,        &
                                       Crsoil,      &
                                       Nmindep,     &
                                       Nminloss,    &
                                       Nminleach,   &
                                       Nupland,     &
                                       Nlittermin,  &
                                       Nsmin,       &
                                       Nsimm,       &
                                       Nsnet,       &
                                       fNminloss,   &
                                       fNminleach,  &
                                       Pdep,        &
                                       Pwea,        &
                                       Pleach,      &
                                       Ploss,       &
                                       Pupland,     &
                                       Plittermin,  &
                                       Psmin,       &
                                       Psimm,       &
                                       Psnet,       &
                                       fPleach,     &
                                       kplab,       &
                                       kpsorb,      &
                                       kpocc,       &
                                       kmlabp,      &
                                       Psorbmax
    REAL(r_2), DIMENSION(:,:),POINTER    :: klitter
    REAL(r_2), DIMENSION(:,:),POINTER    :: ksoil
    REAL(r_2), DIMENSION(:,:,:),POINTER  :: fromLtoS
    REAL(r_2), DIMENSION(:,:,:),POINTER  :: fromStoS
    REAL(r_2), DIMENSION(:,:),POINTER    :: fromLtoCO2
    REAL(r_2), DIMENSION(:,:),POINTER    :: fromStoCO2
    REAL(r_2), DIMENSION(:,:),POINTER    :: FluxCtolitter
    REAL(r_2), DIMENSION(:,:),POINTER    :: FluxNtolitter
    REAL(r_2), DIMENSION(:,:),POINTER    :: FluxPtolitter
    REAL(r_2), DIMENSION(:,:),POINTER    :: FluxCtosoil
    REAL(r_2), DIMENSION(:,:),POINTER    :: FluxNtosoil
    REAL(r_2), DIMENSION(:,:),POINTER    :: FluxPtosoil
    REAL(r_2), DIMENSION(:),POINTER      :: FluxCtoCO2
  END TYPE casa_flux

  TYPE casa_met
    REAL(r_2), DIMENSION(:),POINTER    :: glai,     &
                                          Tairk,    &
                                          precip,   &
                                          tsoilavg, &
                                          moistavg, &
                                          btran
    INTEGER, DIMENSION(:), POINTER     :: lnonwood
    REAL(r_2), DIMENSION(:,:), POINTER :: Tsoil,    &
                                          moist
    INTEGER, DIMENSION(:), POINTER     :: iveg2,    &  
                                          ijgcm,    &
                                          isorder
    REAL(r_2), DIMENSION(:), POINTER   :: lat,      &
                                          lon,      &
                                          areacell
  END TYPE casa_met

  TYPE casa_balance
    REAL(r_2), DIMENSION(:),POINTER   :: FCgppyear,FCnppyear,                 &
            FCrmleafyear,FCrmwoodyear,FCrmrootyear,FCrgrowyear,               &
            FCrpyear, FCrsyear,FCneeyear,                                     &
            FNdepyear,FNfixyear, FNsnetyear,FNupyear, FNleachyear,FNlossyear, &
            FPweayear,FPdustyear,FPsnetyear,FPupyear, FPleachyear,FPlossyear
    REAL(r_2), DIMENSION(:,:),POINTER :: glaimon,glaimonx
    REAL(r_2), DIMENSION(:,:),POINTER :: cplantlast,nplantlast,pplantlast
    REAL(r_2), DIMENSION(:,:),POINTER :: clitterlast,nlitterlast,plitterlast
    REAL(r_2), DIMENSION(:,:),POINTER :: csoillast,nsoillast,psoillast
    REAL(r_2), DIMENSION(:),  POINTER :: nsoilminlast,psoillablast,  &
                                         psoilsorblast,psoilocclast, &
                                         cbalance,nbalance,pbalance, &
                                         sumcbal,sumnbal,sumpbal
    REAL(r_2), DIMENSION(:),POINTER   :: clabilelast
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
    CHARACTER(LEN=99) :: cnpbiome    ! file for biome-specific BGC parameters
    CHARACTER(LEN=99) :: cnppoint    ! file for point-specific BGC inputs
    CHARACTER(LEN=99) :: cnpepool    ! file for end-of-run pool sizes
    CHARACTER(LEN=99) :: cnpipool    ! file for inital pool sizes
    CHARACTER(LEN=99) :: cnpmetin      ! met file for spin up 
    CHARACTER(LEN=99) :: cnpmetout     ! met file for spin up 
    CHARACTER(LEN=99) :: phen        ! leaf phenology datafile
    CHARACTER(LEN=99) :: cnpflux     ! modelled mean yearly CNP fluxes
  END TYPE casafiles_type
  TYPE(casafiles_type) :: casafile

Contains

SUBROUTINE alloc_casavariable(casabiome,casapool,casaflux,casamet, &
                              casabal,arraysize)
!SUBROUTINE alloc_casavariable(casabiome,casapool,casaflux,casamet, &
!                              casabal,arraysize,mvt)
  IMPLICIT NONE
  TYPE (casa_biome)  , INTENT(INOUT) :: casabiome
  TYPE (casa_pool)   , INTENT(INOUT) :: casapool
  TYPE (casa_flux)   , INTENT(INOUT) :: casaflux
  TYPE (casa_met)    , INTENT(INOUT) :: casamet
  TYPE (casa_balance), INTENT(INOUT) :: casabal
  INTEGER,             INTENT(IN) :: arraysize
!  INTEGER,        INTENT(IN) :: mvt
!  INTEGER :: mvt

!  mvt = mvtype
  ALLOCATE(casabiome%ivt2(mvtype),                   &
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
           casabiome%ratioPCplantmin(mvtype,mplant), &
           casabiome%ratioPCplantmax(mvtype,mplant), &
           casabiome%fracLigninplant(mvtype,mplant), &
           casabiome%ftransNPtoL(mvtype,mplant),     &
           casabiome%ftransPPtoL(mvtype,mplant),     &
           casabiome%litterrate(mvtype,mlitter),     &
           casabiome%soilrate(mvtype,msoil))

  ALLOCATE(casapool%Clabile(arraysize),               &
           casapool%dClabiledt(arraysize),            &
           casapool%Cplant(arraysize,mplant),         &
           casapool%Nplant(arraysize,mplant),         &
           casapool%Pplant(arraysize,mplant),         &
           casapool%dCplantdt(arraysize,mplant),      &
           casapool%dNplantdt(arraysize,mplant),      &
           casapool%dPplantdt(arraysize,mplant),      &
           casapool%ratioNCplant(arraysize,mplant),   &
           casapool%ratioPCplant(arraysize,mplant),   &
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
           casapool%ratioPClitter(arraysize,mlitter), &
           casapool%Csoil(arraysize,msoil),           &
           casapool%Nsoil(arraysize,msoil),           &
           casapool%Psoil(arraysize,msoil),           &
           casapool%dCsoildt(arraysize,msoil),        &
           casapool%dNsoildt(arraysize,msoil),        &
           casapool%dPsoildt(arraysize,msoil),        &
           casapool%ratioNCsoil(arraysize,msoil),     &
           casapool%ratioPCsoil(arraysize,msoil),     &
           casapool%ratioNCsoilnew(arraysize,msoil),  &
           casapool%ratioNCsoilmin(arraysize,msoil),  &
           casapool%ratioNCsoilmax(arraysize,msoil))

  ALLOCATE(casaflux%Cgpp(arraysize),                     &
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
           casaflux%fromStoCO2(arraysize,msoil))

  ALLOCATE(casaflux%FluxCtolitter(arraysize,mlitter),    &
           casaflux%FluxNtolitter(arraysize,mlitter),    &
           casaflux%FluxPtolitter(arraysize,mlitter))

  ALLOCATE(casaflux%FluxCtosoil(arraysize,msoil),        &
           casaflux%FluxNtosoil(arraysize,msoil),        &
           casaflux%FluxPtosoil(arraysize,msoil))

  ALLOCATE(casaflux%FluxCtoco2(arraysize))

  ALLOCATE(casamet%glai(arraysize),                &
           casamet%lnonwood(arraysize),            &
           casamet%Tairk(arraysize),               &
           casamet%precip(arraysize),              &
           casamet%tsoilavg(arraysize),            &
           casamet%moistavg(arraysize),            &
           casamet%btran(arraysize),               &
           casamet%Tsoil(arraysize,ms),            &
           casamet%moist(arraysize,ms),            &
           casamet%iveg2(arraysize),               &  
           casamet%ijgcm(arraysize),               &
           casamet%isorder(arraysize),             &
           casamet%lat(arraysize),                 &
           casamet%lon(arraysize),                 &
           casamet%areacell(arraysize))

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
           casabal%FPlossyear(arraysize))
            
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

END MODULE casavariable


MODULE phenvariable
  USE casadimension
  IMPLICIT NONE
  TYPE phen_variable
    INTEGER,   DIMENSION(:),  POINTER :: phase        
    REAL(r_2), DIMENSION(:),  POINTER :: TKshed
    INTEGER,   DIMENSION(:,:),POINTER :: doyphase
  END type phen_variable

CONTAINS

SUBROUTINE alloc_phenvariable(phen,arraysize)
!SUBROUTINE alloc_phenvariable(phen,arraysize,mvt)
  IMPLICIT NONE
  TYPE(phen_variable), INTENT(INOUT) :: phen
  INTEGER,             INTENT(IN) :: arraysize
!  INTEGER,        INTENT(IN) :: mvt

  ALLOCATE(phen%Tkshed(mvtype))
!  ALLOCATE(phen%Tkshed(mvt))
  ALLOCATE(phen%phase(arraysize),         &
           phen%doyphase(arraysize,mphase))
END SUBROUTINE alloc_phenvariable

End MODULE phenvariable

