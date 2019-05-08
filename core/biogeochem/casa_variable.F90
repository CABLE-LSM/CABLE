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

MODULE casavariable
  USE casadimension
  IMPLICIT NONE

  SAVE

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
          ratioNPplantmin, &
          ratioNPplantmax, &
          fracLigninplant, &
          ftransNPtoL,     &
          ftransPPtoL,     &
          litterrate,      &
          ratioPcplantmin, &
          ratioPcplantmax
     REAL(r_2), DIMENSION(:,:),POINTER :: soilrate
  END TYPE casa_biome

  TYPE casa_pool
     REAL(r_2), DIMENSION(:),POINTER :: Clabile,       &
          dClabiledt,    &
          Ctot ,         &          !! vh_js !!
          Ctot_0
     REAL(r_2), DIMENSION(:,:),POINTER :: Cplant,      &
          Nplant,        &
          Pplant,        &
          dCplantdt,     &
          dNplantdt,     &
          dPplantdt,     &
          ratioNCplant,  &
          ratioNPplant
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
          ratioNPlitter
     REAL(r_2), DIMENSION(:,:),POINTER :: Csoil,       &
          Nsoil,         &
          Psoil,         &
          dCsoildt,      &
          dNsoildt,      &
          dPsoildt,      &
          ratioNCsoil,   &
          ratioNCsoilnew,&
          ratioNPsoil,   &
          ratioNCsoilmin,&
          ratioNCsoilmax,&
          ratioPcsoil,   &
          ratioPcplant,  &
          ratioPclitter
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
          fracClabile, &
                                !! vh_js !! the 3 variables below are needed for POP coupling to CASA
          stemnpp, &
          frac_sapwood, &
          sapwood_area
     REAL(r_2), DIMENSION(:,:),POINTER :: fracCalloc,  &
          fracNalloc,    &
          fracPalloc,    &
          Crmplant,      &
          kplant, &
                                !! vh_js !! additional diagnostic
          Cplant_turnover
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
          Psorbmax,    &
                                !! additional diagnostics for partitioning biomass turnover
          Cplant_turnover_disturbance, &
          Cplant_turnover_crowding , &
          Cplant_turnover_resource_limitation

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
     REAL(r_2), DIMENSION(:),POINTER    :: FluxCtohwp
     REAL(r_2), DIMENSION(:),POINTER    :: FluxNtohwp
     REAL(r_2), DIMENSION(:),POINTER    :: FluxPtohwp
     REAL(r_2), DIMENSION(:),POINTER    :: FluxCtoclear
     REAL(r_2), DIMENSION(:),POINTER    :: FluxNtoclear
     REAL(r_2), DIMENSION(:),POINTER    :: FluxPtoclear
     REAL(r_2), DIMENSION(:),POINTER    :: CtransferLUC


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
     ! added yp wang 5/nov/2012
     REAL(r_2), DIMENSION(:,:), POINTER :: Tairkspin,&
          cgppspin,&
          crmplantspin_1,&
          crmplantspin_2,&
          crmplantspin_3,&
          Tsoilspin_1,&
          Tsoilspin_2,&
          Tsoilspin_3,&
          Tsoilspin_4,&
          Tsoilspin_5,&
          Tsoilspin_6,&
          moistspin_1,&
          moistspin_2,&
          moistspin_3,&
          moistspin_4,&
          moistspin_5,&
          moistspin_6, &
          mtempspin

  END TYPE casa_met

  TYPE casa_balance
     REAL(r_2), DIMENSION(:),POINTER   :: FCgppyear,FCnppyear,                 &
          FCrmleafyear,FCrmwoodyear,FCrmrootyear,FCrgrowyear,               &
          FCrpyear, FCrsyear,FCneeyear,  dCdtyear ,                         &
          LAImax, Cleafmean, Crootmean,                          &
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
     CHARACTER(LEN=99) :: cnpipool=''    ! file for inital pool sizes
     CHARACTER(LEN=99) :: cnpmetin      ! met file for spin up
     CHARACTER(LEN=99) :: cnpmetout     ! met file for spin up
     CHARACTER(LEN=99) :: ndep          ! N deposition input file
     ! added yp wang
     CHARACTER(LEN=99) :: cnpspin       ! input file for spin up
     CHARACTER(LEN=99) :: dump_cnpspin  ! name of dump file for spinning casa-cnp

     CHARACTER(LEN=99) :: phen        ! leaf phenology datafile
     CHARACTER(LEN=99) :: cnpflux     ! modelled mean yearly CNP fluxes
     LOGICAL           :: l_ndep
     ! added vh
     CHARACTER(LEN=99) :: c2cdumppath='' ! cable2casa dump for casa spinup
  END TYPE casafiles_type
  TYPE(casafiles_type) :: casafile

CONTAINS

  SUBROUTINE alloc_casavariable(casabiome,casapool,casaflux, &
       casamet,casabal,arraysize)

    USE casaparm, ONLY : leaf
    IMPLICIT NONE
    TYPE (casa_biome)  , INTENT(INOUT) :: casabiome
    TYPE (casa_pool)   , INTENT(INOUT) :: casapool
    TYPE (casa_flux)   , INTENT(INOUT) :: casaflux
    TYPE (casa_met)    , INTENT(INOUT) :: casamet
    TYPE (casa_balance), INTENT(INOUT) :: casabal
    INTEGER,             INTENT(IN) :: arraysize

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
         casabiome%ratioPcplantmin(mvtype,mplant)    &
         )

    ALLOCATE(casapool%Clabile(arraysize),               &
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
         casapool%ratioPcsoil(arraysize,msoil),     &
         casapool%ratioPcplant(arraysize,mplant),   &
         casapool%ratioPclitter(arraysize,mlitter), &
         casapool%Ctot_0(arraysize),                &
         casapool%Ctot(arraysize)   )

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
         casaflux%fromStoCO2(arraysize,msoil),         &
         casaflux%stemnpp(arraysize),                  &
         casaflux%frac_sapwood(arraysize),             &
         casaflux%sapwood_area(arraysize), &
         casaflux%Cplant_turnover(arraysize,mplant) , &
         casaflux%Cplant_turnover_disturbance(arraysize) , &
         casaflux%Cplant_turnover_crowding(arraysize) , &
         casaflux%Cplant_turnover_resource_limitation(arraysize))

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
         casabal%Crootmean(arraysize)            )


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

  SUBROUTINE alloc_sum_casavariable(  sum_casapool, sum_casaflux &
       ,arraysize)

    USE casaparm, ONLY : leaf
    IMPLICIT NONE
    TYPE (casa_pool)   , INTENT(INOUT) :: sum_casapool
    TYPE (casa_flux)   , INTENT(INOUT) :: sum_casaflux
    INTEGER,             INTENT(IN) :: arraysize


    ALLOCATE(sum_casapool%Clabile(arraysize),               &
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
         sum_casapool%ratioPcsoil(arraysize,msoil),     &
         sum_casapool%ratioPcplant(arraysize,mplant),   &
         sum_casapool%ratioPclitter(arraysize,mlitter)  &
         )

    ALLOCATE(sum_casaflux%Cgpp(arraysize),                     &
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

    ALLOCATE(sum_casaflux%FluxCtoco2(arraysize))

  END SUBROUTINE alloc_sum_casavariable

  SUBROUTINE zero_sum_casa(sum_casapool, sum_casaflux)

    IMPLICIT NONE
    TYPE (casa_pool)   , INTENT(INOUT) :: sum_casapool
    TYPE (casa_flux)   , INTENT(INOUT) :: sum_casaflux

    sum_casapool%Clabile = 0
    sum_casapool%dClabiledt = 0
    sum_casapool%Cplant = 0
    sum_casapool%Nplant = 0
    sum_casapool%Pplant = 0
    sum_casapool%dCplantdt = 0
    sum_casapool%dNplantdt = 0
    sum_casapool%dPplantdt = 0
    sum_casapool%ratioNCplant = 0
    sum_casapool%ratioNPplant = 0
    sum_casapool%Nsoilmin = 0
    sum_casapool%Psoillab = 0
    sum_casapool%Psoilsorb  = 0
    sum_casapool%Psoilocc = 0
    sum_casapool%dNsoilmindt = 0
    sum_casapool%dPsoillabdt = 0
    sum_casapool%dPsoilsorbdt = 0
    sum_casapool%dPsoiloccdt = 0
    sum_casapool%Clitter = 0
    sum_casapool%Nlitter = 0
    sum_casapool%Plitter = 0
    sum_casapool%dClitterdt = 0
    sum_casapool%dNlitterdt = 0
    sum_casapool%dPlitterdt = 0
    sum_casapool%ratioNClitter = 0
    sum_casapool%ratioNPlitter = 0
    sum_casapool%Csoil = 0
    sum_casapool%Nsoil = 0
    sum_casapool%Psoil = 0
    sum_casapool%dCsoildt = 0
    sum_casapool%dNsoildt = 0
    sum_casapool%dPsoildt = 0
    sum_casapool%ratioNCsoil = 0
    sum_casapool%ratioNPsoil = 0
    sum_casapool%ratioNCsoilnew = 0
    sum_casapool%ratioNCsoilmin = 0
    sum_casapool%ratioNCsoilmax = 0
    sum_casapool%ratioPcsoil = 0
    sum_casapool%ratioPcplant = 0
    sum_casapool%ratioPclitter = 0

    sum_casaflux%Cgpp = 0
    sum_casaflux%Cnpp = 0
    sum_casaflux%Crp = 0
    sum_casaflux%Crgplant = 0
    sum_casaflux%Nminfix = 0
    sum_casaflux%Nminuptake = 0
    sum_casaflux%Plabuptake = 0
    sum_casaflux%Clabloss = 0
    sum_casaflux%fracClabile = 0
    sum_casaflux%fracCalloc = 0
    sum_casaflux%fracNalloc = 0
    sum_casaflux%fracPalloc = 0
    sum_casaflux%kplant = 0
    sum_casaflux%Crmplant = 0
    sum_casaflux%fromPtoL = 0
    sum_casaflux%Cnep = 0
    sum_casaflux%Crsoil = 0
    sum_casaflux%Nmindep = 0
    sum_casaflux%Nminloss = 0
    sum_casaflux%Nminleach = 0
    sum_casaflux%Nupland = 0
    sum_casaflux%Nlittermin = 0
    sum_casaflux%Nsmin = 0
    sum_casaflux%Nsimm = 0
    sum_casaflux%Nsnet = 0
    sum_casaflux%fNminloss = 0
    sum_casaflux%fNminleach = 0
    sum_casaflux%Pdep = 0
    sum_casaflux%Pwea = 0
    sum_casaflux%Pleach = 0
    sum_casaflux%Ploss = 0
    sum_casaflux%Pupland = 0
    sum_casaflux%Plittermin = 0
    sum_casaflux%Psmin = 0
    sum_casaflux%Psimm = 0
    sum_casaflux%Psnet = 0
    sum_casaflux%fPleach = 0
    sum_casaflux%kplab = 0
    sum_casaflux%kpsorb = 0
    sum_casaflux%kpocc = 0
    sum_casaflux%kmlabP = 0
    sum_casaflux%Psorbmax = 0
    sum_casaflux%klitter = 0
    sum_casaflux%ksoil = 0
    sum_casaflux%fromLtoS = 0
    sum_casaflux%fromStoS = 0
    sum_casaflux%fromLtoCO2 = 0
    sum_casaflux%fromStoCO2 = 0
    sum_casaflux%stemnpp = 0
    sum_casaflux%frac_sapwood = 0
    sum_casaflux%sapwood_area = 0
    sum_casaflux%Cplant_turnover = 0
    sum_casaflux%Cplant_turnover_disturbance = 0
    sum_casaflux% Cplant_turnover_crowding = 0
    sum_casaflux%Cplant_turnover_resource_limitation = 0

    sum_casaflux%FluxCtolitter = 0
    sum_casaflux%FluxNtolitter = 0
    sum_casaflux%FluxPtolitter = 0

    sum_casaflux%FluxCtosoil = 0
    sum_casaflux%FluxNtosoil = 0
    sum_casaflux%FluxPtosoil = 0

    sum_casaflux%FluxCtoco2 = 0




  END SUBROUTINE zero_sum_casa


  SUBROUTINE update_sum_casa(sum_casapool, sum_casaflux, casapool, casaflux, &
       sum_now, average_now, nsteps)

    IMPLICIT NONE
    TYPE (casa_pool)   , INTENT(INOUT) :: sum_casapool
    TYPE (casa_flux)   , INTENT(INOUT) :: sum_casaflux
    TYPE (casa_pool)   , INTENT(IN) :: casapool
    TYPE (casa_flux)   , INTENT(IN) :: casaflux
    LOGICAL, INTENT(IN) :: sum_now, average_now
    INTEGER, INTENT(IN) :: nsteps

    IF (sum_now) THEN

       sum_casapool%Clabile = sum_casapool%Clabile + casapool%Clabile
       sum_casapool%dClabiledt = sum_casapool%Clabile + casapool%Clabile
       sum_casapool%Cplant =sum_casapool%Cplant + casapool%Cplant
       sum_casapool%Nplant =  sum_casapool%Nplant + casapool%Nplant
       sum_casapool%Pplant =  sum_casapool%Pplant + casapool%Pplant
       sum_casapool%dCplantdt = sum_casapool%dCplantdt + casapool%dCplantdt
       sum_casapool%dNplantdt = sum_casapool%dNplantdt + casapool%dNplantdt
       sum_casapool%dPplantdt = sum_casapool%dPplantdt + casapool%dPplantdt
       sum_casapool%ratioNCplant = sum_casapool%ratioNCplant + casapool%ratioNCplant
       sum_casapool%ratioNPplant = sum_casapool%ratioNPplant + casapool%ratioNPplant
       sum_casapool%Nsoilmin =  sum_casapool%Nsoilmin + casapool%Nsoilmin
       sum_casapool%Psoillab = sum_casapool%Psoillab + casapool%Psoillab
       sum_casapool%Psoilsorb  = sum_casapool%Psoilsorb + casapool%Psoilsorb
       sum_casapool%Psoilocc = sum_casapool%Psoilocc + casapool%Psoilocc
       sum_casapool%dNsoilmindt =  sum_casapool%dNsoilmindt + casapool%dNsoilmindt
       sum_casapool%dPsoillabdt = sum_casapool%dPsoillabdt + casapool%dPsoillabdt
       sum_casapool%dPsoilsorbdt =sum_casapool%dPsoilsorbdt + casapool%dPsoilsorbdt
       sum_casapool%dPsoiloccdt = sum_casapool%dPsoiloccdt +casapool%dPsoiloccdt
       sum_casapool%Clitter = sum_casapool%Clitter + casapool%Clitter
       sum_casapool%Nlitter = sum_casapool%Nlitter  + casapool%Nlitter
       sum_casapool%Plitter =  sum_casapool%Plitter + casapool%Plitter
       sum_casapool%dClitterdt = sum_casapool%dClitterdt  + casapool%dClitterdt
       sum_casapool%dNlitterdt = sum_casapool%dNlitterdt  + casapool%dNlitterdt
       sum_casapool%dPlitterdt = sum_casapool%dPlitterdt + casapool%dPlitterdt
       sum_casapool%ratioNClitter = sum_casapool%ratioNClitter + casapool%ratioNClitter
       sum_casapool%ratioNPlitter = sum_casapool%ratioNPlitter + casapool%ratioNPlitter
       sum_casapool%Csoil =  sum_casapool%Csoil + casapool%Csoil
       sum_casapool%Nsoil =  sum_casapool%Nsoil + casapool%Nsoil
       sum_casapool%Psoil =  sum_casapool%Psoil + casapool%Psoil
       sum_casapool%dCsoildt = sum_casapool%dCsoildt + casapool%dCsoildt
       sum_casapool%dNsoildt = sum_casapool%dNsoildt + casapool%dNsoildt
       sum_casapool%dPsoildt = sum_casapool%dPsoildt + casapool%dPsoildt
       sum_casapool%ratioNCsoil = sum_casapool%ratioNCsoil + casapool%ratioNCsoil
       sum_casapool%ratioNPsoil = sum_casapool%ratioNPsoil + casapool%ratioNPsoil
       sum_casapool%ratioNCsoilnew = sum_casapool%ratioNCsoilnew + casapool%ratioNCsoilnew
       sum_casapool%ratioNCsoilmin =  sum_casapool%ratioNCsoilmin + casapool%ratioNCsoilmin
       sum_casapool%ratioNCsoilmax = sum_casapool%ratioNCsoilmax + casapool%ratioNCsoilmax
       sum_casapool%ratioPcsoil =  sum_casapool%ratioPcsoil  + casapool%ratioPcsoil
       sum_casapool%ratioPcplant =  sum_casapool%ratioPcplant + casapool%ratioPcplant
       sum_casapool%ratioPclitter =   sum_casapool%ratioPclitter + casapool%ratioPclitter

       sum_casaflux%Cgpp = sum_casaflux%Cgpp  + casaflux%Cgpp
       sum_casaflux%Cnpp = sum_casaflux%Cnpp + casaflux%Cnpp
       sum_casaflux%Crp = sum_casaflux%Crp  + casaflux%Crp
       sum_casaflux%Crgplant = sum_casaflux%Crgplant + casaflux%Crgplant
       sum_casaflux%Nminfix =  sum_casaflux%Nminfix + casaflux%Nminfix
       sum_casaflux%Nminuptake =  sum_casaflux%Nminuptake + casaflux%Nminuptake
       sum_casaflux%Plabuptake =  sum_casaflux%Plabuptake + casaflux%Plabuptake
       sum_casaflux%Clabloss =  sum_casaflux%Clabloss + casaflux%Clabloss
       sum_casaflux%fracClabile = sum_casaflux%fracClabile +  casaflux%fracClabile
       sum_casaflux%fracCalloc =  sum_casaflux%fracCalloc + casaflux%fracCalloc*casapool%cplant
       sum_casaflux%fracNalloc = sum_casaflux%fracNalloc + casaflux%fracNalloc
       sum_casaflux%fracPalloc =  sum_casaflux%fracPalloc + casaflux%fracPalloc
       sum_casaflux%kplant =  sum_casaflux%kplant + casaflux%kplant*casapool%cplant
       sum_casaflux%Crmplant =   sum_casaflux%Crmplant + casaflux%Crmplant
       sum_casaflux%fromPtoL =  sum_casaflux%fromPtoL + casaflux%fromPtoL
       sum_casaflux%Cnep =  sum_casaflux%Cnep +  casaflux%Cnep
       sum_casaflux%Crsoil = sum_casaflux%Crsoil + casaflux%Crsoil
       sum_casaflux%Nmindep = sum_casaflux%Nmindep + casaflux%Nmindep
       sum_casaflux%Nminloss = sum_casaflux%Nminloss + casaflux%Nminloss
       sum_casaflux%Nminleach =  sum_casaflux%Nminleach + casaflux%Nminleach
       sum_casaflux%Nupland = sum_casaflux%Nupland + casaflux%Nupland
       sum_casaflux%Nlittermin =  sum_casaflux%Nlittermin +  casaflux%Nlittermin
       sum_casaflux%Nsmin =  sum_casaflux%Nsmin +  casaflux%Nsmin
       sum_casaflux%Nsimm =  sum_casaflux%Nsimm + casaflux%Nsimm
       sum_casaflux%Nsnet = sum_casaflux%Nsnet + casaflux%Nsnet
       sum_casaflux%fNminloss =  sum_casaflux%fNminloss + casaflux%fNminloss
       sum_casaflux%fNminleach =  sum_casaflux%fNminleach + casaflux%fNminleach
       sum_casaflux%Pdep = sum_casaflux%Pdep +  casaflux%Pdep
       sum_casaflux%Pwea = sum_casaflux%Pwea + casaflux%Pwea
       sum_casaflux%Pleach =  sum_casaflux%Pleach + casaflux%Pleach
       sum_casaflux%Ploss =  sum_casaflux%Ploss +  casaflux%Ploss
       sum_casaflux%Pupland =  sum_casaflux%Pupland  + casaflux%Pupland
       sum_casaflux%Plittermin =  sum_casaflux%Plittermin + casaflux%Plittermin
       sum_casaflux%Psmin = sum_casaflux%Psmin +  casaflux%Psmin
       sum_casaflux%Psimm = sum_casaflux%Psimm + casaflux%Psimm
       sum_casaflux%Psnet =  sum_casaflux%Psnet + casaflux%Psnet
       sum_casaflux%fPleach =  sum_casaflux%fPleach + casaflux%fPleach
       sum_casaflux%kplab =  sum_casaflux%kplab + casaflux%kplab
       sum_casaflux%kpsorb = sum_casaflux%kpsorb + casaflux%kpsorb
       sum_casaflux%kpocc =  sum_casaflux%kpocc + casaflux%kpocc
       sum_casaflux%kmlabP =  sum_casaflux%kmlabP + casaflux%kmlabP
       sum_casaflux%Psorbmax =  sum_casaflux%Psorbmax + casaflux%Psorbmax
       sum_casaflux%klitter =  sum_casaflux%klitter + casaflux%klitter
       sum_casaflux%ksoil =  sum_casaflux%ksoil + casaflux%ksoil
       sum_casaflux%fromLtoS =  sum_casaflux%fromLtoS + casaflux%fromLtoS
       sum_casaflux%fromStoS =  sum_casaflux%fromStoS + casaflux%fromStoS
       sum_casaflux%fromLtoCO2 =  sum_casaflux%fromLtoCO2 + casaflux%fromLtoCO2
       sum_casaflux%fromStoCO2 = sum_casaflux%fromStoCO2 + casaflux%fromStoCO2
       sum_casaflux%stemnpp =  sum_casaflux%stemnpp + casaflux%stemnpp
       sum_casaflux%frac_sapwood = sum_casaflux%frac_sapwood + casaflux%frac_sapwood
       sum_casaflux%sapwood_area = sum_casaflux%sapwood_area + casaflux%sapwood_area
       sum_casaflux%Cplant_turnover = sum_casaflux%Cplant_turnover + casaflux%Cplant_turnover
       sum_casaflux%Cplant_turnover_disturbance = sum_casaflux%Cplant_turnover_disturbance + &
            casaflux%Cplant_turnover_disturbance
       sum_casaflux%Cplant_turnover_crowding = sum_casaflux%Cplant_turnover_crowding + &
            casaflux%Cplant_turnover_crowding
       sum_casaflux%Cplant_turnover_resource_limitation =  sum_casaflux%Cplant_turnover_resource_limitation +  &
            casaflux%Cplant_turnover_resource_limitation


       sum_casaflux%FluxCtolitter = sum_casaflux%FluxCtolitter + casaflux%FluxCtolitter
       sum_casaflux%FluxNtolitter =  sum_casaflux%FluxNtolitter + casaflux%FluxNtolitter
       sum_casaflux%FluxPtolitter = sum_casaflux%FluxPtolitter + casaflux%FluxPtolitter

       sum_casaflux%FluxCtosoil = sum_casaflux%FluxCtosoil + casaflux%FluxCtosoil
       sum_casaflux%FluxNtosoil =  sum_casaflux%FluxNtosoil + casaflux%FluxNtosoil
       sum_casaflux%FluxPtosoil = sum_casaflux%FluxPtosoil  + casaflux%FluxPtosoil

       sum_casaflux%FluxCtoco2 =  sum_casaflux%FluxCtoco2 + casaflux%FluxCtoco2
    ENDIF

    IF (average_now) THEN
       sum_casapool%Clabile = sum_casapool%Clabile/REAL(nsteps)
       sum_casapool%dClabiledt = sum_casapool%Clabile/REAL(nsteps)
       WHERE (sum_casapool%Cplant.GT.1.e-12)
          sum_casaflux%fracCalloc =  sum_casaflux%fracCalloc/sum_casapool%Cplant
       ELSEWHERE
          sum_casaflux%fracCalloc = 0.0
       endwhere

       WHERE (sum_casapool%Cplant.GT.1.e-12)
          sum_casaflux%kplant =  sum_casaflux%kplant/sum_casapool%Cplant
       ELSEWHERE
          sum_casaflux%kplant = 0.0
       endwhere
       sum_casapool%Cplant =sum_casapool%Cplant/REAL(nsteps)
       sum_casapool%Nplant =  sum_casapool%Nplant/REAL(nsteps)
       sum_casapool%Pplant =  sum_casapool%Pplant/REAL(nsteps)
       sum_casapool%dCplantdt = sum_casapool%dCplantdt/REAL(nsteps)
       sum_casapool%dNplantdt = sum_casapool%dNplantdt/REAL(nsteps)
       sum_casapool%dPplantdt = sum_casapool%dPplantdt/REAL(nsteps)
       sum_casapool%ratioNCplant = sum_casapool%ratioNCplant/REAL(nsteps)
       sum_casapool%ratioNPplant = sum_casapool%ratioNPplant/REAL(nsteps)
       sum_casapool%Nsoilmin =  sum_casapool%Nsoilmin/REAL(nsteps)
       sum_casapool%Psoillab = sum_casapool%Psoillab/REAL(nsteps)
       sum_casapool%Psoilsorb  = sum_casapool%Psoilsorb/REAL(nsteps)
       sum_casapool%Psoilocc = sum_casapool%Psoilocc/REAL(nsteps)
       sum_casapool%dNsoilmindt =  sum_casapool%dNsoilmindt/REAL(nsteps)
       sum_casapool%dPsoillabdt = sum_casapool%dPsoillabdt/REAL(nsteps)
       sum_casapool%dPsoilsorbdt =sum_casapool%dPsoilsorbdt/REAL(nsteps)
       sum_casapool%dPsoiloccdt = sum_casapool%dPsoiloccdt/REAL(nsteps)
       sum_casapool%Clitter = sum_casapool%Clitter/REAL(nsteps)
       sum_casapool%Nlitter = sum_casapool%Nlitter /REAL(nsteps)
       sum_casapool%Plitter =  sum_casapool%Plitter/REAL(nsteps)
       sum_casapool%dClitterdt = sum_casapool%dClitterdt /REAL(nsteps)
       sum_casapool%dNlitterdt = sum_casapool%dNlitterdt /REAL(nsteps)
       sum_casapool%dPlitterdt = sum_casapool%dPlitterdt/REAL(nsteps)
       sum_casapool%ratioNClitter = sum_casapool%ratioNClitter/REAL(nsteps)
       sum_casapool%ratioNPlitter = sum_casapool%ratioNPlitter/REAL(nsteps)
       sum_casapool%Csoil =  sum_casapool%Csoil/REAL(nsteps)
       sum_casapool%Nsoil =  sum_casapool%Nsoil/REAL(nsteps)
       sum_casapool%Psoil =  sum_casapool%Psoil/REAL(nsteps)
       sum_casapool%dCsoildt = sum_casapool%dCsoildt/REAL(nsteps)
       sum_casapool%dNsoildt = sum_casapool%dNsoildt/REAL(nsteps)
       sum_casapool%dPsoildt = sum_casapool%dPsoildt/REAL(nsteps)
       sum_casapool%ratioNCsoil = sum_casapool%ratioNCsoil/REAL(nsteps)
       sum_casapool%ratioNPsoil = sum_casapool%ratioNPsoil/REAL(nsteps)
       sum_casapool%ratioNCsoilnew = sum_casapool%ratioNCsoilnew/REAL(nsteps)
       sum_casapool%ratioNCsoilmin =  sum_casapool%ratioNCsoilmin/REAL(nsteps)
       sum_casapool%ratioNCsoilmax = sum_casapool%ratioNCsoilmax/REAL(nsteps)
       sum_casapool%ratioPcsoil =  sum_casapool%ratioPcsoil /REAL(nsteps)
       sum_casapool%ratioPcplant =  sum_casapool%ratioPcplant/REAL(nsteps)
       sum_casapool%ratioPclitter =   sum_casapool%ratioPclitter/REAL(nsteps)

       sum_casaflux%Cgpp = sum_casaflux%Cgpp /REAL(nsteps)
       sum_casaflux%Cnpp = sum_casaflux%Cnpp/REAL(nsteps)
       sum_casaflux%Crp = sum_casaflux%Crp /REAL(nsteps)
       sum_casaflux%Crgplant = sum_casaflux%Crgplant/REAL(nsteps)
       sum_casaflux%Nminfix =  sum_casaflux%Nminfix/REAL(nsteps)
       sum_casaflux%Nminuptake =  sum_casaflux%Nminuptake/REAL(nsteps)
       sum_casaflux%Plabuptake =  sum_casaflux%Plabuptake/REAL(nsteps)
       sum_casaflux%Clabloss =  sum_casaflux%Clabloss/REAL(nsteps)
       sum_casaflux%fracClabile = sum_casaflux%fracClabile/REAL(nsteps)
       !  sum_casaflux%fracCalloc =  sum_casaflux%fracCalloc/real(nsteps)
       sum_casaflux%fracNalloc = sum_casaflux%fracNalloc/REAL(nsteps)
       sum_casaflux%fracPalloc =  sum_casaflux%fracPalloc/REAL(nsteps)
       ! sum_casaflux%kplant =  sum_casaflux%kplant/real(nsteps)


       sum_casaflux%Crmplant =   sum_casaflux%Crmplant/REAL(nsteps)
       sum_casaflux%fromPtoL =  sum_casaflux%fromPtoL/REAL(nsteps)
       sum_casaflux%Cnep =  sum_casaflux%Cnep/REAL(nsteps)
       sum_casaflux%Crsoil = sum_casaflux%Crsoil/REAL(nsteps)
       sum_casaflux%Nmindep = sum_casaflux%Nmindep/REAL(nsteps)
       sum_casaflux%Nminloss = sum_casaflux%Nminloss/REAL(nsteps)
       sum_casaflux%Nminleach =  sum_casaflux%Nminleach/REAL(nsteps)
       sum_casaflux%Nupland = sum_casaflux%Nupland/REAL(nsteps)
       sum_casaflux%Nlittermin =  sum_casaflux%Nlittermin/REAL(nsteps)
       sum_casaflux%Nsmin =  sum_casaflux%Nsmin/REAL(nsteps)
       sum_casaflux%Nsimm =  sum_casaflux%Nsimm/REAL(nsteps)
       sum_casaflux%Nsnet = sum_casaflux%Nsnet/REAL(nsteps)
       sum_casaflux%fNminloss =  sum_casaflux%fNminloss/REAL(nsteps)
       sum_casaflux%fNminleach =  sum_casaflux%fNminleach/REAL(nsteps)
       sum_casaflux%Pdep = sum_casaflux%Pdep/REAL(nsteps)
       sum_casaflux%Pwea = sum_casaflux%Pwea/REAL(nsteps)
       sum_casaflux%Pleach =  sum_casaflux%Pleach/REAL(nsteps)
       sum_casaflux%Ploss =  sum_casaflux%Ploss/REAL(nsteps)
       sum_casaflux%Pupland =  sum_casaflux%Pupland /REAL(nsteps)
       sum_casaflux%Plittermin =  sum_casaflux%Plittermin/REAL(nsteps)
       sum_casaflux%Psmin = sum_casaflux%Psmin/REAL(nsteps)
       sum_casaflux%Psimm = sum_casaflux%Psimm/REAL(nsteps)
       sum_casaflux%Psnet =  sum_casaflux%Psnet/REAL(nsteps)
       sum_casaflux%fPleach =  sum_casaflux%fPleach/REAL(nsteps)
       sum_casaflux%kplab =  sum_casaflux%kplab/REAL(nsteps)
       sum_casaflux%kpsorb = sum_casaflux%kpsorb/REAL(nsteps)
       sum_casaflux%kpocc =  sum_casaflux%kpocc/REAL(nsteps)
       sum_casaflux%kmlabP =  sum_casaflux%kmlabP/REAL(nsteps)
       sum_casaflux%Psorbmax =  sum_casaflux%Psorbmax/REAL(nsteps)
       sum_casaflux%klitter =  sum_casaflux%klitter/REAL(nsteps)
       sum_casaflux%ksoil =  sum_casaflux%ksoil/REAL(nsteps)
       sum_casaflux%fromLtoS =  sum_casaflux%fromLtoS/REAL(nsteps)
       sum_casaflux%fromStoS =  sum_casaflux%fromStoS/REAL(nsteps)
       sum_casaflux%fromLtoCO2 =  sum_casaflux%fromLtoCO2/REAL(nsteps)
       sum_casaflux%fromStoCO2 = sum_casaflux%fromStoCO2/REAL(nsteps)
       sum_casaflux%stemnpp =  sum_casaflux%stemnpp/REAL(nsteps)
       sum_casaflux%frac_sapwood = sum_casaflux%frac_sapwood/REAL(nsteps)
       sum_casaflux%sapwood_area = sum_casaflux%sapwood_area/REAL(nsteps)
       sum_casaflux%Cplant_turnover = sum_casaflux%Cplant_turnover/REAL(nsteps)
       sum_casaflux%Cplant_turnover_disturbance = casaflux%Cplant_turnover_disturbance/REAL(nsteps)
       sum_casaflux% Cplant_turnover_crowding = sum_casaflux%Cplant_turnover_crowding/REAL(nsteps)
       sum_casaflux%Cplant_turnover_resource_limitation = &
            sum_casaflux%Cplant_turnover_resource_limitation/REAL(nsteps)
       sum_casaflux%FluxCtolitter = sum_casaflux%FluxCtolitter/REAL(nsteps)
       sum_casaflux%FluxNtolitter =  sum_casaflux%FluxNtolitter/REAL(nsteps)
       sum_casaflux%FluxPtolitter = sum_casaflux%FluxPtolitter/REAL(nsteps)

       sum_casaflux%FluxCtosoil = sum_casaflux%FluxCtosoil/REAL(nsteps)
       sum_casaflux%FluxNtosoil =  sum_casaflux%FluxNtosoil/REAL(nsteps)
       sum_casaflux%FluxPtosoil = sum_casaflux%FluxPtosoil /REAL(nsteps)

       sum_casaflux%FluxCtoco2 =  sum_casaflux%FluxCtoco2/REAL(nsteps)
    ENDIF



  END SUBROUTINE update_sum_casa


END MODULE casavariable
