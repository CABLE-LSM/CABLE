!==============================================================================
! This source code is part of the 
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CABLE Academic User Licence Agreement 
! (the "Licence").
! You may not use this file except in compliance with the Licence.
! A copy of the Licence and registration form can be obtained from 
! http://www.cawcr.gov.au/projects/access/cable
! You need to register and read the Licence agreement before use.
! Please contact cable_help@nf.nci.org.au for any questions on 
! registration and the Licence.
!
! Unless required by applicable law or agreed to in writing, 
! software distributed under the Licence is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the Licence for the specific language governing permissions and 
! limitations under the Licence.
! ==============================================================================
!
! Purpose: defines/allocates variables for CASA-CNP
!
! Contact: Yingping.Wang@csiro.au
!
! History: Developed for offline CASA-CNP, code revision likely to better
!          suit ACCESS and to merge more consistently with CABLE code
!
!
! ==============================================================================
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
                                       kclabrate

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
    REAL(r_2), DIMENSION(:),POINTER      :: meangpp
    REAL(r_2), DIMENSION(:),POINTER      :: meanrleaf 
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
    REAL(r_2), DIMENSION(:),POINTER   :: FCgppyear,FCnppyear,             &
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

  ALLOCATE(casaflux%meangpp(arraysize))
  ALLOCATE(casaflux%meanrleaf(arraysize))

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

