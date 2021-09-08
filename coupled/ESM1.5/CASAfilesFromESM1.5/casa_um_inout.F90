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
! Purpose: CASA-CNP
!
! Contact: Yingping Wang, Lauren Stevens
!
! History: 2012
!
!
! ==============================================================================
!
! Lauren Stevens Jan-Apr 2011
! Casa CNP for UM
!
! Includes:
! init_casacnp
! casa_readpoint_pk
! casa_init_pk 
! pack_cnppool
! casa_poolout_unpk
! unpack_cnppool

MODULE casa_um_inout_mod

  USE landuse_mod
  USE feedback_mod
  USE casa_inout_mod

IMPLICIT NONE 

CONTAINS 

!========================================================================
!========================================================================
!========================================================================

!SUBROUTINE cable_casa_init(sin_theta_latitude,um1,cpool_tile,npool_tile,ppool_tile, &
SUBROUTINE init_casacnp(sin_theta_latitude,cpool_tile,npool_tile,ppool_tile, &
                           soil_order,nidep,nifix,pwea,pdust,&
                           wood_hvest_c,wood_hvest_n,wood_hvest_P,&
                           wood_flux_c,wood_flux_n,wood_flux_P,&
                           wresp_c,wresp_n,wresp_P,thinning,&
                           GLAI,PHENPHASE,PREV_YR_SFRAC,idoy)
! Lest 20 Jan 2011
    USE cable_def_types_mod
    USE cable_um_tech_mod, ONLY : um1, canopy
    USE cable_params_mod, ONLY : veg => veg_cbl 
    USE cable_params_mod, ONLY : soil => soil_cbl 
    USE casavariable
    USE phenvariable
    USE casa_types_mod
    !USE casa_inout

IMPLICIT NONE

!    TYPE (um_dimensions), INTENT(IN)          :: um1
!   LOGICAL, INTENT(INOUT), DIMENSION(um1%land_pts, um1%ntiles) :: L_tile_pts
    REAL   , INTENT(IN)   , DIMENSION(um1%row_length,um1%rows)  :: sin_theta_latitude
    REAL   , INTENT(INOUT) :: cpool_tile(um1%land_pts,um1%ntiles,10)
    REAL   , INTENT(INOUT) :: npool_tile(um1%land_pts,um1%ntiles,10)
    REAL   , INTENT(INOUT) :: ppool_tile(um1%land_pts,um1%ntiles,12)
    REAL   , INTENT(INOUT) :: soil_order(um1%land_pts)
    REAL   , INTENT(INOUT) :: nidep(um1%land_pts) 
    REAL   , INTENT(INOUT) :: nifix(um1%land_pts)
    REAL   , INTENT(INOUT) :: pwea(um1%land_pts)
    REAL   , INTENT(INOUT) :: pdust(um1%land_pts)
    REAL   , INTENT(INOUT) :: GLAI(um1%land_pts,um1%ntiles)
    REAL   , INTENT(INOUT) :: PHENPHASE(um1%land_pts,um1%ntiles)
    !INTEGER, INTENT(INOUT) :: PHENPHASE(um1%land_pts,um1%ntiles)
    REAL   , INTENT(INOUT) :: PREV_YR_SFRAC(um1%land_pts,um1%ntiles)
    REAL   , INTENT(INOUT) :: WOOD_HVEST_C(um1%land_pts,um1%ntiles,3)
    REAL   , INTENT(INOUT) :: WOOD_HVEST_N(um1%land_pts,um1%ntiles,3)
    REAL   , INTENT(INOUT) :: WOOD_HVEST_P(um1%land_pts,um1%ntiles,3)
    REAL   , INTENT(INOUT) :: WRESP_C(um1%land_pts,um1%ntiles,3)
    REAL   , INTENT(INOUT) :: WRESP_N(um1%land_pts,um1%ntiles,3)
    REAL   , INTENT(INOUT) :: WRESP_P(um1%land_pts,um1%ntiles,3)
    REAL   , INTENT(INOUT) :: WOOD_FLUX_C(um1%land_pts,um1%ntiles)
    REAL   , INTENT(INOUT) :: WOOD_FLUX_N(um1%land_pts,um1%ntiles)
    REAL   , INTENT(INOUT) :: WOOD_FLUX_P(um1%land_pts,um1%ntiles)
    REAL   , INTENT(INOUT) :: THINNING(um1%land_pts,um1%ntiles)
    INTEGER :: idoy

!    casafile%cnpbiome ='/home/599/lxs599/surface_data/pftlookup_csiro_v16_17tiles.csv'
!    casafile%phen     ='/home/599/lxs599/surface_data/modis_phenology_v1.txt'
!    casafile%cnppoint =' '
!    casafile%cnpepool =' '
!    casafile%cnpipool =' '
!    casafile%cnpmetin =' '
!    casafile%cnpmetout=' '
!!   casafile%cnpmetin ='/home/599/lxs599/surface_data/casametDH.csv'
!    casafile%cnpflux  =' '

!print *,'Lest um_inout',mp,um1%land_pts,um1%ntiles

    CALL alloc_casavariable(casabiome,casapool,casaflux,casamet,casabal,mp)
    CALL alloc_phenvariable(phen,mp)

    !CALL casa_readpoint_pk(sin_theta_latitude,veg,soil,casaflux,casamet,um1, &
    CALL casa_readpoint_pk(sin_theta_latitude,veg,soil,casaflux,casamet, &
                          nidep,nifix,pwea,pdust,soil_order)
    CALL casa_readbiome(veg,soil,casabiome,casapool,casaflux,casamet,phen)
    CALL casa_readphen(veg,casamet,phen)
    CALL casa_init_pk(casabiome,casaflux,casamet,casapool,casabal,veg,canopy,phen, &
                     cpool_tile,npool_tile,ppool_tile,wood_hvest_c,wood_hvest_n,wood_hvest_p,&
                     wood_flux_c,wood_flux_n,wood_flux_p,wresp_c,wresp_n,wresp_p,thinning,&
                     GLAI,PHENPHASE,PREV_YR_SFRAC,idoy)
                     !um1,cpool_tile,npool_tile,ppool_tile,GLAI,PHENPHASE)
 return
END SUBROUTINE init_casacnp

!========================================================================
!========================================================================
!========================================================================

SUBROUTINE casa_readpoint_pk(sin_theta_latitude,veg,soil,casaflux,casamet, &
                           nidep,nifix,pwea,pdust,soil_order)
!SUBROUTINE casa_readpoint_pk(sin_theta_latitude,veg,soil,casaflux,casamet,um1, &

    USE cable_def_types_mod
    !USE define_dimensions    ! mp, r_1, r_2, i_d
    !USE define_types
    USE cable_um_tech_mod, ONLY : um1
    USE casavariable
    USE casaparm
    USE cable_um_init_subrs_mod, ONLY : um2cable_rr, um2cable_lp
    USE cable_data_module, ONLY : MATH
!    USE math_constants

IMPLICIT NONE

! passed in
    TYPE (veg_parameter_type)           :: veg
    TYPE (soil_parameter_type)          :: soil
    TYPE (casa_flux), INTENT(INOUT)     :: casaflux
    TYPE (casa_met) , INTENT(INOUT)     :: casamet
!    TYPE (um_dimensions), INTENT(IN)   :: um1
!   LOGICAL, INTENT(INOUT),DIMENSION(um1%land_pts, um1%ntiles) :: L_tile_pts
    REAL   , INTENT(IN), DIMENSION(um1%row_length,um1%rows) :: sin_theta_latitude
    REAL, DIMENSION(um1%land_pts)       :: soil_order
    REAL, DIMENSION(um1%land_pts)       :: nidep,nifix,pwea,pdust
! local variables
    INTEGER :: k,p,i
    !REAL,DIMENSION(mp)             :: annNdep,annNfix,annPwea,annPdust
    REAL(r_2),DIMENSION(mp)             :: annNdep,annNfix,annPwea,annPdust
    !REAL(r_2),DIMENSION(:),ALLOCATABLE :: annNdep,annNfix,annPwea,annPdust
    REAL(r_2)                           :: annNfert,annPfert
    !REAL(r_2),ALLOCATABLE               :: annNfert,annPfert
    LOGICAL                             :: skip =.TRUE.
!   REAL,DIMENSION(um1%land_pts,um1%ntiles)     :: ndep,nfix,pweather,pdus
!   INTEGER, DIMENSION(um1%land_pts,um1%ntiles) :: sorder
    INTEGER, DIMENSION(um1%land_pts)            :: sorder
    !INTEGER, DIMENSION(:),ALLOCATABLE           :: sorder

    !allocate(annNdep(mp))
    !allocate(annNfix(mp))
    !allocate(annPwea(mp))
    !allocate(annPdust(mp))
    !allocate(annNfert)
    !allocate(annPfert)
    !allocate(sorder(um1%land_pts))

! initialise 
    sorder(:)  = 0
    annNdep(:) = 0.0; annNfix(:) = 0.0 
    annPwea(:) = 0.0; annPdust(:)= 0.0
    annPfert = 0.7/365.0
    annNfert = 4.3/365.0

! pack variables
    call um2cable_rr((asin(sin_theta_latitude)/math%pi180) ,casamet%lat)
    !call um2cable_rr(um1%latitude ,casamet%lat)
    call um2cable_rr(um1%longitude,casamet%lon)
! Lest Nov2011 - not correct, but areacell not needed/used for UM ?
! Lest Feb2014 - use UM's A_BOXAREAS
    casamet%areacell = pack(um1%tile_frac,um1%l_tile_pts)
    !casamet%areacell = pack(,um1%l_tile_pts)

    sorder = INT(soil_order)
    call um2cable_ilp(sorder,sorder(1:10),casamet%isorder,soil%isoilm,skip)
    !call um2cable_lp(nidep,nidep(1:10),annNdep ,soil%isoilm, skip)
    call um2cable_lp(nifix,nifix(1:10),annNfix ,soil%isoilm, skip)
    call um2cable_lp(pwea ,pwea(1:10) ,annPwea ,soil%isoilm, skip)
    call um2cable_lp(pdust,pdust(1:10),annPdust,soil%isoilm, skip)

    !casaflux%Nmindep = annNdep/365.0      ! gN/m2/day
    casaflux%Nminfix = annNfix/365.0      ! gN/m2/day
    casaflux%Pdep    = annPdust/365.0     ! gP/m2/day
    casaflux%Pwea    = annPwea/365.0      ! gP/m2/day

    do p = 1,mp
       if (veg%iveg(p)==cropland .or. veg%iveg(p)==croplnd2) then
       ! P fertilizer =13 Mt P globally in 1994
         casaflux%Pdep(p) = casaflux%Pdep(p)+annPfert
       ! N fertilizer =86 Mt N globally in 1994
         casaflux%Nminfix(p) = casaflux%Nminfix(p)+annNfert
       endif
    enddo

    !deallocate(annNdep,annNfix,annPwea,annPdust,annNfert,annPfert,sorder)

END SUBROUTINE casa_readpoint_pk

!========================================================================
!========================================================================
!========================================================================

SUBROUTINE casa_init_pk(casabiome,casaflux,casamet,casapool,casabal,veg,canopy,phen, &
                       cpool_tile,npool_tile,ppool_tile,wood_hvest_c,wood_hvest_n,wood_hvest_p,&
                       wood_flux_c,wood_flux_n,wood_flux_p,wresp_c,wresp_n,wresp_p,thinning,&
                       GLAI,PHENPHASE,PREV_YR_SFRAC,idoy)
!                       um1,cpool_tile,npool_tile,ppool_tile,GLAI,PHENPHASE)
!  initialize some values in phenology parameters and leaf growth phase

   USE cable_def_types_mod
!!  USE define_dimensions    ! mp, r_2
  USE casadimension        ! icycle,mplant,mlitter,msoil
!  USE define_types
  USE cable_um_tech_mod, ONLY : um1!, veg
  USE cable_common_module, ONLY : ktau_gl, l_luc
  USE casaparm             !, ONLY : initcasa
  USE casavariable
  USE phenvariable

IMPLICIT NONE

  TYPE (casa_biome),   INTENT(IN)       :: casabiome
  TYPE (casa_flux),    INTENT(INOUT)    :: casaflux
  TYPE (casa_met),     INTENT(INOUT)    :: casamet
  TYPE (casa_pool),    INTENT(INOUT)    :: casapool
  TYPE (casa_balance), INTENT(INOUT)    :: casabal
  TYPE (canopy_type), INTENT(INOUT)     :: canopy
  TYPE (veg_parameter_type), INTENT(INOUT) :: veg
  TYPE (phen_variable),   INTENT(INOUT) :: phen
!  TYPE (um_dimensions), INTENT(IN)      :: um1
! LOGICAL, INTENT(INOUT),DIMENSION(um1%land_pts, um1%ntiles) :: L_tile_pts
  REAL   , INTENT(INOUT) :: cpool_tile(um1%land_pts,um1%ntiles,10)
  REAL   , INTENT(INOUT) :: npool_tile(um1%land_pts,um1%ntiles,10)
  REAL   , INTENT(INOUT) :: ppool_tile(um1%land_pts,um1%ntiles,12)
  REAL   , INTENT(INOUT) :: GLAI(um1%land_pts,um1%ntiles)
  REAL   , INTENT(INOUT) ::PHENPHASE(um1%land_pts,um1%ntiles)
  !INTEGER, INTENT(INOUT) :: PHENPHASE(um1%land_pts,um1%ntiles)
  REAL   , INTENT(INOUT) ::PREV_YR_SFRAC(um1%land_pts,um1%ntiles)
  REAL   , INTENT(INOUT) ::WOOD_HVEST_C(um1%land_pts,um1%ntiles,3)
  REAL   , INTENT(INOUT) ::WOOD_HVEST_N(um1%land_pts,um1%ntiles,3)
  REAL   , INTENT(INOUT) ::WOOD_HVEST_P(um1%land_pts,um1%ntiles,3)
  REAL   , INTENT(INOUT) ::WRESP_C(um1%land_pts,um1%ntiles,3)
  REAL   , INTENT(INOUT) ::WRESP_N(um1%land_pts,um1%ntiles,3)
  REAL   , INTENT(INOUT) ::WRESP_P(um1%land_pts,um1%ntiles,3)
  REAL   , INTENT(INOUT) ::WOOD_FLUX_C(um1%land_pts,um1%ntiles)
  REAL   , INTENT(INOUT) ::WOOD_FLUX_N(um1%land_pts,um1%ntiles)
  REAL   , INTENT(INOUT) ::WOOD_FLUX_P(um1%land_pts,um1%ntiles)
  REAL   , INTENT(INOUT) ::THINNING(um1%land_pts,um1%ntiles)
  INTEGER :: idoy,mtau

  print *, 'LUC1', l_luc
  !l_luc = .TRUE.
  print *, 'LUC2', l_luc

!  PRINT *, 'initcasa,mp ', initcasa,mp

! Lest 6dec11 - moved from bgcdriver
  casamet%tairk(:)   = 0.0
  casamet%tsoil(:,:) = 0.0
  casamet%moist(:,:) = 0.0
  casaflux%cgpp(:)   = 0.0
  casaflux%cnpp(:)   = canopy%fnpp*86400.     !0.0 TZ initilize npp from prognostic
  casaflux%Crsoil(:) = canopy%frs*86400.     !0.0
  casaflux%crgplant(:)   = canopy%frp*86400. !0.0 TZ this might not be corrext crgplant ne frp 
  casaflux%crmplant(:,:) = 0.0
  casaflux%clabloss(:)   = 0.0
!print *,'Lest - Crsoil',casaflux%Crsoil


  ! Lest 19/2/14 - will work for coupled and amip
  ! mtau is the step number of the day (1,2,..47,0)
  mtau = mod(ktau_gl,int(24.*3600./um1%timestep))
  print *, 'Lest mtau', ktau_gl, mtau, idoy, um1%timestep, l_luc

  IF (initcasa==1) THEN
     IF (l_luc .and. idoy == 1 .and. mtau == 1) THEN

        print *, 'Lest LUC', mtau
        CALL casa_reinit_pk(casabiome,casamet,casapool,casabal,veg,phen, &
                          cpool_tile,npool_tile,ppool_tile,&
                          wood_hvest_c,wood_hvest_n,wood_hvest_p,&
                          wood_flux_c,wood_flux_n,wood_flux_p,&
                          wresp_c,wresp_n,wresp_p,thinning,&
                          GLAI,PHENPHASE,PREV_YR_SFRAC)
        print *, 'TEST_LS reinit DONE'

     ELSE

        print *, 'Lest No LUC', mtau
        ! CALL pack_cnppool(casamet,casapool,casabal,phen,um1,cpool_tile,npool_tile, &
        CALL pack_cnppool(casamet,casapool,casabal,phen,cpool_tile,npool_tile, &
                          ppool_tile,GLAI,PHENPHASE)
     ENDIF

  ENDIF

! reset labile C pool,comment out by Q.Zhang 10/09/2011
! casapool%clabile    = 0.0
! check pool sizes
  casapool%cplant     = MAX(0.0,casapool%cplant)
  casapool%clitter    = MAX(0.0,casapool%clitter)
  casapool%csoil      = MAX(0.0,casapool%csoil)
  casabal%cplantlast  = casapool%cplant
  casabal%clitterlast = casapool%clitter
  casabal%csoillast   = casapool%csoil
  casabal%clabilelast = casapool%clabile
  casabal%sumcbal     = 0.0
       casabal%FCgppyear=0.0;casabal%FCrpyear=0.0
       casabal%FCnppyear=0;casabal%FCrsyear=0.0;casabal%FCneeyear=0.0

  IF (icycle==1) THEN
    casapool%nplant(:,:) = casapool%cplant(:,:) * casapool%rationcplant(:,:)
    casapool%Nsoil(:,:)  = casapool%ratioNCsoil(:,:) * casapool%Csoil(:,:)
    casapool%Psoil(:,:)  = casapool%ratioPCsoil(:,:) * casapool%Csoil(:,:)
    casapool%Nsoilmin(:) = 2.5
  ENDIF

  IF (icycle >1) THEN
    casapool%nplant      = MAX(1.e-6,casapool%nplant)
    casapool%nlitter     = MAX(1.e-6,casapool%nlitter)
    casapool%nsoil       = MAX(1.e-6,casapool%nsoil)
    casapool%nsoilmin    = MAX(1.e-6,casapool%nsoilmin)
    casabal%nplantlast   = casapool%nplant
    casabal%nlitterlast  = casapool%nlitter
    casabal%nsoillast    = casapool%nsoil
    casabal%nsoilminlast = casapool%nsoilmin
    casabal%sumnbal      = 0.0
       casabal%FNdepyear=0.0;casabal%FNfixyear=0.0;casabal%FNsnetyear=0.0
       casabal%FNupyear=0.0;casabal%FNleachyear=0.0;casabal%FNlossyear=0.0
  ENDIF

  IF (icycle >2) THEN
    casapool%pplant       = MAX(1.0e-7,casapool%pplant)
    casapool%plitter      = MAX(1.0e-7,casapool%plitter)
    casapool%psoil        = MAX(1.0e-7,casapool%psoil)
    casapool%Psoillab     = MAX(1.0e-7,casapool%psoillab)  ! was 2.0, changed according to  YP
    casapool%psoilsorb    = MAX(1.0e-7,casapool%psoilsorb) ! was 10.0, -
    casapool%psoilocc     = MAX(1.0e-7,casapool%psoilocc)  ! was 50.0, -
    casabal%pplantlast    = casapool%pplant
    casabal%plitterlast   = casapool%plitter
    casabal%psoillast     = casapool%psoil
    casabal%psoillablast  = casapool%psoillab
    casabal%psoilsorblast = casapool%psoilsorb
    casabal%psoilocclast  = casapool%psoilocc
    casabal%sumpbal       = 0.0
       casabal%FPweayear=0.0;casabal%FPdustyear=0.0; casabal%FPsnetyear=0.0
       casabal%FPupyear=0.0;casabal%FPleachyear=0.0;casabal%FPlossyear=0.0
  ENDIF

END SUBROUTINE casa_init_pk

!========================================================================
!========================================================================
!========================================================================

SUBROUTINE casa_reinit_pk(casabiome,casamet,casapool,casabal,veg,phen, &
                          cpool_tile,npool_tile,ppool_tile,&
                          woodhvest_c,woodhvest_n,woodhvest_p,&
                          logc,logn,logp,&
                          wresp_c,wresp_n,wresp_p,thinning,&
                          GLAI,PHENPHASE,PREV_YR_SFRAC)

   USE cable_def_types_mod ! combines def_dimensions (mp,r_2) and define_types (mland)  
   USE casadimension
   USE casaparm
   USE casavariable
   USE phenvariable
   USE cable_common_module, ONLY : ktau_gl, l_thinforest
   
   USE cable_um_tech_mod, ONLY : um1

   IMPLICIT NONE

   TYPE(casa_biome),         INTENT(IN)    :: casabiome
   TYPE(casa_met),           INTENT(INOUT) :: casamet
   TYPE(casa_pool),          INTENT(INOUT) :: casapool
   TYPE(casa_balance),       INTENT(INOUT) :: casabal    
   TYPE(veg_parameter_type), INTENT(INOUT)    :: veg
   TYPE(phen_variable),      INTENT(INOUT) :: phen

   REAL,  INTENT(INOUT)  :: cpool_tile(um1%land_pts,um1%ntiles,10)
   REAL,  INTENT(INOUT)  :: npool_tile(um1%land_pts,um1%ntiles,10)
   REAL,  INTENT(INOUT)  :: ppool_tile(um1%land_pts,um1%ntiles,12)
   REAL,  INTENT(INOUT)  :: GLAI(um1%land_pts,um1%ntiles)
   REAL,  INTENT(INOUT)  :: PHENPHASE(um1%land_pts,um1%ntiles)
   REAL,  INTENT(INOUT)  :: PREV_YR_SFRAC(um1%land_pts,um1%ntiles)

   ! local variables
   REAL(r_2) :: clabile_x(um1%land_pts,um1%ntiles)
   REAL(r_2) :: cplant_x(um1%land_pts,um1%ntiles,mplant)
   REAL(r_2) :: clitter_x(um1%land_pts,um1%ntiles,mlitter)
   REAL(r_2) :: csoil_x(um1%land_pts,um1%ntiles,msoil)
   REAL(r_2) :: nplant_x(um1%land_pts,um1%ntiles,mplant)
   REAL(r_2) :: nlitter_x(um1%land_pts,um1%ntiles,mlitter)
   REAL(r_2) :: nsoil_x(um1%land_pts,um1%ntiles,msoil)
   REAL(r_2) :: nsoilmin_x(um1%land_pts,um1%ntiles)
   REAL(r_2) :: pplant_x(um1%land_pts,um1%ntiles,mplant)
   REAL(r_2) :: plitter_x(um1%land_pts,um1%ntiles,mlitter)
   REAL(r_2) :: psoil_x(um1%land_pts,um1%ntiles,msoil)
   REAL(r_2) :: psoillab_x(um1%land_pts,um1%ntiles)
   REAL(r_2) :: psoilsorb_x(um1%land_pts,um1%ntiles)
   REAL(r_2) :: psoilocc_x(um1%land_pts,um1%ntiles)

   REAL(r_2) :: clabile_y(um1%land_pts,um1%ntiles)
   REAL(r_2) :: cplant_y(um1%land_pts,um1%ntiles,mplant)
   REAL(r_2) :: clitter_y(um1%land_pts,um1%ntiles,mlitter)
   REAL(r_2) :: csoil_y(um1%land_pts,um1%ntiles,msoil)
   REAL(r_2) :: nplant_y(um1%land_pts,um1%ntiles,mplant)
   REAL(r_2) :: nlitter_y(um1%land_pts,um1%ntiles,mlitter)
   REAL(r_2) :: nsoil_y(um1%land_pts,um1%ntiles,msoil)
   REAL(r_2) :: nsoilmin_y(um1%land_pts,um1%ntiles)
   REAL(r_2) :: pplant_y(um1%land_pts,um1%ntiles,mplant)
   REAL(r_2) :: plitter_y(um1%land_pts,um1%ntiles,mlitter)
   REAL(r_2) :: psoil_y(um1%land_pts,um1%ntiles,msoil)
   REAL(r_2) :: psoillab_y(um1%land_pts,um1%ntiles)
   REAL(r_2) :: psoilsorb_y(um1%land_pts,um1%ntiles)
   REAL(r_2) :: psoilocc_y(um1%land_pts,um1%ntiles)

   REAL    :: frac_x(um1%land_pts,um1%ntiles), frac_y(um1%land_pts,um1%ntiles)
   LOGICAL :: ifpre_x(um1%land_pts,um1%ntiles), ifpre_y(um1%land_pts,um1%ntiles)
   ! To be recorded wood log variables.
   REAL(r_2) :: logc(um1%land_pts,um1%ntiles), logn(um1%land_pts,um1%ntiles),logp(um1%land_pts,um1%ntiles) ! wood_flux
   REAL(r_2) :: woodhvest_c(um1%land_pts,um1%ntiles,3),woodhvest_n(um1%land_pts,um1%ntiles,3),woodhvest_p(um1%land_pts,um1%ntiles,3)
   REAL(r_2) :: wresp_c(um1%land_pts,um1%ntiles,3),wresp_n(um1%land_pts,um1%ntiles,3),wresp_p(um1%land_pts,um1%ntiles,3)
   REAL(r_2) :: thinning(um1%land_pts,um1%ntiles)
   !REAL(r_2), DIMENSION(3) :: pool_frac, pool_time
   REAL,PARAMETER:: POOL_FRAC(3) =(/0.33, 0.33, 0.34/)
   REAL,PARAMETER:: POOL_TIME(3) =(/1.00, 0.10, 0.01/)
   REAL(r_2) :: cplant_z(um1%land_pts,um1%ntiles,mplant)
   REAL(r_2) :: nplant_z(um1%land_pts,um1%ntiles,mplant)
   REAL(r_2) :: pplant_z(um1%land_pts,um1%ntiles,mplant)

   ! check if all of them are required
   INTEGER   :: g, p, k, y ! np
   REAL(r_2) :: cbal,nbal,pbal

! Initialize temporary variables

  cplant_x = 0.
  nplant_x = 0.
  pplant_x = 0.

  clitter_x = 0.
  nlitter_x = 0.
  plitter_x = 0.

  csoil_x = 0.
  nsoil_x = 0.
  psoil_x = 0.

  clabile_x = 0.
  nsoilmin_x = 0.
  psoillab_x = 0.
  psoilsorb_x = 0.
  psoilocc_x = 0.

  cplant_y = 0.
  nplant_y = 0.
  pplant_y = 0.

  clitter_y = 0.
  nlitter_y = 0.
  plitter_y = 0.

  csoil_y = 0.
  nsoil_y = 0.
  psoil_y = 0.

  clabile_y = 0.
  nsoilmin_y = 0.
  psoillab_y = 0.
  psoilsorb_y = 0.
  psoilocc_y = 0.

  frac_x = 0.
  frac_y = 0.
  ifpre_x = .FALSE.
  ifpre_y = .FALSE.

  logc = 0.
  logn = 0.
  logp = 0.

  cplant_z = 0.
  nplant_z = 0.
  pplant_z = 0.

  wresp_c = 0.
  wresp_n = 0.
  wresp_p = 0.

  thinning = 0.

  ! assign "old" cnp pool values (from last dump file, initilization)
  clabile_x(:,:)   = cpool_tile(:,:,1)
  cplant_x(:,:,1)  = cpool_tile(:,:,2)
  cplant_x(:,:,2)  = cpool_tile(:,:,3)
  cplant_x(:,:,3)  = cpool_tile(:,:,4)
  clitter_x(:,:,1) = cpool_tile(:,:,5)
  clitter_x(:,:,2) = cpool_tile(:,:,6)
  clitter_x(:,:,3) = cpool_tile(:,:,7)
  csoil_x(:,:,1)   = cpool_tile(:,:,8)
  csoil_x(:,:,2)   = cpool_tile(:,:,9)
  csoil_x(:,:,3)   = cpool_tile(:,:,10)

  IF (icycle>1) THEN
     nplant_x(:,:,1)  = npool_tile(:,:,1)
     nplant_x(:,:,2)  = npool_tile(:,:,2)
     nplant_x(:,:,3)  = npool_tile(:,:,3)
     nlitter_x(:,:,1) = npool_tile(:,:,4)
     nlitter_x(:,:,2) = npool_tile(:,:,5)
     nlitter_x(:,:,3) = npool_tile(:,:,6)
     nsoil_x(:,:,1)   = npool_tile(:,:,7)
     nsoil_x(:,:,2)   = npool_tile(:,:,8)
     nsoil_x(:,:,3)   = npool_tile(:,:,9)
     nsoilmin_x(:,:)  = npool_tile(:,:,10)
  END IF

  IF (icycle>2) THEN
     pplant_x(:,:,1)  = ppool_tile(:,:,1)
     pplant_x(:,:,2)  = ppool_tile(:,:,2)
     pplant_x(:,:,3)  = ppool_tile(:,:,3)
     plitter_x(:,:,1) = ppool_tile(:,:,4)
     plitter_x(:,:,2) = ppool_tile(:,:,5)
     plitter_x(:,:,3) = ppool_tile(:,:,6)
     psoil_x(:,:,1)   = ppool_tile(:,:,7)
     psoil_x(:,:,2)   = ppool_tile(:,:,8)
     psoil_x(:,:,3)   = ppool_tile(:,:,9)
     psoillab_x(:,:)  = ppool_tile(:,:,10)
     psoilsorb_x(:,:) = ppool_tile(:,:,11)
     psoilocc_x(:,:)  = ppool_tile(:,:,12)
  END IF

  ! assign fractions (previous) (need to get fractions from previous year)
  !frac_x(:,:) = um1%tile_frac(:,:)
  frac_x(:,:) = PREV_YR_SFRAC(:,:)

  !!!! DEBUG START !!! 
  !if (mp.ge.4) then
  !   write(6,*)"DEBUG frac first 17 tiles  :", frac_x(1,:)
  !   write(6,*)"DEBUG lon of first point   :", casamet%lon(1:4)
  !   write(6,*)"DEBUG lat of first point   :", casamet%lat(1:4)
  !   write(6,*)"DEBUG frac from metcasa    :", casamet%areacell(1:4)
  !   write(6,*)"DEBUG mp           :", mp
  !   write(6,*)"DEBUG mland        :", mland
  !   write(6,*)"DEBUG um1%land_pts :", um1%land_pts
  !   write(6,*)"DEBUG mvtype       :", mvtype
  !   write(6,*)"DEBUG um1%ntiles   :", um1%ntiles
  !end if
  !!!! DEBUG END !!!!!

  ! assign fractions (current)
  frac_y(:,:) = um1%tile_frac(:,:)

  ! set the ifpre_x and ifpre_y values to true where fractions are .ne. 0
  where(frac_x > 0.) ifpre_x = .TRUE.
  where(frac_y > 0.) ifpre_y = .TRUE.

  !Do k = 1,um1%ntiles
  !DO g = 1, um1%land_pts
  !print *,'Lest frac', g, k, frac_y(g,k), frac_x(g,k)
  !END DO
  !END DO

  ! start main loop
  DO g = 1, um1%land_pts
!DO N=1,um1%NTILES
! DO K=1,um1%TILE_PTS(N)
!    g = um1%TILE_INDEX(K,N)

     !write(6,*)"TEST_TZ start of main loop"

     ! Check all glacier tiles
     IF (ifpre_x(g,iceland).and.ifpre_y(g,iceland)) THEN
         IF (abs(frac_x(g,iceland)-1.0)>0.01 .or. abs(frac_y(g,iceland)-1.0)>0.01) THEN
         print *,'Lest Glacier',g,ifpre_x(g,iceland),ifpre_y(g,iceland),&
                                 frac_x(g,iceland), frac_y(g,iceland)
            STOP "Glacier fraction .ne. 1"
         ELSE
            print *, 'Lest cycle', g,ifpre_x(g,iceland),ifpre_y(g,iceland),&
                                 frac_x(g,iceland), frac_y(g,iceland)
            cycle
         END IF
     ELSEIF (.not.ifpre_x(g,iceland) .and. .not.ifpre_y(g,iceland)) THEN
         IF (abs(frac_x(g,iceland)-0.0)>0.01 .or. abs(frac_y(g,iceland)-0.0)>0.01) THEN
         print *, 'Lest landpt fracs x', frac_x(g,:)
         print *, 'Lest landpt fracs y', frac_y(g,:)
         print *,'Lest Glacier',g,ifpre_x(g,iceland),ifpre_y(g,iceland),&
                                 frac_x(g,iceland), frac_y(g,iceland)
            STOP "Glacier fraction .ne. 0"
         END IF
     ELSE
         print *,'Lest Glacier',g,ifpre_x(g,iceland),ifpre_y(g,iceland),&
                                 frac_x(g,iceland), frac_y(g,iceland)
         STOP "Glacier tiles not consistent"
     END IF

     !write(6,*)"TEST_TZ just before call newplant"

     ! For none glacier tiles
     ! Re-calculate plant C, N, P pools
     CALL newplant(cplant_x(g,:,:),frac_x(g,:),ifpre_x(g,:), &
                   cplant_y(g,:,:),frac_y(g,:),ifpre_y(g,:),logc(g,:))
     IF (icycle > 1) CALL newplant(nplant_x(g,:,:),frac_x(g,:),ifpre_x(g,:), &
                                   nplant_y(g,:,:),frac_y(g,:),ifpre_y(g,:),logn(g,:))
     IF (icycle > 2) CALL newplant(pplant_x(g,:,:),frac_x(g,:),ifpre_x(g,:), &
                                   pplant_y(g,:,:),frac_y(g,:),ifpre_y(g,:),logp(g,:))

     write(6,*)"TEST_TZ just after call newplant"

! Lestevens 24 Nov 2017 - Implement Yingping's flux to state pools
     !DATA pool_frac/0.33,0.33,0.34/
     !DATA pool_time/1.0 ,0.1 ,0.01/
     DO y = 1,3 ! pools NOT leaf/wood/root
                      woodhvest_c(g,:,y) = woodhvest_c(g,:,y) + pool_frac(y)*logc(g,:)     !slogc
      IF (icycle > 1) woodhvest_n(g,:,y) = woodhvest_n(g,:,y) + pool_frac(y)*logn(g,:)     !slogn
      IF (icycle > 2) woodhvest_p(g,:,y) = woodhvest_p(g,:,y) + pool_frac(y)*logp(g,:)     !slogp
     END DO
! Lestevens 24 Nov 2017 - Implement Yingping's flux to state pools

     ! Re-calculate litter C, N, P pools
     CALL newlitter(casabiome,frac_x(g,:),ifpre_x(g,:),frac_y(g,:),ifpre_y(g,:), &
                    cplant_x(g,:,:),nplant_x(g,:,:),pplant_x(g,:,:), &
                    cplant_y(g,:,:),nplant_y(g,:,:),pplant_y(g,:,:), &
                    clitter_x(g,:,:),nlitter_x(g,:,:),plitter_x(g,:,:), &
                    clitter_y(g,:,:),nlitter_y(g,:,:),plitter_y(g,:,:))

     !! Re-calculate soil C, N, P pools
     !CALL newsoil(msoil,csoil_x(g,:,:),frac_x(g,:),ifpre_x(g,:),&
     !             csoil_y(g,:,:),frac_y(g,:),ifpre_y(g,:))
     !CALL newsoil(1,clabile_x(g,:),frac_x(g,:),ifpre_x(g,:),&
     !             clabile_y(g,:),frac_y(g,:),ifpre_y(g,:))

     ! Re-calculate soil C, N, P pools
     CALL newsoil(msoil,csoil_x(g,:,:),frac_x(g,:),ifpre_x(g,:),&
                  csoil_y(g,:,:),frac_y(g,:),ifpre_y(g,:))
     CALL newsoil(1,clabile_x(g,:),frac_x(g,:),ifpre_x(g,:),&
                  clabile_y(g,:),frac_y(g,:),ifpre_y(g,:))
     IF (icycle > 1) THEN
        CALL newsoil(msoil,nsoil_x(g,:,:),frac_x(g,:),ifpre_x(g,:),&
                     nsoil_y(g,:,:),frac_y(g,:),ifpre_y(g,:))
        CALL newsoil(1,nsoilmin_x(g,:),frac_x(g,:),ifpre_x(g,:),&
                     nsoilmin_y(g,:),frac_y(g,:),ifpre_y(g,:))
     ENDIF

     IF (icycle > 2) THEN
        CALL newsoil(msoil,psoil_x(g,:,:),frac_x(g,:),ifpre_x(g,:),&
                     psoil_y(g,:,:),frac_y(g,:),ifpre_y(g,:))
        CALL newsoil(1,psoillab_x(g,:),frac_x(g,:),ifpre_x(g,:),&
                     psoillab_y(g,:),frac_y(g,:),ifpre_y(g,:))
        CALL newsoil(1,psoilsorb_x(g,:),frac_x(g,:),ifpre_x(g,:),&
                     psoilsorb_y(g,:),frac_y(g,:),ifpre_y(g,:))
        CALL newsoil(1,psoilocc_x(g,:),frac_x(g,:),ifpre_x(g,:),&
                     psoilocc_y(g,:),frac_y(g,:),ifpre_y(g,:))
     ENDIF


     ! TEST Lestevens 6june18 - thinning forests after luc ----
     IF (l_thinforest) THEN
                      cplant_z(g,:,:) = cplant_y(g,:,:)
      if (icycle > 1) nplant_z(g,:,:) = nplant_y(g,:,:)
      if (icycle > 2) pplant_z(g,:,:) = pplant_y(g,:,:)
      DO y=1,3 ! pools for whvest
                        woodhvest_c(g,:,y) = woodhvest_c(g,:,y) + &
                                   (1-thinning(g,:)) * pool_frac(y) * cplant_y(g,:,wood)
        if (icycle > 1) woodhvest_n(g,:,y) = woodhvest_n(g,:,y) + &
                                   (1-thinning(g,:)) * pool_frac(y) * nplant_y(g,:,wood)
        if (icycle > 2) woodhvest_p(g,:,y) = woodhvest_p(g,:,y) + &
                                   (1-thinning(g,:)) * pool_frac(y) * pplant_y(g,:,wood)
      END DO
      DO y=1,mplant
                        cplant_z(g,:,y) = thinning(g,:) * cplant_y(g,:,y)
        if (icycle > 1) nplant_z(g,:,y) = thinning(g,:) * nplant_y(g,:,y)
        if (icycle > 2) pplant_z(g,:,y) = thinning(g,:) * pplant_y(g,:,y)
      END DO
      CALL newlitter_thin(casabiome,frac_y(g,:),ifpre_y(g,:),frac_y(g,:),ifpre_y(g,:), &
                     cplant_y(g,:,:),nplant_y(g,:,:),pplant_y(g,:,:), &
                     cplant_z(g,:,:),nplant_z(g,:,:),pplant_z(g,:,:), &
                     clitter_y(g,:,:),nlitter_y(g,:,:),plitter_y(g,:,:), &
                     clitter_y(g,:,:),nlitter_y(g,:,:),plitter_y(g,:,:),thinning(g,:))
                      cplant_y(g,:,:) = cplant_z(g,:,:)
      if (icycle > 1) nplant_y(g,:,:) = nplant_z(g,:,:)
      if (icycle > 2) pplant_y(g,:,:) = pplant_z(g,:,:)
     ENDIF
     ! TEST Lestevens 6june18 - thinning forests after luc ----

! Lestevens 24 Nov 2017 - Implement Yingping's flux to state pools
! Lestevens 7 June 2018 - Moved to after thinning
     !DATA pool_frac/0.33,0.33,0.34/
     !DATA pool_time/1.0 ,0.1 ,0.01/
     DO y = 1,3 ! pools NOT leaf/wood/root
                      wresp_c(g,:,y)     = woodhvest_c(g,:,y) * (1.-exp(-1.*pool_time(y))) !
      IF (icycle > 1) wresp_n(g,:,y)     = woodhvest_n(g,:,y) * (1.-exp(-1.*pool_time(y))) !
      IF (icycle > 2) wresp_p(g,:,y)     = woodhvest_p(g,:,y) * (1.-exp(-1.*pool_time(y))) !
                      woodhvest_c(g,:,y) = woodhvest_c(g,:,y) - wresp_c(g,:,y)
      IF (icycle > 1) woodhvest_n(g,:,y) = woodhvest_n(g,:,y) - wresp_n(g,:,y)
      IF (icycle > 2) woodhvest_p(g,:,y) = woodhvest_p(g,:,y) - wresp_p(g,:,y)
     END DO
! Lestevens 24 Nov 2017 - Implement Yingping's flux to state pools


     ! Balance check
     cbal = sum((sum(cplant_x(g,:,:),2) + sum(clitter_x(g,:,:),2)  &
          + sum(csoil_x(g,:,:),2) + clabile_x(g,:)) * frac_x(g,:))  &
          - (sum((sum(cplant_y(g,:,:),2) + sum(clitter_y(g,:,:),2)  &
          + sum(csoil_y(g,:,:),2) + clabile_y(g,:)) * frac_y(g,:)) + sum(logc(g,:)))

     IF (icycle > 1) nbal = sum((sum(nplant_x(g,:,:),2) + sum(nlitter_x(g,:,:),2)  &
                          + sum(nsoil_x(g,:,:),2) + nsoilmin_x(g,:)) * frac_x(g,:))  &
                          - (sum((sum(nplant_y(g,:,:),2) + sum(nlitter_y(g,:,:),2)  &
                          + sum(nsoil_y(g,:,:),2) + nsoilmin_y(g,:)) * frac_y(g,:)) + sum(logn(g,:)))

     IF (icycle > 2) pbal = sum((sum(pplant_x(g,:,:),2) + sum(plitter_x(g,:,:),2)  &
                          + sum(psoil_x(g,:,:),2) + psoillab_x(g,:)  &
                          + psoilsorb_x(g,:) + psoilocc_x(g,:)) * frac_x(g,:))  &
                          - (sum((sum(pplant_y(g,:,:),2) + sum(plitter_y(g,:,:),2)  &
                          + sum(psoil_y(g,:,:),2) + psoillab_y(g,:) + psoilsorb_y(g,:)  &
                          + psoilocc_y(g,:)) * frac_y(g,:)) + sum(logp(g,:)))


     IF(abs(cbal)>1.e-3 .or.abs(nbal)>1.e-4 .or.abs(pbal)>1.e-4) THEN
        print*, 'imbalance on grid:',g,cbal,nbal,pbal
     END IF

  END DO ! end main loop
  !END DO ! end main loop

  ! write values back into cnp pool files
  cpool_tile(:,:,1)  = clabile_y(:,:)
  cpool_tile(:,:,2)  = cplant_y(:,:,1)
  cpool_tile(:,:,3)  = cplant_y(:,:,2)
  cpool_tile(:,:,4)  = cplant_y(:,:,3)
  cpool_tile(:,:,5)  = clitter_y(:,:,1)
  cpool_tile(:,:,6)  = clitter_y(:,:,2)
  cpool_tile(:,:,7)  = clitter_y(:,:,3)
  cpool_tile(:,:,8)  = csoil_y(:,:,1)
  cpool_tile(:,:,9)  = csoil_y(:,:,2)
  cpool_tile(:,:,10) = csoil_y(:,:,3)

  ! if n is switched on
  IF (icycle > 1) THEN
     npool_tile(:,:,1)  = nplant_y(:,:,1)
     npool_tile(:,:,2)  = nplant_y(:,:,2)
     npool_tile(:,:,3)  = nplant_y(:,:,3)
     npool_tile(:,:,4)  = nlitter_y(:,:,1)
     npool_tile(:,:,5)  = nlitter_y(:,:,2)
     npool_tile(:,:,6)  = nlitter_y(:,:,3)
     npool_tile(:,:,7)  = nsoil_y(:,:,1)
     npool_tile(:,:,8)  = nsoil_y(:,:,2)
     npool_tile(:,:,9)  = nsoil_y(:,:,3)
     npool_tile(:,:,10) = nsoilmin_y(:,:)     
  END IF

  ! if p is switched on
  IF (icycle > 2) THEN
     ppool_tile(:,:,1)  = pplant_y(:,:,1)
     ppool_tile(:,:,2)  = pplant_y(:,:,2)
     ppool_tile(:,:,3)  = pplant_y(:,:,3)
     ppool_tile(:,:,4)  = plitter_y(:,:,1)
     ppool_tile(:,:,5)  = plitter_y(:,:,2)
     ppool_tile(:,:,6)  = plitter_y(:,:,3)
     ppool_tile(:,:,7)  = psoil_y(:,:,1)
     ppool_tile(:,:,8)  = psoil_y(:,:,2)
     ppool_tile(:,:,9)  = psoil_y(:,:,3)
     ppool_tile(:,:,10) = psoillab_y(:,:)
     ppool_tile(:,:,11) = psoilsorb_y(:,:)
     ppool_tile(:,:,12) = psoilocc_y(:,:)
  END IF

  ! pack everything
  CALL pack_cnppool(casamet,casapool,casabal,phen,cpool_tile,npool_tile, &
                    ppool_tile,GLAI,PHENPHASE)


  ! update LAI, Vcmax
  IF (icycle > 1) call casa_feedback(ktau_gl,veg,casabiome,casapool,casamet)

  DO p=1,mp
     IF (casamet%iveg2(p) == icewater) THEN
        casamet%glai(p)   = 0.0
     ELSE
        casamet%glai(p)   = MIN(casabiome%glaimax(veg%iveg(p)), MAX(casabiome%glaimin(veg%iveg(p)), &
                               casabiome%sla(veg%iveg(p)) * casapool%cplant(p,leaf)))
     ENDIF
  ! PRINT *, 'p,ivt,glai,vcmax = ', p,veg%iveg(p),casamet%glai(p),veg%vcmax(p)*1.0e6
  END DO

  

END SUBROUTINE casa_reinit_pk

!========================================================================
!========================================================================
!========================================================================

SUBROUTINE pack_cnppool(casamet,casapool,casabal,phen,cpool_tile,npool_tile, &
!SUBROUTINE pack_cnppool(casamet,casapool,casabal,phen,um1,cpool_tile,npool_tile, &
                        ppool_tile,GLAI,PHENPHASE)

    USE cable_um_tech_mod, ONLY : um1
    USE casadimension        ! icycle,mplant,mlitter,msoil
    USE casavariable
    USE phenvariable

IMPLICIT NONE

! passed in
    TYPE (casa_met)     , INTENT(INOUT) :: casamet
    TYPE (casa_pool)    , INTENT(INOUT) :: casapool
    TYPE (casa_balance) , INTENT(INOUT) :: casabal
    TYPE (phen_variable), INTENT(INOUT) :: phen
!    TYPE (um_dimensions), INTENT(IN)    :: um1
!   LOGICAL, INTENT(INOUT),DIMENSION(um1%land_pts, um1%ntiles) :: L_tile_pts
    REAL   , INTENT(INOUT) :: cpool_tile(um1%land_pts,um1%ntiles,10)
    REAL   , INTENT(INOUT) :: npool_tile(um1%land_pts,um1%ntiles,10)
    REAL   , INTENT(INOUT) :: ppool_tile(um1%land_pts,um1%ntiles,12)
    REAL   , INTENT(INOUT) :: GLAI(um1%land_pts,um1%ntiles)
    REAL   , INTENT(INOUT) ::PHENPHASE(um1%land_pts,um1%ntiles)
    !INTEGER, INTENT(INOUT) :: PHENPHASE(um1%land_pts,um1%ntiles)
! local vars
    INTEGER k
    REAL,DIMENSION(um1%land_pts,um1%ntiles)         :: clabile,nsoilmin,psoillab,&
                                                       psoilsorb,psoilocc
    REAL,DIMENSION(um1%land_pts,um1%ntiles,mplant)  :: cplant,nplant,pplant
    REAL,DIMENSION(um1%land_pts,um1%ntiles,mlitter) :: clitter,nlitter,plitter
    REAL,DIMENSION(um1%land_pts,um1%ntiles,msoil)   :: csoil,nsoil,psoil
    INTEGER,DIMENSION(um1%land_pts,um1%ntiles)      :: phenph

! initialise variables
    clabile(:,:)  = 0.0 
    cplant(:,:,:) = 0.0; clitter(:,:,:) = 0.0; csoil(:,:,:) = 0.0
    nplant(:,:,:) = 0.0; nlitter(:,:,:) = 0.0; nsoil(:,:,:) = 0.0 
    nsoilmin(:,:) = 0.0
    pplant(:,:,:) = 0.0; plitter(:,:,:) = 0.0; psoil(:,:,:) = 0.0 
    psoillab(:,:) = 0.0; psoilsorb(:,:) = 0.0; psoilocc(:,:)= 0.0
    phenph(:,:)   = 0

! set to appropriate pools
    clabile(:,:)   = cpool_tile(:,:,1)
    cplant(:,:,1)  = cpool_tile(:,:,2)
    cplant(:,:,2)  = cpool_tile(:,:,3)
    cplant(:,:,3)  = cpool_tile(:,:,4)
    clitter(:,:,1) = cpool_tile(:,:,5)
    clitter(:,:,2) = cpool_tile(:,:,6)
    clitter(:,:,3) = cpool_tile(:,:,7)
    csoil(:,:,1)   = cpool_tile(:,:,8)
    csoil(:,:,2)   = cpool_tile(:,:,9)
    csoil(:,:,3)   = cpool_tile(:,:,10)

    if (icycle >1) then
     nplant(:,:,1)  = npool_tile(:,:,1)
     nplant(:,:,2)  = npool_tile(:,:,2)
     nplant(:,:,3)  = npool_tile(:,:,3)
     nlitter(:,:,1) = npool_tile(:,:,4)
     nlitter(:,:,2) = npool_tile(:,:,5)
     nlitter(:,:,3) = npool_tile(:,:,6)
     nsoil(:,:,1)   = npool_tile(:,:,7)
     nsoil(:,:,2)   = npool_tile(:,:,8)
     nsoil(:,:,3)   = npool_tile(:,:,9)
     nsoilmin(:,:)  = npool_tile(:,:,10)
    endif

    if (icycle >2) then
     pplant(:,:,1)  = ppool_tile(:,:,1)
     pplant(:,:,2)  = ppool_tile(:,:,2)
     pplant(:,:,3)  = ppool_tile(:,:,3)
     plitter(:,:,1) = ppool_tile(:,:,4)
     plitter(:,:,2) = ppool_tile(:,:,5)
     plitter(:,:,3) = ppool_tile(:,:,6)
     psoil(:,:,1)   = ppool_tile(:,:,7)
     psoil(:,:,2)   = ppool_tile(:,:,8)
     psoil(:,:,3)   = ppool_tile(:,:,9)
     psoillab(:,:)  = ppool_tile(:,:,10)
     psoilsorb(:,:) = ppool_tile(:,:,11)
     psoilocc(:,:)  = ppool_tile(:,:,12)
    endif

! pack variables

!! test initialization of LAI
!    casamet%glai = pack(GLAI(:,:)/2.  ,um1%l_tile_pts)
!! test end

    casamet%glai = pack(GLAI(:,:)  ,um1%l_tile_pts)
    phenph = INT(PHENPHASE)
    phen%phase   = pack(phenph(:,:),um1%l_tile_pts)

    casapool%clabile       = pack(clabile(:,:)  ,um1%l_tile_pts)
    do k=1,3

!! test initilization of plant pool  
!     casapool%cplant(:,k)  = pack(cplant(:,:,k)/2. ,um1%l_tile_pts)
!! test end

     casapool%cplant(:,k)  = pack(cplant(:,:,k) ,um1%l_tile_pts)
     casapool%clitter(:,k) = pack(clitter(:,:,k),um1%l_tile_pts)
     casapool%csoil(:,k)   = pack(csoil(:,:,k)  ,um1%l_tile_pts)
    enddo

    if (icycle>1) then
    do k=1,3
     casapool%nplant(:,k)  = pack(nplant(:,:,k) ,um1%l_tile_pts)
     casapool%nlitter(:,k) = pack(nlitter(:,:,k),um1%l_tile_pts)
     casapool%nsoil(:,k)   = pack(nsoil(:,:,k)  ,um1%l_tile_pts)
    enddo
     casapool%nsoilmin     = pack(nsoilmin(:,:) ,um1%l_tile_pts)
    endif

    if (icycle>2) then
    do k=1,3
     casapool%pplant(:,k)  = pack(pplant(:,:,k) ,um1%l_tile_pts)
     casapool%plitter(:,k) = pack(plitter(:,:,k),um1%l_tile_pts)
     casapool%psoil(:,k)   = pack(psoil(:,:,k)  ,um1%l_tile_pts)
    enddo
     casapool%psoillab     = pack(psoillab(:,:) ,um1%l_tile_pts)
     casapool%psoilsorb    = pack(psoilsorb(:,:),um1%l_tile_pts)
     casapool%psoilocc     = pack(psoilocc(:,:) ,um1%l_tile_pts)
    endif

END SUBROUTINE pack_cnppool

!========================================================================
!========================================================================
!========================================================================

SUBROUTINE casa_poolout_unpk(casapool,casaflux,casamet,casabal,phen,  &
                             cpool_tile,npool_tile,ppool_tile,    &
!                             um1,cpool_tile,npool_tile,ppool_tile,    &
                             GLAI,PHENPHASE)

  USE cable_um_tech_mod, ONLY : um1
  USE cable_def_types_mod
  !USE define_dimensions    ! mp, r_1, r_2, i_d
  USE casadimension        ! icycle,mplant,mlitter,msoil
!!   USE define_types
!   USE casaparm
    USE casavariable
    USE phenvariable


IMPLICIT NONE
    TYPE (casa_pool),         INTENT(INOUT) :: casapool
    TYPE (casa_flux),         INTENT(INOUT) :: casaflux
    TYPE (casa_met),          INTENT(INOUT) :: casamet
    TYPE (casa_balance),      INTENT(INOUT) :: casabal
    TYPE (phen_variable),     INTENT(INOUT) :: phen
!    TYPE (um_dimensions),     INTENT(IN)    :: um1
!   LOGICAL, INTENT(INOUT),DIMENSION(um1%land_pts, um1%ntiles) :: L_tile_pts
    REAL   , INTENT(INOUT) :: cpool_tile(um1%land_pts,um1%ntiles,10)
    REAL   , INTENT(INOUT) :: npool_tile(um1%land_pts,um1%ntiles,10)
    REAL   , INTENT(INOUT) :: ppool_tile(um1%land_pts,um1%ntiles,12)
    REAL   , INTENT(INOUT) :: GLAI(um1%land_pts,um1%ntiles)
    REAL   , INTENT(INOUT) ::PHENPHASE(um1%land_pts,um1%ntiles)
    !INTEGER, INTENT(INOUT) :: PHENPHASE(um1%land_pts,um1%ntiles)
! local variables
    REAL(r_2), DIMENSION(mso) :: Psorder,pweasoil,xpsoil50
    REAL(r_2), DIMENSION(mso) :: fracPlab,fracPsorb,fracPocc,fracPorg
    REAL(r_2), DIMENSION(mp)  :: totpsoil
    INTEGER  npt,nso

  ! Soiltype     soilnumber soil P(g P/m2)
  ! Alfisol   1       61.3
  ! Andisol   2       103.9
  ! Aridisol  3       92.8
  ! Entisol   4       136.9
  ! Gellisol  5       98.2
  ! Histosol  6       107.6
  ! Inceptisol        7      84.1
  ! Mollisol  8       110.1
  ! Oxisol    9       35.4
  ! Spodosol  10      41.0
  ! Ultisol   11      51.5
  ! Vertisol  12      190.6
  DATA psorder/61.3,103.9,92.8,136.9,98.2,107.6,84.1,110.1,35.4,41.0,51.5,190.6/
  DATA pweasoil/0.05,0.04,0.03,0.02,0.01,0.009,0.008,0.007,0.006,0.005,0.004,0.003/
  DATA fracpLab/0.08,0.08,0.10,0.02,0.08,0.08,0.08,0.06,0.02,0.05,0.09,0.05/
  DATA fracPsorb/0.32,0.37,0.57,0.67,0.37,0.37,0.37,0.32,0.24,0.22,0.21,0.38/
  DATA fracPocc/0.36,0.38,0.25,0.26,0.38,0.38,0.38,0.44,0.38,0.38,0.37,0.45/
  DATA fracPorg/0.25,0.17,0.08,0.05,0.17,0.17,0.17,0.18,0.36,0.35,0.34,0.12/
  DATA xpsoil50/7.6,4.1,4.2,3.4,4.1,4.1,4.8,4.1,6.9,6.9,6.9,1.7/

!  WRITE(*,91) nyear,cplantsum,clittersum,csoilsum
  casabal%sumcbal=MIN(9999.0,MAX(-9999.0,casabal%sumcbal))
  casabal%sumnbal=MIN(9999.0,MAX(-9999.0,casabal%sumnbal))
  casabal%sumpbal=MIN(9999.0,MAX(-9999.0,casabal%sumpbal))

  DO npt =1, mp
    nso = casamet%isorder(npt)
    totpsoil(npt) = psorder(nso) *xpsoil50(nso)

    IF (icycle<2) THEN
      casapool%nplant(npt,:) = casapool%rationcplant(npt,:)  &
                             * casapool%cplant(npt,:)
      casapool%nlitter(npt,:)= casapool%rationclitter(npt,:) &
                             * casapool%clitter(npt,:)
      casapool%nsoil(npt,:)  = casapool%ratioNCsoil(npt,:)   &
                             * casapool%Csoil(npt,:)
      casapool%nsoilmin(npt) = 2.0
      casabal%sumnbal(npt)   = 0.0
    ENDIF

    IF (icycle<3) THEN
      casabal%sumpbal(npt)   = 0.0
      casapool%pplant(npt,:) = casapool%ratiopcplant(npt,:)  &
                             * casapool%cplant(npt,:)
      casapool%plitter(npt,:)= casapool%ratiopclitter(npt,:) &
                             * casapool%clitter(npt,:)
      casapool%psoil(npt,:)  = casapool%ratioPCsoil(npt,:)   &
                             * casapool%Csoil(npt,:)
      casapool%psoillab(npt) = totpsoil(npt) *fracpLab(nso)
      casapool%psoilsorb(npt)= casaflux%psorbmax(npt) * casapool%psoillab(npt) &
                                /(casaflux%kmlabp(npt)+casapool%psoillab(npt))
      casapool%psoilocc(npt) = totpsoil(npt) *fracPocc(nso)
    ENDIF
  ENDDO

! Lestevens 28 Jan 2011
!    CALL unpack_cnppool(casamet,casapool,casabal,phen,um1,cpool_tile,npool_tile, &
    CALL unpack_cnppool(casamet,casapool,casabal,phen,cpool_tile,npool_tile, &
                        ppool_tile,GLAI,PHENPHASE)


END SUBROUTINE casa_poolout_unpk

!========================================================================
!========================================================================
!========================================================================

!SUBROUTINE unpack_cnppool(casamet,casapool,casabal,phen,um1,cpool_tile,npool_tile, &
SUBROUTINE unpack_cnppool(casamet,casapool,casabal,phen,cpool_tile,npool_tile, &
                          ppool_tile,GLAI,PHENPHASE)

! Lest 25 Jan 2011 - casa_poolout

  USE cable_um_tech_mod, ONLY : um1
  USE cable_def_types_mod
  !USE define_dimensions    ! mp, r_1, r_2, i_d
  USE casadimension        ! icycle,mplant,mlitter,msoil
  USE casavariable
  USE phenvariable

 IMPLICIT NONE

! passed in
    TYPE (casa_met)                  :: casamet
    TYPE (casa_pool)                 :: casapool
    TYPE (casa_balance)              :: casabal
    TYPE (phen_variable)             :: phen
!    TYPE (um_dimensions), INTENT(IN) :: um1
!   LOGICAL, INTENT(INOUT),DIMENSION(um1%land_pts, um1%ntiles) :: L_tile_pts
    REAL   , INTENT(INOUT) :: cpool_tile(um1%land_pts,um1%ntiles,10)
    REAL   , INTENT(INOUT) :: npool_tile(um1%land_pts,um1%ntiles,10)
    REAL   , INTENT(INOUT) :: ppool_tile(um1%land_pts,um1%ntiles,12)
    REAL   , INTENT(INOUT) :: GLAI(um1%land_pts,um1%ntiles)
    REAL   , INTENT(INOUT) ::PHENPHASE(um1%land_pts,um1%ntiles)
    !INTEGER, INTENT(INOUT) :: PHENPHASE(um1%land_pts,um1%ntiles)
! local vars
    INTEGER k
    REAL(r_2) :: miss   = 0.0
    !INTEGER   :: i_miss = 0
    REAL, DIMENSION(um1%land_pts,um1%ntiles)         :: clab,nsmin,pslab,&
                                                        pssorb,psocc,sumcbal,&
                                                        sumnbal,sumpbal
    REAL, DIMENSION(um1%land_pts,um1%ntiles,mplant)  :: cpl,npl,ppl
    REAL, DIMENSION(um1%land_pts,um1%ntiles,mlitter) :: clit,nlit,plit
    REAL, DIMENSION(um1%land_pts,um1%ntiles,msoil)   :: cs,ns,ps

! initialise variables
    clab(:,:)  = 0.0 
    cpl(:,:,:) = 0.0; clit(:,:,:) = 0.0; cs(:,:,:) = 0.0
    npl(:,:,:) = 0.0; nlit(:,:,:) = 0.0; ns(:,:,:) = 0.0  
    nsmin(:,:) = 0.0
    ppl(:,:,:) = 0.0; plit(:,:,:) = 0.0; ps(:,:,:) = 0.0  
    pslab(:,:) = 0.0; pssorb(:,:) = 0.0; psocc(:,:)= 0.0
    sumcbal(:,:)  = 0.0; sumnbal(:,:)   = 0.0; sumpbal(:,:) = 0.0

! unpack variables

    GLAI      = unpack(casamet%glai,um1%l_tile_pts,miss)
    PHENPHASE = unpack(REAL(phen%phase),um1%l_tile_pts,miss)
!   PHENPHASE = unpack(phen%phase,um1%l_tile_pts,i_miss)

    do k= 1,3
     cpl(:,:,k)  = unpack(casapool%cplant(:,k) ,um1%l_tile_pts,miss)
     clit(:,:,k) = unpack(casapool%clitter(:,k),um1%l_tile_pts,miss)
     cs(:,:,k)   = unpack(casapool%csoil(:,k)  ,um1%l_tile_pts,miss)
     npl(:,:,k)  = unpack(casapool%nplant(:,k) ,um1%l_tile_pts,miss)
     nlit(:,:,k) = unpack(casapool%nlitter(:,k),um1%l_tile_pts,miss)
     ns(:,:,k)   = unpack(casapool%nsoil(:,k)  ,um1%l_tile_pts,miss)
     ppl(:,:,k)  = unpack(casapool%pplant(:,k) ,um1%l_tile_pts,miss)
     plit(:,:,k) = unpack(casapool%plitter(:,k),um1%l_tile_pts,miss)
     ps(:,:,k)   = unpack(casapool%psoil(:,k)  ,um1%l_tile_pts,miss)
    enddo

    clab   = unpack(casapool%clabile  ,um1%l_tile_pts,miss)
    nsmin  = unpack(casapool%nsoilmin ,um1%l_tile_pts,miss)
    pslab  = unpack(casapool%psoillab ,um1%l_tile_pts,miss)
    pssorb = unpack(casapool%psoilsorb,um1%l_tile_pts,miss)
    psocc  = unpack(casapool%psoilocc ,um1%l_tile_pts,miss)
    sumcbal   = unpack(casabal%sumcbal   ,um1%l_tile_pts,miss)
    sumnbal   = unpack(casabal%sumnbal   ,um1%l_tile_pts,miss)
    sumpbal   = unpack(casabal%sumpbal   ,um1%l_tile_pts,miss)

! prognostics
    cpool_tile(:,:,1)  =  clab(:,:)    
    cpool_tile(:,:,2)  =  cpl(:,:,1)   
    cpool_tile(:,:,3)  =  cpl(:,:,2)   
    cpool_tile(:,:,4)  =  cpl(:,:,3)   
    cpool_tile(:,:,5)  =  clit(:,:,1)  
    cpool_tile(:,:,6)  =  clit(:,:,2)  
    cpool_tile(:,:,7)  =  clit(:,:,3)  
    cpool_tile(:,:,8)  =  cs(:,:,1)    
    cpool_tile(:,:,9)  =  cs(:,:,2)    
    cpool_tile(:,:,10) =  cs(:,:,3)    
!   sumcbal(:,:)       =  0.0

!   if (icycle > 1) then
     npool_tile(:,:,1)  =  npl(:,:,1)   
     npool_tile(:,:,2)  =  npl(:,:,2)   
     npool_tile(:,:,3)  =  npl(:,:,3)   
     npool_tile(:,:,4)  =  nlit(:,:,1)  
     npool_tile(:,:,5)  =  nlit(:,:,2)  
     npool_tile(:,:,6)  =  nlit(:,:,3)  
     npool_tile(:,:,7)  =  ns(:,:,1)    
     npool_tile(:,:,8)  =  ns(:,:,2)    
     npool_tile(:,:,9)  =  ns(:,:,3)    
     npool_tile(:,:,10) =  nsmin(:,:)   
!    sumnbal(:,:)       =  0.0
!   endif

!   if (icycle > 2) then
     ppool_tile(:,:,1)  =  ppl(:,:,1)   
     ppool_tile(:,:,2)  =  ppl(:,:,2)   
     ppool_tile(:,:,3)  =  ppl(:,:,3)   
     ppool_tile(:,:,4)  =  plit(:,:,1)  
     ppool_tile(:,:,5)  =  plit(:,:,2)  
     ppool_tile(:,:,6)  =  plit(:,:,3)  
     ppool_tile(:,:,7)  =  ps(:,:,1)    
     ppool_tile(:,:,8)  =  ps(:,:,2)    
     ppool_tile(:,:,9)  =  ps(:,:,3)    
     ppool_tile(:,:,10) =  pslab(:,:)   
     ppool_tile(:,:,11) =  pssorb(:,:)  
     ppool_tile(:,:,12) =  psocc(:,:)   
!    sumpbal(:,:)       =  0.0
!   endif

END SUBROUTINE unpack_cnppool

!========================================================================
!========================================================================
!========================================================================

SUBROUTINE unpack_glai(casamet,phen,GLAI,PHENPHASE)

! Lest 15 Jan 2013 - casa_poolout

  USE cable_um_tech_mod, ONLY : um1
  USE cable_def_types_mod
  USE casadimension
  USE casavariable
  USE phenvariable

 IMPLICIT NONE

! passed in
    TYPE (casa_met), INTENT(INOUT)      :: casamet
    TYPE (phen_variable), INTENT(INOUT) :: phen
    REAL, INTENT(INOUT)                 :: GLAI(um1%land_pts,um1%ntiles)
    REAL, INTENT(INOUT)                 ::PHENPHASE(um1%land_pts,um1%ntiles)
    !INTEGER, INTENT(INOUT)              :: PHENPHASE(um1%land_pts,um1%ntiles)
! local vars
    REAL(r_2) :: miss   = 0.0
    !INTEGER   :: i_miss = 0

! unpack variables

    GLAI      = unpack(casamet%glai,um1%l_tile_pts,miss)
    PHENPHASE = unpack(REAL(phen%phase),um1%l_tile_pts,miss)
!   PHENPHASE = unpack(phen%phase,um1%l_tile_pts,i_miss)

END SUBROUTINE unpack_glai

!========================================================================
!========================================================================
!========================================================================

   !--- Lestevens 23Nov11: Based on Jhan's um2cable_lp but for Integers.
   !--- UM met forcing vars needed by CABLE which have UM dimensions
   !---(land_points)[_lp], which is no good to cable. These have to be
   !--- re-packed in a single vector of active tiles. Hence we use
   !--- conditional "mask" l_tile_pts(land_pts,ntiles) which is .true.
   !--- if the land point is/has an active tile

   SUBROUTINE um2cable_ilp(umvar, defaultin, cablevar, soiltype, skip )
         USE cable_def_types_mod, ONLY : mp
         !USE define_dimensions, ONLY : mp
         USE cable_um_tech_mod, ONLY : um1
      implicit none
      INTEGER, INTENT(IN), DIMENSION(um1%land_pts) :: umvar
      INTEGER, INTENT(IN), DIMENSION(10)           :: defaultin
      INTEGER, INTENT(INOUT), DIMENSION(mp)        :: cablevar
      INTEGER, INTENT(INOUT), DIMENSION(mp)        :: soiltype
      INTEGER, DIMENSION(:,:), allocatable         :: fvar
      LOGICAL, optional                            :: skip
      INTEGER                                      :: n,k,l,i


         allocate( fvar(um1%land_pts,um1%ntiles) )
         fvar = 0.0

         DO N=1,um1%NTILES
            DO K=1,um1%TILE_PTS(N)
               L = um1%TILE_INDEX(K,N)
               fvar(L,N) = umvar(L)
               if(.not. present(skip) ) then
                  if( N == um1%ntiles ) then
                     fvar(L,N) =  defaultin(9)
                  endif
               endif
            ENDDO
         ENDDO

         cablevar     =  pack(fvar,um1%l_tile_pts)

         if(.not. present(skip) ) then
            do i=1,mp
               if(soiltype(i)==9) cablevar(i) =  defaultin(9)
            enddo
         endif

         deallocate(fvar)

         return
   END SUBROUTINE um2cable_ilp

  SUBROUTINE redistr_luc(prev_yr_sfrac,inVar,outVar)
      USE cable_um_tech_mod, ONLY : um1
      IMPLICIT NONE

      REAL, INTENT(IN) ,DIMENSION(um1%land_pts,um1%ntiles) :: prev_yr_sfrac
      REAL, INTENT(IN) ,DIMENSION(um1%land_pts,um1%ntiles) :: inVar
      REAL, INTENT(INOUT),DIMENSION(um1%land_pts,um1%ntiles) :: outVar
      ! local variables
      REAL, DIMENSION(um1%ntiles) :: dfrac
      REAL                        :: tmpVar, Rcount
      INTEGER                     :: L,N
  
      !DO L = 1, um1%LAND_PTS
      !  !dfrac = prev_yr_sfrac(L,:) - um1%tile_frac(L,:)
      !  DO N = 1,um1%NTILES
      !    IF (um1%tile_frac(L,N)>=1.e-6) THEN
      !       outVar(L,N) = inVar(L,N)
      !  !  ELSE
      !  !    itype = maxloc(dfrac)
      !  !    outVar(L,N) = inVar(L,itype)  !check
      !  !    if (prev_yr_sfrac(L,itype)<1.e-6) stop "ERROR: CABLE redistr_luc"
      !    ENDIF
      !  END DO
      !END DO

    DO L = 1, um1%land_pts
      dfrac(:) = um1%tile_frac(L,:) - prev_yr_sfrac(L,:)
      ! Collect all cut out areas from various decreasing tiles for averaging
      tmpVar = 0.0
      Rcount = 0.0
      DO N = 1, um1%ntiles
        IF (dfrac(N) < 0.0) THEN
          tmpVar = tmpVar + ABS(dfrac(N)) * inVar(L,N)
          Rcount = Rcount + ABS(dfrac(N))
        ENDIF
      Enddo
      if (Rcount > 0) tmpVar = tmpVar / Rcount
      
      ! Add the averaged amount to those increasing tiles
      DO N = 1, um1%ntiles
        IF (dfrac(N) > 0.0) THEN   ! those tiles increasing in size
          outVar(L,N) = (dfrac(N) * tmpVar + &
                        prev_yr_sfrac(L,N) * inVar(L,N)) / um1%tile_frac(L,N)
        ELSE   ! those that are decreasing or no change in size
          outVar(L,N) = inVar(L,N)
        ENDIF
      Enddo
    Enddo
  
  END SUBROUTINE redistr_luc

  !SUBROUTINE redistr_luc_i(prev_yr_sfrac,inVar,outVar)
  !   IMPLICIT NONE
  !    REAL, INTENT(IN) ,DIMENSION(um1%land_pts,um1%ntiles) :: prev_yr_sfrac
  !    INTEGER, INTENT(IN) ,DIMENSION(um1%land_pts,um1%ntiles) :: inVar
  !    INTEGER, INTENT(INOUT),DIMENSION(um1%land_pts,um1%ntiles) :: outVar
  !    ! local variables
  !    !REAL, DIMENSION(um1%ntiles) :: dfrac
  !    !INTEGER                     :: itype
  !    INTEGER           :: L,N
  !
  !    DO L = 1, um1%LAND_PTS
  !      !dfrac = prev_yr_sfrac(L,:) - um1%tile_frac(L,:)
  !      DO N = 1,um1%NTILES
  !        IF (um1%tile_frac(L,N)>=1.e-6) THEN
  !           outVar(L,N) = inVar(L,N)
  !      !  ELSE
  !      !    itype = maxloc(dfrac)
  !      !    outVar(L,N) = inVar(L,itype)  !check
  !      !    if (prev_yr_sfrac(L,itype)<1.e-6) stop "ERROR: CABLE redistr_luc"
  !        ENDIF
  !      END DO
  !    END DO
  !
  !END SUBROUTINE redistr_luc_i

SUBROUTINE casa_ndep_pk(nidep)

    USE cable_um_tech_mod, ONLY : um1
    USE cable_params_mod, ONLY : soil => soil_cbl 
    USE cable_def_types_mod
    USE casavariable
    USE casa_types_mod
    USE cable_um_init_subrs_mod, ONLY : um2cable_lp

IMPLICIT NONE

! passed in
    REAL, DIMENSION(um1%land_pts)       :: nidep
! local variables
    LOGICAL                             :: skip =.TRUE.

! pack variables gN/m2/day
    call um2cable_lp(nidep,nidep(1:10),casaflux%Nmindep,soil%isoilm,skip)

END SUBROUTINE casa_ndep_pk

END MODULE casa_um_inout_mod

