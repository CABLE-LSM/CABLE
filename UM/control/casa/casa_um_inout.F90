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
! Purpose: CASA-CNP
!
! Contact: Yingping Wang, Lauren Stevens
!
! History: 2012
!
!
! ==============================================================================

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

!!  USE cable_um_tech_mod    !, ONLY : um1, veg, soil   
!!  USE define_dimensions    ! mp, r_1, r_2, i_d
!!  USE casadimension        ! icycle,mplant,mlitter,msoil
!!  USE cable_common_module  ! ktau_gl, kend_gl

IMPLICIT NONE 

CONTAINS 

!========================================================================
!========================================================================
!========================================================================

SUBROUTINE init_casacnp(sin_theta_latitude,cpool_tile,npool_tile,ppool_tile, &
                           soil_order,nidep,nifix,pwea,pdust,GLAI,PHENPHASE)
! Les 20 Jan 2011
    USE cable_def_types_mod
    !USE define_dimensions
    !USE define_types
    USE cable_um_tech_mod, ONLY : um1, veg, soil, canopy
    USE casavariable
    USE phenvariable
    USE casa_types_mod
    USE casa_inout_module, only : casa_readbiome, casa_readphen

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

!    casafile%cnpbiome ='/home/599/lxs599/surface_data/pftlookup_csiro_v16_17tiles.csv'
!    casafile%phen     ='/home/599/lxs599/surface_data/modis_phenology_v1.txt'
!    casafile%cnppoint =' '
!    casafile%cnpepool =' '
!    casafile%cnpipool =' '
!    casafile%cnpmetin =' '
!    casafile%cnpmetout=' '
!!   casafile%cnpmetin ='/home/599/lxs599/surface_data/casametDH.csv'
!    casafile%cnpflux  =' '

!print *,'Les um_inout',mp,um1%land_pts,um1%ntiles

    CALL alloc_casavariable(casabiome,casapool,casaflux,casamet,casabal,mp)
    CALL alloc_phenvariable(phen,mp)

    CALL casa_readpoint_pk(sin_theta_latitude,veg,soil,casaflux,casamet, &
                          nidep,nifix,pwea,pdust,soil_order)
    CALL casa_readbiome(veg,soil,casabiome,casapool,casaflux,casamet,phen)
    CALL casa_readphen(veg,casamet,phen)
    CALL casa_init_pk(casabiome,casaflux,casamet,casapool,casabal,canopy,phen, &
                     cpool_tile,npool_tile,ppool_tile,GLAI,PHENPHASE)

 return
END SUBROUTINE init_casacnp

!========================================================================
!========================================================================
!========================================================================

SUBROUTINE casa_readpoint_pk(sin_theta_latitude,veg,soil,casaflux,casamet, &
                           nidep,nifix,pwea,pdust,soil_order)

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
! Les Nov2011 - not correct, but areacell not needed/used for UM ?
    casamet%areacell = pack(um1%tile_frac,um1%l_tile_pts)

    sorder = INT(soil_order)
    call um2cable_ilp(sorder,sorder(1:10),casamet%isorder,soil%isoilm,skip)
    call um2cable_lp(nidep,nidep(1:10),annNdep ,soil%isoilm, skip)
    call um2cable_lp(nifix,nifix(1:10),annNfix ,soil%isoilm, skip)
    call um2cable_lp(pwea ,pwea(1:10) ,annPwea ,soil%isoilm, skip)
    call um2cable_lp(pdust,pdust(1:10),annPdust,soil%isoilm, skip)

    casaflux%Nmindep = annNdep/365.0      ! gN/m2/day
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

SUBROUTINE casa_init_pk(casabiome,casaflux,casamet,casapool,casabal,canopy,phen, &
                       cpool_tile,npool_tile,ppool_tile,GLAI,PHENPHASE)
!                       um1,cpool_tile,npool_tile,ppool_tile,GLAI,PHENPHASE)
!  initialize some values in phenology parameters and leaf growth phase

   USE cable_def_types_mod
!!  USE define_dimensions    ! mp, r_1, r_2, i_d
  USE casadimension        ! icycle,mplant,mlitter,msoil
!  USE define_types
  USE cable_um_tech_mod, ONLY : um1!, veg
  USE casaparm             !, ONLY : initcasa
  USE casavariable
  USE phenvariable

IMPLICIT NONE

  TYPE (casa_biome),   INTENT(IN)       :: casabiome
  TYPE (casa_flux),    INTENT(INOUT)    :: casaflux
  TYPE (casa_met),     INTENT(INOUT)    :: casamet
  TYPE (casa_pool),    INTENT(INOUT)    :: casapool
  TYPE (casa_balance), INTENT(INOUT)    :: casabal
  TYPE (canopy_type),  INTENT(INOUT)       :: canopy
  TYPE (phen_variable),   INTENT(INOUT) :: phen
!  TYPE (um_dimensions), INTENT(IN)      :: um1
! LOGICAL, INTENT(INOUT),DIMENSION(um1%land_pts, um1%ntiles) :: L_tile_pts
  REAL   , INTENT(INOUT) :: cpool_tile(um1%land_pts,um1%ntiles,10)
  REAL   , INTENT(INOUT) :: npool_tile(um1%land_pts,um1%ntiles,10)
  REAL   , INTENT(INOUT) :: ppool_tile(um1%land_pts,um1%ntiles,12)
  REAL   , INTENT(INOUT) :: GLAI(um1%land_pts,um1%ntiles)
  REAL   , INTENT(INOUT) ::PHENPHASE(um1%land_pts,um1%ntiles)
  !INTEGER, INTENT(INOUT) :: PHENPHASE(um1%land_pts,um1%ntiles)

!  PRINT *, 'initcasa,mp ', initcasa,mp

! Les 6dec11 - copied from bgcdriver
! implications for cruns if not init. here
  casamet%tairk(:)   = 0.0
  casamet%tsoil(:,:) = 0.0
  casamet%moist(:,:) = 0.0
  casaflux%cgpp(:)   = 0.0
  casaflux%cnpp(:)   = 0.0
  casaflux%Crsoil(:) = canopy%frs*86400.     !0.0
  casaflux%crgplant(:)   = canopy%frp*86400. !0.0
  casaflux%crmplant(:,:) = 0.0
  casaflux%clabloss(:)   = 0.0
!print *,'Les - Crsoil',casaflux%Crsoil

  IF (initcasa==1) THEN
   CALL pack_cnppool(casamet,casapool,casabal,phen,cpool_tile,npool_tile, &
                     ppool_tile,GLAI,PHENPHASE)
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
    casapool%Psoillab     = MAX(1.0e-7,casapool%psoillab)  ! was 2.0, YPW
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

SUBROUTINE pack_cnppool(casamet,casapool,casabal,phen,cpool_tile,npool_tile, &
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
    casamet%glai = pack(GLAI(:,:)  ,um1%l_tile_pts)
    phenph = INT(PHENPHASE)
    phen%phase   = pack(phenph(:,:),um1%l_tile_pts)

    casapool%clabile       = pack(clabile(:,:)  ,um1%l_tile_pts)
    do k=1,3
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
    CALL unpack_cnppool(casamet,casapool,casabal,phen,cpool_tile,npool_tile, &
                        ppool_tile,GLAI,PHENPHASE)


END SUBROUTINE casa_poolout_unpk

!========================================================================
!========================================================================
!========================================================================

SUBROUTINE unpack_cnppool(casamet,casapool,casabal,phen,cpool_tile,npool_tile, &
                          ppool_tile,GLAI,PHENPHASE)

! Les 25 Jan 2011 - casa_poolout

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

!SUBROUTINE unpack_glai(casamet,phen,GLAI,PHENPHASE)
!
!! Les 15 Jan 2013 - casa_poolout
!
!  USE cable_um_tech_mod, ONLY : um1
!  USE cable_def_types_mod
!  USE casadimension
!  USE casavariable
!  USE phenvariable
!
! IMPLICIT NONE
!
!! passed in
!    TYPE (casa_met), INTENT(INOUT)      :: casamet
!    TYPE (phen_variable), INTENT(INOUT) :: phen
!    REAL, INTENT(INOUT)                 :: GLAI(um1%land_pts,um1%ntiles)
!    REAL, INTENT(INOUT)                 ::PHENPHASE(um1%land_pts,um1%ntiles)
!    !INTEGER, INTENT(INOUT)              :: PHENPHASE(um1%land_pts,um1%ntiles)
!! local vars
!    REAL(r_2) :: miss   = 0.0
!    !INTEGER   :: i_miss = 0
!
!! unpack variables
!
!    GLAI      = unpack(casamet%glai,um1%l_tile_pts,miss)
!    PHENPHASE = unpack(REAL(phen%phase),um1%l_tile_pts,miss)
!!   PHENPHASE = unpack(phen%phase,um1%l_tile_pts,i_miss)
!
!END SUBROUTINE unpack_glai

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

END MODULE casa_um_inout_mod

