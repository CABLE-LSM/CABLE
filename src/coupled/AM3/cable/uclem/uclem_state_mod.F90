! UCLEM state and parameter storage for the ACCESS3 CABLE library.

module uclem_state_mod

use uclem_science
implicit none
private
public uclem_newangle, uclem_alb, uclem_fbeam, calculate_solar_angles
public facetparams, facetdata, hydrodata, vegdata, intldata
public fparmdata, pdiagdata
public uclem_type, uclem_data_type
public alloc_uclem_type, dealloc_uclem_type, assoc_uclem_type, nullify_uclem
public update_uclem_geometry, update_uclem_optical
public normalise_uclem_road_volcp, rebuild_uclem_derived_state

integer, save :: soilunder=1

!integer, dimension(:), allocatable, save :: ufull_g
!logical, save :: uclem_active = .false.
!logical, dimension(:,:), allocatable, save :: upack_g


! UCLEM parameters:
TYPE uclem_data_type

    !LOGICAL :: uclem_active = .FALSE.
    !LOGICAL, ALLOCATABLE :: upack_g(:,:)
    !INTEGER, ALLOCATABLE :: ufull_g(:)
    ! 7 outputs
    REAL, ALLOCATABLE :: u_fg(:)
    REAL, ALLOCATABLE :: u_eg(:)
    REAL, ALLOCATABLE :: u_ts(:)
    REAL, ALLOCATABLE :: u_wf(:)
    REAL, ALLOCATABLE :: u_rn(:)
    REAL, ALLOCATABLE :: u_evspsbl(:)
    REAL, ALLOCATABLE :: u_sbl(:)
    REAL, ALLOCATABLE :: qf(:)
    ! Internal temperature outputs for later pseudo-level export:
    ! 1 roof_t0, 2 road_t0, 3 walle_t0, 4 wallw_t0, 5 slab_t0,
    ! 6 intm_t0, 7 room_t, 8 cnveg_t, 9 rfveg_t
    REAL, ALLOCATABLE :: u_temp9(:,:)
    ! Internal sensible heat flux component outputs for later pseudo-level export:
    ! 1 fg_roof, 2 fg_road, 3 fg_walle, 4 fg_wallw, 5 fg_cnveg,
    ! 6 fg_rfveg, 7 ac
    REAL, ALLOCATABLE :: u_shf9(:,:)
    !

    TYPE(facetdata), ALLOCATABLE :: roof
    TYPE(facetdata), ALLOCATABLE :: road
    TYPE(facetdata), ALLOCATABLE :: walle
    TYPE(facetdata), ALLOCATABLE :: wallw
    TYPE(facetdata), ALLOCATABLE :: slab
    TYPE(facetdata), ALLOCATABLE :: intm
    TYPE(facetdata), ALLOCATABLE :: room
    TYPE(hydrodata), ALLOCATABLE :: rfhyd
    TYPE(hydrodata), ALLOCATABLE :: rdhyd
    TYPE(pdiagdata), ALLOCATABLE :: pd
    !
    TYPE(facetparams), ALLOCATABLE :: fp_roof
    TYPE(facetparams), ALLOCATABLE :: fp_road
    TYPE(facetparams), ALLOCATABLE :: fp_wall
    TYPE(facetparams), ALLOCATABLE :: fp_slab
    TYPE(facetparams), ALLOCATABLE :: fp_intm
    TYPE(vegdata), ALLOCATABLE :: cnveg
    TYPE(vegdata), ALLOCATABLE :: rfveg
    TYPE(intldata), ALLOCATABLE :: intl
    TYPE(fparmdata), ALLOCATABLE :: fp


END TYPE uclem_data_type


TYPE uclem_type


    REAL, POINTER :: u_fg(:)
    REAL, POINTER :: u_eg(:)
    REAL, POINTER :: u_ts(:)
    REAL, POINTER :: u_wf(:)
    REAL, POINTER :: u_rn(:)
    REAL, POINTER :: u_evspsbl(:)
    REAL, POINTER :: u_sbl(:)
    REAL, POINTER :: qf(:)
    REAL, POINTER :: u_temp9(:,:)
    REAL, POINTER :: u_shf9(:,:)

    TYPE(facetdata), POINTER :: roof
    TYPE(facetdata), POINTER :: road
    TYPE(facetdata), POINTER :: walle
    TYPE(facetdata), POINTER :: wallw
    TYPE(facetdata), POINTER :: slab
    TYPE(facetdata), POINTER :: intm
    TYPE(facetdata), POINTER :: room

    TYPE(hydrodata), POINTER :: rfhyd
    TYPE(hydrodata), POINTER :: rdhyd
    TYPE(pdiagdata), POINTER :: pd

    TYPE(facetparams), POINTER :: fp_roof
    TYPE(facetparams), POINTER :: fp_road
    TYPE(facetparams), POINTER :: fp_wall
    TYPE(facetparams), POINTER :: fp_slab
    TYPE(facetparams), POINTER :: fp_intm

    TYPE(vegdata), POINTER :: cnveg
    TYPE(vegdata), POINTER :: rfveg
    TYPE(intldata), POINTER :: intl
    TYPE(fparmdata), POINTER :: fp
END TYPE uclem_type


contains

! from MODULE cable_veg_type_mod
SUBROUTINE alloc_uclem_type(uclem, mp)

    IMPLICIT NONE

    integer, parameter :: maxtype = 8

    TYPE(uclem_data_type), INTENT(INOUT) :: uclem


    ! In-canyon vegetation fraction
    real, dimension(maxtype) ::    csigvegc=(/ 0.38, 0.45, 0.38, 0.34, 0.05, 0.40, 0.30, 0.20 /)
    ! Green roof vegetation fraction
    real, dimension(maxtype) ::    csigvegr=(/ 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /)
    ! Area fraction occupied by buildings
    real, dimension(maxtype) ::   csigmabld=(/ 0.45, 0.40, 0.45, 0.46, 0.65, 0.40, 0.45, 0.50 /)
    ! Building height (m)
    real, dimension(maxtype) ::  cbldheight=(/   6.,   4.,   6.,   8.,  18.,   4.,   8.,  12. /)
    ! Building height to width ratio
    real, dimension(maxtype) ::    chwratio=(/  0.4,  0.2,  0.4,  0.6,   2.,  0.5,   1.,  1.5 /)
    ! Industrial sensible heat flux (W m^-2)
    real, dimension(maxtype) :: cindustryfg=(/   0.,   0.,   0.,   0.,   0.,  10.,  20.,  30. /)
    ! Internal gains sensible heat flux [floor] (W m^-2)
    real, dimension(maxtype) ::   cintgains=(/   5.,   5.,   5.,   5.,   5.,   5.,   5.,   5. /)
    ! Daily averaged traffic sensible heat flux (W m^-2)
    real, dimension(maxtype) ::  ctrafficfg=(/  1.5,  1.5,  1.5,  1.5,  1.5,  1.5,  1.5,  1.5 /)
    ! Comfort temperature (K)
    real, dimension(maxtype) :: cbldtemp=(/ 295.16, 295.16, 295.16, 295.16, 295.16, 295.16, 295.16, 295.16 /)
    ! Roof albedo
    real, dimension(maxtype) ::  croofalpha=(/ 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20 /)    ! (Fortuniak 08) Masson = 0.15
    ! Wall albedo
    real, dimension(maxtype) ::  cwallalpha=(/ 0.30, 0.30, 0.30, 0.30, 0.30, 0.30, 0.30, 0.30 /)    ! (Fortuniak 08) Masson = 0.25
    ! Road albedo
    real, dimension(maxtype) ::  croadalpha=(/ 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10 /)    ! (Fortuniak 08) Masson = 0.08
    ! Canyon veg albedo
    real, dimension(maxtype) ::  cvegalphac=(/ 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20 /)
    ! Roof veg albedo
    real, dimension(maxtype) ::  cvegalphar=(/ 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20 /)
    ! Roof emissivity
    real, dimension(maxtype) ::  croofemiss=(/ 0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90 /)
    ! Wall emissivity
    real, dimension(maxtype) ::  cwallemiss=(/ 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85 /)
    ! Road emissivity
    real, dimension(maxtype) ::  croademiss=(/ 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94 /)
    ! Slab emissivity
    real, dimension(maxtype) ::  cslabemiss=(/ 0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90 /)
    ! Canyon veg emissivity
    real, dimension(maxtype) ::  cvegemissc=(/ 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96 /)
    ! Roof veg emissivity
    real, dimension(maxtype) ::  cvegemissr=(/ 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96 /)
    ! Green roof soil depth
    real, dimension(maxtype) ::   cvegdeptr=(/ 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10 /)
    ! Roughness length of in-canyon vegetation (m)
    real, dimension(maxtype) ::    czovegc=(/   0.1,   0.1,   0.1,   0.1,   0.1,   0.1,   0.1,   0.1 /)
    ! In-canyon vegetation LAI
    real, dimension(maxtype) ::  cvegrlaic=(/   2.0,   2.0,   2.0,   2.0,   2.0,   2.0,   2.0,   2.0 /)
    ! Unconstrained canopy stomatal resistance
    real, dimension(maxtype) :: cvegrsminc=(/ 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0 /)
    ! Roughness length of green roof vegetation (m)
    real, dimension(maxtype) ::    czovegr=(/   0.1,   0.1,   0.1,   0.1,   0.1,   0.1,   0.1,   0.1 /)
    ! Green roof vegetation LAI
    real, dimension(maxtype) ::  cvegrlair=(/   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0 /)
    ! Unconstrained canopy stomatal resistance
    real, dimension(maxtype) :: cvegrsminr=(/ 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0 /)
    ! Soil wilting point (m^3 m^-3)
    real, dimension(maxtype) ::     cswilt=(/  0.18,  0.18,  0.18,  0.18,  0.18,  0.18,  0.18,  0.18 /)
    ! Soil field capacity (m^3 m^-3)
    real, dimension(maxtype) ::       csfc=(/  0.26,  0.26,  0.26,  0.26,  0.26,  0.26,  0.26,  0.26 /)
    ! Soil saturation point (m^3 m^-3)
    real, dimension(maxtype) ::      cssat=(/  0.42,  0.42,  0.42,  0.42,  0.42,  0.42,  0.42,  0.42 /)
    ! Infiltration air volume changes per hour (m^3 m^-3)
    real, dimension(maxtype) ::  cinfilach=(/  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,  0.50 /)
    ! Ventilation air volume changes per hour (m^3 m^-3)
    real, dimension(maxtype) :: cventilach=(/  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00 /)
    ! Fraction of spaces with heating devices
    real, dimension(maxtype) ::  cheatprop=(/  0.5,  0.5,  0.5,  0.5,  1.0,  0.00,  0.00,  0.00 /)
    ! Fraction of spaces with cooling devices
    real, dimension(maxtype) ::  ccoolprop=(/  0.5,  0.5,  0.5,  0.5,  1.0,  0.00,  0.00,  0.00 /)





    INTEGER, INTENT(IN) :: mp
    INTEGER :: ntiles = 1
    INTEGER :: itmp = 1 ! using 1 for now
    integer, parameter :: nl=4                 ! Number of layers (default 4, must be a multiple of 4).
    integer ii,j,nlp
    real, save         :: maxsnowalpha=0.85    ! max snow albedo
    real, save         :: minsnowalpha=0.5     ! min snow albedo
    real, save         :: maxsnowden=300.      ! max snow density (kg m^-3)
    real, save         :: minsnowden=100.      ! min snow density (kg m^-3)
    real, parameter    :: vkar=0.4             ! von Karman constant
    real, parameter    :: waterden=1000.       ! water density (kg m^-3)
    ! should be in uclem_type
    real, dimension(mp) :: tsigveg,tsigmabld,fp_coeffbldheight

    real, dimension(maxtype,nl) :: croofdepth
    real, dimension(maxtype,nl) :: cwalldepth
    real, dimension(maxtype,nl) :: croaddepth
    real, dimension(maxtype,nl) :: cslabdepth
    real, dimension(maxtype,nl) :: croofcp
    real, dimension(maxtype,nl) :: cwallcp
    real, dimension(maxtype,nl) :: croadcp
    real, dimension(maxtype,nl) :: cslabcp
    real, dimension(maxtype,nl) :: crooflambda
    real, dimension(maxtype,nl) :: cwalllambda
    real, dimension(maxtype,nl) :: croadlambda
    real, dimension(maxtype,nl) :: cslablambda

    ! facet array where: rows=maxtypes (landtypes) and columns=nl (material layers)
    nlp=nl/4 ! number of layers in each material segment (over 4 material segments)
    ! depths (m)
    croofdepth= reshape((/ ((0.01/real(nlp), ii=1,maxtype),j=1,nlp),    &
                        ((0.09/real(nlp), ii=1,maxtype),j=1,nlp),    &
                        ((0.40/real(nlp), ii=1,maxtype),j=1,nlp),    &
                        ((0.10/real(nlp), ii=1,maxtype),j=1,nlp) /), &
                        (/maxtype,nl/))
    cwalldepth= reshape((/ ((0.01/real(nlp), ii=1,maxtype),j=1,nlp),    &
                        ((0.04/real(nlp), ii=1,maxtype),j=1,nlp),    &
                        ((0.10/real(nlp), ii=1,maxtype),j=1,nlp),    &
                        ((0.05/real(nlp), ii=1,maxtype),j=1,nlp) /), &
                        (/maxtype,nl/))
    croaddepth= reshape((/ ((0.01/real(nlp), ii=1,maxtype),j=1,nlp),    &
                        ((0.04/real(nlp), ii=1,maxtype),j=1,nlp),    &
                        ((0.45/real(nlp), ii=1,maxtype),j=1,nlp),    &
                        ((3.50/real(nlp), ii=1,maxtype),j=1,nlp) /), &
                        (/maxtype,nl/))
    cslabdepth=reshape((/  ((0.05/real(nlp), ii=1,maxtype),j=1,nlp),    &
                        ((0.05/real(nlp), ii=1,maxtype),j=1,nlp),    &
                        ((0.05/real(nlp), ii=1,maxtype),j=1,nlp),    &
                        ((0.05/real(nlp), ii=1,maxtype),j=1,nlp) /), &
                        (/maxtype,nl/))
    ! heat capacity (J m^-3 K^-1)
    croofcp =   reshape((/ ((2.11E6, ii=1,maxtype),j=1,nlp),    & ! dense concrete (Oke 87)
                        ((2.11E6, ii=1,maxtype),j=1,nlp),    & ! dense concrete (Oke 87)
                        ((0.28E6, ii=1,maxtype),j=1,nlp),    & ! light concrete (Oke 87)
                        ((0.29E6, ii=1,maxtype),j=1,nlp) /), & ! insulation (Oke 87)
                        (/maxtype,nl/))
    cwallcp =   reshape((/ ((1.55E6, ii=1,maxtype),j=1,nlp),    & ! concrete (Mills 93)
                        ((1.55E6, ii=1,maxtype),j=1,nlp),    & ! concrete (Mills 93)
                        ((1.55E6, ii=1,maxtype),j=1,nlp),    & ! concrete (Mills 93)
                        ((0.29E6, ii=1,maxtype),j=1,nlp) /), & ! insulation (Oke 87)
                        (/maxtype,nl/))
    croadcp =   reshape((/ ((1.94E6, ii=1,maxtype),j=1,nlp),    & ! asphalt (Mills 93)
                        ((1.94E6, ii=1,maxtype),j=1,nlp),    & ! asphalt (Mills 93)
                        ((1.28E6, ii=1,maxtype),j=1,nlp),    & ! dry soil (Mills 93)
                        ((1.28E6, ii=1,maxtype),j=1,nlp) /), & ! dry soil (Mills 93)
                        (/maxtype,nl/))
    cslabcp=   reshape((/  ((1.55E6, ii=1,maxtype),j=1,nlp),    & ! concrete (Mills 93)
                        ((1.55E6, ii=1,maxtype),j=1,nlp),    & ! concrete (Mills 93)
                        ((1.55E6, ii=1,maxtype),j=1,nlp),    & ! concrete (Mills 93)
                        ((1.55E6, ii=1,maxtype),j=1,nlp) /), & ! concrete (Mills 93)
                        (/maxtype,nl/))
    ! heat conductivity (W m^-1 K^-1)
    crooflambda=reshape((/ ((1.5100, ii=1,maxtype),j=1,nlp),    & ! dense concrete (Oke 87)
                        ((1.5100, ii=1,maxtype),j=1,nlp),    & ! dense concrete (Oke 87)
                        ((0.0800, ii=1,maxtype),j=1,nlp),    & ! light concrete (Oke 87)
                        ((0.0500, ii=1,maxtype),j=1,nlp) /), & ! insulation (Oke 87)
                        (/maxtype,nl/))
    cwalllambda=reshape((/ ((0.9338, ii=1,maxtype),j=1,nlp),    & ! concrete (Mills 93)
                        ((0.9338, ii=1,maxtype),j=1,nlp),    & ! concrete (Mills 93)
                        ((0.9338, ii=1,maxtype),j=1,nlp),    & ! concrete (Mills 93)
                        ((0.0500, ii=1,maxtype),j=1,nlp) /), & ! insulation (Oke 87)
                        (/maxtype,nl/))
    croadlambda=reshape((/ ((0.7454, ii=1,maxtype),j=1,nlp),    & ! asphalt (Mills 93)
                        ((0.7454, ii=1,maxtype),j=1,nlp),    & ! asphalt (Mills 93)
                        ((0.2513, ii=1,maxtype),j=1,nlp),    & ! dry soil (Mills 93)
                        ((0.2513, ii=1,maxtype),j=1,nlp) /), & ! dry soil (Mills 93)
                        (/maxtype,nl/))
    cslablambda=reshape((/ ((0.9338, ii=1,maxtype),j=1,nlp),    & ! concrete (Mills 93)
                        ((0.9338, ii=1,maxtype),j=1,nlp),    & ! concrete (Mills 93)
                        ((0.9338, ii=1,maxtype),j=1,nlp),    & ! concrete (Mills 93)
                        ((0.9338, ii=1,maxtype),j=1,nlp) /), & ! concrete (Mills 93)
                        (/maxtype,nl/))






    ALLOCATE(uclem%u_fg(mp))
    ALLOCATE(uclem%u_eg(mp))
    ALLOCATE(uclem%u_ts(mp))
    ALLOCATE(uclem%u_wf(mp))
    ALLOCATE(uclem%u_rn(mp))
    ALLOCATE(uclem%u_evspsbl(mp))
    ALLOCATE(uclem%u_sbl(mp))
    ALLOCATE(uclem%qf(mp))
    ALLOCATE(uclem%u_temp9(mp,9))
    ALLOCATE(uclem%u_shf9(mp,9))
    !ALLOCATE(uclem%ufull_g(ntiles))
    !ALLOCATE(uclem%upack_g(mp, ntiles))
    ALLOCATE(uclem%roof)
    ALLOCATE(uclem%road)
    ALLOCATE(uclem%walle)
    ALLOCATE(uclem%wallw)
    ALLOCATE(uclem%slab)
    ALLOCATE(uclem%intm)
    ALLOCATE(uclem%room)
    ALLOCATE(uclem%rfhyd)
    ALLOCATE(uclem%rdhyd)
    ALLOCATE(uclem%pd)
    ALLOCATE(uclem%fp_roof)
    ALLOCATE(uclem%fp_road)
    ALLOCATE(uclem%fp_wall)
    ALLOCATE(uclem%fp_slab)
    ALLOCATE(uclem%fp_intm)
    ALLOCATE(uclem%cnveg)
    ALLOCATE(uclem%rfveg)
    ALLOCATE(uclem%intl)
    ALLOCATE(uclem%fp)

    !fg -> fp
    allocate(uclem%fp%rfvegdepth(mp))
    allocate(uclem%fp%ctime(mp))
    allocate(uclem%fp%bldairtemp(mp))
    allocate(uclem%fp%weekdayload(mp))
    allocate(uclem%fp%hangle(mp))
    allocate(uclem%fp%vangle(mp))
    allocate(uclem%fp%fbeam(mp))
    allocate(uclem%fp%hwratio(mp))
    allocate(uclem%fp%effhwratio(mp))
    allocate(uclem%fp%bldheight(mp))
    allocate(uclem%fp%sigmabld(mp))
    allocate(uclem%fp%industryfg(mp))
    allocate(uclem%fp%intgains_flr(mp))
    allocate(uclem%fp%trafficfg(mp))
    allocate(uclem%fp%swilt(mp))
    allocate(uclem%fp%sfc(mp))
    allocate(uclem%fp%ssat(mp))
    allocate(uclem%fp%heatprop(mp))
    allocate(uclem%fp%coolprop(mp))
    allocate(uclem%fp%intmassn(mp))
    allocate(uclem%fp%bldwidth(mp))
    allocate(uclem%fp%infilach(mp))
    allocate(uclem%fp%ventilach(mp))
    allocate(uclem%fp%sigmau(mp))
    !f_roof -> fp_roof
    allocate(uclem%fp_roof%depth(mp,nl))
    allocate(uclem%fp_roof%lambda(mp,nl))
    allocate(uclem%fp_roof%volcp(mp,nl))
    allocate(uclem%fp_roof%emiss(mp))
    allocate(uclem%fp_roof%alpha(mp))
    !!f_wall -> fp_wall
    allocate(uclem%fp_wall%depth(mp,nl))
    allocate(uclem%fp_wall%lambda(mp,nl))
    allocate(uclem%fp_wall%volcp(mp,nl))
    allocate(uclem%fp_wall%emiss(mp))
    allocate(uclem%fp_wall%alpha(mp))
    !!f_road -> fp_road
    allocate(uclem%fp_road%depth(mp,nl))
    allocate(uclem%fp_road%lambda(mp,nl))
    allocate(uclem%fp_road%volcp(mp,nl))
    allocate(uclem%fp_road%emiss(mp))
    allocate(uclem%fp_road%alpha(mp))
    !!f_slab -> fp_slab
    allocate(uclem%fp_slab%depth(mp,nl))
    allocate(uclem%fp_slab%lambda(mp,nl))
    allocate(uclem%fp_slab%volcp(mp,nl))
    allocate(uclem%fp_slab%emiss(mp))
    allocate(uclem%fp_slab%alpha(mp))
    !!f_intm -> fp_intm
    allocate(uclem%fp_intm%depth(mp,nl))
    allocate(uclem%fp_intm%lambda(mp,nl))
    allocate(uclem%fp_intm%volcp(mp,nl))
    allocate(uclem%fp_intm%emiss(mp))
    allocate(uclem%fp_intm%alpha(mp))
    !!cnveg_g -> cnveg
    allocate(uclem%cnveg%emiss(mp))
    allocate(uclem%cnveg%alpha(mp))
    allocate(uclem%cnveg%sigma(mp))
    allocate(uclem%cnveg%rsmin(mp))
    allocate(uclem%cnveg%zo(mp))
    allocate(uclem%cnveg%lai(mp))
    allocate(uclem%cnveg%temp(mp))
    !!rfveg_g -> rfveg
    allocate(uclem%rfveg%emiss(mp))
    allocate(uclem%rfveg%alpha(mp))
    allocate(uclem%rfveg%sigma(mp))
    allocate(uclem%rfveg%rsmin(mp))
    allocate(uclem%rfveg%zo(mp))
    allocate(uclem%rfveg%lai(mp))
    allocate(uclem%rfveg%temp(mp))
    !!intl_g -> intl
    allocate(uclem%intl%viewf(mp,4,4))
    allocate(uclem%intl%psi(mp,4,4))
    !rfhyd
    allocate(uclem%rfhyd%surfwater(mp))
    allocate(uclem%rfhyd%snow(mp))
    allocate(uclem%rfhyd%den(mp))
    allocate(uclem%rfhyd%snowalpha(mp))
    allocate(uclem%rfhyd%leafwater(mp))
    allocate(uclem%rfhyd%soilwater(mp))
    !rdhyd
    allocate(uclem%rdhyd%surfwater(mp))
    allocate(uclem%rdhyd%snow(mp))
    allocate(uclem%rdhyd%den(mp))
    allocate(uclem%rdhyd%snowalpha(mp))
    allocate(uclem%rdhyd%leafwater(mp))
    allocate(uclem%rdhyd%soilwater(mp))
    !roof
    allocate(uclem%roof%nodetemp(mp,0:nl))
    allocate(uclem%roof%storage(mp,nl))
    !road
    allocate(uclem%road%nodetemp(mp,0:nl))
    allocate(uclem%road%storage(mp,nl))
    !wallw
    allocate(uclem%wallw%nodetemp(mp,0:nl))
    allocate(uclem%wallw%storage(mp,nl))
    !walle
    allocate(uclem%walle%nodetemp(mp,0:nl))
    allocate(uclem%walle%storage(mp,nl))
    !slab
    allocate(uclem%slab%nodetemp(mp,0:nl))
    allocate(uclem%slab%storage(mp,nl))
    !intm
    allocate(uclem%intm%nodetemp(mp,0:nl))
    allocate(uclem%intm%storage(mp,nl))
    !room
    allocate(uclem%room%nodetemp(mp,1))
    allocate(uclem%room%storage(mp,1))
    !pd
    allocate(uclem%pd%lzom(mp))
    allocate(uclem%pd%lzoh(mp))
    allocate(uclem%pd%cndzmin(mp))
    allocate(uclem%pd%cduv(mp))
    allocate(uclem%pd%cdtq(mp))
    allocate(uclem%pd%zol(mp))
    allocate(uclem%pd%tscrn(mp))
    allocate(uclem%pd%qscrn(mp))
    allocate(uclem%pd%uscrn(mp))
    allocate(uclem%pd%u10(mp))
    allocate(uclem%pd%emiss(mp))
    allocate(uclem%pd%surferr(mp))
    allocate(uclem%pd%atmoserr(mp))
    allocate(uclem%pd%u_alb(mp))
    allocate(uclem%pd%u_melt(mp))
    allocate(uclem%pd%d_atmosflux(mp))
    allocate(uclem%pd%surferr_bias(mp))
    allocate(uclem%pd%atmoserr_bias(mp))
    allocate(uclem%pd%bldheat(mp))
    allocate(uclem%pd%bldcool(mp))
    allocate(uclem%pd%traf(mp))
    allocate(uclem%pd%intgains_full(mp))
    allocate(uclem%pd%surfrunoff(mp))
    allocate(uclem%pd%irrig(mp))
    allocate(uclem%pd%acond_vegw(mp))
    allocate(uclem%pd%soilwetness(mp))
    allocate(uclem%pd%soilwater(mp))
    allocate(uclem%pd%transveg(mp))
    allocate(uclem%pd%soilmoist(mp))
    allocate(uclem%pd%delsoilmoist(mp))
    allocate(uclem%pd%rootmoistc(mp))
    allocate(uclem%pd%ulai(mp))
    allocate(uclem%pd%storage_flux(mp))
    allocate(uclem%pd%snowmelt(mp))
    allocate(uclem%pd%delintercept(mp))
    allocate(uclem%pd%snowt(mp))
    allocate(uclem%pd%vegt(mp))
    allocate(uclem%pd%swe(mp))
    allocate(uclem%pd%surfstor(mp))
    allocate(uclem%pd%snowfrac(mp))
    allocate(uclem%pd%salbedo(mp))
    allocate(uclem%pd%calbedo(mp))
    allocate(uclem%pd%taircanyon(mp))
    allocate(uclem%pd%delswe(mp))
    allocate(uclem%pd%ac_heat_on(mp))
    allocate(uclem%pd%ac_cool_on(mp))


    allocate(uclem%pd%frac_sigma(mp))
    allocate(uclem%pd%roof_water_runoff(mp))
    allocate(uclem%pd%roof_snow_runoff(mp))
    allocate(uclem%pd%roof_soil_runoff(mp))
    allocate(uclem%pd%road_water_runoff(mp))
    allocate(uclem%pd%road_snow_runoff(mp))
    allocate(uclem%pd%road_soil_runoff(mp))

    allocate(uclem%pd%acond_roof(mp))
    allocate(uclem%pd%acond_road(mp))
    allocate(uclem%pd%acond_walle(mp))
    allocate(uclem%pd%acond_wallw(mp))
    allocate(uclem%pd%d_tempr(mp))
    allocate(uclem%pd%d_canyontemp(mp))


  allocate(uclem%pd%rg_roof(mp))
  allocate(uclem%pd%rg_road(mp))
  allocate(uclem%pd%rg_walle(mp))
  allocate(uclem%pd%rg_wallw(mp))
  allocate(uclem%pd%rg_vegc(mp))
  allocate(uclem%pd%rg_vegr(mp))
  allocate(uclem%pd%rg_rfsn(mp))
  allocate(uclem%pd%rg_rdsn(mp))


  allocate(uclem%pd%sg_roof(mp))
  allocate(uclem%pd%sg_road(mp))
  allocate(uclem%pd%sg_walle(mp))
  allocate(uclem%pd%sg_wallw(mp))
  allocate(uclem%pd%sg_vegc(mp))
  allocate(uclem%pd%sg_vegr(mp))
  allocate(uclem%pd%sg_rfsn(mp))
  allocate(uclem%pd%sg_rdsn(mp))


  allocate(uclem%pd%d_canyonrgout(mp))
  allocate(uclem%pd%d_roofrgout(mp))

  allocate(uclem%pd%acond_vegc(mp))
  allocate(uclem%pd%acond_vegr(mp))
  allocate(uclem%pd%acond_rfsn(mp))
  allocate(uclem%pd%acond_rdsn(mp))
  allocate(uclem%pd%d_rfsndelta(mp))
  allocate(uclem%pd%d_rdsndelta(mp))
  allocate(uclem%pd%rfveg_sigma(mp))
  allocate(uclem%pd%cnveg_sigma(mp))
  allocate(uclem%pd%rdsntemp(mp))
  allocate(uclem%pd%rfsntemp(mp))
  allocate(uclem%pd%rfveg_temp(mp))
  allocate(uclem%pd%d_ac_inside(mp))
  allocate(uclem%pd%int_infilflux(mp))
  allocate(uclem%pd%ac_coeff(mp))

  allocate(uclem%pd%eg_roof(mp))
  allocate(uclem%pd%eg_vegr(mp))
  allocate(uclem%pd%eg_road(mp))
  allocate(uclem%pd%eg_vegc(mp))
  allocate(uclem%pd%eg_rdsn(mp))
  allocate(uclem%pd%eg_rfsn(mp))

  allocate(uclem%pd%d_totdepth(mp))
    ! Initialise state variables - temp
    uclem%roof%nodetemp=1._8  ! + urbtemp
    uclem%roof%storage =0._8
    uclem%road%nodetemp=1._8  ! + urbtemp
    uclem%road%storage =0._8
    uclem%walle%nodetemp=1._8 ! + urbtemp
    uclem%walle%storage=0._8
    uclem%wallw%nodetemp=1._8 ! + urbtemp
    uclem%wallw%storage=0._8
    uclem%slab%nodetemp=1._8 ! + urbtemp
    uclem%slab%storage=0._8
    uclem%intm%nodetemp=1._8 ! + urbtemp
    uclem%intm%storage=0._8
    uclem%room%nodetemp=1._8  ! + urbtemp
    uclem%room%storage=0._8
    uclem%u_temp9=1.
    uclem%u_shf9=1.

    ! Initialise state variables - hydro roof
    uclem%rfhyd%surfwater=0.
    uclem%rfhyd%snow=0.
    uclem%rfhyd%den=minsnowden
    uclem%rfhyd%snowalpha=maxsnowalpha
    uclem%rfhyd%leafwater=0.
    uclem%rfhyd%soilwater=0.
    uclem%rdhyd%surfwater=0.
    uclem%rdhyd%snow=0.
    uclem%rdhyd%den=minsnowden
    uclem%rdhyd%snowalpha=maxsnowalpha
    uclem%rdhyd%leafwater=0.
    uclem%rdhyd%soilwater=0.25



    !  pd

    uclem%pd%zol=0.                                                    ! updated in uclem_calc
    uclem%pd%tscrn=1.      ! + urbtemp                                   ! updated in uclem_calc
    uclem%pd%qscrn=0.                                                    ! updated in uclem_calc
    uclem%pd%uscrn=0.                                                    ! updated in uclem_calc
    uclem%pd%u10=0.                                                      ! updated in uclem_calc
    uclem%pd%emiss=0.97                                                  ! updated in uclem_calc
    uclem%pd%bldheat=0.
    uclem%pd%bldcool=0.
    uclem%pd%surfrunoff=0.
    uclem%pd%soilwetness=0.
    uclem%pd%soilwater=0.
    uclem%pd%transveg=0.
    ! Placeholders only.  Both are recomputed from the final geometry before
    ! they are read: pd%soilmoist by rebuild_uclem_derived_state during
    ! initialisation, pd%rootmoistc by uclem_eval on every call.  The former
    ! hardcoded road depth of 4 m here was wrong whenever the ancillary set a
    ! different croaddepth, and it also omitted the roof term and the vegetated
    ! area normalisation.
    uclem%pd%soilmoist=0.
    uclem%pd%rootmoistc=0.
    uclem%pd%delsoilmoist=0.
    uclem%pd%irrig=0.
    uclem%pd%acond_vegw=0.1
    uclem%pd%traf=0.
    uclem%pd%intgains_full=0.



    uclem%pd%frac_sigma=1.
    uclem%pd%surferr=0._8
    uclem%pd%atmoserr=0._8
    uclem%pd%u_alb=0.
    uclem%pd%u_melt=0.
    uclem%pd%d_atmosflux=0._8
    uclem%pd%surferr_bias=0._8
    uclem%pd%atmoserr_bias=0._8
    uclem%pd%storage_flux=0._8
    uclem%pd%first_call=.true.

    uclem%pd%delintercept=0.
    uclem%pd%snowt=273.
    uclem%pd%vegt=273.
    uclem%pd%swe=0.
    uclem%pd%surfstor=0.
    uclem%pd%snowfrac=0.
    uclem%pd%salbedo=0.
    uclem%pd%calbedo=0.
    uclem%pd%taircanyon=273.
    uclem%pd%delswe=0.
    uclem%pd%ulai=1.


    uclem%pd%roof_water_runoff=0.
    uclem%pd%roof_snow_runoff=0.
    uclem%pd%roof_soil_runoff=0.
    uclem%pd%road_water_runoff=0.
    uclem%pd%road_snow_runoff=0.
    uclem%pd%road_soil_runoff=0.

    uclem%pd%ac_heat_on=1.
    uclem%pd%ac_cool_on=1.






    tsigveg = csigvegc(itmp)
    tsigmabld = csigmabld(itmp)
    !uclem%cnveg%sigma = max(min(tsigveg/(1.-tsigmabld),1.),0.)
    !uclem%cnveg%sigma = max(min(tsigveg/max(1.-tsigmabld,0.001),1.),0.)
    uclem%cnveg%sigma = max(min(tsigveg/(1.-tsigmabld),1.),0.)

    uclem%rfveg%sigma = max(min(csigvegr(itmp),1.),0.)
    uclem%fp%sigmabld = max(min(tsigmabld,1.),0.)
    uclem%fp%hwratio = chwratio(itmp)


    uclem%fp%industryfg = cindustryfg(itmp)
    uclem%fp%intgains_flr = cintgains(itmp)
    uclem%fp%trafficfg = ctrafficfg(itmp)
    uclem%fp%bldheight = cbldheight(itmp)
    ! The canyon aerodynamic state needs the building height, so it is set
    ! only once fp%bldheight holds a value.  Evaluated earlier, it divided by
    ! uninitialised memory.
    uclem%pd%cndzmin=max(10.,0.1*uclem%fp%bldheight+2.)                 ! updated in uclem_calc
    uclem%pd%lzom=log(uclem%pd%cndzmin/(0.1*uclem%fp%bldheight)) ! updated in uclem_calc
    uclem%pd%lzoh=6.+uclem%pd%lzom ! (Kanda et al 2005)           ! updated in uclem_calc
    uclem%pd%cduv=(vkar/uclem%pd%lzom)**2                         ! updated in uclem_calc
    uclem%pd%cdtq=vkar**2/(uclem%pd%lzom*uclem%pd%lzoh)    ! updated in uclem_calc
    uclem%fp_roof%alpha = croofalpha(itmp)
    uclem%fp_wall%alpha = cwallalpha(itmp)
    uclem%fp_road%alpha = croadalpha(itmp)
    uclem%cnveg%alpha = cvegalphac(itmp)
    uclem%rfveg%alpha = cvegalphar(itmp)
    uclem%fp_roof%emiss = croofemiss(itmp)
    uclem%fp_wall%emiss = cwallemiss(itmp)
    uclem%fp_road%emiss = croademiss(itmp)
    uclem%cnveg%emiss = cvegemissc(itmp)
    uclem%rfveg%emiss = cvegemissr(itmp)
    uclem%fp%bldairtemp = cbldtemp(itmp)  - urbtemp
    uclem%fp%rfvegdepth = cvegdeptr(itmp)


    do ii = 1,nl
    uclem%fp_roof%depth(:,ii) = croofdepth(itmp,ii)
    uclem%fp_wall%depth(:,ii) = cwalldepth(itmp,ii)
    uclem%fp_road%depth(:,ii) = croaddepth(itmp,ii)
    uclem%fp_roof%lambda(:,ii) = crooflambda(itmp,ii)
    uclem%fp_wall%lambda(:,ii) = cwalllambda(itmp,ii)
    uclem%fp_road%lambda(:,ii) = croadlambda(itmp,ii)
    uclem%fp_roof%volcp(:,ii) = croofcp(itmp,ii)
    uclem%fp_wall%volcp(:,ii) = cwallcp(itmp,ii)
    ! Bare material heat capacity.  The soilunder normalisation is not applied
    ! here because it depends on cnveg%sigma and sigmabld, both of which the
    ! ancillary overrides can still change.  normalise_uclem_road_volcp applies
    ! it once, after those overrides are final.
    uclem%fp_road%volcp(:,ii) = croadcp(itmp,ii)
    end do


    uclem%cnveg%zo = czovegc(itmp)
    uclem%cnveg%lai = cvegrlaic(itmp)
    uclem%cnveg%rsmin = cvegrsminc(itmp)/max(uclem%cnveg%lai,1.E-8)
    uclem%rfveg%zo = czovegr(itmp)
    uclem%rfveg%lai = cvegrlair(itmp)
    uclem%rfveg%rsmin = cvegrsminr(itmp)/max(uclem%rfveg%lai,1.E-8)

    uclem%fp%swilt = cswilt(itmp)
    uclem%fp%sfc = csfc(itmp)
    uclem%fp%ssat = cssat(itmp)

    ! for varying internal temperature
    uclem%fp_slab%emiss = cslabemiss(itmp)
    uclem%fp%infilach = cinfilach(itmp)
    uclem%fp%ventilach = cventilach(itmp)
    uclem%fp%heatprop = cheatprop(itmp)
    uclem%fp%coolprop = ccoolprop(itmp)

    do ii = 1,nl
    uclem%fp_slab%depth(:,ii) = cslabdepth(itmp,ii)
    uclem%fp_intm%depth(:,ii) = cslabdepth(itmp,ii)  ! internal mass material same as slab
    uclem%fp_slab%lambda(:,ii) = cslablambda(itmp,ii)
    uclem%fp_intm%lambda(:,ii) = cslablambda(itmp,ii)
    uclem%fp_slab%volcp(:,ii) = cslabcp(itmp,ii)
    uclem%fp_intm%volcp(:,ii) = cslabcp(itmp,ii)
    end do



    ! Here we modify the effective canyon geometry to account for in-canyon vegetation tall vegetation
    fp_coeffbldheight = max(uclem%fp%bldheight-6.*uclem%cnveg%zo,0.2)/uclem%fp%bldheight
    uclem%fp%effhwratio     = uclem%fp%hwratio*fp_coeffbldheight


    uclem%cnveg%temp=1. ! + urbtemp             ! updated in uclem_calc
    uclem%rfveg%temp=1. ! + urbtemp             ! updated in uclem_calc


    call init_internal(uclem%fp)
    call init_lwcoeff(uclem%fp,uclem%intl,mp)


    ! calculated values from uclem_type - to be moved to the right place later




    ! define grid arrays
    ! Default value before ACCESS urban fraction packing overwrites sigmau.
    uclem%fp%sigmau = 0.2
    uclem%fp%vangle=0.
    uclem%fp%hangle=0.
    uclem%fp%ctime=0.
    uclem%fp%fbeam=1.
    uclem%fp%weekdayload=1.0


RETURN
END SUBROUTINE alloc_uclem_type


! Override UCLEM urban geometry from UM ancillary fields (hgt, hwr, wrr).
! Must be called after alloc_uclem_type; re-derives all quantities that
! depend on building height, H/W ratio and building plan-area fraction.
SUBROUTINE update_uclem_geometry(uclem, hgt_mp, hwr_mp, wrr_mp, mp)

  IMPLICIT NONE

  TYPE(uclem_data_type), INTENT(INOUT) :: uclem
  INTEGER,               INTENT(IN)    :: mp
  REAL,                  INTENT(IN)    :: hgt_mp(mp)   ! building height (m)
  REAL,                  INTENT(IN)    :: hwr_mp(mp)   ! height-to-canyon-width ratio
  REAL,                  INTENT(IN)    :: wrr_mp(mp)   ! road width fraction (1 - sigmabld)

  REAL, DIMENSION(mp) :: old_sigmabld, new_sigmabld, tsigveg, fp_coeffbldheight

  ! Recover the absolute canyon-veg fraction (fraction of total plan area)
  ! before we change sigmabld, so we can rescale cnveg%sigma afterwards.
  old_sigmabld = uclem%fp%sigmabld
  tsigveg      = uclem%cnveg%sigma * (1. - old_sigmabld)

  ! --- Apply UM geometry overrides where the ancillary provides valid values ---
  WHERE (hgt_mp > 0.)
    uclem%fp%bldheight = hgt_mp
  END WHERE
  WHERE (hwr_mp > 0.)
    uclem%fp%hwratio = hwr_mp
  END WHERE
  ! wrr is the road-width fraction; sigmabld = 1 - wrr
  WHERE (wrr_mp > 0. .AND. wrr_mp < 1.)
    uclem%fp%sigmabld = MAX(1. - wrr_mp, 0.)
  END WHERE
  new_sigmabld = uclem%fp%sigmabld

  ! Re-derive in-canyon veg fraction relative to updated building fraction
  uclem%cnveg%sigma = MAX(MIN(tsigveg / MAX(1. - new_sigmabld, 0.001), 1.), 0.)

  ! fp_road%volcp is deliberately not touched here.  It holds the bare material
  ! heat capacity, and the soilunder normalisation that depends on cnveg%sigma
  ! and sigmabld is applied once by normalise_uclem_road_volcp after
  ! update_uclem_optical, which can still change cnveg%sigma.

  ! Re-derive effective H/W (accounts for tall in-canyon vegetation reducing
  ! the effective canyon depth seen by radiation)
  fp_coeffbldheight = MAX(uclem%fp%bldheight - 6.*uclem%cnveg%zo, 0.2) &
                      / uclem%fp%bldheight
  uclem%fp%effhwratio = uclem%fp%hwratio * fp_coeffbldheight

  ! Re-initialise building width, number of internal-mass floors, and the
  ! precomputed interior longwave view-factor / psi matrices, all geometry-
  ! dependent and only computed once per run.
  CALL init_internal(uclem%fp)
  CALL init_lwcoeff(uclem%fp, uclem%intl, mp)

END SUBROUTINE update_uclem_geometry

!===============================================================================
SUBROUTINE update_uclem_optical(uclem, albroof_mp, albwall_mp, albroad_mp,     &
                                 sigvegc_mp, albvegc_mp,                        &
                                 emisroof_mp, emiswall_mp, emisroad_mp,         &
                                 emisvegc_mp, mp)
! Override UCLEM optical and vegetation properties from ancillary (STASH 923-931).
! sigvegc_mp is the total-area canyon veg fraction; converted to in-canyon
! fraction using sigmabld, which must already be set by update_uclem_geometry.

  IMPLICIT NONE

  TYPE(uclem_data_type), INTENT(INOUT) :: uclem
  INTEGER,               INTENT(IN)    :: mp
  REAL,                  INTENT(IN)    :: albroof_mp(mp)
  REAL,                  INTENT(IN)    :: albwall_mp(mp)
  REAL,                  INTENT(IN)    :: albroad_mp(mp)
  REAL,                  INTENT(IN)    :: sigvegc_mp(mp)
  REAL,                  INTENT(IN)    :: albvegc_mp(mp)
  REAL,                  INTENT(IN)    :: emisroof_mp(mp)
  REAL,                  INTENT(IN)    :: emiswall_mp(mp)
  REAL,                  INTENT(IN)    :: emisroad_mp(mp)
  REAL,                  INTENT(IN)    :: emisvegc_mp(mp)

  WHERE (albroof_mp  > 0.) uclem%fp_roof%alpha = albroof_mp
  WHERE (albwall_mp  > 0.) uclem%fp_wall%alpha  = albwall_mp
  WHERE (albroad_mp  > 0.) uclem%fp_road%alpha  = albroad_mp
  WHERE (albvegc_mp  > 0.) uclem%cnveg%alpha    = albvegc_mp
  WHERE (emisroof_mp > 0.) uclem%fp_roof%emiss  = emisroof_mp
  WHERE (emiswall_mp > 0.) uclem%fp_wall%emiss  = emiswall_mp
  WHERE (emisroad_mp > 0.) uclem%fp_road%emiss  = emisroad_mp
  WHERE (emisvegc_mp > 0.) uclem%cnveg%emiss    = emisvegc_mp
  ! Convert total-area veg fraction to in-canyon fraction; guard on albroof
  ! as proxy for valid ancil coverage (same cells where geometry is set).
  WHERE (albroof_mp > 0.)
    uclem%cnveg%sigma = MAX(MIN(sigvegc_mp / MAX(1. - uclem%fp%sigmabld, 0.001), 1.), 0.)
  END WHERE

END SUBROUTINE update_uclem_optical

!===============================================================================
! Apply the soilunder normalisation to the road volumetric heat capacity.
!
! soilunder selects how much of the gridcell the soil column under the road is
! taken to represent.  Modes 1 to 3 spread the same bulk storage over a smaller
! plan area, so the effective volumetric heat capacity is scaled up by the
! reciprocal of that area.  The scaling therefore depends on cnveg%sigma and
! sigmabld.
!
! This is applied here rather than where croadcp is read, because both
! update_uclem_geometry and update_uclem_optical can still change those two
! fractions.  update_uclem_optical in particular replaces cnveg%sigma outright
! from the ancillary.  Normalising early and rescaling afterwards requires
! recovering the bare material value from the already scaled one, which is what
! the previous version did and which silently produced a road heat capacity
! inconsistent with the final vegetation fraction whenever the sigvegc ancillary
! disagreed with the surface type table.  Working from the bare value once, at
! the end, removes the ordering dependence.
!
! Must be called after update_uclem_geometry and update_uclem_optical, and
! before rebuild_uclem_derived_state, which uses volcp to rebuild the storage.
! Calling it more than once would apply the scaling twice.
!===============================================================================
SUBROUTINE normalise_uclem_road_volcp(uclem, mp)

  IMPLICIT NONE

  TYPE(uclem_data_type), INTENT(INOUT) :: uclem
  INTEGER,               INTENT(IN)    :: mp

  INTEGER :: ii
  REAL, DIMENSION(mp) :: veg_norm, bld_norm

  veg_norm = MAX(1. - uclem%cnveg%sigma,   0.001)
  bld_norm = MAX(1. - uclem%fp%sigmabld,   0.001)

  SELECT CASE(soilunder)
    CASE(0)  ! storage under road only, no dependence on the fractions
      ! nothing to do

    CASE(1)  ! storage under road and canyon vegetation
      DO ii = 1, nl
        uclem%fp_road%volcp(:,ii) = uclem%fp_road%volcp(:,ii) / veg_norm
      END DO

    CASE(2)  ! storage under road and buildings
      DO ii = 1, nl
        uclem%fp_road%volcp(:,ii) = uclem%fp_road%volcp(:,ii)                   &
                                  * ( (1. / veg_norm) * (1. / bld_norm - 1.) + 1. )
      END DO

    CASE(3)  ! storage under road, canyon vegetation and buildings
      DO ii = 1, nl
        uclem%fp_road%volcp(:,ii) = uclem%fp_road%volcp(:,ii)                   &
                                  / veg_norm / bld_norm
      END DO

    CASE DEFAULT
      WRITE(6,*) "ERROR: Unknown soilunder mode ", soilunder
      STOP
  END SELECT

END SUBROUTINE normalise_uclem_road_volcp

!===============================================================================
! Thin wrapper over uclem_rebuild_derived, which lives beside energyclosure in
! uclem_science because it mirrors the storage expressions there and needs the
! module-private constant aircp.
!
! Must be called after update_uclem_geometry and update_uclem_optical.  The
! rebuild uses bldheight, sigmabld, the vegetation fractions, the facet depths
! and the volumetric heat capacities, all of which those two routines set from
! the ancillaries.  Calling it earlier would rebuild against default parameters.
!===============================================================================
SUBROUTINE rebuild_uclem_derived_state(uclem, mp)

  IMPLICIT NONE

  TYPE(uclem_data_type), INTENT(INOUT) :: uclem
  INTEGER,               INTENT(IN)    :: mp

  CALL uclem_rebuild_derived(uclem%roof, uclem%road, uclem%walle, uclem%wallw,  &
                             uclem%slab, uclem%intm, uclem%room,                &
                             uclem%fp, uclem%fp_roof, uclem%fp_road,            &
                             uclem%fp_wall, uclem%fp_slab, uclem%fp_intm,       &
                             uclem%rdhyd, uclem%rfhyd,                          &
                             uclem%cnveg, uclem%rfveg, uclem%pd, mp)

END SUBROUTINE rebuild_uclem_derived_state










SUBROUTINE dealloc_uclem_type(uclem)
    IMPLICIT NONE
    TYPE(uclem_data_type), INTENT(INOUT) :: uclem

    DEALLOCATE(uclem%u_fg)
    DEALLOCATE(uclem%u_eg)
    DEALLOCATE(uclem%u_ts)
    DEALLOCATE(uclem%u_wf)
    DEALLOCATE(uclem%u_rn)
    DEALLOCATE(uclem%u_evspsbl)
    DEALLOCATE(uclem%u_sbl)
    DEALLOCATE(uclem%qf)
    DEALLOCATE(uclem%u_temp9)
    DEALLOCATE(uclem%u_shf9)

    DEALLOCATE(uclem%roof)
    DEALLOCATE(uclem%road)
    DEALLOCATE(uclem%walle)
    DEALLOCATE(uclem%wallw)
    DEALLOCATE(uclem%slab)
    DEALLOCATE(uclem%intm)
    DEALLOCATE(uclem%room)
    DEALLOCATE(uclem%rfhyd)
    DEALLOCATE(uclem%rdhyd)
    DEALLOCATE(uclem%pd)
    DEALLOCATE(uclem%fp_roof)
    DEALLOCATE(uclem%fp_road)
    DEALLOCATE(uclem%fp_wall)
    DEALLOCATE(uclem%fp_slab)
    DEALLOCATE(uclem%fp_intm)
    DEALLOCATE(uclem%cnveg)
    DEALLOCATE(uclem%rfveg)
    DEALLOCATE(uclem%intl)
    DEALLOCATE(uclem%fp)


    deallocate(uclem%fp%rfvegdepth)
    deallocate(uclem%fp%ctime)
    deallocate(uclem%fp%bldairtemp)
    deallocate(uclem%fp%weekdayload)
    deallocate(uclem%fp%hangle)
    deallocate(uclem%fp%vangle)
    deallocate(uclem%fp%fbeam)
    deallocate(uclem%fp%hwratio)
    deallocate(uclem%fp%effhwratio)
    deallocate(uclem%fp%bldheight)
    deallocate(uclem%fp%sigmabld)
    deallocate(uclem%fp%industryfg)
    deallocate(uclem%fp%intgains_flr)
    deallocate(uclem%fp%trafficfg)
    deallocate(uclem%fp%swilt)
    deallocate(uclem%fp%sfc)
    deallocate(uclem%fp%ssat)
    deallocate(uclem%fp%heatprop)
    deallocate(uclem%fp%coolprop)
    deallocate(uclem%fp%intmassn)
    deallocate(uclem%fp%bldwidth)
    deallocate(uclem%fp%infilach)
    deallocate(uclem%fp%ventilach)
    deallocate(uclem%fp%sigmau)
!
    !!f_roof -> fp_roof
    deallocate(uclem%fp_roof%depth)
    deallocate(uclem%fp_roof%lambda)
    deallocate(uclem%fp_roof%volcp)
    deallocate(uclem%fp_roof%emiss)
    deallocate(uclem%fp_roof%alpha)
    !!f_wall -> fp_wall
    deallocate(uclem%fp_wall%depth)
    deallocate(uclem%fp_wall%lambda)
    deallocate(uclem%fp_wall%volcp)
    deallocate(uclem%fp_wall%emiss)
    deallocate(uclem%fp_wall%alpha)
    !!f_road -> fp_road
    deallocate(uclem%fp_road%depth)
    deallocate(uclem%fp_road%lambda)
    deallocate(uclem%fp_road%volcp)
    deallocate(uclem%fp_road%emiss)
    deallocate(uclem%fp_road%alpha)
    !!f_slab -> fp_slab
    deallocate(uclem%fp_slab%depth)
    deallocate(uclem%fp_slab%lambda)
    deallocate(uclem%fp_slab%volcp)
    deallocate(uclem%fp_slab%emiss)
    deallocate(uclem%fp_slab%alpha)
    !!f_intm -> fp_intm
    deallocate(uclem%fp_intm%depth)
    deallocate(uclem%fp_intm%lambda)
    deallocate(uclem%fp_intm%volcp)
    deallocate(uclem%fp_intm%emiss)
    deallocate(uclem%fp_intm%alpha)
    !!cnveg_g -> cnveg
    deallocate(uclem%cnveg%emiss)
    deallocate(uclem%cnveg%alpha)
    deallocate(uclem%cnveg%sigma)
    deallocate(uclem%cnveg%rsmin)
    deallocate(uclem%cnveg%zo)
    deallocate(uclem%cnveg%lai)
    deallocate(uclem%cnveg%temp)
    !!rfveg_g -> rfveg
    deallocate(uclem%rfveg%emiss)
    deallocate(uclem%rfveg%alpha)
    deallocate(uclem%rfveg%sigma)
    deallocate(uclem%rfveg%rsmin)
    deallocate(uclem%rfveg%zo)
    deallocate(uclem%rfveg%lai)
    deallocate(uclem%rfveg%temp)
    !!intl_g -> intl
    deallocate(uclem%intl%viewf)
    deallocate(uclem%intl%psi)

    !rfhyd
    deallocate(uclem%rfhyd%surfwater)
    deallocate(uclem%rfhyd%snow)
    deallocate(uclem%rfhyd%den)
    deallocate(uclem%rfhyd%snowalpha)
    deallocate(uclem%rfhyd%leafwater)
    deallocate(uclem%rfhyd%soilwater)
    !rdhyd
    deallocate(uclem%rdhyd%surfwater)
    deallocate(uclem%rdhyd%snow)
    deallocate(uclem%rdhyd%den)
    deallocate(uclem%rdhyd%snowalpha)
    deallocate(uclem%rdhyd%leafwater)
    deallocate(uclem%rdhyd%soilwater)
    !roof
    deallocate(uclem%roof%nodetemp)
    deallocate(uclem%roof%storage)
    !road
    deallocate(uclem%road%nodetemp)
    deallocate(uclem%road%storage)
    !wallw
    deallocate(uclem%wallw%nodetemp)
    deallocate(uclem%wallw%storage)
    !walle
    deallocate(uclem%walle%nodetemp)
    deallocate(uclem%walle%storage)
    !slab
    deallocate(uclem%slab%nodetemp)
    deallocate(uclem%slab%storage)
    !intm
    deallocate(uclem%intm%nodetemp)
    deallocate(uclem%intm%storage)
    !room
    deallocate(uclem%room%nodetemp)
    deallocate(uclem%room%storage)
    !pd
    deallocate(uclem%pd%lzom)
    deallocate(uclem%pd%lzoh)
    deallocate(uclem%pd%cndzmin)
    deallocate(uclem%pd%cduv)
    deallocate(uclem%pd%cdtq)
    deallocate(uclem%pd%zol)
    deallocate(uclem%pd%tscrn)
    deallocate(uclem%pd%qscrn)
    deallocate(uclem%pd%uscrn)
    deallocate(uclem%pd%u10)
    deallocate(uclem%pd%emiss)
    deallocate(uclem%pd%surferr)
    deallocate(uclem%pd%atmoserr)

    deallocate(uclem%pd%u_alb)
    deallocate(uclem%pd%u_melt)
    deallocate(uclem%pd%d_atmosflux)

    deallocate(uclem%pd%surferr_bias)
    deallocate(uclem%pd%atmoserr_bias)
    deallocate(uclem%pd%bldheat)
    deallocate(uclem%pd%bldcool)
    deallocate(uclem%pd%traf)
    deallocate(uclem%pd%intgains_full)
    deallocate(uclem%pd%surfrunoff)
    deallocate(uclem%pd%irrig)
    deallocate(uclem%pd%acond_vegw)
    deallocate(uclem%pd%soilwetness)
    deallocate(uclem%pd%soilwater)
    deallocate(uclem%pd%transveg)
    deallocate(uclem%pd%soilmoist)
    deallocate(uclem%pd%delsoilmoist)
    deallocate(uclem%pd%rootmoistc)
    deallocate(uclem%pd%ulai)
    deallocate(uclem%pd%storage_flux)
    deallocate(uclem%pd%snowmelt)
    deallocate(uclem%pd%delintercept)
    deallocate(uclem%pd%snowt)
    deallocate(uclem%pd%vegt)
    deallocate(uclem%pd%swe)
    deallocate(uclem%pd%surfstor)
    deallocate(uclem%pd%snowfrac)
    deallocate(uclem%pd%salbedo)
    deallocate(uclem%pd%calbedo)
    deallocate(uclem%pd%taircanyon)
    deallocate(uclem%pd%delswe)
    deallocate(uclem%pd%ac_heat_on)
    deallocate(uclem%pd%ac_cool_on)

    deallocate(uclem%pd%frac_sigma)
    deallocate(uclem%pd%roof_water_runoff)
    deallocate(uclem%pd%roof_snow_runoff)
    deallocate(uclem%pd%roof_soil_runoff)
    deallocate(uclem%pd%road_water_runoff)
    deallocate(uclem%pd%road_snow_runoff)
    deallocate(uclem%pd%road_soil_runoff)

    deallocate(uclem%pd%acond_roof)
    deallocate(uclem%pd%acond_road)
    deallocate(uclem%pd%acond_walle)
    deallocate(uclem%pd%acond_wallw)
    deallocate(uclem%pd%d_tempr)
    deallocate(uclem%pd%d_canyontemp)

  deallocate(uclem%pd%rg_roof)
  deallocate(uclem%pd%rg_road)
  deallocate(uclem%pd%rg_walle)
  deallocate(uclem%pd%rg_wallw)
  deallocate(uclem%pd%rg_vegc)
  deallocate(uclem%pd%rg_vegr)
  deallocate(uclem%pd%rg_rfsn)
  deallocate(uclem%pd%rg_rdsn)

  deallocate(uclem%pd%sg_roof)
  deallocate(uclem%pd%sg_road)
  deallocate(uclem%pd%sg_walle)
  deallocate(uclem%pd%sg_wallw)
  deallocate(uclem%pd%sg_vegc)
  deallocate(uclem%pd%sg_vegr)
  deallocate(uclem%pd%sg_rfsn)
  deallocate(uclem%pd%sg_rdsn)

  deallocate(uclem%pd%d_canyonrgout)
  deallocate(uclem%pd%d_roofrgout)


  deallocate(uclem%pd%acond_vegc)
  deallocate(uclem%pd%acond_vegr)
  deallocate(uclem%pd%acond_rfsn)
  deallocate(uclem%pd%acond_rdsn)
  deallocate(uclem%pd%d_rfsndelta)
  deallocate(uclem%pd%d_rdsndelta)
  deallocate(uclem%pd%rfveg_sigma)
  deallocate(uclem%pd%cnveg_sigma)
  deallocate(uclem%pd%rdsntemp)
  deallocate(uclem%pd%rfsntemp)
  deallocate(uclem%pd%rfveg_temp)
  deallocate(uclem%pd%d_ac_inside)
  deallocate(uclem%pd%int_infilflux)
  deallocate(uclem%pd%ac_coeff)


  deallocate(uclem%pd%eg_roof)
  deallocate(uclem%pd%eg_vegr)
  deallocate(uclem%pd%eg_road)
  deallocate(uclem%pd%eg_vegc)
  deallocate(uclem%pd%eg_rdsn)
  deallocate(uclem%pd%eg_rfsn)

  deallocate(uclem%pd%d_totdepth)
RETURN
END SUBROUTINE dealloc_uclem_type



SUBROUTINE assoc_uclem_type(uclem, uclem_data)
    IMPLICIT NONE

    TYPE(uclem_type), INTENT(INOUT) :: uclem
    TYPE(uclem_data_type), INTENT(INOUT), TARGET :: uclem_data
    CALL nullify_uclem(uclem)

    uclem%u_fg => uclem_data%u_fg
    uclem%u_eg => uclem_data%u_eg
    uclem%u_ts => uclem_data%u_ts
    uclem%u_wf => uclem_data%u_wf
    uclem%u_rn => uclem_data%u_rn
    uclem%u_evspsbl => uclem_data%u_evspsbl
    uclem%u_sbl => uclem_data%u_sbl
    uclem%qf    => uclem_data%qf
    uclem%u_temp9 => uclem_data%u_temp9
    uclem%u_shf9 => uclem_data%u_shf9

    uclem%roof => uclem_data%roof
    uclem%road => uclem_data%road
    uclem%walle => uclem_data%walle
    uclem%wallw => uclem_data%wallw
    uclem%slab => uclem_data%slab
    uclem%intm => uclem_data%intm
    uclem%room => uclem_data%room
    uclem%rfhyd => uclem_data%rfhyd
    uclem%rdhyd => uclem_data%rdhyd
    uclem%pd => uclem_data%pd
    uclem%fp => uclem_data%fp
    uclem%fp_roof => uclem_data%fp_roof
    uclem%fp_road => uclem_data%fp_road
    uclem%fp_wall => uclem_data%fp_wall
    uclem%fp_slab => uclem_data%fp_slab
    uclem%cnveg => uclem_data%cnveg
    uclem%fp_intm => uclem_data%fp_intm
    uclem%rfveg => uclem_data%rfveg
    uclem%intl => uclem_data%intl



RETURN
END SUBROUTINE assoc_uclem_type



SUBROUTINE nullify_uclem(uclem)

    IMPLICIT NONE
    TYPE(uclem_type), INTENT(INOUT) :: uclem

    ! output
    NULLIFY(uclem%u_fg)
    NULLIFY(uclem%u_eg)
    NULLIFY(uclem%u_ts)
    NULLIFY(uclem%u_wf)
    NULLIFY(uclem%u_rn)
    NULLIFY(uclem%u_evspsbl)
    NULLIFY(uclem%u_sbl)
    NULLIFY(uclem%qf)
    NULLIFY(uclem%u_temp9)
    NULLIFY(uclem%u_shf9)

    NULLIFY(uclem%roof)
    NULLIFY(uclem%road)
    NULLIFY(uclem%walle)
    NULLIFY(uclem%wallw)
    NULLIFY(uclem%slab)
    NULLIFY(uclem%intm)
    NULLIFY(uclem%room)
    NULLIFY(uclem%rfhyd)
    NULLIFY(uclem%rdhyd)
    NULLIFY(uclem%pd)
    NULLIFY(uclem%fp)
    NULLIFY(uclem%fp_roof)
    NULLIFY(uclem%fp_wall)
    NULLIFY(uclem%fp_road)
    NULLIFY(uclem%fp_slab)
    NULLIFY(uclem%fp_intm)
    NULLIFY(uclem%rfveg)
    NULLIFY(uclem%cnveg)
    NULLIFY(uclem%intl)

RETURN
END SUBROUTINE nullify_uclem





!===============================================================================
! Set UCLEM solar angles and local time from host-provided forcing.
!===============================================================================
subroutine uclem_newangle(uclem,mp,cosin,ctimein,azimuthin)

implicit none

type(uclem_type), intent(inout) :: uclem
integer, intent(in) :: mp
real, dimension(mp), intent(in) :: cosin     ! cosine of zenith angle
real, dimension(mp), intent(in) :: azimuthin ! azimuthal angle
real, dimension(mp), intent(in) :: ctimein   ! local hour (0<=ctime<=1)


    uclem%fp%hangle = 0.5*pi-azimuthin
    uclem%fp%vangle=acos(cosin)
    uclem%fp%ctime=ctimein

return
end subroutine uclem_newangle




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculates the urban contribution to albedo.
! (selected grid points only)

! raw   (.false.=blend, .true.=output only)
! split (0=net albedo, 1=direct albedo, 2=diffuse albedo)

subroutine uclem_alb(uclem, mp, alb, split)

implicit none

integer, intent(in) :: mp
integer albmode
integer, intent(in), optional :: split
real, dimension(mp), intent(inout) :: alb


real, dimension(mp) :: ualb
real, dimension(mp) :: dumfbeam,snowdeltac,snowdeltar
real, dimension(mp) :: tmpalb,wallpsi,roadpsi
real, dimension(mp) :: sg_roof,sg_vegr,sg_road,sg_walle,sg_wallw,sg_vegc,sg_rfsn,sg_rdsn
logical outmode
type(uclem_type), intent(inout) :: uclem


! if ( .not.uclem_active ) return

albmode=0 ! net albedo
if (present(split)) albmode=split


select case(albmode)
    case default ! net albedo
    dumfbeam=uclem%fp%fbeam
    case(1)      ! direct albedo
    dumfbeam=1.
    case(2)      ! diffuse albedo
    dumfbeam=0.
end select


ualb = 0.

call uclem_calc_alb(tmpalb,uclem%fp%hwratio,uclem%fp%effhwratio, &
                    uclem%fp%vangle,uclem%fp%hangle,dumfbeam,    &
                    uclem%cnveg%sigma,uclem%rfveg%sigma,              &
                    uclem%fp%sigmabld,uclem%fp_road%alpha,                &
                    uclem%fp_wall%alpha,uclem%fp_roof%alpha,                &
                    uclem%cnveg%alpha,uclem%rfveg%alpha,              &
                    uclem%rdhyd%snowalpha,                               &
                    uclem%rfhyd%snowalpha,                               &
                    uclem%rdhyd%snow,                                    &
                    uclem%rfhyd%snow)

ualb = tmpalb
alb = ualb

return
end subroutine uclem_alb



subroutine uclem_fbeam(uclem,mp,fbeam)

implicit none

integer, intent(in) :: mp
real, dimension(mp), intent(in) :: fbeam
type(uclem_type), intent(inout) :: uclem


uclem%fp%fbeam=fbeam

return
end subroutine uclem_fbeam



SUBROUTINE calculate_solar_angles(lat_deg, lon_deg, doy, hod_utc, azimuth)
      REAL, INTENT(IN)  :: lat_deg, lon_deg
      REAL, INTENT(IN)  :: doy, hod_utc
      REAL, INTENT(OUT) :: azimuth

      REAL, PARAMETER :: pi = 3.14159265359, deg2rad = pi/180.0
      REAL :: lat_rad, lon_rad, dec, sindec, mean_omega, local_hod


      lat_rad = lat_deg * deg2rad
      lon_rad = lon_deg * deg2rad


      dec = -23.45 * deg2rad * COS(2.0*pi * (doy + 10.0)/365.25)
      sindec = SIN(dec)


      local_hod = hod_utc + lon_deg/15.0


      mean_omega = (local_hod - 12.0) * pi / 12.0


      azimuth = ATAN2(-COS(dec) * SIN(mean_omega), &
                      COS(lat_rad) * sindec - SIN(lat_rad) * COS(dec) * COS(mean_omega))

      !
        IF (azimuth >= pi) THEN
          azimuth = azimuth - 2.0*pi
        ELSEIF (azimuth < -pi) THEN
          azimuth = azimuth + 2.0*pi
        END IF

END SUBROUTINE calculate_solar_angles


end module uclem_state_mod
