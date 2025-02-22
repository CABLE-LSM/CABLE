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
! Purpose: 
!          sumcflux  - accumulating carbon fluxes (not required for UM)
!
! Called from: cable_serial for offline version
!              Not currently called/available for ACCESS version
!
! Contact: Yingping.Wang@csiro.au
!
! History: Model development by Yingping Wang, coupling to Mk3L by Bernard Pak
!          ssoil changed to ssnow
!
! ==============================================================================
!CABLE_LSM:This has to be commented for offline

module casa_cable

contains

SUBROUTINE POPdriver(casaflux,casabal,veg, POP)

  USE cable_def_types_mod
  USE cable_common_module, only: cable_runtime
  USE casadimension
  USE casaparm
  USE casavariable
  USE phenvariable
  USE cable_common_module,  ONLY: CurYear, CABLE_USER
  USE TypeDef,              ONLY: i4b, dp
  USE POPMODULE,            ONLY: POPStep
  USE POP_TYPES,            ONLY: POP_TYPE


  IMPLICIT NONE


  TYPE (veg_parameter_type),  INTENT(IN) :: veg  ! vegetation parameters
  TYPE (casa_flux),           INTENT(IN) :: casaflux
  TYPE (casa_balance),        INTENT(IN) :: casabal
  TYPE(POP_TYPE),             INTENT(INOUT) :: POP

  INTEGER                                   :: it, nit
  REAL(dp)                               :: StemNPP(mp,2)
  REAL(dp), allocatable :: NPPtoGPP(:)
  REAL(dp), allocatable ::  LAImax(:)  , Cleafmean(:),  Crootmean(:)
  CHARACTER                                 :: cyear*4
  CHARACTER                                 :: ncfile*99
  !! vh_js !!
  INTEGER, allocatable :: Iw(:) ! array of indices corresponding to woody (shrub or forest) tiles

  if (.NOT.Allocated(LAIMax)) allocate(LAIMax(mp))
  if (.NOT.Allocated(Cleafmean))  allocate(Cleafmean(mp))
  if (.NOT.Allocated(Crootmean)) allocate(Crootmean(mp))
  if (.NOT.Allocated(NPPtoGPP)) allocate(NPPtoGPP(mp))
  if (.NOT.Allocated(Iw)) allocate(Iw(POP%np))

  IF (cable_user%CALL_POP .and. POP%np.gt.0) THEN ! CALL_POP
     Iw = POP%Iwood

     StemNPP(:,1) = casaflux%stemnpp
     StemNPP(:,2) = 0.0
     WHERE (casabal%FCgppyear > 1.e-5 .and. casabal%FCnppyear > 1.e-5  )
        NPPtoGPP = casabal%FCnppyear/casabal%FCgppyear
     ELSEWHERE
        NPPtoGPP = 0.5
     ENDWHERE
     LAImax = casabal%LAImax
     Cleafmean = casabal%cleafmean
     Crootmean = casabal%Crootmean

     CALL POPStep(pop, max(StemNPP(Iw,:)/1000.,0.0001), int(veg%disturbance_interval(Iw,:), i4b),&
          real(veg%disturbance_intensity(Iw,:),dp)      ,&
          max(LAImax(Iw),0.001), Cleafmean(Iw), Crootmean(Iw), NPPtoGPP(Iw))


  ENDIF ! CALL_POP

END SUBROUTINE POPdriver
! ==============================================================================
SUBROUTINE read_casa_dump(  ncfile, casamet, casaflux,phen, climate, ncall, kend, allATonce )
      USE netcdf
      USE cable_def_types_mod,   ONLY : r_2,ms,mp, climate_type
      USE casadimension,         ONLY : mplant,mdyear
      USE casavariable,          ONLY : casa_met, casa_flux
      USE phenvariable
      USE cable_common_module,  ONLY:  CABLE_USER
      USE casa_ncdf_module,     ONLY : get_var_ncr2, &
                                        get_var_ncr3, stderr_nc
      IMPLICIT NONE

      TYPE (casa_flux), INTENT(INOUT) :: casaflux
      TYPE (casa_met), INTENT(inout)  :: casamet
      TYPE (phen_variable),         INTENT(INOUT) :: phen
      TYPE (climate_type), INTENT(INOUT)       :: climate  ! climate variables
      INTEGER, INTENT(in)             :: kend, ncall
      CHARACTER(len=*), INTENT(in)    :: ncfile
      LOGICAL, INTENT(in)             :: allATonce

      !netcdf IDs/ names
      INTEGER            :: num_vars
      INTEGER, PARAMETER :: num_dims=3
      INTEGER, SAVE                        :: ncrid  ! netcdf file ID
      !INTEGER , DIMENSION(num_vars)        :: varrID ! (1) tvair, (2) pmb


      !vars
      CHARACTER, DIMENSION(:), POINTER :: var_name*15

      REAL     , DIMENSION(mp)        :: lat, lon
      REAL(r_2), DIMENSION(mp)        :: tairk,  cgpp, mtemp, Ndep, Pdep
      REAL(r_2), DIMENSION(mp,ms)     :: tsoil, moist
      REAL(r_2), DIMENSION(mp,mplant) :: crmplant
      REAL(r_2), DIMENSION(mp)        :: phenphase, phendoyphase1, &
           phendoyphase2,  phendoyphase3,  phendoyphase4
      INTEGER :: ncok,  idoy


      num_vars=14

      !Add extra mtemp variable when running with climate
      IF(cable_user%CALL_climate) THEN
          num_vars=num_vars+1
      ENDIF

      allocate(var_name(num_vars))

      !Variable names
      var_name =  (/"lat          ", &
                    "lon          ", &
                    "casamet_tairk", &
                    "tsoil        ", &
                    "moist        ", &
                    "cgpp         ", &
                    "crmplant     ", &
                    "phenphase    ", &
                    "phendoyphase1", &
                    "phendoyphase2", &
                    "phendoyphase3", &
                    "phendoyphase4", &
                    "Ndep         ", &
                    "Pdep         "/)

      !Add extra mtemp variable when running with climate
      IF (cable_user%CALL_climate) THEN
        var_name(num_vars)="mtemp"
      ENDIF

 IF ( allATonce .OR. ncall .EQ. 1 ) THEN
         ncok = NF90_OPEN(TRIM(ncfile), nf90_nowrite, ncrid)
         IF (ncok /= nf90_noerr ) CALL stderr_nc(ncok,'re-opening ', ncfile)
      ENDIF
      IF ( allATonce ) THEN
         DO idoy=1,mdyear

            CALL get_var_ncr2(ncrid, var_name(3), tairk   , idoy )
            CALL get_var_ncr3(ncrid, var_name(4), tsoil   , idoy ,ms)
            CALL get_var_ncr3(ncrid, var_name(5), moist   , idoy ,ms)
            CALL get_var_ncr2(ncrid, var_name(6), cgpp    , idoy )
            CALL get_var_ncr3(ncrid, var_name(7), crmplant, idoy ,mplant)
            CALL get_var_ncr2(ncrid, var_name(8), phenphase, idoy)
            CALL get_var_ncr2(ncrid, var_name(9), phendoyphase1, idoy)
            CALL get_var_ncr2(ncrid, var_name(10), phendoyphase2, idoy)
            CALL get_var_ncr2(ncrid, var_name(11), phendoyphase3, idoy)
            CALL get_var_ncr2(ncrid, var_name(12), phendoyphase4, idoy)
            CALL get_var_ncr2(ncrid, var_name(13), Ndep   , idoy )
            CALL get_var_ncr2(ncrid, var_name(14), Pdep   , idoy )
            !amu561 this need to be in if-block
            IF (cable_user%CALL_climate ) THEN
               CALL get_var_ncr2(ncrid, var_name(15), mtemp   , idoy )
            ENDIF

            casamet%Tairkspin(:,idoy) = tairk
            casamet%cgppspin (:,idoy) = cgpp
            casamet%crmplantspin_1(:,idoy) = crmplant(:,1)
            casamet%crmplantspin_2(:,idoy) = crmplant(:,2)
            casamet%crmplantspin_3(:,idoy) = crmplant(:,3)
            casamet%Tsoilspin_1(:,idoy)    = tsoil(:,1)
            casamet%Tsoilspin_2(:,idoy)    = tsoil(:,2)
            casamet%Tsoilspin_3(:,idoy)    = tsoil(:,3)
            casamet%Tsoilspin_4(:,idoy)    = tsoil(:,4)
            casamet%Tsoilspin_5(:,idoy)    = tsoil(:,5)
            casamet%Tsoilspin_6(:,idoy)    = tsoil(:,6)
            casamet%moistspin_1(:,idoy)    = moist(:,1)
            casamet%moistspin_2(:,idoy)    = moist(:,2)
            casamet%moistspin_3(:,idoy)    = moist(:,3)
            casamet%moistspin_4(:,idoy)    = moist(:,4)
            casamet%moistspin_5(:,idoy)    = moist(:,5)
            casamet%moistspin_6(:,idoy)    = moist(:,6)
            phen%phasespin(:,idoy) = int(phenphase)
            phen%doyphasespin_1(:,idoy) = int(phendoyphase1)
            phen%doyphasespin_2(:,idoy) = int(phendoyphase2)
            phen%doyphasespin_3(:,idoy) = int(phendoyphase3)
            phen%doyphasespin_4(:,idoy) = int(phendoyphase4)
            casaflux%Nmindep = Ndep
            casaflux%Pdep = Pdep
            !amu561 this need to be in if-block
            IF(cable_user%CALL_climate) THEN
               casamet%mtempspin(:,idoy) = mtemp
            ENDIF
         END DO
      ELSE

         CALL get_var_ncr2(ncrid, var_name(3), tairk   ,ncall )
         CALL get_var_ncr3(ncrid, var_name(4), tsoil   ,ncall , ms)
         CALL get_var_ncr3(ncrid, var_name(5), moist   ,ncall , ms)
         CALL get_var_ncr2(ncrid, var_name(6), cgpp    ,ncall )
         CALL get_var_ncr3(ncrid, var_name(7), crmplant,ncall , mplant)
         CALL get_var_ncr2(ncrid, var_name(8), phenphase    ,ncall )
         CALL get_var_ncr2(ncrid, var_name(9), phendoyphase1    ,ncall )
         CALL get_var_ncr2(ncrid, var_name(10), phendoyphase2    ,ncall )
         CALL get_var_ncr2(ncrid, var_name(11), phendoyphase3    ,ncall )
         CALL get_var_ncr2(ncrid, var_name(12), phendoyphase4    ,ncall )
         CALL get_var_ncr2(ncrid, var_name(13), Ndep   , ncall )
         CALL get_var_ncr2(ncrid, var_name(14), Pdep   , ncall )
         IF(cable_user%CALL_climate) THEN
            CALL get_var_ncr2(ncrid, var_name(15), mtemp   , ncall )
         ENDIF

         casamet%tairk     = tairk
         casamet%tsoil     = tsoil
         casamet%moist     = moist
         casaflux%cgpp     = cgpp
         casaflux%crmplant = crmplant
         phen%phase = int(phenphase)
         phen%doyphase(:,1) = int(phendoyphase1)
         phen%doyphase(:,2) = int(phendoyphase2)
         phen%doyphase(:,3) = int(phendoyphase3)
         phen%doyphase(:,4) = int(phendoyphase4)
         casaflux%Nmindep = Ndep
         casaflux%Pdep = Pdep
         IF (cable_user%CALL_climate) THEN
            climate%mtemp_max = mtemp
         ENDIF

      ENDIF

      IF ( allATonce .OR. ncall .EQ. kend ) THEN
         ncok = NF90_CLOSE(ncrid)
         IF (ncok /= nf90_noerr ) CALL stderr_nc(ncok,'closing ', ncfile)
      ENDIF

      deallocate(var_name)

   END SUBROUTINE read_casa_dump


SUBROUTINE write_casa_dump( ncfile, casamet, casaflux, phen, climate, n_call, kend )
  USE netcdf
  USE cable_def_types_mod,   ONLY : r_2,ms,mp, climate_type
  USE cable_common_module,   ONLY : kend_gl
  USE casa_ncdf_module,     ONLY : def_dims, def_vars, def_var_atts, &
       put_var_ncr1, put_var_ncr2,       &
       put_var_ncr3, stderr_nc
  USE casavariable,          ONLY : CASA_MET, CASA_FLUX
  USE casadimension,         ONLY : mplant
  USE phenvariable
  USE cable_common_module,  ONLY:  CABLE_USER

  IMPLICIT NONE

  INTEGER, INTENT(in) :: &
       n_call, &         ! this timestep #
       kend              ! final timestep of run

  TYPE (casa_flux),             INTENT(IN) :: casaflux
  TYPE (casa_met),              INTENT(IN) :: casamet
  TYPE (phen_variable),         INTENT(IN) :: phen
  TYPE (climate_type), INTENT(IN)       :: climate  ! climate variables

  !number of instances. dummied here and so=1
  !integer :: inst =1

  !netcdf IDs/ names
  CHARACTER(len=*)   :: ncfile
  INTEGER            :: num_vars
  INTEGER, PARAMETER :: num_dims=3
  INTEGER, SAVE :: ncid       ! netcdf file ID

  !vars
  CHARACTER, DIMENSION(:), POINTER :: var_name*15


  INTEGER, DIMENSION(:), POINTER :: varID ! (1) tvair, (2) pmb

  !dims
  CHARACTER(len=*), DIMENSION(num_dims), PARAMETER :: &
       dim_name =  (/ "pnt ", &
       "soil", &
       "time" /)

  INTEGER, PARAMETER :: soil_dim = 6

  INTEGER, DIMENSION(soil_dim), PARAMETER  :: soil = (/ 1,2,3,4,5,6 /)

  INTEGER, DIMENSION(num_dims)  :: &
       dimID   ! (1) x, (2) y, (3) time

  INTEGER, DIMENSION(num_dims)  :: &
                                !x,y generally lat/lon BUT for single site = 1,1
       dim_len = (/-1,soil_dim,-1/)  ! (1) mp, (2) soil, (3) time [re-set]



  !local only
  INTEGER :: ncok      !ncdf return status

  ! END header
  dim_len(1)        = mp
  dim_len(num_dims) = NF90_unlimited

  num_vars = 14

  !Add extra mtemp variable when running with climate
  IF (cable_user%CALL_climate) THEN
    num_vars=num_vars+1
  ENDIF

  allocate(var_name(num_vars))
  allocate(varID(num_vars))

  var_name =  (/"lat          ", &
                "lon          ", &
                "casamet_tairk", &
                "tsoil        ", &
                "moist        ", &
                "cgpp         ", &
                "crmplant     ", &
                "phenphase    ", &
                "phendoyphase1", &
                "phendoyphase2", &
                "phendoyphase3", &
                "phendoyphase4", &
                "Ndep         ", &
                "Pdep         "/)

  !Add extra mtemp variable when running with climate
  IF (cable_user%CALL_climate) THEN
    var_name(num_vars)="mtemp"
  ENDIF

  IF (n_call == 1) THEN

     ! create netCDF dataset: enter define mode
     ncok = nf90_create(path = TRIM(ncfile), cmode = nf90_clobber, ncid = ncid)
     IF (ncok /= nf90_noerr) CALL stderr_nc(ncok,'ncdf creating ', ncfile)

     !ncok = nf90_redef(ncid)
     !if (ncok /= nf90_noerr) call stderr_nc(ncok,'enter def mode', ncfile)

     ! define dimensions: from name and length
     CALL def_dims(num_dims, ncid, dimID, dim_len, dim_name )

     ! define variables: from name, type, dims
     CALL def_vars(num_vars, ncid,  nf90_float, dimID, var_name, varID )

     ! define variable attributes
     !CLN LATER!             CALL def_var_atts( ncfile, ncid, varID )

     ncok = nf90_enddef(ncid)
     if (ncok /= nf90_noerr) call stderr_nc(ncok,'end def mode', ncfile)


     CALL put_var_ncr1(ncid, var_name(1), REAL(casamet%lat)  )
     CALL put_var_ncr1(ncid, var_name(2), REAL(casamet%lon)  )


  ENDIF


  CALL put_var_ncr2(ncid, var_name(3), casamet%tairk    ,n_call )
  CALL put_var_ncr3(ncid, var_name(4), casamet%tsoil    ,n_call, ms )
  CALL put_var_ncr3(ncid, var_name(5), casamet%moist    ,n_call, ms )
  CALL put_var_ncr2(ncid, var_name(6), casaflux%cgpp    ,n_call )
  CALL put_var_ncr3(ncid, var_name(7), casaflux%crmplant,n_call, mplant )
  CALL put_var_ncr2(ncid, var_name(8), real(phen%phase , r_2)    ,n_call )
  CALL put_var_ncr2(ncid, var_name(9), real(phen%doyphase(:,1), r_2)    ,n_call )
  CALL put_var_ncr2(ncid, var_name(10), real(phen%doyphase(:,2), r_2)    ,n_call )
  CALL put_var_ncr2(ncid, var_name(11), real(phen%doyphase(:,3), r_2)    ,n_call )
  CALL put_var_ncr2(ncid, var_name(12), real(phen%doyphase(:,4), r_2)    ,n_call )
  CALL put_var_ncr2(ncid, var_name(13), real(casaflux%Nmindep,r_2)    ,n_call )
  CALL put_var_ncr2(ncid, var_name(14), real(casaflux%Pdep,r_2)    ,n_call )
  if (cable_user%CALL_climate) then
     CALL put_var_ncr2(ncid, var_name(15), real(climate%mtemp_max,r_2)    ,n_call )
  endif

  deallocate(var_name)

  IF (n_call == kend ) &
       ncok = nf90_close(ncid)            ! close: save new netCDF dataset

END SUBROUTINE write_casa_dump

 SUBROUTINE casa_feedback(ktau,veg,casabiome,casapool,casamet)
  USE cable_def_types_mod
  USE casadimension
  USE casaparm
  USE casavariable
  USE casa_cnp_module, ONLY: vcmax_np
  USE cable_common_module,  ONLY:  CABLE_USER
  IMPLICIT NONE
  INTEGER,      INTENT(IN) :: ktau ! integration step number
  TYPE (veg_parameter_type),  INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (casa_biome),          INTENT(INOUT) :: casabiome
  TYPE (casa_pool),           INTENT(IN) :: casapool
  TYPE (casa_met),            INTENT(IN) :: casamet

  integer np,ivt
  real, dimension(mp)  :: ncleafx,npleafx, pleafx, nleafx ! local variables
  real, dimension(17)                   ::  xnslope
  data xnslope/0.80,1.00,2.00,1.00,1.00,1.00,0.50,1.00,0.34,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00/

  ! first initialize
  ncleafx(:) = casabiome%ratioNCplantmax(veg%iveg(:),leaf)
  npleafx = 14.2

  DO np=1,mp
    ivt=veg%iveg(np)
    IF (casamet%iveg2(np)/=icewater &
        .AND. casamet%glai(np)>casabiome%glaimin(ivt)  &
        .AND. casapool%cplant(np,leaf)>0.0) THEN
      ncleafx(np) = MIN(casabiome%ratioNCplantmax(ivt,leaf), &
                        MAX(casabiome%ratioNCplantmin(ivt,leaf), &
                            casapool%nplant(np,leaf)/casapool%cplant(np,leaf)))
      IF (icycle>2 .AND. casapool%pplant(np,leaf)>0.0) THEN
        npleafx(np) = MIN(30.0,MAX(8.0,real(casapool%nplant(np,leaf) &
                /casapool%pplant(np,leaf))))
      ENDIF
    ENDIF

    IF (TRIM(cable_user%vcmax).eq.'standard') then
       IF (casamet%glai(np) > casabiome%glaimin(ivt)) THEN
          IF (ivt/=2) THEN
             veg%vcmax(np) = ( casabiome%nintercept(ivt) &
                  + casabiome%nslope(ivt)*ncleafx(np)/casabiome%sla(ivt) ) * 1.0e-6
          ELSE
             IF (casapool%nplant(np,leaf)>0.0.AND.casapool%pplant(np,leaf)>0.0) THEN
                veg%vcmax(np) = ( casabiome%nintercept(ivt)  &
                     + casabiome%nslope(ivt)*(0.4+9.0/npleafx(np)) &
                     * ncleafx(np)/casabiome%sla(ivt) ) * 1.0e-6
             ELSE
                veg%vcmax(np) = ( casabiome%nintercept(ivt) &
                     + casabiome%nslope(ivt)*ncleafx(np)/casabiome%sla(ivt) )*1.0e-6
             ENDIF
          ENDIF
          veg%vcmax(np) =veg%vcmax(np)* xnslope(ivt)
       ENDIF

    elseif (TRIM(cable_user%vcmax).eq.'Walker2014') then
       !Walker, A. P. et al.: The relationship of leaf photosynthetic traits – Vcmax and Jmax –
       !to leaf nitrogen, leaf phosphorus, and specific leaf area:
       !a meta-analysis and modeling study, Ecology and Evolution, 4, 3218-3235, 2014.
       ! veg%vcmax(np) = exp(3.946 + 0.921*log(nleafx(np)) + 0.121*log(pleafx(np)) + &
       !      0.282*log(pleafx(np))*log(nleafx(np))) * 1.0e-6
       nleafx(np) = ncleafx(np)/casabiome%sla(ivt) ! leaf N in g N m-2 leaf
       pleafx(np) = nleafx(np)/npleafx(np) ! leaf P in g P m-2 leaf
       if (ivt .EQ. 7) then
          veg%vcmax(np) = 1.0e-5 ! special for C4 grass: set here to value from  parameter file
       else
          veg%vcmax(np) = vcmax_np(nleafx(np), pleafx(np))
       endif
    else
       stop ('invalid vcmax flag')
    endif

  ENDDO

  veg%ejmax = 2.0 * veg%vcmax
!991 format(i6,2x,i4,2x,2(f9.3,2x))
 END SUBROUTINE casa_feedback


SUBROUTINE sumcflux(ktau, kstart, kend, dels, bgc, canopy,  &
     soil, ssnow, sum_flux, veg, met, casaflux, l_vcmaxFeedbk)

  USE cable_def_types_mod
  USE cable_carbon_module
  USE casadimension
  USE casaparm
  USE casavariable
  IMPLICIT NONE
  INTEGER, INTENT(IN)    :: ktau ! integration step number
  INTEGER, INTENT(IN)    :: kstart ! starting value of ktau
  INTEGER, INTENT(IN)    :: kend ! total # timesteps in run
!  INTEGER, INTENT(IN)    :: mvtype  ! Number of veg types
!  INTEGER, INTENT(IN)    :: mstype ! Number of soil types
  REAL,    INTENT(IN)    :: dels ! time setp size (s)
  TYPE (bgc_pool_type),       INTENT(INOUT) :: bgc
  TYPE (canopy_type),         INTENT(INOUT) :: canopy
  TYPE (soil_parameter_type), INTENT(INOUT) :: soil
  TYPE (soil_snow_type),      INTENT(INOUT) :: ssnow
  TYPE (sum_flux_type),       INTENT(INOUT) :: sum_flux
  TYPE (met_type),            INTENT(IN)    :: met
  TYPE (veg_parameter_type),  INTENT(INOUT) :: veg
  TYPE (casa_flux),           INTENT(INOUT) :: casaflux
  LOGICAL, INTENT(IN)   :: l_vcmaxFeedbk ! using prognostic Vcmax

!   if(icycle<=0) then
!     these are executed in cbm
!      CALL soilcarb(soil, ssoil, veg, bgc, met, canopy)
!      CALL carbon_pl(dels, soil, ssoil, veg, canopy, bgc)
!   else
    if(icycle>0) then
       canopy%frp(:) = (casaflux%crmplant(:,wood)+casaflux%crmplant(:,froot) &
                        +casaflux%crgplant(:))/86400.0
       canopy%frs(:) = casaflux%Crsoil(:)/86400.0
       canopy%frpw(:)= casaflux%crmplant(:,wood)/86400.0
       canopy%frpr(:)= casaflux%crmplant(:,froot)/86400.0
    endif
    if(ktau == kstart) then
       sum_flux%sumpn  = canopy%fpn*dels
       sum_flux%sumrd  = canopy%frday*dels
       sum_flux%dsumpn = canopy%fpn*dels
       sum_flux%dsumrd = canopy%frday*dels
       sum_flux%sumrpw = canopy%frpw*dels
       sum_flux%sumrpr = canopy%frpr*dels
       sum_flux%sumrp  = canopy%frp*dels
       sum_flux%dsumrp = canopy%frp*dels
    ! canopy%frs set in soilcarb
       sum_flux%sumrs = canopy%frs*dels
    else
       sum_flux%sumpn  = sum_flux%sumpn  + canopy%fpn*dels
       sum_flux%sumrd  = sum_flux%sumrd  + canopy%frday*dels
       sum_flux%dsumpn = sum_flux%dsumpn + canopy%fpn*dels
       sum_flux%dsumrd = sum_flux%dsumrd + canopy%frday*dels
       sum_flux%sumrpw = sum_flux%sumrpw + canopy%frpw*dels
       sum_flux%sumrpr = sum_flux%sumrpr + canopy%frpr*dels
       sum_flux%sumrp  = sum_flux%sumrp  + canopy%frp*dels
       sum_flux%dsumrp = sum_flux%dsumrp + canopy%frp*dels
    ! canopy%frs set in soilcarb
       sum_flux%sumrs = sum_flux%sumrs+canopy%frs*dels
    endif
    ! Set net ecosystem exchange after adjustments to frs:
    canopy%fnpp = -1.0* canopy%fpn - canopy%frp
    IF (icycle <= 1) THEN
      canopy%fnee = canopy%fpn + canopy%frs + canopy%frp
    ELSE
      IF (l_vcmaxFeedbk) THEN
        canopy%fnee = canopy%fpn + canopy%frs + canopy%frp &
                    + casaflux%clabloss(:)/86400.0
      ELSE
        canopy%fnee = (casaflux%Crsoil-casaflux%cnpp+casaflux%clabloss)/86400.0
      ENDIF
    ENDIF

END SUBROUTINE sumcflux

  SUBROUTINE totcnppools(kloop,veg,casamet,casapool,bmcplant,bmnplant,bmpplant,bmclitter,bmnlitter,bmplitter, &
                               bmcsoil,bmnsoil,bmpsoil,bmnsoilmin,bmpsoillab,bmpsoilsorb,bmpsoilocc,bmarea)
  ! this subroutine is temporary, and its needs to be modified for multiple tiles within a cell
  USE cable_def_types_mod
  USE casadimension
  USE casaparm
  USE casavariable
  IMPLICIT NONE
  INTEGER,                     INTENT(IN)  :: kloop
  TYPE (veg_parameter_type),   INTENT(IN)  :: veg  ! vegetation parameters
  TYPE(casa_pool),             INTENT(IN)  :: casapool
  TYPE(casa_met),              INTENT(IN)  :: casamet
  real,      dimension(5,mvtype,mplant)    :: bmcplant,  bmnplant,  bmpplant
  real,      dimension(5,mvtype,mlitter)   :: bmclitter, bmnlitter, bmplitter
  real,      dimension(5,mvtype,msoil)     :: bmcsoil,   bmnsoil,   bmpsoil
  real,      dimension(5,mvtype)           :: bmnsoilmin,bmpsoillab,bmpsoilsorb, bmpsoilocc
  real,      dimension(mvtype)             :: bmarea
  ! local variables
  INTEGER  npt,nvt


      bmcplant(kloop,:,:)  = 0.0;  bmnplant(kloop,:,:)  = 0.0; bmpplant(kloop,:,:)  = 0.0
      bmclitter(kloop,:,:) = 0.0;  bmnlitter(kloop,:,:) = 0.0; bmplitter(kloop,:,:) = 0.0
      bmcsoil(kloop,:,:)   = 0.0;  bmnsoil(kloop,:,:)   = 0.0; bmpsoil(kloop,:,:)   = 0.0
      bmnsoilmin(kloop,:)  = 0.0;  bmpsoillab(kloop,:)  = 0.0; bmpsoilsorb(kloop,:) = 0.0;  bmpsoilocc(kloop,:) = 0.0

      bmarea(:) = 0.0

      do npt=1,mp
         nvt=veg%iveg(npt)
         bmcplant(kloop,nvt,:) = bmcplant(kloop,nvt,:)   + casapool%cplant(npt,:) * casamet%areacell(npt)
         bmnplant(kloop,nvt,:) = bmnplant(kloop,nvt,:)   + casapool%nplant(npt,:) * casamet%areacell(npt)
         bmpplant(kloop,nvt,:) = bmpplant(kloop,nvt,:)   + casapool%pplant(npt,:) * casamet%areacell(npt)

         bmclitter(kloop,nvt,:) = bmclitter(kloop,nvt,:) + casapool%clitter(npt,:) * casamet%areacell(npt)
         bmnlitter(kloop,nvt,:) = bmnlitter(kloop,nvt,:) + casapool%nlitter(npt,:) * casamet%areacell(npt)
         bmplitter(kloop,nvt,:) = bmplitter(kloop,nvt,:) + casapool%plitter(npt,:) * casamet%areacell(npt)

         bmcsoil(kloop,nvt,:) = bmcsoil(kloop,nvt,:)     + casapool%csoil(npt,:) * casamet%areacell(npt)
         bmnsoil(kloop,nvt,:) = bmnsoil(kloop,nvt,:)     + casapool%nsoil(npt,:) * casamet%areacell(npt)
         bmpsoil(kloop,nvt,:) = bmpsoil(kloop,nvt,:)     + casapool%psoil(npt,:) * casamet%areacell(npt)

         bmnsoilmin(kloop,nvt)  = bmnsoilmin(kloop,nvt)   + casapool%nsoilmin(npt) * casamet%areacell(npt)
         bmpsoillab(kloop,nvt)  = bmpsoillab(kloop,nvt)   + casapool%psoillab(npt) * casamet%areacell(npt)
         bmpsoilsorb(kloop,nvt) = bmpsoilsorb(kloop,nvt)  + casapool%psoilsorb(npt) * casamet%areacell(npt)
         bmpsoilocc(kloop,nvt)  = bmpsoilocc(kloop,nvt)   + casapool%psoilocc(npt) * casamet%areacell(npt)
         bmarea(nvt)  = bmarea(nvt) + casamet%areacell(npt)
      enddo

      do nvt=1,mvtype
         bmcplant(kloop,nvt,:) = bmcplant(kloop,nvt,:)/bmarea(nvt)
         bmnplant(kloop,nvt,:) = bmnplant(kloop,nvt,:)/bmarea(nvt)
         bmpplant(kloop,nvt,:) = bmpplant(kloop,nvt,:)/bmarea(nvt)

         bmclitter(kloop,nvt,:) = bmclitter(kloop,nvt,:)/bmarea(nvt)
         bmnlitter(kloop,nvt,:) = bmnlitter(kloop,nvt,:)/bmarea(nvt)
         bmplitter(kloop,nvt,:) = bmplitter(kloop,nvt,:)/bmarea(nvt)

         bmcsoil(kloop,nvt,:) = bmcsoil(kloop,nvt,:)/bmarea(nvt)
         bmnsoil(kloop,nvt,:) = bmnsoil(kloop,nvt,:)/bmarea(nvt)
         bmpsoil(kloop,nvt,:) = bmpsoil(kloop,nvt,:)/bmarea(nvt)

         bmnsoilmin(kloop,nvt)  = bmnsoilmin(kloop,nvt)/bmarea(nvt)
         bmpsoillab(kloop,nvt)  = bmpsoillab(kloop,nvt)/bmarea(nvt)
         bmpsoilsorb(kloop,nvt) = bmpsoilsorb(kloop,nvt)/bmarea(nvt)
         bmpsoilocc(kloop,nvt)  = bmpsoilocc(kloop,nvt)/bmarea(nvt)
      enddo

  END SUBROUTINE totcnppools

  SUBROUTINE analyticpool(kend,veg,soil,casabiome,casapool,                                 &
                          casaflux,casamet,casabal,phen,                                    &
                          avgcleaf2met,avgcleaf2str,avgcroot2met,avgcroot2str,avgcwood2cwd, &
                          avgnleaf2met,avgnleaf2str,avgnroot2met,avgnroot2str,avgnwood2cwd, &
                          avgpleaf2met,avgpleaf2str,avgproot2met,avgproot2str,avgpwood2cwd, &
                          avgcgpp, avgcnpp, avgnuptake, avgpuptake,                         &
                          avgxnplimit,avgxkNlimiting,avgxklitter,avgxksoil,                 &
                          avgratioNCsoilmic,avgratioNCsoilslow,avgratioNCsoilpass,          &
                          avgnsoilmin,avgpsoillab,avgpsoilsorb,avgpsoilocc,                 &
                          avg_af, avg_aw, avg_ar, avg_lf, avg_lw, avg_lr,                   &
                          avg_annual_cnpp)
  USE cable_def_types_mod
  USE cable_carbon_module
  USE casadimension
  USE casaparm
  USE casavariable
  USE phenvariable
  IMPLICIT NONE
  INTEGER, INTENT(IN)    :: kend
  TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (soil_parameter_type),   INTENT(INOUT) :: soil ! soil parameters
  TYPE (casa_biome),            INTENT(INOUT) :: casabiome
  TYPE (casa_pool),             INTENT(INOUT) :: casapool
  TYPE (casa_flux),             INTENT(INOUT) :: casaflux
  TYPE (casa_met),              INTENT(INOUT) :: casamet
  TYPE (casa_balance),          INTENT(INOUT) :: casabal
  TYPE (phen_variable),         INTENT(INOUT) :: phen


  ! local variables
  real, dimension(mp)      :: avgcleaf2met,avgcleaf2str,avgcroot2met,avgcroot2str,avgcwood2cwd
  real, dimension(mp)      :: avgnleaf2met,avgnleaf2str,avgnroot2met,avgnroot2str,avgnwood2cwd
  real, dimension(mp)      :: avgpleaf2met,avgpleaf2str,avgproot2met,avgproot2str,avgpwood2cwd
  real, dimension(mp)      :: avgcgpp, avgcnpp, avgnuptake, avgpuptake
  real(r_2), dimension(mp) :: avgxnplimit,avgxkNlimiting,avgxklitter,avgxksoil
  real, dimension(mp)      :: avgratioNCsoilmic,  avgratioNCsoilslow,  avgratioNCsoilpass
  real,      dimension(mp) :: avgnsoilmin,avgpsoillab,avgpsoilsorb,avgpsoilocc
  real,      dimension(mp) :: avg_af, avg_aw, avg_ar, avg_lf, avg_lw, avg_lr
  real,      dimension(mp) :: avg_annual_cnpp

  ! local variables
  REAL(r_2), DIMENSION(mso) :: Psorder,pweasoil,xpsoil50
  REAL(r_2), DIMENSION(mso) :: fracPlab,fracPsorb,fracPocc,fracPorg
  REAL(r_2), DIMENSION(mp)  :: totpsoil
  INTEGER  npt,nout,nso

  ! Soiltype     soilnumber soil P(g P/m2)
  ! Alfisol     1       61.3
  ! Andisol     2       103.9
  ! Aridisol    3       92.8
  ! Entisol     4       136.9
  ! Gellisol    5       98.2
  ! Histosol    6       107.6
  ! Inceptisol  7       84.1
  ! Mollisol    8       110.1
  ! Oxisol      9       35.4
  ! Spodosol    10      41.0
  ! Ultisol     11      51.5
  ! Vertisol    12      190.6
  DATA psorder/61.3,103.9,92.8,136.9,98.2,107.6,84.1,110.1,35.4,41.0,51.5,190.6/
  DATA pweasoil/0.05,0.04,0.03,0.02,0.01,0.009,0.008,0.007,0.006,0.005,0.004,0.003/
  DATA fracpLab/0.08,0.08,0.10,0.02,0.08,0.08,0.08,0.06,0.02,0.05,0.09,0.05/
  DATA fracPsorb/0.32,0.37,0.57,0.67,0.37,0.37,0.37,0.32,0.24,0.22,0.21,0.38/
  DATA fracPocc/0.36,0.38,0.25,0.26,0.38,0.38,0.38,0.44,0.38,0.38,0.37,0.45/
  DATA fracPorg/0.25,0.17,0.08,0.05,0.17,0.17,0.17,0.18,0.36,0.35,0.34,0.12/
  DATA xpsoil50/7.6,4.1,4.2,3.4,4.1,4.1,4.8,4.1,6.9,6.9,6.9,1.7/

  INTEGER :: nyear,iyear
  real year

  ! compute the mean litter input in g(C, N and P)/day from plant pools
  casaflux%fromLtoS = 0.0
  casaflux%fromStoS = 0.0

  casabal%sumcbal(:)   = 0.0
  casabal%sumnbal(:)   = 0.0
  casabal%sumpbal(:)   = 0.0

  do npt=1,mp
  if(casamet%iveg2(npt)/=icewater.and.avgcnpp(npt) > 0.0) THEN
    casaflux%fromLtoS(npt,mic,metb)   = 0.45
                                          ! metb -> mic
    casaflux%fromLtoS(npt,mic,str)   = 0.45*(1.0-casabiome%fracLigninplant(veg%iveg(npt),leaf))
                                          ! str -> mic
    casaflux%fromLtoS(npt,slow,str)  = 0.7 * casabiome%fracLigninplant(veg%iveg(npt),leaf)
                                          ! str -> slow
    casaflux%fromLtoS(npt,mic,cwd)   = 0.40*(1.0 - casabiome%fracLigninplant(veg%iveg(npt),wood))
                                          ! CWD -> fmic
    casaflux%fromLtoS(npt,slow,cwd)  = 0.7 * casabiome%fracLigninplant(veg%iveg(npt),wood)
                                          ! CWD -> slow
!! set the following two backflow to set (see Bolker 199x)
!    casaflux%fromStoS(npt,mic,slow)  = 0.45 * (0.997 - 0.009 *soil%clay(npt))
!    casaflux%fromStoS(npt,mic,pass)  = 0.45

    casaflux%fromStoS(npt,slow,mic)  = (0.85 - 0.68 * (soil%clay(npt)+soil%silt(npt))) &
                                     * (0.997 - 0.032*soil%clay(npt))
    casaflux%fromStoS(npt,pass,mic)  = (0.85 - 0.68 * (soil%clay(npt)+soil%silt(npt))) &
                                     * (0.003 + 0.032*soil%clay(npt))
    casaflux%fromStoS(npt,pass,slow) = 0.45 * (0.003 + 0.009 * soil%clay(npt) )

    casaflux%klitter(npt,metb) = avgxkNlimiting(npt) * avgxklitter(npt)*casabiome%litterrate(veg%iveg(npt),metb)
    casaflux%klitter(npt,str)  = avgxkNlimiting(npt) * avgxklitter(npt)*casabiome%litterrate(veg%iveg(npt),str)&
                               * exp(-3.0*casabiome%fracLigninplant(veg%iveg(npt),leaf))
    casaflux%klitter(npt,cwd)  = avgxkNlimiting(npt) * avgxklitter(npt)*casabiome%litterrate(veg%iveg(npt),cwd)

    casaflux%ksoil(npt,mic)    = avgxksoil(npt)*casabiome%soilrate(veg%iveg(npt),mic)  &
                               * (1.0 - 0.75 *(soil%silt(npt)+soil%clay(npt)))
    casaflux%ksoil(npt,slow)   = avgxksoil(npt) * casabiome%soilrate(veg%iveg(npt),slow)
    casaflux%ksoil(npt,pass)   = avgxksoil(npt) * casabiome%soilrate(veg%iveg(npt),pass)

    if(veg%iveg(npt)==cropland) THEN     ! for cultivated land type
       casaflux%ksoil(npt,mic)  = casaflux%ksoil(npt,mic) * 1.25
       casaflux%ksoil(npt,slow) = casaflux%ksoil(npt,slow)* 1.5
       casaflux%ksoil(npt,pass) = casaflux%ksoil(npt,pass)* 1.5
    endif
  endif
  enddo


  do npt=1,mp
   if(casamet%iveg2(npt)/=icewater.and.avgcnpp(npt) > 0.0) then

      ! Compute steady-state plant C pool sizes
      casapool%cplant(npt,leaf) = (avg_annual_cnpp(npt) * avg_af(npt)) / &
                                   avg_lf(npt)
      casapool%cplant(npt,wood) = (avg_annual_cnpp(npt) * avg_aw(npt)) / &
                                   avg_lw(npt)
      casapool%cplant(npt,froot) = (avg_annual_cnpp(npt) * avg_ar(npt)) / &
                                    avg_lr(npt)

    ! compute steady-state litter and soil C pool sizes
     casapool%clitter(npt,metb) = (avgcleaf2met(npt)+avgcroot2met(npt))/casaflux%klitter(npt,metb)
     casapool%clitter(npt,str) = (avgcleaf2str(npt)+avgcroot2str(npt))/casaflux%klitter(npt,str)
     casapool%clitter(npt,cwd) = (avgcwood2cwd(npt))/casaflux%klitter(npt,cwd)
     casapool%csoil(npt,mic)   = (casaflux%fromLtoS(npt,mic,metb)*casaflux%klitter(npt,metb)*casapool%clitter(npt,metb)   &
                                 +casaflux%fromLtoS(npt,mic,str) *casaflux%klitter(npt,str)*casapool%clitter(npt,str)  &
                                 +casaflux%fromLtoS(npt,mic,cwd) *casaflux%klitter(npt,cwd)*casapool%clitter(npt,cwd) ) &
                               /casaflux%ksoil(npt,mic)
      casapool%csoil(npt,slow)  = (casaflux%fromLtoS(npt,slow,metb)*casaflux%klitter(npt,metb)*casapool%clitter(npt,metb) &
                                 + casaflux%fromLtoS(npt,slow,str)*casaflux%klitter(npt,str)*casapool%clitter(npt,str) &
                                 + casaflux%fromLtoS(npt,slow,cwd)*casaflux%klitter(npt,cwd)*casapool%clitter(npt,cwd) &
                                 + casaflux%fromStoS(npt,slow,mic) *casaflux%ksoil(npt,mic) *casapool%csoil(npt,mic)  ) &
                                /casaflux%ksoil(npt,slow)
      casapool%csoil(npt,pass)  = (casaflux%fromStoS(npt,pass,mic) *casaflux%ksoil(npt,mic) *casapool%csoil(npt,mic)    &
                                  +casaflux%fromStoS(npt,pass,slow)*casaflux%ksoil(npt,slow)*casapool%csoil(npt,slow) ) &
                                /casaflux%ksoil(npt,pass)

      casabal%clitterlast = casapool%clitter
      casabal%csoillast   = casapool%csoil

      if(icycle <=1) then
         casapool%nlitter(npt,:)= casapool%rationclitter(npt,:) * casapool%clitter(npt,:)
         casapool%nsoil(npt,:)  = casapool%ratioNCsoil(npt,:)   * casapool%Csoil(npt,:)
         casapool%nsoilmin(npt) = 2.0
         casabal%sumnbal(npt)   = 0.0
       ELSE
       ! compute steady-state litter and soil N pool sizes
         casapool%nlitter(npt,metb) = (avgnleaf2met(npt)+avgnroot2met(npt))/casaflux%klitter(npt,metb)
         casapool%nlitter(npt,str) = (avgnleaf2str(npt)+avgnroot2str(npt))/casaflux%klitter(npt,str)
         casapool%nlitter(npt,cwd) = (avgnwood2cwd(npt))/casaflux%klitter(npt,cwd)

         casapool%nsoil(npt,mic)   = (casaflux%fromLtoS(npt,mic,metb)*casaflux%klitter(npt,metb)*casapool%clitter(npt,metb)   &
                                     +casaflux%fromLtoS(npt,mic,str) *casaflux%klitter(npt,str)*casapool%clitter(npt,str)  &
                                     +casaflux%fromLtoS(npt,mic,cwd) *casaflux%klitter(npt,cwd)*casapool%clitter(npt,cwd) ) &
                                   * avgratioNCsoilmic(npt)/casaflux%ksoil(npt,mic)
         casapool%nsoil(npt,slow)  = (casaflux%fromLtoS(npt,slow,metb)*casaflux%klitter(npt,metb)*casapool%clitter(npt,metb) &
                                     + casaflux%fromLtoS(npt,slow,str)*casaflux%klitter(npt,str)*casapool%clitter(npt,str) &
                                     + casaflux%fromLtoS(npt,slow,cwd)*casaflux%klitter(npt,cwd)*casapool%clitter(npt,cwd) &
                                     + casaflux%fromStoS(npt,slow,mic) *casaflux%ksoil(npt,mic) *casapool%csoil(npt,mic)  ) &
                                   * avgratioNCsoilslow(npt)/casaflux%ksoil(npt,slow)
         casapool%nsoil(npt,pass)  = (casaflux%fromStoS(npt,pass,mic) *casaflux%ksoil(npt,mic) *casapool%csoil(npt,mic)    &
                                     +casaflux%fromStoS(npt,pass,slow)*casaflux%ksoil(npt,slow)*casapool%csoil(npt,slow) ) &
                                   * avgratioNCsoilpass(npt)/casaflux%ksoil(npt,pass)
          casapool%Nsoilmin(npt)    = avgnsoilmin(npt)

        ENDIF

        IF (icycle<=2) THEN
            totpsoil(npt)          = psorder(casamet%isorder(npt)) *xpsoil50(casamet%isorder(npt))
           casapool%plitter(npt,:)= casapool%Nlitter(npt,:)/casapool%ratioNPlitter(npt,:)
            casapool%psoil(npt,:)  = casapool%Nsoil(npt,:)/casapool%ratioNPsoil(npt,:)
            ! why is this commented here but used in UM
            ! casapool%plitter(npt,:)= casapool%ratiopclitter(npt,:)  * casapool%clitter(npt,:)
            ! casapool%psoil(npt,:)  = casapool%ratioPCsoil(npt,:)    * casapool%Csoil(npt,:)
            casapool%psoillab(npt) = totpsoil(npt) *fracpLab(casamet%isorder(npt))
            casapool%psoilsorb(npt)= casaflux%psorbmax(npt) * casapool%psoillab(npt) &
                                    /(casaflux%kmlabp(npt)+casapool%psoillab(npt))
            casapool%psoilocc(npt) = totpsoil(npt) *fracPocc(casamet%isorder(npt))
        ELSE
        ! compute the steady-state litter and soil P pools
          casapool%plitter(npt,metb) = (avgpleaf2met(npt)+avgproot2met(npt))/casaflux%klitter(npt,metb)
          casapool%plitter(npt,str) = (avgpleaf2str(npt)+avgproot2str(npt))/casaflux%klitter(npt,str)
          casapool%plitter(npt,cwd) = (avgpwood2cwd(npt))/casaflux%klitter(npt,cwd)

          casapool%psoil(npt,mic)   = (casaflux%fromLtoS(npt,mic,metb)*casaflux%klitter(npt,metb)*casapool%clitter(npt,metb)   &
                                     +casaflux%fromLtoS(npt,mic,str) *casaflux%klitter(npt,str)*casapool%clitter(npt,str)  &
                                     +casaflux%fromLtoS(npt,mic,cwd) *casaflux%klitter(npt,cwd)*casapool%clitter(npt,cwd) ) &
                                   * (casapool%ratioNCsoil(npt,mic)/casapool%ratioNPsoil(npt,mic))/casaflux%ksoil(npt,mic)
          casapool%psoil(npt,slow)  = (casaflux%fromLtoS(npt,slow,metb)*casaflux%klitter(npt,metb)*casapool%clitter(npt,metb) &
                                     + casaflux%fromLtoS(npt,slow,str)*casaflux%klitter(npt,str)*casapool%clitter(npt,str) &
                                     + casaflux%fromLtoS(npt,slow,cwd)*casaflux%klitter(npt,cwd)*casapool%clitter(npt,cwd) &
                                     + casaflux%fromStoS(npt,slow,mic) *casaflux%ksoil(npt,mic) *casapool%csoil(npt,mic)  ) &
                                   * (casapool%ratioNCsoil(npt,slow)/casapool%ratioNPsoil(npt,slow))/casaflux%ksoil(npt,slow)
          casapool%psoil(npt,pass)  = (casaflux%fromStoS(npt,pass,mic) *casaflux%ksoil(npt,mic) *casapool%csoil(npt,mic)    &
                                     +casaflux%fromStoS(npt,pass,slow)*casaflux%ksoil(npt,slow)*casapool%csoil(npt,slow) ) &
                                   *  (casapool%ratioNCsoil(npt,pass)/casapool%ratioNPsoil(npt,pass))/casaflux%ksoil(npt,pass)
          ! assign the mineral pools
          casapool%psoillab(npt)      = avgpsoillab(npt)
          casapool%psoilsorb(npt)     = avgPsoilsorb(npt)
          casapool%psoilocc(npt)      = avgPsoilocc(npt)
        ENDIF
  ENDIF
  ENDDO

  END SUBROUTINE analyticpool


End module casa_cable
