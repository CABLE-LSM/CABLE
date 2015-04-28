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
! Purpose: bgcdriver - interface between casacnp and cable
!          sumcflux  - accumulating carbon fluxes (not required for UM)
!
! Called from: cable_driver for offline version
!              Not currently called/available for ACCESS version
!
! Contact: Yingping.Wang@csiro.au
!
! History: Model development by Yingping Wang, coupling to Mk3L by Bernard Pak
!          ssoil changed to ssnow
!
! ==============================================================================

SUBROUTINE bgcdriver(ktau,kstart,kend,dels,met,ssnow,canopy,veg,soil, &
                     casabiome,casapool,casaflux,casamet,casabal,phen, &
                     spinConv, spinup, ktauday, idoy, dump_read, dump_write )

   USE cable_def_types_mod
   USE cable_common_module, only: cable_runtime
   USE casadimension
   USE casaparm
   USE casavariable
   USE phenvariable
   IMPLICIT NONE
 
   INTEGER,      INTENT(IN) :: ktau ! integration step number
   INTEGER,      INTENT(IN) :: kstart ! starting value of ktau
   INTEGER,      INTENT(IN) :: kend ! total # timesteps in run
   
   INTEGER,      INTENT(IN)                  :: idoy ! day of year (1-365)
   INTEGER,      INTENT(IN)                  :: ktauday
   logical,      INTENT(IN) :: spinConv, spinup
   logical,      INTENT(IN) :: dump_read, dump_write 
        
   REAL,         INTENT(IN) :: dels ! time setp size (s)
   TYPE (met_type), INTENT(INOUT)       :: met  ! met input variables
   TYPE (soil_snow_type), INTENT(INOUT) :: ssnow ! soil and snow variables
   TYPE (canopy_type), INTENT(INOUT) :: canopy ! vegetation variables
   TYPE (veg_parameter_type),  INTENT(INOUT) :: veg  ! vegetation parameters
   TYPE (soil_parameter_type), INTENT(INOUT) :: soil ! soil parameters  
   TYPE (casa_biome),          INTENT(INOUT) :: casabiome
   TYPE (casa_pool),           INTENT(INOUT) :: casapool
   TYPE (casa_flux),           INTENT(INOUT) :: casaflux
   TYPE (casa_met),            INTENT(INOUT) :: casamet
   TYPE (casa_balance),        INTENT(INOUT) :: casabal
   TYPE (phen_variable),       INTENT(INOUT) :: phen

   !    phen%phase = 2

   if ( .NOT. dump_read ) then
   ! Lest 13may13: will require loop when prog resp are init nonzero
   ! need for mk3l ?
   IF( .NOT. cable_runtime%UM ) THEN
      if(ktau == kstart) then
         casamet%tairk  = 0.0
         casamet%tsoil  = 0.0
         casamet%moist  = 0.0
         casaflux%cgpp  = 0.0
         ! add initializations (BP jul2010)
         !! Les 10jan13 - init cnpp ?
         !casaflux%cnpp  = 0.0
         casaflux%Crsoil   = 0.0
         casaflux%crgplant = 0.0
         casaflux%crmplant = 0.0
         ! Lest 13may13 ---
         casaflux%clabloss = 0.0
         ! casaflux%crmplant(:,leaf) = 0.0
         ! end changes (BP jul2010)
      ENDIF
   ENDIF
      IF(mod(ktau,ktauday)==1) THEN
         casamet%tairk = met%tk
         casamet%tsoil = ssnow%tgg
         casamet%moist = ssnow%wb
         casaflux%cgpp = (-canopy%fpn+canopy%frday)*dels
         casaflux%crmplant(:,leaf) = canopy%frday*dels
      ELSE
         Casamet%tairk  =casamet%tairk + met%tk
         casamet%tsoil = casamet%tsoil + ssnow%tgg
         casamet%moist = casamet%moist + ssnow%wb
         casaflux%cgpp = casaflux%cgpp + (-canopy%fpn+canopy%frday)*dels
         casaflux%crmplant(:,leaf) = casaflux%crmplant(:,leaf) + &
                                       canopy%frday*dels
      ENDIF

      IF(mod((ktau-kstart+1),ktauday)==0) THEN

         casamet%tairk  =casamet%tairk/FLOAT(ktauday)
         casamet%tsoil=casamet%tsoil/FLOAT(ktauday)
         casamet%moist=casamet%moist/FLOAT(ktauday)
   
         CALL biogeochem(ktau,dels,idoy,veg,soil,casabiome,casapool,casaflux, &
                    casamet,casabal,phen)
   
         IF((.NOT.spinup).OR.(spinup.AND.spinConv)) THEN 
            IF ( dump_write ) &
               call ncdf_dump( casamet%tairk, casamet%tsoil, casamet%moist, &
                               casaflux%cgpp, casaflux%crmplant, idoy, &
                               kend/ktauday )
         ENDIF

      ENDIF

   ELSE 



      IF( mod((ktau-kstart+1),ktauday) == 0 ) & 
         CALL biogeochem(ktau,dels,idoy,veg,soil,casabiome,casapool,casaflux, &
                    casamet,casabal,phen)


   ENDIF




END SUBROUTINE bgcdriver






!! DOES THIS NEED TO BE DELETED FOR NOW - REPLACED WITH BP CODE (LATER?)

   subroutine ncdf_dump(tairk, tsoil, moist, &
                        cgpp, crmplant, &
                        n_call, kend)
!      use netcdf
!      use cable_common_module, only : kend_gl
!      use cable_diag_module, only : def_dims, def_vars, def_var_atts, & 
!                                    put_var_nc, stderr_nc
!
!      implicit none  
!      !var (type) to write 
!      real(r_2), dimension(mp), intent(in) :: & 
!         tairk, &
!         cgpp, &
!         crmplant
!
!      real(r_2), dimension(mp,ms), intent(in) :: & 
!         tsoil, &
!         moist
!      
!      integer, intent(in) :: &
!         n_call, &         ! this timestep # 
!         kend              ! final timestep of run
!
!     
!      !number of instances. dummied here and so=1 
!      !integer :: inst =1
!
!      !netcdf IDs/ names 
!      character(len=*), parameter :: ncfile = "CASA_dump.nc"
!      integer, parameter :: num_vars=5 
!      integer, parameter :: num_dims=4 
!      integer, save :: ncid       ! netcdf file ID
!      
!      !vars 
!      character(len=*), dimension(num_vars), parameter :: &
!            var_name =  (/  "casamet_tairk", & 
!                            "tsoil        ", &
!                            "moist        ", &
!                            "cgpp         ", &
!                            "crmplant     " /)
!
!      integer, dimension(num_vars) :: varID ! (1) tvair, (2) pmb 
!      
!      !dims 
!      character(len=*), dimension(num_dims), parameter :: & 
!            dim_name =  (/ "lat ", &
!                           "lon ", &
!                           "soil", &
!                           "time" /)
!      
!      integer, parameter :: soil_dim = 6
!            
!      integer, dimension(soil_dim), parameter  :: soil = (/ 1,2,3,4,5,6 /)
!      
!      integer, dimension(num_dims)  :: &
!            dimID   ! (1) x, (2) y, (3) time
!      
!      integer, dimension(num_dims)  :: &
!            !x,y generally lat/lon BUT for single site = 1,1       
!            dim_len = (/1,1,soil_dim,-1/)  ! (1) x, (2) y, (3) soil, (4) time [re-set] 
!      
!      
!      !local only
!      integer :: ncok      !ncdf return status
!      
!      ! END header
!
!      dim_len(num_dims) = kend
!
!      if (n_call == 1) then
!         ! create netCDF dataset: enter define mode
!         ncok = nf90_create(path = ncfile, cmode = nf90_noclobber, ncid = ncid)
!            if (ncok /= nf90_noerr) call stderr_nc('ncdf creating ', ncfile) 
!      
!            ! define dimensions: from name and length
!            call def_dims(num_dims, ncid, dimID, dim_len, dim_name )
!     
!            ! define variables: from name, type, dims
!            call def_vars(num_vars, ncid,  nf90_float, dimID, var_name, varID )
!      
!            ! define variable attributes
!            call def_var_atts(ncfile, ncid, varID )
!               
!            ncok = nf90_enddef(ncid) 
!         
!      endif 
!      
!      call put_var_nc(ncid, var_name(1), tairk, n_call )
!      call put_var_nc(ncid, var_name(2), tsoil, n_call )
!      call put_var_nc(ncid, var_name(3), moist, n_call )
!      call put_var_nc(ncid, var_name(4), cgpp, n_call )
!      call put_var_nc(ncid, var_name(5), crmplant, n_call )
!      
!      if (n_call == kend ) & 
!         ncok = nf90_close(ncid)            ! close: save new netCDF dataset
     
   end subroutine ncdf_dump


!! DOES THIS NEED TO BE DELETED FOR NOW - REPLACED WITH BP CODE (LATER?)
   subroutine read_casa_dump( casamet, casaflux, ktau, kend )
!      use netcdf
!      USE casa_cnp_module  
!      use cable_diag_module, only : get_var_nc, stderr_nc
!
!      TYPE (casa_flux), intent(inout) :: casaflux
!      TYPE (casa_met), intent(inout)  :: casamet
!      integer, intent(in) :: kend, ktau 
!
!
!      !netcdf IDs/ names 
!      character(len=*), parameter :: ncfile = "CASA_dump.nc"
!      integer, parameter :: num_vars=5
!      integer, parameter :: num_dims=4
!      integer:: ncid       ! netcdf file ID
! 
!      !vars 
!      character(len=*), dimension(num_vars), parameter :: &
!            var_name =  (/  "casamet_tairk", & 
!                            "tsoil        ", &
!                            "moist        ", &
!                            "cgpp         ", &
!                            "crmplant     " /)
!
!      integer, dimension(num_vars) :: varID ! (1) tvair, (2) pmb 
!      
!      ncok = NF90_OPEN(ncfile, nf90_nowrite, ncid)           
!         if (ncok /= nf90_noerr ) call stderr_nc('re-opening ', ncfile)      
!
!      call get_var_nc(ncid, var_name(1), casamet%tairk,ktau, kend )
!      call get_var_nc(ncid, var_name(2), casamet%tsoil,ktau, kend )
!      call get_var_nc(ncid, var_name(3), casamet%moist,ktau, kend )
!      call get_var_nc(ncid, var_name(4), casaflux%cgpp,ktau, kend )
!      call get_var_nc(ncid, var_name(5), casaflux%crmplant, ktau, kend )
!      
!      ncok = NF90_CLOSE(ncid)            
!         if (ncok /= nf90_noerr ) call stderr_nc('closing ', ncfile)      
!      
   end subroutine read_casa_dump



 SUBROUTINE casa_feedback(ktau,veg,casabiome,casapool,casamet)
  USE cable_def_types_mod
  USE casadimension
  USE casaparm
  USE casavariable
  IMPLICIT NONE
  INTEGER,      INTENT(IN) :: ktau ! integration step number
  TYPE (veg_parameter_type),  INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (casa_biome),          INTENT(INOUT) :: casabiome
  TYPE (casa_pool),           INTENT(IN) :: casapool
  TYPE (casa_met),            INTENT(IN) :: casamet

  integer np,ivt
  real, dimension(mp)  :: ncleafx,npleafx  ! local variables
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
        npleafx(np) = MIN(30.0,MAX(8.0,casapool%nplant(np,leaf) &
                                      /casapool%pplant(np,leaf)))
      ENDIF
    ENDIF

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

!    veg%vcmax(np) = ( nintercept(ivt)  &
!                  + nslope(ivt)*(0.4+8.5/npleafx(np)) &
!                  * ncleafx(np)/casabiome%sla(ivt))*(1.0e-6)
!    veg%vcmax(np) =veg%vcmax(np)* xnslope(ivt)

!write(*,991) np, ivt,veg%vlai(np),veg%vcmax(np)*1.0e6
!write(*,891) np,ivt,casapool%cplant(np,leaf),casapool%nplant(np,leaf),casapool%pplant(np,leaf)
!891 format(2(i6),3(f9.3,2x))
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

!    canopy%fnee = canopy%fpn + canopy%frs + canopy%frp 
!!   For prognostic Vcmax, and NEE should include clabloss under nutrient limitation
!!   Q.Zhang 12/09/2011
!    canopy%fnee(:) = canopy%fnee(:) + casaflux%clabloss(:)/86400.0
!    ! Q.Zhang 08/06/2011. return NEE from casaflux when N/NP mode is activated.
!    ! NPP of CABLE's output is "potential" NPP, not "real" C input to casacnp
!    ! To derive nutrient limited NPP from CABLE's standard output, use NEE+Crsoil
!!    if (icycle>1) then
!!     canopy%fnee = (casaflux%Crsoil-casaflux%cnpp)/86400.0
!!    else
!!     canopy%fnee = canopy%fpn + canopy%frs + canopy%frp
!!    end if

!    write(*,101) ktau,casaflux%Crsoil(:)
!    write(*,101) ktau,dels,veg%vlai,veg%vcmax*1.0e6,casaflux%cgpp,canopy%fpn*1.0e6/12.0,canopy%frp*1.0e6/12.0,canopy%frs*1.0e6/12.0,canopy%fnee*1.0e6/12.0
!101  format(i6,2x,100(f12.5,2x))
!    if(ktau==kend) then
!       PRINT *, 'carbon fluxes'
!       PRINT *, 'sumpn', sum_flux%sumpn
!       PRINT *, 'sumrd', sum_flux%sumrd
!       PRINT *, 'sumrp', sum_flux%sumrp
!       PRINT *, 'sumrs', sum_flux%sumrs
!       PRINT *, 'npp', sum_flux%sumpn+sum_flux%sumrp
!       PRINT *, 'nee', sum_flux%sumpn+sum_flux%sumrp+sum_flux%sumrs
!     !  PRINT *, 'carbon pools', leaf,wood,froot
!     !  PRINT *,  casapool%cplant(1,2),casaflux%crmplant(1,wood),casaflux%Crsoil(1)
!     !  PRINT *, 'respiration rate'
!     !  PRINT *,  casabiome%rmplant(1,2)*365.0
!    endif

END SUBROUTINE sumcflux



