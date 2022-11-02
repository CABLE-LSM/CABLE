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

module bgcdriver_mod

USE casa_inout_mod, ONLY : biogeochem
contains

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
         ! Les 10jan13 - init cnpp ?
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
!         casaflux%cgpp = (-canopy%fpn+canopy%frday)*dels
!         casaflux%crmplant(:,leaf) = canopy%frday*dels
         casaflux%meangpp = (-canopy%fpn+canopy%frday)*dels
         casaflux%meanrleaf = canopy%frday*dels
      ELSE
         Casamet%tairk  =casamet%tairk + met%tk
         casamet%tsoil = casamet%tsoil + ssnow%tgg
         casamet%moist = casamet%moist + ssnow%wb
!         casaflux%cgpp = casaflux%cgpp + (-canopy%fpn+canopy%frday)*dels
!         casaflux%crmplant(:,leaf) = casaflux%crmplant(:,leaf) + canopy%frday*dels
         casaflux%meangpp = casaflux%meangpp + (-canopy%fpn+canopy%frday)*dels
         casaflux%meanrleaf = casaflux%meanrleaf + canopy%frday*dels
      ENDIF

      IF(mod((ktau-kstart+1),ktauday)==0) THEN

         casamet%tairk  =casamet%tairk/FLOAT(ktauday)
         casamet%tsoil=casamet%tsoil/FLOAT(ktauday)
         casamet%moist=casamet%moist/FLOAT(ktauday)

         casaflux%cgpp = casaflux%meangpp
         casaflux%crmplant(:,leaf) = casaflux%meanrleaf
   
         CALL biogeochem(ktau,dels,idoy,veg,soil,casabiome,casapool,casaflux, &
                    casamet,casabal,phen)
   
         !jh!IF((.NOT.spinup).OR.(spinup.AND.spinConv)) THEN 
         !jh!   IF ( dump_write ) &
         !jh!      call ncdf_dump( casamet%tairk, casamet%tsoil, casamet%moist, &
         !jh!                      casaflux%cgpp, casaflux%crmplant, idoy, &
         !jh!                      kend/ktauday )
         !jh!ENDIF

      ENDIF

   ELSE 



      IF( mod((ktau-kstart+1),ktauday) == 0 ) & 
         CALL biogeochem(ktau,dels,idoy,veg,soil,casabiome,casapool,casaflux, &
                    casamet,casabal,phen)


   ENDIF




END SUBROUTINE bgcdriver

End module bgcdriver_mod


 SUBROUTINE casa_feedback(ktau,veg,casabiome,casapool,casamet)
  USE cable_def_types_mod
  USE casadimension
  USE casaparm
  USE casavariable
  IMPLICIT NONE
  INTEGER,      INTENT(IN) :: ktau ! integration step number
  TYPE (veg_parameter_type),  INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (casa_biome),          INTENT(IN) :: casabiome
  TYPE (casa_pool),           INTENT(IN) :: casapool
  TYPE (casa_met),            INTENT(IN) :: casamet
  real, dimension(17)                   ::  nintercept,nslope,xnslope
! Will move to look-up table in later version
  data nintercept/6.32,4.19,6.32,5.73,14.71,6.42,2.00,14.71,4.71,14.71,14.71,7.00,14.71,14.71,14.71,14.71,14.71/
! r29
  data nslope/9.00,28.00,12.00,40.00,10.00,25.00,20.00,4.00,45.00,23.15,23.15,10.00,23.15,23.15,23.15,23.15,23.15/
! Modified further in ACCESS-1.4
  data xnslope/0.80,1.00,2.00,1.00,1.00,1.00,0.50,1.00,0.34,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00/

  integer np,ivt
  real, dimension(mp)  :: ncleafx,npleafx  ! local variables

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
        veg%vcmax(np) = ( nintercept(ivt) &
                        + nslope(ivt)*ncleafx(np)/casabiome%sla(ivt) ) * 1.0e-6
      ELSE
        IF (casapool%nplant(np,leaf)>0.0.AND.casapool%pplant(np,leaf)>0.0) THEN
          veg%vcmax(np) = ( nintercept(ivt)  &
                          + nslope(ivt)*(0.4+9.0/npleafx(np)) &
                          * ncleafx(np)/casabiome%sla(ivt) ) * 1.0e-6
        ELSE
          veg%vcmax(np) = ( nintercept(ivt) &
                          + nslope(ivt)*ncleafx(np)/casabiome%sla(ivt) )*1.0e-6
        ENDIF
      ENDIF
      veg%vcmax(np) = veg%vcmax(np) * xnslope(ivt)
    ENDIF

  ENDDO

  veg%ejmax = 2.0 * veg%vcmax
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



