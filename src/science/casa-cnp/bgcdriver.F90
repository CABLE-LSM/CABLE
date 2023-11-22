!# define UM_BUILD YES

MODULE bgcdriver_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE bgcdriver(ktau,kstart,kend,dels,met,ssnow,canopy,veg,soil, &
                     climate,casabiome,casapool,casaflux,casamet,casabal,phen, &
                     pop, spinConv, spinup, ktauday, idoy,loy, dump_read,   &
                     dump_write, LALLOC)

   USE cable_def_types_mod
   USE cable_common_module, only: cable_runtime
   USE casadimension
   USE casaparm
   USE casavariable
   USE phenvariable
   USE cable_common_module,  ONLY: CurYear, CABLE_USER
   USE TypeDef,              ONLY: i4b, dp
#  ifndef UM_BUILD
   USE POPMODULE,            ONLY: POPStep
#  endif
   USE POP_TYPES,            ONLY: POP_TYPE
#  ifndef UM_BUILD
   USE cable_phenology_module, ONLY: cable_phenology_clim
#  endif
  USE biogeochem_mod, ONLY : biogeochem 
   IMPLICIT NONE

   INTEGER,      INTENT(IN) :: ktau ! integration step number
   INTEGER,      INTENT(IN) :: kstart ! starting value of ktau
   INTEGER,      INTENT(IN) :: kend ! total # timesteps in run

   INTEGER,      INTENT(IN)                  :: idoy ,LOY ! day of year (1-365) , Length oy
   INTEGER,      INTENT(IN)                  :: ktauday
   logical,      INTENT(IN) :: spinConv, spinup
   logical,      INTENT(IN) :: dump_read, dump_write
   INTEGER,      INTENT(IN)                  :: LALLOC

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
   TYPE(POP_TYPE),             INTENT(INOUT) :: POP
   TYPE (climate_type), INTENT(IN)       :: climate  ! climate variables

   ! local variables added ypwang 5/nov/2012
   real,      dimension(mp)  :: cleaf2met, cleaf2str, croot2met, croot2str, cwood2cwd
   real,      dimension(mp)  :: nleaf2met, nleaf2str, nroot2met, nroot2str, nwood2cwd
   real,      dimension(mp)  :: pleaf2met, pleaf2str, proot2met, proot2str, pwood2cwd
   real(r_2), dimension(mp)  :: xnplimit,  xkNlimiting, xklitter, xksoil ,xkleaf,xkleafcold,xkleafdry

   INTEGER                                   :: it, nit
   REAL(dp)                               :: StemNPP(mp,2)
   CHARACTER                                 :: cyear*4
   CHARACTER                                 :: ncfile*99

   INTEGER , parameter :: wlogn=6

IF ( .NOT. dump_read ) THEN  ! construct casa met and flux inputs from current CABLE run

      IF ( TRIM(cable_user%MetType) .EQ. 'cru' ) THEN
         casaflux%Nmindep = met%Ndep
      ENDIF

      IF(ktau == kstart) THEN
         casamet%tairk  = 0.0
         casamet%tsoil  = 0.0
         casamet%moist  = 0.0

      ENDIF

      IF(MOD(ktau,ktauday)==1) THEN
         casamet%tairk = met%tk
         casamet%tsoil = ssnow%tgg
         casamet%moist = ssnow%wb
         casaflux%meangpp = (-canopy%fpn+canopy%frday)*dels
         casaflux%meanrleaf = canopy%frday*dels
      ELSE
         Casamet%tairk  =casamet%tairk + met%tk
         casamet%tsoil = casamet%tsoil + ssnow%tgg
         casamet%moist = casamet%moist + ssnow%wb
         casaflux%meangpp = casaflux%meangpp + (-canopy%fpn+canopy%frday)*dels
         casaflux%meanrleaf = casaflux%meanrleaf + canopy%frday*dels
      ENDIF

      IF(MOD((ktau-kstart+1),ktauday)==0) THEN  ! end of day
  
         casamet%tairk  =casamet%tairk/FLOAT(ktauday)
         casamet%tsoil=casamet%tsoil/FLOAT(ktauday)
         casamet%moist=casamet%moist/FLOAT(ktauday)
         casaflux%cgpp = casaflux%meangpp
         casaflux%crmplant(:,leaf) = casaflux%meanrleaf

         IF ( icycle .GT. 0 ) THEN

#     ifndef UM_BUILD
            IF (trim(cable_user%PHENOLOGY_SWITCH)=='climate') THEN
               ! get climate_dependent phenology
               call cable_phenology_clim(veg, climate, phen)
            ENDIF
#     endif

            CALL biogeochem(ktau,dels,idoy,LALLOC,veg,soil,casabiome,casapool,casaflux, &
                casamet,casabal,phen,POP,climate, xnplimit,xkNlimiting,xklitter,xksoil, &
                xkleaf,xkleafcold,xkleafdry,&
                cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,         &
                nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,         &
                pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd)

            IF (cable_user%CALL_POP) THEN ! accumulate input variables for POP
        
               ! accumulate annual variables for use in POP
               IF(MOD(ktau/ktauday,LOY)==1 ) THEN
          
                  casaflux%stemnpp =  casaflux%cnpp * casaflux%fracCalloc(:,2) * 0.7
                  ! (assumes 70% of wood NPP is allocated above ground)
                  casabal%LAImax = casamet%glai
                  casabal%Cleafmean = casapool%cplant(:,1)/real(LOY)/1000.
                  casabal%Crootmean = casapool%cplant(:,3)/real(LOY)/1000.
        
               ELSE
        
                  casaflux%stemnpp = casaflux%stemnpp + casaflux%cnpp * casaflux%fracCalloc(:,2) * 0.7
                  casabal%LAImax = max(casamet%glai, casabal%LAImax)
                  casabal%Cleafmean = casabal%Cleafmean + casapool%cplant(:,1)/real(LOY)/1000.
                  casabal%Crootmean = casabal%Crootmean + casapool%cplant(:,3)/real(LOY)/1000.
        
               ENDIF
      
            ELSE
      
               casaflux%stemnpp = 0.
      
            ENDIF ! CALL_POP

         ENDIF  ! icycle .gt. 0

      ENDIF  ! end of day

ELSE ! dump_read: ! use casa met and flux inputs from dumpfile

      IF( MOD((ktau-kstart+1),ktauday) == 0 ) THEN  ! end of day

         CALL biogeochem(ktau,dels,idoy,LALLOC,veg,soil,casabiome,casapool,casaflux, &
              casamet,casabal,phen,POP,climate,xnplimit,xkNlimiting,xklitter,xksoil,xkleaf, &
              xkleafcold,xkleafdry,&
              cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,         &
              nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,         &
              pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd)

         IF (cable_user%CALL_POP) THEN ! accumulate input variables for POP

            ! accumulate annual variables for use in POP
            IF(MOD(ktau/ktauday,LOY)==1) THEN
        
               casaflux%stemnpp =  casaflux%cnpp * casaflux%fracCalloc(:,2) * 0.7
               ! (assumes 70% of wood NPP is allocated above ground)
               casabal%LAImax = casamet%glai
               casabal%Cleafmean = casapool%cplant(:,1)/real(LOY)/1000.
               casabal%Crootmean = casapool%cplant(:,3)/real(LOY)/1000.
      
            ELSE
      
               casaflux%stemnpp = casaflux%stemnpp + casaflux%cnpp * casaflux%fracCalloc(:,2) * 0.7
               casabal%LAImax = max(casamet%glai, casabal%LAImax)
               casabal%Cleafmean = casabal%Cleafmean + casapool%cplant(:,1)/real(LOY)/1000.
               casabal%Crootmean = casabal%Crootmean +casapool%cplant(:,3)/real(LOY)/1000.
      
            ENDIF

         ENDIF ! CALL_POP

      ENDIF ! end of day

ENDIF ! dump_read

END SUBROUTINE bgcdriver

END MODULE bgcdriver_mod
