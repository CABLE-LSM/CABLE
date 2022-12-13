
!#define UM_BUILD YES

MODULE biogeochem_mod

IMPLICIT NONE

CONTAINS

  SUBROUTINE biogeochem(ktau,dels,idoY,LALLOC,veg,soil,casabiome,casapool,casaflux, &
       casamet,casabal,phen,POP,climate,xnplimit,xkNlimiting,xklitter,xksoil,xkleaf,xkleafcold,xkleafdry,&
       cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,         &
       nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,         &
       pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd)
USE cable_def_types_mod
USE cable_common_module, ONLY : cable_runtime, cable_user
USE casadimension
USE casa_cnp_module
USE casa_inout_module, ONLY : casa_cnpflux
USE POP_TYPES,            ONLY: POP_TYPE
#ifdef UM_BUILD 
USE casa_rplant_module, ONLY: casa_rplant 
#else
USE casa_rplant_module, ONLY: casa_rplant1 
#endif

IMPLICIT NONE
INTEGER, INTENT(IN)    :: ktau
REAL,    INTENT(IN)    :: dels
INTEGER, INTENT(IN)    :: idoy
INTEGER, INTENT(IN)    :: LALLOC
TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
TYPE (soil_parameter_type),   INTENT(INOUT) :: soil ! soil parameters
TYPE (casa_biome),            INTENT(INOUT) :: casabiome
TYPE (casa_pool),             INTENT(INOUT) :: casapool
TYPE (casa_flux),             INTENT(INOUT) :: casaflux
TYPE (casa_met),              INTENT(INOUT) :: casamet
TYPE (casa_balance),          INTENT(INOUT) :: casabal
TYPE (phen_variable),         INTENT(INOUT) :: phen
TYPE(POP_TYPE),             INTENT(IN) :: POP
TYPE(climate_TYPE),             INTENT(IN) :: climate

REAL, DIMENSION(mp), INTENT(OUT)   :: cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,         &
     nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,         &
     pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd

! local variables
REAL(r_2),    DIMENSION(mp) :: xnplimit,xNPuptake
REAL(r_2),    DIMENSION(mp) :: xklitter,xksoil,xkNlimiting
REAL(r_2),    DIMENSION(mp) :: xkleafcold,xkleafdry,xkleaf
INTEGER  npt,j
REAL, ALLOCATABLE :: tmp(:)
real latx,lonx
integer ivegx,isox


xKNlimiting = 1.0

! zero annual sums
IF (idoy==1) CALL casa_cnpflux(casaflux,casapool,casabal,.TRUE.)

IF (cable_user%PHENOLOGY_SWITCH.EQ.'MODIS' .OR. cable_runtime%esm15 ) THEN
       CALL phenology(idoy,veg,phen)
ENDIF
CALL avgsoil(veg,soil,casamet)
#ifdef UM_BUILD 
CALL casa_rplant(veg,casabiome,casapool,casaflux,casamet,climate)
#else
CALL casa_rplant1(veg,casabiome,casapool,casaflux,casamet)
#endif
IF (.NOT.cable_user%CALL_POP) THEN
   CALL casa_allocation(veg,soil,casabiome,casaflux,casapool,casamet,phen,LALLOC)
ENDIF

CALL casa_xrateplant(xkleafcold,xkleafdry,xkleaf,veg,casabiome, &
         casamet,phen)
CALL casa_coeffplant(xkleafcold,xkleafdry,xkleaf,veg,casabiome,casapool, &
         casaflux,casamet,phen)

CALL casa_xnp(xnplimit,xNPuptake,veg,casabiome,casapool,casaflux,casamet)

IF (cable_user%CALL_POP) THEN

       CALL casa_allocation(veg,soil,casabiome,casaflux,casapool,casamet,phen,LALLOC)
  
       WHERE (pop%pop_grid(:)%cmass_sum_old.GT.0.001 .AND. pop%pop_grid(:)%cmass_sum.GT.0.001 )

          casaflux%frac_sapwood(POP%Iwood) = POP%pop_grid(:)%csapwood_sum/ POP%pop_grid(:)%cmass_sum
          casaflux%sapwood_area(POP%Iwood) = MAX(POP%pop_grid(:)%sapwood_area/10000., 1e-6)
          veg%hc(POP%Iwood) = POP%pop_grid(:)%height_max

          WHERE (pop%pop_grid(:)%LU ==2)

             casaflux%kplant(POP%Iwood,2) =  1.0 -  &
                  (1.0-  MAX( MIN((POP%pop_grid(:)%stress_mortality + &
                  POP%pop_grid(:)%crowding_mortality+ &
                  + POP%pop_grid(:)%fire_mortality ) &
                  /(POP%pop_grid(:)%cmass_sum+POP%pop_grid(:)%growth) + &
                  1.0/veg%disturbance_interval(POP%Iwood,1), 0.99), 0.0))**(1.0/365.0)

          ELSEWHERE
             casaflux%kplant(POP%Iwood,2) =  1.0 -  &
                  (1.0-  MAX( MIN((POP%pop_grid(:)%stress_mortality + &
                  POP%pop_grid(:)%crowding_mortality+ &
                  + POP%pop_grid(:)%fire_mortality+POP%pop_grid(:)%cat_mortality  ) &
                  /(POP%pop_grid(:)%cmass_sum+POP%pop_grid(:)%growth), 0.99), 0.0))**(1.0/365.0)

          ENDWHERE

          veg%hc(POP%Iwood) = POP%pop_grid(:)%height_max
  
       ELSEWHERE
    
          casaflux%frac_sapwood(POP%Iwood) = 1.0
          casaflux%sapwood_area(POP%Iwood) = MAX(POP%pop_grid(:)%sapwood_area/10000., 1e-6)
          casaflux%kplant(POP%Iwood,2) = 0.0
          veg%hc(POP%Iwood) = POP%pop_grid(:)%height_max
  
       ENDWHERE

ENDIF

CALL casa_xratesoil(xklitter,xksoil,veg,soil,casamet,casabiome)

CALL casa_coeffsoil(xklitter,xksoil,veg,soil,casabiome,casaflux,casamet)

IF (icycle>1) THEN
       CALL casa_xkN(xkNlimiting,casapool,casaflux,casamet,casabiome,veg)
       DO j=1,mlitter
          casaflux%klitter(:,j) = casaflux%klitter(:,j)* xkNlimiting(:)
       ENDDO
       CALL casa_nuptake(veg,xkNlimiting,casabiome,casapool,casaflux,casamet)
       IF (icycle >2) CALL casa_puptake(veg,xkNlimiting,casabiome, &
            casapool,casaflux,casamet)
ENDIF

CALL casa_delplant(veg,casabiome,casapool,casaflux,casamet,                &
                       cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,  &
                       nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,  &
                       pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd)

casaflux%Cplant_turnover_disturbance = 0
casaflux%Cplant_turnover_crowding = 0
casaflux%Cplant_turnover_resource_limitation = 0

IF (cable_user%CALL_POP) THEN
       IF (.NOT.ALLOCATED(tmp)) ALLOCATE(tmp(SIZE(POP%pop_grid)))
       tmp = (POP%pop_grid(:)%stress_mortality + POP%pop_grid(:)%crowding_mortality &
            +POP%pop_grid(:)%cat_mortality &
            + POP%pop_grid(:)%fire_mortality  )
       WHERE (tmp.GT. 1.e-12)
          casaflux%Cplant_turnover_disturbance(POP%Iwood) =  &
               casaflux%Cplant_turnover(POP%Iwood,2)*(POP%pop_grid(:)%cat_mortality &
               + POP%pop_grid(:)%fire_mortality  )/tmp
          casaflux%Cplant_turnover_crowding(POP%Iwood) =  &
               casaflux%Cplant_turnover(POP%Iwood,2)*POP%pop_grid(:)%crowding_mortality/tmp
          casaflux%Cplant_turnover_resource_limitation(POP%Iwood) = &
               casaflux%Cplant_turnover(POP%Iwood,2)*POP%pop_grid(:)%stress_mortality/tmp
       endwhere
ENDIF

CALL casa_delsoil(veg,casapool,casaflux,casamet,casabiome)

! add "casabal" to store pool sizes before updating
CALL casa_cnpcycle(veg,casabiome,casapool,casaflux,casamet,casabal, LALLOC)

!! vh_js !!
!CLN ndummy must be before pdummy!!!!
IF (icycle<3) THEN
    IF (icycle<2) CALL casa_ndummy(casamet,casabal,casapool)
    CALL casa_pdummy(casamet,casabal,casaflux,casapool)
ENDIF

CALL casa_cnpbal(veg,casamet,casapool,casaflux,casabal)

CALL casa_cnpflux(casaflux,casapool,casabal,.FALSE.)

! Limit labile for spinning up only
IF(cable_user%l_limit_labile .AND. icycle > 1) THEN
  casapool%Nsoilmin = max(casapool%Nsoilmin,0.5)
  casapool%Psoillab = max(casapool%Psoillab,0.1)
ENDIF

END SUBROUTINE biogeochem

END MODULE biogeochem_mod
