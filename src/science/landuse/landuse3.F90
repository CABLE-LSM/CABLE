MODULE landuse_variable
  use landuse_constant
  IMPLICIT NONE

  SAVE

  TYPE landuse_mland
    ! patch generic
    INTEGER,   DIMENSION(:,:),       ALLOCATABLE :: iveg_x
    INTEGER,   DIMENSION(:,:),       ALLOCATABLE :: isoil_x
    INTEGER,   DIMENSION(:,:),       ALLOCATABLE :: soilorder_x
    INTEGER,   DIMENSION(:,:),       ALLOCATABLE :: phase_x
    real(r_2), DIMENSION(:,:),       ALLOCATABLE :: phen_x
    real(r_2), DIMENSION(:,:),       ALLOCATABLE :: aphen_x
    integer,   DIMENSION(:,:),       ALLOCATABLE :: doyphase3_x
    real(r_2), DIMENSION(:,:),       ALLOCATABLE :: frac_sapwood_x
    real(r_2), DIMENSION(:,:),       ALLOCATABLE :: sapwood_area_x
    INTEGER,   DIMENSION(:,:),       ALLOCATABLE :: isflag_x
    REAL(r_2), DIMENSION(:,:),       ALLOCATABLE :: patchfrac_x
    REAL(r_2), DIMENSION(:,:),       ALLOCATABLE :: lai_x
    REAL(r_2), DIMENSION(:,:),       ALLOCATABLE :: sla_x
    ! biophysical
    real(r_2), dimension(:,:,:), allocatable :: albsoilsn_x
    real(r_2), dimension(:,:,:), allocatable :: albedo_x
    real(r_2), dimension(:,:,:), allocatable :: albsoil_x
    real(r_2),dimension(:,:),   allocatable :: dgdtg_x
    real(r_2),dimension(:,:,:), allocatable :: gammzz_x
    real(r_2), dimension(:,:,:), allocatable :: tgg_x
    real(r_2), dimension(:,:,:), allocatable :: wb_x
    real(r_2), dimension(:,:,:), allocatable :: wbice_x
    real(r_2), dimension(:,:,:), allocatable :: tggsn_x
    real(r_2), dimension(:,:,:), allocatable :: ssdn_x
    real(r_2), dimension(:,:,:), allocatable :: smass_x
    real(r_2), dimension(:,:,:), allocatable :: sdepth_x
    real(r_2), dimension(:,:),   allocatable :: tss_x
    real(r_2), dimension(:,:),   allocatable :: rtsoil_x
    real(r_2), dimension(:,:),   allocatable :: runoff_x
    real(r_2), dimension(:,:),   allocatable :: rnof1_x
    real(r_2), dimension(:,:),   allocatable :: rnof2_x
    real(r_2), dimension(:,:),   allocatable :: ssdnn_x
    real(r_2), dimension(:,:),   allocatable :: snowd_x
    real(r_2), dimension(:,:),   allocatable :: snage_x
    real(r_2), dimension(:,:),   allocatable :: osnowd_x

    real(r_2), dimension(:,:),   allocatable :: cansto_x
    real(r_2), dimension(:,:),   allocatable :: ghflux_x
    real(r_2), dimension(:,:),   allocatable :: sghflux_x
    real(r_2), dimension(:,:),   allocatable :: ga_x
    real(r_2), dimension(:,:),   allocatable :: fev_x
    real(r_2), dimension(:,:),   allocatable :: fes_x
    real(r_2), dimension(:,:),   allocatable :: fhs_x
    real(r_2), dimension(:,:),   allocatable :: wbtot0_x
    real(r_2), dimension(:,:),   allocatable :: osnowd0_x
    real(r_2), dimension(:,:),   allocatable :: trad_x
    real(r_2), dimension(:,:),   allocatable :: GWwb_x
    real(r_2), dimension(:,:,:), allocatable :: cplantx_x
    real(r_2), dimension(:,:,:), allocatable :: csoilx_x

    REAL(r_2), DIMENSION(:,:),   ALLOCATABLE :: clabile_x
    REAL(r_2), DIMENSION(:,:,:), ALLOCATABLE :: cplant_x
    REAL(r_2), DIMENSION(:,:,:), ALLOCATABLE :: clitter_x
    REAL(r_2), DIMENSION(:,:,:), ALLOCATABLE :: csoil_x
    REAL(r_2), DIMENSION(:,:,:), ALLOCATABLE :: cwoodprod_x
    REAL(r_2), DIMENSION(:,:,:), ALLOCATABLE :: nplant_x
    REAL(r_2), DIMENSION(:,:,:), ALLOCATABLE :: nlitter_x
    REAL(r_2), DIMENSION(:,:,:), ALLOCATABLE :: nsoil_x
    REAL(r_2), DIMENSION(:,:),   ALLOCATABLE :: nsoilmin_x
    REAL(r_2), DIMENSION(:,:,:), ALLOCATABLE :: nwoodprod_x
    REAL(r_2), DIMENSION(:,:,:), ALLOCATABLE :: pplant_x
    REAL(r_2), DIMENSION(:,:,:), ALLOCATABLE :: plitter_x
    REAL(r_2), DIMENSION(:,:,:), ALLOCATABLE :: psoil_x
    REAL(r_2), DIMENSION(:,:),   ALLOCATABLE :: psoillab_x
    REAL(r_2), DIMENSION(:,:),   ALLOCATABLE :: psoilsorb_x
    REAL(r_2), DIMENSION(:,:),   ALLOCATABLE :: psoilocc_x
    REAL(r_2), DIMENSION(:,:,:), ALLOCATABLE :: pwoodprod_x

    INTEGER,   DIMENSION(:,:),       ALLOCATABLE :: iveg_y
    INTEGER,   DIMENSION(:,:),       ALLOCATABLE :: isoil_y
    INTEGER,   DIMENSION(:,:),       ALLOCATABLE :: soilorder_y
    INTEGER,   DIMENSION(:,:),       ALLOCATABLE :: phase_y
    real(r_2), DIMENSION(:,:),       ALLOCATABLE :: phen_y
    real(r_2), DIMENSION(:,:),       ALLOCATABLE :: aphen_y
    integer,   DIMENSION(:,:),       ALLOCATABLE :: doyphase3_y
    real(r_2), DIMENSION(:,:),       ALLOCATABLE :: frac_sapwood_y
    real(r_2), DIMENSION(:,:),       ALLOCATABLE :: sapwood_area_y
    INTEGER,   DIMENSION(:,:),       ALLOCATABLE :: isflag_y
    REAL(r_2), DIMENSION(:,:),   ALLOCATABLE :: patchfrac_y
    REAL(r_2), DIMENSION(:,:),   ALLOCATABLE :: lai_y
    REAL(r_2), DIMENSION(:,:),   ALLOCATABLE :: sla_y
    ! biophysical
    real(r_2), dimension(:,:,:), allocatable :: albsoilsn_y
    real(r_2), dimension(:,:,:), allocatable :: albedo_y
    real(r_2), dimension(:,:,:), allocatable :: albsoil_y
    real(r_2),dimension(:,:),   allocatable :: dgdtg_y
    real(r_2),dimension(:,:,:), allocatable :: gammzz_y
    real(r_2), dimension(:,:,:), allocatable :: tgg_y
    real(r_2), dimension(:,:,:), allocatable :: wb_y
    real(r_2), dimension(:,:,:), allocatable :: wbice_y
    real(r_2), dimension(:,:,:), allocatable :: tggsn_y
    real(r_2), dimension(:,:,:), allocatable :: ssdn_y
    real(r_2), dimension(:,:,:), allocatable :: smass_y
    real(r_2), dimension(:,:,:), allocatable :: sdepth_y
    real(r_2), dimension(:,:),   allocatable :: tss_y
    real(r_2), dimension(:,:),   allocatable :: rtsoil_y
    real(r_2), dimension(:,:),   allocatable :: runoff_y
    real(r_2), dimension(:,:),   allocatable :: rnof1_y
    real(r_2), dimension(:,:),   allocatable :: rnof2_y
    real(r_2), dimension(:,:),   allocatable :: ssdnn_y
    real(r_2), dimension(:,:),   allocatable :: snowd_y
    real(r_2), dimension(:,:),   allocatable :: snage_y
    real(r_2), dimension(:,:),   allocatable :: osnowd_y

    real(r_2), dimension(:,:),   allocatable :: cansto_y
    real(r_2), dimension(:,:),   allocatable :: ghflux_y
    real(r_2), dimension(:,:),   allocatable :: sghflux_y
    real(r_2), dimension(:,:),   allocatable :: ga_y
    real(r_2), dimension(:,:),   allocatable :: fev_y
    real(r_2), dimension(:,:),   allocatable :: fes_y
    real(r_2), dimension(:,:),   allocatable :: fhs_y
    real(r_2), dimension(:,:),   allocatable :: wbtot0_y
    real(r_2), dimension(:,:),   allocatable :: osnowd0_y
    real(r_2), dimension(:,:),   allocatable :: trad_y
    real(r_2), dimension(:,:),   allocatable :: GWwb_y
    real(r_2), dimension(:,:,:), allocatable :: cplantx_y
    real(r_2), dimension(:,:,:), allocatable :: csoilx_y

    REAL(r_2), DIMENSION(:,:),   ALLOCATABLE :: clabile_y
    REAL(r_2), DIMENSION(:,:,:), ALLOCATABLE :: cplant_y
    REAL(r_2), DIMENSION(:,:,:), ALLOCATABLE :: clitter_y
    REAL(r_2), DIMENSION(:,:,:), ALLOCATABLE :: csoil_y
    REAL(r_2), DIMENSION(:,:,:), ALLOCATABLE :: cwoodprod_y
    REAL(r_2), DIMENSION(:,:,:), ALLOCATABLE :: nplant_y
    REAL(r_2), DIMENSION(:,:,:), ALLOCATABLE :: nlitter_y
    REAL(r_2), DIMENSION(:,:,:), ALLOCATABLE :: nsoil_y 
    REAL(r_2), DIMENSION(:,:),   ALLOCATABLE :: nsoilmin_y
    REAL(r_2), DIMENSION(:,:,:), ALLOCATABLE :: nwoodprod_y
    REAL(r_2), DIMENSION(:,:,:), ALLOCATABLE :: pplant_y
    REAL(r_2), DIMENSION(:,:,:), ALLOCATABLE :: plitter_y
    REAL(r_2), DIMENSION(:,:,:), ALLOCATABLE :: psoil_y
    REAL(r_2), DIMENSION(:,:),   ALLOCATABLE :: psoillab_y
    REAL(r_2), DIMENSION(:,:),   ALLOCATABLE :: psoilsorb_y
    REAL(r_2), DIMENSION(:,:),   ALLOCATABLE :: psoilocc_y
    REAL(r_2), DIMENSION(:,:,:), ALLOCATABLE :: pwoodprod_y

  ! landuse data
    REAL(r_2), DIMENSION(:,:),   ALLOCATABLE :: pftfrac
    REAL(r_2), DIMENSION(:,:),   ALLOCATABLE :: fharvw 
    REAL(r_2), DIMENSION(:,:,:), ALLOCATABLE :: xluh2cable 
    REAL(r_2), DIMENSION(:,:,:), ALLOCATABLE :: atransit
  END TYPE landuse_mland

  TYPE landuse_mp

   ! generic patch properties
   integer,    dimension(:),        allocatable   :: iveg,isoil,soilorder,phase,isflag           !(mp)
   integer,    dimension(:),        allocatable   :: doyphase3                                   !(mp)
   real(r_2),  dimension(:),        allocatable   :: lat,lon                                     !(mp)
   real(r_2),  dimension(:),        allocatable   :: phen,aphen,frac_sapwood, sapwood_area       !(mp)
   real(r_2),  dimension(:),        allocatable   :: patchfrac,areacell,lai,sla               !float(mp)

   ! biophysical variables
   real(r_2),  dimension(:,:),      allocatable   :: albsoilsn,albedo,albsoil      !float(mp,rad)
   real(r_2),  dimension(:),        allocatable   :: dgdtg                         !double(mp)
   real(r_2),  dimension(:,:),      allocatable   :: gammzz                        !double(mp,soil)
   real(r_2),  dimension(:,:),      allocatable   :: tgg,wb,wbice                  !float(mp,soil)
   real(r_2),  dimension(:,:),      allocatable   :: tggsn,ssdn,smass,sdepth       !float(mp,snow)
   real(r_2),  dimension(:),        allocatable   :: tss,rtsoil,runoff,rnof1,rnof2, &
                                                         ssdnn,snowd,snage,osnowd,      &
                                                         cansto,ghflux,sghflux,ga,      &
                                                         fev,fes,fhs,wbtot0,osnowd0,    &
                                                         trad,GWwb                     !float(mp)
   real(r_2),  dimension(:,:),      allocatable   :: cplantx, csoilx               !float(mp,plant_carbon_pools/soil_carbon_pools)

   ! biogeochemical variables
   real(r_2),  dimension(:),        allocatable   :: sumcbal,sumnbal,sumpbal                              !float(mp)                   
   real(r_2),  dimension(:),        allocatable   :: clabile,nsoilmin,psoillab,psoilsorb,psoilocc         !float(mp)
   real(r_2),  dimension(:,:),      allocatable   :: cplant,nplant,pplant                                 !float(mp,mplant)
   real(r_2),  dimension(:,:),      allocatable   :: clitter,nlitter,plitter                              !float(mp,mlitter)
   real(r_2),  dimension(:,:),      allocatable   :: csoil,nsoil,psoil                                    !float(mp,msoil)
   real(r_2),  dimension(:,:),      allocatable   :: cwoodprod,nwoodprod,pwoodprod                        !float(mp,mwood)

 END TYPE landuse_mp
 
 CONTAINS

  SUBROUTINE landuse_allocate_mland(mland,luc)
   use landuse_constant
   IMPLICIT NONE
   TYPE(landuse_mland), INTENT(INOUT)  :: luc
   integer  mland
   ! patch-genric variables
   ALLOCATE(luc%iveg_x(mland,mvmax),             &
            luc%isoil_x(mland,mvmax),            &
            luc%soilorder_x(mland,mvmax),        &
            luc%phase_x(mland,mvmax),            &
            luc%phen_x(mland,mvmax),             &
            luc%aphen_x(mland,mvmax),            &
            luc%doyphase3_x(mland,mvmax),        &
            luc%frac_sapwood_x(mland,mvmax),     &
            luc%sapwood_area_x(mland,mvmax),     &
            luc%isflag_x(mland,mvmax),           &
            luc%patchfrac_x(mland,mvmax),        &
            luc%lai_x(mland,mvmax),              &
            luc%sla_x(mland,mvmax))

   ALLOCATE(luc%iveg_y(mland,mvmax),             &
            luc%isoil_y(mland,mvmax),            &
            luc%soilorder_y(mland,mvmax),        &
            luc%phase_y(mland,mvmax),            &
            luc%phen_y(mland,mvmax),             &
            luc%aphen_y(mland,mvmax),            &
            luc%doyphase3_y(mland,mvmax),        &
            luc%frac_sapwood_y(mland,mvmax),     &
            luc%sapwood_area_y(mland,mvmax),     &
            luc%isflag_y(mland,mvmax),           &
            luc%patchfrac_y(mland,mvmax),        &
            luc%lai_y(mland,mvmax),              &
            luc%sla_y(mland,mvmax))

    ! biophysical
   ALLOCATE(luc%albsoilsn_x(mland,mvmax,nrb),   &
            luc%albedo_x(mland,mvmax,nrb),      &
            luc%albsoil_x(mland,mvmax,nrb),     &
            luc%dgdtg_x(mland,mvmax),           &
            luc%gammzz_x(mland,mvmax,ms),       &
            luc%tgg_x(mland,mvmax,ms),          &
            luc%wb_x(mland,mvmax,ms),           &
            luc%wbice_x(mland,mvmax,ms),        &
            luc%tggsn_x(mland,mvmax,msn),       &
            luc%ssdn_x(mland,mvmax,msn),        &
            luc%smass_x(mland,mvmax,msn),       &
            luc%sdepth_x(mland,mvmax,msn),      &
            luc%tss_x(mland,mvmax),             &
            luc%rtsoil_x(mland,mvmax),          &
            luc%runoff_x(mland,mvmax),          &
            luc%rnof1_x(mland,mvmax),           &
            luc%rnof2_x(mland,mvmax),           &
            luc%ssdnn_x(mland,mvmax),           &
            luc%snowd_x(mland,mvmax),           &
            luc%snage_x(mland,mvmax),           &
            luc%osnowd_x(mland,mvmax),          &
            luc%cansto_x(mland,mvmax),          &
            luc%ghflux_x(mland,mvmax),          &
            luc%sghflux_x(mland,mvmax),         &
            luc%ga_x(mland,mvmax),              &
            luc%fev_x(mland,mvmax),             &
            luc%fes_x(mland,mvmax),             &
            luc%fhs_x(mland,mvmax),             &
            luc%wbtot0_x(mland,mvmax),          &
            luc%osnowd0_x(mland,mvmax),         & 
            luc%trad_x(mland,mvmax),            &
            luc%GWwb_x(mland,mvmax),            &
            luc%cplantx_x(mland,mvmax,ncp),     &
            luc%csoilx_x(mland,mvmax,ncs))

   ALLOCATE(luc%albsoilsn_y(mland,mvmax,nrb),   &
            luc%albedo_y(mland,mvmax,nrb),      &
            luc%albsoil_y(mland,mvmax,nrb),     &
            luc%dgdtg_y(mland,mvmax),           &
            luc%gammzz_y(mland,mvmax,ms),       &
            luc%tgg_y(mland,mvmax,ms),          &
            luc%wb_y(mland,mvmax,ms),           &
            luc%wbice_y(mland,mvmax,ms),        &
            luc%tggsn_y(mland,mvmax,msn),       &
            luc%ssdn_y(mland,mvmax,msn),        &
            luc%smass_y(mland,mvmax,msn),       &
            luc%sdepth_y(mland,mvmax,msn),      &
            luc%tss_y(mland,mvmax),             &
            luc%rtsoil_y(mland,mvmax),          &
            luc%runoff_y(mland,mvmax),          &
            luc%rnof1_y(mland,mvmax),           &
            luc%rnof2_y(mland,mvmax),           &
            luc%ssdnn_y(mland,mvmax),           &
            luc%snowd_y(mland,mvmax),           &
            luc%snage_y(mland,mvmax),           &
            luc%osnowd_y(mland,mvmax),          &
            luc%cansto_y(mland,mvmax),          &
            luc%ghflux_y(mland,mvmax),          &
            luc%sghflux_y(mland,mvmax),         &
            luc%ga_y(mland,mvmax),              &
            luc%fev_y(mland,mvmax),             &
            luc%fes_y(mland,mvmax),             &
            luc%fhs_y(mland,mvmax),             &
            luc%wbtot0_y(mland,mvmax),          &
            luc%osnowd0_y(mland,mvmax),         & 
            luc%trad_y(mland,mvmax),            &
            luc%GWwb_y(mland,mvmax),            &
            luc%cplantx_y(mland,mvmax,ncp),     &
            luc%csoilx_y(mland,mvmax,ncs))

    ! biogeochemical variables
   ALLOCATE(luc%cplant_x(mland,mvmax,mplant),    &
            luc%nplant_x(mland,mvmax,mplant),    &
            luc%pplant_x(mland,mvmax,mplant),    &
            luc%clitter_x(mland,mvmax,mlitter),  &
            luc%nlitter_x(mland,mvmax,mlitter),  &
            luc%plitter_x(mland,mvmax,mlitter),  &
            luc%csoil_x(mland,mvmax,msoil),      &
            luc%nsoil_x(mland,mvmax,msoil),      &
            luc%psoil_x(mland,mvmax,msoil),      &
            luc%clabile_x(mland,mvmax),          &
            luc%nsoilmin_x(mland,mvmax),         &
            luc%psoillab_x(mland,mvmax),         &
            luc%psoilsorb_x(mland,mvmax),        &
            luc%psoilocc_x(mland,mvmax),         &
            luc%cwoodprod_x(mland,mvmax,mwood),  &
            luc%nwoodprod_x(mland,mvmax,mwood),  &
            luc%pwoodprod_x(mland,mvmax,mwood),  &
            luc%cplant_y(mland,mvmax,mplant),    &
            luc%nplant_y(mland,mvmax,mplant),    &
            luc%pplant_y(mland,mvmax,mplant),    &
            luc%clitter_y(mland,mvmax,mlitter),  &
            luc%nlitter_y(mland,mvmax,mlitter),  &
            luc%plitter_y(mland,mvmax,mlitter),  &
            luc%csoil_y(mland,mvmax,msoil),      &
            luc%nsoil_y(mland,mvmax,msoil),      &
            luc%psoil_y(mland,mvmax,msoil),      &
            luc%clabile_y(mland,mvmax),          &
            luc%nsoilmin_y(mland,mvmax),         &
            luc%psoillab_y(mland,mvmax),         &
            luc%psoilsorb_y(mland,mvmax),        &
            luc%psoilocc_y(mland,mvmax),         &
            luc%cwoodprod_y(mland,mvmax,mwood),  &
            luc%nwoodprod_y(mland,mvmax,mwood),  &
            luc%pwoodprod_y(mland,mvmax,mwood))

    ! land use variables
   ALLOCATE(luc%pftfrac(mland,mvtype),           &
            luc%fharvw(mland,mharvw),            &
            luc%xluh2cable(mland,mvmax,mstate),  &
            luc%atransit(mland,mvmax,mvmax))

    !        luc%phen_y(mland,mvmax),             &
    !        luc%aphen_y(mland,mvmax),            &
    !        luc%doyphase3_y(mland,mvmax),        &
    !        luc%frac_sapwood_y(mland,mvmax),     &
    !        luc%sapwood_area_y(mland,mvmax))        
   ! Initialize temporary variables
           ! patch-genric variables
           luc%iveg_x   = -1;     luc%isoil_x=-1;          luc%soilorder_x=-1;     luc%phase_x=0;    luc%isflag_x=0
           luc%phen_x   = 0.0;    luc%aphen_x=0.0;         luc%doyphase3_x=-1;     luc%frac_sapwood_x=1.0;  luc%sapwood_area_x=0.0
           luc%patchfrac_x=0.0;   luc%lai_x=0.0;           luc%sla_x=0.0
           luc%iveg_y   = -1;     luc%isoil_y=-1;          luc%soilorder_y=-1;     luc%phase_y=0;    luc%isflag_y=0
           luc%phen_y   = 0.0;    luc%aphen_y=0.0;         luc%doyphase3_y=-1;     luc%frac_sapwood_y=1.0;  luc%sapwood_area_y=0.0
           luc%patchfrac_y=0.0;   luc%lai_y=0.0;           luc%sla_y=0.0

           ! biophysical
           luc%albsoilsn_x=0.0
           luc%albedo_x=0.0
           luc%albsoil_x=0.0
           luc%dgdtg_x=0.0
           luc%gammzz_x=0.0
           luc%tgg_x=0.0
           luc%wb_x=0.0
           luc%wbice_x=0.0
           luc%tggsn_x=0.0
           luc%ssdn_x=0.0
           luc%smass_x=0.0
           luc%sdepth_x=0.0
           luc%tss_x=0.0
           luc%rtsoil_x=0.0
           luc%runoff_x=0.0
           luc%rnof1_x=0.0
           luc%rnof2_x=0.0
           luc%ssdnn_x=0.0
           luc%snowd_x=0.0
           luc%snage_x=0.0
           luc%osnowd_x=0.0
           luc%cansto_x=0.0
           luc%ghflux_x=0.0
           luc%sghflux_x=0.0
           luc%ga_x=0.0
           luc%fev_x=0.0
           luc%fes_x=0.0
           luc%fhs_x=0.0
           luc%wbtot0_x=0.0
           luc%osnowd0_x=0.0
           luc%trad_x=0.0
           luc%GWwb_x=0.0
           luc%cplantx_x=0.0
           luc%csoilx_x=0.0

           luc%albsoilsn_y=0.0
           luc%albedo_y=0.0
           luc%albsoil_y=0.0
           luc%dgdtg_y=0.0
           luc%gammzz_y=0.0
           luc%tgg_y=0.0
           luc%wb_y=0.0
           luc%wbice_y=0.0
           luc%tggsn_y=0.0
           luc%ssdn_y=0.0
           luc%smass_y=0.0
           luc%sdepth_y=0.0
           luc%tss_y=0.0
           luc%rtsoil_y=0.0
           luc%runoff_y=0.0
           luc%rnof1_y=0.0
           luc%rnof2_y=0.0
           luc%ssdnn_y=0.0
           luc%snowd_y=0.0
           luc%snage_y=0.0
           luc%osnowd_y=0.0
           luc%cansto_y=0.0
           luc%ghflux_y=0.0
           luc%sghflux_y=0.0
           luc%ga_y=0.0
           luc%fev_y=0.0
           luc%fes_y=0.0
           luc%fhs_y=0.0
           luc%wbtot0_y=0.0
           luc%osnowd0_y=0.0
           luc%trad_y=0.0
           luc%GWwb_y=0.0
           luc%cplantx_y=0.0
           luc%csoilx_y=0.0

           ! biogeochemical
           luc%cplant_x  = 0.0;   luc%nplant_x   = 0.0;    luc%pplant_x   = 0.0
           luc%clitter_x = 0.0;   luc%nlitter_x  = 0.0;    luc%plitter_x  = 0.0
           luc%csoil_x   = 0.0;   luc%nsoil_x    = 0.0;    luc%psoil_x    = 0.0
           luc%clabile_x = 0.0;   luc%nsoilmin_x = 0.0;    luc%psoillab_x = 0.0
           luc%psoilsorb_x = 0.0; luc%psoilocc_x = 0.
           luc%cwoodprod_x=0.0;   luc%nwoodprod_x=0.0;     luc%pwoodprod_x=0.0

           luc%phase_y   = -1;    luc%patchfrac_y=0.0
           luc%cplant_y  = 0.0;   luc%nplant_y   = 0.0;    luc%pplant_y   = 0.0
           luc%clitter_y = 0.0;   luc%nlitter_y  = 0.0;    luc%plitter_y  = 0.0
           luc%csoil_y   = 0.0;   luc%nsoil_y    = 0.0;    luc%psoil_y    = 0.0
           luc%clabile_y = 0.0;   luc%nsoilmin_y = 0.0;    luc%psoillab_y = 0.0
           luc%psoilsorb_y = 0.0; luc%psoilocc_y = 0.0
           luc%cwoodprod_y=0.0;   luc%nwoodprod_y=0.0;     luc%pwoodprod_y=0.0

           luc%pftfrac = 0.0;     luc%fharvw=0.0;          luc%xluh2cable=0.0;    luc%atransit=0.0
   END SUBROUTINE landuse_allocate_mland

   SUBROUTINE landuse_deallocate_mland(luc)
   IMPLICIT NONE
   TYPE(landuse_mland), INTENT(INOUT)  :: luc

   !patch-generic variables
   DEALLOCATE(luc%iveg_x,      luc%isoil_x, luc%soilorder_x, luc%phase_x, luc%isflag_x)
   DEALLOCATE(luc%phen_x,      luc%aphen_x, luc%doyphase3_x, luc%frac_sapwood_x,  luc%sapwood_area_x)
   DEALLOCATE(luc%patchfrac_x, luc%lai_x,   luc%sla_x)

   DEALLOCATE(luc%iveg_y,      luc%isoil_y, luc%soilorder_y, luc%phase_y, luc%isflag_y)
   DEALLOCATE(luc%phen_y,      luc%aphen_y, luc%doyphase3_y, luc%frac_sapwood_y,  luc%sapwood_area_y)
   DEALLOCATE(luc%patchfrac_y, luc%lai_y,   luc%sla_y)


   ! biophysical variables
   DEALLOCATE(luc%albsoilsn_x,    &
              luc%albedo_x,       &
              luc%albsoil_x,      &
              luc%dgdtg_x,        &
              luc%gammzz_x,       &
              luc%tgg_x,          &
              luc%wb_x,           &
              luc%wbice_x,        &
              luc%tggsn_x,        &
              luc%ssdn_x,         &
              luc%smass_x,        &
              luc%sdepth_x,       &
              luc%tss_x,          &
              luc%rtsoil_x,       &
              luc%runoff_x,       &
              luc%rnof1_x,        &
              luc%rnof2_x,        &
              luc%ssdnn_x,        &
              luc%snowd_x,        &
              luc%snage_x,        &
              luc%osnowd_x,       &
              luc%cansto_x,       &
              luc%ghflux_x,       &
              luc%sghflux_x,      &
              luc%ga_x,           &
              luc%fev_x,          &
              luc%fes_x,          &
              luc%fhs_x,          &
              luc%wbtot0_x,       &
              luc%osnowd0_x,      &
              luc%trad_x,         &
              luc%GWwb_x,         &
              luc%cplantx_x,      &
              luc%csoilx_x)

   DEALLOCATE(luc%albsoilsn_y,    &
              luc%albedo_y,       &
              luc%albsoil_y,      &
              luc%dgdtg_y,        &
              luc%gammzz_y,       &
              luc%tgg_y,          &
              luc%wb_y,           &
              luc%wbice_y,        &
              luc%tggsn_y,        &
              luc%ssdn_y,         &
              luc%smass_y,        &
              luc%sdepth_y,       &
              luc%tss_y,          &
              luc%rtsoil_y,       &
              luc%runoff_y,       &
              luc%rnof1_y,        &
              luc%rnof2_y,        &
              luc%ssdnn_y,        &
              luc%snowd_y,        &
              luc%snage_y,        &
              luc%osnowd_y,       &
              luc%cansto_y,       &
              luc%ghflux_y,       &
              luc%sghflux_y,      &
              luc%ga_y,           &
              luc%fev_y,          &
              luc%fes_y,          &
              luc%fhs_y,          &
              luc%wbtot0_y,       &
              luc%osnowd0_y,      &
              luc%trad_y,         &
              luc%GWwb_y,         &
              luc%cplantx_y,      &
              luc%csoilx_y)

   ! biogeochemical variables
   DEALLOCATE(luc%cplant_x,          &
              luc%nplant_x,          &
              luc%pplant_x,          &
              luc%clitter_x,         &
              luc%nlitter_x,         &
              luc%plitter_x,         &
              luc%csoil_x,           &
              luc%nsoil_x,           &
              luc%psoil_x,           &
              luc%clabile_x,         &
              luc%nsoilmin_x,        &
              luc%psoillab_x,        &
              luc%psoilsorb_x,       &
              luc%psoilocc_x,        &
              luc%cwoodprod_x,       &
              luc%nwoodprod_x,       &
              luc%pwoodprod_x,       &
              luc%cplant_y,          &
              luc%nplant_y,          &
              luc%pplant_y,          &
              luc%clitter_y,         &
              luc%nlitter_y,         &
              luc%plitter_y,         &
              luc%csoil_y,           &
              luc%nsoil_y,           &
              luc%psoil_y,           &
              luc%clabile_y,         &
              luc%nsoilmin_y,        &
              luc%psoillab_y,        &
              luc%psoilsorb_y,       &
              luc%psoilocc_y,        &
              luc%cwoodprod_y,       &
              luc%nwoodprod_y,       &
              luc%pwoodprod_y)
   !land use variables    
   DEALLOCATE(luc%pftfrac,           &
              luc%fharvw,            &
              luc%xluh2cable,        &
              luc%atransit)

   END SUBROUTINE landuse_deallocate_mland

   SUBROUTINE landuse_allocate_mp(mpx,ms,msn,nrb,mplant,mlitter,msoil,mwood,ncp,ncs,lucmp)
   integer    mpx,ms,msn,nrb,mplant,mlitter,msoil,mwood,ncp,ncs
   TYPE(landuse_mp), INTENT(INOUT)  :: lucmp

   ! generic patch properties
     allocate(lucmp%iveg(mpx),lucmp%isoil(mpx),lucmp%soilorder(mpx),lucmp%phase(mpx),lucmp%isflag(mpx))

     allocate(lucmp%lat(mpx),lucmp%lon(mpx))
     allocate(lucmp%doyphase3(mpx))
     allocate(lucmp%phen(mpx),lucmp%aphen(mpx),lucmp%frac_sapwood(mpx), lucmp%sapwood_area(mpx))       !(mp)

     allocate(lucmp%patchfrac(mpx),lucmp%areacell(mpx),lucmp%lai(mpx),lucmp%sla(mpx))

   ! biophysical variables
     allocate(lucmp%albsoilsn(mpx,nrb),lucmp%albedo(mpx,nrb),lucmp%albsoil(mpx,nrb))                !float(mp,rad)
     allocate(lucmp%dgdtg(mpx))                                                                     !double(mp)
     allocate(lucmp%gammzz(mpx,ms))                                                                    !double(mp,soil)
     allocate(lucmp%tgg(mpx,ms),lucmp%wb(mpx,ms),lucmp%wbice(mpx,ms))                               !float(mp,soil)
     allocate(lucmp%tggsn(mpx,msn),lucmp%ssdn(mpx,msn),lucmp%smass(mpx,msn),lucmp%sdepth(mpx,msn)) !float(mp,snow)
 
     allocate(lucmp%tss(mpx),lucmp%rtsoil(mpx),lucmp%runoff(mpx),lucmp%rnof1(mpx),lucmp%rnof2(mpx), &
              lucmp%ssdnn(mpx),lucmp%snowd(mpx),lucmp%snage(mpx),lucmp%osnowd(mpx),                 &
              lucmp%cansto(mpx),lucmp%ghflux(mpx),lucmp%sghflux(mpx),lucmp%ga(mpx),                 &
              lucmp%fev(mpx),lucmp%fes(mpx),lucmp%fhs(mpx),lucmp%wbtot0(mpx),lucmp%osnowd0(mpx),    &
              lucmp%trad(mpx),lucmp%GWwb(mpx))                                                  !float(mp)
     allocate(lucmp%cplantx(mpx,ncp), lucmp%csoilx(mpx,ncs))                                    !float(mp,ncp/ncs)

    ! biogeochemical variables 
     allocate(lucmp%sumcbal(mpx),lucmp%sumnbal(mpx),lucmp%sumpbal(mpx))
     allocate(lucmp%clabile(mpx))
     allocate(lucmp%cplant(mpx,mplant),lucmp%nplant(mpx,mplant),lucmp%pplant(mpx,mplant))
     allocate(lucmp%clitter(mpx,mlitter),lucmp%nlitter(mpx,mlitter),lucmp%plitter(mpx,mlitter))
     allocate(lucmp%csoil(mpx,msoil),lucmp%nsoil(mpx,msoil),lucmp%psoil(mpx,msoil))
     allocate(lucmp%nsoilmin(mpx))
     allocate(lucmp%psoillab(mpx),lucmp%psoilsorb(mpx),lucmp%psoilocc(mpx))
     allocate(lucmp%cwoodprod(mpx,mwood),lucmp%nwoodprod(mpx,mwood),lucmp%pwoodprod(mpx,mwood))

     ! initialization
     lucmp%iveg=-1;lucmp%isoil=-1;lucmp%soilorder=-1;lucmp%phase=0;lucmp%isflag=0
     lucmp%doyphase3=-1;lucmp%phen=0.0;lucmp%aphen=0.0;lucmp%frac_sapwood=1.0;lucmp%sapwood_area=0.0       !(mp)

     lucmp%patchfrac=0.0;lucmp%areacell=0.0;lucmp%lai=0.0;lucmp%sla=0.0

     ! biophysical variables
     lucmp%albsoilsn(:,:)=0.0; lucmp%albedo(:,:)=0.0; lucmp%albsoil(:,:)=0.0                 !float(mp,rad)
     lucmp%dgdtg(:)=0.0                                                                      !double(mp)
     lucmp%gammzz(:,:)=0.0                                                                   !double(mp,soil)
     lucmp%tgg(:,:)=0.0;lucmp%wb(:,:)=0.0; lucmp%wbice(:,:)=0.0                              !float(mp,soil)
     lucmp%tggsn(:,:)=0.0; lucmp%ssdn(:,:)=0.0; lucmp%smass(:,:)=0.0; lucmp%sdepth(:,:)=0.0  !float(mp,snow)
 
     lucmp%tss(:)=0.0; lucmp%rtsoil(:)=0.0; lucmp%runoff(:)=0.0; lucmp%rnof1(:)=0.0; lucmp%rnof2(:)=0.0
     lucmp%ssdnn(:)=0.0; lucmp%snowd(:)=0.0; lucmp%snage(:)=0.0; lucmp%osnowd(:)=0.0
     lucmp%cansto(:)=0.0; lucmp%ghflux(:)=0.0; lucmp%sghflux(:)=0.0; lucmp%ga(:)=0.0
     lucmp%fev(:)=0.0; lucmp%fes(:)=0.0; lucmp%fhs(:)=0.0; lucmp%wbtot0(:)=0.0; lucmp%osnowd0(:)=0.0
     lucmp%trad(:)=0.0; lucmp%GWwb(:)=0.0                                              !float(mp)
     lucmp%cplantx(:,:)=0.0; lucmp%csoilx(:,:)=0.0                                     !float(mp,plant_carbon_pools/soil_carbon_pools)

     ! biogeochemical variables
     lucmp%sumcbal=0.0;lucmp%sumnbal=0.0;lucmp%sumpbal=0.0
     lucmp%clabile = 0.0
     lucmp%cplant=0.0;  lucmp%nplant=0.0;  lucmp%pplant=0.0
     lucmp%clitter=0.0; lucmp%nlitter=0.0; lucmp%plitter=0.0
     lucmp%csoil=0.0;   lucmp%nsoil=0.0;   lucmp%psoil=0.0
     lucmp%nsoilmin=0.0
     lucmp%psoillab=0.0;lucmp%psoilsorb=0.0;lucmp%psoilocc=0.0
     lucmp%cwoodprod=0.0;lucmp%nwoodprod=0.0;lucmp%pwoodprod=0.0
   END SUBROUTINE landuse_allocate_mp

   SUBROUTINE landuse_deallocate_mp(mpx,ms,msn,nrb,mplant,mlitter,msoil,mwood,lucmp)
   integer     mpx,ms,msn,nrb,mplant,mlitter,msoil,mwood
   TYPE(landuse_mp), INTENT(INOUT)  :: lucmp
     ! patch-generic variables
     deallocate(lucmp%iveg,lucmp%isoil,lucmp%soilorder,lucmp%phase,lucmp%isflag)
     deallocate(lucmp%lat,lucmp%lon)
     deallocate(lucmp%doyphase3)
     deallocate(lucmp%phen,lucmp%aphen,lucmp%frac_sapwood, lucmp%sapwood_area)       !(mp)

     deallocate(lucmp%patchfrac,lucmp%areacell,lucmp%lai,lucmp%sla)

     ! biophysical variables
     deallocate(lucmp%albsoilsn,lucmp%albedo,lucmp%albsoil)                   !float(mp,rad)
     deallocate(lucmp%dgdtg)                                                  !double(mp)
     deallocate(lucmp%gammzz)                                                 !double(mp,soil)
     deallocate(lucmp%tgg,lucmp%wb,lucmp%wbice)                               !float(mp,soil)
     deallocate(lucmp%tggsn,lucmp%ssdn,lucmp%smass,lucmp%sdepth)             !float(mp,snow)
 
     deallocate(lucmp%tss,lucmp%rtsoil,lucmp%runoff,lucmp%rnof1,lucmp%rnof2,    &
              lucmp%ssdnn,lucmp%snowd,lucmp%snage,lucmp%osnowd,                 &
              lucmp%cansto,lucmp%ghflux,lucmp%sghflux,lucmp%ga,                 &
              lucmp%fev,lucmp%fes,lucmp%fhs,lucmp%wbtot0,lucmp%osnowd0,         &
              lucmp%trad,lucmp%GWwb)                        !float(mp)
     deallocate(lucmp%cplantx, lucmp%csoilx)                !float(mp,plant_carbon_pools/soil_carbon_pools)


     !biogeochemical variables
     deallocate(lucmp%sumcbal,lucmp%sumnbal,lucmp%sumpbal)
     deallocate(lucmp%clabile)
     deallocate(lucmp%cplant,lucmp%nplant,lucmp%pplant)
     deallocate(lucmp%clitter,lucmp%nlitter,lucmp%plitter)
     deallocate(lucmp%csoil,lucmp%nsoil,lucmp%psoil)
     deallocate(lucmp%nsoilmin)
     deallocate(lucmp%psoillab,lucmp%psoilsorb,lucmp%psoilocc)
     deallocate(lucmp%cwoodprod,lucmp%nwoodprod,lucmp%pwoodprod)

   END SUBROUTINE landuse_deallocate_mp

END MODULE landuse_variable

  subroutine landuse_driver(mlon,mlat,landmask,arealand,ssnow,soil,veg,bal,canopy,  &
                            phen,casapool,casabal,casamet,casabiome,casaflux,bgc,rad, &
                            cstart,cend,nap,lucmp)
  USE cable_IO_vars_module, ONLY: mask,patch,landpt, latitude, longitude
  USE cable_def_types_mod,  ONLY: mp,mvtype,mstype,mland,r_2,ms,msn,nrb,ncp,ncs,           &
                                  soil_parameter_type, soil_snow_type, veg_parameter_type, &
                                  balances_type, canopy_type, bgc_pool_type, radiation_type
  USE casadimension,        ONLY: icycle,mplant,mlitter,msoil,mwood,mso
  USE casavariable,         ONLY: casa_pool,casa_balance,casa_met,casa_biome,casa_flux
  USE phenvariable,         ONLY: phen_variable
  USE landuse_variable
  IMPLICIT NONE
  TYPE (soil_snow_type)          :: ssnow   ! soil and snow variables
  TYPE (soil_parameter_type)     :: soil    ! soil parameters
  TYPE (veg_parameter_type)      :: veg     ! vegetation parameters
  TYPE (balances_type)           :: bal
  TYPE (canopy_type)             :: canopy
  TYPE (phen_variable)           :: phen
  TYPE (casa_pool)               :: casapool
  TYPE (casa_biome)              :: casabiome
  TYPE (casa_balance)            :: casabal
  TYPE (casa_met)                :: casamet
  TYPE (casa_flux)               :: casaflux
  TYPE (bgc_pool_type)           :: bgc
  TYPE (radiation_type)          :: rad    ! met data

  TYPE (landuse_mland)           :: luc
  TYPE (landuse_mp)              :: lucmp
  ! input
  integer mlon,mlat
  integer,       dimension(mlon,mlat)         :: landmask
  real(r_2),     dimension(mland)             :: arealand
  ! output
  ! "mland" variables
  integer,       dimension(mland)             :: cstart,cend,nap

  character*500   fxpft,fxluh2cable
  integer ivt,ee,hh,np,p,q,np1
  integer ncid,ok,xID,yID,varID,i,j,m,mpx

     print *, 'calling allocate mp: landuse'
     call landuse_allocate_mp(mp,ms,msn,nrb,mplant,mlitter,msoil,mwood,ncp,ncs,lucmp)  
     print *, 'calling allocate mland: landuse'
     call landuse_allocate_mland(mland,luc)                                                     !setup "varx(mland,:)"        
     print *, 'exiting  allocating mland: landuse'

     ! get the mapping matrix from state to PFT
     ! call landuse_getxluh2(mlat,mlon,landmask,luc,filename%fxluh2cable)    !"xluh2cable"
     ! call landuse_getdata(mlat,mlon,landmask,filename%fxpft,luc)     !"luc(t-1)" and "xpft(t-1)"

     ! get pool sizes and other states in the "restart", "gridinfo" and "poolout" file
     ! patch-generic variables
     do p=1,mp
!        print *, 'p', p, veg%iveg(p),soil%isoilm(p),ssnow%isflag(p)
        lucmp%iveg(p)      = veg%iveg(p)
        lucmp%isoil(p)     = soil%isoilm(p)          
        lucmp%soilorder(p) = casamet%isorder(p)          
        lucmp%isflag(p)    = ssnow%isflag(p)

     enddo   
     print *, 'point A: landuse'
     !
     print *, 'patchfraC',size(patch%frac)
     print *, 'veglai= ',size(veg%vlai)
     print *, 'landuse: casabiome:sla', size(casabiome%sla),  casabiome%sla(:) 

     do p=1,mp
     !   print *, 'landuse b', p, veg%iveg(p),veg%vlai(p),patch(p)%frac

        lucmp%patchfrac(p) = patch(p)%frac             ! maybe we should create another variable for "primary%patch"
        lucmp%lai(p)       = veg%vlai(p)
        lucmp%sla(p)       = casabiome%sla(veg%iveg(p)) 
     enddo
     print *, 'point b: landuse'

     ! biophysical variables 
     do p=1,mp
        lucmp%albsoilsn(p,:) = ssnow%albsoilsn(p,:)
        lucmp%albedo(p,:)    = rad%albedo(p,:)
        lucmp%albsoil(p,:)   = soil%albsoil(p,:)
        lucmp%gammzz(p,:)    = ssnow%gammzz(p,:)
        lucmp%tgg(p,:)       = ssnow%tgg(p,:)
        lucmp%wb(p,:)        = ssnow%wb(p,:)
        lucmp%wbice(p,:)     = ssnow%wbice(p,:)
        lucmp%tggsn(p,:)     = ssnow%tggsn(p,:)
        lucmp%ssdn(p,:)      = ssnow%ssdn(p,:)
        lucmp%smass(p,:)     = ssnow%smass(p,:)
        lucmp%sdepth(p,:)    = ssnow%sdepth(p,:)
        lucmp%tss(p)         = ssnow%tss(p)
     enddo

     print *, 'point C: landuse'
     lucmp%runoff(:)      = ssnow%runoff(:)
     lucmp%rnof1(:)       = ssnow%rnof1(:)
     lucmp%rnof2(:)       = ssnow%rnof2(:)
     lucmp%ssdnn(:)       = ssnow%ssdnn(:)
     lucmp%snowd(:)       = ssnow%snowd(:)
     lucmp%snage(:)       = ssnow%snage(:)
     lucmp%osnowd(:)      = ssnow%osnowd(:)
     lucmp%cansto(:)      = canopy%cansto(:)
     lucmp%ghflux(:)      = canopy%ghflux(:)
     lucmp%sghflux(:)     = canopy%sghflux(:)
     print *, 'point D: landuse'
     lucmp%ga(:)          = canopy%ga(:)
     lucmp%dgdtg(:)       = canopy%dgdtg(:)
     lucmp%fev(:)         = canopy%fev(:)
     lucmp%fes(:)         = canopy%fes(:)
     lucmp%fhs(:)         = canopy%fhs(:) 
     lucmp%wbtot0(:)      = bal%wbtot0(:)
     lucmp%osnowd0(:)     = bal%osnowd0(:)
     lucmp%trad(:)        = rad%trad(:)
     lucmp%GWwb(:)        = ssnow%GWwb(:)
     lucmp%cplantx(:,:)   = bgc%cplant(:,:)
     lucmp%csoilx(:,:)    = bgc%csoil(:,:)

     print *, 'point E: landuse'
     ! biogeochemical variables    
     do m=1,mland
       do np=cstart(m),cend(m)
          ivt = lucmp%iveg(np)
          if(ivt <=mvtype) then
             luc%pftfrac(m,ivt) = patch(np)%frac
          endif
       enddo
     enddo

     print *, 'point F: landuse'
     if(icycle>0) then 
     do p=1,mp        
  !      print *, 'landuse F: ', p, phen%phase(p),phen%doyphase(p,3),phen%phen(p),phen%aphen(p)
 !       print *, 'landuse F2: ',   casaflux%frac_sapwood(p),casaflux%sapwood_area(p)
 !       print *, 'landuse F3: ', casapool%clabile(p),casapool%cplant(p,:),casapool%clitter(p,:),casapool%csoil(p,:),        &
 !                               casapool%cwoodprod(p,:)
!        print *, 'landuse F4: ',casabal%sumcbal(p)

        lucmp%phase(p)       = phen%phase(p)
        lucmp%doyphase3(p)   = phen%doyphase(p,3)
        lucmp%phen(p)        = phen%phen(p)
        lucmp%aphen(p)       = phen%aphen(p)
        lucmp%frac_sapwood(p)= casaflux%frac_sapwood(p)
        lucmp%sapwood_area(p)= casaflux%sapwood_area(p)
        lucmp%clabile(p)     = casapool%clabile(p)
        lucmp%cplant(p,:)    = casapool%cplant(p,:)
        lucmp%clitter(p,:)   = casapool%clitter(p,:)
        lucmp%csoil(p,:)     = casapool%csoil(p,:)
        lucmp%cwoodprod(p,:) = casapool%cwoodprod(p,:)
        lucmp%sumcbal(p)     = casabal%sumcbal(p)
     enddo
     endif
     print *, 'point G: landuse'
     if(icycle>1) then
     do p=1,mp        
        lucmp%nplant(p,:)    = casapool%nplant(p,:)
        lucmp%nlitter(p,:)   = casapool%nlitter(p,:)
        lucmp%nsoil(p,:)     = casapool%nsoil(p,:)
        lucmp%nwoodprod(p,:) = casapool%nwoodprod(p,:)
        lucmp%nsoilmin(p)    = casapool%nsoilmin(p)
        lucmp%sumnbal(p)     = casabal%sumnbal(p)
     enddo
     endif
     print *, 'point H: landuse'
     if(icycle >2) then
     do p=1,mp        
        lucmp%pplant(p,:)    = casapool%pplant(p,:)
        lucmp%plitter(p,:)   = casapool%plitter(p,:)
        lucmp%psoil(p,:)     = casapool%psoil(p,:)
        lucmp%pwoodprod(p,:) = casapool%pwoodprod(p,:)
        lucmp%psoillab(p)    = casapool%psoillab(p)
        lucmp%psoilsorb(p)   = casapool%psoilsorb(p)
        lucmp%psoilocc(p)    = casapool%psoilocc(p)
        lucmp%sumpbal(p)     = casabal%sumpbal(p)
     enddo
     endif
      
     ! assign variables var(mp,:) to luc%var_x(mland,mvmax,:)
     print *, 'calling mp2land: landuse'
     call landuse_mp2land(luc,lucmp,mp,cstart,cend)

     ! we need to deallocate "lucmp" because "mp" will be updated after land use change
     print *, 'calling deallocate mp: landuse'
     call landuse_deallocate_mp(mp,ms,msn,nrb,mplant,mlitter,msoil,mwood,lucmp)

     print *, 'calling transitx: landuse'
     call landuse_transitx(luc,casabiome)

     print *, 'calling checks: landuse'
     call landuse_checks(mlon,mlat,landmask,luc)

     print *, 'calling update mland: landuse'
     call landuse_update_mland(luc)                    ! assign "var_y" to "var_x"

     ! update "cstart", "cend","nap" and "mp=>mpx"
      cstart=0;cend=0;nap=0
      np =0; cstart(:) = 0; cend(:) =0; nap(:) = 0
      do p=1,mland
         np1 =0
         if(sum(luc%patchfrac_y(p,:))<thresh_frac) then
            print *, 'WARNING! patch area sum too low',p,luc%patchfrac_y(p,:)
         else
            do q=1,mvmax
               if(luc%patchfrac_y(p,q) >thresh_frac) then
                  np  = np + 1
                  np1 = np1 + 1
                  if(np1==1) cstart(p) = np
              endif
            enddo
            cend(p) = np
            nap(p)  = max(0,cend(p)-cstart(p) +1)
         endif
      enddo
      mpx = np
     ! allocate "lucmp" with "mpx"
     print *, 'calling allocate mp: landuse'
     call landuse_allocate_mp(mpx,ms,msn,nrb,mplant,mlitter,msoil,mwood,ncp,ncs,lucmp)
 
     ! assign lucmp%lat lucmp%lon
     do p=1,mland
        do q=cstart(p),cend(p)
           lucmp%lat(q) = latitude(p)
           lucmp%lon(q) = longitude(p)
        enddo
     enddo      

     print *, 'calling land2mpx: landuse'
     call landuse_land2mpx(luc,lucmp,mpx)
  !   call landuse_land2mpx(luc,lucmp,mpx,cstart,cend,nap)

     print *, 'calling deallocate mland: landuse'
     call landuse_deallocate_mland(luc)

     print *, 'landuse: exit landuse_driver mpx', mpx

     close(21)
211  format(i4,a120)
end subroutine landuse_driver

 SUBROUTINE landuse_mp2land(luc,lucmp,mp,cstart,cend)
 use landuse_variable
 USE cable_def_types_mod,  ONLY: mvtype,mstype,mland,r_2,ms,msn,nrb,ncp,ncs
 USE casadimension,        ONLY: icycle,mplant,mlitter,msoil,mwood,mso
 IMPLICIT NONE
 integer mp
 type(landuse_mland)   :: luc
 type(landuse_mp)      :: lucmp
 integer g,np,ivt,i
 integer,        dimension(mland)        :: cstart,cend

  do g=1,mland
     do ivt=1,mvmax
        luc%iveg_x(g,ivt) = ivt
     enddo

  do np= cstart(g),cend(g)

     ivt = lucmp%iveg(np)
     if(ivt<1.or.ivt>17) then
        print *, 'at landuse_mp2land: vegtype outy of range!',g,np,ivt
        print *, 'stop!'
     endif

     ! patch-genric variables 
     luc%isoil_x(g,ivt)       = lucmp%isoil(np)
     luc%soilorder_x(g,ivt)   = lucmp%soilorder(np) 
     luc%phase_x(g,ivt)       = lucmp%phase(np)
  
     luc%doyphase3_x(g,ivt)   = lucmp%doyphase3(np)
     luc%phen_x(g,ivt)        = lucmp%phen(np)
     luc%aphen_x(g,ivt)       = lucmp%aphen(np)
     luc%frac_sapwood_x(g,ivt)= lucmp%frac_sapwood(np)
     luc%sapwood_area_x(g,ivt)= lucmp%sapwood_area(np)

     luc%isflag_x(g,ivt)      = lucmp%isflag(np)
     luc%patchfrac_x(g,ivt)   = lucmp%patchfrac(np)   
     luc%lai_x(g,ivt)         = lucmp%lai(np)
     luc%sla_x(g,ivt)         = lucmp%sla(np)  
    
     ! biophysical variables
     do i=1,nrb
        luc%albsoilsn_x(g,ivt,i) = lucmp%albsoilsn(np,i)
        luc%albedo_x(g,ivt,i)    = lucmp%albedo(np,i)
        luc%albsoil_x(g,ivt,i)   = lucmp%albsoil(np,i)
     enddo 

     luc%dgdtg_x(g,ivt)        = lucmp%dgdtg(np)

     do i=1,ms
        luc%gammzz_x(g,ivt,i)  = lucmp%gammzz(np,i)
        luc%tgg_x(g,ivt,i)     = lucmp%tgg(np,i)
        luc%wb_x(g,ivt,i)      = lucmp%wb(np,i)
        luc%wbice_x(g,ivt,i)   = lucmp%wbice(np,i)
     enddo

     do i=1,msn
        luc%tggsn_x(g,ivt,i)   = lucmp%tggsn(np,i)
        luc%ssdn_x(g,ivt,i)    = lucmp%ssdn(np,i)
        luc%smass_x(g,ivt,i)   = lucmp%smass(np,i)
        luc%sdepth_x(g,ivt,i)  = lucmp%sdepth(np,i)
     enddo

     luc%tss_x(g,ivt)          = lucmp%tss(np)
     luc%rtsoil_x(g,ivt)       = lucmp%rtsoil(np)
     luc%runoff_x(g,ivt)       = lucmp%runoff(np)
     luc%rnof1_x(g,ivt)        = lucmp%rnof1(np)
     luc%rnof2_x(g,ivt)        = lucmp%rnof2(np)
     luc%ssdnn_x(g,ivt)        = lucmp%ssdnn(np)
     luc%snowd_x(g,ivt)        = lucmp%snowd(np)
     luc%snage_x(g,ivt)        = lucmp%snage(np)
     luc%osnowd_x(g,ivt)       = lucmp%osnowd(np)
     luc%cansto_x(g,ivt)       = lucmp%cansto(np)
     luc%ghflux_x(g,ivt)       = lucmp%ghflux(np)
     luc%sghflux_x(g,ivt)      = lucmp%sghflux(np)
     luc%ga_x(g,ivt)           = lucmp%ga(np)
     luc%fev_x(g,ivt)          = lucmp%fev(np)
     luc%fes_x(g,ivt)          = lucmp%fes(np)
     luc%fhs_x(g,ivt)          = lucmp%fhs(np)
     luc%wbtot0_x(g,ivt)       = lucmp%wbtot0(np)
     luc%osnowd0_x(g,ivt)      = lucmp%osnowd0(np) 
     luc%trad_x(g,ivt)         = lucmp%trad(np)
     luc%GWwb_x(g,ivt)         = lucmp%GWwb(np)

     do i =1,ncp
        luc%cplantx_x(g,ivt,i) = lucmp%cplantx(np,i)
     enddo
     do i=1,ncs
        luc%csoilx_x(g,ivt,i)  = lucmp%csoilx(np,i)
     enddo

     ! biogeochemical variables
     luc%clabile_x(g,ivt)     = lucmp%clabile(np)
     luc%cplant_x(g,ivt,:)    = lucmp%cplant(np,:)
     luc%clitter_x(g,ivt,:)   = lucmp%clitter(np,:)
     luc%csoil_x(g,ivt,:)     = lucmp%csoil(np,:)
     luc%cwoodprod_x(g,ivt,:) = lucmp%cwoodprod(np,:)

    IF(icycle>1) THEN
      luc%nplant_x(g,ivt,:)    = lucmp%nplant(np,:)
      luc%nlitter_x(g,ivt,:)   = lucmp%nlitter(np,:)
      luc%nsoil_x(g,ivt,:)     = lucmp%nsoil(np,:)
      luc%nsoilmin_x(g,ivt)    = lucmp%nsoilmin(np)
      luc%nwoodprod_x(g,ivt,:) = lucmp%nwoodprod(np,:)
    END IF
    IF(icycle>2) THEN
      luc%pplant_x(g,ivt,:)    = lucmp%pplant(np,:)
      luc%plitter_x(g,ivt,:)   = lucmp%plitter(np,:)
      luc%psoil_x(g,ivt,:)     = lucmp%psoil(np,:)
      luc%psoillab_x(g,ivt)    = lucmp%psoillab(np)
      luc%psoilsorb_x(g,ivt)   = lucmp%psoilsorb(np)
      luc%psoilocc_x(g,ivt)    = lucmp%psoilocc(np)
      luc%pwoodprod_x(g,ivt,:) = lucmp%pwoodprod(np,:)
    END IF
    
  enddo
  enddo
  
  ! patch-genric variables
  luc%isoil_y       = luc%isoil_x
  luc%soilorder_y   = luc%soilorder_x
  luc%phase_y       = luc%phase_x 
  luc%isflag_y      = luc%isflag_x 
  luc%patchfrac_y   = luc%patchfrac_x  
  luc%lai_y         = luc%lai_x 
  luc%sla_y         = luc%sla_x
  
  luc%doyphase3_y   = luc%doyphase3_x
  luc%phen_y        = luc%phen_x
  luc%aphen_y       = luc%aphen_x
  luc%frac_sapwood_y= luc%frac_sapwood_x
  luc%sapwood_area_y= luc%sapwood_area_x

  ! biophysical variables
  luc%albsoilsn_y   = luc%albsoilsn_x
  luc%albedo_y      = luc%albedo_x 
  luc%albsoil_y     = luc%albsoil_x
  luc%dgdtg_y       = luc%dgdtg_x 
  luc%gammzz_y      = luc%gammzz_x
  luc%tgg_y         = luc%tgg_x 
  luc%wb_y          = luc%wb_x
  luc%wbice_y       = luc%wbice_x 
  luc%tggsn_y       = luc%tggsn_x
  luc%ssdn_y        = luc%ssdn_x 
  luc%smass_y       = luc%smass_x
  luc%sdepth_y      = luc%sdepth_x
  luc%tss_y         = luc%tss_x
  luc%rtsoil_y      = luc%rtsoil_x 
  luc%runoff_y      = luc%runoff_x
  luc%rnof1_y       = luc%rnof1_x
  luc%rnof2_y       = luc%rnof2_x
  luc%ssdnn_y       = luc%ssdnn_x  
  luc%snowd_y       = luc%snowd_x 
  luc%snage_y       = luc%snage_x
  luc%osnowd_y      = luc%osnowd_x 
  luc%cansto_y      = luc%cansto_x
  luc%ghflux_y      = luc%ghflux_x
  luc%sghflux_y     = luc%sghflux_x
  luc%ga_y          = luc%ga_x 
  luc%fev_y         = luc%fev_x
  luc%fes_y         = luc%fes_x
  luc%fhs_y         = luc%fhs_x    
  luc%wbtot0_y      = luc%wbtot0_x
  luc%osnowd0_y     = luc%osnowd0_x 
  luc%trad_y        = luc%trad_x
  luc%GWwb_y        = luc%GWwb_x 
  luc%cplantx_y     = luc%cplantx_x
  luc%csoilx_y      = luc%csoilx_x 

  ! biogeochemical variables
  luc%clabile_y   = luc%clabile_x
  luc%cplant_y    = luc%cplant_x
  luc%clitter_y   = luc%clitter_x
  luc%csoil_y     = luc%csoil_x
  luc%cwoodprod_y = luc%cwoodprod_x

  IF(icycle>1) THEN
     luc%nplant_y    = luc%nplant_x
     luc%nlitter_y   = luc%nlitter_x
     luc%nsoil_y     = luc%nsoil_x
     luc%nsoilmin_y  = luc%nsoilmin_x
     luc%nwoodprod_y = luc%nwoodprod_x
  END IF
  IF(icycle>2) THEN
     luc%pplant_y    = luc%pplant_x
     luc%plitter_y   = luc%plitter_x
     luc%psoil_y     = luc%psoil_x
     luc%psoillab_y  = luc%psoillab_x
     luc%psoilsorb_y = luc%psoilsorb_x
     luc%psoilocc_y  = luc%psoilocc_x
     luc%pwoodprod_y = luc%pwoodprod_x
  END IF

END SUBROUTINE landuse_mp2land
  
SUBROUTINE landuse_transitx(luc,casabiome)
   USE casaparm
   USE landuse_constant
   USE casavariable,        ONLY: casa_biome
   USE landuse_variable,    ONLY: landuse_mland

   USE cable_def_types_mod,  ONLY: mland,mvtype,r_2,nrb,ncp,ncs
   USE casadimension,        ONLY: icycle,mplant,mlitter,msoil,mwood

   IMPLICIT NONE
   TYPE(casa_biome)    :: casabiome
   TYPE(landuse_mland) :: luc
   integer,   dimension(mvtype)               :: ivt2
   real(r_2), dimension(mland,mvmax)          :: dclabile
   real(r_2), dimension(mland,mvmax,mplant)   :: dcplant,dnplant,dpplant
   real(r_2), dimension(mland,mvmax,mlitter)  :: dclitter,dnlitter,dplitter
   real(r_2), dimension(mland,mvmax,msoil)    :: dcsoil,dnsoil,dpsoil
   real(r_2), dimension(mland,mvmax)          :: dnsoilmin,dpsoillab,dpsoilsorb,dpsoilocc
   real(r_2), dimension(mland,mvmax,mwood)    :: dcwoodprod,dnwoodprod,dpwoodprod
   real(r_2), DIMENSION(3)                    :: ratioLignintoN
   real(r_2), DIMENSION(mvmax,3,3)            :: fromPtoL
   real(r_2), DIMENSION(mland,mvmax)          :: delarea
   real(r_2), DIMENSION(mvmax,mvmax)          :: afwhpri, afwhsec, transitx
   real(r_2), DIMENSION(mvmax)                :: delfwhpri,delfwhsec
   real(r_2), DIMENSION(3)                    :: totcwoodprod, totclitter, totcsoil
   real(r_2), DIMENSION(3)                    :: totnwoodprod, totnlitter, totnsoil
   real(r_2), DIMENSION(3)                    :: totpwoodprod, totplitter, totpsoil
   real(r_2)                                     totclabile, totnsoilmin, totpsoillab, totpsoilsorb, totpsoilocc
   real(r_2)                                     tempx
   integer p,d,r,q,r1,r2,r3,r4,ierror,ivt,k
   integer irb,is,icp,ics

   ivt2=(/3,3,3,3,2,1,1,2,1,1,3,3,3,1,0,0,0/)
   delarea(:,:)     = 0.0
   dcplant(:,:,:)   = 0.0; dnplant(:,:,:)   = 0.0; dpplant(:,:,:)    = 0.0; dclabile(:,:) = 0.0
   dclitter(:,:,:)  = 0.0; dnlitter(:,:,:)  = 0.0; dplitter(:,:,:)   = 0.0
   dcsoil(:,:,:)    = 0.0; dnsoil(:,:,:)    = 0.0; dpsoil(:,:,:)     = 0.0
   dnsoilmin(:,:)   = 0.0; dpsoillab(:,:)   = 0.0; dpsoilsorb(:,:)   = 0.0; dpsoilocc(:,:) = 0.0
   dcwoodprod(:,:,:) =0.0; dnwoodprod(:,:,:)=0.0;  dpwoodprod(:,:,:) = 0.0

   do p = 1,mland
      fromPtoL(:,1,:)=1.0; fromPtoL(:,2,:) = 0.0
      do d=1,mvmax
         if(luc%cplant_x(p,d,leaf) > 0.001) then
            ! calculate the fraction of litter or root litter into metabolic litter pool
            ivt=mvtype

            ratioLignintoN(leaf) = (luc%cplant_x(p,d,leaf) &
                                 /(max(1.0e-10,luc%nplant_x(p,d,leaf)) * casabiome%ftransNPtoL(ivt,leaf))) &
                                 * casabiome%fracLigninplant(ivt,leaf)
            ratioLignintoN(froot)= (luc%cplant_x(p,d,froot)&
                                 /(max(1.0e-10,luc%nplant_x(p,d,froot))* casabiome%ftransNPtoL(ivt,froot))) &
                                 * casabiome%fracLigninplant(ivt,froot)
            fromPtoL(d,metb,leaf)  = max(0.001, 0.85 - 0.018 *ratioLignintoN(leaf))
            fromPtoL(d,metb,froot) = max(0.001, 0.85 - 0.018 *ratioLignintoN(froot))
            fromPtoL(d,str,leaf)   = 1.0 - fromPtoL(d,metb,leaf)
            fromPtoL(d,str,froot)  = 1.0 - fromPtoL(d,metb,froot)
         endif
      enddo

      ! compute the transition matrix relating to primary forest harvest
      delfwhpri(1:mvmax) = 0.0; afwhpri(1:mvmax,1:mvmax) = 0.0      
      if(luc%fharvw(p,1) >0.0) then
         ! calculate the transition matrix for primary land 
          do d=1,mvtype
            if(d<11.or.d>13) then
               delfwhpri(d)  = luc%fharvw(p,1) * luc%xluh2cable(p,d,1)    ! donor    (positive)
            else
               delfwhpri(d) = -luc%fharvw(p,1) * luc%xluh2cable(p,r,3)    ! receiver (negative)
            endif          
         enddo  ! of "d"
         call landuse_redistribution(p,mvmax,delfwhpri,afwhpri)  
      endif

      transitx(1:mvmax,1:mvmax) = luc%atransit(p,1:mvmax,1:mvmax)+afwhpri(1:mvmax,1:mvmax)

      do d = 1,mvmax
      do r = 1,mvmax
         ! transfer leaf and root into litter (metabolic and structural litter),
         ! transfer wood into wood prodoct pool (three wood product pools)
         ! then calculate the delpool(p,d,:), and delpool(p,r,:) for C, N and P

         if(transitx(r,d) > 0.0.and.d/=r.and.luc%patchfrac_x(p,d)>0.0) then
            ! transfer the area from donor (d) to receiver (r)
            delarea(p,d) = delarea(p,d) - transitx(r,d)
            delarea(p,r) = delarea(p,r) + transitx(r,d)
            ! donor pool changes
            dcplant(p,d,:)    = dcplant(p,d,:)    - transitx(r,d) * luc%cplant_x(p,d,:)
            dclitter(p,d,:)   = dclitter(p,d,:)   - transitx(r,d) * luc%clitter_x(p,d,:)
            dcsoil(p,d,:)     = dcsoil(p,d,:)     - transitx(r,d) * luc%csoil_x(p,d,:)
            dclabile(p,d)     = dclabile(p,d)     - transitx(r,d) * luc%clabile_x(p,d)
            dcwoodprod(p,d,:) = dcwoodprod(p,d,:) - transitx(r,d) * luc%cwoodprod_x(p,d,:)

            ! receiver pool changes
            ! move the donor leaf and root biomass to receiver structural litter pool
            !(to avoid C:N imbalance)
            ! using max function to avoid dividing by zero, ypw 14/may/2008

            ! calculate the fraction of litter or root litter into metabolic litter pool

            dclitter(p,r,1) = dclitter(p,r,1) + transitx(r,d)  &
                            * (luc%clitter_x(p,d,1) + luc%cplant_x(p,d,leaf) * fromPtoL(d,metb,leaf) &
                                                    + luc%cplant_x(p,d,froot)* fromPtoL(d,metb,froot))
            dclitter(p,r,2) = dclitter(p,r,2) + transitx(r,d)  &
                            * (luc%clitter_x(p,d,2) + luc%cplant_x(p,d,leaf) * fromPtoL(d,str,leaf)  &
                                                    + luc%cplant_x(p,d,froot) *fromPtoL(d,str,froot))
            dcsoil(p,r,:)   = dcsoil(p,r,:)   + transitx(r,d) * luc%csoil_x(p,d,:)

            !move the labile carbon to receiving tile (maybe better to fast-decomposing wood product pool)
            dclabile(p,r)   = dclabile(p,r)   + transitx(r,d) * luc%clabile_x(p,d)

            !move the donor wood to receiver wood product pool including labile C
            dcwoodprod(p,r,1) = dcwoodprod(p,r,1) + transitx(r,d)   &
                                                  *(luc%cplant_x(p,d,wood)*fwoodprod(1) + luc%cwoodprod_x(p,d,1)) 

            dcwoodprod(p,r,2) = dcwoodprod(p,r,2) + transitx(r,d)   &
                                                  *(luc%cplant_x(p,d,wood)*fwoodprod(2) + luc%cwoodprod_x(p,d,2))

            dcwoodprod(p,r,3) = dcwoodprod(p,r,3) + transitx(r,d)   &
                                                  *(luc%cplant_x(p,d,wood)*fwoodprod(3) + luc%cwoodprod_x(p,d,3))

             if(icycle >1) then
                dnplant(p,d,:)  = dnplant(p,d,:)  - transitx(r,d) * luc%nplant_x(p,d,:)
                dnlitter(p,d,:) = dnlitter(p,d,:) - transitx(r,d) * luc%nlitter_x(p,d,:)
                dnsoil(p,d,:)   = dnsoil(p,d,:)   - transitx(r,d) * luc%nsoil_x(p,d,:)
                dnsoilmin(p,d)  = dnsoilmin(p,d)  - transitx(r,d) * luc%nsoilmin_x(p,d)

                dnlitter(p,r,1) = dnlitter(p,r,1) + transitx(r,d)  &
                                * (luc%nlitter_x(p,d,1) + luc%nplant_x(p,d,leaf) * fromPtoL(d,metb,leaf) &
                                                        + luc%nplant_x(p,d,froot)* fromPtoL(d,metb,froot))
                dnlitter(p,r,2) = dnlitter(p,r,2) + transitx(r,d)  &
                                * (luc%nlitter_x(p,d,2) + luc%nplant_x(p,d,leaf) * fromPtoL(d,str,leaf) &
                                                        + luc%nplant_x(p,d,froot) *fromPtoL(d,str,froot))
                dnsoil(p,r,:)   = dnsoil(p,r,:)   + transitx(r,d) * luc%nsoil_x(p,d,:)
                dnsoilmin(p,r)  = dnsoilmin(p,r)  + transitx(r,d) * luc%nsoilmin_x(p,d)
                dnwoodprod(p,r,:) = dnwoodprod(p,r,:) + transitx(r,d) &
                                                      *(fwoodprod(:)*luc%nplant_x(p,d,wood) + luc%nwoodprod_x(p,d,:))
             endif
             if(icycle >2) then
                dpplant(p,d,:)  = dpplant(p,d,:)  - transitx(r,d) * luc%pplant_x(p,d,:)
                dplitter(p,d,:) = dplitter(p,d,:) - transitx(r,d) * luc%nlitter_x(p,d,:)
                dpsoil(p,d,:)   = dpsoil(p,d,:)   - transitx(r,d) * luc%psoil_x(p,d,:)
                dpsoillab(p,d)  = dpsoillab(p,d)  - transitx(r,d) * luc%psoillab_x(p,d)
                dpsoilsorb(p,d) = dpsoilsorb(p,d) - transitx(r,d) * luc%psoilsorb_x(p,d)
                dpsoilocc(p,d)  = dpsoilocc(p,d)  - transitx(r,d) * luc%psoilocc_x(p,d)

                dplitter(p,r,1) = dplitter(p,r,1) + transitx(r,d)  &
                                * (luc%plitter_x(p,d,1) + luc%pplant_x(p,d,leaf) * fromPtoL(d,metb,leaf) &
                                                        + luc%pplant_x(p,d,froot)* fromPtoL(d,metb,froot))
                dplitter(p,r,2) = dplitter(p,r,2) + transitx(r,d)  &
                                * (luc%plitter_x(p,d,2) + luc%pplant_x(p,d,leaf) * fromPtoL(d,str,leaf) &
                                                        + luc%pplant_x(p,d,froot) *fromPtoL(d,str,froot))
                dpsoil(p,r,:)   = dpsoil(p,r,:)   + transitx(r,d) * luc%psoil_x(p,d,:)
                dpsoillab(p,r)  = dpsoillab(p,r)  + transitx(r,d) * luc%psoillab_x(p,d)
                dpsoilsorb(p,r) = dpsoilsorb(p,r) + transitx(r,d) * luc%psoilsorb_x(p,d)
                dpsoilocc(p,r)  = dpsoilocc(p,r)  + transitx(r,d) * luc%psoilocc_x(p,d)
                dpwoodprod(p,r,:) = dpwoodprod(p,r,:) + transitx(r,d)   &
                                                      * (fwoodprod(:) * luc%pplant_x(p,d,wood) +luc%pwoodprod_x(p,d,:))
             endif
         endif  ! of "atransit" >0.0, "luctype" >0 etc.
      enddo   ! of "d_tile"
      enddo   ! of "r_tile"

      !@@@@ here we deal with wood harvest from secondary forest (not done yet)

      luc%patchfrac_y(p,:)   = luc%patchfrac_x(p,:) + delarea(p,:)

      do d=1,mvmax
      luc%patchfrac_y(p,d)   = luc%patchfrac_x(p,d) + delarea(p,d)
      if(luc%patchfrac_y(p,d)>0.0) then
         luc%cplant_y(p,d,leaf)  = (luc%cplant_x(p,d,leaf) * luc%patchfrac_x(p,d) + dcplant(p,d,leaf))  &
                                    /luc%patchfrac_y(p,d)
         luc%cplant_y(p,d,wood)  = (luc%cplant_x(p,d,wood) * luc%patchfrac_x(p,d) + dcplant(p,d,wood))  &
                                    /luc%patchfrac_y(p,d)
         luc%cplant_y(p,d,froot) = (luc%cplant_x(p,d,froot)* luc%patchfrac_x(p,d) + dcplant(p,d,froot)) &
                                    /luc%patchfrac_y(p,d)
         luc%clabile_y(p,d)      = (luc%clabile_x(p,d)     * luc%patchfrac_x(p,d) + dclabile(p,d))      &
                                    /luc%patchfrac_y(p,d)
         luc%clitter_y(p,d,metb)  = (luc%clitter_x(p,d,metb) * luc%patchfrac_x(p,d) + dclitter(p,d,metb))  &
                                    /luc%patchfrac_y(p,d)
         luc%clitter_y(p,d,str)  = (luc%clitter_x(p,d,str) * luc%patchfrac_x(p,d) + dclitter(p,d,str))  &
                                    /luc%patchfrac_y(p,d)
         luc%clitter_y(p,d,cwd)  = (luc%clitter_x(p,d,cwd) * luc%patchfrac_x(p,d) + dclitter(p,d,cwd))  &
                                    /luc%patchfrac_y(p,d)
         luc%csoil_y(p,d,mic)    = (luc%csoil_x(p,d,mic)   * luc%patchfrac_x(p,d) + dcsoil(p,d,mic))    &
                                    /luc%patchfrac_y(p,d)
         luc%csoil_y(p,d,slow)   = (luc%csoil_x(p,d,slow)  * luc%patchfrac_x(p,d) + dcsoil(p,d,slow))   &
                                    /luc%patchfrac_y(p,d)
         luc%csoil_y(p,d,3)      = (luc%csoil_x(p,d,3)     * luc%patchfrac_x(p,d) + dcsoil(p,d,3))      &
                                    /luc%patchfrac_y(p,d)
         luc%cwoodprod_y(p,d,1)  = (luc%cwoodprod_x(p,d,1) * luc%patchfrac_x(p,d) + dcwoodprod(p,d,1))  &
                                    /luc%patchfrac_y(p,d)
         luc%cwoodprod_y(p,d,2)  = (luc%cwoodprod_x(p,d,2) * luc%patchfrac_x(p,d) + dcwoodprod(p,d,2))  &
                                    /luc%patchfrac_y(p,d)
         luc%cwoodprod_y(p,d,3)  = (luc%cwoodprod_x(p,d,3) * luc%patchfrac_x(p,d) + dcwoodprod(p,d,3))  &
                                    /luc%patchfrac_y(p,d)

         if(icycle >1) then
            luc%nplant_y(p,d,leaf)  = (luc%nplant_x(p,d,leaf) * luc%patchfrac_x(p,d) + dnplant(p,d,leaf))  &
                                    /luc%patchfrac_y(p,d)
            luc%nplant_y(p,d,wood)  = (luc%nplant_x(p,d,wood) * luc%patchfrac_x(p,d) + dnplant(p,d,wood))  &
                                    /luc%patchfrac_y(p,d)
            luc%nplant_y(p,d,froot) = (luc%nplant_x(p,d,froot)* luc%patchfrac_x(p,d) + dnplant(p,d,froot)) &
                                    /luc%patchfrac_y(p,d)
            luc%nlitter_y(p,d,metb) = (luc%nlitter_x(p,d,metb) * luc%patchfrac_x(p,d) + dnlitter(p,d,metb))  &
                                    /luc%patchfrac_y(p,d)
            luc%nlitter_y(p,d,str)  = (luc%nlitter_x(p,d,str) * luc%patchfrac_x(p,d) + dnlitter(p,d,str))  &
                                    /luc%patchfrac_y(p,d)
            luc%nlitter_y(p,d,cwd)  = (luc%nlitter_x(p,d,cwd) * luc%patchfrac_x(p,d) + dnlitter(p,d,cwd))  &
                                    /luc%patchfrac_y(p,d)
            luc%nsoil_y(p,d,mic)    = (luc%nsoil_x(p,d,mic)   * luc%patchfrac_x(p,d) + dnsoil(p,d,mic))    &
                                    /luc%patchfrac_y(p,d)
            luc%nsoil_y(p,d,slow)   = (luc%nsoil_x(p,d,slow)  * luc%patchfrac_x(p,d) + dnsoil(p,d,slow))   &
                                    /luc%patchfrac_y(p,d)
            luc%nsoil_y(p,d,3)      = (luc%nsoil_x(p,d,3)     * luc%patchfrac_x(p,d) + dnsoil(p,d,3))      &
                                    /luc%patchfrac_y(p,d)
            luc%nsoilmin_y(p,d)     = (luc%nsoilmin_x(p,d)    * luc%patchfrac_x(p,d) + dnsoilmin(p,d))      &
                                    /luc%patchfrac_y(p,d)

            luc%nwoodprod_y(p,d,1)  = (luc%nwoodprod_x(p,d,1) * luc%patchfrac_x(p,d) + dnwoodprod(p,d,1))  &
                                    /luc%patchfrac_y(p,d)
            luc%nwoodprod_y(p,d,2)  = (luc%nwoodprod_x(p,d,2) * luc%patchfrac_x(p,d) + dnwoodprod(p,d,2))  &
                                    /luc%patchfrac_y(p,d)
            luc%nwoodprod_y(p,d,3)  = (luc%nwoodprod_x(p,d,3) * luc%patchfrac_x(p,d) + dnwoodprod(p,d,3))  &
                                    /luc%patchfrac_y(p,d)
         endif
         if(icycle >2) then
            luc%pplant_y(p,d,leaf)  = (luc%pplant_x(p,d,leaf) * luc%patchfrac_x(p,d) + dpplant(p,d,leaf))  &
                                    /luc%patchfrac_y(p,d)
            luc%pplant_y(p,d,wood)  = (luc%pplant_x(p,d,wood) * luc%patchfrac_x(p,d) + dpplant(p,d,wood))  &
                                    /luc%patchfrac_y(p,d)
            luc%pplant_y(p,d,froot) = (luc%pplant_x(p,d,froot)* luc%patchfrac_x(p,d) + dpplant(p,d,froot)) &
                                    /luc%patchfrac_y(p,d)
            luc%plitter_y(p,d,metb)  = (luc%plitter_x(p,d,metb) * luc%patchfrac_x(p,d) + dplitter(p,d,metb))  &
                                    /luc%patchfrac_y(p,d)
            luc%plitter_y(p,d,str)  = (luc%plitter_x(p,d,str) * luc%patchfrac_x(p,d) + dplitter(p,d,str))  &
                                    /luc%patchfrac_y(p,d)
            luc%plitter_y(p,d,cwd)  = (luc%plitter_x(p,d,cwd) * luc%patchfrac_x(p,d) + dplitter(p,d,cwd))  &
                                    /luc%patchfrac_y(p,d)
            luc%psoil_y(p,d,mic)    = (luc%psoil_x(p,d,mic)   * luc%patchfrac_x(p,d) + dpsoil(p,d,mic))    &
                                    /luc%patchfrac_y(p,d)
            luc%psoil_y(p,d,slow)   = (luc%psoil_x(p,d,slow)  * luc%patchfrac_x(p,d) + dpsoil(p,d,slow))   &
                                    /luc%patchfrac_y(p,d)
            luc%psoil_y(p,d,3)      = (luc%psoil_x(p,d,3)     * luc%patchfrac_x(p,d) + dpsoil(p,d,3))      &
                                    /luc%patchfrac_y(p,d)
            luc%psoillab_y(p,d)     = (luc%psoillab_x(p,d)    * luc%patchfrac_x(p,d) + dpsoillab(p,d))     &
                                    /luc%patchfrac_y(p,d)
            luc%psoilsorb_y(p,d)    = (luc%psoilsorb_x(p,d)   * luc%patchfrac_x(p,d) + dpsoilsorb(p,d))    &
                                    /luc%patchfrac_y(p,d)
            luc%psoilocc_y(p,d)     = (luc%psoilocc_x(p,d)    * luc%patchfrac_x(p,d) + dpsoilocc(p,d))     &
                                    /luc%patchfrac_y(p,d)

            luc%pwoodprod_y(p,d,1)  = (luc%pwoodprod_x(p,d,1) * luc%patchfrac_x(p,d) + dpwoodprod(p,d,1))  &
                                    /luc%patchfrac_y(p,d)
            luc%pwoodprod_y(p,d,2)  = (luc%pwoodprod_x(p,d,2) * luc%patchfrac_x(p,d) + dpwoodprod(p,d,2))  &
                                    /luc%patchfrac_y(p,d)
            luc%pwoodprod_y(p,d,3)  = (luc%pwoodprod_x(p,d,3) * luc%patchfrac_x(p,d) + dpwoodprod(p,d,3))  &
                                    /luc%patchfrac_y(p,d)
         endif

      endif   ! luc%patchfrac_y(p,d)>0.0
      enddo   ! of "d"


      !update the biophysical variables here
      ! "isoil" and "soilorder" are contant within a landcell (see gridinfo)
      luc%isoil_y(p,1)     = dominantx(1,12,luc%patchfrac_x(p,1:mvmax),luc%isoil_x(p,1:mvmax))
      luc%soilorder_y(p,1) = dominantx(1,12,luc%patchfrac_x(p,1:mvmax),luc%soilorder_x(p,1:mvmax))
      ! may need to differentiate tree, grass and crop
      luc%phase_y(p,1)     = dominantx(0,3,luc%patchfrac_x(p,1:mvmax),luc%phase_x(p,1:mvmax))
      luc%doyphase3_y(p,1) = dominantx(0,365,luc%patchfrac_x(p,1:mvmax),luc%doyphase3_x(p,1:mvmax))
      do d=2,mvmax
         luc%isoil_y(p,d)     = luc%isoil_y(p,1)
         luc%soilorder_y(p,d) = luc%soilorder_y(p,1)
         luc%phase_y(p,d)     = luc%phase_y(p,1)
         luc%doyphase3_y(p,d) = luc%doyphase3_y(p,1)
      enddo

      do d=1,mvmax

         luc%iveg_y(p,d) = d

         luc%phen_y(p,d)  = avgpatchr2(d,luc%patchfrac_x(p,d),transitx(1:mvmax,1:mvmax),luc%phen_x(p,1:mvmax))
         luc%aphen_y(p,d) = avgpatchr2(d,luc%patchfrac_x(p,d),transitx(1:mvmax,1:mvmax),luc%aphen_x(p,1:mvmax))
         luc%frac_sapwood_y(p,d) =avgpatchr2(d,luc%patchfrac_x(p,d),transitx(1:mvmax,1:mvmax),luc%frac_sapwood_x(p,1:mvmax))
         luc%sapwood_area_y(p,d) =avgpatchr2(d,luc%patchfrac_x(p,d),transitx(1:mvmax,1:mvmax),luc%sapwood_area_x(p,1:mvmax))

         luc%sla_y(p,d) = luc%sla_x(p,d)
         luc%lai_y(p,d) = luc%sla_y(p,d) * max(0.0,luc%cplant_y(p,d,leaf))

         do irb=1,nrb
            luc%albsoilsn_y(p,d,irb)= avgpatchr2(d,luc%patchfrac_x(p,d),transitx(1:mvmax,1:mvmax),luc%albsoilsn_x(p,1:mvmax,irb))
            luc%albedo_y(p,d,irb)   = avgpatchr2(d,luc%patchfrac_x(p,d),transitx(1:mvmax,1:mvmax),luc%albedo_x(p,1:mvmax,irb))
            luc%albsoil_y(p,d,irb)  = avgpatchr2(d,luc%patchfrac_x(p,d),transitx(1:mvmax,1:mvmax),luc%albsoil_x(p,1:mvmax,irb))
         enddo

         luc%dgdtg_y(p,d) = avgpatchr2(d,luc%patchfrac_x(p,d),transitx(1:mvmax,1:mvmax),luc%dgdtg_x(p,1:mvmax))

         do is=1,ms
            luc%gammzz_y(p,d,is) = avgpatchr2(d,luc%patchfrac_x(p,d),transitx(1:mvmax,1:mvmax),luc%gammzz_x(p,1:mvmax,is))
            luc%tgg_y(p,d,is)    = avgpatchr2(d,luc%patchfrac_x(p,d),transitx(1:mvmax,1:mvmax),luc%tgg_x(p,1:mvmax,is))
            luc%wb_y(p,d,is)     = avgpatchr2(d,luc%patchfrac_x(p,d),transitx(1:mvmax,1:mvmax),luc%wb_x(p,1:mvmax,is))
            luc%wbice_y(p,d,is)  = avgpatchr2(d,luc%patchfrac_x(p,d),transitx(1:mvmax,1:mvmax),luc%wbice_x(p,1:mvmax,is))
         enddo

         do is=1,msn
            luc%tggsn_y(p,d,is)  = avgpatchr2(d,luc%patchfrac_x(p,d),transitx(1:mvmax,1:mvmax),luc%tggsn_x(p,1:mvmax,is))
            luc%ssdn_y(p,d,is)   = avgpatchr2(d,luc%patchfrac_x(p,d),transitx(1:mvmax,1:mvmax),luc%ssdn_x(p,1:mvmax,is))
            luc%smass_y(p,d,is)  = avgpatchr2(d,luc%patchfrac_x(p,d),transitx(1:mvmax,1:mvmax),luc%smass_x(p,1:mvmax,is))
            luc%sdepth_y(p,d,is) = avgpatchr2(d,luc%patchfrac_x(p,d),transitx(1:mvmax,1:mvmax),luc%sdepth_x(p,1:mvmax,is))
         enddo

         if(luc%smass_y(p,d,1)<=0.0) then
            luc%isflag_y(p,d) = 0
         else
            luc%isflag_y(p,d) = 1
         endif

         luc%tss_y(p,d)       = avgpatchr2(d,luc%patchfrac_x(p,d),transitx(1:mvmax,1:mvmax),luc%tss_x(p,1:mvmax))
         luc%rtsoil_y(p,d)    = avgpatchr2(d,luc%patchfrac_x(p,d),transitx(1:mvmax,1:mvmax),luc%rtsoil_x(p,1:mvmax))
         luc%runoff_y(p,d)    = avgpatchr2(d,luc%patchfrac_x(p,d),transitx(1:mvmax,1:mvmax),luc%runoff_x(p,1:mvmax))
         luc%rnof1_y(p,d)     = avgpatchr2(d,luc%patchfrac_x(p,d),transitx(1:mvmax,1:mvmax),luc%rnof1_x(p,1:mvmax))
         luc%rnof2_y(p,d)     = avgpatchr2(d,luc%patchfrac_x(p,d),transitx(1:mvmax,1:mvmax),luc%rnof2_x(p,1:mvmax))
         luc%ssdnn_y(p,d)     = avgpatchr2(d,luc%patchfrac_x(p,d),transitx(1:mvmax,1:mvmax),luc%ssdnn_x(p,1:mvmax))
         luc%snowd_y(p,d)     = avgpatchr2(d,luc%patchfrac_x(p,d),transitx(1:mvmax,1:mvmax),luc%snowd_x(p,1:mvmax))
         luc%snage_y(p,d)     = avgpatchr2(d,luc%patchfrac_x(p,d),transitx(1:mvmax,1:mvmax),luc%snage_x(p,1:mvmax))
         luc%osnowd_y(p,d)    = avgpatchr2(d,luc%patchfrac_x(p,d),transitx(1:mvmax,1:mvmax),luc%osnowd_x(p,1:mvmax))
         luc%cansto_y(p,d)    = avgpatchr2(d,luc%patchfrac_x(p,d),transitx(1:mvmax,1:mvmax),luc%cansto_x(p,1:mvmax))
         luc%ghflux_y(p,d)    = avgpatchr2(d,luc%patchfrac_x(p,d),transitx(1:mvmax,1:mvmax),luc%ghflux_x(p,1:mvmax))
         luc%sghflux_y(p,d)   = avgpatchr2(d,luc%patchfrac_x(p,d),transitx(1:mvmax,1:mvmax),luc%sghflux_x(p,1:mvmax))
         luc%ga_y(p,d)        = avgpatchr2(d,luc%patchfrac_x(p,d),transitx(1:mvmax,1:mvmax),luc%ga_x(p,1:mvmax))
         luc%fev_y(p,d)       = avgpatchr2(d,luc%patchfrac_x(p,d),transitx(1:mvmax,1:mvmax),luc%fev_x(p,1:mvmax))
         luc%fes_y(p,d)       = avgpatchr2(d,luc%patchfrac_x(p,d),transitx(1:mvmax,1:mvmax),luc%fes_x(p,1:mvmax))
         luc%fhs_y(p,d)       = avgpatchr2(d,luc%patchfrac_x(p,d),transitx(1:mvmax,1:mvmax),luc%fhs_x(p,1:mvmax))
         luc%wbtot0_y(p,d)    = avgpatchr2(d,luc%patchfrac_x(p,d),transitx(1:mvmax,1:mvmax),luc%wbtot0_x(p,1:mvmax))
         luc%osnowd0_y(p,d)   = avgpatchr2(d,luc%patchfrac_x(p,d),transitx(1:mvmax,1:mvmax),luc%osnowd0_x(p,1:mvmax)) 
         luc%trad_y(p,d)      = avgpatchr2(d,luc%patchfrac_x(p,d),transitx(1:mvmax,1:mvmax),luc%trad_x(p,1:mvmax))
         luc%GWwb_y(p,d)      = avgpatchr2(d,luc%patchfrac_x(p,d),transitx(1:mvmax,1:mvmax),luc%GWwb_x(p,1:mvmax))

         do icp=1,ncp
            luc%cplantx_y(p,d,icp) = avgpatchr2(d,luc%patchfrac_x(p,d),transitx(1:mvmax,1:mvmax),luc%cplantx_x(p,1:mvmax,icp))
         enddo

         do ics=1,ncs
            luc%csoilx_y(p,d,ics)  = avgpatchr2(d,luc%patchfrac_x(p,d),transitx(1:mvmax,1:mvmax),luc%csoilx_x(p,1:mvmax,ics))
         enddo

      enddo   ! of "d"


      ! adding seeding biomass if area > critical value but biomass is too low
      do d=1,mvmax
         ivt = d
         if(luc%patchfrac_x(p,d)<thresh_frac.and.luc%patchfrac_y(p,d)>thresh_frac) then  !newly-born patch
            if(d==1.or.d==2) then
               luc%phase_y(p,d) = 2
            else
               luc%phase_y(p,d) = 0
            endif
         else
            luc%phase_y(p,d) = luc%phase_x(p,d)
         endif

         if(luc%patchfrac_y(p,d) > thresh_frac.and.sum(luc%cplant_y(p,d,1:3))<fseedling) then
            if(ivt2(ivt)==1) then
               luc%cplant_y(p,d,:) = fseedling * fracgrassseed(:)
            else
               luc%cplant_y(p,d,:) = fseedling * fracwoodseed(:)
            endif
            if(icycle >1) then
               luc%nplant_y(p,d,:) = luc%cplant_y(p,d,:) * casabiome%ratioNCplantmax(ivt,:)
            endif
            if(icycle >2) then
               luc%pplant_y(p,d,:) = luc%nplant_y(p,d,:) /casabiome%ratioNPplantmin(ivt,:)
            endif
         endif
      enddo  ! of "d"

    enddo    ! of "p"

    CONTAINS
      function avgpatchr2(q,areax,x2y,x) result(avgr2)
      ! check note: workbook, 2017, p121
      integer q, k1, k2
      real(r_2)                         :: avgr2
      real(r_2)                         :: areax
      real(r_2)                         :: xloss,xgain,delareax
      real(r_2), dimension(mvmax,mvmax) :: x2y
      real(r_2), dimension(mvmax)       :: x

        avgr2=x(q)
        do k1=1,mvmax   ! loss
           if(k1/=q) then
              xloss = xloss + x2y(k1,q) * x(q)
              delareax = delareax - x2y(k1,q)
           endif
        enddo
 
        do k2=1,mvmax
           if(k2/=q) then
              xgain = xgain +x2y(q,k2)  * x(k2)
              delareax = delareax + x2y(q,k2)
           endif
        enddo
        if((areax+delareax) > thresh_frac) then
           avgr2= (x(q)*areax +xgain-xloss)/(areax+delareax)
        endif
      end function avgpatchr2


      function dominantx(xmin,xmax,fracx,xk) result(dominantint)
      ! take the value for the patch with maximum area fraction
      ! check note: workbook, 2017, p121
      integer dominantint
      integer xmin,xmax,k
      integer,   dimension(mvmax)       :: xk
      real(r_2), dimension(mvmax)       :: fracx
      real(r_2)  fracxmax

        dominantint = min(xmax,max(xmin,xk(1)))
        fracxmax  = fracx(1)
        do k=2,mvmax   
           if(fracx(k)>fracxmax) then
              fracxmax = fracx(k)
              dominantint = xk(k)
           endif
        enddo
      end function dominantx


END SUBROUTINE landuse_transitx

   SUBROUTINE landuse_redistribution(p,mvmax,delfracx,atransx)
      ! redistribution the PFT atrsnition to ensure that the change in PFT fractions from the state
      ! data is consistent with the estimates from previous states and transition
      USE cable_def_types_mod,    ONLY: r_2
      implicit none
      real,    parameter                          :: thresh_frac=1.0e-6
      integer p,mvmax
      real(r_2),    dimension(mvmax,mvmax)  :: atransx
      ! local variables
      integer itemp,i,j,k,ndonor,nreceive,vi,vj,np
      real(r_2)     temp
      real(r_2),    dimension(mvmax,mvmax)  :: transx,transy
      real(r_2),    dimension(mvmax)        :: delx,delfracx,donor,receive
      integer,          dimension(mvmax)    :: ivt,ivtdonor,ivtreceive

       atransx(:,:)   = 0.0;     delx=delfracx

       ! check the sum is zero
       if(abs(sum(delx(1:mvmax)))>thresh_frac) then
          print *, 'unbalanced ', p,delx(1:mvmax),sum(delx(1:mvmax))
       endif

       do i=1,mvmax
          ivt(i) = i
       enddo

       ! sort data from the smallest to the largest
       ! receiveer: negative; donor: positive
       do j=mvmax-1,1,-1
          do i=1,j
             if(delx(i) < delx(i+1)) then
                temp      = delx(i)
                itemp     = ivt(i)
                delx(i)   = delx(i+1)
                ivt(i)    = ivt(i+1) 
                delx(i+1) = temp
                ivt(i+1)  = itemp
             endif
          enddo
       enddo 

       ! determine number of negative (receive) and positives (donor)
       ndonor=0;nreceive=0
       do i=1,mvmax
          if(delx(i) >0.0) ndonor = ndonor+1
          if(delx(i) <0.0) nreceive=nreceive+1
       enddo

       donor(1:ndonor)    = delx(1:ndonor)
       ivtdonor(1:ndonor) = ivt(1:ndonor)

       do i=mvmax,mvmax-nreceive+1,-1
          j=mvmax-i+1
          receive(j)    = delx(i)
          ivtreceive(j) = ivt(i)
       enddo
         
       ! donor to receive
       do i=1,ndonor
          if(donor(i) > 0.0) then
             do j=1,nreceive
                if(receive(j) <0.0) then
                   vi = ivtdonor(i)
                   vj = ivtreceive(j)
                   if(donor(i)>=-receive(j)) then
                      atransx(vj,vi) = atransx(vj,vi) - receive(j)
                      donor(i)       = donor(i) + receive(j)        
                      receive(j)     = 0.0 
                   else
                      atransx(vj,vi) = atransx(vj,vi) + donor(i)
                      receive(j)     = receive(j) + donor(i)
                      donor(i)       = 0.0
                   endif
                endif
             enddo  ! "j"
          endif
       enddo    ! "i"


       ! verify the results
       delx = delfracx
       do j=1,mvmax    ! primary land
       do i=11,13,1    ! secondary forest
          delx(j) = delx(j) - atransx(i,j)
          delx(i) = delx(i) + atransx(i,j)
       enddo
       enddo

       do i=1,mvmax
          if(delx(i) > thresh_frac) then
             print *, 'warning: landuse_redistribution faile ', p,i,delx(i)
          endif
       enddo

 END SUBROUTINE landuse_redistribution

 SUBROUTINE landuse_update_mland(luc)                     ! assign "var_y" to "var_x"
 USE landuse_variable,   ONLY: landuse_mland
 IMPLICIT NONE
 TYPE(landuse_mland) :: luc

    ! general patch variables
    luc%iveg_x      = luc%iveg_y
    luc%isoil_x     = luc%isoil_y
    luc%soilorder_x = luc%soilorder_y
    luc%phase_x     = luc%phase_y

    luc%doyphase3_x = luc%doyphase3_y
    luc%phen_x      = luc%phen_y
    luc%aphen_x     = luc%aphen_y
    luc%frac_sapwood_x=luc%frac_sapwood_y
    luc%sapwood_area_x=luc%sapwood_area_y

    luc%isflag_x    = luc%isflag_y
    luc%patchfrac_x = luc%patchfrac_y
    luc%lai_x       = luc%lai_y
    luc%sla_x       = luc%sla_y

    ! biophysical
    luc%albsoilsn_x  = luc%albsoilsn_y
    luc%albedo_x     = luc%albedo_y
    luc%albsoil_x    = luc%albsoil_y
    luc%dgdtg_x      = luc%dgdtg_y
    luc%gammzz_x     = luc%gammzz_y
    luc%tgg_x        = luc%tgg_y
    luc%wb_x         = luc%wb_y
    luc%wbice_x      = luc%wbice_y
    luc%tggsn_x      = luc%tggsn_y
    luc%ssdn_x       = luc%ssdn_y
    luc%smass_x      = luc%smass_y
    luc%sdepth_x     = luc%sdepth_y
    luc%tss_x        = luc%tss_y
    luc%rtsoil_x     = luc%rtsoil_y
    luc%runoff_x     = luc%runoff_y
    luc%rnof1_x      = luc%rnof1_y
    luc%rnof2_x      = luc%rnof2_y
    luc%ssdnn_x      = luc%ssdnn_y
    luc%snowd_x      = luc%snowd_y
    luc%snage_x      = luc%snage_y
    luc%osnowd_x     = luc%osnowd_y
    luc%cansto_x     = luc%cansto_y
    luc%ghflux_x     = luc%ghflux_y
    luc%sghflux_x    = luc%sghflux_y
    luc%ga_x         = luc%ga_y
    luc%fev_x        = luc%fev_y
    luc%fes_x        = luc%fes_y
    luc%fhs_x        = luc%fhs_y
    luc%wbtot0_x     = luc%wbtot0_y
    luc%osnowd0_x    = luc%osnowd0_y
    luc%trad_x       = luc%trad_y
    luc%GWwb_x       = luc%GWwb_y
    luc%cplantx_x    = luc%cplantx_y
    luc%csoilx_x     = luc%csoilx_y

    ! biogeochemical variables
    luc%cplant_x    = luc%cplant_y
    luc%nplant_x    = luc%nplant_y
    luc%pplant_x    = luc%pplant_y
    luc%clitter_x   = luc%clitter_y
    luc%nlitter_x   = luc%nlitter_y
    luc%plitter_x   = luc%plitter_y
    luc%csoil_x     = luc%csoil_y
    luc%nsoil_x     = luc%nsoil_y
    luc%psoil_x     = luc%psoil_y
    luc%clabile_x   = luc%clabile_y
    luc%nsoilmin_x  = luc%nsoilmin_y
    luc%psoillab_x  = luc%psoillab_y
    luc%psoilsorb_x = luc%psoilsorb_y
    luc%psoilocc_x  = luc%psoilocc_y
    luc%cwoodprod_x = luc%cwoodprod_y
    luc%nwoodprod_x = luc%nwoodprod_y
    luc%pwoodprod_x = luc%pwoodprod_y

 END SUBROUTINE landuse_update_mland

 SUBROUTINE landuse_land2mpx(luc,lucmp,mpx)
 USE landuse_constant,     ONLY: mvmax
 USE landuse_variable
 USE cable_def_types_mod,  ONLY: mland
 IMPLICIT NONE
 TYPE(landuse_mland)         :: luc
 TYPE(landuse_mp)            :: lucmp
! integer, dimension(mland)   :: cstart,cend,nap
 integer mpx
 integer np,np1,p,q,n,npnew,npold

    npnew=0; npold=0
    do p=1,mland
       do q=1,mvmax
          if(luc%patchfrac_x(p,q)>thresh_frac) then
             npold=npold +1
          endif
          if(luc%patchfrac_y(p,q)>thresh_frac) then
             npnew = npnew +1
             lucmp%iveg(npnew)      = q
             lucmp%isoil(npnew)     = luc%isoil_y(p,q)
             lucmp%soilorder(npnew) = luc%soilorder_y(p,q)
             lucmp%phase(npnew)      = luc%phase_y(p,q)
             lucmp%isflag(npnew)     = luc%isflag_y(p,q)
             lucmp%patchfrac(npnew)  = luc%patchfrac_y(p,q)
             lucmp%lai(npnew)        = luc%lai_y(p,q)
             lucmp%sla(npnew)        = luc%sla_y(p,q)

             lucmp%doyphase3(npnew) = luc%doyphase3_y(p,q)
             lucmp%phen(npnew)      = luc%phen_y(p,q)
             lucmp%aphen(npnew)     = luc%aphen_y(p,q)
             lucmp%frac_sapwood(npnew) =luc%frac_sapwood_y(p,q)
             lucmp%sapwood_area(npnew) =luc%sapwood_area_y(p,q)

             ! biophysical
             lucmp%albsoilsn(npnew,:)= luc%albsoilsn_y(p,q,:)
             lucmp%albedo(npnew,:)   = luc%albedo_y(p,q,:)
             lucmp%albsoil(npnew,:)  = luc%albsoil_y(p,q,:)
             lucmp%dgdtg(npnew)      = luc%dgdtg_y(p,q)
             lucmp%gammzz(npnew,:)   = luc%gammzz_y(p,q,:)
             lucmp%tgg(npnew,:)      = luc%tgg_y(p,q,:)
             lucmp%wb(npnew,:)       = luc%wb_y(p,q,:)
             lucmp%wbice(npnew,:)    = luc%wbice_y(p,q,:)
             lucmp%tggsn(npnew,:)    = luc%tggsn_y(p,q,:)
             lucmp%ssdn(npnew,:)     = luc%ssdn_y(p,q,:)
             lucmp%smass(npnew,:)    = luc%smass_y(p,q,:)
             lucmp%sdepth(npnew,:)   = luc%sdepth_y(p,q,:)
             lucmp%tss(npnew)        = luc%tss_y(p,q)
             lucmp%rtsoil(npnew)     = luc%rtsoil_y(p,q)
             lucmp%runoff(npnew)     = luc%runoff_y(p,q)
             lucmp%rnof1(npnew)      = luc%rnof1_y(p,q)
             lucmp%rnof2(npnew)      = luc%rnof2_y(p,q)
             lucmp%ssdnn(npnew)      = luc%ssdnn_y(p,q)
             lucmp%snowd(npnew)      = luc%snowd_y(p,q)
             lucmp%snage(npnew)      = luc%snage_y(p,q)
             lucmp%osnowd(npnew)     = luc%osnowd_y(p,q)
             lucmp%cansto(npnew)     = luc%cansto_y(p,q)
             lucmp%ghflux(npnew)     = luc%ghflux_y(p,q)
             lucmp%sghflux(npnew)    = luc%sghflux_y(p,q)
             lucmp%ga(npnew)         = luc%ga_y(p,q)
             lucmp%fev(npnew)        = luc%fev_y(p,q)
             lucmp%fes(npnew)        = luc%fes_y(p,q)
             lucmp%fhs(npnew)        = luc%fhs_y(p,q)
             lucmp%wbtot0(npnew)     = luc%wbtot0_y(p,q)
             lucmp%osnowd0(npnew)    = luc%osnowd0_y(p,q)
             lucmp%trad(npnew)       = luc%trad_y(p,q)
             lucmp%GWwb(npnew)       = luc%GWwb_y(p,q)
             lucmp%cplantx(npnew,:)  = luc%cplantx_y(p,q,:)
             lucmp%csoilx(npnew,:)   = luc%csoilx_y(p,q,:)

             ! assign the new biogeochemocal state variables
             if(icycle > 0) then
                lucmp%cplant(npnew,:)    = luc%cplant_y(p,q,:)
                lucmp%clitter(npnew,:)   = luc%clitter_y(p,q,:)
                lucmp%csoil(npnew,:)     = luc%csoil_y(p,q,:)
                lucmp%clabile(npnew)     = luc%clabile_y(p,q)
                lucmp%cwoodprod(npnew,:) = luc%cwoodprod_y(p,q,:)
             endif
             if(icycle >1) then
                lucmp%nplant(npnew,:)    = luc%nplant_y(p,q,:)
                lucmp%nlitter(npnew,:)   = luc%nlitter_y(p,q,:)
                lucmp%nsoil(npnew,:)     = luc%nsoil_y(p,q,:)
                lucmp%nsoilmin(npnew)    = luc%nsoilmin_y(p,q)
                lucmp%nwoodprod(npnew,:) = luc%nwoodprod_y(p,q,:)
             endif
             if(icycle >2) then
                lucmp%pplant(npnew,:)    = luc%pplant_y(p,q,:)
                lucmp%plitter(npnew,:)   = luc%plitter_y(p,q,:)
                lucmp%psoil(npnew,:)     = luc%psoil_y(p,q,:)
                lucmp%psoillab(npnew)    = luc%psoillab_y(p,q)
                lucmp%psoilsorb(npnew)   = luc%psoilsorb_y(p,q)
                lucmp%psoilocc(npnew)    = luc%psoilocc_y(p,q)
                lucmp%pwoodprod(npnew,:) = luc%pwoodprod_y(p,q,:)
             endif
          endif
          !update patch_type
       enddo   ! end of "q"
    enddo  ! end of "p"

    print *, 'npnew npold', npnew,npold
     
 END SUBROUTINE landuse_land2mpx

 SUBROUTINE landuse_checks(mlon,mlat,landmask,luc)
 ! check mass balance and write output CNP pool sizes for each PFT
 use landuse_constant,     ONLY: mvmax
 use landuse_variable,     ONLY: landuse_mland
 USE cable_def_types_mod,  ONLY: mland,r_2

 IMPLICIT NONE
 integer mlon,mlat
 real, parameter      :: xunit = 1.0e-15
 TYPE(landuse_mland)  :: luc
 integer,       dimension(mlon,mlat) :: landmask
 real(r_2), dimension(mvmax)     :: areapft                    
 real(r_2), dimension(mvmax)     :: cpland,npland,ppland
 real(r_2), dimension(mvmax)     :: clland,nlland,plland
 real(r_2), dimension(mvmax)     :: csland,nsland,psland
 real(r_2), dimension(mvmax)     :: clabland,nsminland,pslabland,pssorbland,psoccland
 real(r_2), dimension(mvmax)     :: cwoodland,nwoodland,pwoodland
 integer n,v
 real(r_2)  totalc,totaln,totalp,totarea

    areapft=0.0
    totalc=0.0;      totaln=0.0;      totalp=0.0;     totarea=0.0
    cpland = 0.0;    npland = 0.0;    ppland=0.0
    clland = 0.0;    nlland = 0.0;    plland = 0.0
    csland = 0.0;    nsland = 0.0;    psland = 0.0
    nsminland=0.0;   pslabland=0.0;   pssorbland=0.0; psoccland = 0.0
    cwoodland = 0.0; nwoodland = 0.0; pwoodland = 0.0

    do n=1,mland
    do v=1,mvmax
    if(luc%patchfrac_y(n,v) > 0.0) then

       areapft(v)= areapft(v)+ luc%patchfrac_y(n,v)
       cpland(v) = cpland(v) + luc%patchfrac_y(n,v) * sum(luc%cplant_y(n,v,1:3))
       npland(v) = npland(v) + luc%patchfrac_y(n,v) * sum(luc%nplant_y(n,v,1:3))
       ppland(v) = ppland(v) + luc%patchfrac_y(n,v) * sum(luc%pplant_y(n,v,1:3))
       clland(v) = clland(v) + luc%patchfrac_y(n,v) * sum(luc%clitter_y(n,v,1:3))
       nlland(v) = nlland(v) + luc%patchfrac_y(n,v) * sum(luc%nlitter_y(n,v,1:3))
       plland(v) = plland(v) + luc%patchfrac_y(n,v) * sum(luc%plitter_y(n,v,1:3))
       csland(v) = csland(v) + luc%patchfrac_y(n,v) * sum(luc%csoil_y(n,v,1:3))
       nsland(v) = nsland(v) + luc%patchfrac_y(n,v) * sum(luc%nsoil_y(n,v,1:3))
       psland(v) = psland(v) + luc%patchfrac_y(n,v) * sum(luc%psoil_y(n,v,1:3))
 
       clabland(v)   = clabland(v)   + luc%patchfrac_y(n,v) * luc%clabile_y(n,v)
       nsminland(v)  = nsminland(v)  + luc%patchfrac_y(n,v) * luc%nsoilmin_y(n,v)
       pslabland(v)  = pslabland(v)  + luc%patchfrac_y(n,v) * luc%psoillab_y(n,v)
       pssorbland(v) = pssorbland(v) + luc%patchfrac_y(n,v) * luc%psoilsorb_y(n,v)
       psoccland(v)  = psoccland(v)  + luc%patchfrac_y(n,v) * luc%psoilocc_y(n,v)
       cwoodland(v)  = cwoodland(v)  + luc%patchfrac_y(n,v) * sum(luc%cwoodprod_y(n,v,1:3))
       nwoodland(v)  = nwoodland(v)  + luc%patchfrac_y(n,v) * sum(luc%nwoodprod_y(n,v,1:3))
       pwoodland(v)  = pwoodland(v)  + luc%patchfrac_y(n,v) * sum(luc%pwoodprod_y(n,v,1:3))
    endif
    enddo
    enddo

    totalc = sum(cpland+clland+csland+clabland+cwoodland)
    totaln = sum(npland+nlland+nsland+nsminland+nwoodland)
    totalp = sum(ppland+plland+psland+pslabland+pssorbland+psoccland+pwoodland)

    do v=1,mvmax
       write(21,201) v, areapft(v),cpland(v)+clland(v)+csland(v)+clabland(v)+cwoodland(v),     &
                        npland(v)+nlland(v)+nsland(v)+nsminland(v)+nwoodland(v),               &
                        ppland(v)+plland(v)+psland(v)+pslabland(v)+pssorbland(v)+psoccland(v)+pwoodland(v)
    enddo
    write(21,202) sum(areapft),totalc,totaln,totalp

201 format(i3,2x,20(f10.5,2x))
202 format(20(f10.5,2x))
   
 END SUBROUTINE landuse_checks
