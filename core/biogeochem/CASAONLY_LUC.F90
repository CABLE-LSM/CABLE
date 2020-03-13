SUBROUTINE CASAONLY_LUC( dels,kstart,kend,veg,soil,casabiome,casapool, &
     casaflux,casamet,casabal,phen,POP,climate,LALLOC,LUC_EXPT, POPLUC, &
     sum_casapool, sum_casaflux, c13o2flux, c13o2pools, sum_c13o2pools, c13o2luc )


  USE cable_def_types_mod
  USE cable_carbon_module
  USE cable_common_module, ONLY: CABLE_USER, is_casa_time
  USE cable_IO_vars_module, ONLY: logn, landpt
  USE casadimension
  USE casaparm
  USE casavariable
  USE phenvariable
  USE POP_Types,  Only: POP_TYPE
  USE POPMODULE,            ONLY: POPStep, POP_init_single
  USE TypeDef,              ONLY: dp
  USE CABLE_LUC_EXPT, ONLY: LUC_EXPT_TYPE, read_LUH2, &
       ptos, ptog, stog, gtos, pharv, smharv, syharv, &
       ptoc, ptoq, stoc, stoq, ctos, qtos
  USE POPLUC_Types
  USE POPLUC_Module, ONLY: POPLUCStep, POPLUC_weights_Transfer, WRITE_LUC_OUTPUT_NC, &
       POP_LUC_CASA_transfer,  WRITE_LUC_RESTART_NC, READ_LUC_RESTART_NC, &
       POPLUC_set_patchfrac, WRITE_LUC_OUTPUT_GRID_NC
  ! 13C
  use cable_c13o2_def, only: c13o2_flux, c13o2_pool, c13o2_luc, c13o2_update_sum_pools, c13o2_zero_sum_pools
  use cable_c13o2,     only: c13o2_save_casapool, c13o2_update_pools, c13o2_save_luc, c13o2_update_luc, &
       c13o2_create_output, c13o2_write_output, c13o2_close_output, &
       c13o2_nvars_output, &
       c13o2_print_delta_pools, c13o2_print_delta_luc

  IMPLICIT NONE
  
  !!CLN  CHARACTER(LEN=99), INTENT(IN)  :: fcnpspin
  REAL,    INTENT(IN)    :: dels
  INTEGER, INTENT(IN)    :: kstart
  INTEGER, INTENT(IN)    :: kend
  INTEGER, INTENT(IN)    :: LALLOC
  TYPE (veg_parameter_type),  INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (soil_parameter_type), INTENT(INOUT) :: soil ! soil parameters
  TYPE (casa_biome),          INTENT(INOUT) :: casabiome
  TYPE (casa_pool),           INTENT(INOUT) :: casapool
  TYPE (casa_flux),           INTENT(INOUT) :: casaflux
  TYPE (casa_met),            INTENT(INOUT) :: casamet
  TYPE (casa_balance),        INTENT(INOUT) :: casabal
  TYPE (phen_variable),       INTENT(INOUT) :: phen
  TYPE (POP_TYPE),            INTENT(INOUT) :: POP
  TYPE (climate_TYPE),        INTENT(INOUT) :: climate
  TYPE (LUC_EXPT_TYPE),       INTENT(INOUT) :: LUC_EXPT
  TYPE(POPLUC_TYPE),          INTENT(INOUT) :: POPLUC
  TYPE (casa_pool),           INTENT(INOUT) :: sum_casapool
  TYPE (casa_flux),           INTENT(INOUT) :: sum_casaflux
  ! 13C
  type(c13o2_flux),           intent(inout) :: c13o2flux
  type(c13o2_pool),           intent(inout) :: c13o2pools
  type(c13o2_pool),           intent(inout) :: sum_c13o2pools
  type(c13o2_luc),            intent(inout) :: c13o2luc

  ! local variables
  INTEGER                  :: myearspin,nyear, yyyy, nyear_dump
  CHARACTER(LEN=99)        :: ncfile
  CHARACTER(LEN=4)         :: cyear
  INTEGER                  :: ktau,ktauday,nday,idoy
  real(r_2), dimension(mp)      :: cleaf2met, cleaf2str, croot2met, croot2str, cwood2cwd
  real(r_2), dimension(mp)      :: nleaf2met, nleaf2str, nroot2met, nroot2str, nwood2cwd
  real(r_2), dimension(mp)      :: pleaf2met, pleaf2str, proot2met, proot2str, pwood2cwd
  real(r_2), dimension(mp)      :: xnplimit,  xkNlimiting, xklitter, xksoil,xkleaf, xkleafcold, xkleafdry

  ! more variables to store the spinup pool size over the last 10 loops. Added by Yp Wang 30 Nov 2012
  integer :: ctime, k, j, l
  INTEGER, allocatable :: Iw(:) ! array of indices corresponding to woody (shrub or forest) tiles
  INTEGER :: count_sum_casa ! number of time steps over which casa pools &
  !and fluxes are aggregated (for output)

  ! 13C
  real(dp), dimension(c13o2pools%ntile,c13o2pools%npools) :: casasave
  real(dp), dimension(c13o2luc%nland,c13o2luc%npools)     :: lucsave
  integer :: c13o2_file_id
  character(len=40), dimension(c13o2_nvars_output) :: c13o2_vars
  integer,           dimension(c13o2_nvars_output) :: c13o2_var_ids

  if (.NOT.Allocated(Iw)) allocate(Iw(POP%np))

  !! vh_js !!
  IF (cable_user%CALL_POP) THEN
     Iw = POP%Iwood
  ENDIF

  CALL zero_sum_casa(sum_casapool, sum_casaflux)
  ! 13C
  if (cable_user%c13o2) call c13o2_zero_sum_pools(sum_c13o2pools)

  ktauday = int(24.0*3600.0/dels)
  nday    = (kend-kstart+1)/ktauday
  ctime   = 0
  count_sum_casa = 0
  myearspin = cable_user%yearend - cable_user%yearstart + 1
  yyyy      = cable_user%yearstart - 1
  do nyear=1, myearspin
     yyyy = yyyy+1
     write(*,*) 'casaonly_LUC', yyyy

     nyear_dump = mod(nyear, cable_user%casa_spin_endyear - cable_user%casa_spin_startyear + 1)
     if (nyear_dump == 0) &
          nyear_dump = cable_user%casa_spin_endyear - cable_user%casa_spin_startyear + 1
     write(cyear,fmt="(i4)") cable_user%casa_spin_startyear + nyear_dump - 1
     ncfile = trim(casafile%c2cdumppath)//'c2c_'//cyear//'_dump.nc'
     call read_casa_dump( ncfile, casamet, casaflux, phen, climate, c13o2flux, 1, 1, .true. )
    
     !!CLN901  format(A99)
     do idoy=1, mdyear
        ktau = (idoy-1)*ktauday + ktauday

        casamet%tairk(:)       = casamet%Tairkspin(:,idoy)
        casamet%tsoil(:,1)     = casamet%Tsoilspin_1(:,idoy)
        casamet%tsoil(:,2)     = casamet%Tsoilspin_2(:,idoy)
        casamet%tsoil(:,3)     = casamet%Tsoilspin_3(:,idoy)
        casamet%tsoil(:,4)     = casamet%Tsoilspin_4(:,idoy)
        casamet%tsoil(:,5)     = casamet%Tsoilspin_5(:,idoy)
        casamet%tsoil(:,6)     = casamet%Tsoilspin_6(:,idoy)
        casamet%moist(:,1)     = casamet%moistspin_1(:,idoy)
        casamet%moist(:,2)     = casamet%moistspin_2(:,idoy)
        casamet%moist(:,3)     = casamet%moistspin_3(:,idoy)
        casamet%moist(:,4)     = casamet%moistspin_4(:,idoy)
        casamet%moist(:,5)     = casamet%moistspin_5(:,idoy)
        casamet%moist(:,6)     = casamet%moistspin_6(:,idoy)
        casaflux%cgpp(:)       = casamet%cgppspin(:,idoy)
        casaflux%crmplant(:,1) = casamet%crmplantspin_1(:,idoy)
        casaflux%crmplant(:,2) = casamet%crmplantspin_2(:,idoy)
        casaflux%crmplant(:,3) = casamet%crmplantspin_3(:,idoy)
        phen%phase(:)      = phen%phasespin(:,idoy)
        phen%doyphase(:,1) = phen%doyphasespin_1(:,idoy)
        phen%doyphase(:,2) = phen%doyphasespin_2(:,idoy)
        phen%doyphase(:,3) = phen%doyphasespin_3(:,idoy)
        phen%doyphase(:,4) = phen%doyphasespin_4(:,idoy)
        climate%qtemp_max_last_year(:) = real(casamet%mtempspin(:,idoy))
        ! 13C
        if (cable_user%c13o2) then
           c13o2flux%cAn12(:) = casamet%cAn12spin(:,idoy)
           c13o2flux%cAn(:)   = casamet%cAn13spin(:,idoy)
        endif
        
        ! 13C
        if (cable_user%c13o2) then
           print*, '13C in casaonly_luc - 00 ', nyear, idoy
           call c13o2_print_delta_pools(casapool, casaflux, c13o2pools)
           call c13o2_print_delta_luc(popluc, c13o2luc)
        endif
        if (cable_user%c13o2) call c13o2_save_casapool(casapool, casasave)
        CALL biogeochem(ktau,dels,idoy,LALLOC,veg,soil,casabiome,casapool,casaflux, &
             casamet,casabal,phen,POP,climate,xnplimit,xkNlimiting,xklitter, &
             xksoil,xkleaf,xkleafcold,xkleafdry,&
             cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,         &
             nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,         &
             pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd)
        ! 13C
        if (cable_user%c13o2) call c13o2_update_pools(casasave, casaflux, c13o2flux, c13o2pools)
        if (cable_user%c13o2) then
           print*, '13C in casaonly_luc - 01 ', nyear, idoy
           call c13o2_print_delta_pools(casapool, casaflux, c13o2pools)
           call c13o2_print_delta_luc(popluc, c13o2luc)
        endif
 
        ! update time-aggregates of casa pools and fluxes
        CALL update_sum_casa(sum_casapool, sum_casaflux, casapool, casaflux, &
             & .TRUE. , .FALSE., 1)
        ! 13C
        if (cable_user%c13o2) &
             call c13o2_update_sum_pools(sum_c13o2pools, c13o2pools, .true., .false., 1)
        count_sum_casa = count_sum_casa + 1

        ! accumulate annual variables for use in POP
        IF (idoy==1) THEN
           ! (assumes 70% of wood NPP is allocated above ground)
           casaflux%stemnpp  = casaflux%cnpp * casaflux%fracCalloc(:,2) * 0.7_dp
           casabal%LAImax    = casamet%glai
           casabal%Cleafmean = casapool%cplant(:,1) / real(mdyear,dp) / 1000._dp
           casabal%Crootmean = casapool%cplant(:,3) / real(mdyear,dp) / 1000._dp
        ELSE
           casaflux%stemnpp  = casaflux%stemnpp + casaflux%cnpp * casaflux%fracCalloc(:,2) * 0.7_dp
           casabal%LAImax    = max(casamet%glai, casabal%LAImax)
           casabal%Cleafmean = casabal%Cleafmean + casapool%cplant(:,1) / real(mdyear,dp) / 1000._dp
           casabal%Crootmean = casabal%Crootmean + casapool%cplant(:,3) / real(mdyear,dp) / 1000._dp
        ENDIF
      IF (CABLE_USER%POPLUC) THEN
        IF(idoy==mdyear) THEN ! end of year
           LUC_EXPT%CTSTEP = yyyy -  LUC_EXPT%FirstYear + 1

           CALL READ_LUH2(LUC_EXPT)

           DO k=1,mland
              POPLUC%ptos(k)   = real(LUC_EXPT%INPUT(ptos)%VAL(k), dp)
              POPLUC%ptog(k)   = real(LUC_EXPT%INPUT(ptog)%VAL(k), dp)
              POPLUC%stog(k)   = real(LUC_EXPT%INPUT(stog)%VAL(k), dp)
              POPLUC%gtop(k)   = 0.0_dp
              POPLUC%gtos(k)   = real(LUC_EXPT%INPUT(gtos)%VAL(k), dp)
              POPLUC%pharv(k)  = real(LUC_EXPT%INPUT(pharv)%VAL(k), dp)
              POPLUC%smharv(k) = real(LUC_EXPT%INPUT(smharv)%VAL(k), dp)
              POPLUC%syharv(k) = real(LUC_EXPT%INPUT(syharv)%VAL(k), dp)

              POPLUC%ptoc(k) = real(LUC_EXPT%INPUT(ptoc)%VAL(k), dp)
              POPLUC%ptoq(k) = real(LUC_EXPT%INPUT(ptoq)%VAL(k), dp)
              POPLUC%stoc(k) = real(LUC_EXPT%INPUT(stoc)%VAL(k), dp)
              POPLUC%stoq(k) = real(LUC_EXPT%INPUT(stoq)%VAL(k), dp)
              POPLUC%ctos(k) = real(LUC_EXPT%INPUT(ctos)%VAL(k), dp)
              POPLUC%qtos(k) = real(LUC_EXPT%INPUT(qtos)%VAL(k), dp)
             
              POPLUC%thisyear = yyyy
           ENDDO
           !stop
           ! set landuse index for secondary forest POP landscapes
           DO k=1,POP%np
              IF (yyyy.eq.LUC_EXPT%YearStart) THEN
                 if (veg%iLU(POP%Iwood(k)).eq.2) then
                    POP%pop_grid(k)%LU = 2
                 endif
              endif
           ENDDO

           ! zero secondary forest tiles in POP where secondary forest area is zero
           DO k=1,mland
              if ((POPLUC%primf(k)-POPLUC%frac_forest(k))==0.0_dp &
                   .and. (.not.LUC_EXPT%prim_only(k))) then

                 j = landpt(k)%cstart+1
                 do l=1,size(POP%Iwood)
                    if ( POP%Iwood(l) == j) then
                       CALL POP_init_single(POP,veg%disturbance_interval,l)
                       exit
                    endif
                 enddo

                 casapool%cplant(j,leaf)  = 0.01_dp
                 casapool%nplant(j,leaf)  = casabiome%ratioNCplantmin(veg%iveg(j),leaf) * casapool%cplant(j,leaf)
                 casapool%pplant(j,leaf)  = casabiome%ratioPCplantmin(veg%iveg(j),leaf) * casapool%cplant(j,leaf)

                 casapool%cplant(j,froot) = 0.01_dp
                 casapool%nplant(j,froot) = casabiome%ratioNCplantmin(veg%iveg(j),froot) * casapool%cplant(j,froot)
                 casapool%pplant(j,froot) = casabiome%ratioPCplantmin(veg%iveg(j),froot) * casapool%cplant(j,froot)

                 casapool%cplant(j,wood)  = 0.01_dp
                 casapool%nplant(j,wood)  = casabiome%ratioNCplantmin(veg%iveg(j),wood) * casapool%cplant(j,wood)
                 casapool%pplant(j,wood)  = casabiome%ratioPCplantmin(veg%iveg(j),wood) * casapool%cplant(j,wood)
                 casaflux%frac_sapwood(j) = 1.0_dp

                 ! 13C
                 if (cable_user%c13o2) then
                    c13o2pools%cplant(j,leaf)  = 0.01_dp ! * vpdbc13 / vpdbc13 ! Divide by 13C
                    c13o2pools%cplant(j,wood)  = 0.01_dp ! * vpdbc13 / vpdbc13 ! so that about same numerical precision as 12C
                    c13o2pools%cplant(j,froot) = 0.01_dp ! * vpdbc13 / vpdbc13 !
                 endif
                 if (cable_user%c13o2) then
                    print*, '13C in casaonly_luc - 11'
                    call c13o2_print_delta_pools(casapool, casaflux, c13o2pools)
                    call c13o2_print_delta_luc(popluc, c13o2luc)
                 endif
              endif
           ENDDO

           CALL POPLUCStep(POPLUC,yyyy)

           CALL POPLUC_weights_transfer(POPLUC,POP,LUC_EXPT)

           CALL POPdriver(casaflux,casabal,veg, POP)

           CALL POP_IO( pop, casamet, YYYY, 'WRITE_EPI', &
                ( YYYY.EQ.cable_user%YearEnd ) )

           ! 13C
           if (cable_user%c13o2) then
              print*, '13C in casaonly_luc - 03'
              call c13o2_print_delta_pools(casapool, casaflux, c13o2pools)
              call c13o2_print_delta_luc(popluc, c13o2luc)
           endif
           if (cable_user%c13o2) call c13o2_save_luc(casapool, popluc, casasave, lucsave)
           CALL POP_LUC_CASA_transfer(POPLUC,POP,LUC_EXPT,casapool,casabal,casaflux,ktauday)
           ! 13C
#ifdef __C13DEBUG__
           if (cable_user%c13o2) &
                call c13o2_update_luc(casasave, lucsave, popluc, luc_expt%prim_only, c13o2pools, c13o2luc, casapool)
#else
           if (cable_user%c13o2) call c13o2_update_luc(casasave, lucsave, popluc, luc_expt%prim_only, c13o2pools, c13o2luc)
#endif
           if (cable_user%c13o2) then
              print*, '13C in casaonly_luc - 04'
              call c13o2_print_delta_pools(casapool, casaflux, c13o2pools)
              call c13o2_print_delta_luc(popluc, c13o2luc)
           endif
           CALL WRITE_LUC_OUTPUT_NC(POPLUC, YYYY, (YYYY.EQ.cable_user%YearEnd))
           CALL POPLUC_set_patchfrac(POPLUC, LUC_EXPT)

        ENDIF  ! end of year
     ELSE
        IF(idoy==mdyear) THEN ! end of year          
           CALL POPdriver(casaflux, casabal, veg, POP)
        endif
       ! CALL POP_IO( pop, casamet, YYYY, 'WRITE_EPI', &
       !         ( YYYY.EQ.cable_user%YearEnd ) )
     ENDIF ! IF (CABLE_USER%POPLUC) 

        IF ( IS_CASA_TIME("write", yyyy, ktau, kstart, 0, kend, ktauday, logn) ) THEN
           ctime = ctime + 1

           CALL update_sum_casa(sum_casapool, sum_casaflux, casapool, casaflux, &
                .FALSE. , .TRUE. , count_sum_casa)
           ! 13C
           if (cable_user%c13o2) then
              call c13o2_update_sum_pools(sum_c13o2pools, c13o2pools, &
                   .false., .true., count_sum_casa)
           endif

           CALL WRITE_CASA_OUTPUT_NC( veg, casamet, sum_casapool, casabal, sum_casaflux, &
                .true., ctime, ( idoy.eq.mdyear .AND. YYYY .EQ. cable_user%YearEnd ) )
           ! 13C
           if (cable_user%c13o2) then
              if (ctime == 1) then
                 call c13o2_create_output(casamet, sum_c13o2pools, c13o2_file_id, c13o2_vars, c13o2_var_ids)
              endif
              call c13o2_write_output(c13o2_file_id, c13o2_vars, c13o2_var_ids, ctime, sum_c13o2pools)
              if ( (idoy == mdyear) .and. (YYYY == cable_user%YearEnd) ) &
                   call c13o2_close_output(c13o2_file_id)
           end if
           count_sum_casa = 0
           CALL zero_sum_casa(sum_casapool, sum_casaflux)
           if (cable_user%c13o2) call c13o2_zero_sum_pools(sum_c13o2pools)

        ENDIF
     enddo
     
  enddo

  IF (CABLE_USER%POPLUC) CALL WRITE_LUC_RESTART_NC( POPLUC, YYYY )

END SUBROUTINE CASAONLY_LUC
