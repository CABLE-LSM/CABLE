SUBROUTINE spincasacnp(dels,kstart,kend,mloop,veg,soil,casabiome,casapool, &
     casaflux,casamet,casabal,phen,POP,climate,LALLOC, c13o2flux, c13o2pools, &
     BLAZE, SIMFIRE)

  USE cable_def_types_mod
  USE cable_carbon_module
  USE cable_common_module, ONLY: CABLE_USER
  USE casadimension
  USE casaparm ! e.g. leaf
  USE casavariable
  USE phenvariable
  USE POP_Types,           Only: POP_TYPE
  USE POP_Constants, ONLY: shootfrac
  use TypeDef,             only: dp
  ! 13C
  use cable_c13o2_def,     only: c13o2_pool, c13o2_flux
  use cable_c13o2,         only: c13o2_save_casapool, c13o2_update_pools, &
       c13o2_create_output, c13o2_write_output, c13o2_close_output, &
       c13o2_nvars_output
  use mo_isotope,          only: isoratio
  USE BLAZE_MOD,     ONLY: TYPE_BLAZE,  BLAZE_ACCOUNTING,  WRITE_BLAZE_OUTPUT_NC
  USE SIMFIRE_MOD,   ONLY: TYPE_SIMFIRE
  implicit none

  !!CLN  character(len=99), intent(in)  :: fcnpspin
  real,                      intent(in)    :: dels
  integer,                   intent(in)    :: kstart
  integer,                   intent(in)    :: kend
  integer,                   intent(in)    :: mloop
  integer,                   intent(in)    :: lalloc
  type(veg_parameter_type),  intent(inout) :: veg       ! vegetation parameters
  type(soil_parameter_type), intent(inout) :: soil      ! soil parameters
  type(casa_biome),          intent(inout) :: casabiome
  type(casa_pool),           intent(inout) :: casapool
  type(casa_flux),           intent(inout) :: casaflux
  type(casa_met),            intent(inout) :: casamet
  type(casa_balance),        intent(inout) :: casabal
  type(phen_variable),       intent(inout) :: phen
  type(pop_type),            intent(inout) :: pop
  type(climate_type),        intent(inout) :: climate
  ! 13C
  type(c13o2_flux),          intent(inout) :: c13o2flux
  type(c13o2_pool),          intent(inout) :: c13o2pools
  ! BLAZE
  type(TYPE_BLAZE),          intent(inout) :: BLAZE
  type(TYPE_SIMFIRE),        intent(inout) :: SIMFIRE

  ! type(casa_met) :: casaspin

  ! local variables
  real(r_2), dimension(:), allocatable, save  :: avg_cleaf2met, avg_cleaf2str, avg_croot2met, avg_croot2str, avg_cwood2cwd
  real(r_2), dimension(:), allocatable, save  :: avg_nleaf2met, avg_nleaf2str, avg_nroot2met, avg_nroot2str, avg_nwood2cwd
  real(r_2), dimension(:), allocatable, save  :: avg_pleaf2met, avg_pleaf2str, avg_proot2met, avg_proot2str, avg_pwood2cwd
  real(r_2), dimension(:), allocatable, save  :: avg_cgpp,      avg_cnpp,      avg_nuptake,   avg_puptake
  real(r_2), dimension(:), allocatable, save  :: avg_nsoilmin,  avg_psoillab,  avg_psoilsorb, avg_psoilocc
  !chris 12/oct/2012 for spin up casa
  real(r_2), dimension(:), allocatable, save  :: avg_ratioNCsoilmic,  avg_ratioNCsoilslow,  avg_ratioNCsoilpass
  real(r_2), dimension(:), allocatable, save  :: avg_xnplimit,  avg_xkNlimiting,avg_xklitter, avg_xksoil

  ! local variables
  integer                  :: myearspin, nyear, nloop1
  character(len=99)        :: ncfile
  character(len=4)         :: cyear
  integer                  :: ktau, ktauday, nday, idoy, ktauy, nloop, LOY, YYYY
  real(r_2), dimension(mp) :: cleaf2met, cleaf2str, croot2met, croot2str, cwood2cwd
  real(r_2), dimension(mp) :: nleaf2met, nleaf2str, nroot2met, nroot2str, nwood2cwd
  real(r_2), dimension(mp) :: pleaf2met, pleaf2str, proot2met, proot2str, pwood2cwd
  real(r_2), dimension(mp) :: xnplimit,  xkNlimiting, xklitter, xksoil, xkleaf, xkleafcold, xkleafdry

  ! more variables to store the spinup pool size over the last 10 loops. Added by Yp Wang 30 Nov 2012
  real(dp), allocatable, save :: LAImax(:), Cleafmean(:), Crootmean(:)
  real(dp), allocatable :: NPPtoGPP(:)
  integer,  allocatable :: Iw(:) ! array of indices corresponding to woody (shrub or forest) tiles
  integer(i_d) :: ctime
  real(dp) :: rday

  ! 13C
  real(dp), dimension(c13o2pools%ntile,c13o2pools%npools) :: casasave
  integer :: c13o2_file_id
  character(len=40), dimension(c13o2_nvars_output) :: c13o2_vars
  integer,           dimension(c13o2_nvars_output) :: c13o2_var_ids
  real(r_2), dimension(:), allocatable :: avg_c13leaf2met, avg_c13leaf2str, avg_c13root2met, &
       avg_c13root2str, avg_c13wood2cwd

   if (.not. allocated(LAIMax))    allocate(LAIMax(mp))
   if (.not. allocated(Cleafmean)) allocate(Cleafmean(mp))
   if (.not. allocated(Crootmean)) allocate(Crootmean(mp))
   if (.not. allocated(NPPtoGPP))  allocate(NPPtoGPP(mp))
   if (.not. allocated(Iw))        allocate(Iw(POP%np))

   ctime = 0
   LOY = 365
   !! vh_js !!
   if (cable_user%CALL_POP) Iw = POP%Iwood

   ktauday = int(24.0*3600.0/dels)
   nday    = (kend-kstart+1)/ktauday

   !chris 12/oct/2012 for spin up casa
   if (.not.(allocated(avg_cleaf2met))) &
        allocate(avg_cleaf2met(mp), avg_cleaf2str(mp), avg_croot2met(mp), avg_croot2str(mp), avg_cwood2cwd(mp), &
        avg_nleaf2met(mp), avg_nleaf2str(mp), avg_nroot2met(mp), avg_nroot2str(mp), avg_nwood2cwd(mp), &
        avg_pleaf2met(mp), avg_pleaf2str(mp), avg_proot2met(mp), avg_proot2str(mp), avg_pwood2cwd(mp), &
        avg_cgpp(mp), avg_cnpp(mp), avg_nuptake(mp), avg_puptake(mp), &
        avg_xnplimit(mp), avg_xkNlimiting(mp), avg_xklitter(mp), avg_xksoil(mp), &
        avg_rationcsoilmic(mp), avg_rationcsoilslow(mp), avg_rationcsoilpass(mp), &
        avg_nsoilmin(mp), avg_psoillab(mp), avg_psoilsorb(mp), avg_psoilocc(mp))
   ! 13C - allocate in any case even if cable_user%c13o2==.false. to pass to analytic soil and litter pools
   allocate(avg_c13leaf2met(mp))
   allocate(avg_c13leaf2str(mp))
   allocate(avg_c13root2met(mp))
   allocate(avg_c13root2str(mp))
   allocate(avg_c13wood2cwd(mp))

   !!CLN  OPEN(91, file=fcnpspin)
   !!CLN  read(91,*) myearspin
   myearspin = CABLE_USER%CASA_SPIN_ENDYEAR - CABLE_USER%CASA_SPIN_STARTYEAR + 1
   
   ! compute the mean fluxes and residence time of each carbon pool
   avg_cleaf2met       = 0.0_dp
   avg_cleaf2str       = 0.0_dp
   avg_croot2met       = 0.0_dp
   avg_croot2str       = 0.0_dp
   avg_cwood2cwd       = 0.0_dp
   avg_nleaf2met       = 0.0_dp
   avg_nleaf2str       = 0.0_dp
   avg_nroot2met       = 0.0_dp
   avg_nroot2str       = 0.0_dp
   avg_nwood2cwd       = 0.0_dp
   avg_pleaf2met       = 0.0_dp
   avg_pleaf2str       = 0.0_dp
   avg_proot2met       = 0.0_dp
   avg_proot2str       = 0.0_dp
   avg_pwood2cwd       = 0.0_dp
   avg_cgpp            = 0.0_dp
   avg_cnpp            = 0.0_dp
   avg_nuptake         = 0.0_dp
   avg_puptake         = 0.0_dp
   avg_xnplimit        = 0.0_dp
   avg_xkNlimiting     = 0.0_dp
   avg_xklitter        = 0.0_dp
   avg_xksoil          = 0.0_dp
   avg_nsoilmin        = 0.0_dp
   avg_psoillab        = 0.0_dp
   avg_psoilsorb       = 0.0_dp
   avg_psoilocc        = 0.0_dp
   avg_rationcsoilmic  = 0.0_dp
   avg_rationcsoilslow = 0.0_dp
   avg_rationcsoilpass = 0.0_dp
   ! 13C
   if (cable_user%c13o2) then
      avg_c13leaf2met = 0.0_dp
      avg_c13leaf2str = 0.0_dp
      avg_c13root2met = 0.0_dp
      avg_c13root2str = 0.0_dp
      avg_c13wood2cwd = 0.0_dp
   endif

!   write(600,*) 'csoil3 init: ', casapool%csoil(3,:)
!   write(600,*) 'csoil1 init: ', casapool%csoil(1,:)

   do nyear=1, myearspin
      write(cyear,FMT="(I4)") CABLE_USER%CASA_SPIN_STARTYEAR + nyear - 1
      YYYY = CABLE_USER%CASA_SPIN_STARTYEAR + nyear - 1
     ncfile = trim(casafile%c2cdumppath)//'c2c_'//cyear//'_dump.nc'
     call read_casa_dump(ncfile, casamet, casaflux, phen, climate, c13o2flux, ktau, kend, .true.)

     do idoy=1,mdyear
        ktau = (idoy-1)*ktauday + 1
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
        phen%phase(:)          = phen%phasespin(:,idoy)
        phen%doyphase(:,1)     = phen%doyphasespin_1(:,idoy)
        phen%doyphase(:,2)     = phen%doyphasespin_2(:,idoy)
        phen%doyphase(:,3)     = phen%doyphasespin_3(:,idoy)
        phen%doyphase(:,4)     = phen%doyphasespin_4(:,idoy)
        climate%qtemp_max_last_year(:) = real(casamet%mtempspin(:,idoy))
        climate%frec(:)        = real(casamet%frecspin(:,idoy))
        ! casaflux%Nmindep and casaflux%Pdep set in read_casa_dump
        ! 13C
        if (cable_user%c13o2) then
           c13o2flux%cAn12(:) = casamet%cAn12spin(:,idoy)
           c13o2flux%cAn(:)   = casamet%cAn13spin(:,idoy)
        endif

        ! 13C
        if (cable_user%c13o2) call c13o2_save_casapool(casapool, casasave)
        call biogeochem(ktau,dels,idoy,LALLOC,veg,soil,casabiome,casapool,casaflux, &
             casamet,casabal,phen,POP,climate,xnplimit,xkNlimiting,xklitter, &
             xksoil,xkleaf,xkleafcold,xkleafdry, &
             cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd, &
             nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd, &
             pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd)
        ! 13C
        if (cable_user%c13o2) then
           avg_c13leaf2met(:) = avg_c13leaf2met(:) + &
                cleaf2met(:) * isoratio(c13o2pools%cplant(:,leaf), casasave(:,leaf), 0.0_dp, tiny(1.0_dp)) ! 1.0_dp
           avg_c13leaf2str(:) = avg_c13leaf2str(:) + &
                cleaf2str(:) * isoratio(c13o2pools%cplant(:,leaf), casasave(:,leaf), 0.0_dp, tiny(1.0_dp))
           avg_c13root2met(:) = avg_c13root2met(:) + &
                croot2met(:) * isoratio(c13o2pools%cplant(:,froot), casasave(:,froot), 0.0_dp, tiny(1.0_dp))
           avg_c13root2str(:) = avg_c13root2str(:) + &
                croot2str(:) * isoratio(c13o2pools%cplant(:,froot), casasave(:,froot), 0.0_dp, tiny(1.0_dp))
           avg_c13wood2cwd(:) = avg_c13wood2cwd(:) + &
                cwood2cwd(:) * isoratio(c13o2pools%cplant(:,wood), casasave(:,wood), 0.0_dp, tiny(1.0_dp))
           call c13o2_update_pools(casasave, casaflux, c13o2flux, c13o2pools)
        endif

        if (cable_user%CALL_BLAZE) then
           CALL BLAZE_ACCOUNTING(BLAZE, climate, ktau, dels, YYYY , idoy)

           call blaze_driver(blaze%ncells, blaze, simfire, casapool, casaflux, &
                casamet, climate, real(shootfrac), idoy, YYYY, 1, POP, veg)

        endif
           

        if (cable_user%CALL_POP .and. (POP%np .gt. 0)) then ! CALL_POP

           if (cable_user%CALL_POP) then ! accumulate input variables for POP
              ! accumulate annual variables for use in POP
              if (mod(ktau/ktauday,LOY) == 1) then
                 ! (assumes 70% of wood NPP is allocated above ground)
                 casaflux%stemnpp  = casaflux%cnpp * casaflux%fracCalloc(:,2) * 0.7_dp
                 casabal%LAImax    = casamet%glai
                 casabal%Cleafmean = casapool%cplant(:,1) / real(LOY,dp) / 1000._dp
                 casabal%Crootmean = casapool%cplant(:,3) / real(LOY,dp) / 1000._dp
              else
                 casaflux%stemnpp  = casaflux%stemnpp + casaflux%cnpp * casaflux%fracCalloc(:,2) * 0.7_dp
                 casabal%LAImax    = max(casamet%glai, casabal%LAImax)
                 casabal%Cleafmean = casabal%Cleafmean + casapool%cplant(:,1) / real(LOY,dp) / 1000.0_dp
                 casabal%Crootmean = casabal%Crootmean + casapool%cplant(:,3) / real(LOY,dp) / 1000.0_dp
              endif
           else
              casaflux%stemnpp = 0.0_dp
           endif ! CALL_POP

           !CALL WRITE_CASA_OUTPUT_NC (veg, casamet, casapool, casabal, casaflux, &
           !            .true., ctime, .FALSE.  )
           !            ctime = ctime+1
           if (idoy==mdyear) then ! end of year
              call POPdriver(casaflux, casabal, veg, POP)
              !MC - ToDo - update 13CO2 harvest
              ! CALL POP_IO( pop, casamet, nyear, 'WRITE_EPI', &
              !     (.FALSE.))
              ! CALL WRITE_CASA_OUTPUT_NC (veg, casamet, casapool, casabal, casaflux, &
              !     .true., ctime, .FALSE.  )
              ! ctime = ctime+1
           endif  ! end of year
            
        else ! (cable_user%CALL_POP .and. (POP%np.gt.0))

           ! if (idoy==mdyear) then ! end of year
           !    ! call write_casa_output_nc(veg, casamet, casapool, casabal, casaflux, .true., ctime, .false.)
           !    ctime = ctime + 1
           ! endif  ! end of year
           casaflux%stemnpp = 0.0_dp
           
        endif ! (cable_user%CALL_POP .and. (POP%np.gt.0))

        ! WHERE(xkNlimiting .eq. 0)  !Chris Lu 4/June/2012
        !    xkNlimiting = 0.001
        ! END WHERE
        ! nptx = 8173

        avg_cleaf2met       = avg_cleaf2met       + cleaf2met
        avg_cleaf2str       = avg_cleaf2str       + cleaf2str
        avg_croot2met       = avg_croot2met       + croot2met
        avg_croot2str       = avg_croot2str       + croot2str
        avg_cwood2cwd       = avg_cwood2cwd       + cwood2cwd

        avg_nleaf2met       = avg_nleaf2met       + nleaf2met
        avg_nleaf2str       = avg_nleaf2str       + nleaf2str
        avg_nroot2met       = avg_nroot2met       + nroot2met
        avg_nroot2str       = avg_nroot2str       + nroot2str
        avg_nwood2cwd       = avg_nwood2cwd       + nwood2cwd

        avg_pleaf2met       = avg_pleaf2met       + pleaf2met
        avg_pleaf2str       = avg_pleaf2str       + pleaf2str
        avg_proot2met       = avg_proot2met       + proot2met
        avg_proot2str       = avg_proot2str       + proot2str
        avg_pwood2cwd       = avg_pwood2cwd       + pwood2cwd

        avg_cgpp            = avg_cgpp            + casaflux%cgpp
        avg_cnpp            = avg_cnpp            + casaflux%cnpp
        avg_nuptake         = avg_nuptake         + casaflux%Nminuptake
        avg_puptake         = avg_puptake         + casaflux%Plabuptake

        avg_xnplimit        = avg_xnplimit        + xnplimit
        avg_xkNlimiting     = avg_xkNlimiting     + xkNlimiting
        avg_xklitter        = avg_xklitter        + xklitter
        avg_xksoil          = avg_xksoil          + xksoil

        avg_nsoilmin        = avg_nsoilmin        + casapool%nsoilmin
        avg_psoillab        = avg_psoillab        + casapool%psoillab
        avg_psoilsorb       = avg_psoilsorb       + casapool%psoilsorb
        avg_psoilocc        = avg_psoilocc        + casapool%psoilocc

        avg_rationcsoilmic  = avg_rationcsoilmic  + casapool%ratioNCsoilnew(:,mic)
        avg_rationcsoilslow = avg_rationcsoilslow + casapool%ratioNCsoilnew(:,slow)
        avg_rationcsoilpass = avg_rationcsoilpass + casapool%ratioNCsoilnew(:,pass)
     enddo ! idoy
  enddo ! nyear

  !!CLN    CLOSE(91)
  ! average
  rday = 1.0_dp / real(nday*myearspin, dp)
  
  avg_cleaf2met       = avg_cleaf2met       * rday
  avg_cleaf2str       = avg_cleaf2str       * rday
  avg_croot2met       = avg_croot2met       * rday
  avg_croot2str       = avg_croot2str       * rday
  avg_cwood2cwd       = avg_cwood2cwd       * rday

  avg_nleaf2met       = avg_nleaf2met       * rday
  avg_nleaf2str       = avg_nleaf2str       * rday
  avg_nroot2met       = avg_nroot2met       * rday
  avg_nroot2str       = avg_nroot2str       * rday
  avg_nwood2cwd       = avg_nwood2cwd       * rday

  avg_pleaf2met       = avg_pleaf2met       * rday
  avg_pleaf2str       = avg_pleaf2str       * rday
  avg_proot2met       = avg_proot2met       * rday
  avg_proot2str       = avg_proot2str       * rday
  avg_pwood2cwd       = avg_pwood2cwd       * rday

  avg_cgpp            = avg_cgpp            * rday
  avg_cnpp            = avg_cnpp            * rday
  avg_nuptake         = avg_nuptake         * rday
  avg_puptake         = avg_puptake         * rday

  avg_xnplimit        = avg_xnplimit        * rday
  avg_xkNlimiting     = avg_xkNlimiting     * rday
  avg_xklitter        = avg_xklitter        * rday
  avg_xksoil          = avg_xksoil          * rday

  avg_nsoilmin        = avg_nsoilmin        * rday
  avg_psoillab        = avg_psoillab        * rday
  avg_psoilsorb       = avg_psoilsorb       * rday
  avg_psoilocc        = avg_psoilocc        * rday

  avg_rationcsoilmic  = avg_rationcsoilmic  * rday
  avg_rationcsoilslow = avg_rationcsoilslow * rday
  avg_rationcsoilpass = avg_rationcsoilpass * rday

  ! 13C
  if (cable_user%c13o2) then
     avg_c13leaf2met = avg_c13leaf2met * rday
     avg_c13leaf2str = avg_c13leaf2str * rday
     avg_c13root2met = avg_c13root2met * rday
     avg_c13root2str = avg_c13root2str * rday
     avg_c13wood2cwd = avg_c13wood2cwd * rday
  endif

  !write(600,*) 'pmet pre-analytic: ' ,  casapool%plitter(1,metb)
  !write(600,*) 'nmet pre-analytic: ' ,  casapool%nlitter(1,metb)
  !write(600,*) 'csoil3 pre-analytic: ', casapool%csoil(3,:)
  !write(600,*) 'csoil1 pre-analytic: ', casapool%csoil(1,:)
  call analyticpool(kend,veg,soil,casabiome,casapool,                                          &
       casaflux,casamet,casabal,phen,                                         &
       avg_cleaf2met,avg_cleaf2str,avg_croot2met,avg_croot2str,avg_cwood2cwd, &
       avg_nleaf2met,avg_nleaf2str,avg_nroot2met,avg_nroot2str,avg_nwood2cwd, &
       avg_pleaf2met,avg_pleaf2str,avg_proot2met,avg_proot2str,avg_pwood2cwd, &
       avg_cgpp, avg_cnpp, avg_nuptake, avg_puptake,                          &
       avg_xnplimit,avg_xkNlimiting,avg_xklitter,avg_xksoil,                  &
       avg_ratioNCsoilmic,avg_ratioNCsoilslow,avg_ratioNCsoilpass,            &
       avg_nsoilmin,avg_psoillab,avg_psoilsorb,avg_psoilocc, &
       avg_c13leaf2met, avg_c13leaf2str, avg_c13root2met, &
       avg_c13root2str, avg_c13wood2cwd, c13o2pools)
  !write(600,*) 'csoil3 post-analytic: ', casapool%csoil(3,:)
  !write(600,*) 'csoil1 post-analytic: ', casapool%csoil(1,:)

  !write(600,*) 'pmet post analytic: ', avg_pleaf2met, avg_proot2met, casaflux%klitter(1,metb), casapool%plitter(1,metb)
  ! write(600,*) 'nmet post analytic: ', avg_nleaf2met, avg_nroot2met, casaflux%klitter(1,metb), casapool%nlitter(1,metb)

  !!$  call totcnppools(1,veg,casamet,casapool,bmcplant,bmnplant,bmpplant,bmclitter,bmnlitter,bmplitter, &
  !!$       bmcsoil,bmnsoil,bmpsoil,bmnsoilmin,bmpsoillab,bmpsoilsorb,bmpsoilocc,bmarea)

  nloop1 = max(1,mloop-3)
  DO nloop=1, mloop
     !!CLN  OPEN(91,file=fcnpspin)
     !!CLN  read(91,*)
     ctime = 0
     DO nyear=1, myearspin
        !!CLN      read(91,901) ncfile
        !write(*,*) 'spincasa CYEAR', CYEAR, ncfile
        
        
        WRITE(CYEAR,FMT="(I4)") CABLE_USER%CASA_SPIN_STARTYEAR + nyear - 1
        ncfile = TRIM(casafile%c2cdumppath)//'c2c_'//CYEAR//'_dump.nc'
        YYYY = CABLE_USER%CASA_SPIN_STARTYEAR + nyear - 1
        call read_casa_dump(ncfile, casamet, casaflux, phen, climate, c13o2flux, ktau, kend, .TRUE.)

        DO idoy=1,mdyear
           ctime = ctime + 1
           ktauy = idoy*ktauday
           ktau  = (idoy-1)*ktauday +1
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
           phen%phase(:)          = phen%phasespin(:,idoy)
           phen%doyphase(:,1)     = phen%doyphasespin_1(:,idoy)
           phen%doyphase(:,2)     = phen%doyphasespin_2(:,idoy)
           phen%doyphase(:,3)     = phen%doyphasespin_3(:,idoy)
           phen%doyphase(:,4)     = phen%doyphasespin_4(:,idoy)
           climate%qtemp_max_last_year(:) = real(casamet%mtempspin(:,idoy))
           climate%frec(:)                = real(casamet%frecspin(:,idoy))
           ! 13C
           if (cable_user%c13o2) then
              c13o2flux%cAn12(:) = casamet%cAn12spin(:,idoy)
              c13o2flux%cAn(:)   = casamet%cAn13spin(:,idoy)
           endif
           ! BLAZE
           if (cable_user%call_blaze) then
              climate%dprecip(:)      = real(casamet%dprecip_spin(:,idoy))
              climate%aprecip_av20(:) = real(casamet%aprecip_av20_spin(:,idoy))
              climate%du10_max(:)     = real(casamet%du10_max_spin(:,idoy))
              climate%drhum(:)        = real(casamet%drhum_spin(:,idoy))
              climate%dtemp_max(:)    = real(casamet%dtemp_max_spin(:,idoy))
              climate%dtemp_min(:)    = real(casamet%dtemp_max_spin(:,idoy))
              climate%KBDI(:)         = real(casamet%KBDI_spin(:,idoy))
              climate%D_MacArthur(:)  = real(casamet%D_MacArthur_spin(:,idoy))
              climate%FFDI(:)         = real(casamet%FFDI_spin(:,idoy))
              climate%DSLR(:)         = casamet%DSLR_spin(:,idoy)
              climate%last_precip(:)  = real(casamet%last_precip_spin(:,idoy))
           endif

           ! if (nloop==1 .and. nyear==1) then
           !    write(6002, "(200e16.6)") casamet%tairk(3), casamet%tsoil(3,:),  casamet%moist(3,:), &
           !         casaflux%cgpp(3) ,casaflux%crmplant(3,1), real(phen%phase(3)) ,  &
           !         real(phen%doyphase(3,:)), climate%qtemp_max_last_year(3)
           ! endif

           ! 13C
           if (cable_user%c13o2) call c13o2_save_casapool(casapool, casasave)
           call biogeochem(ktauy,dels,idoy,LALLOC,veg,soil,casabiome,casapool,casaflux, &
                casamet,casabal,phen,POP,climate,xnplimit,xkNlimiting,xklitter,xksoil,xkleaf,&
                xkleafcold,xkleafdry,&
                cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,         &
                nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,         &
                pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd)
           ! 13C
           if (cable_user%c13o2) call c13o2_update_pools(casasave, casaflux, c13o2flux, c13o2pools)

           if (cable_user%CALL_BLAZE) then
              CALL BLAZE_ACCOUNTING(BLAZE, climate, ktau, dels, YYYY, idoy)
              
              call blaze_driver(blaze%ncells, blaze, simfire, casapool, casaflux, &
                   casamet, climate, real(shootfrac), idoy, YYYY, 1, POP, veg)

              if (nloop==mloop) then
                 if  (nyear ==myearspin .and. idoy == mdyear  ) then
                     call write_blaze_output_nc( BLAZE, .true.)
                 else
                    call write_blaze_output_nc( BLAZE, .false.)
                 endif
              endif
              
           endif
           
           
           !MC - Should this be nyear==myearspin instead of nyear==1?
           if (nloop==mloop .and. nyear==myearspin) then
              !MC - Should ctime be replaced by idoy?
              CALL WRITE_CASA_OUTPUT_NC( veg, casamet, casapool, casabal, casaflux, &
                   .true., ctime,  idoy.eq.mdyear )
              ! 13C
              if (cable_user%c13o2) then
                 if (idoy == 1) then
                    call c13o2_create_output(casamet, c13o2pools, c13o2_file_id, c13o2_vars, c13o2_var_ids)
                 endif
                 call c13o2_write_output(c13o2_file_id, c13o2_vars, c13o2_var_ids, ctime, c13o2pools)
              end if
           endif
           ! 13C
           if (cable_user%c13o2) then
              if ( (nloop.eq.mloop) .and. (nyear.eq.myearspin) .and. (idoy.eq.mdyear) ) &
                   call c13o2_close_output(c13o2_file_id)
           end if

           IF (cable_user%CALL_POP .and. POP%np.gt.0) THEN ! CALL_POP
              IF (cable_user%CALL_POP) THEN ! accumulate input variables for POP
                 ! accumulate annual variables for use in POP
                 IF (MOD(ktau/ktauday,LOY) == 1) THEN
                    ! (assumes 70% of wood NPP is allocated above ground)
                    casaflux%stemnpp  = casaflux%cnpp * casaflux%fracCalloc(:,2) * 0.7_dp
                    casabal%LAImax    = casamet%glai
                    casabal%Cleafmean = casapool%cplant(:,1) / real(LOY,dp) / 1000.0_dp
                    casabal%Crootmean = casapool%cplant(:,3) / real(LOY,dp) / 1000.0_dp
                 ELSE
                    casaflux%stemnpp  = casaflux%stemnpp + casaflux%cnpp * casaflux%fracCalloc(:,2) * 0.7_dp
                    casabal%LAImax    = max(casamet%glai, casabal%LAImax)
                    casabal%Cleafmean = casabal%Cleafmean + casapool%cplant(:,1) / real(LOY,dp) / 1000.0_dp
                    casabal%Crootmean = casabal%Crootmean + casapool%cplant(:,3) / real(LOY,dp) / 1000.0_dp
                 ENDIF
              ELSE
                 casaflux%stemnpp = 0.0_dp
              ENDIF ! CALL_POP

              IF (idoy==mdyear) THEN ! end of year
                 !MC - ToDo - update 13CO2 harvest
                 CALL POPdriver(casaflux, casabal, veg, POP)
                 ! CALL POP_IO( pop, casamet, NYEAR, 'WRITE_EPI', &
                 !     (nloop.eq.mloop .and. nyear.eq.myearspin) )
                 ! CALL WRITE_CASA_OUTPUT_NC(veg, casamet, casapool, casabal, casaflux, &
                 !     .TRUE., ctime, &
                 !     (nloop.eq.mloop .and. nyear.eq.myearspin.and.idoy.eq.mdyear))
                 ! ctime = ctime + 1
              ENDIF ! end of year
           ELSE
              ! IF(idoy==mdyear) THEN ! end of year
              !     CALL WRITE_CASA_OUTPUT_NC (veg, casamet, casapool, casabal, casaflux, &
              !         .TRUE., ctime, &
              !         (nloop.eq.mloop .and. nyear.eq.myearspin.and.idoy.eq.mdyear)  )
              !     ctime = ctime+1
              ! ENDIF ! end of year
              casaflux%stemnpp = 0.0_dp
           ENDIF ! CALL_POP

        ENDDO ! end of idoy

     ENDDO ! end of nyear

  ENDDO ! end of nloop

  CALL casa_fluxout(CABLE_USER%CASA_SPIN_STARTYEAR+myearspin-1, veg, soil, casabal, casamet)

END SUBROUTINE spincasacnp
