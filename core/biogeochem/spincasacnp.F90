SUBROUTINE spincasacnp(dels,kstart,kend,mloop,veg,soil,casabiome,casapool, &
     casaflux,casamet,casabal,phen,POP,climate,LALLOC, c13o2flux, c13o2pools)


  USE cable_def_types_mod
  USE cable_carbon_module
  USE cable_common_module, ONLY: CABLE_USER
  USE casadimension
  USE casaparm
  USE casavariable
  USE phenvariable
  USE POP_Types,  Only: POP_TYPE
  USE POPMODULE,            ONLY: POPStep
  USE TypeDef,              ONLY: i4b, dp
  use cable_c13o2_def, only: c13o2_pool, c13o2_flux
  use cable_c13o2,     only: c13o2_save_casapool, c13o2_update_pools, &
       c13o2_print_delta_pools

  IMPLICIT NONE
  !!CLN  CHARACTER(LEN=99), INTENT(IN)  :: fcnpspin
  REAL,    INTENT(IN)    :: dels
  INTEGER, INTENT(IN)    :: kstart
  INTEGER, INTENT(IN)    :: kend
  INTEGER, INTENT(IN)    :: mloop
  INTEGER, INTENT(IN)    :: LALLOC
  TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (soil_parameter_type),   INTENT(INOUT) :: soil ! soil parameters
  TYPE (casa_biome),            INTENT(INOUT) :: casabiome
  TYPE (casa_pool),             INTENT(INOUT) :: casapool
  TYPE (casa_flux),             INTENT(INOUT) :: casaflux
  TYPE (casa_met),              INTENT(INOUT) :: casamet
  TYPE (casa_balance),          INTENT(INOUT) :: casabal
  TYPE (phen_variable),         INTENT(INOUT) :: phen
  TYPE (POP_TYPE),              INTENT(INOUT) :: POP
  TYPE (climate_TYPE),          INTENT(INOUT) :: climate
  type(c13o2_flux),             intent(in)    :: c13o2flux
  type(c13o2_pool),             intent(inout) :: c13o2pools

  ! TYPE(casa_met) :: casaspin

  ! local variables
  real,      dimension(:), allocatable, save  :: avg_cleaf2met, avg_cleaf2str, avg_croot2met, avg_croot2str, avg_cwood2cwd
  real,      dimension(:), allocatable, save  :: avg_nleaf2met, avg_nleaf2str, avg_nroot2met, avg_nroot2str, avg_nwood2cwd
  real,      dimension(:), allocatable, save  :: avg_pleaf2met, avg_pleaf2str, avg_proot2met, avg_proot2str, avg_pwood2cwd
  real,      dimension(:), allocatable, save  :: avg_cgpp,      avg_cnpp,      avg_nuptake,   avg_puptake
  real,      dimension(:), allocatable, save  :: avg_nsoilmin,  avg_psoillab,  avg_psoilsorb, avg_psoilocc
  !chris 12/oct/2012 for spin up casa
  real,      dimension(:), allocatable, save  :: avg_ratioNCsoilmic,  avg_ratioNCsoilslow,  avg_ratioNCsoilpass
  real(r_2), dimension(:), allocatable, save  :: avg_xnplimit,  avg_xkNlimiting,avg_xklitter, avg_xksoil

  ! local variables
  INTEGER                  :: myearspin,nyear, nloop1
  CHARACTER(LEN=99)        :: ncfile
  CHARACTER(LEN=4)         :: cyear
  INTEGER                  :: ktau,ktauday,nday,idoy,ktaux,ktauy,nloop, LOY
  INTEGER, save            :: ndays
  real,      dimension(mp)      :: cleaf2met, cleaf2str, croot2met, croot2str, cwood2cwd
  real,      dimension(mp)      :: nleaf2met, nleaf2str, nroot2met, nroot2str, nwood2cwd
  real,      dimension(mp)      :: pleaf2met, pleaf2str, proot2met, proot2str, pwood2cwd
  real,      dimension(mp)      :: xcgpp,     xcnpp,     xnuptake,  xpuptake
  real,      dimension(mp)      :: xnsoilmin, xpsoillab, xpsoilsorb,xpsoilocc
  real(r_2), dimension(mp)      :: xnplimit,  xkNlimiting, xklitter, xksoil,xkleaf, xkleafcold, xkleafdry

  ! more variables to store the spinup pool size over the last 10 loops. Added by Yp Wang 30 Nov 2012
  real,      dimension(5,mvtype,mplant)  :: bmcplant,  bmnplant,  bmpplant
  real,      dimension(5,mvtype,mlitter) :: bmclitter, bmnlitter, bmplitter
  real,      dimension(5,mvtype,msoil)   :: bmcsoil,   bmnsoil,   bmpsoil
  real,      dimension(5,mvtype)         :: bmnsoilmin,bmpsoillab,bmpsoilsorb, bmpsoilocc
  real,      dimension(mvtype)           :: bmarea
  integer nptx,nvt,kloop

   REAL(dp)                               :: StemNPP(mp,2)
   REAL(dp), allocatable, save ::  LAImax(:)    , Cleafmean(:),  Crootmean(:)
   REAL(dp), allocatable :: NPPtoGPP(:)
   INTEGER, allocatable :: Iw(:) ! array of indices corresponding to woody (shrub or forest) tiles
   INTEGER ::ctime

   ! 13C
   real(dp), dimension(c13o2pools%ntile,c13o2pools%npools) :: casasave

   
   if (.NOT.Allocated(LAIMax)) allocate(LAIMax(mp))
   if (.NOT.Allocated(Cleafmean))  allocate(Cleafmean(mp))
   if (.NOT.Allocated(Crootmean)) allocate(Crootmean(mp))
   if (.NOT.Allocated(NPPtoGPP)) allocate(NPPtoGPP(mp))
   if (.NOT.Allocated(Iw)) allocate(Iw(POP%np))

   ctime = 1
   LOY = 365
   !! vh_js !!
    IF (cable_user%CALL_POP) THEN

       Iw = POP%Iwood

    ENDIF

  ktauday=int(24.0*3600.0/dels)
  nday=(kend-kstart+1)/ktauday

  !chris 12/oct/2012 for spin up casa
  IF (.not.(allocated(avg_cleaf2met)))  allocate(avg_cleaf2met(mp), avg_cleaf2str(mp), avg_croot2met(mp), avg_croot2str(mp), &
       avg_cwood2cwd(mp), &
       avg_nleaf2met(mp), avg_nleaf2str(mp), avg_nroot2met(mp), avg_nroot2str(mp), avg_nwood2cwd(mp), &
       avg_pleaf2met(mp), avg_pleaf2str(mp), avg_proot2met(mp), avg_proot2str(mp), avg_pwood2cwd(mp), &
       avg_cgpp(mp),      avg_cnpp(mp),      avg_nuptake(mp),   avg_puptake(mp),                     &
       avg_xnplimit(mp),  avg_xkNlimiting(mp), avg_xklitter(mp), avg_xksoil(mp),                      &
       avg_rationcsoilmic(mp),avg_rationcsoilslow(mp),avg_rationcsoilpass(mp),                        &
       avg_nsoilmin(mp),  avg_psoillab(mp),    avg_psoilsorb(mp), avg_psoilocc(mp))

  !!CLN  OPEN(91, file=fcnpspin)
  !!CLN  read(91,*) myearspin
  myearspin = CABLE_USER%CASA_SPIN_ENDYEAR - CABLE_USER%CASA_SPIN_STARTYEAR + 1
  ! compute the mean fluxes and residence time of each carbon pool
  avg_cleaf2met=0.0; avg_cleaf2str=0.0; avg_croot2met=0.0; avg_croot2str=0.0; avg_cwood2cwd=0.0
  avg_nleaf2met=0.0; avg_nleaf2str=0.0; avg_nroot2met=0.0; avg_nroot2str=0.0; avg_nwood2cwd=0.0
  avg_pleaf2met=0.0; avg_pleaf2str=0.0; avg_proot2met=0.0; avg_proot2str=0.0; avg_pwood2cwd=0.0
  avg_cgpp=0.0;      avg_cnpp=0.0;      avg_nuptake=0.0;   avg_puptake=0.0
  avg_xnplimit=0.0;  avg_xkNlimiting=0.0; avg_xklitter=0.0; avg_xksoil=0.0
  avg_nsoilmin=0.0;  avg_psoillab=0.0;    avg_psoilsorb=0.0; avg_psoilocc=0.0
  avg_rationcsoilmic=0.0;avg_rationcsoilslow=0.0;avg_rationcsoilpass=0.0

write(600,*) 'csoil3 init: ', casapool%csoil(3,:)
  write(600,*) 'csoil1 init: ', casapool%csoil(1,:)
  
  do nyear=1,myearspin
     !     read(91,901) ncfile
     WRITE(CYEAR,FMT="(I4)") CABLE_USER%CASA_SPIN_STARTYEAR + nyear - 1
     ncfile = TRIM(casafile%c2cdumppath)//'c2c_'//CYEAR//'_dump.nc'


     call read_casa_dump( ncfile, casamet, casaflux, phen, climate, c13o2flux, ktau, kend, .TRUE. )
     !!CLN901  format(A99)
     do idoy=1,mdyear
        ktau=(idoy-1)*ktauday +1

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
        phen%phase(:) = phen%phasespin(:,idoy)
        phen%doyphase(:,1) = phen%doyphasespin_1(:,idoy)
        phen%doyphase(:,2) =  phen%doyphasespin_2(:,idoy)
        phen%doyphase(:,3) =  phen%doyphasespin_3(:,idoy)
        phen%doyphase(:,4) =  phen%doyphasespin_4(:,idoy)
        climate%qtemp_max_last_year(:) =  casamet%mtempspin(:,idoy)
        if (cable_user%c13o2) then
           c13o2flux%cAn12(:) = casamet%cAn12spin(:,idoy)
           c13o2flux%cAn13(:) = casamet%cAn13spin(:,idoy)
        endif

        if (cable_user%c13o2) call c13o2_save_casapool(casapool, casasave)
        if (cable_user%c13o2) then
           write(*,*) '13C in spincasacnp - 01'
           call c13o2_print_delta_pools(casapool, casaflux, c13o2pools)
        endif
        CALL biogeochem(ktau,dels,idoy,LALLOC,veg,soil,casabiome,casapool,casaflux, &
             casamet,casabal,phen,POP,climate,xnplimit,xkNlimiting,xklitter, &
             xksoil,xkleaf,xkleafcold,xkleafdry,&
             cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,         &
             nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,         &
             pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd)
        if (cable_user%c13o2) call c13o2_update_pools(casasave, casaflux, c13o2flux, c13o2pools)
        if (cable_user%c13o2) then
           write(*,*) '13C in spincasacnp - 02'
           call c13o2_print_delta_pools(casapool, casaflux, c13o2pools)
        endif
         
        IF (cable_user%CALL_POP .and. POP%np.gt.0) THEN ! CALL_POP

           IF (cable_user%CALL_POP) THEN ! accumulate input variables for POP
               ! accumulate annual variables for use in POP
               IF(MOD(ktau/ktauday,LOY)==1 ) THEN
                  casaflux%stemnpp =  casaflux%cnpp * casaflux%fracCalloc(:,2) * 0.7 ! (assumes 70% of wood NPP is allocated above ground)
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

 !CALL WRITE_CASA_OUTPUT_NC (veg, casamet, casapool, casabal, casaflux, &
 !            .true., ctime, .FALSE.  )
 !            ctime = ctime+1
           IF(idoy==mdyear) THEN ! end of year

              CALL POPdriver(casaflux,casabal,veg, POP)

             ! CALL POP_IO( pop, casamet, nyear, 'WRITE_EPI', &
             ! 		 (.FALSE.))
             !CALL WRITE_CASA_OUTPUT_NC (veg, casamet, casapool, casabal, casaflux, &
             ! .true., ctime, .FALSE.  )
             !ctime = ctime+1

           ENDIF  ! end of year
        
        ELSE
           IF(idoy==mdyear) THEN ! end of year

             !CALL WRITE_CASA_OUTPUT_NC (veg, casamet, casapool, casabal, casaflux, &
             ! .true., ctime, .FALSE.  )
             ctime = ctime+1

           ENDIF  ! end of year
           
           casaflux%stemnpp = 0.
        ENDIF ! CALL_POP
      
!!$        WHERE(xkNlimiting .eq. 0)  !Chris Lu 4/June/2012
!!$           xkNlimiting = 0.001
!!$        END WHERE
        nptx=8173

        avg_cleaf2met = avg_cleaf2met + cleaf2met
        avg_cleaf2str = avg_cleaf2str + cleaf2str
        avg_croot2met = avg_croot2met + croot2met
        avg_croot2str = avg_croot2str + croot2str
        avg_cwood2cwd = avg_cwood2cwd + cwood2cwd

        avg_nleaf2met = avg_nleaf2met + nleaf2met
        avg_nleaf2str = avg_nleaf2str + nleaf2str
        avg_nroot2met = avg_nroot2met + nroot2met
        avg_nroot2str = avg_nroot2str + nroot2str
        avg_nwood2cwd = avg_nwood2cwd + nwood2cwd

        avg_pleaf2met = avg_pleaf2met + pleaf2met
        avg_pleaf2str = avg_pleaf2str + pleaf2str
        avg_proot2met = avg_proot2met + proot2met
        avg_proot2str = avg_proot2str + proot2str
        avg_pwood2cwd = avg_pwood2cwd + pwood2cwd

        avg_cgpp      = avg_cgpp      + casaflux%cgpp
        avg_cnpp      = avg_cnpp      + casaflux%cnpp
        avg_nuptake   = avg_nuptake   + casaflux%Nminuptake
        avg_puptake   = avg_puptake   + casaflux%Plabuptake

        avg_xnplimit    = avg_xnplimit    + xnplimit
        avg_xkNlimiting = avg_xkNlimiting + xkNlimiting
        avg_xklitter    = avg_xklitter    + xklitter
        avg_xksoil      = avg_xksoil      + xksoil

        avg_nsoilmin    = avg_nsoilmin    + casapool%nsoilmin
        avg_psoillab    = avg_psoillab    + casapool%psoillab
        avg_psoilsorb   = avg_psoilsorb   + casapool%psoilsorb
        avg_psoilocc    = avg_psoilocc    + casapool%psoilocc

        avg_rationcsoilmic  = avg_rationcsoilmic  + casapool%ratioNCsoilnew(:,mic)
        avg_rationcsoilslow = avg_rationcsoilslow + casapool%ratioNCsoilnew(:,slow)
        avg_rationcsoilpass = avg_rationcsoilpass + casapool%ratioNCsoilnew(:,pass)
     enddo
  enddo

  !!CLN    CLOSE(91)

  avg_cleaf2met = avg_cleaf2met/real(nday*myearspin)
  avg_cleaf2str = avg_cleaf2str/real(nday*myearspin)
  avg_croot2met = avg_croot2met/real(nday*myearspin)
  avg_croot2str = avg_croot2str/real(nday*myearspin)
  avg_cwood2cwd = avg_cwood2cwd/real(nday*myearspin)

  avg_nleaf2met = avg_nleaf2met/real(nday*myearspin)
  avg_nleaf2str = avg_nleaf2str/real(nday*myearspin)
  avg_nroot2met = avg_nroot2met/real(nday*myearspin)
  avg_nroot2str = avg_nroot2str/real(nday*myearspin)
  avg_nwood2cwd = avg_nwood2cwd/real(nday*myearspin)

  avg_pleaf2met = avg_pleaf2met/real(nday*myearspin)
  avg_pleaf2str = avg_pleaf2str/real(nday*myearspin)
  avg_proot2met = avg_proot2met/real(nday*myearspin)
  avg_proot2str = avg_proot2str/real(nday*myearspin)
  avg_pwood2cwd = avg_pwood2cwd/real(nday*myearspin)

  avg_cgpp      = avg_cgpp/real(nday*myearspin)
  avg_cnpp      = avg_cnpp/real(nday*myearspin)
  avg_nuptake   = avg_nuptake/real(nday*myearspin)
  avg_puptake   = avg_puptake/real(nday*myearspin)

  avg_xnplimit    = avg_xnplimit/real(nday*myearspin)
  avg_xkNlimiting = avg_xkNlimiting/real(nday*myearspin)
  avg_xklitter    = avg_xklitter/real(nday*myearspin)
  avg_xksoil      = avg_xksoil/real(nday*myearspin)

  avg_nsoilmin    = avg_nsoilmin/real(nday*myearspin)
  avg_psoillab    = avg_psoillab/real(nday*myearspin)
  avg_psoilsorb   = avg_psoilsorb/real(nday*myearspin)
  avg_psoilocc    = avg_psoilocc/real(nday*myearspin)

  avg_rationcsoilmic  = avg_rationcsoilmic  /real(nday*myearspin)
  avg_rationcsoilslow = avg_rationcsoilslow /real(nday*myearspin)
  avg_rationcsoilpass = avg_rationcsoilpass /real(nday*myearspin)

  !write(600,*) 'pmet pre-analytic: ' ,  casapool%plitter(1,metb)
  !write(600,*) 'nmet pre-analytic: ' ,  casapool%nlitter(1,metb)
  write(600,*) 'csoil3 pre-analytic: ', casapool%csoil(3,:)
  write(600,*) 'csoil1 pre-analytic: ', casapool%csoil(1,:)
  call analyticpool(kend,veg,soil,casabiome,casapool,                                          &
       casaflux,casamet,casabal,phen,                                         &
       avg_cleaf2met,avg_cleaf2str,avg_croot2met,avg_croot2str,avg_cwood2cwd, &
       avg_nleaf2met,avg_nleaf2str,avg_nroot2met,avg_nroot2str,avg_nwood2cwd, &
       avg_pleaf2met,avg_pleaf2str,avg_proot2met,avg_proot2str,avg_pwood2cwd, &
       avg_cgpp, avg_cnpp, avg_nuptake, avg_puptake,                          &
       avg_xnplimit,avg_xkNlimiting,avg_xklitter,avg_xksoil,                  &
       avg_ratioNCsoilmic,avg_ratioNCsoilslow,avg_ratioNCsoilpass,            &
       avg_nsoilmin,avg_psoillab,avg_psoilsorb,avg_psoilocc)
  write(600,*) 'csoil3 post-analytic: ', casapool%csoil(3,:)
  write(600,*) 'csoil1 post-analytic: ', casapool%csoil(1,:)

  !write(600,*) 'pmet post analytic: ', avg_pleaf2met, avg_proot2met, casaflux%klitter(1,metb), casapool%plitter(1,metb)
  ! write(600,*) 'nmet post analytic: ', avg_nleaf2met, avg_nroot2met, casaflux%klitter(1,metb), casapool%nlitter(1,metb)

!!$  call totcnppools(1,veg,casamet,casapool,bmcplant,bmnplant,bmpplant,bmclitter,bmnlitter,bmplitter, &
!!$       bmcsoil,bmnsoil,bmpsoil,bmnsoilmin,bmpsoillab,bmpsoilsorb,bmpsoilocc,bmarea)

  nloop1= max(1,mloop-3)

  DO nloop=1,mloop

     !!CLN  OPEN(91,file=fcnpspin)
     !!CLN  read(91,*)
     DO nyear=1,myearspin
        !!CLN      read(91,901) ncfile
        !write(*,*) 'spincasa CYEAR', CYEAR, ncfile
        WRITE(CYEAR,FMT="(I4)") CABLE_USER%CASA_SPIN_STARTYEAR + nyear - 1
        ncfile = TRIM(casafile%c2cdumppath)//'c2c_'//CYEAR//'_dump.nc'
        call read_casa_dump( ncfile, casamet, casaflux, phen, climate, c13o2flux, ktau, kend, .TRUE. )

        DO idoy=1,mdyear
           ktauy=idoy*ktauday
           ktau=(idoy-1)*ktauday +1
      
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
           phen%phase(:) = phen%phasespin(:,idoy)
           phen%doyphase(:,1) = phen%doyphasespin_1(:,idoy)
           phen%doyphase(:,2) =  phen%doyphasespin_2(:,idoy)
           phen%doyphase(:,3) =  phen%doyphasespin_3(:,idoy)
           phen%doyphase(:,4) =  phen%doyphasespin_4(:,idoy)
           climate%qtemp_max_last_year(:) =  casamet%mtempspin(:,idoy)
           if (cable_user%c13o2) then
              c13o2flux%cAn12(:) = casamet%cAn12spin(:,idoy)
              c13o2flux%cAn13(:) = casamet%cAn13spin(:,idoy)
           endif

           if (nloop==1 .and. nyear==1) then
              write(6002, "(200e16.6)") casamet%tairk(3), casamet%tsoil(3,:),  casamet%moist(3,:), &
                   casaflux%cgpp(3) ,casaflux%crmplant(3,1), real(phen%phase(3)) ,  &
                    real(phen%doyphase(3,:)), climate%qtemp_max_last_year(3)
                   

           endif

           if (cable_user%c13o2) call c13o2_save_casapool(casapool, casasave)
           if (cable_user%c13o2) then
              write(*,*) '13C in spincasacnp - 03'
              call c13o2_print_delta_pools(casapool, casaflux, c13o2pools)
           endif
           call biogeochem(ktauy,dels,idoy,LALLOC,veg,soil,casabiome,casapool,casaflux, &
                casamet,casabal,phen,POP,climate,xnplimit,xkNlimiting,xklitter,xksoil,xkleaf,&
                xkleafcold,xkleafdry,&
                cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,         &
                nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,         &
                pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd)
           if (cable_user%c13o2) call c13o2_update_pools(casasave, casaflux, c13o2flux, c13o2pools)
           if (cable_user%c13o2) then
              write(*,*) '13C in spincasacnp - 04'
              call c13o2_print_delta_pools(casapool, casaflux, c13o2pools)
           endif

            if (nloop==mloop .and. nyear==1) then
              CALL WRITE_CASA_OUTPUT_NC (veg, casamet, casapool, casabal, casaflux, &
                                             .TRUE., ctime, &
                             (nloop.eq.mloop .and. nyear.eq.myearspin.and.idoy.eq.mdyear)  )
              ctime = ctime+1
            endif

           IF (cable_user%CALL_POP .and. POP%np.gt.0) THEN ! CALL_POP

              IF (cable_user%CALL_POP) THEN ! accumulate input variables for POP
                 ! accumulate annual variables for use in POP
                 IF(MOD(ktau/ktauday,LOY)==1 ) THEN
                    casaflux%stemnpp =  casaflux%cnpp * casaflux%fracCalloc(:,2) * 0.7 ! (assumes 70% of wood NPP is allocated above ground)
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
              
           
              IF(idoy==mdyear) THEN ! end of year

                 CALL POPdriver(casaflux,casabal,veg, POP)
                 !CALL POP_IO( pop, casamet, NYEAR, 'WRITE_EPI', &
                ! 	 (nloop.eq.mloop .and. nyear.eq.myearspin) )
                 !CALL WRITE_CASA_OUTPUT_NC (veg, casamet, casapool, casabal, casaflux, &
                 !                            .TRUE., ctime, &
                  !           (nloop.eq.mloop .and. nyear.eq.myearspin.and.idoy.eq.mdyear)  )
                 ctime = ctime+1               
                 
              ENDIF  ! end of year

           ELSE

           !IF(idoy==mdyear) THEN ! end of year
               
               !  CALL WRITE_CASA_OUTPUT_NC (veg, casamet, casapool, casabal, casaflux, &
               !                              .TRUE., ctime, &
               !              (nloop.eq.mloop .and. nyear.eq.myearspin.and.idoy.eq.mdyear)  )
               !  ctime = ctime+1               
                 
           ! ENDIF  ! end of year   
           casaflux%stemnpp = 0.
        ENDIF ! CALL_POP
        
     
     ENDDO   ! end of idoy
  ENDDO   ! end of nyear

  
!!$  if(nloop>=nloop1) &
!!$       call totcnppools(2+nloop-nloop1,veg,casamet,casapool,bmcplant,bmnplant,bmpplant,bmclitter,bmnlitter,bmplitter, &
!!$       bmcsoil,bmnsoil,bmpsoil,bmnsoilmin,bmpsoillab,bmpsoilsorb,bmpsoilocc,bmarea)

ENDDO     ! end of nloop

CALL casa_fluxout(CABLE_USER%CASA_SPIN_STARTYEAR + myearspin - 1 , veg, soil, casabal, casamet)

!STOP
write(600,*) 'csoil3 end: ', casapool%csoil(3,:)
  write(600,*) 'csoil1 end: ', casapool%csoil(1,:)
  
!!$! write the last five loop pool size by PFT type
!!$open(92,file='cnpspinlast5.txt')
!!$write(92,921)
!!$921 format('PFT total area in 10**12 m2', f12.4)
!!$  do nvt=1,mvtype
!!$     write(92,*) bmarea(nvt)
!!$  enddo
!!$
!!$  do nvt=1,mvtype
!!$     if(bmarea(nvt) >0.0) then
!!$        do kloop=1,5
!!$           write(92,922) nvt, bmcplant(kloop,nvt,:),bmclitter(kloop,nvt,:),bmcsoil(kloop,nvt,:)
!!$        enddo
!!$        if (icycle >1) then
!!$           do kloop=1,5
!!$              write(92,922) nvt, bmnplant(kloop,nvt,:),bmnlitter(kloop,nvt,:),bmnsoil(kloop,nvt,:), bmnsoilmin(kloop,nvt)
!!$           enddo
!!$        endif
!!$
!!$        if(icycle >2) then
!!$           do kloop=1,5
!!$              write(92,922) nvt, bmpplant(kloop,nvt,:),bmplitter(kloop,nvt,:),bmpsoil(kloop,nvt,:),  &
!!$                   bmpsoillab(kloop,nvt), bmpsoilsorb(kloop,nvt), bmpsoilocc(kloop,nvt)
!!$           enddo
!!$        endif
!!$     endif
!!$  enddo
!!$922 format(i4,20(f10.4,2x))
!!$  CLOSE(92)


151 FORMAT(i6,100(f12.5,2x))
END SUBROUTINE spincasacnp
