SUBROUTINE spincasacnp( dels,kstart,kend,mloop,veg,soil,casabiome,casapool, &
     casaflux,casamet,casabal,phen,POP,climate,LALLOC )


  USE cable_def_types_mod
  USE cable_carbon_module
  USE cable_common_module, ONLY: CABLE_USER
  USE casadimension
  USE casaparm
  USE casa_cable !jhan:also put this in mod
  USE casa_inout_module
  USE casavariable
  USE phenvariable
  USE POP_Types,  ONLY: POP_TYPE
  USE POPMODULE,            ONLY: POPStep
  USE TypeDef,              ONLY: i4b, dp

  IMPLICIT NONE
  !CLN  CHARACTER(LEN=99), INTENT(IN)  :: fcnpspin
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
  TYPE (POP_TYPE), INTENT(INOUT)     :: POP
  TYPE (climate_TYPE), INTENT(INOUT)     :: climate


  TYPE (casa_met)  :: casaspin

  ! local variables
  REAL,      DIMENSION(:), ALLOCATABLE, SAVE  :: avg_cleaf2met, avg_cleaf2str, avg_croot2met, avg_croot2str, avg_cwood2cwd
  REAL,      DIMENSION(:), ALLOCATABLE, SAVE  :: avg_nleaf2met, avg_nleaf2str, avg_nroot2met, avg_nroot2str, avg_nwood2cwd
  REAL,      DIMENSION(:), ALLOCATABLE, SAVE  :: avg_pleaf2met, avg_pleaf2str, avg_proot2met, avg_proot2str, avg_pwood2cwd
  REAL,      DIMENSION(:), ALLOCATABLE, SAVE  :: avg_cgpp,      avg_cnpp,      avg_nuptake,   avg_puptake
  REAL,      DIMENSION(:), ALLOCATABLE, SAVE  :: avg_nsoilmin,  avg_psoillab,  avg_psoilsorb, avg_psoilocc
  !chris 12/oct/2012 for spin up casa
  REAL,      DIMENSION(:), ALLOCATABLE, SAVE  :: avg_ratioNCsoilmic,  avg_ratioNCsoilslow,  avg_ratioNCsoilpass
  REAL(r_2), DIMENSION(:), ALLOCATABLE, SAVE  :: avg_xnplimit,  avg_xkNlimiting,avg_xklitter, avg_xksoil

  REAL,      DIMENSION(:), ALLOCATABLE, SAVE  :: avg_af
  REAL,      DIMENSION(:), ALLOCATABLE, SAVE  :: avg_aw
  REAL,      DIMENSION(:), ALLOCATABLE, SAVE  :: avg_ar
  REAL,      DIMENSION(:), ALLOCATABLE, SAVE  :: avg_lf
  REAL,      DIMENSION(:), ALLOCATABLE, SAVE  :: avg_lw
  REAL,      DIMENSION(:), ALLOCATABLE, SAVE  :: avg_lr
  REAL,      DIMENSION(:), ALLOCATABLE, SAVE  :: avg_annual_cnpp

  ! local variables
  INTEGER                  :: myearspin,nyear, nloop1
  CHARACTER(LEN=99)        :: ncfile
  CHARACTER(LEN=4)         :: cyear
  INTEGER                  :: ktau,ktauday,nday,idoy,ktaux,ktauy,nloop, LOY
  INTEGER, SAVE            :: ndays
  REAL,      DIMENSION(mp)      :: cleaf2met, cleaf2str, croot2met, croot2str, cwood2cwd
  REAL,      DIMENSION(mp)      :: nleaf2met, nleaf2str, nroot2met, nroot2str, nwood2cwd
  REAL,      DIMENSION(mp)      :: pleaf2met, pleaf2str, proot2met, proot2str, pwood2cwd
  REAL,      DIMENSION(mp)      :: xcgpp,     xcnpp,     xnuptake,  xpuptake
  REAL,      DIMENSION(mp)      :: xnsoilmin, xpsoillab, xpsoilsorb,xpsoilocc
  REAL(r_2), DIMENSION(mp)      :: xnplimit,  xkNlimiting, xklitter, xksoil,xkleaf, xkleafcold, xkleafdry

  ! more variables to store the spinup pool size over the last 10 loops. Added by Yp Wang 30 Nov 2012
  REAL,      DIMENSION(5,mvtype,mplant)  :: bmcplant,  bmnplant,  bmpplant
  REAL,      DIMENSION(5,mvtype,mlitter) :: bmclitter, bmnlitter, bmplitter
  REAL,      DIMENSION(5,mvtype,msoil)   :: bmcsoil,   bmnsoil,   bmpsoil
  REAL,      DIMENSION(5,mvtype)         :: bmnsoilmin,bmpsoillab,bmpsoilsorb, bmpsoilocc
  REAL,      DIMENSION(mvtype)           :: bmarea
  INTEGER nptx,nvt,kloop

  REAL(dp)                               :: StemNPP(mp,2)
  REAL(dp), ALLOCATABLE, SAVE ::  LAImax(:)    , Cleafmean(:),  Crootmean(:)
  REAL(dp), ALLOCATABLE :: NPPtoGPP(:)
  INTEGER, ALLOCATABLE :: Iw(:) ! array of indices corresponding to woody (shrub or forest) tiles

  IF (.NOT.ALLOCATED(LAIMax)) ALLOCATE(LAIMax(mp))
  IF (.NOT.ALLOCATED(Cleafmean))  ALLOCATE(Cleafmean(mp))
  IF (.NOT.ALLOCATED(Crootmean)) ALLOCATE(Crootmean(mp))
  IF (.NOT.ALLOCATED(NPPtoGPP)) ALLOCATE(NPPtoGPP(mp))
  IF (.NOT.ALLOCATED(Iw)) ALLOCATE(Iw(POP%np))

  LOY = 365
  ! vh_js !
  IF (cable_user%CALL_POP) THEN

     Iw = POP%Iwood

  ENDIF

  ktauday=INT(24.0*3600.0/dels)
  nday=(kend-kstart+1)/ktauday

  !chris 12/oct/2012 for spin up casa
  IF (.NOT.(ALLOCATED(avg_cleaf2met)))  ALLOCATE(avg_cleaf2met(mp), avg_cleaf2str(mp), avg_croot2met(mp), avg_croot2str(mp), &
       avg_cwood2cwd(mp), &
       avg_nleaf2met(mp), avg_nleaf2str(mp), avg_nroot2met(mp), avg_nroot2str(mp), avg_nwood2cwd(mp), &
       avg_pleaf2met(mp), avg_pleaf2str(mp), avg_proot2met(mp), avg_proot2str(mp), avg_pwood2cwd(mp), &
       avg_cgpp(mp),      avg_cnpp(mp),      avg_nuptake(mp),   avg_puptake(mp),                     &
       avg_xnplimit(mp),  avg_xkNlimiting(mp), avg_xklitter(mp), avg_xksoil(mp),                      &
       avg_rationcsoilmic(mp),avg_rationcsoilslow(mp),avg_rationcsoilpass(mp),                        &
       avg_nsoilmin(mp),  avg_psoillab(mp),    avg_psoilsorb(mp), avg_psoilocc(mp))

  ALLOCATE(avg_af(mp))
  ALLOCATE(avg_aw(mp))
  ALLOCATE(avg_ar(mp))
  ALLOCATE(avg_lf(mp))
  ALLOCATE(avg_lw(mp))
  ALLOCATE(avg_lr(mp))
  ALLOCATE(avg_annual_cnpp(mp))
  avg_af = 0.0
  avg_aw = 0.0
  avg_ar = 0.0
  avg_lf = 0.0
  avg_lw = 0.0
  avg_lr = 0.0
  avg_annual_cnpp = 0.0

  !CLN  OPEN(91, file=fcnpspin)
  !CLN  read(91,*) myearspin
  myearspin = CABLE_USER%CASA_SPIN_ENDYEAR - CABLE_USER%CASA_SPIN_STARTYEAR + 1
  ! compute the mean fluxes and residence time of each carbon pool
  avg_cleaf2met=0.0; avg_cleaf2str=0.0; avg_croot2met=0.0; avg_croot2str=0.0; avg_cwood2cwd=0.0
  avg_nleaf2met=0.0; avg_nleaf2str=0.0; avg_nroot2met=0.0; avg_nroot2str=0.0; avg_nwood2cwd=0.0
  avg_pleaf2met=0.0; avg_pleaf2str=0.0; avg_proot2met=0.0; avg_proot2str=0.0; avg_pwood2cwd=0.0
  avg_cgpp=0.0;      avg_cnpp=0.0;      avg_nuptake=0.0;   avg_puptake=0.0
  avg_xnplimit=0.0;  avg_xkNlimiting=0.0; avg_xklitter=0.0; avg_xksoil=0.0
  avg_nsoilmin=0.0;  avg_psoillab=0.0;    avg_psoilsorb=0.0; avg_psoilocc=0.0
  avg_rationcsoilmic=0.0;avg_rationcsoilslow=0.0;avg_rationcsoilpass=0.0

  DO nyear=1,myearspin
     !     read(91,901) ncfile
     WRITE(CYEAR,FMT="(I4)") CABLE_USER%CASA_SPIN_STARTYEAR + nyear - 1
     ncfile = TRIM(casafile%c2cdumppath)//'c2c_'//CYEAR//'_dump.nc'


     CALL read_casa_dump( ncfile,casamet, casaflux, phen,climate, ktau ,kend,.TRUE. )
     !CLN901  format(A99)
     DO idoy=1,mdyear
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
        IF (cable_user%call_climate) THEN
           climate%qtemp_max_last_year(:) =  casamet%mtempspin(:,idoy)
        ENDIF
        ! write(6699,*) casaflux%cgpp(1), climate%mtemp(1),  casaflux%crmplant(1,1)

        CALL biogeochem(ktau,dels,idoy,LALLOC,veg,soil,casabiome,casapool,casaflux, &
             casamet,casabal,phen,POP,climate,xnplimit,xkNlimiting,xklitter, &
             xksoil,xkleaf,xkleafcold,xkleafdry,&
             cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,         &
             nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,         &
             pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd)


        IF (cable_user%CALL_POP .AND. POP%np.GT.0) THEN ! CALL_POP

           IF (cable_user%CALL_POP) THEN ! accumulate input variables for POP
              ! accumulate annual variables for use in POP
              IF(MOD(ktau/ktauday,LOY)==1 ) THEN
                 casaflux%stemnpp =  casaflux%cnpp * casaflux%fracCalloc(:,2) * 0.7 ! (assumes 70% of wood NPP is allocated above ground)
                 casabal%LAImax = casamet%glai
                 casabal%Cleafmean = casapool%cplant(:,1)/REAL(LOY)/1000.
                 casabal%Crootmean = casapool%cplant(:,3)/REAL(LOY)/1000.
              ELSE
                 casaflux%stemnpp = casaflux%stemnpp + casaflux%cnpp * casaflux%fracCalloc(:,2) * 0.7
                 casabal%LAImax = MAX(casamet%glai, casabal%LAImax)
                 casabal%Cleafmean = casabal%Cleafmean + casapool%cplant(:,1)/REAL(LOY)/1000.
                 casabal%Crootmean = casabal%Crootmean + casapool%cplant(:,3)/REAL(LOY)/1000.
              ENDIF
           ELSE
              casaflux%stemnpp = 0.
           ENDIF ! CALL_POP


           IF(idoy==mdyear) THEN ! end of year

              CALL POPdriver(casaflux,casabal,veg, POP)

           ENDIF  ! end of year
        ELSE
           casaflux%stemnpp = 0.
        ENDIF ! CALL_POP


!$        WHERE(xkNlimiting .eq. 0)  !Chris Lu 4/June/2012
!$           xkNlimiting = 0.001
!$        END WHERE
        nptx=8173

        ! Calculate average allocation fractions  (-) for the plant pools
        avg_af = avg_af + casaflux%fracCalloc(:,leaf)
        avg_aw = avg_aw + casaflux%fracCalloc(:,wood)
        avg_ar = avg_ar + casaflux%fracCalloc(:,froot)

        ! Calculate average turnover rates for the plant pools (yr-1)
        avg_lf = avg_lf + (casaflux%kplant(:,leaf) * REAL(LOY))
        avg_lw = avg_lw + (casaflux%kplant(:,wood) * REAL(LOY))
        avg_lr = avg_lr + (casaflux%kplant(:,froot) * REAL(LOY))

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
     ENDDO
  ENDDO

  ! Average the plant allocation fraction
  avg_af = avg_af / REAL(nday)
  avg_aw = avg_aw / REAL(nday)
  avg_ar = avg_ar / REAL(nday)

  ! Average the plant turnover fraction
  avg_lf = avg_lf / REAL(nday)
  avg_lw = avg_lw / REAL(nday)
  avg_lr = avg_lr / REAL(nday)

  ! Need the annual NPP to solve plant pools g C m-2 y-1
  avg_annual_cnpp = avg_cnpp / REAL(myearspin)

  avg_cleaf2met = avg_cleaf2met/REAL(nday)
  avg_cleaf2str = avg_cleaf2str/REAL(nday)
  avg_croot2met = avg_croot2met/REAL(nday)
  avg_croot2str = avg_croot2str/REAL(nday)
  avg_cwood2cwd = avg_cwood2cwd/REAL(nday)

  avg_nleaf2met = avg_nleaf2met/REAL(nday)
  avg_nleaf2str = avg_nleaf2str/REAL(nday)
  avg_nroot2met = avg_nroot2met/REAL(nday)
  avg_nroot2str = avg_nroot2str/REAL(nday)
  avg_nwood2cwd = avg_nwood2cwd/REAL(nday)

  avg_pleaf2met = avg_pleaf2met/REAL(nday)
  avg_pleaf2str = avg_pleaf2str/REAL(nday)
  avg_proot2met = avg_proot2met/REAL(nday)
  avg_proot2str = avg_proot2str/REAL(nday)
  avg_pwood2cwd = avg_pwood2cwd/REAL(nday)

  avg_cgpp      = avg_cgpp/REAL(nday)
  avg_cnpp      = avg_cnpp/REAL(nday)

  avg_nuptake   = avg_nuptake/REAL(nday)
  avg_puptake   = avg_puptake/REAL(nday)

  avg_xnplimit    = avg_xnplimit/REAL(nday)
  avg_xkNlimiting = avg_xkNlimiting/REAL(nday)
  avg_xklitter    = avg_xklitter/REAL(nday)

  avg_xksoil      = avg_xksoil/REAL(nday)

  avg_nsoilmin    = avg_nsoilmin/REAL(nday)
  avg_psoillab    = avg_psoillab/REAL(nday)
  avg_psoilsorb   = avg_psoilsorb/REAL(nday)
  avg_psoilocc    = avg_psoilocc/REAL(nday)

  avg_rationcsoilmic  = avg_rationcsoilmic  /REAL(nday)
  avg_rationcsoilslow = avg_rationcsoilslow /REAL(nday)
  avg_rationcsoilpass = avg_rationcsoilpass /REAL(nday)

  CALL analyticpool(kend,veg,soil,casabiome,casapool,                                          &
       casaflux,casamet,casabal,phen,                                         &
       avg_cleaf2met,avg_cleaf2str,avg_croot2met,avg_croot2str,avg_cwood2cwd, &
       avg_nleaf2met,avg_nleaf2str,avg_nroot2met,avg_nroot2str,avg_nwood2cwd, &
       avg_pleaf2met,avg_pleaf2str,avg_proot2met,avg_proot2str,avg_pwood2cwd, &
       avg_cgpp, avg_cnpp, avg_nuptake, avg_puptake,                          &
       avg_xnplimit,avg_xkNlimiting,avg_xklitter,avg_xksoil,                  &
       avg_ratioNCsoilmic,avg_ratioNCsoilslow,avg_ratioNCsoilpass,            &
       avg_nsoilmin,avg_psoillab,avg_psoilsorb,avg_psoilocc,                  &
       avg_af, avg_aw, avg_ar, avg_lf, avg_lw, avg_lr, avg_annual_cnpp)

!$  call totcnppools(1,veg,casamet,casapool,bmcplant,bmnplant,bmpplant,bmclitter,bmnlitter,bmplitter, &
!$       bmcsoil,bmnsoil,bmpsoil,bmnsoilmin,bmpsoillab,bmpsoilsorb,bmpsoilocc,bmarea)

  nloop1= MAX(1,mloop-3)

  DO nloop=1,mloop

     !CLN  OPEN(91,file=fcnpspin)
     !CLN  read(91,*)
     DO nyear=1,myearspin
        !CLN      read(91,901) ncfile
        !write(*,*) 'spincasa CYEAR', CYEAR, ncfile
        WRITE(CYEAR,FMT="(I4)") CABLE_USER%CASA_SPIN_STARTYEAR + nyear - 1
        ncfile = TRIM(casafile%c2cdumppath)//'c2c_'//CYEAR//'_dump.nc'
        CALL read_casa_dump( ncfile, casamet, casaflux, phen,climate, ktau, kend, .TRUE. )

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
           IF (cable_user%call_climate) THEN
              climate%qtemp_max_last_year(:) =  casamet%mtempspin(:,idoy)
           ENDIF



           CALL biogeochem(ktauy,dels,idoy,LALLOC,veg,soil,casabiome,casapool,casaflux, &
                casamet,casabal,phen,POP,climate,xnplimit,xkNlimiting,xklitter,xksoil,xkleaf,&
                xkleafcold,xkleafdry,&
                cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,         &
                nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,         &
                pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd)



           IF (cable_user%CALL_POP .AND. POP%np.GT.0) THEN ! CALL_POP

              IF (cable_user%CALL_POP) THEN ! accumulate input variables for POP
                 ! accumulate annual variables for use in POP
                 IF(MOD(ktau/ktauday,LOY)==1 ) THEN
                    casaflux%stemnpp =  casaflux%cnpp * casaflux%fracCalloc(:,2) * 0.7 ! (assumes 70% of wood NPP is allocated above ground)
                    casabal%LAImax = casamet%glai
                    casabal%Cleafmean = casapool%cplant(:,1)/REAL(LOY)/1000.
                    casabal%Crootmean = casapool%cplant(:,3)/REAL(LOY)/1000.
                 ELSE
                    casaflux%stemnpp = casaflux%stemnpp + casaflux%cnpp * casaflux%fracCalloc(:,2) * 0.7
                    casabal%LAImax = MAX(casamet%glai, casabal%LAImax)
                    casabal%Cleafmean = casabal%Cleafmean + casapool%cplant(:,1)/REAL(LOY)/1000.
                    casabal%Crootmean = casabal%Crootmean + casapool%cplant(:,3)/REAL(LOY)/1000.
                 ENDIF
              ELSE
                 casaflux%stemnpp = 0.
              ENDIF ! CALL_POP


              IF(idoy==mdyear) THEN ! end of year

                 CALL POPdriver(casaflux,casabal,veg, POP)


              ENDIF  ! end of year
           ELSE
              casaflux%stemnpp = 0.
           ENDIF ! CALL_POP


        ENDDO   ! end of idoy
     ENDDO   ! end of nyear


!$  if(nloop>=nloop1) &
!$       call totcnppools(2+nloop-nloop1,veg,casamet,casapool,bmcplant,bmnplant,bmpplant,bmclitter,bmnlitter,bmplitter, &
!$       bmcsoil,bmnsoil,bmpsoil,bmnsoilmin,bmpsoillab,bmpsoilsorb,bmpsoilocc,bmarea)

  ENDDO     ! end of nloop

  CALL casa_fluxout(CABLE_USER%CASA_SPIN_STARTYEAR + myearspin - 1 , veg, soil, casabal, casamet)

  !STOP

  ! write the last five loop pool size by PFT type
  OPEN(92,file='cnpspinlast5.txt')
  WRITE(92,921)
921 FORMAT('PFT total area in 10**12 m2', f12.4)
  DO nvt=1,mvtype
     WRITE(92,*) bmarea(nvt)
  ENDDO

  DO nvt=1,mvtype
     IF(bmarea(nvt) >0.0) THEN
        DO kloop=1,5
           WRITE(92,922) nvt, bmcplant(kloop,nvt,:),bmclitter(kloop,nvt,:),bmcsoil(kloop,nvt,:)
        ENDDO
        IF (icycle >1) THEN
           DO kloop=1,5
              WRITE(92,922) nvt, bmnplant(kloop,nvt,:),bmnlitter(kloop,nvt,:),bmnsoil(kloop,nvt,:), bmnsoilmin(kloop,nvt)
           ENDDO
        ENDIF

        IF(icycle >2) THEN
           DO kloop=1,5
              WRITE(92,922) nvt, bmpplant(kloop,nvt,:),bmplitter(kloop,nvt,:),bmpsoil(kloop,nvt,:),  &
                   bmpsoillab(kloop,nvt), bmpsoilsorb(kloop,nvt), bmpsoilocc(kloop,nvt)
           ENDDO
        ENDIF
     ENDIF
  ENDDO
922 FORMAT(i4,20(f10.4,2x))
  CLOSE(92)


151 FORMAT(i6,100(f12.5,2x))
END SUBROUTINE spincasacnp
