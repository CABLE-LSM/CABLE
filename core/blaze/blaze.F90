MODULE BLAZE_MOD

TYPE TYPE_BLAZE
   INTEGER,  DIMENSION(:),  ALLOCATABLE :: DSLR,ilon, jlat, Flix
   REAL,     DIMENSION(:),  ALLOCATABLE :: RAINF, KBDI, LR, U10,RH,TMAX,TMIN,AREA
   REAL,     DIMENSION(:),  ALLOCATABLE :: FFDI,FLI,ROS,Z,D,w, LAT, LON,DFLI,AB,CAvgAnnRainf
   REAL,     DIMENSION(:),  ALLOCATABLE :: DEADWOOD,POP_TO, POP_CWD, POP_STR,shootfrac
   REAL,     DIMENSION(:,:),ALLOCATABLE :: AnnRAINF, ABM, TO, AGC_g, AGC_w
   REAL,     DIMENSION(:,:),ALLOCATABLE :: AGLit_w, AGLit_g, BGLit_w, BGLit_g
   REAL,     DIMENSION(:,:),ALLOCATABLE :: AvgAnnRAINF, FLUXES
   CHARACTER,DIMENSION(:),  ALLOCATABLE :: FTYPE*6
   INTEGER                              :: T_AVG, YEAR, MONTH, DAY, DOY, NCELLS
   INTEGER                              :: BURNMODE ! 0=off, 1=BLAZE only, 2=BLAZE with POP
!CRM   INTEGER                              :: IGNITION ! 0=GFED3, 1=SIMFIRE
   REAL                                 :: FT,tstp
   LOGICAL                              :: USE_POP = .FALSE., ERR=.FALSE.
   CHARACTER                            :: GFEDP*80, FSTEP*7, BURNT_AREA_SRC*10
   CHARACTER                            :: OUTMODE = "std" ! "full" for diagnostical purposes
END TYPE TYPE_BLAZE

TYPE TYPE_TURNOVER
   REAL :: TO_ATM, TO_CWD, TO_STR
END TYPE TYPE_TURNOVER

REAL, DIMENSION(:,:), ALLOCATABLE  :: BLAZEFLX ! To BLAZE!!!

INTEGER, PARAMETER :: NTO     = 7 ! Number of TurnOver Parameters ,i.e. #lines below
INTEGER, PARAMETER :: LEAF    = 1
INTEGER, PARAMETER :: WOOD    = 2
INTEGER, PARAMETER :: FROOT   = 3
INTEGER, PARAMETER :: DEAD    = 4
INTEGER, PARAMETER :: MLIT    = 5
INTEGER, PARAMETER :: SLIT    = 6
INTEGER, PARAMETER :: CLIT    = 7

INTEGER, PARAMETER :: METB    = 1
INTEGER, PARAMETER :: STR     = 2
INTEGER, PARAMETER :: CWD     = 3

INTEGER, PARAMETER :: r_2     = 4
INTEGER, PARAMETER :: T_AVG   = 5

REAL,    PARAMETER :: fbranch = 0.05 ! fraction of CPLANT_w(wood) attributed to branches
REAL,    PARAMETER :: fbark   = 0.01 ! fraction of CPLANT_w(wood) attributed to bark

! Controlflag to determine whether full blaze or blaze/POP turn-overs used
INTEGER, PARAMETER :: DO_BLAZE_TO = 1
INTEGER, PARAMETER :: DO_POP_TO   = -1
INTEGER, PARAMETER :: NFLUX       = 12

!CLN !WRONG VARNAME! REAL, DIMENSION(5),PARAMETER :: SF = (/ 50.,50.,30.,20.,15./) !Euc*2, pinus123

! Turn Over Factors Suravski et al. 2012
REAL, DIMENSION(5,12),PARAMETER   ::   &
     TOF = (/ .0 , .0 , .05, .2 , .2 , &  !  1 Stems       -> ATM
              .0 , .0 , .15, .2 , .2 , &  !  2 Branches    -> ATM
              .03, .13, .25, .5 , .5 , &  !  3 Bark        -> ATM
              .02, .05, .1 , .6 , .6 , &  !  4 Leaves      -> ATM
              .0 , .0 , .05, .2 , .8 , &  !  5 Stems       -> Litter (CWD) !corrected*
              .0 , .02, .07, .2 , .8 , &  !  6 Branches    -> Litter (CWD) !corrected*
              .03, .13, .25, .5 , .5 , &  !  7 Bark        -> Litter (str)
              .05, .1 , .15, .3 , .4 , &  !  8 Leaves      -> Litter (str)
              .0 , .02, .02, .04, .04, &  !  9 FDEAD roots -> ATM
              .5 , .75, .75, .8 , .8 , &  ! 10 Deadwood    -> ATM
              .6 , .65, .85, 1. , 1. , &  ! 11 Bark Litter -> ATM
              .6 , .65, .85, 1. , 1. /)   ! 12 Leaf Litter -> ATM
!CLN DEADWOOD!!!
!CLN              .7 , .75, .8, .8 , .8 , &  ! 10 Deadwood    -> ATM
!CLN              .9 , .95, .95, 1. , 1. , &  ! 11 Bark Litter -> ATM
!CLN              .9 , .95, .95, 1. , 1. /)   ! 12 Leaf Litter -> ATM

REAL, PARAMETER :: MIN_FUEL = 120. ! Min fuel to spark a fire [g(C)/m2]

CONTAINS

SUBROUTINE INI_BLAZE (POPFLAG, BURNT_AREA_SOURCE, TSTEP, np, LAT, LON, BLAZE)

  !! Called from cable_input now

  IMPLICIT NONE
  
  LOGICAL            , INTENT(IN)    :: POPFLAG
  CHARACTER(LEN=*)   , INTENT(IN)    :: BURNT_AREA_SOURCE, TSTEP
  INTEGER            , INTENT(IN)    :: np
  REAL, DIMENSION(np), INTENT(IN)    :: LAT, LON
  TYPE(TYPE_BLAZE)   , INTENT(INOUT) :: BLAZE

  ! READ ini-nml
  BLAZE%NCELLS = np
  ALLOCATE ( BLAZE%DSLR    ( np ) )
  ALLOCATE ( BLAZE%RAINF   ( np ) )
  ALLOCATE ( BLAZE%LR      ( np ) )
  ALLOCATE ( BLAZE%KBDI    ( np ) )
  ALLOCATE ( BLAZE%FTYPE   ( np ) )
  BLAZE%FTYPE(:) = "Seeder"
  ALLOCATE ( BLAZE%AB      ( np ) )
  ALLOCATE ( BLAZE%AREA    ( np ) )
  ALLOCATE ( BLAZE%U10     ( np ) )
  ALLOCATE ( BLAZE%RH      ( np ) )
  ALLOCATE ( BLAZE%LAT     ( np ) )
  ALLOCATE ( BLAZE%LON     ( np ) )
  ALLOCATE ( BLAZE%JLAT    ( np ) )
  ALLOCATE ( BLAZE%ILON    ( np ) )
  ALLOCATE ( BLAZE%TMAX    ( np ) )
  ALLOCATE ( BLAZE%TMIN    ( np ) )
!  ALLOCATE ( BLAZE%F       ( np ) )
  ALLOCATE ( BLAZE%FLI     ( np ) )
  ALLOCATE ( BLAZE%FFDI     ( np ) )
  ALLOCATE ( BLAZE%FLIx    ( np ) )
  ALLOCATE ( BLAZE%DFLI    ( np ) )
  ALLOCATE ( BLAZE%ROS     ( np ) )
  ALLOCATE ( BLAZE%Z       ( np ) )
  ALLOCATE ( BLAZE%D       ( np ) )
  ALLOCATE ( BLAZE%w       ( np ) )
  ALLOCATE ( BLAZE%TO      ( np, NTO ) )
  ALLOCATE ( BLAZE%AnnRainf( np, 366 ) )
  ALLOCATE ( BLAZE%CAvgAnnRainf(np))
  ALLOCATE ( BLAZE%DEADWOOD( np ) )
  ALLOCATE ( BLAZE%SHOOTFRAC( np ) ) 
  ! POP related vars
  ALLOCATE ( BLAZE%POP_TO  ( np ) )
  ALLOCATE ( BLAZE%POP_CWD ( np ) )
  ALLOCATE ( BLAZE%POP_STR ( np ) )
  ALLOCATE ( BLAZE%FLUXES  ( np, 13 ) )
  ! SETTINGS FOR BLAZE (BLAZEFLAG)
  ! bit value:               0            | 1
  ! 0th bit(1), general    : off          | on
  ! 1st bit(2): BLAZE-mode : FLI-only     | Full-BLAZE
  ! 2nd bit(4): Ignition   : other (GFED) | SIMFIRE 
  ! more possible...for other ignition sources
  ! examples: 
  ! BLAZE(+1) with FLI-only(+0) and SIMFIRE(+4) -> BLAZEFLAG=5
  ! BLAZE(+1) full (+2) with GFED (+0)          -> BLAZEFLAG=3
   
  IF ( POPFLAG ) THEN
     BLAZE%BURNMODE = 2    ! POP-MODE, FLI generated, Fluxes according to POP
     WRITE(*,*) " BLAZE: POP-MODE (FLI and POP related fluxes)"
     PRINT*,"CLN request tstep = 1yr here!"
  ELSE
     BLAZE%BURNMODE = 1    ! Full. All Fluxes copmuted by BLAZE
     WRITE(*,*) " BLAZE: Full Mode"
  END IF
!CRM  IF ( TRIM(BURNT_AREA_SOURCE) == "SIMFIRE" ) THEN
!CRM     BLAZE%IGNITION = 1    
!CRM  ELSE IF ( TRIM(BURNT_AREA_SOURCE) == "PRESCRIBED" ) THEN
!CRM     BLAZE%IGNITION = 2    
!CRM  ELSE IF ( TRIM(BURNT_AREA_SOURCE) == "GFED3.1" ) THEN
!CRM     BLAZE%IGNITION = 3    
!CRM  ELSE 
!CRM     WRITE(*,*) " ERROR: Invalid cable_user%BURNT_AREA : ",TRIM(BURNT_AREA_SOURCE)
!CRM!!!CLN hier MPI-tauglichen abbruch!!!
!CRM     BLAZE%ERR = .TRUE.
!CRM     RETURN
!CRM  END IF
  WRITE(*,*) " Burnt-area source: ", TRIM(BURNT_AREA_SOURCE)

  BLAZE%DSLR      = 0
  BLAZE%RAINF     = 0.
  BLAZE%LR        = 0.
  BLAZE%KBDI      = 0.
  BLAZE%AB        = 0.
  BLAZE%DEADWOOD  = 0. 
  BLAZE%FSTEP     = TRIM(TSTEP)

  BLAZE%LAT       = LAT
  BLAZE%LON       = LON

  ! time for averaging annual rainfall [a]
  BLAZE%T_AVG = 5            
  ! deadwood decay scale time [a]  
  BLAZE%FT    = 30               
  ! Read initial Annual Rainfall Data
  ! CLN check for validity and more options...
  ! vh ! needs a BLAZE restart file?
  
END SUBROUTINE INI_BLAZE

SUBROUTINE BLAZE_ACCOUNTING(BLAZE, met, ktau, dels, year, doy)

  USE CABLE_DEF_TYPES_MOD, ONLY: MET_TYPE
  USE CABLE_COMMON_MODULE, ONLY: IS_LEAPYEAR
  USE cable_IO_vars_module, ONLY:  landpt

  IMPLICIT NONE 

  TYPE (TYPE_BLAZE)   :: BLAZE
  TYPE (MET_TYPE)     :: met
  INTEGER, INTENT(IN) :: ktau, year, doy
  REAL,    INTENT(IN) :: dels
  LOGICAL, SAVE       :: CALL1 = .TRUE.
  INTEGER             :: mp, x, i
  REAL                :: t_fac
  INTEGER,PARAMETER   :: sod = 86400
  LOGICAL             :: is_new_day, is_end_of_day
  REAL                :: v, rh, t, prec, avgAnnRf, dkbdi, FFDI

  mp = BLAZE%NCELLS
  
  is_new_day    = (MOD((ktau-1) * NINT(dels), sod) == 0 )
  is_end_of_day = (MOD( ktau    * NINT(dels), sod) == 0 )
  
  ! factor to get to daily data
  t_fac = 86400. / dels
 
  ! 1.Jan
  ! Update avg ann rainfall 

  ! Add last years 
  x = MOD(year, T_AVG) + 1
  IF ( is_new_day) THEN
     IF ( doy == 1 ) THEN
        BLAZE%AvgAnnRainf(:,x) = SUM(BLAZE%AnnRAINF(:,:),dim=2) 
        BLAZE%CAvgAnnRainf(:)  = SUM(BLAZE%AvgAnnRainf(:,:),dim=2) /REAL(T_AVG)
        BLAZE%AnnRAINF(:,:) = 0.
    END IF
     BLAZE%U10 (:) = 0.
     BLAZE%RH  (:) = 0.
     BLAZE%TMAX(:) = -999.
     BLAZE%TMIN(:) =  999.
  END IF

  DO i = 1, mp
     BLAZE%AnnRAINF(i,doy) = BLAZE%AnnRAINF(i,doy) + met%precip(landpt(i)%cstart)    
     BLAZE%U10 (i) = MAX(met%u10(landpt(i)%cstart),BLAZE%U10(i))        ! m/s -> km/s
     BLAZE%RH  (i) = BLAZE%RH(i) + met%rhum(landpt(i)%cstart)/t_fac     ! daily average rel. humidity
     BLAZE%TMAX(i) = MAX(met%tvair(landpt(i)%cstart),BLAZE%TMAX(i))
     BLAZE%TMIN(i) = MIN(met%tvair(landpt(i)%cstart),BLAZE%TMIN(i))     
  END DO

  if ( .NOT. is_leapyear(year) .AND. doy .EQ. 365 ) BLAZE%AnnRAINF(:,366) = 0.

  ! End of the day prepare daily data 
  IF ( is_end_of_day ) THEN
  
     do i = 1, mp
        v        = BLAZE%U10(i) * 3.6 ! m/s->km/h
        ! Gust parameterisation following ???
        v        = (214.7* (v+10)**(-1.6968) + 1 ) * v  
        rh       = BLAZE%RH(i)
        t        = BLAZE%TMAX(i)
        prec     = BLAZE%AnnRAINF(i,doy)
        avgAnnRF = BLAZE%CAvgAnnRainf(i)
        
        ! Days-since-last-rain and Last-Rainfall
        IF ( prec > 0.01 ) THEN
           IF ( BLAZE%DSLR(i) > 0 ) THEN
              BLAZE%LR(i) = prec
           ELSE
              BLAZE%LR(i) = BLAZE%LR(i) + prec
           END IF
           BLAZE%DSLR(i) = 0
        ELSE
           BLAZE%DSLR(i) = BLAZE%DSLR(i) + 1
        END IF

        ! Keetch-Byram Drought Index
        IF ( BLAZE%DSLR(i) == 0 ) THEN
           IF ( BLAZE%LR(i) > 5. ) THEN
              dkbdi = 5. - BLAZE%LR(i)
           ELSE
              dkbdi = 0.
           END IF

        ELSE
           dkbdi = (( 800. - BLAZE%KBDI(i) ) * (.968 * EXP(.0486 * (t * 9./5. + 32.)) &
                - 8.3) / 1000. / (1. + 10.88 * EXP(-.0441 * avgAnnRF/25.4)) * .254)
        END IF
        BLAZE%KBDI(i) = MAX(0., BLAZE%KBDI(i) + dkbdi)

        ! MacArthur Drought-Factor D
        BLAZE%D(i) = .191 * (BLAZE%KBDI(i) + 104.) * (BLAZE%DSLR(i) + 1.)**1.5 / &
             ( 3.52 * (BLAZE%DSLR(i) + 1.)**1.5 + BLAZE%LR(i) - 1. )
        BLAZE%D(i) = MAX(0.,MIN(10.,BLAZE%D(i)))

        ! MacArthur FFDI
        FFDI = 2. * EXP( -.45 + .987 * LOG(BLAZE%D(i)+.001) &
             - .03456 * RH + .0338 * T + .0234 * V )

        BLAZE%FFDI(i) = MAX(BLAZE%FFDI(i),FFDI)

     END do
  
  END IF
  
END SUBROUTINE BLAZE_ACCOUNTING


FUNCTION AVAIL_FUEL(FLIx, CPLANT_w, CPLANT_g, AGL_w, AGL_g)

  IMPLICIT NONE

  INTEGER, INTENT(IN)               :: FLIx
  REAL, DIMENSION(3), INTENT(IN)    :: CPLANT_w, CPLANT_g, AGL_w, AGL_g
  REAL                              :: AVAIL_FUEL

  AVAIL_FUEL = TOF(FLIx, 3) * fbark  * CPLANT_w(WOOD)             + &
               TOF(FLIx,10)          * AGL_w(CWD)                 + &
               TOF(FLIx,12)          * (AGL_w(METB) + AGL_w(STR)) + &
               TOF(FLIx, 9)          * CPLANT_w(FROOT)            + &
               CPLANT_g(LEAF) + AGL_g(METB) + AGL_g(STR)         
               
  ! CLN maybe take fbark out again

END FUNCTION AVAIL_FUEL

SUBROUTINE BLAZE_TURNOVER(AB, CPLANT_g, CPLANT_w, AGL_g, AGL_w, &
     BGL_g, BGL_w, shootfrac, TO, BT, BURNMODE, POP_TO)
!CRM SUBROUTINE BLAZE_TURNOVER(AB, CPLANT_g, CPLANT_w, AGL_g, AGL_w, &
!CRM     BGL_g, BGL_w, shootfrac, TO, POPFLAG, BT, POP_TO)
 
  IMPLICIT NONE

  !get index of relativ life-biomass loss in TOF table plus fract between 2bins 
  !-> apply idx and factor on other TOs 

  TYPE(TYPE_TURNOVER),INTENT(INOUT) :: TO(7)
  INTEGER,            INTENT(IN)    :: BURNMODE
  REAL,               INTENT(IN)    :: AB, shootfrac
  REAL,               INTENT(INOUT) :: CPLANT_w(3) , CPLANT_g(3)
  REAL, DIMENSION(3), INTENT(INOUT) :: AGL_w, AGL_g, BGL_w, BGL_g
  REAL,               INTENT(OUT)   :: BT(13)
  REAL,     OPTIONAL, INTENT(IN)    :: POP_TO
  TYPE(TYPE_TURNOVER)               :: MTO(7)
  REAL    :: rd, BLAZE_EQUIV, fAB, totw

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! TOF Tuning:

!CRM  IF ( POPFLAG .EQ. -1 ) THEN ! woody TO done by POP if active
  IF ( BURNMODE .EQ. 2 ) THEN ! woody TO done by POP if active
     IF ( .NOT. PRESENT( POP_TO ) ) STOP  "Please provide POP TO!"

     ! Leave TO is adjusted to fraction of total leaf loss = fraction of total wood loos
     ! and then distributed to Atm ans metb litter by relative distribution as computed
     ! by BLAZE

     IF ( POP_TO .GT. 0. ) THEN

        rd   = POP_TO / CPLANT_w(WOOD) 
        totw = TO(WOOD)%TO_ATM + TO(WOOD)%TO_CWD + TO(WOOD)%TO_STR
        IF ( totw .GT. 0. ) THEN      
           TO(WOOD)%TO_ATM = MIN(rd / totw * TO(WOOD)%TO_ATM , 0.999)
           TO(WOOD)%TO_STR = MIN(rd / totw * TO(WOOD)%TO_STR , 0.999)
           TO(WOOD)%TO_CWD = MIN(rd / totw * TO(WOOD)%TO_CWD , 0.999)
        ELSE
           TO(WOOD)%TO_ATM = rd * 0.48  
           TO(WOOD)%TO_STR = rd * 0.04  
           TO(WOOD)%TO_CWD = rd * 0.48  
        ENDIF
        ! fine roots removed prportional to a.g. woody biomass
        TO(FROOT)%TO_STR = MAX( 0., rd - TO(FROOT)%TO_ATM )

        totw = TO(LEAF)%TO_ATM + TO(LEAF)%TO_STR
        IF ( totw .GT. 0. ) THEN
           TO(LEAF)%TO_ATM = MIN(rd / totw * TO(LEAF)%TO_ATM , 0.999)
           TO(LEAF)%TO_STR = MIN(rd / totw * TO(LEAF)%TO_STR , 0.999) 
        ELSE
           TO(LEAF)%TO_ATM = rd * 1. / 3.
           TO(LEAF)%TO_STR = rd * 2. / 3. 
        ENDIF
           
        fAB = 1.
     ELSE
        fAB = 0.
     ENDIF
  ELSE     
     fAB = AB
  ENDIF

  ! Mass Fluxes
  MTO(LEAF )%TO_ATM = fAB * TO(LEAF )%TO_ATM * CPLANT_w (LEAF )
  MTO(WOOD )%TO_ATM = fAB * TO(WOOD )%TO_ATM * CPLANT_w (WOOD )
  MTO(LEAF )%TO_STR = fAB * TO(LEAF )%TO_STR * CPLANT_w (LEAF )
  !CVH should the line below use only the above-ground wood pool?
  MTO(WOOD )%TO_STR = fAB * TO(WOOD )%TO_STR * CPLANT_w (WOOD )
  MTO(WOOD )%TO_CWD = fAB * TO(WOOD )%TO_CWD * CPLANT_w (WOOD )
  MTO(FROOT)%TO_ATM = fAB * TO(FROOT)%TO_ATM * CPLANT_w (FROOT)
  MTO(FROOT)%TO_STR = fAB * TO(FROOT)%TO_STR * CPLANT_w (FROOT)
  MTO(MLIT )%TO_ATM =  AB * TO(MLIT )%TO_ATM * AGL_w(METB ) 
  MTO(SLIT )%TO_ATM =  AB * TO(SLIT )%TO_ATM * AGL_w(STR  ) 
  MTO(CLIT )%TO_ATM =  AB * TO(CLIT )%TO_ATM * AGL_w(CWD  ) 

  ! Diagnostics
  BT( 1) = MTO(LEAF )%TO_ATM 
  BT( 2) = MTO(WOOD )%TO_ATM
  BT( 3) = MTO(LEAF )%TO_STR
  BT( 4) = MTO(WOOD )%TO_STR
  BT( 5) = MTO(WOOD )%TO_CWD
  BT( 6) = MTO(MLIT )%TO_ATM 
  BT( 7) = MTO(SLIT )%TO_ATM 
  BT( 8) = MTO(CLIT )%TO_ATM  
  BT( 9) = MTO(FROOT)%TO_ATM  
  BT(10) = MTO(FROOT)%TO_STR  
  BT(11) = AB * CPLANT_g (LEAF)
  BT(12) = AB * AGL_g(METB)
  BT(13) = AB * AGL_g(STR )

!CLN  PRINT*,'BT( 1) =  MTO(LEAF )%TO_ATM    ',MTO(LEAF)%TO_ATM    
!CLN  PRINT*,'BT( 2) =  MTO(WOOD )%TO_ATM    ',MTO(WOOD)%TO_ATM    
!CLN  PRINT*,'BT( 3) =  MTO(LEAF )%TO_STR    ',MTO(LEAF)%TO_STR    
!CLN  PRINT*,'BT( 4) =  MTO(WOOD )%TO_STR    ',MTO(WOOD)%TO_STR    
!CLN  PRINT*,'BT( 5) =  MTO(WOOD )%TO_CWD    ',MTO(WOOD)%TO_CWD    
!CLN  PRINT*,'BT( 6) =  MTO(MLIT )%TO_ATM    ',MTO(MLIT)%TO_ATM    
!CLN  PRINT*,'BT( 7) =  MTO(SLIT )%TO_ATM    ',MTO(SLIT)%TO_ATM    
!CLN  PRINT*,'BT( 8) =  MTO(CLIT )%TO_ATM    ',MTO(CLIT)%TO_ATM    
!CLN  PRINT*,'BT( 9) =  MTO(FROOT)%TO_ATM    ',MTO(FROOT)%TO_ATM
!CLN  PRINT*,'BT(10) =  MTO(FROOT)%TO_STR    ',MTO(FROOT)%TO_STR 
!CLN  PRINT*,'BT(11) =  AB * CPLANT_g (LEAF )',AB * CPLANT_g (LEAF)
!CLN  PRINT*,'BT(12) =  AB * CLITTER_g(METB )',AB * AGL_g(METB)
!CLN  PRINT*,'BT(13) =  AB * CLITTER_g(STR  )',AB * AGL_g(STR )

  ! Update Pools
!CVH what about BGL_g?
  ! Grass fuels burned by area
  CPLANT_g (LEAF)  = CPLANT_g (LEAF) * (1. - AB)
  AGL_g(METB)      = AGL_g(METB)     * (1. - AB)
  AGL_g(STR )      = AGL_g(STR )     * (1. - AB)
  
  ! Leaves
  CPLANT_w (LEAF)  = &
       CPLANT_w (LEAF) - (MTO(LEAF)%TO_ATM + MTO(LEAF)%TO_STR)
  ! Wood
  CPLANT_w (WOOD)  = &
       CPLANT_w (WOOD) - (MTO(WOOD)%TO_ATM + MTO(WOOD)%TO_STR + MTO(WOOD)%TO_CWD) 
  ! Fine roots
  CPLANT_w (FROOT) = &
       CPLANT_w (FROOT)- (MTO(FROOT)%TO_ATM+ MTO(FROOT)%TO_STR)

  ! Litter above ground
  AGL_w(METB)      = AGL_w(METB) - MTO(MLIT)%TO_ATM 
  AGL_w(STR )      = AGL_w(STR ) - MTO(SLIT)%TO_ATM + MTO(WOOD)%TO_STR + &
                     MTO(LEAF)%TO_STR 
  AGL_w(CWD )      = AGL_w(CWD ) - MTO(CLIT)%TO_ATM + MTO(WOOD)%TO_CWD * shootfrac

  ! Litter below ground
  BGL_w(STR ) = BGL_w(STR ) + MTO(FROOT)%TO_STR 
  BGL_w(CWD ) = BGL_w(CWD ) + MTO(WOOD)%TO_CWD * (1. - shootfrac)

  IF ( ANY ( AGL_w .LT. 0. )) THEN
     WRITE(*,*)"CLIITE_W < 0 ", AGL_w
     STOP -1
  ENDIF

END SUBROUTINE BLAZE_TURNOVER

FUNCTION BURNTIME( YEAR, DOY, FSTEP )

  USE CABLE_COMMON_MODULE, ONLY : DOYSOD2YMDHMS, IS_LEAPYEAR

  INTEGER,  INTENT(IN):: YEAR, DOY
  CHARACTER(LEN=*),INTENT(IN):: FSTEP
  LOGICAL             :: BURNTIME 
  INTEGER             :: MM, DD, EOM
  INTEGER,DIMENSION(12),PARAMETER :: LOM = (/31,28,31,30,31,30,31,31,30,31,30,31/)

  BURNTIME = .FALSE.
  
  CALL DOYSOD2YMDHMS( YEAR, DOY, 1, MM, DD )
  IF ( TRIM(FSTEP) .EQ. "annual" ) THEN
     IF ( (.NOT. (IS_LEAPYEAR(YEAR)) .AND. DOY .EQ. 365 ) .OR. DOY .EQ. 366 ) & 
          BURNTIME = .TRUE.
  ELSE IF ( TRIM(FSTEP) .EQ. "monthly" ) THEN
     IF ( IS_LEAPYEAR(YEAR) .AND. MM.EQ.2 ) THEN
        EOM = 29
     ELSE
        EOM = LOM(MM)
     ENDIF
     IF ( DD .EQ. EOM ) BURNTIME = .TRUE.

  ELSE IF ( TRIM(FSTEP) .EQ. "daily" ) THEN
     BURNTIME = .TRUE.
  ELSE IF ( TRIM(FSTEP) .EQ. "none" ) THEN
     WRITE(6,FMT='(A36,I2.2,x,I2.2,x,I4)') &
          "No fire info available from GFED on ",DD, MM, YEAR
  ELSE 
     STOP "Wrong fire time-step!"
  END IF

END FUNCTION BURNTIME

SUBROUTINE COMBUST (BLAZE, np, CPLANT_g, CPLANT_w, TO, BURN )

  IMPLICIT NONE

  TYPE(TYPE_BLAZE),   INTENT(INOUT) :: BLAZE
  INTEGER,            INTENT(IN)    :: NP
  REAL,               INTENT(IN)    :: CPLANT_w(3), CPLANT_g(3)
!CRM  REAL,               INTENT(IN)    :: AGL_w(3),AGL_g(3), DEADWOOD
  LOGICAL,            INTENT(IN)    :: BURN
  TYPE(TYPE_TURNOVER),INTENT(INOUT) :: TO(NTO)

  REAL, PARAMETER       :: H = 20. ! Heat Yield, [MJ/kg]
  REAL      :: FT, LR, dKBDI
  REAL      :: F,w, D, ROS, RH, T, V, Z
  REAL      :: AvAnnRain, FLI, TOT_g
  INTEGER   :: DSLR, FLIx, i

! Fuel Moisture driving Heat yield or Available Fuel...?
! CLN Dry matter? 

!  IF ( .NOT. ANY (TRIM(BLAZE%FSTEP) .EQ. "none" ) ) &
!     STOP "Fire time-step not specified. Abort combust"
    

  ft  = BLAZE%FT
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Get FLI (Byram's Law)  = H * w * r
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Get Rate-Of-Spread following McArthur fire-danger meter Mark 5 Version (Noble et al. 1980) 

  !  V  from Wind [m/s] 
  !  Rule of thumb: U2 in Forest ~ 1/4*U10 to 1/11*U10 [height of canopy small->large]
  !  T  from Tair [deg C]
  !  RH from Qair [%]
  !  DSLR days since rainfall [d]
  !  P      precipitation (24 hours) [mm]
  !  from SIMFIRE or Rainf 

  ! Keetch-Byram Drought Index [ ] 
  ! Keetch, j.j. and G.M. Byram, 1968, A drought index for forest fire control, 
  ! US Dept, Agr. Forest Service Res. Paper SE-38
  IF ( BLAZE%RAINF(np) .GT. 0. ) THEN
     IF ( BLAZE%DSLR(np) .GT. 0 ) THEN
        BLAZE%LR(np) = BLAZE%RAINF(np) 
     ELSE
        BLAZE%LR(np) = BLAZE%LR(np) + BLAZE%RAINF(np)
     ENDIF
     BLAZE%DSLR(np) = 0
  ELSE
     BLAZE%DSLR(np) = BLAZE%DSLR(np) + 1
  ENDIF
  
  DSLR      = BLAZE%DSLR(np)
  AvAnnRain = SUM(BLAZE%AnnRAINF(np,:)) 
  V         = BLAZE%U10(np)
  RH        = BLAZE%RH(np)
  T         = BLAZE%TMAX(np)

  IF ( DSLR .EQ. 0 ) THEN
     IF ( BLAZE%LR(np) .GE. 5 ) THEN
        dKBDI = 5. - BLAZE%LR(np)
     ELSE
        dKBDI = 0.
     ENDIF
  ELSE
     dKBDI = ((800. - BLAZE%KBDI(np)) * (.968 * EXP(.0486 * (T * 9./5. &
          + 32.)) - 8.3) / 1000. / (1. + 10.88 * EXP(-.0441 * AvAnnRain/25.4)) * .254)
  ENDIF

  BLAZE%KBDI(np) = MAX(0.,BLAZE%KBDI(np) + dKBDI)

!  IF ( .NOT. BURN ) RETURN

  ! McArthur drought factor [ ]
  D   = .191 * ( BLAZE%KBDI(np) + 104. ) * ( REAL(DSLR) + 1. )**1.5 / &
       ( 3.52 * ( REAL(DSLR) + 1. )**1.5 + BLAZE%LR(np) - 1. )
  D   = MAX(0.,MIN(10.,D))

  ! McArthur Fire-Danger index [ ]
  !                             Forest  Grass
  ! Catastrophic (Code Red)	100 +   150 +
  ! Extreme	                75 – 99 100 – 149
  ! Severe	                50 – 74 50 – 99
  ! Very High	                25 – 49 25 – 49
  ! High	                12 – 24 12 – 24
  ! Low–Moderate            	 0 – 11  0 – 11

  F   = 2. * EXP( -.45 + .987 * LOG(D+.001) - .03456 * RH + .0338 * T + .0234 * V )
  F   = MAX(0.,F) 

  ! available Fuels (Litter, deadwood and under/midstorey[when is it distinguished?])) [kg/m^2]
  ! Assume totel Grass is burnt
  ! Factor 0.6 / 0.5  taken from fires < 750 kW/m of Suravski et al. 2012
  FLIx = 1

  w = AVAIL_FUEL(FLIx, CPLANT_w, CPLANT_g, BLAZE%AGLit_w(np,:), BLAZE%AGLit_g(np,:) )
  
  DO i = 1, 4 ! THE LADDER

     ! Rate of Spread [m/s]
     ! original ROS = 0.0012 * F * w 
     ROS = 3.3333E-05 * F * w 
  
     ! Flame Height  [m]
     ! original Z   = 13. * ROS + 0.24 * w - 2.
     Z   = 46.8 * ROS + 0.024 * w - 2.
     Z   = MAX(0.,Z)

     FLI = H * w * ROS

     IF ( w .LT. MIN_FUEL ) THEN
        FLIx = 0
        EXIT
     ENDIF
     ! The LADDER
     IF ( FLI .GT. 7000. ) THEN
        FLIx = 4
        IF ( TRIM(BLAZE%FTYPE(np)) .EQ. "Seeder" ) FLIx = 5
     ELSE IF ( FLI .GT. 3000. ) THEN
        FLIx = 3
     ELSE IF ( FLI .GT. 750. ) THEN
        FLIx = 2
     ENDIF
     ! End here if available fuel doesn't get you into next regime
     IF ( i .GE. FLIx ) EXIT
     w = AVAIL_FUEL(FLIx, CPLANT_w, CPLANT_g, BLAZE%AGLit_w(np,:), BLAZE%AGLit_g(np,:) )
  END DO

  BLAZE%FLIX(np) = flix
  BLAZE%FFDI(np) = F
  BLAZE%FLI(np)  = FLI
  BLAZE%ROS(np)  = ROS
  BLAZE%Z(np)    = Z
  BLAZE%D(np)    = D
  BLAZE%w(np)    = w

  ! END McArthur
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! ROS

  ! Speer 2005, modelling fire weather and fire spread rates for two bushfires near sydney
  ! Aust. Met. Mag. 50, 2001 [grass and shrubs only]
!CLN  ROS = 0.049 * U2^1.21 * VegHeight^0.54

  ! Cheney 2012 
! USE Cheney for Grassla...
!CLN  ROS = 
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Get combustion factors from FLI (Suravski 2012) AS IS
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! TURN-OVER RATES
  ! 

  ! Live to Atmosphere
  TO(LEAF)%TO_ATM    =                          TOF(FLIx,  4) 
  TO(WOOD)%TO_ATM    = (1. - fbranch - fbark) * TOF(FLIx,  1) + &
                       fbranch                * TOF(FLIx,  2) + &
                       fbark                  * TOF(FLIx,  3) 
  TO(FROOT)%TO_ATM   =                          TOF(FLIx,  9) 
  
  ! Live to Litter                         
  TO(LEAF)%TO_STR    =                          TOF(FLIx,  8) 
  TO(WOOD)%TO_STR    = fbark                  * TOF(FLIx,  7) 
  TO(WOOD)%TO_CWD    = (1. - fbranch - fbark) * TOF(FLIx,  5) + &
                       fbranch                * TOF(FLIx,  6) 
  TO(FROOT)%TO_STR   = MAX( TO(WOOD)%TO_ATM + TO(WOOD)%TO_STR + TO(WOOD)%TO_CWD - &
       TO(FROOT)%TO_ATM , 0. )

  ! Litter to Atmosphere                
  TO(MLIT)%TO_ATM    = TOF(FLIx, 11) 
  TO(SLIT)%TO_ATM    = TOF(FLIx, 11) 
  TO(CLIT)%TO_ATM    = TOF(FLIx, 10) 

END SUBROUTINE COMBUST


SUBROUTINE RUN_BLAZE(BLAZE, SF, CPLANT_g, CPLANT_w, tstp, YYYY, doy, TO ) !, CTRL)
!CRMSUBROUTINE RUN_BLAZE(ncells, LAT, LON, shootfrac,CPLANT_g, CPLANT_w, AGL_g, AGL_w, &
!CRM     BGL_g, BGL_w, RAINF, TMIN, TMAX, RH, U10,AvgAnnMaxFAPAR, modis_igbp, &
!CRM     AvgAnnRainf, AB, FLI, DFLI, FFDI, TO, tstp, YYYY, doy, POPFLAG,popd,mnest,BLAZE_FSTEP ) !, CTRL)

  USE CABLE_COMMON_MODULE, ONLY: IS_LEAPYEAR, DOYSOD2YMDHMS
  USE SIMFIRE_MOD 

  IMPLICIT NONE

!CRM  INTEGER, INTENT(IN) :: POPFLAG

  TYPE(TYPE_BLAZE)    :: BLAZE
  TYPE(TYPE_SIMFIRE)  :: SF
  TYPE(TYPE_TURNOVER) :: TO(BLAZE%NCELLS,7)

  TYPE(TYPE_TURNOVER),ALLOCATABLE,SAVE :: TOsav(:,:)
  REAL               ,ALLOCATABLE,SAVE :: FLIsav(:)

  INTEGER          :: np, doy, YYYY, CTRL, MM, DD, DOM(12)
  INTEGER          :: i
  REAL             :: CPLANT_g(BLAZE%NCELLS,3), CPLANT_w(BLAZE%NCELLS,3), tstp
  REAL, DIMENSION(BLAZE%NCELLS,3) :: AGL_g, AGL_w, BGL_g, BGL_w
  REAL, DIMENSION(BLAZE%NCELLS)   :: &
       RAINF,             & ! [mm/d]
       TMIN,              & ! [deg C]
       TMAX,              & ! [deg C]
       RH,                & ! [%]
       U10,               & ! [km/h]
!CRM       AvgAnnMaxFAPAR,    & ! [frac.]
!CRM       AvgAnnRainf,       & ! [mm/a]
       AB,                & ! [frac.]
       FLI,               & ! [kW/m]
       DFLI,              & ! [kW/m]
       FFDI,              &  
       popd,              & ! [kW/m]
       mnest,             &
       shootfrac
               
! CLN to simfire!!! 
!CRM  INTEGER,DIMENSION(BLAZE%NCELLS)  :: modis_igbp  ! [0,17] landcover index

  LOGICAL, SAVE :: BCALL1 = .TRUE.

!CRM  REAL          :: DEADWOOD(ncells)

!CRM  CHARACTER     :: BLAZE_FSTEP*7

  DOM = (/ 31,28,31,30,31,30,31,31,30,31,30,31 /)
  IF ( IS_LEAPYEAR(YYYY) ) DOM(2) = 29

  BLAZE%YEAR = YYYY
  BLAZE%DOY  = DOY
  CALL DOYSOD2YMDHMS( YYYY, doy,10, MM, DD )

  BLAZE%DAY   = DD
  BLAZE%MONTH = MM

!!CLN Go to init part

!!CLN  IF ( BCALL1 ) THEN 
!!CLN     ALLOCATE( TOsav(NCELLS,NTO), FLIsav(NCELLS) )
!!CLN     CALL INI_BLAZE( BLAZE, ncells, npatch)
!!CLN     BLAZE%tstp = tstp
!!CLN
!!CLN     BLAZE%LAT = LAT
!!CLN     BLAZE%LON = LON
!!CLN
!!CLN     ! SET STARTING VALS ON AnnRainf (avg. 1990 - 2012)
!!CLN     DO i = 1, ncells
!!CLN        BLAZE%AnnRAINF(i,1:365) = AvgAnnRainf(i)/ 365.25
!!CLN        BLAZE%AnnRAINF(i,366)   = AvgAnnRainf(i)/(365.25*4.)
!!CLN     END DO
!!CLN
!!CLN!CLN     IF ( BLAZE%IGNITION .EQ. 1 ) THEN
!!CLN     CALL INI_SIMFIRE( NCELLS, SF, modis_igbp, LAT, LON )
!!CLN     SF%FAPAR = AvgAnnMaxFAPAR
!!CLN     IF ( BLAZE%IGNITION .EQ. 1 ) BLAZE%FSTEP = "annual"
!!CLN
!!CLN     BCALL1 = .FALSE.
!!CLN  END IF

  IF ( DD .EQ. 1 ) BLAZE%FLI(:) = 0.

  ! READ GFED BA DATA
  AB(:) = 0.
  IF ( TRIM(BLAZE%BURNT_AREA_SRC) .EQ. "GFED3.1" ) THEN 
!     CALL GET_GFED( BLAZE )
!     CALL GET_GFED4_BA( BLAZE )
!     CALL GET_GFED41s_BA( BLAZE )
     WRITE(*,*)'GFED4 BA not available. Set cable_user%BURNT_AREA == "SIMFIRE"'
     STOP(-1)
     IF ( TRIM(BLAZE%FSTEP) .EQ. "none" ) THEN
        CALL SIMFIRE ( SF, RAINF, TMAX, TMIN, DOY, YYYY, BLAZE%AB )
        popd = SF%POPD
        mnest= SF%MAX_NESTEROV
        BLAZE%FSTEP = "annual"
     ENDIF
  ELSEIF ( TRIM(BLAZE%BURNT_AREA_SRC) .EQ. "SIMFIRE" ) THEN
     ! CALL SIMFIRE DAILY FOR ACOUNTING OF PARAMETERS
     CALL SIMFIRE ( SF, RAINF, TMAX, TMIN, DOY, YYYY, BLAZE%AB )

     DO np = 1, BLAZE%NCELLS
        IF ( AVAIL_FUEL(1, CPLANT_w(np,:), CPLANT_g(np,:), AGL_w(np,:), AGL_g(np,:)) .LE. MIN_FUEL ) &
             BLAZE%AB(np) = 0.
     END DO
     popd = SF%POPD
     mnest= SF%MAX_NESTEROV
  ELSE
     WRITE(*,*) "Wrong ignition type chosen: ",BLAZE%BURNT_AREA_SRC
     STOP -1
  ENDIF
  
!CRM  DEADWOOD    = BLAZE%DEADWOOD
!CRM  BLAZE%Rainf = RAINF
!CRM  BLAZE%TMAX  = TMAX ! deg C
!CRM  BLAZE%TMIN  = TMIN ! deg C
!CRM  BLAZE%U10   = U10 
!CRM  BLAZE%RH    = RH 
 
  ! Apply half of former deadwood to atm now How to distribut (str
  ! set following Fraver 2013 pinus rosinosa (hardwood/decid. wood to be added
   
  DO np = 1, BLAZE%NCELLS

!CLN     AGL(np,CWD) = AGL(np,CWD) + &
!CLN          !       (1.-exp(-0.5*tstp/SF(ft))) * SUM(DEADWOOD,dim=2) 
!CLN          (1.-exp(-0.5*tstp/15.)) * SUM( DEADWOOD(np,:)) 
!CLN      DEADWOOD(np,:)  = DEADWOOD(np,:) * exp(-0.5*tstp/15.)
     
     IF ( BLAZE%FSTEP .NE. "daily" ) THEN
        IF ( DD.EQ.1 ) THEN
           FLIsav(np) = 0
        ELSE
           FLIsav(np) = BLAZE%FLI(np)
        ENDIF
     ENDIF

     CALL COMBUST( BLAZE, np, CPLANT_g(np,:), CPLANT_w(np,:),TO(np,:),BLAZE%AB(np).GT.0 )
!CRM     CALL COMBUST( BLAZE, np, CPLANT_g(np,:), CPLANT_w(np,:), AGL_g(np,:), &
!CRM          AGL_w(np,:),DEADWOOD(np), TO(np,:), BLAZE%AB(np).GT.0 )

     BLAZE%DFLI(np) = BLAZE%FLI(np)
     IF ( BLAZE%FSTEP .NE. "daily" ) THEN
        IF ( BLAZE%FLI(np).GT.FLIsav(np)) THEN
           TOsav(np,:)   = TO(np,:)
        ELSE IF ( BLAZE%FSTEP.EQ."monthly" .AND. DD.EQ.DOM(MM) ) THEN
           TO(np,:)      = TOsav(np,:)
        ELSE 
           BLAZE%FLI(np) = FLIsav(np)
        ENDIF
     ENDIF

     
     BLAZE%TO(np,:) = 0.
     IF ( BLAZE%BURNMODE .EQ. 1 .AND. BURNTIME(BLAZE%YEAR, BLAZE%DOY, BLAZE%FSTEP) &
          .AND. BLAZE%AB(np) .GT. 0. ) THEN
        !Apply to patches or only cell
        CALL BLAZE_TURNOVER( BLAZE%AB(np), CPLANT_g(np,:), CPLANT_w(np,:), &
             AGL_g(np,:), AGL_w(np,:),BGL_g(np,:), BGL_w(np,:), &
             BLAZE%shootfrac(np), TO(np,:), BLAZE%FLUXES(np,:), BLAZE%BURNMODE )
!CRM        CALL BLAZE_TURNOVER( BLAZE%AB(np), CPLANT_g(np,:), CPLANT_w(np,:), &
!CRM             AGL_g(np,:), AGL_w(np,:),BGL_g(np,:), BGL_w(np,:), &
!CRM             shootrfac(np), TO(np,:), POPFLAG, BLAZEFLX(np,:) )
     ENDIF

     ! FLI Accounting for POP
     if ( BLAZE%BURNMODE .EQ. 2 ) THEN
        IF ( DOY.EQ. 1 )FLI(np) = 0.
        FLI(np) = MAX(BLAZE%FLI(np),FLI(np))
     ENDIF

     DFLI(np) = BLAZE%DFLI(np)

     ! Rainfall accounting
     IF ( doy .EQ. 366 ) THEN
        BLAZE%AnnRAINF(np,doy) = ( REAL(BLAZE%T_AVG - 1) *  BLAZE%AnnRAINF(np,doy) + &
             Rainf(np)/4. ) / REAL(BLAZE%T_AVG)
     ELSE
        BLAZE%AnnRAINF(np,doy) = ( REAL(BLAZE%T_AVG - 1) *  BLAZE%AnnRAINF(np,doy) + &
             Rainf(np) ) / REAL(BLAZE%T_AVG)
     ENDIF
     BLAZE%RAINF(np) = Rainf(np)

     AB(np) = BLAZE%AB(np)

     ! Apply other half of former deadwood to litter now How to distribut (str
!CLN     AGL(np,CWD) = AGL(np,CWD) + &
!CLN          !       (1.-exp(-0.5*tstp/SF(ft))) * SUM(DEADWOOD(np,:)) 
!CLN          (1.-exp(-0.5*tstp/15.)) * SUM(DEADWOOD(np,:)) 

!CLN      BLAZE%DEADWOOD(np,:) = DEADWOOD(np,:) * exp(-0.5*tstp/15.)

  END DO

!CRM  FFDI        = BLAZE%FFDI
!CRM  BLAZE_FSTEP = BLAZE%FSTEP

!  WHERE ( BLAZE%AB .GT. 0. ) 
!  FLI = SF%FLI
!  END WHERE
!  IF ( BURNTIME(BLAZE%YEAR, BLAZE%DOY, BLAZE%FSTEP) )&
!  IF ( BLAZE%DAY .EQ. 1 ) &
!       CALL BLAZE_DIAG( NCELLS, BLAZE, CPLANT_g, CPLANT_w, AGL_g, AGL_w, TO, "dref")

END SUBROUTINE RUN_BLAZE

END MODULE BLAZE_MOD

