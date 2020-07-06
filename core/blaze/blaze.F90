MODULE BLAZE_MOD

TYPE TYPE_BLAZE
   INTEGER,  DIMENSION(:),  ALLOCATABLE :: DSLR,ilon, jlat, Flix
   REAL,     DIMENSION(:),  ALLOCATABLE :: RAINF, KBDI, LR, U10,RH,TMAX,TMIN,AREA, w_prior, FDI
   REAL,     DIMENSION(:),  ALLOCATABLE :: FFDI,FLI,ROS,Z,D,w, LAT, LON,DFLI,AB,CAvgAnnRainf
   REAL,     DIMENSION(:),  ALLOCATABLE :: DEADWOOD,POP_TO, POP_CWD, POP_STR,shootfrac
   REAL,     DIMENSION(:,:),ALLOCATABLE :: AnnRAINF, ABM, TO, AGC_g, AGC_w
   REAL,     DIMENSION(:,:),ALLOCATABLE :: AGLit_w, AGLit_g, BGLit_w, BGLit_g
   REAL,     DIMENSION(:,:),ALLOCATABLE :: CPLANT_g, CPLANT_w
   REAL,     DIMENSION(:,:),ALLOCATABLE :: AvgAnnRAINF, FLUXES
   CHARACTER,DIMENSION(:),  ALLOCATABLE :: FTYPE*6
   INTEGER                              :: T_AVG, YEAR, MONTH, DAY, DOY, NCELLS, time
   INTEGER                              :: BURNMODE ! 0=off, 1=BLAZE only, 2=BLAZE with POP
   !CRM INTEGER                              :: IGNITION ! 0=GFED3, 1=SIMFIRE
   REAL                                 :: FT,tstp
   LOGICAL                              :: USE_POP = .FALSE., ERR=.FALSE.
   CHARACTER                            :: GFEDP*80, FSTEP*7, BURNT_AREA_SRC*10
   CHARACTER(LEN=4)                     :: OUTMODE = "full" !"std" ! "full" for diagnostical purposes
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

! Turn Over Factors Surawski et al. 2012
REAL, DIMENSION(5,12),PARAMETER   ::   &
     TOF = reshape( &
           (/ .0 , .0 , .05, .2 , .2 , &   !  1 Stems       -> ATM
              .0 , .0 , .15, .2 , .2 , &   !  2 Branches    -> ATM
              .03, .13, .25, .5 , .5 , &   !  3 Bark        -> ATM
              .02, .05, .1 , .6 , .6 , &   !  4 Leaves      -> ATM
              .0 , .0 , .05, .2 , .8 , &   !  5 Stems       -> Litter (CWD) !corrected*
              .0 , .02, .07, .2 , .8 , &   !  6 Branches    -> Litter (CWD) !corrected*
              .03, .13, .25, .5 , .5 , &   !  7 Bark        -> Litter (str)
              .05, .1 , .15, .3 , .4 , &   !  8 Leaves      -> Litter (str)
              .0 , .02, .02, .04, .04, &   !  9 FDEAD roots -> ATM
              .5 , .75, .75, .8 , .8 , &   ! 10 Deadwood    -> ATM
              .6 , .65, .85, 1. , 1. , &   ! 11 Bark Litter -> ATM
              .6 , .65, .85, 1. , 1. /), & ! 12 Leaf Litter -> ATM
              (/5,12/) )
!CLN DEADWOOD!!!
!CLN              .7 , .75, .8, .8 , .8 , &  ! 10 Deadwood    -> ATM
!CLN              .9 , .95, .95, 1. , 1. , &  ! 11 Bark Litter -> ATM
!CLN              .9 , .95, .95, 1. , 1. /)   ! 12 Leaf Litter -> ATM

REAL, PARAMETER :: MIN_FUEL = 120. ! Min fuel to spark a fire [g(C)/m2]


CONTAINS

SUBROUTINE INI_BLAZE ( np, LAT, LON, BLAZE)

  !! Called from cable_driver now
  USE cable_common_module, ONLY: get_unit

  IMPLICIT NONE

  INTEGER            , INTENT(IN)    :: np
  REAL, DIMENSION(np), INTENT(IN)    :: LAT, LON
  TYPE(TYPE_BLAZE)   , INTENT(INOUT) :: BLAZE
  INTEGER, PARAMETER        :: NPOOLS = 3
  CHARACTER(len=400)   :: HydePath, BurnedAreaFile, &
       BurnedAreaClimatologyFile, SIMFIRE_REGION
  CHARACTER(len=10)   :: BurnedAreaSource
  INTEGER :: iu

  NAMELIST /BLAZENML/ HydePath,  BurnedAreaSource, BurnedAreaFile, BurnedAreaClimatologyFile, &
       SIMFIRE_REGION

  ! READ BLAZE settings
  CALL GET_UNIT(iu)
  OPEN (iu,FILE="BLAZE.nml",STATUS='OLD',ACTION='READ')
  READ (iu,NML=BLAZENML)
  CLOSE(iu)

  ! READ ini-nml
  BLAZE%NCELLS = np
  BLAZE%time = 0
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
  ALLOCATE ( BLAZE%FFDI    ( np ) )
  ALLOCATE ( BLAZE%FLIx    ( np ) )
  ALLOCATE ( BLAZE%DFLI    ( np ) )
  ALLOCATE ( BLAZE%ROS     ( np ) )
  ALLOCATE ( BLAZE%Z       ( np ) )
  ALLOCATE ( BLAZE%D       ( np ) )
  ALLOCATE ( BLAZE%w       ( np ) )
  ALLOCATE ( BLAZE%w_prior       ( np ) )
  ALLOCATE ( BLAZE%TO      ( np, NTO ) )
  ALLOCATE ( BLAZE%AnnRainf( np, 366 ) )
  ALLOCATE ( BLAZE%AvgAnnRAINF(np, 5))
  ALLOCATE ( BLAZE%CAvgAnnRainf(np))
  ALLOCATE ( BLAZE%DEADWOOD( np ) )
  ALLOCATE ( BLAZE%SHOOTFRAC( np ) )
  ! POP related vars
  ALLOCATE ( BLAZE%POP_TO  ( np ) )
  ALLOCATE ( BLAZE%POP_CWD ( np ) )
  ALLOCATE ( BLAZE%POP_STR ( np ) )
  ALLOCATE ( BLAZE%FLUXES  ( np, 13 ) )

  ALLOCATE( BLAZE%AGLit_g(np,NPOOLS),BLAZE%AGLit_w(np,NPOOLS) )
  ALLOCATE( BLAZE%BGLit_g(np,NPOOLS),BLAZE%BGLit_w(np,NPOOLS) )
  ALLOCATE( BLAZE%CPLANT_g(np,NPOOLS),BLAZE%CPLANT_w(np,NPOOLS) )

  call zero_blaze(BLAZE)

  ! SETTINGS FOR BLAZE (BLAZEFLAG)
  ! bit value:               0            | 1
  ! 0th bit(1), general    : off          | on
  ! 1st bit(2): BLAZE-mode : FLI-only     | Full-BLAZE
  ! 2nd bit(4): Ignition   : other (GFED) | SIMFIRE
  ! more possible...for other ignition sources
  ! examples:
  ! BLAZE(+1) with FLI-only(+0) and SIMFIRE(+4) -> BLAZEFLAG=5
  ! BLAZE(+1) full (+2) with GFED (+0)          -> BLAZEFLAG=3


  BLAZE%BURNMODE = 1    ! Full. All Fluxes computed by BLAZE
  WRITE(*,*) " BLAZE: Full Mode"
  WRITE(*,*) " Burnt-area source: ", TRIM(BurnedAreaSource)

  BLAZE%DSLR      = 0
  BLAZE%RAINF     = 0.
  BLAZE%LR        = 0.
  BLAZE%KBDI      = 0.
  BLAZE%AB        = 0.
  BLAZE%DEADWOOD  = 0.
  BLAZE%FSTEP     = 'DAILY'

  BLAZE%LAT       = LAT
  BLAZE%LON       = LON

!!$  ! time for averaging annual rainfall [a]
!!$  BLAZE%T_AVG = 5
!!$  ! deadwood decay scale time [a]
!!$  BLAZE%FT    = 30

  BLAZE%BURNT_AREA_SRC = trim(BurnedAreaSource)

  BLAZE%CPLANT_g = 0.0
  BLAZE%CPLANT_w = 0.0

END SUBROUTINE INI_BLAZE


subroutine zero_blaze(blaze)

  type(type_blaze), intent(inout) :: blaze

  blaze%DSLR         = 0
  blaze%ilon         = 0
  blaze%jlat         = 0
  blaze%Flix         = 0
  blaze%RAINF        = 0
  blaze%KBDI         = 0
  blaze%LR           = 0
  blaze%U10          = 0
  blaze%RH           = 0
  blaze%TMAX         = 0
  blaze%TMIN         = 0
  blaze%AREA         = 0
  blaze%w_prior      = 0
  !blaze%FDI          = 0
  blaze%FFDI         = 0
  blaze%FLI          = 0
  blaze%ROS          = 0
  blaze%Z            = 0
  blaze%D            = 0
  blaze%w            = 0
  blaze%LAT          = 0
  blaze%LON          = 0
  blaze%DFLI         = 0
  blaze%AB           = 0
  Blaze%CAvgAnnRainf = 0
  blaze%DEADWOOD     = 0
  blaze%POP_TO       = 0
  blaze%POP_CWD      = 0
  blaze%POP_STR      = 0
  blaze%shootfrac    = 0
  Blaze%AnnRAINF     = 0
  !blaze%ABM          = 0
  blaze%TO           = 0
  !blaze%AGC_g        = 0
  !blaze%AGC_w        = 0
  blaze%AGLit_w      = 0
  blaze%AGLit_g      = 0
  blaze%BGLit_w      = 0
  blaze%BGLit_g      = 0
  blaze%CPLANT_g     = 0
  blaze%CPLANT_w     = 0
  blaze%AvgAnnRAINF  = 0
  blaze%FLUXES       = 0

end subroutine zero_blaze


subroutine print_blaze(blaze)

  type(type_blaze), intent(in) :: blaze

  write(*,*) 'DSLR ', blaze%DSLR
  write(*,*) 'ilon ', blaze%ilon
  write(*,*) 'jlat ', blaze%jlat
  write(*,*) 'Flix ', blaze%Flix
  write(*,*) 'RAINF ', blaze%RAINF
  write(*,*) 'KBDI ', blaze%KBDI
  write(*,*) 'LR ', blaze%LR
  write(*,*) 'U10 ', blaze%U10
  write(*,*) 'RH ', blaze%RH
  write(*,*) 'TMAX ', blaze%TMAX
  write(*,*) 'TMIN ', blaze%TMIN
  write(*,*) 'AREA ', blaze%AREA
  write(*,*) 'w_prior ', blaze%w_prior
  write(*,*) 'FFDI ', blaze%FFDI
  write(*,*) 'FLI ', blaze%FLI
  write(*,*) 'ROS ', blaze%ROS
  write(*,*) 'Z ', blaze%Z
  write(*,*) 'D ', blaze%D
  write(*,*) 'w ', blaze%w
  write(*,*) 'LAT ', blaze%LAT
  write(*,*) 'LON ', blaze%LON
  write(*,*) 'DFLI ', blaze%DFLI
  write(*,*) 'AB ', blaze%AB
  write(*,*) 'CAvgAnnRainf ', blaze%CAvgAnnRainf
  write(*,*) 'DEADWOOD ', blaze%DEADWOOD
  write(*,*) 'POP_TO ', blaze%POP_TO
  write(*,*) 'POP_CWD ', blaze%POP_CWD
  write(*,*) 'POP_STR ', blaze%POP_STR
  write(*,*) 'shootfrac ', blaze%shootfrac
  write(*,*) 'AnnRAINF ', blaze%AnnRAINF
  write(*,*) 'TO ', blaze%TO
  write(*,*) 'AGLit_w ', blaze%AGLit_w
  write(*,*) 'AGLit_g ', blaze%AGLit_g
  write(*,*) 'BGLit_w ', blaze%BGLit_w
  write(*,*) 'BGLit_g ', blaze%BGLit_g
  write(*,*) 'CPLANT_g ', blaze%CPLANT_g
  write(*,*) 'CPLANT_w ', blaze%CPLANT_w
  write(*,*) 'AvgAnnRAINF ', blaze%AvgAnnRAINF
  write(*,*) 'FLUXES ', blaze%FLUXES

end subroutine print_blaze


SUBROUTINE BLAZE_ACCOUNTING(BLAZE, climate,  ktau, dels, year, doy)

  USE CABLE_DEF_TYPES_MOD, ONLY: climate_type
  USE cable_IO_vars_module, ONLY:  landpt

  IMPLICIT NONE

  TYPE (TYPE_BLAZE)   :: BLAZE
  TYPE (CLIMATE_TYPE), INTENT(IN)     :: climate
  INTEGER, INTENT(IN) :: ktau, year, doy
  REAL,    INTENT(IN) :: dels
  INTEGER             :: mp
  REAL                :: t_fac
  INTEGER,PARAMETER   :: sod = 86400
  LOGICAL             :: is_new_day, is_end_of_day

  mp = BLAZE%NCELLS

  is_new_day    = (MOD((ktau-1) * NINT(dels), sod) == 0 )
  is_end_of_day = (MOD( ktau    * NINT(dels), sod) == 0 )

  ! reset values for each BLAZE time step
  IF ( is_new_day ) THEN
     IF ( ( TRIM(BLAZE%FSTEP) .EQ. "annual" .AND. doy == 1 ) .OR. &
          ( TRIM(BLAZE%FSTEP) .EQ. "daily" ) ) THEN
        BLAZE%FFDI = 0.
        ! ??? Maybe more ...
     ENDIF
  ENDIF

  ! factor to get to daily data
  t_fac = 86400. / dels

  ! assign values of BLAZE variables to those held in the climate structure
  BLAZE%Rainf(:) = climate%dprecip(landpt(:)%cstart)
  BLAZE%CAvgAnnRainf(:) = climate%aprecip_av20(landpt(:)%cstart)
  BLAZE%U10(:) = climate%du10_max(landpt(:)%cstart)
  BLAZE%U10 = BLAZE%u10 * 3.6 ! m/s -> km/h
   ! Wind T. McVicar 201?
   ! Conversion to Windmax following S. Matthews, 2014 (pers. comm. so far)
   ! in km/h
  BLAZE%U10 = ( 214.7 * ( BLAZE%U10 + 10. ) ** (-1.6968 ) + 1. ) * BLAZE%U10

  BLAZE%RH(:) = climate%drhum(landpt(:)%cstart)
  BLAZE%TMAX(:) = climate%dtemp_max(landpt(:)%cstart)
  BLAZE%TMIN(:) = climate%dtemp_min(landpt(:)%cstart)
  BLAZE%KBDI(:) = climate%KBDI(landpt(:)%cstart)
  BLAZE%D(:)    = climate%D_MacArthur(landpt(:)%cstart)
  BLAZE%FFDI(:) = climate%FFDI(landpt(:)%cstart)
  BLAZE%DSLR(:) = climate%DSLR(landpt(:)%cstart)
  BLAZE%LR(:) =  climate%last_precip(landpt(:)%cstart)

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

  IMPLICIT NONE

  !get index of relative live-biomass loss in TOF table plus fract between 2bins
  !-> apply idx and factor on other TOs

  TYPE(TYPE_TURNOVER),INTENT(INOUT) :: TO(7)
  INTEGER,            INTENT(IN)    :: BURNMODE
  REAL,               INTENT(IN)    :: AB, shootfrac
  REAL,               INTENT(INOUT) :: CPLANT_w(3) , CPLANT_g(3)
  REAL, DIMENSION(3), INTENT(INOUT) :: AGL_w, AGL_g, BGL_w, BGL_g
  REAL,               INTENT(OUT)   :: BT(13)
  REAL,     OPTIONAL, INTENT(IN)    :: POP_TO
  TYPE(TYPE_TURNOVER)               :: MTO(7)
  REAL    :: fAB

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  fAB = AB

  ! Mass Fluxes
  MTO(LEAF )%TO_ATM = fAB * TO(LEAF )%TO_ATM * CPLANT_w (LEAF )
  MTO(WOOD )%TO_ATM = fAB * TO(WOOD )%TO_ATM * CPLANT_w (WOOD )
  MTO(LEAF )%TO_STR = fAB * TO(LEAF )%TO_STR * CPLANT_w (LEAF )
  !CVH should the line below use only the above-ground wood pool?
  MTO(WOOD )%TO_STR = fAB * TO(WOOD )%TO_STR * CPLANT_w (WOOD )
  MTO(WOOD )%TO_CWD = fAB * TO(WOOD )%TO_CWD * CPLANT_w (WOOD )
  MTO(FROOT)%TO_ATM = fAB * TO(FROOT)%TO_ATM * CPLANT_w (FROOT)
  MTO(FROOT)%TO_STR = fAB * TO(FROOT)%TO_STR * CPLANT_w (FROOT)
  MTO(MLIT )%TO_ATM = fAB * TO(MLIT )%TO_ATM * AGL_w(METB )
  MTO(SLIT )%TO_ATM = fAB * TO(SLIT )%TO_ATM * AGL_w(STR  )
  MTO(CLIT )%TO_ATM = fAB * TO(CLIT )%TO_ATM * AGL_w(CWD  )

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
!CVH omit: pools are re-contructed each time from casa pools
!!$  ! Grass fuels burned by area
!!$  CPLANT_g (LEAF)  = CPLANT_g (LEAF) * (1. - AB)
!!$  AGL_g(METB)      = AGL_g(METB)     * (1. - AB)
!!$  AGL_g(STR )      = AGL_g(STR )     * (1. - AB)
!!$
!!$  ! Leaves
!!$  CPLANT_w (LEAF)  = &
!!$       CPLANT_w (LEAF) - (MTO(LEAF)%TO_ATM + MTO(LEAF)%TO_STR)
!!$  ! Wood
!!$  CPLANT_w (WOOD)  = &
!!$       CPLANT_w (WOOD) - (MTO(WOOD)%TO_ATM + MTO(WOOD)%TO_STR + MTO(WOOD)%TO_CWD)
!!$  ! Fine roots
!!$  CPLANT_w (FROOT) = &
!!$       CPLANT_w (FROOT)- (MTO(FROOT)%TO_ATM+ MTO(FROOT)%TO_STR)
!!$
!!$  ! Litter above ground
!!$  AGL_w(METB)      = AGL_w(METB) - MTO(MLIT)%TO_ATM
!!$  AGL_w(STR )      = AGL_w(STR ) - MTO(SLIT)%TO_ATM + MTO(WOOD)%TO_STR + &
!!$                     MTO(LEAF)%TO_STR
!!$  AGL_w(CWD )      = AGL_w(CWD ) - MTO(CLIT)%TO_ATM + MTO(WOOD)%TO_CWD * shootfrac
!!$
!!$  ! Litter below ground
!!$  BGL_w(STR ) = BGL_w(STR ) + MTO(FROOT)%TO_STR
!!$  BGL_w(CWD ) = BGL_w(CWD ) + MTO(WOOD)%TO_CWD * (1. - shootfrac)
!!$
!!$  IF ( ANY ( AGL_w .LT. 0. )) THEN
!!$     WRITE(*,*)"CLIITE_W < 0 ", AGL_w
!!$     STOP 1
!!$  ENDIF


END SUBROUTINE BLAZE_TURNOVER


FUNCTION p_surv_OzSavanna(hgt, fli)

  implicit none

  real, intent(in) :: hgt, fli
  real             :: p_surv_OzSavanna

  real :: max_prob_hgt, intensity, min_hgt

  !  Cook fire-mortality
  max_prob_hgt = 8.5
  intensity    = fli / 1000.
  min_hgt      = 3.7 * (1.-exp(-0.19 * intensity))

  if ( hgt > max_prob_hgt .and. hgt > min_hgt ) then
     p_surv_OzSavanna = ( -.0011 * intensity - .00002) * hgt + .0075 * intensity + 1.
  elseif (hgt > min_hgt) then
     p_surv_OzSavanna = (  .0178 * intensity + .0144 ) * &
          hgt + (-.1174 * intensity + 0.9158)
  else
     p_surv_OzSavanna = 0.001
  endif

  p_surv_OzSavanna = max(0.001, min(1.,p_surv_OzSavanna))

END FUNCTION p_surv_OzSavanna



FUNCTION BURNTIME( YEAR, DOY, FSTEP )

  USE CABLE_COMMON_MODULE, ONLY : DOYSOD2YMDHMS, IS_LEAPYEAR

  implicit none

  INTEGER,          INTENT(IN) :: YEAR, DOY
  CHARACTER(LEN=*), INTENT(IN) :: FSTEP
  LOGICAL                      :: BURNTIME

  INTEGER :: MM, DD, EOM
  INTEGER, DIMENSION(12), PARAMETER :: LOM = (/31,28,31,30,31,30,31,31,30,31,30,31/)

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
     WRITE(6,FMT='(A36,I2.2,1x,I2.2,1x,I4)') &
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

  REAL, PARAMETER :: H = 20. ! Heat Yield, [MJ/kg]
  REAL      :: FT
  REAL      :: F, w, ROS, Z
  REAL      :: FLI
  INTEGER   :: FLIx, i

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

  ! vh KBDI now calculated in cable_climate.F90

  ! Keetch-Byram Drought Index [ ]
!!$  ! Keetch, j.j. and G.M. Byram, 1968, A drought index for forest fire control,
!!$  ! US Dept, Agr. Forest Service Res. Paper SE-38
!!$  IF ( BLAZE%RAINF(np) .GT. 0. ) THEN
!!$     IF ( BLAZE%DSLR(np) .GT. 0 ) THEN
!!$        BLAZE%LR(np) = BLAZE%RAINF(np)
!!$     ELSE
!!$        BLAZE%LR(np) = BLAZE%LR(np) + BLAZE%RAINF(np)
!!$     ENDIF
!!$     BLAZE%DSLR(np) = 0
!!$  ELSE
!!$     BLAZE%DSLR(np) = BLAZE%DSLR(np) + 1
!!$  ENDIF
!!$
!!$  DSLR      = BLAZE%DSLR(np)
!!$  AvAnnRain = BLAZE%CAvgAnnRainf(np)
!!$  V         = BLAZE%U10(np)
!!$  RH        = BLAZE%RH(np)
!!$  T         = BLAZE%TMAX(np)
!!$
!!$  IF ( DSLR .EQ. 0 ) THEN
!!$     IF ( BLAZE%LR(np) .GE. 5 ) THEN
!!$        dKBDI = 5. - BLAZE%LR(np)
!!$     ELSE
!!$        dKBDI = 0.
!!$     ENDIF
!!$  ELSE
!!$     dKBDI = ((800. - BLAZE%KBDI(np)) * (.968 * EXP(.0486 * (T * 9./5. &
!!$          + 32.)) - 8.3) / 1000. / (1. + 10.88 * EXP(-.0441 * AvAnnRain/25.4)) * .254)
!!$  ENDIF
!!$
!!$  print*, 'KBDI', BLAZE%KBDI(np) , dKBDI, BLAZE%LR(np), T, AvAnnRain
!!$  BLAZE%KBDI(np) = MAX(0.,BLAZE%KBDI(np) + dKBDI)

!  IF ( .NOT. BURN ) RETURN

!!$  ! McArthur drought factor [ ]
!!$  D   = .191 * ( BLAZE%KBDI(np) + 104. ) * ( REAL(DSLR) + 1. )**1.5 / &
!!$       ( 3.52 * ( REAL(DSLR) + 1. )**1.5 + BLAZE%LR(np) - 1. )
!!$  D   = MAX(0.,MIN(10.,D))

  ! McArthur Fire-Danger index [ ] (FFDI)
  !                             Forest  Grass
  ! Catastrophic (Code Red)	100 +   150 +
  ! Extreme	                75 - 99 100 - 149
  ! Severe	                50 - 74 50 - 99
  ! Very High	                25 - 49 25 - 49
  ! High	                12 - 24 12 - 24
  ! Low-Moderate            	 0 - 11  0 - 11

  F = BLAZE%FFDI(np)

!!$  F   = 2. * EXP( -.45 + .987 * LOG(D+.001) - .03456 * RH + .0338 * T + .0234 * V ) ! V in km/h
!!$  F   = MAX(0.,F)


  ! available Fuels (Litter, deadwood and under/midstorey[when is it distinguished?])) [kg/m^2]
  ! Assume totel Grass is burnt
  ! Factor 0.6 / 0.5  taken from fires < 750 kW/m of Suravski et al. 2012
  FLIx = 1

  w = AVAIL_FUEL(FLIx, CPLANT_w, CPLANT_g, BLAZE%AGLit_w(np,:), BLAZE%AGLit_g(np,:) )

  BLAZE%w_prior(np) = w

  !write(59,*) w, CPLANT_w, CPLANT_g, BLAZE%AGLit_w, BLAZE%AGLit_g

  DO i = 1, 4 ! THE LADDER

     ! Rate of Spread [m/s]
     ! original ROS = 0.0012 * F * w  (Noble 1980: w in tonnes ha-1, ROS in km/h)
     ROS = 2.4e-5 * F* W / 3.6 ! (vh: conversions factor of 1/50 applied to go from gCm-2 to T dm ha-1 )
     !    and 1/3.6 to go from km/h to ms-1)


     !ROS = 3.3333E-05 * F * w  (LN)

     ! Flame Height  [m]
     ! original Z   = 13. * ROS + 0.24 * w - 2.  (Z in m; ROS in km/h; w in t ha-1)
     Z = 13.0 * ROS * 3.6 + 0.24 * w / 50.0 -2. ! vh (conversion factors of 3.6 for ROS &
                                                !   50.0 for w)

     !Z   = 46.8 * ROS + 0.024 * w - 2. (LN)
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
  !BLAZE%FFDI(np) = F
  BLAZE%FLI(np)  = FLI
  BLAZE%ROS(np)  = ROS
  BLAZE%Z(np)    = Z
  !BLAZE%D(np)    = D
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
  ! Get combustion factors from FLI (Surawski 2012) AS IS
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! TURN-OVER RATES
  !
  !print*, FLIx, CPLANT_w, CPLANT_g, BLAZE%AGLit_w(1,:), BLAZE%AGLit_g(1,:)
  !print*, 'FLI, FLIX, w, MIN_FUEL:', FLI, FLIX, w, MIN_FUEL

  IF (FLIx .gt.0) then
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

ELSE
   TO(LEAF)%TO_ATM    = 0.0
   TO(WOOD)%TO_ATM    = 0.0
   TO(FROOT)%TO_ATM   = 0.0
   TO(LEAF)%TO_STR    = 0.0
   TO(WOOD)%TO_CWD    = 0.0
   TO(FROOT)%TO_STR   = 0.0
   TO(MLIT)%TO_ATM    = 0.0
   TO(SLIT)%TO_ATM    = 0.0
   TO(CLIT)%TO_ATM    = 0.0
ENDIF
END SUBROUTINE COMBUST


SUBROUTINE RUN_BLAZE(BLAZE, SF, CPLANT_g, CPLANT_w, tstp, YYYY, doy, TO , climate)

  USE CABLE_COMMON_MODULE, ONLY: IS_LEAPYEAR, DOYSOD2YMDHMS
  USE SIMFIRE_MOD
  USE CABLE_DEF_TYPES_MOD, ONLY:  climate_type

  IMPLICIT NONE

!CRM  INTEGER, INTENT(IN) :: POPFLAG
  TYPE(TYPE_BLAZE)    :: BLAZE
  TYPE(TYPE_SIMFIRE)  :: SF
  TYPE(TYPE_TURNOVER) :: TO(BLAZE%NCELLS,7)
  TYPE (climate_type ), INTENT(IN)         :: climate

  INTEGER          :: np, doy, YYYY, MM, DD, DOM(12)
  REAL             :: CPLANT_g(BLAZE%NCELLS,3), CPLANT_w(BLAZE%NCELLS,3), tstp
  REAL, DIMENSION(BLAZE%NCELLS,3) :: AGL_g, AGL_w, BGL_g, BGL_w
  REAL, DIMENSION(BLAZE%NCELLS)   :: &
       RAINF,             & ! [mm/d]
       TMIN,              & ! [deg C]
       TMAX,              & ! [deg C]
       AB,                & ! [frac.]
       popd,              & ! [kW/m]
       mnest

  ! CLN to simfire!!!
  !CRM  INTEGER,DIMENSION(BLAZE%NCELLS)  :: modis_igbp  ! [0,17] landcover index

  DOM = (/ 31,28,31,30,31,30,31,31,30,31,30,31 /)
  IF ( IS_LEAPYEAR(YYYY) ) DOM(2) = 29

  BLAZE%YEAR = YYYY
  BLAZE%DOY  = DOY
  CALL DOYSOD2YMDHMS( YYYY, doy,10, MM, DD )

  BLAZE%DAY   = DD
  BLAZE%MONTH = MM
  BLAZE%tstp = tstp



!! CRM  IF ( DD .EQ. 1 ) BLAZE%FLI(:) = 0. now in accounting

  ! READ GFED BA DATA
  AB(:) = 0.
  IF ( TRIM(BLAZE%BURNT_AREA_SRC) .EQ. "GFED3.1" ) THEN
!     CALL GET_GFED( BLAZE )
!     CALL GET_GFED4_BA( BLAZE )
!     CALL GET_GFED41s_BA( BLAZE )
     WRITE(*,*)'GFED4 BA not available. Set cable_user%BURNT_AREA == "SIMFIRE"'
     STOP -1
     IF ( TRIM(BLAZE%FSTEP) .EQ. "none" ) THEN
        CALL SIMFIRE ( SF, RAINF, TMAX, TMIN, DOY,MM, YYYY, BLAZE%AB, climate )
        popd = SF%POPD
        mnest= SF%MAX_NESTEROV
        BLAZE%FSTEP = "annual"
     ENDIF
  ELSEIF ( TRIM(BLAZE%BURNT_AREA_SRC) .EQ. "SIMFIRE" ) THEN
     ! CALL SIMFIRE DAILY FOR ACOUNTING OF PARAMETERS
     CALL SIMFIRE ( SF, RAINF, TMAX, TMIN, DOY,MM, YYYY, BLAZE%AB , climate)


     DO np = 1, BLAZE%NCELLS
        IF ( AVAIL_FUEL(1, CPLANT_w(np,:), CPLANT_g(np,:),BLAZE%AGLit_w(np,:),BLAZE%AGLit_g(np,:) ) .LE. MIN_FUEL ) &
             BLAZE%AB(np) = 0.
     END DO
     popd = SF%POPD
     mnest= SF%MAX_NESTEROV
  ELSE
     WRITE(*,*) "Wrong ignition type chosen: ", BLAZE%BURNT_AREA_SRC
     STOP -1
  ENDIF

  ! Apply half of former deadwood to atm now How to distribut (str
  ! set following Fraver 2013 pinus rosinosa (hardwood/decid. wood to be added

  DO np = 1, BLAZE%NCELLS



     CALL COMBUST( BLAZE, np, CPLANT_g(np,:), CPLANT_w(np,:),TO(np,:),BLAZE%AB(np).GT.0 )



     BLAZE%DFLI(np) = BLAZE%FLI(np)
     BLAZE%TO(np,:) = 0.

     IF (BLAZE%AB(np) .GT. 0. ) THEN
        CALL BLAZE_TURNOVER( BLAZE%AB(np), CPLANT_g(np,:), CPLANT_w(np,:), &
             AGL_g(np,:), AGL_w(np,:),BGL_g(np,:), BGL_w(np,:), &
             BLAZE%shootfrac(np), TO(np,:), BLAZE%FLUXES(np,:), BLAZE%BURNMODE )
     ENDIF



     AB(np) = BLAZE%AB(np)

     ! Apply other half of former deadwood to litter now How to distribut (str
!CLN     AGL(np,CWD) = AGL(np,CWD) + &
!CLN          !       (1.-exp(-0.5*tstp/SF(ft))) * SUM(DEADWOOD(np,:))
!CLN          (1.-exp(-0.5*tstp/15.)) * SUM(DEADWOOD(np,:))

!CLN      BLAZE%DEADWOOD(np,:) = DEADWOOD(np,:) * exp(-0.5*tstp/15.)

  END DO


!  IF ( BLAZE%DAY .EQ. 1 ) &
!       CALL BLAZE_DIAG( NCELLS, BLAZE, CPLANT_g, CPLANT_w, AGL_g, AGL_w, TO, "dref")

END SUBROUTINE RUN_BLAZE
!*******************************************************************************
  SUBROUTINE WRITE_BLAZE_OUTPUT_NC ( BLAZE, FINAL )

    USE CABLE_COMMON_MODULE, ONLY:  filename, cable_user, HANDLE_ERR
    USE cable_IO_vars_module, ONLY: timeunits, calendar
    USE netcdf

    IMPLICIT NONE

    TYPE(TYPE_BLAZE), INTENT(IN) :: BLAZE
    LOGICAL, INTENT(IN)    :: FINAL
    !INTEGER, INTENT(IN)    :: ctime

    INTEGER   :: STATUS
    INTEGER   :: land_ID, t_ID
    INTEGER   :: i, mp
    CHARACTER :: FNAME*99,dum*50
    LOGICAL, SAVE :: CALL1 = .TRUE.

     ! 1 dim arrays (mp )
    CHARACTER(len=20),DIMENSION(2) :: A0
    ! 2 dim real arrays (mp,t)
    CHARACTER(len=20),DIMENSION(25):: A1
    ! 2 dim integer arrays (mp,t)
    CHARACTER(len=20),DIMENSION(1):: AI1

    INTEGER, SAVE :: VIDtime, VID0(SIZE(A0)),VID1(SIZE(A1)),VIDI1(SIZE(AI1))
    INTEGER, SAVE :: FILE_ID, CNT = 0

    A0(1) = 'latitude'
    A0(2) = 'longitude'

    A1(1) = 'BurnedArea'
    A1(2) = 'Rainf'
    A1(3) = 'KBDI'
    A1(4) = 'LastRain'
    A1(5) = 'U10'
    A1(6) = 'RelHum'
    A1(7) = 'TMAX'
    A1(8) = 'TMIN'
    A1(9) = 'FFDI'
    A1(10) = 'FLI'
    A1(11) = 'ROS'
    A1(12) = 'FlameHeight'
    A1(13) = 'DMacArthur'
    A1(14) = 'AvailFuel'
    A1(15) = 'AvailFuelPrior'

    A1(16) = 'CPLANT_w_froot'
    A1(17) = 'CPLANT_w_leaf'
    A1(18) = 'CPLANT_w_wood'
    A1(19) = 'CPLANT_g_froot'
    A1(20) = 'CPLANT_g_leaf'
    A1(21) = 'AGL_w_metb'
    A1(22) = 'AGL_w_STR'
    A1(23) = 'AGL_w_CWD'
    A1(24) = 'AGL_g_metb'
    A1(25) = 'AGL_g_str'

    AI1(1) = 'DaysSinceLastRain'

    CNT = CNT + 1

    mp = BLAZE%ncells


     IF ( CALL1 ) THEN
        ! Get File-Name

       WRITE( dum, FMT="(I4,'_',I4)")CABLE_USER%YEARSTART,CABLE_USER%YEAREND
       IF (CABLE_USER%YEARSTART.lt.1000.and.CABLE_USER%YEAREND.lt.1000) THEN
          WRITE( dum, FMT="(I3,'_',I3)")CABLE_USER%YEARSTART,CABLE_USER%YEAREND
       ELSEIF (CABLE_USER%YEARSTART.lt.1000) THEN
          WRITE( dum, FMT="(I3,'_',I4)")CABLE_USER%YEARSTART,CABLE_USER%YEAREND
       ENDIF

       IF  (LEN_TRIM( TRIM(cable_user%BLAZE_outfile) ) .gt. 0 ) THEN
          fname = TRIM(cable_user%BLAZE_outfile)
       ELSE
          fname = TRIM(filename%path)//'/'//TRIM(cable_user%RunIden)//'_'//&
               TRIM(dum)//'_BLAZE_out.nc'
       ENDIF

       ! Create NetCDF file:
       STATUS = NF90_create(fname, ior(nf90_clobber,nf90_64bit_offset), FILE_ID)
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "StartYear", CABLE_USER%YEARSTART )
       STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "EndYear"  , CABLE_USER%YEAREND   )
       STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "RunIden"  , CABLE_USER%RunIden   )

       ! Define dimensions:
       ! Land (number of points)
       STATUS = NF90_def_dim(FILE_ID, 'land'   , mp     , land_ID)
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

       STATUS = NF90_def_dim(FILE_ID, 'time'   , NF90_UNLIMITED, t_ID)
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

        ! Define variables
       STATUS = NF90_def_var(FILE_ID,'time' ,NF90_INT,(/t_ID/),VIDtime )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

       !write(*,*) 'timeunits', TRIM(timeunits), t_ID, FILE_ID

       STATUS = NF90_PUT_ATT(FILE_ID, VIDtime, 'units', TRIM(timeunits))
       IF (STATUS /= NF90_NOERR)  CALL handle_err(STATUS)

       STATUS = NF90_PUT_ATT(FILE_ID,VIDtime, 'calendar', calendar)
       IF (STATUS /= NF90_NOERR)  CALL handle_err(STATUS)



       DO i = 1, SIZE(A0)
          STATUS = NF90_def_var(FILE_ID,TRIM(A0(i)) ,NF90_FLOAT,(/land_ID/),VID0(i))
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
           write(*,*) 'def A0'
       END DO

       DO i = 1, SIZE(A1)
          STATUS = NF90_def_var(FILE_ID,TRIM(A1(i)) ,NF90_FLOAT,(/land_ID,t_ID/),VID1(i))
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       END DO

        DO i = 1, SIZE(AI1)
           STATUS = NF90_def_var(FILE_ID,TRIM(AI1(i)) ,NF90_INT,(/land_ID,t_ID/),VIDI1(i))
           IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
           write(*,*) 'def AI1'
        END DO

       ! End define mode:
       STATUS = NF90_enddef(FILE_ID)
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)


       ! PUT LAT / LON ( mp )
       STATUS = NF90_PUT_VAR(FILE_ID, VID0(1), BLAZE%LAT )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID0(2), BLAZE%LON )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       CALL1 = .FALSE.

    ENDIF ! CALL1



    ! TIME  ( t )
    STATUS = NF90_PUT_VAR(FILE_ID, VIDtime, BLAZE%time, start=(/ CNT /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)


     ! PUT 2D VARS ( mp, t )
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 1), BLAZE%AB, start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 2), BLAZE%Rainf, start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 3), BLAZE%KBDI, start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)


    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 4), BLAZE%LR, start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 5), BLAZE%U10, start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 6), BLAZE%RH, start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 7), BLAZE%TMAX, start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 8), BLAZE%TMIN, start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 9), BLAZE%FFDI, start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 10), BLAZE%FLI, start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 11), BLAZE%ROS, start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 12), BLAZE%Z, start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 13), BLAZE%D, start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 14), BLAZE%w, start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 15), BLAZE%w_prior, start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

!!$    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 15), BLAZE%CAvgAnnRainf, start=(/ 1, CNT /), &
!!$         count=(/ mp, 1 /) )


    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 16), BLAZE%CPLANT_w(:,FROOT), start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 17), BLAZE%CPLANT_w(:,leaf), start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 18), BLAZE%CPLANT_w(:,WOOD), start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 19), BLAZE%CPLANT_g(:,froot), start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 20), BLAZE%CPLANT_g(:,leaf), start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 21), BLAZE%AGLit_w(:,metb), start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 22), BLAZE%AGlit_w(:,str), start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 23), BLAZE%AGlit_w(:,CWD), start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 24), BLAZE%AGlit_g(:,metb), start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 25), BLAZE%AGlit_g(:,str), start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)




    STATUS = NF90_PUT_VAR(FILE_ID, VIDI1(1), BLAZE%DSLR, start=(/ 1, CNT /), count=(/ mp, 1 /)  )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    IF ( FINAL ) THEN
       ! Close NetCDF file:
       STATUS = NF90_close(FILE_ID)
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       WRITE(*,*) " BLAZE Output written to ",fname
    ENDIF

 END SUBROUTINE WRITE_BLAZE_OUTPUT_NC


END MODULE BLAZE_MOD
