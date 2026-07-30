MODULE casa_readbiome_module

CONTAINS

! jh:is ntiles used here actually ntiles?
SUBROUTINE casa_readbiome( mp, nsl, ntiles, veg_iveg, soil_zse, casabiome,     &
                           casapool, casaflux, casamet, phen ) 

!USE cable_def_types_mod !jh? for ms, ntiles
    USE casaparm

!subrs
!jh!USE pft_lookup_parms_csv1_mod_cnp, ONLY: casa_pft_params_csv1 

!data
USE casa_files_type_mod, ONLY: casafile   ! TYPE instance 
USE casavariable,        ONLY: casa_biome ! TYPE declaration
USE casavariable,        ONLY: casa_pool  ! TYPE declaration
USE casavariable,        ONLY: casa_flux  ! TYPE declaration
USE casavariable,        ONLY: casa_met   ! TYPE declaration
   
USE cable_common_module, ONLY: cable_user
USE casadimension, ONLY: mso, msoil, mplant
USE casadimension, ONLY: deltcasa 
USE casadimension, ONLY: icycle
USE phenology_type_mod, ONLY: phenology_type 

IMPLICIT NONE

INTEGER, INTENT(IN) :: mp                ! # active tile points  
INTEGER, INTENT(IN) :: nsl               ! # soil layers              
INTEGER, INTENT(IN) :: ntiles            ! # SurfaceTypes
INTEGER, INTENT(IN) :: veg_iveg(mp)      ! SurfaceType index
REAL,    INTENT(IN) :: soil_zse(nsl)     ! soil depth per layer

TYPE (casa_biome),          INTENT(OUT) :: casabiome
TYPE (casa_pool),           INTENT(OUT) :: casapool
TYPE (casa_flux),           INTENT(OUT) :: casaflux
TYPE (casa_met),            INTENT(OUT) :: casamet
TYPE (phenology_type),      INTENT(OUT) :: phen

    ! local variables
    REAL(r_2), DIMENSION(ntiles)       :: leafage,frootage,woodage
    REAL(r_2), DIMENSION(ntiles)       :: totroot
    REAL(r_2), DIMENSION(ntiles)       :: cwdage,metage,strage
    REAL(r_2), DIMENSION(ntiles)       :: micage,slowage,passage,clabileage,slax
    REAL(r_2), DIMENSION(ntiles,mplant):: ratioCNplant
    REAL(r_2), DIMENSION(ntiles,msoil) :: ratioCNsoil,ratioCNsoilmin,ratioCNsoilmax
    REAL(r_2), DIMENSION(nsl)           :: depthsoila,depthsoilb
    REAL(r_2), DIMENSION(ntiles)       :: xfNminloss, xfNminleach, xnfixrate
    REAL(r_2), DIMENSION(ntiles)       :: cleaf,cwood,cfroot,      &
         cmet,cstr,ccwd,          &
         cmic,cslow,cpass
    REAL(r_2), DIMENSION(ntiles)       :: nleaf,nwood,nfroot,      &
         nmet,nstr,ncwd,          &
         nmic,nslow,npass,xnsoilmin
    REAL(r_2), DIMENSION(ntiles)       :: xpleaf, xpwood, xpfroot, &
         xpmet, xpstr, xpcwd,     &
         xpmic,xpslow,xppass,xplab,xpsorb,xpocc
    REAL(r_2), DIMENSION(mso)       :: xkmlabp,xpsorbmax,xfPleach
    REAL(r_2), DIMENSION(mso,msoil) :: ratioNPsoil
    REAL(r_2), DIMENSION(ntiles)       :: xfherbivore,xxkleafcoldmax, xxkleafdrymax
    REAL(r_2), DIMENSION(ntiles)       :: xkuplabp
    REAL(r_2), DIMENSION(ntiles,nsl)    :: fracroot
    REAL(r_2) ::  xratioNPleafmin,xratioNPleafmax,         &
         xratioNPwoodmin,xratioNPwoodmax,         &
         xratioNPfrootmin,xratioNPfrootmax
    INTEGER :: i,iv1,nv,ns,npt,iv,is,iso
    INTEGER :: nv0,nv1,nv2,nv3,nv4,nv5,nv6,nv7,nv8,nv9,nv10,nv11,nv12
    REAL(r_2), DIMENSION(ntiles)       :: xxnpmax,xq10soil,xxkoptlitter,xxkoptsoil,xprodptase, &
         xcostnpup,xmaxfinelitter,xmaxcwd,xnintercept,xnslope
    REAL(r_2), DIMENSION(mso)          :: xxkplab,xxkpsorb,xxkpocc


casabiome%ivt2( 1 )  = 3
casabiome%ivt2( 2 )  = 3
casabiome%ivt2( 3 )  = 3
casabiome%ivt2( 4 )  = 3
casabiome%ivt2( 5 )  = 2
casabiome%ivt2( 6 )  = 1
casabiome%ivt2( 7 )  = 1
casabiome%ivt2( 8 )  = 2
casabiome%ivt2( 9 )  = 1
casabiome%ivt2( 10 ) = 1
casabiome%ivt2( 11 ) = 0
casabiome%ivt2( 12 ) = 0
casabiome%ivt2( 13 ) = 0
casabiome%ivt2( 14 ) = 0
casabiome%ivt2( 15 ) = 0
casabiome%ivt2( 16 ) = 0
casabiome%ivt2( 17 ) = 0

! where ice,water,urban=0)/grass=1/shrub=2/woody=3

! Tranche 1 
!===============================================================
casabiome % kroot    (  1 ) =    5.500
casabiome %rootdepth (  1 ) =    1.500
casabiome % kuptake  (  1 ) =    2.000
casabiome % krootlen (  1 ) =   14.878
casabiome % kminN    (  1 ) =    2.000
casabiome % kuplabP  (  1 ) =    0.500
xfherbivore          (  1 ) =    0.068
leafage              (  1 ) =    3.000
woodage              (  1 ) =   70.000
frootage             (  1 ) =   18.000
metage               (  1 ) =    0.040
strage               (  1 ) =    0.230
cwdage               (  1 ) =    0.824
micage               (  1 ) =    0.137
slowage              (  1 ) =    5.000
passage              (  1 ) =  222.220
clabileage           (  1 ) =    0.200
slax                 (  1 ) =    0.007
      
casabiome % kroot    (  2 ) =    3.900
casabiome %rootdepth (  2 ) =    1.500
casabiome % kuptake  (  2 ) =    1.900
casabiome % krootlen (  2 ) =   14.386
casabiome % kminN    (  2 ) =    2.000
casabiome % kuplabP  (  2 ) =    0.500
xfherbivore          (  2 ) =    0.406
leafage              (  2 ) =    2.000
woodage              (  2 ) =   60.000
frootage             (  2 ) =   10.000
metage               (  2 ) =    0.040
strage               (  2 ) =    0.230
cwdage               (  2 ) =    0.824
micage               (  2 ) =    0.137
slowage              (  2 ) =    5.000
passage              (  2 ) =  222.220
clabileage           (  2 ) =    0.200
slax                 (  2 ) =    0.015
      
casabiome % kroot    (  3 ) =    5.500
casabiome %rootdepth (  3 ) =    1.500
casabiome % kuptake  (  3 ) =    2.000
casabiome % krootlen (  3 ) =   14.026
casabiome % kminN    (  3 ) =    2.000
casabiome % kuplabP  (  3 ) =    0.500
xfherbivore          (  3 ) =    0.068
leafage              (  3 ) =    0.500
woodage              (  3 ) =   80.000
frootage             (  3 ) =   10.000
metage               (  3 ) =    0.040
strage               (  3 ) =    0.230
cwdage               (  3 ) =    0.824
micage               (  3 ) =    0.137
slowage              (  3 ) =    5.000
passage              (  3 ) =  222.220
clabileage           (  3 ) =    0.200
slax                 (  3 ) =    0.023
      
casabiome % kroot    (  4 ) =    3.900
casabiome %rootdepth (  4 ) =    1.500
casabiome % kuptake  (  4 ) =    2.000
casabiome % krootlen (  4 ) =   18.947
casabiome % kminN    (  4 ) =    2.000
casabiome % kuplabP  (  4 ) =    0.500
xfherbivore          (  4 ) =    0.134
leafage              (  4 ) =    0.518
woodage              (  4 ) =   40.000
frootage             (  4 ) =   10.000
metage               (  4 ) =    0.040
strage               (  4 ) =    0.230
cwdage               (  4 ) =    0.824
micage               (  4 ) =    0.137
slowage              (  4 ) =    5.000
passage              (  4 ) =  222.220
clabileage           (  4 ) =    0.200
slax                 (  4 ) =    0.026
      
casabiome % kroot    (  5 ) =    2.000
casabiome %rootdepth (  5 ) =    0.500
casabiome % kuptake  (  5 ) =    1.800
casabiome % krootlen (  5 ) =   32.308
casabiome % kminN    (  5 ) =    2.000
casabiome % kuplabP  (  5 ) =    0.500
xfherbivore          (  5 ) =    0.022
leafage              (  5 ) =    1.440
woodage              (  5 ) =   40.000
frootage             (  5 ) =    5.000
metage               (  5 ) =    0.040
strage               (  5 ) =    0.230
cwdage               (  5 ) =    0.824
micage               (  5 ) =    0.137
slowage              (  5 ) =    5.000
passage              (  5 ) =  222.220
clabileage           (  5 ) =    0.200
slax                 (  5 ) =    0.010
      
casabiome % kroot    (  6 ) =    5.500
casabiome %rootdepth (  6 ) =    0.500
casabiome % kuptake  (  6 ) =    2.000
casabiome % krootlen (  6 ) =   84.000
casabiome % kminN    (  6 ) =    2.000
casabiome % kuplabP  (  6 ) =    0.500
xfherbivore          (  6 ) =    0.109
leafage              (  6 ) =    0.500
woodage              (  6 ) =    1.000
frootage             (  6 ) =    3.000
metage               (  6 ) =    0.040
strage               (  6 ) =    0.230
cwdage               (  6 ) =    0.824
micage               (  6 ) =    0.137
slowage              (  6 ) =    5.000
passage              (  6 ) =  222.220
clabileage           (  6 ) =    0.200
slax                 (  6 ) =    0.030
      
casabiome % kroot    (  7 ) =    5.500
casabiome %rootdepth (  7 ) =    0.500
casabiome % kuptake  (  7 ) =    2.000
casabiome % krootlen (  7 ) =   84.000
casabiome % kminN    (  7 ) =    2.000
casabiome % kuplabP  (  7 ) =    0.500
xfherbivore          (  7 ) =    0.109
leafage              (  7 ) =    0.500
woodage              (  7 ) =    1.000
frootage             (  7 ) =    3.000
metage               (  7 ) =    0.040
strage               (  7 ) =    0.230
cwdage               (  7 ) =    0.824
micage               (  7 ) =    0.137
slowage              (  7 ) =    5.000
passage              (  7 ) =  222.220
clabileage           (  7 ) =    0.200
slax                 (  7 ) =    0.022
      
casabiome % kroot    (  8 ) =    5.500
casabiome %rootdepth (  8 ) =    0.500
casabiome % kuptake  (  8 ) =    2.000
casabiome % krootlen (  8 ) =   84.000
casabiome % kminN    (  8 ) =    2.000
casabiome % kuplabP  (  8 ) =    0.500
xfherbivore          (  8 ) =    0.109
leafage              (  8 ) =    0.750
woodage              (  8 ) =    1.000
frootage             (  8 ) =    3.000
metage               (  8 ) =    0.040
strage               (  8 ) =    0.230
cwdage               (  8 ) =    0.824
micage               (  8 ) =    0.137
slowage              (  8 ) =    5.000
passage              (  8 ) =  222.220
clabileage           (  8 ) =    0.200
slax                 (  8 ) =    0.027
      
casabiome % kroot    (  9 ) =    5.500
casabiome %rootdepth (  9 ) =    0.500
casabiome % kuptake  (  9 ) =    1.600
casabiome % krootlen (  9 ) =  120.500
casabiome % kminN    (  9 ) =    2.000
casabiome % kuplabP  (  9 ) =    0.500
xfherbivore          (  9 ) =    0.140
leafage              (  9 ) =    0.370
woodage              (  9 ) =    1.000
frootage             (  9 ) =    0.884
metage               (  9 ) =    0.040
strage               (  9 ) =    0.230
cwdage               (  9 ) =    0.824
micage               (  9 ) =    0.137
slowage              (  9 ) =    5.000
passage              (  9 ) =  222.220
clabileage           (  9 ) =    0.200
slax                 (  9 ) =    0.030
      
casabiome % kroot    ( 10 ) =    5.500
casabiome %rootdepth ( 10 ) =    0.500
casabiome % kuptake  ( 10 ) =    1.600
casabiome % krootlen ( 10 ) =  120.500
casabiome % kminN    ( 10 ) =    2.000
casabiome % kuplabP  ( 10 ) =    0.500
xfherbivore          ( 10 ) =    0.140
leafage              ( 10 ) =    0.370
woodage              ( 10 ) =    1.000
frootage             ( 10 ) =    0.884
metage               ( 10 ) =    0.040
strage               ( 10 ) =    0.230
cwdage               ( 10 ) =    0.824
micage               ( 10 ) =    0.137
slowage              ( 10 ) =    5.000
passage              ( 10 ) =  222.220
clabileage           ( 10 ) =    0.200
slax                 ( 10 ) =    0.022
      
casabiome % kroot    ( 11 ) =    5.500
casabiome %rootdepth ( 11 ) =    0.500
casabiome % kuptake  ( 11 ) =    1.600
casabiome % krootlen ( 11 ) =    0.000
casabiome % kminN    ( 11 ) =    2.000
casabiome % kuplabP  ( 11 ) =    0.500
xfherbivore          ( 11 ) =    0.000
leafage              ( 11 ) =    1.000
woodage              ( 11 ) =    1.000
frootage             ( 11 ) =    1.000
metage               ( 11 ) =    0.040
strage               ( 11 ) =    0.230
cwdage               ( 11 ) =    0.824
micage               ( 11 ) =    0.137
slowage              ( 11 ) =    5.000
passage              ( 11 ) =  222.220
clabileage           ( 11 ) =    0.200
slax                 ( 11 ) =    0.020
      
casabiome % kroot    ( 12 ) =    5.500
casabiome %rootdepth ( 12 ) =    0.500
casabiome % kuptake  ( 12 ) =    1.800
casabiome % krootlen ( 12 ) =    0.000
casabiome % kminN    ( 12 ) =    2.000
casabiome % kuplabP  ( 12 ) =    0.500
xfherbivore          ( 12 ) =    0.000
leafage              ( 12 ) =    1.000
woodage              ( 12 ) =    1.000
frootage             ( 12 ) =    1.000
metage               ( 12 ) =    0.040
strage               ( 12 ) =    0.230
cwdage               ( 12 ) =    0.824
micage               ( 12 ) =    0.137
slowage              ( 12 ) =    5.000
passage              ( 12 ) =  222.220
clabileage           ( 12 ) =    0.200
slax                 ( 12 ) =    0.020
      
casabiome % kroot    ( 13 ) =    5.500
casabiome %rootdepth ( 13 ) =    0.500
casabiome % kuptake  ( 13 ) =    1.800
casabiome % krootlen ( 13 ) =    0.000
casabiome % kminN    ( 13 ) =    2.000
casabiome % kuplabP  ( 13 ) =    0.500
xfherbivore          ( 13 ) =    0.000
leafage              ( 13 ) =    1.000
woodage              ( 13 ) =    1.000
frootage             ( 13 ) =    1.000
metage               ( 13 ) =    0.040
strage               ( 13 ) =    0.230
cwdage               ( 13 ) =    0.824
micage               ( 13 ) =    0.137
slowage              ( 13 ) =    5.000
passage              ( 13 ) =  222.220
clabileage           ( 13 ) =    0.200
slax                 ( 13 ) =    0.020
      
casabiome % kroot    ( 14 ) =    2.000
casabiome %rootdepth ( 14 ) =    0.500
casabiome % kuptake  ( 14 ) =    1.800
casabiome % krootlen ( 14 ) =   30.769
casabiome % kminN    ( 14 ) =    2.000
casabiome % kuplabP  ( 14 ) =    0.500
xfherbivore          ( 14 ) =    0.010
leafage              ( 14 ) =    0.461
woodage              ( 14 ) =    5.000
frootage             ( 14 ) =    4.000
metage               ( 14 ) =    0.040
strage               ( 14 ) =    0.230
cwdage               ( 14 ) =    0.824
micage               ( 14 ) =    0.137
slowage              ( 14 ) =    5.000
passage              ( 14 ) =  222.220
clabileage           ( 14 ) =    0.200
slax                 ( 14 ) =    0.024
      
casabiome % kroot    ( 15 ) =    2.000
casabiome %rootdepth ( 15 ) =    0.500
casabiome % kuptake  ( 15 ) =    1.800
casabiome % krootlen ( 15 ) =    0.000
casabiome % kminN    ( 15 ) =    2.000
casabiome % kuplabP  ( 15 ) =    0.500
xfherbivore          ( 15 ) =    0.000
leafage              ( 15 ) =    1.000
woodage              ( 15 ) =    1.000
frootage             ( 15 ) =    1.000
metage               ( 15 ) =    0.040
strage               ( 15 ) =    0.230
cwdage               ( 15 ) =    0.824
micage               ( 15 ) =    0.137
slowage              ( 15 ) =    5.000
passage              ( 15 ) =  222.220
clabileage           ( 15 ) =    0.200
slax                 ( 15 ) =    0.020
      
casabiome % kroot    ( 16 ) =    5.500
casabiome %rootdepth ( 16 ) =    1.500
casabiome % kuptake  ( 16 ) =    1.800
casabiome % krootlen ( 16 ) =    0.000
casabiome % kminN    ( 16 ) =    2.000
casabiome % kuplabP  ( 16 ) =    0.500
xfherbivore          ( 16 ) =    0.000
leafage              ( 16 ) =    1.000
woodage              ( 16 ) =    1.000
frootage             ( 16 ) =    1.000
metage               ( 16 ) =    0.040
strage               ( 16 ) =    0.230
cwdage               ( 16 ) =    0.824
micage               ( 16 ) =    0.137
slowage              ( 16 ) =    5.000
passage              ( 16 ) =  222.220
clabileage           ( 16 ) =    0.200
slax                 ( 16 ) =    0.020
      
casabiome % kroot    ( 17 ) =    5.500
casabiome %rootdepth ( 17 ) =    0.500
casabiome % kuptake  ( 17 ) =    1.800
casabiome % krootlen ( 17 ) =    0.000
casabiome % kminN    ( 17 ) =    2.000
casabiome % kuplabP  ( 17 ) =    0.500
xfherbivore          ( 17 ) =    0.000
leafage              ( 17 ) =    1.000
woodage              ( 17 ) =    1.000
frootage             ( 17 ) =    1.000
metage               ( 17 ) =    0.040
strage               ( 17 ) =    0.230
cwdage               ( 17 ) =    0.824
micage               ( 17 ) =    0.137
slowage              ( 17 ) =    5.000
passage              ( 17 ) =  222.220
clabileage           ( 17 ) =    0.200
slax                 ( 17 ) =    0.020

! Tranche *2* should be 2a
!===============================================================
casabiome%fracnpptoP (  1,  1 ) =    0.250
casabiome%fracnpptoP (  1,  2 ) =    0.400
casabiome%fracnpptoP (  1,  3 ) =    0.350
casabiome%rmplant    (  1,  1 ) =    0.100
casabiome%rmplant    (  1,  2 ) =    2.000
casabiome%rmplant    (  1,  3 ) =   10.000
      
casabiome%fracnpptoP (  2,  1 ) =    0.250
casabiome%fracnpptoP (  2,  2 ) =    0.400
casabiome%fracnpptoP (  2,  3 ) =    0.350
casabiome%rmplant    (  2,  1 ) =    0.100
casabiome%rmplant    (  2,  2 ) =    2.000
casabiome%rmplant    (  2,  3 ) =   10.000
      
casabiome%fracnpptoP (  3,  1 ) =    0.250
casabiome%fracnpptoP (  3,  2 ) =    0.400
casabiome%fracnpptoP (  3,  3 ) =    0.350
casabiome%rmplant    (  3,  1 ) =    0.100
casabiome%rmplant    (  3,  2 ) =    2.000
casabiome%rmplant    (  3,  3 ) =   10.000
      
casabiome%fracnpptoP (  4,  1 ) =    0.250
casabiome%fracnpptoP (  4,  2 ) =    0.400
casabiome%fracnpptoP (  4,  3 ) =    0.350
casabiome%rmplant    (  4,  1 ) =    0.100
casabiome%rmplant    (  4,  2 ) =    2.000
casabiome%rmplant    (  4,  3 ) =   10.000
      
casabiome%fracnpptoP (  5,  1 ) =    0.250
casabiome%fracnpptoP (  5,  2 ) =    0.400
casabiome%fracnpptoP (  5,  3 ) =    0.350
casabiome%rmplant    (  5,  1 ) =    0.100
casabiome%rmplant    (  5,  2 ) =    2.000
casabiome%rmplant    (  5,  3 ) =   10.000
      
casabiome%fracnpptoP (  6,  1 ) =    0.250
casabiome%fracnpptoP (  6,  2 ) =    0.400
casabiome%fracnpptoP (  6,  3 ) =    0.350
casabiome%rmplant    (  6,  1 ) =    0.100
casabiome%rmplant    (  6,  2 ) =    2.000
casabiome%rmplant    (  6,  3 ) =   10.000
      
casabiome%fracnpptoP (  7,  1 ) =    0.250
casabiome%fracnpptoP (  7,  2 ) =    0.400
casabiome%fracnpptoP (  7,  3 ) =    0.350
casabiome%rmplant    (  7,  1 ) =    0.100
casabiome%rmplant    (  7,  2 ) =    2.000
casabiome%rmplant    (  7,  3 ) =   10.000
      
casabiome%fracnpptoP (  8,  1 ) =    0.250
casabiome%fracnpptoP (  8,  2 ) =    0.400
casabiome%fracnpptoP (  8,  3 ) =    0.350
casabiome%rmplant    (  8,  1 ) =    0.100
casabiome%rmplant    (  8,  2 ) =    2.000
casabiome%rmplant    (  8,  3 ) =   10.000
      
casabiome%fracnpptoP (  9,  1 ) =    0.250
casabiome%fracnpptoP (  9,  2 ) =    0.400
casabiome%fracnpptoP (  9,  3 ) =    0.350
casabiome%rmplant    (  9,  1 ) =    0.100
casabiome%rmplant    (  9,  2 ) =    2.000
casabiome%rmplant    (  9,  3 ) =   10.000
      
casabiome%fracnpptoP ( 10,  1 ) =    0.250
casabiome%fracnpptoP ( 10,  2 ) =    0.400
casabiome%fracnpptoP ( 10,  3 ) =    0.350
casabiome%rmplant    ( 10,  1 ) =    0.100
casabiome%rmplant    ( 10,  2 ) =    2.000
casabiome%rmplant    ( 10,  3 ) =   10.000
      
casabiome%fracnpptoP ( 11,  1 ) =    0.250
casabiome%fracnpptoP ( 11,  2 ) =    0.400
casabiome%fracnpptoP ( 11,  3 ) =    0.350
casabiome%rmplant    ( 11,  1 ) =    0.100
casabiome%rmplant    ( 11,  2 ) =    2.000
casabiome%rmplant    ( 11,  3 ) =   10.000
      
casabiome%fracnpptoP ( 12,  1 ) =    0.250
casabiome%fracnpptoP ( 12,  2 ) =    0.400
casabiome%fracnpptoP ( 12,  3 ) =    0.350
casabiome%rmplant    ( 12,  1 ) =    0.100
casabiome%rmplant    ( 12,  2 ) =    2.000
casabiome%rmplant    ( 12,  3 ) =   10.000
      
casabiome%fracnpptoP ( 13,  1 ) =    0.250
casabiome%fracnpptoP ( 13,  2 ) =    0.400
casabiome%fracnpptoP ( 13,  3 ) =    0.350
casabiome%rmplant    ( 13,  1 ) =    0.100
casabiome%rmplant    ( 13,  2 ) =    2.000
casabiome%rmplant    ( 13,  3 ) =   10.000
      
casabiome%fracnpptoP ( 14,  1 ) =    0.250
casabiome%fracnpptoP ( 14,  2 ) =    0.400
casabiome%fracnpptoP ( 14,  3 ) =    0.350
casabiome%rmplant    ( 14,  1 ) =    0.100
casabiome%rmplant    ( 14,  2 ) =    2.000
casabiome%rmplant    ( 14,  3 ) =   10.000
      
casabiome%fracnpptoP ( 15,  1 ) =    0.250
casabiome%fracnpptoP ( 15,  2 ) =    0.400
casabiome%fracnpptoP ( 15,  3 ) =    0.350
casabiome%rmplant    ( 15,  1 ) =    0.100
casabiome%rmplant    ( 15,  2 ) =    2.000
casabiome%rmplant    ( 15,  3 ) =   10.000
      
casabiome%fracnpptoP ( 16,  1 ) =    0.250
casabiome%fracnpptoP ( 16,  2 ) =    0.400
casabiome%fracnpptoP ( 16,  3 ) =    0.350
casabiome%rmplant    ( 16,  1 ) =    0.100
casabiome%rmplant    ( 16,  2 ) =    2.000
casabiome%rmplant    ( 16,  3 ) =   10.000
      
casabiome%fracnpptoP ( 17,  1 ) =    0.250
casabiome%fracnpptoP ( 17,  2 ) =    0.400
casabiome%fracnpptoP ( 17,  3 ) =    0.350
casabiome%rmplant    ( 17,  1 ) =    0.100
casabiome%rmplant    ( 17,  2 ) =    2.000
casabiome%rmplant    ( 17,  3 ) =   10.000
    
! Tranche 2b
!===============================================================
ratioCNplant                (  1,  1 ) =   49.800
ratioCNplant                (  1,  2 ) =  238.100
ratioCNplant                (  1,  3 ) =   73.700
casabiome % ftransNPtoL     (  1,  1 ) =    0.500
casabiome % ftransNPtoL     (  1,  2 ) =    0.950
casabiome % ftransNPtoL     (  1,  3 ) =    0.900
casabiome % fracligninplant (  1,  1 ) =    0.250
casabiome % fracligninplant (  1,  2 ) =    0.400
casabiome % fracligninplant (  1,  3 ) =    0.250
ratioCNsoil                 (  1,  1 ) =    8.000
ratioCNsoil                 (  1,  2 ) =   16.100
ratioCNsoil                 (  1,  3 ) =   16.100
ratioCNsoilmin              (  1,  1 ) =    3.000
ratioCNsoilmin              (  1,  2 ) =   12.000
ratioCNsoilmin              (  1,  3 ) =    7.000
ratioCNsoilmax              (  1,  1 ) =   15.000
ratioCNsoilmax              (  1,  2 ) =   30.000
ratioCNsoilmax              (  1,  3 ) =   15.000
casabiome % glaimax         (  1 )     =    7.000
casabiome % glaimin         (  1 )     =    1.000
      
ratioCNplant                (  2,  1 ) =   23.100
ratioCNplant                (  2,  2 ) =  134.900
ratioCNplant                (  2,  3 ) =   61.200
casabiome % ftransNPtoL     (  2,  1 ) =    0.500
casabiome % ftransNPtoL     (  2,  2 ) =    0.950
casabiome % ftransNPtoL     (  2,  3 ) =    0.900
casabiome % fracligninplant (  2,  1 ) =    0.200
casabiome % fracligninplant (  2,  2 ) =    0.400
casabiome % fracligninplant (  2,  3 ) =    0.200
ratioCNsoil                 (  2,  1 ) =    8.000
ratioCNsoil                 (  2,  2 ) =   12.800
ratioCNsoil                 (  2,  3 ) =   12.800
ratioCNsoilmin              (  2,  1 ) =    3.000
ratioCNsoilmin              (  2,  2 ) =   12.000
ratioCNsoilmin              (  2,  3 ) =    7.000
ratioCNsoilmax              (  2,  1 ) =   15.000
ratioCNsoilmax              (  2,  2 ) =   30.000
ratioCNsoilmax              (  2,  3 ) =   15.000
casabiome % glaimax         (  2 ) =    7.000
casabiome % glaimin         (  2 ) =    1.000
      
ratioCNplant                (  3,  1 ) =   59.300
ratioCNplant                (  3,  2 ) =  243.800
ratioCNplant                (  3,  3 ) =   75.000
casabiome % ftransNPtoL     (  3,  1 ) =    0.500
casabiome % ftransNPtoL     (  3,  2 ) =    0.950
casabiome % ftransNPtoL     (  3,  3 ) =    0.900
casabiome % fracligninplant (  3,  1 ) =    0.200
casabiome % fracligninplant (  3,  2 ) =    0.400
casabiome % fracligninplant (  3,  3 ) =    0.200
ratioCNsoil                 (  3,  1 ) =    8.000
ratioCNsoil                 (  3,  2 ) =   24.800
ratioCNsoil                 (  3,  3 ) =   24.800
ratioCNsoilmin              (  3,  1 ) =    3.000
ratioCNsoilmin              (  3,  2 ) =   12.000
ratioCNsoilmin              (  3,  3 ) =    7.000
ratioCNsoilmax              (  3,  1 ) =   15.000
ratioCNsoilmax              (  3,  2 ) =   30.000
ratioCNsoilmax              (  3,  3 ) =   15.000
casabiome % glaimax         (  3 ) =    7.000
casabiome % glaimin         (  3 ) =    0.500
      
ratioCNplant                (  4,  1 ) =   31.400
ratioCNplant                (  4,  2 ) =  156.200
ratioCNplant                (  4,  3 ) =   63.200
casabiome % ftransNPtoL     (  4,  1 ) =    0.500
casabiome % ftransNPtoL     (  4,  2 ) =    0.950
casabiome % ftransNPtoL     (  4,  3 ) =    0.900
casabiome % fracligninplant (  4,  1 ) =    0.200
casabiome % fracligninplant (  4,  2 ) =    0.400
casabiome % fracligninplant (  4,  3 ) =    0.200
ratioCNsoil                 (  4,  1 ) =    8.000
ratioCNsoil                 (  4,  2 ) =   30.000
ratioCNsoil                 (  4,  3 ) =   30.000
ratioCNsoilmin              (  4,  1 ) =    3.000
ratioCNsoilmin              (  4,  2 ) =   12.000
ratioCNsoilmin              (  4,  3 ) =    7.000
ratioCNsoilmax              (  4,  1 ) =   15.000
ratioCNsoilmax              (  4,  2 ) =   30.000
ratioCNsoilmax              (  4,  3 ) =   15.000
casabiome % glaimax         (  4 ) =    7.000
casabiome % glaimin         (  4 ) =    0.500
      
ratioCNplant                (  5,  1 ) =   37.600
ratioCNplant                (  5,  2 ) =  142.100
ratioCNplant                (  5,  3 ) =   67.100
casabiome % ftransNPtoL     (  5,  1 ) =    0.500
casabiome % ftransNPtoL     (  5,  2 ) =    0.950
casabiome % ftransNPtoL     (  5,  3 ) =    0.900
casabiome % fracligninplant (  5,  1 ) =    0.200
casabiome % fracligninplant (  5,  2 ) =    0.400
casabiome % fracligninplant (  5,  3 ) =    0.200
ratioCNsoil                 (  5,  1 ) =    8.000
ratioCNsoil                 (  5,  2 ) =   19.300
ratioCNsoil                 (  5,  3 ) =   19.300
ratioCNsoilmin              (  5,  1 ) =    3.000
ratioCNsoilmin              (  5,  2 ) =   12.000
ratioCNsoilmin              (  5,  3 ) =    7.000
ratioCNsoilmax              (  5,  1 ) =   15.000
ratioCNsoilmax              (  5,  2 ) =   30.000
ratioCNsoilmax              (  5,  3 ) =   15.000
casabiome % glaimax         (  5 ) =    3.000
casabiome % glaimin         (  5 ) =    0.100
      
ratioCNplant                (  6,  1 ) =   34.800
ratioCNplant                (  6,  2 ) =  150.000
ratioCNplant                (  6,  3 ) =   64.500
casabiome % ftransNPtoL     (  6,  1 ) =    0.500
casabiome % ftransNPtoL     (  6,  2 ) =    0.950
casabiome % ftransNPtoL     (  6,  3 ) =    0.900
casabiome % fracligninplant (  6,  1 ) =    0.100
casabiome % fracligninplant (  6,  2 ) =    0.400
casabiome % fracligninplant (  6,  3 ) =    0.100
ratioCNsoil                 (  6,  1 ) =    8.000
ratioCNsoil                 (  6,  2 ) =   13.100
ratioCNsoil                 (  6,  3 ) =   13.100
ratioCNsoilmin              (  6,  1 ) =    3.000
ratioCNsoilmin              (  6,  2 ) =   12.000
ratioCNsoilmin              (  6,  3 ) =    7.000
ratioCNsoilmax              (  6,  1 ) =   15.000
ratioCNsoilmax              (  6,  2 ) =   30.000
ratioCNsoilmax              (  6,  3 ) =   15.000
casabiome % glaimax         (  6 ) =    3.000
casabiome % glaimin         (  6 ) =    0.100
      
ratioCNplant                (  7,  1 ) =   44.000
ratioCNplant                (  7,  2 ) =  150.000
ratioCNplant                (  7,  3 ) =   62.700
casabiome % ftransNPtoL     (  7,  1 ) =    0.500
casabiome % ftransNPtoL     (  7,  2 ) =    0.950
casabiome % ftransNPtoL     (  7,  3 ) =    0.900
casabiome % fracligninplant (  7,  1 ) =    0.100
casabiome % fracligninplant (  7,  2 ) =    0.400
casabiome % fracligninplant (  7,  3 ) =    0.100
ratioCNsoil                 (  7,  1 ) =    8.000
ratioCNsoil                 (  7,  2 ) =   13.100
ratioCNsoil                 (  7,  3 ) =   13.100
ratioCNsoilmin              (  7,  1 ) =    3.000
ratioCNsoilmin              (  7,  2 ) =   12.000
ratioCNsoilmin              (  7,  3 ) =    7.000
ratioCNsoilmax              (  7,  1 ) =   15.000
ratioCNsoilmax              (  7,  2 ) =   30.000
ratioCNsoilmax              (  7,  3 ) =   15.000
casabiome % glaimax         (  7 ) =    3.000
casabiome % glaimin         (  7 ) =    0.100
      
ratioCNplant                (  8,  1 ) =   49.200
ratioCNplant                (  8,  2 ) =  147.300
ratioCNplant                (  8,  3 ) =   69.000
casabiome % ftransNPtoL     (  8,  1 ) =    0.500
casabiome % ftransNPtoL     (  8,  2 ) =    0.950
casabiome % ftransNPtoL     (  8,  3 ) =    0.900
casabiome % fracligninplant (  8,  1 ) =    0.100
casabiome % fracligninplant (  8,  2 ) =    0.400
casabiome % fracligninplant (  8,  3 ) =    0.100
ratioCNsoil                 (  8,  1 ) =    8.000
ratioCNsoil                 (  8,  2 ) =   13.100
ratioCNsoil                 (  8,  3 ) =   13.100
ratioCNsoilmin              (  8,  1 ) =    3.000
ratioCNsoilmin              (  8,  2 ) =   12.000
ratioCNsoilmin              (  8,  3 ) =    7.000
ratioCNsoilmax              (  8,  1 ) =   15.000
ratioCNsoilmax              (  8,  2 ) =   30.000
ratioCNsoilmax              (  8,  3 ) =   15.000
casabiome % glaimax         (  8 ) =    3.000
casabiome % glaimin         (  8 ) =    0.100
      
ratioCNplant                (  9,  1 ) =   21.600
ratioCNplant                (  9,  2 ) =  150.000
ratioCNplant                (  9,  3 ) =   60.700
casabiome % ftransNPtoL     (  9,  1 ) =    0.500
casabiome % ftransNPtoL     (  9,  2 ) =    0.950
casabiome % ftransNPtoL     (  9,  3 ) =    0.900
casabiome % fracligninplant (  9,  1 ) =    0.100
casabiome % fracligninplant (  9,  2 ) =    0.400
casabiome % fracligninplant (  9,  3 ) =    0.100
ratioCNsoil                 (  9,  1 ) =    8.000
ratioCNsoil                 (  9,  2 ) =   13.200
ratioCNsoil                 (  9,  3 ) =   13.200
ratioCNsoilmin              (  9,  1 ) =    3.000
ratioCNsoilmin              (  9,  2 ) =   12.000
ratioCNsoilmin              (  9,  3 ) =    7.000
ratioCNsoilmax              (  9,  1 ) =   15.000
ratioCNsoilmax              (  9,  2 ) =   30.000
ratioCNsoilmax              (  9,  3 ) =   15.000
casabiome % glaimax         (  9 ) =    6.000
casabiome % glaimin         (  9 ) =    0.100
      
ratioCNplant                ( 10,  1 ) =   25.000
ratioCNplant                ( 10,  2 ) =  125.000
ratioCNplant                ( 10,  3 ) =   71.000
casabiome % ftransNPtoL     ( 10,  1 ) =    0.500
casabiome % ftransNPtoL     ( 10,  2 ) =    0.950
casabiome % ftransNPtoL     ( 10,  3 ) =    0.900
casabiome % fracligninplant ( 10,  1 ) =    0.100
casabiome % fracligninplant ( 10,  2 ) =    0.400
casabiome % fracligninplant ( 10,  3 ) =    0.100
ratioCNsoil                 ( 10,  1 ) =    8.000
ratioCNsoil                 ( 10,  2 ) =   13.200
ratioCNsoil                 ( 10,  3 ) =   13.200
ratioCNsoilmin              ( 10,  1 ) =    3.000
ratioCNsoilmin              ( 10,  2 ) =   12.000
ratioCNsoilmin              ( 10,  3 ) =    7.000
ratioCNsoilmax              ( 10,  1 ) =   15.000
ratioCNsoilmax              ( 10,  2 ) =   30.000
ratioCNsoilmax              ( 10,  3 ) =   15.000
casabiome % glaimax         ( 10 ) =    6.000
casabiome % glaimin         ( 10 ) =    0.100
      
ratioCNplant                ( 11,  1 ) =   30.000
ratioCNplant                ( 11,  2 ) =  150.000
ratioCNplant                ( 11,  3 ) =   71.000
casabiome % ftransNPtoL     ( 11,  1 ) =    0.500
casabiome % ftransNPtoL     ( 11,  2 ) =    0.950
casabiome % ftransNPtoL     ( 11,  3 ) =    0.900
casabiome % fracligninplant ( 11,  1 ) =    0.150
casabiome % fracligninplant ( 11,  2 ) =    0.400
casabiome % fracligninplant ( 11,  3 ) =    0.150
ratioCNsoil                 ( 11,  1 ) =    8.000
ratioCNsoil                 ( 11,  2 ) =   13.100
ratioCNsoil                 ( 11,  3 ) =   13.100
ratioCNsoilmin              ( 11,  1 ) =    3.000
ratioCNsoilmin              ( 11,  2 ) =   12.000
ratioCNsoilmin              ( 11,  3 ) =    7.000
ratioCNsoilmax              ( 11,  1 ) =   15.000
ratioCNsoilmax              ( 11,  2 ) =   30.000
ratioCNsoilmax              ( 11,  3 ) =   15.000
casabiome % glaimax         ( 11 ) =    5.000
casabiome % glaimin         ( 11 ) =    0.050
      
ratioCNplant                ( 12,  1 ) =   30.000
ratioCNplant                ( 12,  2 ) =  150.000
ratioCNplant                ( 12,  3 ) =   71.000
casabiome % ftransNPtoL     ( 12,  1 ) =    0.500
casabiome % ftransNPtoL     ( 12,  2 ) =    0.950
casabiome % ftransNPtoL     ( 12,  3 ) =    0.900
casabiome % fracligninplant ( 12,  1 ) =    0.150
casabiome % fracligninplant ( 12,  2 ) =    0.400
casabiome % fracligninplant ( 12,  3 ) =    0.150
ratioCNsoil                 ( 12,  1 ) =    8.000
ratioCNsoil                 ( 12,  2 ) =   13.100
ratioCNsoil                 ( 12,  3 ) =   13.100
ratioCNsoilmin              ( 12,  1 ) =    3.000
ratioCNsoilmin              ( 12,  2 ) =   12.000
ratioCNsoilmin              ( 12,  3 ) =    7.000
ratioCNsoilmax              ( 12,  1 ) =   15.000
ratioCNsoilmax              ( 12,  2 ) =   30.000
ratioCNsoilmax              ( 12,  3 ) =   15.000
casabiome % glaimax         ( 12 ) =    5.000
casabiome % glaimin         ( 12 ) =    0.050
      
ratioCNplant                ( 13,  1 ) =   30.000
ratioCNplant                ( 13,  2 ) =  150.000
ratioCNplant                ( 13,  3 ) =   71.000
casabiome % ftransNPtoL     ( 13,  1 ) =    0.500
casabiome % ftransNPtoL     ( 13,  2 ) =    0.950
casabiome % ftransNPtoL     ( 13,  3 ) =    0.900
casabiome % fracligninplant ( 13,  1 ) =    0.150
casabiome % fracligninplant ( 13,  2 ) =    0.400
casabiome % fracligninplant ( 13,  3 ) =    0.150
ratioCNsoil                 ( 13,  1 ) =    8.000
ratioCNsoil                 ( 13,  2 ) =   13.100
ratioCNsoil                 ( 13,  3 ) =   13.100
ratioCNsoilmin              ( 13,  1 ) =    3.000
ratioCNsoilmin              ( 13,  2 ) =   12.000
ratioCNsoilmin              ( 13,  3 ) =    7.000
ratioCNsoilmax              ( 13,  1 ) =   15.000
ratioCNsoilmax              ( 13,  2 ) =   30.000
ratioCNsoilmax              ( 13,  3 ) =   15.000
casabiome % glaimax         ( 13 ) =    5.000
casabiome % glaimin         ( 13 ) =    0.050
      
ratioCNplant                ( 14,  1 ) =   50.000
ratioCNplant                ( 14,  2 ) =  150.000
ratioCNplant                ( 14,  3 ) =   71.000
casabiome % ftransNPtoL     ( 14,  1 ) =    0.500
casabiome % ftransNPtoL     ( 14,  2 ) =    0.950
casabiome % ftransNPtoL     ( 14,  3 ) =    0.900
casabiome % fracligninplant ( 14,  1 ) =    0.150
casabiome % fracligninplant ( 14,  2 ) =    0.400
casabiome % fracligninplant ( 14,  3 ) =    0.150
ratioCNsoil                 ( 14,  1 ) =    8.000
ratioCNsoil                 ( 14,  2 ) =   26.800
ratioCNsoil                 ( 14,  3 ) =   26.800
ratioCNsoilmin              ( 14,  1 ) =    3.000
ratioCNsoilmin              ( 14,  2 ) =   12.000
ratioCNsoilmin              ( 14,  3 ) =    7.000
ratioCNsoilmax              ( 14,  1 ) =   15.000
ratioCNsoilmax              ( 14,  2 ) =   30.000
ratioCNsoilmax              ( 14,  3 ) =   15.000
casabiome % glaimax         ( 14 ) =    1.000
casabiome % glaimin         ( 14 ) =    0.050
      
ratioCNplant                ( 15,  1 ) =   40.000
ratioCNplant                ( 15,  2 ) =  150.000
ratioCNplant                ( 15,  3 ) =   71.000
casabiome % ftransNPtoL     ( 15,  1 ) =    0.500
casabiome % ftransNPtoL     ( 15,  2 ) =    0.950
casabiome % ftransNPtoL     ( 15,  3 ) =    0.900
casabiome % fracligninplant ( 15,  1 ) =    0.150
casabiome % fracligninplant ( 15,  2 ) =    0.400
casabiome % fracligninplant ( 15,  3 ) =    0.150
ratioCNsoil                 ( 15,  1 ) =    8.000
ratioCNsoil                 ( 15,  2 ) =   20.000
ratioCNsoil                 ( 15,  3 ) =   20.000
ratioCNsoilmin              ( 15,  1 ) =    3.000
ratioCNsoilmin              ( 15,  2 ) =   12.000
ratioCNsoilmin              ( 15,  3 ) =    7.000
ratioCNsoilmax              ( 15,  1 ) =   15.000
ratioCNsoilmax              ( 15,  2 ) =   30.000
ratioCNsoilmax              ( 15,  3 ) =   15.000
casabiome % glaimax         ( 15 ) =    6.000
casabiome % glaimin         ( 15 ) =    0.050
      
ratioCNplant                ( 16,  1 ) =   40.000
ratioCNplant                ( 16,  2 ) =  135.000
ratioCNplant                ( 16,  3 ) =   71.000
casabiome % ftransNPtoL     ( 16,  1 ) =    0.500
casabiome % ftransNPtoL     ( 16,  2 ) =    0.950
casabiome % ftransNPtoL     ( 16,  3 ) =    0.900
casabiome % fracligninplant ( 16,  1 ) =    0.250
casabiome % fracligninplant ( 16,  2 ) =    0.400
casabiome % fracligninplant ( 16,  3 ) =    0.250
ratioCNsoil                 ( 16,  1 ) =    8.000
ratioCNsoil                 ( 16,  2 ) =   20.000
ratioCNsoil                 ( 16,  3 ) =   20.000
ratioCNsoilmin              ( 16,  1 ) =    3.000
ratioCNsoilmin              ( 16,  2 ) =   12.000
ratioCNsoilmin              ( 16,  3 ) =    7.000
ratioCNsoilmax              ( 16,  1 ) =   15.000
ratioCNsoilmax              ( 16,  2 ) =   30.000
ratioCNsoilmax              ( 16,  3 ) =   15.000
casabiome % glaimax         ( 16 ) =    1.000
casabiome % glaimin         ( 16 ) =    0.050
      
ratioCNplant                ( 17,  1 ) =   40.000
ratioCNplant                ( 17,  2 ) =  150.000
ratioCNplant                ( 17,  3 ) =   71.000
casabiome % ftransNPtoL     ( 17,  1 ) =    0.500
casabiome % ftransNPtoL     ( 17,  2 ) =    0.950
casabiome % ftransNPtoL     ( 17,  3 ) =    0.900
casabiome % fracligninplant ( 17,  1 ) =    0.100
casabiome % fracligninplant ( 17,  2 ) =    0.400
casabiome % fracligninplant ( 17,  3 ) =    0.100
ratioCNsoil                 ( 17,  1 ) =    8.000
ratioCNsoil                 ( 17,  2 ) =   20.000
ratioCNsoil                 ( 17,  3 ) =   20.000
ratioCNsoilmin              ( 17,  1 ) =    3.000
ratioCNsoilmin              ( 17,  2 ) =   12.000
ratioCNsoilmin              ( 17,  3 ) =    7.000
ratioCNsoilmax              ( 17,  1 ) =   15.000
ratioCNsoilmax              ( 17,  2 ) =   30.000
ratioCNsoilmax              ( 17,  3 ) =   15.000
casabiome % glaimax         ( 17 ) =    0.000
casabiome % glaimin         ( 17 ) =    0.000
      
 ! Tranche 3 
 !===============================================================
cleaf  (  1 ) =    5.500
cwood  (  1 ) =    5.500
cfroot (  1 ) =    5.500
cmet   (  1 ) =    5.500
cstr   (  1 ) =    5.500
ccwd   (  1 ) =    5.500
cmic   (  1 ) =    5.500
cslow  (  1 ) =    5.500
cpass  (  1 ) =    5.500
      
cleaf  (  2 ) =    3.900
cwood  (  2 ) =    3.900
cfroot (  2 ) =    3.900
cmet   (  2 ) =    3.900
cstr   (  2 ) =    3.900
ccwd   (  2 ) =    3.900
cmic   (  2 ) =    3.900
cslow  (  2 ) =    3.900
cpass  (  2 ) =    3.900
      
cleaf  (  3 ) =    5.500
cwood  (  3 ) =    5.500
cfroot (  3 ) =    5.500
cmet   (  3 ) =    5.500
cstr   (  3 ) =    5.500
ccwd   (  3 ) =    5.500
cmic   (  3 ) =    5.500
cslow  (  3 ) =    5.500
cpass  (  3 ) =    5.500
      
cleaf  (  4 ) =    3.900
cwood  (  4 ) =    3.900
cfroot (  4 ) =    3.900
cmet   (  4 ) =    3.900
cstr   (  4 ) =    3.900
ccwd   (  4 ) =    3.900
cmic   (  4 ) =    3.900
cslow  (  4 ) =    3.900
cpass  (  4 ) =    3.900
      
cleaf  (  5 ) =    2.000
cwood  (  5 ) =    2.000
cfroot (  5 ) =    2.000
cmet   (  5 ) =    2.000
cstr   (  5 ) =    2.000
ccwd   (  5 ) =    2.000
cmic   (  5 ) =    2.000
cslow  (  5 ) =    2.000
cpass  (  5 ) =    2.000
      
cleaf  (  6 ) =    5.500
cwood  (  6 ) =    5.500
cfroot (  6 ) =    5.500
cmet   (  6 ) =    5.500
cstr   (  6 ) =    5.500
ccwd   (  6 ) =    5.500
cmic   (  6 ) =    5.500
cslow  (  6 ) =    5.500
cpass  (  6 ) =    5.500
      
cleaf  (  7 ) =    5.500
cwood  (  7 ) =    5.500
cfroot (  7 ) =    5.500
cmet   (  7 ) =    5.500
cstr   (  7 ) =    5.500
ccwd   (  7 ) =    5.500
cmic   (  7 ) =    5.500
cslow  (  7 ) =    5.500
cpass  (  7 ) =    5.500
      
cleaf  (  8 ) =    5.500
cwood  (  8 ) =    5.500
cfroot (  8 ) =    5.500
cmet   (  8 ) =    5.500
cstr   (  8 ) =    5.500
ccwd   (  8 ) =    5.500
cmic   (  8 ) =    5.500
cslow  (  8 ) =    5.500
cpass  (  8 ) =    5.500
      
cleaf  (  9 ) =    5.500
cwood  (  9 ) =    5.500
cfroot (  9 ) =    5.500
cmet   (  9 ) =    5.500
cstr   (  9 ) =    5.500
ccwd   (  9 ) =    5.500
cmic   (  9 ) =    5.500
cslow  (  9 ) =    5.500
cpass  (  9 ) =    5.500
      
cleaf  ( 10 ) =    5.500
cwood  ( 10 ) =    5.500
cfroot ( 10 ) =    5.500
cmet   ( 10 ) =    5.500
cstr   ( 10 ) =    5.500
ccwd   ( 10 ) =    5.500
cmic   ( 10 ) =    5.500
cslow  ( 10 ) =    5.500
cpass  ( 10 ) =    5.500
      
cleaf  ( 11 ) =    5.500
cwood  ( 11 ) =    5.500
cfroot ( 11 ) =    5.500
cmet   ( 11 ) =    5.500
cstr   ( 11 ) =    5.500
ccwd   ( 11 ) =    5.500
cmic   ( 11 ) =    5.500
cslow  ( 11 ) =    5.500
cpass  ( 11 ) =    5.500
      
cleaf  ( 12 ) =    5.500
cwood  ( 12 ) =    5.500
cfroot ( 12 ) =    5.500
cmet   ( 12 ) =    5.500
cstr   ( 12 ) =    5.500
ccwd   ( 12 ) =    5.500
cmic   ( 12 ) =    5.500
cslow  ( 12 ) =    5.500
cpass  ( 12 ) =    5.500
      
cleaf  ( 13 ) =    5.500
cwood  ( 13 ) =    5.500
cfroot ( 13 ) =    5.500
cmet   ( 13 ) =    5.500
cstr   ( 13 ) =    5.500
ccwd   ( 13 ) =    5.500
cmic   ( 13 ) =    5.500
cslow  ( 13 ) =    5.500
cpass  ( 13 ) =    5.500
      
cleaf  ( 14 ) =    2.000
cwood  ( 14 ) =    2.000
cfroot ( 14 ) =    2.000
cmet   ( 14 ) =    2.000
cstr   ( 14 ) =    2.000
ccwd   ( 14 ) =    2.000
cmic   ( 14 ) =    2.000
cslow  ( 14 ) =    2.000
cpass  ( 14 ) =    2.000
      
cleaf  ( 15 ) =    2.000
cwood  ( 15 ) =    2.000
cfroot ( 15 ) =    2.000
cmet   ( 15 ) =    2.000
cstr   ( 15 ) =    2.000
ccwd   ( 15 ) =    2.000
cmic   ( 15 ) =    2.000
cslow  ( 15 ) =    2.000
cpass  ( 15 ) =    2.000
      
cleaf  ( 16 ) =    5.500
cwood  ( 16 ) =    5.500
cfroot ( 16 ) =    5.500
cmet   ( 16 ) =    5.500
cstr   ( 16 ) =    5.500
ccwd   ( 16 ) =    5.500
cmic   ( 16 ) =    5.500
cslow  ( 16 ) =    5.500
cpass  ( 16 ) =    5.500
      
cleaf  ( 17 ) =    5.500
cwood  ( 17 ) =    5.500
cfroot ( 17 ) =    5.500
cmet   ( 17 ) =    5.500
cstr   ( 17 ) =    5.500
ccwd   ( 17 ) =    5.500
cmic   ( 17 ) =    5.500
cslow  ( 17 ) =    5.500
cpass  ( 17 ) =    5.500
      
 ! Tranche 4 
 !===============================================================
phen%TKshed              (  1 ) =  268.000
xxkleafcoldmax           (  1 ) =    0.200
casabiome%xkleafcoldexp  (  1 ) =    3.000
xxkleafdrymax            (  1 ) =    0.100
casabiome%xkleafdryexp   (  1 ) =    3.000
      
phen%TKshed              (  2 ) =  260.000
xxkleafcoldmax           (  2 ) =    0.100
casabiome%xkleafcoldexp  (  2 ) =    3.000
xxkleafdrymax            (  2 ) =    0.100
casabiome%xkleafdryexp   (  2 ) =    3.000
      
phen%TKshed              (  3 ) =  263.150
xxkleafcoldmax           (  3 ) =    0.100
casabiome%xkleafcoldexp  (  3 ) =    3.000
xxkleafdrymax            (  3 ) =    0.100
casabiome%xkleafdryexp   (  3 ) =    3.000
      
phen%TKshed              (  4 ) =  268.150
xxkleafcoldmax           (  4 ) =    0.600
casabiome%xkleafcoldexp  (  4 ) =    3.000
xxkleafdrymax            (  4 ) =    1.000
casabiome%xkleafdryexp   (  4 ) =    3.000
      
phen%TKshed              (  5 ) =  277.150
xxkleafcoldmax           (  5 ) =    1.000
casabiome%xkleafcoldexp  (  5 ) =    3.000
xxkleafdrymax            (  5 ) =    0.100
casabiome%xkleafdryexp   (  5 ) =    3.000
      
phen%TKshed              (  6 ) =  275.150
xxkleafcoldmax           (  6 ) =    0.200
casabiome%xkleafcoldexp  (  6 ) =    3.000
xxkleafdrymax            (  6 ) =    0.100
casabiome%xkleafdryexp   (  6 ) =    3.000
      
phen%TKshed              (  7 ) =  275.150
xxkleafcoldmax           (  7 ) =    0.200
casabiome%xkleafcoldexp  (  7 ) =    3.000
xxkleafdrymax            (  7 ) =    0.100
casabiome%xkleafdryexp   (  7 ) =    3.000
      
phen%TKshed              (  8 ) =  275.150
xxkleafcoldmax           (  8 ) =    0.200
casabiome%xkleafcoldexp  (  8 ) =    3.000
xxkleafdrymax            (  8 ) =    0.100
casabiome%xkleafdryexp   (  8 ) =    3.000
      
phen%TKshed              (  9 ) =  278.150
xxkleafcoldmax           (  9 ) =    0.300
casabiome%xkleafcoldexp  (  9 ) =    3.000
xxkleafdrymax            (  9 ) =    0.100
casabiome%xkleafdryexp   (  9 ) =    3.000
      
phen%TKshed              ( 10 ) =  278.150
xxkleafcoldmax           ( 10 ) =    0.300
casabiome%xkleafcoldexp  ( 10 ) =    3.000
xxkleafdrymax            ( 10 ) =    0.100
casabiome%xkleafdryexp   ( 10 ) =    3.000
      
phen%TKshed              ( 11 ) =  277.150
xxkleafcoldmax           ( 11 ) =    0.100
casabiome%xkleafcoldexp  ( 11 ) =    3.000
xxkleafdrymax            ( 11 ) =    0.100
casabiome%xkleafdryexp   ( 11 ) =    3.000
      
phen%TKshed              ( 12 ) =  277.150
xxkleafcoldmax           ( 12 ) =    0.100
casabiome%xkleafcoldexp  ( 12 ) =    3.000
xxkleafdrymax            ( 12 ) =    0.100
casabiome%xkleafdryexp   ( 12 ) =    3.000
      
phen%TKshed              ( 13 ) =  277.150
xxkleafcoldmax           ( 13 ) =    0.100
casabiome%xkleafcoldexp  ( 13 ) =    3.000
xxkleafdrymax            ( 13 ) =    0.100
casabiome%xkleafdryexp   ( 13 ) =    3.000
      
phen%TKshed              ( 14 ) =  277.150
xxkleafcoldmax           ( 14 ) =    0.100
casabiome%xkleafcoldexp  ( 14 ) =    3.000
xxkleafdrymax            ( 14 ) =    0.100
casabiome%xkleafdryexp   ( 14 ) =    3.000
      
phen%TKshed              ( 15 ) =  277.150
xxkleafcoldmax           ( 15 ) =    0.100
casabiome%xkleafcoldexp  ( 15 ) =    3.000
xxkleafdrymax            ( 15 ) =    0.100
casabiome%xkleafdryexp   ( 15 ) =    3.000
      
phen%TKshed              ( 16 ) =  277.150
xxkleafcoldmax           ( 16 ) =    0.100
casabiome%xkleafcoldexp  ( 16 ) =    3.000
xxkleafdrymax            ( 16 ) =    0.100
casabiome%xkleafdryexp   ( 16 ) =    3.000
      
phen%TKshed              ( 17 ) =  283.150
xxkleafcoldmax           ( 17 ) =    0.100
casabiome%xkleafcoldexp  ( 17 ) =    3.000
xxkleafdrymax            ( 17 ) =    0.100
casabiome%xkleafdryexp   ( 17 ) =    3.000
      
 ! Tranche 5 gets skipped
 !===============================================================
      
 ! Tranche 6 
 !===============================================================
casabiome%ratioNCplantmin (  1,  1 ) =    0.020
casabiome%ratioNCplantmax (  1,  1 ) =    0.024
casabiome%ratioNCplantmin (  1,  2 ) =    0.004
casabiome%ratioNCplantmax (  1,  2 ) =    0.005
casabiome%ratioNCplantmin (  1,  3 ) =    0.013
casabiome%ratioNCplantmax (  1,  3 ) =    0.015
xfNminloss                (  1 ) =    0.050
xfNminleach               (  1 ) =    0.050
xnfixrate                 (  1 ) =    0.080
      
casabiome%ratioNCplantmin (  2,  1 ) =    0.040
casabiome%ratioNCplantmax (  2,  1 ) =    0.048
casabiome%ratioNCplantmin (  2,  2 ) =    0.007
casabiome%ratioNCplantmax (  2,  2 ) =    0.008
casabiome%ratioNCplantmin (  2,  3 ) =    0.015
casabiome%ratioNCplantmax (  2,  3 ) =    0.018
xfNminloss                (  2 ) =    0.050
xfNminleach               (  2 ) =    0.050
xnfixrate                 (  2 ) =    2.600
      
casabiome%ratioNCplantmin (  3,  1 ) =    0.017
casabiome%ratioNCplantmax (  3,  1 ) =    0.020
casabiome%ratioNCplantmin (  3,  2 ) =    0.004
casabiome%ratioNCplantmax (  3,  2 ) =    0.005
casabiome%ratioNCplantmin (  3,  3 ) =    0.013
casabiome%ratioNCplantmax (  3,  3 ) =    0.015
xfNminloss                (  3 ) =    0.050
xfNminleach               (  3 ) =    0.050
xnfixrate                 (  3 ) =    0.210
      
casabiome%ratioNCplantmin (  4,  1 ) =    0.029
casabiome%ratioNCplantmax (  4,  1 ) =    0.034
casabiome%ratioNCplantmin (  4,  2 ) =    0.006
casabiome%ratioNCplantmax (  4,  2 ) =    0.007
casabiome%ratioNCplantmin (  4,  3 ) =    0.014
casabiome%ratioNCplantmax (  4,  3 ) =    0.017
xfNminloss                (  4 ) =    0.050
xfNminleach               (  4 ) =    0.050
xnfixrate                 (  4 ) =    1.640
      
casabiome%ratioNCplantmin (  5,  1 ) =    0.025
casabiome%ratioNCplantmax (  5,  1 ) =    0.030
casabiome%ratioNCplantmin (  5,  2 ) =    0.007
casabiome%ratioNCplantmax (  5,  2 ) =    0.008
casabiome%ratioNCplantmin (  5,  3 ) =    0.014
casabiome%ratioNCplantmax (  5,  3 ) =    0.017
xfNminloss                (  5 ) =    0.050
xfNminleach               (  5 ) =    0.050
xnfixrate                 (  5 ) =    0.370
      
casabiome%ratioNCplantmin (  6,  1 ) =    0.026
casabiome%ratioNCplantmax (  6,  1 ) =    0.032
casabiome%ratioNCplantmin (  6,  2 ) =    0.007
casabiome%ratioNCplantmax (  6,  2 ) =    0.008
casabiome%ratioNCplantmin (  6,  3 ) =    0.014
casabiome%ratioNCplantmax (  6,  3 ) =    0.017
xfNminloss                (  6 ) =    0.050
xfNminleach               (  6 ) =    0.050
xnfixrate                 (  6 ) =    0.950
      
casabiome%ratioNCplantmin (  7,  1 ) =    0.020
casabiome%ratioNCplantmax (  7,  1 ) =    0.024
casabiome%ratioNCplantmin (  7,  2 ) =    0.007
casabiome%ratioNCplantmax (  7,  2 ) =    0.008
casabiome%ratioNCplantmin (  7,  3 ) =    0.014
casabiome%ratioNCplantmax (  7,  3 ) =    0.017
xfNminloss                (  7 ) =    0.050
xfNminleach               (  7 ) =    0.050
xnfixrate                 (  7 ) =    0.950
      
casabiome%ratioNCplantmin (  8,  1 ) =    0.020
casabiome%ratioNCplantmax (  8,  1 ) =    0.024
casabiome%ratioNCplantmin (  8,  2 ) =    0.007
casabiome%ratioNCplantmax (  8,  2 ) =    0.008
casabiome%ratioNCplantmin (  8,  3 ) =    0.014
casabiome%ratioNCplantmax (  8,  3 ) =    0.017
xfNminloss                (  8 ) =    0.050
xfNminleach               (  8 ) =    0.050
xnfixrate                 (  8 ) =    0.950
      
casabiome%ratioNCplantmin (  9,  1 ) =    0.040
casabiome%ratioNCplantmax (  9,  1 ) =    0.048
casabiome%ratioNCplantmin (  9,  2 ) =    0.008
casabiome%ratioNCplantmax (  9,  2 ) =    0.010
casabiome%ratioNCplantmin (  9,  3 ) =    0.014
casabiome%ratioNCplantmax (  9,  3 ) =    0.017
xfNminloss                (  9 ) =    0.050
xfNminleach               (  9 ) =    0.050
xnfixrate                 (  9 ) =    4.000
      
casabiome%ratioNCplantmin ( 10,  1 ) =    0.040
casabiome%ratioNCplantmax ( 10,  1 ) =    0.048
casabiome%ratioNCplantmin ( 10,  2 ) =    0.008
casabiome%ratioNCplantmax ( 10,  2 ) =    0.010
casabiome%ratioNCplantmin ( 10,  3 ) =    0.014
casabiome%ratioNCplantmax ( 10,  3 ) =    0.017
xfNminloss                ( 10 ) =    0.050
xfNminleach               ( 10 ) =    0.050
xnfixrate                 ( 10 ) =    4.000
      
casabiome%ratioNCplantmin ( 11,  1 ) =    0.033
casabiome%ratioNCplantmax ( 11,  1 ) =    0.040
casabiome%ratioNCplantmin ( 11,  2 ) =    0.007
casabiome%ratioNCplantmax ( 11,  2 ) =    0.008
casabiome%ratioNCplantmin ( 11,  3 ) =    0.014
casabiome%ratioNCplantmax ( 11,  3 ) =    0.017
xfNminloss                ( 11 ) =    0.050
xfNminleach               ( 11 ) =    0.050
xnfixrate                 ( 11 ) =    0.000
      
casabiome%ratioNCplantmin ( 12,  1 ) =    0.025
casabiome%ratioNCplantmax ( 12,  1 ) =    0.030
casabiome%ratioNCplantmin ( 12,  2 ) =    0.007
casabiome%ratioNCplantmax ( 12,  2 ) =    0.008
casabiome%ratioNCplantmin ( 12,  3 ) =    0.014
casabiome%ratioNCplantmax ( 12,  3 ) =    0.017
xfNminloss                ( 12 ) =    0.050
xfNminleach               ( 12 ) =    0.050
xnfixrate                 ( 12 ) =    0.000
      
casabiome%ratioNCplantmin ( 13,  1 ) =    0.025
casabiome%ratioNCplantmax ( 13,  1 ) =    0.030
casabiome%ratioNCplantmin ( 13,  2 ) =    0.007
casabiome%ratioNCplantmax ( 13,  2 ) =    0.008
casabiome%ratioNCplantmin ( 13,  3 ) =    0.014
casabiome%ratioNCplantmax ( 13,  3 ) =    0.017
xfNminloss                ( 13 ) =    0.050
xfNminleach               ( 13 ) =    0.050
xnfixrate                 ( 13 ) =    0.000
      
casabiome%ratioNCplantmin ( 14,  1 ) =    0.018
casabiome%ratioNCplantmax ( 14,  1 ) =    0.022
casabiome%ratioNCplantmin ( 14,  2 ) =    0.007
casabiome%ratioNCplantmax ( 14,  2 ) =    0.008
casabiome%ratioNCplantmin ( 14,  3 ) =    0.014
casabiome%ratioNCplantmax ( 14,  3 ) =    0.017
xfNminloss                ( 14 ) =    0.050
xfNminleach               ( 14 ) =    0.050
xnfixrate                 ( 14 ) =    0.350
      
casabiome%ratioNCplantmin ( 15,  1 ) =    0.025
casabiome%ratioNCplantmax ( 15,  1 ) =    0.030
casabiome%ratioNCplantmin ( 15,  2 ) =    0.007
casabiome%ratioNCplantmax ( 15,  2 ) =    0.008
casabiome%ratioNCplantmin ( 15,  3 ) =    0.014
casabiome%ratioNCplantmax ( 15,  3 ) =    0.017
xfNminloss                ( 15 ) =    0.050
xfNminleach               ( 15 ) =    0.050
xnfixrate                 ( 15 ) =    0.000
      
casabiome%ratioNCplantmin ( 16,  1 ) =    0.025
casabiome%ratioNCplantmax ( 16,  1 ) =    0.030
casabiome%ratioNCplantmin ( 16,  2 ) =    0.007
casabiome%ratioNCplantmax ( 16,  2 ) =    0.009
casabiome%ratioNCplantmin ( 16,  3 ) =    0.014
casabiome%ratioNCplantmax ( 16,  3 ) =    0.017
xfNminloss                ( 16 ) =    0.050
xfNminleach               ( 16 ) =    0.050
xnfixrate                 ( 16 ) =    0.000
      
casabiome%ratioNCplantmin ( 17,  1 ) =    0.025
casabiome%ratioNCplantmax ( 17,  1 ) =    0.030
casabiome%ratioNCplantmin ( 17,  2 ) =    0.007
casabiome%ratioNCplantmax ( 17,  2 ) =    0.008
casabiome%ratioNCplantmin ( 17,  3 ) =    0.014
casabiome%ratioNCplantmax ( 17,  3 ) =    0.017
xfNminloss                ( 17 ) =    0.050
xfNminleach               ( 17 ) =    0.050
xnfixrate                 ( 17 ) =    0.000
      
 ! Tranche 7 
 !===============================================================
Nleaf      (  1 ) =    7.541
Nwood      (  1 ) =   31.462
Nfroot     (  1 ) =    6.098
Nmet       (  1 ) =    0.064
Nstr       (  1 ) =    1.394
Ncwd       (  1 ) =    2.424
Nmic       (  1 ) =   52.866
Nslow      (  1 ) =  919.729
Npass      (  1 ) =  295.026
xnsoilmin  (  1 ) = 1000.000
      
Nleaf      (  2 ) =    9.900
Nwood      (  2 ) =  102.000
Nfroot     (  2 ) =   38.000
Nmet       (  2 ) =    0.744
Nstr       (  2 ) =    2.892
Ncwd       (  2 ) =    8.524
Nmic       (  2 ) =    1.138
Nslow      (  2 ) =   20.787
Npass      (  2 ) =  880.121
xnsoilmin  (  2 ) = 1000.000
      
Nleaf      (  3 ) =    1.610
Nwood      (  3 ) =   22.734
Nfroot     (  3 ) =    5.366
Nmet       (  3 ) =    0.059
Nstr       (  3 ) =    1.852
Ncwd       (  3 ) =    3.107
Nmic       (  3 ) =   59.708
Nslow      (  3 ) = 1074.741
Npass      (  3 ) =  338.787
xnsoilmin  (  3 ) = 1000.000
      
Nleaf      (  4 ) =    3.757
Nwood      (  4 ) =   80.250
Nfroot     (  4 ) =    5.366
Nmet       (  4 ) =    0.137
Nstr       (  4 ) =    2.084
Ncwd       (  4 ) =    6.582
Nmic       (  4 ) =   40.556
Nslow      (  4 ) =  743.550
Npass      (  4 ) =  336.079
xnsoilmin  (  4 ) = 1000.000
      
Nleaf      (  5 ) =    2.933
Nwood      (  5 ) =    2.756
Nfroot     (  5 ) =    3.415
Nmet       (  5 ) =    0.054
Nstr       (  5 ) =    0.263
Ncwd       (  5 ) =    0.827
Nmic       (  5 ) =   16.805
Nslow      (  5 ) =  297.699
Npass      (  5 ) =   92.432
xnsoilmin  (  5 ) = 1000.000
      
Nleaf      (  6 ) =    4.572
Nwood      (  6 ) =    0.000
Nfroot     (  6 ) =    6.415
Nmet       (  6 ) =    0.476
Nstr       (  6 ) =    0.339
Ncwd       (  6 ) =    0.000
Nmic       (  6 ) =   42.564
Nslow      (  6 ) =  379.629
Npass      (  6 ) =  278.661
xnsoilmin  (  6 ) = 1000.000
      
Nleaf      (  7 ) =    4.572
Nwood      (  7 ) =    0.000
Nfroot     (  7 ) =    6.415
Nmet       (  7 ) =    0.476
Nstr       (  7 ) =    0.339
Ncwd       (  7 ) =    0.000
Nmic       (  7 ) =   42.564
Nslow      (  7 ) =  379.629
Npass      (  7 ) =  278.661
xnsoilmin  (  7 ) = 1000.000
      
Nleaf      (  8 ) =    4.572
Nwood      (  8 ) =    0.000
Nfroot     (  8 ) =    6.415
Nmet       (  8 ) =    0.476
Nstr       (  8 ) =    0.339
Ncwd       (  8 ) =    0.000
Nmic       (  8 ) =   42.564
Nslow      (  8 ) =  379.629
Npass      (  8 ) =  278.661
xnsoilmin  (  8 ) = 1000.000
      
Nleaf      (  9 ) =    5.333
Nwood      (  9 ) =    0.000
Nfroot     (  9 ) =    5.854
Nmet       (  9 ) =    0.476
Nstr       (  9 ) =    0.339
Ncwd       (  9 ) =    0.000
Nmic       (  9 ) =   51.242
Nslow      (  9 ) =  457.029
Npass      (  9 ) =  335.476
xnsoilmin  (  9 ) = 1000.000
      
Nleaf      ( 10 ) =    5.333
Nwood      ( 10 ) =    0.000
Nfroot     ( 10 ) =    5.854
Nmet       ( 10 ) =    0.476
Nstr       ( 10 ) =    0.339
Ncwd       ( 10 ) =    0.000
Nmic       ( 10 ) =   51.242
Nslow      ( 10 ) =  457.029
Npass      ( 10 ) =  335.476
xnsoilmin  ( 10 ) = 1000.000
      
Nleaf      ( 11 ) =    0.000
Nwood      ( 11 ) =    0.000
Nfroot     ( 11 ) =    0.000
Nmet       ( 11 ) =    0.000
Nstr       ( 11 ) =    0.000
Ncwd       ( 11 ) =    0.000
Nmic       ( 11 ) =    0.000
Nslow      ( 11 ) =    0.000
Npass      ( 11 ) =    0.000
xnsoilmin  ( 11 ) = 1000.000
      
Nleaf      ( 12 ) =    0.000
Nwood      ( 12 ) =    0.000
Nfroot     ( 12 ) =    0.000
Nmet       ( 12 ) =    0.000
Nstr       ( 12 ) =    0.000
Ncwd       ( 12 ) =    0.000
Nmic       ( 12 ) =    0.000
Nslow      ( 12 ) =    0.000
Npass      ( 12 ) =    0.000
xnsoilmin  ( 12 ) = 1000.000
      
Nleaf      ( 13 ) =    0.000
Nwood      ( 13 ) =    0.000
Nfroot     ( 13 ) =    0.000
Nmet       ( 13 ) =    0.000
Nstr       ( 13 ) =    0.000
Ncwd       ( 13 ) =    0.000
Nmic       ( 13 ) =    0.000
Nslow      ( 13 ) =    0.000
Npass      ( 13 ) =    0.000
xnsoilmin  ( 13 ) = 1000.000
      
Nleaf      ( 14 ) =    0.500
Nwood      ( 14 ) =    0.126
Nfroot     ( 14 ) =    1.537
Nmet       ( 14 ) =    0.018
Nstr       ( 14 ) =    0.033
Ncwd       ( 14 ) =    0.211
Nmic       ( 14 ) =    5.778
Nslow      ( 14 ) =   88.337
Npass      ( 14 ) =   34.478
xnsoilmin  ( 14 ) = 1000.000
      
Nleaf      ( 15 ) =    0.000
Nwood      ( 15 ) =    0.000
Nfroot     ( 15 ) =    0.000
Nmet       ( 15 ) =    0.000
Nstr       ( 15 ) =    0.000
Ncwd       ( 15 ) =    0.000
Nmic       ( 15 ) =    0.000
Nslow      ( 15 ) =    0.000
Npass      ( 15 ) =    0.000
xnsoilmin  ( 15 ) = 1000.000
      
Nleaf      ( 16 ) =    0.000
Nwood      ( 16 ) =    0.000
Nfroot     ( 16 ) =    0.000
Nmet       ( 16 ) =    0.000
Nstr       ( 16 ) =    0.000
Ncwd       ( 16 ) =    0.000
Nmic       ( 16 ) =    0.000
Nslow      ( 16 ) =    0.000
Npass      ( 16 ) =    0.000
xnsoilmin  ( 16 ) = 1000.000
      
Nleaf      ( 17 ) =    0.000
Nwood      ( 17 ) =    0.000
Nfroot     ( 17 ) =    0.000
Nmet       ( 17 ) =    0.000
Nstr       ( 17 ) =    0.000
Ncwd       ( 17 ) =    0.000
Nmic       ( 17 ) =    0.000
Nslow      ( 17 ) =    0.000
Npass      ( 17 ) =    0.000
xnsoilmin  ( 17 ) = 1000.000
      
 ! Tranche 8 
 !===============================================================
xratioNPleafmin                            =   10.000
xratioNPleafmax                            =   10.000
xratioNPwoodmin                            =   15.000
xratioNPwoodmax                            =   15.000
xratioNPfrootmin                           =   15.000
xratioNPfrootmax                           =   15.000
casabiome%ftransPPtoL   (  1,  1 ) =    0.500
casabiome%ftransPPtoL   (  1,  2 ) =    0.950
casabiome%ftransPPtoL   (  1,  3 ) =    0.900
      
xratioNPleafmin                            =   10.000
xratioNPleafmax                            =   10.000
xratioNPwoodmin                            =   15.000
xratioNPwoodmax                            =   15.000
xratioNPfrootmin                           =   15.000
xratioNPfrootmax                           =   15.000
casabiome%ftransPPtoL   (  2,  1 ) =    0.500
casabiome%ftransPPtoL   (  2,  2 ) =    0.950
casabiome%ftransPPtoL   (  2,  3 ) =    0.900
      
xratioNPleafmin                            =   10.000
xratioNPleafmax                            =   10.000
xratioNPwoodmin                            =   15.000
xratioNPwoodmax                            =   15.000
xratioNPfrootmin                           =   15.000
xratioNPfrootmax                           =   15.000
casabiome%ftransPPtoL   (  3,  1 ) =    0.500
casabiome%ftransPPtoL   (  3,  2 ) =    0.950
casabiome%ftransPPtoL   (  3,  3 ) =    0.900
      
xratioNPleafmin                            =   10.000
xratioNPleafmax                            =   10.000
xratioNPwoodmin                            =   15.000
xratioNPwoodmax                            =   15.000
xratioNPfrootmin                           =   15.000
xratioNPfrootmax                           =   15.000
casabiome%ftransPPtoL   (  4,  1 ) =    0.500
casabiome%ftransPPtoL   (  4,  2 ) =    0.950
casabiome%ftransPPtoL   (  4,  3 ) =    0.900
      
xratioNPleafmin                            =   10.000
xratioNPleafmax                            =   10.000
xratioNPwoodmin                            =   15.000
xratioNPwoodmax                            =   15.000
xratioNPfrootmin                           =   15.000
xratioNPfrootmax                           =   15.000
casabiome%ftransPPtoL   (  5,  1 ) =    0.500
casabiome%ftransPPtoL   (  5,  2 ) =    0.950
casabiome%ftransPPtoL   (  5,  3 ) =    0.900
      
xratioNPleafmin                            =   10.000
xratioNPleafmax                            =   10.000
xratioNPwoodmin                            =   15.000
xratioNPwoodmax                            =   15.000
xratioNPfrootmin                           =   15.000
xratioNPfrootmax                           =   15.000
casabiome%ftransPPtoL   (  6,  1 ) =    0.500
casabiome%ftransPPtoL   (  6,  2 ) =    0.950
casabiome%ftransPPtoL   (  6,  3 ) =    0.900
      
xratioNPleafmin                            =   10.000
xratioNPleafmax                            =   10.000
xratioNPwoodmin                            =   15.000
xratioNPwoodmax                            =   15.000
xratioNPfrootmin                           =   15.000
xratioNPfrootmax                           =   15.000
casabiome%ftransPPtoL   (  7,  1 ) =    0.500
casabiome%ftransPPtoL   (  7,  2 ) =    0.950
casabiome%ftransPPtoL   (  7,  3 ) =    0.900
      
xratioNPleafmin                            =   10.000
xratioNPleafmax                            =   10.000
xratioNPwoodmin                            =   15.000
xratioNPwoodmax                            =   15.000
xratioNPfrootmin                           =   15.000
xratioNPfrootmax                           =   15.000
casabiome%ftransPPtoL   (  8,  1 ) =    0.500
casabiome%ftransPPtoL   (  8,  2 ) =    0.950
casabiome%ftransPPtoL   (  8,  3 ) =    0.900
      
xratioNPleafmin                            =   10.000
xratioNPleafmax                            =   10.000
xratioNPwoodmin                            =   15.000
xratioNPwoodmax                            =   15.000
xratioNPfrootmin                           =   15.000
xratioNPfrootmax                           =   15.000
casabiome%ftransPPtoL   (  9,  1 ) =    0.500
casabiome%ftransPPtoL   (  9,  2 ) =    0.950
casabiome%ftransPPtoL   (  9,  3 ) =    0.900
      
xratioNPleafmin                            =   10.000
xratioNPleafmax                            =   10.000
xratioNPwoodmin                            =   15.000
xratioNPwoodmax                            =   15.000
xratioNPfrootmin                           =   15.000
xratioNPfrootmax                           =   15.000
casabiome%ftransPPtoL   ( 10,  1 ) =    0.500
casabiome%ftransPPtoL   ( 10,  2 ) =    0.950
casabiome%ftransPPtoL   ( 10,  3 ) =    0.900
      
xratioNPleafmin                            =   10.000
xratioNPleafmax                            =   10.000
xratioNPwoodmin                            =   15.000
xratioNPwoodmax                            =   15.000
xratioNPfrootmin                           =   15.000
xratioNPfrootmax                           =   15.000
casabiome%ftransPPtoL   ( 11,  1 ) =    0.500
casabiome%ftransPPtoL   ( 11,  2 ) =    0.950
casabiome%ftransPPtoL   ( 11,  3 ) =    0.900
      
xratioNPleafmin                            =   10.000
xratioNPleafmax                            =   10.000
xratioNPwoodmin                            =   15.000
xratioNPwoodmax                            =   15.000
xratioNPfrootmin                           =   15.000
xratioNPfrootmax                           =   15.000
casabiome%ftransPPtoL   ( 12,  1 ) =    0.500
casabiome%ftransPPtoL   ( 12,  2 ) =    0.950
casabiome%ftransPPtoL   ( 12,  3 ) =    0.900
      
xratioNPleafmin                            =   10.000
xratioNPleafmax                            =   10.000
xratioNPwoodmin                            =   15.000
xratioNPwoodmax                            =   15.000
xratioNPfrootmin                           =   15.000
xratioNPfrootmax                           =   15.000
casabiome%ftransPPtoL   ( 13,  1 ) =    0.500
casabiome%ftransPPtoL   ( 13,  2 ) =    0.950
casabiome%ftransPPtoL   ( 13,  3 ) =    0.900
      
xratioNPleafmin                            =   10.000
xratioNPleafmax                            =   10.000
xratioNPwoodmin                            =   15.000
xratioNPwoodmax                            =   15.000
xratioNPfrootmin                           =   15.000
xratioNPfrootmax                           =   15.000
casabiome%ftransPPtoL   ( 14,  1 ) =    0.500
casabiome%ftransPPtoL   ( 14,  2 ) =    0.950
casabiome%ftransPPtoL   ( 14,  3 ) =    0.900
      
xratioNPleafmin                            =   10.000
xratioNPleafmax                            =   10.000
xratioNPwoodmin                            =   15.000
xratioNPwoodmax                            =   15.000
xratioNPfrootmin                           =   15.000
xratioNPfrootmax                           =   15.000
casabiome%ftransPPtoL   ( 15,  1 ) =    0.500
casabiome%ftransPPtoL   ( 15,  2 ) =    0.950
casabiome%ftransPPtoL   ( 15,  3 ) =    0.900
      
xratioNPleafmin                            =   10.000
xratioNPleafmax                            =   10.000
xratioNPwoodmin                            =   15.000
xratioNPwoodmax                            =   15.000
xratioNPfrootmin                           =   15.000
xratioNPfrootmax                           =   15.000
casabiome%ftransPPtoL   ( 16,  1 ) =    0.500
casabiome%ftransPPtoL   ( 16,  2 ) =    0.950
casabiome%ftransPPtoL   ( 16,  3 ) =    0.900
      
xratioNPleafmin                            =   10.000
xratioNPleafmax                            =   10.000
xratioNPwoodmin                            =   15.000
xratioNPwoodmax                            =   15.000
xratioNPfrootmin                           =   15.000
xratioNPfrootmax                           =   15.000
casabiome%ftransPPtoL   ( 17,  1 ) =    0.500
casabiome%ftransPPtoL   ( 17,  2 ) =    0.950
casabiome%ftransPPtoL   ( 17,  3 ) =    0.900
      
 ! Tranche 9 
 !===============================================================
xkmlabp                    (  1 ) =   74.541
xpsorbmax                  (  1 ) =  745.408
xfPleach                   (  1 ) =    0.001
ratioNPsoil                (  1,  1 ) =    4.000
ratioNPsoil                (  1,  2 ) =    5.000
ratioNPsoil                (  1,  3 ) =    5.000
xxkplab                    (  1 ) =    0.001
xxkpsorb                   (  1 ) =    0.000
xxkpocc                    (  1 ) =    0.000
      
xkmlabp                    (  2 ) =   68.158
xpsorbmax                  (  2 ) =  788.082
xfPleach                   (  2 ) =    0.001
ratioNPsoil                (  2,  1 ) =    4.000
ratioNPsoil                (  2,  2 ) =    5.000
ratioNPsoil                (  2,  3 ) =    5.000
xxkplab                    (  2 ) =    0.001
xxkpsorb                   (  2 ) =    0.000
xxkpocc                    (  2 ) =    0.000
      
xkmlabp                    (  3 ) =   77.952
xpsorbmax                  (  3 ) = 1110.816
xfPleach                   (  3 ) =    0.001
ratioNPsoil                (  3,  1 ) =    4.000
ratioNPsoil                (  3,  2 ) =    5.000
ratioNPsoil                (  3,  3 ) =    5.000
xxkplab                    (  3 ) =    0.001
xxkpsorb                   (  3 ) =    0.000
xxkpocc                    (  3 ) =    0.000
      
xkmlabp                    (  4 ) =   64.419
xpsorbmax                  (  4 ) =  744.847
xfPleach                   (  4 ) =    0.001
ratioNPsoil                (  4,  1 ) =    4.000
ratioNPsoil                (  4,  2 ) =   15.000
ratioNPsoil                (  4,  3 ) =   15.000
xxkplab                    (  4 ) =    0.001
xxkpsorb                   (  4 ) =    0.000
xxkpocc                    (  4 ) =    0.000
      
xkmlabp                    (  5 ) =   64.419
xpsorbmax                  (  5 ) =  744.847
xfPleach                   (  5 ) =    0.001
ratioNPsoil                (  5,  1 ) =    4.000
ratioNPsoil                (  5,  2 ) =    5.000
ratioNPsoil                (  5,  3 ) =    5.000
xxkplab                    (  5 ) =    0.001
xxkpsorb                   (  5 ) =    0.000
xxkpocc                    (  5 ) =    0.000
      
xkmlabp                    (  6 ) =   70.586
xpsorbmax                  (  6 ) =  816.146
xfPleach                   (  6 ) =    0.001
ratioNPsoil                (  6,  1 ) =    4.000
ratioNPsoil                (  6,  2 ) =    5.000
ratioNPsoil                (  6,  3 ) =    5.000
xxkplab                    (  6 ) =    0.001
xxkpsorb                   (  6 ) =    0.000
xxkpocc                    (  6 ) =    0.000
      
xkmlabp                    (  7 ) =   64.589
xpsorbmax                  (  7 ) =  746.808
xfPleach                   (  7 ) =    0.001
ratioNPsoil                (  7,  1 ) =    4.000
ratioNPsoil                (  7,  2 ) =    5.000
ratioNPsoil                (  7,  3 ) =    5.000
xxkplab                    (  7 ) =    0.001
xxkpsorb                   (  7 ) =    0.000
xxkpocc                    (  7 ) =    0.000
      
xkmlabp                    (  8 ) =   54.169
xpsorbmax                  (  8 ) =  722.256
xfPleach                   (  8 ) =    0.001
ratioNPsoil                (  8,  1 ) =    4.000
ratioNPsoil                (  8,  2 ) =    5.000
ratioNPsoil                (  8,  3 ) =    5.000
xxkplab                    (  8 ) =    0.001
xxkpsorb                   (  8 ) =    0.000
xxkpocc                    (  8 ) =    0.000
      
xkmlabp                    (  9 ) =    9.770
xpsorbmax                  (  9 ) =  293.112
xfPleach                   (  9 ) =    0.001
ratioNPsoil                (  9,  1 ) =    4.000
ratioNPsoil                (  9,  2 ) =    7.000
ratioNPsoil                (  9,  3 ) =    7.000
xxkplab                    (  9 ) =    0.001
xxkpsorb                   (  9 ) =    0.000
xxkpocc                    (  9 ) =    0.000
      
xkmlabp                    ( 10 ) =   28.290
xpsorbmax                  ( 10 ) =  311.190
xfPleach                   ( 10 ) =    0.001
ratioNPsoil                ( 10,  1 ) =    4.000
ratioNPsoil                ( 10,  2 ) =    7.000
ratioNPsoil                ( 10,  3 ) =    7.000
xxkplab                    ( 10 ) =    0.001
xxkpsorb                   ( 10 ) =    0.000
xxkpocc                    ( 10 ) =    0.000
      
xkmlabp                    ( 11 ) =   63.963
xpsorbmax                  ( 11 ) =  373.118
xfPleach                   ( 11 ) =    0.001
ratioNPsoil                ( 11,  1 ) =    4.000
ratioNPsoil                ( 11,  2 ) =    7.000
ratioNPsoil                ( 11,  3 ) =    7.000
xxkplab                    ( 11 ) =    0.001
xxkpsorb                   ( 11 ) =    0.000
xxkpocc                    ( 11 ) =    0.000
      
xkmlabp                    ( 12 ) =   32.402
xpsorbmax                  ( 12 ) =  615.638
xfPleach                   ( 12 ) =    0.001
ratioNPsoil                ( 12,  1 ) =    4.000
ratioNPsoil                ( 12,  2 ) =    7.000
ratioNPsoil                ( 12,  3 ) =    7.000
xxkplab                    ( 12 ) =    0.001
xxkpsorb                   ( 12 ) =    0.000
xxkpocc                    ( 12 ) =    0.000
      
 ! Tranche 10
 !===============================================================
xpleaf  (  1 ) =    0.192
xpwood  (  1 ) =    0.954
xpfroot (  1 ) =    0.077
xpmet   (  1 ) =    0.004
xpstr   (  1 ) =    0.070
xpcwd   (  1 ) =    0.101
xpmic   (  1 ) =    6.873
xpslow  (  1 ) =  119.565
xppass  (  1 ) =   38.353
xplab   (  1 ) =   26.737
xpsorb  (  1 ) =  126.730
xpocc   (  1 ) =  138.571
      
xpleaf  (  2 ) =    0.415
xpwood  (  2 ) =    5.880
xpfroot (  2 ) =    1.950
xpmet   (  2 ) =    0.030
xpstr   (  2 ) =    0.145
xpcwd   (  2 ) =    0.192
xpmic   (  2 ) =    0.148
xpslow  (  2 ) =    2.702
xppass  (  2 ) =  114.416
xplab   (  2 ) =   19.947
xpsorb  (  2 ) =   92.263
xpocc   (  2 ) =  120.374
      
xpleaf  (  3 ) =    0.116
xpwood  (  3 ) =    0.644
xpfroot (  3 ) =    0.081
xpmet   (  3 ) =    0.005
xpstr   (  3 ) =    0.093
xpcwd   (  3 ) =    0.129
xpmic   (  3 ) =    7.762
xpslow  (  3 ) =  139.716
xppass  (  3 ) =   44.042
xplab   (  3 ) =   29.107
xpsorb  (  3 ) =  134.639
xpocc   (  3 ) =  138.220
      
xpleaf  (  4 ) =    0.135
xpwood  (  4 ) =    2.425
xpfroot (  4 ) =    0.141
xpmet   (  4 ) =    0.007
xpstr   (  4 ) =    0.104
xpcwd   (  4 ) =    0.148
xpmic   (  4 ) =    5.272
xpslow  (  4 ) =   96.662
xppass  (  4 ) =   43.690
xplab   (  4 ) =   30.509
xpsorb  (  4 ) =  132.012
xpocc   (  4 ) =  148.083
      
xpleaf  (  5 ) =    0.023
xpwood  (  5 ) =    0.000
xpfroot (  5 ) =    0.037
xpmet   (  5 ) =    0.002
xpstr   (  5 ) =    0.013
xpcwd   (  5 ) =    0.019
xpmic   (  5 ) =    2.185
xpslow  (  5 ) =   38.701
xppass  (  5 ) =   12.016
xplab   (  5 ) =   23.206
xpsorb  (  5 ) =  173.470
xpocc   (  5 ) =  114.496
      
xpleaf  (  6 ) =    0.151
xpwood  (  6 ) =    0.000
xpfroot (  6 ) =    0.151
xpmet   (  6 ) =    0.019
xpstr   (  6 ) =    0.017
xpcwd   (  6 ) =    0.000
xpmic   (  6 ) =    5.533
xpslow  (  6 ) =   49.352
xppass  (  6 ) =   36.226
xplab   (  6 ) =   25.538
xpsorb  (  6 ) =  186.207
xpocc   (  6 ) =  145.163
      
xpleaf  (  7 ) =    0.151
xpwood  (  7 ) =    0.000
xpfroot (  7 ) =    0.151
xpmet   (  7 ) =    0.019
xpstr   (  7 ) =    0.017
xpcwd   (  7 ) =    0.000
xpmic   (  7 ) =    5.533
xpslow  (  7 ) =   49.352
xppass  (  7 ) =   36.226
xplab   (  7 ) =   25.538
xpsorb  (  7 ) =  186.207
xpocc   (  7 ) =  145.163
      
xpleaf  (  8 ) =    0.151
xpwood  (  8 ) =    0.000
xpfroot (  8 ) =    0.151
xpmet   (  8 ) =    0.019
xpstr   (  8 ) =    0.017
xpcwd   (  8 ) =    0.000
xpmic   (  8 ) =    5.533
xpslow  (  8 ) =   49.352
xppass  (  8 ) =   36.226
xplab   (  8 ) =   25.538
xpsorb  (  8 ) =  186.207
xpocc   (  8 ) =  145.163
      
xpleaf  (  9 ) =    0.151
xpwood  (  9 ) =    0.000
xpfroot (  9 ) =    0.151
xpmet   (  9 ) =    0.019
xpstr   (  9 ) =    0.017
xpcwd   (  9 ) =    0.000
xpmic   (  9 ) =    6.662
xpslow  (  9 ) =   59.414
xppass  (  9 ) =   43.612
xplab   (  9 ) =   27.729
xpsorb  (  9 ) =  155.518
xpocc   (  9 ) =  158.884
      
xpleaf  ( 10 ) =    0.151
xpwood  ( 10 ) =    0.000
xpfroot ( 10 ) =    0.151
xpmet   ( 10 ) =    0.019
xpstr   ( 10 ) =    0.017
xpcwd   ( 10 ) =    0.000
xpmic   ( 10 ) =    6.662
xpslow  ( 10 ) =   59.414
xppass  ( 10 ) =   43.612
xplab   ( 10 ) =   27.729
xpsorb  ( 10 ) =  155.518
xpocc   ( 10 ) =  158.884
      
xpleaf  ( 11 ) =    0.000
xpwood  ( 11 ) =    0.000
xpfroot ( 11 ) =    0.000
xpmet   ( 11 ) =    0.000
xpstr   ( 11 ) =    0.000
xpcwd   ( 11 ) =    0.000
xpmic   ( 11 ) =    0.000
xpslow  ( 11 ) =    0.000
xppass  ( 11 ) =    0.000
xplab   ( 11 ) =    0.000
xpsorb  ( 11 ) =    0.000
xpocc   ( 11 ) =    0.000
      
xpleaf  ( 12 ) =    0.000
xpwood  ( 12 ) =    0.000
xpfroot ( 12 ) =    0.000
xpmet   ( 12 ) =    0.000
xpstr   ( 12 ) =    0.000
xpcwd   ( 12 ) =    0.000
xpmic   ( 12 ) =    0.000
xpslow  ( 12 ) =    0.000
xppass  ( 12 ) =    0.000
xplab   ( 12 ) =    0.000
xpsorb  ( 12 ) =    0.000
xpocc   ( 12 ) =    0.000
      
xpleaf  ( 13 ) =    0.000
xpwood  ( 13 ) =    0.000
xpfroot ( 13 ) =    0.000
xpmet   ( 13 ) =    0.000
xpstr   ( 13 ) =    0.000
xpcwd   ( 13 ) =    0.000
xpmic   ( 13 ) =    0.000
xpslow  ( 13 ) =    0.000
xppass  ( 13 ) =    0.000
xplab   ( 13 ) =    0.000
xpsorb  ( 13 ) =    0.000
xpocc   ( 13 ) =    0.000
      
xpleaf  ( 14 ) =    0.007
xpwood  ( 14 ) =    0.000
xpfroot ( 14 ) =    0.009
xpmet   ( 14 ) =    0.001
xpstr   ( 14 ) =    0.002
xpcwd   ( 14 ) =    0.000
xpmic   ( 14 ) =    0.751
xpslow  ( 14 ) =   11.484
xppass  ( 14 ) =    4.482
xplab   ( 14 ) =   21.038
xpsorb  ( 14 ) =  255.790
xpocc   ( 14 ) =  108.897
      
xpleaf  ( 15 ) =    0.000
xpwood  ( 15 ) =    0.000
xpfroot ( 15 ) =    0.000
xpmet   ( 15 ) =    0.000
xpstr   ( 15 ) =    0.000
xpcwd   ( 15 ) =    0.000
xpmic   ( 15 ) =    0.000
xpslow  ( 15 ) =    0.000
xppass  ( 15 ) =    0.000
xplab   ( 15 ) =    0.000
xpsorb  ( 15 ) =    0.000
xpocc   ( 15 ) =    0.000
      
xpleaf  ( 16 ) =    0.000
xpwood  ( 16 ) =    0.000
xpfroot ( 16 ) =    0.000
xpmet   ( 16 ) =    0.000
xpstr   ( 16 ) =    0.000
xpcwd   ( 16 ) =    0.000
xpmic   ( 16 ) =    0.000
xpslow  ( 16 ) =    0.000
xppass  ( 16 ) =    0.000
xplab   ( 16 ) =    0.000
xpsorb  ( 16 ) =    0.000
xpocc   ( 16 ) =    0.000
      
xpleaf  ( 17 ) =    0.000
xpwood  ( 17 ) =    0.000
xpfroot ( 17 ) =    0.000
xpmet   ( 17 ) =    0.000
xpstr   ( 17 ) =    0.000
xpcwd   ( 17 ) =    0.000
xpmic   ( 17 ) =    0.000
xpslow  ( 17 ) =    0.000
xppass  ( 17 ) =    0.000
xplab   ( 17 ) =    0.103
xpsorb  ( 17 ) =    1.176
xpocc   ( 17 ) =    0.688
      
 ! Tranche 11
 !===============================================================
xxnpmax        (  1 ) =    1.511
xq10soil       (  1 ) =    1.720
xxkoptlitter   (  1 ) =    0.400
xxkoptsoil     (  1 ) =    0.330
xprodptase     (  1 ) =    0.500
xcostnpup      (  1 ) =   40.000
xmaxfinelitter (  1 ) = 1524.000
xmaxcwd        (  1 ) = 1795.000
xnintercept    (  1 ) =    6.320
xnslope        (  1 ) =   18.150
      
xxnpmax        (  2 ) =    1.279
xq10soil       (  2 ) =    1.720
xxkoptlitter   (  2 ) =    0.400
xxkoptsoil     (  2 ) =    0.600
xprodptase     (  2 ) =    0.200
xcostnpup      (  2 ) =   25.000
xmaxfinelitter (  2 ) =  384.000
xmaxcwd        (  2 ) =  613.000
xnintercept    (  2 ) =    4.190
xnslope        (  2 ) =   26.190
      
xxnpmax        (  3 ) =    1.591
xq10soil       (  3 ) =    1.720
xxkoptlitter   (  3 ) =    0.400
xxkoptsoil     (  3 ) =    0.150
xprodptase     (  3 ) =    0.500
xcostnpup      (  3 ) =   40.000
xmaxfinelitter (  3 ) = 1527.000
xmaxcwd        (  3 ) = 1918.000
xnintercept    (  3 ) =    6.320
xnslope        (  3 ) =   18.150
      
xxnpmax        (  4 ) =    1.186
xq10soil       (  4 ) =    1.720
xxkoptlitter   (  4 ) =    0.400
xxkoptsoil     (  4 ) =    0.600
xprodptase     (  4 ) =    0.500
xcostnpup      (  4 ) =   40.000
xmaxfinelitter (  4 ) =  887.000
xmaxcwd        (  4 ) = 1164.000
xnintercept    (  4 ) =    5.730
xnslope        (  4 ) =   29.810
      
xxnpmax        (  5 ) =    1.358
xq10soil       (  5 ) =    1.720
xxkoptlitter   (  5 ) =    0.400
xxkoptsoil     (  5 ) =    0.160
xprodptase     (  5 ) =    0.500
xcostnpup      (  5 ) =   40.000
xmaxfinelitter (  5 ) =  157.000
xmaxcwd        (  5 ) =  107.000
xnintercept    (  5 ) =   14.710
xnslope        (  5 ) =   23.150
      
xxnpmax        (  6 ) =    1.456
xq10soil       (  6 ) =    1.720
xxkoptlitter   (  6 ) =    0.400
xxkoptsoil     (  6 ) =    0.400
xprodptase     (  6 ) =    0.500
xcostnpup      (  6 ) =   40.000
xmaxfinelitter (  6 ) =  361.000
xmaxcwd        (  6 ) =  420.000
xnintercept    (  6 ) =    6.420
xnslope        (  6 ) =   40.960
      
xxnpmax        (  7 ) =    1.456
xq10soil       (  7 ) =    1.720
xxkoptlitter   (  7 ) =    0.400
xxkoptsoil     (  7 ) =    0.300
xprodptase     (  7 ) =    0.500
xcostnpup      (  7 ) =   40.000
xmaxfinelitter (  7 ) =  225.000
xmaxcwd        (  7 ) =  228.000
xnintercept    (  7 ) =    2.000
xnslope        (  7 ) =    8.000
      
xxnpmax        (  8 ) =    1.456
xq10soil       (  8 ) =    1.720
xxkoptlitter   (  8 ) =    0.400
xxkoptsoil     (  8 ) =    0.200
xprodptase     (  8 ) =    0.500
xcostnpup      (  8 ) =   40.000
xmaxfinelitter (  8 ) =  913.000
xmaxcwd        (  8 ) =  573.000
xnintercept    (  8 ) =   14.710
xnslope        (  8 ) =   23.150
      
xxnpmax        (  9 ) =    1.210
xq10soil       (  9 ) =    1.720
xxkoptlitter   (  9 ) =    0.400
xxkoptsoil     (  9 ) =    0.200
xprodptase     (  9 ) =    0.500
xcostnpup      (  9 ) =   40.000
xmaxfinelitter (  9 ) =  660.000
xmaxcwd        (  9 ) =  811.000
xnintercept    (  9 ) =    4.710
xnslope        (  9 ) =   59.230
      
xxnpmax        ( 10 ) =    1.210
xq10soil       ( 10 ) =    1.720
xxkoptlitter   ( 10 ) =    0.400
xxkoptsoil     ( 10 ) =    0.250
xprodptase     ( 10 ) =    0.500
xcostnpup      ( 10 ) =   40.000
xmaxfinelitter ( 10 ) =  100.000
xmaxcwd        ( 10 ) =  100.000
xnintercept    ( 10 ) =   14.710
xnslope        ( 10 ) =   23.150
      
xxnpmax        ( 11 ) =    1.456
xq10soil       ( 11 ) =    1.720
xxkoptlitter   ( 11 ) =    0.400
xxkoptsoil     ( 11 ) =    1.000
xprodptase     ( 11 ) =    0.500
xcostnpup      ( 11 ) =   40.000
xmaxfinelitter ( 11 ) =  100.000
xmaxcwd        ( 11 ) =  100.000
xnintercept    ( 11 ) =   14.710
xnslope        ( 11 ) =   23.150
      
xxnpmax        ( 12 ) =    1.366
xq10soil       ( 12 ) =    1.720
xxkoptlitter   ( 12 ) =    0.400
xxkoptsoil     ( 12 ) =    0.650
xprodptase     ( 12 ) =    4.000
xcostnpup      ( 12 ) =   40.000
xmaxfinelitter ( 12 ) =  100.000
xmaxcwd        ( 12 ) =  100.000
xnintercept    ( 12 ) =    7.000
xnslope        ( 12 ) =   10.000
      
xxnpmax        ( 13 ) =    1.210
xq10soil       ( 13 ) =    1.720
xxkoptlitter   ( 13 ) =    0.400
xxkoptsoil     ( 13 ) =    0.500
xprodptase     ( 13 ) =    0.500
xcostnpup      ( 13 ) =   40.000
xmaxfinelitter ( 13 ) =  100.000
xmaxcwd        ( 13 ) =  100.000
xnintercept    ( 13 ) =   14.710
xnslope        ( 13 ) =   23.150
      
xxnpmax        ( 14 ) =    1.000
xq10soil       ( 14 ) =    1.720
xxkoptlitter   ( 14 ) =    0.400
xxkoptsoil     ( 14 ) =    2.000
xprodptase     ( 14 ) =    0.500
xcostnpup      ( 14 ) =   40.000
xmaxfinelitter ( 14 ) =   83.000
xmaxcwd        ( 14 ) =   23.000
xnintercept    ( 14 ) =   14.710
xnslope        ( 14 ) =   23.150
      
xxnpmax        ( 15 ) =    1.400
xq10soil       ( 15 ) =    1.720
xxkoptlitter   ( 15 ) =    0.400
xxkoptsoil     ( 15 ) =    0.500
xprodptase     ( 15 ) =    0.500
xcostnpup      ( 15 ) =   40.000
xmaxfinelitter ( 15 ) =  100.000
xmaxcwd        ( 15 ) =  100.000
xnintercept    ( 15 ) =   14.710
xnslope        ( 15 ) =   23.150
      
xxnpmax        ( 16 ) =    1.000
xq10soil       ( 16 ) =    1.720
xxkoptlitter   ( 16 ) =    0.400
xxkoptsoil     ( 16 ) =    1.000
xprodptase     ( 16 ) =    0.500
xcostnpup      ( 16 ) =   40.000
xmaxfinelitter ( 16 ) =  100.000
xmaxcwd        ( 16 ) =  100.000
xnintercept    ( 16 ) =   14.710
xnslope        ( 16 ) =   23.150
      
xxnpmax        ( 17 ) =    1.000
xq10soil       ( 17 ) =    1.720
xxkoptlitter   ( 17 ) =    0.400
xxkoptsoil     ( 17 ) =    1.000
xprodptase     ( 17 ) =    0.500
xcostnpup      ( 17 ) =   40.000
xmaxfinelitter ( 17 ) =  100.000
xmaxcwd        ( 17 ) =  100.000
xnintercept    ( 17 ) =   14.710
xnslope        ( 17 ) =   23.150
      





    DO nv=1,ntiles
       casabiome%ratioPcplantmin(nv,leaf)  = 1.0/(xratioNPleafmax*ratioCNplant(nv,leaf))
       casabiome%ratioPcplantmax(nv,leaf)  = 1.0/(xratioNPleafmin*ratioCNplant(nv,leaf))
       casabiome%ratioPcplantmin(nv,wood)  = 1.0/(xratioNPwoodmax*ratioCNplant(nv,wood))
       casabiome%ratioPcplantmax(nv,wood)  = 1.0/(xratioNPwoodmin*ratioCNplant(nv,wood))
       casabiome%ratioPcplantmin(nv,froot) = 1.0/(xratioNPfrootmax*ratioCNplant(nv,froot))
       casabiome%ratioPcplantmax(nv,froot) = 1.0/(xratioNPfrootmin*ratioCNplant(nv,froot))


       casabiome%ratioNPplantmin(nv,leaf)  = xratioNPleafmin
       casabiome%ratioNPplantmax(nv,leaf)  = xratioNPleafmax
       casabiome%ratioNPplantmin(nv,wood)  = xratioNPwoodmin
       casabiome%ratioNPplantmax(nv,wood)  = xratioNPwoodmax
       casabiome%ratioNPplantmin(nv,froot) = xratioNPfrootmin
       casabiome%ratioNPplantmax(nv,froot) = xratioNPfrootmax

    ENDDO

    fracroot   = 0.0
    depthsoila = 0.0
    depthsoilb = 0.0
    DO ns=1, nsl
       depthsoilb(ns) = depthsoilb(ns) + soil_zse(ns)
       IF (ns==1) THEN
          depthsoila(ns) = 0.0
       ELSE
          depthsoila(ns) = depthsoilb(ns-1)
       ENDIF
    ENDDO

    DO nv=1,ntiles
       casabiome%sla(nv)             = slax(nv)
       casabiome%fraclabile(nv,leaf) = deltcasa*0.6    !1/day
       casabiome%fraclabile(nv,froot)= deltcasa*0.4    !1/day
       casabiome%fraclabile(nv,wood) = deltcasa*0.0
       casabiome%plantrate(nv,leaf)  = deltcasa/(leafage(nv)*(1.0-xfherbivore(nv)))
       casabiome%plantrate(nv,froot) = deltcasa/frootage(nv)
       casabiome%plantrate(nv,wood)  = deltcasa/woodage(nv)
       casabiome%litterrate(nv,metb) = deltcasa/metage(nv)
       casabiome%litterrate(nv,str)  = deltcasa/strage(nv)
       casabiome%litterrate(nv,cwd)  = deltcasa/cwdage(nv)
       casabiome%soilrate(nv,mic)    = deltcasa/micage(nv)
       casabiome%soilrate(nv,slow)   = deltcasa/slowage(nv)
       casabiome%soilrate(nv,pass)   = deltcasa/passage(nv)
       casabiome%xkleafcoldmax(nv)   = deltcasa * xxkleafcoldmax(nv)
       casabiome%xkleafdrymax(nv)    = deltcasa * xxkleafdrymax(nv)
       !    casabiome%kuplabp(nv)         = xkuplabp(nv)
       casabiome%rmplant(nv,:)       = casabiome%rmplant(nv,:)*deltcasa
       casabiome%kclabrate(nv)       = deltcasa/clabileage(nv)

       !@@@@@@@@@@@@@@@@@
       casabiome%xnpmax(nv)          = xxnpmax(nv)
       casabiome%q10soil(nv)         = xq10soil(nv)
       casabiome%xkoptlitter(nv)     = xxkoptlitter(nv)
       casabiome%xkoptsoil(nv)       = xxkoptsoil(nv)
       casabiome%prodptase(nv)       = xprodptase(nv)/365.0   ! convert from yearly to daily
       casabiome%costnpup(nv)        = xcostnpup(nv)
       casabiome%maxfinelitter(nv)   = xmaxfinelitter(nv)
       casabiome%maxcwd(nv)          = xmaxcwd(nv)
       casabiome%nintercept(nv)      = xnintercept(nv)
       casabiome%nslope(nv)          = xnslope(nv)
       !@@@@@@@@@@@@@@
    ENDDO

    !@@@@@@@@@@@@@@
    DO ns=1,mso
       casabiome%xkplab(ns)          =  xxkplab(ns)
       casabiome%xkpsorb(ns)         =  xxkpsorb(ns)
       casabiome%xkpocc(ns)          =  xxkpocc(ns)
    ENDDO

    !@@@@@@@@@@@@@@

    ! PRINT *, 'casabiome%xkoptsoil = ', casabiome%xkoptsoil(2)

    DO npt = 1, mp
       iv1=veg_iveg(npt)
       iso=casamet%isorder(npt)
       ! The following to be commented out when coupled to CABLE
       !    veg%froot(npt,:) =fracroot(iv1,:)
       !    PRINT *, 'npt,iv1,iso ', npt,iv1, iso
       casamet%iveg2(npt) =casabiome%ivt2(iv1)
       casamet%lnonwood(npt) = 1
       casapool%cplant(npt,wood)  = 0.0
       casapool%clitter(npt,cwd)  = 0.0
       casapool%nplant(npt,wood)  = 0.0
       casapool%nlitter(npt,cwd)  = 0.0
       casapool%pplant(npt,wood)  = 0.0
       casapool%plitter(npt,cwd)  = 0.0
       IF (casamet%iveg2(npt)==forest.OR.casamet%iveg2(npt)==shrub) THEN
          casamet%lnonwood(npt) = 0
          casapool%cplant(npt,wood)  = Cwood(iv1)
          casapool%clitter(npt,cwd)  = ccwd(iv1)
          casapool%nplant(npt,wood)  = nwood(iv1)
          casapool%nlitter(npt,cwd)  = ncwd(iv1)
          casapool%pplant(npt,wood)  = xpwood(iv1)
          casapool%plitter(npt,cwd)  = xpcwd(iv1)
          !! vh_js !!
          IF (cable_user%CALL_POP) THEN  ! initialise very small wood pool, so POP can start from zero.
             casapool%cplant(npt,wood) = 0.01
             casapool%nplant(npt,wood)= casabiome%ratioNCplantmin(iv1,wood)* casapool%cplant(npt,wood)
             casapool%pplant(npt,wood)= casabiome%ratioPCplantmin(iv1,wood)* casapool%cplant(npt,wood)
          ENDIF
          !! vh_js

       ENDIF
       casapool%cplant(npt,leaf)     = cleaf(iv1)
       casapool%cplant(npt,froot)    = cfroot(iv1)
       casapool%clabile(npt)         = 0.0
       casapool%clitter(npt,metb)     = cmet(iv1)
       casapool%clitter(npt,str)     = cstr(iv1)
       casapool%csoil(npt,mic)       = cmic(iv1)
       casapool%csoil(npt,slow)      = cslow(iv1)
       casapool%csoil(npt,pass)      = cpass(iv1)
       IF (icycle==1) THEN
          casapool%ratioNCplant(npt,:)  = 1.0/ratioCNplant(iv1,:)
       ENDIF

       ! initializing glai in case not reading pool file (eg. during spin)
       casamet%glai(npt) = MAX(casabiome%glaimin(iv1), &
            casabiome%sla(iv1) * casapool%cplant(npt,leaf))

       casaflux%fNminloss(npt)   = xfNminloss(iv1)
       ! comment out by ypw 12/07/2009
       casaflux%fNminleach(npt)  = 10.0*xfNminleach(iv1) * deltcasa
       !    casaflux%fNminleach(npt)  = xfNminleach(iv1)
       casapool%nplant(npt,leaf) = nleaf(iv1)
       casapool%nplant(npt,froot)= nfroot(iv1)
       casapool%nlitter(npt,metb) = nmet(iv1)
       !    casapool%nlitter(npt,str) = nstr(iv1)
       casapool%nlitter(npt,str) = cstr(iv1)*ratioNCstrfix
       casapool%nsoil(npt,mic)   = nmic(iv1)
       casapool%nsoil(npt,slow)  = nslow(iv1)
       casapool%nsoil(npt,pass)  = npass(iv1)
       casapool%nsoilmin(npt)    = xnsoilmin(iv1)
       casapool%pplant(npt,leaf) = xpleaf(iv1)
       casapool%pplant(npt,froot)= xpfroot(iv1)
       casapool%plitter(npt,metb) = xpmet(iv1)
       !    casapool%plitter(npt,str) = xpstr(iv1)
       casapool%plitter(npt,str) = casapool%nlitter(npt,str)/ratioNPstrfix
       casapool%psoil(npt,mic)   = xpmic(iv1)
       casapool%psoil(npt,slow)  = xpslow(iv1)
       casapool%psoil(npt,pass)  = xppass(iv1)
       casapool%psoillab(npt)    = xplab(iv1)
       casapool%psoilsorb(npt)   = xpsorb(iv1)
       casapool%psoilocc(npt)    = xpocc(iv1)
       casaflux%kmlabp(npt)      = xkmlabp(iso)
       casaflux%psorbmax(npt)    = xpsorbmax(iso)
       casaflux%fpleach(npt)     = xfPleach(iso) 
       !   we used the spatially explicit estimate N fixation by Wang and Houlton (GRL)
       !    casaflux%Nminfix(npt)     = xnfixrate(iv1)/365.0

       ! use the PFT-specific C:N:P stoichiometry
       casapool%ratioNCplant(npt,:)  = 1.0/ratioCNplant(iv1,:)
       casapool%ratioNPplant(npt,:)  = casabiome%ratioNPplantmin(iv1,:)
       casapool%ratioPCplant(npt,:)  = 1.0/(ratioCNplant(iv1,:) *casabiome%ratioNPplantmin(iv1,:) )

       casapool%ratioNClitter(npt,metb) = casapool%ratioNCplant(npt,leaf)  * casabiome%ftransNPtoL(iv1,leaf)
       casapool%ratioNClitter(npt,str)  = casapool%ratioNCplant(npt,froot) * casabiome%ftransNPtoL(iv1,froot)
       casapool%ratioNClitter(npt,cwd)  = casapool%ratioNCplant(npt,wood)  * casabiome%ftransNPtoL(iv1,wood)

       casapool%ratioPClitter(npt,metb) = casapool%ratioPCplant(npt,leaf)  * casabiome%ftransPPtoL(iv1,leaf)
       casapool%ratioPClitter(npt,str)  = casapool%ratioPCplant(npt,froot) * casabiome%ftransPPtoL(iv1,froot)
       casapool%ratioPClitter(npt,cwd)  = casapool%ratioPCplant(npt,wood)  * casabiome%ftransPPtoL(iv1,wood)

       casapool%ratioNPlitter(npt,metb) = casapool%ratioNClitter(npt,metb)/(casapool%ratioPClitter(npt,metb) +1.0e-10)
       casapool%ratioNPlitter(npt,str)  = casapool%ratioNClitter(npt,str)/(casapool%ratioPClitter(npt,str) +1.0e-10)
       casapool%ratioNPlitter(npt,cwd)  = casapool%ratioNClitter(npt,cwd)/(casapool%ratioPClitter(npt,cwd) +1.0e-10)

       casapool%ratioNCsoil(npt,:)   = 1.0/ratioCNsoil(iv1,:)
       casapool%ratioNPsoil(npt,:)   = ratioNPsoil(iso,:)
       casapool%ratioPCsoil(npt,:)   = 1.0/(ratioCNsoil(iv1,:)*ratioNPsoil(iso,:))

       casapool%ratioNCsoilmin(npt,:)   = 1.0/ratioCNsoilmax(iv1,:)
       casapool%ratioNCsoilmax(npt,:)   = 1.0/ratioCNsoilmin(iv1,:)
       casapool%ratioNCsoilnew(npt,:)   = casapool%ratioNCsoilmax(npt,:)
    ENDDO

    IF(icycle<2) THEN
       casapool%Nplant(:,:)  = casapool%Cplant(:,:)  * casapool%ratioNCplant(:,:)
       casapool%Nlitter(:,:) = casapool%Clitter(:,:) * casapool%ratioNClitter(:,:)
       casapool%Nsoil(:,:)   = casapool%Csoil(:,:)   * casapool%ratioNCsoil(:,:) 
    ENDIF
    IF(icycle<3) THEN
       casapool%Pplant(:,:)  = casapool%Cplant(:,:)  * casapool%ratioPCplant(:,:)
       casapool%Plitter(:,:) = casapool%Clitter(:,:) * casapool%ratioPClitter(:,:)
       casapool%Psoil(:,:)   = casapool%Csoil(:,:)   * casapool%ratioPCsoil(:,:) 
       casapool%Psoilsorb(:) = casaflux%psorbmax(:) * casapool%psoillab(:) &
            /(casaflux%kmlabp(:)+casapool%psoillab(:))
    ENDIF

RETURN
END SUBROUTINE casa_readbiome

END MODULE casa_readbiome_module
