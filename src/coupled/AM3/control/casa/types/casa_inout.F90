!==============================================================================
! This source code is part of the
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CSIRO Open Source Software License
! Agreement (variation of the BSD / MIT License).
!
! You may not use this file except in compliance with this License.
! A copy of the License (CSIRO_BSD_MIT_License_v2.0_CABLE.txt) is located
! in each directory containing CABLE code.
!
! ==============================================================================
! Purpose: Input and output code for CASA-CNP when run offline
!          ACCESS version may use some of this code but split into different files?
!
! Contact: Yingping.Wang@csiro.au and Bernard.Pak@csiro.au
!
! History: Developed for offline code.  Expect to re-write for MPI and ACCESS
!          versions
!
!
! ==============================================================================
! casa_inout.f90
!
! the following routines are used when "casacnp" is coupled to "cable"
!   casa_readbiome
!   casa_readphen
!   casa_readpoint   (removed, now done in parameter_module)
!   casa_init
!   casa_poolout
!   casa_cnpflux  (zeros casabal quantites on doy 1 and updates casabal at end of biogeochem)
!   biogeochem

MODULE casa_inout_module

USE casa_files_type_mod, ONLY: casafile   ! TYPE instance 

CONTAINS

  SUBROUTINE casa_readphen(veg,casamet,phen)
    !SUBROUTINE casa_readphen(mvt,veg,casamet,phen)
    ! read in the tabulated modis-derived leaf phenology data
    ! for latitude bands of 79.75 to -55.25
    USE cable_def_types_mod
    USE casadimension
    USE casaparm
    USE casavariable
USE phenology_type_mod,          ONLY: phenology_type 

IMPLICIT NONE

    TYPE (veg_parameter_type), INTENT(IN)    :: veg  ! vegetation parameters
    TYPE (casa_met),           INTENT(IN)    :: casamet
TYPE (phenology_type),      INTENT(INOUT) :: phen

    ! local variables
    INTEGER, PARAMETER            :: nphen=8! was 10(IGBP). changed by Q.Zhang @01/12/2011
    INTEGER np,nx,ilat
    INTEGER, DIMENSION(271,mvtype) :: greenup, fall,  phendoy1
    INTEGER, DIMENSION(nphen)     :: greenupx,fallx,xphendoy1
    INTEGER, DIMENSION(271,nphen)     :: greenupxy,fallxy,xyphendoy1
    INTEGER, DIMENSION(nphen)     :: ivtx
    REAL(r_2), DIMENSION(271)     :: xlat

    ! initilize for evergreen PFTs
    greenup(:,:) = -50
    fall(:,:)    = 367
    phendoy1(:,:)= 2

    !jh!OPEN(101,file=casafile%phen)
    !jh!READ(101,*)
    !jh!READ(101,*) (ivtx(nx),nx=1,nphen) ! fixed at 10, as only 10 of 17 IGBP PFT
    !jh!! have seasonal leaf phenology
    !jh!DO ilat=271,1,-1
    !jh!   READ(101,*) xlat(ilat),(greenupx(nx),nx=1,nphen), &
    !jh!        (fallx(nx),nx=1,nphen),(xphendoy1(nx),nx=1,nphen)
    !jh!   DO nx=1,nphen
    !jh!      greenup(ilat,ivtx(nx)) = greenupx(nx)
    !jh!      fall(ilat,ivtx(nx))    = fallx(nx)
    !jh!      phendoy1(ilat,ivtx(nx))= xphendoy1(nx)
    !jh!   ENDDO
    !jh!ENDDO

 ! xlat         79.7500000000000     
greenupxy  ( 271, 1 ) =  167
greenupxy  ( 271, 2 ) =  156
greenupxy  ( 271, 3 ) =  184
greenupxy  ( 271, 4 ) =  184
greenupxy  ( 271, 5 ) =  184
greenupxy  ( 271, 6 ) =  184
greenupxy  ( 271, 7 ) =  145
greenupxy  ( 271, 8 ) =  145
fallxy     ( 271, 1 ) =  231
fallxy     ( 271, 2 ) =  231
fallxy     ( 271, 3 ) =  229
fallxy     ( 271, 4 ) =  229
fallxy     ( 271, 5 ) =  229
fallxy     ( 271, 6 ) =  229
fallxy     ( 271, 7 ) =  230
fallxy     ( 271, 8 ) =  230
xyphendoy1 ( 271, 1 ) =    0
xyphendoy1 ( 271, 2 ) =    0
xyphendoy1 ( 271, 3 ) =    0
xyphendoy1 ( 271, 4 ) =    0
xyphendoy1 ( 271, 5 ) =    0
xyphendoy1 ( 271, 6 ) =    0
xyphendoy1 ( 271, 7 ) =    0
xyphendoy1 ( 271, 8 ) =    0
      
 ! xlat         79.2500000000000     
greenupxy  ( 270, 1 ) =  167
greenupxy  ( 270, 2 ) =  156
greenupxy  ( 270, 3 ) =  184
greenupxy  ( 270, 4 ) =  184
greenupxy  ( 270, 5 ) =  184
greenupxy  ( 270, 6 ) =  184
greenupxy  ( 270, 7 ) =  145
greenupxy  ( 270, 8 ) =  145
fallxy     ( 270, 1 ) =  231
fallxy     ( 270, 2 ) =  231
fallxy     ( 270, 3 ) =  229
fallxy     ( 270, 4 ) =  229
fallxy     ( 270, 5 ) =  229
fallxy     ( 270, 6 ) =  229
fallxy     ( 270, 7 ) =  230
fallxy     ( 270, 8 ) =  230
xyphendoy1 ( 270, 1 ) =    0
xyphendoy1 ( 270, 2 ) =    0
xyphendoy1 ( 270, 3 ) =    0
xyphendoy1 ( 270, 4 ) =    0
xyphendoy1 ( 270, 5 ) =    0
xyphendoy1 ( 270, 6 ) =    0
xyphendoy1 ( 270, 7 ) =    0
xyphendoy1 ( 270, 8 ) =    0
      
 ! xlat         78.7500000000000     
greenupxy  ( 269, 1 ) =  167
greenupxy  ( 269, 2 ) =  156
greenupxy  ( 269, 3 ) =  184
greenupxy  ( 269, 4 ) =  184
greenupxy  ( 269, 5 ) =  184
greenupxy  ( 269, 6 ) =  184
greenupxy  ( 269, 7 ) =  145
greenupxy  ( 269, 8 ) =  145
fallxy     ( 269, 1 ) =  231
fallxy     ( 269, 2 ) =  231
fallxy     ( 269, 3 ) =  229
fallxy     ( 269, 4 ) =  229
fallxy     ( 269, 5 ) =  229
fallxy     ( 269, 6 ) =  229
fallxy     ( 269, 7 ) =  230
fallxy     ( 269, 8 ) =  230
xyphendoy1 ( 269, 1 ) =    0
xyphendoy1 ( 269, 2 ) =    0
xyphendoy1 ( 269, 3 ) =    0
xyphendoy1 ( 269, 4 ) =    0
xyphendoy1 ( 269, 5 ) =    0
xyphendoy1 ( 269, 6 ) =    0
xyphendoy1 ( 269, 7 ) =    0
xyphendoy1 ( 269, 8 ) =    0
      
 ! xlat         78.2500000000000     
greenupxy  ( 268, 1 ) =  167
greenupxy  ( 268, 2 ) =  156
greenupxy  ( 268, 3 ) =  184
greenupxy  ( 268, 4 ) =  184
greenupxy  ( 268, 5 ) =  184
greenupxy  ( 268, 6 ) =  184
greenupxy  ( 268, 7 ) =  145
greenupxy  ( 268, 8 ) =  145
fallxy     ( 268, 1 ) =  231
fallxy     ( 268, 2 ) =  231
fallxy     ( 268, 3 ) =  229
fallxy     ( 268, 4 ) =  229
fallxy     ( 268, 5 ) =  229
fallxy     ( 268, 6 ) =  229
fallxy     ( 268, 7 ) =  230
fallxy     ( 268, 8 ) =  230
xyphendoy1 ( 268, 1 ) =    0
xyphendoy1 ( 268, 2 ) =    0
xyphendoy1 ( 268, 3 ) =    0
xyphendoy1 ( 268, 4 ) =    0
xyphendoy1 ( 268, 5 ) =    0
xyphendoy1 ( 268, 6 ) =    0
xyphendoy1 ( 268, 7 ) =    0
xyphendoy1 ( 268, 8 ) =    0
      
 ! xlat         77.7500000000000     
greenupxy  ( 267, 1 ) =  167
greenupxy  ( 267, 2 ) =  156
greenupxy  ( 267, 3 ) =  184
greenupxy  ( 267, 4 ) =  184
greenupxy  ( 267, 5 ) =  184
greenupxy  ( 267, 6 ) =  184
greenupxy  ( 267, 7 ) =  145
greenupxy  ( 267, 8 ) =  145
fallxy     ( 267, 1 ) =  231
fallxy     ( 267, 2 ) =  231
fallxy     ( 267, 3 ) =  229
fallxy     ( 267, 4 ) =  229
fallxy     ( 267, 5 ) =  229
fallxy     ( 267, 6 ) =  229
fallxy     ( 267, 7 ) =  230
fallxy     ( 267, 8 ) =  230
xyphendoy1 ( 267, 1 ) =    0
xyphendoy1 ( 267, 2 ) =    0
xyphendoy1 ( 267, 3 ) =    0
xyphendoy1 ( 267, 4 ) =    0
xyphendoy1 ( 267, 5 ) =    0
xyphendoy1 ( 267, 6 ) =    0
xyphendoy1 ( 267, 7 ) =    0
xyphendoy1 ( 267, 8 ) =    0
      
 ! xlat         77.2500000000000     
greenupxy  ( 266, 1 ) =  167
greenupxy  ( 266, 2 ) =  156
greenupxy  ( 266, 3 ) =  184
greenupxy  ( 266, 4 ) =  184
greenupxy  ( 266, 5 ) =  184
greenupxy  ( 266, 6 ) =  184
greenupxy  ( 266, 7 ) =  145
greenupxy  ( 266, 8 ) =  145
fallxy     ( 266, 1 ) =  231
fallxy     ( 266, 2 ) =  231
fallxy     ( 266, 3 ) =  229
fallxy     ( 266, 4 ) =  229
fallxy     ( 266, 5 ) =  229
fallxy     ( 266, 6 ) =  229
fallxy     ( 266, 7 ) =  230
fallxy     ( 266, 8 ) =  230
xyphendoy1 ( 266, 1 ) =    0
xyphendoy1 ( 266, 2 ) =    0
xyphendoy1 ( 266, 3 ) =    0
xyphendoy1 ( 266, 4 ) =    0
xyphendoy1 ( 266, 5 ) =    0
xyphendoy1 ( 266, 6 ) =    0
xyphendoy1 ( 266, 7 ) =    0
xyphendoy1 ( 266, 8 ) =    0
      
 ! xlat         76.7500000000000     
greenupxy  ( 265, 1 ) =  167
greenupxy  ( 265, 2 ) =  156
greenupxy  ( 265, 3 ) =  184
greenupxy  ( 265, 4 ) =  184
greenupxy  ( 265, 5 ) =  184
greenupxy  ( 265, 6 ) =  184
greenupxy  ( 265, 7 ) =  145
greenupxy  ( 265, 8 ) =  145
fallxy     ( 265, 1 ) =  231
fallxy     ( 265, 2 ) =  231
fallxy     ( 265, 3 ) =  229
fallxy     ( 265, 4 ) =  230
fallxy     ( 265, 5 ) =  230
fallxy     ( 265, 6 ) =  230
fallxy     ( 265, 7 ) =  230
fallxy     ( 265, 8 ) =  230
xyphendoy1 ( 265, 1 ) =    0
xyphendoy1 ( 265, 2 ) =    0
xyphendoy1 ( 265, 3 ) =    0
xyphendoy1 ( 265, 4 ) =    0
xyphendoy1 ( 265, 5 ) =    0
xyphendoy1 ( 265, 6 ) =    0
xyphendoy1 ( 265, 7 ) =    0
xyphendoy1 ( 265, 8 ) =    0
      
 ! xlat         76.2500000000000     
greenupxy  ( 264, 1 ) =  167
greenupxy  ( 264, 2 ) =  156
greenupxy  ( 264, 3 ) =  184
greenupxy  ( 264, 4 ) =  184
greenupxy  ( 264, 5 ) =  184
greenupxy  ( 264, 6 ) =  184
greenupxy  ( 264, 7 ) =  145
greenupxy  ( 264, 8 ) =  145
fallxy     ( 264, 1 ) =  231
fallxy     ( 264, 2 ) =  231
fallxy     ( 264, 3 ) =  230
fallxy     ( 264, 4 ) =  230
fallxy     ( 264, 5 ) =  230
fallxy     ( 264, 6 ) =  230
fallxy     ( 264, 7 ) =  230
fallxy     ( 264, 8 ) =  230
xyphendoy1 ( 264, 1 ) =    0
xyphendoy1 ( 264, 2 ) =    0
xyphendoy1 ( 264, 3 ) =    0
xyphendoy1 ( 264, 4 ) =    0
xyphendoy1 ( 264, 5 ) =    0
xyphendoy1 ( 264, 6 ) =    0
xyphendoy1 ( 264, 7 ) =    0
xyphendoy1 ( 264, 8 ) =    0
      
 ! xlat         75.7500000000000     
greenupxy  ( 263, 1 ) =  167
greenupxy  ( 263, 2 ) =  156
greenupxy  ( 263, 3 ) =  184
greenupxy  ( 263, 4 ) =  184
greenupxy  ( 263, 5 ) =  184
greenupxy  ( 263, 6 ) =  184
greenupxy  ( 263, 7 ) =  145
greenupxy  ( 263, 8 ) =  145
fallxy     ( 263, 1 ) =  231
fallxy     ( 263, 2 ) =  231
fallxy     ( 263, 3 ) =  230
fallxy     ( 263, 4 ) =  231
fallxy     ( 263, 5 ) =  231
fallxy     ( 263, 6 ) =  231
fallxy     ( 263, 7 ) =  230
fallxy     ( 263, 8 ) =  230
xyphendoy1 ( 263, 1 ) =    0
xyphendoy1 ( 263, 2 ) =    0
xyphendoy1 ( 263, 3 ) =    0
xyphendoy1 ( 263, 4 ) =    0
xyphendoy1 ( 263, 5 ) =    0
xyphendoy1 ( 263, 6 ) =    0
xyphendoy1 ( 263, 7 ) =    0
xyphendoy1 ( 263, 8 ) =    0
      
 ! xlat         75.2500000000000     
greenupxy  ( 262, 1 ) =  167
greenupxy  ( 262, 2 ) =  156
greenupxy  ( 262, 3 ) =  185
greenupxy  ( 262, 4 ) =  185
greenupxy  ( 262, 5 ) =  185
greenupxy  ( 262, 6 ) =  185
greenupxy  ( 262, 7 ) =  145
greenupxy  ( 262, 8 ) =  145
fallxy     ( 262, 1 ) =  231
fallxy     ( 262, 2 ) =  231
fallxy     ( 262, 3 ) =  231
fallxy     ( 262, 4 ) =  232
fallxy     ( 262, 5 ) =  232
fallxy     ( 262, 6 ) =  232
fallxy     ( 262, 7 ) =  230
fallxy     ( 262, 8 ) =  230
xyphendoy1 ( 262, 1 ) =    0
xyphendoy1 ( 262, 2 ) =    0
xyphendoy1 ( 262, 3 ) =    0
xyphendoy1 ( 262, 4 ) =    0
xyphendoy1 ( 262, 5 ) =    0
xyphendoy1 ( 262, 6 ) =    0
xyphendoy1 ( 262, 7 ) =    0
xyphendoy1 ( 262, 8 ) =    0
      
 ! xlat         74.7500000000000     
greenupxy  ( 261, 1 ) =  167
greenupxy  ( 261, 2 ) =  156
greenupxy  ( 261, 3 ) =  185
greenupxy  ( 261, 4 ) =  185
greenupxy  ( 261, 5 ) =  185
greenupxy  ( 261, 6 ) =  185
greenupxy  ( 261, 7 ) =  145
greenupxy  ( 261, 8 ) =  145
fallxy     ( 261, 1 ) =  231
fallxy     ( 261, 2 ) =  231
fallxy     ( 261, 3 ) =  231
fallxy     ( 261, 4 ) =  232
fallxy     ( 261, 5 ) =  232
fallxy     ( 261, 6 ) =  232
fallxy     ( 261, 7 ) =  230
fallxy     ( 261, 8 ) =  230
xyphendoy1 ( 261, 1 ) =    0
xyphendoy1 ( 261, 2 ) =    0
xyphendoy1 ( 261, 3 ) =    0
xyphendoy1 ( 261, 4 ) =    0
xyphendoy1 ( 261, 5 ) =    0
xyphendoy1 ( 261, 6 ) =    0
xyphendoy1 ( 261, 7 ) =    0
xyphendoy1 ( 261, 8 ) =    0
      
 ! xlat         74.2500000000000     
greenupxy  ( 260, 1 ) =  167
greenupxy  ( 260, 2 ) =  156
greenupxy  ( 260, 3 ) =  185
greenupxy  ( 260, 4 ) =  184
greenupxy  ( 260, 5 ) =  184
greenupxy  ( 260, 6 ) =  184
greenupxy  ( 260, 7 ) =  145
greenupxy  ( 260, 8 ) =  145
fallxy     ( 260, 1 ) =  231
fallxy     ( 260, 2 ) =  231
fallxy     ( 260, 3 ) =  232
fallxy     ( 260, 4 ) =  232
fallxy     ( 260, 5 ) =  232
fallxy     ( 260, 6 ) =  232
fallxy     ( 260, 7 ) =  230
fallxy     ( 260, 8 ) =  230
xyphendoy1 ( 260, 1 ) =    0
xyphendoy1 ( 260, 2 ) =    0
xyphendoy1 ( 260, 3 ) =    0
xyphendoy1 ( 260, 4 ) =    0
xyphendoy1 ( 260, 5 ) =    0
xyphendoy1 ( 260, 6 ) =    0
xyphendoy1 ( 260, 7 ) =    0
xyphendoy1 ( 260, 8 ) =    0
      
 ! xlat         73.7500000000000     
greenupxy  ( 259, 1 ) =  167
greenupxy  ( 259, 2 ) =  156
greenupxy  ( 259, 3 ) =  185
greenupxy  ( 259, 4 ) =  184
greenupxy  ( 259, 5 ) =  184
greenupxy  ( 259, 6 ) =  184
greenupxy  ( 259, 7 ) =  145
greenupxy  ( 259, 8 ) =  145
fallxy     ( 259, 1 ) =  231
fallxy     ( 259, 2 ) =  231
fallxy     ( 259, 3 ) =  233
fallxy     ( 259, 4 ) =  233
fallxy     ( 259, 5 ) =  233
fallxy     ( 259, 6 ) =  233
fallxy     ( 259, 7 ) =  230
fallxy     ( 259, 8 ) =  230
xyphendoy1 ( 259, 1 ) =    0
xyphendoy1 ( 259, 2 ) =    0
xyphendoy1 ( 259, 3 ) =    0
xyphendoy1 ( 259, 4 ) =    0
xyphendoy1 ( 259, 5 ) =    0
xyphendoy1 ( 259, 6 ) =    0
xyphendoy1 ( 259, 7 ) =    0
xyphendoy1 ( 259, 8 ) =    0
      
 ! xlat         73.2500000000000     
greenupxy  ( 258, 1 ) =  167
greenupxy  ( 258, 2 ) =  156
greenupxy  ( 258, 3 ) =  186
greenupxy  ( 258, 4 ) =  183
greenupxy  ( 258, 5 ) =  183
greenupxy  ( 258, 6 ) =  183
greenupxy  ( 258, 7 ) =  145
greenupxy  ( 258, 8 ) =  145
fallxy     ( 258, 1 ) =  231
fallxy     ( 258, 2 ) =  231
fallxy     ( 258, 3 ) =  232
fallxy     ( 258, 4 ) =  233
fallxy     ( 258, 5 ) =  233
fallxy     ( 258, 6 ) =  233
fallxy     ( 258, 7 ) =  230
fallxy     ( 258, 8 ) =  230
xyphendoy1 ( 258, 1 ) =    0
xyphendoy1 ( 258, 2 ) =    0
xyphendoy1 ( 258, 3 ) =    0
xyphendoy1 ( 258, 4 ) =    0
xyphendoy1 ( 258, 5 ) =    0
xyphendoy1 ( 258, 6 ) =    0
xyphendoy1 ( 258, 7 ) =    0
xyphendoy1 ( 258, 8 ) =    0
      
 ! xlat         72.7500000000000     
greenupxy  ( 257, 1 ) =  167
greenupxy  ( 257, 2 ) =  156
greenupxy  ( 257, 3 ) =  186
greenupxy  ( 257, 4 ) =  184
greenupxy  ( 257, 5 ) =  184
greenupxy  ( 257, 6 ) =  184
greenupxy  ( 257, 7 ) =  145
greenupxy  ( 257, 8 ) =  145
fallxy     ( 257, 1 ) =  231
fallxy     ( 257, 2 ) =  231
fallxy     ( 257, 3 ) =  233
fallxy     ( 257, 4 ) =  234
fallxy     ( 257, 5 ) =  234
fallxy     ( 257, 6 ) =  234
fallxy     ( 257, 7 ) =  230
fallxy     ( 257, 8 ) =  230
xyphendoy1 ( 257, 1 ) =    0
xyphendoy1 ( 257, 2 ) =    0
xyphendoy1 ( 257, 3 ) =    0
xyphendoy1 ( 257, 4 ) =    0
xyphendoy1 ( 257, 5 ) =    0
xyphendoy1 ( 257, 6 ) =    0
xyphendoy1 ( 257, 7 ) =    0
xyphendoy1 ( 257, 8 ) =    0
      
 ! xlat         72.2500000000000     
greenupxy  ( 256, 1 ) =  167
greenupxy  ( 256, 2 ) =  156
greenupxy  ( 256, 3 ) =  186
greenupxy  ( 256, 4 ) =  183
greenupxy  ( 256, 5 ) =  183
greenupxy  ( 256, 6 ) =  183
greenupxy  ( 256, 7 ) =  145
greenupxy  ( 256, 8 ) =  145
fallxy     ( 256, 1 ) =  231
fallxy     ( 256, 2 ) =  231
fallxy     ( 256, 3 ) =  234
fallxy     ( 256, 4 ) =  234
fallxy     ( 256, 5 ) =  234
fallxy     ( 256, 6 ) =  234
fallxy     ( 256, 7 ) =  230
fallxy     ( 256, 8 ) =  230
xyphendoy1 ( 256, 1 ) =    0
xyphendoy1 ( 256, 2 ) =    0
xyphendoy1 ( 256, 3 ) =    0
xyphendoy1 ( 256, 4 ) =    0
xyphendoy1 ( 256, 5 ) =    0
xyphendoy1 ( 256, 6 ) =    0
xyphendoy1 ( 256, 7 ) =    0
xyphendoy1 ( 256, 8 ) =    0
      
 ! xlat         71.7500000000000     
greenupxy  ( 255, 1 ) =  167
greenupxy  ( 255, 2 ) =  156
greenupxy  ( 255, 3 ) =  185
greenupxy  ( 255, 4 ) =  183
greenupxy  ( 255, 5 ) =  183
greenupxy  ( 255, 6 ) =  183
greenupxy  ( 255, 7 ) =  145
greenupxy  ( 255, 8 ) =  145
fallxy     ( 255, 1 ) =  231
fallxy     ( 255, 2 ) =  231
fallxy     ( 255, 3 ) =  234
fallxy     ( 255, 4 ) =  234
fallxy     ( 255, 5 ) =  234
fallxy     ( 255, 6 ) =  234
fallxy     ( 255, 7 ) =  230
fallxy     ( 255, 8 ) =  230
xyphendoy1 ( 255, 1 ) =    0
xyphendoy1 ( 255, 2 ) =    0
xyphendoy1 ( 255, 3 ) =    0
xyphendoy1 ( 255, 4 ) =    0
xyphendoy1 ( 255, 5 ) =    0
xyphendoy1 ( 255, 6 ) =    0
xyphendoy1 ( 255, 7 ) =    0
xyphendoy1 ( 255, 8 ) =    0
      
 ! xlat         71.2500000000000     
greenupxy  ( 254, 1 ) =  167
greenupxy  ( 254, 2 ) =  156
greenupxy  ( 254, 3 ) =  184
greenupxy  ( 254, 4 ) =  182
greenupxy  ( 254, 5 ) =  182
greenupxy  ( 254, 6 ) =  182
greenupxy  ( 254, 7 ) =  145
greenupxy  ( 254, 8 ) =  145
fallxy     ( 254, 1 ) =  231
fallxy     ( 254, 2 ) =  231
fallxy     ( 254, 3 ) =  234
fallxy     ( 254, 4 ) =  234
fallxy     ( 254, 5 ) =  234
fallxy     ( 254, 6 ) =  234
fallxy     ( 254, 7 ) =  230
fallxy     ( 254, 8 ) =  230
xyphendoy1 ( 254, 1 ) =    0
xyphendoy1 ( 254, 2 ) =    0
xyphendoy1 ( 254, 3 ) =    0
xyphendoy1 ( 254, 4 ) =    0
xyphendoy1 ( 254, 5 ) =    0
xyphendoy1 ( 254, 6 ) =    0
xyphendoy1 ( 254, 7 ) =    0
xyphendoy1 ( 254, 8 ) =    0
      
 ! xlat         70.7500000000000     
greenupxy  ( 253, 1 ) =  167
greenupxy  ( 253, 2 ) =  156
greenupxy  ( 253, 3 ) =  183
greenupxy  ( 253, 4 ) =  181
greenupxy  ( 253, 5 ) =  181
greenupxy  ( 253, 6 ) =  181
greenupxy  ( 253, 7 ) =  145
greenupxy  ( 253, 8 ) =  145
fallxy     ( 253, 1 ) =  231
fallxy     ( 253, 2 ) =  231
fallxy     ( 253, 3 ) =  234
fallxy     ( 253, 4 ) =  233
fallxy     ( 253, 5 ) =  233
fallxy     ( 253, 6 ) =  233
fallxy     ( 253, 7 ) =  230
fallxy     ( 253, 8 ) =  230
xyphendoy1 ( 253, 1 ) =    0
xyphendoy1 ( 253, 2 ) =    0
xyphendoy1 ( 253, 3 ) =    0
xyphendoy1 ( 253, 4 ) =    0
xyphendoy1 ( 253, 5 ) =    0
xyphendoy1 ( 253, 6 ) =    0
xyphendoy1 ( 253, 7 ) =    0
xyphendoy1 ( 253, 8 ) =    0
      
 ! xlat         70.2500000000000     
greenupxy  ( 252, 1 ) =  167
greenupxy  ( 252, 2 ) =  156
greenupxy  ( 252, 3 ) =  182
greenupxy  ( 252, 4 ) =  180
greenupxy  ( 252, 5 ) =  180
greenupxy  ( 252, 6 ) =  180
greenupxy  ( 252, 7 ) =  145
greenupxy  ( 252, 8 ) =  145
fallxy     ( 252, 1 ) =  231
fallxy     ( 252, 2 ) =  231
fallxy     ( 252, 3 ) =  234
fallxy     ( 252, 4 ) =  234
fallxy     ( 252, 5 ) =  234
fallxy     ( 252, 6 ) =  234
fallxy     ( 252, 7 ) =  230
fallxy     ( 252, 8 ) =  230
xyphendoy1 ( 252, 1 ) =    0
xyphendoy1 ( 252, 2 ) =    0
xyphendoy1 ( 252, 3 ) =    0
xyphendoy1 ( 252, 4 ) =    0
xyphendoy1 ( 252, 5 ) =    0
xyphendoy1 ( 252, 6 ) =    0
xyphendoy1 ( 252, 7 ) =    0
xyphendoy1 ( 252, 8 ) =    0
      
 ! xlat         69.7500000000000     
greenupxy  ( 251, 1 ) =  167
greenupxy  ( 251, 2 ) =  156
greenupxy  ( 251, 3 ) =  181
greenupxy  ( 251, 4 ) =  178
greenupxy  ( 251, 5 ) =  178
greenupxy  ( 251, 6 ) =  178
greenupxy  ( 251, 7 ) =  145
greenupxy  ( 251, 8 ) =  145
fallxy     ( 251, 1 ) =  231
fallxy     ( 251, 2 ) =  231
fallxy     ( 251, 3 ) =  234
fallxy     ( 251, 4 ) =  234
fallxy     ( 251, 5 ) =  234
fallxy     ( 251, 6 ) =  234
fallxy     ( 251, 7 ) =  230
fallxy     ( 251, 8 ) =  230
xyphendoy1 ( 251, 1 ) =    0
xyphendoy1 ( 251, 2 ) =    0
xyphendoy1 ( 251, 3 ) =    0
xyphendoy1 ( 251, 4 ) =    0
xyphendoy1 ( 251, 5 ) =    0
xyphendoy1 ( 251, 6 ) =    0
xyphendoy1 ( 251, 7 ) =    0
xyphendoy1 ( 251, 8 ) =    0
      
 ! xlat         69.2500000000000     
greenupxy  ( 250, 1 ) =  167
greenupxy  ( 250, 2 ) =  156
greenupxy  ( 250, 3 ) =  179
greenupxy  ( 250, 4 ) =  177
greenupxy  ( 250, 5 ) =  177
greenupxy  ( 250, 6 ) =  177
greenupxy  ( 250, 7 ) =  145
greenupxy  ( 250, 8 ) =  145
fallxy     ( 250, 1 ) =  231
fallxy     ( 250, 2 ) =  231
fallxy     ( 250, 3 ) =  234
fallxy     ( 250, 4 ) =  234
fallxy     ( 250, 5 ) =  234
fallxy     ( 250, 6 ) =  234
fallxy     ( 250, 7 ) =  230
fallxy     ( 250, 8 ) =  230
xyphendoy1 ( 250, 1 ) =    0
xyphendoy1 ( 250, 2 ) =    0
xyphendoy1 ( 250, 3 ) =    0
xyphendoy1 ( 250, 4 ) =    0
xyphendoy1 ( 250, 5 ) =    0
xyphendoy1 ( 250, 6 ) =    0
xyphendoy1 ( 250, 7 ) =    0
xyphendoy1 ( 250, 8 ) =    0
      
 ! xlat         68.7500000000000     
greenupxy  ( 249, 1 ) =  167
greenupxy  ( 249, 2 ) =  156
greenupxy  ( 249, 3 ) =  177
greenupxy  ( 249, 4 ) =  176
greenupxy  ( 249, 5 ) =  176
greenupxy  ( 249, 6 ) =  176
greenupxy  ( 249, 7 ) =  145
greenupxy  ( 249, 8 ) =  145
fallxy     ( 249, 1 ) =  231
fallxy     ( 249, 2 ) =  231
fallxy     ( 249, 3 ) =  234
fallxy     ( 249, 4 ) =  233
fallxy     ( 249, 5 ) =  233
fallxy     ( 249, 6 ) =  233
fallxy     ( 249, 7 ) =  230
fallxy     ( 249, 8 ) =  230
xyphendoy1 ( 249, 1 ) =    0
xyphendoy1 ( 249, 2 ) =    0
xyphendoy1 ( 249, 3 ) =    0
xyphendoy1 ( 249, 4 ) =    0
xyphendoy1 ( 249, 5 ) =    0
xyphendoy1 ( 249, 6 ) =    0
xyphendoy1 ( 249, 7 ) =    0
xyphendoy1 ( 249, 8 ) =    0
      
 ! xlat         68.2500000000000     
greenupxy  ( 248, 1 ) =  167
greenupxy  ( 248, 2 ) =  156
greenupxy  ( 248, 3 ) =  176
greenupxy  ( 248, 4 ) =  174
greenupxy  ( 248, 5 ) =  174
greenupxy  ( 248, 6 ) =  174
greenupxy  ( 248, 7 ) =  145
greenupxy  ( 248, 8 ) =  145
fallxy     ( 248, 1 ) =  231
fallxy     ( 248, 2 ) =  231
fallxy     ( 248, 3 ) =  233
fallxy     ( 248, 4 ) =  233
fallxy     ( 248, 5 ) =  233
fallxy     ( 248, 6 ) =  233
fallxy     ( 248, 7 ) =  230
fallxy     ( 248, 8 ) =  230
xyphendoy1 ( 248, 1 ) =    0
xyphendoy1 ( 248, 2 ) =    0
xyphendoy1 ( 248, 3 ) =    0
xyphendoy1 ( 248, 4 ) =    0
xyphendoy1 ( 248, 5 ) =    0
xyphendoy1 ( 248, 6 ) =    0
xyphendoy1 ( 248, 7 ) =    0
xyphendoy1 ( 248, 8 ) =    0
      
 ! xlat         67.7500000000000     
greenupxy  ( 247, 1 ) =  167
greenupxy  ( 247, 2 ) =  156
greenupxy  ( 247, 3 ) =  175
greenupxy  ( 247, 4 ) =  173
greenupxy  ( 247, 5 ) =  173
greenupxy  ( 247, 6 ) =  173
greenupxy  ( 247, 7 ) =  145
greenupxy  ( 247, 8 ) =  145
fallxy     ( 247, 1 ) =  231
fallxy     ( 247, 2 ) =  231
fallxy     ( 247, 3 ) =  233
fallxy     ( 247, 4 ) =  233
fallxy     ( 247, 5 ) =  233
fallxy     ( 247, 6 ) =  233
fallxy     ( 247, 7 ) =  230
fallxy     ( 247, 8 ) =  230
xyphendoy1 ( 247, 1 ) =    0
xyphendoy1 ( 247, 2 ) =    0
xyphendoy1 ( 247, 3 ) =    0
xyphendoy1 ( 247, 4 ) =    0
xyphendoy1 ( 247, 5 ) =    0
xyphendoy1 ( 247, 6 ) =    0
xyphendoy1 ( 247, 7 ) =    0
xyphendoy1 ( 247, 8 ) =    0
      
 ! xlat         67.2500000000000     
greenupxy  ( 246, 1 ) =  167
greenupxy  ( 246, 2 ) =  156
greenupxy  ( 246, 3 ) =  173
greenupxy  ( 246, 4 ) =  171
greenupxy  ( 246, 5 ) =  171
greenupxy  ( 246, 6 ) =  171
greenupxy  ( 246, 7 ) =  145
greenupxy  ( 246, 8 ) =  145
fallxy     ( 246, 1 ) =  231
fallxy     ( 246, 2 ) =  231
fallxy     ( 246, 3 ) =  233
fallxy     ( 246, 4 ) =  233
fallxy     ( 246, 5 ) =  233
fallxy     ( 246, 6 ) =  233
fallxy     ( 246, 7 ) =  230
fallxy     ( 246, 8 ) =  230
xyphendoy1 ( 246, 1 ) =    0
xyphendoy1 ( 246, 2 ) =    0
xyphendoy1 ( 246, 3 ) =    0
xyphendoy1 ( 246, 4 ) =    0
xyphendoy1 ( 246, 5 ) =    0
xyphendoy1 ( 246, 6 ) =    0
xyphendoy1 ( 246, 7 ) =    0
xyphendoy1 ( 246, 8 ) =    0
      
 ! xlat         66.7500000000000     
greenupxy  ( 245, 1 ) =  166
greenupxy  ( 245, 2 ) =  156
greenupxy  ( 245, 3 ) =  172
greenupxy  ( 245, 4 ) =  170
greenupxy  ( 245, 5 ) =  170
greenupxy  ( 245, 6 ) =  170
greenupxy  ( 245, 7 ) =  145
greenupxy  ( 245, 8 ) =  145
fallxy     ( 245, 1 ) =  231
fallxy     ( 245, 2 ) =  230
fallxy     ( 245, 3 ) =  233
fallxy     ( 245, 4 ) =  233
fallxy     ( 245, 5 ) =  233
fallxy     ( 245, 6 ) =  233
fallxy     ( 245, 7 ) =  230
fallxy     ( 245, 8 ) =  230
xyphendoy1 ( 245, 1 ) =    0
xyphendoy1 ( 245, 2 ) =    0
xyphendoy1 ( 245, 3 ) =    0
xyphendoy1 ( 245, 4 ) =    0
xyphendoy1 ( 245, 5 ) =    0
xyphendoy1 ( 245, 6 ) =    0
xyphendoy1 ( 245, 7 ) =    0
xyphendoy1 ( 245, 8 ) =    0
      
 ! xlat         66.2500000000000     
greenupxy  ( 244, 1 ) =  166
greenupxy  ( 244, 2 ) =  155
greenupxy  ( 244, 3 ) =  171
greenupxy  ( 244, 4 ) =  168
greenupxy  ( 244, 5 ) =  168
greenupxy  ( 244, 6 ) =  168
greenupxy  ( 244, 7 ) =  145
greenupxy  ( 244, 8 ) =  145
fallxy     ( 244, 1 ) =  231
fallxy     ( 244, 2 ) =  230
fallxy     ( 244, 3 ) =  233
fallxy     ( 244, 4 ) =  233
fallxy     ( 244, 5 ) =  233
fallxy     ( 244, 6 ) =  233
fallxy     ( 244, 7 ) =  230
fallxy     ( 244, 8 ) =  230
xyphendoy1 ( 244, 1 ) =    0
xyphendoy1 ( 244, 2 ) =    0
xyphendoy1 ( 244, 3 ) =    0
xyphendoy1 ( 244, 4 ) =    0
xyphendoy1 ( 244, 5 ) =    0
xyphendoy1 ( 244, 6 ) =    0
xyphendoy1 ( 244, 7 ) =    0
xyphendoy1 ( 244, 8 ) =    0
      
 ! xlat         65.7500000000000     
greenupxy  ( 243, 1 ) =  164
greenupxy  ( 243, 2 ) =  154
greenupxy  ( 243, 3 ) =  170
greenupxy  ( 243, 4 ) =  167
greenupxy  ( 243, 5 ) =  167
greenupxy  ( 243, 6 ) =  167
greenupxy  ( 243, 7 ) =  145
greenupxy  ( 243, 8 ) =  145
fallxy     ( 243, 1 ) =  231
fallxy     ( 243, 2 ) =  230
fallxy     ( 243, 3 ) =  233
fallxy     ( 243, 4 ) =  233
fallxy     ( 243, 5 ) =  233
fallxy     ( 243, 6 ) =  233
fallxy     ( 243, 7 ) =  230
fallxy     ( 243, 8 ) =  230
xyphendoy1 ( 243, 1 ) =    0
xyphendoy1 ( 243, 2 ) =    0
xyphendoy1 ( 243, 3 ) =    0
xyphendoy1 ( 243, 4 ) =    0
xyphendoy1 ( 243, 5 ) =    0
xyphendoy1 ( 243, 6 ) =    0
xyphendoy1 ( 243, 7 ) =    0
xyphendoy1 ( 243, 8 ) =    0
      
 ! xlat         65.2500000000000     
greenupxy  ( 242, 1 ) =  163
greenupxy  ( 242, 2 ) =  153
greenupxy  ( 242, 3 ) =  169
greenupxy  ( 242, 4 ) =  166
greenupxy  ( 242, 5 ) =  166
greenupxy  ( 242, 6 ) =  166
greenupxy  ( 242, 7 ) =  145
greenupxy  ( 242, 8 ) =  145
fallxy     ( 242, 1 ) =  231
fallxy     ( 242, 2 ) =  230
fallxy     ( 242, 3 ) =  233
fallxy     ( 242, 4 ) =  233
fallxy     ( 242, 5 ) =  233
fallxy     ( 242, 6 ) =  233
fallxy     ( 242, 7 ) =  230
fallxy     ( 242, 8 ) =  230
xyphendoy1 ( 242, 1 ) =    0
xyphendoy1 ( 242, 2 ) =    0
xyphendoy1 ( 242, 3 ) =    0
xyphendoy1 ( 242, 4 ) =    0
xyphendoy1 ( 242, 5 ) =    0
xyphendoy1 ( 242, 6 ) =    0
xyphendoy1 ( 242, 7 ) =    0
xyphendoy1 ( 242, 8 ) =    0
      
 ! xlat         64.7500000000000     
greenupxy  ( 241, 1 ) =  161
greenupxy  ( 241, 2 ) =  153
greenupxy  ( 241, 3 ) =  168
greenupxy  ( 241, 4 ) =  165
greenupxy  ( 241, 5 ) =  165
greenupxy  ( 241, 6 ) =  165
greenupxy  ( 241, 7 ) =  145
greenupxy  ( 241, 8 ) =  145
fallxy     ( 241, 1 ) =  231
fallxy     ( 241, 2 ) =  230
fallxy     ( 241, 3 ) =  233
fallxy     ( 241, 4 ) =  232
fallxy     ( 241, 5 ) =  232
fallxy     ( 241, 6 ) =  232
fallxy     ( 241, 7 ) =  230
fallxy     ( 241, 8 ) =  230
xyphendoy1 ( 241, 1 ) =    0
xyphendoy1 ( 241, 2 ) =    0
xyphendoy1 ( 241, 3 ) =    0
xyphendoy1 ( 241, 4 ) =    0
xyphendoy1 ( 241, 5 ) =    0
xyphendoy1 ( 241, 6 ) =    0
xyphendoy1 ( 241, 7 ) =    0
xyphendoy1 ( 241, 8 ) =    0
      
 ! xlat         64.2500000000000     
greenupxy  ( 240, 1 ) =  160
greenupxy  ( 240, 2 ) =  152
greenupxy  ( 240, 3 ) =  167
greenupxy  ( 240, 4 ) =  164
greenupxy  ( 240, 5 ) =  164
greenupxy  ( 240, 6 ) =  164
greenupxy  ( 240, 7 ) =  145
greenupxy  ( 240, 8 ) =  145
fallxy     ( 240, 1 ) =  231
fallxy     ( 240, 2 ) =  230
fallxy     ( 240, 3 ) =  233
fallxy     ( 240, 4 ) =  233
fallxy     ( 240, 5 ) =  233
fallxy     ( 240, 6 ) =  233
fallxy     ( 240, 7 ) =  230
fallxy     ( 240, 8 ) =  230
xyphendoy1 ( 240, 1 ) =    0
xyphendoy1 ( 240, 2 ) =    0
xyphendoy1 ( 240, 3 ) =    0
xyphendoy1 ( 240, 4 ) =    0
xyphendoy1 ( 240, 5 ) =    0
xyphendoy1 ( 240, 6 ) =    0
xyphendoy1 ( 240, 7 ) =    0
xyphendoy1 ( 240, 8 ) =    0
      
 ! xlat         63.7500000000000     
greenupxy  ( 239, 1 ) =  159
greenupxy  ( 239, 2 ) =  151
greenupxy  ( 239, 3 ) =  167
greenupxy  ( 239, 4 ) =  163
greenupxy  ( 239, 5 ) =  163
greenupxy  ( 239, 6 ) =  163
greenupxy  ( 239, 7 ) =  144
greenupxy  ( 239, 8 ) =  144
fallxy     ( 239, 1 ) =  231
fallxy     ( 239, 2 ) =  229
fallxy     ( 239, 3 ) =  233
fallxy     ( 239, 4 ) =  233
fallxy     ( 239, 5 ) =  233
fallxy     ( 239, 6 ) =  233
fallxy     ( 239, 7 ) =  229
fallxy     ( 239, 8 ) =  229
xyphendoy1 ( 239, 1 ) =    0
xyphendoy1 ( 239, 2 ) =    0
xyphendoy1 ( 239, 3 ) =    0
xyphendoy1 ( 239, 4 ) =    0
xyphendoy1 ( 239, 5 ) =    0
xyphendoy1 ( 239, 6 ) =    0
xyphendoy1 ( 239, 7 ) =    0
xyphendoy1 ( 239, 8 ) =    0
      
 ! xlat         63.2500000000000     
greenupxy  ( 238, 1 ) =  158
greenupxy  ( 238, 2 ) =  150
greenupxy  ( 238, 3 ) =  166
greenupxy  ( 238, 4 ) =  163
greenupxy  ( 238, 5 ) =  163
greenupxy  ( 238, 6 ) =  163
greenupxy  ( 238, 7 ) =  144
greenupxy  ( 238, 8 ) =  144
fallxy     ( 238, 1 ) =  232
fallxy     ( 238, 2 ) =  229
fallxy     ( 238, 3 ) =  233
fallxy     ( 238, 4 ) =  233
fallxy     ( 238, 5 ) =  233
fallxy     ( 238, 6 ) =  233
fallxy     ( 238, 7 ) =  227
fallxy     ( 238, 8 ) =  227
xyphendoy1 ( 238, 1 ) =    0
xyphendoy1 ( 238, 2 ) =    0
xyphendoy1 ( 238, 3 ) =    0
xyphendoy1 ( 238, 4 ) =    0
xyphendoy1 ( 238, 5 ) =    0
xyphendoy1 ( 238, 6 ) =    0
xyphendoy1 ( 238, 7 ) =    0
xyphendoy1 ( 238, 8 ) =    0
      
 ! xlat         62.7500000000000     
greenupxy  ( 237, 1 ) =  157
greenupxy  ( 237, 2 ) =  148
greenupxy  ( 237, 3 ) =  166
greenupxy  ( 237, 4 ) =  162
greenupxy  ( 237, 5 ) =  162
greenupxy  ( 237, 6 ) =  162
greenupxy  ( 237, 7 ) =  144
greenupxy  ( 237, 8 ) =  144
fallxy     ( 237, 1 ) =  232
fallxy     ( 237, 2 ) =  229
fallxy     ( 237, 3 ) =  233
fallxy     ( 237, 4 ) =  233
fallxy     ( 237, 5 ) =  233
fallxy     ( 237, 6 ) =  233
fallxy     ( 237, 7 ) =  227
fallxy     ( 237, 8 ) =  227
xyphendoy1 ( 237, 1 ) =    0
xyphendoy1 ( 237, 2 ) =    0
xyphendoy1 ( 237, 3 ) =    0
xyphendoy1 ( 237, 4 ) =    0
xyphendoy1 ( 237, 5 ) =    0
xyphendoy1 ( 237, 6 ) =    0
xyphendoy1 ( 237, 7 ) =    0
xyphendoy1 ( 237, 8 ) =    0
      
 ! xlat         62.2500000000000     
greenupxy  ( 236, 1 ) =  156
greenupxy  ( 236, 2 ) =  147
greenupxy  ( 236, 3 ) =  165
greenupxy  ( 236, 4 ) =  161
greenupxy  ( 236, 5 ) =  161
greenupxy  ( 236, 6 ) =  161
greenupxy  ( 236, 7 ) =  143
greenupxy  ( 236, 8 ) =  143
fallxy     ( 236, 1 ) =  232
fallxy     ( 236, 2 ) =  229
fallxy     ( 236, 3 ) =  234
fallxy     ( 236, 4 ) =  233
fallxy     ( 236, 5 ) =  233
fallxy     ( 236, 6 ) =  233
fallxy     ( 236, 7 ) =  226
fallxy     ( 236, 8 ) =  226
xyphendoy1 ( 236, 1 ) =    0
xyphendoy1 ( 236, 2 ) =    0
xyphendoy1 ( 236, 3 ) =    0
xyphendoy1 ( 236, 4 ) =    0
xyphendoy1 ( 236, 5 ) =    0
xyphendoy1 ( 236, 6 ) =    0
xyphendoy1 ( 236, 7 ) =    0
xyphendoy1 ( 236, 8 ) =    0
      
 ! xlat         61.7500000000000     
greenupxy  ( 235, 1 ) =  155
greenupxy  ( 235, 2 ) =  146
greenupxy  ( 235, 3 ) =  165
greenupxy  ( 235, 4 ) =  161
greenupxy  ( 235, 5 ) =  161
greenupxy  ( 235, 6 ) =  161
greenupxy  ( 235, 7 ) =  142
greenupxy  ( 235, 8 ) =  142
fallxy     ( 235, 1 ) =  232
fallxy     ( 235, 2 ) =  228
fallxy     ( 235, 3 ) =  234
fallxy     ( 235, 4 ) =  233
fallxy     ( 235, 5 ) =  233
fallxy     ( 235, 6 ) =  233
fallxy     ( 235, 7 ) =  224
fallxy     ( 235, 8 ) =  224
xyphendoy1 ( 235, 1 ) =    0
xyphendoy1 ( 235, 2 ) =    0
xyphendoy1 ( 235, 3 ) =    0
xyphendoy1 ( 235, 4 ) =    0
xyphendoy1 ( 235, 5 ) =    0
xyphendoy1 ( 235, 6 ) =    0
xyphendoy1 ( 235, 7 ) =    0
xyphendoy1 ( 235, 8 ) =    0
      
 ! xlat         61.2500000000000     
greenupxy  ( 234, 1 ) =  155
greenupxy  ( 234, 2 ) =  145
greenupxy  ( 234, 3 ) =  165
greenupxy  ( 234, 4 ) =  160
greenupxy  ( 234, 5 ) =  160
greenupxy  ( 234, 6 ) =  160
greenupxy  ( 234, 7 ) =  141
greenupxy  ( 234, 8 ) =  141
fallxy     ( 234, 1 ) =  232
fallxy     ( 234, 2 ) =  228
fallxy     ( 234, 3 ) =  234
fallxy     ( 234, 4 ) =  233
fallxy     ( 234, 5 ) =  233
fallxy     ( 234, 6 ) =  233
fallxy     ( 234, 7 ) =  222
fallxy     ( 234, 8 ) =  222
xyphendoy1 ( 234, 1 ) =    0
xyphendoy1 ( 234, 2 ) =    0
xyphendoy1 ( 234, 3 ) =    0
xyphendoy1 ( 234, 4 ) =    0
xyphendoy1 ( 234, 5 ) =    0
xyphendoy1 ( 234, 6 ) =    0
xyphendoy1 ( 234, 7 ) =    0
xyphendoy1 ( 234, 8 ) =    0
      
 ! xlat         60.7500000000000     
greenupxy  ( 233, 1 ) =  154
greenupxy  ( 233, 2 ) =  144
greenupxy  ( 233, 3 ) =  165
greenupxy  ( 233, 4 ) =  159
greenupxy  ( 233, 5 ) =  159
greenupxy  ( 233, 6 ) =  159
greenupxy  ( 233, 7 ) =  140
greenupxy  ( 233, 8 ) =  140
fallxy     ( 233, 1 ) =  232
fallxy     ( 233, 2 ) =  227
fallxy     ( 233, 3 ) =  234
fallxy     ( 233, 4 ) =  232
fallxy     ( 233, 5 ) =  232
fallxy     ( 233, 6 ) =  232
fallxy     ( 233, 7 ) =  220
fallxy     ( 233, 8 ) =  220
xyphendoy1 ( 233, 1 ) =    0
xyphendoy1 ( 233, 2 ) =    0
xyphendoy1 ( 233, 3 ) =    0
xyphendoy1 ( 233, 4 ) =    0
xyphendoy1 ( 233, 5 ) =    0
xyphendoy1 ( 233, 6 ) =    0
xyphendoy1 ( 233, 7 ) =    0
xyphendoy1 ( 233, 8 ) =    0
      
 ! xlat         60.2500000000000     
greenupxy  ( 232, 1 ) =  153
greenupxy  ( 232, 2 ) =  144
greenupxy  ( 232, 3 ) =  165
greenupxy  ( 232, 4 ) =  159
greenupxy  ( 232, 5 ) =  159
greenupxy  ( 232, 6 ) =  159
greenupxy  ( 232, 7 ) =  139
greenupxy  ( 232, 8 ) =  139
fallxy     ( 232, 1 ) =  232
fallxy     ( 232, 2 ) =  226
fallxy     ( 232, 3 ) =  235
fallxy     ( 232, 4 ) =  232
fallxy     ( 232, 5 ) =  232
fallxy     ( 232, 6 ) =  232
fallxy     ( 232, 7 ) =  220
fallxy     ( 232, 8 ) =  220
xyphendoy1 ( 232, 1 ) =    0
xyphendoy1 ( 232, 2 ) =    0
xyphendoy1 ( 232, 3 ) =    0
xyphendoy1 ( 232, 4 ) =    0
xyphendoy1 ( 232, 5 ) =    0
xyphendoy1 ( 232, 6 ) =    0
xyphendoy1 ( 232, 7 ) =    0
xyphendoy1 ( 232, 8 ) =    0
      
 ! xlat         59.7500000000000     
greenupxy  ( 231, 1 ) =  153
greenupxy  ( 231, 2 ) =  143
greenupxy  ( 231, 3 ) =  165
greenupxy  ( 231, 4 ) =  158
greenupxy  ( 231, 5 ) =  158
greenupxy  ( 231, 6 ) =  158
greenupxy  ( 231, 7 ) =  138
greenupxy  ( 231, 8 ) =  138
fallxy     ( 231, 1 ) =  232
fallxy     ( 231, 2 ) =  226
fallxy     ( 231, 3 ) =  235
fallxy     ( 231, 4 ) =  232
fallxy     ( 231, 5 ) =  232
fallxy     ( 231, 6 ) =  232
fallxy     ( 231, 7 ) =  218
fallxy     ( 231, 8 ) =  218
xyphendoy1 ( 231, 1 ) =    0
xyphendoy1 ( 231, 2 ) =    0
xyphendoy1 ( 231, 3 ) =    0
xyphendoy1 ( 231, 4 ) =    0
xyphendoy1 ( 231, 5 ) =    0
xyphendoy1 ( 231, 6 ) =    0
xyphendoy1 ( 231, 7 ) =    0
xyphendoy1 ( 231, 8 ) =    0
      
 ! xlat         59.2500000000000     
greenupxy  ( 230, 1 ) =  152
greenupxy  ( 230, 2 ) =  143
greenupxy  ( 230, 3 ) =  165
greenupxy  ( 230, 4 ) =  157
greenupxy  ( 230, 5 ) =  157
greenupxy  ( 230, 6 ) =  157
greenupxy  ( 230, 7 ) =  138
greenupxy  ( 230, 8 ) =  138
fallxy     ( 230, 1 ) =  233
fallxy     ( 230, 2 ) =  226
fallxy     ( 230, 3 ) =  235
fallxy     ( 230, 4 ) =  232
fallxy     ( 230, 5 ) =  232
fallxy     ( 230, 6 ) =  232
fallxy     ( 230, 7 ) =  217
fallxy     ( 230, 8 ) =  217
xyphendoy1 ( 230, 1 ) =    0
xyphendoy1 ( 230, 2 ) =    0
xyphendoy1 ( 230, 3 ) =    0
xyphendoy1 ( 230, 4 ) =    0
xyphendoy1 ( 230, 5 ) =    0
xyphendoy1 ( 230, 6 ) =    0
xyphendoy1 ( 230, 7 ) =    0
xyphendoy1 ( 230, 8 ) =    0
      
 ! xlat         58.7500000000000     
greenupxy  ( 229, 1 ) =  151
greenupxy  ( 229, 2 ) =  143
greenupxy  ( 229, 3 ) =  165
greenupxy  ( 229, 4 ) =  156
greenupxy  ( 229, 5 ) =  156
greenupxy  ( 229, 6 ) =  156
greenupxy  ( 229, 7 ) =  138
greenupxy  ( 229, 8 ) =  138
fallxy     ( 229, 1 ) =  233
fallxy     ( 229, 2 ) =  226
fallxy     ( 229, 3 ) =  235
fallxy     ( 229, 4 ) =  232
fallxy     ( 229, 5 ) =  232
fallxy     ( 229, 6 ) =  232
fallxy     ( 229, 7 ) =  217
fallxy     ( 229, 8 ) =  217
xyphendoy1 ( 229, 1 ) =    0
xyphendoy1 ( 229, 2 ) =    0
xyphendoy1 ( 229, 3 ) =    0
xyphendoy1 ( 229, 4 ) =    0
xyphendoy1 ( 229, 5 ) =    0
xyphendoy1 ( 229, 6 ) =    0
xyphendoy1 ( 229, 7 ) =    0
xyphendoy1 ( 229, 8 ) =    0
      
 ! xlat         58.2500000000000     
greenupxy  ( 228, 1 ) =  151
greenupxy  ( 228, 2 ) =  142
greenupxy  ( 228, 3 ) =  164
greenupxy  ( 228, 4 ) =  156
greenupxy  ( 228, 5 ) =  156
greenupxy  ( 228, 6 ) =  156
greenupxy  ( 228, 7 ) =  137
greenupxy  ( 228, 8 ) =  137
fallxy     ( 228, 1 ) =  233
fallxy     ( 228, 2 ) =  225
fallxy     ( 228, 3 ) =  235
fallxy     ( 228, 4 ) =  232
fallxy     ( 228, 5 ) =  232
fallxy     ( 228, 6 ) =  232
fallxy     ( 228, 7 ) =  217
fallxy     ( 228, 8 ) =  217
xyphendoy1 ( 228, 1 ) =    0
xyphendoy1 ( 228, 2 ) =    0
xyphendoy1 ( 228, 3 ) =    0
xyphendoy1 ( 228, 4 ) =    0
xyphendoy1 ( 228, 5 ) =    0
xyphendoy1 ( 228, 6 ) =    0
xyphendoy1 ( 228, 7 ) =    0
xyphendoy1 ( 228, 8 ) =    0
      
 ! xlat         57.7500000000000     
greenupxy  ( 227, 1 ) =  150
greenupxy  ( 227, 2 ) =  141
greenupxy  ( 227, 3 ) =  164
greenupxy  ( 227, 4 ) =  155
greenupxy  ( 227, 5 ) =  155
greenupxy  ( 227, 6 ) =  155
greenupxy  ( 227, 7 ) =  136
greenupxy  ( 227, 8 ) =  136
fallxy     ( 227, 1 ) =  233
fallxy     ( 227, 2 ) =  225
fallxy     ( 227, 3 ) =  236
fallxy     ( 227, 4 ) =  231
fallxy     ( 227, 5 ) =  231
fallxy     ( 227, 6 ) =  231
fallxy     ( 227, 7 ) =  216
fallxy     ( 227, 8 ) =  216
xyphendoy1 ( 227, 1 ) =    0
xyphendoy1 ( 227, 2 ) =    0
xyphendoy1 ( 227, 3 ) =    0
xyphendoy1 ( 227, 4 ) =    0
xyphendoy1 ( 227, 5 ) =    0
xyphendoy1 ( 227, 6 ) =    0
xyphendoy1 ( 227, 7 ) =    0
xyphendoy1 ( 227, 8 ) =    0
      
 ! xlat         57.2500000000000     
greenupxy  ( 226, 1 ) =  149
greenupxy  ( 226, 2 ) =  140
greenupxy  ( 226, 3 ) =  164
greenupxy  ( 226, 4 ) =  154
greenupxy  ( 226, 5 ) =  154
greenupxy  ( 226, 6 ) =  154
greenupxy  ( 226, 7 ) =  135
greenupxy  ( 226, 8 ) =  135
fallxy     ( 226, 1 ) =  233
fallxy     ( 226, 2 ) =  225
fallxy     ( 226, 3 ) =  236
fallxy     ( 226, 4 ) =  231
fallxy     ( 226, 5 ) =  231
fallxy     ( 226, 6 ) =  231
fallxy     ( 226, 7 ) =  216
fallxy     ( 226, 8 ) =  216
xyphendoy1 ( 226, 1 ) =    0
xyphendoy1 ( 226, 2 ) =    0
xyphendoy1 ( 226, 3 ) =    0
xyphendoy1 ( 226, 4 ) =    0
xyphendoy1 ( 226, 5 ) =    0
xyphendoy1 ( 226, 6 ) =    0
xyphendoy1 ( 226, 7 ) =    0
xyphendoy1 ( 226, 8 ) =    0
      
 ! xlat         56.7500000000000     
greenupxy  ( 225, 1 ) =  149
greenupxy  ( 225, 2 ) =  139
greenupxy  ( 225, 3 ) =  164
greenupxy  ( 225, 4 ) =  153
greenupxy  ( 225, 5 ) =  153
greenupxy  ( 225, 6 ) =  153
greenupxy  ( 225, 7 ) =  134
greenupxy  ( 225, 8 ) =  134
fallxy     ( 225, 1 ) =  233
fallxy     ( 225, 2 ) =  225
fallxy     ( 225, 3 ) =  236
fallxy     ( 225, 4 ) =  231
fallxy     ( 225, 5 ) =  231
fallxy     ( 225, 6 ) =  231
fallxy     ( 225, 7 ) =  215
fallxy     ( 225, 8 ) =  215
xyphendoy1 ( 225, 1 ) =    0
xyphendoy1 ( 225, 2 ) =    0
xyphendoy1 ( 225, 3 ) =    0
xyphendoy1 ( 225, 4 ) =    0
xyphendoy1 ( 225, 5 ) =    0
xyphendoy1 ( 225, 6 ) =    0
xyphendoy1 ( 225, 7 ) =    0
xyphendoy1 ( 225, 8 ) =    0
      
 ! xlat         56.2500000000000     
greenupxy  ( 224, 1 ) =  148
greenupxy  ( 224, 2 ) =  137
greenupxy  ( 224, 3 ) =  163
greenupxy  ( 224, 4 ) =  152
greenupxy  ( 224, 5 ) =  152
greenupxy  ( 224, 6 ) =  152
greenupxy  ( 224, 7 ) =  133
greenupxy  ( 224, 8 ) =  133
fallxy     ( 224, 1 ) =  234
fallxy     ( 224, 2 ) =  225
fallxy     ( 224, 3 ) =  236
fallxy     ( 224, 4 ) =  230
fallxy     ( 224, 5 ) =  230
fallxy     ( 224, 6 ) =  230
fallxy     ( 224, 7 ) =  215
fallxy     ( 224, 8 ) =  215
xyphendoy1 ( 224, 1 ) =    0
xyphendoy1 ( 224, 2 ) =    0
xyphendoy1 ( 224, 3 ) =    0
xyphendoy1 ( 224, 4 ) =    0
xyphendoy1 ( 224, 5 ) =    0
xyphendoy1 ( 224, 6 ) =    0
xyphendoy1 ( 224, 7 ) =    0
xyphendoy1 ( 224, 8 ) =    0
      
 ! xlat         55.7500000000000     
greenupxy  ( 223, 1 ) =  147
greenupxy  ( 223, 2 ) =  136
greenupxy  ( 223, 3 ) =  163
greenupxy  ( 223, 4 ) =  150
greenupxy  ( 223, 5 ) =  150
greenupxy  ( 223, 6 ) =  150
greenupxy  ( 223, 7 ) =  132
greenupxy  ( 223, 8 ) =  132
fallxy     ( 223, 1 ) =  234
fallxy     ( 223, 2 ) =  224
fallxy     ( 223, 3 ) =  236
fallxy     ( 223, 4 ) =  229
fallxy     ( 223, 5 ) =  229
fallxy     ( 223, 6 ) =  229
fallxy     ( 223, 7 ) =  214
fallxy     ( 223, 8 ) =  214
xyphendoy1 ( 223, 1 ) =    0
xyphendoy1 ( 223, 2 ) =    0
xyphendoy1 ( 223, 3 ) =    0
xyphendoy1 ( 223, 4 ) =    0
xyphendoy1 ( 223, 5 ) =    0
xyphendoy1 ( 223, 6 ) =    0
xyphendoy1 ( 223, 7 ) =    0
xyphendoy1 ( 223, 8 ) =    0
      
 ! xlat         55.2500000000000     
greenupxy  ( 222, 1 ) =  147
greenupxy  ( 222, 2 ) =  134
greenupxy  ( 222, 3 ) =  162
greenupxy  ( 222, 4 ) =  149
greenupxy  ( 222, 5 ) =  149
greenupxy  ( 222, 6 ) =  149
greenupxy  ( 222, 7 ) =  132
greenupxy  ( 222, 8 ) =  132
fallxy     ( 222, 1 ) =  234
fallxy     ( 222, 2 ) =  224
fallxy     ( 222, 3 ) =  236
fallxy     ( 222, 4 ) =  228
fallxy     ( 222, 5 ) =  228
fallxy     ( 222, 6 ) =  228
fallxy     ( 222, 7 ) =  214
fallxy     ( 222, 8 ) =  214
xyphendoy1 ( 222, 1 ) =    0
xyphendoy1 ( 222, 2 ) =    0
xyphendoy1 ( 222, 3 ) =    0
xyphendoy1 ( 222, 4 ) =    0
xyphendoy1 ( 222, 5 ) =    0
xyphendoy1 ( 222, 6 ) =    0
xyphendoy1 ( 222, 7 ) =    0
xyphendoy1 ( 222, 8 ) =    0
      
 ! xlat         54.7500000000000     
greenupxy  ( 221, 1 ) =  146
greenupxy  ( 221, 2 ) =  133
greenupxy  ( 221, 3 ) =  161
greenupxy  ( 221, 4 ) =  147
greenupxy  ( 221, 5 ) =  147
greenupxy  ( 221, 6 ) =  147
greenupxy  ( 221, 7 ) =  131
greenupxy  ( 221, 8 ) =  131
fallxy     ( 221, 1 ) =  234
fallxy     ( 221, 2 ) =  223
fallxy     ( 221, 3 ) =  235
fallxy     ( 221, 4 ) =  227
fallxy     ( 221, 5 ) =  227
fallxy     ( 221, 6 ) =  227
fallxy     ( 221, 7 ) =  214
fallxy     ( 221, 8 ) =  214
xyphendoy1 ( 221, 1 ) =    0
xyphendoy1 ( 221, 2 ) =    0
xyphendoy1 ( 221, 3 ) =    0
xyphendoy1 ( 221, 4 ) =    0
xyphendoy1 ( 221, 5 ) =    0
xyphendoy1 ( 221, 6 ) =    0
xyphendoy1 ( 221, 7 ) =    0
xyphendoy1 ( 221, 8 ) =    0
      
 ! xlat         54.2500000000000     
greenupxy  ( 220, 1 ) =  146
greenupxy  ( 220, 2 ) =  132
greenupxy  ( 220, 3 ) =  161
greenupxy  ( 220, 4 ) =  145
greenupxy  ( 220, 5 ) =  145
greenupxy  ( 220, 6 ) =  145
greenupxy  ( 220, 7 ) =  130
greenupxy  ( 220, 8 ) =  130
fallxy     ( 220, 1 ) =  234
fallxy     ( 220, 2 ) =  223
fallxy     ( 220, 3 ) =  235
fallxy     ( 220, 4 ) =  225
fallxy     ( 220, 5 ) =  225
fallxy     ( 220, 6 ) =  225
fallxy     ( 220, 7 ) =  214
fallxy     ( 220, 8 ) =  214
xyphendoy1 ( 220, 1 ) =    0
xyphendoy1 ( 220, 2 ) =    0
xyphendoy1 ( 220, 3 ) =    0
xyphendoy1 ( 220, 4 ) =    0
xyphendoy1 ( 220, 5 ) =    0
xyphendoy1 ( 220, 6 ) =    0
xyphendoy1 ( 220, 7 ) =    0
xyphendoy1 ( 220, 8 ) =    0
      
 ! xlat         53.7500000000000     
greenupxy  ( 219, 1 ) =  145
greenupxy  ( 219, 2 ) =  130
greenupxy  ( 219, 3 ) =  160
greenupxy  ( 219, 4 ) =  143
greenupxy  ( 219, 5 ) =  143
greenupxy  ( 219, 6 ) =  143
greenupxy  ( 219, 7 ) =  129
greenupxy  ( 219, 8 ) =  129
fallxy     ( 219, 1 ) =  234
fallxy     ( 219, 2 ) =  222
fallxy     ( 219, 3 ) =  235
fallxy     ( 219, 4 ) =  224
fallxy     ( 219, 5 ) =  224
fallxy     ( 219, 6 ) =  224
fallxy     ( 219, 7 ) =  214
fallxy     ( 219, 8 ) =  214
xyphendoy1 ( 219, 1 ) =    0
xyphendoy1 ( 219, 2 ) =    0
xyphendoy1 ( 219, 3 ) =    0
xyphendoy1 ( 219, 4 ) =    0
xyphendoy1 ( 219, 5 ) =    0
xyphendoy1 ( 219, 6 ) =    0
xyphendoy1 ( 219, 7 ) =    0
xyphendoy1 ( 219, 8 ) =    0
      
 ! xlat         53.2500000000000     
greenupxy  ( 218, 1 ) =  144
greenupxy  ( 218, 2 ) =  129
greenupxy  ( 218, 3 ) =  158
greenupxy  ( 218, 4 ) =  141
greenupxy  ( 218, 5 ) =  141
greenupxy  ( 218, 6 ) =  141
greenupxy  ( 218, 7 ) =  128
greenupxy  ( 218, 8 ) =  128
fallxy     ( 218, 1 ) =  234
fallxy     ( 218, 2 ) =  222
fallxy     ( 218, 3 ) =  234
fallxy     ( 218, 4 ) =  222
fallxy     ( 218, 5 ) =  222
fallxy     ( 218, 6 ) =  222
fallxy     ( 218, 7 ) =  214
fallxy     ( 218, 8 ) =  214
xyphendoy1 ( 218, 1 ) =    0
xyphendoy1 ( 218, 2 ) =    0
xyphendoy1 ( 218, 3 ) =    0
xyphendoy1 ( 218, 4 ) =    0
xyphendoy1 ( 218, 5 ) =    0
xyphendoy1 ( 218, 6 ) =    0
xyphendoy1 ( 218, 7 ) =    0
xyphendoy1 ( 218, 8 ) =    0
      
 ! xlat         52.7500000000000     
greenupxy  ( 217, 1 ) =  144
greenupxy  ( 217, 2 ) =  127
greenupxy  ( 217, 3 ) =  157
greenupxy  ( 217, 4 ) =  139
greenupxy  ( 217, 5 ) =  139
greenupxy  ( 217, 6 ) =  139
greenupxy  ( 217, 7 ) =  127
greenupxy  ( 217, 8 ) =  127
fallxy     ( 217, 1 ) =  234
fallxy     ( 217, 2 ) =  221
fallxy     ( 217, 3 ) =  234
fallxy     ( 217, 4 ) =  221
fallxy     ( 217, 5 ) =  221
fallxy     ( 217, 6 ) =  221
fallxy     ( 217, 7 ) =  214
fallxy     ( 217, 8 ) =  214
xyphendoy1 ( 217, 1 ) =    0
xyphendoy1 ( 217, 2 ) =    0
xyphendoy1 ( 217, 3 ) =    0
xyphendoy1 ( 217, 4 ) =    0
xyphendoy1 ( 217, 5 ) =    0
xyphendoy1 ( 217, 6 ) =    0
xyphendoy1 ( 217, 7 ) =    0
xyphendoy1 ( 217, 8 ) =    0
      
 ! xlat         52.2500000000000     
greenupxy  ( 216, 1 ) =  143
greenupxy  ( 216, 2 ) =  126
greenupxy  ( 216, 3 ) =  155
greenupxy  ( 216, 4 ) =  137
greenupxy  ( 216, 5 ) =  137
greenupxy  ( 216, 6 ) =  137
greenupxy  ( 216, 7 ) =  126
greenupxy  ( 216, 8 ) =  126
fallxy     ( 216, 1 ) =  234
fallxy     ( 216, 2 ) =  221
fallxy     ( 216, 3 ) =  233
fallxy     ( 216, 4 ) =  219
fallxy     ( 216, 5 ) =  219
fallxy     ( 216, 6 ) =  219
fallxy     ( 216, 7 ) =  214
fallxy     ( 216, 8 ) =  214
xyphendoy1 ( 216, 1 ) =    0
xyphendoy1 ( 216, 2 ) =    0
xyphendoy1 ( 216, 3 ) =    0
xyphendoy1 ( 216, 4 ) =    0
xyphendoy1 ( 216, 5 ) =    0
xyphendoy1 ( 216, 6 ) =    0
xyphendoy1 ( 216, 7 ) =    0
xyphendoy1 ( 216, 8 ) =    0
      
 ! xlat         51.7500000000000     
greenupxy  ( 215, 1 ) =  142
greenupxy  ( 215, 2 ) =  124
greenupxy  ( 215, 3 ) =  153
greenupxy  ( 215, 4 ) =  135
greenupxy  ( 215, 5 ) =  135
greenupxy  ( 215, 6 ) =  135
greenupxy  ( 215, 7 ) =  125
greenupxy  ( 215, 8 ) =  125
fallxy     ( 215, 1 ) =  234
fallxy     ( 215, 2 ) =  221
fallxy     ( 215, 3 ) =  233
fallxy     ( 215, 4 ) =  218
fallxy     ( 215, 5 ) =  218
fallxy     ( 215, 6 ) =  218
fallxy     ( 215, 7 ) =  214
fallxy     ( 215, 8 ) =  214
xyphendoy1 ( 215, 1 ) =    0
xyphendoy1 ( 215, 2 ) =    0
xyphendoy1 ( 215, 3 ) =    0
xyphendoy1 ( 215, 4 ) =    0
xyphendoy1 ( 215, 5 ) =    0
xyphendoy1 ( 215, 6 ) =    0
xyphendoy1 ( 215, 7 ) =    0
xyphendoy1 ( 215, 8 ) =    0
      
 ! xlat         51.2500000000000     
greenupxy  ( 214, 1 ) =  142
greenupxy  ( 214, 2 ) =  123
greenupxy  ( 214, 3 ) =  151
greenupxy  ( 214, 4 ) =  133
greenupxy  ( 214, 5 ) =  133
greenupxy  ( 214, 6 ) =  133
greenupxy  ( 214, 7 ) =  124
greenupxy  ( 214, 8 ) =  124
fallxy     ( 214, 1 ) =  234
fallxy     ( 214, 2 ) =  221
fallxy     ( 214, 3 ) =  232
fallxy     ( 214, 4 ) =  217
fallxy     ( 214, 5 ) =  217
fallxy     ( 214, 6 ) =  217
fallxy     ( 214, 7 ) =  215
fallxy     ( 214, 8 ) =  215
xyphendoy1 ( 214, 1 ) =    0
xyphendoy1 ( 214, 2 ) =    0
xyphendoy1 ( 214, 3 ) =    0
xyphendoy1 ( 214, 4 ) =    0
xyphendoy1 ( 214, 5 ) =    0
xyphendoy1 ( 214, 6 ) =    0
xyphendoy1 ( 214, 7 ) =    0
xyphendoy1 ( 214, 8 ) =    0
      
 ! xlat         50.7500000000000     
greenupxy  ( 213, 1 ) =  142
greenupxy  ( 213, 2 ) =  122
greenupxy  ( 213, 3 ) =  149
greenupxy  ( 213, 4 ) =  131
greenupxy  ( 213, 5 ) =  131
greenupxy  ( 213, 6 ) =  131
greenupxy  ( 213, 7 ) =  123
greenupxy  ( 213, 8 ) =  123
fallxy     ( 213, 1 ) =  234
fallxy     ( 213, 2 ) =  221
fallxy     ( 213, 3 ) =  230
fallxy     ( 213, 4 ) =  215
fallxy     ( 213, 5 ) =  215
fallxy     ( 213, 6 ) =  215
fallxy     ( 213, 7 ) =  215
fallxy     ( 213, 8 ) =  215
xyphendoy1 ( 213, 1 ) =    0
xyphendoy1 ( 213, 2 ) =    0
xyphendoy1 ( 213, 3 ) =    0
xyphendoy1 ( 213, 4 ) =    0
xyphendoy1 ( 213, 5 ) =    0
xyphendoy1 ( 213, 6 ) =    0
xyphendoy1 ( 213, 7 ) =    0
xyphendoy1 ( 213, 8 ) =    0
      
 ! xlat         50.2500000000000     
greenupxy  ( 212, 1 ) =  141
greenupxy  ( 212, 2 ) =  121
greenupxy  ( 212, 3 ) =  146
greenupxy  ( 212, 4 ) =  129
greenupxy  ( 212, 5 ) =  129
greenupxy  ( 212, 6 ) =  129
greenupxy  ( 212, 7 ) =  121
greenupxy  ( 212, 8 ) =  121
fallxy     ( 212, 1 ) =  234
fallxy     ( 212, 2 ) =  222
fallxy     ( 212, 3 ) =  228
fallxy     ( 212, 4 ) =  213
fallxy     ( 212, 5 ) =  213
fallxy     ( 212, 6 ) =  213
fallxy     ( 212, 7 ) =  216
fallxy     ( 212, 8 ) =  216
xyphendoy1 ( 212, 1 ) =    0
xyphendoy1 ( 212, 2 ) =    0
xyphendoy1 ( 212, 3 ) =    0
xyphendoy1 ( 212, 4 ) =    0
xyphendoy1 ( 212, 5 ) =    0
xyphendoy1 ( 212, 6 ) =    0
xyphendoy1 ( 212, 7 ) =    0
xyphendoy1 ( 212, 8 ) =    0
      
 ! xlat         49.7500000000000     
greenupxy  ( 211, 1 ) =  141
greenupxy  ( 211, 2 ) =  120
greenupxy  ( 211, 3 ) =  143
greenupxy  ( 211, 4 ) =  127
greenupxy  ( 211, 5 ) =  127
greenupxy  ( 211, 6 ) =  127
greenupxy  ( 211, 7 ) =  120
greenupxy  ( 211, 8 ) =  120
fallxy     ( 211, 1 ) =  233
fallxy     ( 211, 2 ) =  222
fallxy     ( 211, 3 ) =  226
fallxy     ( 211, 4 ) =  212
fallxy     ( 211, 5 ) =  212
fallxy     ( 211, 6 ) =  212
fallxy     ( 211, 7 ) =  217
fallxy     ( 211, 8 ) =  217
xyphendoy1 ( 211, 1 ) =    0
xyphendoy1 ( 211, 2 ) =    0
xyphendoy1 ( 211, 3 ) =    0
xyphendoy1 ( 211, 4 ) =    0
xyphendoy1 ( 211, 5 ) =    0
xyphendoy1 ( 211, 6 ) =    0
xyphendoy1 ( 211, 7 ) =    0
xyphendoy1 ( 211, 8 ) =    0
      
 ! xlat         49.2500000000000     
greenupxy  ( 210, 1 ) =  140
greenupxy  ( 210, 2 ) =  119
greenupxy  ( 210, 3 ) =  141
greenupxy  ( 210, 4 ) =  125
greenupxy  ( 210, 5 ) =  125
greenupxy  ( 210, 6 ) =  125
greenupxy  ( 210, 7 ) =  119
greenupxy  ( 210, 8 ) =  119
fallxy     ( 210, 1 ) =  233
fallxy     ( 210, 2 ) =  222
fallxy     ( 210, 3 ) =  225
fallxy     ( 210, 4 ) =  211
fallxy     ( 210, 5 ) =  211
fallxy     ( 210, 6 ) =  211
fallxy     ( 210, 7 ) =  218
fallxy     ( 210, 8 ) =  218
xyphendoy1 ( 210, 1 ) =    0
xyphendoy1 ( 210, 2 ) =    0
xyphendoy1 ( 210, 3 ) =    0
xyphendoy1 ( 210, 4 ) =    0
xyphendoy1 ( 210, 5 ) =    0
xyphendoy1 ( 210, 6 ) =    0
xyphendoy1 ( 210, 7 ) =    0
xyphendoy1 ( 210, 8 ) =    0
      
 ! xlat         48.7500000000000     
greenupxy  ( 209, 1 ) =  139
greenupxy  ( 209, 2 ) =  119
greenupxy  ( 209, 3 ) =  138
greenupxy  ( 209, 4 ) =  124
greenupxy  ( 209, 5 ) =  124
greenupxy  ( 209, 6 ) =  124
greenupxy  ( 209, 7 ) =  119
greenupxy  ( 209, 8 ) =  119
fallxy     ( 209, 1 ) =  233
fallxy     ( 209, 2 ) =  222
fallxy     ( 209, 3 ) =  222
fallxy     ( 209, 4 ) =  209
fallxy     ( 209, 5 ) =  209
fallxy     ( 209, 6 ) =  209
fallxy     ( 209, 7 ) =  219
fallxy     ( 209, 8 ) =  219
xyphendoy1 ( 209, 1 ) =    0
xyphendoy1 ( 209, 2 ) =    0
xyphendoy1 ( 209, 3 ) =    0
xyphendoy1 ( 209, 4 ) =    0
xyphendoy1 ( 209, 5 ) =    0
xyphendoy1 ( 209, 6 ) =    0
xyphendoy1 ( 209, 7 ) =    0
xyphendoy1 ( 209, 8 ) =    0
      
 ! xlat         48.2500000000000     
greenupxy  ( 208, 1 ) =  139
greenupxy  ( 208, 2 ) =  119
greenupxy  ( 208, 3 ) =  135
greenupxy  ( 208, 4 ) =  122
greenupxy  ( 208, 5 ) =  122
greenupxy  ( 208, 6 ) =  122
greenupxy  ( 208, 7 ) =  118
greenupxy  ( 208, 8 ) =  118
fallxy     ( 208, 1 ) =  232
fallxy     ( 208, 2 ) =  222
fallxy     ( 208, 3 ) =  220
fallxy     ( 208, 4 ) =  208
fallxy     ( 208, 5 ) =  208
fallxy     ( 208, 6 ) =  208
fallxy     ( 208, 7 ) =  220
fallxy     ( 208, 8 ) =  220
xyphendoy1 ( 208, 1 ) =    0
xyphendoy1 ( 208, 2 ) =    0
xyphendoy1 ( 208, 3 ) =    0
xyphendoy1 ( 208, 4 ) =    0
xyphendoy1 ( 208, 5 ) =    0
xyphendoy1 ( 208, 6 ) =    0
xyphendoy1 ( 208, 7 ) =    0
xyphendoy1 ( 208, 8 ) =    0
      
 ! xlat         47.7500000000000     
greenupxy  ( 207, 1 ) =  138
greenupxy  ( 207, 2 ) =  119
greenupxy  ( 207, 3 ) =  132
greenupxy  ( 207, 4 ) =  120
greenupxy  ( 207, 5 ) =  120
greenupxy  ( 207, 6 ) =  120
greenupxy  ( 207, 7 ) =  118
greenupxy  ( 207, 8 ) =  118
fallxy     ( 207, 1 ) =  232
fallxy     ( 207, 2 ) =  223
fallxy     ( 207, 3 ) =  217
fallxy     ( 207, 4 ) =  207
fallxy     ( 207, 5 ) =  207
fallxy     ( 207, 6 ) =  207
fallxy     ( 207, 7 ) =  222
fallxy     ( 207, 8 ) =  222
xyphendoy1 ( 207, 1 ) =    0
xyphendoy1 ( 207, 2 ) =    0
xyphendoy1 ( 207, 3 ) =    0
xyphendoy1 ( 207, 4 ) =    0
xyphendoy1 ( 207, 5 ) =    0
xyphendoy1 ( 207, 6 ) =    0
xyphendoy1 ( 207, 7 ) =    0
xyphendoy1 ( 207, 8 ) =    0
      
 ! xlat         47.2500000000000     
greenupxy  ( 206, 1 ) =  138
greenupxy  ( 206, 2 ) =  119
greenupxy  ( 206, 3 ) =  129
greenupxy  ( 206, 4 ) =  118
greenupxy  ( 206, 5 ) =  118
greenupxy  ( 206, 6 ) =  118
greenupxy  ( 206, 7 ) =  118
greenupxy  ( 206, 8 ) =  118
fallxy     ( 206, 1 ) =  232
fallxy     ( 206, 2 ) =  223
fallxy     ( 206, 3 ) =  215
fallxy     ( 206, 4 ) =  205
fallxy     ( 206, 5 ) =  205
fallxy     ( 206, 6 ) =  205
fallxy     ( 206, 7 ) =  223
fallxy     ( 206, 8 ) =  223
xyphendoy1 ( 206, 1 ) =    0
xyphendoy1 ( 206, 2 ) =    0
xyphendoy1 ( 206, 3 ) =    0
xyphendoy1 ( 206, 4 ) =    0
xyphendoy1 ( 206, 5 ) =    0
xyphendoy1 ( 206, 6 ) =    0
xyphendoy1 ( 206, 7 ) =    0
xyphendoy1 ( 206, 8 ) =    0
      
 ! xlat         46.7500000000000     
greenupxy  ( 205, 1 ) =  138
greenupxy  ( 205, 2 ) =  119
greenupxy  ( 205, 3 ) =  127
greenupxy  ( 205, 4 ) =  118
greenupxy  ( 205, 5 ) =  118
greenupxy  ( 205, 6 ) =  118
greenupxy  ( 205, 7 ) =  118
greenupxy  ( 205, 8 ) =  118
fallxy     ( 205, 1 ) =  231
fallxy     ( 205, 2 ) =  224
fallxy     ( 205, 3 ) =  213
fallxy     ( 205, 4 ) =  205
fallxy     ( 205, 5 ) =  205
fallxy     ( 205, 6 ) =  205
fallxy     ( 205, 7 ) =  224
fallxy     ( 205, 8 ) =  224
xyphendoy1 ( 205, 1 ) =    0
xyphendoy1 ( 205, 2 ) =    0
xyphendoy1 ( 205, 3 ) =    0
xyphendoy1 ( 205, 4 ) =    0
xyphendoy1 ( 205, 5 ) =    0
xyphendoy1 ( 205, 6 ) =    0
xyphendoy1 ( 205, 7 ) =    0
xyphendoy1 ( 205, 8 ) =    0
      
 ! xlat         46.2500000000000     
greenupxy  ( 204, 1 ) =  137
greenupxy  ( 204, 2 ) =  119
greenupxy  ( 204, 3 ) =  124
greenupxy  ( 204, 4 ) =  117
greenupxy  ( 204, 5 ) =  117
greenupxy  ( 204, 6 ) =  117
greenupxy  ( 204, 7 ) =  118
greenupxy  ( 204, 8 ) =  118
fallxy     ( 204, 1 ) =  230
fallxy     ( 204, 2 ) =  224
fallxy     ( 204, 3 ) =  210
fallxy     ( 204, 4 ) =  205
fallxy     ( 204, 5 ) =  205
fallxy     ( 204, 6 ) =  205
fallxy     ( 204, 7 ) =  226
fallxy     ( 204, 8 ) =  226
xyphendoy1 ( 204, 1 ) =    0
xyphendoy1 ( 204, 2 ) =    0
xyphendoy1 ( 204, 3 ) =    0
xyphendoy1 ( 204, 4 ) =    0
xyphendoy1 ( 204, 5 ) =    0
xyphendoy1 ( 204, 6 ) =    0
xyphendoy1 ( 204, 7 ) =    0
xyphendoy1 ( 204, 8 ) =    0
      
 ! xlat         45.7500000000000     
greenupxy  ( 203, 1 ) =  136
greenupxy  ( 203, 2 ) =  119
greenupxy  ( 203, 3 ) =  121
greenupxy  ( 203, 4 ) =  116
greenupxy  ( 203, 5 ) =  116
greenupxy  ( 203, 6 ) =  116
greenupxy  ( 203, 7 ) =  118
greenupxy  ( 203, 8 ) =  118
fallxy     ( 203, 1 ) =  229
fallxy     ( 203, 2 ) =  226
fallxy     ( 203, 3 ) =  206
fallxy     ( 203, 4 ) =  205
fallxy     ( 203, 5 ) =  205
fallxy     ( 203, 6 ) =  205
fallxy     ( 203, 7 ) =  227
fallxy     ( 203, 8 ) =  227
xyphendoy1 ( 203, 1 ) =    0
xyphendoy1 ( 203, 2 ) =    0
xyphendoy1 ( 203, 3 ) =    0
xyphendoy1 ( 203, 4 ) =    0
xyphendoy1 ( 203, 5 ) =    0
xyphendoy1 ( 203, 6 ) =    0
xyphendoy1 ( 203, 7 ) =    0
xyphendoy1 ( 203, 8 ) =    0
      
 ! xlat         45.2500000000000     
greenupxy  ( 202, 1 ) =  136
greenupxy  ( 202, 2 ) =  120
greenupxy  ( 202, 3 ) =  117
greenupxy  ( 202, 4 ) =  116
greenupxy  ( 202, 5 ) =  116
greenupxy  ( 202, 6 ) =  116
greenupxy  ( 202, 7 ) =  117
greenupxy  ( 202, 8 ) =  117
fallxy     ( 202, 1 ) =  229
fallxy     ( 202, 2 ) =  227
fallxy     ( 202, 3 ) =  203
fallxy     ( 202, 4 ) =  205
fallxy     ( 202, 5 ) =  205
fallxy     ( 202, 6 ) =  205
fallxy     ( 202, 7 ) =  229
fallxy     ( 202, 8 ) =  229
xyphendoy1 ( 202, 1 ) =    0
xyphendoy1 ( 202, 2 ) =    0
xyphendoy1 ( 202, 3 ) =    0
xyphendoy1 ( 202, 4 ) =    0
xyphendoy1 ( 202, 5 ) =    0
xyphendoy1 ( 202, 6 ) =    0
xyphendoy1 ( 202, 7 ) =    0
xyphendoy1 ( 202, 8 ) =    0
      
 ! xlat         44.7500000000000     
greenupxy  ( 201, 1 ) =  136
greenupxy  ( 201, 2 ) =  120
greenupxy  ( 201, 3 ) =  115
greenupxy  ( 201, 4 ) =  116
greenupxy  ( 201, 5 ) =  116
greenupxy  ( 201, 6 ) =  116
greenupxy  ( 201, 7 ) =  117
greenupxy  ( 201, 8 ) =  117
fallxy     ( 201, 1 ) =  228
fallxy     ( 201, 2 ) =  227
fallxy     ( 201, 3 ) =  200
fallxy     ( 201, 4 ) =  206
fallxy     ( 201, 5 ) =  206
fallxy     ( 201, 6 ) =  206
fallxy     ( 201, 7 ) =  230
fallxy     ( 201, 8 ) =  230
xyphendoy1 ( 201, 1 ) =    0
xyphendoy1 ( 201, 2 ) =    0
xyphendoy1 ( 201, 3 ) =    0
xyphendoy1 ( 201, 4 ) =    0
xyphendoy1 ( 201, 5 ) =    0
xyphendoy1 ( 201, 6 ) =    0
xyphendoy1 ( 201, 7 ) =    0
xyphendoy1 ( 201, 8 ) =    0
      
 ! xlat         44.2500000000000     
greenupxy  ( 200, 1 ) =  135
greenupxy  ( 200, 2 ) =  119
greenupxy  ( 200, 3 ) =  112
greenupxy  ( 200, 4 ) =  115
greenupxy  ( 200, 5 ) =  115
greenupxy  ( 200, 6 ) =  115
greenupxy  ( 200, 7 ) =  117
greenupxy  ( 200, 8 ) =  117
fallxy     ( 200, 1 ) =  227
fallxy     ( 200, 2 ) =  228
fallxy     ( 200, 3 ) =  197
fallxy     ( 200, 4 ) =  206
fallxy     ( 200, 5 ) =  206
fallxy     ( 200, 6 ) =  206
fallxy     ( 200, 7 ) =  230
fallxy     ( 200, 8 ) =  230
xyphendoy1 ( 200, 1 ) =    0
xyphendoy1 ( 200, 2 ) =    0
xyphendoy1 ( 200, 3 ) =    0
xyphendoy1 ( 200, 4 ) =    0
xyphendoy1 ( 200, 5 ) =    0
xyphendoy1 ( 200, 6 ) =    0
xyphendoy1 ( 200, 7 ) =    0
xyphendoy1 ( 200, 8 ) =    0
      
 ! xlat         43.7500000000000     
greenupxy  ( 199, 1 ) =  135
greenupxy  ( 199, 2 ) =  119
greenupxy  ( 199, 3 ) =  109
greenupxy  ( 199, 4 ) =  115
greenupxy  ( 199, 5 ) =  115
greenupxy  ( 199, 6 ) =  115
greenupxy  ( 199, 7 ) =  116
greenupxy  ( 199, 8 ) =  116
fallxy     ( 199, 1 ) =  226
fallxy     ( 199, 2 ) =  228
fallxy     ( 199, 3 ) =  194
fallxy     ( 199, 4 ) =  206
fallxy     ( 199, 5 ) =  206
fallxy     ( 199, 6 ) =  206
fallxy     ( 199, 7 ) =  231
fallxy     ( 199, 8 ) =  231
xyphendoy1 ( 199, 1 ) =    0
xyphendoy1 ( 199, 2 ) =    0
xyphendoy1 ( 199, 3 ) =    0
xyphendoy1 ( 199, 4 ) =    0
xyphendoy1 ( 199, 5 ) =    0
xyphendoy1 ( 199, 6 ) =    0
xyphendoy1 ( 199, 7 ) =    0
xyphendoy1 ( 199, 8 ) =    0
      
 ! xlat         43.2500000000000     
greenupxy  ( 198, 1 ) =  134
greenupxy  ( 198, 2 ) =  119
greenupxy  ( 198, 3 ) =  107
greenupxy  ( 198, 4 ) =  115
greenupxy  ( 198, 5 ) =  115
greenupxy  ( 198, 6 ) =  115
greenupxy  ( 198, 7 ) =  116
greenupxy  ( 198, 8 ) =  116
fallxy     ( 198, 1 ) =  226
fallxy     ( 198, 2 ) =  228
fallxy     ( 198, 3 ) =  192
fallxy     ( 198, 4 ) =  206
fallxy     ( 198, 5 ) =  206
fallxy     ( 198, 6 ) =  206
fallxy     ( 198, 7 ) =  232
fallxy     ( 198, 8 ) =  232
xyphendoy1 ( 198, 1 ) =    0
xyphendoy1 ( 198, 2 ) =    0
xyphendoy1 ( 198, 3 ) =    0
xyphendoy1 ( 198, 4 ) =    0
xyphendoy1 ( 198, 5 ) =    0
xyphendoy1 ( 198, 6 ) =    0
xyphendoy1 ( 198, 7 ) =    0
xyphendoy1 ( 198, 8 ) =    0
      
 ! xlat         42.7500000000000     
greenupxy  ( 197, 1 ) =  133
greenupxy  ( 197, 2 ) =  118
greenupxy  ( 197, 3 ) =  106
greenupxy  ( 197, 4 ) =  114
greenupxy  ( 197, 5 ) =  114
greenupxy  ( 197, 6 ) =  114
greenupxy  ( 197, 7 ) =  115
greenupxy  ( 197, 8 ) =  115
fallxy     ( 197, 1 ) =  225
fallxy     ( 197, 2 ) =  228
fallxy     ( 197, 3 ) =  191
fallxy     ( 197, 4 ) =  207
fallxy     ( 197, 5 ) =  207
fallxy     ( 197, 6 ) =  207
fallxy     ( 197, 7 ) =  233
fallxy     ( 197, 8 ) =  233
xyphendoy1 ( 197, 1 ) =    0
xyphendoy1 ( 197, 2 ) =    0
xyphendoy1 ( 197, 3 ) =    0
xyphendoy1 ( 197, 4 ) =    0
xyphendoy1 ( 197, 5 ) =    0
xyphendoy1 ( 197, 6 ) =    0
xyphendoy1 ( 197, 7 ) =    0
xyphendoy1 ( 197, 8 ) =    0
      
 ! xlat         42.2500000000000     
greenupxy  ( 196, 1 ) =  133
greenupxy  ( 196, 2 ) =  117
greenupxy  ( 196, 3 ) =  104
greenupxy  ( 196, 4 ) =  114
greenupxy  ( 196, 5 ) =  114
greenupxy  ( 196, 6 ) =  114
greenupxy  ( 196, 7 ) =  114
greenupxy  ( 196, 8 ) =  114
fallxy     ( 196, 1 ) =  225
fallxy     ( 196, 2 ) =  228
fallxy     ( 196, 3 ) =  187
fallxy     ( 196, 4 ) =  207
fallxy     ( 196, 5 ) =  207
fallxy     ( 196, 6 ) =  207
fallxy     ( 196, 7 ) =  233
fallxy     ( 196, 8 ) =  233
xyphendoy1 ( 196, 1 ) =    0
xyphendoy1 ( 196, 2 ) =    0
xyphendoy1 ( 196, 3 ) =    0
xyphendoy1 ( 196, 4 ) =    0
xyphendoy1 ( 196, 5 ) =    0
xyphendoy1 ( 196, 6 ) =    0
xyphendoy1 ( 196, 7 ) =    0
xyphendoy1 ( 196, 8 ) =    0
      
 ! xlat         41.7500000000000     
greenupxy  ( 195, 1 ) =  132
greenupxy  ( 195, 2 ) =  116
greenupxy  ( 195, 3 ) =  102
greenupxy  ( 195, 4 ) =  113
greenupxy  ( 195, 5 ) =  113
greenupxy  ( 195, 6 ) =  113
greenupxy  ( 195, 7 ) =  112
greenupxy  ( 195, 8 ) =  112
fallxy     ( 195, 1 ) =  224
fallxy     ( 195, 2 ) =  227
fallxy     ( 195, 3 ) =  186
fallxy     ( 195, 4 ) =  208
fallxy     ( 195, 5 ) =  208
fallxy     ( 195, 6 ) =  208
fallxy     ( 195, 7 ) =  232
fallxy     ( 195, 8 ) =  232
xyphendoy1 ( 195, 1 ) =    0
xyphendoy1 ( 195, 2 ) =    0
xyphendoy1 ( 195, 3 ) =    0
xyphendoy1 ( 195, 4 ) =    0
xyphendoy1 ( 195, 5 ) =    0
xyphendoy1 ( 195, 6 ) =    0
xyphendoy1 ( 195, 7 ) =    0
xyphendoy1 ( 195, 8 ) =    0
      
 ! xlat         41.2500000000000     
greenupxy  ( 194, 1 ) =  132
greenupxy  ( 194, 2 ) =  115
greenupxy  ( 194, 3 ) =  101
greenupxy  ( 194, 4 ) =  113
greenupxy  ( 194, 5 ) =  113
greenupxy  ( 194, 6 ) =  113
greenupxy  ( 194, 7 ) =  111
greenupxy  ( 194, 8 ) =  111
fallxy     ( 194, 1 ) =  225
fallxy     ( 194, 2 ) =  227
fallxy     ( 194, 3 ) =  184
fallxy     ( 194, 4 ) =  209
fallxy     ( 194, 5 ) =  209
fallxy     ( 194, 6 ) =  209
fallxy     ( 194, 7 ) =  233
fallxy     ( 194, 8 ) =  233
xyphendoy1 ( 194, 1 ) =    0
xyphendoy1 ( 194, 2 ) =    0
xyphendoy1 ( 194, 3 ) =    0
xyphendoy1 ( 194, 4 ) =    0
xyphendoy1 ( 194, 5 ) =    0
xyphendoy1 ( 194, 6 ) =    0
xyphendoy1 ( 194, 7 ) =    0
xyphendoy1 ( 194, 8 ) =    0
      
 ! xlat         40.7500000000000     
greenupxy  ( 193, 1 ) =  131
greenupxy  ( 193, 2 ) =  113
greenupxy  ( 193, 3 ) =  100
greenupxy  ( 193, 4 ) =  112
greenupxy  ( 193, 5 ) =  112
greenupxy  ( 193, 6 ) =  112
greenupxy  ( 193, 7 ) =  109
greenupxy  ( 193, 8 ) =  109
fallxy     ( 193, 1 ) =  224
fallxy     ( 193, 2 ) =  226
fallxy     ( 193, 3 ) =  180
fallxy     ( 193, 4 ) =  210
fallxy     ( 193, 5 ) =  210
fallxy     ( 193, 6 ) =  210
fallxy     ( 193, 7 ) =  234
fallxy     ( 193, 8 ) =  234
xyphendoy1 ( 193, 1 ) =    0
xyphendoy1 ( 193, 2 ) =    0
xyphendoy1 ( 193, 3 ) =    0
xyphendoy1 ( 193, 4 ) =    0
xyphendoy1 ( 193, 5 ) =    0
xyphendoy1 ( 193, 6 ) =    0
xyphendoy1 ( 193, 7 ) =    0
xyphendoy1 ( 193, 8 ) =    0
      
 ! xlat         40.2500000000000     
greenupxy  ( 192, 1 ) =  131
greenupxy  ( 192, 2 ) =  112
greenupxy  ( 192, 3 ) =   98
greenupxy  ( 192, 4 ) =  111
greenupxy  ( 192, 5 ) =  111
greenupxy  ( 192, 6 ) =  111
greenupxy  ( 192, 7 ) =  108
greenupxy  ( 192, 8 ) =  108
fallxy     ( 192, 1 ) =  224
fallxy     ( 192, 2 ) =  226
fallxy     ( 192, 3 ) =  176
fallxy     ( 192, 4 ) =  210
fallxy     ( 192, 5 ) =  210
fallxy     ( 192, 6 ) =  210
fallxy     ( 192, 7 ) =  234
fallxy     ( 192, 8 ) =  234
xyphendoy1 ( 192, 1 ) =    0
xyphendoy1 ( 192, 2 ) =    0
xyphendoy1 ( 192, 3 ) =    0
xyphendoy1 ( 192, 4 ) =    0
xyphendoy1 ( 192, 5 ) =    0
xyphendoy1 ( 192, 6 ) =    0
xyphendoy1 ( 192, 7 ) =    0
xyphendoy1 ( 192, 8 ) =    0
      
 ! xlat         39.7500000000000     
greenupxy  ( 191, 1 ) =  131
greenupxy  ( 191, 2 ) =  111
greenupxy  ( 191, 3 ) =   95
greenupxy  ( 191, 4 ) =  110
greenupxy  ( 191, 5 ) =  110
greenupxy  ( 191, 6 ) =  110
greenupxy  ( 191, 7 ) =  107
greenupxy  ( 191, 8 ) =  107
fallxy     ( 191, 1 ) =  225
fallxy     ( 191, 2 ) =  226
fallxy     ( 191, 3 ) =  172
fallxy     ( 191, 4 ) =  211
fallxy     ( 191, 5 ) =  211
fallxy     ( 191, 6 ) =  211
fallxy     ( 191, 7 ) =  235
fallxy     ( 191, 8 ) =  235
xyphendoy1 ( 191, 1 ) =    0
xyphendoy1 ( 191, 2 ) =    0
xyphendoy1 ( 191, 3 ) =    0
xyphendoy1 ( 191, 4 ) =    0
xyphendoy1 ( 191, 5 ) =    0
xyphendoy1 ( 191, 6 ) =    0
xyphendoy1 ( 191, 7 ) =    0
xyphendoy1 ( 191, 8 ) =    0
      
 ! xlat         39.2500000000000     
greenupxy  ( 190, 1 ) =  131
greenupxy  ( 190, 2 ) =  109
greenupxy  ( 190, 3 ) =   92
greenupxy  ( 190, 4 ) =  110
greenupxy  ( 190, 5 ) =  110
greenupxy  ( 190, 6 ) =  110
greenupxy  ( 190, 7 ) =  105
greenupxy  ( 190, 8 ) =  105
fallxy     ( 190, 1 ) =  226
fallxy     ( 190, 2 ) =  225
fallxy     ( 190, 3 ) =  167
fallxy     ( 190, 4 ) =  213
fallxy     ( 190, 5 ) =  213
fallxy     ( 190, 6 ) =  213
fallxy     ( 190, 7 ) =  235
fallxy     ( 190, 8 ) =  235
xyphendoy1 ( 190, 1 ) =    0
xyphendoy1 ( 190, 2 ) =    0
xyphendoy1 ( 190, 3 ) =    0
xyphendoy1 ( 190, 4 ) =    0
xyphendoy1 ( 190, 5 ) =    0
xyphendoy1 ( 190, 6 ) =    0
xyphendoy1 ( 190, 7 ) =    0
xyphendoy1 ( 190, 8 ) =    0
      
 ! xlat         38.7500000000000     
greenupxy  ( 189, 1 ) =  131
greenupxy  ( 189, 2 ) =  107
greenupxy  ( 189, 3 ) =   89
greenupxy  ( 189, 4 ) =  110
greenupxy  ( 189, 5 ) =  110
greenupxy  ( 189, 6 ) =  110
greenupxy  ( 189, 7 ) =  103
greenupxy  ( 189, 8 ) =  103
fallxy     ( 189, 1 ) =  226
fallxy     ( 189, 2 ) =  225
fallxy     ( 189, 3 ) =  163
fallxy     ( 189, 4 ) =  215
fallxy     ( 189, 5 ) =  215
fallxy     ( 189, 6 ) =  215
fallxy     ( 189, 7 ) =  236
fallxy     ( 189, 8 ) =  236
xyphendoy1 ( 189, 1 ) =    0
xyphendoy1 ( 189, 2 ) =    0
xyphendoy1 ( 189, 3 ) =    0
xyphendoy1 ( 189, 4 ) =    0
xyphendoy1 ( 189, 5 ) =    0
xyphendoy1 ( 189, 6 ) =    0
xyphendoy1 ( 189, 7 ) =    0
xyphendoy1 ( 189, 8 ) =    0
      
 ! xlat         38.2500000000000     
greenupxy  ( 188, 1 ) =  131
greenupxy  ( 188, 2 ) =  106
greenupxy  ( 188, 3 ) =   86
greenupxy  ( 188, 4 ) =  111
greenupxy  ( 188, 5 ) =  111
greenupxy  ( 188, 6 ) =  111
greenupxy  ( 188, 7 ) =  101
greenupxy  ( 188, 8 ) =  101
fallxy     ( 188, 1 ) =  227
fallxy     ( 188, 2 ) =  224
fallxy     ( 188, 3 ) =  166
fallxy     ( 188, 4 ) =  216
fallxy     ( 188, 5 ) =  216
fallxy     ( 188, 6 ) =  216
fallxy     ( 188, 7 ) =  236
fallxy     ( 188, 8 ) =  236
xyphendoy1 ( 188, 1 ) =    0
xyphendoy1 ( 188, 2 ) =    0
xyphendoy1 ( 188, 3 ) =    0
xyphendoy1 ( 188, 4 ) =    0
xyphendoy1 ( 188, 5 ) =    0
xyphendoy1 ( 188, 6 ) =    0
xyphendoy1 ( 188, 7 ) =    0
xyphendoy1 ( 188, 8 ) =    0
      
 ! xlat         37.7500000000000     
greenupxy  ( 187, 1 ) =  130
greenupxy  ( 187, 2 ) =  104
greenupxy  ( 187, 3 ) =   88
greenupxy  ( 187, 4 ) =  111
greenupxy  ( 187, 5 ) =  111
greenupxy  ( 187, 6 ) =  111
greenupxy  ( 187, 7 ) =   99
greenupxy  ( 187, 8 ) =   99
fallxy     ( 187, 1 ) =  227
fallxy     ( 187, 2 ) =  223
fallxy     ( 187, 3 ) =  170
fallxy     ( 187, 4 ) =  217
fallxy     ( 187, 5 ) =  217
fallxy     ( 187, 6 ) =  217
fallxy     ( 187, 7 ) =  235
fallxy     ( 187, 8 ) =  235
xyphendoy1 ( 187, 1 ) =    0
xyphendoy1 ( 187, 2 ) =    0
xyphendoy1 ( 187, 3 ) =    0
xyphendoy1 ( 187, 4 ) =    0
xyphendoy1 ( 187, 5 ) =    0
xyphendoy1 ( 187, 6 ) =    0
xyphendoy1 ( 187, 7 ) =    0
xyphendoy1 ( 187, 8 ) =    0
      
 ! xlat         37.2500000000000     
greenupxy  ( 186, 1 ) =  130
greenupxy  ( 186, 2 ) =  102
greenupxy  ( 186, 3 ) =   89
greenupxy  ( 186, 4 ) =  111
greenupxy  ( 186, 5 ) =  111
greenupxy  ( 186, 6 ) =  111
greenupxy  ( 186, 7 ) =   96
greenupxy  ( 186, 8 ) =   96
fallxy     ( 186, 1 ) =  229
fallxy     ( 186, 2 ) =  222
fallxy     ( 186, 3 ) =  175
fallxy     ( 186, 4 ) =  219
fallxy     ( 186, 5 ) =  219
fallxy     ( 186, 6 ) =  219
fallxy     ( 186, 7 ) =  235
fallxy     ( 186, 8 ) =  235
xyphendoy1 ( 186, 1 ) =    0
xyphendoy1 ( 186, 2 ) =    0
xyphendoy1 ( 186, 3 ) =    0
xyphendoy1 ( 186, 4 ) =    0
xyphendoy1 ( 186, 5 ) =    0
xyphendoy1 ( 186, 6 ) =    0
xyphendoy1 ( 186, 7 ) =    0
xyphendoy1 ( 186, 8 ) =    0
      
 ! xlat         36.7500000000000     
greenupxy  ( 185, 1 ) =  130
greenupxy  ( 185, 2 ) =  100
greenupxy  ( 185, 3 ) =   93
greenupxy  ( 185, 4 ) =  111
greenupxy  ( 185, 5 ) =  111
greenupxy  ( 185, 6 ) =  111
greenupxy  ( 185, 7 ) =   94
greenupxy  ( 185, 8 ) =   94
fallxy     ( 185, 1 ) =  230
fallxy     ( 185, 2 ) =  221
fallxy     ( 185, 3 ) =  180
fallxy     ( 185, 4 ) =  220
fallxy     ( 185, 5 ) =  220
fallxy     ( 185, 6 ) =  220
fallxy     ( 185, 7 ) =  235
fallxy     ( 185, 8 ) =  235
xyphendoy1 ( 185, 1 ) =    0
xyphendoy1 ( 185, 2 ) =    0
xyphendoy1 ( 185, 3 ) =    0
xyphendoy1 ( 185, 4 ) =    0
xyphendoy1 ( 185, 5 ) =    0
xyphendoy1 ( 185, 6 ) =    0
xyphendoy1 ( 185, 7 ) =    0
xyphendoy1 ( 185, 8 ) =    0
      
 ! xlat         36.2500000000000     
greenupxy  ( 184, 1 ) =  129
greenupxy  ( 184, 2 ) =   98
greenupxy  ( 184, 3 ) =   97
greenupxy  ( 184, 4 ) =  110
greenupxy  ( 184, 5 ) =  110
greenupxy  ( 184, 6 ) =  110
greenupxy  ( 184, 7 ) =   92
greenupxy  ( 184, 8 ) =   92
fallxy     ( 184, 1 ) =  231
fallxy     ( 184, 2 ) =  221
fallxy     ( 184, 3 ) =  184
fallxy     ( 184, 4 ) =  222
fallxy     ( 184, 5 ) =  222
fallxy     ( 184, 6 ) =  222
fallxy     ( 184, 7 ) =  235
fallxy     ( 184, 8 ) =  235
xyphendoy1 ( 184, 1 ) =    0
xyphendoy1 ( 184, 2 ) =    0
xyphendoy1 ( 184, 3 ) =    0
xyphendoy1 ( 184, 4 ) =    0
xyphendoy1 ( 184, 5 ) =    0
xyphendoy1 ( 184, 6 ) =    0
xyphendoy1 ( 184, 7 ) =    0
xyphendoy1 ( 184, 8 ) =    0
      
 ! xlat         35.7500000000000     
greenupxy  ( 183, 1 ) =  128
greenupxy  ( 183, 2 ) =   96
greenupxy  ( 183, 3 ) =  100
greenupxy  ( 183, 4 ) =  110
greenupxy  ( 183, 5 ) =  110
greenupxy  ( 183, 6 ) =  110
greenupxy  ( 183, 7 ) =   91
greenupxy  ( 183, 8 ) =   91
fallxy     ( 183, 1 ) =  233
fallxy     ( 183, 2 ) =  220
fallxy     ( 183, 3 ) =  188
fallxy     ( 183, 4 ) =  224
fallxy     ( 183, 5 ) =  224
fallxy     ( 183, 6 ) =  224
fallxy     ( 183, 7 ) =  235
fallxy     ( 183, 8 ) =  235
xyphendoy1 ( 183, 1 ) =    0
xyphendoy1 ( 183, 2 ) =    0
xyphendoy1 ( 183, 3 ) =    0
xyphendoy1 ( 183, 4 ) =    0
xyphendoy1 ( 183, 5 ) =    0
xyphendoy1 ( 183, 6 ) =    0
xyphendoy1 ( 183, 7 ) =    0
xyphendoy1 ( 183, 8 ) =    0
      
 ! xlat         35.2500000000000     
greenupxy  ( 182, 1 ) =  127
greenupxy  ( 182, 2 ) =   94
greenupxy  ( 182, 3 ) =  104
greenupxy  ( 182, 4 ) =  110
greenupxy  ( 182, 5 ) =  110
greenupxy  ( 182, 6 ) =  110
greenupxy  ( 182, 7 ) =   90
greenupxy  ( 182, 8 ) =   90
fallxy     ( 182, 1 ) =  234
fallxy     ( 182, 2 ) =  220
fallxy     ( 182, 3 ) =  192
fallxy     ( 182, 4 ) =  226
fallxy     ( 182, 5 ) =  226
fallxy     ( 182, 6 ) =  226
fallxy     ( 182, 7 ) =  234
fallxy     ( 182, 8 ) =  234
xyphendoy1 ( 182, 1 ) =    0
xyphendoy1 ( 182, 2 ) =    0
xyphendoy1 ( 182, 3 ) =    0
xyphendoy1 ( 182, 4 ) =    0
xyphendoy1 ( 182, 5 ) =    0
xyphendoy1 ( 182, 6 ) =    0
xyphendoy1 ( 182, 7 ) =    0
xyphendoy1 ( 182, 8 ) =    0
      
 ! xlat         34.7500000000000     
greenupxy  ( 181, 1 ) =  126
greenupxy  ( 181, 2 ) =   93
greenupxy  ( 181, 3 ) =  107
greenupxy  ( 181, 4 ) =  110
greenupxy  ( 181, 5 ) =  110
greenupxy  ( 181, 6 ) =  110
greenupxy  ( 181, 7 ) =   88
greenupxy  ( 181, 8 ) =   88
fallxy     ( 181, 1 ) =  235
fallxy     ( 181, 2 ) =  220
fallxy     ( 181, 3 ) =  197
fallxy     ( 181, 4 ) =  228
fallxy     ( 181, 5 ) =  228
fallxy     ( 181, 6 ) =  228
fallxy     ( 181, 7 ) =  233
fallxy     ( 181, 8 ) =  233
xyphendoy1 ( 181, 1 ) =    0
xyphendoy1 ( 181, 2 ) =    0
xyphendoy1 ( 181, 3 ) =    0
xyphendoy1 ( 181, 4 ) =    0
xyphendoy1 ( 181, 5 ) =    0
xyphendoy1 ( 181, 6 ) =    0
xyphendoy1 ( 181, 7 ) =    0
xyphendoy1 ( 181, 8 ) =    0
      
 ! xlat         34.2500000000000     
greenupxy  ( 180, 1 ) =  125
greenupxy  ( 180, 2 ) =   91
greenupxy  ( 180, 3 ) =  109
greenupxy  ( 180, 4 ) =  111
greenupxy  ( 180, 5 ) =  111
greenupxy  ( 180, 6 ) =  111
greenupxy  ( 180, 7 ) =   87
greenupxy  ( 180, 8 ) =   87
fallxy     ( 180, 1 ) =  236
fallxy     ( 180, 2 ) =  221
fallxy     ( 180, 3 ) =  201
fallxy     ( 180, 4 ) =  231
fallxy     ( 180, 5 ) =  231
fallxy     ( 180, 6 ) =  231
fallxy     ( 180, 7 ) =  233
fallxy     ( 180, 8 ) =  233
xyphendoy1 ( 180, 1 ) =    0
xyphendoy1 ( 180, 2 ) =    0
xyphendoy1 ( 180, 3 ) =    0
xyphendoy1 ( 180, 4 ) =    0
xyphendoy1 ( 180, 5 ) =    0
xyphendoy1 ( 180, 6 ) =    0
xyphendoy1 ( 180, 7 ) =    0
xyphendoy1 ( 180, 8 ) =    0
      
 ! xlat         33.7500000000000     
greenupxy  ( 179, 1 ) =  123
greenupxy  ( 179, 2 ) =   90
greenupxy  ( 179, 3 ) =  113
greenupxy  ( 179, 4 ) =  114
greenupxy  ( 179, 5 ) =  114
greenupxy  ( 179, 6 ) =  114
greenupxy  ( 179, 7 ) =   86
greenupxy  ( 179, 8 ) =   86
fallxy     ( 179, 1 ) =  237
fallxy     ( 179, 2 ) =  222
fallxy     ( 179, 3 ) =  208
fallxy     ( 179, 4 ) =  233
fallxy     ( 179, 5 ) =  233
fallxy     ( 179, 6 ) =  233
fallxy     ( 179, 7 ) =  233
fallxy     ( 179, 8 ) =  233
xyphendoy1 ( 179, 1 ) =    0
xyphendoy1 ( 179, 2 ) =    0
xyphendoy1 ( 179, 3 ) =    0
xyphendoy1 ( 179, 4 ) =    0
xyphendoy1 ( 179, 5 ) =    0
xyphendoy1 ( 179, 6 ) =    0
xyphendoy1 ( 179, 7 ) =    0
xyphendoy1 ( 179, 8 ) =    0
      
 ! xlat         33.2500000000000     
greenupxy  ( 178, 1 ) =  120
greenupxy  ( 178, 2 ) =   89
greenupxy  ( 178, 3 ) =  118
greenupxy  ( 178, 4 ) =  115
greenupxy  ( 178, 5 ) =  115
greenupxy  ( 178, 6 ) =  115
greenupxy  ( 178, 7 ) =   85
greenupxy  ( 178, 8 ) =   85
fallxy     ( 178, 1 ) =  238
fallxy     ( 178, 2 ) =  223
fallxy     ( 178, 3 ) =  213
fallxy     ( 178, 4 ) =  235
fallxy     ( 178, 5 ) =  235
fallxy     ( 178, 6 ) =  235
fallxy     ( 178, 7 ) =  236
fallxy     ( 178, 8 ) =  236
xyphendoy1 ( 178, 1 ) =    0
xyphendoy1 ( 178, 2 ) =    0
xyphendoy1 ( 178, 3 ) =    0
xyphendoy1 ( 178, 4 ) =    0
xyphendoy1 ( 178, 5 ) =    0
xyphendoy1 ( 178, 6 ) =    0
xyphendoy1 ( 178, 7 ) =    0
xyphendoy1 ( 178, 8 ) =    0
      
 ! xlat         32.7500000000000     
greenupxy  ( 177, 1 ) =  118
greenupxy  ( 177, 2 ) =   88
greenupxy  ( 177, 3 ) =  122
greenupxy  ( 177, 4 ) =  117
greenupxy  ( 177, 5 ) =  117
greenupxy  ( 177, 6 ) =  117
greenupxy  ( 177, 7 ) =   84
greenupxy  ( 177, 8 ) =   84
fallxy     ( 177, 1 ) =  238
fallxy     ( 177, 2 ) =  225
fallxy     ( 177, 3 ) =  218
fallxy     ( 177, 4 ) =  238
fallxy     ( 177, 5 ) =  238
fallxy     ( 177, 6 ) =  238
fallxy     ( 177, 7 ) =  235
fallxy     ( 177, 8 ) =  235
xyphendoy1 ( 177, 1 ) =    0
xyphendoy1 ( 177, 2 ) =    0
xyphendoy1 ( 177, 3 ) =    0
xyphendoy1 ( 177, 4 ) =    0
xyphendoy1 ( 177, 5 ) =    0
xyphendoy1 ( 177, 6 ) =    0
xyphendoy1 ( 177, 7 ) =    0
xyphendoy1 ( 177, 8 ) =    0
      
 ! xlat         32.2500000000000     
greenupxy  ( 176, 1 ) =  114
greenupxy  ( 176, 2 ) =   88
greenupxy  ( 176, 3 ) =  126
greenupxy  ( 176, 4 ) =  120
greenupxy  ( 176, 5 ) =  120
greenupxy  ( 176, 6 ) =  120
greenupxy  ( 176, 7 ) =   85
greenupxy  ( 176, 8 ) =   85
fallxy     ( 176, 1 ) =  235
fallxy     ( 176, 2 ) =  226
fallxy     ( 176, 3 ) =  227
fallxy     ( 176, 4 ) =  241
fallxy     ( 176, 5 ) =  241
fallxy     ( 176, 6 ) =  241
fallxy     ( 176, 7 ) =  236
fallxy     ( 176, 8 ) =  236
xyphendoy1 ( 176, 1 ) =    0
xyphendoy1 ( 176, 2 ) =    0
xyphendoy1 ( 176, 3 ) =    0
xyphendoy1 ( 176, 4 ) =    0
xyphendoy1 ( 176, 5 ) =    0
xyphendoy1 ( 176, 6 ) =    0
xyphendoy1 ( 176, 7 ) =    0
xyphendoy1 ( 176, 8 ) =    0
      
 ! xlat         31.7500000000000     
greenupxy  ( 175, 1 ) =  113
greenupxy  ( 175, 2 ) =   88
greenupxy  ( 175, 3 ) =  133
greenupxy  ( 175, 4 ) =  122
greenupxy  ( 175, 5 ) =  122
greenupxy  ( 175, 6 ) =  122
greenupxy  ( 175, 7 ) =   85
greenupxy  ( 175, 8 ) =   85
fallxy     ( 175, 1 ) =  236
fallxy     ( 175, 2 ) =  228
fallxy     ( 175, 3 ) =  236
fallxy     ( 175, 4 ) =  244
fallxy     ( 175, 5 ) =  244
fallxy     ( 175, 6 ) =  244
fallxy     ( 175, 7 ) =  236
fallxy     ( 175, 8 ) =  236
xyphendoy1 ( 175, 1 ) =    0
xyphendoy1 ( 175, 2 ) =    0
xyphendoy1 ( 175, 3 ) =    0
xyphendoy1 ( 175, 4 ) =    0
xyphendoy1 ( 175, 5 ) =    0
xyphendoy1 ( 175, 6 ) =    0
xyphendoy1 ( 175, 7 ) =    0
xyphendoy1 ( 175, 8 ) =    0
      
 ! xlat         31.2500000000000     
greenupxy  ( 174, 1 ) =  110
greenupxy  ( 174, 2 ) =   89
greenupxy  ( 174, 3 ) =  140
greenupxy  ( 174, 4 ) =  124
greenupxy  ( 174, 5 ) =  124
greenupxy  ( 174, 6 ) =  124
greenupxy  ( 174, 7 ) =   85
greenupxy  ( 174, 8 ) =   85
fallxy     ( 174, 1 ) =  236
fallxy     ( 174, 2 ) =  230
fallxy     ( 174, 3 ) =  245
fallxy     ( 174, 4 ) =  247
fallxy     ( 174, 5 ) =  247
fallxy     ( 174, 6 ) =  247
fallxy     ( 174, 7 ) =  236
fallxy     ( 174, 8 ) =  236
xyphendoy1 ( 174, 1 ) =    0
xyphendoy1 ( 174, 2 ) =    0
xyphendoy1 ( 174, 3 ) =    0
xyphendoy1 ( 174, 4 ) =    0
xyphendoy1 ( 174, 5 ) =    0
xyphendoy1 ( 174, 6 ) =    0
xyphendoy1 ( 174, 7 ) =    0
xyphendoy1 ( 174, 8 ) =    0
      
 ! xlat         30.7500000000000     
greenupxy  ( 173, 1 ) =  108
greenupxy  ( 173, 2 ) =   90
greenupxy  ( 173, 3 ) =  148
greenupxy  ( 173, 4 ) =  126
greenupxy  ( 173, 5 ) =  126
greenupxy  ( 173, 6 ) =  126
greenupxy  ( 173, 7 ) =   86
greenupxy  ( 173, 8 ) =   86
fallxy     ( 173, 1 ) =  236
fallxy     ( 173, 2 ) =  232
fallxy     ( 173, 3 ) =  255
fallxy     ( 173, 4 ) =  251
fallxy     ( 173, 5 ) =  251
fallxy     ( 173, 6 ) =  251
fallxy     ( 173, 7 ) =  237
fallxy     ( 173, 8 ) =  237
xyphendoy1 ( 173, 1 ) =    0
xyphendoy1 ( 173, 2 ) =    0
xyphendoy1 ( 173, 3 ) =    0
xyphendoy1 ( 173, 4 ) =    0
xyphendoy1 ( 173, 5 ) =    0
xyphendoy1 ( 173, 6 ) =    0
xyphendoy1 ( 173, 7 ) =    0
xyphendoy1 ( 173, 8 ) =    0
      
 ! xlat         30.2500000000000     
greenupxy  ( 172, 1 ) =  107
greenupxy  ( 172, 2 ) =   93
greenupxy  ( 172, 3 ) =  156
greenupxy  ( 172, 4 ) =  128
greenupxy  ( 172, 5 ) =  128
greenupxy  ( 172, 6 ) =  128
greenupxy  ( 172, 7 ) =   88
greenupxy  ( 172, 8 ) =   88
fallxy     ( 172, 1 ) =  237
fallxy     ( 172, 2 ) =  234
fallxy     ( 172, 3 ) =  265
fallxy     ( 172, 4 ) =  253
fallxy     ( 172, 5 ) =  253
fallxy     ( 172, 6 ) =  253
fallxy     ( 172, 7 ) =  237
fallxy     ( 172, 8 ) =  237
xyphendoy1 ( 172, 1 ) =    0
xyphendoy1 ( 172, 2 ) =    0
xyphendoy1 ( 172, 3 ) =    0
xyphendoy1 ( 172, 4 ) =    0
xyphendoy1 ( 172, 5 ) =    0
xyphendoy1 ( 172, 6 ) =    0
xyphendoy1 ( 172, 7 ) =    0
xyphendoy1 ( 172, 8 ) =    0
      
 ! xlat         29.7500000000000     
greenupxy  ( 171, 1 ) =  108
greenupxy  ( 171, 2 ) =   95
greenupxy  ( 171, 3 ) =  163
greenupxy  ( 171, 4 ) =  129
greenupxy  ( 171, 5 ) =  129
greenupxy  ( 171, 6 ) =  129
greenupxy  ( 171, 7 ) =   90
greenupxy  ( 171, 8 ) =   90
fallxy     ( 171, 1 ) =  238
fallxy     ( 171, 2 ) =  236
fallxy     ( 171, 3 ) =  268
fallxy     ( 171, 4 ) =  256
fallxy     ( 171, 5 ) =  256
fallxy     ( 171, 6 ) =  256
fallxy     ( 171, 7 ) =  238
fallxy     ( 171, 8 ) =  238
xyphendoy1 ( 171, 1 ) =    0
xyphendoy1 ( 171, 2 ) =    0
xyphendoy1 ( 171, 3 ) =    0
xyphendoy1 ( 171, 4 ) =    0
xyphendoy1 ( 171, 5 ) =    0
xyphendoy1 ( 171, 6 ) =    0
xyphendoy1 ( 171, 7 ) =    0
xyphendoy1 ( 171, 8 ) =    0
      
 ! xlat         29.2500000000000     
greenupxy  ( 170, 1 ) =  110
greenupxy  ( 170, 2 ) =  100
greenupxy  ( 170, 3 ) =  167
greenupxy  ( 170, 4 ) =  132
greenupxy  ( 170, 5 ) =  132
greenupxy  ( 170, 6 ) =  132
greenupxy  ( 170, 7 ) =   93
greenupxy  ( 170, 8 ) =   93
fallxy     ( 170, 1 ) =  238
fallxy     ( 170, 2 ) =  238
fallxy     ( 170, 3 ) =  270
fallxy     ( 170, 4 ) =  259
fallxy     ( 170, 5 ) =  259
fallxy     ( 170, 6 ) =  259
fallxy     ( 170, 7 ) =  241
fallxy     ( 170, 8 ) =  241
xyphendoy1 ( 170, 1 ) =    0
xyphendoy1 ( 170, 2 ) =    0
xyphendoy1 ( 170, 3 ) =    0
xyphendoy1 ( 170, 4 ) =    0
xyphendoy1 ( 170, 5 ) =    0
xyphendoy1 ( 170, 6 ) =    0
xyphendoy1 ( 170, 7 ) =    0
xyphendoy1 ( 170, 8 ) =    0
      
 ! xlat         28.7500000000000     
greenupxy  ( 169, 1 ) =  109
greenupxy  ( 169, 2 ) =  105
greenupxy  ( 169, 3 ) =  170
greenupxy  ( 169, 4 ) =  134
greenupxy  ( 169, 5 ) =  134
greenupxy  ( 169, 6 ) =  134
greenupxy  ( 169, 7 ) =   97
greenupxy  ( 169, 8 ) =   97
fallxy     ( 169, 1 ) =  238
fallxy     ( 169, 2 ) =  241
fallxy     ( 169, 3 ) =  272
fallxy     ( 169, 4 ) =  262
fallxy     ( 169, 5 ) =  262
fallxy     ( 169, 6 ) =  262
fallxy     ( 169, 7 ) =  244
fallxy     ( 169, 8 ) =  244
xyphendoy1 ( 169, 1 ) =    0
xyphendoy1 ( 169, 2 ) =    0
xyphendoy1 ( 169, 3 ) =    0
xyphendoy1 ( 169, 4 ) =    0
xyphendoy1 ( 169, 5 ) =    0
xyphendoy1 ( 169, 6 ) =    0
xyphendoy1 ( 169, 7 ) =    0
xyphendoy1 ( 169, 8 ) =    0
      
 ! xlat         28.2500000000000     
greenupxy  ( 168, 1 ) =  112
greenupxy  ( 168, 2 ) =  109
greenupxy  ( 168, 3 ) =  170
greenupxy  ( 168, 4 ) =  136
greenupxy  ( 168, 5 ) =  136
greenupxy  ( 168, 6 ) =  136
greenupxy  ( 168, 7 ) =   99
greenupxy  ( 168, 8 ) =   99
fallxy     ( 168, 1 ) =  239
fallxy     ( 168, 2 ) =  243
fallxy     ( 168, 3 ) =  273
fallxy     ( 168, 4 ) =  264
fallxy     ( 168, 5 ) =  264
fallxy     ( 168, 6 ) =  264
fallxy     ( 168, 7 ) =  247
fallxy     ( 168, 8 ) =  247
xyphendoy1 ( 168, 1 ) =    0
xyphendoy1 ( 168, 2 ) =    0
xyphendoy1 ( 168, 3 ) =    0
xyphendoy1 ( 168, 4 ) =    0
xyphendoy1 ( 168, 5 ) =    0
xyphendoy1 ( 168, 6 ) =    0
xyphendoy1 ( 168, 7 ) =    0
xyphendoy1 ( 168, 8 ) =    0
      
 ! xlat         27.7500000000000     
greenupxy  ( 167, 1 ) =  111
greenupxy  ( 167, 2 ) =  114
greenupxy  ( 167, 3 ) =  170
greenupxy  ( 167, 4 ) =  138
greenupxy  ( 167, 5 ) =  138
greenupxy  ( 167, 6 ) =  138
greenupxy  ( 167, 7 ) =  101
greenupxy  ( 167, 8 ) =  101
fallxy     ( 167, 1 ) =  239
fallxy     ( 167, 2 ) =  246
fallxy     ( 167, 3 ) =  275
fallxy     ( 167, 4 ) =  266
fallxy     ( 167, 5 ) =  266
fallxy     ( 167, 6 ) =  266
fallxy     ( 167, 7 ) =  251
fallxy     ( 167, 8 ) =  251
xyphendoy1 ( 167, 1 ) =    0
xyphendoy1 ( 167, 2 ) =    0
xyphendoy1 ( 167, 3 ) =    0
xyphendoy1 ( 167, 4 ) =    0
xyphendoy1 ( 167, 5 ) =    0
xyphendoy1 ( 167, 6 ) =    0
xyphendoy1 ( 167, 7 ) =    0
xyphendoy1 ( 167, 8 ) =    0
      
 ! xlat         27.2500000000000     
greenupxy  ( 166, 1 ) =  116
greenupxy  ( 166, 2 ) =  119
greenupxy  ( 166, 3 ) =  169
greenupxy  ( 166, 4 ) =  140
greenupxy  ( 166, 5 ) =  140
greenupxy  ( 166, 6 ) =  140
greenupxy  ( 166, 7 ) =  104
greenupxy  ( 166, 8 ) =  104
fallxy     ( 166, 1 ) =  242
fallxy     ( 166, 2 ) =  249
fallxy     ( 166, 3 ) =  276
fallxy     ( 166, 4 ) =  268
fallxy     ( 166, 5 ) =  268
fallxy     ( 166, 6 ) =  268
fallxy     ( 166, 7 ) =  254
fallxy     ( 166, 8 ) =  254
xyphendoy1 ( 166, 1 ) =    0
xyphendoy1 ( 166, 2 ) =    0
xyphendoy1 ( 166, 3 ) =    0
xyphendoy1 ( 166, 4 ) =    0
xyphendoy1 ( 166, 5 ) =    0
xyphendoy1 ( 166, 6 ) =    0
xyphendoy1 ( 166, 7 ) =    0
xyphendoy1 ( 166, 8 ) =    0
      
 ! xlat         26.7500000000000     
greenupxy  ( 165, 1 ) =  115
greenupxy  ( 165, 2 ) =  123
greenupxy  ( 165, 3 ) =  169
greenupxy  ( 165, 4 ) =  141
greenupxy  ( 165, 5 ) =  141
greenupxy  ( 165, 6 ) =  141
greenupxy  ( 165, 7 ) =  106
greenupxy  ( 165, 8 ) =  106
fallxy     ( 165, 1 ) =  242
fallxy     ( 165, 2 ) =  251
fallxy     ( 165, 3 ) =  278
fallxy     ( 165, 4 ) =  270
fallxy     ( 165, 5 ) =  270
fallxy     ( 165, 6 ) =  270
fallxy     ( 165, 7 ) =  260
fallxy     ( 165, 8 ) =  260
xyphendoy1 ( 165, 1 ) =    0
xyphendoy1 ( 165, 2 ) =    0
xyphendoy1 ( 165, 3 ) =    0
xyphendoy1 ( 165, 4 ) =    0
xyphendoy1 ( 165, 5 ) =    0
xyphendoy1 ( 165, 6 ) =    0
xyphendoy1 ( 165, 7 ) =    0
xyphendoy1 ( 165, 8 ) =    0
      
 ! xlat         26.2500000000000     
greenupxy  ( 164, 1 ) =  115
greenupxy  ( 164, 2 ) =  127
greenupxy  ( 164, 3 ) =  168
greenupxy  ( 164, 4 ) =  144
greenupxy  ( 164, 5 ) =  144
greenupxy  ( 164, 6 ) =  144
greenupxy  ( 164, 7 ) =  109
greenupxy  ( 164, 8 ) =  109
fallxy     ( 164, 1 ) =  245
fallxy     ( 164, 2 ) =  254
fallxy     ( 164, 3 ) =  279
fallxy     ( 164, 4 ) =  272
fallxy     ( 164, 5 ) =  272
fallxy     ( 164, 6 ) =  272
fallxy     ( 164, 7 ) =  262
fallxy     ( 164, 8 ) =  262
xyphendoy1 ( 164, 1 ) =    0
xyphendoy1 ( 164, 2 ) =    0
xyphendoy1 ( 164, 3 ) =    0
xyphendoy1 ( 164, 4 ) =    0
xyphendoy1 ( 164, 5 ) =    0
xyphendoy1 ( 164, 6 ) =    0
xyphendoy1 ( 164, 7 ) =    0
xyphendoy1 ( 164, 8 ) =    0
      
 ! xlat         25.7500000000000     
greenupxy  ( 163, 1 ) =  115
greenupxy  ( 163, 2 ) =  131
greenupxy  ( 163, 3 ) =  169
greenupxy  ( 163, 4 ) =  145
greenupxy  ( 163, 5 ) =  145
greenupxy  ( 163, 6 ) =  145
greenupxy  ( 163, 7 ) =  112
greenupxy  ( 163, 8 ) =  112
fallxy     ( 163, 1 ) =  245
fallxy     ( 163, 2 ) =  256
fallxy     ( 163, 3 ) =  279
fallxy     ( 163, 4 ) =  273
fallxy     ( 163, 5 ) =  273
fallxy     ( 163, 6 ) =  273
fallxy     ( 163, 7 ) =  263
fallxy     ( 163, 8 ) =  263
xyphendoy1 ( 163, 1 ) =    0
xyphendoy1 ( 163, 2 ) =    0
xyphendoy1 ( 163, 3 ) =    0
xyphendoy1 ( 163, 4 ) =    0
xyphendoy1 ( 163, 5 ) =    0
xyphendoy1 ( 163, 6 ) =    0
xyphendoy1 ( 163, 7 ) =    0
xyphendoy1 ( 163, 8 ) =    0
      
 ! xlat         25.2500000000000     
greenupxy  ( 162, 1 ) =  116
greenupxy  ( 162, 2 ) =  134
greenupxy  ( 162, 3 ) =  168
greenupxy  ( 162, 4 ) =  145
greenupxy  ( 162, 5 ) =  145
greenupxy  ( 162, 6 ) =  145
greenupxy  ( 162, 7 ) =  115
greenupxy  ( 162, 8 ) =  115
fallxy     ( 162, 1 ) =  245
fallxy     ( 162, 2 ) =  257
fallxy     ( 162, 3 ) =  279
fallxy     ( 162, 4 ) =  273
fallxy     ( 162, 5 ) =  273
fallxy     ( 162, 6 ) =  273
fallxy     ( 162, 7 ) =  265
fallxy     ( 162, 8 ) =  265
xyphendoy1 ( 162, 1 ) =    0
xyphendoy1 ( 162, 2 ) =    0
xyphendoy1 ( 162, 3 ) =    0
xyphendoy1 ( 162, 4 ) =    0
xyphendoy1 ( 162, 5 ) =    0
xyphendoy1 ( 162, 6 ) =    0
xyphendoy1 ( 162, 7 ) =    0
xyphendoy1 ( 162, 8 ) =    0
      
 ! xlat         24.7500000000000     
greenupxy  ( 161, 1 ) =  118
greenupxy  ( 161, 2 ) =  137
greenupxy  ( 161, 3 ) =  167
greenupxy  ( 161, 4 ) =  146
greenupxy  ( 161, 5 ) =  146
greenupxy  ( 161, 6 ) =  146
greenupxy  ( 161, 7 ) =  118
greenupxy  ( 161, 8 ) =  118
fallxy     ( 161, 1 ) =  240
fallxy     ( 161, 2 ) =  258
fallxy     ( 161, 3 ) =  279
fallxy     ( 161, 4 ) =  274
fallxy     ( 161, 5 ) =  274
fallxy     ( 161, 6 ) =  274
fallxy     ( 161, 7 ) =  267
fallxy     ( 161, 8 ) =  267
xyphendoy1 ( 161, 1 ) =    0
xyphendoy1 ( 161, 2 ) =    0
xyphendoy1 ( 161, 3 ) =    0
xyphendoy1 ( 161, 4 ) =    0
xyphendoy1 ( 161, 5 ) =    0
xyphendoy1 ( 161, 6 ) =    0
xyphendoy1 ( 161, 7 ) =    0
xyphendoy1 ( 161, 8 ) =    0
      
 ! xlat         24.2500000000000     
greenupxy  ( 160, 1 ) =  121
greenupxy  ( 160, 2 ) =  140
greenupxy  ( 160, 3 ) =  165
greenupxy  ( 160, 4 ) =  146
greenupxy  ( 160, 5 ) =  146
greenupxy  ( 160, 6 ) =  146
greenupxy  ( 160, 7 ) =  120
greenupxy  ( 160, 8 ) =  120
fallxy     ( 160, 1 ) =  241
fallxy     ( 160, 2 ) =  259
fallxy     ( 160, 3 ) =  279
fallxy     ( 160, 4 ) =  274
fallxy     ( 160, 5 ) =  274
fallxy     ( 160, 6 ) =  274
fallxy     ( 160, 7 ) =  269
fallxy     ( 160, 8 ) =  269
xyphendoy1 ( 160, 1 ) =    0
xyphendoy1 ( 160, 2 ) =    0
xyphendoy1 ( 160, 3 ) =    0
xyphendoy1 ( 160, 4 ) =    0
xyphendoy1 ( 160, 5 ) =    0
xyphendoy1 ( 160, 6 ) =    0
xyphendoy1 ( 160, 7 ) =    0
xyphendoy1 ( 160, 8 ) =    0
      
 ! xlat         23.7500000000000     
greenupxy  ( 159, 1 ) =  127
greenupxy  ( 159, 2 ) =  142
greenupxy  ( 159, 3 ) =  164
greenupxy  ( 159, 4 ) =  146
greenupxy  ( 159, 5 ) =  146
greenupxy  ( 159, 6 ) =  146
greenupxy  ( 159, 7 ) =  121
greenupxy  ( 159, 8 ) =  121
fallxy     ( 159, 1 ) =  249
fallxy     ( 159, 2 ) =  260
fallxy     ( 159, 3 ) =  279
fallxy     ( 159, 4 ) =  274
fallxy     ( 159, 5 ) =  274
fallxy     ( 159, 6 ) =  274
fallxy     ( 159, 7 ) =  271
fallxy     ( 159, 8 ) =  271
xyphendoy1 ( 159, 1 ) =    0
xyphendoy1 ( 159, 2 ) =    0
xyphendoy1 ( 159, 3 ) =    0
xyphendoy1 ( 159, 4 ) =    0
xyphendoy1 ( 159, 5 ) =    0
xyphendoy1 ( 159, 6 ) =    0
xyphendoy1 ( 159, 7 ) =    0
xyphendoy1 ( 159, 8 ) =    0
      
 ! xlat         23.2500000000000     
greenupxy  ( 158, 1 ) =  131
greenupxy  ( 158, 2 ) =  145
greenupxy  ( 158, 3 ) =  162
greenupxy  ( 158, 4 ) =  145
greenupxy  ( 158, 5 ) =  145
greenupxy  ( 158, 6 ) =  145
greenupxy  ( 158, 7 ) =  122
greenupxy  ( 158, 8 ) =  122
fallxy     ( 158, 1 ) =  251
fallxy     ( 158, 2 ) =  262
fallxy     ( 158, 3 ) =  278
fallxy     ( 158, 4 ) =  274
fallxy     ( 158, 5 ) =  274
fallxy     ( 158, 6 ) =  274
fallxy     ( 158, 7 ) =  273
fallxy     ( 158, 8 ) =  273
xyphendoy1 ( 158, 1 ) =    0
xyphendoy1 ( 158, 2 ) =    0
xyphendoy1 ( 158, 3 ) =    0
xyphendoy1 ( 158, 4 ) =    0
xyphendoy1 ( 158, 5 ) =    0
xyphendoy1 ( 158, 6 ) =    0
xyphendoy1 ( 158, 7 ) =    0
xyphendoy1 ( 158, 8 ) =    0
      
 ! xlat         22.7500000000000     
greenupxy  ( 157, 1 ) =  134
greenupxy  ( 157, 2 ) =  147
greenupxy  ( 157, 3 ) =  159
greenupxy  ( 157, 4 ) =  144
greenupxy  ( 157, 5 ) =  144
greenupxy  ( 157, 6 ) =  144
greenupxy  ( 157, 7 ) =  121
greenupxy  ( 157, 8 ) =  121
fallxy     ( 157, 1 ) =  254
fallxy     ( 157, 2 ) =  262
fallxy     ( 157, 3 ) =  278
fallxy     ( 157, 4 ) =  274
fallxy     ( 157, 5 ) =  274
fallxy     ( 157, 6 ) =  274
fallxy     ( 157, 7 ) =  276
fallxy     ( 157, 8 ) =  276
xyphendoy1 ( 157, 1 ) =    0
xyphendoy1 ( 157, 2 ) =    0
xyphendoy1 ( 157, 3 ) =    0
xyphendoy1 ( 157, 4 ) =    0
xyphendoy1 ( 157, 5 ) =    0
xyphendoy1 ( 157, 6 ) =    0
xyphendoy1 ( 157, 7 ) =    0
xyphendoy1 ( 157, 8 ) =    0
      
 ! xlat         22.2500000000000     
greenupxy  ( 156, 1 ) =  137
greenupxy  ( 156, 2 ) =  148
greenupxy  ( 156, 3 ) =  157
greenupxy  ( 156, 4 ) =  144
greenupxy  ( 156, 5 ) =  144
greenupxy  ( 156, 6 ) =  144
greenupxy  ( 156, 7 ) =  122
greenupxy  ( 156, 8 ) =  122
fallxy     ( 156, 1 ) =  256
fallxy     ( 156, 2 ) =  263
fallxy     ( 156, 3 ) =  278
fallxy     ( 156, 4 ) =  274
fallxy     ( 156, 5 ) =  274
fallxy     ( 156, 6 ) =  274
fallxy     ( 156, 7 ) =  279
fallxy     ( 156, 8 ) =  279
xyphendoy1 ( 156, 1 ) =    0
xyphendoy1 ( 156, 2 ) =    0
xyphendoy1 ( 156, 3 ) =    0
xyphendoy1 ( 156, 4 ) =    0
xyphendoy1 ( 156, 5 ) =    0
xyphendoy1 ( 156, 6 ) =    0
xyphendoy1 ( 156, 7 ) =    0
xyphendoy1 ( 156, 8 ) =    0
      
 ! xlat         21.7500000000000     
greenupxy  ( 155, 1 ) =  140
greenupxy  ( 155, 2 ) =  148
greenupxy  ( 155, 3 ) =  159
greenupxy  ( 155, 4 ) =  143
greenupxy  ( 155, 5 ) =  143
greenupxy  ( 155, 6 ) =  143
greenupxy  ( 155, 7 ) =  123
greenupxy  ( 155, 8 ) =  123
fallxy     ( 155, 1 ) =  256
fallxy     ( 155, 2 ) =  264
fallxy     ( 155, 3 ) =  277
fallxy     ( 155, 4 ) =  274
fallxy     ( 155, 5 ) =  274
fallxy     ( 155, 6 ) =  274
fallxy     ( 155, 7 ) =  281
fallxy     ( 155, 8 ) =  281
xyphendoy1 ( 155, 1 ) =    0
xyphendoy1 ( 155, 2 ) =    0
xyphendoy1 ( 155, 3 ) =    0
xyphendoy1 ( 155, 4 ) =    0
xyphendoy1 ( 155, 5 ) =    0
xyphendoy1 ( 155, 6 ) =    0
xyphendoy1 ( 155, 7 ) =    0
xyphendoy1 ( 155, 8 ) =    0
      
 ! xlat         21.2500000000000     
greenupxy  ( 154, 1 ) =  141
greenupxy  ( 154, 2 ) =  148
greenupxy  ( 154, 3 ) =  160
greenupxy  ( 154, 4 ) =  142
greenupxy  ( 154, 5 ) =  142
greenupxy  ( 154, 6 ) =  142
greenupxy  ( 154, 7 ) =  124
greenupxy  ( 154, 8 ) =  124
fallxy     ( 154, 1 ) =  257
fallxy     ( 154, 2 ) =  265
fallxy     ( 154, 3 ) =  277
fallxy     ( 154, 4 ) =  274
fallxy     ( 154, 5 ) =  274
fallxy     ( 154, 6 ) =  274
fallxy     ( 154, 7 ) =  283
fallxy     ( 154, 8 ) =  283
xyphendoy1 ( 154, 1 ) =    0
xyphendoy1 ( 154, 2 ) =    0
xyphendoy1 ( 154, 3 ) =    0
xyphendoy1 ( 154, 4 ) =    0
xyphendoy1 ( 154, 5 ) =    0
xyphendoy1 ( 154, 6 ) =    0
xyphendoy1 ( 154, 7 ) =    0
xyphendoy1 ( 154, 8 ) =    0
      
 ! xlat         20.7500000000000     
greenupxy  ( 153, 1 ) =  141
greenupxy  ( 153, 2 ) =  146
greenupxy  ( 153, 3 ) =  161
greenupxy  ( 153, 4 ) =  146
greenupxy  ( 153, 5 ) =  146
greenupxy  ( 153, 6 ) =  146
greenupxy  ( 153, 7 ) =  124
greenupxy  ( 153, 8 ) =  124
fallxy     ( 153, 1 ) =  261
fallxy     ( 153, 2 ) =  265
fallxy     ( 153, 3 ) =  276
fallxy     ( 153, 4 ) =  274
fallxy     ( 153, 5 ) =  274
fallxy     ( 153, 6 ) =  274
fallxy     ( 153, 7 ) =  284
fallxy     ( 153, 8 ) =  284
xyphendoy1 ( 153, 1 ) =    0
xyphendoy1 ( 153, 2 ) =    0
xyphendoy1 ( 153, 3 ) =    0
xyphendoy1 ( 153, 4 ) =    0
xyphendoy1 ( 153, 5 ) =    0
xyphendoy1 ( 153, 6 ) =    0
xyphendoy1 ( 153, 7 ) =    0
xyphendoy1 ( 153, 8 ) =    0
      
 ! xlat         20.2500000000000     
greenupxy  ( 152, 1 ) =  141
greenupxy  ( 152, 2 ) =  144
greenupxy  ( 152, 3 ) =  163
greenupxy  ( 152, 4 ) =  149
greenupxy  ( 152, 5 ) =  149
greenupxy  ( 152, 6 ) =  149
greenupxy  ( 152, 7 ) =  125
greenupxy  ( 152, 8 ) =  125
fallxy     ( 152, 1 ) =  262
fallxy     ( 152, 2 ) =  266
fallxy     ( 152, 3 ) =  274
fallxy     ( 152, 4 ) =  273
fallxy     ( 152, 5 ) =  273
fallxy     ( 152, 6 ) =  273
fallxy     ( 152, 7 ) =  285
fallxy     ( 152, 8 ) =  285
xyphendoy1 ( 152, 1 ) =    0
xyphendoy1 ( 152, 2 ) =    0
xyphendoy1 ( 152, 3 ) =    0
xyphendoy1 ( 152, 4 ) =    0
xyphendoy1 ( 152, 5 ) =    0
xyphendoy1 ( 152, 6 ) =    0
xyphendoy1 ( 152, 7 ) =    0
xyphendoy1 ( 152, 8 ) =    0
      
 ! xlat         19.7500000000000     
greenupxy  ( 151, 1 ) =  130
greenupxy  ( 151, 2 ) =  142
greenupxy  ( 151, 3 ) =  166
greenupxy  ( 151, 4 ) =  152
greenupxy  ( 151, 5 ) =  152
greenupxy  ( 151, 6 ) =  152
greenupxy  ( 151, 7 ) =  124
greenupxy  ( 151, 8 ) =  124
fallxy     ( 151, 1 ) =  263
fallxy     ( 151, 2 ) =  267
fallxy     ( 151, 3 ) =  273
fallxy     ( 151, 4 ) =  272
fallxy     ( 151, 5 ) =  272
fallxy     ( 151, 6 ) =  272
fallxy     ( 151, 7 ) =  285
fallxy     ( 151, 8 ) =  285
xyphendoy1 ( 151, 1 ) =    0
xyphendoy1 ( 151, 2 ) =    0
xyphendoy1 ( 151, 3 ) =    0
xyphendoy1 ( 151, 4 ) =    0
xyphendoy1 ( 151, 5 ) =    0
xyphendoy1 ( 151, 6 ) =    0
xyphendoy1 ( 151, 7 ) =    0
xyphendoy1 ( 151, 8 ) =    0
      
 ! xlat         19.2500000000000     
greenupxy  ( 150, 1 ) =  130
greenupxy  ( 150, 2 ) =  140
greenupxy  ( 150, 3 ) =  169
greenupxy  ( 150, 4 ) =  155
greenupxy  ( 150, 5 ) =  155
greenupxy  ( 150, 6 ) =  155
greenupxy  ( 150, 7 ) =  125
greenupxy  ( 150, 8 ) =  125
fallxy     ( 150, 1 ) =  263
fallxy     ( 150, 2 ) =  268
fallxy     ( 150, 3 ) =  273
fallxy     ( 150, 4 ) =  271
fallxy     ( 150, 5 ) =  271
fallxy     ( 150, 6 ) =  271
fallxy     ( 150, 7 ) =  285
fallxy     ( 150, 8 ) =  285
xyphendoy1 ( 150, 1 ) =    0
xyphendoy1 ( 150, 2 ) =    0
xyphendoy1 ( 150, 3 ) =    0
xyphendoy1 ( 150, 4 ) =    0
xyphendoy1 ( 150, 5 ) =    0
xyphendoy1 ( 150, 6 ) =    0
xyphendoy1 ( 150, 7 ) =    0
xyphendoy1 ( 150, 8 ) =    0
      
 ! xlat         18.7500000000000     
greenupxy  ( 149, 1 ) =  123
greenupxy  ( 149, 2 ) =  138
greenupxy  ( 149, 3 ) =  172
greenupxy  ( 149, 4 ) =  158
greenupxy  ( 149, 5 ) =  158
greenupxy  ( 149, 6 ) =  158
greenupxy  ( 149, 7 ) =  128
greenupxy  ( 149, 8 ) =  128
fallxy     ( 149, 1 ) =  263
fallxy     ( 149, 2 ) =  269
fallxy     ( 149, 3 ) =  272
fallxy     ( 149, 4 ) =  270
fallxy     ( 149, 5 ) =  270
fallxy     ( 149, 6 ) =  270
fallxy     ( 149, 7 ) =  285
fallxy     ( 149, 8 ) =  285
xyphendoy1 ( 149, 1 ) =    0
xyphendoy1 ( 149, 2 ) =    0
xyphendoy1 ( 149, 3 ) =    0
xyphendoy1 ( 149, 4 ) =    0
xyphendoy1 ( 149, 5 ) =    0
xyphendoy1 ( 149, 6 ) =    0
xyphendoy1 ( 149, 7 ) =    0
xyphendoy1 ( 149, 8 ) =    0
      
 ! xlat         18.2500000000000     
greenupxy  ( 148, 1 ) =  123
greenupxy  ( 148, 2 ) =  137
greenupxy  ( 148, 3 ) =  175
greenupxy  ( 148, 4 ) =  160
greenupxy  ( 148, 5 ) =  160
greenupxy  ( 148, 6 ) =  160
greenupxy  ( 148, 7 ) =  131
greenupxy  ( 148, 8 ) =  131
fallxy     ( 148, 1 ) =  263
fallxy     ( 148, 2 ) =  269
fallxy     ( 148, 3 ) =  272
fallxy     ( 148, 4 ) =  270
fallxy     ( 148, 5 ) =  270
fallxy     ( 148, 6 ) =  270
fallxy     ( 148, 7 ) =  284
fallxy     ( 148, 8 ) =  284
xyphendoy1 ( 148, 1 ) =    0
xyphendoy1 ( 148, 2 ) =    0
xyphendoy1 ( 148, 3 ) =    0
xyphendoy1 ( 148, 4 ) =    0
xyphendoy1 ( 148, 5 ) =    0
xyphendoy1 ( 148, 6 ) =    0
xyphendoy1 ( 148, 7 ) =    0
xyphendoy1 ( 148, 8 ) =    0
      
 ! xlat         17.7500000000000     
greenupxy  ( 147, 1 ) =  126
greenupxy  ( 147, 2 ) =  138
greenupxy  ( 147, 3 ) =  178
greenupxy  ( 147, 4 ) =  163
greenupxy  ( 147, 5 ) =  163
greenupxy  ( 147, 6 ) =  163
greenupxy  ( 147, 7 ) =  134
greenupxy  ( 147, 8 ) =  134
fallxy     ( 147, 1 ) =  262
fallxy     ( 147, 2 ) =  270
fallxy     ( 147, 3 ) =  271
fallxy     ( 147, 4 ) =  269
fallxy     ( 147, 5 ) =  269
fallxy     ( 147, 6 ) =  269
fallxy     ( 147, 7 ) =  284
fallxy     ( 147, 8 ) =  284
xyphendoy1 ( 147, 1 ) =    0
xyphendoy1 ( 147, 2 ) =    0
xyphendoy1 ( 147, 3 ) =    0
xyphendoy1 ( 147, 4 ) =    0
xyphendoy1 ( 147, 5 ) =    0
xyphendoy1 ( 147, 6 ) =    0
xyphendoy1 ( 147, 7 ) =    0
xyphendoy1 ( 147, 8 ) =    0
      
 ! xlat         17.2500000000000     
greenupxy  ( 146, 1 ) =  126
greenupxy  ( 146, 2 ) =  138
greenupxy  ( 146, 3 ) =  180
greenupxy  ( 146, 4 ) =  165
greenupxy  ( 146, 5 ) =  165
greenupxy  ( 146, 6 ) =  165
greenupxy  ( 146, 7 ) =  135
greenupxy  ( 146, 8 ) =  135
fallxy     ( 146, 1 ) =  265
fallxy     ( 146, 2 ) =  271
fallxy     ( 146, 3 ) =  271
fallxy     ( 146, 4 ) =  269
fallxy     ( 146, 5 ) =  269
fallxy     ( 146, 6 ) =  269
fallxy     ( 146, 7 ) =  283
fallxy     ( 146, 8 ) =  283
xyphendoy1 ( 146, 1 ) =    0
xyphendoy1 ( 146, 2 ) =    0
xyphendoy1 ( 146, 3 ) =    0
xyphendoy1 ( 146, 4 ) =    0
xyphendoy1 ( 146, 5 ) =    0
xyphendoy1 ( 146, 6 ) =    0
xyphendoy1 ( 146, 7 ) =    0
xyphendoy1 ( 146, 8 ) =    0
      
 ! xlat         16.7500000000000     
greenupxy  ( 145, 1 ) =  131
greenupxy  ( 145, 2 ) =  139
greenupxy  ( 145, 3 ) =  183
greenupxy  ( 145, 4 ) =  167
greenupxy  ( 145, 5 ) =  167
greenupxy  ( 145, 6 ) =  167
greenupxy  ( 145, 7 ) =  138
greenupxy  ( 145, 8 ) =  138
fallxy     ( 145, 1 ) =  266
fallxy     ( 145, 2 ) =  271
fallxy     ( 145, 3 ) =  271
fallxy     ( 145, 4 ) =  269
fallxy     ( 145, 5 ) =  269
fallxy     ( 145, 6 ) =  269
fallxy     ( 145, 7 ) =  283
fallxy     ( 145, 8 ) =  283
xyphendoy1 ( 145, 1 ) =    0
xyphendoy1 ( 145, 2 ) =    0
xyphendoy1 ( 145, 3 ) =    0
xyphendoy1 ( 145, 4 ) =    0
xyphendoy1 ( 145, 5 ) =    0
xyphendoy1 ( 145, 6 ) =    0
xyphendoy1 ( 145, 7 ) =    0
xyphendoy1 ( 145, 8 ) =    0
      
 ! xlat         16.2500000000000     
greenupxy  ( 144, 1 ) =  135
greenupxy  ( 144, 2 ) =  140
greenupxy  ( 144, 3 ) =  185
greenupxy  ( 144, 4 ) =  168
greenupxy  ( 144, 5 ) =  168
greenupxy  ( 144, 6 ) =  168
greenupxy  ( 144, 7 ) =  141
greenupxy  ( 144, 8 ) =  141
fallxy     ( 144, 1 ) =  275
fallxy     ( 144, 2 ) =  271
fallxy     ( 144, 3 ) =  271
fallxy     ( 144, 4 ) =  269
fallxy     ( 144, 5 ) =  269
fallxy     ( 144, 6 ) =  269
fallxy     ( 144, 7 ) =  282
fallxy     ( 144, 8 ) =  282
xyphendoy1 ( 144, 1 ) =    0
xyphendoy1 ( 144, 2 ) =    0
xyphendoy1 ( 144, 3 ) =    0
xyphendoy1 ( 144, 4 ) =    0
xyphendoy1 ( 144, 5 ) =    0
xyphendoy1 ( 144, 6 ) =    0
xyphendoy1 ( 144, 7 ) =    0
xyphendoy1 ( 144, 8 ) =    0
      
 ! xlat         15.7500000000000     
greenupxy  ( 143, 1 ) =  137
greenupxy  ( 143, 2 ) =  139
greenupxy  ( 143, 3 ) =  186
greenupxy  ( 143, 4 ) =  169
greenupxy  ( 143, 5 ) =  169
greenupxy  ( 143, 6 ) =  169
greenupxy  ( 143, 7 ) =  144
greenupxy  ( 143, 8 ) =  144
fallxy     ( 143, 1 ) =  277
fallxy     ( 143, 2 ) =  272
fallxy     ( 143, 3 ) =  271
fallxy     ( 143, 4 ) =  269
fallxy     ( 143, 5 ) =  269
fallxy     ( 143, 6 ) =  269
fallxy     ( 143, 7 ) =  281
fallxy     ( 143, 8 ) =  281
xyphendoy1 ( 143, 1 ) =    0
xyphendoy1 ( 143, 2 ) =    0
xyphendoy1 ( 143, 3 ) =    0
xyphendoy1 ( 143, 4 ) =    0
xyphendoy1 ( 143, 5 ) =    0
xyphendoy1 ( 143, 6 ) =    0
xyphendoy1 ( 143, 7 ) =    0
xyphendoy1 ( 143, 8 ) =    0
      
 ! xlat         15.2500000000000     
greenupxy  ( 142, 1 ) =  137
greenupxy  ( 142, 2 ) =  139
greenupxy  ( 142, 3 ) =  187
greenupxy  ( 142, 4 ) =  170
greenupxy  ( 142, 5 ) =  170
greenupxy  ( 142, 6 ) =  170
greenupxy  ( 142, 7 ) =  147
greenupxy  ( 142, 8 ) =  147
fallxy     ( 142, 1 ) =  280
fallxy     ( 142, 2 ) =  272
fallxy     ( 142, 3 ) =  271
fallxy     ( 142, 4 ) =  270
fallxy     ( 142, 5 ) =  270
fallxy     ( 142, 6 ) =  270
fallxy     ( 142, 7 ) =  281
fallxy     ( 142, 8 ) =  281
xyphendoy1 ( 142, 1 ) =    0
xyphendoy1 ( 142, 2 ) =    0
xyphendoy1 ( 142, 3 ) =    0
xyphendoy1 ( 142, 4 ) =    0
xyphendoy1 ( 142, 5 ) =    0
xyphendoy1 ( 142, 6 ) =    0
xyphendoy1 ( 142, 7 ) =    0
xyphendoy1 ( 142, 8 ) =    0
      
 ! xlat         14.7500000000000     
greenupxy  ( 141, 1 ) =  137
greenupxy  ( 141, 2 ) =  139
greenupxy  ( 141, 3 ) =  188
greenupxy  ( 141, 4 ) =  172
greenupxy  ( 141, 5 ) =  172
greenupxy  ( 141, 6 ) =  172
greenupxy  ( 141, 7 ) =  150
greenupxy  ( 141, 8 ) =  150
fallxy     ( 141, 1 ) =  282
fallxy     ( 141, 2 ) =  272
fallxy     ( 141, 3 ) =  271
fallxy     ( 141, 4 ) =  270
fallxy     ( 141, 5 ) =  270
fallxy     ( 141, 6 ) =  270
fallxy     ( 141, 7 ) =  280
fallxy     ( 141, 8 ) =  280
xyphendoy1 ( 141, 1 ) =    0
xyphendoy1 ( 141, 2 ) =    0
xyphendoy1 ( 141, 3 ) =    0
xyphendoy1 ( 141, 4 ) =    0
xyphendoy1 ( 141, 5 ) =    0
xyphendoy1 ( 141, 6 ) =    0
xyphendoy1 ( 141, 7 ) =    0
xyphendoy1 ( 141, 8 ) =    0
      
 ! xlat         14.2500000000000     
greenupxy  ( 140, 1 ) =  135
greenupxy  ( 140, 2 ) =  138
greenupxy  ( 140, 3 ) =  189
greenupxy  ( 140, 4 ) =  173
greenupxy  ( 140, 5 ) =  173
greenupxy  ( 140, 6 ) =  173
greenupxy  ( 140, 7 ) =  154
greenupxy  ( 140, 8 ) =  154
fallxy     ( 140, 1 ) =  283
fallxy     ( 140, 2 ) =  272
fallxy     ( 140, 3 ) =  271
fallxy     ( 140, 4 ) =  269
fallxy     ( 140, 5 ) =  269
fallxy     ( 140, 6 ) =  269
fallxy     ( 140, 7 ) =  279
fallxy     ( 140, 8 ) =  279
xyphendoy1 ( 140, 1 ) =    0
xyphendoy1 ( 140, 2 ) =    0
xyphendoy1 ( 140, 3 ) =    0
xyphendoy1 ( 140, 4 ) =    0
xyphendoy1 ( 140, 5 ) =    0
xyphendoy1 ( 140, 6 ) =    0
xyphendoy1 ( 140, 7 ) =    0
xyphendoy1 ( 140, 8 ) =    0
      
 ! xlat         13.7500000000000     
greenupxy  ( 139, 1 ) =  136
greenupxy  ( 139, 2 ) =  138
greenupxy  ( 139, 3 ) =  189
greenupxy  ( 139, 4 ) =  174
greenupxy  ( 139, 5 ) =  174
greenupxy  ( 139, 6 ) =  174
greenupxy  ( 139, 7 ) =  157
greenupxy  ( 139, 8 ) =  157
fallxy     ( 139, 1 ) =  284
fallxy     ( 139, 2 ) =  272
fallxy     ( 139, 3 ) =  269
fallxy     ( 139, 4 ) =  269
fallxy     ( 139, 5 ) =  269
fallxy     ( 139, 6 ) =  269
fallxy     ( 139, 7 ) =  278
fallxy     ( 139, 8 ) =  278
xyphendoy1 ( 139, 1 ) =    0
xyphendoy1 ( 139, 2 ) =    0
xyphendoy1 ( 139, 3 ) =    0
xyphendoy1 ( 139, 4 ) =    0
xyphendoy1 ( 139, 5 ) =    0
xyphendoy1 ( 139, 6 ) =    0
xyphendoy1 ( 139, 7 ) =    0
xyphendoy1 ( 139, 8 ) =    0
      
 ! xlat         13.2500000000000     
greenupxy  ( 138, 1 ) =  134
greenupxy  ( 138, 2 ) =  138
greenupxy  ( 138, 3 ) =  185
greenupxy  ( 138, 4 ) =  174
greenupxy  ( 138, 5 ) =  174
greenupxy  ( 138, 6 ) =  174
greenupxy  ( 138, 7 ) =  159
greenupxy  ( 138, 8 ) =  159
fallxy     ( 138, 1 ) =  284
fallxy     ( 138, 2 ) =  273
fallxy     ( 138, 3 ) =  268
fallxy     ( 138, 4 ) =  269
fallxy     ( 138, 5 ) =  269
fallxy     ( 138, 6 ) =  269
fallxy     ( 138, 7 ) =  278
fallxy     ( 138, 8 ) =  278
xyphendoy1 ( 138, 1 ) =    0
xyphendoy1 ( 138, 2 ) =    0
xyphendoy1 ( 138, 3 ) =    0
xyphendoy1 ( 138, 4 ) =    0
xyphendoy1 ( 138, 5 ) =    0
xyphendoy1 ( 138, 6 ) =    0
xyphendoy1 ( 138, 7 ) =    0
xyphendoy1 ( 138, 8 ) =    0
      
 ! xlat         12.7500000000000     
greenupxy  ( 137, 1 ) =  134
greenupxy  ( 137, 2 ) =  138
greenupxy  ( 137, 3 ) =  178
greenupxy  ( 137, 4 ) =  173
greenupxy  ( 137, 5 ) =  173
greenupxy  ( 137, 6 ) =  173
greenupxy  ( 137, 7 ) =  163
greenupxy  ( 137, 8 ) =  163
fallxy     ( 137, 1 ) =  284
fallxy     ( 137, 2 ) =  273
fallxy     ( 137, 3 ) =  266
fallxy     ( 137, 4 ) =  269
fallxy     ( 137, 5 ) =  269
fallxy     ( 137, 6 ) =  269
fallxy     ( 137, 7 ) =  278
fallxy     ( 137, 8 ) =  278
xyphendoy1 ( 137, 1 ) =    0
xyphendoy1 ( 137, 2 ) =    0
xyphendoy1 ( 137, 3 ) =    0
xyphendoy1 ( 137, 4 ) =    0
xyphendoy1 ( 137, 5 ) =    0
xyphendoy1 ( 137, 6 ) =    0
xyphendoy1 ( 137, 7 ) =    0
xyphendoy1 ( 137, 8 ) =    0
      
 ! xlat         12.2500000000000     
greenupxy  ( 136, 1 ) =  131
greenupxy  ( 136, 2 ) =  137
greenupxy  ( 136, 3 ) =  176
greenupxy  ( 136, 4 ) =  167
greenupxy  ( 136, 5 ) =  167
greenupxy  ( 136, 6 ) =  167
greenupxy  ( 136, 7 ) =  166
greenupxy  ( 136, 8 ) =  166
fallxy     ( 136, 1 ) =  284
fallxy     ( 136, 2 ) =  273
fallxy     ( 136, 3 ) =  264
fallxy     ( 136, 4 ) =  268
fallxy     ( 136, 5 ) =  268
fallxy     ( 136, 6 ) =  268
fallxy     ( 136, 7 ) =  278
fallxy     ( 136, 8 ) =  278
xyphendoy1 ( 136, 1 ) =    0
xyphendoy1 ( 136, 2 ) =    0
xyphendoy1 ( 136, 3 ) =    0
xyphendoy1 ( 136, 4 ) =    0
xyphendoy1 ( 136, 5 ) =    0
xyphendoy1 ( 136, 6 ) =    0
xyphendoy1 ( 136, 7 ) =    0
xyphendoy1 ( 136, 8 ) =    0
      
 ! xlat         11.7500000000000     
greenupxy  ( 135, 1 ) =  130
greenupxy  ( 135, 2 ) =  136
greenupxy  ( 135, 3 ) =  176
greenupxy  ( 135, 4 ) =  161
greenupxy  ( 135, 5 ) =  161
greenupxy  ( 135, 6 ) =  161
greenupxy  ( 135, 7 ) =  169
greenupxy  ( 135, 8 ) =  169
fallxy     ( 135, 1 ) =  284
fallxy     ( 135, 2 ) =  274
fallxy     ( 135, 3 ) =  263
fallxy     ( 135, 4 ) =  267
fallxy     ( 135, 5 ) =  267
fallxy     ( 135, 6 ) =  267
fallxy     ( 135, 7 ) =  277
fallxy     ( 135, 8 ) =  277
xyphendoy1 ( 135, 1 ) =    0
xyphendoy1 ( 135, 2 ) =    0
xyphendoy1 ( 135, 3 ) =    0
xyphendoy1 ( 135, 4 ) =    0
xyphendoy1 ( 135, 5 ) =    0
xyphendoy1 ( 135, 6 ) =    0
xyphendoy1 ( 135, 7 ) =    0
xyphendoy1 ( 135, 8 ) =    0
      
 ! xlat         11.2500000000000     
greenupxy  ( 134, 1 ) =  138
greenupxy  ( 134, 2 ) =  134
greenupxy  ( 134, 3 ) =  174
greenupxy  ( 134, 4 ) =  154
greenupxy  ( 134, 5 ) =  154
greenupxy  ( 134, 6 ) =  154
greenupxy  ( 134, 7 ) =  171
greenupxy  ( 134, 8 ) =  171
fallxy     ( 134, 1 ) =  287
fallxy     ( 134, 2 ) =  274
fallxy     ( 134, 3 ) =  263
fallxy     ( 134, 4 ) =  268
fallxy     ( 134, 5 ) =  268
fallxy     ( 134, 6 ) =  268
fallxy     ( 134, 7 ) =  277
fallxy     ( 134, 8 ) =  277
xyphendoy1 ( 134, 1 ) =    0
xyphendoy1 ( 134, 2 ) =    0
xyphendoy1 ( 134, 3 ) =    0
xyphendoy1 ( 134, 4 ) =    0
xyphendoy1 ( 134, 5 ) =    0
xyphendoy1 ( 134, 6 ) =    0
xyphendoy1 ( 134, 7 ) =    0
xyphendoy1 ( 134, 8 ) =    0
      
 ! xlat         10.7500000000000     
greenupxy  ( 133, 1 ) =  135
greenupxy  ( 133, 2 ) =  132
greenupxy  ( 133, 3 ) =  172
greenupxy  ( 133, 4 ) =  147
greenupxy  ( 133, 5 ) =  147
greenupxy  ( 133, 6 ) =  147
greenupxy  ( 133, 7 ) =  169
greenupxy  ( 133, 8 ) =  169
fallxy     ( 133, 1 ) =  286
fallxy     ( 133, 2 ) =  273
fallxy     ( 133, 3 ) =  263
fallxy     ( 133, 4 ) =  268
fallxy     ( 133, 5 ) =  268
fallxy     ( 133, 6 ) =  268
fallxy     ( 133, 7 ) =  276
fallxy     ( 133, 8 ) =  276
xyphendoy1 ( 133, 1 ) =    0
xyphendoy1 ( 133, 2 ) =    0
xyphendoy1 ( 133, 3 ) =    0
xyphendoy1 ( 133, 4 ) =    0
xyphendoy1 ( 133, 5 ) =    0
xyphendoy1 ( 133, 6 ) =    0
xyphendoy1 ( 133, 7 ) =    0
xyphendoy1 ( 133, 8 ) =    0
      
 ! xlat         10.2500000000000     
greenupxy  ( 132, 1 ) =  137
greenupxy  ( 132, 2 ) =  128
greenupxy  ( 132, 3 ) =  158
greenupxy  ( 132, 4 ) =  139
greenupxy  ( 132, 5 ) =  139
greenupxy  ( 132, 6 ) =  139
greenupxy  ( 132, 7 ) =  165
greenupxy  ( 132, 8 ) =  165
fallxy     ( 132, 1 ) =  286
fallxy     ( 132, 2 ) =  273
fallxy     ( 132, 3 ) =  249
fallxy     ( 132, 4 ) =  268
fallxy     ( 132, 5 ) =  268
fallxy     ( 132, 6 ) =  268
fallxy     ( 132, 7 ) =  274
fallxy     ( 132, 8 ) =  274
xyphendoy1 ( 132, 1 ) =    0
xyphendoy1 ( 132, 2 ) =    0
xyphendoy1 ( 132, 3 ) =    0
xyphendoy1 ( 132, 4 ) =    0
xyphendoy1 ( 132, 5 ) =    0
xyphendoy1 ( 132, 6 ) =    0
xyphendoy1 ( 132, 7 ) =    0
xyphendoy1 ( 132, 8 ) =    0
      
 ! xlat         9.75000000000000     
greenupxy  ( 131, 1 ) =  129
greenupxy  ( 131, 2 ) =  123
greenupxy  ( 131, 3 ) =  145
greenupxy  ( 131, 4 ) =  133
greenupxy  ( 131, 5 ) =  133
greenupxy  ( 131, 6 ) =  133
greenupxy  ( 131, 7 ) =  158
greenupxy  ( 131, 8 ) =  158
fallxy     ( 131, 1 ) =  286
fallxy     ( 131, 2 ) =  272
fallxy     ( 131, 3 ) =  235
fallxy     ( 131, 4 ) =  268
fallxy     ( 131, 5 ) =  268
fallxy     ( 131, 6 ) =  268
fallxy     ( 131, 7 ) =  273
fallxy     ( 131, 8 ) =  273
xyphendoy1 ( 131, 1 ) =    0
xyphendoy1 ( 131, 2 ) =    0
xyphendoy1 ( 131, 3 ) =    0
xyphendoy1 ( 131, 4 ) =    0
xyphendoy1 ( 131, 5 ) =    0
xyphendoy1 ( 131, 6 ) =    0
xyphendoy1 ( 131, 7 ) =    0
xyphendoy1 ( 131, 8 ) =    0
      
 ! xlat         9.25000000000000     
greenupxy  ( 130, 1 ) =  122
greenupxy  ( 130, 2 ) =  118
greenupxy  ( 130, 3 ) =  132
greenupxy  ( 130, 4 ) =  126
greenupxy  ( 130, 5 ) =  126
greenupxy  ( 130, 6 ) =  126
greenupxy  ( 130, 7 ) =  152
greenupxy  ( 130, 8 ) =  152
fallxy     ( 130, 1 ) =  285
fallxy     ( 130, 2 ) =  272
fallxy     ( 130, 3 ) =  220
fallxy     ( 130, 4 ) =  268
fallxy     ( 130, 5 ) =  268
fallxy     ( 130, 6 ) =  268
fallxy     ( 130, 7 ) =  274
fallxy     ( 130, 8 ) =  274
xyphendoy1 ( 130, 1 ) =    0
xyphendoy1 ( 130, 2 ) =    0
xyphendoy1 ( 130, 3 ) =    0
xyphendoy1 ( 130, 4 ) =    0
xyphendoy1 ( 130, 5 ) =    0
xyphendoy1 ( 130, 6 ) =    0
xyphendoy1 ( 130, 7 ) =    0
xyphendoy1 ( 130, 8 ) =    0
      
 ! xlat         8.75000000000000     
greenupxy  ( 129, 1 ) =  122
greenupxy  ( 129, 2 ) =  112
greenupxy  ( 129, 3 ) =  127
greenupxy  ( 129, 4 ) =  119
greenupxy  ( 129, 5 ) =  119
greenupxy  ( 129, 6 ) =  119
greenupxy  ( 129, 7 ) =  145
greenupxy  ( 129, 8 ) =  145
fallxy     ( 129, 1 ) =  285
fallxy     ( 129, 2 ) =  272
fallxy     ( 129, 3 ) =  218
fallxy     ( 129, 4 ) =  268
fallxy     ( 129, 5 ) =  268
fallxy     ( 129, 6 ) =  268
fallxy     ( 129, 7 ) =  273
fallxy     ( 129, 8 ) =  273
xyphendoy1 ( 129, 1 ) =    0
xyphendoy1 ( 129, 2 ) =    0
xyphendoy1 ( 129, 3 ) =    0
xyphendoy1 ( 129, 4 ) =    0
xyphendoy1 ( 129, 5 ) =    0
xyphendoy1 ( 129, 6 ) =    0
xyphendoy1 ( 129, 7 ) =    0
xyphendoy1 ( 129, 8 ) =    0
      
 ! xlat         8.25000000000000     
greenupxy  ( 128, 1 ) =  116
greenupxy  ( 128, 2 ) =  108
greenupxy  ( 128, 3 ) =  114
greenupxy  ( 128, 4 ) =  111
greenupxy  ( 128, 5 ) =  111
greenupxy  ( 128, 6 ) =  111
greenupxy  ( 128, 7 ) =  134
greenupxy  ( 128, 8 ) =  134
fallxy     ( 128, 1 ) =  285
fallxy     ( 128, 2 ) =  271
fallxy     ( 128, 3 ) =  205
fallxy     ( 128, 4 ) =  268
fallxy     ( 128, 5 ) =  268
fallxy     ( 128, 6 ) =  268
fallxy     ( 128, 7 ) =  272
fallxy     ( 128, 8 ) =  272
xyphendoy1 ( 128, 1 ) =    0
xyphendoy1 ( 128, 2 ) =    0
xyphendoy1 ( 128, 3 ) =    0
xyphendoy1 ( 128, 4 ) =    0
xyphendoy1 ( 128, 5 ) =    0
xyphendoy1 ( 128, 6 ) =    0
xyphendoy1 ( 128, 7 ) =    0
xyphendoy1 ( 128, 8 ) =    0
      
 ! xlat         7.75000000000000     
greenupxy  ( 127, 1 ) =  114
greenupxy  ( 127, 2 ) =  103
greenupxy  ( 127, 3 ) =  102
greenupxy  ( 127, 4 ) =  102
greenupxy  ( 127, 5 ) =  102
greenupxy  ( 127, 6 ) =  102
greenupxy  ( 127, 7 ) =  129
greenupxy  ( 127, 8 ) =  129
fallxy     ( 127, 1 ) =  286
fallxy     ( 127, 2 ) =  271
fallxy     ( 127, 3 ) =  191
fallxy     ( 127, 4 ) =  267
fallxy     ( 127, 5 ) =  267
fallxy     ( 127, 6 ) =  267
fallxy     ( 127, 7 ) =  272
fallxy     ( 127, 8 ) =  272
xyphendoy1 ( 127, 1 ) =    0
xyphendoy1 ( 127, 2 ) =    0
xyphendoy1 ( 127, 3 ) =    0
xyphendoy1 ( 127, 4 ) =    0
xyphendoy1 ( 127, 5 ) =    0
xyphendoy1 ( 127, 6 ) =    0
xyphendoy1 ( 127, 7 ) =    0
xyphendoy1 ( 127, 8 ) =    0
      
 ! xlat         7.25000000000000     
greenupxy  ( 126, 1 ) =  110
greenupxy  ( 126, 2 ) =  100
greenupxy  ( 126, 3 ) =   90
greenupxy  ( 126, 4 ) =   93
greenupxy  ( 126, 5 ) =   93
greenupxy  ( 126, 6 ) =   93
greenupxy  ( 126, 7 ) =  124
greenupxy  ( 126, 8 ) =  124
fallxy     ( 126, 1 ) =  286
fallxy     ( 126, 2 ) =  271
fallxy     ( 126, 3 ) =  185
fallxy     ( 126, 4 ) =  267
fallxy     ( 126, 5 ) =  267
fallxy     ( 126, 6 ) =  267
fallxy     ( 126, 7 ) =  269
fallxy     ( 126, 8 ) =  269
xyphendoy1 ( 126, 1 ) =    0
xyphendoy1 ( 126, 2 ) =    0
xyphendoy1 ( 126, 3 ) =    0
xyphendoy1 ( 126, 4 ) =    0
xyphendoy1 ( 126, 5 ) =    0
xyphendoy1 ( 126, 6 ) =    0
xyphendoy1 ( 126, 7 ) =    0
xyphendoy1 ( 126, 8 ) =    0
      
 ! xlat         6.75000000000000     
greenupxy  ( 125, 1 ) =  107
greenupxy  ( 125, 2 ) =   96
greenupxy  ( 125, 3 ) =   79
greenupxy  ( 125, 4 ) =   85
greenupxy  ( 125, 5 ) =   85
greenupxy  ( 125, 6 ) =   85
greenupxy  ( 125, 7 ) =  119
greenupxy  ( 125, 8 ) =  119
fallxy     ( 125, 1 ) =  284
fallxy     ( 125, 2 ) =  271
fallxy     ( 125, 3 ) =  178
fallxy     ( 125, 4 ) =  267
fallxy     ( 125, 5 ) =  267
fallxy     ( 125, 6 ) =  267
fallxy     ( 125, 7 ) =  270
fallxy     ( 125, 8 ) =  270
xyphendoy1 ( 125, 1 ) =    0
xyphendoy1 ( 125, 2 ) =    0
xyphendoy1 ( 125, 3 ) =    0
xyphendoy1 ( 125, 4 ) =    0
xyphendoy1 ( 125, 5 ) =    0
xyphendoy1 ( 125, 6 ) =    0
xyphendoy1 ( 125, 7 ) =    0
xyphendoy1 ( 125, 8 ) =    0
      
 ! xlat         6.25000000000000     
greenupxy  ( 124, 1 ) =  103
greenupxy  ( 124, 2 ) =   93
greenupxy  ( 124, 3 ) =   70
greenupxy  ( 124, 4 ) =   80
greenupxy  ( 124, 5 ) =   80
greenupxy  ( 124, 6 ) =   80
greenupxy  ( 124, 7 ) =  115
greenupxy  ( 124, 8 ) =  115
fallxy     ( 124, 1 ) =  283
fallxy     ( 124, 2 ) =  272
fallxy     ( 124, 3 ) =  162
fallxy     ( 124, 4 ) =  266
fallxy     ( 124, 5 ) =  266
fallxy     ( 124, 6 ) =  266
fallxy     ( 124, 7 ) =  271
fallxy     ( 124, 8 ) =  271
xyphendoy1 ( 124, 1 ) =    0
xyphendoy1 ( 124, 2 ) =    0
xyphendoy1 ( 124, 3 ) =    0
xyphendoy1 ( 124, 4 ) =    0
xyphendoy1 ( 124, 5 ) =    0
xyphendoy1 ( 124, 6 ) =    0
xyphendoy1 ( 124, 7 ) =    0
xyphendoy1 ( 124, 8 ) =    0
      
 ! xlat         5.75000000000000     
greenupxy  ( 123, 1 ) =  102
greenupxy  ( 123, 2 ) =   91
greenupxy  ( 123, 3 ) =   62
greenupxy  ( 123, 4 ) =   75
greenupxy  ( 123, 5 ) =   75
greenupxy  ( 123, 6 ) =   75
greenupxy  ( 123, 7 ) =  110
greenupxy  ( 123, 8 ) =  110
fallxy     ( 123, 1 ) =  278
fallxy     ( 123, 2 ) =  272
fallxy     ( 123, 3 ) =  148
fallxy     ( 123, 4 ) =  266
fallxy     ( 123, 5 ) =  266
fallxy     ( 123, 6 ) =  266
fallxy     ( 123, 7 ) =  272
fallxy     ( 123, 8 ) =  272
xyphendoy1 ( 123, 1 ) =    0
xyphendoy1 ( 123, 2 ) =    0
xyphendoy1 ( 123, 3 ) =    0
xyphendoy1 ( 123, 4 ) =    0
xyphendoy1 ( 123, 5 ) =    0
xyphendoy1 ( 123, 6 ) =    0
xyphendoy1 ( 123, 7 ) =    0
xyphendoy1 ( 123, 8 ) =    0
      
 ! xlat         5.25000000000000     
greenupxy  ( 122, 1 ) =   96
greenupxy  ( 122, 2 ) =   88
greenupxy  ( 122, 3 ) =   52
greenupxy  ( 122, 4 ) =   69
greenupxy  ( 122, 5 ) =   69
greenupxy  ( 122, 6 ) =   69
greenupxy  ( 122, 7 ) =  105
greenupxy  ( 122, 8 ) =  105
fallxy     ( 122, 1 ) =  276
fallxy     ( 122, 2 ) =  269
fallxy     ( 122, 3 ) =  133
fallxy     ( 122, 4 ) =  265
fallxy     ( 122, 5 ) =  265
fallxy     ( 122, 6 ) =  265
fallxy     ( 122, 7 ) =  271
fallxy     ( 122, 8 ) =  271
xyphendoy1 ( 122, 1 ) =    0
xyphendoy1 ( 122, 2 ) =    0
xyphendoy1 ( 122, 3 ) =    0
xyphendoy1 ( 122, 4 ) =    0
xyphendoy1 ( 122, 5 ) =    0
xyphendoy1 ( 122, 6 ) =    0
xyphendoy1 ( 122, 7 ) =    0
xyphendoy1 ( 122, 8 ) =    0
      
 ! xlat         4.75000000000000     
greenupxy  ( 121, 1 ) =   93
greenupxy  ( 121, 2 ) =   84
greenupxy  ( 121, 3 ) =   41
greenupxy  ( 121, 4 ) =   61
greenupxy  ( 121, 5 ) =   61
greenupxy  ( 121, 6 ) =   61
greenupxy  ( 121, 7 ) =   96
greenupxy  ( 121, 8 ) =   96
fallxy     ( 121, 1 ) =  276
fallxy     ( 121, 2 ) =  267
fallxy     ( 121, 3 ) =  117
fallxy     ( 121, 4 ) =  265
fallxy     ( 121, 5 ) =  265
fallxy     ( 121, 6 ) =  265
fallxy     ( 121, 7 ) =  269
fallxy     ( 121, 8 ) =  269
xyphendoy1 ( 121, 1 ) =    0
xyphendoy1 ( 121, 2 ) =    0
xyphendoy1 ( 121, 3 ) =    0
xyphendoy1 ( 121, 4 ) =    0
xyphendoy1 ( 121, 5 ) =    0
xyphendoy1 ( 121, 6 ) =    0
xyphendoy1 ( 121, 7 ) =    0
xyphendoy1 ( 121, 8 ) =    0
      
 ! xlat         4.25000000000000     
greenupxy  ( 120, 1 ) =   93
greenupxy  ( 120, 2 ) =   71
greenupxy  ( 120, 3 ) =   34
greenupxy  ( 120, 4 ) =   55
greenupxy  ( 120, 5 ) =   55
greenupxy  ( 120, 6 ) =   55
greenupxy  ( 120, 7 ) =   85
greenupxy  ( 120, 8 ) =   85
fallxy     ( 120, 1 ) =  276
fallxy     ( 120, 2 ) =  266
fallxy     ( 120, 3 ) =  102
fallxy     ( 120, 4 ) =  264
fallxy     ( 120, 5 ) =  264
fallxy     ( 120, 6 ) =  264
fallxy     ( 120, 7 ) =  267
fallxy     ( 120, 8 ) =  267
xyphendoy1 ( 120, 1 ) =    0
xyphendoy1 ( 120, 2 ) =    0
xyphendoy1 ( 120, 3 ) =    0
xyphendoy1 ( 120, 4 ) =    0
xyphendoy1 ( 120, 5 ) =    0
xyphendoy1 ( 120, 6 ) =    0
xyphendoy1 ( 120, 7 ) =    0
xyphendoy1 ( 120, 8 ) =    0
      
 ! xlat         3.75000000000000     
greenupxy  ( 119, 1 ) =   90
greenupxy  ( 119, 2 ) =   58
greenupxy  ( 119, 3 ) =   32
greenupxy  ( 119, 4 ) =   48
greenupxy  ( 119, 5 ) =   48
greenupxy  ( 119, 6 ) =   48
greenupxy  ( 119, 7 ) =   79
greenupxy  ( 119, 8 ) =   79
fallxy     ( 119, 1 ) =  259
fallxy     ( 119, 2 ) =  266
fallxy     ( 119, 3 ) =   87
fallxy     ( 119, 4 ) =  265
fallxy     ( 119, 5 ) =  265
fallxy     ( 119, 6 ) =  265
fallxy     ( 119, 7 ) =  267
fallxy     ( 119, 8 ) =  267
xyphendoy1 ( 119, 1 ) =    0
xyphendoy1 ( 119, 2 ) =    0
xyphendoy1 ( 119, 3 ) =    0
xyphendoy1 ( 119, 4 ) =    0
xyphendoy1 ( 119, 5 ) =    0
xyphendoy1 ( 119, 6 ) =    0
xyphendoy1 ( 119, 7 ) =    0
xyphendoy1 ( 119, 8 ) =    0
      
 ! xlat         3.25000000000000     
greenupxy  ( 118, 1 ) =   90
greenupxy  ( 118, 2 ) =   51
greenupxy  ( 118, 3 ) =   13
greenupxy  ( 118, 4 ) =   38
greenupxy  ( 118, 5 ) =   38
greenupxy  ( 118, 6 ) =   38
greenupxy  ( 118, 7 ) =   67
greenupxy  ( 118, 8 ) =   67
fallxy     ( 118, 1 ) =  259
fallxy     ( 118, 2 ) =  264
fallxy     ( 118, 3 ) =   84
fallxy     ( 118, 4 ) =  265
fallxy     ( 118, 5 ) =  265
fallxy     ( 118, 6 ) =  265
fallxy     ( 118, 7 ) =  266
fallxy     ( 118, 8 ) =  266
xyphendoy1 ( 118, 1 ) =    0
xyphendoy1 ( 118, 2 ) =    0
xyphendoy1 ( 118, 3 ) =    0
xyphendoy1 ( 118, 4 ) =    0
xyphendoy1 ( 118, 5 ) =    0
xyphendoy1 ( 118, 6 ) =    0
xyphendoy1 ( 118, 7 ) =    0
xyphendoy1 ( 118, 8 ) =    0
      
 ! xlat         2.75000000000000     
greenupxy  ( 117, 1 ) =   81
greenupxy  ( 117, 2 ) =   42
greenupxy  ( 117, 3 ) =    8
greenupxy  ( 117, 4 ) =   34
greenupxy  ( 117, 5 ) =   34
greenupxy  ( 117, 6 ) =   34
greenupxy  ( 117, 7 ) =   57
greenupxy  ( 117, 8 ) =   57
fallxy     ( 117, 1 ) =  251
fallxy     ( 117, 2 ) =  253
fallxy     ( 117, 3 ) =   83
fallxy     ( 117, 4 ) =  241
fallxy     ( 117, 5 ) =  241
fallxy     ( 117, 6 ) =  241
fallxy     ( 117, 7 ) =  258
fallxy     ( 117, 8 ) =  258
xyphendoy1 ( 117, 1 ) =    0
xyphendoy1 ( 117, 2 ) =    0
xyphendoy1 ( 117, 3 ) =    0
xyphendoy1 ( 117, 4 ) =    0
xyphendoy1 ( 117, 5 ) =    0
xyphendoy1 ( 117, 6 ) =    0
xyphendoy1 ( 117, 7 ) =    0
xyphendoy1 ( 117, 8 ) =    0
      
 ! xlat         2.25000000000000     
greenupxy  ( 116, 1 ) =   67
greenupxy  ( 116, 2 ) =   35
greenupxy  ( 116, 3 ) =    3
greenupxy  ( 116, 4 ) =   26
greenupxy  ( 116, 5 ) =   26
greenupxy  ( 116, 6 ) =   26
greenupxy  ( 116, 7 ) =   48
greenupxy  ( 116, 8 ) =   48
fallxy     ( 116, 1 ) =  245
fallxy     ( 116, 2 ) =  240
fallxy     ( 116, 3 ) =   88
fallxy     ( 116, 4 ) =  226
fallxy     ( 116, 5 ) =  226
fallxy     ( 116, 6 ) =  226
fallxy     ( 116, 7 ) =  249
fallxy     ( 116, 8 ) =  249
xyphendoy1 ( 116, 1 ) =    0
xyphendoy1 ( 116, 2 ) =    0
xyphendoy1 ( 116, 3 ) =    0
xyphendoy1 ( 116, 4 ) =    0
xyphendoy1 ( 116, 5 ) =    0
xyphendoy1 ( 116, 6 ) =    0
xyphendoy1 ( 116, 7 ) =    0
xyphendoy1 ( 116, 8 ) =    0
      
 ! xlat         1.75000000000000     
greenupxy  ( 115, 1 ) =   51
greenupxy  ( 115, 2 ) =   24
greenupxy  ( 115, 3 ) =  363
greenupxy  ( 115, 4 ) =   19
greenupxy  ( 115, 5 ) =   19
greenupxy  ( 115, 6 ) =   19
greenupxy  ( 115, 7 ) =   38
greenupxy  ( 115, 8 ) =   38
fallxy     ( 115, 1 ) =  231
fallxy     ( 115, 2 ) =  226
fallxy     ( 115, 3 ) =   87
fallxy     ( 115, 4 ) =  210
fallxy     ( 115, 5 ) =  210
fallxy     ( 115, 6 ) =  210
fallxy     ( 115, 7 ) =  242
fallxy     ( 115, 8 ) =  242
xyphendoy1 ( 115, 1 ) =    0
xyphendoy1 ( 115, 2 ) =    0
xyphendoy1 ( 115, 3 ) =    2
xyphendoy1 ( 115, 4 ) =    0
xyphendoy1 ( 115, 5 ) =    0
xyphendoy1 ( 115, 6 ) =    0
xyphendoy1 ( 115, 7 ) =    0
xyphendoy1 ( 115, 8 ) =    0
      
 ! xlat         1.25000000000000     
greenupxy  ( 114, 1 ) =   58
greenupxy  ( 114, 2 ) =   15
greenupxy  ( 114, 3 ) =  360
greenupxy  ( 114, 4 ) =   10
greenupxy  ( 114, 5 ) =   10
greenupxy  ( 114, 6 ) =   10
greenupxy  ( 114, 7 ) =   30
greenupxy  ( 114, 8 ) =   30
fallxy     ( 114, 1 ) =  231
fallxy     ( 114, 2 ) =  211
fallxy     ( 114, 3 ) =   91
fallxy     ( 114, 4 ) =  194
fallxy     ( 114, 5 ) =  194
fallxy     ( 114, 6 ) =  194
fallxy     ( 114, 7 ) =  233
fallxy     ( 114, 8 ) =  233
xyphendoy1 ( 114, 1 ) =    0
xyphendoy1 ( 114, 2 ) =    0
xyphendoy1 ( 114, 3 ) =    2
xyphendoy1 ( 114, 4 ) =    0
xyphendoy1 ( 114, 5 ) =    0
xyphendoy1 ( 114, 6 ) =    0
xyphendoy1 ( 114, 7 ) =    0
xyphendoy1 ( 114, 8 ) =    0
      
 ! xlat        0.750000000000000     
greenupxy  ( 113, 1 ) =   67
greenupxy  ( 113, 2 ) =    5
greenupxy  ( 113, 3 ) =  358
greenupxy  ( 113, 4 ) =    4
greenupxy  ( 113, 5 ) =    4
greenupxy  ( 113, 6 ) =    4
greenupxy  ( 113, 7 ) =   21
greenupxy  ( 113, 8 ) =   21
fallxy     ( 113, 1 ) =  207
fallxy     ( 113, 2 ) =  196
fallxy     ( 113, 3 ) =   95
fallxy     ( 113, 4 ) =  176
fallxy     ( 113, 5 ) =  176
fallxy     ( 113, 6 ) =  176
fallxy     ( 113, 7 ) =  222
fallxy     ( 113, 8 ) =  222
xyphendoy1 ( 113, 1 ) =    0
xyphendoy1 ( 113, 2 ) =    0
xyphendoy1 ( 113, 3 ) =    2
xyphendoy1 ( 113, 4 ) =    0
xyphendoy1 ( 113, 5 ) =    0
xyphendoy1 ( 113, 6 ) =    0
xyphendoy1 ( 113, 7 ) =    0
xyphendoy1 ( 113, 8 ) =    0
      
 ! xlat        0.250000000000000     
greenupxy  ( 112, 1 ) =   67
greenupxy  ( 112, 2 ) =  365
greenupxy  ( 112, 3 ) =  357
greenupxy  ( 112, 4 ) =  364
greenupxy  ( 112, 5 ) =  364
greenupxy  ( 112, 6 ) =  364
greenupxy  ( 112, 7 ) =   15
greenupxy  ( 112, 8 ) =   15
fallxy     ( 112, 1 ) =  207
fallxy     ( 112, 2 ) =  179
fallxy     ( 112, 3 ) =   85
fallxy     ( 112, 4 ) =  157
fallxy     ( 112, 5 ) =  157
fallxy     ( 112, 6 ) =  157
fallxy     ( 112, 7 ) =  213
fallxy     ( 112, 8 ) =  213
xyphendoy1 ( 112, 1 ) =    0
xyphendoy1 ( 112, 2 ) =    2
xyphendoy1 ( 112, 3 ) =    2
xyphendoy1 ( 112, 4 ) =    2
xyphendoy1 ( 112, 5 ) =    2
xyphendoy1 ( 112, 6 ) =    2
xyphendoy1 ( 112, 7 ) =    0
xyphendoy1 ( 112, 8 ) =    0
      
 ! xlat       -0.250000000000000     
greenupxy  ( 111, 1 ) =   67
greenupxy  ( 111, 2 ) =  359
greenupxy  ( 111, 3 ) =  355
greenupxy  ( 111, 4 ) =  360
greenupxy  ( 111, 5 ) =  360
greenupxy  ( 111, 6 ) =  360
greenupxy  ( 111, 7 ) =   13
greenupxy  ( 111, 8 ) =   13
fallxy     ( 111, 1 ) =  207
fallxy     ( 111, 2 ) =  164
fallxy     ( 111, 3 ) =   86
fallxy     ( 111, 4 ) =  140
fallxy     ( 111, 5 ) =  140
fallxy     ( 111, 6 ) =  140
fallxy     ( 111, 7 ) =  202
fallxy     ( 111, 8 ) =  202
xyphendoy1 ( 111, 1 ) =    0
xyphendoy1 ( 111, 2 ) =    2
xyphendoy1 ( 111, 3 ) =    2
xyphendoy1 ( 111, 4 ) =    2
xyphendoy1 ( 111, 5 ) =    2
xyphendoy1 ( 111, 6 ) =    2
xyphendoy1 ( 111, 7 ) =    0
xyphendoy1 ( 111, 8 ) =    0
      
 ! xlat       -0.750000000000000     
greenupxy  ( 110, 1 ) =   67
greenupxy  ( 110, 2 ) =  345
greenupxy  ( 110, 3 ) =  353
greenupxy  ( 110, 4 ) =  358
greenupxy  ( 110, 5 ) =  358
greenupxy  ( 110, 6 ) =  358
greenupxy  ( 110, 7 ) =    5
greenupxy  ( 110, 8 ) =    5
fallxy     ( 110, 1 ) =  207
fallxy     ( 110, 2 ) =  148
fallxy     ( 110, 3 ) =   90
fallxy     ( 110, 4 ) =  124
fallxy     ( 110, 5 ) =  124
fallxy     ( 110, 6 ) =  124
fallxy     ( 110, 7 ) =  192
fallxy     ( 110, 8 ) =  192
xyphendoy1 ( 110, 1 ) =    0
xyphendoy1 ( 110, 2 ) =    2
xyphendoy1 ( 110, 3 ) =    2
xyphendoy1 ( 110, 4 ) =    2
xyphendoy1 ( 110, 5 ) =    2
xyphendoy1 ( 110, 6 ) =    2
xyphendoy1 ( 110, 7 ) =    0
xyphendoy1 ( 110, 8 ) =    0
      
 ! xlat        -1.25000000000000     
greenupxy  ( 109, 1 ) =   67
greenupxy  ( 109, 2 ) =  329
greenupxy  ( 109, 3 ) =  353
greenupxy  ( 109, 4 ) =  356
greenupxy  ( 109, 5 ) =  356
greenupxy  ( 109, 6 ) =  356
greenupxy  ( 109, 7 ) =  363
greenupxy  ( 109, 8 ) =  363
fallxy     ( 109, 1 ) =  207
fallxy     ( 109, 2 ) =  143
fallxy     ( 109, 3 ) =   93
fallxy     ( 109, 4 ) =  124
fallxy     ( 109, 5 ) =  124
fallxy     ( 109, 6 ) =  124
fallxy     ( 109, 7 ) =  187
fallxy     ( 109, 8 ) =  187
xyphendoy1 ( 109, 1 ) =    0
xyphendoy1 ( 109, 2 ) =    2
xyphendoy1 ( 109, 3 ) =    2
xyphendoy1 ( 109, 4 ) =    2
xyphendoy1 ( 109, 5 ) =    2
xyphendoy1 ( 109, 6 ) =    2
xyphendoy1 ( 109, 7 ) =    2
xyphendoy1 ( 109, 8 ) =    2
      
 ! xlat        -1.75000000000000     
greenupxy  ( 108, 1 ) =   48
greenupxy  ( 108, 2 ) =  314
greenupxy  ( 108, 3 ) =  353
greenupxy  ( 108, 4 ) =  355
greenupxy  ( 108, 5 ) =  355
greenupxy  ( 108, 6 ) =  355
greenupxy  ( 108, 7 ) =  354
greenupxy  ( 108, 8 ) =  354
fallxy     ( 108, 1 ) =  187
fallxy     ( 108, 2 ) =  139
fallxy     ( 108, 3 ) =   95
fallxy     ( 108, 4 ) =  124
fallxy     ( 108, 5 ) =  124
fallxy     ( 108, 6 ) =  124
fallxy     ( 108, 7 ) =  177
fallxy     ( 108, 8 ) =  177
xyphendoy1 ( 108, 1 ) =    0
xyphendoy1 ( 108, 2 ) =    2
xyphendoy1 ( 108, 3 ) =    2
xyphendoy1 ( 108, 4 ) =    2
xyphendoy1 ( 108, 5 ) =    2
xyphendoy1 ( 108, 6 ) =    2
xyphendoy1 ( 108, 7 ) =    2
xyphendoy1 ( 108, 8 ) =    2
      
 ! xlat        -2.25000000000000     
greenupxy  ( 107, 1 ) =   48
greenupxy  ( 107, 2 ) =  297
greenupxy  ( 107, 3 ) =  351
greenupxy  ( 107, 4 ) =  355
greenupxy  ( 107, 5 ) =  355
greenupxy  ( 107, 6 ) =  355
greenupxy  ( 107, 7 ) =  345
greenupxy  ( 107, 8 ) =  345
fallxy     ( 107, 1 ) =  187
fallxy     ( 107, 2 ) =  135
fallxy     ( 107, 3 ) =   97
fallxy     ( 107, 4 ) =  124
fallxy     ( 107, 5 ) =  124
fallxy     ( 107, 6 ) =  124
fallxy     ( 107, 7 ) =  164
fallxy     ( 107, 8 ) =  164
xyphendoy1 ( 107, 1 ) =    0
xyphendoy1 ( 107, 2 ) =    2
xyphendoy1 ( 107, 3 ) =    2
xyphendoy1 ( 107, 4 ) =    2
xyphendoy1 ( 107, 5 ) =    2
xyphendoy1 ( 107, 6 ) =    2
xyphendoy1 ( 107, 7 ) =    2
xyphendoy1 ( 107, 8 ) =    2
      
 ! xlat        -2.75000000000000     
greenupxy  ( 106, 1 ) =   33
greenupxy  ( 106, 2 ) =  293
greenupxy  ( 106, 3 ) =  349
greenupxy  ( 106, 4 ) =  351
greenupxy  ( 106, 5 ) =  351
greenupxy  ( 106, 6 ) =  351
greenupxy  ( 106, 7 ) =  336
greenupxy  ( 106, 8 ) =  336
fallxy     ( 106, 1 ) =  162
fallxy     ( 106, 2 ) =  131
fallxy     ( 106, 3 ) =   98
fallxy     ( 106, 4 ) =  123
fallxy     ( 106, 5 ) =  123
fallxy     ( 106, 6 ) =  123
fallxy     ( 106, 7 ) =  152
fallxy     ( 106, 8 ) =  152
xyphendoy1 ( 106, 1 ) =    0
xyphendoy1 ( 106, 2 ) =    2
xyphendoy1 ( 106, 3 ) =    2
xyphendoy1 ( 106, 4 ) =    2
xyphendoy1 ( 106, 5 ) =    2
xyphendoy1 ( 106, 6 ) =    2
xyphendoy1 ( 106, 7 ) =    2
xyphendoy1 ( 106, 8 ) =    2
      
 ! xlat        -3.25000000000000     
greenupxy  ( 105, 1 ) =   26
greenupxy  ( 105, 2 ) =  289
greenupxy  ( 105, 3 ) =  347
greenupxy  ( 105, 4 ) =  347
greenupxy  ( 105, 5 ) =  347
greenupxy  ( 105, 6 ) =  347
greenupxy  ( 105, 7 ) =  327
greenupxy  ( 105, 8 ) =  327
fallxy     ( 105, 1 ) =  151
fallxy     ( 105, 2 ) =  119
fallxy     ( 105, 3 ) =  100
fallxy     ( 105, 4 ) =  122
fallxy     ( 105, 5 ) =  122
fallxy     ( 105, 6 ) =  122
fallxy     ( 105, 7 ) =  142
fallxy     ( 105, 8 ) =  142
xyphendoy1 ( 105, 1 ) =    0
xyphendoy1 ( 105, 2 ) =    2
xyphendoy1 ( 105, 3 ) =    2
xyphendoy1 ( 105, 4 ) =    2
xyphendoy1 ( 105, 5 ) =    2
xyphendoy1 ( 105, 6 ) =    2
xyphendoy1 ( 105, 7 ) =    2
xyphendoy1 ( 105, 8 ) =    2
      
 ! xlat        -3.75000000000000     
greenupxy  ( 104, 1 ) =   26
greenupxy  ( 104, 2 ) =  285
greenupxy  ( 104, 3 ) =  347
greenupxy  ( 104, 4 ) =  345
greenupxy  ( 104, 5 ) =  345
greenupxy  ( 104, 6 ) =  345
greenupxy  ( 104, 7 ) =  323
greenupxy  ( 104, 8 ) =  323
fallxy     ( 104, 1 ) =  151
fallxy     ( 104, 2 ) =  108
fallxy     ( 104, 3 ) =  102
fallxy     ( 104, 4 ) =  121
fallxy     ( 104, 5 ) =  121
fallxy     ( 104, 6 ) =  121
fallxy     ( 104, 7 ) =  132
fallxy     ( 104, 8 ) =  132
xyphendoy1 ( 104, 1 ) =    0
xyphendoy1 ( 104, 2 ) =    2
xyphendoy1 ( 104, 3 ) =    2
xyphendoy1 ( 104, 4 ) =    2
xyphendoy1 ( 104, 5 ) =    2
xyphendoy1 ( 104, 6 ) =    2
xyphendoy1 ( 104, 7 ) =    2
xyphendoy1 ( 104, 8 ) =    2
      
 ! xlat        -4.25000000000000     
greenupxy  ( 103, 1 ) =   26
greenupxy  ( 103, 2 ) =  281
greenupxy  ( 103, 3 ) =  345
greenupxy  ( 103, 4 ) =  343
greenupxy  ( 103, 5 ) =  343
greenupxy  ( 103, 6 ) =  343
greenupxy  ( 103, 7 ) =  320
greenupxy  ( 103, 8 ) =  320
fallxy     ( 103, 1 ) =  151
fallxy     ( 103, 2 ) =  106
fallxy     ( 103, 3 ) =  104
fallxy     ( 103, 4 ) =  120
fallxy     ( 103, 5 ) =  120
fallxy     ( 103, 6 ) =  120
fallxy     ( 103, 7 ) =  123
fallxy     ( 103, 8 ) =  123
xyphendoy1 ( 103, 1 ) =    0
xyphendoy1 ( 103, 2 ) =    2
xyphendoy1 ( 103, 3 ) =    2
xyphendoy1 ( 103, 4 ) =    2
xyphendoy1 ( 103, 5 ) =    2
xyphendoy1 ( 103, 6 ) =    2
xyphendoy1 ( 103, 7 ) =    2
xyphendoy1 ( 103, 8 ) =    2
      
 ! xlat        -4.75000000000000     
greenupxy  ( 102, 1 ) =  362
greenupxy  ( 102, 2 ) =  278
greenupxy  ( 102, 3 ) =  344
greenupxy  ( 102, 4 ) =  341
greenupxy  ( 102, 5 ) =  341
greenupxy  ( 102, 6 ) =  341
greenupxy  ( 102, 7 ) =  320
greenupxy  ( 102, 8 ) =  320
fallxy     ( 102, 1 ) =  133
fallxy     ( 102, 2 ) =  104
fallxy     ( 102, 3 ) =  107
fallxy     ( 102, 4 ) =  119
fallxy     ( 102, 5 ) =  119
fallxy     ( 102, 6 ) =  119
fallxy     ( 102, 7 ) =  122
fallxy     ( 102, 8 ) =  122
xyphendoy1 ( 102, 1 ) =    2
xyphendoy1 ( 102, 2 ) =    2
xyphendoy1 ( 102, 3 ) =    2
xyphendoy1 ( 102, 4 ) =    2
xyphendoy1 ( 102, 5 ) =    2
xyphendoy1 ( 102, 6 ) =    2
xyphendoy1 ( 102, 7 ) =    2
xyphendoy1 ( 102, 8 ) =    2
      
 ! xlat        -5.25000000000000     
greenupxy  ( 101, 1 ) =  362
greenupxy  ( 101, 2 ) =  271
greenupxy  ( 101, 3 ) =  347
greenupxy  ( 101, 4 ) =  343
greenupxy  ( 101, 5 ) =  343
greenupxy  ( 101, 6 ) =  343
greenupxy  ( 101, 7 ) =  318
greenupxy  ( 101, 8 ) =  318
fallxy     ( 101, 1 ) =  133
fallxy     ( 101, 2 ) =  102
fallxy     ( 101, 3 ) =  110
fallxy     ( 101, 4 ) =  119
fallxy     ( 101, 5 ) =  119
fallxy     ( 101, 6 ) =  119
fallxy     ( 101, 7 ) =  122
fallxy     ( 101, 8 ) =  122
xyphendoy1 ( 101, 1 ) =    2
xyphendoy1 ( 101, 2 ) =    2
xyphendoy1 ( 101, 3 ) =    2
xyphendoy1 ( 101, 4 ) =    2
xyphendoy1 ( 101, 5 ) =    2
xyphendoy1 ( 101, 6 ) =    2
xyphendoy1 ( 101, 7 ) =    2
xyphendoy1 ( 101, 8 ) =    2
      
 ! xlat        -5.75000000000000     
greenupxy  ( 100, 1 ) =  355
greenupxy  ( 100, 2 ) =  265
greenupxy  ( 100, 3 ) =  348
greenupxy  ( 100, 4 ) =  339
greenupxy  ( 100, 5 ) =  339
greenupxy  ( 100, 6 ) =  339
greenupxy  ( 100, 7 ) =  316
greenupxy  ( 100, 8 ) =  316
fallxy     ( 100, 1 ) =  131
fallxy     ( 100, 2 ) =   97
fallxy     ( 100, 3 ) =  112
fallxy     ( 100, 4 ) =  120
fallxy     ( 100, 5 ) =  120
fallxy     ( 100, 6 ) =  120
fallxy     ( 100, 7 ) =  118
fallxy     ( 100, 8 ) =  118
xyphendoy1 ( 100, 1 ) =    2
xyphendoy1 ( 100, 2 ) =    2
xyphendoy1 ( 100, 3 ) =    2
xyphendoy1 ( 100, 4 ) =    2
xyphendoy1 ( 100, 5 ) =    2
xyphendoy1 ( 100, 6 ) =    2
xyphendoy1 ( 100, 7 ) =    2
xyphendoy1 ( 100, 8 ) =    2
      
 ! xlat        -6.25000000000000     
greenupxy  (  99, 1 ) =  349
greenupxy  (  99, 2 ) =  260
greenupxy  (  99, 3 ) =  349
greenupxy  (  99, 4 ) =  338
greenupxy  (  99, 5 ) =  338
greenupxy  (  99, 6 ) =  338
greenupxy  (  99, 7 ) =  316
greenupxy  (  99, 8 ) =  316
fallxy     (  99, 1 ) =  131
fallxy     (  99, 2 ) =   92
fallxy     (  99, 3 ) =  109
fallxy     (  99, 4 ) =  117
fallxy     (  99, 5 ) =  117
fallxy     (  99, 6 ) =  117
fallxy     (  99, 7 ) =  116
fallxy     (  99, 8 ) =  116
xyphendoy1 (  99, 1 ) =    2
xyphendoy1 (  99, 2 ) =    2
xyphendoy1 (  99, 3 ) =    2
xyphendoy1 (  99, 4 ) =    2
xyphendoy1 (  99, 5 ) =    2
xyphendoy1 (  99, 6 ) =    2
xyphendoy1 (  99, 7 ) =    2
xyphendoy1 (  99, 8 ) =    2
      
 ! xlat        -6.75000000000000     
greenupxy  (  98, 1 ) =  349
greenupxy  (  98, 2 ) =  257
greenupxy  (  98, 3 ) =  350
greenupxy  (  98, 4 ) =  338
greenupxy  (  98, 5 ) =  338
greenupxy  (  98, 6 ) =  338
greenupxy  (  98, 7 ) =  316
greenupxy  (  98, 8 ) =  316
fallxy     (  98, 1 ) =  131
fallxy     (  98, 2 ) =   89
fallxy     (  98, 3 ) =  111
fallxy     (  98, 4 ) =  115
fallxy     (  98, 5 ) =  115
fallxy     (  98, 6 ) =  115
fallxy     (  98, 7 ) =  114
fallxy     (  98, 8 ) =  114
xyphendoy1 (  98, 1 ) =    2
xyphendoy1 (  98, 2 ) =    2
xyphendoy1 (  98, 3 ) =    2
xyphendoy1 (  98, 4 ) =    2
xyphendoy1 (  98, 5 ) =    2
xyphendoy1 (  98, 6 ) =    2
xyphendoy1 (  98, 7 ) =    2
xyphendoy1 (  98, 8 ) =    2
      
 ! xlat        -7.25000000000000     
greenupxy  (  97, 1 ) =  349
greenupxy  (  97, 2 ) =  255
greenupxy  (  97, 3 ) =  349
greenupxy  (  97, 4 ) =  337
greenupxy  (  97, 5 ) =  337
greenupxy  (  97, 6 ) =  337
greenupxy  (  97, 7 ) =  316
greenupxy  (  97, 8 ) =  316
fallxy     (  97, 1 ) =  131
fallxy     (  97, 2 ) =   86
fallxy     (  97, 3 ) =  109
fallxy     (  97, 4 ) =  113
fallxy     (  97, 5 ) =  113
fallxy     (  97, 6 ) =  113
fallxy     (  97, 7 ) =  112
fallxy     (  97, 8 ) =  112
xyphendoy1 (  97, 1 ) =    2
xyphendoy1 (  97, 2 ) =    2
xyphendoy1 (  97, 3 ) =    2
xyphendoy1 (  97, 4 ) =    2
xyphendoy1 (  97, 5 ) =    2
xyphendoy1 (  97, 6 ) =    2
xyphendoy1 (  97, 7 ) =    2
xyphendoy1 (  97, 8 ) =    2
      
 ! xlat        -7.75000000000000     
greenupxy  (  96, 1 ) =  344
greenupxy  (  96, 2 ) =  253
greenupxy  (  96, 3 ) =  346
greenupxy  (  96, 4 ) =  334
greenupxy  (  96, 5 ) =  334
greenupxy  (  96, 6 ) =  334
greenupxy  (  96, 7 ) =  315
greenupxy  (  96, 8 ) =  315
fallxy     (  96, 1 ) =  127
fallxy     (  96, 2 ) =   84
fallxy     (  96, 3 ) =  107
fallxy     (  96, 4 ) =  111
fallxy     (  96, 5 ) =  111
fallxy     (  96, 6 ) =  111
fallxy     (  96, 7 ) =  110
fallxy     (  96, 8 ) =  110
xyphendoy1 (  96, 1 ) =    2
xyphendoy1 (  96, 2 ) =    2
xyphendoy1 (  96, 3 ) =    2
xyphendoy1 (  96, 4 ) =    2
xyphendoy1 (  96, 5 ) =    2
xyphendoy1 (  96, 6 ) =    2
xyphendoy1 (  96, 7 ) =    2
xyphendoy1 (  96, 8 ) =    2
      
 ! xlat        -8.25000000000000     
greenupxy  (  95, 1 ) =  344
greenupxy  (  95, 2 ) =  254
greenupxy  (  95, 3 ) =  342
greenupxy  (  95, 4 ) =  331
greenupxy  (  95, 5 ) =  331
greenupxy  (  95, 6 ) =  331
greenupxy  (  95, 7 ) =  315
greenupxy  (  95, 8 ) =  315
fallxy     (  95, 1 ) =  127
fallxy     (  95, 2 ) =   81
fallxy     (  95, 3 ) =  107
fallxy     (  95, 4 ) =  110
fallxy     (  95, 5 ) =  110
fallxy     (  95, 6 ) =  110
fallxy     (  95, 7 ) =  106
fallxy     (  95, 8 ) =  106
xyphendoy1 (  95, 1 ) =    2
xyphendoy1 (  95, 2 ) =    2
xyphendoy1 (  95, 3 ) =    2
xyphendoy1 (  95, 4 ) =    2
xyphendoy1 (  95, 5 ) =    2
xyphendoy1 (  95, 6 ) =    2
xyphendoy1 (  95, 7 ) =    2
xyphendoy1 (  95, 8 ) =    2
      
 ! xlat        -8.75000000000000     
greenupxy  (  94, 1 ) =  343
greenupxy  (  94, 2 ) =  254
greenupxy  (  94, 3 ) =  339
greenupxy  (  94, 4 ) =  327
greenupxy  (  94, 5 ) =  327
greenupxy  (  94, 6 ) =  327
greenupxy  (  94, 7 ) =  314
greenupxy  (  94, 8 ) =  314
fallxy     (  94, 1 ) =  125
fallxy     (  94, 2 ) =   78
fallxy     (  94, 3 ) =  106
fallxy     (  94, 4 ) =  108
fallxy     (  94, 5 ) =  108
fallxy     (  94, 6 ) =  108
fallxy     (  94, 7 ) =  105
fallxy     (  94, 8 ) =  105
xyphendoy1 (  94, 1 ) =    2
xyphendoy1 (  94, 2 ) =    2
xyphendoy1 (  94, 3 ) =    2
xyphendoy1 (  94, 4 ) =    2
xyphendoy1 (  94, 5 ) =    2
xyphendoy1 (  94, 6 ) =    2
xyphendoy1 (  94, 7 ) =    2
xyphendoy1 (  94, 8 ) =    2
      
 ! xlat        -9.25000000000000     
greenupxy  (  93, 1 ) =  343
greenupxy  (  93, 2 ) =  254
greenupxy  (  93, 3 ) =  335
greenupxy  (  93, 4 ) =  324
greenupxy  (  93, 5 ) =  324
greenupxy  (  93, 6 ) =  324
greenupxy  (  93, 7 ) =  313
greenupxy  (  93, 8 ) =  313
fallxy     (  93, 1 ) =  125
fallxy     (  93, 2 ) =   76
fallxy     (  93, 3 ) =  103
fallxy     (  93, 4 ) =  106
fallxy     (  93, 5 ) =  106
fallxy     (  93, 6 ) =  106
fallxy     (  93, 7 ) =  103
fallxy     (  93, 8 ) =  103
xyphendoy1 (  93, 1 ) =    2
xyphendoy1 (  93, 2 ) =    2
xyphendoy1 (  93, 3 ) =    2
xyphendoy1 (  93, 4 ) =    2
xyphendoy1 (  93, 5 ) =    2
xyphendoy1 (  93, 6 ) =    2
xyphendoy1 (  93, 7 ) =    2
xyphendoy1 (  93, 8 ) =    2
      
 ! xlat        -9.75000000000000     
greenupxy  (  92, 1 ) =  343
greenupxy  (  92, 2 ) =  255
greenupxy  (  92, 3 ) =  331
greenupxy  (  92, 4 ) =  321
greenupxy  (  92, 5 ) =  321
greenupxy  (  92, 6 ) =  321
greenupxy  (  92, 7 ) =  311
greenupxy  (  92, 8 ) =  311
fallxy     (  92, 1 ) =  125
fallxy     (  92, 2 ) =   74
fallxy     (  92, 3 ) =  101
fallxy     (  92, 4 ) =  104
fallxy     (  92, 5 ) =  104
fallxy     (  92, 6 ) =  104
fallxy     (  92, 7 ) =   98
fallxy     (  92, 8 ) =   98
xyphendoy1 (  92, 1 ) =    2
xyphendoy1 (  92, 2 ) =    2
xyphendoy1 (  92, 3 ) =    2
xyphendoy1 (  92, 4 ) =    2
xyphendoy1 (  92, 5 ) =    2
xyphendoy1 (  92, 6 ) =    2
xyphendoy1 (  92, 7 ) =    2
xyphendoy1 (  92, 8 ) =    2
      
 ! xlat        -10.2500000000000     
greenupxy  (  91, 1 ) =  337
greenupxy  (  91, 2 ) =  255
greenupxy  (  91, 3 ) =  327
greenupxy  (  91, 4 ) =  318
greenupxy  (  91, 5 ) =  318
greenupxy  (  91, 6 ) =  318
greenupxy  (  91, 7 ) =  311
greenupxy  (  91, 8 ) =  311
fallxy     (  91, 1 ) =  118
fallxy     (  91, 2 ) =   72
fallxy     (  91, 3 ) =   99
fallxy     (  91, 4 ) =  102
fallxy     (  91, 5 ) =  102
fallxy     (  91, 6 ) =  102
fallxy     (  91, 7 ) =   95
fallxy     (  91, 8 ) =   95
xyphendoy1 (  91, 1 ) =    2
xyphendoy1 (  91, 2 ) =    2
xyphendoy1 (  91, 3 ) =    2
xyphendoy1 (  91, 4 ) =    2
xyphendoy1 (  91, 5 ) =    2
xyphendoy1 (  91, 6 ) =    2
xyphendoy1 (  91, 7 ) =    2
xyphendoy1 (  91, 8 ) =    2
      
 ! xlat        -10.7500000000000     
greenupxy  (  90, 1 ) =  333
greenupxy  (  90, 2 ) =  258
greenupxy  (  90, 3 ) =  323
greenupxy  (  90, 4 ) =  315
greenupxy  (  90, 5 ) =  315
greenupxy  (  90, 6 ) =  315
greenupxy  (  90, 7 ) =  311
greenupxy  (  90, 8 ) =  311
fallxy     (  90, 1 ) =  112
fallxy     (  90, 2 ) =   71
fallxy     (  90, 3 ) =   97
fallxy     (  90, 4 ) =  101
fallxy     (  90, 5 ) =  101
fallxy     (  90, 6 ) =  101
fallxy     (  90, 7 ) =   92
fallxy     (  90, 8 ) =   92
xyphendoy1 (  90, 1 ) =    2
xyphendoy1 (  90, 2 ) =    2
xyphendoy1 (  90, 3 ) =    2
xyphendoy1 (  90, 4 ) =    2
xyphendoy1 (  90, 5 ) =    2
xyphendoy1 (  90, 6 ) =    2
xyphendoy1 (  90, 7 ) =    2
xyphendoy1 (  90, 8 ) =    2
      
 ! xlat        -11.2500000000000     
greenupxy  (  89, 1 ) =  333
greenupxy  (  89, 2 ) =  262
greenupxy  (  89, 3 ) =  321
greenupxy  (  89, 4 ) =  313
greenupxy  (  89, 5 ) =  313
greenupxy  (  89, 6 ) =  313
greenupxy  (  89, 7 ) =  310
greenupxy  (  89, 8 ) =  310
fallxy     (  89, 1 ) =  112
fallxy     (  89, 2 ) =   71
fallxy     (  89, 3 ) =   96
fallxy     (  89, 4 ) =   99
fallxy     (  89, 5 ) =   99
fallxy     (  89, 6 ) =   99
fallxy     (  89, 7 ) =   90
fallxy     (  89, 8 ) =   90
xyphendoy1 (  89, 1 ) =    2
xyphendoy1 (  89, 2 ) =    2
xyphendoy1 (  89, 3 ) =    2
xyphendoy1 (  89, 4 ) =    2
xyphendoy1 (  89, 5 ) =    2
xyphendoy1 (  89, 6 ) =    2
xyphendoy1 (  89, 7 ) =    2
xyphendoy1 (  89, 8 ) =    2
      
 ! xlat        -11.7500000000000     
greenupxy  (  88, 1 ) =  328
greenupxy  (  88, 2 ) =  266
greenupxy  (  88, 3 ) =  319
greenupxy  (  88, 4 ) =  311
greenupxy  (  88, 5 ) =  311
greenupxy  (  88, 6 ) =  311
greenupxy  (  88, 7 ) =  309
greenupxy  (  88, 8 ) =  309
fallxy     (  88, 1 ) =  104
fallxy     (  88, 2 ) =   72
fallxy     (  88, 3 ) =   95
fallxy     (  88, 4 ) =   98
fallxy     (  88, 5 ) =   98
fallxy     (  88, 6 ) =   98
fallxy     (  88, 7 ) =   89
fallxy     (  88, 8 ) =   89
xyphendoy1 (  88, 1 ) =    2
xyphendoy1 (  88, 2 ) =    2
xyphendoy1 (  88, 3 ) =    2
xyphendoy1 (  88, 4 ) =    2
xyphendoy1 (  88, 5 ) =    2
xyphendoy1 (  88, 6 ) =    2
xyphendoy1 (  88, 7 ) =    2
xyphendoy1 (  88, 8 ) =    2
      
 ! xlat        -12.2500000000000     
greenupxy  (  87, 1 ) =  328
greenupxy  (  87, 2 ) =  269
greenupxy  (  87, 3 ) =  317
greenupxy  (  87, 4 ) =  310
greenupxy  (  87, 5 ) =  310
greenupxy  (  87, 6 ) =  310
greenupxy  (  87, 7 ) =  308
greenupxy  (  87, 8 ) =  308
fallxy     (  87, 1 ) =  104
fallxy     (  87, 2 ) =   72
fallxy     (  87, 3 ) =   94
fallxy     (  87, 4 ) =   97
fallxy     (  87, 5 ) =   97
fallxy     (  87, 6 ) =   97
fallxy     (  87, 7 ) =   87
fallxy     (  87, 8 ) =   87
xyphendoy1 (  87, 1 ) =    2
xyphendoy1 (  87, 2 ) =    2
xyphendoy1 (  87, 3 ) =    2
xyphendoy1 (  87, 4 ) =    2
xyphendoy1 (  87, 5 ) =    2
xyphendoy1 (  87, 6 ) =    2
xyphendoy1 (  87, 7 ) =    2
xyphendoy1 (  87, 8 ) =    2
      
 ! xlat        -12.7500000000000     
greenupxy  (  86, 1 ) =  328
greenupxy  (  86, 2 ) =  272
greenupxy  (  86, 3 ) =  316
greenupxy  (  86, 4 ) =  309
greenupxy  (  86, 5 ) =  309
greenupxy  (  86, 6 ) =  309
greenupxy  (  86, 7 ) =  307
greenupxy  (  86, 8 ) =  307
fallxy     (  86, 1 ) =  104
fallxy     (  86, 2 ) =   72
fallxy     (  86, 3 ) =   94
fallxy     (  86, 4 ) =   96
fallxy     (  86, 5 ) =   96
fallxy     (  86, 6 ) =   96
fallxy     (  86, 7 ) =   85
fallxy     (  86, 8 ) =   85
xyphendoy1 (  86, 1 ) =    2
xyphendoy1 (  86, 2 ) =    2
xyphendoy1 (  86, 3 ) =    2
xyphendoy1 (  86, 4 ) =    2
xyphendoy1 (  86, 5 ) =    2
xyphendoy1 (  86, 6 ) =    2
xyphendoy1 (  86, 7 ) =    2
xyphendoy1 (  86, 8 ) =    2
      
 ! xlat        -13.2500000000000     
greenupxy  (  85, 1 ) =  326
greenupxy  (  85, 2 ) =  275
greenupxy  (  85, 3 ) =  316
greenupxy  (  85, 4 ) =  309
greenupxy  (  85, 5 ) =  309
greenupxy  (  85, 6 ) =  309
greenupxy  (  85, 7 ) =  306
greenupxy  (  85, 8 ) =  306
fallxy     (  85, 1 ) =  105
fallxy     (  85, 2 ) =   72
fallxy     (  85, 3 ) =   94
fallxy     (  85, 4 ) =   95
fallxy     (  85, 5 ) =   95
fallxy     (  85, 6 ) =   95
fallxy     (  85, 7 ) =   83
fallxy     (  85, 8 ) =   83
xyphendoy1 (  85, 1 ) =    2
xyphendoy1 (  85, 2 ) =    2
xyphendoy1 (  85, 3 ) =    2
xyphendoy1 (  85, 4 ) =    2
xyphendoy1 (  85, 5 ) =    2
xyphendoy1 (  85, 6 ) =    2
xyphendoy1 (  85, 7 ) =    2
xyphendoy1 (  85, 8 ) =    2
      
 ! xlat        -13.7500000000000     
greenupxy  (  84, 1 ) =  321
greenupxy  (  84, 2 ) =  277
greenupxy  (  84, 3 ) =  317
greenupxy  (  84, 4 ) =  311
greenupxy  (  84, 5 ) =  311
greenupxy  (  84, 6 ) =  311
greenupxy  (  84, 7 ) =  305
greenupxy  (  84, 8 ) =  305
fallxy     (  84, 1 ) =   97
fallxy     (  84, 2 ) =   73
fallxy     (  84, 3 ) =   93
fallxy     (  84, 4 ) =   94
fallxy     (  84, 5 ) =   94
fallxy     (  84, 6 ) =   94
fallxy     (  84, 7 ) =   80
fallxy     (  84, 8 ) =   80
xyphendoy1 (  84, 1 ) =    2
xyphendoy1 (  84, 2 ) =    2
xyphendoy1 (  84, 3 ) =    2
xyphendoy1 (  84, 4 ) =    2
xyphendoy1 (  84, 5 ) =    2
xyphendoy1 (  84, 6 ) =    2
xyphendoy1 (  84, 7 ) =    2
xyphendoy1 (  84, 8 ) =    2
      
 ! xlat        -14.2500000000000     
greenupxy  (  83, 1 ) =  320
greenupxy  (  83, 2 ) =  279
greenupxy  (  83, 3 ) =  318
greenupxy  (  83, 4 ) =  313
greenupxy  (  83, 5 ) =  313
greenupxy  (  83, 6 ) =  313
greenupxy  (  83, 7 ) =  304
greenupxy  (  83, 8 ) =  304
fallxy     (  83, 1 ) =   93
fallxy     (  83, 2 ) =   74
fallxy     (  83, 3 ) =   93
fallxy     (  83, 4 ) =   93
fallxy     (  83, 5 ) =   93
fallxy     (  83, 6 ) =   93
fallxy     (  83, 7 ) =   79
fallxy     (  83, 8 ) =   79
xyphendoy1 (  83, 1 ) =    2
xyphendoy1 (  83, 2 ) =    2
xyphendoy1 (  83, 3 ) =    2
xyphendoy1 (  83, 4 ) =    2
xyphendoy1 (  83, 5 ) =    2
xyphendoy1 (  83, 6 ) =    2
xyphendoy1 (  83, 7 ) =    2
xyphendoy1 (  83, 8 ) =    2
      
 ! xlat        -14.7500000000000     
greenupxy  (  82, 1 ) =  320
greenupxy  (  82, 2 ) =  282
greenupxy  (  82, 3 ) =  321
greenupxy  (  82, 4 ) =  316
greenupxy  (  82, 5 ) =  316
greenupxy  (  82, 6 ) =  316
greenupxy  (  82, 7 ) =  303
greenupxy  (  82, 8 ) =  303
fallxy     (  82, 1 ) =   93
fallxy     (  82, 2 ) =   75
fallxy     (  82, 3 ) =   92
fallxy     (  82, 4 ) =   93
fallxy     (  82, 5 ) =   93
fallxy     (  82, 6 ) =   93
fallxy     (  82, 7 ) =   78
fallxy     (  82, 8 ) =   78
xyphendoy1 (  82, 1 ) =    2
xyphendoy1 (  82, 2 ) =    2
xyphendoy1 (  82, 3 ) =    2
xyphendoy1 (  82, 4 ) =    2
xyphendoy1 (  82, 5 ) =    2
xyphendoy1 (  82, 6 ) =    2
xyphendoy1 (  82, 7 ) =    2
xyphendoy1 (  82, 8 ) =    2
      
 ! xlat        -15.2500000000000     
greenupxy  (  81, 1 ) =  320
greenupxy  (  81, 2 ) =  284
greenupxy  (  81, 3 ) =  324
greenupxy  (  81, 4 ) =  318
greenupxy  (  81, 5 ) =  318
greenupxy  (  81, 6 ) =  318
greenupxy  (  81, 7 ) =  302
greenupxy  (  81, 8 ) =  302
fallxy     (  81, 1 ) =   93
fallxy     (  81, 2 ) =   75
fallxy     (  81, 3 ) =   92
fallxy     (  81, 4 ) =   92
fallxy     (  81, 5 ) =   92
fallxy     (  81, 6 ) =   92
fallxy     (  81, 7 ) =   77
fallxy     (  81, 8 ) =   77
xyphendoy1 (  81, 1 ) =    2
xyphendoy1 (  81, 2 ) =    2
xyphendoy1 (  81, 3 ) =    2
xyphendoy1 (  81, 4 ) =    2
xyphendoy1 (  81, 5 ) =    2
xyphendoy1 (  81, 6 ) =    2
xyphendoy1 (  81, 7 ) =    2
xyphendoy1 (  81, 8 ) =    2
      
 ! xlat        -15.7500000000000     
greenupxy  (  80, 1 ) =  320
greenupxy  (  80, 2 ) =  287
greenupxy  (  80, 3 ) =  328
greenupxy  (  80, 4 ) =  321
greenupxy  (  80, 5 ) =  321
greenupxy  (  80, 6 ) =  321
greenupxy  (  80, 7 ) =  301
greenupxy  (  80, 8 ) =  301
fallxy     (  80, 1 ) =   93
fallxy     (  80, 2 ) =   76
fallxy     (  80, 3 ) =   92
fallxy     (  80, 4 ) =   91
fallxy     (  80, 5 ) =   91
fallxy     (  80, 6 ) =   91
fallxy     (  80, 7 ) =   76
fallxy     (  80, 8 ) =   76
xyphendoy1 (  80, 1 ) =    2
xyphendoy1 (  80, 2 ) =    2
xyphendoy1 (  80, 3 ) =    2
xyphendoy1 (  80, 4 ) =    2
xyphendoy1 (  80, 5 ) =    2
xyphendoy1 (  80, 6 ) =    2
xyphendoy1 (  80, 7 ) =    2
xyphendoy1 (  80, 8 ) =    2
      
 ! xlat        -16.2500000000000     
greenupxy  (  79, 1 ) =  322
greenupxy  (  79, 2 ) =  289
greenupxy  (  79, 3 ) =  332
greenupxy  (  79, 4 ) =  325
greenupxy  (  79, 5 ) =  325
greenupxy  (  79, 6 ) =  325
greenupxy  (  79, 7 ) =  300
greenupxy  (  79, 8 ) =  300
fallxy     (  79, 1 ) =   94
fallxy     (  79, 2 ) =   77
fallxy     (  79, 3 ) =   92
fallxy     (  79, 4 ) =   90
fallxy     (  79, 5 ) =   90
fallxy     (  79, 6 ) =   90
fallxy     (  79, 7 ) =   77
fallxy     (  79, 8 ) =   77
xyphendoy1 (  79, 1 ) =    2
xyphendoy1 (  79, 2 ) =    2
xyphendoy1 (  79, 3 ) =    2
xyphendoy1 (  79, 4 ) =    2
xyphendoy1 (  79, 5 ) =    2
xyphendoy1 (  79, 6 ) =    2
xyphendoy1 (  79, 7 ) =    2
xyphendoy1 (  79, 8 ) =    2
      
 ! xlat        -16.7500000000000     
greenupxy  (  78, 1 ) =  325
greenupxy  (  78, 2 ) =  292
greenupxy  (  78, 3 ) =  337
greenupxy  (  78, 4 ) =  329
greenupxy  (  78, 5 ) =  329
greenupxy  (  78, 6 ) =  329
greenupxy  (  78, 7 ) =  298
greenupxy  (  78, 8 ) =  298
fallxy     (  78, 1 ) =   96
fallxy     (  78, 2 ) =   78
fallxy     (  78, 3 ) =   91
fallxy     (  78, 4 ) =   90
fallxy     (  78, 5 ) =   90
fallxy     (  78, 6 ) =   90
fallxy     (  78, 7 ) =   77
fallxy     (  78, 8 ) =   77
xyphendoy1 (  78, 1 ) =    2
xyphendoy1 (  78, 2 ) =    2
xyphendoy1 (  78, 3 ) =    2
xyphendoy1 (  78, 4 ) =    2
xyphendoy1 (  78, 5 ) =    2
xyphendoy1 (  78, 6 ) =    2
xyphendoy1 (  78, 7 ) =    2
xyphendoy1 (  78, 8 ) =    2
      
 ! xlat        -17.2500000000000     
greenupxy  (  77, 1 ) =  328
greenupxy  (  77, 2 ) =  294
greenupxy  (  77, 3 ) =  343
greenupxy  (  77, 4 ) =  331
greenupxy  (  77, 5 ) =  331
greenupxy  (  77, 6 ) =  331
greenupxy  (  77, 7 ) =  297
greenupxy  (  77, 8 ) =  297
fallxy     (  77, 1 ) =   93
fallxy     (  77, 2 ) =   79
fallxy     (  77, 3 ) =   91
fallxy     (  77, 4 ) =   90
fallxy     (  77, 5 ) =   90
fallxy     (  77, 6 ) =   90
fallxy     (  77, 7 ) =   78
fallxy     (  77, 8 ) =   78
xyphendoy1 (  77, 1 ) =    2
xyphendoy1 (  77, 2 ) =    2
xyphendoy1 (  77, 3 ) =    2
xyphendoy1 (  77, 4 ) =    2
xyphendoy1 (  77, 5 ) =    2
xyphendoy1 (  77, 6 ) =    2
xyphendoy1 (  77, 7 ) =    2
xyphendoy1 (  77, 8 ) =    2
      
 ! xlat        -17.7500000000000     
greenupxy  (  76, 1 ) =  329
greenupxy  (  76, 2 ) =  296
greenupxy  (  76, 3 ) =  348
greenupxy  (  76, 4 ) =  333
greenupxy  (  76, 5 ) =  333
greenupxy  (  76, 6 ) =  333
greenupxy  (  76, 7 ) =  296
greenupxy  (  76, 8 ) =  296
fallxy     (  76, 1 ) =   97
fallxy     (  76, 2 ) =   80
fallxy     (  76, 3 ) =   91
fallxy     (  76, 4 ) =   90
fallxy     (  76, 5 ) =   90
fallxy     (  76, 6 ) =   90
fallxy     (  76, 7 ) =   79
fallxy     (  76, 8 ) =   79
xyphendoy1 (  76, 1 ) =    2
xyphendoy1 (  76, 2 ) =    2
xyphendoy1 (  76, 3 ) =    2
xyphendoy1 (  76, 4 ) =    2
xyphendoy1 (  76, 5 ) =    2
xyphendoy1 (  76, 6 ) =    2
xyphendoy1 (  76, 7 ) =    2
xyphendoy1 (  76, 8 ) =    2
      
 ! xlat        -18.2500000000000     
greenupxy  (  75, 1 ) =  326
greenupxy  (  75, 2 ) =  298
greenupxy  (  75, 3 ) =  352
greenupxy  (  75, 4 ) =  333
greenupxy  (  75, 5 ) =  333
greenupxy  (  75, 6 ) =  333
greenupxy  (  75, 7 ) =  295
greenupxy  (  75, 8 ) =  295
fallxy     (  75, 1 ) =   94
fallxy     (  75, 2 ) =   81
fallxy     (  75, 3 ) =   91
fallxy     (  75, 4 ) =   90
fallxy     (  75, 5 ) =   90
fallxy     (  75, 6 ) =   90
fallxy     (  75, 7 ) =   79
fallxy     (  75, 8 ) =   79
xyphendoy1 (  75, 1 ) =    2
xyphendoy1 (  75, 2 ) =    2
xyphendoy1 (  75, 3 ) =    2
xyphendoy1 (  75, 4 ) =    2
xyphendoy1 (  75, 5 ) =    2
xyphendoy1 (  75, 6 ) =    2
xyphendoy1 (  75, 7 ) =    2
xyphendoy1 (  75, 8 ) =    2
      
 ! xlat        -18.7500000000000     
greenupxy  (  74, 1 ) =  327
greenupxy  (  74, 2 ) =  299
greenupxy  (  74, 3 ) =  357
greenupxy  (  74, 4 ) =  333
greenupxy  (  74, 5 ) =  333
greenupxy  (  74, 6 ) =  333
greenupxy  (  74, 7 ) =  294
greenupxy  (  74, 8 ) =  294
fallxy     (  74, 1 ) =   92
fallxy     (  74, 2 ) =   82
fallxy     (  74, 3 ) =   93
fallxy     (  74, 4 ) =   90
fallxy     (  74, 5 ) =   90
fallxy     (  74, 6 ) =   90
fallxy     (  74, 7 ) =   80
fallxy     (  74, 8 ) =   80
xyphendoy1 (  74, 1 ) =    2
xyphendoy1 (  74, 2 ) =    2
xyphendoy1 (  74, 3 ) =    2
xyphendoy1 (  74, 4 ) =    2
xyphendoy1 (  74, 5 ) =    2
xyphendoy1 (  74, 6 ) =    2
xyphendoy1 (  74, 7 ) =    2
xyphendoy1 (  74, 8 ) =    2
      
 ! xlat        -19.2500000000000     
greenupxy  (  73, 1 ) =  330
greenupxy  (  73, 2 ) =  299
greenupxy  (  73, 3 ) =  361
greenupxy  (  73, 4 ) =  332
greenupxy  (  73, 5 ) =  332
greenupxy  (  73, 6 ) =  332
greenupxy  (  73, 7 ) =  294
greenupxy  (  73, 8 ) =  294
fallxy     (  73, 1 ) =   94
fallxy     (  73, 2 ) =   82
fallxy     (  73, 3 ) =   93
fallxy     (  73, 4 ) =   90
fallxy     (  73, 5 ) =   90
fallxy     (  73, 6 ) =   90
fallxy     (  73, 7 ) =   80
fallxy     (  73, 8 ) =   80
xyphendoy1 (  73, 1 ) =    2
xyphendoy1 (  73, 2 ) =    2
xyphendoy1 (  73, 3 ) =    2
xyphendoy1 (  73, 4 ) =    2
xyphendoy1 (  73, 5 ) =    2
xyphendoy1 (  73, 6 ) =    2
xyphendoy1 (  73, 7 ) =    2
xyphendoy1 (  73, 8 ) =    2
      
 ! xlat        -19.7500000000000     
greenupxy  (  72, 1 ) =  330
greenupxy  (  72, 2 ) =  298
greenupxy  (  72, 3 ) =  364
greenupxy  (  72, 4 ) =  331
greenupxy  (  72, 5 ) =  331
greenupxy  (  72, 6 ) =  331
greenupxy  (  72, 7 ) =  294
greenupxy  (  72, 8 ) =  294
fallxy     (  72, 1 ) =   99
fallxy     (  72, 2 ) =   82
fallxy     (  72, 3 ) =   94
fallxy     (  72, 4 ) =   90
fallxy     (  72, 5 ) =   90
fallxy     (  72, 6 ) =   90
fallxy     (  72, 7 ) =   80
fallxy     (  72, 8 ) =   80
xyphendoy1 (  72, 1 ) =    2
xyphendoy1 (  72, 2 ) =    2
xyphendoy1 (  72, 3 ) =    2
xyphendoy1 (  72, 4 ) =    2
xyphendoy1 (  72, 5 ) =    2
xyphendoy1 (  72, 6 ) =    2
xyphendoy1 (  72, 7 ) =    2
xyphendoy1 (  72, 8 ) =    2
      
 ! xlat        -20.2500000000000     
greenupxy  (  71, 1 ) =  328
greenupxy  (  71, 2 ) =  297
greenupxy  (  71, 3 ) =    2
greenupxy  (  71, 4 ) =  330
greenupxy  (  71, 5 ) =  330
greenupxy  (  71, 6 ) =  330
greenupxy  (  71, 7 ) =  293
greenupxy  (  71, 8 ) =  293
fallxy     (  71, 1 ) =  101
fallxy     (  71, 2 ) =   82
fallxy     (  71, 3 ) =   95
fallxy     (  71, 4 ) =   90
fallxy     (  71, 5 ) =   90
fallxy     (  71, 6 ) =   90
fallxy     (  71, 7 ) =   80
fallxy     (  71, 8 ) =   80
xyphendoy1 (  71, 1 ) =    2
xyphendoy1 (  71, 2 ) =    2
xyphendoy1 (  71, 3 ) =    0
xyphendoy1 (  71, 4 ) =    2
xyphendoy1 (  71, 5 ) =    2
xyphendoy1 (  71, 6 ) =    2
xyphendoy1 (  71, 7 ) =    2
xyphendoy1 (  71, 8 ) =    2
      
 ! xlat        -20.7500000000000     
greenupxy  (  70, 1 ) =  328
greenupxy  (  70, 2 ) =  296
greenupxy  (  70, 3 ) =    5
greenupxy  (  70, 4 ) =  330
greenupxy  (  70, 5 ) =  330
greenupxy  (  70, 6 ) =  330
greenupxy  (  70, 7 ) =  292
greenupxy  (  70, 8 ) =  292
fallxy     (  70, 1 ) =  101
fallxy     (  70, 2 ) =   82
fallxy     (  70, 3 ) =   96
fallxy     (  70, 4 ) =   90
fallxy     (  70, 5 ) =   90
fallxy     (  70, 6 ) =   90
fallxy     (  70, 7 ) =   80
fallxy     (  70, 8 ) =   80
xyphendoy1 (  70, 1 ) =    2
xyphendoy1 (  70, 2 ) =    2
xyphendoy1 (  70, 3 ) =    0
xyphendoy1 (  70, 4 ) =    2
xyphendoy1 (  70, 5 ) =    2
xyphendoy1 (  70, 6 ) =    2
xyphendoy1 (  70, 7 ) =    2
xyphendoy1 (  70, 8 ) =    2
      
 ! xlat        -21.2500000000000     
greenupxy  (  69, 1 ) =  321
greenupxy  (  69, 2 ) =  296
greenupxy  (  69, 3 ) =    8
greenupxy  (  69, 4 ) =  328
greenupxy  (  69, 5 ) =  328
greenupxy  (  69, 6 ) =  328
greenupxy  (  69, 7 ) =  293
greenupxy  (  69, 8 ) =  293
fallxy     (  69, 1 ) =   99
fallxy     (  69, 2 ) =   82
fallxy     (  69, 3 ) =   98
fallxy     (  69, 4 ) =   89
fallxy     (  69, 5 ) =   89
fallxy     (  69, 6 ) =   89
fallxy     (  69, 7 ) =   81
fallxy     (  69, 8 ) =   81
xyphendoy1 (  69, 1 ) =    2
xyphendoy1 (  69, 2 ) =    2
xyphendoy1 (  69, 3 ) =    0
xyphendoy1 (  69, 4 ) =    2
xyphendoy1 (  69, 5 ) =    2
xyphendoy1 (  69, 6 ) =    2
xyphendoy1 (  69, 7 ) =    2
xyphendoy1 (  69, 8 ) =    2
      
 ! xlat        -21.7500000000000     
greenupxy  (  68, 1 ) =  321
greenupxy  (  68, 2 ) =  295
greenupxy  (  68, 3 ) =   11
greenupxy  (  68, 4 ) =  325
greenupxy  (  68, 5 ) =  325
greenupxy  (  68, 6 ) =  325
greenupxy  (  68, 7 ) =  292
greenupxy  (  68, 8 ) =  292
fallxy     (  68, 1 ) =   96
fallxy     (  68, 2 ) =   82
fallxy     (  68, 3 ) =  101
fallxy     (  68, 4 ) =   89
fallxy     (  68, 5 ) =   89
fallxy     (  68, 6 ) =   89
fallxy     (  68, 7 ) =   82
fallxy     (  68, 8 ) =   82
xyphendoy1 (  68, 1 ) =    2
xyphendoy1 (  68, 2 ) =    2
xyphendoy1 (  68, 3 ) =    0
xyphendoy1 (  68, 4 ) =    2
xyphendoy1 (  68, 5 ) =    2
xyphendoy1 (  68, 6 ) =    2
xyphendoy1 (  68, 7 ) =    2
xyphendoy1 (  68, 8 ) =    2
      
 ! xlat        -22.2500000000000     
greenupxy  (  67, 1 ) =  326
greenupxy  (  67, 2 ) =  293
greenupxy  (  67, 3 ) =   14
greenupxy  (  67, 4 ) =  322
greenupxy  (  67, 5 ) =  322
greenupxy  (  67, 6 ) =  322
greenupxy  (  67, 7 ) =  292
greenupxy  (  67, 8 ) =  292
fallxy     (  67, 1 ) =  106
fallxy     (  67, 2 ) =   81
fallxy     (  67, 3 ) =  104
fallxy     (  67, 4 ) =   88
fallxy     (  67, 5 ) =   88
fallxy     (  67, 6 ) =   88
fallxy     (  67, 7 ) =   82
fallxy     (  67, 8 ) =   82
xyphendoy1 (  67, 1 ) =    2
xyphendoy1 (  67, 2 ) =    2
xyphendoy1 (  67, 3 ) =    0
xyphendoy1 (  67, 4 ) =    2
xyphendoy1 (  67, 5 ) =    2
xyphendoy1 (  67, 6 ) =    2
xyphendoy1 (  67, 7 ) =    2
xyphendoy1 (  67, 8 ) =    2
      
 ! xlat        -22.7500000000000     
greenupxy  (  66, 1 ) =  326
greenupxy  (  66, 2 ) =  292
greenupxy  (  66, 3 ) =   17
greenupxy  (  66, 4 ) =  319
greenupxy  (  66, 5 ) =  319
greenupxy  (  66, 6 ) =  319
greenupxy  (  66, 7 ) =  292
greenupxy  (  66, 8 ) =  292
fallxy     (  66, 1 ) =  106
fallxy     (  66, 2 ) =   81
fallxy     (  66, 3 ) =  107
fallxy     (  66, 4 ) =   88
fallxy     (  66, 5 ) =   88
fallxy     (  66, 6 ) =   88
fallxy     (  66, 7 ) =   82
fallxy     (  66, 8 ) =   82
xyphendoy1 (  66, 1 ) =    2
xyphendoy1 (  66, 2 ) =    2
xyphendoy1 (  66, 3 ) =    0
xyphendoy1 (  66, 4 ) =    2
xyphendoy1 (  66, 5 ) =    2
xyphendoy1 (  66, 6 ) =    2
xyphendoy1 (  66, 7 ) =    2
xyphendoy1 (  66, 8 ) =    2
      
 ! xlat        -23.2500000000000     
greenupxy  (  65, 1 ) =  321
greenupxy  (  65, 2 ) =  290
greenupxy  (  65, 3 ) =   20
greenupxy  (  65, 4 ) =  316
greenupxy  (  65, 5 ) =  316
greenupxy  (  65, 6 ) =  316
greenupxy  (  65, 7 ) =  293
greenupxy  (  65, 8 ) =  293
fallxy     (  65, 1 ) =  107
fallxy     (  65, 2 ) =   81
fallxy     (  65, 3 ) =  111
fallxy     (  65, 4 ) =   88
fallxy     (  65, 5 ) =   88
fallxy     (  65, 6 ) =   88
fallxy     (  65, 7 ) =   82
fallxy     (  65, 8 ) =   82
xyphendoy1 (  65, 1 ) =    2
xyphendoy1 (  65, 2 ) =    2
xyphendoy1 (  65, 3 ) =    0
xyphendoy1 (  65, 4 ) =    2
xyphendoy1 (  65, 5 ) =    2
xyphendoy1 (  65, 6 ) =    2
xyphendoy1 (  65, 7 ) =    2
xyphendoy1 (  65, 8 ) =    2
      
 ! xlat        -23.7500000000000     
greenupxy  (  64, 1 ) =  320
greenupxy  (  64, 2 ) =  289
greenupxy  (  64, 3 ) =   23
greenupxy  (  64, 4 ) =  313
greenupxy  (  64, 5 ) =  313
greenupxy  (  64, 6 ) =  313
greenupxy  (  64, 7 ) =  291
greenupxy  (  64, 8 ) =  291
fallxy     (  64, 1 ) =  108
fallxy     (  64, 2 ) =   81
fallxy     (  64, 3 ) =  114
fallxy     (  64, 4 ) =   89
fallxy     (  64, 5 ) =   89
fallxy     (  64, 6 ) =   89
fallxy     (  64, 7 ) =   82
fallxy     (  64, 8 ) =   82
xyphendoy1 (  64, 1 ) =    2
xyphendoy1 (  64, 2 ) =    2
xyphendoy1 (  64, 3 ) =    0
xyphendoy1 (  64, 4 ) =    2
xyphendoy1 (  64, 5 ) =    2
xyphendoy1 (  64, 6 ) =    2
xyphendoy1 (  64, 7 ) =    2
xyphendoy1 (  64, 8 ) =    2
      
 ! xlat        -24.2500000000000     
greenupxy  (  63, 1 ) =  321
greenupxy  (  63, 2 ) =  288
greenupxy  (  63, 3 ) =   26
greenupxy  (  63, 4 ) =  310
greenupxy  (  63, 5 ) =  310
greenupxy  (  63, 6 ) =  310
greenupxy  (  63, 7 ) =  290
greenupxy  (  63, 8 ) =  290
fallxy     (  63, 1 ) =  109
fallxy     (  63, 2 ) =   81
fallxy     (  63, 3 ) =  119
fallxy     (  63, 4 ) =   90
fallxy     (  63, 5 ) =   90
fallxy     (  63, 6 ) =   90
fallxy     (  63, 7 ) =   82
fallxy     (  63, 8 ) =   82
xyphendoy1 (  63, 1 ) =    2
xyphendoy1 (  63, 2 ) =    2
xyphendoy1 (  63, 3 ) =    0
xyphendoy1 (  63, 4 ) =    2
xyphendoy1 (  63, 5 ) =    2
xyphendoy1 (  63, 6 ) =    2
xyphendoy1 (  63, 7 ) =    2
xyphendoy1 (  63, 8 ) =    2
      
 ! xlat        -24.7500000000000     
greenupxy  (  62, 1 ) =  313
greenupxy  (  62, 2 ) =  286
greenupxy  (  62, 3 ) =   29
greenupxy  (  62, 4 ) =  307
greenupxy  (  62, 5 ) =  307
greenupxy  (  62, 6 ) =  307
greenupxy  (  62, 7 ) =  289
greenupxy  (  62, 8 ) =  289
fallxy     (  62, 1 ) =  104
fallxy     (  62, 2 ) =   80
fallxy     (  62, 3 ) =  125
fallxy     (  62, 4 ) =   90
fallxy     (  62, 5 ) =   90
fallxy     (  62, 6 ) =   90
fallxy     (  62, 7 ) =   82
fallxy     (  62, 8 ) =   82
xyphendoy1 (  62, 1 ) =    2
xyphendoy1 (  62, 2 ) =    2
xyphendoy1 (  62, 3 ) =    0
xyphendoy1 (  62, 4 ) =    2
xyphendoy1 (  62, 5 ) =    2
xyphendoy1 (  62, 6 ) =    2
xyphendoy1 (  62, 7 ) =    2
xyphendoy1 (  62, 8 ) =    2
      
 ! xlat        -25.2500000000000     
greenupxy  (  61, 1 ) =  305
greenupxy  (  61, 2 ) =  285
greenupxy  (  61, 3 ) =   33
greenupxy  (  61, 4 ) =  304
greenupxy  (  61, 5 ) =  304
greenupxy  (  61, 6 ) =  304
greenupxy  (  61, 7 ) =  281
greenupxy  (  61, 8 ) =  281
fallxy     (  61, 1 ) =  106
fallxy     (  61, 2 ) =   80
fallxy     (  61, 3 ) =  131
fallxy     (  61, 4 ) =   91
fallxy     (  61, 5 ) =   91
fallxy     (  61, 6 ) =   91
fallxy     (  61, 7 ) =   82
fallxy     (  61, 8 ) =   82
xyphendoy1 (  61, 1 ) =    2
xyphendoy1 (  61, 2 ) =    2
xyphendoy1 (  61, 3 ) =    0
xyphendoy1 (  61, 4 ) =    2
xyphendoy1 (  61, 5 ) =    2
xyphendoy1 (  61, 6 ) =    2
xyphendoy1 (  61, 7 ) =    2
xyphendoy1 (  61, 8 ) =    2
      
 ! xlat        -25.7500000000000     
greenupxy  (  60, 1 ) =  291
greenupxy  (  60, 2 ) =  283
greenupxy  (  60, 3 ) =   37
greenupxy  (  60, 4 ) =  301
greenupxy  (  60, 5 ) =  301
greenupxy  (  60, 6 ) =  301
greenupxy  (  60, 7 ) =  274
greenupxy  (  60, 8 ) =  274
fallxy     (  60, 1 ) =  106
fallxy     (  60, 2 ) =   80
fallxy     (  60, 3 ) =  138
fallxy     (  60, 4 ) =   91
fallxy     (  60, 5 ) =   91
fallxy     (  60, 6 ) =   91
fallxy     (  60, 7 ) =   81
fallxy     (  60, 8 ) =   81
xyphendoy1 (  60, 1 ) =    2
xyphendoy1 (  60, 2 ) =    2
xyphendoy1 (  60, 3 ) =    0
xyphendoy1 (  60, 4 ) =    2
xyphendoy1 (  60, 5 ) =    2
xyphendoy1 (  60, 6 ) =    2
xyphendoy1 (  60, 7 ) =    2
xyphendoy1 (  60, 8 ) =    2
      
 ! xlat        -26.2500000000000     
greenupxy  (  59, 1 ) =  291
greenupxy  (  59, 2 ) =  282
greenupxy  (  59, 3 ) =   41
greenupxy  (  59, 4 ) =  298
greenupxy  (  59, 5 ) =  298
greenupxy  (  59, 6 ) =  298
greenupxy  (  59, 7 ) =  266
greenupxy  (  59, 8 ) =  266
fallxy     (  59, 1 ) =  102
fallxy     (  59, 2 ) =   79
fallxy     (  59, 3 ) =  145
fallxy     (  59, 4 ) =   92
fallxy     (  59, 5 ) =   92
fallxy     (  59, 6 ) =   92
fallxy     (  59, 7 ) =   81
fallxy     (  59, 8 ) =   81
xyphendoy1 (  59, 1 ) =    2
xyphendoy1 (  59, 2 ) =    2
xyphendoy1 (  59, 3 ) =    0
xyphendoy1 (  59, 4 ) =    2
xyphendoy1 (  59, 5 ) =    2
xyphendoy1 (  59, 6 ) =    2
xyphendoy1 (  59, 7 ) =    2
xyphendoy1 (  59, 8 ) =    2
      
 ! xlat        -26.7500000000000     
greenupxy  (  58, 1 ) =  289
greenupxy  (  58, 2 ) =  281
greenupxy  (  58, 3 ) =   47
greenupxy  (  58, 4 ) =  296
greenupxy  (  58, 5 ) =  296
greenupxy  (  58, 6 ) =  296
greenupxy  (  58, 7 ) =  257
greenupxy  (  58, 8 ) =  257
fallxy     (  58, 1 ) =  105
fallxy     (  58, 2 ) =   79
fallxy     (  58, 3 ) =  154
fallxy     (  58, 4 ) =   93
fallxy     (  58, 5 ) =   93
fallxy     (  58, 6 ) =   93
fallxy     (  58, 7 ) =   81
fallxy     (  58, 8 ) =   81
xyphendoy1 (  58, 1 ) =    2
xyphendoy1 (  58, 2 ) =    2
xyphendoy1 (  58, 3 ) =    0
xyphendoy1 (  58, 4 ) =    2
xyphendoy1 (  58, 5 ) =    2
xyphendoy1 (  58, 6 ) =    2
xyphendoy1 (  58, 7 ) =    2
xyphendoy1 (  58, 8 ) =    2
      
 ! xlat        -27.2500000000000     
greenupxy  (  57, 1 ) =  292
greenupxy  (  57, 2 ) =  280
greenupxy  (  57, 3 ) =   52
greenupxy  (  57, 4 ) =  294
greenupxy  (  57, 5 ) =  294
greenupxy  (  57, 6 ) =  294
greenupxy  (  57, 7 ) =  249
greenupxy  (  57, 8 ) =  249
fallxy     (  57, 1 ) =  105
fallxy     (  57, 2 ) =   79
fallxy     (  57, 3 ) =  161
fallxy     (  57, 4 ) =   93
fallxy     (  57, 5 ) =   93
fallxy     (  57, 6 ) =   93
fallxy     (  57, 7 ) =   81
fallxy     (  57, 8 ) =   81
xyphendoy1 (  57, 1 ) =    2
xyphendoy1 (  57, 2 ) =    2
xyphendoy1 (  57, 3 ) =    0
xyphendoy1 (  57, 4 ) =    2
xyphendoy1 (  57, 5 ) =    2
xyphendoy1 (  57, 6 ) =    2
xyphendoy1 (  57, 7 ) =    2
xyphendoy1 (  57, 8 ) =    2
      
 ! xlat        -27.7500000000000     
greenupxy  (  56, 1 ) =  292
greenupxy  (  56, 2 ) =  279
greenupxy  (  56, 3 ) =   58
greenupxy  (  56, 4 ) =  292
greenupxy  (  56, 5 ) =  292
greenupxy  (  56, 6 ) =  292
greenupxy  (  56, 7 ) =  241
greenupxy  (  56, 8 ) =  241
fallxy     (  56, 1 ) =  105
fallxy     (  56, 2 ) =   79
fallxy     (  56, 3 ) =  169
fallxy     (  56, 4 ) =   93
fallxy     (  56, 5 ) =   93
fallxy     (  56, 6 ) =   93
fallxy     (  56, 7 ) =   81
fallxy     (  56, 8 ) =   81
xyphendoy1 (  56, 1 ) =    2
xyphendoy1 (  56, 2 ) =    2
xyphendoy1 (  56, 3 ) =    0
xyphendoy1 (  56, 4 ) =    2
xyphendoy1 (  56, 5 ) =    2
xyphendoy1 (  56, 6 ) =    2
xyphendoy1 (  56, 7 ) =    2
xyphendoy1 (  56, 8 ) =    2
      
 ! xlat        -28.2500000000000     
greenupxy  (  55, 1 ) =  292
greenupxy  (  55, 2 ) =  279
greenupxy  (  55, 3 ) =   64
greenupxy  (  55, 4 ) =  287
greenupxy  (  55, 5 ) =  287
greenupxy  (  55, 6 ) =  287
greenupxy  (  55, 7 ) =  233
greenupxy  (  55, 8 ) =  233
fallxy     (  55, 1 ) =  102
fallxy     (  55, 2 ) =   78
fallxy     (  55, 3 ) =  179
fallxy     (  55, 4 ) =   93
fallxy     (  55, 5 ) =   93
fallxy     (  55, 6 ) =   93
fallxy     (  55, 7 ) =   82
fallxy     (  55, 8 ) =   82
xyphendoy1 (  55, 1 ) =    2
xyphendoy1 (  55, 2 ) =    2
xyphendoy1 (  55, 3 ) =    0
xyphendoy1 (  55, 4 ) =    2
xyphendoy1 (  55, 5 ) =    2
xyphendoy1 (  55, 6 ) =    2
xyphendoy1 (  55, 7 ) =    2
xyphendoy1 (  55, 8 ) =    2
      
 ! xlat        -28.7500000000000     
greenupxy  (  54, 1 ) =  292
greenupxy  (  54, 2 ) =  278
greenupxy  (  54, 3 ) =   70
greenupxy  (  54, 4 ) =  286
greenupxy  (  54, 5 ) =  286
greenupxy  (  54, 6 ) =  286
greenupxy  (  54, 7 ) =  225
greenupxy  (  54, 8 ) =  225
fallxy     (  54, 1 ) =  102
fallxy     (  54, 2 ) =   78
fallxy     (  54, 3 ) =  188
fallxy     (  54, 4 ) =   94
fallxy     (  54, 5 ) =   94
fallxy     (  54, 6 ) =   94
fallxy     (  54, 7 ) =   83
fallxy     (  54, 8 ) =   83
xyphendoy1 (  54, 1 ) =    2
xyphendoy1 (  54, 2 ) =    2
xyphendoy1 (  54, 3 ) =    0
xyphendoy1 (  54, 4 ) =    2
xyphendoy1 (  54, 5 ) =    2
xyphendoy1 (  54, 6 ) =    2
xyphendoy1 (  54, 7 ) =    2
xyphendoy1 (  54, 8 ) =    2
      
 ! xlat        -29.2500000000000     
greenupxy  (  53, 1 ) =  297
greenupxy  (  53, 2 ) =  278
greenupxy  (  53, 3 ) =   75
greenupxy  (  53, 4 ) =  284
greenupxy  (  53, 5 ) =  284
greenupxy  (  53, 6 ) =  284
greenupxy  (  53, 7 ) =  218
greenupxy  (  53, 8 ) =  218
fallxy     (  53, 1 ) =   97
fallxy     (  53, 2 ) =   78
fallxy     (  53, 3 ) =  196
fallxy     (  53, 4 ) =   93
fallxy     (  53, 5 ) =   93
fallxy     (  53, 6 ) =   93
fallxy     (  53, 7 ) =   83
fallxy     (  53, 8 ) =   83
xyphendoy1 (  53, 1 ) =    2
xyphendoy1 (  53, 2 ) =    2
xyphendoy1 (  53, 3 ) =    0
xyphendoy1 (  53, 4 ) =    2
xyphendoy1 (  53, 5 ) =    2
xyphendoy1 (  53, 6 ) =    2
xyphendoy1 (  53, 7 ) =    2
xyphendoy1 (  53, 8 ) =    2
      
 ! xlat        -29.7500000000000     
greenupxy  (  52, 1 ) =  290
greenupxy  (  52, 2 ) =  277
greenupxy  (  52, 3 ) =   79
greenupxy  (  52, 4 ) =  281
greenupxy  (  52, 5 ) =  281
greenupxy  (  52, 6 ) =  281
greenupxy  (  52, 7 ) =  210
greenupxy  (  52, 8 ) =  210
fallxy     (  52, 1 ) =   96
fallxy     (  52, 2 ) =   77
fallxy     (  52, 3 ) =  204
fallxy     (  52, 4 ) =   93
fallxy     (  52, 5 ) =   93
fallxy     (  52, 6 ) =   93
fallxy     (  52, 7 ) =   82
fallxy     (  52, 8 ) =   82
xyphendoy1 (  52, 1 ) =    2
xyphendoy1 (  52, 2 ) =    2
xyphendoy1 (  52, 3 ) =    0
xyphendoy1 (  52, 4 ) =    2
xyphendoy1 (  52, 5 ) =    2
xyphendoy1 (  52, 6 ) =    2
xyphendoy1 (  52, 7 ) =    2
xyphendoy1 (  52, 8 ) =    2
      
 ! xlat        -30.2500000000000     
greenupxy  (  51, 1 ) =  290
greenupxy  (  51, 2 ) =  276
greenupxy  (  51, 3 ) =   82
greenupxy  (  51, 4 ) =  280
greenupxy  (  51, 5 ) =  280
greenupxy  (  51, 6 ) =  280
greenupxy  (  51, 7 ) =  203
greenupxy  (  51, 8 ) =  203
fallxy     (  51, 1 ) =   94
fallxy     (  51, 2 ) =   77
fallxy     (  51, 3 ) =  211
fallxy     (  51, 4 ) =   92
fallxy     (  51, 5 ) =   92
fallxy     (  51, 6 ) =   92
fallxy     (  51, 7 ) =   80
fallxy     (  51, 8 ) =   80
xyphendoy1 (  51, 1 ) =    2
xyphendoy1 (  51, 2 ) =    2
xyphendoy1 (  51, 3 ) =    0
xyphendoy1 (  51, 4 ) =    2
xyphendoy1 (  51, 5 ) =    2
xyphendoy1 (  51, 6 ) =    2
xyphendoy1 (  51, 7 ) =    2
xyphendoy1 (  51, 8 ) =    2
      
 ! xlat        -30.7500000000000     
greenupxy  (  50, 1 ) =  285
greenupxy  (  50, 2 ) =  275
greenupxy  (  50, 3 ) =   85
greenupxy  (  50, 4 ) =  280
greenupxy  (  50, 5 ) =  280
greenupxy  (  50, 6 ) =  280
greenupxy  (  50, 7 ) =  196
greenupxy  (  50, 8 ) =  196
fallxy     (  50, 1 ) =   90
fallxy     (  50, 2 ) =   77
fallxy     (  50, 3 ) =  216
fallxy     (  50, 4 ) =   91
fallxy     (  50, 5 ) =   91
fallxy     (  50, 6 ) =   91
fallxy     (  50, 7 ) =   83
fallxy     (  50, 8 ) =   83
xyphendoy1 (  50, 1 ) =    2
xyphendoy1 (  50, 2 ) =    2
xyphendoy1 (  50, 3 ) =    0
xyphendoy1 (  50, 4 ) =    2
xyphendoy1 (  50, 5 ) =    2
xyphendoy1 (  50, 6 ) =    2
xyphendoy1 (  50, 7 ) =    2
xyphendoy1 (  50, 8 ) =    2
      
 ! xlat        -31.2500000000000     
greenupxy  (  49, 1 ) =  284
greenupxy  (  49, 2 ) =  274
greenupxy  (  49, 3 ) =   81
greenupxy  (  49, 4 ) =  281
greenupxy  (  49, 5 ) =  281
greenupxy  (  49, 6 ) =  281
greenupxy  (  49, 7 ) =  188
greenupxy  (  49, 8 ) =  188
fallxy     (  49, 1 ) =   88
fallxy     (  49, 2 ) =   76
fallxy     (  49, 3 ) =  213
fallxy     (  49, 4 ) =   91
fallxy     (  49, 5 ) =   91
fallxy     (  49, 6 ) =   91
fallxy     (  49, 7 ) =   86
fallxy     (  49, 8 ) =   86
xyphendoy1 (  49, 1 ) =    2
xyphendoy1 (  49, 2 ) =    2
xyphendoy1 (  49, 3 ) =    0
xyphendoy1 (  49, 4 ) =    2
xyphendoy1 (  49, 5 ) =    2
xyphendoy1 (  49, 6 ) =    2
xyphendoy1 (  49, 7 ) =    2
xyphendoy1 (  49, 8 ) =    2
      
 ! xlat        -31.7500000000000     
greenupxy  (  48, 1 ) =  282
greenupxy  (  48, 2 ) =  273
greenupxy  (  48, 3 ) =   75
greenupxy  (  48, 4 ) =  280
greenupxy  (  48, 5 ) =  280
greenupxy  (  48, 6 ) =  280
greenupxy  (  48, 7 ) =  180
greenupxy  (  48, 8 ) =  180
fallxy     (  48, 1 ) =   85
fallxy     (  48, 2 ) =   76
fallxy     (  48, 3 ) =  212
fallxy     (  48, 4 ) =   90
fallxy     (  48, 5 ) =   90
fallxy     (  48, 6 ) =   90
fallxy     (  48, 7 ) =   83
fallxy     (  48, 8 ) =   83
xyphendoy1 (  48, 1 ) =    2
xyphendoy1 (  48, 2 ) =    2
xyphendoy1 (  48, 3 ) =    0
xyphendoy1 (  48, 4 ) =    2
xyphendoy1 (  48, 5 ) =    2
xyphendoy1 (  48, 6 ) =    2
xyphendoy1 (  48, 7 ) =    2
xyphendoy1 (  48, 8 ) =    2
      
 ! xlat        -32.2500000000000     
greenupxy  (  47, 1 ) =  279
greenupxy  (  47, 2 ) =  272
greenupxy  (  47, 3 ) =   68
greenupxy  (  47, 4 ) =  279
greenupxy  (  47, 5 ) =  279
greenupxy  (  47, 6 ) =  279
greenupxy  (  47, 7 ) =  174
greenupxy  (  47, 8 ) =  174
fallxy     (  47, 1 ) =   81
fallxy     (  47, 2 ) =   75
fallxy     (  47, 3 ) =  211
fallxy     (  47, 4 ) =   89
fallxy     (  47, 5 ) =   89
fallxy     (  47, 6 ) =   89
fallxy     (  47, 7 ) =   83
fallxy     (  47, 8 ) =   83
xyphendoy1 (  47, 1 ) =    2
xyphendoy1 (  47, 2 ) =    2
xyphendoy1 (  47, 3 ) =    0
xyphendoy1 (  47, 4 ) =    2
xyphendoy1 (  47, 5 ) =    2
xyphendoy1 (  47, 6 ) =    2
xyphendoy1 (  47, 7 ) =    2
xyphendoy1 (  47, 8 ) =    2
      
 ! xlat        -32.7500000000000     
greenupxy  (  46, 1 ) =  273
greenupxy  (  46, 2 ) =  270
greenupxy  (  46, 3 ) =   61
greenupxy  (  46, 4 ) =  276
greenupxy  (  46, 5 ) =  276
greenupxy  (  46, 6 ) =  276
greenupxy  (  46, 7 ) =  173
greenupxy  (  46, 8 ) =  173
fallxy     (  46, 1 ) =   77
fallxy     (  46, 2 ) =   73
fallxy     (  46, 3 ) =  208
fallxy     (  46, 4 ) =   88
fallxy     (  46, 5 ) =   88
fallxy     (  46, 6 ) =   88
fallxy     (  46, 7 ) =   83
fallxy     (  46, 8 ) =   83
xyphendoy1 (  46, 1 ) =    2
xyphendoy1 (  46, 2 ) =    2
xyphendoy1 (  46, 3 ) =    0
xyphendoy1 (  46, 4 ) =    2
xyphendoy1 (  46, 5 ) =    2
xyphendoy1 (  46, 6 ) =    2
xyphendoy1 (  46, 7 ) =    2
xyphendoy1 (  46, 8 ) =    2
      
 ! xlat        -33.2500000000000     
greenupxy  (  45, 1 ) =  273
greenupxy  (  45, 2 ) =  269
greenupxy  (  45, 3 ) =   52
greenupxy  (  45, 4 ) =  274
greenupxy  (  45, 5 ) =  274
greenupxy  (  45, 6 ) =  274
greenupxy  (  45, 7 ) =  172
greenupxy  (  45, 8 ) =  172
fallxy     (  45, 1 ) =   79
fallxy     (  45, 2 ) =   73
fallxy     (  45, 3 ) =  207
fallxy     (  45, 4 ) =   87
fallxy     (  45, 5 ) =   87
fallxy     (  45, 6 ) =   87
fallxy     (  45, 7 ) =   83
fallxy     (  45, 8 ) =   83
xyphendoy1 (  45, 1 ) =    2
xyphendoy1 (  45, 2 ) =    2
xyphendoy1 (  45, 3 ) =    0
xyphendoy1 (  45, 4 ) =    2
xyphendoy1 (  45, 5 ) =    2
xyphendoy1 (  45, 6 ) =    2
xyphendoy1 (  45, 7 ) =    2
xyphendoy1 (  45, 8 ) =    2
      
 ! xlat        -33.7500000000000     
greenupxy  (  44, 1 ) =  274
greenupxy  (  44, 2 ) =  269
greenupxy  (  44, 3 ) =   41
greenupxy  (  44, 4 ) =  271
greenupxy  (  44, 5 ) =  271
greenupxy  (  44, 6 ) =  271
greenupxy  (  44, 7 ) =  177
greenupxy  (  44, 8 ) =  177
fallxy     (  44, 1 ) =   73
fallxy     (  44, 2 ) =   73
fallxy     (  44, 3 ) =  199
fallxy     (  44, 4 ) =   82
fallxy     (  44, 5 ) =   82
fallxy     (  44, 6 ) =   82
fallxy     (  44, 7 ) =   83
fallxy     (  44, 8 ) =   83
xyphendoy1 (  44, 1 ) =    2
xyphendoy1 (  44, 2 ) =    2
xyphendoy1 (  44, 3 ) =    0
xyphendoy1 (  44, 4 ) =    2
xyphendoy1 (  44, 5 ) =    2
xyphendoy1 (  44, 6 ) =    2
xyphendoy1 (  44, 7 ) =    2
xyphendoy1 (  44, 8 ) =    2
      
 ! xlat        -34.2500000000000     
greenupxy  (  43, 1 ) =  277
greenupxy  (  43, 2 ) =  268
greenupxy  (  43, 3 ) =   33
greenupxy  (  43, 4 ) =  270
greenupxy  (  43, 5 ) =  270
greenupxy  (  43, 6 ) =  270
greenupxy  (  43, 7 ) =  182
greenupxy  (  43, 8 ) =  182
fallxy     (  43, 1 ) =   73
fallxy     (  43, 2 ) =   73
fallxy     (  43, 3 ) =  198
fallxy     (  43, 4 ) =   81
fallxy     (  43, 5 ) =   81
fallxy     (  43, 6 ) =   81
fallxy     (  43, 7 ) =   83
fallxy     (  43, 8 ) =   83
xyphendoy1 (  43, 1 ) =    2
xyphendoy1 (  43, 2 ) =    2
xyphendoy1 (  43, 3 ) =    0
xyphendoy1 (  43, 4 ) =    2
xyphendoy1 (  43, 5 ) =    2
xyphendoy1 (  43, 6 ) =    2
xyphendoy1 (  43, 7 ) =    2
xyphendoy1 (  43, 8 ) =    2
      
 ! xlat        -34.7500000000000     
greenupxy  (  42, 1 ) =  270
greenupxy  (  42, 2 ) =  267
greenupxy  (  42, 3 ) =   23
greenupxy  (  42, 4 ) =  262
greenupxy  (  42, 5 ) =  262
greenupxy  (  42, 6 ) =  262
greenupxy  (  42, 7 ) =  187
greenupxy  (  42, 8 ) =  187
fallxy     (  42, 1 ) =   69
fallxy     (  42, 2 ) =   72
fallxy     (  42, 3 ) =  197
fallxy     (  42, 4 ) =   79
fallxy     (  42, 5 ) =   79
fallxy     (  42, 6 ) =   79
fallxy     (  42, 7 ) =   83
fallxy     (  42, 8 ) =   83
xyphendoy1 (  42, 1 ) =    2
xyphendoy1 (  42, 2 ) =    2
xyphendoy1 (  42, 3 ) =    0
xyphendoy1 (  42, 4 ) =    2
xyphendoy1 (  42, 5 ) =    2
xyphendoy1 (  42, 6 ) =    2
xyphendoy1 (  42, 7 ) =    2
xyphendoy1 (  42, 8 ) =    2
      
 ! xlat        -35.2500000000000     
greenupxy  (  41, 1 ) =  273
greenupxy  (  41, 2 ) =  265
greenupxy  (  41, 3 ) =   14
greenupxy  (  41, 4 ) =  261
greenupxy  (  41, 5 ) =  261
greenupxy  (  41, 6 ) =  261
greenupxy  (  41, 7 ) =  193
greenupxy  (  41, 8 ) =  193
fallxy     (  41, 1 ) =   68
fallxy     (  41, 2 ) =   69
fallxy     (  41, 3 ) =  193
fallxy     (  41, 4 ) =   75
fallxy     (  41, 5 ) =   75
fallxy     (  41, 6 ) =   75
fallxy     (  41, 7 ) =   21
fallxy     (  41, 8 ) =   21
xyphendoy1 (  41, 1 ) =    2
xyphendoy1 (  41, 2 ) =    2
xyphendoy1 (  41, 3 ) =    0
xyphendoy1 (  41, 4 ) =    2
xyphendoy1 (  41, 5 ) =    2
xyphendoy1 (  41, 6 ) =    2
xyphendoy1 (  41, 7 ) =    2
xyphendoy1 (  41, 8 ) =    2
      
 ! xlat        -35.7500000000000     
greenupxy  (  40, 1 ) =  272
greenupxy  (  40, 2 ) =  264
greenupxy  (  40, 3 ) =    7
greenupxy  (  40, 4 ) =  261
greenupxy  (  40, 5 ) =  261
greenupxy  (  40, 6 ) =  261
greenupxy  (  40, 7 ) =  194
greenupxy  (  40, 8 ) =  194
fallxy     (  40, 1 ) =   68
fallxy     (  40, 2 ) =   64
fallxy     (  40, 3 ) =  191
fallxy     (  40, 4 ) =   74
fallxy     (  40, 5 ) =   74
fallxy     (  40, 6 ) =   74
fallxy     (  40, 7 ) =   21
fallxy     (  40, 8 ) =   21
xyphendoy1 (  40, 1 ) =    2
xyphendoy1 (  40, 2 ) =    2
xyphendoy1 (  40, 3 ) =    0
xyphendoy1 (  40, 4 ) =    2
xyphendoy1 (  40, 5 ) =    2
xyphendoy1 (  40, 6 ) =    2
xyphendoy1 (  40, 7 ) =    2
xyphendoy1 (  40, 8 ) =    2
      
 ! xlat        -36.2500000000000     
greenupxy  (  39, 1 ) =  276
greenupxy  (  39, 2 ) =  262
greenupxy  (  39, 3 ) =  363
greenupxy  (  39, 4 ) =  259
greenupxy  (  39, 5 ) =  259
greenupxy  (  39, 6 ) =  259
greenupxy  (  39, 7 ) =  196
greenupxy  (  39, 8 ) =  196
fallxy     (  39, 1 ) =   68
fallxy     (  39, 2 ) =   59
fallxy     (  39, 3 ) =  187
fallxy     (  39, 4 ) =   72
fallxy     (  39, 5 ) =   72
fallxy     (  39, 6 ) =   72
fallxy     (  39, 7 ) =   21
fallxy     (  39, 8 ) =   21
xyphendoy1 (  39, 1 ) =    2
xyphendoy1 (  39, 2 ) =    2
xyphendoy1 (  39, 3 ) =    2
xyphendoy1 (  39, 4 ) =    2
xyphendoy1 (  39, 5 ) =    2
xyphendoy1 (  39, 6 ) =    2
xyphendoy1 (  39, 7 ) =    2
xyphendoy1 (  39, 8 ) =    2
      
 ! xlat        -36.7500000000000     
greenupxy  (  38, 1 ) =  272
greenupxy  (  38, 2 ) =  259
greenupxy  (  38, 3 ) =  347
greenupxy  (  38, 4 ) =  262
greenupxy  (  38, 5 ) =  262
greenupxy  (  38, 6 ) =  262
greenupxy  (  38, 7 ) =  198
greenupxy  (  38, 8 ) =  198
fallxy     (  38, 1 ) =   68
fallxy     (  38, 2 ) =   58
fallxy     (  38, 3 ) =  181
fallxy     (  38, 4 ) =   70
fallxy     (  38, 5 ) =   70
fallxy     (  38, 6 ) =   70
fallxy     (  38, 7 ) =   21
fallxy     (  38, 8 ) =   21
xyphendoy1 (  38, 1 ) =    2
xyphendoy1 (  38, 2 ) =    2
xyphendoy1 (  38, 3 ) =    2
xyphendoy1 (  38, 4 ) =    2
xyphendoy1 (  38, 5 ) =    2
xyphendoy1 (  38, 6 ) =    2
xyphendoy1 (  38, 7 ) =    2
xyphendoy1 (  38, 8 ) =    2
      
 ! xlat        -37.2500000000000     
greenupxy  (  37, 1 ) =  273
greenupxy  (  37, 2 ) =  257
greenupxy  (  37, 3 ) =  329
greenupxy  (  37, 4 ) =  260
greenupxy  (  37, 5 ) =  260
greenupxy  (  37, 6 ) =  260
greenupxy  (  37, 7 ) =  199
greenupxy  (  37, 8 ) =  199
fallxy     (  37, 1 ) =   67
fallxy     (  37, 2 ) =   53
fallxy     (  37, 3 ) =  171
fallxy     (  37, 4 ) =   68
fallxy     (  37, 5 ) =   68
fallxy     (  37, 6 ) =   68
fallxy     (  37, 7 ) =   21
fallxy     (  37, 8 ) =   21
xyphendoy1 (  37, 1 ) =    2
xyphendoy1 (  37, 2 ) =    2
xyphendoy1 (  37, 3 ) =    2
xyphendoy1 (  37, 4 ) =    2
xyphendoy1 (  37, 5 ) =    2
xyphendoy1 (  37, 6 ) =    2
xyphendoy1 (  37, 7 ) =    2
xyphendoy1 (  37, 8 ) =    2
      
 ! xlat        -37.7500000000000     
greenupxy  (  36, 1 ) =  270
greenupxy  (  36, 2 ) =  255
greenupxy  (  36, 3 ) =  310
greenupxy  (  36, 4 ) =  258
greenupxy  (  36, 5 ) =  258
greenupxy  (  36, 6 ) =  258
greenupxy  (  36, 7 ) =  199
greenupxy  (  36, 8 ) =  199
fallxy     (  36, 1 ) =   69
fallxy     (  36, 2 ) =   48
fallxy     (  36, 3 ) =  160
fallxy     (  36, 4 ) =   61
fallxy     (  36, 5 ) =   61
fallxy     (  36, 6 ) =   61
fallxy     (  36, 7 ) =   21
fallxy     (  36, 8 ) =   21
xyphendoy1 (  36, 1 ) =    2
xyphendoy1 (  36, 2 ) =    2
xyphendoy1 (  36, 3 ) =    2
xyphendoy1 (  36, 4 ) =    2
xyphendoy1 (  36, 5 ) =    2
xyphendoy1 (  36, 6 ) =    2
xyphendoy1 (  36, 7 ) =    2
xyphendoy1 (  36, 8 ) =    2
      
 ! xlat        -38.2500000000000     
greenupxy  (  35, 1 ) =  274
greenupxy  (  35, 2 ) =  255
greenupxy  (  35, 3 ) =  307
greenupxy  (  35, 4 ) =  260
greenupxy  (  35, 5 ) =  260
greenupxy  (  35, 6 ) =  260
greenupxy  (  35, 7 ) =  202
greenupxy  (  35, 8 ) =  202
fallxy     (  35, 1 ) =   67
fallxy     (  35, 2 ) =   45
fallxy     (  35, 3 ) =  151
fallxy     (  35, 4 ) =   59
fallxy     (  35, 5 ) =   59
fallxy     (  35, 6 ) =   59
fallxy     (  35, 7 ) =   21
fallxy     (  35, 8 ) =   21
xyphendoy1 (  35, 1 ) =    2
xyphendoy1 (  35, 2 ) =    2
xyphendoy1 (  35, 3 ) =    2
xyphendoy1 (  35, 4 ) =    2
xyphendoy1 (  35, 5 ) =    2
xyphendoy1 (  35, 6 ) =    2
xyphendoy1 (  35, 7 ) =    2
xyphendoy1 (  35, 8 ) =    2
      
 ! xlat        -38.7500000000000     
greenupxy  (  34, 1 ) =  273
greenupxy  (  34, 2 ) =  253
greenupxy  (  34, 3 ) =  304
greenupxy  (  34, 4 ) =  260
greenupxy  (  34, 5 ) =  260
greenupxy  (  34, 6 ) =  260
greenupxy  (  34, 7 ) =  202
greenupxy  (  34, 8 ) =  202
fallxy     (  34, 1 ) =   67
fallxy     (  34, 2 ) =   43
fallxy     (  34, 3 ) =  141
fallxy     (  34, 4 ) =   58
fallxy     (  34, 5 ) =   58
fallxy     (  34, 6 ) =   58
fallxy     (  34, 7 ) =   14
fallxy     (  34, 8 ) =   14
xyphendoy1 (  34, 1 ) =    2
xyphendoy1 (  34, 2 ) =    2
xyphendoy1 (  34, 3 ) =    2
xyphendoy1 (  34, 4 ) =    2
xyphendoy1 (  34, 5 ) =    2
xyphendoy1 (  34, 6 ) =    2
xyphendoy1 (  34, 7 ) =    2
xyphendoy1 (  34, 8 ) =    2
      
 ! xlat        -39.2500000000000     
greenupxy  (  33, 1 ) =  274
greenupxy  (  33, 2 ) =  254
greenupxy  (  33, 3 ) =  299
greenupxy  (  33, 4 ) =  255
greenupxy  (  33, 5 ) =  255
greenupxy  (  33, 6 ) =  255
greenupxy  (  33, 7 ) =  201
greenupxy  (  33, 8 ) =  201
fallxy     (  33, 1 ) =   66
fallxy     (  33, 2 ) =   40
fallxy     (  33, 3 ) =  141
fallxy     (  33, 4 ) =   56
fallxy     (  33, 5 ) =   56
fallxy     (  33, 6 ) =   56
fallxy     (  33, 7 ) =   14
fallxy     (  33, 8 ) =   14
xyphendoy1 (  33, 1 ) =    2
xyphendoy1 (  33, 2 ) =    2
xyphendoy1 (  33, 3 ) =    2
xyphendoy1 (  33, 4 ) =    2
xyphendoy1 (  33, 5 ) =    2
xyphendoy1 (  33, 6 ) =    2
xyphendoy1 (  33, 7 ) =    2
xyphendoy1 (  33, 8 ) =    2
      
 ! xlat        -39.7500000000000     
greenupxy  (  32, 1 ) =  275
greenupxy  (  32, 2 ) =  253
greenupxy  (  32, 3 ) =  293
greenupxy  (  32, 4 ) =  250
greenupxy  (  32, 5 ) =  250
greenupxy  (  32, 6 ) =  250
greenupxy  (  32, 7 ) =  205
greenupxy  (  32, 8 ) =  205
fallxy     (  32, 1 ) =   65
fallxy     (  32, 2 ) =   38
fallxy     (  32, 3 ) =  136
fallxy     (  32, 4 ) =   53
fallxy     (  32, 5 ) =   53
fallxy     (  32, 6 ) =   53
fallxy     (  32, 7 ) =   14
fallxy     (  32, 8 ) =   14
xyphendoy1 (  32, 1 ) =    2
xyphendoy1 (  32, 2 ) =    2
xyphendoy1 (  32, 3 ) =    2
xyphendoy1 (  32, 4 ) =    2
xyphendoy1 (  32, 5 ) =    2
xyphendoy1 (  32, 6 ) =    2
xyphendoy1 (  32, 7 ) =    2
xyphendoy1 (  32, 8 ) =    2
      
 ! xlat        -40.2500000000000     
greenupxy  (  31, 1 ) =  276
greenupxy  (  31, 2 ) =  254
greenupxy  (  31, 3 ) =  288
greenupxy  (  31, 4 ) =  246
greenupxy  (  31, 5 ) =  246
greenupxy  (  31, 6 ) =  246
greenupxy  (  31, 7 ) =  211
greenupxy  (  31, 8 ) =  211
fallxy     (  31, 1 ) =   64
fallxy     (  31, 2 ) =   36
fallxy     (  31, 3 ) =  133
fallxy     (  31, 4 ) =   53
fallxy     (  31, 5 ) =   53
fallxy     (  31, 6 ) =   53
fallxy     (  31, 7 ) =   11
fallxy     (  31, 8 ) =   11
xyphendoy1 (  31, 1 ) =    2
xyphendoy1 (  31, 2 ) =    2
xyphendoy1 (  31, 3 ) =    2
xyphendoy1 (  31, 4 ) =    2
xyphendoy1 (  31, 5 ) =    2
xyphendoy1 (  31, 6 ) =    2
xyphendoy1 (  31, 7 ) =    2
xyphendoy1 (  31, 8 ) =    2
      
 ! xlat        -40.7500000000000     
greenupxy  (  30, 1 ) =  277
greenupxy  (  30, 2 ) =  255
greenupxy  (  30, 3 ) =  285
greenupxy  (  30, 4 ) =  245
greenupxy  (  30, 5 ) =  245
greenupxy  (  30, 6 ) =  245
greenupxy  (  30, 7 ) =  217
greenupxy  (  30, 8 ) =  217
fallxy     (  30, 1 ) =   64
fallxy     (  30, 2 ) =   34
fallxy     (  30, 3 ) =  128
fallxy     (  30, 4 ) =   53
fallxy     (  30, 5 ) =   53
fallxy     (  30, 6 ) =   53
fallxy     (  30, 7 ) =   11
fallxy     (  30, 8 ) =   11
xyphendoy1 (  30, 1 ) =    2
xyphendoy1 (  30, 2 ) =    2
xyphendoy1 (  30, 3 ) =    2
xyphendoy1 (  30, 4 ) =    2
xyphendoy1 (  30, 5 ) =    2
xyphendoy1 (  30, 6 ) =    2
xyphendoy1 (  30, 7 ) =    2
xyphendoy1 (  30, 8 ) =    2
      
 ! xlat        -41.2500000000000     
greenupxy  (  29, 1 ) =  279
greenupxy  (  29, 2 ) =  256
greenupxy  (  29, 3 ) =  280
greenupxy  (  29, 4 ) =  244
greenupxy  (  29, 5 ) =  244
greenupxy  (  29, 6 ) =  244
greenupxy  (  29, 7 ) =  217
greenupxy  (  29, 8 ) =  217
fallxy     (  29, 1 ) =   65
fallxy     (  29, 2 ) =   33
fallxy     (  29, 3 ) =  116
fallxy     (  29, 4 ) =   52
fallxy     (  29, 5 ) =   52
fallxy     (  29, 6 ) =   52
fallxy     (  29, 7 ) =   15
fallxy     (  29, 8 ) =   15
xyphendoy1 (  29, 1 ) =    2
xyphendoy1 (  29, 2 ) =    2
xyphendoy1 (  29, 3 ) =    2
xyphendoy1 (  29, 4 ) =    2
xyphendoy1 (  29, 5 ) =    2
xyphendoy1 (  29, 6 ) =    2
xyphendoy1 (  29, 7 ) =    2
xyphendoy1 (  29, 8 ) =    2
      
 ! xlat        -41.7500000000000     
greenupxy  (  28, 1 ) =  280
greenupxy  (  28, 2 ) =  256
greenupxy  (  28, 3 ) =  278
greenupxy  (  28, 4 ) =  244
greenupxy  (  28, 5 ) =  244
greenupxy  (  28, 6 ) =  244
greenupxy  (  28, 7 ) =  217
greenupxy  (  28, 8 ) =  217
fallxy     (  28, 1 ) =   62
fallxy     (  28, 2 ) =   34
fallxy     (  28, 3 ) =  108
fallxy     (  28, 4 ) =   52
fallxy     (  28, 5 ) =   52
fallxy     (  28, 6 ) =   52
fallxy     (  28, 7 ) =   18
fallxy     (  28, 8 ) =   18
xyphendoy1 (  28, 1 ) =    2
xyphendoy1 (  28, 2 ) =    2
xyphendoy1 (  28, 3 ) =    2
xyphendoy1 (  28, 4 ) =    2
xyphendoy1 (  28, 5 ) =    2
xyphendoy1 (  28, 6 ) =    2
xyphendoy1 (  28, 7 ) =    2
xyphendoy1 (  28, 8 ) =    2
      
 ! xlat        -42.2500000000000     
greenupxy  (  27, 1 ) =  281
greenupxy  (  27, 2 ) =  257
greenupxy  (  27, 3 ) =  276
greenupxy  (  27, 4 ) =  245
greenupxy  (  27, 5 ) =  245
greenupxy  (  27, 6 ) =  245
greenupxy  (  27, 7 ) =  218
greenupxy  (  27, 8 ) =  218
fallxy     (  27, 1 ) =   60
fallxy     (  27, 2 ) =   35
fallxy     (  27, 3 ) =  107
fallxy     (  27, 4 ) =   55
fallxy     (  27, 5 ) =   55
fallxy     (  27, 6 ) =   55
fallxy     (  27, 7 ) =   21
fallxy     (  27, 8 ) =   21
xyphendoy1 (  27, 1 ) =    2
xyphendoy1 (  27, 2 ) =    2
xyphendoy1 (  27, 3 ) =    2
xyphendoy1 (  27, 4 ) =    2
xyphendoy1 (  27, 5 ) =    2
xyphendoy1 (  27, 6 ) =    2
xyphendoy1 (  27, 7 ) =    2
xyphendoy1 (  27, 8 ) =    2
      
 ! xlat        -42.7500000000000     
greenupxy  (  26, 1 ) =  281
greenupxy  (  26, 2 ) =  258
greenupxy  (  26, 3 ) =  272
greenupxy  (  26, 4 ) =  247
greenupxy  (  26, 5 ) =  247
greenupxy  (  26, 6 ) =  247
greenupxy  (  26, 7 ) =  220
greenupxy  (  26, 8 ) =  220
fallxy     (  26, 1 ) =   60
fallxy     (  26, 2 ) =   36
fallxy     (  26, 3 ) =  107
fallxy     (  26, 4 ) =   56
fallxy     (  26, 5 ) =   56
fallxy     (  26, 6 ) =   56
fallxy     (  26, 7 ) =   27
fallxy     (  26, 8 ) =   27
xyphendoy1 (  26, 1 ) =    2
xyphendoy1 (  26, 2 ) =    2
xyphendoy1 (  26, 3 ) =    2
xyphendoy1 (  26, 4 ) =    2
xyphendoy1 (  26, 5 ) =    2
xyphendoy1 (  26, 6 ) =    2
xyphendoy1 (  26, 7 ) =    2
xyphendoy1 (  26, 8 ) =    2
      
 ! xlat        -43.2500000000000     
greenupxy  (  25, 1 ) =  284
greenupxy  (  25, 2 ) =  260
greenupxy  (  25, 3 ) =  268
greenupxy  (  25, 4 ) =  256
greenupxy  (  25, 5 ) =  256
greenupxy  (  25, 6 ) =  256
greenupxy  (  25, 7 ) =  224
greenupxy  (  25, 8 ) =  224
fallxy     (  25, 1 ) =   62
fallxy     (  25, 2 ) =   37
fallxy     (  25, 3 ) =   99
fallxy     (  25, 4 ) =   54
fallxy     (  25, 5 ) =   54
fallxy     (  25, 6 ) =   54
fallxy     (  25, 7 ) =   30
fallxy     (  25, 8 ) =   30
xyphendoy1 (  25, 1 ) =    2
xyphendoy1 (  25, 2 ) =    2
xyphendoy1 (  25, 3 ) =    2
xyphendoy1 (  25, 4 ) =    2
xyphendoy1 (  25, 5 ) =    2
xyphendoy1 (  25, 6 ) =    2
xyphendoy1 (  25, 7 ) =    2
xyphendoy1 (  25, 8 ) =    2
      
 ! xlat        -43.7500000000000     
greenupxy  (  24, 1 ) =  282
greenupxy  (  24, 2 ) =  262
greenupxy  (  24, 3 ) =  263
greenupxy  (  24, 4 ) =  258
greenupxy  (  24, 5 ) =  258
greenupxy  (  24, 6 ) =  258
greenupxy  (  24, 7 ) =  229
greenupxy  (  24, 8 ) =  229
fallxy     (  24, 1 ) =   62
fallxy     (  24, 2 ) =   38
fallxy     (  24, 3 ) =   83
fallxy     (  24, 4 ) =   55
fallxy     (  24, 5 ) =   55
fallxy     (  24, 6 ) =   55
fallxy     (  24, 7 ) =   35
fallxy     (  24, 8 ) =   35
xyphendoy1 (  24, 1 ) =    2
xyphendoy1 (  24, 2 ) =    2
xyphendoy1 (  24, 3 ) =    2
xyphendoy1 (  24, 4 ) =    2
xyphendoy1 (  24, 5 ) =    2
xyphendoy1 (  24, 6 ) =    2
xyphendoy1 (  24, 7 ) =    2
xyphendoy1 (  24, 8 ) =    2
      
 ! xlat        -44.2500000000000     
greenupxy  (  23, 1 ) =  279
greenupxy  (  23, 2 ) =  263
greenupxy  (  23, 3 ) =  264
greenupxy  (  23, 4 ) =  260
greenupxy  (  23, 5 ) =  260
greenupxy  (  23, 6 ) =  260
greenupxy  (  23, 7 ) =  234
greenupxy  (  23, 8 ) =  234
fallxy     (  23, 1 ) =   61
fallxy     (  23, 2 ) =   41
fallxy     (  23, 3 ) =   64
fallxy     (  23, 4 ) =   55
fallxy     (  23, 5 ) =   55
fallxy     (  23, 6 ) =   55
fallxy     (  23, 7 ) =   35
fallxy     (  23, 8 ) =   35
xyphendoy1 (  23, 1 ) =    2
xyphendoy1 (  23, 2 ) =    2
xyphendoy1 (  23, 3 ) =    2
xyphendoy1 (  23, 4 ) =    2
xyphendoy1 (  23, 5 ) =    2
xyphendoy1 (  23, 6 ) =    2
xyphendoy1 (  23, 7 ) =    2
xyphendoy1 (  23, 8 ) =    2
      
 ! xlat        -44.7500000000000     
greenupxy  (  22, 1 ) =  274
greenupxy  (  22, 2 ) =  266
greenupxy  (  22, 3 ) =  265
greenupxy  (  22, 4 ) =  263
greenupxy  (  22, 5 ) =  263
greenupxy  (  22, 6 ) =  263
greenupxy  (  22, 7 ) =  237
greenupxy  (  22, 8 ) =  237
fallxy     (  22, 1 ) =   59
fallxy     (  22, 2 ) =   44
fallxy     (  22, 3 ) =   64
fallxy     (  22, 4 ) =   53
fallxy     (  22, 5 ) =   53
fallxy     (  22, 6 ) =   53
fallxy     (  22, 7 ) =   39
fallxy     (  22, 8 ) =   39
xyphendoy1 (  22, 1 ) =    2
xyphendoy1 (  22, 2 ) =    2
xyphendoy1 (  22, 3 ) =    2
xyphendoy1 (  22, 4 ) =    2
xyphendoy1 (  22, 5 ) =    2
xyphendoy1 (  22, 6 ) =    2
xyphendoy1 (  22, 7 ) =    2
xyphendoy1 (  22, 8 ) =    2
      
 ! xlat        -45.2500000000000     
greenupxy  (  21, 1 ) =  282
greenupxy  (  21, 2 ) =  271
greenupxy  (  21, 3 ) =  263
greenupxy  (  21, 4 ) =  265
greenupxy  (  21, 5 ) =  265
greenupxy  (  21, 6 ) =  265
greenupxy  (  21, 7 ) =  237
greenupxy  (  21, 8 ) =  237
fallxy     (  21, 1 ) =   59
fallxy     (  21, 2 ) =   46
fallxy     (  21, 3 ) =   64
fallxy     (  21, 4 ) =   53
fallxy     (  21, 5 ) =   53
fallxy     (  21, 6 ) =   53
fallxy     (  21, 7 ) =   39
fallxy     (  21, 8 ) =   39
xyphendoy1 (  21, 1 ) =    2
xyphendoy1 (  21, 2 ) =    2
xyphendoy1 (  21, 3 ) =    2
xyphendoy1 (  21, 4 ) =    2
xyphendoy1 (  21, 5 ) =    2
xyphendoy1 (  21, 6 ) =    2
xyphendoy1 (  21, 7 ) =    2
xyphendoy1 (  21, 8 ) =    2
      
 ! xlat        -45.7500000000000     
greenupxy  (  20, 1 ) =  280
greenupxy  (  20, 2 ) =  276
greenupxy  (  20, 3 ) =  263
greenupxy  (  20, 4 ) =  265
greenupxy  (  20, 5 ) =  265
greenupxy  (  20, 6 ) =  265
greenupxy  (  20, 7 ) =  244
greenupxy  (  20, 8 ) =  244
fallxy     (  20, 1 ) =   59
fallxy     (  20, 2 ) =   50
fallxy     (  20, 3 ) =   56
fallxy     (  20, 4 ) =   53
fallxy     (  20, 5 ) =   53
fallxy     (  20, 6 ) =   53
fallxy     (  20, 7 ) =   39
fallxy     (  20, 8 ) =   39
xyphendoy1 (  20, 1 ) =    2
xyphendoy1 (  20, 2 ) =    2
xyphendoy1 (  20, 3 ) =    2
xyphendoy1 (  20, 4 ) =    2
xyphendoy1 (  20, 5 ) =    2
xyphendoy1 (  20, 6 ) =    2
xyphendoy1 (  20, 7 ) =    2
xyphendoy1 (  20, 8 ) =    2
      
 ! xlat        -46.2500000000000     
greenupxy  (  19, 1 ) =  280
greenupxy  (  19, 2 ) =  276
greenupxy  (  19, 3 ) =  262
greenupxy  (  19, 4 ) =  266
greenupxy  (  19, 5 ) =  266
greenupxy  (  19, 6 ) =  266
greenupxy  (  19, 7 ) =  251
greenupxy  (  19, 8 ) =  251
fallxy     (  19, 1 ) =   59
fallxy     (  19, 2 ) =   52
fallxy     (  19, 3 ) =   44
fallxy     (  19, 4 ) =   53
fallxy     (  19, 5 ) =   53
fallxy     (  19, 6 ) =   53
fallxy     (  19, 7 ) =   39
fallxy     (  19, 8 ) =   39
xyphendoy1 (  19, 1 ) =    2
xyphendoy1 (  19, 2 ) =    2
xyphendoy1 (  19, 3 ) =    2
xyphendoy1 (  19, 4 ) =    2
xyphendoy1 (  19, 5 ) =    2
xyphendoy1 (  19, 6 ) =    2
xyphendoy1 (  19, 7 ) =    2
xyphendoy1 (  19, 8 ) =    2
      
 ! xlat        -46.7500000000000     
greenupxy  (  18, 1 ) =  280
greenupxy  (  18, 2 ) =  277
greenupxy  (  18, 3 ) =  262
greenupxy  (  18, 4 ) =  266
greenupxy  (  18, 5 ) =  266
greenupxy  (  18, 6 ) =  266
greenupxy  (  18, 7 ) =  253
greenupxy  (  18, 8 ) =  253
fallxy     (  18, 1 ) =   59
fallxy     (  18, 2 ) =   51
fallxy     (  18, 3 ) =   41
fallxy     (  18, 4 ) =   50
fallxy     (  18, 5 ) =   50
fallxy     (  18, 6 ) =   50
fallxy     (  18, 7 ) =   41
fallxy     (  18, 8 ) =   41
xyphendoy1 (  18, 1 ) =    2
xyphendoy1 (  18, 2 ) =    2
xyphendoy1 (  18, 3 ) =    2
xyphendoy1 (  18, 4 ) =    2
xyphendoy1 (  18, 5 ) =    2
xyphendoy1 (  18, 6 ) =    2
xyphendoy1 (  18, 7 ) =    2
xyphendoy1 (  18, 8 ) =    2
      
 ! xlat        -47.2500000000000     
greenupxy  (  17, 1 ) =  280
greenupxy  (  17, 2 ) =  279
greenupxy  (  17, 3 ) =  263
greenupxy  (  17, 4 ) =  266
greenupxy  (  17, 5 ) =  266
greenupxy  (  17, 6 ) =  266
greenupxy  (  17, 7 ) =  260
greenupxy  (  17, 8 ) =  260
fallxy     (  17, 1 ) =   59
fallxy     (  17, 2 ) =   51
fallxy     (  17, 3 ) =   38
fallxy     (  17, 4 ) =   48
fallxy     (  17, 5 ) =   48
fallxy     (  17, 6 ) =   48
fallxy     (  17, 7 ) =   44
fallxy     (  17, 8 ) =   44
xyphendoy1 (  17, 1 ) =    2
xyphendoy1 (  17, 2 ) =    2
xyphendoy1 (  17, 3 ) =    2
xyphendoy1 (  17, 4 ) =    2
xyphendoy1 (  17, 5 ) =    2
xyphendoy1 (  17, 6 ) =    2
xyphendoy1 (  17, 7 ) =    2
xyphendoy1 (  17, 8 ) =    2
      
 ! xlat        -47.7500000000000     
greenupxy  (  16, 1 ) =  280
greenupxy  (  16, 2 ) =  280
greenupxy  (  16, 3 ) =  265
greenupxy  (  16, 4 ) =  268
greenupxy  (  16, 5 ) =  268
greenupxy  (  16, 6 ) =  268
greenupxy  (  16, 7 ) =  267
greenupxy  (  16, 8 ) =  267
fallxy     (  16, 1 ) =   59
fallxy     (  16, 2 ) =   52
fallxy     (  16, 3 ) =   35
fallxy     (  16, 4 ) =   45
fallxy     (  16, 5 ) =   45
fallxy     (  16, 6 ) =   45
fallxy     (  16, 7 ) =   43
fallxy     (  16, 8 ) =   43
xyphendoy1 (  16, 1 ) =    2
xyphendoy1 (  16, 2 ) =    2
xyphendoy1 (  16, 3 ) =    2
xyphendoy1 (  16, 4 ) =    2
xyphendoy1 (  16, 5 ) =    2
xyphendoy1 (  16, 6 ) =    2
xyphendoy1 (  16, 7 ) =    2
xyphendoy1 (  16, 8 ) =    2
      
 ! xlat        -48.2500000000000     
greenupxy  (  15, 1 ) =  280
greenupxy  (  15, 2 ) =  282
greenupxy  (  15, 3 ) =  267
greenupxy  (  15, 4 ) =  269
greenupxy  (  15, 5 ) =  269
greenupxy  (  15, 6 ) =  269
greenupxy  (  15, 7 ) =  269
greenupxy  (  15, 8 ) =  269
fallxy     (  15, 1 ) =   59
fallxy     (  15, 2 ) =   51
fallxy     (  15, 3 ) =   34
fallxy     (  15, 4 ) =   46
fallxy     (  15, 5 ) =   46
fallxy     (  15, 6 ) =   46
fallxy     (  15, 7 ) =   43
fallxy     (  15, 8 ) =   43
xyphendoy1 (  15, 1 ) =    2
xyphendoy1 (  15, 2 ) =    2
xyphendoy1 (  15, 3 ) =    2
xyphendoy1 (  15, 4 ) =    2
xyphendoy1 (  15, 5 ) =    2
xyphendoy1 (  15, 6 ) =    2
xyphendoy1 (  15, 7 ) =    2
xyphendoy1 (  15, 8 ) =    2
      
 ! xlat        -48.7500000000000     
greenupxy  (  14, 1 ) =  280
greenupxy  (  14, 2 ) =  282
greenupxy  (  14, 3 ) =  268
greenupxy  (  14, 4 ) =  271
greenupxy  (  14, 5 ) =  271
greenupxy  (  14, 6 ) =  271
greenupxy  (  14, 7 ) =  271
greenupxy  (  14, 8 ) =  271
fallxy     (  14, 1 ) =   59
fallxy     (  14, 2 ) =   51
fallxy     (  14, 3 ) =   30
fallxy     (  14, 4 ) =   44
fallxy     (  14, 5 ) =   44
fallxy     (  14, 6 ) =   44
fallxy     (  14, 7 ) =   46
fallxy     (  14, 8 ) =   46
xyphendoy1 (  14, 1 ) =    2
xyphendoy1 (  14, 2 ) =    2
xyphendoy1 (  14, 3 ) =    2
xyphendoy1 (  14, 4 ) =    2
xyphendoy1 (  14, 5 ) =    2
xyphendoy1 (  14, 6 ) =    2
xyphendoy1 (  14, 7 ) =    2
xyphendoy1 (  14, 8 ) =    2
      
 ! xlat        -49.2500000000000     
greenupxy  (  13, 1 ) =  280
greenupxy  (  13, 2 ) =  283
greenupxy  (  13, 3 ) =  269
greenupxy  (  13, 4 ) =  272
greenupxy  (  13, 5 ) =  272
greenupxy  (  13, 6 ) =  272
greenupxy  (  13, 7 ) =  273
greenupxy  (  13, 8 ) =  273
fallxy     (  13, 1 ) =   59
fallxy     (  13, 2 ) =   51
fallxy     (  13, 3 ) =   28
fallxy     (  13, 4 ) =   42
fallxy     (  13, 5 ) =   42
fallxy     (  13, 6 ) =   42
fallxy     (  13, 7 ) =   47
fallxy     (  13, 8 ) =   47
xyphendoy1 (  13, 1 ) =    2
xyphendoy1 (  13, 2 ) =    2
xyphendoy1 (  13, 3 ) =    2
xyphendoy1 (  13, 4 ) =    2
xyphendoy1 (  13, 5 ) =    2
xyphendoy1 (  13, 6 ) =    2
xyphendoy1 (  13, 7 ) =    2
xyphendoy1 (  13, 8 ) =    2
      
 ! xlat        -49.7500000000000     
greenupxy  (  12, 1 ) =  280
greenupxy  (  12, 2 ) =  284
greenupxy  (  12, 3 ) =  270
greenupxy  (  12, 4 ) =  274
greenupxy  (  12, 5 ) =  274
greenupxy  (  12, 6 ) =  274
greenupxy  (  12, 7 ) =  275
greenupxy  (  12, 8 ) =  275
fallxy     (  12, 1 ) =   59
fallxy     (  12, 2 ) =   51
fallxy     (  12, 3 ) =   32
fallxy     (  12, 4 ) =   41
fallxy     (  12, 5 ) =   41
fallxy     (  12, 6 ) =   41
fallxy     (  12, 7 ) =   48
fallxy     (  12, 8 ) =   48
xyphendoy1 (  12, 1 ) =    2
xyphendoy1 (  12, 2 ) =    2
xyphendoy1 (  12, 3 ) =    2
xyphendoy1 (  12, 4 ) =    2
xyphendoy1 (  12, 5 ) =    2
xyphendoy1 (  12, 6 ) =    2
xyphendoy1 (  12, 7 ) =    2
xyphendoy1 (  12, 8 ) =    2
      
 ! xlat        -50.2500000000000     
greenupxy  (  11, 1 ) =  280
greenupxy  (  11, 2 ) =  286
greenupxy  (  11, 3 ) =  272
greenupxy  (  11, 4 ) =  276
greenupxy  (  11, 5 ) =  276
greenupxy  (  11, 6 ) =  276
greenupxy  (  11, 7 ) =  277
greenupxy  (  11, 8 ) =  277
fallxy     (  11, 1 ) =   59
fallxy     (  11, 2 ) =   52
fallxy     (  11, 3 ) =   32
fallxy     (  11, 4 ) =   41
fallxy     (  11, 5 ) =   41
fallxy     (  11, 6 ) =   41
fallxy     (  11, 7 ) =   50
fallxy     (  11, 8 ) =   50
xyphendoy1 (  11, 1 ) =    2
xyphendoy1 (  11, 2 ) =    2
xyphendoy1 (  11, 3 ) =    2
xyphendoy1 (  11, 4 ) =    2
xyphendoy1 (  11, 5 ) =    2
xyphendoy1 (  11, 6 ) =    2
xyphendoy1 (  11, 7 ) =    2
xyphendoy1 (  11, 8 ) =    2
      
 ! xlat        -50.7500000000000     
greenupxy  (  10, 1 ) =  280
greenupxy  (  10, 2 ) =  288
greenupxy  (  10, 3 ) =  273
greenupxy  (  10, 4 ) =  278
greenupxy  (  10, 5 ) =  278
greenupxy  (  10, 6 ) =  278
greenupxy  (  10, 7 ) =  278
greenupxy  (  10, 8 ) =  278
fallxy     (  10, 1 ) =   59
fallxy     (  10, 2 ) =   53
fallxy     (  10, 3 ) =   31
fallxy     (  10, 4 ) =   42
fallxy     (  10, 5 ) =   42
fallxy     (  10, 6 ) =   42
fallxy     (  10, 7 ) =   50
fallxy     (  10, 8 ) =   50
xyphendoy1 (  10, 1 ) =    2
xyphendoy1 (  10, 2 ) =    2
xyphendoy1 (  10, 3 ) =    2
xyphendoy1 (  10, 4 ) =    2
xyphendoy1 (  10, 5 ) =    2
xyphendoy1 (  10, 6 ) =    2
xyphendoy1 (  10, 7 ) =    2
xyphendoy1 (  10, 8 ) =    2
      
 ! xlat        -51.2500000000000     
greenupxy  (   9, 1 ) =  280
greenupxy  (   9, 2 ) =  289
greenupxy  (   9, 3 ) =  275
greenupxy  (   9, 4 ) =  279
greenupxy  (   9, 5 ) =  279
greenupxy  (   9, 6 ) =  279
greenupxy  (   9, 7 ) =  278
greenupxy  (   9, 8 ) =  278
fallxy     (   9, 1 ) =   59
fallxy     (   9, 2 ) =   53
fallxy     (   9, 3 ) =   34
fallxy     (   9, 4 ) =   42
fallxy     (   9, 5 ) =   42
fallxy     (   9, 6 ) =   42
fallxy     (   9, 7 ) =   50
fallxy     (   9, 8 ) =   50
xyphendoy1 (   9, 1 ) =    2
xyphendoy1 (   9, 2 ) =    2
xyphendoy1 (   9, 3 ) =    2
xyphendoy1 (   9, 4 ) =    2
xyphendoy1 (   9, 5 ) =    2
xyphendoy1 (   9, 6 ) =    2
xyphendoy1 (   9, 7 ) =    2
xyphendoy1 (   9, 8 ) =    2
      
 ! xlat        -51.7500000000000     
greenupxy  (   8, 1 ) =  280
greenupxy  (   8, 2 ) =  290
greenupxy  (   8, 3 ) =  276
greenupxy  (   8, 4 ) =  280
greenupxy  (   8, 5 ) =  280
greenupxy  (   8, 6 ) =  280
greenupxy  (   8, 7 ) =  278
greenupxy  (   8, 8 ) =  278
fallxy     (   8, 1 ) =   59
fallxy     (   8, 2 ) =   54
fallxy     (   8, 3 ) =   34
fallxy     (   8, 4 ) =   42
fallxy     (   8, 5 ) =   42
fallxy     (   8, 6 ) =   42
fallxy     (   8, 7 ) =   50
fallxy     (   8, 8 ) =   50
xyphendoy1 (   8, 1 ) =    2
xyphendoy1 (   8, 2 ) =    2
xyphendoy1 (   8, 3 ) =    2
xyphendoy1 (   8, 4 ) =    2
xyphendoy1 (   8, 5 ) =    2
xyphendoy1 (   8, 6 ) =    2
xyphendoy1 (   8, 7 ) =    2
xyphendoy1 (   8, 8 ) =    2
      
 ! xlat        -52.2500000000000     
greenupxy  (   7, 1 ) =  280
greenupxy  (   7, 2 ) =  290
greenupxy  (   7, 3 ) =  276
greenupxy  (   7, 4 ) =  280
greenupxy  (   7, 5 ) =  280
greenupxy  (   7, 6 ) =  280
greenupxy  (   7, 7 ) =  278
greenupxy  (   7, 8 ) =  278
fallxy     (   7, 1 ) =   59
fallxy     (   7, 2 ) =   54
fallxy     (   7, 3 ) =   35
fallxy     (   7, 4 ) =   42
fallxy     (   7, 5 ) =   42
fallxy     (   7, 6 ) =   42
fallxy     (   7, 7 ) =   50
fallxy     (   7, 8 ) =   50
xyphendoy1 (   7, 1 ) =    2
xyphendoy1 (   7, 2 ) =    2
xyphendoy1 (   7, 3 ) =    2
xyphendoy1 (   7, 4 ) =    2
xyphendoy1 (   7, 5 ) =    2
xyphendoy1 (   7, 6 ) =    2
xyphendoy1 (   7, 7 ) =    2
xyphendoy1 (   7, 8 ) =    2
      
 ! xlat        -52.7500000000000     
greenupxy  (   6, 1 ) =  280
greenupxy  (   6, 2 ) =  290
greenupxy  (   6, 3 ) =  276
greenupxy  (   6, 4 ) =  280
greenupxy  (   6, 5 ) =  280
greenupxy  (   6, 6 ) =  280
greenupxy  (   6, 7 ) =  278
greenupxy  (   6, 8 ) =  278
fallxy     (   6, 1 ) =   59
fallxy     (   6, 2 ) =   54
fallxy     (   6, 3 ) =   35
fallxy     (   6, 4 ) =   42
fallxy     (   6, 5 ) =   42
fallxy     (   6, 6 ) =   42
fallxy     (   6, 7 ) =   50
fallxy     (   6, 8 ) =   50
xyphendoy1 (   6, 1 ) =    2
xyphendoy1 (   6, 2 ) =    2
xyphendoy1 (   6, 3 ) =    2
xyphendoy1 (   6, 4 ) =    2
xyphendoy1 (   6, 5 ) =    2
xyphendoy1 (   6, 6 ) =    2
xyphendoy1 (   6, 7 ) =    2
xyphendoy1 (   6, 8 ) =    2
      
 ! xlat        -53.2500000000000     
greenupxy  (   5, 1 ) =  280
greenupxy  (   5, 2 ) =  290
greenupxy  (   5, 3 ) =  276
greenupxy  (   5, 4 ) =  280
greenupxy  (   5, 5 ) =  280
greenupxy  (   5, 6 ) =  280
greenupxy  (   5, 7 ) =  278
greenupxy  (   5, 8 ) =  278
fallxy     (   5, 1 ) =   59
fallxy     (   5, 2 ) =   54
fallxy     (   5, 3 ) =   35
fallxy     (   5, 4 ) =   42
fallxy     (   5, 5 ) =   42
fallxy     (   5, 6 ) =   42
fallxy     (   5, 7 ) =   50
fallxy     (   5, 8 ) =   50
xyphendoy1 (   5, 1 ) =    2
xyphendoy1 (   5, 2 ) =    2
xyphendoy1 (   5, 3 ) =    2
xyphendoy1 (   5, 4 ) =    2
xyphendoy1 (   5, 5 ) =    2
xyphendoy1 (   5, 6 ) =    2
xyphendoy1 (   5, 7 ) =    2
xyphendoy1 (   5, 8 ) =    2
      
 ! xlat        -53.7500000000000     
greenupxy  (   4, 1 ) =  280
greenupxy  (   4, 2 ) =  290
greenupxy  (   4, 3 ) =  276
greenupxy  (   4, 4 ) =  280
greenupxy  (   4, 5 ) =  280
greenupxy  (   4, 6 ) =  280
greenupxy  (   4, 7 ) =  278
greenupxy  (   4, 8 ) =  278
fallxy     (   4, 1 ) =   59
fallxy     (   4, 2 ) =   54
fallxy     (   4, 3 ) =   35
fallxy     (   4, 4 ) =   42
fallxy     (   4, 5 ) =   42
fallxy     (   4, 6 ) =   42
fallxy     (   4, 7 ) =   50
fallxy     (   4, 8 ) =   50
xyphendoy1 (   4, 1 ) =    2
xyphendoy1 (   4, 2 ) =    2
xyphendoy1 (   4, 3 ) =    2
xyphendoy1 (   4, 4 ) =    2
xyphendoy1 (   4, 5 ) =    2
xyphendoy1 (   4, 6 ) =    2
xyphendoy1 (   4, 7 ) =    2
xyphendoy1 (   4, 8 ) =    2
      
 ! xlat        -54.2500000000000     
greenupxy  (   3, 1 ) =  280
greenupxy  (   3, 2 ) =  290
greenupxy  (   3, 3 ) =  276
greenupxy  (   3, 4 ) =  280
greenupxy  (   3, 5 ) =  280
greenupxy  (   3, 6 ) =  280
greenupxy  (   3, 7 ) =  278
greenupxy  (   3, 8 ) =  278
fallxy     (   3, 1 ) =   59
fallxy     (   3, 2 ) =   54
fallxy     (   3, 3 ) =   35
fallxy     (   3, 4 ) =   42
fallxy     (   3, 5 ) =   42
fallxy     (   3, 6 ) =   42
fallxy     (   3, 7 ) =   50
fallxy     (   3, 8 ) =   50
xyphendoy1 (   3, 1 ) =    2
xyphendoy1 (   3, 2 ) =    2
xyphendoy1 (   3, 3 ) =    2
xyphendoy1 (   3, 4 ) =    2
xyphendoy1 (   3, 5 ) =    2
xyphendoy1 (   3, 6 ) =    2
xyphendoy1 (   3, 7 ) =    2
xyphendoy1 (   3, 8 ) =    2
      
 ! xlat        -54.7500000000000     
greenupxy  (   2, 1 ) =  280
greenupxy  (   2, 2 ) =  290
greenupxy  (   2, 3 ) =  276
greenupxy  (   2, 4 ) =  280
greenupxy  (   2, 5 ) =  280
greenupxy  (   2, 6 ) =  280
greenupxy  (   2, 7 ) =  278
greenupxy  (   2, 8 ) =  278
fallxy     (   2, 1 ) =   59
fallxy     (   2, 2 ) =   54
fallxy     (   2, 3 ) =   35
fallxy     (   2, 4 ) =   42
fallxy     (   2, 5 ) =   42
fallxy     (   2, 6 ) =   42
fallxy     (   2, 7 ) =   50
fallxy     (   2, 8 ) =   50
xyphendoy1 (   2, 1 ) =    2
xyphendoy1 (   2, 2 ) =    2
xyphendoy1 (   2, 3 ) =    2
xyphendoy1 (   2, 4 ) =    2
xyphendoy1 (   2, 5 ) =    2
xyphendoy1 (   2, 6 ) =    2
xyphendoy1 (   2, 7 ) =    2
xyphendoy1 (   2, 8 ) =    2
      
 ! xlat        -55.2500000000000     
greenupxy  (   1, 1 ) =  280
greenupxy  (   1, 2 ) =  290
greenupxy  (   1, 3 ) =  276
greenupxy  (   1, 4 ) =  280
greenupxy  (   1, 5 ) =  280
greenupxy  (   1, 6 ) =  280
greenupxy  (   1, 7 ) =  278
greenupxy  (   1, 8 ) =  278
fallxy     (   1, 1 ) =   59
fallxy     (   1, 2 ) =   54
fallxy     (   1, 3 ) =   35
fallxy     (   1, 4 ) =   42
fallxy     (   1, 5 ) =   42
fallxy     (   1, 6 ) =   42
fallxy     (   1, 7 ) =   50
fallxy     (   1, 8 ) =   50
xyphendoy1 (   1, 1 ) =    2
xyphendoy1 (   1, 2 ) =    2
xyphendoy1 (   1, 3 ) =    2
xyphendoy1 (   1, 4 ) =    2
xyphendoy1 (   1, 5 ) =    2
xyphendoy1 (   1, 6 ) =    2
xyphendoy1 (   1, 7 ) =    2
xyphendoy1 (   1, 8 ) =    2
      
!jhan: alternative to reading file
DO nx=1,nphen
   greenup(ilat,ivtx(nx)) = greenupxy (ilat, nx)
   fall(ilat,ivtx(nx))    = fallxy    (ilat, nx)
   phendoy1(ilat,ivtx(nx))= xyphendoy1(ilat, nx)
ENDDO

    DO np=1,mp
       ilat=(casamet%lat(np)+55.25)/0.5+1
       ilat= MIN(271,MAX(1,ilat))
       phen%phase(np) = phendoy1(ilat,veg%iveg(np))
       phen%doyphase(np,1) = greenup(ilat,veg%iveg(np)) ! DOY for greenup
       phen%doyphase(np,2) = phen%doyphase(np,1) +14    ! DOY for steady LAI
       phen%doyphase(np,3) = fall(ilat,veg%iveg(np))    ! DOY for leaf senescence
       phen%doyphase(np,4) = phen%doyphase(np,3) +14    ! DOY for minimal LAI season
       IF (phen%doyphase(np,2) > 365) phen%doyphase(np,2)=phen%doyphase(np,2)-365
       IF (phen%doyphase(np,4) > 365) phen%doyphase(np,4)=phen%doyphase(np,4)-365

    ENDDO

  END SUBROUTINE casa_readphen

  SUBROUTINE casa_init(casabiome,casamet,casaflux,casapool,casabal,veg,phen)
    ! mst not used (BP sep2010)
    ! ! for first time reading file *_1220.csv  (BP may2010)
    !SUBROUTINE casa_init(mst,casapool,casabal,veg)
    ! !SUBROUTINE casa_init(mst,casapool,casabal)
    ! ! end addition (BP may2010)
    !  initialize some values in phenology parameters and leaf growth phase
    USE casadimension
    USE casaparm

USE phenology_type_mod,    ONLY: phenology_type 
USE casa_biome_type_mod,   ONLY: casa_biome   => casa_biome_type
USE casa_pool_type_mod,    ONLY: casa_pool    => casa_pool_type
USE casa_flux_type_mod,    ONLY: casa_flux    => casa_flux_type
USE casa_met_type_mod,     ONLY: casa_met     => casa_met_type
USE casa_balance_type_mod, ONLY: casa_balance => casa_bal_type             

    ! for first time reading file *_1220.csv  (BP may2010)
    USE cable_def_types_mod
!jh: we get lat/lon from UM anyway and I want to get rid of this dependency
!jh!    USE cable_io_vars_module, ONLY: patch
    USE cable_common_module, ONLY: cable_user
#ifndef UM_CBL 
USE casa_offline_inout_module, ONLY : READ_CASA_RESTART_NC
#endif

    ! end addition (BP may2010)
    IMPLICIT NONE
    !  INTEGER,        INTENT(IN)    :: mst
    TYPE (casa_biome),   INTENT(IN)    :: casabiome
    TYPE (casa_met),     INTENT(INOUT) :: casamet
    TYPE (casa_flux),    INTENT(INOUT) :: casaflux
    TYPE (casa_pool),    INTENT(INOUT) :: casapool
    TYPE (casa_balance), INTENT(INOUT) :: casabal
    ! for first time reading file *_1220.csv  (BP may2010)
    TYPE (veg_parameter_type), INTENT(IN) :: veg
    TYPE (phenology_type),   INTENT(INOUT) :: phen
    REAL(r_2) :: clabile,cplant(3),clitter(3),csoil(3)
    REAL(r_2) :: nplant(3),nlitter(3),nsoil(3),nsoilmin,pplant(3)
    REAL(r_2) :: plitter(3),psoil(3),psoillab,psoilsorb,psoilocc
    ! end addition (BP may2010)

    ! local variables
    INTEGER   :: np,npt,npz
    INTEGER   :: nyearz,ivtz,istz,isoz
    REAL(r_2) :: latz,lonz,areacellz,glaiz,slaz
    LOGICAL   :: EXRST


    IF (.NOT.cable_user%casa_fromzero) THEN
       PRINT *, 'initial pool from restart file'
    ENDIF
    PRINT *, 'icycle,initcasa,mp ', icycle,initcasa,mp
    !phen%phase = 2

    !CLN initialise all ! THIS NEEDS FIXING because of e.g. ICE-WATER
    casaflux%Cgpp         = 0.
    casaflux%Cnpp         = 0.
    casaflux%Crp          = 0.
    casaflux%Crgplant     = 0.
    ! casaflux%Nminfix      = 0.
    casaflux%Nminuptake   = 0.
    casaflux%Plabuptake   = 0.
    casaflux%Clabloss     = 0.
    casaflux%fracClabile  = 0.
    casaflux%stemnpp      = 0.
    casaflux%frac_sapwood = 0.
    casaflux%sapwood_area = 0.
    casaflux%FluxCtohwp = 0.
    casaflux%FluxCtoClear = 0.
    casaflux%fracCalloc   = 0.
    casaflux%fracNalloc   = 0.
    casaflux%fracPalloc   = 0.
    casaflux%Crmplant     = 0.
    casaflux%kplant       = 0.

    casaflux%fromPtoL     = 0.

    casaflux%Cnep         = 0.
    casaflux%Crsoil       = 0.
    casapool%dClabiledt = 0.0
    !casaflux%Nmindep      =  casaflux%Nmindep /2.0
    !casaflux%Nmindep      = 0.
    casaflux%Nminloss     = 0.
    casaflux%Nminleach    = 0.
    casaflux%Nupland      = 0.
    casaflux%Nlittermin   = 0.
    casaflux%Nsmin        = 0.
    casaflux%Nsimm        = 0.
    casaflux%Nsnet        = 0.
    !casaflux%fNminloss    = 0.
    !casaflux%fNminleach   = 0.
    !casaflux%Pdep         = 0.
    !casaflux%Pwea         = 0.
    casaflux%Pleach       = 0.
    casaflux%Ploss        = 0.
    casaflux%Pupland      = 0.
    casaflux%Plittermin   = 0.
    casaflux%Psmin        = 0.
    casaflux%Psimm        = 0.
    casaflux%Psnet        = 0.
    !  casaflux%fPleach      = 0. !vh ! this should be a parameter, not a flux variable
    casaflux%kplab        = 0.
    casaflux%kpsorb       = 0.
    casaflux%kpocc        = 0.
    !  casaflux%Psorbmax     = 0. !vh ! this should be a paramter, not a flux variable

    casaflux%klitter      = 0.
    casaflux%ksoil        = 0.
    casaflux%fromLtoS     = 0.
    casaflux%fromStoS     = 0.
    casaflux%fromLtoCO2   = 0.
    casaflux%fromStoCO2   = 0.
    casaflux%FluxCtolitter= 0.
    casaflux%FluxNtolitter= 0.
    casaflux%FluxPtolitter= 0.
    casaflux%FluxCtosoil  = 0.
    casaflux%FluxNtosoil  = 0.
    casaflux%FluxPtosoil  = 0.
    casaflux%FluxCtoCO2   = 0.

    casaflux%Cplant_turnover = 0.
!Ticket #204 - rm phen% clobbing here AND incorrectly so anyway

    IF (initcasa==1) THEN
       IF (.NOT.cable_user%casa_fromzero) THEN
#ifndef UM_CBL 
          CALL READ_CASA_RESTART_NC (  casamet, casapool, casaflux, phen )
#endif
       ELSE
          WRITE(*,*)'casa_init: not using restart file!'
          WRITE(*,*)'Using input from readbiome.!'
          WRITE(*,*) 'initialising frac_sapwood=1 and sapwood_area = 0)'
          casaflux%frac_sapwood(:) = 1.0
          casaflux%sapwood_area(:) = 0.0
       ENDIF
    ENDIF
    WHERE(casamet%lnonwood==1) casapool%cplant(:,WOOD) = 0.0
    WHERE(casamet%lnonwood==1) casapool%nplant(:,WOOD) = 0.0
    WHERE(casamet%lnonwood==1) casapool%pplant(:,WOOD) = 0.0

!$ENDIF
    !92 format(5(i6,2x),5(f18.6,3x),2(i6,',',2x),',',2x,100(f18.6,3x))
92  FORMAT(5(i6,',',2x),5(f18.6,',',2x),2(i6,',',2x),',',2x,100(f18.6,',',2x))

!jh: we get lat/lon from UM anyway and I want to get rid of this dependency
!jh!    IF(initcasa==0) THEN
!jh!       nyearz = 1
!jh!       DO npt=1,mp
!jh!          casamet%lon(npt) = patch(npt)%longitude
!jh!          casamet%lat(npt) = patch(npt)%latitude
!jh!       ENDDO
!jh!    ENDIF
!jh!
    ! reset labile C pool,comment out by Q.Zhang 10/09/2011
    !  casapool%clabile    = 0.0
    ! check pool sizes
    casapool%cplant     = MAX(0.0,casapool%cplant)
    casapool%clitter    = MAX(0.0,casapool%clitter)
    casapool%csoil      = MAX(0.0,casapool%csoil)
    casabal%cplantlast  = casapool%cplant
    casabal%clitterlast = casapool%clitter
    casabal%csoillast   = casapool%csoil
    casabal%clabilelast = casapool%clabile
    casabal%sumcbal     = 0.0
    casabal%FCgppyear=0.0;casabal%FCrpyear=0.0
    casabal%FCnppyear=0;casabal%FCrsyear=0.0;casabal%FCneeyear=0.0
    !vh !
    WHERE(casamet%lnonwood==1) casapool%cplant(:,WOOD) = 0.0
    IF (icycle==1) THEN
       casapool%Nplant(:,:) = casapool%cplant(:,:) * casapool%ratioNCplant(:,:)
       casapool%Nsoil(:,:)  = casapool%ratioNCsoil(:,:) * casapool%Csoil(:,:)
       casapool%Psoil(:,:)  = casapool%Nsoil(:,:)/casapool%ratioNPsoil(:,:)
       casapool%Nsoilmin(:) = 2.5
    ENDIF

    IF (icycle >=1) THEN
       casapool%nplant     = MAX(1.e-6,casapool%nplant)
       casapool%nlitter    = MAX(1.e-6,casapool%nlitter)
       casapool%nsoil      = MAX(1.e-6,casapool%nsoil)
       casapool%nsoilmin   = MAX(1.e-6,casapool%nsoilmin)
       casabal%nplantlast  = casapool%nplant
       casabal%nlitterlast = casapool%nlitter
       casabal%nsoillast   = casapool%nsoil
       casabal%nsoilminlast= casapool%nsoilmin
       casabal%sumnbal     = 0.0
       casabal%FNdepyear=0.0;casabal%FNfixyear=0.0;casabal%FNsnetyear=0.0
       casabal%FNupyear=0.0;casabal%FNleachyear=0.0;casabal%FNlossyear=0.0
       !vh !
       WHERE(casamet%lnonwood==1) casapool%nplant(:,WOOD) = 0.0
    ENDIF

    IF (icycle >=1) THEN
       casapool%pplant       = MAX(1.0e-7,casapool%pplant)
       casapool%plitter      = MAX(1.0e-7,casapool%plitter)
       casapool%psoil        = MAX(1.0e-7,casapool%psoil)
       casapool%Psoillab     = MAX(1.0e-7,casapool%psoillab)  ! was 2.0, changed according to  YP
       casapool%psoilsorb    = MAX(1.0e-7,casapool%psoilsorb) ! was 10.0, -
       casapool%psoilocc     = MAX(1.0e-7,casapool%psoilocc)  ! was 50.0, -
       casabal%pplantlast    = casapool%pplant
       casabal%plitterlast   = casapool%plitter
       casabal%psoillast     = casapool%psoil
       casabal%psoillablast  = casapool%psoillab
       casabal%psoilsorblast = casapool%psoilsorb
       casabal%psoilocclast  = casapool%psoilocc
       casabal%sumpbal       = 0.0
       casabal%FPweayear=0.0;casabal%FPdustyear=0.0; casabal%FPsnetyear=0.0
       casabal%FPupyear=0.0;casabal%FPleachyear=0.0;casabal%FPlossyear=0.0
       !vh !
       WHERE(casamet%lnonwood==1) casapool%pplant(:,WOOD) = 0.0
    ENDIF
    
    casapool%cwoodprod=0.0; casapool%nwoodprod=0.0;casapool%pwoodprod=0.0

  END SUBROUTINE casa_init


  SUBROUTINE casa_poolout(ktau,veg,soil,casabiome,casapool,casaflux,casamet, &
       casabal,phen)
    USE cable_def_types_mod
    USE casadimension
    USE casaparm
    USE casavariable
    USE cable_common_module, ONLY: cable_user
USE phenology_type_mod,    ONLY: phenology_type 
USE casa_pool_type_mod,    ONLY: casa_pool    => casa_pool_type
USE casa_flux_type_mod,    ONLY: casa_flux    => casa_flux_type
USE casa_biome_type_mod,   ONLY: casa_biome   => casa_biome_type
USE casa_met_type_mod,     ONLY: casa_met     => casa_met_type
USE casa_balance_type_mod, ONLY: casa_balance => casa_bal_type             

    IMPLICIT NONE

    INTEGER,               INTENT(IN)    :: ktau
    TYPE (veg_parameter_type),  INTENT(INOUT) :: veg  ! vegetation parameters
    TYPE (soil_parameter_type), INTENT(INOUT) :: soil ! soil parameters
    TYPE (casa_biome),          INTENT(INOUT) :: casabiome
    TYPE (casa_pool),           INTENT(INOUT) :: casapool
    TYPE (casa_flux),           INTENT(INOUT) :: casaflux
    TYPE (casa_met),            INTENT(INOUT) :: casamet
    TYPE (casa_balance),        INTENT(INOUT) :: casabal
    TYPE (phenology_type),       INTENT(INOUT) :: phen

    ! local variables
    REAL(r_2), DIMENSION(mso) :: Psorder,pweasoil,xpsoil50
    REAL(r_2), DIMENSION(mso) :: fracPlab,fracPsorb,fracPocc,fracPorg
    REAL(r_2), DIMENSION(mp)  :: totpsoil
    INTEGER  npt,nout,nso

    ! Soiltype     soilnumber soil P(g P/m2)
    ! Alfisol     1       61.3
    ! Andisol     2       103.9
    ! Aridisol    3       92.8
    ! Entisol     4       136.9
    ! Gellisol    5       98.2
    ! Histosol    6       107.6
    ! Inceptisol  7       84.1
    ! Mollisol    8       110.1
    ! Oxisol      9       35.4
    ! Spodosol    10      41.0
    ! Ultisol     11      51.5
    ! Vertisol    12      190.6
    DATA psorder/61.3,103.9,92.8,136.9,98.2,107.6,84.1,110.1,35.4,41.0,51.5,190.6/
    DATA pweasoil/0.05,0.04,0.03,0.02,0.01,0.009,0.008,0.007,0.006,0.005,0.004,0.003/
    DATA fracpLab/0.08,0.08,0.10,0.02,0.08,0.08,0.08,0.06,0.02,0.05,0.09,0.05/
    DATA fracPsorb/0.32,0.37,0.57,0.67,0.37,0.37,0.37,0.32,0.24,0.22,0.21,0.38/
    DATA fracPocc/0.36,0.38,0.25,0.26,0.38,0.38,0.38,0.44,0.38,0.38,0.37,0.45/
    DATA fracPorg/0.25,0.17,0.08,0.05,0.17,0.17,0.17,0.18,0.36,0.35,0.34,0.12/
    DATA xpsoil50/7.6,4.1,4.2,3.4,4.1,4.1,4.8,4.1,6.9,6.9,6.9,1.7/
    !
    ! estimated based on Yang, Post and Jain (2013)
    !   Soiltype     soilnumber soil P(g P/m2  top 50 cm)
    !   Alfisol     1       400
    !   Andisol     2       426
    !   Aridisol    3       352
    !   Entisol     4       490
    !   Gellisol    5       403
    !   Histosol    6       441
    !   Inceptisol  7       501
    !   Mollisol    8       358
    !   Oxisol      9       96
    !   Spodosol    10      364
    !   Ultisol     11      272
    !   Vertisol    12      430
    !  DATA psorder/400.0,426.0,352.0,490.0,403.0,441.0,501.0,358.0,96.0,364.0,272.0,430.0/
    !  DATA pweasoil/0.05,0.04,0.03,0.02,0.01,0.009,0.008,0.007,0.006,0.005,0.004,0.003/
    !  DATA fracpLab/0.07,0.04,0.08,0.10,0.08,0.10,0.12,0.05,0.05,0.06,0.06,0.05/
    !  DATA fracPsorb/0.30,0.44,0.69,0.53,0.37,0.14,0.24,0.32,0.15,0.21,0.17,0.35/
    !  DATA fracPocc/0.38,0.22,0.18,0.22,0.38,0.42,0.23,0.44,0.60,0.30,0.51,0.48/
    !  DATA fracPorg/0.25,0.30,0.05,0.15,0.17,0.34,0.41,0.19,0.20,0.43,0.26,0.12/
    !  DATA xpsoil50/1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0/

    PRINT *, 'Within casa_poolout, mp = ', mp
    nout=103
    OPEN(nout,file=casafile%cnpepool)
    PRINT *, 'Opened file ', casafile%cnpepool

    casabal%sumcbal=MIN(9999.0,MAX(-9999.0,casabal%sumcbal))
    casabal%sumnbal=MIN(9999.0,MAX(-9999.0,casabal%sumnbal))
    casabal%sumpbal=MIN(9999.0,MAX(-9999.0,casabal%sumpbal))

    DO npt =1, mp
       nso = casamet%isorder(npt)
       totpsoil(npt) = psorder(nso) *xpsoil50(nso)
       IF(casamet%iveg2(npt)>0 ) THEN
          IF (icycle<2) THEN
             casapool%Nplant(npt,:) = casapool%ratioNCplant(npt,:)  &
                  * casapool%cplant(npt,:)
             casapool%Nlitter(npt,:)= casapool%ratioNClitter(npt,:) &
                  * casapool%clitter(npt,:)
             casapool%Nsoil(npt,:)  = casapool%ratioNCsoil(npt,:)   &
                  * casapool%Csoil(npt,:)
             casapool%nsoilmin(npt) = 2.0
             casabal%sumnbal(npt)   = 0.0
             IF(casamet%iveg2(npt)==grass) THEN
                casapool%nplant(npt,wood) = 0.0
                casapool%nlitter(npt,cwd) = 0.0
             ENDIF
          ENDIF

          IF (icycle<3) THEN
             casabal%sumpbal(npt)   = 0.0
             casapool%pplant(npt,:)  = casapool%Nplant(npt,:)/casapool%ratioNPplant(npt,:)
             casapool%plitter(npt,:) = casapool%Nlitter(npt,:)/(casapool%ratioNPlitter(npt,:)+1.0e-10)
             casapool%psoil(npt,:)   = casapool%Nsoil(npt,:)/casapool%ratioNPsoil(npt,:)
             casapool%psoillab(npt) = totpsoil(npt) *fracpLab(nso)
             casapool%psoilsorb(npt)= casaflux%psorbmax(npt) * casapool%psoillab(npt) &
                  /(casaflux%kmlabp(npt)+casapool%psoillab(npt))
             casapool%psoilocc(npt) = totpsoil(npt) *fracPocc(nso)
             IF(casamet%iveg2(npt)==grass) THEN
                casapool%pplant(npt,wood) = 0.0
                casapool%plitter(npt,cwd) = 0.0
             ENDIF
          ENDIF
       ELSE
          casapool%cplant(npt,:)=0.0; casapool%clitter(npt,:)=0.0; casapool%csoil(npt,:) = 0.0; casapool%clabile(npt) = 0.0
          casapool%nplant(npt,:)=0.0; casapool%nlitter(npt,:)=0.0; casapool%nsoil(npt,:) = 0.0; casapool%nsoilmin(npt) = 0.0
          casapool%pplant(npt,:)=0.0; casapool%plitter(npt,:)=0.0; casapool%psoil(npt,:) = 0.0
          casapool%psoillab(npt) = 0.0; casapool%psoilsorb(npt) = 0.0; casapool%psoilocc(npt) = 0.0
          casabal%sumcbal(npt) =0.0; casabal%sumnbal(npt) =0.0; casabal%sumpbal(npt) = 0.0
       ENDIF

       ! vh_js  !
       IF (cable_user%CALL_POP) THEN

          WRITE(nout,92) ktau,npt,veg%iveg(npt),soil%isoilm(npt) ,     &
               casamet%isorder(npt),casamet%lat(npt),casamet%lon(npt), &
               casamet%areacell(npt)*(1.0e-9),casamet%glai(npt),       &
               casabiome%sla(veg%iveg(npt)), phen%phase(npt), &
               phen%doyphase(npt,3), phen%phen(npt), phen%aphen(npt), &
               casapool%clabile(npt), &
               casapool%cplant(npt,:),casapool%clitter(npt,:),casapool%csoil(npt,:), &
               casaflux%frac_sapwood(npt), casaflux%sapwood_area(npt), &
               casapool%nplant(npt,:),casapool%nlitter(npt,:),casapool%nsoil(npt,:), &
               casapool%nsoilmin(npt),casapool%pplant(npt,:),          &
               casapool%plitter(npt,:), casapool%psoil(npt,:),         &
               casapool%psoillab(npt),casapool%psoilsorb(npt),casapool%psoilocc(npt), &
               casabal%sumcbal(npt),casabal%sumnbal(npt),casabal%sumpbal(npt)


       ELSE
          WRITE(nout,92) ktau,npt,veg%iveg(npt),soil%isoilm(npt),     &
               casamet%isorder(npt),casamet%lat(npt),casamet%lon(npt), &
               casamet%areacell(npt)*(1.0e-9),casamet%glai(npt),       &
               casabiome%sla(veg%iveg(npt)), phen%phase(npt), &
               phen%doyphase(npt,3), phen%phen(npt), phen%aphen(npt), &
               casapool%clabile(npt), &
               casapool%cplant(npt,:),casapool%clitter(npt,:),casapool%csoil(npt,:), &
               casapool%nplant(npt,:),casapool%nlitter(npt,:),casapool%nsoil(npt,:), &
               casapool%nsoilmin(npt),casapool%pplant(npt,:),          &
               casapool%plitter(npt,:), casapool%psoil(npt,:),         &
               casapool%psoillab(npt),casapool%psoilsorb(npt),casapool%psoilocc(npt), &
               casabal%sumcbal(npt),casabal%sumnbal(npt),casabal%sumpbal(npt)
       ENDIF


    ENDDO

    CLOSE(nout)

92  FORMAT(5(i6,',',2x),5(f18.6,',',2x),2(i6,',',2x),100(f18.6,',',2x))
  END SUBROUTINE casa_poolout

SUBROUTINE casa_fluxout(myear,veg,soil,casabal,casamet)

    USE cable_def_types_mod
    USE casadimension
    USE casaparm
USE casa_met_type_mod,     ONLY: casa_met     => casa_met_type
USE casa_balance_type_mod, ONLY: casa_balance => casa_bal_type             
USE phenology_type_mod,          ONLY: phenology_type
 
    IMPLICIT NONE

    TYPE (veg_parameter_type),  INTENT(INOUT) :: veg  ! vegetation parameters
    TYPE (soil_parameter_type), INTENT(INOUT) :: soil ! soil parameters
    TYPE (casa_met),            INTENT(INOUT) :: casamet
    TYPE (casa_balance),        INTENT(INOUT) :: casabal
    INTEGER,               INTENT(IN)    :: myear
    !  REAL(r_2),    INTENT(IN) :: clitterinput(mp,3),csoilinput(mp,3)

    ! local variables
    INTEGER  npt,nout
    REAL(r_2) xyear, totGPP, totNPP

    totGPP =0.0
    totNPP =0.0
    nout=104
    xyear=1.0/FLOAT(myear)
    casabal%FCgppyear=casabal%FCgppyear * xyear
    casabal%FCnppyear=casabal%FCnppyear * xyear
    casabal%FCrmleafyear=casabal%FCrmleafyear * xyear
    casabal%FCrmwoodyear=casabal%FCrmwoodyear * xyear
    casabal%FCrmrootyear=casabal%FCrmrootyear * xyear
    casabal%FCrgrowyear=casabal%FCrgrowyear * xyear
    casabal%FCrsyear=casabal%FCrsyear * xyear
    casabal%FCneeyear=casabal%FCneeyear * xyear
    casabal%FNdepyear=casabal%FNdepyear * xyear
    casabal%FNfixyear=casabal%FNfixyear * xyear
    casabal%FNsnetyear=casabal%FNsnetyear * xyear
    casabal%FNupyear=casabal%FNupyear * xyear
    casabal%FNleachyear=casabal%FNleachyear * xyear
    casabal%FNlossyear=casabal%FNlossyear * xyear
    casabal%FPweayear=casabal%FPweayear * xyear
    casabal%FPdustyear=casabal%FPdustyear * xyear
    casabal%FPsnetyear=casabal%FPsnetyear * xyear
    casabal%FPupyear=casabal%FPupyear * xyear
    casabal%FPleachyear=casabal%FPleachyear * xyear
    casabal%FPlossyear=casabal%FPlossyear * xyear
    !  clitterinput = clitterinput * xyear
    !  csoilinput   = csoilinput   * xyear

    PRINT *, 'writing CNP fluxes out to file ', casafile%cnpflux
    OPEN(nout,file=casafile%cnpflux)
    DO npt =1,mp
       SELECT CASE(icycle)
       CASE(1)

          WRITE(nout,*) myear,npt,veg%iveg(npt),soil%isoilm(npt),    &
               casamet%isorder(npt),casamet%lat(npt),casamet%lon(npt), &
               casamet%areacell(npt)*(1.0e-9),casabal%Fcgppyear(npt),  &
               casabal%Fcnppyear(npt),  &
               casabal%Fcrmleafyear(npt),casabal%Fcrmwoodyear(npt),     &
               casabal%Fcrmrootyear(npt),casabal%Fcrgrowyear(npt),     &
               casabal%Fcrsyear(npt),casabal%Fcneeyear(npt)  ! ,           &
          !            clitterinput(npt,:),csoilinput(npt,:)

       CASE(2)
          WRITE(nout,*) myear,npt,veg%iveg(npt),soil%isoilm(npt),    &
               casamet%isorder(npt),casamet%lat(npt),casamet%lon(npt), &
               casamet%areacell(npt)*(1.0e-9),casabal%Fcgppyear(npt),  &
               casabal%FCnppyear(npt),                                 &
               casabal%FCrmleafyear(npt),casabal%FCrmwoodyear(npt),     &
               casabal%FCrmrootyear(npt),casabal%FCrgrowyear(npt),     &
               casabal%FCrsyear(npt), casabal%FCneeyear(npt),          &
                                !        clitterinput(npt,:),csoilinput(npt,:), &
               casabal%FNdepyear(npt),casabal%FNfixyear(npt),casabal%FNsnetyear(npt), &
               casabal%FNupyear(npt), casabal%FNleachyear(npt),casabal%FNlossyear(npt)

       CASE(3)
          WRITE(nout,*) myear,npt,veg%iveg(npt),soil%isoilm(npt), &
               casamet%isorder(npt),casamet%lat(npt),casamet%lon(npt),  &
               casamet%areacell(npt)*(1.0e-9),casabal%Fcgppyear(npt), &
               casabal%FCnppyear(npt),                                  &
               casabal%Fcrmleafyear(npt),casabal%Fcrmwoodyear(npt),     &
               casabal%Fcrmrootyear(npt),casabal%Fcrgrowyear(npt),     &
               casabal%FCrsyear(npt),   casabal%FCneeyear(npt),         &
                                !        clitterinput(npt,:),csoilinput(npt,:), &
               casabal%FNdepyear(npt),casabal%FNfixyear(npt),  casabal%FNsnetyear(npt),&
               casabal%FNupyear(npt), casabal%FNleachyear(npt),casabal%FNlossyear(npt),&
               casabal%FPweayear(npt),casabal%FPdustyear(npt), casabal%FPsnetyear(npt),&
               casabal%FPupyear(npt), casabal%FPleachyear(npt),casabal%FPlossyear(npt)

       END SELECT
       totGPP = totGPP+casabal%Fcgppyear(npt)* casamet%areacell(npt)


       totNPP = totNPP+casabal%Fcnppyear(npt)* casamet%areacell(npt)
    ENDDO

    PRINT *, 'totGPP global = ', totGPP*(1.0e-15)
    PRINT *, 'totNPP global = ', totNPP*(1.0e-15)
    CLOSE(nout)
92  FORMAT(5(i6,',',2x),100(f15.6,',',2x))
  END SUBROUTINE casa_fluxout

SUBROUTINE casa_cnpflux(casaflux,casapool,casabal,zeroflux)
    USE cable_def_types_mod
    USE casadimension
    USE casaparm
USE casa_pool_type_mod,    ONLY: casa_pool    => casa_pool_type
USE casa_flux_type_mod,    ONLY: casa_flux    => casa_flux_type
USE casa_balance_type_mod, ONLY: casa_balance => casa_bal_type             

IMPLICIT NONE

TYPE (casa_flux),    INTENT(INOUT) :: casaflux
TYPE (casa_pool),    INTENT(INOUT) :: casapool
TYPE (casa_balance), INTENT(INOUT) :: casabal

    LOGICAL :: zeroflux
    !  REAL(r_2), INTENT(INOUT) :: clitterinput(mp,3),csoilinput(mp,3)
    INTEGER n

    IF(zeroflux) THEN
       casabal%FCgppyear    = 0.0
       casabal%FCrpyear     = 0.0
       casabal%FCrmleafyear = 0.0
       casabal%FCrmwoodyear = 0.0
       casabal%FCrmrootyear = 0.0
       casabal%FCrgrowyear  = 0.0
       casabal%FCnppyear    = 0.0
       casabal%FCrsyear     = 0.0
       casabal%FCneeyear    = 0.0
       casabal%dCdtyear    = 0.0


       casabal%FNdepyear    = 0.0
       casabal%FNfixyear    = 0.0
       casabal%FNsnetyear   = 0.0
       casabal%FNupyear     = 0.0
       casabal%FNleachyear  = 0.0
       casabal%FNlossyear   = 0.0

       casabal%FPweayear   = 0.0
       casabal%FPdustyear  = 0.0
       casabal%FPsnetyear  = 0.0
       casabal%FPupyear    = 0.0
       casabal%FPleachyear = 0.0
       casabal%FPlossyear  = 0.0

       casaflux%FluxCtohwp = 0.0
       casaflux%FluxNtohwp = 0.0
       casaflux%FluxPtohwp = 0.0
       casaflux%FluxCtoclear = 0.0
       casaflux%FluxNtoclear = 0.0
       casaflux%FluxPtoclear = 0.0
       casaflux%CtransferLUC = 0.02
    ELSE

       casaflux%Crp(:)   = casaflux%Crmplant(:,leaf) + casaflux%Crmplant(:,wood) + casaflux%Crmplant(:,froot) + casaflux%Crgplant(:)
       casabal%FCgppyear = casabal%FCgppyear + casaflux%Cgpp   * deltpool
       casabal%FCrpyear  = casabal%FCrpyear  + casaflux%Crp    * deltpool
       casabal%FCrmleafyear(:)  = casabal%FCrmleafyear(:)  + casaflux%Crmplant(:,leaf)    * deltpool
       casabal%FCrmwoodyear(:)  = casabal%FCrmwoodyear(:)  + casaflux%Crmplant(:,wood)    * deltpool
       casabal%FCrmrootyear(:)  = casabal%FCrmrootyear(:)  + casaflux%Crmplant(:,froot)   * deltpool
       casabal%FCrgrowyear      = casabal%FCrgrowyear      + casaflux%Crgplant            * deltpool
       ! change made ypwang 17-nov-2013 to accoutn for change in labile carbon pool  size
       casabal%FCnppyear        = casabal%FCnppyear + (casaflux%Cnpp+casapool%dClabiledt)   * deltpool
       casabal%FCrsyear  = casabal%FCrsyear  + casaflux%Crsoil * deltpool
       casabal%FCneeyear = casabal%FCneeyear &
            + (casaflux%Cnpp+casapool%dClabiledt-casaflux%Crsoil) * deltpool
       casabal%dCdtyear =  casabal%dCdtyear + (casapool%Ctot-casapool%Ctot_0)*deltpool

       !  DO n=1,3
       !    clitterinput(:,n)= clitterinput(:,n) + casaflux%kplant(:,n) * casapool%cplant(:,n) * deltpool
       !    csoilinput(:,n) = csoilinput(:,n) + casaflux%fluxCtosoil(:,n) * deltpool
       !    !csoilinput(:,n) = csoilinput(:,n)+casaflux%fluxCtolitter(:,n)*deltpool
       !  ENDDO

       IF (icycle >1) THEN
          casabal%FNdepyear   = casabal%FNdepyear   + casaflux%Nmindep    * deltpool
          casabal%FNfixyear   = casabal%FNfixyear   + casaflux%Nminfix    * deltpool
          casabal%FNsnetyear  = casabal%FNsnetyear  + casaflux%Nsnet      * deltpool
          casabal%FNupyear    = casabal%FNupyear    + casaflux%Nminuptake * deltpool
          casabal%FNleachyear = casabal%FNleachyear + casaflux%Nminleach  * deltpool
          casabal%FNlossyear  = casabal%FNlossyear  + casaflux%Nminloss   * deltpool
       ENDIF

       IF (icycle >2) THEN
          casabal%FPweayear   = casabal%FPweayear   + casaflux%Pwea       * deltpool
          casabal%FPdustyear  = casabal%FPdustyear  + casaflux%Pdep       * deltpool
          casabal%FPsnetyear  = casabal%FPsnetyear  + casaflux%Psnet      * deltpool
          casabal%FPupyear    = casabal%FPupyear    + casaflux%Plabuptake * deltpool
          casabal%FPleachyear = casabal%FPleachyear + casaflux%Pleach     * deltpool
          casabal%FPlossyear  = casabal%FPlossyear  + casaflux%Ploss      * deltpool
       ENDIF
    ENDIF
  END SUBROUTINE casa_cnpflux

END MODULE casa_inout_module
