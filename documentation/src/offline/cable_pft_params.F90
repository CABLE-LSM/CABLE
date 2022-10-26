MODULE cable_pft_params_mod
USE grid_constants_mod_cbl, ONLY: ntype_max

IMPLICIT NONE 

TYPE vegin_type

CHARACTER(LEN=70) :: desc(ntype_max)         ! decriptions of veg type

REAL ::                                                                        &
   canst1(ntype_max),                                                          &
   length(ntype_max),                                                          &
   width(ntype_max),                                                           &
   vcmax(ntype_max),                                                           &
   ejmax(ntype_max),                                                           &
   hc(ntype_max),                                                              &
   xfang(ntype_max),                                                           &
   rp20(ntype_max),                                                            &
   rpcoef(ntype_max),                                                          &
   rs20(ntype_max),                                                            &
   wai(ntype_max),                                                             &
   rootbeta(ntype_max),                                                        &
   shelrb(ntype_max),                                                          &
   vegcf(ntype_max),                                                           &
   frac4(ntype_max),                                                           &
   xalbnir(ntype_max),                                                         &
   extkn(ntype_max),                                                           &
   tminvj(ntype_max),                                                          &
   tmaxvj(ntype_max),                                                          &
   vbeta(ntype_max),                                                           &
   a1gs(ntype_max),                                                            &
   d0gs(ntype_max),                                                            &
   alpha(ntype_max),                                                           &
   convex(ntype_max),                                                          &
   cfrd(ntype_max),                                                            &
   gswmin(ntype_max),                                                          &
   conkc0(ntype_max),                                                          &
   conko0(ntype_max),                                                          &
   ekc(ntype_max),                                                             &
   eko(ntype_max),                                                             &
   g0(ntype_max),                                                              &
   g1(ntype_max),                                                              &
   zr(ntype_max),                                                              &
   clitt(ntype_max),                                                           &
   froot1(ntype_max),                                                          &
   froot2(ntype_max),                                                          &
   froot3(ntype_max),                                                          &
   froot4(ntype_max),                                                          &
   froot5(ntype_max),                                                          &
   froot6(ntype_max),                                                          &
   csoil1(ntype_max),                                                          &
   csoil2(ntype_max),                                                          &
   ratecs1(ntype_max),                                                         &
   ratecs2(ntype_max),                                                         &
   cplant1(ntype_max),                                                         &
   cplant2(ntype_max),                                                         &
   cplant3(ntype_max),                                                         &
   ratecp1(ntype_max),                                                         &
   ratecp2(ntype_max),                                                         &
   ratecp3(ntype_max),                                                         &
   refl1(ntype_max),                                                           &
   refl2(ntype_max),                                                           &
   refl3(ntype_max),                                                           &
   taul1(ntype_max),                                                           &
   taul2(ntype_max),                                                           &
   taul3(ntype_max),                                                           &
   dleaf(ntype_max),                                                           &
   lai(ntype_max)

END TYPE vegin_type


TYPE(vegin_type) :: vegin

CHARACTER(LEN=70) :: veg_desc(ntype_max)         ! decriptions of veg type

CONTAINS

subroutine cable_pft_params()

! Gets parameter values for each vegetation type 
USE cable_def_types_mod, ONLY : mvtype

integer :: ERROR
integer, parameter :: namelist_unit=711179
integer :: j
CHARACTER(LEN=*), parameter :: nml_dir='./' 
CHARACTER(LEN=*), parameter :: iomessage='something wrong with your PFT params file' 
CHARACTER(LEN=*), PARAMETER :: routinename='cable_pft_params'

NAMELIST / cable_pftparm/ vegin

!HACK:offline(cable_parameters) checks mvtype
mvtype=ntype_max    
!-----------------------------------------------------------------------------
! Read namelist
!-----------------------------------------------------------------------------
write (6,*) "Reading CABLE_PFTPARM namelist..."

OPEN( namelist_unit, FILE=(TRIM(nml_dir) // '/' // 'pft_params.nml'),          &
      STATUS='old', POSITION='rewind', ACTION='read', IOSTAT  = ERROR )
IF ( ERROR /= 0 ) write (6,*) "Error opening  CABLE_PFTPARM namelist..."

READ(namelist_unit, NML = cable_pftparm, IOSTAT = ERROR )
IF ( ERROR /= 0 ) write (6,*) "Error reading  CABLE_PFTPARM namelist..."
                 
CLOSE(namelist_unit, IOSTAT = ERROR)

! new calculation dleaf since April 2012 (cable v1.8 did not use width)
vegin%dleaf = SQRT(vegin%width * vegin%length)

veg_desc = vegin%desc
    
End subroutine cable_pft_params

END MODULE cable_pft_params_mod

