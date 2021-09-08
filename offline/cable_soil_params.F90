MODULE cable_soil_params_mod
   IMPLICIT NONE 

   TYPE soilin_type

      REAL, DIMENSION(:),ALLOCATABLE ::                                        &
         silt,    & !
         clay,    & !
         sand,    & !
         swilt,   & !
         sfc,     & !
         ssat,    & !
         bch,     & !
         hyds,    & !
         sucs,    & !
         rhosoil, & !
         css,     & !
         c3         !
   
   END TYPE soilin_type
 
   CHARACTER(LEN=70), DIMENSION(:), POINTER ::                                 &
      soil_desc     ! decriptns of soil type 

   TYPE(soilin_type), SAVE  :: soilin

CONTAINS


subroutine cable_soil_params()

   ! Gets parameter values for each vegetation type and soil type.
   
   USE cable_def_types_mod, ONLY : mstype
  implicit none
  logical, save :: first_call = .true.

  if( first_call ) then
  
  !hard-wired #  of soil types, promote to nml
  mstype = 9 

      ALLOCATE ( soil_desc(mstype) )
      ALLOCATE ( soilin%silt(mstype), soilin%clay(mstype), soilin%sand(mstype) )
      ALLOCATE ( soilin%swilt(mstype), soilin%sfc(mstype), soilin%ssat(mstype) )
      ALLOCATE ( soilin%bch(mstype), soilin%hyds(mstype), soilin%sucs(mstype) )
      ALLOCATE ( soilin%rhosoil(mstype), soilin%css(mstype) )
 
 !SOIL parameters: description and corresponding variable name in code. 
 !SOIL parameters are assigned as TYPE soilin% but later used as soil%
 
!SOIL: Coarse sand/Loamy sand                                                
! =========================================================
    soilin%silt( 1) =        0.080000
    soilin%clay( 1) =        0.090000
    soilin%sand( 1) =        0.830000
   soilin%swilt( 1) =        0.072000
     soilin%sfc( 1) =        0.143000
    soilin%ssat( 1) =        0.398000
     soilin%bch( 1) =        4.200000
    soilin%hyds( 1) =        0.000166
    soilin%sucs( 1) =       -0.106000
 soilin%rhosoil( 1) =     1600.000000
     soilin%css( 1) =      850.000000
 
 !SOIL: Medium clay loam/silty clay loam/silt loam                            
 !=========================================================
    soilin%silt( 2) =        0.330000
    soilin%clay( 2) =        0.300000
    soilin%sand( 2) =        0.370000
   soilin%swilt( 2) =        0.216000
     soilin%sfc( 2) =        0.301000
    soilin%ssat( 2) =        0.479000
     soilin%bch( 2) =        7.100000
    soilin%hyds( 2) =        0.000004
    soilin%sucs( 2) =       -0.591000
 soilin%rhosoil( 2) =     1600.000000
     soilin%css( 2) =      850.000000
 
 !SOIL: Fine clay                                                             
 !=========================================================
    soilin%silt( 3) =        0.170000
    soilin%clay( 3) =        0.670000
    soilin%sand( 3) =        0.160000
   soilin%swilt( 3) =        0.286000
     soilin%sfc( 3) =        0.367000
    soilin%ssat( 3) =        0.482000
     soilin%bch( 3) =       11.400000
    soilin%hyds( 3) =        0.000001
    soilin%sucs( 3) =       -0.405000
 soilin%rhosoil( 3) =     1381.000000
     soilin%css( 3) =      850.000000
 
 !SOIL: Coarse-medium sandy loam/loam                                         
 !=========================================================
    soilin%silt( 4) =        0.200000
    soilin%clay( 4) =        0.200000
    soilin%sand( 4) =        0.600000
   soilin%swilt( 4) =        0.135000
     soilin%sfc( 4) =        0.218000
    soilin%ssat( 4) =        0.443000
     soilin%bch( 4) =        5.150000
    soilin%hyds( 4) =        0.000021
    soilin%sucs( 4) =       -0.348000
 soilin%rhosoil( 4) =     1373.000000
     soilin%css( 4) =      850.000000
 
 !SOIL: Coarse-fine sandy clay                                                
 !=========================================================
    soilin%silt( 5) =        0.060000
    soilin%clay( 5) =        0.420000
    soilin%sand( 5) =        0.520000
   soilin%swilt( 5) =        0.219000
     soilin%sfc( 5) =        0.310000
    soilin%ssat( 5) =        0.426000
     soilin%bch( 5) =       10.400000
    soilin%hyds( 5) =        0.000002
    soilin%sucs( 5) =       -0.153000
 soilin%rhosoil( 5) =     1476.000000
     soilin%css( 5) =      850.000000
 
 !SOIL: Medium-fine silty clay                                                
 !=========================================================
    soilin%silt( 6) =        0.250000
    soilin%clay( 6) =        0.480000
    soilin%sand( 6) =        0.270000
   soilin%swilt( 6) =        0.283000
     soilin%sfc( 6) =        0.370000
    soilin%ssat( 6) =        0.482000
     soilin%bch( 6) =       10.400000
    soilin%hyds( 6) =        0.000001
    soilin%sucs( 6) =       -0.490000
 soilin%rhosoil( 6) =     1521.000000
     soilin%css( 6) =      850.000000
 
 !SOIL: Coarse-medium-fine sandy clay loam                                    
 !=========================================================
    soilin%silt( 7) =        0.150000
    soilin%clay( 7) =        0.270000
    soilin%sand( 7) =        0.580000
   soilin%swilt( 7) =        0.175000
     soilin%sfc( 7) =        0.255000
    soilin%ssat( 7) =        0.420000
     soilin%bch( 7) =        7.120000
    soilin%hyds( 7) =        0.000006
    soilin%sucs( 7) =       -0.299000
 soilin%rhosoil( 7) =     1373.000000
     soilin%css( 7) =      850.000000
 
 !SOIL: Organic peat                                                          
 !=========================================================
    soilin%silt( 8) =        0.700000
    soilin%clay( 8) =        0.170000
    soilin%sand( 8) =        0.130000
   soilin%swilt( 8) =        0.395000
     soilin%sfc( 8) =        0.450000
    soilin%ssat( 8) =        0.451000
     soilin%bch( 8) =        5.830000
    soilin%hyds( 8) =        0.000800
    soilin%sucs( 8) =       -0.356000
 soilin%rhosoil( 8) =     1537.000000
     soilin%css( 8) =     1920.000000
 
 !SOIL: Permanent ice                                                         
 !=========================================================
    soilin%silt( 9) =        0.330000
    soilin%clay( 9) =        0.300000
    soilin%sand( 9) =        0.370000
   soilin%swilt( 9) =        0.216000
     soilin%sfc( 9) =        0.301000
    soilin%ssat( 9) =        0.479000
     soilin%bch( 9) =        7.100000
    soilin%hyds( 9) =        0.000001
    soilin%sucs( 9) =       -0.153000
  soilin%rhosoil( 9) =      910.000000
     soilin%css( 9) =     2100.000000    
  endif

  first_call = .false.

End subroutine cable_soil_params

END MODULE cable_soil_params_mod

