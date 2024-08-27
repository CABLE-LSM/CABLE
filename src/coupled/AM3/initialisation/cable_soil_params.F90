MODULE cable_soil_params_mod

   IMPLICIT NONE 

CONTAINS

subroutine cable_soil_params(pars_io_cbl)

   ! Gets parameter values for each vegetation type and soil type.
   
USE params_io_mod_cbl,        ONLY: params_io_data_type

  IMPLICIT NONE
  
  TYPE(params_io_data_type),      TARGET :: pars_io_cbl

  
 !SOIL parameters: description and corresponding variable name in code. 
 !SOIL parameters are assigned as TYPE pars_io_cbl%soilin_ but later used as soil%
 
!SOIL: Coarse sand/Loamy sand                                                
! =========================================================
    pars_io_cbl%soilin_silt( 1) =        0.080000
    pars_io_cbl%soilin_clay( 1) =        0.090000
    pars_io_cbl%soilin_sand( 1) =        0.830000
   pars_io_cbl%soilin_swilt( 1) =        0.072000
     pars_io_cbl%soilin_sfc( 1) =        0.143000
    pars_io_cbl%soilin_ssat( 1) =        0.398000
     pars_io_cbl%soilin_bch( 1) =        4.200000
    pars_io_cbl%soilin_hyds( 1) =        0.000166
    pars_io_cbl%soilin_sucs( 1) =       -0.106000
 pars_io_cbl%soilin_rhosoil( 1) =     1600.000000
     pars_io_cbl%soilin_css( 1) =      850.000000
 
 !SOIL: Medium clay loam/silty clay loam/silt loam                            
 !=========================================================
    pars_io_cbl%soilin_silt( 2) =        0.330000
    pars_io_cbl%soilin_clay( 2) =        0.300000
    pars_io_cbl%soilin_sand( 2) =        0.370000
   pars_io_cbl%soilin_swilt( 2) =        0.216000
     pars_io_cbl%soilin_sfc( 2) =        0.301000
    pars_io_cbl%soilin_ssat( 2) =        0.479000
     pars_io_cbl%soilin_bch( 2) =        7.100000
    pars_io_cbl%soilin_hyds( 2) =        0.000004
    pars_io_cbl%soilin_sucs( 2) =       -0.591000
 pars_io_cbl%soilin_rhosoil( 2) =     1600.000000
     pars_io_cbl%soilin_css( 2) =      850.000000
 
 !SOIL: Fine clay                                                             
 !=========================================================
    pars_io_cbl%soilin_silt( 3) =        0.170000
    pars_io_cbl%soilin_clay( 3) =        0.670000
    pars_io_cbl%soilin_sand( 3) =        0.160000
   pars_io_cbl%soilin_swilt( 3) =        0.286000
     pars_io_cbl%soilin_sfc( 3) =        0.367000
    pars_io_cbl%soilin_ssat( 3) =        0.482000
     pars_io_cbl%soilin_bch( 3) =       11.400000
    pars_io_cbl%soilin_hyds( 3) =        0.000001
    pars_io_cbl%soilin_sucs( 3) =       -0.405000
 pars_io_cbl%soilin_rhosoil( 3) =     1381.000000
     pars_io_cbl%soilin_css( 3) =      850.000000
 
 !SOIL: Coarse-medium sandy loam/loam                                         
 !=========================================================
    pars_io_cbl%soilin_silt( 4) =        0.200000
    pars_io_cbl%soilin_clay( 4) =        0.200000
    pars_io_cbl%soilin_sand( 4) =        0.600000
   pars_io_cbl%soilin_swilt( 4) =        0.135000
     pars_io_cbl%soilin_sfc( 4) =        0.218000
    pars_io_cbl%soilin_ssat( 4) =        0.443000
     pars_io_cbl%soilin_bch( 4) =        5.150000
    pars_io_cbl%soilin_hyds( 4) =        0.000021
    pars_io_cbl%soilin_sucs( 4) =       -0.348000
 pars_io_cbl%soilin_rhosoil( 4) =     1373.000000
     pars_io_cbl%soilin_css( 4) =      850.000000
 
 !SOIL: Coarse-fine sandy clay                                                
 !=========================================================
    pars_io_cbl%soilin_silt( 5) =        0.060000
    pars_io_cbl%soilin_clay( 5) =        0.420000
    pars_io_cbl%soilin_sand( 5) =        0.520000
   pars_io_cbl%soilin_swilt( 5) =        0.219000
     pars_io_cbl%soilin_sfc( 5) =        0.310000
    pars_io_cbl%soilin_ssat( 5) =        0.426000
     pars_io_cbl%soilin_bch( 5) =       10.400000
    pars_io_cbl%soilin_hyds( 5) =        0.000002
    pars_io_cbl%soilin_sucs( 5) =       -0.153000
 pars_io_cbl%soilin_rhosoil( 5) =     1476.000000
     pars_io_cbl%soilin_css( 5) =      850.000000
 
 !SOIL: Medium-fine silty clay                                                
 !=========================================================
    pars_io_cbl%soilin_silt( 6) =        0.250000
    pars_io_cbl%soilin_clay( 6) =        0.480000
    pars_io_cbl%soilin_sand( 6) =        0.270000
   pars_io_cbl%soilin_swilt( 6) =        0.283000
     pars_io_cbl%soilin_sfc( 6) =        0.370000
    pars_io_cbl%soilin_ssat( 6) =        0.482000
     pars_io_cbl%soilin_bch( 6) =       10.400000
    pars_io_cbl%soilin_hyds( 6) =        0.000001
    pars_io_cbl%soilin_sucs( 6) =       -0.490000
 pars_io_cbl%soilin_rhosoil( 6) =     1521.000000
     pars_io_cbl%soilin_css( 6) =      850.000000
 
 !SOIL: Coarse-medium-fine sandy clay loam                                    
 !=========================================================
    pars_io_cbl%soilin_silt( 7) =        0.150000
    pars_io_cbl%soilin_clay( 7) =        0.270000
    pars_io_cbl%soilin_sand( 7) =        0.580000
   pars_io_cbl%soilin_swilt( 7) =        0.175000
     pars_io_cbl%soilin_sfc( 7) =        0.255000
    pars_io_cbl%soilin_ssat( 7) =        0.420000
     pars_io_cbl%soilin_bch( 7) =        7.120000
    pars_io_cbl%soilin_hyds( 7) =        0.000006
    pars_io_cbl%soilin_sucs( 7) =       -0.299000
 pars_io_cbl%soilin_rhosoil( 7) =     1373.000000
     pars_io_cbl%soilin_css( 7) =      850.000000
 
 !SOIL: Organic peat                                                          
 !=========================================================
    pars_io_cbl%soilin_silt( 8) =        0.700000
    pars_io_cbl%soilin_clay( 8) =        0.170000
    pars_io_cbl%soilin_sand( 8) =        0.130000
   pars_io_cbl%soilin_swilt( 8) =        0.395000
     pars_io_cbl%soilin_sfc( 8) =        0.450000
    pars_io_cbl%soilin_ssat( 8) =        0.451000
     pars_io_cbl%soilin_bch( 8) =        5.830000
    pars_io_cbl%soilin_hyds( 8) =        0.000800
    pars_io_cbl%soilin_sucs( 8) =       -0.356000
 pars_io_cbl%soilin_rhosoil( 8) =     1537.000000
     pars_io_cbl%soilin_css( 8) =     1920.000000
 
 !SOIL: Permanent ice                                                         
 !=========================================================
    pars_io_cbl%soilin_silt( 9) =        0.330000
    pars_io_cbl%soilin_clay( 9) =        0.300000
    pars_io_cbl%soilin_sand( 9) =        0.370000
   pars_io_cbl%soilin_swilt( 9) =        0.216000
     pars_io_cbl%soilin_sfc( 9) =        0.301000
    pars_io_cbl%soilin_ssat( 9) =        0.479000
     pars_io_cbl%soilin_bch( 9) =        7.100000
    pars_io_cbl%soilin_hyds( 9) =        0.000001
    pars_io_cbl%soilin_sucs( 9) =       -0.153000
  pars_io_cbl%soilin_rhosoil( 9) =      917.000000
     pars_io_cbl%soilin_css( 9) =     2100.000000    

End subroutine cable_soil_params

END MODULE cable_soil_params_mod

