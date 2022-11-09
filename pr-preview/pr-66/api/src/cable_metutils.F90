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
! Purpose: lookup tables for possible names for the met. forcing variables in
!          netcdf files
!
! Contact: jhan.srbinovsky@csiro.au
!
! History: First version by Claire Carouge (09/2021)
!
! ==============================================================================
!
module CABLE_METUTILS_MODULE
  USE netcdf
  
  implicit none
    
  PUBLIC ::            &
    possible_varnames, &  
    ! Subroutines
    find_metvarid 

  TYPE input_netcdf_names
    ! Latitude variable
    CHARACTER(LEN=8),DIMENSION(3) :: LatNames=(/                                &
                          'latitude', 'nav_lat ', 'lat     '                    &
                          /)

    ! Longitude variable
    CHARACTER(LEN=9),DIMENSION(3) :: LonNames=(/                                &
                          'longitude', 'nav_lon  ', 'lon      '                 &
                          /)

    ! Mask variable
    CHARACTER(LEN=7),DIMENSION(2) :: MaskNames=(/                               &
                          'mask   ', 'landsea'                                  &
                          /)

    ! Time variable
    CHARACTER(LEN=4),DIMENSION(1) :: TimeNames=(/                               &
                          'time'                                                &
                          /)

    ! Downward Shortwave radiation variable
    CHARACTER(LEN=6),DIMENSION(4) :: SWdownNames=(/                             &
                          'dswrf ', 'rsds  ', 'FSDS  ', 'SWdown'                &
                          /)

    ! Air temperature variable
    CHARACTER(LEN=4),DIMENSION(3) :: TairNames=(/                               &
                          'tas ', 'TBOT', 'Tair'                                &
                          /)

    ! Air humidity variable
    CHARACTER(LEN=4),DIMENSION(4) :: QairNames=(/                               &
                          'shum', 'huss', 'QBOT', 'Qair'                        &
                          /)

    ! Wind variable
    CHARACTER(LEN=4),DIMENSION(3) :: WindNames=(/                            &
                          'wind', 'Wind', 'WIND'                                &
                          /)

    ! Rain variable
    CHARACTER(LEN=8),DIMENSION(5) :: RainNames=(/                            &
                          'prcp    ', 'pr      ', 'RAIN    ', 'Rainf   ',       &
                          'Precip  '                                            &
                          /)

    ! Downward longwave radiation variable
    CHARACTER(LEN=6),DIMENSION(4) :: LWdownNames=(/                          &
                          'dlwrf ', 'rlds  ', 'FLDS  ', 'LWdown'                &
                          /)
      
    ! Surface pressure variable
    CHARACTER(LEN=5),DIMENSION(5) :: PSurfNames=(/                           &
                          'pres ', 'ps   ', 'PBOT ', 'PSurf', 'Psurf'        &
                          /)
      
    ! Elevation variable
    CHARACTER(LEN=9),DIMENSION(1) :: ElevNames=(/                            &
                          'Elevation'                                           &
                          /)
      
    ! CO2 concentration variable
    CHARACTER(LEN=6),DIMENSION(1) :: CO2Names=(/                             &
                          'CO2air'                                              &
                          /)
      
    ! Snow variable
    CHARACTER(LEN=5),DIMENSION(1) :: SnowNames=(/                            &
                          'Snowf'                                               &
                          /)
      
    ! LAI variable
    CHARACTER(LEN=3),DIMENSION(1) :: LAINames=(/                             &
                          'LAI'                                                 &
                          /)
      
    ! avPrecip variable
    CHARACTER(LEN=8),DIMENSION(1) :: APrecipNames=(/                         &
                          'avPrecip'                                            &
                          /)
      
    ! Vegetation types variable
    CHARACTER(LEN=4),DIMENSION(1) :: IVegNames=(/                            &
                          'iveg'                                                &
                          /)
      
    ! Patch fraction variable
    CHARACTER(LEN=9),DIMENSION(1) :: PFracNames=(/                           &
                          'patchfrac'                                           &
                          /)
      
    ! Soil type variable
    CHARACTER(LEN=5),DIMENSION(1) :: ISoilNames=(/                           &
                          'isoil'                                               &
                          /)

  END TYPE input_netcdf_names

  type(input_netcdf_names) :: possible_varnames

CONTAINS
  !==============================================================================
  !
  ! Name: find_metvarid
  !
  ! Purpose: Find the var ID in the met. file. It will use a lookup table of 
  !   possible variable names to find the variable in the netcdf file. 
  !
  ! CALLed from: open_met_file
  !
  !==============================================================================
  SUBROUTINE find_metvarid(file_id, possible_names, varid, status)
    
    IMPLICIT NONE

    ! Arguments
    ! ID of the netcdf with the met. forcing
    INTEGER, INTENT(IN) :: file_id  
    ! Lookup table for names of the variable in netcdf files
    CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: possible_names

    ! Output: variable id in netcdf file
    INTEGER, INTENT(OUT) :: varid 
    ! Output: error status if variable not found.
    INTEGER, INTENT(OUT) :: status
    
    ! Local variables
    INTEGER :: i  ! Loop index

    ! Initialise varid
    varid = -10
    DO i = 1, size(possible_names)
      status = NF90_INQ_VARID(file_id,trim(possible_names(i)),varid)
      ! Stop the loop if variable found
      IF (status == NF90_NOERR) exit
    END DO

  END SUBROUTINE find_metvarid

end module CABLE_METUTILS_MODULE