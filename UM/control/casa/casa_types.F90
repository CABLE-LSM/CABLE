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
! Purpose:
!
! Contact:
!
! History:
!
!
! ==============================================================================

module casa_types_mod

 use casavariable
 use phenvariable

 implicit none

 private

      ! Les 19 Jan 2011
!      TYPE (casa_biome)    , save   :: casabiome  !
!      TYPE (casa_pool)     , save   :: casapool   !
!      TYPE (casa_flux)     , save   :: casaflux   !
!      TYPE (casa_met)      , save   :: casamet    !
!      TYPE (casa_balance)  , save   :: casabal    !
!      TYPE (phen_variable) , save   :: phen       !
      TYPE (casa_biome)       :: casabiome  !
      TYPE (casa_pool)        :: casapool   !
      TYPE (casa_flux)        :: casaflux   !
      TYPE (casa_met)         :: casamet    !
      TYPE (casa_balance)     :: casabal    !
      TYPE (phen_variable)    :: phen       !
     !TYPE (casafiles_type)   :: casafile   !

      ! Save these so only have to do the allocation once.
 save casabiome, casapool, &
      casaflux, casamet, casabal, phen!, casafile
 public casabiome, casapool, &
        casaflux, casamet, casabal, phen!, casafile

end module casa_types_mod
   
