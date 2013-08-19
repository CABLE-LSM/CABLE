!==============================================================================
! This source code is part of the
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CABLE Academic User Licence Agreement
! (the "Licence").
! You may not use this file except in compliance with the Licence.
! A copy of the Licence and registration form can be obtained from
! http://www.accessimulator.org.au/cable
! You need to register and read the Licence agreement before use.
! Please contact cable_help@nf.nci.org.au for any questions on
! registration and the Licence.
!
! Unless required by applicable law or agreed to in writing,
! software distributed under the Licence is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the Licence for the specific language governing permissions and
! limitations under the Licence.
! ==============================================================================
!
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
   
