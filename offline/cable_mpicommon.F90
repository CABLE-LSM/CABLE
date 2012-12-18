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
! Purpose: common modules for MPI wrapper for CABLE
!
! Contact: Bernard.Pak@csiro.au
!
! History: MPI wrapper developed by Maciej Golebiewski (2012)
!
! ==============================================================================
!
MODULE cable_mpicommon

  USE cable_def_types_mod

  IMPLICIT NONE

  SAVE

  PUBLIC

  ! MPI: base number of input fields sent to workers as start up
  ! parameters
  !INTEGER, PARAMETER :: nparam = 68
  !INTEGER, PARAMETER :: nparam = 219
  ! MPI: Bernard commented out two canopy params (potev_c and rwater)
  ! when porting to CABLE_r491
  ! INTEGER, PARAMETER :: nparam = 217
  ! MPI: CABLE_r491, after following up with Bernard on the new variables
  !INTEGER, PARAMETER :: nparam = 260
  ! added 23 params when trying to fix the bug in MPI
  INTEGER, PARAMETER :: nparam = 283

  ! MPI: extra params sent only if nsoilparmnew is true
  INTEGER, PARAMETER :: nsoilnew = 1

  ! MPI: number of casa parameters sent to workers as
  ! start up parameters
  !INTEGER, PARAMETER :: ncasaparam = 68
  !INTEGER, PARAMETER :: ncasaparam = 176
  ! MPI: added casapool fields ratioNCsoilnew, ratioNCsoilmin and ratioNCsoilmax
  INTEGER, PARAMETER :: ncasaparam = 179

  ! MPI: base number of casa_init parameters sent to the workers
  INTEGER, PARAMETER :: ncinit = 18

  ! MPI: number of casa_init parameters sent to the workers only if
  ! icycle = 2 or 3
  INTEGER, PARAMETER :: ncinit2 = 13

  ! MPI: number of casa_init parameters sent to the workers only if
  ! icycle = 3
  INTEGER, PARAMETER :: ncinit3 = 18

  ! MPI: number of input fields sent to workers at the start of each
  ! timestep
  !INTEGER, PARAMETER :: ninput = 11
  ! added 4 time fields in met: year, moy, doy, hod
  INTEGER, PARAMETER :: ninput = 15

  ! MPI: number of 3D array slices / worker (results)
  INTEGER, PARAMETER :: n3d = 1

  ! MPI: number of matrix slices / worker (results)
  !INTEGER, PARAMETER :: nmat = 29
  ! MPI: 2011-07-08 - removed dtmlt from data exchange
  !INTEGER, PARAMETER :: nmat = 28
  ! MPI: gol124: net +1 when Bernard ported to CABLE_r491
  !INTEGER, PARAMETER :: nmat = 29
  ! MPI: CABLE_r491, after following up with Bernard on the new variables
  INTEGER, PARAMETER :: nmat = 36

  ! MPI: number of contig vector parts / worker (results)
  !INTEGER, PARAMETER :: nvec = 149
  ! MPI: 2011-06-28 - removed ebal, ebal_tot, seb, seb_tot from data exchange
  !INTEGER, PARAMETER :: nvec = 145
  ! MPI: 2011-07-08 - removed otss from data exchange
  ! INTEGER, PARAMETER :: nvec = 144
  ! MPI: 2012-02-14 - removed year, moy, doy, hod
  !INTEGER, PARAMETER :: nvec = 140
  ! MPI: gol124: net -3 (removed or changed to 2D) when Bernard
  ! ported to CABLE_r491
  !INTEGER, PARAMETER :: nvec = 137
  ! MPI: CABLE_r491, after following up with Bernard on the new variables
  INTEGER, PARAMETER :: nvec = 161

  ! MPI: number of final casa result matrices and vectors to receive
  ! by the master for casa_poolout and casa_fluxout
  INTEGER, PARAMETER :: ncasa_mat = 15
  INTEGER, PARAMETER :: ncasa_vec = 27

  ! MPI: number of fields included in restart_t type for data
  ! that is returned only for creating a restart file at the end of the run
  !INTEGER, PARAMETER :: nrestart = 16
  ! MPI: gol124: canopy%rwater removed when Bernard ported to CABLE_r491
  INTEGER, PARAMETER :: nrestart = 15

  ! MPI: type to hold landpoint decomposition info
  TYPE lpdecomp_t
          INTEGER :: landp0      ! starting land point index
          INTEGER :: nland       ! number of landpoints

          INTEGER :: patch0      ! starting patch index in global CABLE vars
          INTEGER :: npatch      ! sum of patches for all landpoints of this
                                 ! worker
  END TYPE

  ! MPI: worker's local landpoints and patches
  TYPE(lpdecomp_t) :: wpatch

  ! MPI: Fortran types extents
  INTEGER :: extr1, extr2, extid, extl

CONTAINS

! calculates extents of the Fortran types used by CABLE
SUBROUTINE find_extents

  USE mpi
  USE cable_def_types_mod

  IMPLICIT NONE

  INTEGER, DIMENSION(2) :: itmp
  REAL, DIMENSION(2) :: r1tmp
  REAL(r_2), DIMENSION(2) :: r2tmp
  LOGICAL, DIMENSION(2) :: ltmp

  INTEGER(KIND=MPI_ADDRESS_KIND), DIMENSION(2) :: a

  INTEGER :: ierr

  CALL MPI_Get_address (itmp(1), a(1), ierr)
  CALL MPI_Get_address (itmp(2), a(2), ierr)
  extid = INT(a(2)-a(1))

  CALL MPI_Get_address (r1tmp(1), a(1), ierr)
  CALL MPI_Get_address (r1tmp(2), a(2), ierr)
  extr1 = INT(a(2)-a(1))

  CALL MPI_Get_address (r2tmp(1), a(1), ierr)
  CALL MPI_Get_address (r2tmp(2), a(2), ierr)
  extr2 = INT(a(2)-a(1))

  CALL MPI_Get_address (ltmp(1), a(1), ierr)
  CALL MPI_Get_address (ltmp(2), a(2), ierr)
  extl = INT(a(2)-a(1))

END SUBROUTINE find_extents

! creates MPI derived datatypes for exchanging landpt and patch arrays
SUBROUTINE decomp_types (landpt_t, patch_t)

  USE mpi

  USE cable_IO_vars_module

  IMPLICIT NONE

  INTEGER, INTENT(OUT) :: landpt_t, patch_t

  ! dummy vars to calculate field offsets
  TYPE(land_type) :: dlandpt
  TYPE(patch_type) :: dpatch

  INTEGER(KIND=MPI_ADDRESS_KIND) :: base_d

  INTEGER, PARAMETER :: fields = 5
  INTEGER, DIMENSION(fields) :: blocks, types
  INTEGER(KIND=MPI_ADDRESS_KIND), DIMENSION(fields) :: displs

  INTEGER :: ierr

  blocks = 1

  ! create MPI type to exchange landpt records
  types = MPI_INTEGER
  
  CALL MPI_Get_address (dlandpt, base_d, ierr)

  CALL MPI_Get_address (dlandpt%nap, displs(1), ierr)
  CALL MPI_Get_address (dlandpt%cstart, displs(2), ierr)
  CALL MPI_Get_address (dlandpt%cend, displs(3), ierr)
  CALL MPI_Get_address (dlandpt%ilat, displs(4), ierr)
  CALL MPI_Get_address (dlandpt%ilon, displs(5), ierr)

  displs = displs - base_d

  CALL MPI_Type_create_struct (5, blocks, displs, types, landpt_t, ierr)
  CALL MPI_Type_commit (landpt_t, ierr)

  ! create MPI type to exchange patch records
  types = MPI_REAL
  
  CALL MPI_Get_address (dpatch, base_d, ierr)

  CALL MPI_Get_address (dpatch%frac, displs(1), ierr)
  CALL MPI_Get_address (dpatch%latitude, displs(2), ierr)
  CALL MPI_Get_address (dpatch%longitude, displs(3), ierr)

  displs = displs - base_d

  CALL MPI_Type_create_struct (3, blocks, displs, types, patch_t, ierr)
  CALL MPI_Type_commit (patch_t, ierr)

  RETURN

END SUBROUTINE decomp_types

SUBROUTINE bcast_start_time (comm)

  USE mpi

  USE cable_IO_vars_module

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: comm

  INTEGER :: ierr

  CALL MPI_Bcast (shod, 1, MPI_REAL, 0, comm, ierr)
  CALL MPI_Bcast (sdoy, 1, MPI_INTEGER, 0, comm, ierr)
  CALL MPI_Bcast (smoy, 1, MPI_INTEGER, 0, comm, ierr)
  CALL MPI_Bcast (syear, 1, MPI_INTEGER, 0, comm, ierr)
  CALL MPI_Bcast (time_coord, 3, MPI_CHARACTER, 0, comm, ierr)

  RETURN

END SUBROUTINE bcast_start_time

END MODULE cable_mpicommon
