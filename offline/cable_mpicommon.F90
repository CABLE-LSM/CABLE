! ==============================================================================
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
! Purpose: common modules for MPI wrapper for CABLE
!
! Contact: Bernard.Pak@csiro.au
!
! History: MPI wrapper developed by Maciej Golebiewski (2012)
!
! ==============================================================================
!
module cable_mpicommon

  use cable_def_types_mod

  implicit none

  private

  ! type
  public :: lpdecomp_t          ! type for landpoint decomposition
  ! functions
  public :: add_address_1block  ! passing arrays for worker
  public :: add_address_hvector ! passing arrays for master
  public :: decomp_types        ! derived type for landpt and patch arrays
  public :: bcast_start_time    ! broadcast day, month, ... of datetime
  public :: find_extents        ! byte size of kinds

  ! parameters
  ! MPI: base number of input fields sent to workers as start up parameters
  !   INTEGER, PARAMETER :: nparam = 68
  ! MPI: Bernard commented out two canopy params (potev_c and rwater)
  !   when porting to CABLE_r491 ! nparam = 219 -> nparam = 217
  ! MPI: CABLE_r491, after following up with Bernard on the new variables
  !   INTEGER, PARAMETER :: nparam = 260
  ! Added 23 params when trying to fix the bug in MPI ! nparam -> 283
  ! Add 10 variables to veg% param -> 293
  ! Ticket #56, add 2 new params for the Medlyns Stom Cond model 293 -> 295
  ! Vanessa Haverd: add 4 new params 295 -> 299
  ! VH add 9 params for sli 299 -> 308
  ! ? 3 extra params -> 311
  ! Matthias Cuntz: add 9 canopy params for 13C -> 320
  ! Paul Ryan: add 1 soil params for SLI -> 321
  ! Jurgen Knauer: add 10 veg params -> 331
  ! Matthias Cuntz: add 4 soil params for SLI -> 335
  ! Matthias Cuntz: send all parameters of types -> 418
  integer, parameter, public :: nparam = 418

  ! MPI: number of casa parameters sent to workers as start up parameters
  !   INTEGER, PARAMETER :: ncasaparam = 68
  !   INTEGER, PARAMETER :: ncasaparam = 176
  ! MPI: added casapool fields ratioNCsoilnew, ratioNCsoilmin and ratioNCsoilmax
  !   INTEGER, PARAMETER :: ncasaparam = 212  ! changed lpn added 9 variables
  !   (casaflux%frac_sapwood/sapwood_area,casabiome,casabiome%ratioNPplantmin,%ratioNPplantmax)
  !   casapool%ratioNPplant,%ratioNPlitter,ratioNPsoil
  ! vh added 3 variables la_tosa, vcmax_scalar, disturbance_interval
  !   INTEGER, PARAMETER :: ncasaparam = 215
  ! vh added 5 variables DAMM_EnzPool, DAMM_KMO2, DAMM_KMcp, DAMM_Ea, DAMM_alpha
  ! ? 2 extra params -> 222
  ! Matthias Cuntz: add 7 canopy params for 13C -> 229
  ! Matthias Cuntz: add 3 for fire -> 232
  integer, parameter, public :: ncasaparam = 232

  ! MPI: number of casa_dump parameters sent/rec'd to/from the workers every timestep
  ! Matthias Cuntz: add 2 casamet params for 13C -> 11
  ! Matthias Cuntz: calc as ncdumprw+icycle-1 and add 2 in case of 13C -> 8
  ! Vanessa Haverd: 12 extra climate dump variables (frec + 11 needed for BLAZE) -> 20
  ! Matthias Cuntz: add 11 in case of call_blaze -> 9
  integer, parameter, public :: ncdumprw = 9
  ! MPI: number of casa_LUC parameters sent/rec'd to/from the workers every year
  ! Matthias Cuntz: add 5 1D and 3 2D used in luc_casa_transfer -> 24
  integer, parameter, public :: nLUCrw = 24

  ! MPI: number of POP parameters sent/rec'd to/from the workers every
  ! timestep or at start, end. Here, with POP the dimensions are separate!
  integer, parameter, public :: npop = 988

  ! MPI: number of input fields sent to workers at the start of each timestep
  ! added 4 time fields in met: year, moy, doy, hod
  integer, parameter, public :: ninput = 18

  ! MPI: number of 3D array slices / worker (results)
  integer, parameter, public :: n3d = 1

  ! MPI: number of matrix slices / worker (results)
  !   INTEGER, PARAMETER :: nmat = 29
  ! MPI: 2011-07-08 - removed dtmlt from data exchange
  !   INTEGER, PARAMETER :: nmat = 28
  ! MPI: gol124: net +1 when Bernard ported to CABLE_r491
  !   INTEGER, PARAMETER :: nmat = 29
  ! MPI: CABLE_r491, after following up with Bernard on the new variables
  ! vh sli nmat + 4 36 -> 40
  ! Matthias Cuntz: add 8 2D canopy params for 13C -> 48
  ! Matthias Cuntz: add 1 2D ssnow params for SLI -> 49
  integer, parameter, public :: nmat = 49

  ! MPI: number of contig vector parts / worker (results)
  !   INTEGER, PARAMETER :: nvec = 149
  ! MPI: 2011-06-28 - removed ebal, ebal_tot, seb, seb_tot from data exchange
  !   INTEGER, PARAMETER :: nvec = 145
  ! MPI: 2011-07-08 - removed otss from data exchange
  !   INTEGER, PARAMETER :: nvec = 144
  ! MPI: 2012-02-14 - removed year, moy, doy, hod
  !   INTEGER, PARAMETER :: nvec = 140
  ! MPI: gol124: net -3 (removed or changed to 2D) when Bernard ported to CABLE_r491
  !   INTEGER, PARAMETER :: nvec = 137
  ! MPI: CABLE_r491, after following up with Bernard on the new variables
  ! vh sli nvec + 6 162 -> 168
  ! ? 15 extra params -> 183
  ! Matthias Cuntz: add 1 1D canopy param for 13C -> 184
  ! Matthias Cuntz: add 3 1D canopy params for elasticities -> 187
  integer, parameter, public :: nvec = 187

  ! MPI: number of fields included in restart_t type for data
  ! that is returned only for creating a restart file at the end of the run
  !   INTEGER, PARAMETER :: nrestart = 16
  ! MPI: gol124: canopy%rwater removed when Bernard ported to CABLE_r491 -> 15
  ! Matthias Cuntz: add 7 vars of canopy, veg, ssnow -> 22
  integer, parameter, public :: nrestart = 22
  integer, parameter, public :: nsumcasaflux = 62
  integer, parameter, public :: nsumcasapool = 40
  ! Matthias Cuntz: add 1 2D and 20 1D vars for restart file -> 72
  ! Paul Ryan: add 2 2D vars -> 74
  ! Matthias Cuntz: add 4 vars to pass all climate variables -> 78
  integer, parameter, public :: nclimate = 78

  ! 13C
  ! MPI: number of variables in c13o2_flux type
  integer, parameter, public :: nc13o2_flux = 10
  ! MPI: number of variables in c13o2_pool type
  integer, parameter, public :: nc13o2_pool = 5
  ! MPI: number of variables in c13o2_luc type
  integer, parameter, public :: nc13o2_luc = 3

  ! MPI: type to hold landpoint decomposition info
  type lpdecomp_t
     integer :: landp0      ! starting land point index
     integer :: nland       ! number of landpoints
     integer :: patch0      ! starting patch index in global CABLE vars
     integer :: npatch      ! sum of patches for all landpoints of this worker
     integer :: npop_iwood  ! number of pop-patches for each worker
     integer, allocatable :: iwood(:)  ! number of pop-patches for each worker
  end type lpdecomp_t

  ! MPI: worker's local landpoints and patches
  type(lpdecomp_t), public :: wpatch

  ! MPI: Fortran types extents
  integer, public :: extr1, extr2, extid, extl

  ! Routines for passing arrays in Fortran types with MPI.
  ! The master uses add_address_hvector and the worker receives
  ! with add_address_1block.
  ! add_address_hvector and add_address_1block are equivalent
  ! for 1D-arrays.

  ! add address, block length, and type of variable
  interface add_address_1block
     module procedure add_address_1block_1d_i1, add_address_1block_1d_l, &
          add_address_1block_1d_r1, add_address_1block_1d_r2, &
          add_address_1block_2d_i1, add_address_1block_2d_l, &
          add_address_1block_2d_r1, add_address_1block_2d_r2, &
          add_address_1block_3d_i1, add_address_1block_3d_l, &
          add_address_1block_3d_r1, add_address_1block_3d_r2
  end interface add_address_1block

  ! add address, block length, and type of variable,
  ! creating MPI types for arrays
  interface add_address_hvector
     module procedure add_address_hvector_1d_i1, add_address_hvector_1d_l, &
          add_address_hvector_1d_r1, add_address_hvector_1d_r2, &
          add_address_hvector_2d_i1, add_address_hvector_2d_l, &
          add_address_hvector_2d_r1, add_address_hvector_2d_r2, &
          add_address_hvector_3d_i1, add_address_hvector_3d_l, &
          add_address_hvector_3d_r1, add_address_hvector_3d_r2
  end interface add_address_hvector

contains

  ! calculates extents of the Fortran types used by CABLE
  subroutine find_extents()

    use mpi
    use cable_def_types_mod

    implicit none

    integer,   dimension(2) :: itmp
    real,      dimension(2) :: r1tmp
    real(r_2), dimension(2) :: r2tmp
    logical,   dimension(2) :: ltmp

    integer(KIND=MPI_ADDRESS_KIND), dimension(2) :: a

    integer :: ierr

    call MPI_Get_address(itmp(1), a(1), ierr)
    call MPI_Get_address(itmp(2), a(2), ierr)
    extid = int(a(2)-a(1))

    call MPI_Get_address(r1tmp(1), a(1), ierr)
    call MPI_Get_address(r1tmp(2), a(2), ierr)
    extr1 = int(a(2)-a(1))

    call MPI_Get_address(r2tmp(1), a(1), ierr)
    call MPI_Get_address(r2tmp(2), a(2), ierr)
    extr2 = int(a(2)-a(1))

    call MPI_Get_address(ltmp(1), a(1), ierr)
    call MPI_Get_address(ltmp(2), a(2), ierr)
    extl = int(a(2)-a(1))

  end subroutine find_extents


  ! creates MPI derived datatypes for exchanging landpt and patch arrays
  subroutine decomp_types(landpt_t, patch_t)

    use mpi
    use cable_IO_vars_module

    implicit none

    integer, intent(OUT) :: landpt_t, patch_t

    ! dummy vars to calculate field offsets
    type(land_type)  :: dlandpt(2)
    type(patch_type) :: dpatch(2)

    integer(KIND=MPI_ADDRESS_KIND) :: base_d, el2, text

    integer, parameter :: fields = 5
    integer, dimension(fields) :: blocks, types
    integer(KIND=MPI_ADDRESS_KIND), dimension(fields) :: displs

    ! temp variable for lower bound parameter when setting extent
    integer(KIND=MPI_ADDRESS_KIND) :: lb

    integer :: tmp_t, ierr

    lb = 0
    blocks = 1

    ! create MPI type to exchange landpt records
    types = MPI_INTEGER

    call MPI_Get_address(dlandpt(1), base_d, ierr)

    call MPI_Get_address(dlandpt(1)%nap, displs(1), ierr)
    call MPI_Get_address(dlandpt(1)%cstart, displs(2), ierr)
    call MPI_Get_address(dlandpt(1)%cend, displs(3), ierr)
    call MPI_Get_address(dlandpt(1)%ilat, displs(4), ierr)
    call MPI_Get_address(dlandpt(1)%ilon, displs(5), ierr)

    displs = displs - base_d

    call MPI_Type_create_struct(5, blocks, displs, types, tmp_t, ierr)
    call MPI_Type_commit(tmp_t, ierr)

    ! make sure the type has correct extent for use in arrays
    call MPI_Get_Address(dlandpt(2), el2, ierr)
    text = el2 - base_d
    call MPI_Type_create_resized(tmp_t, lb, text, landpt_t, ierr)
    call MPI_Type_commit(landpt_t, ierr)

    ! create MPI type to exchange patch records
    types = MPI_REAL
    types(1) = MPI_DOUBLE_PRECISION

    call MPI_Get_address(dpatch(1), base_d, ierr)

    call MPI_Get_address(dpatch(1)%frac, displs(1), ierr)
    call MPI_Get_address(dpatch(1)%latitude, displs(2), ierr)
    call MPI_Get_address(dpatch(1)%longitude, displs(3), ierr)

    displs = displs - base_d

    call MPI_Type_create_struct(3, blocks, displs, types, tmp_t, ierr)
    call MPI_Type_commit(tmp_t, ierr)

    ! make sure the type has correct extent for use in arrays
    call MPI_Get_Address(dpatch(2), el2, ierr)
    text = el2 - base_d
    call MPI_Type_create_resized(tmp_t, lb, text, patch_t, ierr)
    call MPI_Type_commit(patch_t, ierr)

    return

  end subroutine decomp_types


  subroutine bcast_start_time(comm)

    use mpi
    use cable_IO_vars_module

    implicit none

    integer, intent(IN) :: comm

    integer :: ierr

    call MPI_Bcast(shod, 1, MPI_REAL, 0, comm, ierr)
    call MPI_Bcast(sdoy, 1, MPI_INTEGER, 0, comm, ierr)
    call MPI_Bcast(smoy, 1, MPI_INTEGER, 0, comm, ierr)
    call MPI_Bcast(syear, 1, MPI_INTEGER, 0, comm, ierr)
    call MPI_Bcast(time_coord, 3, MPI_CHARACTER, 0, comm, ierr)

    return

  end subroutine bcast_start_time


  ! ------------------------------------------------------------------


  ! Routines for passing arrays in Fortran types with MPI.
  ! The master uses add_address_hvector and the worker receives
  ! with add_address_1block.
  ! add_address_hvector and add_address_1block are equivalent
  ! for 1D-arrays.

  subroutine add_address_1block_1d_i1(var, off, cnt, addr, blen, typ, idx)

#ifndef __INTEL__
    use mpi, only: MPI_Address_kind, MPI_Get_address, MPI_Byte
#else
    use mpi
#endif

    implicit none

    ! variable to add from offset *off* with length *cnt*
    integer, dimension(:), intent(in) :: var
    integer,               intent(in) :: off
    integer,               intent(in) :: cnt
    ! MPI adress
    integer(MPI_Address_kind), dimension(:), intent(inout) :: addr
    ! block length
    integer,                   dimension(:), intent(inout) :: blen
    ! MPI data type
    integer,                   dimension(:), intent(inout) :: typ
    ! index in *addr*, *blen*, *typ*
    integer, intent(inout) :: idx

    integer :: ierr

    idx = idx + 1
    call MPI_Get_address(var(off), addr(idx), ierr)
    blen(idx) = cnt * extid
    typ(idx)  = MPI_Byte

  end subroutine add_address_1block_1d_i1


  subroutine add_address_1block_1d_l(var, off, cnt, addr, blen, typ, idx)

#ifndef __INTEL__
    use mpi, only: MPI_Address_kind, MPI_Get_address, MPI_Byte
#else
    use mpi
#endif

    implicit none

    ! variable to add from offset *off* with length *cnt*
    logical, dimension(:), intent(in) :: var
    integer,               intent(in) :: off
    integer,               intent(in) :: cnt
    ! MPI adress
    integer(MPI_Address_kind), dimension(:), intent(inout) :: addr
    ! block length
    integer,                   dimension(:), intent(inout) :: blen
    ! MPI data type
    integer,                   dimension(:), intent(inout) :: typ
    ! index in *addr*, *blen*, *typ*
    integer, intent(inout) :: idx

    integer :: ierr

    idx = idx + 1
    call MPI_Get_address(var(off), addr(idx), ierr)
    blen(idx) = cnt * extl
    typ(idx)  = MPI_Byte

  end subroutine add_address_1block_1d_l


  subroutine add_address_1block_1d_r1(var, off, cnt, addr, blen, typ, idx)

#ifndef __INTEL__
    use mpi, only: MPI_Address_kind, MPI_Get_address, MPI_Byte
#else
    use mpi
#endif

    implicit none

    ! variable to add from offset *off* with length *cnt*
    real, dimension(:), intent(in) :: var
    integer,            intent(in) :: off
    integer,            intent(in) :: cnt
    ! MPI adress
    integer(MPI_Address_kind), dimension(:), intent(inout) :: addr
    ! block length
    integer,                   dimension(:), intent(inout) :: blen
    ! MPI data type
    integer,                   dimension(:), intent(inout) :: typ
    ! index in *addr*, *blen*, *typ*
    integer, intent(inout) :: idx

    integer :: ierr

    idx = idx + 1
    call MPI_Get_address(var(off), addr(idx), ierr)
    blen(idx) = cnt * extr1
    typ(idx)  = MPI_Byte

  end subroutine add_address_1block_1d_r1


  subroutine add_address_1block_1d_r2(var, off, cnt, addr, blen, typ, idx)

#ifndef __INTEL__
    use mpi, only: MPI_Address_kind, MPI_Get_address, MPI_Byte
#else
    use mpi
#endif
    use cable_def_types_mod, only: r_2

    implicit none

    ! variable to add from offset *off* with length *cnt*
    real(r_2), dimension(:), intent(in) :: var
    integer,                 intent(in) :: off
    integer,                 intent(in) :: cnt
    ! MPI adress
    integer(MPI_Address_kind), dimension(:), intent(inout) :: addr
    ! block length
    integer,                   dimension(:), intent(inout) :: blen
    ! MPI data type
    integer,                   dimension(:), intent(inout) :: typ
    ! index in *addr*, *blen*, *typ*
    integer, intent(inout) :: idx

    integer :: ierr

    idx = idx + 1
    call MPI_Get_address(var(off), addr(idx), ierr)
    blen(idx) = cnt * extr2
    typ(idx)  = MPI_Byte

  end subroutine add_address_1block_1d_r2


  subroutine add_address_1block_2d_i1(var, off, cnt, addr, blen, typ, idx)

#ifndef __INTEL__
    use mpi, only: MPI_Address_kind, MPI_Get_address, MPI_Byte
#else
    use mpi
#endif

    implicit none

    ! variable to add from offset *off* with length *cnt*
    integer, dimension(:,:), intent(in) :: var
    integer,                 intent(in) :: off
    integer,                 intent(in) :: cnt
    ! MPI adress
    integer(MPI_Address_kind), dimension(:), intent(inout) :: addr
    ! block length
    integer,                   dimension(:), intent(inout) :: blen
    ! MPI data type
    integer,                   dimension(:), intent(inout) :: typ
    ! index in *addr*, *blen*, *typ*
    integer, intent(inout) :: idx

    integer :: ierr

    idx = idx + 1
    call MPI_Get_address(var(off, 1), addr(idx), ierr)
    blen(idx) = size(var, 2) * cnt * extid
    typ(idx)  = MPI_Byte

  end subroutine add_address_1block_2d_i1


  subroutine add_address_1block_2d_l(var, off, cnt, addr, blen, typ, idx)

#ifndef __INTEL__
    use mpi, only: MPI_Address_kind, MPI_Get_address, MPI_Byte
#else
    use mpi
#endif

    implicit none

    ! variable to add from offset *off* with length *cnt*
    logical, dimension(:,:), intent(in) :: var
    integer,                 intent(in) :: off
    integer,                 intent(in) :: cnt
    ! MPI adress
    integer(MPI_Address_kind), dimension(:), intent(inout) :: addr
    ! block length
    integer,                   dimension(:), intent(inout) :: blen
    ! MPI data type
    integer,                   dimension(:), intent(inout) :: typ
    ! index in *addr*, *blen*, *typ*
    integer, intent(inout) :: idx

    integer :: ierr

    idx = idx + 1
    call MPI_Get_address(var(off, 1), addr(idx), ierr)
    blen(idx) = size(var, 2) * cnt * extl
    typ(idx)  = MPI_Byte

  end subroutine add_address_1block_2d_l


  subroutine add_address_1block_2d_r1(var, off, cnt, addr, blen, typ, idx)

#ifndef __INTEL__
    use mpi, only: MPI_Address_kind, MPI_Get_address, MPI_Byte
#else
    use mpi
#endif

    implicit none

    ! variable to add from offset *off* with length *cnt*
    real, dimension(:,:), intent(in) :: var
    integer,              intent(in) :: off
    integer,              intent(in) :: cnt
    ! MPI adress
    integer(MPI_Address_kind), dimension(:), intent(inout) :: addr
    ! block length
    integer,                   dimension(:), intent(inout) :: blen
    ! MPI data type
    integer,                   dimension(:), intent(inout) :: typ
    ! index in *addr*, *blen*, *typ*
    integer, intent(inout) :: idx

    integer :: ierr

    idx = idx + 1
    call MPI_Get_address(var(off, 1), addr(idx), ierr)
    blen(idx) = size(var, 2) * cnt * extr1
    typ(idx)  = MPI_Byte

  end subroutine add_address_1block_2d_r1


  subroutine add_address_1block_2d_r2(var, off, cnt, addr, blen, typ, idx)

#ifndef __INTEL__
    use mpi, only: MPI_Address_kind, MPI_Get_address, MPI_Byte
#else
    use mpi
#endif
    use cable_def_types_mod, only: r_2

    implicit none

    ! variable to add from offset *off* with length *cnt*
    real(r_2), dimension(:,:), intent(in) :: var
    integer,                   intent(in) :: off
    integer,                   intent(in) :: cnt
    ! MPI adress
    integer(MPI_Address_kind), dimension(:), intent(inout) :: addr
    ! block length
    integer,                   dimension(:), intent(inout) :: blen
    ! MPI data type
    integer,                   dimension(:), intent(inout) :: typ
    ! index in *addr*, *blen*, *typ*
    integer, intent(inout) :: idx

    integer :: ierr

    idx = idx + 1
    call MPI_Get_address(var(off, 1), addr(idx), ierr)
    blen(idx) = size(var, 2) * cnt * extr2
    typ(idx)  = MPI_Byte

  end subroutine add_address_1block_2d_r2


  subroutine add_address_1block_3d_i1(var, off, cnt, addr, blen, typ, idx)

#ifndef __INTEL__
    use mpi, only: MPI_Address_kind, MPI_Get_address, MPI_Byte
#else
    use mpi
#endif

    implicit none

    ! variable to add from offset *off* with length *cnt*
    integer, dimension(:,:,:), intent(in) :: var
    integer,                   intent(in) :: off
    integer,                   intent(in) :: cnt
    ! MPI adress
    integer(MPI_Address_kind), dimension(:), intent(inout) :: addr
    ! block length
    integer,                   dimension(:), intent(inout) :: blen
    ! MPI data type
    integer,                   dimension(:), intent(inout) :: typ
    ! index in *addr*, *blen*, *typ*
    integer, intent(inout) :: idx

    integer :: ierr

    idx = idx + 1
    call MPI_Get_address(var(off, 1, 1), addr(idx), ierr)
    blen(idx) = size(var, 3) * size(var, 2) * cnt * extid
    typ(idx)  = MPI_Byte

  end subroutine add_address_1block_3d_i1


  subroutine add_address_1block_3d_l(var, off, cnt, addr, blen, typ, idx)

#ifndef __INTEL__
    use mpi, only: MPI_Address_kind, MPI_Get_address, MPI_Byte
#else
    use mpi
#endif

    implicit none

    ! variable to add from offset *off* with length *cnt*
    logical, dimension(:,:,:), intent(in) :: var
    integer,                   intent(in) :: off
    integer,                   intent(in) :: cnt
    ! MPI adress
    integer(MPI_Address_kind), dimension(:), intent(inout) :: addr
    ! block length
    integer,                   dimension(:), intent(inout) :: blen
    ! MPI data type
    integer,                   dimension(:), intent(inout) :: typ
    ! index in *addr*, *blen*, *typ*
    integer, intent(inout) :: idx

    integer :: ierr

    idx = idx + 1
    call MPI_Get_address(var(off, 1, 1), addr(idx), ierr)
    blen(idx) = size(var, 3) * size(var, 2) * cnt * extl
    typ(idx)  = MPI_Byte

  end subroutine add_address_1block_3d_l


  subroutine add_address_1block_3d_r1(var, off, cnt, addr, blen, typ, idx)

#ifndef __INTEL__
    use mpi, only: MPI_Address_kind, MPI_Get_address, MPI_Byte
#else
    use mpi
#endif

    implicit none

    ! variable to add from offset *off* with length *cnt*
    real, dimension(:,:,:), intent(in) :: var
    integer,                intent(in) :: off
    integer,                intent(in) :: cnt
    ! MPI adress
    integer(MPI_Address_kind), dimension(:), intent(inout) :: addr
    ! block length
    integer,                   dimension(:), intent(inout) :: blen
    ! MPI data type
    integer,                   dimension(:), intent(inout) :: typ
    ! index in *addr*, *blen*, *typ*
    integer, intent(inout) :: idx

    integer :: ierr

    idx = idx + 1
    call MPI_Get_address(var(off, 1, 1), addr(idx), ierr)
    blen(idx) = size(var, 3) * size(var, 2) * cnt * extr1
    typ(idx)  = MPI_Byte

  end subroutine add_address_1block_3d_r1


  subroutine add_address_1block_3d_r2(var, off, cnt, addr, blen, typ, idx)

#ifndef __INTEL__
    use mpi, only: MPI_Address_kind, MPI_Get_address, MPI_Byte
#else
    use mpi
#endif
    use cable_def_types_mod, only: r_2

    implicit none

    ! variable to add from offset *off* with length *cnt*
    real(r_2), dimension(:,:,:), intent(in) :: var
    integer,                     intent(in) :: off
    integer,                     intent(in) :: cnt
    ! MPI adress
    integer(MPI_Address_kind), dimension(:), intent(inout) :: addr
    ! block length
    integer,                   dimension(:), intent(inout) :: blen
    ! MPI data type
    integer,                   dimension(:), intent(inout) :: typ
    ! index in *addr*, *blen*, *typ*
    integer, intent(inout) :: idx

    integer :: ierr

    idx = idx + 1
    call MPI_Get_address(var(off, 1, 1), addr(idx), ierr)
    blen(idx) = size(var, 3) * size(var, 2) * cnt * extr2
    typ(idx)  = MPI_Byte

  end subroutine add_address_1block_3d_r2


  ! ------------------------------------------------------------------


  subroutine add_address_hvector_1d_i1(var, off, cnt, stride, &
       addr, blen, typ, idx)

#ifndef __INTEL__
    use mpi, only: MPI_Address_kind, MPI_Get_address, MPI_Byte
#else
    use mpi
#endif

    implicit none

    ! variable to add from offset *off* with length *cnt*
    integer, dimension(:), intent(in) :: var
    integer,               intent(in) :: off
    integer,               intent(in) :: cnt
    ! length of 1st dim of original array
    integer,               intent(in) :: stride  ! not used
    ! MPI adress
    integer(MPI_Address_kind), dimension(:), intent(inout) :: addr
    ! block length
    integer,                   dimension(:), intent(inout) :: blen
    ! MPI data type
    integer,                   dimension(:), intent(inout) :: typ
    ! index in *addr*, *blen*, *typ*
    integer, intent(inout) :: idx

    integer :: ierr

    idx = idx + 1
    call MPI_Get_address(var(off), addr(idx), ierr)
    blen(idx) = cnt * extid
    typ(idx)  = MPI_Byte

  end subroutine add_address_hvector_1d_i1


  subroutine add_address_hvector_1d_l(var, off, cnt, stride, &
       addr, blen, typ, idx)

#ifndef __INTEL__
    use mpi, only: MPI_Address_kind, MPI_Get_address, MPI_Byte
#else
    use mpi
#endif

    implicit none

    ! variable to add from offset *off* with length *cnt*
    logical, dimension(:), intent(in) :: var
    integer,               intent(in) :: off
    integer,               intent(in) :: cnt
    ! length of 1st dim of original array
    integer,               intent(in) :: stride  ! not used
    ! MPI adress
    integer(MPI_Address_kind), dimension(:), intent(inout) :: addr
    ! block length
    integer,                   dimension(:), intent(inout) :: blen
    ! MPI data type
    integer,                   dimension(:), intent(inout) :: typ
    ! index in *addr*, *blen*, *typ*
    integer, intent(inout) :: idx

    integer :: ierr

    idx = idx + 1
    call MPI_Get_address(var(off), addr(idx), ierr)
    blen(idx) = cnt * extl
    typ(idx)  = MPI_Byte

  end subroutine add_address_hvector_1d_l


  subroutine add_address_hvector_1d_r1(var, off, cnt, stride, &
       addr, blen, typ, idx)

#ifndef __INTEL__
    use mpi, only: MPI_Address_kind, MPI_Get_address, MPI_Byte
#else
    use mpi
#endif

    implicit none

    ! variable to add from offset *off* with length *cnt*
    real, dimension(:), intent(in) :: var
    integer,            intent(in) :: off
    integer,            intent(in) :: cnt
    ! length of 1st dim of original array
    integer,            intent(in) :: stride  ! not used
    ! MPI adress
    integer(MPI_Address_kind), dimension(:), intent(inout) :: addr
    ! block length
    integer,                   dimension(:), intent(inout) :: blen
    ! MPI data type
    integer,                   dimension(:), intent(inout) :: typ
    ! index in *addr*, *blen*, *typ*
    integer, intent(inout) :: idx

    integer :: ierr

    idx = idx + 1
    call MPI_Get_address(var(off), addr(idx), ierr)
    blen(idx) = cnt * extr1
    typ(idx)  = MPI_Byte

  end subroutine add_address_hvector_1d_r1


  subroutine add_address_hvector_1d_r2(var, off, cnt, stride, &
       addr, blen, typ, idx)

#ifndef __INTEL__
    use mpi, only: MPI_Address_kind, MPI_Get_address, MPI_Byte
#else
    use mpi
#endif
    use cable_def_types_mod, only: r_2

    implicit none

    ! variable to add from offset *off* with length *cnt*
    real(r_2), dimension(:), intent(in) :: var
    integer,                 intent(in) :: off
    integer,                 intent(in) :: cnt
    ! length of 1st dim of original array
    integer,                 intent(in) :: stride  ! not used
    ! MPI adress
    integer(MPI_Address_kind), dimension(:), intent(inout) :: addr
    ! block length
    integer,                   dimension(:), intent(inout) :: blen
    ! MPI data type
    integer,                   dimension(:), intent(inout) :: typ
    ! index in *addr*, *blen*, *typ*
    integer, intent(inout) :: idx

    integer :: ierr

    idx = idx + 1
    call MPI_Get_address(var(off), addr(idx), ierr)
    blen(idx) = cnt * extr2
    typ(idx)  = MPI_Byte

  end subroutine add_address_hvector_1d_r2


  subroutine add_address_hvector_2d_i1(var, off, cnt, stride, &
       addr, blen, typ, idx)

#ifndef __INTEL__
    use mpi, only: MPI_Address_kind, MPI_Get_address, &
         MPI_Type_create_hvector, MPI_Byte
#else
    use mpi
#endif

    implicit none

    ! variable to add from offset *off* with length *cnt*
    integer, dimension(:,:), intent(in) :: var
    integer,                 intent(in) :: off
    integer,                 intent(in) :: cnt
    ! length of 1st dim of original array
    integer,                 intent(in) :: stride
    ! MPI adress
    integer(MPI_Address_kind), dimension(:), intent(inout) :: addr
    ! block length
    integer,                   dimension(:), intent(inout) :: blen
    ! MPI data type
    integer,                   dimension(:), intent(inout) :: typ
    ! index in *addr*, *blen*, *typ*
    integer, intent(inout) :: idx

    integer(MPI_Address_kind) :: istride
    integer :: ierr

    idx = idx + 1
    call MPI_Get_address(var(off, 1), addr(idx), ierr)
    istride = stride * extid
    call MPI_Type_create_hvector(size(var, 2), cnt * extid, istride, &
         MPI_BYTE, typ(idx), ierr)
    blen(idx) = 1

  end subroutine add_address_hvector_2d_i1


  subroutine add_address_hvector_2d_l(var, off, cnt, stride, &
       addr, blen, typ, idx)

#ifndef __INTEL__
    use mpi, only: MPI_Address_kind, MPI_Get_address, &
         MPI_Type_create_hvector, MPI_Byte
#else
    use mpi
#endif

    implicit none

    ! variable to add from offset *off* with length *cnt*
    logical, dimension(:,:), intent(in) :: var
    integer,                 intent(in) :: off
    integer,                 intent(in) :: cnt
    ! length of 1st dim of original array
    integer,                 intent(in) :: stride
    ! MPI adress
    integer(MPI_Address_kind), dimension(:), intent(inout) :: addr
    ! block length
    integer,                   dimension(:), intent(inout) :: blen
    ! MPI data type
    integer,                   dimension(:), intent(inout) :: typ
    ! index in *addr*, *blen*, *typ*
    integer, intent(inout) :: idx

    integer(MPI_Address_kind) :: istride
    integer :: ierr

    idx = idx + 1
    call MPI_Get_address(var(off, 1), addr(idx), ierr)
    istride = stride * extl
    call MPI_Type_create_hvector(size(var, 2), cnt * extl, istride, &
         MPI_BYTE, typ(idx), ierr)
    blen(idx) = 1

  end subroutine add_address_hvector_2d_l


  subroutine add_address_hvector_2d_r1(var, off, cnt, stride, &
       addr, blen, typ, idx)

#ifndef __INTEL__
    use mpi, only: MPI_Address_kind, MPI_Get_address, &
         MPI_Type_create_hvector, MPI_Byte
#else
    use mpi
#endif

    implicit none

    ! variable to add from offset *off* with length *cnt*
    real, dimension(:,:), intent(in) :: var
    integer,              intent(in) :: off
    integer,              intent(in) :: cnt
    ! length of 1st dim of original array
    integer,              intent(in) :: stride
    ! MPI adress
    integer(MPI_Address_kind), dimension(:), intent(inout) :: addr
    ! block length
    integer,                   dimension(:), intent(inout) :: blen
    ! MPI data type
    integer,                   dimension(:), intent(inout) :: typ
    ! index in *addr*, *blen*, *typ*
    integer, intent(inout) :: idx

    integer(MPI_Address_kind) :: istride
    integer :: ierr

    idx = idx + 1
    call MPI_Get_address(var(off, 1), addr(idx), ierr)
    istride = stride * extr1
    call MPI_Type_create_hvector(size(var, 2), cnt * extr1, istride, &
         MPI_BYTE, typ(idx), ierr)
    blen(idx) = 1

  end subroutine add_address_hvector_2d_r1


  subroutine add_address_hvector_2d_r2(var, off, cnt, stride, &
       addr, blen, typ, idx)

#ifndef __INTEL__
    use mpi, only: MPI_Address_kind, MPI_Get_address, &
         MPI_Type_create_hvector, MPI_Byte
#else
    use mpi
#endif
    use cable_def_types_mod, only: r_2

    implicit none

    ! variable to add from offset *off* with length *cnt*
    real(r_2), dimension(:,:), intent(in) :: var
    integer,                   intent(in) :: off
    integer,                   intent(in) :: cnt
    ! length of 1st dim of original array
    integer,                   intent(in) :: stride
    ! MPI adress
    integer(MPI_Address_kind), dimension(:), intent(inout) :: addr
    ! block length
    integer,                   dimension(:), intent(inout) :: blen
    ! MPI data type
    integer,                   dimension(:), intent(inout) :: typ
    ! index in *addr*, *blen*, *typ*
    integer, intent(inout) :: idx

    integer(MPI_Address_kind) :: istride
    integer :: ierr

    idx = idx + 1
    call MPI_Get_address(var(off, 1), addr(idx), ierr)
    istride = stride * extr2
    call MPI_Type_create_hvector(size(var, 2), cnt * extr2, istride, &
         MPI_BYTE, typ(idx), ierr)
    blen(idx) = 1

  end subroutine add_address_hvector_2d_r2


  subroutine add_address_hvector_3d_i1(var, off, cnt, stride, &
       addr, blen, typ, idx)

#ifndef __INTEL__
    use mpi, only: MPI_Address_kind, MPI_Get_address, &
         MPI_Type_create_hvector, MPI_Byte
#else
    use mpi
#endif

    implicit none

    ! variable to add from offset *off* with length *cnt*
    integer, dimension(:,:,:), intent(in) :: var
    integer,                   intent(in) :: off
    integer,                   intent(in) :: cnt
    ! length of 1st dim of original array
    integer,                   intent(in) :: stride
    ! MPI adress
    integer(MPI_Address_kind), dimension(:), intent(inout) :: addr
    ! block length
    integer,                   dimension(:), intent(inout) :: blen
    ! MPI data type
    integer,                   dimension(:), intent(inout) :: typ
    ! index in *addr*, *blen*, *typ*
    integer, intent(inout) :: idx

    integer(MPI_Address_kind) :: istride
    integer :: ierr

    idx = idx + 1
    call MPI_Get_address(var(off, 1, 1), addr(idx), ierr)
    istride = stride * extid
    call MPI_Type_create_hvector(size(var, 3) * size(var, 2), cnt * extid, &
         istride, MPI_BYTE, typ(idx), ierr)
    blen(idx) = 1

  end subroutine add_address_hvector_3d_i1


  subroutine add_address_hvector_3d_l(var, off, cnt, stride, &
       addr, blen, typ, idx)

#ifndef __INTEL__
    use mpi, only: MPI_Address_kind, MPI_Get_address, &
         MPI_Type_create_hvector, MPI_Byte
#else
    use mpi
#endif

    implicit none

    ! variable to add from offset *off* with length *cnt*
    logical, dimension(:,:,:), intent(in) :: var
    integer,                   intent(in) :: off
    integer,                   intent(in) :: cnt
    ! length of 1st dim of original array
    integer,                   intent(in) :: stride
    ! MPI adress
    integer(MPI_Address_kind), dimension(:), intent(inout) :: addr
    ! block length
    integer,                   dimension(:), intent(inout) :: blen
    ! MPI data type
    integer,                   dimension(:), intent(inout) :: typ
    ! index in *addr*, *blen*, *typ*
    integer, intent(inout) :: idx

    integer(MPI_Address_kind) :: istride
    integer :: ierr

    idx = idx + 1
    call MPI_Get_address(var(off, 1, 1), addr(idx), ierr)
    istride = stride * extl
    call MPI_Type_create_hvector(size(var, 3) * size(var, 2), cnt * extl, &
         istride, MPI_BYTE, typ(idx), ierr)
    blen(idx) = 1

  end subroutine add_address_hvector_3d_l


  subroutine add_address_hvector_3d_r1(var, off, cnt, stride, &
       addr, blen, typ, idx)

#ifndef __INTEL__
    use mpi, only: MPI_Address_kind, MPI_Get_address, &
         MPI_Type_create_hvector, MPI_Byte
#else
    use mpi
#endif

    implicit none

    ! variable to add from offset *off* with length *cnt*
    real, dimension(:,:,:), intent(in) :: var
    integer,                intent(in) :: off
    integer,                intent(in) :: cnt
    ! length of 1st dim of original array
    integer,                intent(in) :: stride
    ! MPI adress
    integer(MPI_Address_kind), dimension(:), intent(inout) :: addr
    ! block length
    integer,                   dimension(:), intent(inout) :: blen
    ! MPI data type
    integer,                   dimension(:), intent(inout) :: typ
    ! index in *addr*, *blen*, *typ*
    integer, intent(inout) :: idx

    integer(MPI_Address_kind) :: istride
    integer :: ierr

    idx = idx + 1
    call MPI_Get_address(var(off, 1, 1), addr(idx), ierr)
    istride = stride * extr1
    call MPI_Type_create_hvector(size(var, 3) * size(var, 2), cnt * extr1, &
         istride, MPI_BYTE, typ(idx), ierr)
    blen(idx) = 1

  end subroutine add_address_hvector_3d_r1


  subroutine add_address_hvector_3d_r2(var, off, cnt, stride, &
       addr, blen, typ, idx)

#ifndef __INTEL__
    use mpi, only: MPI_Address_kind, MPI_Get_address, &
         MPI_Type_create_hvector, MPI_Byte
#else
    use mpi
#endif
    use cable_def_types_mod, only: r_2

    implicit none

    ! variable to add from offset *off* with length *cnt*
    real(r_2), dimension(:,:,:), intent(in) :: var
    integer,                     intent(in) :: off
    integer,                     intent(in) :: cnt
    ! length of 1st dim of original array
    integer,                     intent(in) :: stride
    ! MPI adress
    integer(MPI_Address_kind), dimension(:), intent(inout) :: addr
    ! block length
    integer,                   dimension(:), intent(inout) :: blen
    ! MPI data type
    integer,                   dimension(:), intent(inout) :: typ
    ! index in *addr*, *blen*, *typ*
    integer, intent(inout) :: idx

    integer(MPI_Address_kind) :: istride
    integer :: ierr

    idx = idx + 1
    call MPI_Get_address(var(off, 1, 1), addr(idx), ierr)
    istride = stride * extr2
    call MPI_Type_create_hvector(size(var, 3) * size(var, 2), cnt * extr2, &
         istride, MPI_BYTE, typ(idx), ierr)
    blen(idx) = 1

  end subroutine add_address_hvector_3d_r2

end module cable_mpicommon
