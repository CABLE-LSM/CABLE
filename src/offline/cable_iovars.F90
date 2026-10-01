!==============================================================================
! This source code is part of the
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CSIRO Open Source Software License
! Agreement (variation of the BSD / MIT License).
!
! You may not use this file except in compliance with this License.
! A copy of the License (CSIRO_BSD_MIT_License_v2.0_CABLE.txt) is located
! in each directory cTYPE(casa_flux_type), INTENT(IN) :: casaflux ! casa fluxesontaining CABLE code.
!
! ==============================================================================
! Purpose: Defines input/output related variables for CABLE offline
!
! Contact: Bernard.Pak@csiro.au
!
! History: Development by Gab Abramowitz
!          Additional code to use multiple vegetation types per grid-cell (patches)
!
! ==============================================================================
MODULE cable_IO_vars_module

  USE cable_def_types_mod, ONLY : r_2, mvtype, mstype

  IMPLICIT NONE

  PUBLIC
  PRIVATE r_2, mvtype, mstype

  ! ============ Timing variables =====================
  REAL :: shod ! start time hour-of-day

  INTEGER :: sdoy,smoy,syear ! start time day-of-year month and year

  CHARACTER(LEN=200) :: timeunits ! timing info read from nc file

  CHARACTER(LEN=10) :: calendar ! 'standard' if using leap years (set by 
                                ! leaps namelist option), else 'noleap'

  CHARACTER(LEN=3) :: time_coord ! GMT or LOCal time variables

  REAL(r_2),POINTER,DIMENSION(:) :: timevar ! time variable from file

  INTEGER,DIMENSION(12) ::                                                    &
       daysm = (/31,28,31,30,31,30,31,31,30,31,30,31/),                         &
       daysml = (/31,29,31,30,31,30,31,31,30,31,30,31/),                        &
       lastday = (/31,59,90,120,151,181,212,243,273,304,334,365/),              &
       lastdayl = (/31,60,91,121,152,182,213,244,274,305,335,366/)

  LOGICAL :: leaps   ! use leap year timing?

  ! ============ Structure variables ===================
  REAL, POINTER,DIMENSION(:) :: latitude, longitude

  REAL,POINTER, DIMENSION(:,:) :: lat_all, lon_all ! lat and lon

  CHARACTER(LEN=4) :: metGrid ! Either 'land' or 'mask'

  INTEGER,POINTER,DIMENSION(:,:) :: mask ! land/sea mask from met file

  INTEGER, POINTER, DIMENSION(:) :: land_x
    !* An array mapping local land indexes on the current MPI rank to global x
    ! (longitude) indexes.
  INTEGER, POINTER, DIMENSION(:) :: land_y
    !* An array mapping local land indexes on the current MPI rank to global y
    ! (latitude) indexes.

  INTEGER, POINTER, DIMENSION(:) :: land_x_global
    !! An array mapping global land indexes to global x (longitude) indexes.
  INTEGER, POINTER, DIMENSION(:) :: land_y_global
    !! An array mapping global land indexes to global y (latitude) indexes.

  INTEGER ::                                                                  &
       xdimsize,ydimsize,   & ! sizes of x and y dimensions
       ngridcells             ! number of gridcells in simulation

  ! For vegetated surface type
  TYPE patch_type

     REAL ::                                                                  &
          frac,    &  ! fractional cover of each veg patch
          latitude,&
          longitude

  END TYPE patch_type


  TYPE land_type

     INTEGER ::                                                               &
          nap,     & ! number of active (>0%) patches (<=max_vegpatches)
          cstart,  & ! pos of 1st gridcell veg patch in main arrays
          cend,    & ! pos of last gridcell veg patch in main arrays
          ilat,    & ! replacing land_y  ! ??
          ilon       ! replacing land_x  ! ??

  END TYPE land_type


  TYPE(land_type), DIMENSION(:), POINTER :: landpt
    !! Land information for each land point in the local grid of this MPI rank
  TYPE(land_type), DIMENSION(:), POINTER :: landpt_global
    !! Land information for each land point in the global grid
  TYPE(patch_type), DIMENSION(:), POINTER :: patch

  INTEGER ::                                                                  &
       max_vegpatches,   & ! The maximum # of patches in any grid cell
       nmetpatches         ! size of patch dimension in met file, if exists

  INTEGER :: land_decomp_start
    !! Starting land point index of this MPI rank in global land array
  INTEGER :: land_decomp_end
    !! Ending land point index of this MPI rank in global land array
  INTEGER :: patch_decomp_start
    !! Starting patch index of this MPI rank in global patch array
  INTEGER :: patch_decomp_end
    !! Ending patch index of this MPI rank in global patch array

  ! =============== File details ==========================
   TYPE globalMet_type
     LOGICAL           ::                                                     &
       l_gpcc,&! = .FALSE., &         ! ypwang following Chris Lu (30/oct/2012)
       l_gswp,&!= .FALSE. , &         ! BP May 2013
       l_ncar,&! = .FALSE., &         ! BP Dec 2013
       l_access ! = .FALSE.          ! BP May 2013

      CHARACTER(LEN=99) ::                                                     &
         rainf, &
         snowf, &
         LWdown, &
         SWdown, &
         PSurf, &
         Qair, &
         Tair, &
         wind

   END TYPE globalMet_type
   
   TYPE(globalMet_type) :: globalMetfile
 
  TYPE gswp_type

     CHARACTER(LEN=200) ::                                                     &
          rainf, &
          snowf, &
          LWdown, &
          SWdown, &
          PSurf, &
          Qair, &
          Tair, &
          wind, &
          mask

  END TYPE gswp_type

  TYPE(gswp_type)      :: gswpfile


  INTEGER ::                                                                  &
       ncciy,      & ! year number (& switch) for gswp run
       ncid_rin,   & ! input netcdf restart file ID
       logn          ! log file unit number

  LOGICAL ::                                                                  &
       verbose,    & ! print init and param details of all grid cells?
       soilparmnew   ! read IGBP new soil map. Q.Zhang @ 12/20/2010

  ! ================ Veg and soil type variables ============================
  INTEGER, POINTER ::                                                         &
       soiltype_metfile(:,:),  & ! user defined soil type (from met file)
       vegtype_metfile(:,:)      ! user-def veg type (from met file)

   REAL, POINTER :: vegpatch_metfile(:,:) ! Anna: patchfrac for user-def vegtype


  TYPE parID_type ! model parameter IDs in netcdf file

     INTEGER :: bch,latitude,clay,css,rhosoil,hyds,rs20,sand,sfc,silt,        &
          ssat,sucs,swilt,froot,zse,canst1,dleaf,meth,za_tq,za_uv,             &
          ejmax,frac4,hc,lai,rp20,rpcoef,shelrb, vbeta, xalbnir,               &
          vcmax,xfang,ratecp,ratecs,area,refsbare,isoil,iveg,albsoil,          &
          taul,refl,tauw,refw,wai,vegcf,extkn,tminvj,tmaxvj,                   &
          veg_class,soil_class,mvtype,mstype,patchfrac,                        &
          WatSat,GWWatSat,SoilMatPotSat,GWSoilMatPotSat,                       &
          HkSat,GWHkSat,FrcSand,FrcClay,Clappb,Watr,GWWatr,sfc_vec,forg,swilt_vec, &
          slope,slope_std,GWdz,SatFracmax,Qhmax,QhmaxEfold,HKefold,HKdepth
     INTEGER :: ishorizon,nhorizons,clitt, &
          zeta,fsatmax, &
          gamma,ZR,F10

     INTEGER :: g0,g1 ! Ticket #56

  END TYPE parID_type

  ! =============== Logical  variables ============================
  TYPE input_details_type

     LOGICAL ::                                                               &
          Wind,    & ! T => 'Wind' is present; F => use vector component wind
          LWdown,  & ! T=> downward longwave is present in met file
          CO2air,  & ! T=> air CO2 concentration is present in met file
          PSurf,   & ! T=> surface air pressure is present in met file
          Snowf,   & ! T=> snowfall variable is present in met file
          avPrecip,& ! T=> ave rainfall present in met file (use for spinup)
          LAI,     & ! T=> LAI is present in the met file
          LAI_T,   & ! T=> LAI is time dependent, for each time step
          LAI_M,   & ! T=> LAI is time dependent, for each month
          LAI_P,   & ! T=> LAI is patch dependent
          parameters,&! TRUE if non-default parameters are found
          initial, & ! switched to TRUE when initialisation data are loaded
          patch,   & ! T=> met file have a subgrid veg/soil patch dimension
          laiPatch   ! T=> LAI file have a subgrid veg patch dimension

  END TYPE input_details_type

  TYPE(input_details_type) :: exists

  TYPE output_settings_type
    !! Output settings read from the `&cable` namelist. What is written, and how
    !! often, is set by the output configuration file named in `filename%output_config`.

    LOGICAL :: restart = .FALSE.
      !! Write a restart file at the end of the run.

    CHARACTER(LEN=7) :: grid = 'default'
      !! Layout of the output files: 'default' follows the meteorological forcing;
      !! 'land' or 'mask' force the compressed land-point or the lat/lon layout.

  END TYPE output_settings_type

  TYPE(output_settings_type), SAVE :: output

  ENUM, BIND(C)
    ENUMERATOR :: NO_CHECK = 0
    ENUMERATOR :: ON_TIMESTEP = 1
    ENUMERATOR :: ON_WRITE = 2
    ENUMERATOR :: RANGE_CHECK
  END ENUM
  TYPE checks_type
    LOGICAL :: energy_bal, mass_bal
    INTEGER(KIND(RANGE_CHECK)) :: ranges  ! 0 = NO , 1 = TIMESTEP , 2 = WRITE
    LOGICAL :: exit
  END TYPE checks_type

  TYPE(checks_type) :: check ! what types of checks to perform

  ! ============== Proxy input variables ================================
  REAL,POINTER,DIMENSION(:)  :: PrecipScale! precip scaling per site for spinup
  REAL,POINTER,DIMENSION(:,:)  :: defaultLAI ! in case met file/host model
  ! has no LAI
  REAL :: fixedCO2 ! CO2 level if CO2air not in met file

  ! For threading:
  !$OMP THREADPRIVATE(landpt,patch)
CONTAINS

  SUBROUTINE set_decomp_for_coupled()
    ! In coupled mode, we leave the decomposition up to the controlling UM-
    ! treat every process as it's own simulation
    land_decomp_start = 1
  END SUBROUTINE set_decomp_for_coupled

  FUNCTION to_land_index_global(land_index_local) RESULT(land_index_global)
    !! Translate local land index on current MPI rank to global land index
    INTEGER, INTENT(IN) :: land_index_local
    INTEGER :: land_index_global
    land_index_global = land_decomp_start + land_index_local - 1
  END FUNCTION to_land_index_global

END MODULE cable_IO_vars_module
