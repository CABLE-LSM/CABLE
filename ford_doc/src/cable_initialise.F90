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
! Purpose: default initialisation module for CABLE offline
!
! Contact: Bernard.Pak@csiro.au
!
! History: Developed by Gab Abramowitz
!          Since 1.4b fes split into fess and fesp
!          Significant changes: new routine 'extraRestart' for land-use change
!
!
! ==============================================================================
!
! MODULEs used: cable_abort_module
!               cable_def_types_mod
!               cable_IO_vars_module
!               cable_read_module
!               physical_constants
!               netcdf
!
!==============================================================================

MODULE cable_init_module

  USE cable_abort_module,       ONLY: abort, nc_abort
  USE cable_def_types_mod
  USE cable_IO_vars_module,       ONLY: latitude,longitude, patch,            &
       landpt,smoy,ncid_rin,max_vegpatches,          &
       soilparmnew,ncciy, vegtype_metfile,           &
       soiltype_metfile
  USE cable_read_module
  USE netcdf
  USE cable_common_module, ONLY : filename, cable_user

  IMPLICIT NONE

  PRIVATE
  PUBLIC get_default_inits, get_restart_data

  INTEGER :: ok ! netcdf status

CONTAINS
  !==============================================================================
  !
  ! Name: get_default_inits
  !
  ! Purpose: Loads initialisations based on Mk3L 50 year monthly
  !          climatology file
  !
  ! CALLed from: load_parameters
  !
  ! CALLs: abort
  !
  !==============================================================================


  !==============================================================================
  ! changes since version release on
  ! changes made by who on date
  !
  !==============================================================================

  SUBROUTINE get_default_inits(met,soil,ssnow,canopy,logn, EMSOIL)

    IMPLICIT NONE

    ! Input arguments
    TYPE (met_type), INTENT(IN)            :: met
    TYPE (soil_parameter_type), INTENT(IN) :: soil
    TYPE (soil_snow_type), INTENT(INOUT)   :: ssnow
    TYPE (canopy_type), INTENT(OUT)        :: canopy
    INTEGER,INTENT(IN)                :: logn     ! log file unit number
    REAL, INTENT(IN) :: EMSOIL
    ! Local variables
    INTEGER :: e,i,j  ! do loop counter

    WRITE(logn,*) ' Initializing variables.'

    DO e=1,mp ! over all patches

       ! The following write statements are redundant in online runs
       !       ! Write to log file:
       !       WRITE(logn,'(A21,I8,2(A15,1X,F9.4,1X))') '     Land grid point:',e, &
       !            '      Latitude ',latitude(e),'Longitude',longitude(e)
       !       WRITE(logn,'(A46,2(1X,F8.3,1X,A3))') &
       !            '        is closest to default gridcell centred', &
       !            REAL(lat_inits(final_y)),'lat', REAL(lon_inits(final_x)),'lon'
       ! Only the following snow inits are necessary,
       ! soilsnow will update other variables.
       IF(ssnow%snowd(e)>0.0) THEN  ! in cm
          ssnow%ssdnn(e)  = 120.0   ! overall snow density (kg/m3)
          ssnow%ssdn(e,:)   = 120.0 ! snow density per layer (kg/m3)
          ssnow%snage(e)  = 0.0     ! snow age (fresh)
          ssnow%isflag(e) = 0
       ELSE
          ssnow%ssdnn(e)  = 140.0   ! overall snow density (kg/m3)
          ssnow%osnowd(e) = 0.0     ! snow depth prev timestep (mm or kg/m2)
          ssnow%snage(e)  = 0.0     ! snow age
          ssnow%isflag(e) = 0       ! snow layer scheme flag
          ! (0 = no/little snow, 1=snow)
          ssnow%tggsn(e,:)  = 273.1 ! snow temperature per layer (K)
          ssnow%ssdn(e,:)   = 140.0 ! snow density per layer (kg/m3)
          ssnow%smass(e,:)  = 0.0   ! snow mass per layer (kg/m^2)
       END IF
       ! Soil ice:
       WHERE(ssnow%tgg(e,:)<273.15)
          ssnow%wbice(e,:)  = ssnow%wb(e,:)*0.8
       ELSEWHERE
          ssnow%wbice(e,:) = 0.0
       END WHERE

    END DO

    IF(ANY(ssnow%tgg>350.0).OR.ANY(ssnow%tgg<180.0)) CALL abort('Soil temps nuts')
    IF(ANY(ssnow%albsoilsn>1.0).OR.ANY(ssnow%albsoilsn<0.0)) CALL abort('Albedo nuts')

    ! Site independent initialisations (all gridcells):
    ! soil+snow albedo for infrared (other values read in below):
    ssnow%albsoilsn(:,3) = 1.0 - emsoil
    !   ssnow%albsoilsn(:,3) = 0.05  ! YP Nov2009 (fix cold bias)
    canopy%cansto  = 0.0   ! canopy water storage (mm or kg/m2)
    canopy%sghflux = 0.0
    canopy%ghflux  = 0.0
    ssnow%runoff   = 0.0   ! runoff total = subsurface + surface runoff
    ssnow%rnof1    = 0.0   ! surface runoff (mm/timestepsize)
    ssnow%rnof2    = 0.0   ! deep drainage (mm/timestepsize)
    ssnow%rtsoil   = 100.0 ! turbulent resistance for soil
    canopy%ga      = 0.0   ! ground heat flux (W/m2)
    canopy%dgdtg   = 0.0   ! derivative of ground heat flux wrt soil temp
    canopy%fev     = 0.0   ! latent heat flux from vegetation (W/m2)
    canopy%fes     = 0.0   ! latent heat flux from soil (W/m2)
    canopy%fhs     = 0.0   ! sensible heat flux from soil (W/m2)
    canopy%us = 0.1 ! friction velocity (needed in roughness before first call to canopy: should in be in restart?)

  END SUBROUTINE get_default_inits

  !==============================================================================
  !
  ! Name: get_restart_data
  !
  ! Purpose: Reads initialisations and parameters from restart file
  !
  ! CALLed from: load_parameters
  !
  ! CALLs: nc_abort
  !        extraRestart
  !        readpar
  !        abort
  !
  ! Input file: [restart].nc
  !
  !==============================================================================

  SUBROUTINE get_restart_data(logn,ssnow,canopy,rough,bgc,                       &
       bal,veg,soil,rad,vegparmnew, EMSOIL)

    IMPLICIT NONE

    ! Input arguments
    INTEGER, INTENT(IN)                       :: logn       ! log file number
    TYPE (soil_snow_type),INTENT(INOUT)       :: ssnow      ! soil and snow variables
    TYPE (bgc_pool_type),INTENT(INOUT)        :: bgc        ! carbon pool variables
    TYPE (canopy_type),INTENT(INOUT)          :: canopy     ! vegetation variables
    TYPE (roughness_type),INTENT(INOUT)       :: rough      ! roughness varibles
    TYPE (balances_type),INTENT(INOUT)        :: bal        ! energy + water balance variables
    TYPE (veg_parameter_type), INTENT(INOUT)  :: veg        ! vegetation parameters
    TYPE (soil_parameter_type), INTENT(INOUT) :: soil       ! soil parameters
    TYPE (radiation_type),INTENT(INOUT)       :: rad
    LOGICAL,INTENT(IN)                        :: vegparmnew ! are we using the new format?

    REAL, INTENT(IN) :: EMSOIL
    ! Local variables
    REAL, POINTER,DIMENSION(:)           ::                                &
         lat_restart,                                                           &
         lon_restart
    INTEGER,POINTER,DIMENSION(:)         :: INvar
    !    REAL, POINTER,DIMENSION(:,:) :: surffrac ! fraction of each surf type
    INTEGER ::                                                             &
         mland_restart,         & ! number of land points in restart file
         INvegt,                &
         INsoilt,               &
         INpatch,               &
         mpatchID,              &
                                !     surftype_restart,      & ! number of surface types in restart file
         latID, lonID,          & ! lat,lon variable ID
         mvtypeID,              & ! veg type variable ID
         mstypeID,              & ! soil type variable ID
         mlandID,               & ! netcdf ID for land points
         i,                     & ! do loop counter
                                !     jj,                    & ! do loop counter
         parID                    ! parameter's netcdf ID
    LOGICAL ::                                                                  &
         from_restart = .TRUE., & ! insist variables/params load
         dummy                    ! To replace completeSet in parameter read; unused
    REAL,    ALLOCATABLE :: var_r(:)
    REAL,    ALLOCATABLE :: var_r2(:,:)

    ! Write to screen the restart file is found:
    WRITE(*,*) 'Reading restart data from: ' ,TRIM(filename%restart_in)

    ! Check number of gridpoints in restart file is correct:
    ok = NF90_INQ_DIMID(ncid_rin,'mland',mlandID)
    IF(ok /= NF90_NOERR) THEN
       ok = NF90_INQ_DIMID(ncid_rin,'mp',mlandID) ! name used before sep2010
       IF(ok /= NF90_NOERR) CALL nc_abort                                       &
            (ok,'Error finding mland dimension in restart file '                &
            //TRIM(filename%restart_in)//' (SUBROUTINE get_restart)')
    END IF
    ok = NF90_INQUIRE_DIMENSION(ncid_rin,mlandID,len=mland_restart)
    PRINT *, 'number of land point in restart file: ', mland_restart
    IF(ok /= NF90_NOERR) CALL nc_abort                                          &
         (ok,'Error finding number of land points in restart file '             &
         //TRIM(filename%restart_in)//' (SUBROUTINE get_restart)')
    IF(mland_restart /= mland) CALL abort('Number of land points in '//         &
         'restart file '//TRIM(filename%restart_in)//                           &
         ' differs from number in met file '//TRIM(filename%met))

    ! Added the checking of mp; if not equal, redirect to another
    ! subroutine to get grid-based info (BP may2010) for LULUC
    ok = NF90_INQ_DIMID(ncid_rin,'mp_patch',mpatchID)
    IF(ok /= NF90_NOERR) THEN
       ok = NF90_INQ_DIMID(ncid_rin,'mp',mpatchID) ! old file before sep2010
       IF(ok /= NF90_NOERR)  CALL nc_abort                                      &
            (ok,'Error finding mp_patch dimension in restart file '             &
            //TRIM(filename%restart_in)//' (SUBROUTINE get_restart)')
    ENDIF
    ok = NF90_INQUIRE_DIMENSION(ncid_rin,mpatchID,len=INpatch)
    PRINT *, 'total number of patches in restart file: ', INpatch
    IF (INpatch /= mp) THEN
       CALL extraRestart(INpatch,ssnow,canopy,rough,bgc,                        &
            bal,veg,soil,rad, EMSOIL)
       RETURN
    ENDIF

    ! removed the following because already in IGBP types (BP apr08)
    !    ! Check number of surface types is correct:
    !    ok = NF90_INQ_DIMID(ncid_rin,'surftype',surftypeID)
    !    IF(ok /= NF90_NOERR) CALL nc_abort &
    !         (ok,'Error finding surftype dimension in restart file ' &
    !         //TRIM(filename%restart_in)//' (SUBROUTINE get_restart)')
    !    ok = NF90_INQUIRE_DIMENSION(ncid_rin,surftypeID,len=surftype_restart)
    !    IF(ok /= NF90_NOERR) CALL nc_abort &
    !         (ok,'Error finding number of surface types in restart file ' &
    !         //TRIM(filename%restart_in)//' (SUBROUTINE get_restart)')
    !    IF(surftype_restart /= 4) CALL &
    !         abort('Number of surface types per grid cell in '// &
    !         'restart file '//TRIM(filename%restart_in)// &
    !         ' differs from number in cable_variables.f90 ')
    !    ! Get surffrac variable:
    !    ALLOCATE(surffrac(mland,4))
    !    ok = NF90_INQ_VARID(ncid_rin,'surffrac',surffracID)
    !    IF(ok /= NF90_NOERR) CALL nc_abort &
    !         (ok,'Error finding surffrac in restart file ' &
    !         //TRIM(filename%restart_in)//' (SUBROUTINE get_restart)')
    !    ok=NF90_GET_VAR(ncid_rin,surffracID,surffrac)
    !    IF(ok/=NF90_NOERR) CALL nc_abort(ok,'Error reading surffrac in file ' &
    !         //TRIM(filename%restart_in)// '(SUBROUTINE get_restart)')
    !    landpt(:)%veg%frac =  surffrac(:,1)
    !    landpt(:)%urban%frac = surffrac(:,2)
    !    landpt(:)%lake%frac = surffrac(:,3)
    !    landpt(:)%ice%frac = surffrac(:,4)
    !    DEALLOCATE(surffrac)


    ! check that lat/lon b/w run and restart are compatible:
    ALLOCATE(lat_restart(mland),lon_restart(mland))
    ok = NF90_INQ_VARID(ncid_rin,'latitude',latID)
    IF(ok /= NF90_NOERR) CALL nc_abort                                          &
         (ok,'Error finding latitude in restart file '                          &
         //TRIM(filename%restart_in)//' (SUBROUTINE get_restart)')
    ok = NF90_INQ_VARID(ncid_rin,'longitude',lonID)
    IF(ok /= NF90_NOERR) CALL nc_abort                                          &
         (ok,'Error finding longitude in restart file '                         &
         //TRIM(filename%restart_in)//' (SUBROUTINE get_restart)')
    ok=NF90_GET_VAR(ncid_rin,latID,lat_restart)
    IF(ok/=NF90_NOERR) CALL nc_abort(ok,'Error reading latitude in file '       &
         //TRIM(filename%restart_in)// '(SUBROUTINE get_restart)')
    ! Removed rad%latitude from here as it is already done in write_default_params
    ! (BP may2010)
    !    ! Set rad%latitude parameter
    !    DO i=1,mland
    !       ! All patches in a single grid cell have the same latitude:
    !       rad%latitude(landpt(i)%cstart:landpt(i)%cend)=lat_restart(i)
    !    END DO
    ok=NF90_GET_VAR(ncid_rin,lonID,lon_restart)
    IF(ok/=NF90_NOERR) CALL nc_abort(ok,'Error reading longitude in file '      &
         //TRIM(filename%restart_in)// '(SUBROUTINE get_restart)')
    IF(ANY(ABS(lat_restart-latitude)>0.01))                                     &
         CALL abort('Latitude of land points in '//                             &
         'restart file '//TRIM(filename%restart_in)//                           &
         ' differs from met file '//TRIM(filename%met))
    IF(ANY(ABS(lon_restart-longitude)>0.01))                                    &
         CALL abort('Longitude of land points in '//                            &
         'restart file '//TRIM(filename%restart_in)//                           &
         ' differs from met file '//TRIM(filename%met))
    DEALLOCATE(lat_restart,lon_restart)

    ! Check that the number of vegetation types is present in restart file:
    ok = NF90_INQ_VARID(ncid_rin,'mvtype',mvtypeID)
    IF(ok /= NF90_NOERR) THEN
       ok = NF90_INQ_VARID(ncid_rin,'nvegt',mvtypeID)
       IF(ok == NF90_NOERR) THEN
          ok=NF90_GET_VAR(ncid_rin,mvtypeID,INvegt)
          IF(ok/=NF90_NOERR) CALL nc_abort(ok,'Error reading nvegt in file '    &
               //TRIM(filename%restart_in)// '(SUBROUTINE get_restart)')
          IF(INvegt > 17) CALL nc_abort(ok,'Error: nvegt value in file '        &
               //TRIM(filename%restart_in)// ' out of range')
          IF (INvegt /= mvtype) PRINT *, 'Warning: INvegt, nvegt = ', INvegt, mvtype
       ENDIF
       ! Removed the following as mvtype is determined earlier from
       ! reading in def_veg_params_xx.txt (BP may2010)
       !       IF(vegparmnew) THEN
       !          mvtype = 17
       !       ELSE
       !          mvtype = 13
       !       ENDIF
    ELSE
       ! Changed the read-in variable name so that mvtype would not be overwritten
       ! and added some more checking (BP may2010)
       ok=NF90_GET_VAR(ncid_rin,mvtypeID,INvegt)
       IF(ok/=NF90_NOERR) CALL nc_abort(ok,'Error reading mvtype in file '      &
            //TRIM(filename%restart_in)// '(SUBROUTINE get_restart)')
       IF(INvegt > 17) CALL nc_abort(ok,'Error: mvtype value in file '          &
            //TRIM(filename%restart_in)// ' out of range')
       IF (INvegt /= mvtype) PRINT *, 'Warning: INvegt, mvtype = ', INvegt, mvtype
    ENDIF
    ! Check that the number of soil types is present in restart file:
    ok = NF90_INQ_VARID(ncid_rin,'mstype',mstypeID)
    IF(ok /= NF90_NOERR) THEN
       ok = NF90_INQ_VARID(ncid_rin,'nsoilt',mstypeID)
       IF(ok == NF90_NOERR) THEN
          ok=NF90_GET_VAR(ncid_rin,mstypeID,INsoilt)
          IF(ok/=NF90_NOERR) CALL nc_abort(ok,'Error reading nsoilt in file '   &
               //TRIM(filename%restart_in)// '(SUBROUTINE get_restart)')
          IF(INsoilt /= mstype) CALL nc_abort(ok,'Error: nsoilt value in file ' &
               //TRIM(filename%restart_in)// ' is wrong')
       ENDIF
       ! Removed the following as mstype is determined earlier from
       ! reading in def_soil_params.txt (BP may2010)
       !       mstype = 9
    ELSE
       ! Changed the read-in variable name so that mstype would not be overwritten
       ! (BP may2010)
       ok=NF90_GET_VAR(ncid_rin,mstypeID,INsoilt)
       IF(ok/=NF90_NOERR) CALL nc_abort(ok,'Error reading mstype in file '      &
            //TRIM(filename%restart_in)// '(SUBROUTINE get_restart)')
       IF(INsoilt /= mstype) CALL nc_abort(ok,'Error: mstype value in file '    &
            //TRIM(filename%restart_in)// ' is wrong')
    ENDIF

    dummy=.TRUE. ! initialise for completeness only - not used

    ! Get variable initialisations =============================
    ! Arguments are: netcdf file ID; parameter name;
    !   complete set check; parameter value; filename for error messages;
    !   number of veg/soil patches in met file; switch to indicate
    !   size of dimensions of the parameter; an indicator to show
    !   we're reading from the restart file.
    ! Use 'defd' for single dim double precision.
    ! Use, e.g., 'msd' to fetch double precision 2D soil varible
    CALL readpar(ncid_rin,'tgg',dummy,ssnow%tgg,filename%restart_in,            &
         max_vegpatches,'ms',from_restart,mp)
    CALL readpar(ncid_rin,'wb',dummy,ssnow%wb,filename%restart_in,              &
         max_vegpatches,'msd',from_restart,mp)
    CALL readpar(ncid_rin,'wbice',dummy,ssnow%wbice,filename%restart_in,        &
         max_vegpatches,'msd',from_restart,mp)
    !    WHERE (ssnow%tgg > 273.2 .AND. ssnow%wbice >0.0) ssnow%wbice=0.0
    CALL readpar(ncid_rin,'gammzz',dummy,ssnow%gammzz,filename%restart_in,      &
         max_vegpatches,'msd',from_restart,mp)
    CALL readpar(ncid_rin,'tss',dummy,ssnow%tss,filename%restart_in,            &
         max_vegpatches,'def',from_restart,mp)
    CALL readpar(ncid_rin,'ssdnn',dummy,ssnow%ssdnn,filename%restart_in,        &
         max_vegpatches,'def',from_restart,mp)
    CALL readpar(ncid_rin,'ssdn',dummy,ssnow%ssdn,filename%restart_in,          &
         max_vegpatches,'snow',from_restart,mp)
    CALL readpar(ncid_rin,'osnowd',dummy,ssnow%osnowd,filename%restart_in,      &
         max_vegpatches,'def',from_restart,mp)
    CALL readpar(ncid_rin,'smass',dummy,ssnow%smass,filename%restart_in,        &
         max_vegpatches,'snow',from_restart,mp)
    CALL readpar(ncid_rin,'sdepth',dummy,ssnow%sdepth,filename%restart_in,      &
         max_vegpatches,'snow',from_restart,mp)
    CALL readpar(ncid_rin,'tggsn',dummy,ssnow%tggsn,filename%restart_in,        &
         max_vegpatches,'snow',from_restart,mp)
    CALL readpar(ncid_rin,'snage',dummy,ssnow%snage,filename%restart_in,        &
         max_vegpatches,'def',from_restart,mp)
    CALL readpar(ncid_rin,'snowd',dummy,ssnow%snowd,filename%restart_in,        &
         max_vegpatches,'def',from_restart,mp)
    CALL readpar(ncid_rin,'rtsoil',dummy,ssnow%rtsoil,filename%restart_in,      &
         max_vegpatches,'def',from_restart,mp)
    CALL readpar(ncid_rin,'isflag',dummy,ssnow%isflag,filename%restart_in,      &
         max_vegpatches,'def',from_restart,mp)
    CALL readpar(ncid_rin,'albsoilsn',dummy,ssnow%albsoilsn,                    &
         filename%restart_in,max_vegpatches,'nrb',from_restart,mp)
    ssnow%albsoilsn(:,3) = 1.0 - emsoil  !! (BP Nov 2009)
    CALL readpar(ncid_rin,'rnof1',dummy,ssnow%rnof1,filename%restart_in,        &
         max_vegpatches,'def',from_restart,mp)
    CALL readpar(ncid_rin,'rnof2',dummy,ssnow%rnof2,filename%restart_in,        &
         max_vegpatches,'def',from_restart,mp)
    CALL readpar(ncid_rin,'runoff',dummy,ssnow%runoff,filename%restart_in,      &
         max_vegpatches,'def',from_restart,mp)

    !MD
    ok = NF90_INQ_VARID(ncid_rin,'GWwb',parID)
    IF(ok == NF90_NOERR) THEN
       CALL readpar(ncid_rin,'GWwb',dummy,ssnow%GWwb,filename%restart_in,            &
            max_vegpatches,'def',from_restart,mp)
    ELSE
       ssnow%GWwb = 0.95*soil%ssat
    END IF

!!$   IF(cable_user%SOIL_STRUC=='sli'.or.cable_user%FWSOIL_SWITCH=='Haverd2013') THEN
!!$      CALL readpar(ncid_rin,'gamma',dummy,veg%gamma,filename%restart_in,           &
!!$           max_vegpatches,'def',from_restart,mp)
!!$   ENDIF

    IF(cable_user%SOIL_STRUC=='sli') THEN
       CALL readpar(ncid_rin,'S',dummy,ssnow%S,filename%restart_in, &
            max_vegpatches,'ms',from_restart,mp)
       CALL readpar(ncid_rin,'Tsoil',dummy,ssnow%Tsoil,filename%restart_in, &
            max_vegpatches,'ms',from_restart,mp)
       CALL readpar(ncid_rin,'h0',dummy,ssnow%h0,filename%restart_in, &
            max_vegpatches,'def',from_restart,mp)
       CALL readpar(ncid_rin,'nsnow',dummy,ssnow%nsnow,filename%restart_in, &
            max_vegpatches,'def',from_restart,mp)
       CALL readpar(ncid_rin,'Tsurface',dummy,ssnow%Tsurface,filename%restart_in, &
            max_vegpatches,'def',from_restart,mp)
       CALL readpar(ncid_rin,'snowliq',dummy,ssnow%snowliq,filename%restart_in, &
            max_vegpatches,'snow',from_restart,mp)
       CALL readpar(ncid_rin,'sconds',dummy,ssnow%sconds,filename%restart_in, &
            max_vegpatches,'snow',from_restart,mp)
!!$       CALL readpar(ncid_rin,'ZR',dummy,veg%ZR, &
!!$            filename%restart_in,max_vegpatches,'def',from_restart,mp)
!!$       CALL readpar(ncid_rin,'F10',dummy,veg%F10, &
!!$            filename%restart_in,max_vegpatches,'def',from_restart,mp)
!!$       CALL readpar(ncid_rin,'zeta',dummy,soil%zeta,filename%restart_in,           &
!!$            max_vegpatches,'def',from_restart,mp)
!!$       CALL readpar(ncid_rin,'fsatmax',dummy,soil%fsatmax,filename%restart_in,           &
!!$            max_vegpatches,'def',from_restart,mp)
!!$       CALL readpar(ncid_rin,'nhorizons',dummy,soil%nhorizons,filename%restart_in,           &
!!$            max_vegpatches,'def',from_restart,mp)
!!$       ALLOCATE(var_r2(mp,ms))
!!$       CALL readpar(ncid_rin,'ishorizon',dummy,var_r2,filename%restart_in,           &
!!$            max_vegpatches,'ms',from_restart,mp)
!!$       soil%ishorizon = int(var_r2)
!!$       DEALLOCATE(var_r2)
!!$       CALL readpar(ncid_rin,'clitt',dummy,veg%clitt,filename%restart_in,           &
!!$            max_vegpatches,'def',from_restart,mp)
    ENDIF
    CALL readpar(ncid_rin,'cansto',dummy,canopy%cansto,filename%restart_in,     &
         max_vegpatches,'def',from_restart,mp)
    CALL readpar(ncid_rin,'sghflux',dummy,canopy%sghflux,filename%restart_in,   &
         max_vegpatches,'def',from_restart,mp)
    CALL readpar(ncid_rin,'ghflux',dummy,canopy%ghflux,filename%restart_in,     &
         max_vegpatches,'def',from_restart,mp)
    CALL readpar(ncid_rin,'ga',dummy,canopy%ga,filename%restart_in,             &
         max_vegpatches,'def',from_restart,mp)
    CALL readpar(ncid_rin,'dgdtg',dummy,canopy%dgdtg,filename%restart_in,       &
         max_vegpatches,'def',from_restart,mp)
    CALL readpar(ncid_rin,'fev',dummy,canopy%fev,filename%restart_in,           &
         max_vegpatches,'def',from_restart,mp)
    !jhan:hack - elimiinate call as r_2 now
    !    CALL readpar(ncid_rin,'fes',dummy,canopy%fes,filename%restart_in,          &
    !         max_vegpatches,'def',from_restart,mp)
    CALL readpar(ncid_rin,'fhs',dummy,canopy%fhs,filename%restart_in,           &
         max_vegpatches,'def',from_restart,mp)
    CALL readpar(ncid_rin,'cplant',dummy,bgc%cplant,filename%restart_in,        &
         max_vegpatches,'ncp',from_restart,mp)
    CALL readpar(ncid_rin,'csoil',dummy,bgc%csoil,filename%restart_in,          &
         max_vegpatches,'ncs',from_restart,mp)
    CALL readpar(ncid_rin,'wbtot0',dummy,bal%wbtot0,filename%restart_in,        &
         max_vegpatches,'def',from_restart,mp)
    CALL readpar(ncid_rin,'osnowd0',dummy,bal%osnowd0,filename%restart_in,      &
         max_vegpatches,'def',from_restart,mp)
    ! The following two restart file additions are to initialise Mk3L:
    CALL readpar(ncid_rin,'albedo',dummy,rad%albedo,filename%restart_in,        &
         max_vegpatches,'nrb',from_restart,mp)
    CALL readpar(ncid_rin,'trad',dummy,rad%trad,filename%restart_in,            &
         max_vegpatches,'def',from_restart,mp)

    ! Get model parameters =============================================
    ! rad%latitude set above in lat/lon checking section
    ALLOCATE(INvar(mp))
    CALL readpar(ncid_rin,'iveg',dummy,INvar,filename%restart_in,               &
         max_vegpatches,'def',from_restart,mp)
    IF (ASSOCIATED(vegtype_metfile)) THEN
       ! met file iveg info is now in veg%iveg
       IF (ANY(INvar /= veg%iveg)) THEN
          PRINT *, 'Error: veg type in restart file different from met input'
          PRINT *, 'Recommend not using this restart file as parameters have changed.'
          CALL abort('Check iveg in '//filename%restart_in)
       ENDIF
    ELSE
       ! no problem with overwriting default values
       veg%iveg = INvar
    ENDIF
    !    CALL readpar(ncid_rin,'iveg',dummy,veg%iveg,filename%restart_in,           &
    !         max_vegpatches,'def',from_restart,mp)
    IF (.NOT.CABLE_USER%POPLUC) THEN
       CALL readpar(ncid_rin,'patchfrac',dummy,patch(:)%frac,filename%restart_in,  &
            max_vegpatches,'def',from_restart,mp)
    ENDIF
    !    DO i=1, mland
    !    DO jj = landpt(i)%cstart, landpt(i)%cend
    !      IF (INvar(jj) /= veg%iveg(jj)) THEN
    !        PRINT *, 'veg type in restart file is weird.'
    !        PRINT *, 'mland and mp #: ', i, jj
    !        PRINT *, 'INvar, veg%iveg: ', INvar(jj), veg%iveg(jj)
    !        PRINT *, 'lon and lat: ', longitude(i), latitude(i)
    !      END IF
    !    END DO
    !    END DO
    !! getting rid of spurious veg types in Antarctica from the CCAM2Mk3L process
    !! Doing it once will fix the problem in the restart file in subsequent runs
    !    DO i=1, mland
    !      IF ( rad%latitude(landpt(i)%cstart) < -60.0 .AND. &
    !           patch(landpt(i)%cstart)%frac < 1.0 ) THEN
    !        IF ( veg%iveg(landpt(i)%cstart) <= 15 ) THEN
    !          patch(landpt(i)%cstart:landpt(i)%cend)%frac = 0.0
    !          patch(landpt(i)%cstart)%frac = 1.0
    !          veg%iveg(landpt(i)%cstart:landpt(i)%cend) = 15
    !        END IF
    !      END IF
    !    END DO
    !! end of fix to spurious veg types

    CALL readpar(ncid_rin,'isoil',dummy,INvar,filename%restart_in,              &
         max_vegpatches,'def',from_restart,mp)
    IF (ASSOCIATED(soiltype_metfile)) THEN
       ! met file isoil info is now in soil%isoilm
       IF (ANY(INvar /= soil%isoilm)) THEN
          PRINT *, 'Error: soil type in restart file different from met input'
          PRINT *, 'Recommend not using this restart file as parameters have changed.'
          CALL abort('Check isoil in '//filename%restart_in)
       ENDIF
    ELSE
       ! no problem with overwriting default values
       soil%isoilm = INvar
    ENDIF
    !    CALL readpar(ncid_rin,'isoil',dummy,soil%isoilm,filename%restart_in,       &
    !         max_vegpatches,'def',from_restart,mp)
    !   CALL readpar(ncid_rin,'clay',dummy,soil%clay,filename%restart_in,           &
    !                max_vegpatches,'def',from_restart,mp)
    !   CALL readpar(ncid_rin,'sand',dummy,soil%sand,filename%restart_in,           &
    !                max_vegpatches,'def',from_restart,mp)
    !   CALL readpar(ncid_rin,'silt',dummy,soil%silt,filename%restart_in,           &
    !                max_vegpatches,'def',from_restart,mp)
    !   IF ( .NOT. soilparmnew) THEN  ! Q.Zhang @12/20/2010
    !      CALL readpar(ncid_rin,'ssat',dummy,soil%ssat,filename%restart_in,        &
    !                   max_vegpatches,'def',from_restart,mp)
    !      CALL readpar(ncid_rin,'sfc',dummy,soil%sfc,filename%restart_in,          &
    !                   max_vegpatches,'def',from_restart,mp)
    !      CALL readpar(ncid_rin,'swilt',dummy,soil%swilt,filename%restart_in,      &
    !                   max_vegpatches,'def',from_restart,mp)
    !      CALL readpar(ncid_rin,'bch',dummy,soil%bch,filename%restart_in,          &
    !                   max_vegpatches,'def',from_restart,mp)
    !      CALL readpar(ncid_rin,'hyds',dummy,soil%hyds,filename%restart_in,        &
    !                  max_vegpatches,'def',from_restart,mp)
    !      CALL readpar(ncid_rin,'sucs',dummy,soil%sucs,filename%restart_in,        &
    !                   max_vegpatches,'def',from_restart,mp)
    !      CALL readpar(ncid_rin,'css',dummy,soil%css,filename%restart_in,          &
    !                   max_vegpatches,'def',from_restart,mp)
    !      CALL readpar(ncid_rin,'rhosoil',dummy,soil%rhosoil,filename%restart_in,  &
    !                   max_vegpatches,'def',from_restart,mp)
    !      IF (ncciy > 0 .AND. filename%restart_in(22:23) == 'HQ') THEN
    !         ok = NF90_INQ_VARID(ncid_rin,'albsoil',parID)
    !         IF(ok == NF90_NOERR) THEN
    !            ALLOCATE(var_r(mp))
    !            ok= NF90_GET_VAR(ncid_rin,parID,var_r,start=(/1/),count=(/mp/))
    !            IF(ok /= NF90_NOERR) CALL nc_abort                                 &
    !                 (ok,'Error reading '//'albsoil'//' in file '                  &
    !                 //TRIM(filename%restart_in)//' (SUBROUTINE get_restart_data)')
    !            soil%albsoil(:,1) = var_r / 3.0
    !            soil%albsoil(:,2) = var_r * 2.0 / 3.0
    !            soil%albsoil(:,3) = 0.005
    !            DEALLOCATE(var_r)
    !         END IF
    !      ELSE
    !         CALL readpar(ncid_rin,'albsoil',dummy,soil%albsoil,filename%restart_in, &
    !                      max_vegpatches,'nrb',from_restart,mp)
    !      END IF
    !   END IF
    !    CALL readpar(ncid_rin,'rs20',dummy,soil%rs20,filename%restart_in,          &
    !         max_vegpatches,'def',from_restart,mp)
    !   CALL readpar(ncid_rin,'rs20',dummy,veg%rs20,filename%restart_in,            &
    !                max_vegpatches,'def',from_restart,mp)
    !   CALL readpar(ncid_rin,'froot',dummy,veg%froot,filename%restart_in,          &
    !                max_vegpatches,'ms',from_restart,mp)
    !   CALL readpar(ncid_rin,'hc',dummy,veg%hc,filename%restart_in,                &
    !                max_vegpatches,'def',from_restart,mp)
    !   CALL readpar(ncid_rin,'canst1',dummy,veg%canst1,filename%restart_in,        &
    !                max_vegpatches,'def',from_restart,mp)
    !   CALL readpar(ncid_rin,'dleaf',dummy,veg%dleaf,filename%restart_in,          &
    !                max_vegpatches,'def',from_restart,mp)
    !   CALL readpar(ncid_rin,'frac4',dummy,veg%frac4,filename%restart_in,          &
    !                max_vegpatches,'def',from_restart,mp)
    !   CALL readpar(ncid_rin,'ejmax',dummy,veg%ejmax,filename%restart_in,          &
    !                max_vegpatches,'def',from_restart,mp)
    !   CALL readpar(ncid_rin,'vcmax',dummy,veg%vcmax,filename%restart_in,          &
    !                max_vegpatches,'def',from_restart,mp)
    !   CALL readpar(ncid_rin,'rp20',dummy,veg%rp20,filename%restart_in,            &
    !                max_vegpatches,'def',from_restart,mp)
    !   CALL readpar(ncid_rin,'rpcoef',dummy,veg%rpcoef,filename%restart_in,        &
    !                max_vegpatches,'def',from_restart,mp)
    !   CALL readpar(ncid_rin,'shelrb',dummy,veg%shelrb,filename%restart_in,        &
    !                max_vegpatches,'def',from_restart,mp)
    !   CALL readpar(ncid_rin,'xfang',dummy,veg%xfang,filename%restart_in,          &
    !                max_vegpatches,'def',from_restart,mp)
    !   CALL readpar(ncid_rin,'wai',dummy,veg%wai,filename%restart_in,              &
    !                max_vegpatches,'def',from_restart,mp)
    !   CALL readpar(ncid_rin,'vegcf',dummy,veg%vegcf,filename%restart_in,          &
    !                max_vegpatches,'def',from_restart,mp)
    !   CALL readpar(ncid_rin,'extkn',dummy,veg%extkn,filename%restart_in,          &
    !                max_vegpatches,'def',from_restart,mp)
    !   CALL readpar(ncid_rin,'tminvj',dummy,veg%tminvj,filename%restart_in,        &
    !                max_vegpatches,'def',from_restart,mp)
    !   CALL readpar(ncid_rin,'tmaxvj',dummy,veg%tmaxvj,filename%restart_in,        &
    !                max_vegpatches,'def',from_restart,mp)
    !   CALL readpar(ncid_rin,'vbeta',dummy,veg%vbeta,filename%restart_in,          &
    !                max_vegpatches,'def',from_restart,mp)
    !   CALL readpar(ncid_rin,'xalbnir',dummy,veg%xalbnir,filename%restart_in,      &
    !                max_vegpatches,'def',from_restart,mp)
    !   veg%xalbnir = 1.0   ! xalbnir will soon be removed totally
    !   CALL readpar(ncid_rin,'g0',dummy,veg%g0,filename%restart_in,            &
    !                max_vegpatches,'def',from_restart,mp) ! Ticket #56
    !   CALL readpar(ncid_rin,'g1',dummy,veg%g1,filename%restart_in,            &
    !                max_vegpatches,'def',from_restart,mp) ! Ticket #56
    !   CALL readpar(ncid_rin,'meth',dummy,veg%meth,filename%restart_in,            &
    !                max_vegpatches,'def',from_restart,mp)
    !   ! special treatment of za with the introduction of za_uv and za_tq
    ! in case an old restart file is used
    !   ok = NF90_INQ_VARID(ncid_rin,'za',parID)
    !   IF(ok == NF90_NOERR) THEN ! if it does exist
    !      CALL readpar(ncid_rin,'za',dummy,rough%za_uv,filename%restart_in,        &
    !                   max_vegpatches,'def',from_restart,mp)
    !      CALL readpar(ncid_rin,'za',dummy,rough%za_tq,filename%restart_in,        &
    !                   max_vegpatches,'def',from_restart,mp)
    !   ELSE
    !      CALL readpar(ncid_rin,'za_uv',dummy,rough%za_uv,filename%restart_in,     &
    !                   max_vegpatches,'def',from_restart,mp)
    !      CALL readpar(ncid_rin,'za_tq',dummy,rough%za_tq,filename%restart_in,     &
    !                   max_vegpatches,'def',from_restart,mp)
    !   ENDIF
    CALL readpar(ncid_rin,'zse',dummy,soil%zse,filename%restart_in,             &
         max_vegpatches,'ms',from_restart,mp)
    !   CALL readpar(ncid_rin,'ratecp',dummy,bgc%ratecp,filename%restart_in,        &
    !                max_vegpatches,'ncp',from_restart,mp)
    !   CALL readpar(ncid_rin,'ratecs',dummy,bgc%ratecs,filename%restart_in,        &
    !                max_vegpatches,'ncs',from_restart,mp)
    !
    ! Close restart file:
    ok = NF90_CLOSE(ncid_rin)
    IF(ok/=NF90_NOERR) CALL nc_abort(ok,'Error closing restart file '           &
         //TRIM(filename%restart_in)// '(SUBROUTINE get_restart)')

  END SUBROUTINE get_restart_data

  !==============================================================================
  !
  ! Name: extraRestart
  !
  ! Purpose: Redistribute the patches if restart file does not match
  !
  ! CALLed from: get_restart_data
  !
  ! CALLs: nc_abort
  !        readpar
  !        redistr_r2d
  !        redistr_rd
  !        redistr_r
  !        redistr_r2
  !        redistr_i
  !
  ! Input file: [restart].nc
  !
  !==============================================================================

  SUBROUTINE extraRestart(INpatch,ssnow,canopy,rough,bgc,                        &
       bal,veg,soil,rad, EMSOIL)
    ! Assume this subroutine to be used for first simulation year only,
    ! so do not need to read in tgg, wb, iveg, patchfrac and frac4.
    IMPLICIT NONE
    INTEGER, INTENT(IN)                  :: INpatch
    TYPE (soil_snow_type),INTENT(INOUT)       :: ssnow  ! soil and snow variables
    TYPE (bgc_pool_type),INTENT(INOUT)        :: bgc    ! carbon pool variables
    TYPE (canopy_type),INTENT(INOUT)          :: canopy ! vegetation variables
    TYPE (roughness_type),INTENT(INOUT)       :: rough  ! roughness varibles
    TYPE (balances_type),INTENT(INOUT)        :: bal ! energy + water balance variables
    REAL, INTENT(IN) :: EMSOIL
    TYPE (veg_parameter_type), INTENT(INOUT)  :: veg  ! vegetation parameters
    TYPE (soil_parameter_type), INTENT(INOUT) :: soil ! soil parameters
    TYPE (radiation_type),INTENT(INOUT)       :: rad

    ! local variables
    INTEGER, ALLOCATABLE, DIMENSION(:)   ::                                &
         nap,                                                                   &
         var_i
    REAL,    ALLOCATABLE, DIMENSION(:)   :: var_r
    REAL(r_2),    ALLOCATABLE, DIMENSION(:)   :: var_rd
    REAL,    ALLOCATABLE, DIMENSION(:,:) :: var_r2
    REAL(r_2),    ALLOCATABLE, DIMENSION(:,:) :: var_r2d
    LOGICAL                                   ::                                &
         from_restart = .TRUE.,    & ! insist variables/params load
         dummy = .TRUE.              ! To replace completeSet in parameter read; unused
    INTEGER                              :: napID

    PRINT *, '***** NOTE: now in extraRestart. *****'
    ALLOCATE(nap(mland))
    ok = NF90_INQ_VARID(ncid_rin,'nap',napID)
    IF(ok /= NF90_NOERR) CALL nc_abort                                          &
         (ok,'Error finding number of active patches in restart file '          &
         //TRIM(filename%restart_in)//' (SUBROUTINE extraRestart)')
    ok=NF90_GET_VAR(ncid_rin,napID,nap)
    IF(ok/=NF90_NOERR) CALL nc_abort(ok,'Error reading nap in file '            &
         //TRIM(filename%restart_in)// '(SUBROUTINE extraRestart)')

    ALLOCATE(var_i(INpatch))
    ALLOCATE(var_r(INpatch))
    ALLOCATE(var_rd(INpatch))
    ALLOCATE(var_r2(INpatch,msn))
    ALLOCATE(var_r2d(INpatch,ms))
    CALL readpar(ncid_rin,'wbice',dummy,var_r2d,filename%restart_in,            &
         max_vegpatches,'msd',from_restart,INpatch)
    CALL redistr_r2d(INpatch,nap,var_r2d,ssnow%wbice,'wbice',ms)
    CALL readpar(ncid_rin,'gammzz',dummy,var_r2d,filename%restart_in,           &
         max_vegpatches,'msd',from_restart,INpatch)
    CALL redistr_r2d(INpatch,nap,var_r2d,ssnow%gammzz,'gammzz',ms)
    CALL readpar(ncid_rin,'tss',dummy,var_r,filename%restart_in,                &
         max_vegpatches,'def',from_restart,INpatch)
    CALL redistr_r(INpatch,nap,var_r,ssnow%tss,'tss')
    CALL readpar(ncid_rin,'ssdnn',dummy,var_r,filename%restart_in,              &
         max_vegpatches,'def',from_restart,INpatch)
    CALL redistr_r(INpatch,nap,var_r,ssnow%ssdnn,'ssdnn')
    CALL readpar(ncid_rin,'ssdn',dummy,var_r2,filename%restart_in,              &
         max_vegpatches,'snow',from_restart,INpatch)
    CALL redistr_r2(INpatch,nap,var_r2,ssnow%ssdn,'ssdn',msn)
    CALL readpar(ncid_rin,'osnowd',dummy,var_r,filename%restart_in,             &
         max_vegpatches,'def',from_restart,INpatch)
    CALL redistr_r(INpatch,nap,var_r,ssnow%osnowd,'osnowd')
    CALL readpar(ncid_rin,'smass',dummy,var_r2,filename%restart_in,             &
         max_vegpatches,'snow',from_restart,INpatch)
    CALL redistr_r2(INpatch,nap,var_r2,ssnow%smass,'smass',msn)
    CALL readpar(ncid_rin,'sdepth',dummy,var_r2,filename%restart_in,            &
         max_vegpatches,'snow',from_restart,INpatch)
    CALL redistr_r2(INpatch,nap,var_r2,ssnow%sdepth,'sdepth',msn)
    CALL readpar(ncid_rin,'tggsn',dummy,var_r2,filename%restart_in,             &
         max_vegpatches,'snow',from_restart,INpatch)
    CALL redistr_r2(INpatch,nap,var_r2,ssnow%tggsn,'tggsn',msn)
    CALL readpar(ncid_rin,'snage',dummy,var_r,filename%restart_in,              &
         max_vegpatches,'def',from_restart,INpatch)
    CALL redistr_r(INpatch,nap,var_r,ssnow%snage,'snage')
    CALL readpar(ncid_rin,'snowd',dummy,ssnow%snowd,filename%restart_in,        &
         max_vegpatches,'def',from_restart,INpatch)
    CALL redistr_r(INpatch,nap,var_r,ssnow%snowd,'snowd')
    CALL readpar(ncid_rin,'rtsoil',dummy,var_r,filename%restart_in,             &
         max_vegpatches,'def',from_restart,INpatch)
    CALL redistr_r(INpatch,nap,var_r,ssnow%rtsoil,'rtsoil')
    CALL readpar(ncid_rin,'isflag',dummy,var_i,filename%restart_in,             &
         max_vegpatches,'def',from_restart,INpatch)
    CALL redistr_i(INpatch,nap,var_i,ssnow%isflag,'isflag')

    DEALLOCATE(var_r2)
    ALLOCATE(var_r2(INpatch,nrb))
    CALL readpar(ncid_rin,'albsoilsn',dummy,var_r2,filename%restart_in,         &
         max_vegpatches,'nrb',from_restart,INpatch)
    CALL redistr_r2(INpatch,nap,var_r2,ssnow%albsoilsn,'albsoilsn',nrb)
    ssnow%albsoilsn(:,3) = 1.0 - emsoil  !! (BP Nov 2009)
    ! The following two restart file additions are to initialise Mk3L:
    CALL readpar(ncid_rin,'albedo',dummy,var_r2,filename%restart_in,            &
         max_vegpatches,'nrb',from_restart,INpatch)
    CALL redistr_r2(INpatch,nap,var_r2,rad%albedo,'albedo',nrb)
    CALL readpar(ncid_rin,'trad',dummy,var_r,filename%restart_in,               &
         max_vegpatches,'def',from_restart,INpatch)
    CALL redistr_r(INpatch,nap,var_r,rad%trad,'trad')

    CALL readpar(ncid_rin,'rnof1',dummy,var_r,filename%restart_in,              &
         max_vegpatches,'def',from_restart,INpatch)
    CALL redistr_r(INpatch,nap,var_r,ssnow%rnof1,'rnof1')
    CALL readpar(ncid_rin,'rnof2',dummy,var_r,filename%restart_in,              &
         max_vegpatches,'def',from_restart,INpatch)
    CALL redistr_r(INpatch,nap,var_r,ssnow%rnof2,'rnof2')
    CALL readpar(ncid_rin,'runoff',dummy,var_r,filename%restart_in,             &
         max_vegpatches,'def',from_restart,INpatch)
    CALL redistr_r(INpatch,nap,var_r,ssnow%runoff,'runoff')
    CALL readpar(ncid_rin,'cansto',dummy,var_r,filename%restart_in,             &
         max_vegpatches,'def',from_restart,INpatch)
    CALL redistr_r(INpatch,nap,var_r,canopy%cansto,'cansto')
    CALL readpar(ncid_rin,'sghflux',dummy,var_r,filename%restart_in,            &
         max_vegpatches,'def',from_restart,INpatch)
    CALL redistr_r(INpatch,nap,var_r,canopy%sghflux,'sghflux')
    CALL readpar(ncid_rin,'ghflux',dummy,var_r,filename%restart_in,             &
         max_vegpatches,'def',from_restart,INpatch)
    CALL redistr_r(INpatch,nap,var_r,canopy%ghflux,'ghflux')
    CALL readpar(ncid_rin,'ga',dummy,var_r,filename%restart_in,                 &
         max_vegpatches,'def',from_restart,INpatch)
    CALL redistr_r(INpatch,nap,var_r,canopy%ga,'ga')
    CALL readpar(ncid_rin,'dgdtg',dummy,var_rd,filename%restart_in,             &
         max_vegpatches,'def',from_restart,INpatch)
    CALL redistr_rd(INpatch,nap,var_rd,canopy%dgdtg,'dgdtg')
    CALL readpar(ncid_rin,'fev',dummy,var_r,filename%restart_in,                &
         max_vegpatches,'def',from_restart,INpatch)
    CALL redistr_r(INpatch,nap,var_r,canopy%fev,'fev')
    !jhan:hack - elimiinate call as r_2 now
    !    CALL readpar(ncid_rin,'fes',dummy,var_r,filename%restart_in,               &
    !         max_vegpatches,'def',from_restart,INpatch)
    !    CALL redistr_r(INpatch,nap,var_r,canopy%fes,'fes')
    CALL readpar(ncid_rin,'fhs',dummy,var_r,filename%restart_in,                &
         max_vegpatches,'def',from_restart,INpatch)
    CALL redistr_r(INpatch,nap,var_r,canopy%fhs,'fhs')

    DEALLOCATE(var_r2)
    ALLOCATE(var_r2(INpatch,ncp))
    CALL readpar(ncid_rin,'cplant',dummy,var_r2,filename%restart_in,            &
         max_vegpatches,'ncp',from_restart,INpatch)
    CALL redistr_r2(INpatch,nap,var_r2,bgc%cplant,'cplant',ncp)

    DEALLOCATE(var_r2)
    ALLOCATE(var_r2(INpatch,ncs))
    CALL readpar(ncid_rin,'csoil',dummy,var_r2,filename%restart_in,             &
         max_vegpatches,'ncs',from_restart,INpatch)
    CALL redistr_r2(INpatch,nap,var_r2,bgc%csoil,'csoil',ncs)
    CALL readpar(ncid_rin,'wbtot0',dummy,var_r,filename%restart_in,             &
         max_vegpatches,'def',from_restart,INpatch)
    CALL redistr_r(INpatch,nap,var_r,bal%wbtot0,'wbtot0')
    CALL readpar(ncid_rin,'osnowd0',dummy,bal%osnowd0,filename%restart_in,      &
         max_vegpatches,'def',from_restart,INpatch)
    CALL redistr_r(INpatch,nap,var_r,bal%osnowd0,'osnowd0')

    ! assume all soil and veg parameters are done in default_parameters
    ! therefore, no need to do it here again

    veg%xalbnir = 1.0   ! xalbnir will soon be removed totally

    PRINT *, 'Finished extraRestart'

    DEALLOCATE(var_i)
    DEALLOCATE(var_r)
    DEALLOCATE(var_rd)
    DEALLOCATE(var_r2)
    DEALLOCATE(var_r2d)

  END SUBROUTINE extraRestart

END MODULE cable_init_module
