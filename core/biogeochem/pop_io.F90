SUBROUTINE POP_IO( POP, casamet, YEAR, ACTION, CF )
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! POP        : POP structure containing all specific parameter
  ! casamet    : structure containing met and grid specific parameters from CASA
  ! YEAR       : Current year <YYYY>
  ! ACTION     : What do you want?
  !              "READ_RST"  : Read a restart file (will be looking for a file
  !                            either with given name or from YEAR-1
  !              "WRITE_RST" : Write a restart file for YEAR+1
  !              "WRITE_EPI" : Write data at the end of each year
  ! CLOSE_FILE : Flag to close file at the end of Episode (Episode only)
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  USE netcdf
  USE TypeDef, only: dp
  USE POP_constants
  USE POP_types
  USE CASAVARIABLE, only: casa_met, icycle
  USE CABLE_COMMON_MODULE

  IMPLICIT NONE

  TYPE(POP_TYPE),   INTENT(INOUT) :: POP
  TYPE(casa_met),   INTENT(IN)    :: casamet
  INTEGER       ,   INTENT(IN)    :: YEAR
  CHARACTER(LEN=9), INTENT(IN)    :: ACTION
  LOGICAL,          INTENT(IN)    :: CF

  INTEGER          :: STATUS, i, m, p, l, land_ID, patch_ID, ndis_ID
  !CRM  INTEGER    :: ndis1_ID,nlay_ID,hgtb_ID,ncoh_ID,t_ID
  INTEGER          :: nlay_ID, hgtb_ID, ncoh_ID, t_ID
  INTEGER          :: nlayer_dim, ndisturb_dim, land_dim
  INTEGER          :: HEIGHT_BINS_dim, npatch2d_dim, NCOHORT_MAX_dim
  INTEGER          :: dID, t_dim, tx=-1, ntile, mp, CNT
  CHARACTER(len=3) :: typ='rst'
  CHARACTER        :: dum*9, fname*120
  LOGICAL          :: CLOSE_FILE, EXISTFILE

  !   ! 1 dim arrays (np)
  !   CHARACTER(len=40),DIMENSION( 2), PARAMETER :: AR0 = (/'latitude','longitude'/)

  !   ! LANDSCAPE STRUCTURE
  !   ! 2 dim arrays (np,t)
  !   CHARACTER(len=40),DIMENSION( 1), PARAMETER :: AI1 = (/ 'npatch_active' /)
  !   CHARACTER(len=40),DIMENSION(11), PARAMETER :: AR1 = (/ 'cmass_sum',           &
  !        'densindiv','height_mean','height_max','basal_area','stress_mortality',  &
  !        'fire_mortality','growth','crown_cover','crown_area','crown_volume' /)
  !   ! 3 dim arrays (np,nlayer,t)
  !   CHARACTER(len=40),DIMENSION( 4), PARAMETER :: AR2 = (/ 'biomass','density',   &
  !        'hmean','hmax' /)
  !   ! 3 dim arrays (np,height_bins,t)
  !   CHARACTER(len=40),DIMENSION( 4), PARAMETER :: AR3 = (/ 'cmass_stem_bin',      &
  !        'densindiv_bin','height_bin','diameter_bin' /)
  !   ! 3 dim arrays (np,ndisturb,t)
  !   CHARACTER(len=40),DIMENSION( 1), PARAMETER :: AI4 = (/ 'n_age' /)

  !   ! PATCH STRUCTURE
  !   ! 3 dim arrays (np,npatch2d,t)
  !   CHARACTER(len=40),DIMENSION( 1), PARAMETER :: AI5 = (/ 'patch_id' /)
  !   CHARACTER(len=40),DIMENSION(10), PARAMETER :: AR5 = (/ 'patch_freq',          &
  !        'patch_freq_old','patch_freq_old2','patch_factor_recruit',               &
  !        'patch_biomass','patch_biomass_old','patch_biomass_old2',                &
  !        'patch_stress_mortality','patch_fire_mortality','patch_growth' /)
  !   ! 4 dim arrays (np,npatch2d,ndisturb+1,t)
  ! !CRM  CHARACTER(len=40),DIMENSION( 1), PARAMETER :: AI6 =                           &
  ! !CRM       (/ 'patch_ranked_age_unique' /)
  ! !CRM  CHARACTER(len=40),DIMENSION( 1), PARAMETER :: AR6 =                           &
  ! !CRM       (/ 'patch_freq_ranked_age_unique' /)
  !   ! 4 dim arrays (np,npatch2d,ndisturb,t)
  !   CHARACTER(len=40),DIMENSION( 5), PARAMETER :: AI7 =                           &
  !        (/ 'patch_disturbance_interval','patch_first_disturbance_year',          &
  !        'patch_age','patch_age_old','patch_ranked_age_unique' /)
  !   CHARACTER(len=40),DIMENSION( 1), PARAMETER :: AR7 =                           &
  !        (/ 'patch_freq_ranked_age_unique' /)

  !   ! LAYER STRUCTURE
  !   ! 4 dim arrays (np,npatch2d,nlayer,t)
  !   CHARACTER(len=40),DIMENSION( 4), PARAMETER :: AR8 = (/ 'layer_biomass',       &
  !        'layer_density','layer_hmean','layer_hmax' /)
  !   CHARACTER(len=40),DIMENSION( 1), PARAMETER :: AI8 = (/ 'layer_ncohort' /)

  !   ! COHORT STRUCTURE
  !   ! 5 dim arrays (np,npatch2d,nlayer,ncohort_max,t)
  !   CHARACTER(len=40),DIMENSION( 5), PARAMETER :: AR9 = (/ 'cohort_biomass',      &
  !        'cohort_density','cohort_frac_resource_uptake','cohort_height',          &
  !        'cohort_diameter' /)
  !   CHARACTER(len=40),DIMENSION( 2), PARAMETER :: AI9 = (/'cohort_age','cohort_id'/)

  ! 1 dim arrays (np)
  CHARACTER(len=40), DIMENSION( 2) :: AR0
  CHARACTER(len=40), DIMENSION( 2) :: AI0

  ! LANDSCAPE STRUCTURE
  ! 2 dim arrays (np,t)
  CHARACTER(len=40), DIMENSION( 1) :: AI1
  CHARACTER(len=40), DIMENSION(24) :: AR1
  ! 3 dim arrays (np,nlayer,t)
  CHARACTER(len=40), DIMENSION( 4) :: AR2
  ! 3 dim arrays (np,height_bins,t)
  CHARACTER(len=40), DIMENSION( 4) :: AR3
  ! 3 dim arrays (np,ndisturb,t)
  CHARACTER(len=40), DIMENSION( 1) :: AI4

  ! PATCH STRUCTURE
  ! 3 dim arrays (np,npatch2d,t)
  CHARACTER(len=40), DIMENSION( 1) :: AI5
  CHARACTER(len=40), DIMENSION(24) :: AR5
  ! 4 dim arrays (np,npatch2d,ndisturb,t)
  CHARACTER(len=40), DIMENSION( 4) :: AI7
  CHARACTER(len=40), DIMENSION( 1) :: AR7

  ! LAYER STRUCTURE
  ! 4 dim arrays (np,npatch2d,nlayer,t)
  CHARACTER(len=40), DIMENSION( 4) :: AR8
  CHARACTER(len=40), DIMENSION( 1) :: AI8

  ! COHORT STRUCTURE
  ! 5 dim arrays (np,npatch2d,nlayer,ncohort_max,t)
  CHARACTER(len=40), DIMENSION(19) :: AR9
  CHARACTER(len=40), DIMENSION( 2) :: AI9

  INTEGER, SAVE :: VIDtime, &
       VIDR0(SIZE(AR0)), VIDI0(SIZE(AI0)), VIDR1(SIZE(AR1)), VIDI1(SIZE(AI1)), &
       VIDR2(SIZE(AR2)), VIDR3(SIZE(AR3)), VIDI4(SIZE(AI4)), VIDR5(SIZE(AR5)), &
       VIDI5(SIZE(AI5)), VIDI7(SIZE(AI7)), VIDR7(SIZE(AR7)), VIDR8(SIZE(AR8)), &
       VIDI8(SIZE(AI8)), VIDR9(SIZE(AR9)), VIDI9(SIZE(AI9))
  INTEGER, SAVE :: FILE_ID, EPI_CNT = 0

  ! TEMPORARY ARRAYS
  INTEGER,  ALLOCATABLE :: I1(:), I2(:,:), I3(:,:,:), I4(:,:,:,:)
  REAL(dp), ALLOCATABLE :: R1(:), R2(:,:), R3(:,:,:), R4(:,:,:,:)

  mp = POP%np

  AR0(1)  = 'latitude'
  AR0(2)  = 'longitude'

  AI0(1)  = 'Iwood'
  AI0(2)  = 'it_pop' ! Scalar value !!!

  AI1(1)  = 'npatch_active'

  AR1(1)  = 'cmass_sum'
  AR1(2)  = 'cmass_sum_old'
  AR1(3)  = 'cheartwood_sum'
  AR1(4)  = 'csapwood_sum'
  AR1(5)  = 'csapwood_sum_old'
  AR1(6)  = 'densindiv'
  AR1(7)  = 'height_mean'
  AR1(8)  = 'height_max'
  AR1(9)  = 'basal_area'
  AR1(10) = 'sapwood_loss'
  AR1(11) = 'sapwood_area_loss'
  AR1(12) = 'stress_mortality'
  AR1(13) = 'crowding_mortality'
  AR1(14) = 'fire_mortality'
  AR1(15) = 'cat_mortality'
  AR1(16) = 'res_mortality'
  AR1(17) = 'growth'
  AR1(18) = 'area_growth'
  AR1(19) = 'crown_cover'
  AR1(20) = 'crown_area'
  AR1(21) = 'crown_volume'
  AR1(22) = 'sapwood_area'
  AR1(23) = 'sapwood_area_old'
  AR1(24) = 'Kclump'

  AR2(1)  = 'biomass'
  AR2(2)  = 'density'
  AR2(3)  = 'hmean'
  AR2(4)  = 'hmax'

  AR3(1)  = 'cmass_stem_bin'
  AR3(2)  = 'densindiv_bin'
  AR3(3)  = 'height_bin'
  AR3(4)  = 'diameter_bin'

  AI4(1)  = 'n_age'

  AI5(1)  = 'patch_id'

  AR5(1)  = 'patch_freq'
  AR5(2)  = 'patch_freq_old'
  AR5(3)  = 'patch_factor_recruit'
  AR5(4)  = 'patch_pgap'
  AR5(5)  = 'patch_lai'
  AR5(6)  = 'patch_biomass'
  AR5(7)  = 'patch_biomass_old'
  AR5(8)  = 'patch_sapwood'
  AR5(9)  = 'patch_heartwood'
  AR5(10) = 'patch_sapwood_old'
  AR5(11) = 'patch_sapwood_area'
  AR5(12) = 'patch_sapwood_area_old'
  AR5(13) = 'patch_stress_mortality'
  AR5(14) = 'patch_fire_mortality'
  AR5(15) = 'patch_cat_mortality'
  AR5(16) = 'patch_crowding_mortality'
  AR5(17) = 'patch_cpc'
  AR5(18) = 'patch_sapwood_loss'
  AR5(19) = 'patch_sapwood_area_loss'
  AR5(20) = 'patch_growth'
  AR5(21) = 'patch_area_growth'
  AR5(22) = 'patch_frac_NPP'
  AR5(23) = 'patch_frac_respiration'
  AR5(24) = 'patch_frac_GPP'

  AI7(1)  = 'patch_disturbance_interval'
  AI7(2)  = 'patch_first_disturbance_year'
  AI7(3)  = 'patch_age'
  AI7(4)  = 'patch_ranked_age_unique'

  AR7(1)  = 'patch_freq_ranked_age_unique'

  AR8(1)  = 'layer_biomass'
  AR8(2)  = 'layer_density'
  AR8(3)  = 'layer_hmean'
  AR8(4)  = 'layer_hmax'


  AI8(1)  = 'layer_ncohort'

  AR9(1)  = 'cohort_biomass'
  AR9(2)  = 'cohort_density'
  AR9(3)  = 'cohort_frac_resource_uptake'
  AR9(4)  = 'cohort_frac_light_uptake'
  AR9(5)  = 'cohort_frac_interception'
  AR9(6)  = 'cohort_frac_respiration'
  AR9(7)  = 'cohort_frac_NPP'
  AR9(8)  = 'cohort_respiration_scalar'
  AR9(9)  = 'cohort_crown_area'
  AR9(10) = 'cohort_Pgap'
  AR9(11) = 'cohort_height'
  AR9(12) = 'cohort_diameter'
  AR9(13) = 'cohort_sapwood'
  AR9(14) = 'cohort_heartwood'
  AR9(15) = 'cohort_sapwood_area'
  AR9(16) = 'cohort_basal_area'
  AR9(17) = 'cohort_LAI'
  AR9(18) = 'cohort_Cleaf'
  AR9(19) = 'cohort_Croot'

  AI9(1)  = 'cohort_age'
  AI9(2)  = 'cohort_id'

  ntile  = mp

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! IF ( PRESENT( CF ) ) THEN
  CLOSE_FILE = CF
  ! ELSE
  !    CLOSE_FILE = .FALSE.
  ! END IF

  ! Check for valid ACTION
  IF ( INDEX(ACTION,"WRITE_EPI") .GT. 0 ) THEN
     typ = 'out'
  ELSE IF ( INDEX(ACTION,"WRITE_RST") .GT. 0 ) THEN
     typ = 'rst'
  ELSE IF ( INDEX(ACTION,"WRITE_INI") .GT. 0 ) THEN
     typ = 'ini'
  ELSE IF ( INDEX(ACTION,"READ_rst") .GT. 0 ) THEN
     typ = 'rst'
  ELSE
     WRITE(*,*)  "WRONG ACTION:'", TRIM(ACTION), "' in call to pop_io!"
     STOP -1
  ENDIF

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! WRITE POP VALUES TO OUTPUT FILE
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Get File-Name
  IF ( typ .EQ. 'out' ) THEN
     WRITE(dum, FMT="(I4,'_',I4)") CABLE_USER%YEARSTART, CABLE_USER%YEAREND
     IF (CABLE_USER%YEARSTART.lt.1000.and.CABLE_USER%YEAREND.lt.1000) THEN
        WRITE(dum, FMT="(I3,'_',I3)") CABLE_USER%YEARSTART, CABLE_USER%YEAREND
     ELSEIF (CABLE_USER%YEARSTART.lt.1000) THEN
        WRITE(dum, FMT="(I3,'_',I4)") CABLE_USER%YEARSTART, CABLE_USER%YEAREND
     ENDIF
  ELSE
     WRITE(dum, FMT="(I4)") YEAR
  ENDIF

  ! IF (typ.eq.'ini') THEN
  !    fname = TRIM(cable_user%POP_rst)//'/'//'pop_'//TRIM(cable_user%RunIDEN)&
  !         //'_'//typ//'.nc'
  ! ELSE
  !    fname = TRIM(filename%path)//'/'//TRIM(cable_user%RunIden)//'_'//&
  !         TRIM(dum)//'_pop_'//typ//'.nc'
  ! ENDIF
  IF ((typ.eq.'ini') .OR. (typ.eq.'rst')) THEN
     IF (LEN_TRIM( TRIM(cable_user%POP_restart_out) ) .gt. 0) THEN
        fname = TRIM(cable_user%POP_restart_out)
     ELSE
        fname = TRIM(filename%path)//'/'//'pop_'//TRIM(cable_user%RunIDEN)//'_'//typ//'.nc'
     ENDIF
  ELSE
     IF (LEN_TRIM( TRIM(cable_user%POP_outfile) ) .gt. 0) THEN
        fname = TRIM(cable_user%POP_outfile)
     ELSE
        fname = TRIM(filename%path)//'/'//TRIM(cable_user%RunIden)//'_'//TRIM(dum)//'_pop_'//typ//'.nc'
     ENDIF
  ENDIF

  IF ( INDEX(ACTION,"WRITE") .GT. 0 ) THEN

     IF ( typ .EQ. 'out' ) THEN
        EPI_CNT = EPI_CNT + 1
        CNT     = EPI_CNT
     ELSE
        CNT = 1
     ENDIF

     IF ( CNT .EQ. 1 ) THEN

        INQUIRE( FILE=TRIM( fname ), EXIST=EXISTFILE )
        EXISTFILE = .FALSE.

        IF ( EXISTFILE .and. (typ.ne.'ini') .and. (typ.ne.'rst') ) THEN  ! file exists

           STATUS = NF90_open(trim(fname), mode=nf90_write, ncid=FILE_ID)
           ! print*, 'OOpen70.1 ', file_id, trim(fname)
           IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

           STATUS = nf90_inq_dimid(FILE_ID, 'time', t_id)
           IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

           !CLN status = nf90_inquire_dimension(FILE_ID, t_id,name = RecordDimName, len = CNT)
           !CLN if (status /= nf90_noerr) call handle_err(status)
           !CLN CNT = CNT+1
           ! Enquire variable IDs
           STATUS = nf90_inq_varid(FILE_ID, 'Time', VIDTime)
           IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

           DO i = 1, SIZE(AR0)
              STATUS = nf90_inq_varid(FILE_ID, TRIM(AR0(i)), VIDR0(i))
              IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
           END DO
           DO i = 1, SIZE(AI0)
              STATUS = nf90_inq_varid(FILE_ID, TRIM(AI0(i)), VIDI0(i))
              IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
           END DO
           DO i = 1, SIZE(AI1)
              STATUS = nf90_inq_varid(FILE_ID, TRIM(AI1(i)), VIDI1(i))
              IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
           END DO
           DO i = 1, SIZE(AR1)
              STATUS = nf90_inq_varid(FILE_ID, TRIM(AR1(i)), VIDR1(i))
              IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
           END DO
           DO i = 1, SIZE(AR2)
              STATUS = nf90_inq_varid(FILE_ID, TRIM(AR2(i)), VIDR2(i))
              IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
           END DO
           DO i = 1, SIZE(AR3)
              STATUS = nf90_inq_varid(FILE_ID, TRIM(AR3(i)), VIDR3(i))
              IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
           END DO
           DO i = 1, SIZE(AI4)
              STATUS = nf90_inq_varid(FILE_ID, TRIM(AI4(i)), VIDI4(i))
              IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
           END DO
           DO i = 1, SIZE(AI5)
              STATUS = nf90_inq_varid(FILE_ID, TRIM(AI5(i)), VIDI5(i))
              IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
           END DO
           DO i = 1, SIZE(AR5)
              STATUS = nf90_inq_varid(FILE_ID, TRIM(AR5(i)), VIDR5(i))
              IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
           END DO
           DO i = 1, SIZE(AI7)
              STATUS = nf90_inq_varid(FILE_ID, TRIM(AI7(i)), VIDI7(i))
              IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
           END DO
           DO i = 1, SIZE(AR7)
              STATUS = nf90_inq_varid(FILE_ID, TRIM(AR7(i)), VIDR7(i))
              IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
           END DO
           DO i = 1, SIZE(AI8)
              STATUS = nf90_inq_varid(FILE_ID, TRIM(AI8(i)), VIDI8(i))
              IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
           END DO
           DO i = 1, SIZE(AR8)
              STATUS = nf90_inq_varid(FILE_ID, TRIM(AR8(i)), VIDR8(i))
              IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
           END DO
           DO i = 1, SIZE(AI9)
              STATUS = nf90_inq_varid(FILE_ID, TRIM(AI9(i)), VIDI9(i))
              IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
           END DO
           DO i = 1, SIZE(AR9)
              STATUS = nf90_inq_varid(FILE_ID, TRIM(AR9(i)), VIDR9(i))
              IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
           END DO

        ELSE  ! file doesn't already exist, or is RST or INI

           ! Create NetCDF file:
           STATUS = NF90_create(trim(fname), cmode=ior(nf90_clobber,nf90_64bit_offset), ncid=FILE_ID)
           ! print*, 'OCreate70 ', file_id, trim(fname)
           IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

           ! GLOBAL ATTRIBUTES
           STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "Icycle" , icycle             )
           STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "Year"   , YEAR               )
           STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "RunIden", CABLE_USER%RunIden )

           ! Define dimensions:
           STATUS = NF90_def_dim(FILE_ID, 'land'       , mp         , land_ID )
           IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
           STATUS = NF90_def_dim(FILE_ID, 'NPATCH2D'   , NPATCH2D   , patch_ID)
           IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
           STATUS = NF90_def_dim(FILE_ID, 'NDISTURB'   , NDISTURB   , ndis_ID )
           IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
           !CRM     STATUS = NF90_def_dim(FILE_ID, 'NDISTURB+1' , NDISTURB+1 , ndis1_ID)
           !CRM     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
           STATUS = NF90_def_dim(FILE_ID, 'NLAYER'     , NLAYER     , nlay_ID )
           IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
           STATUS = NF90_def_dim(FILE_ID, 'HEIGHT_BINS', HEIGHT_BINS, hgtb_ID )
           IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
           STATUS = NF90_def_dim(FILE_ID, 'NCOHORT_MAX', NCOHORT_MAX, ncoh_ID )
           IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
           IF ( typ .EQ. 'out' ) THEN
              STATUS = NF90_def_dim(FILE_ID, 'time'   ,  NF90_UNLIMITED, t_ID    )
              IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
           ELSE
              STATUS = NF90_def_dim(FILE_ID, 'time'   ,  1, t_ID    )
              IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
           ENDIF

           ! Define variables
           STATUS = NF90_def_var(FILE_ID,'Time' , NF90_INT, (/t_ID/), VIDtime )
           IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

           DO i = 1, SIZE(AR0)
              STATUS = NF90_def_var(FILE_ID,TRIM(AR0(i)), NF90_DOUBLE, (/land_ID/), VIDR0(i))
              IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
           END DO
           DO i = 1, SIZE(AI0)
              STATUS = NF90_def_var(FILE_ID,TRIM(AI0(i)), NF90_INT, (/land_ID/), VIDI0(i))
              IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
           END DO
           DO i = 1, SIZE(AI1)
              STATUS = NF90_def_var(FILE_ID,TRIM(AI1(i)), NF90_INT, (/land_ID,t_ID/), VIDI1(i))
              IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
           END DO
           DO i = 1, SIZE(AR1)
              STATUS = NF90_def_var(FILE_ID,TRIM(AR1(i)), NF90_DOUBLE, (/land_ID,t_ID/), VIDR1(i))
              IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
           END DO
           DO i = 1, SIZE(AR2)
              STATUS = NF90_def_var(FILE_ID,TRIM(AR2(i)), NF90_DOUBLE, (/land_ID,nlay_ID,t_ID/), VIDR2(i))
              IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
           END DO

           DO i = 1, SIZE(AR3)
              STATUS = NF90_def_var(FILE_ID,TRIM(AR3(i)), NF90_DOUBLE, (/land_ID,hgtb_ID,t_ID/), VIDR3(i))
              IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
           END DO
           DO i = 1, SIZE(AI4)
              STATUS = NF90_def_var(FILE_ID,TRIM(AI4(i)), NF90_INT, (/land_ID,ndis_ID,t_ID/), VIDI4(i))
              IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
           END DO
           DO i = 1, SIZE(AI5)
              STATUS = NF90_def_var(FILE_ID,TRIM(AI5(i)), NF90_INT, (/land_ID,patch_ID,t_ID/), VIDI5(i))
              IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
           END DO
           DO i = 1, SIZE(AR5)
              STATUS = NF90_def_var(FILE_ID,TRIM(AR5(i)), NF90_DOUBLE, (/land_ID,patch_ID,t_ID/), VIDR5(i))
              IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
           END DO
           DO i = 1, SIZE(AI7)
              STATUS = NF90_def_var(FILE_ID,TRIM(AI7(i)), NF90_INT, &
                   (/land_ID,patch_ID,ndis_ID,t_ID/), VIDI7(i))
              IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
           END DO

           DO i = 1, SIZE(AR7)
              STATUS = NF90_def_var(FILE_ID,TRIM(AR7(i)), NF90_DOUBLE, &
                   (/land_ID,patch_ID,ndis_ID,t_ID/), VIDR7(i))
              IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
           END DO
           DO i = 1, SIZE(AI8)
              STATUS = NF90_def_var(FILE_ID,TRIM(AI8(i)), NF90_INT, &
                   (/land_ID,patch_ID,nlay_ID,t_ID/), VIDI8(i))
              IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
           END DO
           DO i = 1, SIZE(AR8)
              STATUS = NF90_def_var(FILE_ID,TRIM(AR8(i)), NF90_DOUBLE, &
                   (/land_ID,patch_ID,nlay_ID,t_ID/), VIDR8(i))
              IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
           END DO
           DO i = 1, SIZE(AI9)
              STATUS = NF90_def_var(FILE_ID,TRIM(AI9(i)), NF90_INT, &
                   (/land_ID,patch_ID,nlay_ID,ncoh_ID,t_ID/), VIDI9(i))
              IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
           END DO
           DO i = 1, SIZE(AR9)
              STATUS = NF90_def_var(FILE_ID,TRIM(AR9(i)), NF90_DOUBLE, &
                   (/land_ID,patch_ID,nlay_ID,ncoh_ID,t_ID/), VIDR9(i))
              IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
           END DO

           ! End define mode:
           STATUS = NF90_enddef(FILE_ID)
           IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

           ! PUT LAT / LON ( np )
           STATUS = NF90_PUT_VAR(FILE_ID, VIDR0(1), casamet%lat(POP%Iwood), &
                start=(/ 1 /), count=(/ mp /) )
           IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

           STATUS = NF90_PUT_VAR(FILE_ID, VIDR0(2), casamet%lon(POP%Iwood), &
                start=(/ 1 /), count=(/ mp /) )
           IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

           ! PUT Iwood
           STATUS = NF90_PUT_VAR(FILE_ID, VIDI0(1), POP%Iwood, &
                start=(/ 1 /), count=(/ mp /) )
           IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

           ! PUT it_pop
           STATUS = NF90_PUT_VAR(FILE_ID, VIDI0(2), POP%it_pop, &
                start=(/ 1 /), count=(/ mp /) )
           IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

        END IF ! file exists
     END IF

     ! WRITE CURRENT STATE
     ! TIME  ( t )
     ! print*, 'OWrite70 ', file_id
     STATUS = NF90_PUT_VAR(FILE_ID, VIDtime, YEAR, start=(/ CNT /) )
     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

     ! PUT 2D VARS ( mp, t )
     STATUS = NF90_PUT_VAR(FILE_ID, VIDI1( 1), POP%pop_grid(:)%npatch_active, &
          start=(/ 1, CNT /), count=(/ mp, 1 /) )
     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

     STATUS = NF90_PUT_VAR(FILE_ID, VIDR1( 1), POP%pop_grid(:)%cmass_sum, &
          start=(/ 1, CNT /), count=(/ mp, 1 /) )
     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

     STATUS = NF90_PUT_VAR(FILE_ID, VIDR1( 2), POP%pop_grid(:)%cmass_sum_old, &
          start=(/ 1, CNT /), count=(/ mp, 1 /) )
     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

     STATUS = NF90_PUT_VAR(FILE_ID, VIDR1( 3), POP%pop_grid(:)%cheartwood_sum, &
          start=(/ 1, CNT /), count=(/ mp, 1 /) )
     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

     STATUS = NF90_PUT_VAR(FILE_ID, VIDR1( 4), POP%pop_grid(:)%csapwood_sum, &
          start=(/ 1, CNT /), count=(/ mp, 1 /) )
     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

     STATUS = NF90_PUT_VAR(FILE_ID, VIDR1( 5), POP%pop_grid(:)%csapwood_sum_old,  &
          start=(/ 1, CNT /), count=(/ mp, 1 /) )
     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

     STATUS = NF90_PUT_VAR(FILE_ID, VIDR1( 6), POP%pop_grid(:)%densindiv,         &
          start=(/ 1, CNT /), count=(/ mp, 1 /) )
     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

     STATUS = NF90_PUT_VAR(FILE_ID, VIDR1( 7), POP%pop_grid(:)%height_mean,       &
          start=(/ 1, CNT /), count=(/ mp, 1 /) )
     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

     STATUS = NF90_PUT_VAR(FILE_ID, VIDR1( 8), POP%pop_grid(:)%height_max,        &
          start=(/ 1, CNT /), count=(/ mp, 1 /) )
     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

     STATUS = NF90_PUT_VAR(FILE_ID, VIDR1( 9), POP%pop_grid(:)%basal_area,        &
          start=(/ 1, CNT /), count=(/ mp, 1 /) )
     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

     STATUS = NF90_PUT_VAR(FILE_ID, VIDR1( 10), POP%pop_grid(:)%sapwood_loss,     &
          start=(/ 1, CNT /), count=(/ mp, 1 /) )
     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

     STATUS = NF90_PUT_VAR(FILE_ID, VIDR1( 11), POP%pop_grid(:)%sapwood_area_loss,          &
          start=(/ 1, CNT /), count=(/ mp, 1 /) )
     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

     STATUS = NF90_PUT_VAR(FILE_ID, VIDR1( 12), POP%pop_grid(:)%stress_mortality,   &
          start=(/ 1, CNT /), count=(/ mp, 1 /) )
     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

     STATUS = NF90_PUT_VAR(FILE_ID, VIDR1( 13), POP%pop_grid(:)%crowding_mortality,   &
          start=(/ 1, CNT /), count=(/ mp, 1 /) )
     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

     STATUS = NF90_PUT_VAR(FILE_ID, VIDR1( 14), POP%pop_grid(:)%fire_mortality,     &
          start=(/ 1, CNT /), count=(/ mp, 1 /) )
     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

     STATUS = NF90_PUT_VAR(FILE_ID, VIDR1( 15), POP%pop_grid(:)%cat_mortality,     &
          start=(/ 1, CNT /), count=(/ mp, 1 /) )
     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

     STATUS = NF90_PUT_VAR(FILE_ID, VIDR1( 16), POP%pop_grid(:)%res_mortality,     &
          start=(/ 1, CNT /), count=(/ mp, 1 /) )
     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

     STATUS = NF90_PUT_VAR(FILE_ID, VIDR1( 17), POP%pop_grid(:)%growth,             &
          start=(/ 1, CNT /), count=(/ mp, 1 /) )
     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

     STATUS = NF90_PUT_VAR(FILE_ID, VIDR1( 18), POP%pop_grid(:)%area_growth,             &
          start=(/ 1, CNT /), count=(/ mp, 1 /) )
     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

     STATUS = NF90_PUT_VAR(FILE_ID, VIDR1( 19), POP%pop_grid(:)%crown_cover,        &
          start=(/ 1, CNT /), count=(/ mp, 1 /) )
     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

     STATUS = NF90_PUT_VAR(FILE_ID, VIDR1(20), POP%pop_grid(:)%crown_area,         &
          start=(/ 1, CNT /), count=(/ mp, 1 /) )
     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

     STATUS = NF90_PUT_VAR(FILE_ID, VIDR1(21), POP%pop_grid(:)%crown_volume,       &
          start=(/ 1, CNT /), count=(/ mp, 1 /) )
     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

     STATUS = NF90_PUT_VAR(FILE_ID, VIDR1(22), POP%pop_grid(:)%sapwood_area,         &
          start=(/ 1, CNT /), count=(/ mp, 1 /) )
     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

     STATUS = NF90_PUT_VAR(FILE_ID, VIDR1(23), POP%pop_grid(:)%sapwood_area_old,         &
          start=(/ 1, CNT /), count=(/ mp, 1 /) )
     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

     STATUS = NF90_PUT_VAR(FILE_ID, VIDR1(24), POP%pop_grid(:)%KClump,         &
          start=(/ 1, CNT /), count=(/ mp, 1 /) )
     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

     ! PUT 3D VARS ( mp,nlayer, t )

     ALLOCATE( R2( mp, nlayer ) )
     !  DO i = 1, SIZE(AR2)
     DO i = 1, SIZE(VIDR2)
        DO m = 1, mp
           SELECT CASE (i)
           CASE(1)
              R2(m,:) = POP%pop_grid(m)%biomass
           CASE(2)
              R2(m,:) = POP%pop_grid(m)%density
           CASE(3)
              R2(m,:) = POP%pop_grid(m)%hmean
           CASE(4)
              R2(m,:) = POP%pop_grid(m)%hmax
           CASE default
              STOP "Parameter not assigned in pop_bios_io.f90!"
           END SELECT
        END DO
        STATUS = NF90_PUT_VAR(FILE_ID, VIDR2( i), R2,         &
             start=(/ 1, 1, CNT /), count=(/ mp, nlayer, 1 /) )
        IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
     END DO
     DEALLOCATE( R2 )

     ! PUT 3D VARS ( mp,height_bins, t )

     ALLOCATE( R2( mp, height_bins ) )
     DO i = 1, SIZE(VIDR3)
        DO m = 1, mp
           SELECT CASE(i)
           CASE(1)
              R2(m,:) = POP%pop_grid(m)%cmass_stem_bin
           CASE(2)
              R2(m,:) = POP%pop_grid(m)%densindiv_bin
           CASE(3)
              R2(m,:) = POP%pop_grid(m)%height_bin
           CASE(4)
              R2(m,:) = POP%pop_grid(m)%diameter_bin
           CASE default
              STOP "Parameter not assigned in pop_bios_io.f90!"
           END SELECT
        END DO
        STATUS = NF90_PUT_VAR(FILE_ID, VIDR3( i), R2,         &
             start=(/ 1, 1, CNT /), count=(/ mp, height_bins, 1 /) )
        IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
     END DO
     DEALLOCATE( R2)

     ! PUT 3D VARS ( mp,ndisturb, t )

     ALLOCATE ( I2( mp, ndisturb ) )
     DO i = 1, SIZE(VIDI4)
        DO m = 1, mp
           SELECT CASE (i)
           CASE(1)
              I2( m,: ) = POP%pop_grid(m)%n_age
           CASE default
              STOP "Parameter not assigned in pop_bios_io.f90!"
           END SELECT
        END DO
        STATUS = NF90_PUT_VAR(FILE_ID, VIDI4( i), I2,         &
             start=(/ 1, 1, CNT /), count=(/ mp, ndisturb, 1 /) )
        IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
     END DO
     DEALLOCATE( I2 )

     ! PUT 3D VARS ( mp,npatch2d, t )

     ALLOCATE ( I2( mp, npatch2d ) )
     DO i = 1, SIZE(VIDI5)
        DO m = 1, mp
           SELECT CASE (i)
           CASE( 1)
              I2( m,: ) = POP%pop_grid(m)%patch(:)%id
           CASE default
              STOP "Parameter not assigned in pop_bios_io.f90!"
           END SELECT
        END DO
        STATUS = NF90_PUT_VAR(FILE_ID, VIDI5( i), I2,         &
             start=(/ 1, 1, CNT /), count=(/ mp, npatch2d, 1 /) )
        IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
     END DO
     DEALLOCATE( I2 )

     ALLOCATE( R2( mp, npatch2d ) )
     DO i = 1, SIZE(VIDR5)
        DO m = 1, mp
           SELECT CASE ( i )
           CASE( 1)
              R2(m,:)=POP%pop_grid(m)%freq
           CASE( 2)
              R2(m,:)=POP%pop_grid(m)%freq_old
           CASE( 3)
              R2(m,:)=POP%pop_grid(m)%patch(:)%factor_recruit
           CASE( 4)
              R2(m,:)=POP%pop_grid(m)%patch(:)%pgap
           CASE( 5)
              R2(m,:)=POP%pop_grid(m)%patch(:)%lai
           CASE( 6)
              R2(m,:)=POP%pop_grid(m)%patch(:)%biomass
           CASE( 7)
              R2(m,:)=POP%pop_grid(m)%patch(:)%biomass_old
           CASE( 8)
              R2(m,:)=POP%pop_grid(m)%patch(:)%sapwood
           CASE( 9)
              R2(m,:)=POP%pop_grid(m)%patch(:)%heartwood
           CASE( 10)
              R2(m,:)=POP%pop_grid(m)%patch(:)%sapwood_old
           CASE( 11)
              R2(m,:)=POP%pop_grid(m)%patch(:)%sapwood_area
           CASE( 12)
              R2(m,:)=POP%pop_grid(m)%patch(:)%sapwood_area_old
           CASE( 13)
              R2(m,:)=POP%pop_grid(m)%patch(:)%stress_mortality
           CASE( 14)
              R2(m,:)=POP%pop_grid(m)%patch(:)%fire_mortality
           CASE( 15)
              R2(m,:)=POP%pop_grid(m)%patch(:)%cat_mortality
           CASE( 16)
              R2(m,:)=POP%pop_grid(m)%patch(:)%crowding_mortality
           CASE( 17)
              R2(m,:)=POP%pop_grid(m)%patch(:)%cpc
           CASE( 18)
              R2(m,:)=POP%pop_grid(m)%patch(:)%sapwood_loss
           CASE( 19)
              R2(m,:)=POP%pop_grid(m)%patch(:)%sapwood_area_loss
           CASE( 20)
              R2(m,:)=POP%pop_grid(m)%patch(:)%growth
           CASE( 21)
              R2(m,:)=POP%pop_grid(m)%patch(:)%area_growth
           CASE( 22)
              R2(m,:)=POP%pop_grid(m)%patch(:)%frac_NPP
           CASE( 23)
              R2(m,:)=POP%pop_grid(m)%patch(:)%frac_respiration
           CASE( 24)
              R2(m,:)=POP%pop_grid(m)%patch(:)%frac_light_uptake
           CASE default
              STOP "Parameter not assigned in pop_bios_io.f90!"
           END SELECT
        END DO
        STATUS = NF90_PUT_VAR(FILE_ID, VIDR5( i), R2,         &
             start=(/ 1, 1, CNT /), count=(/ mp, npatch2d, 1 /) )
        IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
     END DO
     DEALLOCATE( R2 )

     ! PUT 4D VARS ( mp,npatch2d, ndisturb,t )
     ALLOCATE( I3( mp, npatch2d, ndisturb ) )
     DO i = 1, SIZE(VIDI7)
        DO p = 1, npatch2d
           DO m = 1, mp
              SELECT CASE ( i )
              CASE( 1)
                 I3(m,p,:)= POP%pop_grid(m)%patch(p)%disturbance_interval
              CASE( 2)
                 I3(m,p,:)= POP%pop_grid(m)%patch(p)%first_disturbance_year
              CASE( 3)
                 I3(m,p,:)= POP%pop_grid(m)%patch(p)%age
              CASE( 4)
                 I3(m,p,:)= POP%pop_grid(m)%ranked_age_unique(p,:)
              CASE default
                 STOP "Parameter not assigned in pop_bios_io.f90!"
              END SELECT
           END DO
        END DO
        STATUS = NF90_PUT_VAR(FILE_ID, VIDI7( i), I3,         &
             start=(/1, 1, 1, CNT /), count=(/ mp, npatch2d,ndisturb, 1 /) )
        IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
     END DO
     DEALLOCATE( I3 )

     ALLOCATE( R3( mp, npatch2d, ndisturb ) )
     DO i = 1, SIZE(VIDR7)

        DO m = 1, mp
           SELECT CASE ( i )
           CASE( 1)
              R3(m,:,:)= POP%pop_grid(m)%freq_ranked_age_unique
           CASE default
              STOP "Parameter not assigned in pop_bios_io.f90!"
           END SELECT
        END DO

        STATUS = NF90_PUT_VAR(FILE_ID, VIDR7( i), R3,         &
             start=(/1, 1, 1, CNT /), count=(/ mp, npatch2d,ndisturb, 1 /) )
        IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
     END DO
     DEALLOCATE( R3 )

     !PUT 4D VARS ( mp,npatch2d, nlayer,t )
     ALLOCATE( I3( mp, npatch2d, nlayer ) )
     DO i = 1, SIZE(VIDI8)
        DO p = 1, npatch2d
           DO m = 1, mp
              SELECT CASE ( i )
              CASE( 1)
                 I3(m,p,:) =  POP%pop_grid(m)%patch(p)%layer(:)%ncohort
              CASE default
                 STOP "Parameter not assigned in pop_bios_io.f90!"
              END SELECT
           END DO
        END DO
        STATUS = NF90_PUT_VAR(FILE_ID, VIDI8( i), I3,         &
             start=(/1, 1, 1, CNT /), count=(/ mp, npatch2d,nlayer, 1 /) )
        IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
     END DO
     DEALLOCATE( I3 )

     ALLOCATE( R3( mp, npatch2d, nlayer ) )
     DO i = 1, SIZE(VIDR8)
        DO p = 1, npatch2d
           DO m = 1, mp
              SELECT CASE ( i )
              CASE( 1)
                 R3(m,p,:)=POP%pop_grid(m)%patch(p)%layer(:)%biomass
              CASE( 2)
                 R3(m,p,:)=POP%pop_grid(m)%patch(p)%layer(:)%density
              CASE( 3)
                 R3(m,p,:)=POP%pop_grid(m)%patch(p)%layer(:)%hmean
              CASE( 4)
                 R3(m,p,:)=POP%pop_grid(m)%patch(p)%layer(:)%hmax
              CASE default
                 STOP "Parameter not assigned in pop_bios_io.f90!"
              END SELECT
           END DO
        END DO
        STATUS = NF90_PUT_VAR(FILE_ID, VIDR8( i), R3,         &
             start=(/1, 1, 1, CNT /), count=(/ mp, npatch2d,nlayer, 1 /) )
        IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
     END DO
     DEALLOCATE( R3 )

     ! PUT 5D VARS ( mp,npatch2d, nlayer,ncohort_max,t )
     ALLOCATE( I4( mp, npatch2d, nlayer, ncohort_max ) )
     DO i = 1, SIZE(VIDI9)
        DO l = 1, nlayer
           DO p = 1, npatch2d
              DO m = 1, mp
                 SELECT CASE ( i )
                 CASE( 1)
                    I4(m,p,l,:) = POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%age
                 CASE( 2)
                    I4(m,p,l,:) = POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%id
                 CASE default
                    STOP "Parameter not assigned in pop_bios_io.f90!"
                 END SELECT
              END DO
           END DO
        END DO
        STATUS = NF90_PUT_VAR(FILE_ID, VIDI9( i), I4,         &
             start=(/ 1, 1,1,1, CNT /), count=(/ mp, npatch2d,nlayer,ncohort_max,1 /) )
        IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

     ENDDO
     DEALLOCATE(I4)

     ALLOCATE( R4( mp, npatch2d, nlayer, ncohort_max ) )
     DO i = 1, SIZE(VIDR9)
        DO l = 1, nlayer
           DO p = 1, npatch2d
              DO m = 1, mp
                 SELECT CASE ( i )
                 CASE( 1)
                    R4(m,p,l,:)=POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%biomass
                 CASE( 2)
                    R4(m,p,l,:)=POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%density
                 CASE( 3)
                    R4(m,p,l,:)= &
                         POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%frac_resource_uptake
                 CASE( 4)
                    R4(m,p,l,:)= &
                         POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%frac_light_uptake
                 CASE( 5)
                    R4(m,p,l,:)= &
                         POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%frac_interception
                 CASE( 6)
                    R4(m,p,l,:)= &
                         POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%frac_respiration
                 CASE( 7)
                    R4(m,p,l,:)=POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%frac_NPP
                 CASE( 8)
                    R4(m,p,l,:)  = &
                         POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%respiration_scalar
                 CASE( 9)
                    R4(m,p,l,:) =POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%crown_area
                 CASE( 10)
                    R4(m,p,l,:) =POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%Pgap
                 CASE( 11)
                    R4(m,p,l,:) =POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%height
                 CASE( 12)
                    R4(m,p,l,:) =POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%diameter
                 CASE( 13)
                    R4(m,p,l,:) =POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%sapwood
                 CASE( 14)
                    R4(m,p,l,:) =POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%heartwood
                 CASE( 15)
                    R4(m,p,l,:) =POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%sapwood_area
                 CASE( 16)
                    R4(m,p,l,:) =POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%basal_area
                 CASE( 17)
                    R4(m,p,l,:) =POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%LAI
                 CASE( 18)
                    R4(m,p,l,:) =POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%Cleaf
                 CASE( 19)
                    R4(m,p,l,:) =POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%Croot
                 CASE default
                    STOP "Parameter not assigned in pop_bios_io.f90!"
                 END SELECT
              END DO
           END DO
        END DO
        STATUS = NF90_PUT_VAR(FILE_ID, VIDR9( i), R4,         &
             start=(/ 1, 1,1,1, CNT /), count=(/ mp, npatch2d,nlayer,ncohort_max,1 /) )
        IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

     ENDDO
     DEALLOCATE(R4)

     ! ! PUT 3D VARS ( mp,nlayer, t )
     ! MPS:DO m = 1, mp


     !    PAT:DO p = 1, npatch2d


     !       STATUS = NF90_PUT_VAR(FILE_ID, VIDR7( 1), POP%pop_grid(m)%freq_ranked_age_unique(p,:),&
     !            start=(/ m, p, 1, CNT /), count=(/ 1, 1, NDISTURB, 1 /) )
     !       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

     !       ! LAYER STRUCTURE
     !       ! PUT 4D VARS ( mp,npatch2d, nlayer,t )
     !       STATUS = NF90_PUT_VAR(FILE_ID, VIDI8( 1), POP%pop_grid(m)%patch(p)%layer(:)%ncohort,&
     !            start=(/ m, p, 1, CNT /), count=(/ 1, 1, nlayer, 1 /) )
     !       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)




     !       LAY:DO l = 1, nlayer
     !          ! COHORT STRUCTURE
     !          ! PUT 5D VARS ( mp,npatch2d, nlayer,ncohort_max,t )
     !          STATUS = NF90_PUT_VAR(FILE_ID, VIDI9( 1), POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%age,&
     !               start=(/ m, p, l, 1, CNT /), count=(/ 1, 1, 1, ncohort_max, 1 /) )
     !          IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
     !          STATUS = NF90_PUT_VAR(FILE_ID, VIDI9( 2), POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%id,&
     !               start=(/ m, p, l, 1, CNT /), count=(/ 1, 1, 1, ncohort_max, 1 /) )
     !          IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)


     !       END DO LAY
     !    END DO PAT
     ! END DO MPS

     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! READ POP VALUES AS RESTART VALUES
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ELSE IF ( INDEX(ACTION,'READ') .GT. 0 ) THEN

     IF (LEN_TRIM(TRIM(cable_user%POP_restart_in)).gt.0) THEN
        fname = TRIM(cable_user%POP_restart_in)
     ELSE
        WRITE( dum, FMT="(I4)")YEAR-1
        fname = TRIM(cable_user%POP_rst)//'/'//TRIM(dum)//'_pop_'//TRIM(cable_user%RunIDEN)//'_'//typ//'.nc'
     ENDIF
     INQUIRE( FILE=TRIM(fname), EXIST=EXISTFILE )
     ! If suitable restart-file, try ini-restart
     IF ( .NOT. EXISTFILE ) THEN
        WRITE(*,*) "Restart file not found: ",TRIM(fname)
        WRITE(*,*) "Looking for initialization file..."
        fname = TRIM(cable_user%POP_rst)//'/'//'pop_'//TRIM(cable_user%RunIDEN)//'_ini.nc'
     ENDIF
     INQUIRE( FILE=TRIM(fname), EXIST=EXISTFILE )
     IF (.NOT. EXISTFILE) THEN
        WRITE(*,*) " No ini-restart file found either! ", TRIM(fname)
        STOP -1
     ELSE
        typ = "ini"
     ENDIF

     WRITE(*,*)"Reading POP-rst file: ", TRIM(fname)

     STATUS = NF90_OPEN( TRIM(fname), NF90_NOWRITE, FILE_ID )
     ! print*, 'OOpen70.2 ', file_id, TRIM(fname)
     IF (STATUS /= NF90_noerr)THEN
        WRITE(*,*)"Error opening file (pop_bios_io.f90) ",TRIM(fname)
        CALL handle_err(STATUS)
     ENDIF
     ! DIMS
     STATUS = NF90_INQ_DIMID( FILE_ID, 'land', dID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     STATUS = NF90_INQUIRE_DIMENSION( FILE_ID, dID, LEN=land_dim )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

     STATUS = NF90_INQ_DIMID( FILE_ID, 'NPATCH2D', dID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     STATUS = NF90_INQUIRE_DIMENSION( FILE_ID, dID, LEN=npatch2d_dim )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

     STATUS = NF90_INQ_DIMID( FILE_ID, 'NDISTURB', dID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     STATUS = NF90_INQUIRE_DIMENSION( FILE_ID, dID, LEN=NDISTURB_dim )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

     !CRM   STATUS = NF90_INQ_DIMID( FILE_ID, 'NDISTURB+1', dID )
     !CRM   IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     !CRM   STATUS = NF90_INQUIRE_DIMENSION( FILE_ID, dID, LEN=NDISTURB1_dim )
     !CRM   IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

     STATUS = NF90_INQ_DIMID( FILE_ID, 'NLAYER', dID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     STATUS = NF90_INQUIRE_DIMENSION( FILE_ID, dID, LEN=NLAYER_dim )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

     STATUS = NF90_INQ_DIMID( FILE_ID, 'HEIGHT_BINS', dID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     STATUS = NF90_INQUIRE_DIMENSION( FILE_ID, dID, LEN=HEIGHT_BINS_dim )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

     STATUS = NF90_INQ_DIMID( FILE_ID, 'NCOHORT_MAX', dID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     STATUS = NF90_INQUIRE_DIMENSION( FILE_ID, dID, LEN=NCOHORT_MAX_dim )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

     IF ( land_dim .NE. mp .OR.  npatch2d_dim .NE. NPATCH2D .OR.  &
          HEIGHT_BINS_dim .NE. HEIGHT_BINS .OR. NCOHORT_MAX_dim .NE. NCOHORT_MAX &
          .OR. NLAYER_dim .NE. NLAYER .OR. NDISTURB_dim .NE. NDISTURB ) THEN
        WRITE(*,*)"Dimension misfit in pop_bios_io.f90!"
        WRITE(*,*)"Restart file  | Current Run"
        WRITE(*,*)"# points   ",land_dim,"     ",mp
        WRITE(*,*)"# patches  ",NPATCH2D_dim,"     ",NPATCH2D
        WRITE(*,*)"# HGT_BINS ",HEIGHT_BINS_dim,"     ",HEIGHT_BINS
        WRITE(*,*)"NCOHORT_MAX",NCOHORT_MAX_dim,"     ",NCOHORT_MAX
        WRITE(*,*)"# NLAYER   ",NLAYER_dim,"     ",NLAYER
        WRITE(*,*)"# NDISTURB ",NDISTURB_dim,"     ",NDISTURB
        STOP
     ENDIF

     ! TIME
     STATUS = NF90_INQ_DIMID( FILE_ID, 'time', dID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     STATUS = NF90_INQUIRE_DIMENSION( FILE_ID, dID, LEN=t_dim )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

     STATUS = NF90_INQ_VARID( FILE_ID, 'Time', t_ID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

     ALLOCATE( I1 (t_dim) )
     STATUS = NF90_GET_VAR( FILE_ID, t_ID, I1 )
     IF ( STATUS /= NF90_noerr ) CALL handle_err(STATUS)
     DO i = 1, t_dim
        IF ( YEAR .EQ. I1(i)+1 ) THEN ! DATA FROM PRECEDING YEAR !
           tx = i
           EXIT
        END IF
     END DO
     DEALLOCATE( I1 )

     IF ( tx .LE. 0 ) THEN
        WRITE(*,*) 'FILE '//TRIM(fname)//" doesn't contain specific data for ",YEAR
        WRITE(*,*) 'Resetting  tx to 1!'
        tx = 1
        IF ( typ .NE. "ini" ) THEN
           WRITE(*,*) "Wrong date in input pop restart-file! ",TRIM(fname)
           STOP
        ENDIF
     ENDIF

     ! CHECK LAT 'N LON
     ALLOCATE( R1( mp ) )
     STATUS = NF90_INQ_VARID( FILE_ID, TRIM(AR0(1)), dID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     STATUS = NF90_GET_VAR ( FILE_ID, dID, R1 )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     IF ( ANY ( casamet%lat(POP%Iwood) .NE. R1 ) ) THEN
        WRITE(*,*)"INPUT LATs don't match casamet! pop_bios_io.f90" &
             , TRIM(fname)
        ! STOP
     ENDIF
     STATUS = NF90_INQ_VARID( FILE_ID, TRIM(AR0(2)), dID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     STATUS = NF90_GET_VAR ( FILE_ID, dID, R1 )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     IF ( ANY ( casamet%lon(POP%Iwood) .NE. R1 ) ) THEN
        WRITE(*,*)"INPUT LONs don't match casamet! pop_bios_io.f90" &
             , TRIM(fname)
        !STOP
     ENDIF
     DEALLOCATE ( R1 )

     ! GET 0D VARS ( np )
     STATUS = NF90_INQ_VARID( FILE_ID, TRIM(AI0(1)), dID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     STATUS = NF90_GET_VAR  ( FILE_ID, dID, POP%Iwood, &
          start=(/ 1,tx /), count=(/ mp, 1 /) )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

     STATUS = NF90_INQ_VARID( FILE_ID, TRIM(AI0(2)), dID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     STATUS = NF90_GET_VAR  ( FILE_ID, dID, POP%it_pop, &
          start=(/ 1,tx /), count=(/ mp, 1 /) )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

     ! GET 1D VARS ( np )
     STATUS = NF90_INQ_VARID( FILE_ID, TRIM(AI1(1)), dID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     STATUS = NF90_GET_VAR  ( FILE_ID, dID, POP%pop_grid(:)%npatch_active, &
          start=(/ 1,tx /), count=(/ mp, 1 /) )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

     ALLOCATE ( R1( mp ) )
     DO i = 1, SIZE(AR1)
        STATUS = NF90_INQ_VARID( FILE_ID, TRIM(AR1(i)), dID )
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        STATUS = NF90_GET_VAR  ( FILE_ID, dID, R1, start=(/1,tx/), count=(/mp,1/) )
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        SELECT CASE ( i )
        CASE( 1)
           POP%pop_grid(:)%cmass_sum          = R1
        CASE( 2)
           POP%pop_grid(:)%cmass_sum_old      = R1
        CASE( 3)
           POP%pop_grid(:)%cheartwood_sum     = R1
        CASE( 4)
           POP%pop_grid(:)%csapwood_sum       = R1
        CASE( 5)
           POP%pop_grid(:)%csapwood_sum_old   = R1
        CASE( 6)
           POP%pop_grid(:)%densindiv          = R1
        CASE( 7)
           POP%pop_grid(:)%height_mean        = R1
        CASE( 8)
           POP%pop_grid(:)%height_max         = R1
        CASE( 9)
           POP%pop_grid(:)%basal_area         = R1
        CASE(10)
           POP%pop_grid(:)%sapwood_loss       = R1
        CASE(11)
           POP%pop_grid(:)%sapwood_area_loss  = R1
        CASE(12)
           POP%pop_grid(:)%stress_mortality   = R1
        CASE(13)
           POP%pop_grid(:)%crowding_mortality = R1
        CASE(14)
           POP%pop_grid(:)%fire_mortality     = R1
        CASE(15)
           POP%pop_grid(:)%cat_mortality      = R1
        CASE(16)
           POP%pop_grid(:)%res_mortality      = R1
        CASE(17)
           POP%pop_grid(:)%growth             = R1
        CASE(18)
           POP%pop_grid(:)%area_growth        = R1
        CASE(19)
           POP%pop_grid(:)%crown_cover        = R1
        CASE(20)
           POP%pop_grid(:)%crown_area         = R1
        CASE(21)
           POP%pop_grid(:)%crown_volume       = R1
        CASE(22)
           POP%pop_grid(:)%sapwood_area       = R1
        CASE(23)
           POP%pop_grid(:)%sapwood_area_old   = R1
        CASE(24)
           POP%pop_grid(:)%KClump             = R1
        CASE default
           STOP "Parameter not assigned in pop_bios_io.f90!"
        END SELECT
     END DO
     DEALLOCATE ( R1 )

     ! GET 2D VARS ( mp,nlayer )
     ALLOCATE( R2( mp, nlayer ) )
     DO i = 1, SIZE(AR2)
        STATUS = NF90_INQ_VARID( FILE_ID, TRIM(AR2(i)), dID )
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        STATUS = NF90_GET_VAR  ( FILE_ID, dID, R2, &
             start=(/ 1, 1, tx /), count=(/ mp, nlayer, 1 /) )
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        DO m = 1, mp
           SELECT CASE ( i )
           CASE( 1)
              POP%pop_grid(m)%biomass = R2(m,:)
           CASE( 2)
              POP%pop_grid(m)%density = R2(m,:)
           CASE( 3)
              POP%pop_grid(m)%hmean   = R2(m,:)
           CASE( 4)
              POP%pop_grid(m)%hmax    = R2(m,:)
           CASE default
              STOP "Parameter not assigned in pop_bios_io.f90!"
           END SELECT
        END DO
     END DO
     DEALLOCATE( R2 )

     ! GET 2D VARS ( mp,height_bins )
     ALLOCATE( R2( mp, height_bins ) )
     DO i = 1, SIZE(AR3)
        STATUS = NF90_INQ_VARID( FILE_ID, TRIM(AR3(i)), dID )
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        STATUS = NF90_GET_VAR  ( FILE_ID, dID, R2, &
             start=(/ 1, 1, tx /), count=(/ mp, height_bins, 1 /) )
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        DO m = 1, mp
           SELECT CASE ( i )
           CASE( 1)
              POP%pop_grid(m)%cmass_stem_bin = R2(m,:)
           CASE( 2)
              POP%pop_grid(m)%densindiv_bin  = R2(m,:)
           CASE( 3)
              POP%pop_grid(m)%height_bin     = R2(m,:)
           CASE( 4)
              POP%pop_grid(m)%diameter_bin   = R2(m,:)
           CASE default
              STOP "Parameter not assigned in pop_bios_io.f90!"
           END SELECT
        END DO
     END DO
     DEALLOCATE( R2 )

     ! GET 2D VARS ( mp,ndisturb)
     ALLOCATE ( I2( mp, ndisturb ) )
     STATUS = NF90_INQ_VARID( FILE_ID, TRIM(AI4(1)), dID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     STATUS = NF90_GET_VAR  ( FILE_ID, dID, I2, &
          start=(/ 1, 1, tx /), count=(/ mp, ndisturb, 1 /) )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     DO m = 1, mp
        POP%pop_grid(m)%n_age = I2( m,: )
     END DO
     DEALLOCATE( I2 )

     ! GET 2D VARS ( mp,npatch2d)
     ALLOCATE ( I2( mp,npatch2d ) )
     STATUS = NF90_INQ_VARID( FILE_ID, TRIM(AI5(1)), dID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     STATUS = NF90_GET_VAR  ( FILE_ID, dID, I2, &
          start=(/ 1, 1, tx /), count=(/ mp, npatch2d, 1 /) )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     DO m = 1, mp
        POP%pop_grid(m)%patch(:)%id = I2( m,: )
     END DO
     DEALLOCATE( I2 )

     ! PATCH STRUCTURE
     ! GET 2D VARS ( mp,npatch2d )
     ALLOCATE( R2( mp, npatch2d ) )
     DO i = 1, SIZE(AR5)
        STATUS = NF90_INQ_VARID( FILE_ID, TRIM(AR5(i)), dID )
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        STATUS = NF90_GET_VAR  ( FILE_ID, dID, R2, &
             start=(/ 1, 1, tx /), count=(/ mp, npatch2d, 1 /) )
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        DO m = 1, mp
           SELECT CASE ( i )
           CASE( 1)
              POP%pop_grid(m)%freq                        = R2(m,:)
           CASE( 2)
              POP%pop_grid(m)%freq_old                    = R2(m,:)
           CASE( 3)
              POP%pop_grid(m)%patch(:)%factor_recruit     = R2(m,:)
           CASE( 4)
              POP%pop_grid(m)%patch(:)%pgap               = R2(m,:)
           CASE( 5)
              POP%pop_grid(m)%patch(:)%lai                = R2(m,:)
           CASE( 6)
              POP%pop_grid(m)%patch(:)%biomass            = R2(m,:)
           CASE( 7)
              POP%pop_grid(m)%patch(:)%biomass_old        = R2(m,:)
           CASE( 8)
              POP%pop_grid(m)%patch(:)%sapwood            = R2(m,:)
           CASE( 9)
              POP%pop_grid(m)%patch(:)%heartwood          = R2(m,:)
           CASE( 10)
              POP%pop_grid(m)%patch(:)%sapwood_old        = R2(m,:)
           CASE( 11)
              POP%pop_grid(m)%patch(:)%sapwood_area       = R2(m,:)
           CASE( 12)
              POP%pop_grid(m)%patch(:)%sapwood_area_old   = R2(m,:)
           CASE( 13)
              POP%pop_grid(m)%patch(:)%stress_mortality   = R2(m,:)
           CASE( 14)
              POP%pop_grid(m)%patch(:)%fire_mortality     = R2(m,:)
           CASE( 15)
              POP%pop_grid(m)%patch(:)%cat_mortality      = R2(m,:)
           CASE( 16)
              POP%pop_grid(m)%patch(:)%crowding_mortality = R2(m,:)
           CASE( 17)
              POP%pop_grid(m)%patch(:)%cpc                = R2(m,:)
           CASE( 18)
              POP%pop_grid(m)%patch(:)%sapwood_loss       = R2(m,:)
           CASE( 19)
              POP%pop_grid(m)%patch(:)%sapwood_area_loss  = R2(m,:)
           CASE( 20)
              POP%pop_grid(m)%patch(:)%growth             = R2(m,:)
           CASE( 21)
              POP%pop_grid(m)%patch(:)%area_growth        = R2(m,:)
           CASE( 22)
              POP%pop_grid(m)%patch(:)%frac_NPP           = R2(m,:)
           CASE( 23)
              POP%pop_grid(m)%patch(:)%frac_respiration   = R2(m,:)
           CASE( 24)
              POP%pop_grid(m)%patch(:)%frac_light_uptake  = R2(m,:)
           CASE default
              STOP "Parameter not assigned in pop_bios_io.f90!"
           END SELECT
        END DO
     END DO
     DEALLOCATE( R2 )

     ! GET 3D VARS ( mp,npatch2d,ndisturb )
     ALLOCATE( I3( mp, npatch2d, ndisturb ) )
     DO i = 1, SIZE(AI7)
        STATUS = NF90_INQ_VARID( FILE_ID, TRIM(AI7(i)), dID )
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        STATUS = NF90_GET_VAR  ( FILE_ID, dID, I3, &
             start=(/ 1, 1, 1, tx /), count=(/ mp, npatch2d, ndisturb, 1 /) )
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        DO p = 1, npatch2d
           DO m = 1, mp
              SELECT CASE ( i )
              CASE( 1)
                 POP%pop_grid(m)%patch(p)%disturbance_interval   = I3(m,p,:)
              CASE( 2)
                 POP%pop_grid(m)%patch(p)%first_disturbance_year = I3(m,p,:)
              CASE( 3)
                 POP%pop_grid(m)%patch(p)%age                    = I3(m,p,:)
              CASE( 4)
                 POP%pop_grid(m)%ranked_age_unique(p,:)          = I3(m,p,:)
              CASE default
                 STOP "Parameter not assigned in pop_bios_io.f90!"
              END SELECT
           END DO
        END DO
     END DO
     DEALLOCATE( I3 )

     ALLOCATE( R3( mp, npatch2d,ndisturb ) )
     STATUS = NF90_INQ_VARID( FILE_ID, TRIM(AR7(1)), dID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     STATUS = NF90_GET_VAR  ( FILE_ID, dID, R3, &
          start=(/ 1, 1, 1, tx /), count=(/ mp, npatch2d, ndisturb, 1 /) )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     DO m = 1, mp
        POP%pop_grid(m)%freq_ranked_age_unique = R3(m,:,:)
     END DO
     DEALLOCATE( R3 )

     ! LAYER STRUCTURE
     ! GET 3D VARS ( mp,npatch2d,nlayer )
     ALLOCATE( I3( mp, npatch2d, nlayer ) )
     DO i = 1, SIZE(AI8)
        STATUS = NF90_INQ_VARID( FILE_ID, TRIM(AI8(i)), dID )
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        STATUS = NF90_GET_VAR  ( FILE_ID, dID, I3, &
             start=(/ 1, 1, 1, tx /), count=(/ mp, npatch2d, nlayer, 1 /) )
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        DO p = 1, npatch2d
           DO m = 1, mp
              SELECT CASE ( i )
              CASE( 1)
                 POP%pop_grid(m)%patch(p)%layer(:)%ncohort = I3(m,p,:)
              CASE default
                 STOP "Parameter not assigned in pop_bios_io.f90!"
              END SELECT
           END DO
        END DO
     END DO
     DEALLOCATE( I3 )

     ALLOCATE( R3( mp, npatch2d, nlayer ) )
     DO i = 1, SIZE(AR8)
        STATUS = NF90_INQ_VARID( FILE_ID, TRIM(AR8(i)), dID )
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        STATUS = NF90_GET_VAR  ( FILE_ID, dID, R3, &
             start=(/ 1, 1, 1, tx /), count=(/ mp, npatch2d, nlayer, 1 /) )
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        DO p = 1, npatch2d
           DO m = 1, mp
              SELECT CASE ( i )
              CASE( 1)
                 POP%pop_grid(m)%patch(p)%layer(:)%biomass = R3(m,p,:)
              CASE( 2)
                 POP%pop_grid(m)%patch(p)%layer(:)%density = R3(m,p,:)
              CASE( 3)
                 POP%pop_grid(m)%patch(p)%layer(:)%hmean   = R3(m,p,:)
              CASE( 4)
                 POP%pop_grid(m)%patch(p)%layer(:)%hmax    = R3(m,p,:)
              CASE default
                 STOP "Parameter not assigned in pop_bios_io.f90!"
              END SELECT
           END DO
        END DO
     END DO
     DEALLOCATE( R3 )

     ! COHORT STRUCTURE
     ! GET 4D VARS ( mp,npatch2d,nlayer,ncohort_max )
     ALLOCATE( I4( mp, npatch2d, nlayer, ncohort_max ) )
     DO i = 1, SIZE(AI9)
        STATUS = NF90_INQ_VARID( FILE_ID, TRIM(AI9(i)), dID )
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        STATUS = NF90_GET_VAR  ( FILE_ID, dID, I4, start=(/ 1, 1, 1, 1, tx /), &
             count=(/ mp,npatch2d,nlayer,ncohort_max,1 /) )
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        DO l = 1, nlayer
           DO p = 1, npatch2d
              DO m = 1, mp
                 SELECT CASE ( i )
                 CASE( 1)
                    POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%age = I4(m,p,l,:)
                 CASE( 2)
                    POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%id  = I4(m,p,l,:)
                 CASE default
                    STOP "Parameter not assigned in pop_bios_io.f90!"
                 END SELECT
              END DO
           END DO
        END DO
     END DO
     DEALLOCATE( I4 )

     ALLOCATE( R4( mp, npatch2d, nlayer, ncohort_max ) )
     DO i = 1, SIZE(AR9)
        STATUS = NF90_INQ_VARID( FILE_ID, TRIM(AR9(i)), dID )
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        STATUS = NF90_GET_VAR  ( FILE_ID, dID, R4, start=(/ 1, 1, 1, 1, tx /), &
             count=(/ mp,npatch2d,nlayer,ncohort_max,1 /) )
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        DO l = 1, nlayer
           DO p = 1, npatch2d
              DO m = 1, mp
                 SELECT CASE ( i )
                 CASE( 1)
                    POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%biomass              = R4(m,p,l,:)
                 CASE( 2)
                    POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%density              = R4(m,p,l,:)
                 CASE( 3)
                    POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%frac_resource_uptake = R4(m,p,l,:)
                 CASE( 4)
                    POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%frac_light_uptake    = R4(m,p,l,:)
                 CASE( 5)
                    POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%frac_interception    = R4(m,p,l,:)
                 CASE( 6)
                    POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%frac_respiration     = R4(m,p,l,:)
                 CASE( 7)
                    POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%frac_NPP             = R4(m,p,l,:)
                 CASE( 8)
                    POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%respiration_scalar   = R4(m,p,l,:)
                 CASE( 9)
                    POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%crown_area           = R4(m,p,l,:)
                 CASE( 10)
                    POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%Pgap                 = R4(m,p,l,:)
                 CASE( 11)
                    POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%height               = R4(m,p,l,:)
                 CASE( 12)
                    POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%diameter             = R4(m,p,l,:)
                 CASE( 13)
                    POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%sapwood              = R4(m,p,l,:)
                 CASE( 14)
                    POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%heartwood            = R4(m,p,l,:)
                 CASE( 15)
                    POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%sapwood_area         = R4(m,p,l,:)
                 CASE( 16)
                    POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%basal_area           = R4(m,p,l,:)
                 CASE( 17)
                    POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%LAI                  = R4(m,p,l,:)
                 CASE( 18)
                    POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%Cleaf                = R4(m,p,l,:)
                 CASE( 19)
                    POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%Croot                = R4(m,p,l,:)
                 CASE default
                    STOP "Parameter not assigned in pop_bios_io.f90!"
                 END SELECT
              END DO
           END DO
        END DO
     END DO
     DEALLOCATE( R4 )

  ELSE
     WRITE(*,*) 'ACTION = ',TRIM(ACTION)
     STOP 'Please, enter either "READ" or "WRITE" when calling pop_bios_io.f90!'
  END IF

  IF ( CLOSE_FILE .OR. typ .EQ. 'rst' .OR. typ .EQ. 'ini' ) THEN
     ! Close NetCDF file:
     ! print*, 'OClose70 ', file_id
     STATUS = NF90_close(FILE_ID)
     file_id = -1
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     WRITE(*,*) "Closed POP-file ", TRIM(fname)
  ENDIF

END SUBROUTINE POP_IO
