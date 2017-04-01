MODULE pop_mpi

  USE POP_Types
  USE POP_Constants, ONLY: NCOHORT_MAX, NLAYER, HEIGHT_BINS, NDISTURB, NPATCH, NPATCH2D, &
       NYEAR_HISTORY, AGEMAX

  ! Total number of type_landscape variables to be communicated
  INTEGER, PARAMETER :: n_landscape_types = 48

  ! Total number of type_patch variables to be communicated
  INTEGER, PARAMETER :: n_patch_types     = 28

  ! Total number of type_layer variables to be communicated
  INTEGER, PARAMETER :: n_layer_types     = 6

  ! Total number of type_cohort variables to be communicated
  INTEGER, PARAMETER :: n_cohort_types    = 21

CONTAINS

  ! create MPI datatype that describes a variable of type cohort
  !
  SUBROUTINE create_cohort (cohort_t, comm)

    USE MPI

    IMPLICIT NONE

    ! the new MPI derived datatype:
    INTEGER, INTENT(OUT) :: cohort_t

    ! communicator for error-messages
    INTEGER, INTENT(IN)  :: comm

    ! temp instance of Cohort for computing displacements
    TYPE(Cohort) :: tmp_coh(2)

    ! temp variables for computing displacements and extents
    INTEGER(KIND=MPI_ADDRESS_KIND) :: a1, a2

    ! temp variable for setting a type's extent
    INTEGER(KIND=MPI_ADDRESS_KIND) :: text

    ! temp variable for lower bound parameter when setting extent
    INTEGER(KIND=MPI_ADDRESS_KIND) :: lb

    INTEGER :: tmp_t

    ! displacement, block length and block type arrays for
    ! for all fields in Type(Patch)
    INTEGER(KIND=MPI_ADDRESS_KIND),DIMENSION(:),ALLOCATABLE :: disp
    INTEGER,DIMENSION(:),ALLOCATABLE :: blen, btype

    INTEGER :: ierr, bidx

    ALLOCATE( disp (n_cohort_types) )
    ALLOCATE( blen (n_cohort_types) )
    ALLOCATE( btype(n_cohort_types) )

    bidx = 0
    lb   = 0

    ! all displacements computed relative to start of first cohort
    CALL MPI_Get_Address (tmp_coh(1), a1, ierr)

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_coh(1)%id, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_INTEGER

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_coh(1)%age, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_INTEGER

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_coh(1)%biomass, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_coh(1)%density, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_coh(1)%frac_resource_uptake, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_coh(1)%frac_light_uptake, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_coh(1)%frac_interception, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_coh(1)%frac_respiration, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_coh(1)%frac_NPP, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_coh(1)%respiration_scalar, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_coh(1)%crown_area, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_coh(1)%Pgap, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_coh(1)%height, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_coh(1)%diameter, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_coh(1)%sapwood, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_coh(1)%heartwood, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_coh(1)%sapwood_area, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_coh(1)%basal_area, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_coh(1)%LAI, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_coh(1)%Cleaf, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_coh(1)%Croot, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    IF ( bidx .NE. n_cohort_types ) THEN
       PRINT*,"Error in pop_mpi layer. bidx ",bidx, " != n_cohort_types ",&
            n_cohort_types
       CALL MPI_ABORT(comm, 0, ierr)
    ENDIF

    CALL MPI_Type_create_struct (n_cohort_types, blen, disp, btype, tmp_t, ierr)
    CALL MPI_Type_commit (tmp_t, ierr)

    ! make sure the type has correct extent for use in arrays
    CALL MPI_Get_Address (tmp_coh(2), a2, ierr)
    text = a2 - a1
    CALL MPI_Type_create_resized (tmp_t, lb, text, cohort_t, ierr)
    CALL MPI_Type_commit (cohort_t, ierr)

    DEALLOCATE( disp  )
    DEALLOCATE( blen  )
    DEALLOCATE( btype )

    RETURN

  END SUBROUTINE create_cohort

  ! create MPI datatype that describes a variable of type layer
  !
  SUBROUTINE create_layer (layer_t, comm)

    USE MPI

    IMPLICIT NONE

    ! the new MPI derived datatype:
    INTEGER, INTENT(OUT) :: layer_t

    ! communicator for error-messages
    INTEGER, INTENT(IN)  :: comm

    ! temp instance of Cohort for computing displacements
    TYPE(Layer) :: tmp_layer(2)

    ! temp variables for computing displacements and extents
    INTEGER(KIND=MPI_ADDRESS_KIND) :: a1, a2

    ! temp variable for setting a type's extent
    INTEGER(KIND=MPI_ADDRESS_KIND) :: text

    ! temp variable for lower bound parameter when setting extent
    INTEGER(KIND=MPI_ADDRESS_KIND) :: lb

    INTEGER :: tmp_t, cohort_t

    ! displacement, block length and block type arrays for
    ! for all fields in Type(Layer)
    INTEGER(KIND=MPI_ADDRESS_KIND),DIMENSION(:),ALLOCATABLE :: disp
    INTEGER,DIMENSION(:),ALLOCATABLE :: blen, btype

    INTEGER :: ierr, bidx

    ALLOCATE( disp (n_layer_types) )
    ALLOCATE( blen (n_layer_types) )
    ALLOCATE( btype(n_layer_types) )

    lb   = 0
    bidx = 0

    ! create MPI derived datatype for Type(Cohort)
    CALL create_cohort (cohort_t, comm)

    CALL MPI_Get_Address (tmp_layer(1), a1, ierr)

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_layer(1)%cohort, a2, ierr)
    disp (bidx) = a2 - a1
    ! always send all cohort array, even if ncohort < ncohort_max
    ! it's a bit inefficient, but makes the code much simpler
    blen (bidx) = NCOHORT_MAX
    btype(bidx) = cohort_t

    ! Scalar INTEGER

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_layer(1)%ncohort, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_INTEGER

    ! Scalar REAL

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_layer(1)%biomass, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_layer(1)%density, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_layer(1)%hmean, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_layer(1)%hmax, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    IF ( bidx .NE. n_layer_types ) THEN
       PRINT*,"Error in pop_mpi layer. bidx ",bidx, " != n_layer_types ",&
            n_layer_types
       CALL MPI_ABORT(comm, 0, ierr)
    ENDIF

    CALL MPI_Type_create_struct (n_layer_types, blen, disp, btype, tmp_t, ierr)
    CALL MPI_Type_commit (tmp_t, ierr)

    ! make sure the type has correct extent for use in arrays
    CALL MPI_Get_Address (tmp_layer(2), a2, ierr)
    text = a2 - a1
    CALL MPI_Type_create_resized (tmp_t, lb, text, layer_t, ierr)
    CALL MPI_Type_commit (layer_t, ierr)

    DEALLOCATE( disp  )
    DEALLOCATE( blen  )
    DEALLOCATE( btype )

    RETURN

  END SUBROUTINE create_layer

  ! create MPI datatype that describes a variable of type patch
  !
  SUBROUTINE create_patch (patch_t, comm)

    USE MPI

    IMPLICIT NONE

    ! the new MPI derived datatype:
    INTEGER, INTENT(OUT) :: patch_t

    ! communicator for error-messages
    INTEGER, INTENT(IN)  :: comm

    ! temp instance of Cohort for computing displacements
    TYPE(patch) :: tmp_patch(2)

    ! temp variables for computing displacements and extents
    INTEGER(KIND=MPI_ADDRESS_KIND) :: a1, a2

    ! temp variable for setting a type's extent
    INTEGER(KIND=MPI_ADDRESS_KIND) :: text

    ! temp variable for lower bound parameter when setting extent
    INTEGER(KIND=MPI_ADDRESS_KIND) :: lb

    INTEGER :: tmp_t, layer_t

    ! displacement, block length and block type arrays for
    ! for all fields in Type(Patch)
    INTEGER(KIND=MPI_ADDRESS_KIND),DIMENSION(:),ALLOCATABLE :: disp
    INTEGER,DIMENSION(:),ALLOCATABLE :: blen, btype

    INTEGER :: ierr, bidx

    ALLOCATE( disp (n_patch_types) )
    ALLOCATE( blen (n_patch_types) )
    ALLOCATE( btype(n_patch_types) )

    bidx = 0
    lb   = 0

    ! create MPI derived datatype for Type(Layer)
    CALL create_layer (layer_t, comm)

    CALL MPI_Get_Address (tmp_patch(1), a1, ierr)

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_patch(1)%Layer, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = NLAYER
    btype(bidx) = layer_t

    ! Scalar REAL

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_patch(1)%factor_recruit, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_patch(1)%pgap, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_patch(1)%lai, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_patch(1)%biomass, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_patch(1)%biomass_old, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_patch(1)%sapwood, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_patch(1)%heartwood, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_patch(1)%sapwood_old, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_patch(1)%sapwood_area, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_patch(1)%sapwood_area_old, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_patch(1)%stress_mortality, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_patch(1)%fire_mortality, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_patch(1)%cat_mortality, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_patch(1)%crowding_mortality, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_patch(1)%cpc, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_patch(1)%mortality, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_patch(1)%sapwood_loss, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_patch(1)%sapwood_area_loss, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_patch(1)%growth, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_patch(1)%area_growth, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_patch(1)%frac_NPP, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_patch(1)%frac_respiration, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_patch(1)%frac_light_uptake, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    ! Scalar INTEGER

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_patch(1)%id, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_INTEGER

    ! INTEGER NDISTURB

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_patch(1)%disturbance_interval(1), a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = NDISTURB
    btype(bidx) = MPI_INTEGER

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_patch(1)%first_disturbance_year(1), a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = NDISTURB
    btype(bidx) = MPI_INTEGER

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_patch(1)%age(1), a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = NDISTURB
    btype(bidx) = MPI_INTEGER

    IF ( bidx .NE. n_patch_types ) THEN
       PRINT*,"Error in pop_mpi patch. bidx ",bidx, " != n_patch_types ",&
            n_patch_types
       CALL MPI_ABORT(comm, 0, ierr)
    ENDIF

    CALL MPI_Type_create_struct (n_patch_types, blen, disp, btype, tmp_t, ierr)
    CALL MPI_Type_commit (tmp_t, ierr)

    ! make sure the type has correct extent for use in arrays
    CALL MPI_Get_Address (tmp_patch(2), a2, ierr)
    text = a2 - a1
    CALL MPI_Type_create_resized (tmp_t, lb, text, patch_t, ierr)
    CALL MPI_Type_commit (patch_t, ierr)

    DEALLOCATE( disp  )
    DEALLOCATE( blen  )
    DEALLOCATE( btype )

    RETURN

  END SUBROUTINE create_patch

  ! create MPI datatype that describes a single grid cell variable of type landscape
  !
  SUBROUTINE create_pop_gridcell_type (gcell_t, comm)

    ! Level one
    ! here Landscape-types are addressed as well as patches generated

    USE MPI

    IMPLICIT NONE

    ! the new MPI derived datatype:
    INTEGER, INTENT(OUT) :: gcell_t

    ! communicator for error-messages
    INTEGER, INTENT(IN)  :: comm

    ! temp instance of Cohort for computing displacements
    TYPE(Landscape) :: tmp_grid(2)

    ! temp variables for computing displacements and extents
    INTEGER(KIND=MPI_ADDRESS_KIND) :: a1, a2

    ! temp variable for setting a type's extent
    INTEGER(KIND=MPI_ADDRESS_KIND) :: text

    ! temp variable for lower bound parameter when setting extent
    INTEGER(KIND=MPI_ADDRESS_KIND) :: lb

    INTEGER :: tmp_t, patch_t

    ! displacement, block length and block type arrays for
    ! for all fields in Type(Patch)
    INTEGER(KIND=MPI_ADDRESS_KIND),DIMENSION(:),ALLOCATABLE :: disp
    INTEGER,DIMENSION(:),ALLOCATABLE :: blen, btype

    INTEGER :: ierr, bidx

    ALLOCATE( disp (n_landscape_types) )
    ALLOCATE( blen (n_landscape_types) )
    ALLOCATE( btype(n_landscape_types) )

    lb   = 0
    bidx = 0

    ! create MPI derived datatype for Type(Patch)
    CALL create_patch (patch_t, comm)

    ! compute displacements for the new derived type
    CALL MPI_Get_Address (tmp_grid(1), a1, ierr)

    ! NPATCH2D

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_grid(1)%patch, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = NPATCH2D
    btype(bidx) = patch_t

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_grid(1)%freq, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = NPATCH2D
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_grid(1)%freq_old, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = NPATCH2D
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_grid(1)%fire_freq, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = NPATCH2D
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_grid(1)%fire_freq_old, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = NPATCH2D
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_grid(1)%cat_freq, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = NPATCH2D
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_grid(1)%cat_freq_old, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = NPATCH2D
    btype(bidx) = MPI_DOUBLE

    ! NPATCH2D * NDISTURB

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_grid(1)%freq_ranked_age_unique, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = NPATCH2D * NDISTURB
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_grid(1)%ranked_age_unique, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = NPATCH2D * NDISTURB
    btype(bidx) = MPI_INTEGER

    ! NDISTURB

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_grid(1)%n_age, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = NPATCH2D
    btype(bidx) = MPI_INTEGER

    ! NLAYER

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_grid(1)%biomass, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = NLAYER
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_grid(1)%density, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = NPATCH2D
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_grid(1)%hmean, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = NLAYER
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_grid(1)%hmax, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = NLAYER
    btype(bidx) = MPI_DOUBLE

    ! HEIGHT_BINS

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_grid(1)%cmass_stem_bin, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = HEIGHT_BINS
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_grid(1)%densindiv_bin, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = HEIGHT_BINS
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_grid(1)%height_bin, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = HEIGHT_BINS
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_grid(1)%diameter_bin, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = HEIGHT_BINS
    btype(bidx) = MPI_DOUBLE

    ! NYEAR_HISTORY

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_grid(1)%fire_mortality_history, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = NYEAR_HISTORY
    btype(bidx) = MPI_DOUBLE

    ! AGEMAX
    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_grid(1)%freq_age, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = AGEMAX
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_grid(1)%biomass_age, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = AGEMAX
    btype(bidx) = MPI_DOUBLE

    ! Scalars REAL

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_grid(1)%cmass_sum, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_grid(1)%cmass_sum_old, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_grid(1)%cheartwood_sum, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_grid(1)%csapwood_sum, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_grid(1)%csapwood_sum_old, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_grid(1)%densindiv, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_grid(1)%height_mean, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_grid(1)%height_max, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_grid(1)%basal_area, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_grid(1)%sapwood_loss, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_grid(1)%sapwood_area_loss, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_grid(1)%stress_mortality, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_grid(1)%crowding_mortality, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_grid(1)%fire_mortality, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_grid(1)%cat_mortality, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_grid(1)%res_mortality, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_grid(1)%growth, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_grid(1)%area_growth, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_grid(1)%crown_cover, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_grid(1)%crown_area, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_grid(1)%crown_volume, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_grid(1)%sapwood_area, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_grid(1)%sapwood_area_old, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_grid(1)%Kclump, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_grid(1)%smoothing_buffer, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_grid(1)%fire_mortality_smoothed, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_DOUBLE

    ! Scalar Integer

!!$    bidx = bidx + 1
!!$    CALL MPI_Get_Address (tmp_grid(1)%npatch_active, a2, ierr)
!!$    disp (bidx) = a2 - a1
!!$    blen (bidx) = 1
!!$    btype(bidx) = MPI_INTEGER

    bidx = bidx + 1
    CALL MPI_Get_Address (tmp_grid(1)%LU, a2, ierr)
    disp (bidx) = a2 - a1
    blen (bidx) = 1
    btype(bidx) = MPI_INTEGER


    IF ( bidx .NE. n_landscape_types ) THEN
       PRINT*,"Error in pop_mpi landscape. bidx ",bidx, " != n_landscape_types ",&
            n_landscape_types
       CALL MPI_ABORT(comm, 0, ierr)
    ENDIF

    CALL MPI_Type_create_struct (n_landscape_types, blen, disp, btype, tmp_t, ierr)
    CALL MPI_Type_commit (tmp_t, ierr)

    ! make sure the type has correct extent for use in arrays
    CALL MPI_Get_Address (tmp_grid(2), a2, ierr)
    text = a2 - a1
    CALL MPI_Type_create_resized (tmp_t, lb, text, gcell_t, ierr)
    CALL MPI_Type_commit (gcell_t, ierr)

    DEALLOCATE( disp )
    DEALLOCATE( blen )
    DEALLOCATE( btype)

    RETURN

  END SUBROUTINE create_pop_gridcell_type

END MODULE pop_mpi

