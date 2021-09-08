MODULE zenith_mod_cbl

REAL, ALLOCATABLE :: latitude(:,:) ! CABLE needs lat for zenith calc per dt   

CONTAINS

SUBROUTINE latitudeFix_cbl( latitude )
USE theta_field_sizes,        ONLY: t_i_length, t_j_length
USE model_grid_mod, ONLY: grid_lat => latitude
IMPLICIT NONE
REAL, ALLOCATABLE :: latitude(:,:)
LOGICAL, SAVE :: first_call = .TRUE.

IF ( first_call) THEN 
  ALLOCATE( latitude(t_i_length,t_j_length) )
  WRITE(6,*) "    "
  WRITE(6,*) "To hyper-accurately test CABLE-JULES vs CABLE-CABLE we need to "
  WRITE(6,*) "use the same zenith angle. It seems that JULES does not        "
  WRITE(6,*) "remember *latitude* past the first timestep anyway, making it  "
  WRITE(6,*) "less intrusive just to include CABLE's calc of zenith angle    "
  WRITE(6,*) "    "
  
  latitude = grid_lat
  first_call = .FALSE.
END IF

END SUBROUTINE latitudeFix_cbl

SUBROUTINE zenith_cbl( cos_zenith )
USE datetime_mod,       ONLY: l_360, l_leap
USE datetime_utils_mod, ONLY: day_of_year
USE model_time_mod,     ONLY: current_time
USE theta_field_sizes,  ONLY: row_length => t_i_length,                       &
                              rows       => t_j_length

USE cbl_sinbet_mod, ONLY: sinbet

IMPLICIT NONE
REAL, INTENT(OUT) :: cos_zenith(row_length, rows) ! Cosine of zenith angle
INTEGER ::  curr_day_number
REAL :: current_hour
INTEGER :: i, j

curr_day_number = day_of_year( current_time%year, current_time%month,         &
                               current_time%day, l_360, l_leap )
current_hour = (current_time%time / 3600.0) + .25

  ! Elemental function, in this case over (row_length,rows)
DO i = 1,  row_length
  DO j = 1,  rows
    cos_zenith(i,j) = sinbet( REAL(curr_day_number),latitude(i,j),            &
                              REAL(current_hour) ) 
    WRITE(6,*) "lat ij_cbl(i,j) ", latitude(i,j)
  END DO
END DO

END SUBROUTINE zenith_cbl

END MODULE zenith_mod_cbl
