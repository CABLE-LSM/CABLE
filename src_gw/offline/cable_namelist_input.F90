MODULE  cable_namelist_util
  !module contains a subroutine to get a namelist file names something
  !besides cable.nml from the first argument passed to CABLE on the
  !command line.  If the first argument does not end in nml it is assumed
  !that is specifies the met or cnppool that will be read from GETARG
  !later

  IMPLICIT NONE

  CHARACTER(:), ALLOCATABLE,SAVE :: CABLE_NAMELIST
  LOGICAL, SAVE                  :: arg_not_namelist

CONTAINS

  SUBROUTINE get_namelist_file_name()

    !LOCAL
    INTEGER :: arg_length
    INTEGER :: command_line_stat

    LOGICAL :: no_arguments
    !get namelist name

    arg_not_namelist = .FALSE.  !assume the best!

    IF (COMMAND_ARGUMENT_COUNT() == 0) THEN
       no_arguments=.TRUE.
    ELSE
       no_arguments=.FALSE.
    END IF

    IF (.NOT.no_arguments) THEN
       CALL GET_COMMAND_ARGUMENT(1, LENGTH=arg_length, STATUS=command_line_stat)

       IF ((command_line_stat .GT. 0) .OR. (arg_length .EQ. 0)) THEN

          WRITE(*,*) "Cannot process command line definition of namelist fie "
          WRITE(*,*) "using cable.nml as the namelsit"
          IF (ALLOCATED(CABLE_NAMELIST)) DEALLOCATE(CABLE_NAMELIST)
          ALLOCATE(CHARACTER(9) :: CABLE_NAMELIST)
          CABLE_NAMELIST = "cable.nml"

       END IF

       ALLOCATE(CHARACTER(arg_length) :: CABLE_NAMELIST)

       CALL GET_COMMAND_ARGUMENT(1,VALUE=CABLE_NAMELIST,STATUS=command_line_stat)
       !check it ends in nml
       IF ((SCAN(CABLE_NAMELIST, 'l',back=.TRUE.) .NE. (SCAN(CABLE_NAMELIST, 'm',back=.TRUE.)+1)) .OR. &
            (SCAN(CABLE_NAMELIST, 'l',back=.TRUE.) .NE. (SCAN(CABLE_NAMELIST, 'n',back=.TRUE.)+2)) ) THEN

          WRITE(*,*) 'FIRST ARGUMENT DOES NOT END IN NML'
          WRITE(*,*) 'ASSUMING filename%met and casafile%cnpipool'
          arg_not_namelist = .TRUE.
       END IF

       IF ((command_line_stat .NE. 0) .OR. (arg_not_namelist)) THEN
          WRITE(*,*) "Using cable.nml as the namelsit"
          IF (ALLOCATED(CABLE_NAMELIST)) DEALLOCATE(CABLE_NAMELIST)
          ALLOCATE(CHARACTER(9) :: CABLE_NAMELIST)
          CABLE_NAMELIST = "cable.nml"
       END IF
    ELSE
       WRITE(*,*) "Cannot process command line definition of namelist fie "
       WRITE(*,*) "using cable.nml as the namelsit"
       IF (ALLOCATED(CABLE_NAMELIST)) DEALLOCATE(CABLE_NAMELIST)
       ALLOCATE(CHARACTER(9) :: CABLE_NAMELIST)
       CABLE_NAMELIST = "cable.nml"

    END IF

    IF (no_arguments) arg_not_namelist = .TRUE.


  END SUBROUTINE get_namelist_file_name



END MODULE  cable_namelist_util
