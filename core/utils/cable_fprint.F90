MODULE cable_fprint_module

  INTERFACE cable_fprintf
     MODULE PROCEDURE cable_fprintf0, cable_fprintf1, cable_fprintf2
  END INTERFACE

CONTAINS

  ! Writes to std output stereamm(=6)
  !===============================================================================

  SUBROUTINE cable_fprintf0( basename, L_fprint )
    USE cable_common_module, ONLY : knode_gl, ktau_gl
    IMPLICIT NONE

    CHARACTER(len=*) :: basename     ! subr name
    LOGICAL :: L_fprint

    IF( .NOT. L_fprint ) RETURN

    CALL cable_print0( 6, basename, ktau_gl, knode_gl )

  END SUBROUTINE cable_fprintf0

  SUBROUTINE cable_print0(idiag, basename, ktau, mype )
    IMPLICIT NONE
    INTEGER :: idiag
    CHARACTER(len=*) :: basename     ! subr name
    INTEGER :: ktau, mype

19  FORMAT(  "CABLE_LSM: ", A20, " Finished."   )
    WRITE (idiag, 19)  basename

  END SUBROUTINE cable_print0

  !===============================================================================

  ! Writes to file
  !===============================================================================

  SUBROUTINE cable_fprintf1( iDiag, basename, knode, ktau, L_fprint )
    USE cable_fFile_module, ONLY : open_file_per_node, fprintf_dir
    USE cable_common_module, ONLY :
    IMPLICIT NONE

    INTEGER :: iDiag  ! f^n creates unique unidID to be returned to calling point
    INTEGER :: ktau, knode
    ! writes file per processor (basename+node)
    CHARACTER(len=*) :: basename     ! subr name
    LOGICAL :: L_fprint
    ! LOCAL vars
    INTEGER, SAVE :: pDiag=3713      ! give unique SEED per module procedure

    IF( .NOT. L_fprint ) RETURN

    ! Returns unique unit=iDiag and modified basename
    CALL open_file_per_node( iDiag, pDiag, fprintf_dir, basename, knode )

    CALL cable_print1( idiag, basename, ktau, knode )

  END SUBROUTINE cable_fprintf1

  SUBROUTINE cable_print1(idiag, basename, ktau, mype )
    IMPLICIT NONE
    INTEGER :: idiag
    CHARACTER(len=*) :: basename     ! subr name
    INTEGER :: ktau, mype

19  FORMAT(  "CABLE_LSM: ", A20, " finished @", I8.1, " on ",I3.1  )
    WRITE (idiag, 19)  basename, ktau, mype

  END SUBROUTINE cable_print1

  !===============================================================================

  ! Writes to file
  !===============================================================================
  SUBROUTINE cable_fprintf2( iDiag, basename, msg, L_fprint )
    USE cable_fFile_module, ONLY : open_file_per_node, fprintf_dir
    USE cable_common_module, ONLY : knode_gl, ktau_gl
    IMPLICIT NONE

    INTEGER :: iDiag  ! f^n creates unique unidID to be returned to calling point
    ! writes file per processor (basename+node)
    CHARACTER(len=*) :: basename     ! subr name
    CHARACTER(len=*) :: msg ! message
    LOGICAL :: L_fprint
    ! LOCAL vars
    INTEGER, SAVE :: pDiag=5713      ! give unique SEED per module procedure

    IF( .NOT. L_fprint ) RETURN

    ! Returns unique unit=iDiag and modified basename
    CALL open_file_per_node( iDiag, pDiag, fprintf_dir, basename, knode_gl )

    CALL cable_print2( idiag, basename, msg, ktau_gl, knode_gl )

  END SUBROUTINE cable_fprintf2

  SUBROUTINE cable_print2(idiag, basename, msg, ktau, mype )
    IMPLICIT NONE
    INTEGER :: idiag
    CHARACTER(len=*) :: basename     ! subr name
    CHARACTER(len=*) :: msg ! message
    INTEGER :: ktau, mype

    WRITE (idiag, *) ""
    WRITE (idiag, *) msg
    WRITE (idiag, *) &
         "lsm_id set CABLE. call succeeded however return and call JULES"
    WRITE (idiag, *) ""

  END SUBROUTINE cable_print2



END MODULE cable_fprint_module
