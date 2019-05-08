MODULE cable_Pyfprint_module
  USE cable_fFile_module

  INTERFACE cable_Pyfprintf
     MODULE PROCEDURE cable_Pyfprintf1, cable_Pyfprintf2
  END INTERFACE

CONTAINS

  ! 1-D REAL
  !==========================================================================!

  SUBROUTINE cable_Pyfprintf1( iDiag, basename, var1, dimx, L_fprint )
    USE cable_fFile_module
    USE cable_common_module, ONLY : knode_gl, ktau_gl
    IMPLICIT NONE
    ! IN vars
    INTEGER :: iDiag  ! f^n creates unique unidID to be returned to calling point
    INTEGER :: dimx   ! 1-D length
    REAL, DIMENSION(dimx) :: var1    ! var CALLed
    ! writes file per processor (basename+node)
    CHARACTER(len=*) :: basename     ! filename based on var
    LOGICAL :: L_fprint
    ! LOCAL vars
    INTEGER, SAVE :: pDiag=1713      ! give unique SEED per module procedure

    IF( .NOT. L_fprint ) RETURN

    ! Returns unique unit=iDiag and modified basename
    CALL open_file_per_node( iDiag, pDiag, fprintf_dir, basename, knode_gl )
    WRITE(6,*) "jh:Pytest:iDiag ",  iDiag
    WRITE(6,*) "jh:Pytest:pDiag ", pdiag
    WRITE(6,*) "jh:Pytest:fprintf_dir ",  fprintf_dir
    WRITE(6,*) "jh:Pytest:basename ",  basename
    WRITE(6,*) "jh:Pytest:knode_gl  ",  knode_gl

    CALL cable_Pyprint1( idiag, dimx, var1, ktau_gl )

  END SUBROUTINE cable_Pyfprintf1


  SUBROUTINE cable_Pyprint1(idiag,dimx,field, ktau)
    IMPLICIT NONE
    INTEGER :: idiag
    INTEGER :: dimx
    REAL, DIMENSION(dimx) :: field
    INTEGER :: ktau

    INTEGER :: j

    WRITE(6,*) "jh:2Pytest:iDiag ",  iDiag
    WRITE(6,*) "jh:2Pytest:dimx ", dimx
    WRITE(6,*) "jh:2Pytest:ktau  ",  ktau
    DO j=1, dimx
       WRITE (iDiag,*) field(j)
    ENDDO

  END SUBROUTINE cable_Pyprint1

  ! 2-D REAL
  !==========================================================================!

  SUBROUTINE cable_Pyfprintf2( iDiag, basename, var1, dimx, dimy, L_fprint )
    USE cable_fFile_module
    USE cable_common_module, ONLY : knode_gl, ktau_gl
    IMPLICIT NONE
    ! IN vars
    INTEGER :: iDiag  ! f^n creates unique unidID to be returned to calling point
    INTEGER :: dimx, dimy   ! 2-D length
    REAL, DIMENSION(dimx,dimy) :: var1    ! var CALLed
    ! writes file per processor (basename+node)
    CHARACTER(len=*) :: basename     ! filename based on var
    LOGICAL :: L_fprint
    ! LOCAL vars
    INTEGER, SAVE :: pDiag=2713      ! give unique SEED per module procedure

    IF( .NOT. L_fprint ) RETURN

    ! Returns unique unit=iDiag and modified basename
    CALL open_file_per_node( iDiag, pDiag, fprintf_dir, basename, knode_gl )
    WRITE(6,*) "jh:Pytest:iDiag ",  iDiag
    WRITE(6,*) "jh:Pytest:pDiag ", pdiag
    WRITE(6,*) "jh:Pytest:fprintf_dir ",  fprintf_dir
    WRITE(6,*) "jh:Pytest:basename ",  basename
    WRITE(6,*) "jh:Pytest:knode_gl  ",  knode_gl

    CALL cable_Pyprint2( iDiag, dimx, dimy, var1, ktau_gl )

  END SUBROUTINE cable_Pyfprintf2


  SUBROUTINE cable_Pyprint2(idiag,dimx,dimy, field, ktau)
    IMPLICIT NONE
    INTEGER :: idiag
    INTEGER :: dimx, dimy
    REAL, DIMENSION(dimx, dimy) :: field
    INTEGER :: ktau

    INTEGER :: i,j

    WRITE(6,*) "jh:2Pytest:iDiag ",  iDiag
    WRITE(6,*) "jh:2Pytest:dimx ", dimx
    WRITE(6,*) "jh:2Pytest:dimy ", dimy
    WRITE(6,*) "jh:2Pytest:ktau  ",  ktau
    DO i=1, dimx
       DO j=1, dimy
          WRITE (iDiag,*) field(i,j)
       ENDDO
    ENDDO

  END SUBROUTINE cable_Pyprint2

END MODULE cable_Pyfprint_module
