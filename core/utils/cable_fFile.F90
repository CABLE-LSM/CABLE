MODULE cable_fFile_module
  USE cable_common_module

  IMPLICIT NONE
  INTEGER, PARAMETER :: gok=0
  INTEGER :: galloctest=1

  !CABLE_LSM: intro'd quick writing capbility. remove from here. keep for ref
  CHARACTER(len=200) :: fprintf_dir_root = "/"

  CHARACTER(len=15) :: unique_subdir = "./"

  LOGICAL :: L_cable_fprint   = .FALSE.,    &
       L_cable_Pyfprint = .FALSE.

  CHARACTER(len=300) :: fprintf_dir

CONTAINS

  SUBROUTINE open_file_per_node( iDiag,pDiag, dir, basename, node, fbasename )
    ! use cable_common_module
    INTEGER :: iDiag, pDiag, node
    INTEGER :: gopenstatus = 1
    CHARACTER(len=*)  :: dir
    CHARACTER(len=*)  :: basename
    CHARACTER(len=*), OPTIONAL  :: fbasename
    CHARACTER(len=300) :: infilename
    CHARACTER(len=30) :: chnode

    WRITE(chnode,10) node
10  FORMAT(i3.3)
    infilename=TRIM( TRIM(dir)//TRIM(basename)//TRIM(chnode) )

    IF(iDiag==0) THEN
       pDiag = pDiag+2
       iDiag=pDiag

       CALL open_iDiag( iDiag, infilename, gopenstatus)

    ENDIF

    !jhan:check if file is open
    IF(gopenstatus==gok) THEN
       fbasename = infilename
       RETURN
    ELSE
       WRITE (*,*) infilename,' NOT open for write. Error(open_file_per_node)'
       WRITE (*,*) '???Check that the path exists???'
       STOP
    ENDIF

  END SUBROUTINE open_file_per_node

  SUBROUTINE open_iDiag( iDiag, infilename, gopenstatus)

    !   use cable_common_module
    INTEGER :: iDiag
    INTEGER :: gopenstatus
    CHARACTER(len=*) :: infilename
    CHARACTER(len=300) :: ffilename

    ffilename=TRIM( TRIM(infilename)// '.txt' )

    WRITE(6,*) "jh:Pytest:fname ", ffilename
    OPEN( unit=iDiag, file=ffilename, status="unknown", &
         action="write", iostat=gopenstatus, form="formatted", &
         position='append' )

  END SUBROUTINE open_iDiag

END MODULE cable_fFile_module
