MODULE cable_fprint_module

  interface cable_fprintf
    module procedure cable_fprintf0, cable_fprintf1, cable_fprintf2
  end interface cable_fprintf

CONTAINS

! Writes to std output stereamm(=6) 
!===============================================================================

SUBROUTINE cable_fprintf0( basename, L_fprint )
  USE cable_common_module, only : knode_gl, ktau_gl
  implicit none  
  
  character(len=*) :: basename     ! subr name 
  logical :: L_fprint

  if( .NOT. L_fprint ) return
   
  call cable_print0( 6, basename, ktau_gl, knode_gl )
                             
END SUBROUTINE cable_fprintf0

SUBROUTINE cable_print0(idiag, basename, ktau, mype )
  implicit none  
  integer :: idiag
  character(len=*) :: basename     ! subr name 
  integer :: ktau, mype 

19 format(  "CABLE_LSM: ", A20, " Finished."   )
  write (idiag, 19)  basename

END SUBROUTINE cable_print0

!===============================================================================

! Writes to file
!===============================================================================

SUBROUTINE cable_fprintf1( iDiag, basename, knode, ktau, L_fprint )
  USE cable_fFile_module, ONLY : open_file_per_node, fprintf_dir
  USE cable_common_module, only : 
  implicit none  

  integer :: iDiag  ! f^n creates unique unidID to be returned to calling point
  integer :: ktau, knode
  ! writes file per processor (basename+node)
  character(len=*) :: basename     ! subr name 
  logical :: L_fprint
  ! LOCAL vars
  integer, SAVE :: pDiag=3713      ! give unique SEED per module procedure 

  if( .NOT. L_fprint ) return

  ! Returns unique unit=iDiag and modified basename
  call open_file_per_node( iDiag, pDiag, fprintf_dir, basename, knode )
   
  call cable_print1( idiag, basename, ktau, knode )
                             
END SUBROUTINE cable_fprintf1

SUBROUTINE cable_print1(idiag, basename, ktau, mype )
  implicit none  
  integer :: idiag
  character(len=*) :: basename     ! subr name 
  integer :: ktau, mype 

19 format(  "CABLE_LSM: ", A20, " finished @", I8.1, " on ",I3.1  )
  write (idiag, 19)  basename, ktau, mype 

END SUBROUTINE cable_print1

!===============================================================================

! Writes to file
!===============================================================================
SUBROUTINE cable_fprintf2( iDiag, basename, msg, L_fprint )
  USE cable_fFile_module, ONLY : open_file_per_node, fprintf_dir
  USE cable_common_module, only : knode_gl, ktau_gl
  implicit none  
  
  integer :: iDiag  ! f^n creates unique unidID to be returned to calling point
  ! writes file per processor (basename+node)
  character(len=*) :: basename     ! subr name 
  character(len=*) :: msg ! message 
  logical :: L_fprint
  ! LOCAL vars
  integer, SAVE :: pDiag=5713      ! give unique SEED per module procedure 

  if( .NOT. L_fprint ) return

  ! Returns unique unit=iDiag and modified basename
  call open_file_per_node( iDiag, pDiag, fprintf_dir, basename, knode_gl )
   
  call cable_print2( idiag, basename, msg, ktau_gl, knode_gl )
                             
END SUBROUTINE cable_fprintf2

SUBROUTINE cable_print2(idiag, basename, msg, ktau, mype )
  implicit none  
  integer :: idiag
  character(len=*) :: basename     ! subr name 
  character(len=*) :: msg ! message 
  integer :: ktau, mype 

  write (idiag, *) ""  
  write (idiag, *) msg  
  write (idiag, *) &
    "lsm_id set CABLE. call succeeded however return and call JULES" 
  write (idiag, *) ""  

END SUBROUTINE cable_print2



END MODULE cable_fprint_module


