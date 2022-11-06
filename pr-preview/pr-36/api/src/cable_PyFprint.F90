MODULE cable_Pyfprint_module
  USE cable_fFile_module

   interface cable_Pyfprintf
      module procedure cable_Pyfprintf1, cable_Pyfprintf2
   end interface cable_Pyfprintf

CONTAINS

! 1-D REAL
!==========================================================================!

SUBROUTINE cable_Pyfprintf1( iDiag, basename, var1, dimx, L_fprint )
  USE cable_fFile_module
  USE cable_common_module, only : knode_gl, ktau_gl
  implicit none  
  ! IN vars
  integer :: iDiag  ! f^n creates unique unidID to be returned to calling point
  integer :: dimx   ! 1-D length
  real, dimension(dimx) :: var1    ! var CALLed    
  ! writes file per processor (basename+node)
  character(len=*) :: basename     ! filename based on var
  logical :: L_fprint
  ! LOCAL vars
  integer, SAVE :: pDiag=1713      ! give unique SEED per module procedure 

  if( .NOT. L_fprint ) return

  ! Returns unique unit=iDiag and modified basename
  call open_file_per_node( iDiag, pDiag, fprintf_dir, basename, knode_gl )
   
  call cable_Pyprint1( idiag, dimx, var1, ktau_gl )
                             
END SUBROUTINE cable_Pyfprintf1


SUBROUTINE cable_Pyprint1(idiag,dimx,field, ktau)
  implicit none  
  integer :: idiag
  integer :: dimx
  real, dimension(dimx) :: field
  integer :: ktau 
  
  integer :: j

  do j=1, dimx  
      write (iDiag,*) field(j) 
  enddo

END SUBROUTINE cable_Pyprint1

! 2-D REAL
!==========================================================================!

SUBROUTINE cable_Pyfprintf2( iDiag, basename, var1, dimx, dimy, L_fprint )
  USE cable_fFile_module
  USE cable_common_module, only : knode_gl, ktau_gl
  implicit none  
  ! IN vars
  integer :: iDiag  ! f^n creates unique unidID to be returned to calling point
  integer :: dimx, dimy   ! 2-D length
  real, dimension(dimx,dimy) :: var1    ! var CALLed    
  ! writes file per processor (basename+node)
  character(len=*) :: basename     ! filename based on var
  logical :: L_fprint
  ! LOCAL vars
  integer, SAVE :: pDiag=2713      ! give unique SEED per module procedure 

  if( .NOT. L_fprint ) return

  ! Returns unique unit=iDiag and modified basename
  call open_file_per_node( iDiag, pDiag, fprintf_dir, basename, knode_gl )
   
  call cable_Pyprint2( iDiag, dimx, dimy, var1, ktau_gl )
                             
END SUBROUTINE cable_Pyfprintf2


SUBROUTINE cable_Pyprint2(idiag,dimx,dimy, field, ktau)
  implicit none  
  integer :: idiag
  integer :: dimx, dimy
  real, dimension(dimx, dimy) :: field
  integer :: ktau 
  
  integer :: i,j

  do i=1, dimx  
  do j=1, dimy  
      write (iDiag,*) field(i,j) 
  enddo
  enddo

END SUBROUTINE cable_Pyprint2

END MODULE cable_Pyfprint_module


