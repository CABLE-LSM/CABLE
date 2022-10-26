MODULE cable_fFile_module
  use cable_common_module
  
  IMPLICIT NONE
  INTEGER, PARAMETER :: gok=0
  INTEGER :: galloctest=1

  !CABLE_LSM: intro'd quick writing capbility. remove from here. keep for ref
  character(len=200) :: fprintf_dir_root = "/"
  
  character(len=15) :: unique_subdir = "./"

  logical :: L_cable_fprint   = .FALSE.,    &  
             L_cable_Pyfprint = .FALSE. 
              
  character(len=300) :: fprintf_dir

CONTAINS

SUBROUTINE open_file_per_node( iDiag,pDiag, dir, basename, node, fbasename )
  ! use cable_common_module
   integer :: iDiag, pDiag, node
   integer :: gopenstatus = 1
   character(len=*)  :: dir 
   character(len=*)  :: basename
   character(len=*), optional  :: fbasename
   character(len=300) :: infilename
   character(len=30) :: chnode

   write(chnode,10) node
10 format(i3.3)   
   infilename=trim( trim(dir)//trim(basename)//trim(chnode) )
    
   IF(iDiag==0) tHEN
      pDiag = pDiag+2  
      iDiag=pDiag
         
      call open_iDiag( iDiag, infilename, gopenstatus) 
      
   ENDIF
   
   !jhan:check if file is open      
   if(gopenstatus==gok) then
      fbasename = infilename
      return
   else
      write (*,*) infilename,' NOT open for write. Error(open_file_per_node)'
      write (*,*) '???Check that the path exists???'
      STOP
   endif

END SUBROUTINE open_file_per_node

subroutine open_iDiag( iDiag, infilename, gopenstatus) 

!   use cable_common_module
   integer :: iDiag 
   integer :: gopenstatus
   character(len=*) :: infilename
   character(len=300) :: ffilename
   
   ffilename=trim( trim(infilename)// '.txt' )

   open( unit=iDiag, file=ffilename, status="unknown", &
     action="write", iostat=gopenstatus, form="formatted", &
     position='append' )

End subroutine open_iDiag

subroutine qprint( iDiag, infilename) 

   integer :: iDiag 
   character(len=*) :: infilename
   character(len=300) :: ffilename
   
   ffilename=trim( trim(infilename)// '.txt' )

   open( unit=iDiag, file=ffilename, status="unknown", &
     action="write", form="formatted", &
     position='append' )

End subroutine qprint 


END MODULE cable_fFile_module


