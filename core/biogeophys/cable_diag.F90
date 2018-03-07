!==============================================================================
! This source code is part of the
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CSIRO Open Source Software License
! Agreement (variation of the BSD / MIT License).
!
! You may not use this file except in compliance with this License.
! A copy of the License (CSIRO_BSD_MIT_License_v2.0_CABLE.txt) is located
! in each directory containing CABLE code.
!
! ==============================================================================
! Purpose: handles additional, dynamically decided diagnostic output from model.
!          permanently used for bitwise identical testing. more applications
!          will follow.
!
! Contact: Jhan.Srbinovsky@csiro.au
!
! History: Currently stripped down version of cable_diag here. will be
!          re-implemented in time.
!
! ==============================================================================

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++ USE this module in any subr. you wish to write vars from.             +++!
!+++ x is typically the number of landpoints(tiles). binary file is        +++!
!+++ then appended every timestep with the new foo(x_i)                    +++!
!+++                                                                       +++!
!+++ CALL syntax:                                                          +++!
!+++                                                                       +++!
!+++ cable_diag( Nvars, filename, dimx, dimy, timestep, vname1, var1 )     +++!
!+++                                                                       +++!
!+++ output binaries can be interpreted from the command line              +++!
!+++ using a suite of tools. Currently, only zero_diff.ksh is supported.   +++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!CABLE_LSM:This has to be commented for offline
!#define UM_BUILD YES
MODULE cable_diag_module
  use cable_def_types_mod, only : r_2
   IMPLICIT NONE
   INTEGER, PARAMETER :: gok=0
   INTEGER :: galloctest=1

   !--- subrs overloaded to respond to call cable_diag
   INTERFACE cable_diag
      MODULE PROCEDURE cable_diag1
   END INTERFACE cable_diag
!CABLE_LSM:"A" version of diagnostics along the way. farray builds arrays
   interface cable_farray 
      module procedure cable_farray1, cable_farray2
   end interface cable_farray 
   
   interface cable_NaN
      module procedure cable_NaN1, cable_NaN2
   end interface cable_NaN

#ifndef UM_BUILD
  interface put_var_nc
     module procedure put_var_ncr1, put_var_ncr2, put_var_ncr3
  end interface put_var_nc

  interface get_var_nc
     module procedure get_var_ncr2, get_var_ncr3
  end interface get_var_nc

#endif
!CABLE_LSM: procedures for writing out vars in seperate text files
   interface cable_fprintf
      module procedure cable_fprintf1, cable_fprintf2, cable_fprintf3,  & 
                       cable_iprintf1, cable_iprintf2, cable_Lprintf2
   end interface cable_fprintf

!CABLE_LSM: farray builds array for _diags. keep for ref.
   integer, parameter ::                     & 
      farray_nmax = 50
    
   character(len=30), dimension(farray_nmax) :: &
      farray_names  
       
   real, dimension(:,:), allocatable :: & 
      farray_fields

   real, dimension(:,:,:), allocatable :: & 
      farray_fields2

!CABLE_LSM: make avail for diag. check
   integer :: frow_length, frows, fland_pts, fntiles
   real, dimension(:), allocatable :: fFland 
   real, dimension(:,:), allocatable :: ftile_frac

CONTAINS

!==========================================================================!
!CABLE_LSM: make args avail for diag. check across CABLE
SUBROUTINE cable_fsend( row_length, rows, land_pts, ntiles,Fland, tile_frac )
integer :: row_length, rows, land_pts, ntiles
real, dimension(land_pts) :: Fland 
real, dimension(land_pts, ntiles) :: tile_frac

   if(.NOT. allocated(fFland) ) allocate( fFland(land_pts) )
   if(.NOT. allocated(ftile_frac) ) allocate( ftile_frac(land_pts, ntiles) )
   
   frow_length = row_length 
   frows = rows 
   fland_pts = land_pts 
   fntiles = ntiles
   fFland =  Fland 
   ftile_frac =  tile_frac

END SUBROUTINE cable_fsend

!CABLE_LSM: procedures for writing out vars in seperate text files
SUBROUTINE fprint_fland( iDiag, dir, basename, node )

   use cable_common_module
   ! IN vars
   integer :: iDiag  ! f^n creates unique unidID to be returned to calling point
   
   ! writes file per processor (basename+node)
   integer :: node                  ! processor number 
   character(len=*) :: dir          ! dir name based on var
   character(len=*) :: basename     ! filename based on var
   character(len=29) :: fbasename     ! filename based on var

   ! LOCAL vars
   real :: fsum
   integer, SAVE :: pDiag=713       ! give unique SEED per module procedure 
   integer :: j, jctr, k            ! local counters

   fbasename = basename 
   ! Returns unique unit=iDiag and modified basename
   call open_file_per_node( iDiag, pDiag, dir, basename, node, fbasename=fbasename )
   
  if (node==0) then      
    write(iDiag,*) "This meesage is only written for node==0 " 
    write(iDiag,*) "Checking over identified land points that tile fraction "
    write(iDiag,*) "sums to 1 (within tolerance). And is zero for non-land"
    write(iDiag,*) "Write(s) follow(s) ONLY for points where these "
    write(iDiag,*) "conditions NOT met"
    write(iDiag,*) "Jhan:Check that Total land points = 0 is OK"
    write(iDiag,*) "" 
  endif  
  
  jctr=0
  do j=1, fland_pts
    if( fFland(j)> 0.) then
      jctr = jctr + 1
      fsum = sum(ftile_frac(j,:))
      if(fsum < 0.9999 .OR. fsum > 1.001 ) then
         write(iDiag,*) "Summed tile_frac is: ", fsum
         write(iDiag,*) "    for this land_pt ", j
         write(iDiag,*) "" 
      endif   
    endif   
    if( sum(ftile_frac(j,:)) > 0. .AND. fFland(j) == 0.) then
      write(iDiag,*) "Report: Summed tile_frac is: ", sum( ftile_frac(j,:) )
      write(iDiag,*) "    BUT Fland is zero " 
      write(iDiag,*) "" 
    endif
  enddo
  
  !write(iDiag,*) "Total land points ", jctr
   
  close(iDiag)
  
  call remove_empty_file( iDiag, fbasename )
   
                           
END SUBROUTINE fprint_fland

!==========================================================================!
! writes text files Interfaced Module Procedures 
!==========================================================================!

! 1-D REAL
SUBROUTINE cable_fprintf1( iDiag, basename, var1, dimx, L_fprint )
  use cable_common_module
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
   
  call cable_fextremes1( idiag, dimx, var1, ktau_gl )
                             
END SUBROUTINE cable_fprintf1


!SUBROUTINE cable_fextremes1(fname,field,mype)
SUBROUTINE cable_fextremes1(idiag,dimx,field, ktau)

   integer :: idiag
   integer :: dimx
   real, dimension(dimx) :: field
   integer :: ktau 
   
   integer, parameter :: nbins =3 
   real :: emax, emin, emean, emode
   real :: erange 
   real, dimension(nbins) :: bin
   
   integer :: i,j,k
   integer :: n,m,op
   real :: edbin 
   integer, dimension(nbins) :: ibin
   integer :: ib, ibmax, binmax, maxbin
   !logical,save :: first_call=.true.
      
  !T!if( ktau==1 .OR. mod( ktau,10)==0 ) then  
    write (iDiag,*) ""
    write (iDiag,*) "timestep ", ktau 
  !T!endif
   
  if( dimx > 0) then
    emax =  MAXVAL( field )
    emin =  MINVAL(field )
    emean =  SUM(field(:) ) / ( dimx )
    if( ktau==1 ) &   
      write (iDiag,*) "dimx ", dimx
  else
    emax = 123.123 
    emin = 123.123 
    emean = 123.123 
    if( ktau==1 ) &   
      write (iDiag,*) "dimx is zero here"
  endif 
    
  erange = emax - emin
  edbin = erange / nbins 

  bin(1) = emin          

  ! define bins per field
  do ib=2, nbins
    bin(ib) = bin(ib-1) + edbin
  enddo

  ibin =0
  ! for each Element in field 
  do j=1, dimx  
      
    ! Assignn each Element to a bin 
    do ib=1,(nbins-1) 

      IF( field(j) >= bin(ib) .AND. &
          field(j) < bin(ib+1) ) THEN
               
         !if(ib==1 )print *, "jhan:field1 ", field(i,j)
         ibin(ib) = ibin(ib) + 1
      ENDIF   
          
    enddo ! DO LOOP over fill bins

  enddo ! DO LOOP over elements 

!jhan:reform bin(ib) to value at middle of bin
17 format(  "bin(", I2.1, ") ", E15.6, 4X, "ibin(", I2.1, ") ", I6.1 )
 
  binmax = 0 
  maxbin = 1  
  ! find max bin per field 
  do ib=1, (nbins-1)

    !C!if( ktau==1 ) &   
      !C!write (iDiag,17) ib, bin(ib), ib, ibin(ib) 

    ! ibin(ib) = the # of values in each bin 
    IF( ibin(ib) > binmax ) THEN 
      binmax = ibin(ib) 
      maxbin = ib ! this is actually the starting index of the bin
    ENDIF   
       
  enddo ! DO LOOP over bins
     
  !if( ktau==1 ) &   
  !  write (iDiag,*) "ib: maxbin ", maxbin

  emode = bin(maxbin) 

18 format(  " Min.", 18X, " Max.", 18X, " Mean " )
19 format(  E15.6, 6X, E15.6, 6X, E15.6 )
  if( ktau==1 .OR. mod( ktau,10)==0 ) &   
    write (iDiag,18)
  write (iDiag,19) emin, emax, emean

!18 format(  " Min.", 18X, " Max.", 18X, " Mean ", 18X, " Mode" )
!  write (iDiag,19) emin, emax, emean, emode
!19 format(  E15.6, 6X, E15.6, 6X, E15.6, 6X, E15.6 )


!19 format(  "     ", I6.1, "     ", E15.6, "     ", E15.6, "     ", E15.6 )
!19 format(  "ktau: ", I6.1, " Min.", E15.6, " Max.", E15.6, " Med.", E15.6 )

END SUBROUTINE cable_fextremes1




! 2-D REAL
SUBROUTINE cable_fprintf2( iDiag, basename, dimx, dimy, timestep, node, &
                        dir, var1 )
  use cable_common_module
  ! IN vars
  integer :: iDiag  ! f^n creates unique unidID to be returned to calling point
  integer :: dimx   ! 1-D length
  integer :: dimy   ! 2-D length
  integer :: timestep  
  real, dimension(dimx,dimy) :: var1    ! var CALLed    
  character(len=*) :: dir ! obsolete here - remove 
  
  ! writes file per processor (basename+node)
  integer :: node                  ! processor number 
  character(len=*) :: basename     ! filename based on var

  ! LOCAL vars
  integer, SAVE :: pDiag=2713      ! give unique SEED per module procedure 
  integer :: j, k                  ! local counters
  integer :: r, rl                 ! local counters
  logical, save :: first_call = .true.

  ! Returns unique unit=iDiag and modified basename
  call open_file_per_node( iDiag, pDiag, dir, basename, node )
        
  ! Writes when a new timestep   
  CALL check_timestep( iDiag, timestep )

  if( first_call )  write (iDiag,*) "x, y, z" 
  if( dimx==fland_pts .AND. dimy==fntiles ) then
    if( first_call )  & 
      write (iDiag,*) "Filtered: tile fraction > 0 " 
  elseif( dimx==frow_length .AND. dimy==frows ) then 
    if( first_call )  & 
      write (iDiag,*) "Filtered: NOT as yet " 
  else   
    if( first_call )  & 
      write (iDiag,*) "dimensions not recognized ", dimx, dimy 
    return
  endif

  first_call = .false.

17 format (2I3.1,ES15.6)
  do j=1,dimx      
  do k=1,dimy      
    if( dimx==fland_pts .AND. dimy==fntiles ) then 
      if( ftile_frac(j,k) > 0. ) write (iDiag,17) j,k, var1(j,k)
    elseif( dimx==frow_length .AND. dimy==frows ) then 
      !if( ftile_frac(j,k) > 0. ) 
      write (iDiag,17) j,k, var1(j,k)
    else   
      write (iDiag,*) "Should never get here"
      return
    endif   
  enddo   
  enddo   
                            
END SUBROUTINE cable_fprintf2
 
! 3-D REAL
SUBROUTINE cable_fprintf3( iDiag, basename, dimx, dimy, dimz, timestep, node, &
                        dir, var1 )
   use cable_common_module
   ! IN vars
   integer :: iDiag  ! f^n creates unique unidID to be returned to calling point
   integer :: dimx   ! 1-D length
   integer :: dimy   ! 2-D length
   integer :: dimz   ! 3-D length
   integer :: timestep  
   real, dimension(dimx,dimy,dimz) :: var1    ! var CALLed    
   character(len=*) :: dir ! obsolete here - remove 
   
   ! writes file per processor (basename+node)
   integer :: node                  ! processor number 
   character(len=*) :: basename     ! filename based on var

   ! LOCAL vars
   integer, SAVE :: pDiag=3713      ! give unique SEED per module procedure 
   integer :: j, k, m               ! local counters 
   logical, save :: first_call = .true.

   ! Returns unique unit=iDiag and modified basename
   call open_file_per_node( iDiag,pDiag, dir, basename, node )
         
   ! Writes when a new timestep   
   CALL check_timestep( iDiag, timestep )

   if( first_call )  then
      write (iDiag,*) "x, y, z" 
      if( dimx==fland_pts .AND. dimy==fntiles ) then
         write (iDiag,*) "Filtered: tile fraction > 0 " 
      elseif( dimx==frow_length .AND. dimy==frows ) then 
         write (iDiag,*) "Filtered: NOT as yet " 
      else   
         write (iDiag,*) "dimensions not recognized ", dimx, dimy 
      endif
   endif
   first_call = .false.
         
   do j=1,dimx      
   do k=1,dimy      
   do m=1,dimz      
      if( dimx==fland_pts .AND. dimy==fntiles ) then 
         if( ftile_frac(j,k) > 0. ) then 
            write (iDiag,*) j,k,m
            write (iDiag,*) var1(j,k,m)
         endif   
      elseif( dimx==frow_length .AND. dimy==frows ) then 
         write (iDiag,*) j,k,m
         write (iDiag,*) var1(j,k,m)
      else   
         write (iDiag,*) j,k,m
         write (iDiag,*) var1(j,k,m)
      endif   
   enddo   
   enddo   
   enddo   
                           
END SUBROUTINE cable_fprintf3


! 1-D INTEGER
SUBROUTINE cable_iprintf1( iDiag, basename, dimx, timestep, node, &
                        dir, var1 )
   use cable_common_module
   ! IN vars
   integer :: iDiag  ! f^n creates unique unidID to be returned to calling point
   integer :: dimx   ! 1-D length
   integer :: timestep  
   integer, dimension(dimx) :: var1    ! var CALLed    
   character(len=*) :: dir ! obsolete here - remove 
   
   ! writes file per processor (basename+node)
   integer :: node                  ! processor number 
   character(len=*) :: basename     ! filename based on var

   ! LOCAL vars
   real, dimension(dimx) :: invar1  ! store var read in 

   integer, SAVE :: pDiag=11713      ! give unique SEED per module procedure 
   
   integer :: j, inj            ! local counters
   logical, save :: first_call = .true.

   ! Returns unique unit=iDiag and modified basename
   call open_file_per_node( iDiag,pDiag, dir, basename, node )
         
   ! Writes when a new timestep   
   CALL check_timestep( iDiag, timestep )
         
   if( first_call )  then
   write (iDiag,*) "x, y, z" 
   if( dimx==fland_pts ) then
      write (iDiag,*) "Filtered: land fraction > 0 " 
   else   
      write (iDiag,*) "dimensions not recognized ", dimx
   endif
   endif
   first_call = .false.
    
   do j=1,dimx      
      if( dimx==fland_pts ) then 
         if( fFland(j) > 0. ) then 
            write (iDiag,*) j
            write (iDiag,*) var1(j)
         endif   
      else   
         write (iDiag,*) j
         write (iDiag,*) var1(j)
      endif   
   enddo   
                             
END SUBROUTINE cable_iprintf1


! 2-D INTEGER
SUBROUTINE cable_iprintf2( iDiag, basename, dimx, dimy, timestep, node, &
                        dir, var1 )
   use cable_common_module
   ! IN vars
   integer :: iDiag  ! f^n creates unique unidID to be returned to calling point
   integer :: dimx   ! 1-D length
   integer :: dimy   ! 2-D length
   integer :: timestep  
   integer, dimension(dimx,dimy) :: var1    ! var CALLed    
   character(len=*) :: dir ! obsolete here - remove 
   
   ! writes file per processor (basename+node)
   integer :: node                  ! processor number 
   character(len=*) :: basename     ! filename based on var

   ! LOCAL vars
   integer, SAVE :: pDiag=21713      ! give unique SEED per module procedure 
   integer :: j, k                  ! local counters
   logical, save :: first_call = .true.

   ! Returns unique unit=iDiag and modified basename
   call open_file_per_node( iDiag, pDiag, dir, basename, node )
         
   ! Writes when a new timestep   
   CALL check_timestep( iDiag, timestep )
         
   if( first_call )  then
   write (iDiag,*) "x, y, z" 
   if( dimx==fland_pts .AND. dimy==fntiles ) then
      write (iDiag,*) "Filtered: tile fraction > 0 " 
   elseif( dimx==frow_length .AND. dimy==frows ) then 
      write (iDiag,*) "Filtered: NOT as yet " 
   else   
      write (iDiag,*) "dimensions not recognized ", dimx, dimy 
   endif
   endif
   first_call = .false.
    
   do j=1,dimx      
   do k=1,dimy      
      if( dimx==fland_pts .AND. dimy==fntiles ) then 
         if( ftile_frac(j,k) > 0. ) then 
            write (iDiag,*) j,k
            write (iDiag,*) var1(j,k)
         endif   
      elseif( dimx==frow_length .AND. dimy==frows ) then 
         write (iDiag,*) j,k
         write (iDiag,*) var1(j,k)
      else   
         write (iDiag,*) j,k
         write (iDiag,*) var1(j,k)
      endif   
   enddo   
   enddo   
                            
END SUBROUTINE cable_iprintf2


! 2-D LOGICAL
SUBROUTINE cable_Lprintf2( iDiag, basename, dimx, dimy, timestep, node, &
                        dir, var1 )
   use cable_common_module
   ! IN vars
   integer :: iDiag  ! f^n creates unique unidID to be returned to calling point
   integer :: dimx   ! 1-D length
   integer :: dimy   ! 2-D length
   integer :: timestep  
   logical, dimension(dimx,dimy) :: var1    ! var CALLed    
   character(len=*) :: dir ! obsolete here - remove 
   
   ! writes file per processor (basename+node)
   integer :: node                  ! processor number 
   character(len=*) :: basename     ! filename based on var

   ! LOCAL vars
   integer, SAVE :: pDiag=41713      ! give unique SEED per module procedure 
   integer :: j, k                  ! local counters
   logical, save :: first_call = .true.

   ! Returns unique unit=iDiag and modified basename
   call open_file_per_node( iDiag, pDiag, dir, basename, node )
         
   ! Writes when a new timestep   
   CALL check_timestep( iDiag, timestep )
         
   if( first_call )  then
   write (iDiag,*) "x, y, z" 
   if( dimx==fland_pts .AND. dimy==fntiles ) then
      write (iDiag,*) "Filtered: tile fraction > 0 " 
   elseif( dimx==frow_length .AND. dimy==frows ) then 
      write (iDiag,*) "Filtered: NOT as yet " 
   else   
      write (iDiag,*) "dimensions not recognized ", dimx, dimy 
   endif
   endif
   first_call = .false.
    
   do j=1,dimx      
   do k=1,dimy      
      if( dimx==fland_pts .AND. dimy==fntiles ) then 
         if( ftile_frac(j,k) > 0. ) then 
            write (iDiag,*) j,k
            write (iDiag,*) var1(j,k)
         endif   
      elseif( dimx==frow_length .AND. dimy==frows ) then 
         write (iDiag,*) j,k
         write (iDiag,*) var1(j,k)
      else   
         write (iDiag,*) j,k
         write (iDiag,*) var1(j,k)
      endif   
   enddo   
   enddo   
                            
END SUBROUTINE cable_Lprintf2

!==========================================================================!
!==========================================================================!

SUBROUTINE open_file_per_node( iDiag,pDiag, dir, basename, node, fbasename )
   use cable_common_module
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
      STOP
   endif

END SUBROUTINE open_file_per_node

subroutine open_iDiag( iDiag, infilename, gopenstatus) 

   use cable_common_module
   integer :: iDiag 
   integer :: gopenstatus
   character(len=*) :: infilename
   character(len=300) :: ffilename
   
   ffilename=trim( trim(infilename)// '.txt' )
   !if( cable_user%check_write ) then 
      open( unit=iDiag, file=ffilename, status="replace", &
         action="write", iostat=gopenstatus, form="formatted", &
         position='append' )
   !endif 
   !if( cable_user%check_read ) then 
   !   open(unit=iDiag,file=trim(infilename)//'.txt',status="old", &
   !      action="read", iostat=gopenstatus )
   !endif 

End subroutine open_iDiag

subroutine remove_empty_file( iDiag, fbasename )

   use cable_common_module
   integer :: iDiag 
   character(len=*) :: fbasename
   integer :: gopenstatus, gios
   character(len=99) :: blnk 
   character(len=49) :: cmd 

  open(unit=iDiag,status="unknown", file=fbasename, &
    action="read", iostat=gopenstatus, form="formatted", position='append' )
    read(idiag,*, iostat=gios) blnk 
  close(idiag)
  if( gios <0 ) then
    call unlink(fbasename)       
  endif 

return

End subroutine remove_empty_file

SUBROUTINE check_timestep( iDiag, timestep )
   ! IN vars
   integer :: iDiag  ! f^n creates unique unidID to be returned to calling point
   integer :: timestep  
   
   integer, save ::otimestep = 0  
!otimestep = 0 
!write (iDiag,*) "otimestep ", otimestep 
!write (iDiag,*) "timestep ", timestep 
   if(otimestep .NE. timestep ) then
      otimestep=timestep
      !if( timestep==1 .OR. mod( timestep,10)==0 ) then
      !  write (iDiag,*) ""
      !  write (iDiag,*) "timestep ", timestep 
      !endif           
   else 
      return
   endif           
   
   return

END SUBROUTINE check_timestep

!SUBROUTINE cable_fprintf2( iDiag, basename, dimx, dimy, timestep, node, &
!                        dir, var1 )
!  use cable_common_module
!  ! IN vars
!  integer :: iDiag  ! f^n creates unique unidID to be returned to calling point
!  integer :: dimx   ! 1-D length
!  integer :: dimy   ! 2-D length
!  integer :: timestep  
!  real, dimension(dimx,dimy) :: var1    ! var CALLed    
!  character(len=*) :: dir ! obsolete here - remove 
!  
!  ! writes file per processor (basename+node)
!  integer :: node                  ! processor number 
!  character(len=*) :: basename     ! filename based on var
!
!  ! LOCAL vars
!  integer, SAVE :: pDiag=2713      ! give unique SEED per module procedure 
!  integer :: j, k                  ! local counters
!  integer :: r, rl                 ! local counters
!  logical, save :: first_call = .true.
!
!  ! Returns unique unit=iDiag and modified basename
!  call open_file_per_node( iDiag, pDiag, dir, basename, node )
!        
!  ! Writes when a new timestep   
!  CALL check_timestep( iDiag, timestep )
!
!  if( first_call )  write (iDiag,*) "x, y, z" 
!  if( dimx==fland_pts .AND. dimy==fntiles ) then
!    if( first_call )  & 
!      write (iDiag,*) "Filtered: tile fraction > 0 " 
!  elseif( dimx==frow_length .AND. dimy==frows ) then 
!    if( first_call )  & 
!      write (iDiag,*) "Filtered: NOT as yet " 
!  else   
!    if( first_call )  & 
!      write (iDiag,*) "dimensions not recognized ", dimx, dimy 
!    return
!  endif
!
!  first_call = .false.
!
!17 format (2I3.1,ES15.6)
!  do j=1,dimx      
!  do k=1,dimy      
!    if( dimx==fland_pts .AND. dimy==fntiles ) then 
!      if( ftile_frac(j,k) > 0. ) write (iDiag,17) j,k, var1(j,k)
!    elseif( dimx==frow_length .AND. dimy==frows ) then 
!      !if( ftile_frac(j,k) > 0. ) 
!      write (iDiag,17) j,k, var1(j,k)
!    else   
!      write (iDiag,*) "Should never get here"
!      return
!    endif   
!  enddo   
!  enddo   
!                            
!END SUBROUTINE cable_fprintf2
!
!SUBROUTINE cable_range2( cDiag30, vname, row_length, rows, ktau_gl, knode_gl, &
!                         dir, tl_1)
!
!   real, dimension(:,:) :: field
!   character(len=*), dimension(:) :: fname
!   
!   integer, optional :: mype 
!   integer :: i,j,k
!   integer :: n,m,op
!   real :: emax, emin, emean, emode
!   real :: erange 
!   real :: edbin 
!   real, dimension(100) :: bin
!   integer, dimension(100) :: ibin
!   integer :: ib, ibmax, binmax, maxbin
!
!   n = size(fname)
!   m = size(field,2)
!
!   ! for each field in fname(i)
!   do i=1, n  
!      
!      emax =  MAXVAL( field(i,:) )
!      emin =  MINVAL(field(i,:) )
!      emean =  SUM(field(i,:) ) / ( m )
!    
!      erange = emax - emin
!      edbin = erange / 100. ! for 100 bins
!
!      bin(1) = emin          
!
!      ! define bins per fname(i)
!      do ib=2, 100
!         bin(ib) = bin(ib-1) + edbin
!      enddo
!
!      ibin =0
!      ! for each Element in field 
!      do j=1, m  
!      
!         ! Assignn each Element to a bin 
!         do ib=1, 99
!
!            IF( field(i,j) >= bin(ib) .AND. &
!                field(i,j) < bin(ib+1) ) THEN
!               
!               !if(ib==1 )print *, "jhan:field1 ", field(i,j)
!               ibin(ib) = ibin(ib) + 1
!            
!            ENDIF   
!          
!         enddo ! DO LOOP over fill bins
!
!      enddo ! DO LOOP over elements 
! 
!      binmax = 0 
!      
!      ! find max bin per field 
!      do ib=1, 99
!
!         IF( ibin(ib) > binmax ) THEN 
!            binmax = ibin(ib) 
!            maxbin = ib 
!         ENDIF   
!       
!      enddo ! DO LOOP over bins
!     
!     print *, "jhan:bins1 ", bin 
!     print *, "jhan:bins1 count", ibin 
!     
!      
!      emode = bin(maxbin) 
!
!         print *, ""
!         print *, "CABLE_log: "
!         print *, "   Field ", fname(i)
!         print *, "   Min ", emin
!         print *, "   Max ", emax
!         print *, "   Mean ",emean
!         print *, "   Mode ",emode
!         print *, "End CABLE_log: "
!         print *, ""
!      
!      enddo ! DO LOOP over fname(i)
!
!
!END SUBROUTINE cable_range2


!==========================================================================!
!==========================================================================!


 

!

!=============================================================================!
!=============================================================================!


! writes binary files 
!==========================================================================!
! cable_diag1/2/3 call subrs to write filename.dat which contains description
! of data and format etc., and filename.bin containing the data
!==========================================================================!

SUBROUTINE cable_diag1( iDiag, basename, dimx, dimy, timestep, node, &
                        vname1, var1, once )
   integer, intent(inOUT) :: iDiag
   integer, SAVE :: pDiag=713
   integer, intent(in) :: dimx, dimy, timestep,node
   real, intent(in), dimension(:) :: var1
   integer, optional :: once
   integer :: Nvars=1 !this WAS input
   integer :: i=0
   character(len=*), intent(in) :: basename, vname1
   character(len=30) :: filename, chnode

      IF(iDiag==0) tHEN
         pDiag = pDiag+2
         iDiag=pDiag
      ENDIF

      write(chnode,10) node
   10 format(i3.3)
      filename=trim(trim(basename)//trim(chnode))

      if (timestep == 1) &
         call cable_diag_desc1( iDiag, trim(filename), dimx, dimy, vname1 )

      if( present(once) ) then
         if (timestep == 1) &
         ! write data only on first timestep
         call cable_diag_data1( iDiag, trim(filename), dimx, timestep, dimy, &
                             var1 )
      else
         ! write data every timestep
         call cable_diag_data1( iDiag, trim(filename), dimx, timestep, dimy, &
                                var1 )
      endif

END SUBROUTINE cable_diag1

!=============================================================================!
!=============================================================================!

SUBROUTINE cable_diag_desc1( iDiag, filename, dimx, dimy, vname1 )

   integer, intent(in) :: iDiag,dimx,dimy
   integer, PARAMETER :: Nvars=1
   character(len=*), intent(in) :: filename, vname1
   integer, save :: gopenstatus = 1

     open(unit=iDiag,file=filename//'.dat', status="replace", &
          action="write", iostat=gopenstatus )

      if(gopenstatus==gok) then
            write (iDiag,*) 'Number of var(s): '
            write (iDiag,*) Nvars
            write (iDiag,*) 'Name of var(s): '
            write (iDiag,7139) vname1
 7139       format(a)
            write (iDiag,*) 'dimension of var(s) in x: '
            write (iDiag,*) dimx
            write (iDiag,*) 'dimension of var(s) in y: '
            write (iDiag,*) dimy
      else
         write (*,*) filename//'.dat',' Error: unable to write'
      endif

   close(iDiag)

END SUBROUTINE cable_diag_desc1


SUBROUTINE cable_diag_data1( iDiag, filename, dimx, timestep, kend, var1  )

   integer, intent(in) :: iDiag, dimx, timestep, kend
   integer, PARAMETER :: Nvars=1
   real, intent(in), dimension(:) :: var1
   character(len=*), intent(in) :: filename
   integer, save :: gopenstatus = 1

   if (timestep == 1)  then
      open(unit=iDiag+1,file=filename//'.bin',status="unknown", &
           action="write", iostat=gopenstatus, form="unformatted", &
           position='append' )
   endif

   if(gopenstatus==gok) then
         write (iDiag+1) var1
   else
      write (*,*) filename//'.bin',' NOT open for write. Error'
   endif

   if (timestep == kend) &
      close(iDiag+1)

END SUBROUTINE cable_diag_data1

!==========================================================================!
!--- cable generic print status
!==========================================================================!

SUBROUTINE cable_stat( routname)
   use cable_common_module, only : ktau_gl, knode_gl

   character(len=*), intent(in) :: routname
      if(knode_gl==1) & 
         write(6,*) 'CABLE@  ', routname, ktau_gl

END SUBROUTINE cable_stat

!==========================================================================!
!--- cable status NaN
!==========================================================================!


SUBROUTINE cable_NaN1(fname,field,mype)

real, dimension(:,:) :: field
character(len=*), dimension(:) :: fname

integer, optional :: mype 
integer :: i,j
logical :: NoNaN, check
integer :: n,m

   n = size(fname)
   m = size(field,2)
   check = .FALSE.   

   do i=1, n  
      
      NoNaN = .TRUE.
      
      do j=1, m  

         call isnan(  field(i,j), check )

         if( check ) then 
            print *, ""
            print *, "CABLE_log: "
            if( present(mype) ) print *, "proc # ", mype
            print *, "   Element: ",j
            print *, "   of field ", fname(i)
            print *, "   is NaN"
            print *, "End CABLE_log: "
            print *, ""
            NoNaN = .FALSE.
         end if 

      enddo
      
      if(NoNaN) then 
         print *, ""
         print *, "CABLE_log: "
         if( present(mype) ) print *, "proc # ", mype
         print *, '   Field: ', fname(i)
         print *, "   is clear of NaNs"
         print *, "End CABLE_log: "
         print *, ""
      end if
       
   enddo

END SUBROUTINE cable_NaN1


SUBROUTINE cable_NaN2(fname,field,mype)

real, dimension(:,:,:) :: field
character(len=*), dimension(:) :: fname

integer, optional :: mype 
integer :: i,j,k
logical :: NoNaN, check
integer :: n,m,op

   n = size(fname)
   m = size(field,2)
   op = size(field,3)
   check = .FALSE.   

   do i=1, n  
      
      NoNaN = .TRUE.
      
      do j=1, m  
      
         do k=1, op  

            call isnan(  field(i,j,k), check )
   
            if( check ) then 
               print *, ""
               print *, "CABLE_log: "
               if( present(mype) ) print *, "proc # ", mype
               print *, "   Element: ",j,op
               print *, "   of field ", fname(i)
               print *, "   is NaN"
               print *, "End CABLE_log: "
               print *, ""
               NoNaN = .FALSE.
            end if 

         enddo
      
      enddo
      
      if(NoNaN) then 
         print *, ""
         print *, "CABLE_log: "
         if( present(mype) ) print *, "proc # ", mype
         print *, '   Field: ', fname(i)
         print *, "   is clear of NaNs"
         print *, "End CABLE_log: "
         print *, ""
      end if
       
   enddo

END SUBROUTINE cable_NaN2





subroutine isnan(var, check) 
   real :: var 
   logical :: check

   if (var .ne. var) then 
      check = .true. 
   else 
      check = .false. 
   end if 

end subroutine isnan


!logical function isinf(a) 
!real a 
!
!!if ((a*0).ne.0) then 
!!isinf = .true. 
!!else 
!!isinf = .false. 
!!end if 
!!return 
!!end 


SUBROUTINE cable_farray1( mp, CheckNames, CheckFields, &
                  n1,f1, n2,f2, n3,f3, n4,f4, n5,f5, n6,f6, n7,f7, & 
                  n8,f8, n9,f9, n10,f10, n11,f11, n12,f12, n13,f13, &
                  n14,f14, n15,f15, n16,f16, n17,f17, n18,f18, & 
                  n19,f19, n20,f20, n21,f21, n22,f22, n23,f23, & 
                  n24,f24, n25,f25, n26,f26, n27,f27, n28,f28, & 
                  n29,f29, n30,f30, n31,f31, n32,f32, n33,f33, & 
                  n34,f34, n35,f35, n36,f36, n37,f37, n38,f38, & 
                  n39,f39, n40,f40, n41,f41, n42,f42, n43,f43, & 
                  n44,f44, n45,f45, n46,f46, n47,f47, n48,f48, & 
                  n49,f49, n50,f50 & 
               )

   integer :: mp
    
   character(len=*), optional :: & 
                        n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, & 
                        n11, n12, n13, n14, n15, n16, n17, n18, n19, n20, &
                        n21, n22, n23, n24, n25, n26, n27, n28, n29, n30, &
                        n31, n32, n33, n34, n35, n36, n37, n38, n39, n40, &
                        n41, n42, n43, n44, n45, n46, n47, n48, n49, n50  

   real, dimension(:), optional :: & 
                        f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, & 
                        f11, f12, f13, f14, f15, f16, f17, f18, f19, f20,& 
                        f21, f22, f23, f24, f25, f26, f27, f28, f29, f30,& 
                        f31, f32, f33, f34, f35, f36, f37, f38, f39, f40,& 
                        f41, f42, f43, f44, f45, f46, f47, f48, f49, f50
   
   character(len=30), dimension(:), allocatable :: CheckNames 
   real, dimension(:,:), allocatable :: CheckFields

   integer :: i, k

   allocate( farray_fields( farray_nmax, mp) ) 
   i = 0

   if( present (n1) .AND. present(f1) ) then
      CALL fill_farray( n1, f1, i )   
   else
      print *, "CABLE_log: cable_farray missing dummy args"
      return
   endif
   
   if( present (n2) .AND. present(f2) ) CALL fill_farray( n2, f2, i )   
   if( present (n3) .AND. present(f3) ) CALL fill_farray( n3, f3, i )   
   if( present (n4) .AND. present(f4) ) CALL fill_farray( n4, f4, i )   
   if( present (n5) .AND. present(f5) ) CALL fill_farray( n5, f5, i )   
   if( present (n6) .AND. present(f6) ) CALL fill_farray( n6, f6, i )   
   if( present (n7) .AND. present(f7) ) CALL fill_farray( n7, f7, i )   
   if( present (n8) .AND. present(f8) ) CALL fill_farray( n8, f8, i )   
   if( present (n9) .AND. present(f9) ) CALL fill_farray( n9, f9, i )   
   if( present (n10) .AND. present(f10) ) CALL fill_farray( n10, f10, i )   
   if( present (n11) .AND. present(f11) ) CALL fill_farray( n11, f11, i )   
   if( present (n12) .AND. present(f12) ) CALL fill_farray( n12, f12, i )   
   if( present (n13) .AND. present(f13) ) CALL fill_farray( n13, f13, i )   
   if( present (n14) .AND. present(f14) ) CALL fill_farray( n14, f14, i )   
   if( present (n15) .AND. present(f15) ) CALL fill_farray( n15, f15, i )   
   if( present (n16) .AND. present(f16) ) CALL fill_farray( n16, f16, i )   
   if( present (n17) .AND. present(f17) ) CALL fill_farray( n17, f17, i )   
   if( present (n18) .AND. present(f18) ) CALL fill_farray( n18, f18, i )   
   if( present (n19) .AND. present(f19) ) CALL fill_farray( n19, f19, i )   
   if( present (n20) .AND. present(f20) ) CALL fill_farray( n20, f20, i )   
   if( present (n21) .AND. present(f21) ) CALL fill_farray( n21, f21, i )   
   if( present (n22) .AND. present(f22) ) CALL fill_farray( n22, f22, i )   
   if( present (n23) .AND. present(f23) ) CALL fill_farray( n23, f23, i )   
   if( present (n24) .AND. present(f24) ) CALL fill_farray( n24, f24, i )   
   if( present (n25) .AND. present(f25) ) CALL fill_farray( n25, f25, i )   
   if( present (n26) .AND. present(f26) ) CALL fill_farray( n26, f26, i )   
   if( present (n27) .AND. present(f27) ) CALL fill_farray( n27, f27, i )   
   if( present (n29) .AND. present(f29) ) CALL fill_farray( n29, f29, i )   
   if( present (n30) .AND. present(f30) ) CALL fill_farray( n30, f30, i )   
   if( present (n31) .AND. present(f31) ) CALL fill_farray( n31, f31, i )   
   if( present (n32) .AND. present(f32) ) CALL fill_farray( n32, f32, i )   
   if( present (n33) .AND. present(f33) ) CALL fill_farray( n33, f33, i )   
   if( present (n34) .AND. present(f34) ) CALL fill_farray( n34, f34, i )   
   if( present (n35) .AND. present(f35) ) CALL fill_farray( n35, f35, i )   
   if( present (n36) .AND. present(f36) ) CALL fill_farray( n36, f36, i )   
   if( present (n37) .AND. present(f37) ) CALL fill_farray( n37, f37, i )   
   if( present (n38) .AND. present(f38) ) CALL fill_farray( n38, f38, i )   
   if( present (n39) .AND. present(f39) ) CALL fill_farray( n39, f39, i )   
   if( present (n40) .AND. present(f40) ) CALL fill_farray( n40, f40, i )   
   if( present (n41) .AND. present(f41) ) CALL fill_farray( n41, f41, i )   
   if( present (n42) .AND. present(f42) ) CALL fill_farray( n42, f42, i )   
   if( present (n50) .AND. present(f50) ) CALL fill_farray( n50, f50, i )   
   
   allocate( CheckNames(i) )
   allocate( CheckFields(i,mp) )

   do k=1,i
      CheckNames(k) = farray_names(k) 
      CheckFields(k,:) = farray_fields(k,:)
   enddo        

   deallocate( farray_fields )
    
END SUBROUTINE cable_farray1



SUBROUTINE cable_farray2( mp, np, CheckNames, CheckFields, &
                  n1,f1, n2,f2, n3,f3, n4,f4, n5,f5, n6,f6, n7,f7, & 
                  n8,f8, n9,f9, n10,f10, n11,f11, n12,f12, n13,f13, &
                  n14,f14, n15,f15, n16,f16, n17,f17, n18,f18, & 
                  n19,f19, n20,f20, n21,f21, n22,f22, n23,f23, & 
                  n24,f24, n25,f25, n26,f26, n27,f27, n28,f28, & 
                  n29,f29, n30,f30, n31,f31, n32,f32, n33,f33, & 
                  n34,f34, n35,f35, n36,f36, n37,f37, n38,f38, & 
                  n39,f39, n40,f40, n41,f41, n42,f42, n43,f43, & 
                  n44,f44, n45,f45, n46,f46, n47,f47, n48,f48, & 
                  n49,f49, n50,f50 & 
               )

   integer :: mp, np
    
   character(len=*), optional :: & 
                        n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, & 
                        n11, n12, n13, n14, n15, n16, n17, n18, n19, n20, &
                        n21, n22, n23, n24, n25, n26, n27, n28, n29, n30, &
                        n31, n32, n33, n34, n35, n36, n37, n38, n39, n40, &
                        n41, n42, n43, n44, n45, n46, n47, n48, n49, n50  

   real, dimension(:,:), optional :: & 
                        f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, & 
                        f11, f12, f13, f14, f15, f16, f17, f18, f19, f20,& 
                        f21, f22, f23, f24, f25, f26, f27, f28, f29, f30,& 
                        f31, f32, f33, f34, f35, f36, f37, f38, f39, f40,& 
                        f41, f42, f43, f44, f45, f46, f47, f48, f49, f50
   
   character(len=30), dimension(:), allocatable :: CheckNames 
   real, dimension(:,:,:), allocatable :: CheckFields

   integer :: i, k

   allocate( farray_fields2( farray_nmax, mp, np) ) 
   i = 0

   if( present (n1) .AND. present(f1) ) then
      CALL fill_farray2( n1, f1, i )   
   else
      print *, "CABLE_log: cable_farray missing dummy args"
      return
   endif
   
   if( present (n2) .AND. present(f2) ) CALL fill_farray2( n2, f2, i )   
   if( present (n3) .AND. present(f3) ) CALL fill_farray2( n3, f3, i )   
   if( present (n4) .AND. present(f4) ) CALL fill_farray2( n4, f4, i )   
   if( present (n5) .AND. present(f5) ) CALL fill_farray2( n5, f5, i )   
   if( present (n6) .AND. present(f6) ) CALL fill_farray2( n6, f6, i )   
   if( present (n7) .AND. present(f7) ) CALL fill_farray2( n7, f7, i )   
   if( present (n8) .AND. present(f8) ) CALL fill_farray2( n8, f8, i )   
   if( present (n9) .AND. present(f9) ) CALL fill_farray2( n9, f9, i )   
   if( present (n10) .AND. present(f10) ) CALL fill_farray2( n10, f10, i )   
   if( present (n11) .AND. present(f11) ) CALL fill_farray2( n11, f11, i )   
   if( present (n12) .AND. present(f12) ) CALL fill_farray2( n12, f12, i )   
   if( present (n13) .AND. present(f13) ) CALL fill_farray2( n13, f13, i )   
   if( present (n14) .AND. present(f14) ) CALL fill_farray2( n14, f14, i )   
   if( present (n15) .AND. present(f15) ) CALL fill_farray2( n15, f15, i )   
   if( present (n16) .AND. present(f16) ) CALL fill_farray2( n16, f16, i )   
   if( present (n17) .AND. present(f17) ) CALL fill_farray2( n17, f17, i )   
   if( present (n18) .AND. present(f18) ) CALL fill_farray2( n18, f18, i )   
   if( present (n19) .AND. present(f19) ) CALL fill_farray2( n19, f19, i )   
   if( present (n20) .AND. present(f20) ) CALL fill_farray2( n20, f20, i )   
   if( present (n21) .AND. present(f21) ) CALL fill_farray2( n21, f21, i )   
   if( present (n22) .AND. present(f22) ) CALL fill_farray2( n22, f22, i )   
   if( present (n23) .AND. present(f23) ) CALL fill_farray2( n23, f23, i )   
   if( present (n24) .AND. present(f24) ) CALL fill_farray2( n24, f24, i )   
   if( present (n25) .AND. present(f25) ) CALL fill_farray2( n25, f25, i )   
   if( present (n26) .AND. present(f26) ) CALL fill_farray2( n26, f26, i )   
   if( present (n27) .AND. present(f27) ) CALL fill_farray2( n27, f27, i )   
   if( present (n29) .AND. present(f29) ) CALL fill_farray2( n29, f29, i )   
   if( present (n30) .AND. present(f30) ) CALL fill_farray2( n30, f30, i )   
   if( present (n31) .AND. present(f31) ) CALL fill_farray2( n31, f31, i )   
   if( present (n32) .AND. present(f32) ) CALL fill_farray2( n32, f32, i )   
   if( present (n33) .AND. present(f33) ) CALL fill_farray2( n33, f33, i )   
   if( present (n34) .AND. present(f34) ) CALL fill_farray2( n34, f34, i )   
   if( present (n35) .AND. present(f35) ) CALL fill_farray2( n35, f35, i )   
   if( present (n36) .AND. present(f36) ) CALL fill_farray2( n36, f36, i )   
   if( present (n37) .AND. present(f37) ) CALL fill_farray2( n37, f37, i )   
   if( present (n38) .AND. present(f38) ) CALL fill_farray2( n38, f38, i )   
   if( present (n39) .AND. present(f39) ) CALL fill_farray2( n39, f39, i )   
   if( present (n40) .AND. present(f40) ) CALL fill_farray2( n40, f40, i )   
   if( present (n41) .AND. present(f41) ) CALL fill_farray2( n41, f41, i )   
   if( present (n42) .AND. present(f42) ) CALL fill_farray2( n42, f42, i )   
   if( present (n50) .AND. present(f50) ) CALL fill_farray2( n50, f50, i )   
   
   allocate( CheckNames(i) )
   allocate( CheckFields(i,mp, np) )

   do k=1,i
      CheckNames(k) = farray_names(k) 
      CheckFields(k,:,:) = farray_fields2(k,:,:)
   enddo        

   deallocate( farray_fields2 )
    
END SUBROUTINE cable_farray2

subroutine fill_farray( n, f, i )   
   
   character(len=*) :: n
   real, dimension(:) :: f
   integer :: i
      
      i=i+1
      farray_names(i) = n
      farray_fields(i,:) = f
   
end subroutine fill_farray 


subroutine fill_farray2( n, f, i )   
   
   character(len=*) :: n
   real, dimension(:,:) :: f
   integer :: i
      
      i=i+1
      farray_names(i) = n
      farray_fields2(i,:,:) = f
   
end subroutine fill_farray2 


SUBROUTINE cable_extremes1(fname,field,mype)

   real, dimension(:,:) :: field
   character(len=*), dimension(:) :: fname
   
   integer, optional :: mype 
   integer :: i,j,k
   integer :: n,m,op
   real :: emax, emin, emean, emode
   real :: erange 
   real :: edbin 
   real, dimension(100) :: bin
   integer, dimension(100) :: ibin
   integer :: ib, ibmax, binmax, maxbin

   n = size(fname)
   m = size(field,2)

   ! for each field in fname(i)
   do i=1, n  
      
      emax =  MAXVAL( field(i,:) )
      emin =  MINVAL(field(i,:) )
      emean =  SUM(field(i,:) ) / ( m )
    
      erange = emax - emin
      edbin = erange / 100. ! for 100 bins

      bin(1) = emin          

      ! define bins per fname(i)
      do ib=2, 100
         bin(ib) = bin(ib-1) + edbin
      enddo

      ibin =0
      ! for each Element in field 
      do j=1, m  
      
         ! Assignn each Element to a bin 
         do ib=1, 99

            IF( field(i,j) >= bin(ib) .AND. &
                field(i,j) < bin(ib+1) ) THEN
               
               !if(ib==1 )print *, "jhan:field1 ", field(i,j)
               ibin(ib) = ibin(ib) + 1
            
            ENDIF   
          
         enddo ! DO LOOP over fill bins

      enddo ! DO LOOP over elements 
 
      binmax = 0 
      
      ! find max bin per field 
      do ib=1, 99

         IF( ibin(ib) > binmax ) THEN 
            binmax = ibin(ib) 
            maxbin = ib 
         ENDIF   
       
      enddo ! DO LOOP over bins
     
     print *, "jhan:bins1 ", bin 
     print *, "jhan:bins1 count", ibin 
     
      
      emode = bin(maxbin) 

         print *, ""
         print *, "CABLE_log: "
         print *, "   Field ", fname(i)
         print *, "   Min ", emin
         print *, "   Max ", emax
         print *, "   Mean ",emean
         print *, "   Mode ",emode
         print *, "End CABLE_log: "
         print *, ""
      
      enddo ! DO LOOP over fname(i)


END SUBROUTINE cable_extremes1


SUBROUTINE cable_extremes2(fname,field,mype)

   real, dimension(:,:,:) :: field
   character(len=*), dimension(:) :: fname
   
   integer, optional :: mype 
   integer :: i,j,k
   integer :: n,m,op
   real :: emax, emin, emean, emode
   real :: erange 
   real :: edbin 
   real, dimension(100) :: bin
   integer, dimension(100) :: ibin
   integer :: ib, ibmax, binmax, maxbin

   n = size(fname)
   m = size(field,2)
   op= size(field,3)

   ! for each field in fname(i)
   do i=1, n  
      
      emax =  MAXVAL( field(i,:,:) )
      emin =  MINVAL(field(i,:,:) )
      emean =  SUM(field(i,:,:) ) / ( m*op )
    
      erange = emax - emin
      edbin = erange / 100. ! for 100 bins

      bin(1) = emin          

      ! define bins per fname(i)
      do ib=2, 100
         bin(ib) = bin(ib-1) + edbin
      enddo

      ibin =0
      ! for each Element in field 
      do j=1, m  
         
         do k=1, op  
      
            ! Assignn each Element to a bin 
            do ib=1, 99

               IF( field(i,j,k) >= bin(ib) .AND. &
                   field(i,j,k) < bin(ib+1) ) THEN
                  !if(ib==1) print *, "jhan:field2 ", field(i,j,k)
                  ibin(ib) = ibin(ib) + 1
            
               ENDIF   
          
            enddo ! DO LOOP over fill bins
            
         enddo ! DO LOOP over elements 

      enddo ! DO LOOP over elements 
 
      binmax = 0 
      
      ! find max bin per field 
      do ib=1, 99

         IF( ibin(ib) > binmax ) THEN 
            binmax = ibin(ib) 
            maxbin = ib 
         ENDIF   
       
      enddo ! DO LOOP over bins
      
      emode = bin(maxbin) 

     print *, "jhan:bins2 ", bin 
     print *, "jhan:bins2 count", ibin 
     
         print *, ""
         print *, "CABLE_log: "
         print *, "   Field ", fname(i)
         print *, "   Min ", emin
         print *, "   Max ", emax
         print *, "   Mean ",emean
         print *, "   Mode ",emode
         print *, "End CABLE_log: "
         print *, ""
      
      enddo ! DO LOOP over fname(i)


END SUBROUTINE cable_extremes2

#ifndef UM_BUILD
  subroutine def_dims(nd, ncid, dimID, dim_len, dim_name )
    use netcdf
    implicit none
    integer, intent(in) :: nd, ncid
    character(len=*), dimension(:), intent(in) :: dim_name
    integer, dimension(:), intent(out) :: dimID
    integer, dimension(:), intent(in) :: dim_len
    integer :: j, ncok

    do j=1, nd
       ncok = NF90_DEF_DIM(ncid, trim(dim_name(j)), dim_len(j), dimID(j) )
       if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def dim ', dim_name(j))
    enddo

    return
  end subroutine def_dims




  subroutine def_vars(nv, ncid,  xtype, dimID, var_name,varID )
    use netcdf
    implicit none
    integer, intent(in) :: nv, ncid, xtype
    integer, dimension(:), intent(in) :: dimID
    integer, dimension(:), intent(inout) :: varID
    character(len=*), dimension(:), intent(in) :: var_name
    integer :: j, ncok

    ! lat
    ncok = NF90_DEF_VAR( ncid, trim(var_name(1)), xtype, &
         (/ dimID(1) /), varID(1))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(1))

    ! lon
    ncok = NF90_DEF_VAR(ncid, trim(var_name(2)), xtype, &
         (/ dimID(1) /), varID(2))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(2))

    ! tairk
    ncok = NF90_DEF_VAR(ncid, trim(var_name(3)), xtype, &
         (/ dimID(1), dimID(3) /), varID(3))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(3))

    !tsoil
    ncok = NF90_DEF_VAR(ncid, trim(var_name(4)), xtype, &
         (/ dimID(1), dimID(2),dimID(3)/), varID(4))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(4))

    ! moist
    ncok = NF90_DEF_VAR(ncid, trim(var_name(5)), xtype, &
         (/ dimID(1), dimID(2),dimID(3)/), varID(5))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(5))

    !cgpp
    ncok = NF90_DEF_VAR(ncid, trim(var_name(6)), xtype, &
         (/ dimID(1), dimID(3)/), varID(6))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(6))

    !crmplant
    ncok = NF90_DEF_VAR(ncid, trim(var_name(7)), xtype, &
         (/ dimID(1), dimID(2),dimID(3)/), varID(7))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(7))

    !phenphase
    ncok = NF90_DEF_VAR(ncid, trim(var_name(8)), xtype, &
         (/ dimID(1), dimID(3)/), varID(8))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(8))

    !doyphase1
    ncok = NF90_DEF_VAR(ncid, trim(var_name(9)), xtype, &
         (/ dimID(1), dimID(3)/), varID(9))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(9))

    !doyphase2
    ncok = NF90_DEF_VAR(ncid, trim(var_name(10)), xtype, &
         (/ dimID(1), dimID(3)/), varID(10))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(10))

    !doyphase3
    ncok = NF90_DEF_VAR(ncid, trim(var_name(11)), xtype, &
         (/ dimID(1), dimID(3)/), varID(11))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(11))

    !doyphase4
    ncok = NF90_DEF_VAR(ncid, trim(var_name(12)), xtype, &
         (/ dimID(1), dimID(3)/), varID(12))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(12))


    !mtemp
    ncok = NF90_DEF_VAR(ncid, trim(var_name(13)), xtype, &
         (/ dimID(1),dimID(3)/), varID(13))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(13))

    !Ndep
    ncok = NF90_DEF_VAR(ncid, trim(var_name(14)), xtype, &
         (/ dimID(1),dimID(3)/), varID(14))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(14))

    return
  end subroutine def_vars

  subroutine def_var_atts( ncfile_in, ncid, varID )
    use netcdf
    implicit none
    character(len=*), intent(in) :: ncfile_in
    integer, intent(in):: ncid       ! netcdf file ID
    integer, dimension(:), intent(in) :: varID ! (1) ~ tvair, (2) ~ pmb
    integer :: j, ncok
    character(len=10) dummy

    write(dummy,11) varID(1)
11  format(i2)
    ncok = NF90_PUT_ATT(ncid, nf90_global, "Title", "Forcing for define_air subroutine")
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def att ', ncfile_in)
    ncok = NF90_PUT_ATT(ncid, varID(3), "longname", "air temperature within canopy")
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def att ', dummy)
    ncok = NF90_PUT_ATT(ncid, varID(3), "units", "K")
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def att ', dummy)

    write(dummy,11) varID(2)


    return
  end subroutine def_var_atts


  subroutine put_var_ncr1(ncid, var_name, var )
    use netcdf
    use cable_def_types_mod, only : mp
    implicit none
    character(len=*), intent(in) ::  var_name
    real, dimension(:),intent(in) :: var
    integer, intent(in) :: ncid
    integer :: ncok, varID,j

    ncok = NF90_INQ_VARID(ncid, var_name, varId )
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'inquire var ', var_name)

    ncok = NF90_PUT_VAR(ncid, varId, var, start=(/1/), &
         count=(/mp/) )
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'putting var ', var_name)

  end subroutine put_var_ncr1


  subroutine put_var_ncr2(ncid, var_name, var, n_call )
    use netcdf
    use cable_def_types_mod, only : r_2, mp
    implicit none
    character(len=*), intent(in) ::  var_name
    real(r_2), dimension(:),intent(in) :: var
    integer, intent(in) :: ncid, n_call
    integer :: ncok, varID

    ncok = NF90_INQ_VARID(ncid, var_name, varId )
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'inquire var ', var_name)

    ncok = NF90_PUT_VAR(ncid, varId, var, start=(/1,n_call /), &
         count=(/mp,1/) )

    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'putting var ', var_name)

  end subroutine put_var_ncr2

  !soil vars
  subroutine put_var_ncr3(ncid, var_name, var, n_call, nl)
    use netcdf
    use cable_def_types_mod, only : r_2, mp, ms
    implicit none
    character(len=*), intent(in) :: var_name
    real(r_2), dimension(:,:),intent(in) :: var
    integer, intent(in) :: ncid, n_call, nl
    integer :: ncok, varID

    ncok = NF90_INQ_VARID( ncid, var_name, varId )
    IF( ncok /= nf90_noerr ) call stderr_nc(ncok,'inquire var ', var_name )

    ncok = NF90_PUT_VAR(ncid, varId, var, start=(/1,1,n_call /), &
         count=(/mp,nl,1/))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'putting var ', var_name)

    return
  end subroutine put_var_ncr3



  subroutine get_var_ncr2(ncid, var_name, var, n_call )
    use netcdf
    use cable_def_types_mod, only : r_2,mp
    implicit none
    character(len=*), intent(in) :: var_name
    real(r_2), dimension(:),intent(out) :: var
    integer, intent(in) :: ncid
    integer :: ncok, varID, n_call
    real, dimension(mp) :: temp

    temp = 0.

    ncok = NF90_INQ_VARID(ncid, var_name, varId )
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'inquire var ', var_name)
    ncok = NF90_GET_VAR(ncid, varId, temp, start=(/1,n_call/), &
         count=(/mp,1/) )

    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'getting var ', var_name)

    var = real( temp, r_2 )
  end subroutine get_var_ncr2

  subroutine get_var_ncr3(ncid, var_name, var, n_call, nl )
    use netcdf
    use cable_def_types_mod, only : r_2, mp, ms
    implicit none
    character(len=*), intent(in) :: var_name
    real(r_2), dimension(:,:),intent(out) :: var
    integer, intent(in) :: ncid, n_call, nl
    integer :: ncok, varID
    real, dimension(mp,1:nl) :: temp

    ncok = NF90_INQ_VARID(ncid, var_name, varId )
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'inquire var ', var_name)

    ncok = NF90_GET_VAR(ncid, varId, temp, start=(/1,1,n_call /), &
         count=(/mp, nl, 1/))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'putting var ', var_name)
    var = real( temp, r_2 )
  end subroutine get_var_ncr3



  subroutine stderr_nc(status,message, var)
    use netcdf
    character(len=*), intent(in) :: message, var
    INTEGER, INTENT(IN) :: status
    character(len=7) :: err_mess
    err_mess = 'ERROR:'
    print *, (err_mess//message), var
    PRINT*,NF90_STRERROR(status)
    stop
  end subroutine stderr_nc
#endif


END MODULE cable_diag_module



