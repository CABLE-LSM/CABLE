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
  USE cable_def_types_mod, ONLY : r_2
  USE cable_fFile_module, ONLY : fprintf_dir

  IMPLICIT NONE
  INTEGER, PARAMETER :: gok=0
  INTEGER :: galloctest=1

  !--- subrs overloaded to respond to call cable_diag
  INTERFACE cable_diag
     MODULE PROCEDURE cable_diag1
  END INTERFACE
  !CABLE_LSM:"A" version of diagnostics along the way. farray builds arrays
  INTERFACE cable_farray
     MODULE PROCEDURE cable_farray1, cable_farray2
  END INTERFACE

  INTERFACE cable_NaN
     MODULE PROCEDURE cable_NaN1, cable_NaN2
  END INTERFACE

#ifndef UM_BUILD
  INTERFACE put_var_nc
     MODULE PROCEDURE put_var_ncr1, put_var_ncr2, put_var_ncr3
  END INTERFACE

  INTERFACE get_var_nc
     MODULE PROCEDURE get_var_ncr2, get_var_ncr3
  END INTERFACE

#endif
  !CABLE_LSM: procedures for writing out vars in seperate text files
  INTERFACE cable_fprintf
     MODULE PROCEDURE cable_fprintf1, cable_fprintf2, cable_fprintf3,  &
          cable_iprintf1, cable_iprintf2, cable_Lprintf2
  END INTERFACE

  !CABLE_LSM: farray builds array for _diags. keep for ref.
  INTEGER, PARAMETER ::                     &
       farray_nmax = 50

  CHARACTER(len=30), DIMENSION(farray_nmax) :: &
       farray_names

  REAL, DIMENSION(:,:), ALLOCATABLE :: &
       farray_fields

  REAL, DIMENSION(:,:,:), ALLOCATABLE :: &
       farray_fields2

  !CABLE_LSM: make avail for diag. check
  INTEGER :: frow_length, frows, fland_pts, fntiles
  REAL, DIMENSION(:), ALLOCATABLE :: fFland
  REAL, DIMENSION(:,:), ALLOCATABLE :: ftile_frac

CONTAINS

  !==========================================================================!
  !CABLE_LSM: make args avail for diag. check across CABLE
  SUBROUTINE cable_fsend( row_length, rows, land_pts, ntiles,Fland, tile_frac )
    INTEGER :: row_length, rows, land_pts, ntiles
    REAL, DIMENSION(land_pts) :: Fland
    REAL, DIMENSION(land_pts, ntiles) :: tile_frac

    IF(.NOT. ALLOCATED(fFland) ) ALLOCATE( fFland(land_pts) )
    IF(.NOT. ALLOCATED(ftile_frac) ) ALLOCATE( ftile_frac(land_pts, ntiles) )

    frow_length = row_length
    frows = rows
    fland_pts = land_pts
    fntiles = ntiles
    fFland =  Fland
    ftile_frac =  tile_frac

  END SUBROUTINE cable_fsend

  !CABLE_LSM: procedures for writing out vars in seperate text files
  SUBROUTINE fprint_fland( iDiag, dir, basename, node )

    USE cable_common_module
    ! IN vars
    INTEGER :: iDiag  ! f^n creates unique unidID to be returned to calling point

    ! writes file per processor (basename+node)
    INTEGER :: node                  ! processor number
    CHARACTER(len=*) :: dir          ! dir name based on var
    CHARACTER(len=*) :: basename     ! filename based on var
    CHARACTER(len=29) :: fbasename     ! filename based on var

    ! LOCAL vars
    REAL :: fsum
    INTEGER, SAVE :: pDiag=713       ! give unique SEED per module procedure
    INTEGER :: j, jctr, k            ! local counters

    fbasename = basename
    ! Returns unique unit=iDiag and modified basename
    CALL open_file_per_node( iDiag, pDiag, dir, basename, node, fbasename=fbasename )

    IF (node==0) THEN
       WRITE(iDiag,*) "This meesage is only written for node==0 "
       WRITE(iDiag,*) "Checking over identified land points that tile fraction "
       WRITE(iDiag,*) "sums to 1 (within tolerance). And is zero for non-land"
       WRITE(iDiag,*) "Write(s) follow(s) ONLY for points where these "
       WRITE(iDiag,*) "conditions NOT met"
       WRITE(iDiag,*) "Jhan:Check that Total land points = 0 is OK"
       WRITE(iDiag,*) ""
    ENDIF

    jctr=0
    DO j=1, fland_pts
       IF( fFland(j)> 0.) THEN
          jctr = jctr + 1
          fsum = SUM(ftile_frac(j,:))
          IF(fsum < 0.9999 .OR. fsum > 1.001 ) THEN
             WRITE(iDiag,*) "Summed tile_frac is: ", fsum
             WRITE(iDiag,*) "    for this land_pt ", j
             WRITE(iDiag,*) ""
          ENDIF
       ENDIF
       IF( SUM(ftile_frac(j,:)) > 0. .AND. fFland(j) == 0.) THEN
          WRITE(iDiag,*) "Report: Summed tile_frac is: ", SUM( ftile_frac(j,:) )
          WRITE(iDiag,*) "    BUT Fland is zero "
          WRITE(iDiag,*) ""
       ENDIF
    ENDDO

    !write(iDiag,*) "Total land points ", jctr

    CLOSE(iDiag)

    CALL remove_empty_file( iDiag, fbasename )


  END SUBROUTINE fprint_fland

  !==========================================================================!
  ! writes text files Interfaced Module Procedures
  !==========================================================================!

  ! 1-D REAL
  SUBROUTINE cable_fprintf1( iDiag, basename, var1, dimx, L_fprint )
    USE cable_common_module
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

    CALL cable_fextremes1( idiag, dimx, var1, ktau_gl )

  END SUBROUTINE cable_fprintf1


  !SUBROUTINE cable_fextremes1(fname,field,mype)
  SUBROUTINE cable_fextremes1(idiag,dimx,field, ktau)

    INTEGER :: idiag
    INTEGER :: dimx
    REAL, DIMENSION(dimx) :: field
    INTEGER :: ktau

    INTEGER, PARAMETER :: nbins =3
    REAL :: emax, emin, emean, emode
    REAL :: erange
    REAL, DIMENSION(nbins) :: bin

    INTEGER :: i,j,k
    INTEGER :: n,m,op
    REAL :: edbin
    INTEGER, DIMENSION(nbins) :: ibin
    INTEGER :: ib, ibmax, binmax, maxbin
    !logical,save :: first_call=.true.

    !T!if( ktau==1 .OR. mod( ktau,10)==0 ) then
    WRITE (iDiag,*) ""
    WRITE (iDiag,*) "timestep ", ktau
    !T!endif

    IF( dimx > 0) THEN
       emax =  MAXVAL( field )
       emin =  MINVAL(field )
       emean =  SUM(field(:) ) / ( dimx )
       IF( ktau==1 ) &
            WRITE (iDiag,*) "dimx ", dimx
    ELSE
       emax = 123.123
       emin = 123.123
       emean = 123.123
       IF( ktau==1 ) &
            WRITE (iDiag,*) "dimx is zero here"
    ENDIF

    erange = emax - emin
    edbin = erange / nbins

    bin(1) = emin

    ! define bins per field
    DO ib=2, nbins
       bin(ib) = bin(ib-1) + edbin
    ENDDO

    ibin =0
    ! for each Element in field
    DO j=1, dimx

       ! Assignn each Element to a bin
       DO ib=1,(nbins-1)

          IF( field(j) >= bin(ib) .AND. &
               field(j) < bin(ib+1) ) THEN

             !if(ib==1 )print *, "jhan:field1 ", field(i,j)
             ibin(ib) = ibin(ib) + 1
          ENDIF

       ENDDO ! DO LOOP over fill bins

    ENDDO ! DO LOOP over elements

    !jhan:reform bin(ib) to value at middle of bin
17  FORMAT(  "bin(", I2.1, ") ", E15.6, 4X, "ibin(", I2.1, ") ", I6.1 )

    binmax = 0
    maxbin = 1
    ! find max bin per field
    DO ib=1, (nbins-1)

       !C!if( ktau==1 ) &
       !C!write (iDiag,17) ib, bin(ib), ib, ibin(ib)

       ! ibin(ib) = the # of values in each bin
       IF( ibin(ib) > binmax ) THEN
          binmax = ibin(ib)
          maxbin = ib ! this is actually the starting index of the bin
       ENDIF

    ENDDO ! DO LOOP over bins

    !if( ktau==1 ) &
    !  write (iDiag,*) "ib: maxbin ", maxbin

    emode = bin(maxbin)

18  FORMAT(  " Min.", 18X, " Max.", 18X, " Mean " )
19  FORMAT(  E15.6, 6X, E15.6, 6X, E15.6 )
    IF( ktau==1 .OR. MOD( ktau,10)==0 ) &
         WRITE (iDiag,18)
    WRITE (iDiag,19) emin, emax, emean

    !18 format(  " Min.", 18X, " Max.", 18X, " Mean ", 18X, " Mode" )
    !  write (iDiag,19) emin, emax, emean, emode
    !19 format(  E15.6, 6X, E15.6, 6X, E15.6, 6X, E15.6 )


    !19 format(  "     ", I6.1, "     ", E15.6, "     ", E15.6, "     ", E15.6 )
    !19 format(  "ktau: ", I6.1, " Min.", E15.6, " Max.", E15.6, " Med.", E15.6 )

  END SUBROUTINE cable_fextremes1




  ! 2-D REAL
  SUBROUTINE cable_fprintf2( iDiag, basename, dimx, dimy, timestep, node, &
       dir, var1 )
    USE cable_common_module
    ! IN vars
    INTEGER :: iDiag  ! f^n creates unique unidID to be returned to calling point
    INTEGER :: dimx   ! 1-D length
    INTEGER :: dimy   ! 2-D length
    INTEGER :: timestep
    REAL, DIMENSION(dimx,dimy) :: var1    ! var CALLed
    CHARACTER(len=*) :: dir ! obsolete here - remove

    ! writes file per processor (basename+node)
    INTEGER :: node                  ! processor number
    CHARACTER(len=*) :: basename     ! filename based on var

    ! LOCAL vars
    INTEGER, SAVE :: pDiag=2713      ! give unique SEED per module procedure
    INTEGER :: j, k                  ! local counters
    INTEGER :: r, rl                 ! local counters
    LOGICAL, SAVE :: first_call = .TRUE.

    ! Returns unique unit=iDiag and modified basename
    CALL open_file_per_node( iDiag, pDiag, dir, basename, node )

    ! Writes when a new timestep
    CALL check_timestep( iDiag, timestep )

    IF( first_call )  WRITE (iDiag,*) "x, y, z"
    IF( dimx==fland_pts .AND. dimy==fntiles ) THEN
       IF( first_call )  &
            WRITE (iDiag,*) "Filtered: tile fraction > 0 "
    ELSEIF( dimx==frow_length .AND. dimy==frows ) THEN
       IF( first_call )  &
            WRITE (iDiag,*) "Filtered: NOT as yet "
    ELSE
       IF( first_call )  &
            WRITE (iDiag,*) "dimensions not recognized ", dimx, dimy
       RETURN
    ENDIF

    first_call = .FALSE.

17  FORMAT (2I3.1,ES15.6)
    DO j=1,dimx
       DO k=1,dimy
          IF( dimx==fland_pts .AND. dimy==fntiles ) THEN
             IF( ftile_frac(j,k) > 0. ) WRITE (iDiag,17) j,k, var1(j,k)
          ELSEIF( dimx==frow_length .AND. dimy==frows ) THEN
             !if( ftile_frac(j,k) > 0. )
             WRITE (iDiag,17) j,k, var1(j,k)
          ELSE
             WRITE (iDiag,*) "Should never get here"
             RETURN
          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE cable_fprintf2

  ! 3-D REAL
  SUBROUTINE cable_fprintf3( iDiag, basename, dimx, dimy, dimz, timestep, node, &
       dir, var1 )
    USE cable_common_module
    ! IN vars
    INTEGER :: iDiag  ! f^n creates unique unidID to be returned to calling point
    INTEGER :: dimx   ! 1-D length
    INTEGER :: dimy   ! 2-D length
    INTEGER :: dimz   ! 3-D length
    INTEGER :: timestep
    REAL, DIMENSION(dimx,dimy,dimz) :: var1    ! var CALLed
    CHARACTER(len=*) :: dir ! obsolete here - remove

    ! writes file per processor (basename+node)
    INTEGER :: node                  ! processor number
    CHARACTER(len=*) :: basename     ! filename based on var

    ! LOCAL vars
    INTEGER, SAVE :: pDiag=3713      ! give unique SEED per module procedure
    INTEGER :: j, k, m               ! local counters
    LOGICAL, SAVE :: first_call = .TRUE.

    ! Returns unique unit=iDiag and modified basename
    CALL open_file_per_node( iDiag,pDiag, dir, basename, node )

    ! Writes when a new timestep
    CALL check_timestep( iDiag, timestep )

    IF( first_call )  THEN
       WRITE (iDiag,*) "x, y, z"
       IF( dimx==fland_pts .AND. dimy==fntiles ) THEN
          WRITE (iDiag,*) "Filtered: tile fraction > 0 "
       ELSEIF( dimx==frow_length .AND. dimy==frows ) THEN
          WRITE (iDiag,*) "Filtered: NOT as yet "
       ELSE
          WRITE (iDiag,*) "dimensions not recognized ", dimx, dimy
       ENDIF
    ENDIF
    first_call = .FALSE.

    DO j=1,dimx
       DO k=1,dimy
          DO m=1,dimz
             IF( dimx==fland_pts .AND. dimy==fntiles ) THEN
                IF( ftile_frac(j,k) > 0. ) THEN
                   WRITE (iDiag,*) j,k,m
                   WRITE (iDiag,*) var1(j,k,m)
                ENDIF
             ELSEIF( dimx==frow_length .AND. dimy==frows ) THEN
                WRITE (iDiag,*) j,k,m
                WRITE (iDiag,*) var1(j,k,m)
             ELSE
                WRITE (iDiag,*) j,k,m
                WRITE (iDiag,*) var1(j,k,m)
             ENDIF
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE cable_fprintf3


  ! 1-D INTEGER
  SUBROUTINE cable_iprintf1( iDiag, basename, dimx, timestep, node, &
       dir, var1 )
    USE cable_common_module
    ! IN vars
    INTEGER :: iDiag  ! f^n creates unique unidID to be returned to calling point
    INTEGER :: dimx   ! 1-D length
    INTEGER :: timestep
    INTEGER, DIMENSION(dimx) :: var1    ! var CALLed
    CHARACTER(len=*) :: dir ! obsolete here - remove

    ! writes file per processor (basename+node)
    INTEGER :: node                  ! processor number
    CHARACTER(len=*) :: basename     ! filename based on var

    ! LOCAL vars
    REAL, DIMENSION(dimx) :: invar1  ! store var read in

    INTEGER, SAVE :: pDiag=11713      ! give unique SEED per module procedure

    INTEGER :: j, inj            ! local counters
    LOGICAL, SAVE :: first_call = .TRUE.

    ! Returns unique unit=iDiag and modified basename
    CALL open_file_per_node( iDiag,pDiag, dir, basename, node )

    ! Writes when a new timestep
    CALL check_timestep( iDiag, timestep )

    IF( first_call )  THEN
       WRITE (iDiag,*) "x, y, z"
       IF( dimx==fland_pts ) THEN
          WRITE (iDiag,*) "Filtered: land fraction > 0 "
       ELSE
          WRITE (iDiag,*) "dimensions not recognized ", dimx
       ENDIF
    ENDIF
    first_call = .FALSE.

    DO j=1,dimx
       IF( dimx==fland_pts ) THEN
          IF( fFland(j) > 0. ) THEN
             WRITE (iDiag,*) j
             WRITE (iDiag,*) var1(j)
          ENDIF
       ELSE
          WRITE (iDiag,*) j
          WRITE (iDiag,*) var1(j)
       ENDIF
    ENDDO

  END SUBROUTINE cable_iprintf1


  ! 2-D INTEGER
  SUBROUTINE cable_iprintf2( iDiag, basename, dimx, dimy, timestep, node, &
       dir, var1 )
    USE cable_common_module
    ! IN vars
    INTEGER :: iDiag  ! f^n creates unique unidID to be returned to calling point
    INTEGER :: dimx   ! 1-D length
    INTEGER :: dimy   ! 2-D length
    INTEGER :: timestep
    INTEGER, DIMENSION(dimx,dimy) :: var1    ! var CALLed
    CHARACTER(len=*) :: dir ! obsolete here - remove

    ! writes file per processor (basename+node)
    INTEGER :: node                  ! processor number
    CHARACTER(len=*) :: basename     ! filename based on var

    ! LOCAL vars
    INTEGER, SAVE :: pDiag=21713      ! give unique SEED per module procedure
    INTEGER :: j, k                  ! local counters
    LOGICAL, SAVE :: first_call = .TRUE.

    ! Returns unique unit=iDiag and modified basename
    CALL open_file_per_node( iDiag, pDiag, dir, basename, node )

    ! Writes when a new timestep
    CALL check_timestep( iDiag, timestep )

    IF( first_call )  THEN
       WRITE (iDiag,*) "x, y, z"
       IF( dimx==fland_pts .AND. dimy==fntiles ) THEN
          WRITE (iDiag,*) "Filtered: tile fraction > 0 "
       ELSEIF( dimx==frow_length .AND. dimy==frows ) THEN
          WRITE (iDiag,*) "Filtered: NOT as yet "
       ELSE
          WRITE (iDiag,*) "dimensions not recognized ", dimx, dimy
       ENDIF
    ENDIF
    first_call = .FALSE.

    DO j=1,dimx
       DO k=1,dimy
          IF( dimx==fland_pts .AND. dimy==fntiles ) THEN
             IF( ftile_frac(j,k) > 0. ) THEN
                WRITE (iDiag,*) j,k
                WRITE (iDiag,*) var1(j,k)
             ENDIF
          ELSEIF( dimx==frow_length .AND. dimy==frows ) THEN
             WRITE (iDiag,*) j,k
             WRITE (iDiag,*) var1(j,k)
          ELSE
             WRITE (iDiag,*) j,k
             WRITE (iDiag,*) var1(j,k)
          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE cable_iprintf2


  ! 2-D LOGICAL
  SUBROUTINE cable_Lprintf2( iDiag, basename, dimx, dimy, timestep, node, &
       dir, var1 )
    USE cable_common_module
    ! IN vars
    INTEGER :: iDiag  ! f^n creates unique unidID to be returned to calling point
    INTEGER :: dimx   ! 1-D length
    INTEGER :: dimy   ! 2-D length
    INTEGER :: timestep
    LOGICAL, DIMENSION(dimx,dimy) :: var1    ! var CALLed
    CHARACTER(len=*) :: dir ! obsolete here - remove

    ! writes file per processor (basename+node)
    INTEGER :: node                  ! processor number
    CHARACTER(len=*) :: basename     ! filename based on var

    ! LOCAL vars
    INTEGER, SAVE :: pDiag=41713      ! give unique SEED per module procedure
    INTEGER :: j, k                  ! local counters
    LOGICAL, SAVE :: first_call = .TRUE.

    ! Returns unique unit=iDiag and modified basename
    CALL open_file_per_node( iDiag, pDiag, dir, basename, node )

    ! Writes when a new timestep
    CALL check_timestep( iDiag, timestep )

    IF( first_call )  THEN
       WRITE (iDiag,*) "x, y, z"
       IF( dimx==fland_pts .AND. dimy==fntiles ) THEN
          WRITE (iDiag,*) "Filtered: tile fraction > 0 "
       ELSEIF( dimx==frow_length .AND. dimy==frows ) THEN
          WRITE (iDiag,*) "Filtered: NOT as yet "
       ELSE
          WRITE (iDiag,*) "dimensions not recognized ", dimx, dimy
       ENDIF
    ENDIF
    first_call = .FALSE.

    DO j=1,dimx
       DO k=1,dimy
          IF( dimx==fland_pts .AND. dimy==fntiles ) THEN
             IF( ftile_frac(j,k) > 0. ) THEN
                WRITE (iDiag,*) j,k
                WRITE (iDiag,*) var1(j,k)
             ENDIF
          ELSEIF( dimx==frow_length .AND. dimy==frows ) THEN
             WRITE (iDiag,*) j,k
             WRITE (iDiag,*) var1(j,k)
          ELSE
             WRITE (iDiag,*) j,k
             WRITE (iDiag,*) var1(j,k)
          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE cable_Lprintf2

  !==========================================================================!
  !==========================================================================!

  SUBROUTINE open_file_per_node( iDiag,pDiag, dir, basename, node, fbasename )
    USE cable_common_module
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
       STOP
    ENDIF

  END SUBROUTINE open_file_per_node

  SUBROUTINE open_iDiag( iDiag, infilename, gopenstatus)

    USE cable_common_module
    INTEGER :: iDiag
    INTEGER :: gopenstatus
    CHARACTER(len=*) :: infilename
    CHARACTER(len=300) :: ffilename

    ffilename=TRIM( TRIM(infilename)// '.txt' )
    !if( cable_user%check_write ) then
    OPEN( unit=iDiag, file=ffilename, status="replace", &
         action="write", iostat=gopenstatus, form="formatted", &
         position='append' )
    !endif
    !if( cable_user%check_read ) then
    !   open(unit=iDiag,file=trim(infilename)//'.txt',status="old", &
    !      action="read", iostat=gopenstatus )
    !endif

  END SUBROUTINE open_iDiag

  SUBROUTINE remove_empty_file( iDiag, fbasename )

    USE cable_common_module
    INTEGER :: iDiag
    CHARACTER(len=*) :: fbasename
    INTEGER :: gopenstatus, gios
    CHARACTER(len=99) :: blnk
    CHARACTER(len=49) :: cmd

    OPEN(unit=iDiag,status="unknown", file=fbasename, &
         action="read", iostat=gopenstatus, form="formatted", position='append' )
    READ(idiag,*, iostat=gios) blnk
    CLOSE(idiag)
    IF( gios <0 ) THEN
       CALL unlink(fbasename)
    ENDIF

    RETURN

  END SUBROUTINE remove_empty_file

  SUBROUTINE check_timestep( iDiag, timestep )
    ! IN vars
    INTEGER :: iDiag  ! f^n creates unique unidID to be returned to calling point
    INTEGER :: timestep

    INTEGER, SAVE ::otimestep = 0
    !otimestep = 0
    !write (iDiag,*) "otimestep ", otimestep
    !write (iDiag,*) "timestep ", timestep
    IF(otimestep .NE. timestep ) THEN
       otimestep=timestep
       !if( timestep==1 .OR. mod( timestep,10)==0 ) then
       !  write (iDiag,*) ""
       !  write (iDiag,*) "timestep ", timestep
       !endif
    ELSE
       RETURN
    ENDIF

    RETURN

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
    INTEGER, INTENT(inOUT) :: iDiag
    INTEGER, SAVE :: pDiag=713
    INTEGER, INTENT(in) :: dimx, dimy, timestep,node
    REAL, INTENT(in), DIMENSION(:) :: var1
    INTEGER, OPTIONAL :: once
    INTEGER :: Nvars=1 !this WAS input
    INTEGER :: i=0
    CHARACTER(len=*), INTENT(in) :: basename, vname1
    CHARACTER(len=30) :: filename, chnode

    IF(iDiag==0) THEN
       pDiag = pDiag+2
       iDiag=pDiag
    ENDIF

    WRITE(chnode,10) node
10  FORMAT(i3.3)
    filename=TRIM(TRIM(basename)//TRIM(chnode))

    IF (timestep == 1) &
         CALL cable_diag_desc1( iDiag, TRIM(filename), dimx, dimy, vname1 )

    IF( PRESENT(once) ) THEN
       IF (timestep == 1) &
                                ! write data only on first timestep
            CALL cable_diag_data1( iDiag, TRIM(filename), dimx, timestep, dimy, &
            var1 )
    ELSE
       ! write data every timestep
       CALL cable_diag_data1( iDiag, TRIM(filename), dimx, timestep, dimy, &
            var1 )
    ENDIF

  END SUBROUTINE cable_diag1

  !=============================================================================!
  !=============================================================================!

  SUBROUTINE cable_diag_desc1( iDiag, filename, dimx, dimy, vname1 )

    INTEGER, INTENT(in) :: iDiag,dimx,dimy
    INTEGER, PARAMETER :: Nvars=1
    CHARACTER(len=*), INTENT(in) :: filename, vname1
    INTEGER, SAVE :: gopenstatus = 1

    OPEN(unit=iDiag,file=filename//'.dat', status="replace", &
         action="write", iostat=gopenstatus )

    IF(gopenstatus==gok) THEN
       WRITE (iDiag,*) 'Number of var(s): '
       WRITE (iDiag,*) Nvars
       WRITE (iDiag,*) 'Name of var(s): '
       WRITE (iDiag,7139) vname1
7139   FORMAT(a)
       WRITE (iDiag,*) 'dimension of var(s) in x: '
       WRITE (iDiag,*) dimx
       WRITE (iDiag,*) 'dimension of var(s) in y: '
       WRITE (iDiag,*) dimy
    ELSE
       WRITE (*,*) filename//'.dat',' Error: unable to write'
    ENDIF

    CLOSE(iDiag)

  END SUBROUTINE cable_diag_desc1


  SUBROUTINE cable_diag_data1( iDiag, filename, dimx, timestep, kend, var1  )

    INTEGER, INTENT(in) :: iDiag, dimx, timestep, kend
    INTEGER, PARAMETER :: Nvars=1
    REAL, INTENT(in), DIMENSION(:) :: var1
    CHARACTER(len=*), INTENT(in) :: filename
    INTEGER, SAVE :: gopenstatus = 1

    IF (timestep == 1)  THEN
       OPEN(unit=iDiag+1,file=filename//'.bin',status="unknown", &
            action="write", iostat=gopenstatus, form="unformatted", &
            position='append' )
    ENDIF

    IF(gopenstatus==gok) THEN
       WRITE (iDiag+1) var1
    ELSE
       WRITE (*,*) filename//'.bin',' NOT open for write. Error'
    ENDIF

    IF (timestep == kend) &
         CLOSE(iDiag+1)

  END SUBROUTINE cable_diag_data1

  !==========================================================================!
  !--- cable generic print status
  !==========================================================================!

  SUBROUTINE cable_stat( routname)
    USE cable_common_module, ONLY : ktau_gl, knode_gl

    CHARACTER(len=*), INTENT(in) :: routname
    IF(knode_gl==1) &
         WRITE(6,*) 'CABLE@  ', routname, ktau_gl

  END SUBROUTINE cable_stat

  !==========================================================================!
  !--- cable status NaN
  !==========================================================================!


  SUBROUTINE cable_NaN1(fname,field,mype)

    REAL, DIMENSION(:,:) :: field
    CHARACTER(len=*), DIMENSION(:) :: fname

    INTEGER, OPTIONAL :: mype
    INTEGER :: i,j
    LOGICAL :: NoNaN, check
    INTEGER :: n,m

    n = SIZE(fname)
    m = SIZE(field,2)
    check = .FALSE.

    DO i=1, n

       NoNaN = .TRUE.

       DO j=1, m

          CALL isnan(  field(i,j), check )

          IF( check ) THEN
             PRINT *, ""
             PRINT *, "CABLE_log: "
             IF( PRESENT(mype) ) PRINT *, "proc # ", mype
             PRINT *, "   Element: ",j
             PRINT *, "   of field ", fname(i)
             PRINT *, "   is NaN"
             PRINT *, "End CABLE_log: "
             PRINT *, ""
             NoNaN = .FALSE.
          END IF

       ENDDO

       IF(NoNaN) THEN
          PRINT *, ""
          PRINT *, "CABLE_log: "
          IF( PRESENT(mype) ) PRINT *, "proc # ", mype
          PRINT *, '   Field: ', fname(i)
          PRINT *, "   is clear of NaNs"
          PRINT *, "End CABLE_log: "
          PRINT *, ""
       END IF

    ENDDO

  END SUBROUTINE cable_NaN1


  SUBROUTINE cable_NaN2(fname,field,mype)

    REAL, DIMENSION(:,:,:) :: field
    CHARACTER(len=*), DIMENSION(:) :: fname

    INTEGER, OPTIONAL :: mype
    INTEGER :: i,j,k
    LOGICAL :: NoNaN, check
    INTEGER :: n,m,op

    n = SIZE(fname)
    m = SIZE(field,2)
    op = SIZE(field,3)
    check = .FALSE.

    DO i=1, n

       NoNaN = .TRUE.

       DO j=1, m

          DO k=1, op

             CALL isnan(  field(i,j,k), check )

             IF( check ) THEN
                PRINT *, ""
                PRINT *, "CABLE_log: "
                IF( PRESENT(mype) ) PRINT *, "proc # ", mype
                PRINT *, "   Element: ",j,op
                PRINT *, "   of field ", fname(i)
                PRINT *, "   is NaN"
                PRINT *, "End CABLE_log: "
                PRINT *, ""
                NoNaN = .FALSE.
             END IF

          ENDDO

       ENDDO

       IF(NoNaN) THEN
          PRINT *, ""
          PRINT *, "CABLE_log: "
          IF( PRESENT(mype) ) PRINT *, "proc # ", mype
          PRINT *, '   Field: ', fname(i)
          PRINT *, "   is clear of NaNs"
          PRINT *, "End CABLE_log: "
          PRINT *, ""
       END IF

    ENDDO

  END SUBROUTINE cable_NaN2





  SUBROUTINE isnan(var, check)
    REAL :: var
    LOGICAL :: check

    IF (var .NE. var) THEN
       check = .TRUE.
    ELSE
       check = .FALSE.
    END IF

  END SUBROUTINE isnan


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

    INTEGER :: mp

    CHARACTER(len=*), OPTIONAL :: &
         n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, &
         n11, n12, n13, n14, n15, n16, n17, n18, n19, n20, &
         n21, n22, n23, n24, n25, n26, n27, n28, n29, n30, &
         n31, n32, n33, n34, n35, n36, n37, n38, n39, n40, &
         n41, n42, n43, n44, n45, n46, n47, n48, n49, n50

    REAL, DIMENSION(:), OPTIONAL :: &
         f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, &
         f11, f12, f13, f14, f15, f16, f17, f18, f19, f20,&
         f21, f22, f23, f24, f25, f26, f27, f28, f29, f30,&
         f31, f32, f33, f34, f35, f36, f37, f38, f39, f40,&
         f41, f42, f43, f44, f45, f46, f47, f48, f49, f50

    CHARACTER(len=30), DIMENSION(:), ALLOCATABLE :: CheckNames
    REAL, DIMENSION(:,:), ALLOCATABLE :: CheckFields

    INTEGER :: i, k

    ALLOCATE( farray_fields( farray_nmax, mp) )
    i = 0

    IF( PRESENT (n1) .AND. PRESENT(f1) ) THEN
       CALL fill_farray( n1, f1, i )
    ELSE
       PRINT *, "CABLE_log: cable_farray missing dummy args"
       RETURN
    ENDIF

    IF( PRESENT (n2) .AND. PRESENT(f2) ) CALL fill_farray( n2, f2, i )
    IF( PRESENT (n3) .AND. PRESENT(f3) ) CALL fill_farray( n3, f3, i )
    IF( PRESENT (n4) .AND. PRESENT(f4) ) CALL fill_farray( n4, f4, i )
    IF( PRESENT (n5) .AND. PRESENT(f5) ) CALL fill_farray( n5, f5, i )
    IF( PRESENT (n6) .AND. PRESENT(f6) ) CALL fill_farray( n6, f6, i )
    IF( PRESENT (n7) .AND. PRESENT(f7) ) CALL fill_farray( n7, f7, i )
    IF( PRESENT (n8) .AND. PRESENT(f8) ) CALL fill_farray( n8, f8, i )
    IF( PRESENT (n9) .AND. PRESENT(f9) ) CALL fill_farray( n9, f9, i )
    IF( PRESENT (n10) .AND. PRESENT(f10) ) CALL fill_farray( n10, f10, i )
    IF( PRESENT (n11) .AND. PRESENT(f11) ) CALL fill_farray( n11, f11, i )
    IF( PRESENT (n12) .AND. PRESENT(f12) ) CALL fill_farray( n12, f12, i )
    IF( PRESENT (n13) .AND. PRESENT(f13) ) CALL fill_farray( n13, f13, i )
    IF( PRESENT (n14) .AND. PRESENT(f14) ) CALL fill_farray( n14, f14, i )
    IF( PRESENT (n15) .AND. PRESENT(f15) ) CALL fill_farray( n15, f15, i )
    IF( PRESENT (n16) .AND. PRESENT(f16) ) CALL fill_farray( n16, f16, i )
    IF( PRESENT (n17) .AND. PRESENT(f17) ) CALL fill_farray( n17, f17, i )
    IF( PRESENT (n18) .AND. PRESENT(f18) ) CALL fill_farray( n18, f18, i )
    IF( PRESENT (n19) .AND. PRESENT(f19) ) CALL fill_farray( n19, f19, i )
    IF( PRESENT (n20) .AND. PRESENT(f20) ) CALL fill_farray( n20, f20, i )
    IF( PRESENT (n21) .AND. PRESENT(f21) ) CALL fill_farray( n21, f21, i )
    IF( PRESENT (n22) .AND. PRESENT(f22) ) CALL fill_farray( n22, f22, i )
    IF( PRESENT (n23) .AND. PRESENT(f23) ) CALL fill_farray( n23, f23, i )
    IF( PRESENT (n24) .AND. PRESENT(f24) ) CALL fill_farray( n24, f24, i )
    IF( PRESENT (n25) .AND. PRESENT(f25) ) CALL fill_farray( n25, f25, i )
    IF( PRESENT (n26) .AND. PRESENT(f26) ) CALL fill_farray( n26, f26, i )
    IF( PRESENT (n27) .AND. PRESENT(f27) ) CALL fill_farray( n27, f27, i )
    IF( PRESENT (n29) .AND. PRESENT(f29) ) CALL fill_farray( n29, f29, i )
    IF( PRESENT (n30) .AND. PRESENT(f30) ) CALL fill_farray( n30, f30, i )
    IF( PRESENT (n31) .AND. PRESENT(f31) ) CALL fill_farray( n31, f31, i )
    IF( PRESENT (n32) .AND. PRESENT(f32) ) CALL fill_farray( n32, f32, i )
    IF( PRESENT (n33) .AND. PRESENT(f33) ) CALL fill_farray( n33, f33, i )
    IF( PRESENT (n34) .AND. PRESENT(f34) ) CALL fill_farray( n34, f34, i )
    IF( PRESENT (n35) .AND. PRESENT(f35) ) CALL fill_farray( n35, f35, i )
    IF( PRESENT (n36) .AND. PRESENT(f36) ) CALL fill_farray( n36, f36, i )
    IF( PRESENT (n37) .AND. PRESENT(f37) ) CALL fill_farray( n37, f37, i )
    IF( PRESENT (n38) .AND. PRESENT(f38) ) CALL fill_farray( n38, f38, i )
    IF( PRESENT (n39) .AND. PRESENT(f39) ) CALL fill_farray( n39, f39, i )
    IF( PRESENT (n40) .AND. PRESENT(f40) ) CALL fill_farray( n40, f40, i )
    IF( PRESENT (n41) .AND. PRESENT(f41) ) CALL fill_farray( n41, f41, i )
    IF( PRESENT (n42) .AND. PRESENT(f42) ) CALL fill_farray( n42, f42, i )
    IF( PRESENT (n50) .AND. PRESENT(f50) ) CALL fill_farray( n50, f50, i )

    ALLOCATE( CheckNames(i) )
    ALLOCATE( CheckFields(i,mp) )

    DO k=1,i
       CheckNames(k) = farray_names(k)
       CheckFields(k,:) = farray_fields(k,:)
    ENDDO

    DEALLOCATE( farray_fields )

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

    INTEGER :: mp, np

    CHARACTER(len=*), OPTIONAL :: &
         n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, &
         n11, n12, n13, n14, n15, n16, n17, n18, n19, n20, &
         n21, n22, n23, n24, n25, n26, n27, n28, n29, n30, &
         n31, n32, n33, n34, n35, n36, n37, n38, n39, n40, &
         n41, n42, n43, n44, n45, n46, n47, n48, n49, n50

    REAL, DIMENSION(:,:), OPTIONAL :: &
         f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, &
         f11, f12, f13, f14, f15, f16, f17, f18, f19, f20,&
         f21, f22, f23, f24, f25, f26, f27, f28, f29, f30,&
         f31, f32, f33, f34, f35, f36, f37, f38, f39, f40,&
         f41, f42, f43, f44, f45, f46, f47, f48, f49, f50

    CHARACTER(len=30), DIMENSION(:), ALLOCATABLE :: CheckNames
    REAL, DIMENSION(:,:,:), ALLOCATABLE :: CheckFields

    INTEGER :: i, k

    ALLOCATE( farray_fields2( farray_nmax, mp, np) )
    i = 0

    IF( PRESENT (n1) .AND. PRESENT(f1) ) THEN
       CALL fill_farray2( n1, f1, i )
    ELSE
       PRINT *, "CABLE_log: cable_farray missing dummy args"
       RETURN
    ENDIF

    IF( PRESENT (n2) .AND. PRESENT(f2) ) CALL fill_farray2( n2, f2, i )
    IF( PRESENT (n3) .AND. PRESENT(f3) ) CALL fill_farray2( n3, f3, i )
    IF( PRESENT (n4) .AND. PRESENT(f4) ) CALL fill_farray2( n4, f4, i )
    IF( PRESENT (n5) .AND. PRESENT(f5) ) CALL fill_farray2( n5, f5, i )
    IF( PRESENT (n6) .AND. PRESENT(f6) ) CALL fill_farray2( n6, f6, i )
    IF( PRESENT (n7) .AND. PRESENT(f7) ) CALL fill_farray2( n7, f7, i )
    IF( PRESENT (n8) .AND. PRESENT(f8) ) CALL fill_farray2( n8, f8, i )
    IF( PRESENT (n9) .AND. PRESENT(f9) ) CALL fill_farray2( n9, f9, i )
    IF( PRESENT (n10) .AND. PRESENT(f10) ) CALL fill_farray2( n10, f10, i )
    IF( PRESENT (n11) .AND. PRESENT(f11) ) CALL fill_farray2( n11, f11, i )
    IF( PRESENT (n12) .AND. PRESENT(f12) ) CALL fill_farray2( n12, f12, i )
    IF( PRESENT (n13) .AND. PRESENT(f13) ) CALL fill_farray2( n13, f13, i )
    IF( PRESENT (n14) .AND. PRESENT(f14) ) CALL fill_farray2( n14, f14, i )
    IF( PRESENT (n15) .AND. PRESENT(f15) ) CALL fill_farray2( n15, f15, i )
    IF( PRESENT (n16) .AND. PRESENT(f16) ) CALL fill_farray2( n16, f16, i )
    IF( PRESENT (n17) .AND. PRESENT(f17) ) CALL fill_farray2( n17, f17, i )
    IF( PRESENT (n18) .AND. PRESENT(f18) ) CALL fill_farray2( n18, f18, i )
    IF( PRESENT (n19) .AND. PRESENT(f19) ) CALL fill_farray2( n19, f19, i )
    IF( PRESENT (n20) .AND. PRESENT(f20) ) CALL fill_farray2( n20, f20, i )
    IF( PRESENT (n21) .AND. PRESENT(f21) ) CALL fill_farray2( n21, f21, i )
    IF( PRESENT (n22) .AND. PRESENT(f22) ) CALL fill_farray2( n22, f22, i )
    IF( PRESENT (n23) .AND. PRESENT(f23) ) CALL fill_farray2( n23, f23, i )
    IF( PRESENT (n24) .AND. PRESENT(f24) ) CALL fill_farray2( n24, f24, i )
    IF( PRESENT (n25) .AND. PRESENT(f25) ) CALL fill_farray2( n25, f25, i )
    IF( PRESENT (n26) .AND. PRESENT(f26) ) CALL fill_farray2( n26, f26, i )
    IF( PRESENT (n27) .AND. PRESENT(f27) ) CALL fill_farray2( n27, f27, i )
    IF( PRESENT (n29) .AND. PRESENT(f29) ) CALL fill_farray2( n29, f29, i )
    IF( PRESENT (n30) .AND. PRESENT(f30) ) CALL fill_farray2( n30, f30, i )
    IF( PRESENT (n31) .AND. PRESENT(f31) ) CALL fill_farray2( n31, f31, i )
    IF( PRESENT (n32) .AND. PRESENT(f32) ) CALL fill_farray2( n32, f32, i )
    IF( PRESENT (n33) .AND. PRESENT(f33) ) CALL fill_farray2( n33, f33, i )
    IF( PRESENT (n34) .AND. PRESENT(f34) ) CALL fill_farray2( n34, f34, i )
    IF( PRESENT (n35) .AND. PRESENT(f35) ) CALL fill_farray2( n35, f35, i )
    IF( PRESENT (n36) .AND. PRESENT(f36) ) CALL fill_farray2( n36, f36, i )
    IF( PRESENT (n37) .AND. PRESENT(f37) ) CALL fill_farray2( n37, f37, i )
    IF( PRESENT (n38) .AND. PRESENT(f38) ) CALL fill_farray2( n38, f38, i )
    IF( PRESENT (n39) .AND. PRESENT(f39) ) CALL fill_farray2( n39, f39, i )
    IF( PRESENT (n40) .AND. PRESENT(f40) ) CALL fill_farray2( n40, f40, i )
    IF( PRESENT (n41) .AND. PRESENT(f41) ) CALL fill_farray2( n41, f41, i )
    IF( PRESENT (n42) .AND. PRESENT(f42) ) CALL fill_farray2( n42, f42, i )
    IF( PRESENT (n50) .AND. PRESENT(f50) ) CALL fill_farray2( n50, f50, i )

    ALLOCATE( CheckNames(i) )
    ALLOCATE( CheckFields(i,mp, np) )

    DO k=1,i
       CheckNames(k) = farray_names(k)
       CheckFields(k,:,:) = farray_fields2(k,:,:)
    ENDDO

    DEALLOCATE( farray_fields2 )

  END SUBROUTINE cable_farray2

  SUBROUTINE fill_farray( n, f, i )

    CHARACTER(len=*) :: n
    REAL, DIMENSION(:) :: f
    INTEGER :: i

    i=i+1
    farray_names(i) = n
    farray_fields(i,:) = f

  END SUBROUTINE fill_farray


  SUBROUTINE fill_farray2( n, f, i )

    CHARACTER(len=*) :: n
    REAL, DIMENSION(:,:) :: f
    INTEGER :: i

    i=i+1
    farray_names(i) = n
    farray_fields2(i,:,:) = f

  END SUBROUTINE fill_farray2


  SUBROUTINE cable_extremes1(fname,field,mype)

    REAL, DIMENSION(:,:) :: field
    CHARACTER(len=*), DIMENSION(:) :: fname

    INTEGER, OPTIONAL :: mype
    INTEGER :: i,j,k
    INTEGER :: n,m,op
    REAL :: emax, emin, emean, emode
    REAL :: erange
    REAL :: edbin
    REAL, DIMENSION(100) :: bin
    INTEGER, DIMENSION(100) :: ibin
    INTEGER :: ib, ibmax, binmax, maxbin

    n = SIZE(fname)
    m = SIZE(field,2)

    ! for each field in fname(i)
    DO i=1, n

       emax =  MAXVAL( field(i,:) )
       emin =  MINVAL(field(i,:) )
       emean =  SUM(field(i,:) ) / ( m )

       erange = emax - emin
       edbin = erange / 100. ! for 100 bins

       bin(1) = emin

       ! define bins per fname(i)
       DO ib=2, 100
          bin(ib) = bin(ib-1) + edbin
       ENDDO

       ibin =0
       ! for each Element in field
       DO j=1, m

          ! Assignn each Element to a bin
          DO ib=1, 99

             IF( field(i,j) >= bin(ib) .AND. &
                  field(i,j) < bin(ib+1) ) THEN

                !if(ib==1 )print *, "jhan:field1 ", field(i,j)
                ibin(ib) = ibin(ib) + 1

             ENDIF

          ENDDO ! DO LOOP over fill bins

       ENDDO ! DO LOOP over elements

       binmax = 0

       ! find max bin per field
       DO ib=1, 99

          IF( ibin(ib) > binmax ) THEN
             binmax = ibin(ib)
             maxbin = ib
          ENDIF

       ENDDO ! DO LOOP over bins

       PRINT *, "jhan:bins1 ", bin
       PRINT *, "jhan:bins1 count", ibin


       emode = bin(maxbin)

       PRINT *, ""
       PRINT *, "CABLE_log: "
       PRINT *, "   Field ", fname(i)
       PRINT *, "   Min ", emin
       PRINT *, "   Max ", emax
       PRINT *, "   Mean ",emean
       PRINT *, "   Mode ",emode
       PRINT *, "End CABLE_log: "
       PRINT *, ""

    ENDDO ! DO LOOP over fname(i)


  END SUBROUTINE cable_extremes1


  SUBROUTINE cable_extremes2(fname,field,mype)

    REAL, DIMENSION(:,:,:) :: field
    CHARACTER(len=*), DIMENSION(:) :: fname

    INTEGER, OPTIONAL :: mype
    INTEGER :: i,j,k
    INTEGER :: n,m,op
    REAL :: emax, emin, emean, emode
    REAL :: erange
    REAL :: edbin
    REAL, DIMENSION(100) :: bin
    INTEGER, DIMENSION(100) :: ibin
    INTEGER :: ib, ibmax, binmax, maxbin

    n = SIZE(fname)
    m = SIZE(field,2)
    op= SIZE(field,3)

    ! for each field in fname(i)
    DO i=1, n

       emax =  MAXVAL( field(i,:,:) )
       emin =  MINVAL(field(i,:,:) )
       emean =  SUM(field(i,:,:) ) / ( m*op )

       erange = emax - emin
       edbin = erange / 100. ! for 100 bins

       bin(1) = emin

       ! define bins per fname(i)
       DO ib=2, 100
          bin(ib) = bin(ib-1) + edbin
       ENDDO

       ibin =0
       ! for each Element in field
       DO j=1, m

          DO k=1, op

             ! Assignn each Element to a bin
             DO ib=1, 99

                IF( field(i,j,k) >= bin(ib) .AND. &
                     field(i,j,k) < bin(ib+1) ) THEN
                   !if(ib==1) print *, "jhan:field2 ", field(i,j,k)
                   ibin(ib) = ibin(ib) + 1

                ENDIF

             ENDDO ! DO LOOP over fill bins

          ENDDO ! DO LOOP over elements

       ENDDO ! DO LOOP over elements

       binmax = 0

       ! find max bin per field
       DO ib=1, 99

          IF( ibin(ib) > binmax ) THEN
             binmax = ibin(ib)
             maxbin = ib
          ENDIF

       ENDDO ! DO LOOP over bins

       emode = bin(maxbin)

       PRINT *, "jhan:bins2 ", bin
       PRINT *, "jhan:bins2 count", ibin

       PRINT *, ""
       PRINT *, "CABLE_log: "
       PRINT *, "   Field ", fname(i)
       PRINT *, "   Min ", emin
       PRINT *, "   Max ", emax
       PRINT *, "   Mean ",emean
       PRINT *, "   Mode ",emode
       PRINT *, "End CABLE_log: "
       PRINT *, ""

    ENDDO ! DO LOOP over fname(i)


  END SUBROUTINE cable_extremes2

#ifndef UM_BUILD
  SUBROUTINE def_dims(nd, ncid, dimID, dim_len, dim_name )
    USE netcdf
    IMPLICIT NONE
    INTEGER, INTENT(in) :: nd, ncid
    CHARACTER(len=*), DIMENSION(:), INTENT(in) :: dim_name
    INTEGER, DIMENSION(:), INTENT(out) :: dimID
    INTEGER, DIMENSION(:), INTENT(in) :: dim_len
    INTEGER :: j, ncok

    DO j=1, nd
       ncok = NF90_DEF_DIM(ncid, TRIM(dim_name(j)), dim_len(j), dimID(j) )
       IF (ncok /= nf90_noerr ) CALL stderr_nc(ncok,'def dim ', dim_name(j))
    ENDDO

    RETURN
  END SUBROUTINE def_dims




  SUBROUTINE def_vars(nv, ncid,  xtype, dimID, var_name,varID )
    USE netcdf
    IMPLICIT NONE
    INTEGER, INTENT(in) :: nv, ncid, xtype
    INTEGER, DIMENSION(:), INTENT(in) :: dimID
    INTEGER, DIMENSION(:), INTENT(inout) :: varID
    CHARACTER(len=*), DIMENSION(:), INTENT(in) :: var_name
    INTEGER :: j, ncok

    ! lat
    ncok = NF90_DEF_VAR( ncid, TRIM(var_name(1)), xtype, &
         (/ dimID(1) /), varID(1))
    IF (ncok /= nf90_noerr ) CALL stderr_nc(ncok,'def var ', var_name(1))

    ! lon
    ncok = NF90_DEF_VAR(ncid, TRIM(var_name(2)), xtype, &
         (/ dimID(1) /), varID(2))
    IF (ncok /= nf90_noerr ) CALL stderr_nc(ncok,'def var ', var_name(2))

    ! tairk
    ncok = NF90_DEF_VAR(ncid, TRIM(var_name(3)), xtype, &
         (/ dimID(1), dimID(3) /), varID(3))
    IF (ncok /= nf90_noerr ) CALL stderr_nc(ncok,'def var ', var_name(3))

    !tsoil
    ncok = NF90_DEF_VAR(ncid, TRIM(var_name(4)), xtype, &
         (/ dimID(1), dimID(2),dimID(3)/), varID(4))
    IF (ncok /= nf90_noerr ) CALL stderr_nc(ncok,'def var ', var_name(4))

    ! moist
    ncok = NF90_DEF_VAR(ncid, TRIM(var_name(5)), xtype, &
         (/ dimID(1), dimID(2),dimID(3)/), varID(5))
    IF (ncok /= nf90_noerr ) CALL stderr_nc(ncok,'def var ', var_name(5))

    !cgpp
    ncok = NF90_DEF_VAR(ncid, TRIM(var_name(6)), xtype, &
         (/ dimID(1), dimID(3)/), varID(6))
    IF (ncok /= nf90_noerr ) CALL stderr_nc(ncok,'def var ', var_name(6))

    !crmplant
    ncok = NF90_DEF_VAR(ncid, TRIM(var_name(7)), xtype, &
         (/ dimID(1), dimID(2),dimID(3)/), varID(7))
    IF (ncok /= nf90_noerr ) CALL stderr_nc(ncok,'def var ', var_name(7))

    !phenphase
    ncok = NF90_DEF_VAR(ncid, TRIM(var_name(8)), xtype, &
         (/ dimID(1), dimID(3)/), varID(8))
    IF (ncok /= nf90_noerr ) CALL stderr_nc(ncok,'def var ', var_name(8))

    !doyphase1
    ncok = NF90_DEF_VAR(ncid, TRIM(var_name(9)), xtype, &
         (/ dimID(1), dimID(3)/), varID(9))
    IF (ncok /= nf90_noerr ) CALL stderr_nc(ncok,'def var ', var_name(9))

    !doyphase2
    ncok = NF90_DEF_VAR(ncid, TRIM(var_name(10)), xtype, &
         (/ dimID(1), dimID(3)/), varID(10))
    IF (ncok /= nf90_noerr ) CALL stderr_nc(ncok,'def var ', var_name(10))

    !doyphase3
    ncok = NF90_DEF_VAR(ncid, TRIM(var_name(11)), xtype, &
         (/ dimID(1), dimID(3)/), varID(11))
    IF (ncok /= nf90_noerr ) CALL stderr_nc(ncok,'def var ', var_name(11))

    !doyphase4
    ncok = NF90_DEF_VAR(ncid, TRIM(var_name(12)), xtype, &
         (/ dimID(1), dimID(3)/), varID(12))
    IF (ncok /= nf90_noerr ) CALL stderr_nc(ncok,'def var ', var_name(12))


    !mtemp
    ncok = NF90_DEF_VAR(ncid, TRIM(var_name(13)), xtype, &
         (/ dimID(1),dimID(3)/), varID(13))
    IF (ncok /= nf90_noerr ) CALL stderr_nc(ncok,'def var ', var_name(13))

    !Ndep
    ncok = NF90_DEF_VAR(ncid, TRIM(var_name(14)), xtype, &
         (/ dimID(1),dimID(3)/), varID(14))
    IF (ncok /= nf90_noerr ) CALL stderr_nc(ncok,'def var ', var_name(14))

    RETURN
  END SUBROUTINE def_vars

  SUBROUTINE def_var_atts( ncfile_in, ncid, varID )
    USE netcdf
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in) :: ncfile_in
    INTEGER, INTENT(in):: ncid       ! netcdf file ID
    INTEGER, DIMENSION(:), INTENT(in) :: varID ! (1) ~ tvair, (2) ~ pmb
    INTEGER :: j, ncok
    CHARACTER(len=10) dummy

    WRITE(dummy,11) varID(1)
11  FORMAT(i2)
    ncok = NF90_PUT_ATT(ncid, nf90_global, "Title", "Forcing for define_air subroutine")
    IF (ncok /= nf90_noerr ) CALL stderr_nc(ncok,'def att ', ncfile_in)
    ncok = NF90_PUT_ATT(ncid, varID(3), "longname", "air temperature within canopy")
    IF (ncok /= nf90_noerr ) CALL stderr_nc(ncok,'def att ', dummy)
    ncok = NF90_PUT_ATT(ncid, varID(3), "units", "K")
    IF (ncok /= nf90_noerr ) CALL stderr_nc(ncok,'def att ', dummy)

    WRITE(dummy,11) varID(2)


    RETURN
  END SUBROUTINE def_var_atts


  SUBROUTINE put_var_ncr1(ncid, var_name, var )
    USE netcdf
    USE cable_def_types_mod, ONLY : mp
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in) ::  var_name
    REAL, DIMENSION(:),INTENT(in) :: var
    INTEGER, INTENT(in) :: ncid
    INTEGER :: ncok, varID,j

    ncok = NF90_INQ_VARID(ncid, var_name, varId )
    IF (ncok /= nf90_noerr ) CALL stderr_nc(ncok,'inquire var ', var_name)

    ncok = NF90_PUT_VAR(ncid, varId, var, start=(/1/), &
         count=(/mp/) )
    IF (ncok /= nf90_noerr ) CALL stderr_nc(ncok,'putting var ', var_name)

  END SUBROUTINE put_var_ncr1


  SUBROUTINE put_var_ncr2(ncid, var_name, var, n_call )
    USE netcdf
    USE cable_def_types_mod, ONLY : r_2, mp
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in) ::  var_name
    REAL(r_2), DIMENSION(:),INTENT(in) :: var
    INTEGER, INTENT(in) :: ncid, n_call
    INTEGER :: ncok, varID

    ncok = NF90_INQ_VARID(ncid, var_name, varId )
    IF (ncok /= nf90_noerr ) CALL stderr_nc(ncok,'inquire var ', var_name)

    ncok = NF90_PUT_VAR(ncid, varId, var, start=(/1,n_call /), &
         count=(/mp,1/) )

    IF (ncok /= nf90_noerr ) CALL stderr_nc(ncok,'putting var ', var_name)

  END SUBROUTINE put_var_ncr2

  !soil vars
  SUBROUTINE put_var_ncr3(ncid, var_name, var, n_call, nl)
    USE netcdf
    USE cable_def_types_mod, ONLY : r_2, mp, ms
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in) :: var_name
    REAL(r_2), DIMENSION(:,:),INTENT(in) :: var
    INTEGER, INTENT(in) :: ncid, n_call, nl
    INTEGER :: ncok, varID

    ncok = NF90_INQ_VARID( ncid, var_name, varId )
    IF( ncok /= nf90_noerr ) CALL stderr_nc(ncok,'inquire var ', var_name )

    ncok = NF90_PUT_VAR(ncid, varId, var, start=(/1,1,n_call /), &
         count=(/mp,nl,1/))
    IF (ncok /= nf90_noerr ) CALL stderr_nc(ncok,'putting var ', var_name)

    RETURN
  END SUBROUTINE put_var_ncr3



  SUBROUTINE get_var_ncr2(ncid, var_name, var, n_call )
    USE netcdf
    USE cable_def_types_mod, ONLY : r_2,mp
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in) :: var_name
    REAL(r_2), DIMENSION(:),INTENT(out) :: var
    INTEGER, INTENT(in) :: ncid
    INTEGER :: ncok, varID, n_call
    REAL, DIMENSION(mp) :: temp

    temp = 0.

    ncok = NF90_INQ_VARID(ncid, var_name, varId )
    IF (ncok /= nf90_noerr ) CALL stderr_nc(ncok,'inquire var ', var_name)
    ncok = NF90_GET_VAR(ncid, varId, temp, start=(/1,n_call/), &
         count=(/mp,1/) )

    IF (ncok /= nf90_noerr ) CALL stderr_nc(ncok,'getting var ', var_name)

    var = REAL( temp, r_2 )
  END SUBROUTINE get_var_ncr2

  SUBROUTINE get_var_ncr3(ncid, var_name, var, n_call, nl )
    USE netcdf
    USE cable_def_types_mod, ONLY : r_2, mp, ms
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in) :: var_name
    REAL(r_2), DIMENSION(:,:),INTENT(out) :: var
    INTEGER, INTENT(in) :: ncid, n_call, nl
    INTEGER :: ncok, varID
    REAL, DIMENSION(mp,1:nl) :: temp

    ncok = NF90_INQ_VARID(ncid, var_name, varId )
    IF (ncok /= nf90_noerr ) CALL stderr_nc(ncok,'inquire var ', var_name)

    ncok = NF90_GET_VAR(ncid, varId, temp, start=(/1,1,n_call /), &
         count=(/mp, nl, 1/))
    IF (ncok /= nf90_noerr ) CALL stderr_nc(ncok,'putting var ', var_name)
    var = REAL( temp, r_2 )
  END SUBROUTINE get_var_ncr3



  SUBROUTINE stderr_nc(status,message, var)
    USE netcdf
    CHARACTER(len=*), INTENT(in) :: message, var
    INTEGER, INTENT(IN) :: status
    CHARACTER(len=7) :: err_mess
    err_mess = 'ERROR:'
    PRINT *, (err_mess//message), var
    PRINT*,NF90_STRERROR(status)
    STOP
  END SUBROUTINE stderr_nc
#endif


END MODULE cable_diag_module
