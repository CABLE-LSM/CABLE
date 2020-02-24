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

#ifndef UM_BUILD
  interface put_var_nc
     module procedure put_var_ncr1, put_var_ncr2, put_var_ncr3
  end interface put_var_nc

  interface get_var_nc
     module procedure get_var_ncr2, get_var_ncr3
  end interface get_var_nc

#endif
CONTAINS

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


#ifndef UM_BUILD
  subroutine def_dims(nd, ncid, dimID, dim_len, dim_name )
    
    use netcdf

    implicit none

    integer,                        intent(in)  :: nd, ncid
    integer,          dimension(:), intent(out) :: dimID
    integer,          dimension(:), intent(in)  :: dim_len
    character(len=*), dimension(:), intent(in)  :: dim_name

    integer :: j, ncok

    do j=1, nd
       ncok = NF90_DEF_DIM(ncid, trim(dim_name(j)), dim_len(j), dimID(j))
       if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def dim ', dim_name(j))
    enddo

    return

  end subroutine def_dims


  subroutine def_vars(nv, ncid, xtype, dimID, var_name, varID)

    use netcdf
    
    implicit none
    
    integer,                        intent(in)    :: nv, ncid, xtype
    integer,          dimension(:), intent(in)    :: dimID
    character(len=*), dimension(:), intent(in)    :: var_name
    integer,          dimension(:), intent(inout) :: varID
    
    integer :: j, ncok

    ! lat
    ! print*, 'DV01'
    ncok = NF90_DEF_VAR( ncid, trim(var_name(1)), xtype, &
         (/ dimID(1) /), varID(1))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(1))

    ! lon
    ! print*, 'DV02'
    ncok = NF90_DEF_VAR(ncid, trim(var_name(2)), xtype, &
         (/ dimID(1) /), varID(2))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(2))

    ! tairk
    ! print*, 'DV03'
    ncok = NF90_DEF_VAR(ncid, trim(var_name(3)), xtype, &
         (/ dimID(1), dimID(3) /), varID(3))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(3))

    !tsoil
    ! print*, 'DV04'
    ncok = NF90_DEF_VAR(ncid, trim(var_name(4)), xtype, &
         (/ dimID(1), dimID(2), dimID(3)/), varID(4))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(4))

    ! moist
    ! print*, 'DV05'
    ncok = NF90_DEF_VAR(ncid, trim(var_name(5)), xtype, &
         (/ dimID(1), dimID(2), dimID(3)/), varID(5))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(5))

    !cgpp
    ! print*, 'DV06'
    ncok = NF90_DEF_VAR(ncid, trim(var_name(6)), xtype, &
         (/ dimID(1), dimID(3)/), varID(6))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(6))

    !crmplant
    ! print*, 'DV07'
    ncok = NF90_DEF_VAR(ncid, trim(var_name(7)), xtype, &
         (/ dimID(1), dimID(2), dimID(3)/), varID(7))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(7))

    !phenphase
    ! print*, 'DV08'
    ncok = NF90_DEF_VAR(ncid, trim(var_name(8)), xtype, &
         (/ dimID(1), dimID(3)/), varID(8))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(8))

    !doyphase1
    ! print*, 'DV09'
    ncok = NF90_DEF_VAR(ncid, trim(var_name(9)), xtype, &
         (/ dimID(1), dimID(3)/), varID(9))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(9))

    !doyphase2
    ! print*, 'DV10'
    ncok = NF90_DEF_VAR(ncid, trim(var_name(10)), xtype, &
         (/ dimID(1), dimID(3)/), varID(10))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(10))

    !doyphase3
    ! print*, 'DV11'
    ncok = NF90_DEF_VAR(ncid, trim(var_name(11)), xtype, &
         (/ dimID(1), dimID(3)/), varID(11))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(11))

    !doyphase4
    ! print*, 'DV12'
    ncok = NF90_DEF_VAR(ncid, trim(var_name(12)), xtype, &
         (/ dimID(1), dimID(3)/), varID(12))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(12))

    !mtemp
    ! print*, 'DV13'
    ncok = NF90_DEF_VAR(ncid, trim(var_name(13)), xtype, &
         (/ dimID(1), dimID(3)/), varID(13))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(13))

    !Ndep
    ! print*, 'DV14'
    ncok = NF90_DEF_VAR(ncid, trim(var_name(14)), xtype, &
         (/ dimID(1), dimID(3)/), varID(14))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(14))

    !Pdep
    ! print*, 'DV15'
    ncok = NF90_DEF_VAR(ncid, trim(var_name(15)), xtype, &
         (/ dimID(1), dimID(3)/), varID(15))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(15))

    !cAn12spin
    ! print*, 'DV16'
    ncok = NF90_DEF_VAR(ncid, trim(var_name(16)), xtype, &
         (/ dimID(1), dimID(3)/), varID(16))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(16))

    !cAn13spin
    ! print*, 'DV17'
    ncok = NF90_DEF_VAR(ncid, trim(var_name(17)), xtype, &
         (/ dimID(1), dimID(3)/), varID(17))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(17))

    return

  end subroutine def_vars

  
  subroutine def_var_atts( ncfile_in, ncid, varID )
    
    use netcdf
    
    implicit none
    
    character(len=*),      intent(in) :: ncfile_in
    integer,               intent(in) :: ncid       ! netcdf file ID
    integer, dimension(:), intent(in) :: varID      ! (1) ~ tvair, (2) ~ pmb

    integer :: j, ncok
    character(len=10) :: dummy

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


  subroutine put_var_ncr1(ncid, var_name, var)
    use netcdf
    use cable_def_types_mod, only : mp
    implicit none
    character(len=*), intent(in) ::  var_name
    real, dimension(:),intent(in) :: var
    integer, intent(in) :: ncid
    integer :: ncok, varID,j

    ncok = NF90_INQ_VARID(ncid, var_name, varId)
    if (ncok /= nf90_noerr) call stderr_nc(ncok, 'inquire var ', var_name)

    ncok = NF90_PUT_VAR(ncid, varId, var, start=(/1/), count=(/mp/))
    if (ncok /= nf90_noerr) call stderr_nc(ncok, 'putting var ', var_name)

  end subroutine put_var_ncr1


  subroutine put_var_ncr2(ncid, var_name, var, n_call)
    use netcdf
    use cable_def_types_mod, only : r_2, mp
    implicit none
    character(len=*), intent(in) ::  var_name
    real(r_2), dimension(:),intent(in) :: var
    integer, intent(in) :: ncid, n_call
    integer :: ncok, varID

    ncok = NF90_INQ_VARID(ncid, var_name, varId)
    if (ncok /= nf90_noerr) call stderr_nc(ncok, 'inquire var ', var_name)
   
    ncok = NF90_PUT_VAR(ncid, varId, var, start=(/1,n_call/), count=(/mp,1/))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok, 'putting var ', var_name)

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

    ncok = NF90_INQ_VARID( ncid, var_name, varId)
    IF (ncok /= nf90_noerr) call stderr_nc(ncok, 'inquire var ', var_name)

    ncok = NF90_PUT_VAR(ncid, varId, var, start=(/1,1,n_call/), count=(/mp,nl,1/))
    if (ncok /= nf90_noerr) call stderr_nc(ncok, 'putting var ', var_name)

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
!==========================================================================!
!--- cable generic print status
!==========================================================================!

SUBROUTINE cable_stat( routname)
   use cable_common_module, only : ktau_gl, knode_gl

   character(len=*), intent(in) :: routname
      if(knode_gl==1) &
         write(6,*) 'CABLE@  ', routname, ktau_gl

END SUBROUTINE cable_stat


END MODULE cable_diag_module



