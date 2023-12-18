!==============================================================================
! This source code is part of the 
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CABLE Academic User Licence Agreement 
! (the "Licence").
! You may not use this file except in compliance with the Licence.
! A copy of the Licence and registration form can be obtained from 
! http://www.cawcr.gov.au/projects/access/cable
! You need to register and read the Licence agreement before use.
! Please contact cable_help@nf.nci.org.au for any questions on 
! registration and the Licence.
!
! Unless required by applicable law or agreed to in writing, 
! software distributed under the Licence is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the Licence for the specific language governing permissions and 
! limitations under the Licence.
! ==============================================================================
!
! Purpose: Passes UM variables to CABLE, calls cbm, passes CABLE variables 
!          back to UM. 'Explicit' is the first of two routines that call cbm at 
!          different parts of the UM timestep.
!
! Called from: UM code sf_exch
!
! Contact: Jhan.Srbinovsky@csiro.au
!
! History: Developed for CABLE-JULES coupling in UM 10.5
!
!
! ==============================================================================

module cable_unpack_checks_mod

contains

subroutine cable_unpack_check(var_name, um_var, cable_var,  fLtpts, lpts, &
                              nt, mp, dt, row_length, rows, mype, latitude, &
                              longitude, tols, Lunpack, checks )

  USE cable_um_tech_mod,   ONLY : um1

  implicit none

  character(len=*) :: var_name 
  integer :: lpts, nt, mp, row_length, rows
  real  :: um_var(lpts,nt)
  real  :: cable_var(mp)
  logical :: fLtpts(lpts,nt)
  integer :: dt, mype 
  real, optional ::tols
  logical, optional :: Lunpack, checks
  REAL,  DIMENSION(row_length,rows) ::                             &
      latitude,   &
      longitude

  
  REAL :: cable_var_max, cable_var_min
  REAL :: cable_var_sum, cable_var_avg
  real  :: um_var_old(lpts,nt)
  real  :: um_var_new(lpts,nt)
  REAL :: um_var_sum, um_var_avg
  real  :: um_var_upper, um_var_lower
  REAL :: miss = 0.0
  integer :: funit
  integer :: i,j,k,l,n 
  integer :: fi,fj,fk 

  cable_var_max = maxval(cable_var)
  cable_var_min = minval(cable_var)
 
  cable_var_sum = 0.0 
  do i=1,mp
    cable_var_sum = cable_var_sum + cable_var(i)
  enddo
  cable_var_avg = cable_var_sum / mp
  
  funit = 6         
  write(funit,*) "CABLE_LSM:unpack_check::mype,dt,min,max,avg"
  write(funit,10) mype, dt, cable_var_max, cable_var_min, cable_var_avg
10 format("jh: ", I3.1, 1X, I2.1, 1X, 3F6.1 )     
  
  k=0 
  um_var_sum = 0.0 
  do j=1,nt
    do i=1,lpts
      if(Fltpts(i,j))then 
        um_var_sum = um_var_sum + um_var(i,j)
        k=k+1
      endif  
    enddo
  enddo
  um_var_avg = um_var_sum / k 
  
  um_var_old = um_var
  um_var_new = UNPACK(cable_var,  fLtpts, miss)

 
20 format(A10, ": old:new ", 2F6.1 )     
  k=0 
  !fj=0 
  do j=1,nt
    !fi=0 
    do i=1,lpts
      if(Fltpts(i,j))then 
        um_var_upper = um_var_old(i,j) + ( tols * um_var_old(i,j) )
        um_var_lower = um_var_old(i,j) - ( tols * um_var_old(i,j) )
        if( um_var_new(i,j) .gt. um_var_upper .OR. &
                um_var_new(i,j) .lt. um_var_lower) then
            k=k+1
            write(funit,20) var_name, um_var_old(i,j), um_var_new(i,j)  
         endif  
      endif  
    enddo
  enddo
 
  write(funit,*) "jh:This many bad points: ", k 
  
  fk=0 
DO N=1,um1%ntiles
  DO K=1,um1%TILE_PTS(N)
    L = um1%TILE_INDEX(K,N)
    J=(um1%LAND_INDEX(L)-1)/um1%row_length + 1
    I = um1%LAND_INDEX(L) - (J-1)*um1%row_length
    um_var_upper = um_var_old(l,n) + ( tols * um_var_old(l,n) )
    um_var_lower = um_var_old(l,n) - ( tols * um_var_old(l,n) )
    if( um_var_new(l,n) .gt. um_var_upper .OR. &
      um_var_new(l,n) .lt. um_var_lower) then
      fk=fk+1
      write(funit,20) var_name, um_var_old(l,n), um_var_new(l,n), &
                            latitude(i,j), longitude(i,j)      
    endif  
  enddo
enddo

  write(funit,*) "jh:This many fk bad points: ", fk 


!*U_S(I,J) = U_S(I,J)+FLAND(L)*um1%TILE_FRAC(L,N)*U_S_TILE(L,N)
!
  if(Lunpack) &
    um_var = UNPACK(cable_var,  fLtpts, miss)
  
  return

End subroutine cable_unpack_check

End module cable_unpack_checks_mod


