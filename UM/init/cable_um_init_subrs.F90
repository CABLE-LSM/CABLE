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
! Purpose: Routines to pass UM variables into appropriate CABLE variables and 
!          to map parameters for each surface type to CABLE arrays
!
! Contact: Jhan.Srbinovsky@csiro.au
!
! History: Rewrites of code in v1.8 (ACCESS1.3)
!
!
! ==============================================================================

MODULE cable_um_init_subrs_mod
   
   IMPLICIT NONE

CONTAINS

   subroutine initialize_maps(latitude,longitude, tile_index_mp)
      use cable_data_module, only : cable
      use cable_um_tech_mod, only : um1, veg
      use cable_def_types_mod, only : mp

      use cable_diag_module, only : cable_diag 
      use cable_common_module, only : ktau_gl, knode_gl, cable_user 
         
      implicit none
      real, intent(in), dimension(um1%row_length,um1%rows) :: &
         latitude, longitude
      integer, intent(in), dimension(um1%land_pts, um1%ntiles) :: &
         tile_index_mp  ! index of tile
          
      logical, save :: first_call = .true.
      
      INTEGER :: i, j

      real :: dlon
      real, dimension(um1%row_length) :: tlong,acoslong, acoslon1
      real, dimension(um1%row_length, um1%rows) :: new_longitude

      integer, save :: iDiag0, iDiag1, iDiag2, iDiag3, iDiag4, iDiag5 
      real, dimension(um1%row_length,um1%rows) :: &
         asinlatitude, acoslat, acoslon
      !decs to write text files mapping mp points to lat/lon, i/j 
      character(len=*), parameter :: hcomp="mype"
      character(len=*), parameter :: hcompa ="i"
      character(len=*), parameter :: hcompb ="j"
      character(len=*), parameter :: hcompc ="l"
      character(len=*), parameter :: hcompd ="n"
      character(len=*), parameter :: hcomp1 ="mp "
      character(len=*), parameter :: hcomp2 ="lat "
      character(len=*), parameter :: hcomp3 ="lon"
      character(len=*), parameter :: hcomp4 ="frac"
      character(len=*), parameter :: hcomp5 ="iveg"
      character(len=*), parameter :: footer1 =""
      character(len=*), parameter :: footer2 = &
                  "------------------------------------------------------------"
      character(len=*), parameter :: hfmt1 = &
                  '(A8, 4X, A8, 4X, A8, 4X, A8, 4X, A8, 4X, A8, 4X, A8,   4X, A8,   4X, A8,    4X, A8)'
      character(len=*), parameter :: dfmt1 = &
                  '(I8, 4X, I8, 4X, I8, 4X, I8, 4X, I8, 4X, I8, 4X, F8.3, 4X, F8.3, 4X, ES8.2, 4X, I8)'
      character(len=30) :: chnode
      character(len=12) :: filename
      character(len=9), parameter :: basename="cable_mp_"
      integer, dimension(um1%row_length,um1%rows) ::umi, umj 
      integer, dimension(um1%land_pts, um1%ntiles) :: uml,umn 
      integer, dimension(mp) :: cable_umi,cable_umj, cable_uml,cable_umn 
 
           
            allocate( cable%lat(mp), cable%lon(mp), cable%tile(mp), cable%tile_frac(mp) )

            !-------------------------------------   
            !---make indexes for tile, lat, lon
            !-------------------------------------   
           
            !=== LaTITUDE 
            !form acoslat(:,:) & acoslon(:,:)
            !acoslat = ( latitude ) /cable%const%math%pi180
            
            !--- get latitude index corresponding to cable points
            !call um2cable_rr( (asin(latitude)/cable%const%math%pi180), cable%lat )
            !call um2cable_rr( ((latitude)/cable%const%math%pi180), cable%lat )
            !call um2cable_rr( acoslat, cable%lat )
            call um2cable_rr( latitude, cable%lat )
           
            !==================================================================

            !=== LONGITUDE 
            !acoslong =  ( longitude(:,2) ) /cable%const%math%pi180  
            !!acoslon =  ( longitude ) /cable%const%math%pi180  
            !acoslon1 =  acoslon(:,1) 
            !!acoslong =  acos( longitude(:,1) ) /cable%const%math%pi180  
       
            !!--- get longitude index corresponding to cable points.
            !!--- this is not so straight forward as UM longitude index 
            !!--- contains ambiguity. thus define "new_longitude" first
            !tlong(1) = acoslong(1)
            !do j=2, um1%row_length
            !   if( acoslong(j) < acoslong(j-1) ) then  
            !      dlon = acoslong(j) - acoslong(j-1)
            !      tlong(j) = tlong(j-1) - dlon   
            !   else 
            !      tlong(j) = acoslong(j)
            !   endif           
            !enddo
            !
            !do j=1, um1%row_length
            !   new_longitude(j,:) = tlong(j)
            !enddo
            
            call um2cable_rr( longitude, cable%lon )
            
            !--- get tile index/fraction  corresponding to cable points
            cable%tile = pack(tile_index_mp, um1%l_tile_pts)
            cable%tile_frac = pack(um1%tile_frac, um1%l_tile_pts)

return          
            !--- write all these maps.  cable_user%initialize_mapping can be 
            !--- set in namelist cable.nml
            !if ( cable_user%initialize_mapping ) then
            !write indexes for tile, lat, lon  !fudge
            !asinlatitude = ( latitude ) /cable%const%math%pi180
             
            !call cable_diag( iDiag0, 'latitude', um1%rows, 1, ktau_gl,  & 
            !      knode_gl, 'latitude',asinlatitude(1,:)  ) 


            ! ----------------------------------------------------------------------------------
            write(chnode,10) knode_gl
   10       format(I3.3)   
            filename=trim(trim(basename)//trim(chnode))
            
            umi=0; umj=0; uml=0; umn=0            
            do i=1, um1%row_length      
               do j=1, um1%rows     
                 umi(i,j) = i 
                 umj(i,j) = j 
               enddo   
            enddo   
            
            do i=1, um1%land_pts
               do j=1, um1%ntiles
                 uml(i,j) = i 
                 umn(i,j) = j 
               enddo   
            enddo   

            call um2cable_irr( umi, cable_umi )
            call um2cable_irr( umj, cable_umj )
            
            cable_uml = pack(uml, um1%l_tile_pts)
            cable_umn = pack(umn, um1%l_tile_pts)
            
            !open(unit=12517,file=filename,status="unknown", &
            !      action="write", form="formatted",position='append' )
            !   
            !   write (12517, hfmt1) hcomp, hcompa, hcompb, hcompc, hcompd,hcomp1, hcomp2, hcomp3, hcomp4, hcomp5
            !   write (12517, *) footer2 
            !   do i=1, mp 
            !      WRITE(12517,dfmt1) , knode_gl, cable_umi(i), cable_umj(i), cable_uml(i), cable_umn(i), &
            !                        i, cable%lat(i), cable%lon(i),    &
            !                        cable%tile_frac(i), veg%iveg(i)  
            !   enddo   
            !   write (12517, *) footer1 
            !
            !close(12517)

            ! ----------------------------------------------------------------------------------
             
            call cable_diag( iDiag1, 'longitude', um1%row_length, 1, ktau_gl,  & 
                  knode_gl, 'longitude', ( new_longitude(:,1) ) ) 
        
            !write indexes for tile, lat, lon
            call cable_diag( iDiag2, 'lat_index', mp, 1, ktau_gl,  & 
                  knode_gl, 'lat', cable%lat )
            call cable_diag( iDiag3, 'lon_index', mp, 1, ktau_gl,  & 
                  knode_gl, 'lon', cable%lon )
            
            !this should be integer-ed. typecast for now
            call cable_diag( iDiag4, 'tile_index', mp, 1, ktau_gl,  & 
                  knode_gl, 'tile', real(cable%tile) )
            
            call cable_diag( iDiag5, 'tile_frac', mp, 1, ktau_gl,  & 
                  knode_gl, 'tile_frac', cable%tile_frac )
            
      
      return
   end subroutine initialize_maps
  
  
        
SUBROUTINE initialize_soil( bexp, hcon, satcon, sathh, smvcst, smvcwt,         &
                            smvccl, albsoil, tsoil_tile, sthu, sthu_tile,      &
                            dzsoil, slope_avg, slope_std,&
                             dz_gw,aq_perm,drain_dens ) 

   USE cable_def_types_mod, ONLY : ms, mstype, mp, r_2
   USE cable_um_tech_mod,   ONLY : um1, soil, veg, ssnow 
   USE cable_common_module, ONLY : cable_runtime, cable_user,                  &
                                   soilin, knode_gl, gw_params
   
   REAL, INTENT(IN), DIMENSION(um1%land_pts) :: &
      bexp, &
      hcon, &
      satcon, & 
      sathh, &
      smvcst, &
      smvcwt, &
      smvccl, &
      albsoil 
   
   REAL, INTENT(IN), DIMENSION(um1%land_pts) :: &
      slope_avg, &
      slope_std, &
      dz_gw,aq_perm,drain_dens

   REAL, INTENT(IN), DIMENSION(um1%land_pts, um1%sm_levels) :: sthu
   
   REAL, INTENT(IN), DIMENSION(um1%land_pts, um1%ntiles, um1%sm_levels) :: &
      sthu_tile,     &
      tsoil_tile

   REAL, INTENT(IN), DIMENSION(um1%sm_levels) :: dzsoil

   !___defs 1st call to CABLE in this run
   LOGICAL, SAVE :: first_call= .TRUE.
   INTEGER :: i,j,k,L,n
   REAL, ALLOCATABLE :: tempvar(:), tempvar2(:),fwork(:,:)
   LOGICAL, PARAMETER :: skip =.TRUE. 
   REAL, DIMENSION(mstype) :: dummy 
   REAL :: tmp_clay, tmp_sand
   REAL, ALLOCATABLE :: znode(:), ssat_bounded(:,:),rho_soil_bulk(:,:)

   REAL, PARAMETER :: snow_ccnsw = 2.0,&
                      ssat_lo = 0.15,&
                      ssat_hi = 0.65,&
                      rhob_lo = 810.0,&
                      rhob_hi = 2300.0

   REAL :: sucs_sign_factor,hyds_unit_factor,sucs_min_magnitude


   IF (cable_user%gw_model) THEN
      sucs_sign_factor = 1.0
      hyds_unit_factor = 1.0
      sucs_min_magnitude = 106.0
   ELSE
      sucs_sign_factor = 1.0
      hyds_unit_factor = 1.0/1000.0
      sucs_min_magnitude = 106.0/1000.0
   END IF
      
      dummy=0. 

      IF( first_call ) THEN 

         ssnow%pudsto = 0.0; ssnow%pudsmx = 0.0
      
         !--- soil%isoilm defines soiltype. 
         ! currently is either 2 (arbitrarily) or 9.
         ! type 9 -> permanent ice points which are dealt with by CABLE. 
         ! Spatially explicit soil properties are used by 
         ! the UM anyway, and is only really an issue for soil%css & 
         ! soil%rhosoil, which are set to either 2 or 9. 
         ! dealing with this in CASACNP is another issue.
         !--- %isoilm=10 for Lakes
         soil%isoilm  =  2
            
         ! set soil type for permanent ice based on where permanent ice 
         ! located in vegetation map (in v1.8 set by soil albedo value)
         ! hard-wired numbers to be removed in future release
         WHERE( veg%iveg == 17 ) soil%isoilm = 9
           
         !--- set CABLE-var soil%albsoil from UM var albsoil
         ! (see below ~ um2cable_lp)
         CALL um2cable_lp( albsoil, dummy, soil%albsoil(:,1),                &
                           soil%isoilm, skip )

         !--- defined in soil_thick.h in UM
         soil%zse = dzsoil
         soil%zse_vec = spread(dzsoil,1,mp)
         
         ! distance between consecutive layer midpoints
         soil%zshh(1)=0.5*soil%zse(1) 
         soil%zshh(ms+1)=0.5*soil%zse(ms)
         soil%zshh(2:ms) = 0.5 * (soil%zse(1:ms-1) + soil%zse(2:ms))

         !node depths
         IF (allocated(znode)) deallocate(znode)
         allocate(znode(ms))

         znode(1) = soil%zshh(1)
         do k=2,ms
            znode(k) = znode(k-1) * 0.5*(soil%zse(k-1)+soil%zse(k))
         end do

         !-------------------------------------------------------------------
         !--- UM met forcing vars needed by CABLE which have UM dimensions
         !---(land_pts,ntiles)[_lp], which is no good to cable. These have to be 
         !--- re-packed in a single vector of active tiles. Hence we use 
         !--- conditional "mask" l_tile_pts(land_pts,ntiles) which is .true.
         !--- if the land point is/has an active tile. generic format:
         !---     um2cable_lp( UM var, 
         !---                 default value for snow tile  where 
         !---                 peranent ice point to be treatedas a snowtile, 
         !---                 CABLE var, 
         !---                 mask )
         !--- where mask tells um2cable_lp whether or not to use default value 
         !--- for snow tile 
         !-------------------------------------------------------------------
         
         ! parameter b in Campbell equation 
         CALL um2cable_lp( BEXP, soilin%bch, soil%bch, soil%isoilm)
         
         ALLOCATE( tempvar(mstype), tempvar2(mp) )
   
         tempvar = soilin%sand(9) * 0.3  + soilin%clay(9) *0.25 +              &
                   soilin%silt(9) * 0.265
         
         CALL um2cable_lp( HCON, tempvar, tempvar2, soil%isoilm)
         soil%cnsd = REAL( tempvar2, r_2 )
         DEALLOCATE( tempvar, tempvar2 )
         
         ! hydraulic conductivity @saturation (satcon[mm/s], soilin%hyds[m/s] )
         CALL um2cable_lp( satcon,soilin%hyds*1000.0, soil%hyds, soil%isoilm)

         CALL um2cable_lp( sathh, soilin%sucs, soil%sucs, soil%isoilm)
         CALL um2cable_lp( smvcst, soilin%ssat, soil%ssat, soil%isoilm)
         CALL um2cable_lp( smvcwt, soilin%swilt, soil%swilt, soil%isoilm)
         CALL um2cable_lp( smvccl, soilin%sfc, soil%sfc, soil%isoilm)
   
    
         !mrd561
         if (allocated(fwork)) deallocate(fwork) 

         ALLOCATE( fwork(um1%land_pts,um1%ntiles) )

         fwork(:,:) = 20.0
         DO n=1,um1%NTILES
           do k=1,um1%TILE_PTS(N)
              i = um1%tile_index(k,n)
              fwork(i,n) = dz_gw(i)
            end do
         end do

         soil%GWdz(:) = pack(fwork(:,:),um1%l_tile_pts)

         do i=1,mp 
            if (soil%GWdz(i) .lt. 20.0) soil%GWdz(i) = 20.0
            if (soil%GWdz(i) .gt. 150.0) soil%GWdz(i) = 150.0
            if (veg%iveg(i) .eq. 16) soil%GWdz(i) = 150.0
         end do

         fwork(:,:) = 0.02
         DO n=1,um1%NTILES
           do k=1,um1%TILE_PTS(N)
              i = um1%tile_index(k,n)
              fwork(i,n) = slope_avg(i)
            end do
         end do
         soil%slope(:) = pack(fwork(:,:),um1%l_tile_pts)
         do i=1,mp 
            if (soil%slope(i) .lt. 0.0002) soil%slope(i) = 0.0002
            if (soil%slope(i) .gt. 0.2) soil%slope(i) = 0.2
            if (veg%iveg(i) .eq. 16) soil%slope(i) = 0.0002
         end do
 
         fwork(:,:) = .005 
         DO n=1,um1%NTILES
           do k=1,um1%TILE_PTS(N)
              i = um1%tile_index(k,n)
              fwork(i,n) = slope_std(i)
            end do
         end do
         soil%slope_std(:) = pack(fwork(:,:),um1%l_tile_pts) 
         do i=1,mp 
            if (soil%slope_std(i) .lt. 0.00002) soil%slope_std(i) = 0.00002
            if (soil%slope_std(i) .gt. 0.2) soil%slope_std(i) = 0.2
            if (veg%iveg(i) .eq. 16) soil%slope_std(i) = 0.00002
         end do
         fwork(:,:) = 0.0008
         soil%drain_dens(:) = 0.0
         DO n=1,um1%NTILES
           do k=1,um1%TILE_PTS(N)
              i = um1%tile_index(k,n)
              fwork(i,n) = drain_dens(i)
            end do
         end do
         soil%drain_dens(:) = pack(fwork(:,:),um1%l_tile_pts)
         do i=1,mp 
            if (soil%drain_dens(i) .lt. 1.0e-6) soil%drain_dens(i) = 1.0e-6
            if (soil%drain_dens(i) .gt. 0.02) soil%drain_dens(i) = 0.02
            if (veg%iveg(i) .eq. 16) soil%drain_dens(i) = 0.02
         end do

         WRITE(6,*) 'maxval soil%drain_dens',maxval(soil%drain_dens,dim=1)
         WRITE(6,*) 'minval soil%drain_dens',minval(soil%drain_dens,dim=1)
         IF (any(soil%drain_dens(:) .eq. 0.0)) &
                           write(*,*) 'drain_dens has values of zero'

 

         soil%GWhyds_vec(:) = 0.0
         fwork = 3.0e-6
         DO n=1,um1%NTILES
           do k=1,um1%TILE_PTS(N)
              i = um1%tile_index(k,n)
              fwork(i,n) = aq_perm(i)/10.0
            end do
         end do
         soil%GWhyds_vec(:) = pack(fwork(:,:),um1%l_tile_pts) 
         do i=1,mp 
            if (soil%GWhyds_vec(i) .lt. 1.0e-8) soil%GWhyds_vec(i) =1.0e-8
            if (soil%GWhyds_vec(i) .gt. 1.0e-3) soil%GWhyds_vec(i) =1.0e-3
         end do

         WRITE(6,*) 'maxval soil%GWhyds_vec',maxval(soil%GWhyds_vec,dim=1)
         WRITE(6,*) 'minval soil%GWhyds_vec',minval(soil%GWhyds_vec,dim=1)
         IF (any(soil%GWhyds_vec(:) .eq. 0.0)) &
                           write(*,*) 'hyds_vec has values of zero'

         deallocate(fwork) 
            
         !--- (re)set values for CABLE
         soil%ibp2    =  NINT(soil%bch)+2
         soil%i2bp3   =  2*NINT(soil%bch)+3
         
         ! satcon in UM is in mm/s; Cable needs m/s
         soil%hyds    = hyds_unit_factor * soil%hyds
         soil%sucs    = ABS( soil%sucs) * sucs_sign_factor
         soil%sucs    =  MAX(sucs_min_magnitude,soil%sucs)
         soil%ssat    =  MAX( soil%ssat, soil%sfc + 0.01 )

         
         !jhan:coupled runs 
         soil%hsbh    =  soil%hyds*ABS(soil%sucs)*soil%bch

         WHERE(soil%ssat > 0. )                                                &
            soil%pwb_min =  (soil%swilt / soil%ssat )**soil%ibp2
           
         !--- these are temporary 
         soil%rhosoil =  soilin%rhosoil(soil%isoilm)
         soil%css     =  soilin%css(soil%isoilm)

         do k=1,ms
            soil%ssat_vec(:,k)      = real(soil%ssat(:)   ,r_2)    
            soil%sucs_vec(:,k)      = real(soil%sucs(:)   ,r_2)   
            soil%hyds_vec(:,k)      = real(soil%hyds(:)   ,r_2)  
            soil%swilt_vec(:,k)     = real(soil%swilt(:)  ,r_2)  
            soil%bch_vec(:,k)       = real(soil%bch(:)    ,r_2)
            soil%sfc_vec(:,k)       = real(soil%sfc(:)    ,r_2)
            soil%rhosoil_vec(:,k)   = real(soil%rhosoil(:),r_2)   
            soil%cnsd_vec(:,k)      = real(soil%cnsd      ,r_2)
            soil%css_vec(:,k)       = real(soil%css       ,r_2)
            soil%watr(:,k)          = 0.001_r_2
         end do

         where (soil%ssat_vec .le. 0.0 .and. soil%sfc_vec .gt. 0.0)
              soil%ssat_vec = soil%sfc_vec + 0.05
         end where
         !--- Lestevens 28 Sept 2012 - Fix Init for soil% textures 
         !--- needed for CASA-CNP

         !default values, overwrite if cable_uer%gw_model selected
         soil%clay = soilin%clay(soil%isoilm)
         soil%silt = soilin%silt(soil%isoilm)
         soil%sand = soilin%sand(soil%isoilm)
         do k=1,ms
             !should read texture by layer evantually
              soil%clay_vec(:,k) = soil%clay(:)
              soil%sand_vec(:,k) = soil%sand(:)
              soil%silt_vec(:,k) = soil%silt(:)
         end do
         
         do k=1,ms
            do i=1,mp

               if ( (soil%silt_vec(i,k) .gt. 0.99) .or. &
                    (soil%silt_vec(i,k) .lt. 0.01) .or. &
                    (soil%sand_vec(i,k) .gt. 0.99) .or. &
                    (soil%sand_vec(i,k) .lt. 0.01) .or. &
                    (soil%clay_vec(i,k) .gt. 0.99) .or. &
                    (soil%clay_vec(i,k) .lt. 0.01) ) then

                    !all bad
                    soil%clay_vec(i,k) = 0.3
                    soil%sand_vec(i,k) = 0.3
                    soil%silt_vec(i,k) = 0.4

               end if

            end do

         end do


         IF (cable_user%gw_model) THEN
           
            DO k=1,ms

               do i=1,mp  !from reversing pedotransfer functions
                          !,ay cause io issues because not passed into um

                  if (soil%isoilm(i) .ne. 9) then

                        soil%hyds_vec(i,k) = soil%hyds_vec(i,k) * &   !change in hyds
                                            exp(-gw_params%hkrz*( znode(k)-gw_params%zdepth) )
  
                  end if


               end do
            end do

            k=1
            soil%hyds(:) = soil%hyds_vec(:,k)

         END IF

         IF (cable_user%soil_thermal_fix) then

           if (allocated(ssat_bounded)) deallocate(ssat_bounded)
           if (allocated(rho_soil_bulk)) deallocate(rho_soil_bulk)

           allocate(ssat_bounded(size(soil%ssat_vec,dim=1),&
                                 size(soil%ssat_vec,dim=2) ) )

           ssat_bounded(:,:) = min( ssat_hi, max(ssat_lo, &
                                              soil%ssat_vec(:,:) ) )

           allocate(rho_soil_bulk(size(soil%rhosoil_vec,dim=1),&
                                  size(soil%rhosoil_vec,dim=2) ) )

           rho_soil_bulk(:,:) = min(rhob_hi, max(rhob_lo , &
                                  (2700.0*(1.0 - ssat_bounded(:,:)) ) ) )


            do k=1,ms
               do i=1,mp


                  if (soil%isoilm(i) .ne. 9) then

                     soil%rhosoil_vec(i,k) = 2700.0

                     soil%cnsd_vec(i,k) = ( (0.135*(1.0-ssat_bounded(i,k))) +&
                                         (64.7/rho_soil_bulk(i,k)) ) / &
                                       (1.0 - 0.947*(1.0-ssat_bounded(i,k)))

                  end if

               end do
            end do

            k=1
            do i=1,mp
               if (soil%isoilm(i) .ne. 9) then
                  soil%rhosoil(i) = soil%rhosoil_vec(i,1)
                  soil%cnsd(i)    = soil%cnsd_vec(i,1)
               end if
            end do

           if (allocated(ssat_bounded)) deallocate(ssat_bounded)
           if (allocated(rho_soil_bulk)) deallocate(rho_soil_bulk)

         END IF

         !always set these though not needed unless gw_model - true
         !should read in the values but they need calibration
         !where (soil%GWssat_vec .lt. 0.0) &
         where(soil%isoilm .eq. 9 .or. veg%iveg .eq. 16)
         !CAN leave these read in, not enough testing for now
                soil%GWhyds_vec = soil%hyds_vec(:,ms)
         endwhere
         soil%GWssat_vec = soil%ssat_vec(:,ms)
         soil%GWsucs_vec = soil%sucs_vec(:,ms)
         soil%GWbch_vec  = soil%bch_vec(:,ms)
         soil%GWwatr     = 0.0

         
         !for sli   
         soil%nhorizons = 1 ! use 1 soil horizon globally
            
         first_call= .FALSE.
      ENDIF

   END SUBROUTINE initialize_soil
 
!========================================================================
!========================================================================
!========================================================================
          
SUBROUTINE initialize_veg( canht_ft, lai_ft, soil_zse ) 
   USE cable_um_tech_mod
   USE cable_common_module, ONLY : cable_runtime, cable_user, vegin
   
  REAL, INTENT(IN), DIMENSION(um1%land_pts, um1%npft) :: canht_ft, lai_ft 
  real, dimension(ms) :: soil_zse 
   
   LOGICAL, SAVE :: first_call= .TRUE. ! defs 1st call to CABLE in this run

      !---clobbers veg height, lai and resets ivegt for CABLE tiles
      CALL clobber_height_lai( canht_ft, lai_ft )
      
      !--- veg params were read from initialize_soil() 
      IF(first_call)  THEN
         CALL init_veg_pars_fr_vegin( soil_zse ) 
         ! Fix in-canopy turbulence scheme globally:
         veg%meth = 1
      ENDIF
      first_call= .FALSE.
     
END SUBROUTINE initialize_veg

!========================================================================
!========================================================================
!========================================================================

SUBROUTINE clobber_height_lai( um_htveg, um_lai )
   USE cable_um_tech_mod, ONLY : um1, kblum_veg, veg

   REAL, INTENT(IN), DIMENSION(um1%land_pts, um1%npft) ::                      &
                                                          um_htveg, um_lai
   INTEGER :: i,j,n
    
   DO N=1,um1%NTILES
      DO J=1,um1%TILE_PTS(N)
         
         i = um1%TILE_INDEX(j,N)  ! It must be landpt index

         IF( um1%TILE_FRAC(i,N) .gt. 0.0 ) THEN
            
            ! hard-wired vegetation type numbers need to be removed
            IF(N < 5 ) THEN ! rml changed 4 to 5
               ! trees
               kblum_veg%IVEGT(i,N) = N
               kblum_veg%LAIFT(i,N) = max(0.01,um_lai(i,N)) 
               kblum_veg%HTVEG(i,N) = max(1.,um_htveg(i,N)) 
            ELSE IF(N > 4 .AND. N < 14 ) THEN !rml changed 3 to 4
               ! shrubs/grass
               kblum_veg%IVEGT(i,N) = N
               kblum_veg%LAIFT(i,N) = max(0.01, um_lai(i,N)) 
               kblum_veg%HTVEG(i,N) = max(0.1, um_htveg(i,N)) 
             ELSE IF(N > 13 ) THEN
               ! non-vegetated
               kblum_veg%IVEGT(i,N) = N
               kblum_veg%LAIFT(i,N) = 0. 
               kblum_veg%HTVEG(i,N) = 0.
            ENDIF

         ENDIF

      ENDDO
   ENDDO
  
   veg%iveg   = PACK(kblum_veg%ivegt, um1%L_TILE_PTS)
   veg%vlai   = PACK(kblum_veg%laift, um1%L_TILE_PTS)
   veg%hc     = PACK(kblum_veg%htveg, um1%L_TILE_PTS)

END SUBROUTINE clobber_height_lai

!========================================================================
!========================================================================
!========================================================================

SUBROUTINE init_respiration(NPP_FT_ACC,RESP_W_FT_ACC)
   ! Lestevens 23apr13 - for reading in prog soil & plant resp
   USE cable_um_tech_mod,   ONLY : um1, canopy
   !USE cable_common_module, ONLY : cable_runtime, cable_user

   REAL, INTENT(INOUT),DIMENSION(um1%land_pts, um1%ntiles) :: NPP_FT_ACC
   REAL, INTENT(INOUT),DIMENSION(um1%land_pts, um1%ntiles) :: RESP_W_FT_ACC

!   REAL, ALLOCATABLE :: tempvar(:,:), tempvar2(:,:)
   INTEGER :: l,j,n

!   ALLOCATE( tempvar(um1%land_pts,um1%ntiles) )
!   ALLOCATE( tempvar2(um1%land_pts,um1%ntiles) )
!
!      DO N=1,um1%NTILES
!         DO J=1,um1%TILE_PTS(N)
!
!            L = um1%TILE_INDEX(j,N)  ! It must be landpt index
!
!            !IF( um1%TILE_FRAC(L,N) .gt. 0.0 ) THEN
!
!               !IF(N <= 13 ) THEN
!                  tempvar(L,N)  = NPP_FT_ACC(L,N)
!                  tempvar2(L,N) = RESP_W_FT_ACC(L,N)
!               !ELSE IF(N > 13 ) THEN
!               !   tempvar(L,N)  = 0.
!               !   tempvar2(L,N) = 0.
!               !ENDIF
!
!            !ENDIF
!
!         ENDDO
!      ENDDO

      !---set soil & plant respiration (now in dim(land_pts,ntiles))
      canopy%frs = PACK(NPP_FT_ACC   , um1%L_TILE_PTS)
      canopy%frp = PACK(RESP_W_FT_ACC, um1%L_TILE_PTS)
!      canopy%frs = PACK(tempvar, um1%L_TILE_PTS)
!      canopy%frp = PACK(tempvar2, um1%L_TILE_PTS)

      !---convert units to g C m-2 s-1
      canopy%frs = canopy%frs * 1000.
      canopy%frp = canopy%frp * 1000.

!   DEALLOCATE( tempvar, tempvar2 )

END SUBROUTINE init_respiration

!========================================================================
!========================================================================
!========================================================================

SUBROUTINE init_veg_pars_fr_vegin( soil_zse ) 
  USE cable_common_module, ONLY : vegin, init_veg_from_vegin,cable_user, &
                                  knode_gl, ktau_gl
  USE cable_um_tech_mod,   ONLY : veg, soil 
  USE cable_def_types_mod, ONLY : mp,ms

  real, dimension(ms) :: soil_zse 

  CALL init_veg_from_vegin(1, mp, veg, soil_zse) 

  !For Legacy sake - ACCESS1.3 froot distribution WAS fixed for all veg types
  IF (cable_user%access13roots) THEN
    veg%froot(:,1) = 0.05
    veg%froot(:,2) = 0.20
    veg%froot(:,3) = 0.20
    veg%froot(:,4) = 0.20
    veg%froot(:,5) = 0.20
    veg%froot(:,6) = 0.15
  ENDIF

  veg%ejmax    = 2.*veg%vcmax
  
  !offline set init _parameters
  veg%gamma = 3.e-2 !for Haverd2013 switch 
  veg%F10 = 0.85 

END SUBROUTINE init_veg_pars_fr_vegin

!========================================================================
!========================================================================
!========================================================================
        
SUBROUTINE initialize_radiation( sw_down, lw_down, cos_zenith_angle,           &
                                 surf_down_sw, sin_theta_latitude, ls_rain,    &
                                 ls_snow, tl_1, qw_1, vshr_land, pstar,        &
! rml 2/7/13 pass 3d co2 through to cable if required
                   CO2_MMR,CO2_3D,CO2_DIM_LEN,CO2_DIM_ROW,L_CO2_INTERACTIVE )   

   USE cable_def_types_mod, ONLY : mp
   USE cable_data_module,   ONLY : PHYS, OTHER
   USE cable_um_tech_mod,   ONLY : um1, rad, soil, met,                        &
                                   conv_rain_prevstep, conv_snow_prevstep
   USE cable_common_module, ONLY : cable_runtime, cable_user, ktau_gl, kwidth_gl

   REAL, INTENT(INOUT), DIMENSION(um1%row_length, um1%rows) :: sw_down
   
   REAL, INTENT(IN), DIMENSION(um1%row_length, um1%rows) ::                    & 
      lw_down,           &
      sin_theta_latitude
   
   REAL, INTENT(INOUT), DIMENSION(um1%row_length, um1%rows) :: cos_zenith_angle

   REAL, INTENT(IN), DIMENSION(um1%row_length, um1%rows, 4) :: surf_down_sw 
   
   REAL, INTENT(IN), DIMENSION(um1%row_length, um1%rows) ::                    & 
      ls_rain,    &
      ls_snow,    &   
      tl_1,       &
      qw_1,       &
      vshr_land,  & 
      pstar
   
   REAL, INTENT(IN) :: co2_mmr
! rml 2/7/13 Extra atmospheric co2 variables
   LOGICAL, INTENT(IN) :: L_CO2_INTERACTIVE
   INTEGER, INTENT(IN) ::                              &
      CO2_DIM_LEN                                      &
     ,CO2_DIM_ROW
   REAL, INTENT(IN) :: CO2_3D(:,:)  ! co2 mass mixing ratio
             
   !___defs 1st call to CABLE in this run. OK in UM & coupled
   LOGICAL, SAVE :: first_call= .TRUE.
     
   REAL, POINTER :: TFRZ, RAD_THRESH
      
      TFRZ => PHYS%TFRZ
      RAD_THRESH => OTHER%RAD_THRESH
     
      IF( first_call ) THEN
         rad%albedo_T = soil%albsoil(:,1)
         !first_call = .FALSE.            !second use of first_call later
         ALLOCATE( conv_rain_prevstep(mp), conv_snow_prevstep(mp) )
         conv_rain_prevstep = 0. 
         conv_snow_prevstep = 0.
      ENDIF   
      
      ! re-set UM rad. forcings to suit CABLE. also called in explicit call to 
      ! CABLE from subr cable_um_expl_update() 
      CALL update_kblum_radiation( sw_down, cos_zenith_angle, surf_down_sw )
      
      ! set met. and rad. forcings to CABLE. also called in radiation call to 
      ! CABLE from subr cable_rad_() !jhan?
      ! subr.  um2cable_met_rad_alb() USES CABLE types met%, rad%, soil%
      ! and kblum% rad. calculated in  update_kblum_radiation() above 
      CALL um2cable_met_rad( cos_zenith_angle)
         
      ! UM met forcing vars needed by CABLE which have UM dimensions
      !(row_length,rows)[_rr], which is no good to cable. These have to be 
      ! re-packed in a single vector of active tiles. Hence we use 
      ! conditional "mask" l_tile_pts(land_pts,ntiles) which is .true.
      ! if the land point is/has an active tile
      ! generic format:
      ! um2cable_rr( UM var, CABLE var)
      
      ! CABLE met type forcings, not set by um2cable_met_rad()
      CALL um2cable_rr( LW_DOWN, met%fld)
      CALL um2cable_rr( (LS_RAIN*kwidth_gl), met%precip)
      CALL um2cable_rr( (LS_SNOW*kwidth_gl), met%precip_sn)
      CALL um2cable_rr( TL_1, met%tk)

      CALL um2cable_rr( QW_1, met%qv)
      CALL um2cable_rr( VSHR_LAND, met%ua)
      CALL um2cable_rr( PSTAR*0.01, met%pmb)
      
      !---re-set some of CABLE's forcing variables
      met%precip   =  met%precip + met%precip_sn 
      !met%precip   =  (met%precip + conv_rain_prevstep) &
      !               + (met%precip_sn +  conv_snow_prevstep)
      !               + (met%precip_sn +  conv_rain_prevstep)
      met%tvair =     met%tk
      met%tvrad =     met%tk
      met%coszen =    max(met%coszen,1e-8)

      !initialise rad%trad on first call only
      IF (first_call) THEN
         rad%trad = met%tk
         first_call = .FALSE.
      END IF 

      !---this is necessary clobrring at present 
      WHERE(met%ua < 0.001 ) met%ua = 0.001
      
      ! rml 24/2/11 Set atmospheric CO2 seen by cable to CO2_MMR (value seen 
      ! by radiation scheme).  Option in future to have cable see interactive 
      ! (3d) CO2 field Convert CO2 from kg/kg to mol/mol ( m_air, 
      ! 28.966 taken from include/constant/ccarbon.h file )
      ! r935 rml 2/7/13 Add in co2_interactive option
      !IF (L_CO2_INTERACTIVE) THEN
      !  CALL um2cable_rr(CO2_3D, met%ca)
      !ELSE
        met%ca = CO2_MMR
      !ENDIF
      met%ca = met%ca * 28.966/44.

      WHERE (met%coszen < RAD_THRESH ) 
         rad%fbeam(:,1) = REAL(0) 
         rad%fbeam(:,2) = REAL(0) 
         rad%fbeam(:,3) = REAL(0) 
      ENDWHERE

      !--- CABLE radiation type forcings, not set by um2cable_met_rad(
      !--- kblum_rad% vars are computed in subroutine update_kblum_radiation 
      CALL um2cable_rr( um1%longitude,rad%longitude )

END SUBROUTINE initialize_radiation

!========================================================================
!========================================================================
!========================================================================
          
SUBROUTINE initialize_canopy(canopy_tile,visc_sublayer_dz)
   USE cable_um_tech_mod,   ONLY : um1, canopy 
   USE cable_common_module, ONLY : cable_runtime, cable_user
   
   REAL, INTENT(IN),DIMENSION(um1%land_pts, um1%ntiles) :: canopy_tile
   REAL, INTENT(IN),DIMENSION(um1%land_pts, um1%ntiles) :: visc_sublayer_dz
   
   ! defs 1st call to CABLE in this run. OK in UM & coupled
   LOGICAL, SAVE :: first_call= .TRUE.
   
      !--- %ga is computed (on LHS only) in define_canopy and 
      !--- then used in soilsnow() in implicit call, then unpacked
      IF( first_call ) THEN
         canopy%ga = 0.
         canopy%us = 0.01
         canopy%fes_cor = 0.
         canopy%fhs_cor = 0.
         canopy%fwsoil = 1.
         first_call = .FALSE.
      ENDIF

     !---set canopy storage (already in dim(land_pts,ntiles) ) 
     canopy%cansto = pack(CANOPY_TILE, um1%l_tile_pts)
     canopy%oldcansto=canopy%cansto
     canopy%sublayer_dz(:) = 0.0  !junk
     canopy%sublayer_dz(:) = pack(visc_sublayer_dz(:,:),um1%l_tile_pts) 
     where (canopy%sublayer_dz .lt. 1.0e-8) canopy%sublayer_dz = 1.0e-8
     where (canopy%sublayer_dz .gt. 1.0) canopy%sublayer_dz = 1.0

     IF (first_call ) THEN

        WRITE(6,*) 'maxval canopy%sublayer_dz',maxval(canopy%sublayer_dz,dim=1)
        WRITE(6,*) 'minval canopy%sublayer_dz',minval(canopy%sublayer_dz,dim=1)

     END IF

END SUBROUTINE initialize_canopy

!========================================================================
!========================================================================
!========================================================================
 
SUBROUTINE initialize_soilsnow( smvcst, tsoil_tile, sthf_tile,smcl_tile,smgw_tile, &
                                snow_tile, snow_rho1l, snow_age, isnow_flg3l,&
                                snow_rho3l, snow_cond, snow_depth3l,           &
                                snow_mass3l, snow_tmp3l, fland,                &
                                sin_theta_latitude ) 

   USE cable_def_types_mod,  ONLY : mp, msn, ms, r_2,mstype
   USE cable_data_module,   ONLY : PHYS
   USE cable_um_tech_mod,   ONLY : um1, soil, ssnow, met, bal, veg
   USE cable_common_module, ONLY : cable_runtime, cable_user, ktau_gl
   
   REAL, INTENT(IN), DIMENSION(um1%land_pts) :: smvcst
   
   REAL, INTENT(IN), DIMENSION(um1%land_pts, um1%ntiles, um1%sm_levels) ::    &
      sthf_tile, &   !
      smcl_tile, &   !
      tsoil_tile     !

   REAL, INTENT(IN), DIMENSION(um1%land_pts, um1%ntiles) ::    &
      smgw_tile

   INTEGER, INTENT(IN), DIMENSION(um1%land_pts, um1%ntiles) :: isnow_flg3l 

   REAL, INTENT(INOUT), DIMENSION(um1%land_pts, um1%ntiles) :: snow_tile

   REAL, INTENT(IN), DIMENSION(um1%land_pts, um1%ntiles) ::                    &
      snow_rho1l, &  !
      snow_age     !

   REAL, INTENT(INOUT), DIMENSION(um1%land_pts, um1%ntiles,3) :: snow_cond

   REAL, INTENT(IN), DIMENSION(um1%land_pts, um1%ntiles,3) ::                  & 
      snow_rho3l,    & !
      snow_depth3l,  & !
      snow_mass3l,   & !
      snow_tmp3l       !
   
   REAL, INTENT(IN), DIMENSION(um1%land_pts) :: fland 
   
   REAL, INTENT(IN), DIMENSION(um1%row_length, um1%rows) :: sin_theta_latitude
   
   INTEGER :: i,j,k,L,n
   REAL  :: zsetot, max_snow_depth=50000.
   REAL, ALLOCATABLE:: fwork(:,:,:), sfact(:), fvar(:), rtemp(:),&
                       tot_mass_tmp(:,:,:), ice_vol_tmp(:,:,:)
   REAL, POINTER :: TFRZ
   LOGICAL :: skip =.TRUE. 
   LOGICAL, save :: first_call = .TRUE.
   REAL, DIMENSION(mstype) :: dummy 
      
      dummy=0. 

!     not sure if this is in restart file hence repeated again
      IF( first_call) THEN 
        ssnow%pudsto = 0.0 
      endif
      ssnow%pudsmx = 0.0
      ssnow%wbtot1 = 0
      ssnow%wbtot2 = 0
      ssnow%wb_lake = 0.

      !local pntrs to derived type data
      TFRZ => PHYS%TFRZ
       !why is snow updated from um values every timestep
       !but soil moisture not?
      
      !line removed Jun 2018 to assist with water conservation
      !in coupled runs alongside the ice-berg scheme
      !snow_tile = MIN(max_snow_depth, snow_tile)

      ssnow%snowd  = PACK(SNOW_TILE,um1%l_tile_pts)
      ssnow%ssdnn  = PACK(SNOW_RHO1L,um1%l_tile_pts)  
      ssnow%isflag = PACK(int(ISNOW_FLG3L),um1%l_tile_pts)  
!jhan: clobber
!ssnow%isflag = 0. 
      
      DO J=1, msn
         
         ssnow%sdepth(:,J)= PACK(SNOW_DEPTH3L(:,:,J),um1%l_tile_pts)
         ssnow%smass(:,J) = PACK(SNOW_MASS3L(:,:,J),um1%l_tile_pts)  
         ssnow%ssdn(:,J)  = PACK(SNOW_RHO3L(:,:,J),um1%l_tile_pts)  
         ssnow%tggsn(:,J) = PACK(SNOW_TMP3L(:,:,J),um1%l_tile_pts)  
         ssnow%sconds(:,J)= PACK(SNOW_COND(:,:,J),um1%l_tile_pts)  
         
      ENDDO 
      !ssnow%wb_lake = MAX( ssnow%wbtot2 - ssnow%wbtot1, 0.)
       
      DO J=1,um1%sm_levels
         ssnow%tgg(:,J) = PACK(TSOIL_TILE(:,:,J),um1%l_tile_pts)
      ENDDO 
      ssnow%snage = PACK(SNOW_AGE, um1%l_tile_pts)

      ssnow%GWwb(:) = pack(smgw_tile(:,:),um1%l_tile_pts)
      where (ssnow%GWwb .gt. soil%GWssat_vec) ssnow%GWwb = soil%GWssat_vec

      IF( first_call) THEN 
        
         ssnow%wbtot = 0.
         ssnow%wb_lake = 0.0
         ssnow%totwblake = 0.0  ! wb_lake integrated over river timestep
         ssnow%tggav = 0.
         ssnow%rtsoil = 50.
         ssnow%t_snwlr = 0.05

         ! snow depth from prev timestep 
         ssnow%osnowd  = PACK(SNOW_TILE,um1%l_tile_pts)  

         zsetot = sum(soil%zse)
         DO k = 1, um1%sm_levels
            ssnow%tggav = ssnow%tggav  + soil%zse(k)*ssnow%tgg(:,k)/zsetot
         END DO
     
         ! not updated 
         ALLOCATE( sfact( mp ) )
         sfact = 0.68
         WHERE (soil%albsoil(:,1) <= 0.14) 
            sfact = 0.5
         ELSEWHERE (soil%albsoil(:,1) > 0.14 .and. soil%albsoil(:,1) <= 0.20)
           sfact = 0.62
         END WHERE
         ssnow%albsoilsn(:,2) = 2. * soil%albsoil(:,1) / (1. + sfact)
         ssnow%albsoilsn(:,1) = sfact * ssnow%albsoilsn(:,2)
         DEALLOCATE( sfact )

         ALLOCATE( fvar(um1%land_pts ) )
         DO N=1,um1%NTILES
            DO K=1,um1%TILE_PTS(N)
               L = um1%TILE_INDEX(K,N)
               fvar(L) = real(L)
            ENDDO
         ENDDO
         CALL um2cable_lp( fland,dummy , ssnow%fland, soil%isoilm, skip )
         CALL um2cable_lp( fvar, dummy, ssnow%ifland, soil%isoilm, skip )
         DEALLOCATE( fvar )
         
         !--- updated via smcl,sthf etc 
         !previous code required rho_water == rho_ice
         !GWmodel correctly uses rho_ice ~= 0.92*rho_water
         !code now handles when densities are ice==water and ice!=water
         IF (ALLOCATED(tot_mass_tmp)) DEALLOCATE(tot_mass_tmp)
         IF (ALLOCATED(ice_vol_tmp)) DEALLOCATE(ice_vol_tmp)

         ALLOCATE( ice_vol_tmp(um1%land_pts,um1%ntiles,um1%sm_levels) )
         ALLOCATE( tot_mass_tmp(um1%land_pts,um1%ntiles,um1%sm_levels) )

         ice_vol_tmp(:,:,:) = 0.
         tot_mass_tmp(:,:,:) = 0.

         DO N=1,um1%NTILES                                                   
           DO K=1,um1%TILE_PTS(N)                                           
           I = um1%TILE_INDEX(K,N)                                      
             DO J = 1,um1%SM_LEVELS
               tot_mass_tmp(I,N,J) = SMCL_TILE(I,N,J)
               ice_vol_tmp(I,N,J) = STHF_TILE(I,N,J)*SMVCST(I)
             ENDDO ! J
           ENDDO
         ENDDO
   
         DO J = 1,um1%SM_LEVELS
            !ice volume
            ssnow%wbice(:,J) = pack(ice_vol_tmp(:,:,J),um1%l_tile_pts)
            ssnow%wbice(:,J) = max(0.,ssnow%wbice(:,J))  !should not be needed -- mrd561
            !liq volume  from (tot_mass - ice_mass) / (dz*rho_liq)
            ssnow%wbliq(:,j)= (pack(tot_mass_tmp(:,:,J),um1%l_tile_pts) -  &!total mass
                               ssnow%wbice(:,J)*soil%zse(j)*um1%RHO_ICE &!subtract ice mass 
                               )/(soil%zse(j)*um1%RHO_WATER)            !convert units
            ssnow%wb(:,J)    = ssnow%wbice(:,J) + ssnow%wbliq(:,j) 
            ! lakes: removed hard-wired number in future version
            !WHERE( veg%iveg == 16 ) ssnow%wb(:,J) = 0.95*soil%ssat
            !WHERE( veg%iveg == 16 ) ssnow%wb(:,J) = soil%sfc
         ENDDO
         
         DEALLOCATE( tot_mass_tmp )
         DEALLOCATE( ice_vol_tmp )

         
         ssnow%owetfac = MAX( 0., MIN( 1.0,                                    &
                         ( ssnow%wb(:,1) - soil%swilt ) /                      &
                         ( max(0.083, (soil%sfc - soil%swilt) ) )              &
                         ) )

         ! Temporay fix for accounting for reduction of soil evaporation 
         ! due to freezing
         WHERE( ssnow%wbice(:,1) > 0. )                                        &
            ! Prevents divide by zero at glaciated points where both 
            ! wb and wbice=0.
            ssnow%owetfac = ssnow%owetfac * ( 1.0 - ssnow%wbice(:,1) /         &
                            ssnow%wb(:,1) )**2
      
         !jhan: do we want to do this before %owetfac is set 
         DO J = 1, um1%sm_levels
            !should be removed!!!!!!!! This cannot conserve if there are any
            !dynamics 
            WHERE( soil%isoilm == 9 ) ! permanent ice: remove hard-wired no. in future
               ssnow%wb(:,J) = 0.95*soil%ssat
               ssnow%wbice(:,J) = 0.85*ssnow%wb(:,J)
            ENDWHERE
            !no not force rho_water==rho_ice==1000.0
            ssnow%wbtot = ssnow%wbtot + soil%zse(j)*&
                                       (ssnow%wbliq(:,j)*um1%RHO_WATER+&
                                        ssnow%wbice(:,j)*um1%RHO_ICE )                     
         ENDDO
         IF (cable_user%gw_model) THEN
            ssnow%wbtot = ssnow%wbtot + ssnow%GWwb(:)*soil%GWdz*um1%RHO_WATER
         ENDIF
     
         bal%wbtot0 = ssnow%wbtot

         !---set antartic flag using  sin_theta_latitude(row_length,rows)
         ALLOCATE( fwork(1,um1%land_pts,um1%ntiles) )
         fwork = 0.0
         DO N=1,um1%NTILES                     
            DO K=1,um1%TILE_PTS(N)
               L = um1%TILE_INDEX(K,N)
               J=(um1%LAND_INDEX(L)-1)/um1%row_length + 1
               I = um1%LAND_INDEX(L) - (J-1)*um1%row_length
               IF( sin_theta_latitude(I,J) .LT. -0.91 ) fwork(1,L,N) = 1.0
            ENDDO
         ENDDO
         ssnow%iantrct = pack(fwork(1,:,:),um1%L_TILE_PTS)
        
         DEALLOCATE( fwork )

      !! SLI specific initialisations:
      !IF(cable_user%SOIL_STRUC=='sli') THEN
      !   ssnow%h0(:)        = 0.0
      !   ssnow%S(:,:)       = ssnow%wb(:,:)/SPREAD(soil%ssat,2,ms)
      !   ssnow%snowliq(:,:) = 0.0
      !   ssnow%Tsurface     = 25.0
      !   ssnow%nsnow        = 0
      !   !ssnow%Tsoil        = ssnow%tgg - 273.16
      !   ssnow%kth          = 0.3
      !   ssnow%lE           = 0.
      !   ! vh ! should be calculated from soil moisture or be in restart file
      !   !ssnow%sconds(:,:)  = 0.06_r_2    ! vh snow thermal cond (W m-2 K-1),
      !   ! should be in restart file
      !END IF

         first_call = .FALSE.

      ENDIF ! END: if (first_call)       

      !mrd561
      !should be initialized but not sure if here
      ssnow%qrecharge    = 0.0
      ssnow%wtd          = 1.0
      ssnow%rtevap_sat   = 0.0
      ssnow%rtevap_unsat = 0.0
      ssnow%satfrac      = 0.5
      ssnow%qhz          = 0.0
      ssnow%qhlev        = 0.0
      ssnow%wbliq = ssnow%wb - ssnow%wbice

      ! SLI specific initialisations:
      IF(cable_user%SOIL_STRUC=='sli') THEN
         ssnow%h0(:)        = 0.0
         ssnow%S(:,:)       = ssnow%wb(:,:)/soil%ssat_vec !SPREAD(soil%ssat,2,ms)
         ssnow%snowliq(:,:) = 0.0
         ssnow%Tsurface     = 25.0  !why not ssnow%tgg(:,1) - 273.16
         ssnow%nsnow        = 0
         ssnow%Tsoil        = ssnow%tgg - 273.16
         ssnow%kth          = 0.3
         ssnow%lE           = 0.
         ! vh ! should be calculated from soil moisture or be in restart file
         ssnow%sconds(:,:)  = 0.06_r_2    ! vh snow thermal cond (W m-2 K-1),
         ! should be in restart file
      END IF


      IF (cable_user%gw_model) THEN
         ssnow%wb_lake(:) = 0.0  !already prevent drainage unless fully
                                 !saturated and can take from GWwb
         ssnow%wbtot1 = 0.0
         ssnow%wbtot2 = 0.0
         ssnow%wbtot1 = max(0., &
                            (soil%sfc_vec(:,1) - ssnow%wb(:,1))*soil%zse(1)/soil%GWdz(:) )

         WHERE( veg%iveg == 16 .and. &
                ssnow%wb(:,1) < soil%sfc_vec(:,1) .and. &
                ssnow%GWwb(:) .gt. ssnow%wbtot1)

                ssnow%wb(:,1) = soil%sfc_vec(:,1)
                ssnow%GWwb(:) = ssnow%GWwb(:) -  ssnow%wbtot1(:)

          ENDWHERE

          ssnow%wbtot1 = 0.0
          
      ELSE
!     DO J=1, msn
         ssnow%wbtot1 = 0.0
         ssnow%wbtot2 = 0.0
      DO J=1, 1

            WHERE( veg%iveg == 16 .and. ssnow%wb(:,J) < soil%sfc ) 
                ! lakes: remove hard-wired number in future version
            ssnow%wbtot1 = ssnow%wbtot1 + REAL( ssnow%wb(:,J) ) * 1000.0 *     &
                           soil%zse(J)
            ssnow%wb(:,J) = soil%sfc
            ssnow%wbtot2 = ssnow%wbtot2 + REAL( ssnow%wb(:,J) ) * 1000.0 *     &
                           soil%zse(J)
         ENDWHERE

      ENDDO
      ssnow%wb_lake = MAX( ssnow%wbtot2 - ssnow%wbtot1, 0.)
      END IF 

END SUBROUTINE initialize_soilsnow
 
!========================================================================
!========================================================================
!========================================================================
          
SUBROUTINE initialize_roughness( z1_tq, z1_uv, htveg )  
   USE cable_um_tech_mod,   ONLY : um1, rough, veg
   USE cable_common_module, ONLY : ktau_gl
   USE cable_def_types_mod, ONLY : mp
   USE cable_common_module, ONLY : cable_runtime, cable_user
   
   REAL, INTENT(IN), DIMENSION(um1%row_length, um1%rows) ::  z1_tq, z1_uv
   REAL, INTENT(INOUT), DIMENSION(um1%land_pts, um1%ntiles) :: htveg
   INTEGER :: i,j,k,L,n
   REAL, ALLOCATABLE, DIMENSION(:,:) :: jhruff, jhwork

      !--- CABLE roughness type forcings
      CALL um2cable_rr( Z1_TQ, rough%za_tq)
      CALL um2cable_rr( Z1_UV, rough%za_uv)

      ALLOCATE(jhwork (um1%land_pts,um1%ntiles) ) 
      ALLOCATE(jhruff (um1%land_pts,um1%ntiles) ) 

      !Veg height changes seasonally in MOSES hence no updates here due to snow
      jhwork = 0.
      DO N=1,um1%NTILES
        DO K=1,um1%TILE_PTS(N)
          I = um1%TILE_INDEX(K,N)
          jhWORK(I,N) = MAX(.01,HTVEG(I,N))
        ENDDO
      ENDDO

      jHRUFF= 0.01 
      DO l=1,um1%land_pts
        DO n=1,um1%ntiles     
          IF( jHRUFF(L,N) .lt. jhwork(l,n)) jHRUFF(L,:) =  jhwork(l,n)
        ENDDO
      ENDDO
      
      rough%hruff= MAX(0.01,veg%hc)
      rough%hruff_grmx = pack(jHRUFF, um1%l_tile_pts) 

      DEALLOCATE( jhruff, jhwork ) 

END SUBROUTINE initialize_roughness

!========================================================================
!========================================================================
!========================================================================

SUBROUTINE update_kblum_radiation( sw_down, cos_zenith_angle, surf_down_sw )
   USE cable_um_tech_mod!, only : um1, um_rad, kblum_rad
  
   REAL, INTENT(INOUT), DIMENSION(um1%row_length, um1%rows) :: sw_down
   REAL, INTENT(IN), DIMENSION(um1%row_length, um1%rows) :: cos_zenith_angle
   REAL, INTENT(IN), DIMENSION(um1%row_length, um1%rows, 4) :: surf_down_sw 

      !jhan: do you really want to be changing sw_down            
      SW_DOWN = ( surf_down_sw(:,:,1)                                          &
                        + surf_down_sw(:,:,2)                                  &
                        + surf_down_sw(:,:,3)                                  &
                        + surf_down_sw(:,:,4) )                                &
                        * cos_zenith_angle(:,:)

      kblum_rad%SW_DOWN_DIR = ( surf_down_sw(:,:,1)                            &
                        + surf_down_sw(:,:,3) )                                &
                        * cos_zenith_angle(:,:)

      kblum_rad%SW_DOWN_DIF = ( surf_down_sw(:,:,2)                            & 
                              + surf_down_sw(:,:,4) )                          &
                              *cos_zenith_angle(:,:)

      kblum_rad%SW_DOWN_VIS = (surf_down_sw(:,:,1)                             & 
                              + surf_down_sw(:,:,2) )                          &
                              * cos_zenith_angle(:,:)

      kblum_rad%SW_DOWN_NIR = ( surf_down_sw(:,:,3)                            &
                              + surf_down_sw(:,:,4) )                          &
                              *cos_zenith_angle(:,:)
      ! fbeam for VIS
      kblum_rad%FBEAM(:,:,1) = surf_down_sw(:,:,1)                             &
                              * cos_zenith_angle(:,:)                          &
                                 / max( 0.1, kblum_rad%SW_DOWN_VIS )
      ! fbeam for NIR
      kblum_rad%FBEAM(:,:,2) = surf_down_sw(:,:,3)                             &
                              * cos_zenith_angle(:,:)                          &
                              / max( 0.1, kblum_rad%SW_DOWN_NIR )
      !---fbeam for all solar 
      kblum_rad%FBEAM(:,:,3) = kblum_rad%SW_DOWN_DIR /                         &
                              MAX( 0.1, SW_DOWN )
       
END SUBROUTINE Update_kblum_radiation

!========================================================================
!========================================================================
!========================================================================

SUBROUTINE  um2cable_met_rad( cos_zenith_angle)
   USE cable_um_tech_mod, ONLY :um1, kblum_rad, rad, met

   !___ from UM, cosine zenith angle and soil albedo
   REAL, INTENT(INOUT) :: cos_zenith_angle(um1%row_length, um1%rows)

      !--- CABLE met type forcings
      CALL um2cable_rr( cos_zenith_angle, met%coszen)
      CALL um2cable_rr( kblum_rad%SW_DOWN_VIS, met%fsd(:,1))
      CALL um2cable_rr( kblum_rad%SW_DOWN_NIR, met%fsd(:,2))
      
      !--- CABLE radiation type forcings
      !--- kblum_rad% vars are computed in subroutine update_kblum_radiation 
      CALL um2cable_rr( kblum_rad%FBEAM(:,:,1), rad%fbeam(:,1))
      CALL um2cable_rr( kblum_rad%FBEAM(:,:,2), rad%fbeam(:,2))
      CALL um2cable_rr( kblum_rad%FBEAM(:,:,3), rad%fbeam(:,3))

END SUBROUTINE  um2cable_met_rad

!========================================================================
!========================================================================
!========================================================================

!--- UM met forcing vars needed by CABLE commonly have UM dimensions
!---(row_length,rows), which is no good to CABLE. These have to be 
!--- re-packed in a single vector of active tiles. Hence we use 
!--- conditional "mask" l_tile_pts(land_pts,ntiles) which is .true.
!--- if the land point is/has an active tile
SUBROUTINE um2cable_rr(umvar,cablevar)
   USE cable_def_types_mod, ONLY : mp
   USE cable_um_tech_mod,   ONLY :um1
 
   REAL, INTENT(IN), DIMENSION(um1%row_length, um1%rows) :: umvar   
   REAL, INTENT(INOUT), DIMENSION(mp) :: cablevar
   REAL, DIMENSION(um1%land_pts,um1%ntiles) :: fvar   
   INTEGER :: n,k,l,j,i

      fvar = 0.0
      DO N=1,um1%NTILES                     
         DO K=1,um1%TILE_PTS(N)
            L = um1%TILE_INDEX(K,N)
            J=(um1%LAND_INDEX(L)-1)/um1%row_length + 1
            I = um1%LAND_INDEX(L) - (J-1)*um1%row_length
            fvar(L,N) = umvar(I,J)
         ENDDO
      ENDDO
      cablevar =  pack(fvar,um1%l_tile_pts)

END SUBROUTINE um2cable_rr

!========================================================================

SUBROUTINE um2cable_irr(umvar,cablevar)
   USE cable_def_types_mod, ONLY : mp
   USE cable_um_tech_mod,   ONLY :um1
 
   integer, INTENT(IN), DIMENSION(um1%row_length, um1%rows) :: umvar   
   integer, INTENT(INOUT), DIMENSION(mp) :: cablevar
   integer, DIMENSION(um1%land_pts,um1%ntiles) :: fvar   
   INTEGER :: n,k,l,j,i

      fvar = 0.0
      DO N=1,um1%NTILES                     
         DO K=1,um1%TILE_PTS(N)
            L = um1%TILE_INDEX(K,N)
            J=(um1%LAND_INDEX(L)-1)/um1%row_length + 1
            I = um1%LAND_INDEX(L) - (J-1)*um1%row_length
            fvar(L,N) = umvar(I,J)
         ENDDO
      ENDDO
      cablevar =  pack(fvar,um1%l_tile_pts)

END SUBROUTINE um2cable_irr


!========================================================================
!========================================================================

!--- UM met forcing vars needed by CABLE which have UM dimensions
!---(land_points)[_lp], which is no good to cable. These have to be 
!--- re-packed in a single vector of active tiles. Hence we use 
!--- conditional "mask" l_tile_pts(land_pts,ntiles) which is .true.
!--- if the land point is/has an active tile
SUBROUTINE um2cable_lp(umvar, defaultin, cablevar, soiltype, skip )
   USE cable_def_types_mod, ONLY : mp, mstype
   USE cable_um_tech_mod,   ONLY :um1
  
   REAL, INTENT(IN), DIMENSION(um1%land_pts) :: umvar
   REAL, INTENT(IN), DIMENSION(mstype) :: defaultin    
   REAL, INTENT(INOUT), DIMENSION(mp) :: cablevar
   INTEGER, INTENT(INOUT), DIMENSION(mp) :: soiltype
   REAL, DIMENSION(:,:), ALLOCATABLE:: fvar   
   LOGICAL, OPTIONAL :: skip
   INTEGER :: n,k,l,i

         
      ALLOCATE( fvar(um1%land_pts,um1%ntiles) )
      fvar = 0.0

      ! loop over Ntiles
      DO N=1,um1%NTILES
         ! loop over number of points per tile
         DO K=1,um1%TILE_PTS(N)
            ! index of each point per tile in an array of dim=(land_pts,ntiles)
            L = um1%TILE_INDEX(K,N)
            ! at this point fvar=umvar, ELSE=0.0 
            fvar(L,N) = umvar(L)
            ! unless explicitly SKIPPED by including arg in subr call
            IF(.NOT. PRESENT(skip) ) THEN
               ! on perma frost tile, set fvar=defaultin
               IF( N == um1%ntiles ) THEN
                  fvar(L,N) =  defaultin(9)
               ENDIF
            ENDIF
         ENDDO
      ENDDO
     
      cablevar     =  PACK(fvar,um1%l_tile_pts)
  
      ! unless explicitly SKIPPED by including arg in subr call
      IF(.NOT. PRESENT(skip) ) THEN
         DO i=1,mp
            ! soiltype=9 for perma-frost tiles 
            IF(soiltype(i)==9) cablevar(i) =  defaultin(9)         
         ENDDO        
      ENDIF
   
      DEALLOCATE(fvar)

END SUBROUTINE um2cable_lp
 

!========================================================================
!========================================================================
!========================================================================

SUBROUTINE init_bgc_vars() 
   USE cable_def_types_mod, ONLY : ncs, ncp 
   USE cable_um_tech_mod,   ONLY : bgc, veg   
   USE cable_common_module, ONLY : vegin
   
   INTEGER :: k

   ! note that ratecp and ratecs are the same for all veg at the moment. (BP)
    DO k=1,ncp
       bgc%cplant(:,k) = vegin%cplant(k,veg%iveg)
       bgc%ratecp(k) = vegin%ratecp(k,1)
    ENDDO
    DO k=1,ncs
      bgc%csoil(:,k) = vegin%csoil(k,veg%iveg)
      bgc%ratecs(k) = vegin%ratecs(k,1)
    ENDDO

END SUBROUTINE init_bgc_vars

!========================================================================
!========================================================================
!========================================================================

subroutine init_sumflux_zero() 
   USE cable_um_tech_mod, ONLY : sum_flux
      sum_flux%sumpn = 0.; sum_flux%sumrp = 0.; sum_flux%sumrpw = 0.
      sum_flux%sumrpr = 0.; sum_flux%sumrs = 0.; sum_flux%sumrd = 0.
      sum_flux%dsumpn = 0.; sum_flux%dsumrp = 0.; sum_flux%dsumrs = 0.
      sum_flux%dsumrd = 0.; sum_flux%sumxrp = 0.;  sum_flux%sumxrs = 0.
END SUBROUTINE init_sumflux_zero 

!========================================================================
!========================================================================
!========================================================================

SUBROUTINE alloc_cable_types()
   USE cable_def_types_mod, ONLY : mp, alloc_cbm_var
   USE cable_um_tech_mod,   ONLY : air, canopy, met, bal, rad, rough,          &
                                   soil, ssnow, sum_flux, veg, bgc

      CALL alloc_cbm_var(air, mp)
      CALL alloc_cbm_var(canopy, mp)
      CALL alloc_cbm_var(met, mp)
      CALL alloc_cbm_var(bal, mp)
      CALL alloc_cbm_var(rad, mp)
      CALL alloc_cbm_var(rough, mp)
      CALL alloc_cbm_var(soil, mp)
      CALL alloc_cbm_var(ssnow, mp)
      CALL alloc_cbm_var(sum_flux, mp)
      CALL alloc_cbm_var(veg, mp)
      CALL alloc_cbm_var(bgc, mp)

END SUBROUTINE alloc_cable_types

!========================================================================
!========================================================================
!========================================================================


END MODULE cable_um_init_subrs_mod




