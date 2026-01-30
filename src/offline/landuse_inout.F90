  subroutine landuse_data(mlon,mlat,landmask,arealand,luc_atransit,luc_fharvw,luc_xluh2cable)
  use netcdf
  use cable_abort_module,   ONLY: nc_abort 
  use cable_common_module,  ONLY: filename
  USE cable_def_types_mod,  ONLY: mland,r_2
  use landuse_constant,     ONLY: mstate,mvmax,mharvw
  IMPLICIT NONE

  integer     mlon,mlat
  real(r_2), dimension(mland,mvmax,mvmax)         :: luc_atransit
  real(r_2), dimension(mland,mharvw)              :: luc_fharvw
  real(r_2), dimension(mland,mvmax,mstate)        :: luc_xluh2cable
  integer,    dimension(mlon,mlat)                :: landmask
  real(r_2),  dimension(mland)                    :: arealand
  ! "mland" variables
  real(r_2),  dimension(:,:),      allocatable    :: areax    
  !
  integer ivt,ee,hh,np,p,q,np1
  integer ncid,ok,xID,yID,varID,i,j,m,mpx

    ! get " mlon mlat landmask" from "gridinfo"
    ok = NF90_OPEN(filename%type, 0, ncid)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error opening grid info file.')

    ok = NF90_INQ_DIMID(ncid, 'longitude', xID)
    IF (ok /= NF90_NOERR) ok = NF90_INQ_DIMID(ncid, 'x', xID)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error inquiring x dimension.')
    ok = NF90_INQUIRE_DIMENSION(ncid, xID, LEN=mlon)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error getting x dimension.')
    ok = NF90_INQ_DIMID(ncid, 'latitude', yID)
    IF (ok /= NF90_NOERR) ok = NF90_INQ_DIMID(ncid, 'y', yID)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error inquiring y dimension.')
    ok = NF90_INQUIRE_DIMENSION(ncid, yID, LEN=mlat)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error getting y dimension.')

    allocate(areax(mlon,mlat))
     
    ok = NF90_INQ_VARID(ncid, 'area', varID)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok,                                    &
                  'Error finding variable area')
    ok = NF90_GET_VAR(ncid, varID, areax)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok,                                    &
                  'Error reading variable longitude.')

    ok = NF90_CLOSE(ncid)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error closing grid info file.')

     m=0; arealand= 0.0
     do i=1,mlon
     do j=1,mlat
        if(areax(i,j) >0.01) then
           landmask(i,j) = 1
           m=m+1
           arealand(m) = areax(i,j)
        else
           landmask(i,j) =0
        endif
     enddo
     enddo
     if(m/=mland) then
        print *, 'mland not consistent: check gridinof area'
        stop
     endif   

     ! get the mapping matrix (landuse type to PFT)
     call landuse_getxluh2(mlat,mlon,landmask,filename%fxluh2cable,luc_xluh2cable)    !"xluh2cable"
     call landuse_getdata(mlat,mlon,landmask,filename%fxpft,luc_atransit,luc_fharvw)
  end subroutine landuse_data


SUBROUTINE landuse_getxluh2(mlat,mlon,landmask,fxluh2cable,luc_xluh2cable)
! get data: luc%fprimary; luc%fsecondary
  USE netcdf
  USE cable_def_types_mod,  ONLY: mland,r_2
  use landuse_constant,     ONLY: mstate,mvmax
  IMPLICIT NONE
  character*500 fxluh2cable
  integer   mlat,mlon
  integer,  dimension(mlon,mlat)                 :: landmask
  real(r_2),dimension(mland,mvmax,mstate)        :: luc_xluh2cable
  ! local variables
  real(r_2),   dimension(:,:,:,:), allocatable   :: xluh2cable
  integer ok,ncid2,varxid
  integer i,j,k,m,v,s

    allocate(xluh2cable(mlon,mlat,21,mstate))
    ok = nf90_open(fxluh2cable,nf90_nowrite,ncid2)
    ok = nf90_inq_varid(ncid2,"xluh2cable",varxid)
    ok = nf90_get_var(ncid2,varxid,xluh2cable)
    ok = nf90_close(ncid2)
    ! assig the values of luc%variables
    luc_xluh2cable(:,:,:) = 0.0
    m = 0
    do i=1,mlon
    do j=1,mlat
       if(landmask(i,j) ==1) then
          m= m +1
          do k=1,10
             luc_xluh2cable(m,k,:) = xluh2cable(i,j,k,:)
          enddo
          luc_xluh2cable(m,16,:) = xluh2cable(i,j,11,:)+xluh2cable(i,j,16,:)
          luc_xluh2cable(m,14,:) = xluh2cable(i,j,14,:)
          luc_xluh2cable(m,15,:) = xluh2cable(i,j,15,:)
          luc_xluh2cable(m,17,:) = xluh2cable(i,j,17,:)

          luc_xluh2cable(m,11,:) = xluh2cable(i,j,18,:)
          luc_xluh2cable(m,12,:) = xluh2cable(i,j,19,:)
          luc_xluh2cable(m,13,:) = xluh2cable(i,j,21,:)

          do s=1,mstate
             do v=1,mvmax
                luc_xluh2cable(m,v,s) = luc_xluh2cable(m,v,s)/sum(luc_xluh2cable(m,1:mvmax,s))
             enddo
          enddo
       endif
    enddo
    enddo

    deallocate(xluh2cable)

 END SUBROUTINE landuse_getxluh2

SUBROUTINE landuse_getdata(mlat,mlon,landmask,fxpft,luc_atransit,luc_fharvw)
! get LUC data
  USE netcdf
  USE cable_def_types_mod,  ONLY: mland,r_2
  use landuse_constant,     ONLY: mstate,mvmax,mharvw
  IMPLICIT NONE
  character*500 fxpft
  integer mlat,mlon
  integer,   dimension(mlon,mlat)               :: landmask 
  real(r_2), dimension(mland,mvmax,mvmax)       :: luc_atransit
  real(r_2), dimension(mland,mharvw)            :: luc_fharvw
  ! local variables
  real(r_2),  dimension(:,:,:),   allocatable   :: fracharvw
  real(r_2),  dimension(:,:,:,:), allocatable   :: transitx
  integer  ok,ncid1,varxid
  integer  i,j,m,k,ivt

    allocate(fracharvw(mlon,mlat,mharvw))
    allocate(transitx(mlon,mlat,mvmax,mvmax))

    ok = nf90_open(fxpft,nf90_nowrite,ncid1)
    ok = nf90_inq_varid(ncid1,"harvest",varxid)
    ok = nf90_get_var(ncid1,varxid,fracharvw)
    ok = nf90_inq_varid(ncid1,"transition",varxid)
    ok = nf90_get_var(ncid1,varxid,transitx)
    ok = nf90_close(ncid1)

    ! assig the values of luc%variables
    luc_fharvw(:,:) =0.0; luc_atransit(:,:,:)=0.0
    m = 0
    do i=1,mlon
    do j=1,mlat
       if(landmask(i,j) ==1) then
          m= m +1
          luc_atransit(m,:,:)   = transitx(i,j,:,:)
          luc_fharvw(m,:)       = fracharvw(i,j,:)
       endif
    enddo
    enddo
    
    deallocate(fracharvw)
    deallocate(transitx)
END SUBROUTINE landuse_getdata

  subroutine create_new_gridinfo(fgridold,fgridnew,mlon,mlat,landmask,patchfrac_new)
  use netcdf
  use cable_abort_module,   ONLY: nc_abort 
  use cable_common_module,  ONLY: filename
  USE cable_def_types_mod,  ONLY: r_2,nrb,ms
  USE cable_IO_vars_module, ONLY: logn
  use landuse_constant,     ONLY: mvmax, thresh_frac

  implicit none
  character*500 fgridold,fgridnew
  integer       mlon,mlat
  integer,      dimension(mlon,mlat)            :: landmask
  real(r_2),    dimension(mlon,mlat,mvmax)      :: patchfrac_new
  ! local variables
  integer,      parameter                       :: missint  = -99999
  real,         parameter                       :: missreal = 1.0e23
  integer,      parameter                       :: time12=12 
  integer,      dimension(mlon,mlat)            :: isoil_y,soilorder_y
  real(r_2),    dimension(mlon,mlat)            :: ndep_y,nfix_y,pdust_y,pwea_y
  real(r_2),    dimension(mlon,mlat,nrb)        :: albedo_y
  real(r_2),    dimension(mlon)                 :: longitude_y
  real(r_2),    dimension(mlat)                 :: latitude_y
  real(r_2),    dimension(mlon,mlat,time12)     :: lai_y,snowdepth_y
  real(r_2),    dimension(mlon,mlat,ms,time12)  :: soilmoist_y, soiltemp_y
  real(r_2),    dimension(mlon,mlat)            :: albedo2_y,area_y
  real(r_2),    dimension(mlon,mlat)            :: bch_y,clay_y,cnsd_y,css_y,hyds_y,rhosoil_y, &
                                                   sand_y,sfc_y,silt_y,ssat_y,sucs_y,swilt_y
  real(r_2),    dimension(mlon,mlat,mvmax)      :: patchfrac_y
  integer,      dimension(mlon,mlat,mvmax)      :: iveg_y
  integer,      dimension(mvmax)                :: tmpint
  real(r_2),    dimension(mvmax)                :: tmpx
  !
  integer dim_radid,dim_mlatid,dim_mlonid,dim_timeid,dim_patchid,dim_soilid
  integer var_albedoid,var_laiid,var_ndepid,var_nfixid,var_pdustid,var_pweaid
  integer var_snowdepthid,var_soilmoistid,var_soilorderid,var_soiltempid
  integer var_albedo2id,var_areaid,var_bchid,var_clayid,var_cnsdid,var_cssid
  integer var_hydsid,var_isoilid,var_latitudeid,var_longitudeid,var_rhosoilid
  integer var_sandid,var_sfcid,var_siltid,var_ssatid,var_sucsid,var_swiltid
  integer var_ivegid,var_patchfracid
  integer varxid
  integer ncid0,ncid11,ok
  integer i,j,k

     ! then set non-land patch vegtype to -1
     ! then order patch within each land cell by area fraction from largest to the smallest
     ! sort patch by area fraction
     patchfrac_y = 0.0
     do j=1,mlat
     do i=1,mlon
        if(landmask(i,j) ==1) then
           do k=1,mvmax
              tmpx(k)  = patchfrac_new(i,j,k)
              tmpint(k)= k
           enddo

           call sort(mvmax,tmpx,tmpint)

           do k=1,mvmax
              patchfrac_y(i,j,k) = tmpx(k)
              if(patchfrac_y(i,j,k)>=thresh_frac) then
                 iveg_y(i,j,k)  = tmpint(k)
              else
                 iveg_y(i,j,k) = -1
              endif
           enddo
        else
           iveg_y(i,j,:) = -1
           patchfrac_y(i,j,:) = 0.0
        endif
     enddo
     enddo

     write(logn,*) 'landuse on: create new gridinfo'

     ok = NF90_OPEN(fgridold,0,ncid0)
     if(ok/=nf90_noerr) call nc_abort(ok, 'file opening error')

     ok = NF90_INQ_VARID(ncid0,'latitude',varxid)
     ok = NF90_GET_VAR(ncid0,varxid,latitude_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in reading latitude')

     ok = NF90_INQ_VARID(ncid0,'longitude',varxid)
     ok = NF90_GET_VAR(ncid0,varxid,longitude_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in longitude')

     ok = NF90_INQ_VARID(ncid0,'area',varxid)
     ok = NF90_GET_VAR(ncid0,varxid,area_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in area')

     ok = NF90_INQ_VARID(ncid0,'isoil',varxid)
     ok = NF90_GET_VAR(ncid0,varxid,isoil_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in reading isoil') 

     ok = NF90_INQ_VARID(ncid0,'SoilOrder',varxid)
     ok = NF90_GET_VAR(ncid0,varxid,soilorder_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in reading SoilOrder') 

     ok = NF90_INQ_VARID(ncid0,'SnowDepth',varxid)
     ok = NF90_GET_VAR(ncid0,varxid,snowdepth_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in reading SnowDepth') 

     ok = NF90_INQ_VARID(ncid0,'LAI',varxid)
     ok = NF90_GET_VAR(ncid0,varxid,lai_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in reading LAI') 

     ok = NF90_INQ_VARID(ncid0,'SoilMoist',varxid)
     ok = NF90_GET_VAR(ncid0,varxid,soilmoist_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in reading SoilMoist') 

     ok = NF90_INQ_VARID(ncid0,'SoilTemp',varxid)
     ok = NF90_GET_VAR(ncid0,varxid,soiltemp_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in reading SoilTemp') 

     ok = NF90_INQ_VARID(ncid0,'Albedo',varxid)
     ok = NF90_GET_VAR(ncid0,varxid,albedo_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in reading Albedo') 

     ok = NF90_INQ_VARID(ncid0,'albedo2',varxid)
     ok = NF90_GET_VAR(ncid0,varxid,albedo2_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in reading albedo2') 

     ok = NF90_INQ_VARID(ncid0,'bch',varxid)
     ok = NF90_GET_VAR(ncid0,varxid,bch_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in reading bch') 

     ok = NF90_INQ_VARID(ncid0,'clay',varxid)
     ok = NF90_GET_VAR(ncid0,varxid,clay_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in reading clay') 

     ok = NF90_INQ_VARID(ncid0,'cnsd',varxid)
     ok = NF90_GET_VAR(ncid0,varxid,cnsd_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in reading cnsd') 

     ok = NF90_INQ_VARID(ncid0,'css',varxid)
     ok = NF90_GET_VAR(ncid0,varxid,css_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in reading css') 

     ok = NF90_INQ_VARID(ncid0,'hyds',varxid)
     ok = NF90_GET_VAR(ncid0,varxid,hyds_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in reading hyds') 

     ok = NF90_INQ_VARID(ncid0,'rhosoil',varxid)
     ok = NF90_GET_VAR(ncid0,varxid,rhosoil_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in reading rhosoil')
 
     ok = NF90_INQ_VARID(ncid0,'sand',varxid)
     ok = NF90_GET_VAR(ncid0,varxid,sand_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in reading sand') 

     ok = NF90_INQ_VARID(ncid0,'sfc',varxid)
     ok = NF90_GET_VAR(ncid0,varxid,sfc_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in reading sfc') 

     ok = NF90_INQ_VARID(ncid0,'silt',varxid)
     ok = NF90_GET_VAR(ncid0,varxid,silt_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in reading silt') 

     ok = NF90_INQ_VARID(ncid0,'ssat',varxid)
     ok = NF90_GET_VAR(ncid0,varxid,ssat_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in reading ssat') 

     ok = NF90_INQ_VARID(ncid0,'sucs',varxid)
     ok = NF90_GET_VAR(ncid0,varxid,sucs_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in reading sucs') 

     ok = NF90_INQ_VARID(ncid0,'swilt',varxid)
     ok = NF90_GET_VAR(ncid0,varxid,swilt_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in reading swilt') 

     ok = NF90_INQ_VARID(ncid0,'Ndep',varxid)
     ok = NF90_GET_VAR(ncid0,varxid,ndep_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in reading Ndep') 

     ok = NF90_INQ_VARID(ncid0,'Nfix',varxid)
     ok = NF90_GET_VAR(ncid0,varxid,nfix_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in reading Nfix') 

     ok = NF90_INQ_VARID(ncid0,'Pdust',varxid)
     ok = NF90_GET_VAR(ncid0,varxid,pdust_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in reading Pdust') 

     ok = NF90_INQ_VARID(ncid0,'Pwea',varxid)
     ok = NF90_GET_VAR(ncid0,varxid,pwea_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in reading Pwea') 

     ok = nf90_close(ncid0)

     !
     ! create ACCESS ESM pool size file
     ok = nf90_create(fgridnew,nf90_clobber,ncid11)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in opening fgridnew') 

     ! define dimensions
     ok = nf90_def_dim(ncid11,'rad',nrb,dim_radid)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining rad') 

     ok = nf90_def_dim(ncid11,'latitude',mlat,dim_mlatid)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining latitude') 

     ok = nf90_def_dim(ncid11,'longitude',mlon,dim_mlonid)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining longitude') 

     ok = nf90_def_dim(ncid11,'time',12,dim_timeid)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining time') 

     ok = nf90_def_dim(ncid11,'soil',ms,dim_soilid)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining soil') 

     ok = nf90_def_dim(ncid11,'patch',mvmax,dim_patchid)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining patch') 

     ! define variables
     ok = nf90_def_var(ncid11,'Albedo',nf90_float,(/dim_mlonid,dim_mlatid,dim_radid/),var_albedoid)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining Albedo') 
     ok = nf90_put_att(ncid11,var_albedoid,'units',' ---')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining Albedo')
     ok = nf90_put_att(ncid11,var_albedoid,'long_name',' snow-free bareground albedo')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining Albedo')
     ok = nf90_put_att(ncid11,var_albedoid,'_FillValue',missreal)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining patch')
     ok = nf90_put_att(ncid11,var_albedoid,'missing_value',missreal)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining patch')

     ok = nf90_def_var(ncid11,'LAI',nf90_float,(/dim_mlonid,dim_mlatid,dim_timeid/),var_laiid)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining LAI')
     ok = nf90_put_att(ncid11,var_laiid,'units',' ---')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining LAI')
     ok = nf90_put_att(ncid11,var_laiid,'long_name',' monthly canopy LAI')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining LAI')
     ok = nf90_put_att(ncid11,var_laiid,'_FillValue',missreal)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining LAI')
     ok = nf90_put_att(ncid11,var_laiid,'missing_value',missreal)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining LAI')

     ok = nf90_def_var(ncid11,'Ndep',nf90_float,(/dim_mlonid,dim_mlatid/),var_ndepid)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining Ndep')
     ok = nf90_put_att(ncid11,var_ndepid,'units',' gN m-2 yr-1')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining Ndep')
     ok = nf90_put_att(ncid11,var_ndepid,'long_name',' Annual N deposition')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining Ndep')
     ok = nf90_put_att(ncid11,var_ndepid,'_FillValue',missreal)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining Ndep')
     ok = nf90_put_att(ncid11,var_ndepid,'missing_value',missreal)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining Ndep')

     ok = nf90_def_var(ncid11,'Nfix',nf90_float,(/dim_mlonid,dim_mlatid/),var_nfixid)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining Nfix')
     ok = nf90_put_att(ncid11,var_nfixid,'units',' gN m-2 yr-1')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining Nfix')
     ok = nf90_put_att(ncid11,var_nfixid,'long_name',' Annual N fixation')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining Nfix')
     ok = nf90_put_att(ncid11,var_nfixid,'_FillValue',missreal)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining Nfix')
     ok = nf90_put_att(ncid11,var_nfixid,'missing_value',missreal)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining Nfix')

     ok = nf90_def_var(ncid11,'Pdust',nf90_float,(/dim_mlonid,dim_mlatid/),var_pdustid)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining patch')
     ok = nf90_put_att(ncid11,var_pdustid,'units',' gP m-2 yr-1')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining Pdust')
     ok = nf90_put_att(ncid11,var_pdustid,'long_name',' Annual P deposition')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining Pdust')
     ok = nf90_put_att(ncid11,var_pdustid,'_FillValue',missreal)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining Pdust')
     ok = nf90_put_att(ncid11,var_pdustid,'missing_value',missreal)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining Pdust')

     ok = nf90_def_var(ncid11,'Pwea',nf90_float,(/dim_mlonid,dim_mlatid/),var_pweaid)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining Pwea')
     ok = nf90_put_att(ncid11,var_pweaid,'units',' gP m-2 yr-1')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining Pwea')
     ok = nf90_put_att(ncid11,var_pweaid,'long_name',' P weathering rate')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining Pwea')
     ok = nf90_put_att(ncid11,var_pweaid,'_FillValue',missreal)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining Pwea')
     ok = nf90_put_att(ncid11,var_pweaid,'missing_value',missreal)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining Pwea')

     ok = nf90_def_var(ncid11,'SnowDepth',nf90_float,(/dim_mlonid,dim_mlatid,dim_timeid/),var_snowdepthid)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining SnowDepth')
     ok = nf90_put_att(ncid11,var_snowdepthid,'units',' m')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining SnowDepth')
     ok = nf90_put_att(ncid11,var_snowdepthid,'long_name',' snow depth')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining SnowDepth')
     ok = nf90_put_att(ncid11,var_snowdepthid,'_FillValue',missreal)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining SnowDepth')
     ok = nf90_put_att(ncid11,var_snowdepthid,'missing_value',missreal)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining SnowDepth')

     ok = nf90_def_var(ncid11,'SoilMoist',nf90_float,(/dim_mlonid,dim_mlatid,dim_soilid,dim_timeid/),var_soilmoistid)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining SoilMoist')
     ok = nf90_put_att(ncid11,var_soilmoistid,'units',' m3/m3')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining SoilMoist')
     ok = nf90_put_att(ncid11,var_soilmoistid,'long_name',' soil moisture profile')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining SoilMoist')
     ok = nf90_put_att(ncid11,var_soilmoistid,'_FillValue',missreal)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining SoilMoist')
     ok = nf90_put_att(ncid11,var_soilmoistid,'missing_value',missreal)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining SoilMoist')

     ok = nf90_def_var(ncid11,'SoilOrder',nf90_int,(/dim_mlonid,dim_mlatid/),var_soilorderid)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining SoilOrder')
     ok = nf90_put_att(ncid11,var_soilorderid,'units',' class ')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining SoilOrder')
     ok = nf90_put_att(ncid11,var_soilorderid,'long_name',' soil order class')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining SoilOrder')
     ok = nf90_put_att(ncid11,var_soilorderid,'_FillValue',missint)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining SoilOrder')
     ok = nf90_put_att(ncid11,var_soilorderid,'missing_value',missint)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining SoilOrder')

     ok = nf90_def_var(ncid11,'SoilTemp',nf90_float,(/dim_mlonid,dim_mlatid,dim_soilid,dim_timeid/),var_soiltempid)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining SoilTemp')
     ok = nf90_put_att(ncid11,var_soiltempid,'units',' K')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining SoilTemp')
     ok = nf90_put_att(ncid11,var_soiltempid,'long_name',' soil temperature profile')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining SoilTemp')
     ok = nf90_put_att(ncid11,var_soiltempid,'_FillValue',missreal)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining SoilTemp')
     ok = nf90_put_att(ncid11,var_soiltempid,'missing_value',missreal)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining SoilTemp')

     ok = nf90_def_var(ncid11,'albedo2',nf90_float,(/dim_mlonid,dim_mlatid/),var_albedo2id)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining albedo2')
     ok = nf90_put_att(ncid11,var_albedo2id,'units',' ---')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining albedo2')
     ok = nf90_put_att(ncid11,var_albedo2id,'long_name',' snow-free soil albedo')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining albedo2')
     ok = nf90_put_att(ncid11,var_albedo2id,'_FillValue',missreal)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining albedo2')
     ok = nf90_put_att(ncid11,var_albedo2id,'missing_value',missreal)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining albedo2')

     ok = nf90_def_var(ncid11,'area',nf90_float,(/dim_mlonid,dim_mlatid/),var_areaid)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining area')
     ok = nf90_put_att(ncid11,var_areaid,'units',' m2')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining area')
     ok = nf90_put_att(ncid11,var_areaid,'long_name',' land area')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining area')
     ok = nf90_put_att(ncid11,var_areaid,'_FillValue',missreal)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining area')
     ok = nf90_put_att(ncid11,var_areaid,'missing_value',missreal)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining patch')

     ok = nf90_def_var(ncid11,'bch',nf90_float,(/dim_mlonid,dim_mlatid/),var_bchid)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining bch')
     ok = nf90_put_att(ncid11,var_bchid,'units',' ---')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining bch')
     ok = nf90_put_att(ncid11,var_bchid,'long_name',' Clapp-Hornberger B coefficient')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining bch')
     ok = nf90_put_att(ncid11,var_bchid,'_FillValue',missreal)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining bch')
     ok = nf90_put_att(ncid11,var_bchid,'missing_value',missreal)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining bch')

     ok = nf90_def_var(ncid11,'clay',nf90_float,(/dim_mlonid,dim_mlatid/),var_clayid)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining clay')
     ok = nf90_put_att(ncid11,var_clayid,'units',' fraction')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining clay')
     ok = nf90_put_att(ncid11,var_clayid,'long_name',' clay fraction')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining clay')
     ok = nf90_put_att(ncid11,var_clayid,'_FillValue',missreal)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining clay')
     ok = nf90_put_att(ncid11,var_clayid,'missing_value',missreal)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining clay')

     ok = nf90_def_var(ncid11,'cnsd',nf90_float,(/dim_mlonid,dim_mlatid/),var_cnsdid)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining cnsd')
     ok = nf90_put_att(ncid11,var_cnsdid,'units',' W m-2 K-1')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining cnsd')
     ok = nf90_put_att(ncid11,var_cnsdid,'long_name',' thermal conductivity')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining cnsd')
     ok = nf90_put_att(ncid11,var_cnsdid,'_FillValue',missreal)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining cnsd')
     ok = nf90_put_att(ncid11,var_cnsdid,'missing_value',missreal)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining cnsd')

     ok = nf90_def_var(ncid11,'css',nf90_float,(/dim_mlonid,dim_mlatid/),var_cssid)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining css')
     ok = nf90_put_att(ncid11,var_cssid,'units',' J/kg/K')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining css')
     ok = nf90_put_att(ncid11,var_cssid,'long_name',' soil specific heat capacity')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining css')
     ok = nf90_put_att(ncid11,var_cssid,'_FillValue',missreal)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining css')
     ok = nf90_put_att(ncid11,var_cssid,'missing_value',missreal)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining css')

     ok = nf90_def_var(ncid11,'hyds',nf90_float,(/dim_mlonid,dim_mlatid/),var_hydsid)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining hyds')
     ok = nf90_put_att(ncid11,var_hydsid,'units',' m/s')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining hyds')
     ok = nf90_put_att(ncid11,var_hydsid,'long_name',' saturated hydraulic conductivity')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining hyds')
     ok = nf90_put_att(ncid11,var_hydsid,'_FillValue',missreal)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining hyds')
     ok = nf90_put_att(ncid11,var_hydsid,'missing_value',missreal)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining hyds')

     ok = nf90_def_var(ncid11,'isoil',nf90_int,(/dim_mlonid,dim_mlatid/),var_isoilid)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining isoil')
     ok = nf90_put_att(ncid11,var_isoilid,'units',' --')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining isoil')
     ok = nf90_put_att(ncid11,var_isoilid,'long_name',' Zobler soil texture class')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining isoil')
     ok = nf90_put_att(ncid11,var_isoilid,'_FillValue',missint)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining isoil')
     ok = nf90_put_att(ncid11,var_isoilid,'missing_value',missint)
     if(ok/=nf90_noerr)call nc_abort(ok, 'error in defining isoil')

     ok = nf90_def_var(ncid11,'latitude',nf90_float,(/dim_mlatid/),var_latitudeid)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining latitude')
     ok = nf90_put_att(ncid11,var_latitudeid,'units',' degree')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining latitude')
     ok = nf90_put_att(ncid11,var_latitudeid,'long_name',' latitude')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining latitude')
     ok = nf90_put_att(ncid11,var_latitudeid,'_FillValue',missreal)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining latitude')
     ok = nf90_put_att(ncid11,var_latitudeid,'missing_value',missreal)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining latitude')

     ok = nf90_def_var(ncid11,'longitude',nf90_float,(/dim_mlonid/),var_longitudeid)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining longitude')
     ok = nf90_put_att(ncid11,var_longitudeid,'units',' degree')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining longitude')
     ok = nf90_put_att(ncid11,var_longitudeid,'long_name',' longitude')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining longitude')
     ok = nf90_put_att(ncid11,var_longitudeid,'_FillValue',missreal)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining longitude')
     ok = nf90_put_att(ncid11,var_longitudeid,'missing_value',missreal)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining longitude')

     ok = nf90_def_var(ncid11,'rhosoil',nf90_float,(/dim_mlonid,dim_mlatid/),var_rhosoilid)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining rhosoil')
     ok = nf90_put_att(ncid11,var_rhosoilid,'units',' kg m-3')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining rhosoil')
     ok = nf90_put_att(ncid11,var_rhosoilid,'long_name',' soil bulk density')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining rhosoil')
     ok = nf90_put_att(ncid11,var_rhosoilid,'_FillValue',missreal)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining rhosoil')
     ok = nf90_put_att(ncid11,var_rhosoilid,'missing_value',missreal)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining rhosoil')

     ok = nf90_def_var(ncid11,'sand',nf90_float,(/dim_mlonid,dim_mlatid/),var_sandid)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining sand')
     ok = nf90_put_att(ncid11,var_sandid,'units',' ---')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining sand')
     ok = nf90_put_att(ncid11,var_sandid,'long_name',' sand fraction')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining sand')
     ok = nf90_put_att(ncid11,var_sandid,'_FillValue',missreal)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining sand')
     ok = nf90_put_att(ncid11,var_sandid,'missing_value',missreal)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining sand')

     ok = nf90_def_var(ncid11,'sfc',nf90_float,(/dim_mlonid,dim_mlatid/),var_sfcid)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining sfc')
     ok = nf90_put_att(ncid11,var_sfcid,'units',' m3/m3')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining sfc')
     ok = nf90_put_att(ncid11,var_sfcid,'long_name',' soil water at field capacity')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining sfc')
     ok = nf90_put_att(ncid11,var_sfcid,'_FillValue',missreal)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining sfc')
     ok = nf90_put_att(ncid11,var_sfcid,'missing_value',missreal)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining sfc')

     ok = nf90_def_var(ncid11,'silt',nf90_float,(/dim_mlonid,dim_mlatid/),var_siltid)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining silt')
     ok = nf90_put_att(ncid11,var_siltid,'units',' fraction')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining silt')
     ok = nf90_put_att(ncid11,var_siltid,'long_name',' soil silt fraction')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining silt')
     ok = nf90_put_att(ncid11,var_siltid,'_FillValue',missreal)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining silt')
     ok = nf90_put_att(ncid11,var_siltid,'missing_value',missreal)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining silt')

     ok = nf90_def_var(ncid11,'ssat',nf90_float,(/dim_mlonid,dim_mlatid/),var_ssatid)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining ssat')
     ok = nf90_put_att(ncid11,var_ssatid,'units',' m3/m3')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining ssat')
     ok = nf90_put_att(ncid11,var_ssatid,'long_name',' soil water at saturation')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining ssat')
     ok = nf90_put_att(ncid11,var_ssatid,'_FillValue',missreal)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining ssat')
     ok = nf90_put_att(ncid11,var_ssatid,'missing_value',missreal)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining ssat')

     ok = nf90_def_var(ncid11,'sucs',nf90_float,(/dim_mlonid,dim_mlatid/),var_sucsid)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining sucs')
     ok = nf90_put_att(ncid11,var_sucsid,'units',' m')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining sucs')
     ok = nf90_put_att(ncid11,var_sucsid,'long_name',' soil sunction at saturation')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining sucs')
     ok = nf90_put_att(ncid11,var_sucsid,'_FillValue',missreal)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining sucs')
     ok = nf90_put_att(ncid11,var_sucsid,'missing_value',missreal)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining sucs')

     ok = nf90_def_var(ncid11,'swilt',nf90_float,(/dim_mlonid,dim_mlatid/),var_swiltid)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining swilt')
     ok = nf90_put_att(ncid11,var_swiltid,'units',' m3/m3')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining swilt')
     ok = nf90_put_att(ncid11,var_swiltid,'long_name',' soil water at wilting point')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining swilt')
     ok = nf90_put_att(ncid11,var_swiltid,'_FillValue',missreal)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining swilt')
     ok = nf90_put_att(ncid11,var_swiltid,'missing_value',missreal)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining swilt')

     ok = nf90_def_var(ncid11,'iveg',nf90_int,(/dim_mlonid,dim_mlatid,dim_patchid/),var_ivegid)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining iveg')
     ok = nf90_put_att(ncid11,var_ivegid,'units',' class')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining iveg')
     ok = nf90_put_att(ncid11,var_ivegid,'long_name',' CABLE PFT')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining iveg')
     ok = nf90_put_att(ncid11,var_ivegid,'_FillValue',missint)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining iveg')
     ok = nf90_put_att(ncid11,var_ivegid,'missing_value',missint)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining iveg')

     ok = nf90_def_var(ncid11,'patchfrac',nf90_float,(/dim_mlonid,dim_mlatid,dim_patchid/),var_patchfracid)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining patchfrac')
     ok = nf90_put_att(ncid11,var_patchfracid,'units',' area fraction')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining patchfrac')
     ok = nf90_put_att(ncid11,var_patchfracid,'long_name',' CABLE PFT fraction')
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining patchfrac')
     ok = nf90_put_att(ncid11,var_patchfracid,'_FillValue',missreal)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining patchfrac')
     ok = nf90_put_att(ncid11,var_patchfracid,'missing_value',missreal)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in defining patchfrac')

     ok = nf90_enddef(ncid11)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in enddef')

     ! assign values

     call rangechk3(mlon,mlat,nrb,landmask,albedo_y,0.01,0.9)
     ok = nf90_put_var(ncid11,var_albedoid,albedo_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in put albedo')

     call rangechk3(mlon,mlat,time12,landmask,lai_y,0.0,10.0)
     ok = nf90_put_var(ncid11,var_laiid,lai_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in put lai')

     call rangechk2(mlon,mlat,landmask,ndep_y,0.0,10.0)
     ok = nf90_put_var(ncid11,var_ndepid,ndep_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in put ndep')

     call rangechk2(mlon,mlat,landmask,nfix_y,0.0,15.0)
     ok = nf90_put_var(ncid11,var_nfixid,nfix_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in put nfix')

     call rangechk2(mlon,mlat,landmask,pdust_y,0.0,5.0)
     ok = nf90_put_var(ncid11,var_pdustid,pdust_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in put pdep')

     call rangechk2(mlon,mlat,landmask,pwea_y,0.0,5.0)
     ok = nf90_put_var(ncid11,var_pweaid,pwea_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in put pwea')

     ok = nf90_put_var(ncid11,var_snowdepthid,snowdepth_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in put snowdepth')

     ok = nf90_put_var(ncid11,var_soilmoistid,soilmoist_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in put soilmoist')

     ok = nf90_put_var(ncid11,var_soilorderid,soilorder_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in put soilorder')

     ok = nf90_put_var(ncid11,var_soiltempid,soiltemp_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in put soiltemp')

     call rangechk2(mlon,mlat,landmask,albedo2_y,0.01,0.9)
     ok = nf90_put_var(ncid11,var_albedo2id,albedo2_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in put albedo2')

     ok = nf90_put_var(ncid11,var_areaid,area_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in put area')

     call rangechk2(mlon,mlat,landmask,bch_y,2.0,15.0)
     ok = nf90_put_var(ncid11,var_bchid,bch_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in put bch')

     call rangechk2(mlon,mlat,landmask,clay_y,0.0,1.0)
     ok = nf90_put_var(ncid11,var_clayid,clay_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in put clay')

     ok = nf90_put_var(ncid11,var_cnsdid,cnsd_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in put cnsd')

     call rangechk2(mlon,mlat,landmask,css_y,700.0,2200.0)
     ok = nf90_put_var(ncid11,var_cssid,css_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in put css')

     call rangechk2(mlon,mlat,landmask,hyds_y,5.0e-7,8.5e-3)
     ok = nf90_put_var(ncid11,var_hydsid,hyds_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in put hyds')

     ok = nf90_put_var(ncid11,var_isoilid,isoil_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in put isoil')

     ok = nf90_put_var(ncid11,var_latitudeid,latitude_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in put latitude')

     ok = nf90_put_var(ncid11,var_longitudeid,longitude_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in put longitude')

     call rangechk2(mlon,mlat,landmask,rhosoil_y,300.0,3000.0)
     ok = nf90_put_var(ncid11,var_rhosoilid,rhosoil_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in put rhosoil')

     call rangechk2(mlon,mlat,landmask,sand_y,0.0,1.0)
     ok = nf90_put_var(ncid11,var_sandid,sand_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in put sand')

     call rangechk2(mlon,mlat,landmask,sfc_y,0.1,0.5)
     ok = nf90_put_var(ncid11,var_sfcid,sfc_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in put sfc')

     call rangechk2(mlon,mlat,landmask,silt_y,0.0,1.0)
     ok = nf90_put_var(ncid11,var_siltid,silt_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in put silt')

     call rangechk2(mlon,mlat,landmask,ssat_y,0.35,0.50)
     ok = nf90_put_var(ncid11,var_ssatid,ssat_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in put ssat')

     call rangechk2(mlon,mlat,landmask,sucs_y,-0.8,-0.03)
     ok = nf90_put_var(ncid11,var_sucsid,sucs_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in put sucs')

     call rangechk2(mlon,mlat,landmask,swilt_y,0.05,0.4)
     ok = nf90_put_var(ncid11,var_swiltid,swilt_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in put swilt')

     ok = nf90_put_var(ncid11,var_ivegid,iveg_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in put iveg')

     ok = nf90_put_var(ncid11,var_patchfracid,patchfrac_y)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in put patchfrac')

     ok = nf90_close(ncid11)
     if(ok/=nf90_noerr) call nc_abort(ok, 'error in put albedo')

     write(logn,*) 'landuse on: new gridinfo created', fgridnew
  end subroutine create_new_gridinfo

  subroutine rangechk2(mlon,mlat,landmask,varx2,xmin,xmax)
  USE cable_def_types_mod,  ONLY: r_2
  implicit none
  real xmin,xmax
  integer mlon,mlat
  integer, dimension(mlon,mlat)   :: landmask
  real(r_2),    dimension(mlon,mlat)   :: varx2
  integer i,j

    do i=1,mlon
    do j=1,mlat
       if(landmask(i,j)==1) then
          varx2(i,j) = max(xmin,varx2(i,j))
          varx2(i,j) = min(xmax,varx2(i,j))
       endif
    enddo
    enddo
  end subroutine rangechk2


  subroutine rangechk3(mlon,mlat,nx3,landmask,varx3,xmin,xmax)
  USE cable_def_types_mod,  ONLY: r_2
  implicit none
  real xmin,xmax
  integer mlon,mlat,nx3
  integer, dimension(mlon,mlat)       :: landmask
  real(r_2),    dimension(mlon,mlat,nx3)   :: varx3
  integer i,j,k

    do i=1,mlon
    do j=1,mlat
       if(landmask(i,j)==1) then
          do k=1,nx3
             varx3(i,j,k) = max(xmin,varx3(i,j,k))
             varx3(i,j,k) = min(xmax,varx3(i,j,k))
          enddo
       endif
    enddo
    enddo
  end subroutine rangechk3



  subroutine sort(mvmax,tmpx,tmpint)
  ! based on numerical recipes, straight insertion method  p322
  USE cable_def_types_mod,  ONLY: r_2
  implicit none
  integer mvmax
  integer,      dimension(mvmax)  :: tmpint
  real(r_2),    dimension(mvmax)  :: tmpx
  integer i, j, na
  real(r_2)     xa

    do j=2,mvmax
       xa = tmpx(j)
       na = tmpint(j)
       do i=j-1,1,-1
          if(tmpx(i)>=xa) go to 10
            tmpx(i+1) = tmpx(i)
            tmpint(i+1) = tmpint(i)
       enddo
       i=0
10     tmpx(i+1)   = xa
       tmpint(i+1) = na
    enddo

  end subroutine sort

  SUBROUTINE WRITE_LANDUSE_CASA_RESTART_NC(mpx, lucmp, CASAONLY )
    ! if any additional variables are added to CASA restart file, they should be repeated here
    !

    USE netcdf
    USE casavariable,         ONLY : icycle, mplant, mlitter, msoil, mwood, casafile
    USE cable_IO_vars_module, ONLY : logn
    USE cable_common_module
    USE casa_ncdf_module,     ONLY: HANDLE_ERR
    USE landuse_variable,     ONLY: landuse_mp

    IMPLICIT NONE
    type(landuse_mp)         :: lucmp

    INTEGER,      INTENT(IN) :: mpx
    INTEGER*4                :: mp4
    INTEGER*4, PARAMETER     :: pmp4 =0
    INTEGER, PARAMETER       :: fmp4 = KIND(pmp4)
    INTEGER*4                :: STATUS
    INTEGER*4                :: FILE_ID, land_ID, plnt_ID, litt_ID, soil_ID, wood_ID, i
    LOGICAL                  :: CASAONLY
    CHARACTER                :: CYEAR*4, FNAME*99,dum*50

    ! ! 1 dim arrays (npt )
    ! CHARACTER(len=20),DIMENSION(7), PARAMETER :: A1 = (/ 'latitude', 'longitude', 'glai', &
    !      'clabile', 'psoillab','psoilsorb','psoilocc' /)
    ! ! 2 dim arrays (npt,mplant)
    ! CHARACTER(len=20),DIMENSION(3), PARAMETER :: A2 = (/ 'cplant' , 'nplant' , 'pplantc' /)
    ! ! 2 dim arrays (npt,mlitter)
    ! CHARACTER(len=20),DIMENSION(3), PARAMETER :: A3 = (/ 'clitter', 'nlitter', 'plitter' /)
    ! ! 2 dim arrays (npt,msoil)
    ! CHARACTER(len=20),DIMENSION(3), PARAMETER :: A4 = (/ 'csoil', 'nsoil', 'psoil' /)

    ! 1 dim arrays (npt )
    CHARACTER(len=20),DIMENSION(12) :: A1
    CHARACTER(len=20),DIMENSION(2)  :: AI1
    ! 2 dim arrays (npt,mplant)
    CHARACTER(len=20),DIMENSION(3)  :: A2
    ! 2 dim arrays (npt,mlitter)
    CHARACTER(len=20),DIMENSION(3)  :: A3
    ! 2 dim arrays (npt,msoil)
    CHARACTER(len=20),DIMENSION(3)  :: A4
    ! 2 dim arrays (npt,msoil)
    CHARACTER(len=20),DIMENSION(3)  :: A5

    INTEGER*4 :: VID1(SIZE(A1)), VIDI1(SIZE(AI1)), VID2(SIZE(A2)), &
                 VID3(SIZE(A3)), VID4(SIZE(A4)),VID5(SIZE(A5))

    mp4=INT(mpx,fmp4)
    write(logn,*)  ' landuse on: writing casa pool: patch number=', mpx,fmp4,mp4

    A1(1) = 'latitude'
    A1(2) = 'longitude'
    A1(3) = 'glai'
    A1(4) = 'clabile'
    A1(5) = 'psoillab'
    A1(6) = 'psoilsorb'
    A1(7) = 'psoilocc'
    A1(8) = 'frac_sapwood'
    A1(9) = 'sapwood_area'
    A1(10) = 'phen'
    A1(11) = 'aphen'
    A1(12) = 'nsoilmin'

    AI1(1) = 'phase'
    AI1(2) = 'doyphase3'


    A2(1) = 'cplant'
    A2(2) = 'nplant'
    A2(3) = 'pplant'

    A3(1) = 'clitter'
    A3(2) = 'nlitter'
    A3(3) = 'plitter'

    A4(1) = 'csoil'
    A4(2) = 'nsoil'
    A4(3) = 'psoil'

    A5(1) = 'cwoodprod'
    A5(2) = 'nwoodprod'
    A5(3) = 'pwoodprod'

    ! Get File-Name
    WRITE(CYEAR, FMT='(I4)') CurYear + 1

    IF (LEN( TRIM(casafile%cnpepool) ) .GT. 0) THEN
       fname=TRIM(casafile%cnpepool)
    ELSE
       fname = TRIM(filename%path)//'/'//TRIM( cable_user%RunIden )//&
            '_casa_rst.nc'
    ENDIF
    ! Create NetCDF file:
    STATUS = NF90_create(fname, NF90_CLOBBER, FILE_ID)
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    WRITE(*,*) 'writing casa restart', fname
    ! Put the file in define mode:
    STATUS = NF90_redef(FILE_ID)

    STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "Valid restart date", "01/01/"//CYEAR  )
    STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "Icycle", icycle  )
    IF ( CASAONLY ) THEN
       dum = 'CASA-ONLY run'
    ELSE
       dum = 'CABLE-CASA coupled run'
    ENDIF
    STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "Run-Type", TRIM(dum) )

    ! Define dimensions:
    ! Land (number of points)
    STATUS = NF90_def_dim(FILE_ID, 'land'   , mp4     , land_ID)
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    STATUS = NF90_def_dim(FILE_ID, 'mplant' , mplant , plnt_ID)
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    STATUS = NF90_def_dim(FILE_ID, 'mlitter', mlitter, litt_ID)
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    STATUS = NF90_def_dim(FILE_ID, 'msoil'  , msoil  , soil_ID)
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

    STATUS = NF90_def_dim(FILE_ID, 'mwood'  , mwood  , wood_ID)
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

    DO i = 1, SIZE(A1)
       STATUS = NF90_def_var(FILE_ID,TRIM(A1(i)) ,NF90_FLOAT,(/land_ID/),VID1(i))
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    END DO

    DO i = 1, SIZE(AI1)
       STATUS = NF90_def_var(FILE_ID,TRIM(AI1(i)) ,NF90_INT,(/land_ID/),VIDI1(i))
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    END DO

    DO i = 1, SIZE(A2)
       STATUS = NF90_def_var(FILE_ID,TRIM(A2(i)) ,NF90_FLOAT,(/land_ID,plnt_ID/),VID2(i))
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    END DO

    DO i = 1, SIZE(A3)
       STATUS = NF90_def_var(FILE_ID,TRIM(A3(i)) ,NF90_FLOAT,(/land_ID,litt_ID/),VID3(i))
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    END DO

    DO i = 1, SIZE(A4)
       STATUS = NF90_def_var(FILE_ID,TRIM(A4(i)) ,NF90_FLOAT,(/land_ID,soil_ID/),VID4(i))
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    END DO

    DO i = 1, SIZE(A5)
       STATUS = NF90_def_var(FILE_ID,TRIM(A5(i)) ,NF90_FLOAT,(/land_ID,wood_ID/),VID5(i))
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    END DO

    ! End define mode:
    STATUS = NF90_enddef(FILE_ID)
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

    ! PUT LAT / LON
    STATUS = NF90_PUT_VAR(FILE_ID, VID1(1), lucmp%lat )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1(2), lucmp%lon )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    ! PUT VARS
    STATUS = NF90_PUT_VAR(FILE_ID, VID1(3), lucmp%lai )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1(4), lucmp%clabile )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1(5), lucmp%psoillab )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1(6), lucmp%psoilsorb )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1(7), lucmp%psoilocc )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1(8), lucmp%frac_sapwood )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1(9), lucmp%sapwood_area )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1(10), lucmp%phen )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1(11), lucmp%aphen )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1(12), lucmp%Nsoilmin )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VIDI1(1), lucmp%phase )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VIDI1(2), lucmp%doyphase3 )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID2(1), lucmp%cplant  )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID2(2), lucmp%nplant  )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID2(3), lucmp%pplant  )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID3(1), lucmp%clitter  )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID3(2), lucmp%nlitter )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID3(3), lucmp%plitter )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID4(1), lucmp%csoil )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID4(2), lucmp%nsoil )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID4(3), lucmp%psoil )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID5(1), lucmp%cwoodprod )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID5(2), lucmp%nwoodprod )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID4(3), lucmp%pwoodprod )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)


    ! Close NetCDF file:
    STATUS = NF90_close(FILE_ID)
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

    write(logn, *) 'landuse on: casapool writeen to ', fname
  END SUBROUTINE WRITE_LANDUSE_CASA_RESTART_NC


  SUBROUTINE create_landuse_cable_restart(logn,dels,ktau,soil,mpx,lucmp,cstart,cend,nap, met)
    ! Creates a restart file for CABLE using a land only grid cell area occupied by a '//  &
    ! Creates a restart file for CABLE using a land only grid with mland
    ! land points and max_vegpatches veg/soil patches (some of which may
    ! not be active). It uses CABLE's internal variable names.
    ! be aware
    ! if any additional variables are added to CABLE restart file, they should be repeated here
    !
    use netcdf
    USE cable_def_types_mod,        ONLY : r_2, mland, mvtype, mstype,nrb,ncs,ncp,ms,msn,soil_parameter_type, met_type
    use cable_abort_module,         ONLY : nc_abort
    USE cable_IO_vars_module,       ONLY : latitude,longitude,timeunits,calendar,time_coord, timevar
    USE cable_checks_module,        ONLY : ranges
    USE cable_common_module,        ONLY : filename,CurYear,cable_user
    USE landuse_variable,           ONLY : landuse_mp

    implicit none
    type(landuse_mp)                       :: lucmp
    type(soil_parameter_type)              :: soil

    INTEGER, INTENT(IN)                    :: logn ! log file number
    REAL,    INTENT(IN)                    :: dels ! time step size
    INTEGER, INTENT(IN)                    :: ktau ! timestep number in loop which include spinup
    INTEGER, INTENT(IN)                    :: mpx  ! timestep number in loop which include spinup
    INTEGER, DIMENSION(mland), INTENT(in)  :: cstart,cend,nap

    TYPE(met_type), TARGET, INTENT(IN)                         :: met
!    TYPE (soil_parameter_type),INTENT(IN)  :: soil    ! from "cable_de_types_mod"
!    TYPE (ranges_type),        INTENT(IN)  :: ranges  ! from "cable_checks_module"

    INTEGER           :: ncid_restart ! netcdf restart file ID
    ! REAL, POINTER,DIMENSION(:,:) :: surffrac ! fraction of each surf type
    INTEGER           :: mlandID, mpID, radID, soilID, napID,                       &
                         soilcarbID, plantcarbID, tID, snowID ! dimension IDs

    INTEGER           :: patchfrac_id,mvtype_id,mstype_id
    INTEGER           :: iveg_id, isoil_id, zse_id,albsoil_id
    INTEGER           :: tvarID, latID, lonID !,surffracID ! time,lat,lon variable ID
    INTEGER           :: tggID, wbID, wbiceID, tssID, ssdnnID, ssdnID, osnowdID,    &
                         smassID, sdepthID, snageID, snowdID, rtsoilID, isflagID,   &
                         canstoID, albsoilsnID, gammzzID, tggsnID, sghfluxID,       &
                         ghfluxID, runoffID, rnof1ID, rnof2ID, gaID, dgdtgID,       &
                         fevID, fesID, fhsID, wbtot0ID, osnowd0ID, cplantID,        &
                         csoilID, tradID, albedoID, gwID

    INTEGER           :: h0ID, snowliqID, SID, TsurfaceID, scondsID, nsnowID, TsoilID
    CHARACTER(LEN=10) :: todaydate, nowtime ! used to timestamp netcdf file

    CHARACTER         :: FRST_OUT*200, CYEAR*4

    INTEGER              ok

    WRITE(logn, '(A24)') ' Writing restart file...'
    IF ( TRIM(filename%path) .EQ. '' ) filename%path = './'
    frst_out = TRIM(filename%path)//'/'//TRIM(filename%restart_out)
    ! Look for explicit restart file (netCDF). If not, asssume input is path
    IF ( INDEX(TRIM(frst_out),'.nc',BACK=.TRUE.) .NE. LEN_TRIM(frst_out)-2 ) THEN
       WRITE( CYEAR,FMT="(I4)" ) CurYear + 1
       frst_out = TRIM(filename%path)//'/'//TRIM(cable_user%RunIden)//&
            '_'//CYEAR//'_cable_rst.nc'
    ENDIF

    ! Create output file:
    ok = NF90_CREATE(frst_out, NF90_CLOBBER, ncid_restart)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error creating restart file '      &
         //TRIM(frst_out)// '(SUBROUTINE create_landuse_cable_restart)')

    ! Put the file in define mode:
    ok = NF90_REDEF(ncid_restart)
    ! Define dimensions:
    ok = NF90_DEF_DIM(ncid_restart, 'mland', mland, mlandID)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining mland dimension in restart file. '// &
         '(SUBROUTINE create_landuse_cable_restart)')

    ok = NF90_DEF_DIM(ncid_restart, 'mp', mpx, mpID)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining mp dimension in restart file. '// &
         '(SUBROUTINE create_landuse_cable_restart)')

    ok = NF90_DEF_DIM(ncid_restart, 'soil', ms, soilID) ! number of soil layers
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining vertical soil dimension in restart file. '// &
         '(SUBROUTINE create_landuse_cable_restart)')

    ok = NF90_DEF_DIM(ncid_restart, 'snow', 3, snowID) ! number of snow layers
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining vertical snow dimension in restart file. '// &
         '(SUBROUTINE create_landuse_cable_restart)')

    ok = NF90_DEF_DIM(ncid_restart, 'rad', nrb, radID) ! number of rad. bands
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining radiation dimension in restart file. '// &
         '(SUBROUTINE create_landuse_cable_restart)')

    ok = NF90_DEF_DIM(ncid_restart, 'soil_carbon_pools', ncs, soilcarbID)
    ! number of soil carbon pools
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining soil carbon pool dimension in restart file. '// &
         '(SUBROUTINE create_landuse_cable_restart)')

    ok = NF90_DEF_DIM(ncid_restart, 'plant_carbon_pools', ncp, plantcarbID)
    ! number of plant carbon pools
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining plant carbon pool dimension in restart file. '// &
         '(SUBROUTINE create_landuse_cable_restart)')

    ok = NF90_DEF_DIM(ncid_restart, 'time', 1, tID)
    IF (ok /= NF90_NOERR) CALL nc_abort &
         (ok, 'Error defining time dimension in restart file. '// &
         '(SUBROUTINE create_landuse_cable_restart)')

    ! Define "time" variable and its attributes:
    ok=NF90_DEF_VAR(ncid_restart,'time',NF90_DOUBLE,(/tID/),tvarID)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining time variable in restart file. '// &
         '(SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, tvarID, 'units', timeunits)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining time variable attributes in restart file. '// &
         '(SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, tvarID, 'coordinate', time_coord)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining time variable attributes in restart file. '// &
         '(SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, tvarID, 'calendar', calendar)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining time variable attribute calendar in restart file. '// &
         '(SUBROUTINE create_landuse_cable_restart)')

    ! Define latitude and longitude variable:
    ok=NF90_DEF_VAR(ncid_restart, 'latitude', NF90_FLOAT, (/mlandID/), latID)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining latitude variable in restart file. '// &
         '(SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart,latID,'units','degrees_north')
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining latitude variable attributes in restart file. '// &
         '(SUBROUTINE create_landuse_cable_restart)')

    ok=NF90_DEF_VAR(ncid_restart, 'longitude', NF90_FLOAT, (/mlandID/), lonID)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining longitude variable in restart file. '// &
         '(SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, lonID, 'units', 'degrees_east')
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining longitude variable attributes in restart file. '// &
         '(SUBROUTINE create_landuse_cable_restart)')
    ! Define number of active patches variable:
    ok = NF90_DEF_VAR(ncid_restart, 'nap', NF90_FLOAT, (/mlandID/), napID)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining nap variable in restart file. '// &
         '(SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, napID, 'long_name',                        &
         'Number of active patches')
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok,'Error defining nap variable attributes in restart file. '// &
         '(SUBROUTINE create_landuse_cable_restart)')

    ! Define patch fraction variable:
    ok=NF90_DEF_VAR(ncid_restart, 'patchfrac', NF90_FLOAT, (/mpID/),           &
         patchfrac_id)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining patchfrac variable in restart file. '// &
         '(SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, patchfrac_id, 'long_name',               &
         'Fraction of vegetated grid cell area occupied by a '//  &
         'vegetation/soil patch')
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining patchfrac variable attributes in restart file. '// &
         '(SUBROUTINE create_landuse_cable_restart)')

    ! mvtype (Number of vegetation types):
    ok = NF90_DEF_VAR(ncid_restart, 'mvtype', NF90_INT, mvtype_id)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining mvtype variable in restart file. '// &
         '(SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, mvtype_id, "long_name",                  &
         "Number of vegetation types")
    ! mstype (Number of soil types):
    ok = NF90_DEF_VAR(ncid_restart, 'mstype', NF90_INT, mstype_id)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining mstype variable in restart file. '// &
         '(SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, mstype_id, "long_name",                  &
         "Number of soil types")

    !------------------define soil states---------------------------------------

    ok = NF90_DEF_VAR(ncid_restart, 'tgg', NF90_FLOAT, (/mpID, soilID/), tggID)
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining tgg variable in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, tggID, 'long_name', 'Average layer soil temperature')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining tgg variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, tggID, 'units', 'K')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining tgg variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')

    ok = NF90_DEF_VAR(ncid_restart, 'wb', NF90_DOUBLE, (/mpID, soilID/), wbID)
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining wb variable in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, wbID, 'long_name', 'Average layer volumetric soil moisture')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining wb variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, wbID, 'units', 'vol/vol')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining wb variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')

    ok = NF90_DEF_VAR(ncid_restart, 'wbice', NF90_DOUBLE, (/mpID, soilID/), wbiceID)
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining wbice variable in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, wbiceID, 'long_name', 'Average layer volumetric soil ice')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining wbice variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, wbiceID, 'units', 'vol/vol')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining wbice variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')

    ok = NF90_DEF_VAR(ncid_restart, 'tss', NF90_FLOAT, (/mpID/), tssID)
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining tss variable in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, tssID, 'long_name', 'Combined soil/snow temperature')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining tss variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, tssID, 'units', 'K')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining tss variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')

    ok = NF90_DEF_VAR(ncid_restart, 'albsoilsn', NF90_FLOAT, (/mpID, radID/), albsoilsnID)
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining albsoilsn variable in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, albsoilsnID, 'long_name', 'Combined soil/snow albedo')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining albsoilsn variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, albsoilsnID, 'units', '-')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining albsoilsn variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')

    ok = NF90_DEF_VAR(ncid_restart, 'rtsoil', NF90_FLOAT, (/mpID/), rtsoilID)
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining rtsoil variable in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, rtsoilID, 'long_name', 'Turbulent resistance for soil')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining rtsoil variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, rtsoilID, 'units', '??')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining rtsoil variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')

    ok = NF90_DEF_VAR(ncid_restart, 'gammzz', NF90_DOUBLE, (/mpID, soilID/), gammzzID)
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining gammzz variable in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, gammzzID, 'long_name', 'Heat capacity for each soil layer')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining gammzz variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, gammzzID, 'units', 'J/kg/C')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining gammzz variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')

    ok = NF90_DEF_VAR(ncid_restart, 'runoff', NF90_FLOAT, (/mpID/), runoffID)
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining runoff variable in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, runoffID, 'long_name', 'Total runoff')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining runoff variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, runoffID, 'units', 'mm/timestep')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining runoff variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')

    ok = NF90_DEF_VAR(ncid_restart, 'rnof1', NF90_FLOAT, (/mpID/), rnof1ID)
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining rnof1 variable in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, rnof1ID, 'long_name', 'Surface runoff')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining rnof1 variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, rnof1ID, 'units', 'mm/timestep')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining rnof1 variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')

    ok = NF90_DEF_VAR(ncid_restart, 'rnof2', NF90_FLOAT, (/mpID/), rnof2ID)
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining rnof2 variable in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, rnof2ID, 'long_name', 'Subsurface runoff')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining rnof2 variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, rnof2ID, 'units', 'mm/timestep')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining rnof2 variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')

    !---------------define snow states------------------------------------------

    ok = NF90_DEF_VAR(ncid_restart, 'tggsn', NF90_FLOAT, (/mpID, snowID/), tggsnID)
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining tggsn variable in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, tggsnID, 'long_name', 'Average layer snow temperature')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining tggsn variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, tggsnID, 'units', 'K')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining tggsn variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')

    ok = NF90_DEF_VAR(ncid_restart, 'ssdnn', NF90_FLOAT, (/mpID/), ssdnnID)
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining ssdnn variable in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, ssdnnID, 'long_name', 'Average snow density')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining ssdnn variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, ssdnnID, 'units', 'kg/m^3')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining ssdnn variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')

    ok = NF90_DEF_VAR(ncid_restart, 'ssdn', NF90_FLOAT, (/mpID, snowID/), ssdnID)
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining ssdn variable in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, ssdnID, 'long_name', 'Average layer snow density')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining ssdn variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, ssdnID, 'units', 'kg/m^3')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining ssdn variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')

    ok = NF90_DEF_VAR(ncid_restart, 'snowd', NF90_FLOAT, (/mpID/), snowdID)
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining snowd variable in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, snowdID, 'long_name', 'Liquid water equivalent snow depth')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining snowd variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, snowdID, 'units', 'mm')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining snowd variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')

    ok = NF90_DEF_VAR(ncid_restart, 'snage', NF90_FLOAT, (/mpID/), snageID)
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining snage variable in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, snageID, 'long_name', 'Snow age')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining snage variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, snageID, 'units', '??')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining snage variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')

    ok = NF90_DEF_VAR(ncid_restart, 'smass', NF90_FLOAT, (/mpID, snowID/), smassID)
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining smass variable in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, smassID, 'long_name', 'Average layer snow mass')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining smass variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, smassID, 'units', 'kg/m^2')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining smass variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')

    ok = NF90_DEF_VAR(ncid_restart, 'sdepth', NF90_FLOAT, (/mpID, snowID/), sdepthID)
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining sdepth variable in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, sdepthID, 'long_name', 'Snow layer depth')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining sdepth variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, sdepthID, 'units', 'm')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining sdepth variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')

    ok = NF90_DEF_VAR(ncid_restart, 'osnowd', NF90_FLOAT, (/mpID/), osnowdID)
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining osnowd variable in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, osnowdID, 'long_name', 'Previous time step snow depth in water equivalent')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining osnowd variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, osnowdID, 'units', 'mm')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining osnowd variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')

    ok = NF90_DEF_VAR(ncid_restart, 'isflag', NF90_INT, (/mpID/), isflagID)
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining isflag variable in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, isflagID, 'long_name', 'Snow layer scheme flag')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining isflag variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, isflagID, 'units', '-')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining isflag variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')

    !----------------define canopy states----------------------------------
    ok = NF90_DEF_VAR(ncid_restart, 'cansto', NF90_FLOAT, (/mpID/), canstoID)
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining cansto variable in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, canstoID, 'long_name', 'Canopy surface water storage')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining cansto variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, canstoID, 'units', 'mm')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining cansto variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')

    ok = NF90_DEF_VAR(ncid_restart, 'ghflux', NF90_FLOAT, (/mpID/), ghfluxID)
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining ghflux variable in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, ghfluxID, 'long_name', '????')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining ghflux variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, ghfluxID, 'units', 'W/m^2')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining ghflux variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')

    ok = NF90_DEF_VAR(ncid_restart, 'sghflux', NF90_FLOAT, (/mpID/), sghfluxID)
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining sghflux variable in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, sghfluxID, 'long_name', '????')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining sghflux variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, sghfluxID, 'units', 'W/m^2')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining sghflux variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')

    ok = NF90_DEF_VAR(ncid_restart, 'ga', NF90_FLOAT, (/mpID/), gaID)
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining ga variable in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, gaID, 'long_name', 'Ground heat flux')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining ga variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, gaID, 'units', 'W/m^2')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining ga variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')

    ok = NF90_DEF_VAR(ncid_restart, 'dgdtg', NF90_DOUBLE, (/mpID/), dgdtgID)
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining dgdtg variable in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, dgdtgID, 'long_name', 'Derivative of ground heat flux wrt soil temperature')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining dgdtg variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, dgdtgID, 'units', 'W/m^2/K')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining dgdtg variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')

    ok = NF90_DEF_VAR(ncid_restart, 'fev', NF90_FLOAT, (/mpID/), fevID)
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining fev variable in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, fevID, 'long_name', 'Latent heat flux from vegetation')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining fev variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, fevID, 'units', 'W/m^2')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining fev variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')

    ok = NF90_DEF_VAR(ncid_restart, 'fes', NF90_FLOAT, (/mpID/), fesID)
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining fes variable in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, fesID, 'long_name', 'Latent heat flux from soil')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining fes variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, fesID, 'units', 'W/m^2')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining fes variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')

    ok = NF90_DEF_VAR(ncid_restart, 'fhs', NF90_FLOAT, (/mpID/), fhsID)
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining fhs variable in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, fhsID, 'long_name', 'Sensible heat flux from soil')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining fhs variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, fhsID, 'units', 'W/m^2')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining fhs variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')

    !--------------biogeochemical variables------------------------

    ok = NF90_DEF_VAR(ncid_restart, 'cplant', NF90_FLOAT, (/mpID, plantcarbID/), cplantID)
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining cplant variable in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, cplantID, 'long_name', 'Plant carbon stores')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining cplant variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, cplantID, 'units', 'gC/m^2')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining cplant variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')

    ok = NF90_DEF_VAR(ncid_restart, 'csoil', NF90_FLOAT, (/mpID, soilcarbID/), csoilID)
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining csoil variable in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, csoilID, 'long_name', 'Soil carbon stores')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining csoil variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, csoilID, 'units', 'gC/m^2')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining csoil variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')

    !-------------------others---------------------------------

    ok = NF90_DEF_VAR(ncid_restart, 'wbtot0', NF90_FLOAT, (/mpID/), wbtot0ID)
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining wbtot0 variable in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, wbtot0ID, 'long_name', 'Initial time step soil water total')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining wbtot0 variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, wbtot0ID, 'units', 'mm')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining wbtot0 variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')

    ok = NF90_DEF_VAR(ncid_restart, 'osnowd0', NF90_FLOAT, (/mpID/), osnowd0ID)
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining osnowd0 variable in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, osnowd0ID, 'long_name', 'Initial time step snow water total')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining osnowd0 variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, osnowd0ID, 'units', 'mm')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining osnowd0 variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')

    ok = NF90_DEF_VAR(ncid_restart, 'albedo', NF90_FLOAT, (/mpID, radID/), albedoID)
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining albedo variable in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, albedoID, 'long_name', 'Albedo for shortwave and NIR radiation')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining albedo variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, albedoID, 'units', '-')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining albedo variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')

    ok = NF90_DEF_VAR(ncid_restart, 'trad', NF90_FLOAT, (/mpID/), tradID)
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining trad variable in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, tradID, 'long_name', 'Surface radiative temperature (soil/snow/veg inclusive)')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining trad variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, tradID, 'units', 'K')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining trad variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')

    !---------------------MODEL PARAMETERS---------------------------------
    WRITE(logn,'(A43)') '   Writing model parameters to restart file'

    ok = NF90_DEF_VAR(ncid_restart, 'iveg', NF90_INT, (/mpID/), iveg_id)
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining iveg variable in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, iveg_id, 'long_name', 'Vegetation type')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining iveg variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, iveg_id, 'units', '-')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining iveg variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')

    ok = NF90_DEF_VAR(ncid_restart, 'isoil', NF90_INT, (/mpID/), isoil_id)
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining isoil variable in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, isoil_id, 'long_name', 'Soil type')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining isoil variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, isoil_id, 'units', '-')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining isoil variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')

    ! zse (depth of each soil layer):
    ok = NF90_DEF_VAR(ncid_restart, 'zse', NF90_FLOAT, (/soilID/), zse_id)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining zse variable in restart file. '// &
         '(SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, zse_id, "long_name",                     &
         "Depth of each soil layer")
    ok = NF90_PUT_ATT(ncid_restart, zse_id, "units", "m")

    ok = NF90_DEF_VAR(ncid_restart, 'albsoil', NF90_FLOAT, (/mpID, radID/), albsoil_id)
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining albsoil variable in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, albsoil_id, 'long_name', 'Soil reflectance')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining albsoil variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, albsoil_id, 'units', '-')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining albsoil variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')

    ok = NF90_DEF_VAR(ncid_restart, 'GWwb', NF90_FLOAT, (/mpID/), gwID)
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining GWwb variable in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, gwID, 'long_name', 'GW water content')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining GWwb variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, gwID, 'units', 'mm3/mm3')
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining GWwb variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')

    ! Soil-Litter-Iso soil model
    IF(cable_user%SOIL_STRUC=='sli') THEN
       ! Parameters for SLI:
       ok = NF90_DEF_VAR(ncid_restart, 'S', NF90_FLOAT, (/mpID, soilID/), SID)
       if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining S variable in restart file. (SUBROUTINE create_landuse_cable_restart)')
       ok = NF90_PUT_ATT(ncid_restart, SID, 'long_name', 'Fractional soil moisture content relative to saturated value')
       if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining S variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')
       ok = NF90_PUT_ATT(ncid_restart, SID, 'units', '-')
       if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining S variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')

       ok = NF90_DEF_VAR(ncid_restart, 'Tsoil', NF90_FLOAT, (/mpID, soilID/), TsoilID)
       if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining Tsoil variable in restart file. (SUBROUTINE create_landuse_cable_restart)')
       ok = NF90_PUT_ATT(ncid_restart, TsoilID, 'long_name', 'Tsoil')
       if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining Tsoil variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')
       ok = NF90_PUT_ATT(ncid_restart, TsoilID, 'units', 'degC')
       if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining Tsoil variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')

       ok = NF90_DEF_VAR(ncid_restart, 'snowliq', NF90_FLOAT, (/mpID, snowID/), snowliqID)
       if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining snowliq variable in restart file. (SUBROUTINE create_landuse_cable_restart)')
       ok = NF90_PUT_ATT(ncid_restart, snowliqID, 'long_name', 'liquid water content of snowpack')
       if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining snowliq variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')
       ok = NF90_PUT_ATT(ncid_restart, snowliqID, 'units', 'mm')
       if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining snowliq variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')

       ok = NF90_DEF_VAR(ncid_restart, 'sconds', NF90_FLOAT, (/mpID, snowID/), scondsID)
       if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining sconds variable in restart file. (SUBROUTINE create_landuse_cable_restart)')
       ok = NF90_PUT_ATT(ncid_restart, scondsID, 'long_name', 'thermal cond of snowpack')
       if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining sconds variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')
       ok = NF90_PUT_ATT(ncid_restart, scondsID, 'units', 'Wm-1K-1')
       if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining sconds variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')

       ok = NF90_DEF_VAR(ncid_restart, 'h0', NF90_FLOAT, (/mpID/), h0ID)
       if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining h0 variable in restart file. (SUBROUTINE create_landuse_cable_restart)')
       ok = NF90_PUT_ATT(ncid_restart, h0ID, 'long_name', 'Pond height above soil')
       if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining h0 variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')
       ok = NF90_PUT_ATT(ncid_restart, h0ID, 'units', 'm')
       if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining h0 variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')

       ok = NF90_DEF_VAR(ncid_restart, 'nsnow', NF90_INT, (/mpID/), nsnowID)
       if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining nsnow variable in restart file. (SUBROUTINE create_landuse_cable_restart)')
       ok = NF90_PUT_ATT(ncid_restart, nsnowID, 'long_name', 'number of snow layers')
       if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining nsnow variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')
       ok = NF90_PUT_ATT(ncid_restart, nsnowID, 'units', '-')
       if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining nsnow variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')

       ok = NF90_DEF_VAR(ncid_restart, 'Tsurface', NF90_FLOAT, (/mpID/), TsurfaceID)
       if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining Tsurface variable in restart file. (SUBROUTINE create_landuse_cable_restart)')
       ok = NF90_PUT_ATT(ncid_restart, TsurfaceID, 'long_name', 'soil or snow surface T')
       if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining Tsurface variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')
       ok = NF90_PUT_ATT(ncid_restart, TsurfaceID, 'units', 'degC')
       if (ok /= NF90_NOERR) call nc_abort(ok, 'Error defining Tsurface variable attributes in restart file. (SUBROUTINE create_landuse_cable_restart)')
    END IF ! SLI soil model

    ! Write global attributes for file:
    CALL DATE_AND_TIME(todaydate, nowtime)
    todaydate = todaydate(1:4)//'/'//todaydate(5:6)//'/'//todaydate(7:8)
    nowtime = nowtime(1:2)//':'//nowtime(3:4)//':'//nowtime(5:6)
    ok = NF90_PUT_ATT(ncid_restart, NF90_GLOBAL, "Production",                 &
         TRIM(todaydate)//' at '//TRIM(nowtime))
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing global detail to '   &
         //TRIM(frst_out)// ' (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, NF90_GLOBAL, "Source",                     &
         'CABLE LSM restart file')
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing global detail to '   &
         //TRIM(frst_out)// ' (SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_ATT(ncid_restart, NF90_GLOBAL, "CABLE_input_file",           &
         TRIM(filename%met))
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing global detail to '   &
         //TRIM(frst_out)// ' (SUBROUTINE create_landuse_cable_restart)')

    ! End netcdf define mode:
    ok = NF90_ENDDEF(ncid_restart)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error creating restart file '      &
         //TRIM(frst_out)// '(SUBROUTINE create_landuse_cable_restart)')

    ! Write time variable:
    ok = NF90_PUT_VAR(ncid_restart, tvarID, REAL(REAL(ktau) * dels, r_2))
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error time variable to '           &
         //TRIM(frst_out)// '(SUBROUTINE create_landuse_cable_restart)')

    ! Write latitude and longitude variables:
    ok = NF90_PUT_VAR(ncid_restart, latID, latitude)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok,                                     &
         'Error writing latitude variable to '   &
         //TRIM(frst_out)// '(SUBROUTINE create_landuse_cable_restart)')

    ok = NF90_PUT_VAR(ncid_restart, lonID, longitude)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok,                                     &
         'Error writing longitude variable to '  &
         //TRIM(frst_out)// '(SUBROUTINE create_landuse_cable_restart)')

    ! Write number of active patches for each land grid cell:
    ok = NF90_PUT_VAR(ncid_restart, napID, nap)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok,                                     &
         'Error writing nap variable to '        &
         //TRIM(frst_out)// '(SUBROUTINE create_landuse_cable_restart)')

    ! Write vegetated patch fractions
    ok = NF90_PUT_VAR(ncid_restart, patchfrac_id,                            &
         lucmp%patchfrac, start = (/1/), count = (/mpx/))
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing patchfrac to '       &
         //TRIM(frst_out)// '(SUBROUTINE create_landuse_cable_restart)')

    ! Write number of veg and soil types
    ok = NF90_PUT_VAR(ncid_restart, mvtype_id,mvtype)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok,                                     &
         'Error writing mvtype parameter to '    &
         //TRIM(frst_out)// '(SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_VAR(ncid_restart, mstype_id,mstype)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok,                                     &
         'Error writing mstype parameter to '    &
         //TRIM(frst_out)// '(SUBROUTINE create_landuse_cable_restart)')


    ! Write parameters:
    ok = NF90_PUT_VAR(ncid_restart, iveg_id, lucmp%iveg)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing iveg variable to ' &
                                        //TRIM(frst_out)//'(SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_VAR(ncid_restart, isoil_id, lucmp%isoil)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing isoil variable to ' &
                                        //TRIM(frst_out)//'(SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_VAR(ncid_restart, tggID, lucmp%tgg)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing tgg variable to ' &
                                        //TRIM(frst_out)//'(SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_VAR(ncid_restart, wbID, lucmp%wb)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing wb variable to ' &
                                        //TRIM(frst_out)//'(SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_VAR(ncid_restart, wbiceID, lucmp%wbice)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing wbice variable to ' &
                                        //TRIM(frst_out)//'(SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_VAR(ncid_restart, gammzzID, lucmp%gammzz)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing gammzz variable to ' &
                                        //TRIM(frst_out)//'(SUBROUTINE create_landuse_cable_restart)')
    ! Snow dimensioned variables/parameters:
    ok = NF90_PUT_VAR(ncid_restart, ssdnID, lucmp%ssdn)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing ssdn variable to ' &
                                        //TRIM(frst_out)//'(SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_VAR(ncid_restart, smassID, lucmp%smass)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing smass variable to ' &
                                        //TRIM(frst_out)//'(SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_VAR(ncid_restart, sdepthID, lucmp%sdepth)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing sdepth variable to ' &
                                        //TRIM(frst_out)//'(SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_VAR(ncid_restart, tggsnID, lucmp%tggsn)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing tggsn variable to ' &
                                        //TRIM(frst_out)//'(SUBROUTINE create_landuse_cable_restart)')
    ! Other dims
    ok = NF90_PUT_VAR(ncid_restart, albsoilsnID, lucmp%albsoilsn)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing albsoilsn variable to ' &
                                        //TRIM(frst_out)//'(SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_VAR(ncid_restart, cplantID, lucmp%cplantx)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing cplant variable to ' &
                                        //TRIM(frst_out)//'(SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_VAR(ncid_restart, csoilID, lucmp%csoilx)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing csoil variable to ' &
                                        //TRIM(frst_out)//'(SUBROUTINE create_landuse_cable_restart)')
    ok = NF90_PUT_VAR(ncid_restart, zse_id, REAL(soil%zse, 4))
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing zse parameter to '   &
                                        //TRIM(frst_out)// '(SUBROUTINE create_restart)')
    ! Single dim:
    ok = NF90_PUT_VAR(ncid_restart, albsoil_id, lucmp%albsoil)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing albsoil parameter to '   &
                                        //TRIM(frst_out)// '(SUBROUTINE create_restart)')
    ok = NF90_PUT_VAR(ncid_restart, tssID, lucmp%tss)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing tss parameter to '   &
                                        //TRIM(frst_out)// '(SUBROUTINE create_restart)')
    ok = NF90_PUT_VAR(ncid_restart, ssdnnID, lucmp%ssdnn)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing ssdnn parameter to '   &
                                        //TRIM(frst_out)// '(SUBROUTINE create_restart)')
    ok = NF90_PUT_VAR(ncid_restart, osnowdID, lucmp%osnowd)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing osnowd parameter to '   &
                                        //TRIM(frst_out)// '(SUBROUTINE create_restart)')
    ok = NF90_PUT_VAR(ncid_restart, snageID, lucmp%snage)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing snage parameter to '   &
                                        //TRIM(frst_out)// '(SUBROUTINE create_restart)')
    ok = NF90_PUT_VAR(ncid_restart, snowdID, lucmp%snowd)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing snowd parameter to '   &
                                        //TRIM(frst_out)// '(SUBROUTINE create_restart)')
    ok = NF90_PUT_VAR(ncid_restart, rtsoilID, lucmp%rtsoil)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing rtsoil parameter to '   &
                                        //TRIM(frst_out)// '(SUBROUTINE create_restart)')
    ok = NF90_PUT_VAR(ncid_restart, isflagID, lucmp%isflag)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing isflag parameter to '   &
                                        //TRIM(frst_out)// '(SUBROUTINE create_restart)')
    ok = NF90_PUT_VAR(ncid_restart, canstoID, lucmp%cansto)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing cansto parameter to '   &
                                        //TRIM(frst_out)// '(SUBROUTINE create_restart)')
    ok = NF90_PUT_VAR(ncid_restart, sghfluxID, lucmp%sghflux)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing sghflux parameter to '   &
                                        //TRIM(frst_out)// '(SUBROUTINE create_restart)')
    ok = NF90_PUT_VAR(ncid_restart, ghfluxID, lucmp%ghflux)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing ghflux parameter to '   &
                                        //TRIM(frst_out)// '(SUBROUTINE create_restart)')
    ok = NF90_PUT_VAR(ncid_restart, runoffID, lucmp%runoff)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing runoff parameter to '   &
                                        //TRIM(frst_out)// '(SUBROUTINE create_restart)')
    ok = NF90_PUT_VAR(ncid_restart, rnof1ID, lucmp%rnof1)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing rnof1 parameter to '   &
                                        //TRIM(frst_out)// '(SUBROUTINE create_restart)')
    ok = NF90_PUT_VAR(ncid_restart, rnof2ID, lucmp%rnof2)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing rnof2 parameter to '   &
                                        //TRIM(frst_out)// '(SUBROUTINE create_restart)')
    ok = NF90_PUT_VAR(ncid_restart, gaID, lucmp%ga)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing ga parameter to '   &
                                        //TRIM(frst_out)// '(SUBROUTINE create_restart)')
    ok = NF90_PUT_VAR(ncid_restart, dgdtgID, lucmp%dgdtg)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing dgdtg parameter to '   &
                                        //TRIM(frst_out)// '(SUBROUTINE create_restart)')

    ! Close restart file
    ok = NF90_CLOSE(ncid_restart)

    write (logn, *) ' landuse on'
    WRITE (logn, '(A36)') '   Restart file complete and closed.'

  END SUBROUTINE create_landuse_cable_restart
