MODULE cbl_albedo_mod

  IMPLICIT NONE

  PUBLIC albedo
  PRIVATE

CONTAINS

SUBROUTINE Albedo( AlbSnow, AlbSoil,              & 
mp, nrb,                                          &
jls_radiation ,                                   &
veg_mask, sunlit_mask, sunlit_veg_mask,           &  
Ccoszen_tols, CGAUSS_W,                           & 
surface_type, soil_type, VegRefl, VegTaul,        &
coszen, reducedLAIdue2snow,                       &
SnowDepth, SnowDensity, SoilTemp, SnowAge,        &
xk, c1, rhoch,                                    & 
RadFbeam, RadAlbedo,                              &
ExtCoeff_dif, ExtCoeff_beam,                      &
EffExtCoeff_dif, EffExtCoeff_beam,                &
CanopyRefl_dif,CanopyRefl_beam,                   &
CanopyTransmit_dif, CanopyTransmit_beam,          &
EffSurfRefl_dif, EffSurfRefl_beam                 )

!subrs called
USE cbl_snow_albedo_module, ONLY : surface_albedosn

implicit none

!model dimensions
integer :: mp                       !total number of "tiles"  
integer :: nrb                      !number of radiation bands [per legacy=3, but really=2 VIS,NIR. 3rd dim was for LW]

!This is what we are returning here
REAL :: EffSurfRefl_dif(mp,nrb)     !Effective Surface Relectance as seen by atmosphere [Diffuse SW]  (rad%reffdf)
REAL :: EffSurfRefl_beam(mp,nrb)    !Effective Surface Relectance as seen by atmosphere [Direct Beam SW] (rad%reffbm)

!constants
real :: Ccoszen_tols                !threshold cosine of sun's zenith angle, below which considered SUNLIT
real :: Cgauss_w(nrb)
LOGICAL :: jls_radiation            !runtime switch def. in cable_*main routines 
                                    !signifying this is the radiation pathway 

!masks
LOGICAL :: veg_mask(mp)             ! this "mp" is vegetated (uses minimum LAI) 
LOGICAL :: sunlit_mask(mp)          ! this "mp" is sunlit (uses zenith angle)
LOGICAL :: sunlit_veg_mask(mp)      ! this "mp" is BOTH sunlit AND  vegetated  

!Vegetation parameters
REAL :: VegTaul(mp,nrb)             !PARAMETER leaf transmisivity (veg%taul)
REAL :: VegRefl(mp,nrb)             !PARAMETER leaf reflectivity (veg%refl)
integer:: surface_type(mp)          !Integer index of Surface type (veg%iveg)
integer:: soil_type(mp)          !Integer index of Soil    type (soil%isoilm)

real :: reducedLAIdue2snow(mp)      !Reduced LAI given snow coverage

! Albedos
REAL :: AlbSoil(mp,nrb)             !Bare Soil Albedo - parametrized (soil%albsoil)
REAL :: AlbSnow(mp,nrb)             !Ground Albedo given a snow coverage (ssnow%albsoilsn)
REAL :: RadAlbedo(mp,nrb)           !Total albedo given RadFbeam (rad%albedo)

!Forcing
REAL :: coszen(mp)                  !cosine zenith angle  (met%coszen)
REAL :: SW_down(mp,nrb)             !Downward shortwave "forced" (met%fsd)

!Prognostics
REAL :: SnowDepth(mp)               !Total Snow depth - water eqivalent - packed from snow_surft (SnowDepth)
REAL :: SnowDensity(mp)             !Total Snow density (assumes 1 layer describes snow cover) (SnowDensity)
REAL :: SoilTemp(mp)                !Soil Temperature of top layer - for lake alebdo (ssnow%tgg)
REAL :: SnowAge(mp)                 !Snow age (assumes 1 layer describes snow cover) (SnowAge)

REAL :: RadFbeam(mp,nrb)            !Computed Beam Fraction given total SW (rad%fbeam)

!common radiation scalings - computed  in init_radiation()
REAL :: xk(mp,nrb)
REAL :: c1(mp,nrb)
REAL :: rhoch(mp,nrb)

!Variables shared primarily between radiation and albedo and possibly elsewhere
!Extinction co-efficients computed in init_radiation()
REAL :: ExtCoeff_beam(mp)           !"raw" Extinction co-efficient 
                                    !Direct Beam component of SW radiation (rad%extkb)
REAL :: ExtCoeff_dif(mp)            !"raw"Extinction co-efficient for 
                                    !Diffuse component of SW radiation (rad%extkd)
REAL :: EffExtCoeff_beam(mp,nrb)    !Effective Extinction co-eff 
                                    !Direct Beam component of SW radiation (rad%extkbm)
REAL :: EffExtCoeff_dif(mp,nrb)     !Effective Extinction co-eff 
                                    !Diffuse component of SW radiation (rad%extkdm)

!Canopy reflectance/transmitance compued in albedo() 
REAL :: CanopyRefl_dif(mp,nrb)      !Canopy reflectance  (rad%rhocdf   
REAL :: CanopyRefl_beam(mp,nrb)     !Canopy reflectance  (rad%rhocbm)   
REAL :: CanopyTransmit_dif(mp,nrb)  !Canopy Transmitance (rad%cexpkdm)   
REAL :: CanopyTransmit_beam(mp,nrb) !Canopy Transmitance (rad%cexpkbm)

real :: SumEffSurfRefl_beam(1)
real :: SumEffSurfRefl_dif(1)
integer :: i

    ! END header

AlbSnow(:,:) = 0.0
!CanopyTransmit_beam(:,:) = 0.0
CanopyRefl_beam(:,:) = 0.0
CanopyRefl_dif(:,:) = 0.0        
!CanopyTransmit_dif(:,:) = 0.0  ! MPI (at least inits this = 1.0 at dt=0) 

!Modify parametrised soil albedo based on snow coverage 
!call surface_albedosn( AlbSnow, AlbSoil, mp, nrb, jls_radiation, surface_type, soil_type, &
!                       SnowDepth, SnowODepth, SnowFlag_3L,                      & 
!                       SnowDensity, SoilTemp, SnowTemp, SnowAge,                     & 
!                       MetTk, Coszen )
call surface_albedosn( AlbSnow, AlbSoil, mp, nrb, surface_type, soil_type, &
                       SnowDepth, SnowDensity, SoilTemp, SnowAge, Coszen )

! Update fractional leaf transmittance and reflection
!---1 = visible, 2 = nir radiaition

! Define canopy Reflectance for diffuse/direct radiation
! Formerly rad%rhocbm, rad%rhocdf
call CanopyReflectance( CanopyRefl_beam, CanopyRefl_dif, &
                        mp, nrb, CGauss_w, sunlit_veg_mask, &
                        AlbSnow, xk, rhoch,                  &
                        ExtCoeff_beam, ExtCoeff_dif)

! Define canopy diffuse transmittance 
! Formerly rad%cexpkbm, rad%cexpkdm
call CanopyTransmitance(CanopyTransmit_beam, CanopyTransmit_dif, mp, nrb,&
                              sunlit_veg_mask, reducedLAIdue2snow, &
                              EffExtCoeff_dif, EffExtCoeff_beam)

!---1 = visible, 2 = nir radiaition
! Finally compute Effective 4-band albedo for diffuse/direct radiation- 
! In the UM this is the required variable to be passed back on the rad call
! Formerly rad%reffbm, rad%reffdf

! Even when there is no vegetation, albedo is at least snow modified soil albedo
EffSurfRefl_dif = AlbSnow
EffSurfRefl_beam = AlbSnow

call EffectiveSurfaceReflectance( EffSurfRefl_beam, EffSurfRefl_dif,           &
                                  mp, nrb, veg_mask, sunlit_veg_mask,          &
                                  CanopyRefl_beam, CanopyRefl_dif,             &
                                  CanopyTransmit_beam,CanopyTransmit_dif,      &
                                  AlbSnow )

! Compute total albedo to SW given the Effective Surface Reflectance 
! (considering Canopy/Soil/Snow contributions) 
! we dont need to do this on rad call AND may not haveappropriate RadFbeam
RadAlbedo = AlbSnow
if(.NOT. jls_radiation) &
  call FbeamRadAlbedo( RadAlbedo, mp, nrb, veg_mask, radfbeam, &
                       EffSurfRefl_dif, EffSurfRefl_beam, AlbSnow )
 
END SUBROUTINE albedo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CanopyReflectance( CanopyRefl_beam, CanopyRefl_dif, &
                         mp, nrb, CGauss_w, sunlit_veg_mask, &
                         AlbSnow, xk, rhoch,                  &
                         ExtCoeff_beam, ExtCoeff_dif)
implicit none 
!re-decl in args
integer :: mp                       !total number of "tiles"  
integer :: nrb                      !number of radiation bands [per legacy=3, but really=2 VIS,NIR. 3rd dim was for LW]
REAL :: CanopyRefl_dif(mp,nrb)      !Canopy reflectance (rad%cexpkdm) 
REAL :: CanopyRefl_beam(mp,nrb)     !Canopy reflectance (rad%cexpkbm)   
real :: Cgauss_w(nrb)
LOGICAL :: sunlit_veg_mask(mp)      ! this "mp" is BOTH sunlit AND  vegetated  
REAL :: AlbSnow(mp,nrb)             !Ground Albedo given a snow coverage (ssnow%albsoilsn)
REAL :: xk(mp,nrb)
REAL :: rhoch(mp,nrb)
REAL :: ExtCoeff_beam(mp)           !"raw" Extinction co-efficient for Direct Beam component of SW radiation (rad%extkb)
REAL :: ExtCoeff_dif(mp)            !"raw"Extinction co-efficient for Diffuse component of SW radiation (rad%extkd)

! Initialise canopy beam reflectance:
!HACHvstrunk!CanopyRefl_beam  = AlbSnow !Formerly rad%reffbm
!HACHvstrunk!CanopyRefl_dif   = AlbSnow ! Formerly rad%refdfm

call CanopyReflectance_beam( CanopyRefl_beam, mp, nrb, sunlit_veg_mask, &
                             ExtCoeff_beam, ExtCoeff_dif, rhoch )

call CanopyReflectance_dif( CanopyRefl_dif, mp, nrb, CGauss_w, &
                             ExtCoeff_dif, xk, rhoch )
End subroutine CanopyReflectance

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CanopyReflectance_beam( CanopyRefl_beam, mp, nrb, sunlit_veg_mask, &
              ExtCoeff_beam,ExtCoeff_dif, rhoch )
implicit none
integer :: mp
integer ::nrb 
real :: CanopyRefl_beam(mp,nrb)
real :: ExtCoeff_dif(mp) 
real :: ExtCoeff_beam(mp) 
LOGICAL :: sunlit_veg_mask(mp) 
REAL :: rhoch(mp, nrb)
integer :: i, b

! Canopy reflection (6.21) beam:
DO i = 1,mp
  DO b = 1, 2
    IF( sunlit_veg_mask(i) ) &
      CanopyRefl_beam(i,b) = 2. * ExtCoeff_beam(i) / &
                            ( ExtCoeff_beam(i) + ExtCoeff_dif(i) )          & 
                            * rhoch(i,b)
    END DO
END DO

End subroutine CanopyReflectance_beam

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CanopyReflectance_dif( CanopyRefl_dif, mp, nrb, CGauss_w,  &
                                  ExtCoeff_dif, xk, rhoch )

implicit none
INTEGER :: mp
integer :: nrb
real :: Cgauss_w(nrb)
REAL :: CanopyRefl_dif(mp,nrb)  
real :: ExtCoeff_dif(mp)    
REAL :: xk(mp,nrb)      ! extinct. coef.for beam rad. and black leaves
REAL :: rhoch(mp,nrb)      
!local vars
INTEGER :: ictr

! Canopy REFLection of diffuse radiation for black leaves:
DO ictr=1,2

  CanopyRefl_dif(:,ictr) = rhoch(:,ictr) *  2. *                                &
                       ( CGAUSS_W(1) * xk(:,1) / ( xk(:,1) + ExtCoeff_dif(:) )&
                       + CGAUSS_W(2) * xk(:,2) / ( xk(:,2) + ExtCoeff_dif(:) )&
                       + CGAUSS_W(3) * xk(:,3) / ( xk(:,3) + ExtCoeff_dif(:) ) )

ENDDO

End subroutine CanopyReflectance_dif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CanopyTransmitance(CanopyTransmit_beam, CanopyTransmit_dif, mp, nrb,&
                              mask, reducedLAIdue2snow, &
                              EffExtCoeff_dif, EffExtCoeff_beam)
implicit none
!re-decl in args
integer :: mp                       !total number of "tiles"  
integer :: nrb                      !number of radiation bands [per legacy=3, but really=2 VIS,NIR. 3rd dim was for LW]
REAL :: CanopyTransmit_dif(mp,nrb)      !Canopy Transmitance (rad%cexpkdm) 
REAL :: CanopyTransmit_beam(mp,nrb)     !Canopy Transmitance (rad%cexpkbm)   
LOGICAL :: mask(mp)      ! this "mp" is BOTH sunlit AND  vegetated  
LOGICAL :: dummyMask(mp)
real :: reducedLAIdue2snow(mp)
REAL :: EffExtCoeff_beam(mp,nrb)           !"raw" Extinction co-efficient for Direct Beam component of SW radiation (rad%extkb)
REAL :: EffExtCoeff_dif(mp,nrb)            !"raw"Extinction co-efficient for Diffuse component of SW radiation (rad%extkd)

! For beam, compute canopy trasmitance when sunlit (and vegetated)
call CanopyTransmitance_beam( CanopyTransmit_beam, mp, nrb, EffExtCoeff_beam,  &
                              reducedLAIdue2snow, mask )

!'=1.0' initialization remains the calculated value where "mask"=FALSE
dummyMask(:) = .true. 

! For diffuse rad, always compute canopy trasmitance
call CanopyTransmitance_dif( CanopyTransmit_dif, mp, nrb, EffExtCoeff_dif, &
                             reducedLAIdue2snow, dummyMask )

End subroutine CanopyTransmitance

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CanopyTransmitance_dif(CanopyTransmit, mp, nrb, ExtinctionCoeff, reducedLAIdue2snow, mask )
implicit none
integer :: mp 
integer :: nrb
logical :: mask(mp) 
real :: CanopyTransmit(mp,nrb) 
real :: ExtinctionCoeff(mp,nrb) 
real :: reducedLAIdue2snow(mp)
real :: dummy(mp,nrb) 
integer :: i, b
 
DO i = 1,mp
  DO b = 1, 2 
    dummy(i,b) = ExtinctionCoeff(i,b) * reducedLAIdue2snow(i)
    CanopyTransmit(i,b) = EXP( -1.* dummy(i,b) )
  enddo
enddo

End subroutine  CanopyTransmitance_dif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine CanopyTransmitance_beam(CanopyTransmit, mp, nrb, ExtinctionCoeff, reducedLAIdue2snow, mask )
implicit none
integer :: mp 
integer :: nrb
logical :: mask(mp) 
real :: CanopyTransmit(mp,nrb) 
real :: ExtinctionCoeff(mp,nrb) 
real :: reducedLAIdue2snow(mp)
real :: dummy(mp,nrb) 
integer :: i, b
 
DO i = 1,mp
  DO b = 1, 2 
    if( mask(i) ) then 
      dummy(i,b) = min( ExtinctionCoeff(i,b) * reducedLAIdue2snow(i), 20. )
      CanopyTransmit(i,b) = EXP( -1.* dummy(i,b) )
    endif
  enddo
enddo

End subroutine  CanopyTransmitance_beam


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine EffectiveSurfaceReflectance(EffSurfRefl_beam, EffSurfRefl_dif,      &
                                       mp, nrb, veg_mask, sunlit_veg_mask,     &
                                       CanopyRefl_beam, CanopyRefl_dif,        &
                                       CanopyTransmit_beam,CanopyTransmit_dif, & 
                                       AlbSnow )
implicit none
!re-decl input args 
integer :: mp                       !total number of "tiles"  
integer :: nrb                      !number of radiation bands [per legacy=3, but really=2 VIS,NIR. 3rd dim was for LW]
LOGICAL :: veg_mask(mp)             ! this "mp" is vegetated (uses minimum LAI) 
LOGICAL :: sunlit_veg_mask(mp)      ! this "mp" is vegetated (uses minimum LAI) 
REAL :: EffSurfRefl_dif(mp,nrb)     !Effective Surface Relectance as seen by atmosphere [Diffuse SW]  (rad%reffdf)
REAL :: EffSurfRefl_beam(mp,nrb)    !Effective Surface Relectance as seen by atmosphere [Direct Beam SW] (rad%reffbm)
REAL :: CanopyRefl_beam(mp,nrb)  
REAL :: CanopyRefl_dif(mp,nrb)  
REAL :: CanopyTransmit_dif(mp,nrb)      !Canopy reflectance (rad%cexpkdm) 
REAL :: CanopyTransmit_beam(mp,nrb)     !Canopy reflectance (rad%cexpkbm)   
real :: AlbSnow(mp,nrb)



call EffectiveReflectance( EffSurfRefl_dif, mp, nrb, CanopyRefl_dif, AlbSnow, &
                           CanopyTransmit_dif, veg_mask )

call EffectiveReflectance( EffSurfRefl_beam, mp, nrb, CanopyRefl_beam, AlbSnow,&
                           CanopyTransmit_beam, sunlit_veg_mask )

End subroutine EffectiveSurfaceReflectance

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine EffectiveReflectance( EffRefl, mp, nrb, CanopyRefl, AlbSnow, &
          CanopyTransmit, mask )
implicit none
integer :: mp
integer :: nrb
real :: AlbSnow(mp,nrb)
real :: CanopyRefl(mp,nrb)
real :: CanopyTransmit(mp,nrb) 
real :: EffRefl(mp,nrb) 
logical :: mask(mp) 
integer :: i,b  

DO i = 1,mp
  DO b = 1, 2!ithis is fixed as 2  because nrb=3 due to legacy  
      IF( mask(i) ) then 
      
         ! Calculate effective beam reflectance (fraction):
         EffRefl(i,b) = CanopyRefl(i,b) &
                            + ( AlbSnow(i,b) - CanopyRefl(i,b) ) &
                            * CanopyTransmit(i,b)**2

    endif
  END DO
END DO

End subroutine EffectiveReflectance

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine FbeamRadAlbedo( RadAlbedo, mp, nrb, veg_mask, radfbeam, &
                           EffSurfRefl_dif, EffSurfRefl_beam, AlbSnow )
implicit none
!re-decl input args 
integer :: mp                       !total number of "tiles"  
integer :: nrb                      !number of radiation bands [per legacy=3, but really=2 VIS,NIR. 3rd dim was for LW]
REAL :: RadAlbedo(mp,nrb)           !Total albedo given RadFbeam (rad%albedo)
REAL :: AlbSnow(mp,nrb)            !Computed Beam Fraction given total SW (rad%fbeam)
LOGICAL :: veg_mask(mp)             ! this "mp" is vegetated (uses minimum LAI) 
REAL :: RadFbeam(mp,nrb)            !Computed Beam Fraction given total SW (rad%fbeam)
REAL :: EffSurfRefl_dif(mp,nrb)     !Effective Surface Relectance as seen by atmosphere [Diffuse SW]  (rad%reffdf)
REAL :: EffSurfRefl_beam(mp,nrb)    !Effective Surface Relectance as seen by atmosphere [Direct Beam SW] (rad%reffbm)
!local vars
INTEGER :: b    !rad. band 1=visible, 2=near-infrared, 3=long-wave
INTEGER :: i    

! Initialise total albedo:
RadAlbedo = AlbSnow
DO i = 1,mp
  DO b = 1, 2 !nrb -1 -nrb shouldnt be =3 anyway
    ! Define albedo:
    IF( veg_mask(i) )                                      &
       RadAlbedo(i,b) = ( 1. - radfbeam(i,b) )*EffSurfRefl_dif(i,b) +           &
                         radfbeam(i,b) * EffSurfRefl_beam(i,b)
  END DO
END DO

End subroutine FbeamRadAlbedo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE cbl_albedo_mod
