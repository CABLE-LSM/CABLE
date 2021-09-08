MODULE cable_rad_driv_mod
  
CONTAINS

SUBROUTINE cable_rad_driver(                                                  &
mp,                                                                           &
nrb,                                                                          &
Clai_thresh,                                                                  &
Ccoszen_tols,                                                                 &
CGauss_w,                                                                     &
Cpi,                                                                          &
veg_mask,                                                                     &
sunlit_mask,                                                                  &
sunlit_veg_mask,                                                              &
Cpi180,                                                                       &
jls_standalone,                                                               &
jls_radiation,                                                                &
z0surf_min,                                                                   &
SurfaceType,                                                                  &
LAI_pft_cbl,                                                                  &
HGT_pft_cbl,                                                                  &
SnowDepth,                                                                    &
SnowODepth,                                                                   &
SnowFlag_3L,                                                                  &
SnowDensity,                                                                  &
SoilTemp,                                                                     &
SnowTemp,                                                                     &
SnowAge,                                                                      &
HeightAboveSnow,                                                              &
ExtCoeff_beam,                                                                &
ExtCoeff_dif,                                                                 &
EffExtCoeff_beam,                                                             &
EffExtCoeff_dif,                                                              &
CanopyRefl_dif,                                                               &
CanopyRefl_beam,                                                              &
CanopyTransmit_dif,                                                           &
CanopyTransmit_beam,                                                          &
EffSurfRefl_dif,                                                              &
EffSurfRefl_beam,                                                             &
RadAlbedo,                                                                    &
reducedLAIdue2snow,                                                           &
coszen,                                                                       &
VegXfang,                                                                     &
VegTaul,                                                                      &
VegRefl,                                                                      &
c1,                                                                           &
rhoch,                                                                        &
metDoY,                                                                       &
SW_down,                                                                      &
RadFbeam,                                                                     &
xk,                                                                           &
AlbSnow,                                                                      &
AlbSoil,                                                                      &
metTk                                                                         &
                             )

!subrs:HaC1.3
USE cbl_masks_mod,              ONLY: fveg_mask, fsunlit_mask
USE cbl_masks_mod,              ONLY: fsunlit_veg_mask
USE cbl_hruff_mod,              ONLY: HgtAboveSnow
USE cbl_LAI_eff_mod,            ONLY: LAI_eff
USE cbl_albedo_mod,             ONLY: albedo
USE cbl_init_radiation_module,  ONLY: init_radiation

IMPLICIT NONE

!___ re-decl input args

!model dimensions
!-------------------------------------------------------------------------------
!JaC:todo:ultimatelty get this from JaC~
INTEGER :: mp                       !total number of "tiles"  
INTEGER :: nrb                      !number of radiation bands [per legacy=3, but really=2 VIS,NIR. 3rd dim was for LW]
!-------------------------------------------------------------------------------

!constants
!-------------------------------------------------------------------------------
REAL :: Ccoszen_tols                !threshold cosine of sun's zenith angle, below which considered SUNLIT
REAL :: Cgauss_w(nrb)               !Gaussian integration weights
REAL :: Clai_thresh                 !The minimum LAI below which a "cell" is considred NOT vegetated
REAL :: Cpi                         !PI - describing the ratio of circumference to diameter
REAL :: Cpi180                      !PI in radians
REAL :: z0surf_min                  !the minimum roughness of bare soil
LOGICAL :: jls_standalone           !runtime switch defined in cable_*main routines signifying this is a JULES(/UM) run 
LOGICAL :: jls_radiation            !runtime switch defined in cable_*main routines signifying this is the radiation pathway 
!-------------------------------------------------------------------------------

!masks
!-------------------------------------------------------------------------------
LOGICAL, ALLOCATABLE :: veg_mask(:)             ! this "mp" is vegetated (uses minimum LAI) 
LOGICAL, ALLOCATABLE :: sunlit_mask(:)          ! this "mp" is sunlit (uses zenith angle)
LOGICAL, ALLOCATABLE :: sunlit_veg_mask(:)      ! this "mp" is BOTH sunlit AND  vegetated  
!-------------------------------------------------------------------------------

!recieved as spatial maps from the UM. remapped to "mp"
!-------------------------------------------------------------------------------
INTEGER:: surface_type(mp)          ! Integer index of Surface type (veg%iveg)
REAL :: LAI_pft_cbl(mp)             !LAI -  "limited" and remapped
REAL :: HGT_pft_cbl(mp)             !canopy height -  "limited" and remapped
!-------------------------------------------------------------------------------

REAL :: HeightAboveSnow(mp)         !Height of Canopy above snow (rough%hruff)
                                    !compute from  z0surf_min, HGT_pft_cbl, SnowDepth, SnowDensity
REAL :: reducedLAIdue2snow(mp)      ! Reduced LAI given snow coverage

!Forcing
!-------------------------------------------------------------------------------
REAL :: MetTk(mp)                   !Air Temperture at surface - atmospheric forcing (met%tk)
REAL :: coszen(mp)                  !cosine zenith angle  (met%coszen)
INTEGER :: metDoY(mp)                  !Day of the Year - not always available (met%doy)
REAL :: SW_down(mp,nrb)             !Downward shortwave "forced" (met%fsd)
!-------------------------------------------------------------------------------

!Prognostics
!-------------------------------------------------------------------------------
REAL :: SnowDepth(mp)               !Total Snow depth - water eqivalent - packed from snow_surft (ssnow%snowd)
REAL :: SnowODepth(mp)              !Total Snow depth before any update this timestep (ssnow%Osnowd)
REAL :: SnowDensity(mp)             !Total Snow density (assumes 1 layer describes snow cover) (ssnow%ssdnn)
REAL :: SoilTemp(mp)                !Soil Temperature of top layer (soil%tgg)
REAL :: SnowTemp(mp)                !Snow Temperature of top layer (soil%tgg)
REAL :: SnowAge(mp)                 !Snow age (assumes 1 layer describes snow cover) (ssnow%snage)
INTEGER:: SnowFlag_3L(mp)           !Flag to treat snow as 3 layer  - if enough present. Updated depending on total depth (ssnow%isflag)
!-------------------------------------------------------------------------------
                                                                                           
! Albedos
!-------------------------------------------------------------------------------
REAL :: AlbSoil(mp,nrb)             !Bare Soil Albedo - parametrized (soil%albsoil)
REAL :: AlbSnow(mp,nrb)             !Ground Albedo given a snow coverage (ssnow%albsoilsn)
REAL :: RadAlbedo(mp,nrb)           !Total albedo given RadFbeam (rad%albedo)
REAL :: EffSurfRefl_dif(mp,nrb)     !Effective Surface Relectance as seen by atmosphere [Diffuse SW]  (rad%reffdf)
REAL :: EffSurfRefl_beam(mp,nrb)    !Effective Surface Relectance as seen by atmosphere [Direct Beam SW] (rad%reffbm)
!-------------------------------------------------------------------------------
REAL :: RadFbeam(mp,nrb)            !Computed Beam Fraction given total SW (rad%fbeam)

!common radiation scalings [computed in albedo() ]
!-------------------------------------------------------------------------------
REAL :: xk(mp,nrb)
REAL :: c1(mp,nrb)
REAL :: rhoch(mp,nrb)
!-------------------------------------------------------------------------------

!Variables shared primarily between radiation and albedo and possibly elsewhere
!-------------------------------------------------------------------------------
!Extinction co-efficients compued in init_radiation()
REAL :: ExtCoeff_beam(mp)           !"raw" Extinction co-efficient for Direct Beam component of SW radiation (rad%extkb)
REAL :: ExtCoeff_dif(mp)            !"raw"Extinction co-efficient for Diffuse component of SW radiation (rad%extkd)
REAL :: EffExtCoeff_beam(mp,nrb)    !Effective Extinction co-efficient for Direct Beam component of SW radiation (rad%extkbm)
REAL :: EffExtCoeff_dif(mp,nrb)     !Effective Extinction co-efficient for Diffuse component of SW radiation (rad%extkdm)

!Canopy reflectance/transmitance compued in albedo() 
REAL :: CanopyRefl_dif(mp,nrb)      !Canopy reflectance (rad%cexpkdm) 
REAL :: CanopyRefl_beam(mp,nrb)     !Canopy reflectance (rad%cexpkbm)   
REAL :: CanopyTransmit_dif(mp,nrb)  !Canopy Transmitance (rad%rhodf   
REAL :: CanopyTransmit_beam(mp,nrb) !Canopy Transmitance (rad%rhobm)    
!-------------------------------------------------------------------------------

!Vegetation parameters
!-------------------------------------------------------------------------------
REAL :: VegXfang(mp)                !leaf angle PARAMETER (veg%xfang)
REAL :: VegTaul(mp,nrb)             !PARAMETER leaf transmisivity (veg%taul)
REAL :: VegRefl(mp,nrb)             !PARAMETER leaf reflectivity (veg%refl)
!-------------------------------------------------------------------------------
INTEGER:: SurfaceType(mp) 
 
CHARACTER(LEN=*), PARAMETER :: subr_name = "cable_rad_driver"

LOGICAL :: cbl_standalone = .FALSE.    

!set Height of Canopy above snow (rough%hruff) 
CALL HgtAboveSnow( HeightAboveSnow, mp, z0surf_min, HGT_pft_cbl,              &
                   SnowDepth, SnowDensity )

!set Effective LAI considering ipotentail snow coverage
CALL LAI_eff( mp, LAI_pft_cbl, HGT_pft_cbl, HeightAboveSnow,                  &
                reducedLAIdue2snow)

!define logical masks that are used throughout
CALL fveg_mask( veg_mask, mp, Clai_thresh, reducedLAIdue2snow )
CALL fsunlit_mask( sunlit_mask, mp, Ccoszen_tols, coszen )
CALL fsunlit_veg_mask( sunlit_veg_mask, mp )

!Defines Extinction Coefficients to use in calculation of Canopy 
!Reflectance/Transmitance. 
CALL init_radiation( & !rad%extkb, rad%extkd,                                     &
                     ExtCoeff_beam, ExtCoeff_dif,                             &
                     !rad%extkbm, rad%extkdm, Rad%Fbeam,                        &
                     EffExtCoeff_beam, EffExtCoeff_dif, RadFbeam,             &
                     c1, rhoch, xk,                                           &
                     mp,nrb,                                                  &
                     Clai_thresh, Ccoszen_tols, CGauss_w, Cpi, Cpi180,        &
                     cbl_standalone, jls_standalone, jls_radiation,           &
                     subr_name,                                               &
                     veg_mask, sunlit_mask, sunlit_veg_mask,                  &
                     !veg%Xfang, veg%taul, veg%refl,                            &
                     VegXfang, VegTaul, VegRefl,                              &
                     !met%coszen, int(met%DoY), met%fsd,                        &
                     coszen, metDoY, SW_down,                                 &
                     !canopy%vlaiw                                              &
                     reducedLAIdue2snow )
 
!Finally call albedo to get what we really need to fill contract with JULES
!Defines 4-"band" albedos [VIS/NIR bands. direct beam/diffuse components] from 
!considering albedo of Ground (snow?) and Canopy Reflectance/Transmitance. 

CALL Albedo( &!ssnow%AlbSoilsn, soil%AlbSoil,                                &
             AlbSnow, AlbSoil,                                                &
             mp, nrb,                                                         &
             jls_radiation,                                                   &
             veg_mask, sunlit_mask, sunlit_veg_mask,                          &
             Ccoszen_tols, cgauss_w,                                          &
             !veg%iveg, veg%refl, veg%taul,                                 & 
             SurfaceType, VegRefl, VegTaul,                                   &
             !met%tk, met%coszen, canopy%vlaiw,                             &
             metTk, coszen, reducedLAIdue2snow,                               &
             !ssnow%snowd, ssnow%osnowd, ssnow%isflag,                      & 
             SnowDepth, SnowODepth, SnowFlag_3L,                              &
             !ssnow%ssdnn, ssnow%tgg(:,1), ssnow%snage,                     & 
             SnowDensity, SoilTemp, SnowTemp, SnowAge,                        &
             xk, c1, rhoch,                                                   &
             !rad%fbeam, rad%albedo,                                        &
             RadFbeam, RadAlbedo,                                             &
             !rad%extkd, rad%extkb,                                         & 
             ExtCoeff_dif, ExtCoeff_beam,                                     &
             !rad%extkdm, rad%extkbm,                                       & 
             EffExtCoeff_dif, EffExtCoeff_beam,                               &
             !rad%rhocdf, rad%rhocbm,                                       &
             CanopyRefl_dif,CanopyRefl_beam,                                  &
             !rad%cexpkdm, rad%cexpkbm,                                     & 
             CanopyTransmit_dif, CanopyTransmit_beam,                         &
             !rad%reffdf, rad%reffbm                                        &
             EffSurfRefl_dif, EffSurfRefl_beam )


RETURN
 
END SUBROUTINE cable_rad_driver
 
END MODULE cable_rad_driv_mod

