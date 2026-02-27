MODULE cable_soil_snow_type_mod

USE cable_other_constants_mod, ONLY: r_2

IMPLICIT NONE

PUBLIC :: soil_snow_type
PUBLIC :: soil_snow_data_type
PUBLIC :: alloc_soil_snow_type
PUBLIC :: dealloc_soil_snow_type
PUBLIC :: assoc_soil_snow_type
PUBLIC :: nullify_soil_snow_cbl

! Soil and snow variables:
TYPE soil_snow_data_type

  INTEGER, ALLOCATABLE :: isflag(:) ! 0 => no snow 1 => snow
  REAL, ALLOCATABLE :: iantrct   (:) ! pointer to Antarctic land points                             
  REAL, ALLOCATABLE :: pudsto    (:) ! puddle storage                                               
  REAL, ALLOCATABLE :: pudsmx    (:) ! puddle storage                                               
  REAL, ALLOCATABLE :: cls       (:) ! factor for latent heat                                       
  REAL, ALLOCATABLE :: dfn_dtg   (:) ! d(canopy%fns)/d(ssnow%tgg)                                   
  REAL, ALLOCATABLE :: dfh_dtg   (:) ! d(canopy%fhs)/d(ssnow%tgg)                                   
  REAL, ALLOCATABLE :: dfe_ddq   (:) ! d(canopy%fes)/d(dq)        - REV_CORR: no longer necessary   
  REAL, ALLOCATABLE :: ddq_dtg   (:) ! d(dq)/d(ssnow%tgg)         - REV_CORR: no longer necessary   
  REAL, ALLOCATABLE :: dfe_dtg   (:) ! d(canopy%fes)/d(ssnow%tgg) - REV_CORR: covers above vars     
  REAL, ALLOCATABLE :: evapsn    (:) ! snow evaporation                                             
  REAL, ALLOCATABLE :: fwtop     (:) ! water flux to the soil                                       
  REAL, ALLOCATABLE :: fwtop1    (:) ! water flux to the soil                    
  REAL, ALLOCATABLE :: fwtop2    (:) ! water flux to the soil                    
  REAL, ALLOCATABLE :: fwtop3    (:) ! water flux to the soil                    
  REAL, ALLOCATABLE :: osnowd    (:) ! snow depth from previous time step        
  REAL, ALLOCATABLE :: potev     (:) ! potential evapotranspiration              
  REAL, ALLOCATABLE :: runoff    (:) ! total runoff (mm/dels)                    
  REAL, ALLOCATABLE :: rnof1     (:) ! surface runoff (mm/dels)                  
  REAL, ALLOCATABLE :: rnof2     (:) ! deep drainage (mm/dels)                   
  REAL, ALLOCATABLE :: rtsoil    (:) ! turbulent resistance for soil             
  REAL, ALLOCATABLE :: wbtot1    (:) ! total soil water (mm)                     
  REAL, ALLOCATABLE :: wbtot2    (:) ! total soil water (mm)                     
  REAL, ALLOCATABLE :: wb_lake   (:) 
  REAL, ALLOCATABLE :: totwblake (:) !daily integrated wb_lake: used in ACCESS
  REAL, ALLOCATABLE :: sinfil    (:) 
  REAL, ALLOCATABLE :: qstss     (:) 
  REAL, ALLOCATABLE :: wetfac       (:) ! surface wetness fact. at current time step  
  REAL, ALLOCATABLE :: owetfac      (:) ! surface wetness fact. at previous time step 
  REAL, ALLOCATABLE :: t_snwlr      (:) ! top snow layer depth in 3 layer snowpack    
  REAL, ALLOCATABLE :: tggav        (:) ! mean soil temperature in K                  
  REAL, ALLOCATABLE :: otgg         (:) ! soil temperature in K                       
  REAL, ALLOCATABLE :: otss         (:) ! surface temperature (weighted soil, snow)   
  REAL, ALLOCATABLE :: tprecip      (:) 
  REAL, ALLOCATABLE :: tevap        (:) 
  REAL, ALLOCATABLE :: trnoff       (:) 
  REAL, ALLOCATABLE :: totenbal     (:) 
  REAL, ALLOCATABLE :: totenbal2    (:) 
  REAL, ALLOCATABLE :: fland        (:) ! factor for latent heat
  REAL, ALLOCATABLE :: ifland       (:) ! integer soil type
  REAL, ALLOCATABLE :: qasrf        (:) ! heat advected to the snow by precip.
  REAL, ALLOCATABLE :: qfsrf        (:) ! energy of snowpack phase changes
  REAL, ALLOCATABLE :: qssrf        (:) ! sublimation
  REAL, ALLOCATABLE :: snage        (:) ! snow age
  REAL, ALLOCATABLE :: snowd        (:) ! snow depth (liquid water)
  REAL, ALLOCATABLE :: smelt        (:) ! snow melt
  REAL, ALLOCATABLE :: ssdnn        (:) ! average snow density
  REAL, ALLOCATABLE :: tss          (:) ! surface temperature (weighted soil, snow)
  REAL, ALLOCATABLE :: tss_p        (:) ! surface temperature (weighted soil, snow)
  REAL, ALLOCATABLE :: deltss       (:) ! surface temperature (weighted soil, snow)
  REAL, ALLOCATABLE :: owb1         (:) ! surface temperature (weighted soil, snow)

  REAL, ALLOCATABLE :: sconds       (:,:) !
  REAL, ALLOCATABLE :: sdepth       (:,:) ! snow depth
  REAL, ALLOCATABLE :: smass        (:,:) ! snow mass
  REAL, ALLOCATABLE :: ssdn         (:,:) ! snow densities
  REAL, ALLOCATABLE :: tgg          (:,:) ! soil temperature in K
  REAL, ALLOCATABLE :: tggsn        (:,:) ! snow temperature in K
  REAL, ALLOCATABLE :: dtmlt        (:,:) ! water flux to the soil
  REAL, ALLOCATABLE :: albsoilsn    (:,:) ! soil + snow reflectance
  REAL, ALLOCATABLE :: tilefrac     (:,:) ! factor for latent heat
  
  REAL(r_2), ALLOCATABLE :: wbtot      (:) ! total soil water (mm) 
  
  REAL(r_2), ALLOCATABLE :: evapfbl     (:,:) !
  REAL(r_2), ALLOCATABLE :: gammzz      (:,:) ! heat capacity for each soil layer
  REAL(r_2), ALLOCATABLE :: wb          (:,:) ! volumetric soil moisture (solid+liq)
  REAL(r_2), ALLOCATABLE :: wbice       (:,:) ! soil ice 
  REAL(r_2), ALLOCATABLE :: wblf        (:,:) 
  REAL(r_2), ALLOCATABLE :: wbfice      (:,:) 

  ! variables for the revised soil moisture + GW scheme
  REAL(r_2), ALLOCATABLE :: GWwb          (:) ! water content in aquifer [mm3/mm3]
  REAL(r_2), ALLOCATABLE :: GWhk          (:) ! aquifer hydraulic conductivity  [mm/s]
  REAL(r_2), ALLOCATABLE :: GWdhkdw       (:) ! aquifer d(hk) over d(water content) [(mm/s)/(mm3/mm3)]
  REAL(r_2), ALLOCATABLE :: GWdsmpdw      (:) ! aquifer d(smp) / dw   [(mm)/(mm3/mm3)]
  REAL(r_2), ALLOCATABLE :: wtd           (:) ! water table depth   [mm]
  REAL(r_2), ALLOCATABLE :: GWsmp         (:) ! aquifer soil matric potential [mm]
  REAL(r_2), ALLOCATABLE :: GWwbeq        (:) ! equilibrium aquifer water content [mm3/mm3]
  REAL(r_2), ALLOCATABLE :: GWzq          (:) ! equilibrium aquifer smp   [mm]
  REAL(r_2), ALLOCATABLE :: qhz           (:) ! horizontal hydraulic conductivity in 1D gw model for soil layers  [mm/s] 
  REAL(r_2), ALLOCATABLE :: satfrac       (:) 
  REAL(r_2), ALLOCATABLE :: Qrecharge     (:) 
  REAL(r_2), ALLOCATABLE :: rh_srf        (:) 
  REAL(r_2), ALLOCATABLE :: rtevap_sat    (:) 
  REAL(r_2), ALLOCATABLE :: rtevap_unsat  (:) 
  REAL(r_2), ALLOCATABLE :: rt_qh_sublayer(:)
   
  REAL(r_2), ALLOCATABLE :: wbeq         (:,:) ! equilibrium water content [mm3/mm3]
  REAL(r_2), ALLOCATABLE :: zq           (:,:) ! equilibrium smp       [mm]
  REAL(r_2), ALLOCATABLE :: icefrac      (:,:) ! ice fraction  [none]  -> ice mass / total mass
  REAL(r_2), ALLOCATABLE :: fracice      (:,:) ! alternate ice fraction  [none] - parameterized
  REAL(r_2), ALLOCATABLE :: hk           (:,:) ! hydraulic conductivity for soil layers [mm/s]
  REAL(r_2), ALLOCATABLE :: smp          (:,:) ! soil matric potential for soil layers         [mm]
  REAL(r_2), ALLOCATABLE :: dhkdw        (:,:) ! d(hydraulic conductivity ) d(water) for soil layers [(mm/s)/(mm3/mm3)]
  REAL(r_2), ALLOCATABLE :: dsmpdw       (:,:) ! d(smp)/ d(water) for soil layers   [(mm)/(mm3/mm3)]
  REAL(r_2), ALLOCATABLE :: wbliq        (:,:) ! volumetric liquid water content  [mm3/mm3]
  REAL(r_2), ALLOCATABLE :: wmliq        (:,:) ! water mass [mm] liq
  REAL(r_2), ALLOCATABLE :: wmice        (:,:) ! water mass [mm] ice
  REAL(r_2), ALLOCATABLE :: wmtot        (:,:) ! water mass [mm] liq+ice ->total
  REAL(r_2), ALLOCATABLE :: qhlev        (:,:) 

  ! Additional SLI variables:
  REAL(r_2), ALLOCATABLE :: S                (:,:) ! moisture content relative to sat value    (edit vh 23/01/08)
  REAL(r_2), ALLOCATABLE :: Tsoil            (:,:) !     Tsoil (deg C)
  REAL(r_2), ALLOCATABLE :: SL               (:)   ! litter moisture content relative to sat value (edit vh 23/01/08)
  REAL(r_2), ALLOCATABLE :: TL               (:)   ! litter temperature in K     (edit vh 23/01/08)
  REAL(r_2), ALLOCATABLE :: h0               (:)   ! pond height in m            (edit vh 23/01/08)
  REAL(r_2), ALLOCATABLE :: rex              (:,:) ! root extraction from each layer (mm/dels)
  REAL(r_2), ALLOCATABLE :: wflux            (:,:) ! water flux at layer boundaries (mm s-1)
  REAL(r_2), ALLOCATABLE :: delwcol          (:)   ! change in water column (mm / dels)
  REAL(r_2), ALLOCATABLE :: zdelta           (:)   ! water table depth           (edit vh 23/06/08)
  REAL(r_2), ALLOCATABLE :: kth              (:,:) ! thermal conductivity           (edit vh 29/07/08)
  REAL(r_2), ALLOCATABLE :: Tsurface         (:)   !  tepmerature at surface (soil, pond or litter) (edit vh 22/10/08)
  REAL(r_2), ALLOCATABLE :: lE               (:)   ! soil latent heat flux
  REAL(r_2), ALLOCATABLE :: evap             (:)   ! soil evaporation (mm / dels)
  REAL(r_2), ALLOCATABLE :: ciso             (:,:) ! concentration of minor isotopologue in soil water (kg m-3 water)
  REAL(r_2), ALLOCATABLE :: cisoL            (:)   ! concentration of minor isotopologue in litter water (kg m-3 water)
  REAL(r_2), ALLOCATABLE :: rlitt            (:)   ! resistance to heat/moisture transfer through litter (m-1 s)
  REAL(r_2), ALLOCATABLE :: thetai           (:,:) ! volumetric ice content (MC)
  REAL(r_2), ALLOCATABLE :: snowliq          (:,:) ! liquid snow content (mm H2O)
  REAL(r_2), ALLOCATABLE :: nsteps           (:)   ! number of iterations at each timestep
  REAL(r_2), ALLOCATABLE :: TsurfaceFR       (:)   !  tepmerature at surface (soil, pond or litter) (edit vh 22/10/08)
  REAL(r_2), ALLOCATABLE :: Ta_daily         (:,:) ! air temp averaged over last 24h
  INTEGER,   ALLOCATABLE :: nsnow            (:)   ! number of layers in snow-pack (0-nsnow_max)
  REAL(r_2), ALLOCATABLE :: Qadv_daily       (:)   ! advective heat flux into surface , daily average (W m-2)
  REAL(r_2), ALLOCATABLE :: G0_daily         (:)   ! conductive heat flux into surface , daily average (W m-2)
  REAL(r_2), ALLOCATABLE :: Qevap_daily      (:)   ! evaporative flux at surface, daily average (m s-1)
  REAL(r_2), ALLOCATABLE :: Qprec_daily      (:)   ! liquid precip, daily average (m s-1)
  REAL(r_2), ALLOCATABLE :: Qprec_snow_daily (:)   ! solid precip, daily average (m s-1)

END TYPE soil_snow_data_type

TYPE soil_snow_type

  INTEGER, POINTER :: isflag(:) ! 0 => no snow 1 => snow
  REAL, POINTER :: iantrct   (:) ! pointer to Antarctic land points                             
  REAL, POINTER :: pudsto    (:) ! puddle storage                                               
  REAL, POINTER :: pudsmx    (:) ! puddle storage                                               
  REAL, POINTER :: cls       (:) ! factor for latent heat                                       
  REAL, POINTER :: dfn_dtg   (:) ! d(canopy%fns)/d(ssnow%tgg)                                   
  REAL, POINTER :: dfh_dtg   (:) ! d(canopy%fhs)/d(ssnow%tgg)                                   
  REAL, POINTER :: dfe_ddq   (:) ! d(canopy%fes)/d(dq)        - REV_CORR: no longer necessary   
  REAL, POINTER :: ddq_dtg   (:) ! d(dq)/d(ssnow%tgg)         - REV_CORR: no longer necessary   
  REAL, POINTER :: dfe_dtg   (:) ! d(canopy%fes)/d(ssnow%tgg) - REV_CORR: covers above vars     
  REAL, POINTER :: evapsn    (:) ! snow evaporation                                             
  REAL, POINTER :: fwtop     (:) ! water flux to the soil                                       
  REAL, POINTER :: fwtop1    (:) ! water flux to the soil                    
  REAL, POINTER :: fwtop2    (:) ! water flux to the soil                    
  REAL, POINTER :: fwtop3    (:) ! water flux to the soil                    
  REAL, POINTER :: osnowd    (:) ! snow depth from previous time step        
  REAL, POINTER :: potev     (:) ! potential evapotranspiration              
  REAL, POINTER :: runoff    (:) ! total runoff (mm/dels)                    
  REAL, POINTER :: rnof1     (:) ! surface runoff (mm/dels)                  
  REAL, POINTER :: rnof2     (:) ! deep drainage (mm/dels)                   
  REAL, POINTER :: rtsoil    (:) ! turbulent resistance for soil             
  REAL, POINTER :: wbtot1    (:) ! total soil water (mm)                     
  REAL, POINTER :: wbtot2    (:) ! total soil water (mm)                     
  REAL, POINTER :: wb_lake   (:) 
  REAL, POINTER :: totwblake (:) !daily integrated wb_lake: used in ACCESS
  REAL, POINTER :: sinfil    (:) 
  REAL, POINTER :: qstss     (:) 
  REAL, POINTER :: wetfac       (:) ! surface wetness fact. at current time step  
  REAL, POINTER :: owetfac      (:) ! surface wetness fact. at previous time step 
  REAL, POINTER :: t_snwlr      (:) ! top snow layer depth in 3 layer snowpack    
  REAL, POINTER :: tggav        (:) ! mean soil temperature in K                  
  REAL, POINTER :: otgg         (:) ! soil temperature in K                       
  REAL, POINTER :: otss         (:) ! surface temperature (weighted soil, snow)   
  REAL, POINTER :: tprecip      (:) 
  REAL, POINTER :: tevap        (:) 
  REAL, POINTER :: trnoff       (:) 
  REAL, POINTER :: totenbal     (:) 
  REAL, POINTER :: totenbal2    (:) 
  REAL, POINTER :: fland        (:) ! factor for latent heat
  REAL, POINTER :: ifland       (:) ! integer soil type
  REAL, POINTER :: qasrf        (:) ! heat advected to the snow by precip.
  REAL, POINTER :: qfsrf        (:) ! energy of snowpack phase changes
  REAL, POINTER :: qssrf        (:) ! sublimation
  REAL, POINTER :: snage        (:) ! snow age
  REAL, POINTER :: snowd        (:) ! snow depth (liquid water)
  REAL, POINTER :: smelt        (:) ! snow melt
  REAL, POINTER :: ssdnn        (:) ! average snow density
  REAL, POINTER :: tss          (:) ! surface temperature (weighted soil, snow)
  REAL, POINTER :: tss_p        (:) ! surface temperature (weighted soil, snow)
  REAL, POINTER :: deltss       (:) ! surface temperature (weighted soil, snow)
  REAL, POINTER :: owb1         (:) ! surface temperature (weighted soil, snow)

  REAL, POINTER :: sconds       (:,:) !
  REAL, POINTER :: sdepth       (:,:) ! snow depth
  REAL, POINTER :: smass        (:,:) ! snow mass
  REAL, POINTER :: ssdn         (:,:) ! snow densities
  REAL, POINTER :: tgg          (:,:) ! soil temperature in K
  REAL, POINTER :: tggsn        (:,:) ! snow temperature in K
  REAL, POINTER :: dtmlt        (:,:) ! water flux to the soil
  REAL, POINTER :: albsoilsn    (:,:) ! soil + snow reflectance
  REAL, POINTER :: tilefrac     (:,:) ! factor for latent heat
  
  REAL(r_2), POINTER :: wbtot      (:) ! total soil water (mm) 
  
  REAL(r_2), POINTER :: evapfbl     (:,:) !
  REAL(r_2), POINTER :: gammzz      (:,:) ! heat capacity for each soil layer
  REAL(r_2), POINTER :: wb          (:,:) ! volumetric soil moisture (solid+liq)
  REAL(r_2), POINTER :: wbice       (:,:) ! soil ice 
  REAL(r_2), POINTER :: wblf        (:,:) 
  REAL(r_2), POINTER :: wbfice      (:,:) 

  ! variables for the revised soil moisture + GW scheme
  REAL(r_2), POINTER :: GWwb          (:) ! water content in aquifer [mm3/mm3]
  REAL(r_2), POINTER :: GWhk          (:) ! aquifer hydraulic conductivity  [mm/s]
  REAL(r_2), POINTER :: GWdhkdw       (:) ! aquifer d(hk) over d(water content) [(mm/s)/(mm3/mm3)]
  REAL(r_2), POINTER :: GWdsmpdw      (:) ! aquifer d(smp) / dw   [(mm)/(mm3/mm3)]
  REAL(r_2), POINTER :: wtd           (:) ! water table depth   [mm]
  REAL(r_2), POINTER :: GWsmp         (:) ! aquifer soil matric potential [mm]
  REAL(r_2), POINTER :: GWwbeq        (:) ! equilibrium aquifer water content [mm3/mm3]
  REAL(r_2), POINTER :: GWzq          (:) ! equilibrium aquifer smp   [mm]
  REAL(r_2), POINTER :: qhz           (:) ! horizontal hydraulic conductivity in 1D gw model for soil layers  [mm/s] 
  REAL(r_2), POINTER :: satfrac       (:) 
  REAL(r_2), POINTER :: Qrecharge     (:) 
  REAL(r_2), POINTER :: rh_srf        (:) 
  REAL(r_2), POINTER :: rtevap_sat    (:) 
  REAL(r_2), POINTER :: rtevap_unsat  (:) 
  REAL(r_2), POINTER :: rt_qh_sublayer(:)
   
  REAL(r_2), POINTER :: wbeq         (:,:) ! equilibrium water content [mm3/mm3]
  REAL(r_2), POINTER :: zq           (:,:) ! equilibrium smp       [mm]
  REAL(r_2), POINTER :: icefrac      (:,:) ! ice fraction  [none]  -> ice mass / total mass
  REAL(r_2), POINTER :: fracice      (:,:) ! alternate ice fraction  [none] - parameterized
  REAL(r_2), POINTER :: hk           (:,:) ! hydraulic conductivity for soil layers [mm/s]
  REAL(r_2), POINTER :: smp          (:,:) ! soil matric potential for soil layers         [mm]
  REAL(r_2), POINTER :: dhkdw        (:,:) ! d(hydraulic conductivity ) d(water) for soil layers [(mm/s)/(mm3/mm3)]
  REAL(r_2), POINTER :: dsmpdw       (:,:) ! d(smp)/ d(water) for soil layers   [(mm)/(mm3/mm3)]
  REAL(r_2), POINTER :: wbliq        (:,:) ! volumetric liquid water content  [mm3/mm3]
  REAL(r_2), POINTER :: wmliq        (:,:) ! water mass [mm] liq
  REAL(r_2), POINTER :: wmice        (:,:) ! water mass [mm] ice
  REAL(r_2), POINTER :: wmtot        (:,:) ! water mass [mm] liq+ice ->total
  REAL(r_2), POINTER :: qhlev        (:,:) 

  ! Additional SLI variables:
  REAL(r_2), POINTER :: S                (:,:) ! moisture content relative to sat value    (edit vh 23/01/08)
  REAL(r_2), POINTER :: Tsoil            (:,:) !     Tsoil (deg C)
  REAL(r_2), POINTER :: SL               (:)   ! litter moisture content relative to sat value (edit vh 23/01/08)
  REAL(r_2), POINTER :: TL               (:)   ! litter temperature in K     (edit vh 23/01/08)
  REAL(r_2), POINTER :: h0               (:)   ! pond height in m            (edit vh 23/01/08)
  REAL(r_2), POINTER :: rex              (:,:) ! root extraction from each layer (mm/dels)
  REAL(r_2), POINTER :: wflux            (:,:) ! water flux at layer boundaries (mm s-1)
  REAL(r_2), POINTER :: delwcol          (:)   ! change in water column (mm / dels)
  REAL(r_2), POINTER :: zdelta           (:)   ! water table depth           (edit vh 23/06/08)
  REAL(r_2), POINTER :: kth              (:,:) ! thermal conductivity           (edit vh 29/07/08)
  REAL(r_2), POINTER :: Tsurface         (:)   !  tepmerature at surface (soil, pond or litter) (edit vh 22/10/08)
  REAL(r_2), POINTER :: lE               (:)   ! soil latent heat flux
  REAL(r_2), POINTER :: evap             (:)   ! soil evaporation (mm / dels)
  REAL(r_2), POINTER :: ciso             (:,:) ! concentration of minor isotopologue in soil water (kg m-3 water)
  REAL(r_2), POINTER :: cisoL            (:)   ! concentration of minor isotopologue in litter water (kg m-3 water)
  REAL(r_2), POINTER :: rlitt            (:)   ! resistance to heat/moisture transfer through litter (m-1 s)
  REAL(r_2), POINTER :: thetai           (:,:) ! volumetric ice content (MC)
  REAL(r_2), POINTER :: snowliq          (:,:) ! liquid snow content (mm H2O)
  REAL(r_2), POINTER :: nsteps           (:)   ! number of iterations at each timestep
  REAL(r_2), POINTER :: TsurfaceFR       (:)   !  tepmerature at surface (soil, pond or litter) (edit vh 22/10/08)
  REAL(r_2), POINTER :: Ta_daily         (:,:) ! air temp averaged over last 24h
  INTEGER,   POINTER :: nsnow            (:)   ! number of layers in snow-pack (0-nsnow_max)
  REAL(r_2), POINTER :: Qadv_daily       (:)   ! advective heat flux into surface , daily average (W m-2)
  REAL(r_2), POINTER :: G0_daily         (:)   ! conductive heat flux into surface , daily average (W m-2)
  REAL(r_2), POINTER :: Qevap_daily      (:)   ! evaporative flux at surface, daily average (m s-1)
  REAL(r_2), POINTER :: Qprec_daily      (:)   ! liquid precip, daily average (m s-1)
  REAL(r_2), POINTER :: Qprec_snow_daily (:)   ! solid precip, daily average (m s-1)

END TYPE soil_snow_type

CONTAINS

SUBROUTINE alloc_soil_snow_type(soil_snow, mp)

USE grid_constants_mod_cbl,   ONLY: nsnl             ! # snow layers 
USE grid_constants_mod_cbl,   ONLY: nsl              ! # soil layers                
USE grid_constants_mod_cbl,   ONLY: nrb              ! # radiation bands *2SW, *1LW(legacy) 
USE grid_constants_mod_cbl,   ONLY: ntype_max        ! # PFTs 

IMPLICIT NONE

TYPE(soil_snow_data_type), INTENT(INOUT) :: soil_snow
INTEGER, INTENT(IN) :: mp

ALLOCATE( soil_snow% isflag       (mp) )
ALLOCATE( soil_snow% iantrct      (mp) )
ALLOCATE( soil_snow% pudsto       (mp) )
ALLOCATE( soil_snow% pudsmx       (mp) )
ALLOCATE( soil_snow% cls          (mp) )
ALLOCATE( soil_snow% dfn_dtg      (mp) )
ALLOCATE( soil_snow% dfh_dtg      (mp) )
ALLOCATE( soil_snow% dfe_ddq      (mp) )
ALLOCATE( soil_snow% ddq_dtg      (mp) )
ALLOCATE( soil_snow% dfe_dtg      (mp) )
ALLOCATE( soil_snow% evapsn       (mp) )
ALLOCATE( soil_snow% fwtop        (mp) )
ALLOCATE( soil_snow% fwtop1       (mp) )
ALLOCATE( soil_snow% fwtop2       (mp) )
ALLOCATE( soil_snow% fwtop3       (mp) )
ALLOCATE( soil_snow% osnowd       (mp) )
ALLOCATE( soil_snow% potev        (mp) )
ALLOCATE( soil_snow% runoff       (mp) )
ALLOCATE( soil_snow% rnof1        (mp) )
ALLOCATE( soil_snow% rnof2        (mp) )
ALLOCATE( soil_snow% rtsoil       (mp) )
ALLOCATE( soil_snow% wbtot1       (mp) )
ALLOCATE( soil_snow% wbtot2       (mp) )
ALLOCATE( soil_snow% wb_lake      (mp) )
ALLOCATE( soil_snow% totwblake    (mp) )
ALLOCATE( soil_snow% sinfil       (mp) )
ALLOCATE( soil_snow% qstss        (mp) )
ALLOCATE( soil_snow% wetfac       (mp) )
ALLOCATE( soil_snow% owetfac      (mp) )
ALLOCATE( soil_snow% t_snwlr      (mp) )
ALLOCATE( soil_snow% tggav        (mp) )
ALLOCATE( soil_snow% otgg         (mp) )
ALLOCATE( soil_snow% otss         (mp) )
ALLOCATE( soil_snow% tprecip      (mp) )
ALLOCATE( soil_snow% tevap        (mp) )
ALLOCATE( soil_snow% trnoff       (mp) )
ALLOCATE( soil_snow% totenbal     (mp) )
ALLOCATE( soil_snow% totenbal2    (mp) )
ALLOCATE( soil_snow% fland        (mp) )
ALLOCATE( soil_snow% ifland       (mp) )
ALLOCATE( soil_snow% qasrf        (mp) )
ALLOCATE( soil_snow% qfsrf        (mp) )
ALLOCATE( soil_snow% qssrf        (mp) )
ALLOCATE( soil_snow% snage        (mp) )
ALLOCATE( soil_snow% snowd        (mp) )
ALLOCATE( soil_snow% smelt        (mp) )
ALLOCATE( soil_snow% ssdnn        (mp) )
ALLOCATE( soil_snow% tss          (mp) )
ALLOCATE( soil_snow% tss_p        (mp) )
ALLOCATE( soil_snow% deltss       (mp) )
ALLOCATE( soil_snow% owb1         (mp) )
ALLOCATE( soil_snow% sconds       (mp,nsnl) )
ALLOCATE( soil_snow% sdepth       (mp,nsnl) )
ALLOCATE( soil_snow% smass        (mp,nsnl) )
ALLOCATE( soil_snow% ssdn         (mp,nsnl) )
ALLOCATE( soil_snow% tgg          (mp,nsl) )
ALLOCATE( soil_snow% tggsn        (mp,nsnl) )
ALLOCATE( soil_snow% dtmlt        (mp,nsnl) )
ALLOCATE( soil_snow% albsoilsn    (mp,nrb) )
ALLOCATE( soil_snow% evapfbl      (mp,nsl) )
ALLOCATE( soil_snow% tilefrac     (mp,ntype_max) )
ALLOCATE( soil_snow% wbtot        (mp) )
ALLOCATE( soil_snow% gammzz       (mp,nsl) )
ALLOCATE( soil_snow% wb           (mp,nsl) )
ALLOCATE( soil_snow% wbice        (mp,nsl) )
ALLOCATE( soil_snow% wblf         (mp,nsl) )
ALLOCATE( soil_snow% wbfice       (mp,nsl) )
ALLOCATE( soil_snow% GWwb          (mp) )
ALLOCATE( soil_snow% GWhk          (mp) )
ALLOCATE( soil_snow% GWdhkdw       (mp) )
ALLOCATE( soil_snow% GWdsmpdw      (mp) )
ALLOCATE( soil_snow% wtd           (mp) )
ALLOCATE( soil_snow% GWsmp         (mp) )
ALLOCATE( soil_snow% GWwbeq        (mp) )
ALLOCATE( soil_snow% GWzq          (mp) )
ALLOCATE( soil_snow% qhz           (mp) )
ALLOCATE( soil_snow% satfrac       (mp) )
ALLOCATE( soil_snow% Qrecharge     (mp) )
ALLOCATE( soil_snow% rh_srf        (mp) )
ALLOCATE( soil_snow% rtevap_sat    (mp) )
ALLOCATE( soil_snow% rtevap_unsat  (mp) )
ALLOCATE( soil_snow% rt_qh_sublayer(mp) )
ALLOCATE( soil_snow% wbeq           (mp,nsl) )
ALLOCATE( soil_snow% zq             (mp,nsl) )
ALLOCATE( soil_snow% icefrac        (mp,nsl) )
ALLOCATE( soil_snow% fracice        (mp,nsl) )
ALLOCATE( soil_snow% hk             (mp,nsl) )
ALLOCATE( soil_snow% smp            (mp,nsl) )
ALLOCATE( soil_snow% dhkdw          (mp,nsl) )
ALLOCATE( soil_snow% dsmpdw         (mp,nsl) )
ALLOCATE( soil_snow% wbliq          (mp,nsl) )
ALLOCATE( soil_snow% wmliq          (mp,nsl) )
ALLOCATE( soil_snow% wmice          (mp,nsl) )
ALLOCATE( soil_snow% wmtot          (mp,nsl) )
ALLOCATE( soil_snow% qhlev          (mp,nsl+1) )
ALLOCATE( soil_snow% S              (mp,nsl) )
ALLOCATE( soil_snow% Tsoil          (mp,nsl) )
ALLOCATE( soil_snow% SL             (mp) )
ALLOCATE( soil_snow% TL             (mp) )
ALLOCATE( soil_snow% h0             (mp) )
ALLOCATE( soil_snow% rex            (mp,nsl) )
ALLOCATE( soil_snow% wflux          (mp,0:nsl) )
ALLOCATE( soil_snow% delwcol        (mp) )
ALLOCATE( soil_snow% zdelta         (mp) )
ALLOCATE( soil_snow% kth            (mp,nsl) )
ALLOCATE( soil_snow% Tsurface       (mp) )
ALLOCATE( soil_snow% lE             (mp) )
ALLOCATE( soil_snow% evap           (mp) )
ALLOCATE( soil_snow% ciso           (mp,nsl+1) )
ALLOCATE( soil_snow% cisoL          (mp) )
ALLOCATE( soil_snow% rlitt          (mp) )
ALLOCATE( soil_snow% thetai         (mp,nsl) )
ALLOCATE( soil_snow% snowliq        (mp,nsnl) )
ALLOCATE( soil_snow% nsteps         (mp) )
ALLOCATE( soil_snow% TsurfaceFR     (mp) )
ALLOCATE( soil_snow% Ta_daily       (mp,100) )
ALLOCATE( soil_snow% nsnow          (mp) )
ALLOCATE( soil_snow% Qadv_daily     (mp) )
ALLOCATE( soil_snow% G0_daily       (mp) )
ALLOCATE( soil_snow% Qevap_daily    (mp) )
ALLOCATE( soil_snow% Qprec_daily    (mp) )
ALLOCATE( soil_snow% Qprec_snow_daily (mp) )

soil_snow % isflag          (:)    = 0.0      
soil_snow % iantrct         (:)    = 0.0      
soil_snow % pudsto          (:)    = 0.0      
soil_snow % pudsmx          (:)    = 0.0      
soil_snow % cls             (:)    = 0.0      
soil_snow % dfn_dtg         (:)    = 0.0      
soil_snow % dfh_dtg         (:)    = 0.0      
soil_snow % dfe_ddq         (:)    = 0.0      
soil_snow % ddq_dtg         (:)    = 0.0      
soil_snow % dfe_dtg         (:)    = 0.0      
soil_snow % evapsn          (:)    = 0.0      
soil_snow % fwtop           (:)    = 0.0      
soil_snow % fwtop1          (:)    = 0.0      
soil_snow % fwtop2          (:)    = 0.0      
soil_snow % fwtop3          (:)    = 0.0      
soil_snow % osnowd          (:)    = 0.0      
soil_snow % potev           (:)    = 0.0      
soil_snow % runoff          (:)    = 0.0      
soil_snow % rnof1           (:)    = 0.0      
soil_snow % rnof2           (:)    = 0.0      
soil_snow % rtsoil          (:)    = 0.0      
soil_snow % wbtot1          (:)    = 0.0      
soil_snow % wbtot2          (:)    = 0.0      
soil_snow % wb_lake         (:)    = 0.0      
soil_snow % totwblake       (:)    = 0.0      
soil_snow % sinfil          (:)    = 0.0      
soil_snow % qstss           (:)    = 0.0      
soil_snow % wetfac          (:)    = 0.0      
soil_snow % owetfac         (:)    = 0.0      
soil_snow % t_snwlr         (:)    = 0.0      
soil_snow % tggav           (:)    = 0.0      
soil_snow % otgg            (:)    = 0.0      
soil_snow % otss            (:)    = 0.0      
soil_snow % tprecip         (:)    = 0.0      
soil_snow % tevap           (:)    = 0.0      
soil_snow % trnoff          (:)    = 0.0      
soil_snow % totenbal        (:)    = 0.0      
soil_snow % totenbal2       (:)    = 0.0      
soil_snow % fland           (:)    = 0.0      
soil_snow % ifland          (:)    = 0.0      
soil_snow % qasrf           (:)    = 0.0      
soil_snow % qfsrf           (:)    = 0.0      
soil_snow % qssrf           (:)    = 0.0      
soil_snow % snage           (:)    = 0.0      
soil_snow % snowd           (:)    = 0.0      
soil_snow % smelt           (:)    = 0.0      
soil_snow % ssdnn           (:)    = 0.0      
soil_snow % tss             (:)    = 0.0      
soil_snow % tss_p           (:)    = 0.0      
soil_snow % deltss          (:)    = 0.0      
soil_snow % owb1            (:)    = 0.0      
soil_snow % sconds          (:,:)    = 0.0      
soil_snow % sdepth          (:,:)    = 0.0      
soil_snow % smass           (:,:)    = 0.0      
soil_snow % ssdn            (:,:)    = 0.0      
soil_snow % tgg             (:,:)    = 0.0      
soil_snow % tggsn           (:,:)    = 0.0      
soil_snow % dtmlt           (:,:)    = 0.0      
soil_snow % albsoilsn       (:,:)    = 0.0      
soil_snow % evapfbl         (:,:)    = 0.0      
soil_snow % tilefrac        (:,:)    = 0.0      
soil_snow % wbtot           (:)    = 0.0      
soil_snow % gammzz          (:,:)    = 0.0      
soil_snow % wb              (:,:)    = 0.0      
soil_snow % wbice           (:,:)    = 0.0      
soil_snow % wblf            (:,:)    = 0.0      
soil_snow % wbfice          (:,:)    = 0.0      
soil_snow % GWwb            (:)    = 0.0      
soil_snow % GWhk            (:)    = 0.0      
soil_snow % GWdhkdw         (:)    = 0.0      
soil_snow % GWdsmpdw        (:)    = 0.0      
soil_snow % wtd             (:)    = 0.0      
soil_snow % GWsmp           (:)    = 0.0      
soil_snow % GWwbeq          (:)    = 0.0      
soil_snow % GWzq            (:)    = 0.0      
soil_snow % qhz             (:)    = 0.0      
soil_snow % satfrac         (:)    = 0.0      
soil_snow % Qrecharge       (:)    = 0.0      
soil_snow % rh_srf          (:)    = 0.0      
soil_snow % rtevap_sat      (:)    = 0.0      
soil_snow % rtevap_unsat    (:)    = 0.0      
soil_snow % rt_qh_sublayer  (:)    = 0.0      
soil_snow % wbeq            (:,:)    = 0.0      
soil_snow % zq              (:,:)    = 0.0      
soil_snow % icefrac         (:,:)    = 0.0      
soil_snow % fracice         (:,:)    = 0.0      
soil_snow % hk              (:,:)    = 0.0      
soil_snow % smp             (:,:)    = 0.0      
soil_snow % dhkdw           (:,:)    = 0.0      
soil_snow % dsmpdw          (:,:)    = 0.0      
soil_snow % wbliq           (:,:)    = 0.0      
soil_snow % wmliq           (:,:)    = 0.0      
soil_snow % wmice           (:,:)    = 0.0      
soil_snow % wmtot           (:,:)    = 0.0      
soil_snow % qhlev           (:,:)    = 0.0      
soil_snow % S               (:,:)    = 0.0      
soil_snow % Tsoil           (:,:)    = 0.0      
soil_snow % SL              (:)    = 0.0      
soil_snow % TL              (:)    = 0.0      
soil_snow % h0              (:)    = 0.0      
soil_snow % rex             (:,:)    = 0.0      
soil_snow % wflux           (:,:)    = 0.0      
soil_snow % delwcol         (:)    = 0.0      
soil_snow % zdelta          (:)    = 0.0      
soil_snow % kth             (:,:)    = 0.0      
soil_snow % Tsurface        (:)    = 0.0      
soil_snow % lE              (:)    = 0.0      
soil_snow % evap            (:)    = 0.0      
soil_snow % ciso            (:,:)    = 0.0      
soil_snow % cisoL           (:)    = 0.0      
soil_snow % rlitt           (:)    = 0.0      
soil_snow % thetai          (:,:)    = 0.0      
soil_snow % snowliq         (:,:)    = 0.0      
soil_snow % nsteps          (:)    = 0.0      
soil_snow % TsurfaceFR      (:)    = 0.0      
soil_snow % Ta_daily        (:,:)    = 0.0      
soil_snow % nsnow           (:)    = 0.0      
soil_snow % Qadv_daily      (:)    = 0.0      
soil_snow % G0_daily        (:)    = 0.0      
soil_snow % Qevap_daily     (:)    = 0.0      
soil_snow % Qprec_daily     (:)    = 0.0      
soil_snow % Qprec_snow_daily(:)    = 0.0      

RETURN
END SUBROUTINE alloc_soil_snow_type

SUBROUTINE dealloc_soil_snow_type(soil_snow)

TYPE(soil_snow_type), INTENT(inout) :: soil_snow

DEALLOCATE ( soil_snow % isflag          )
DEALLOCATE ( soil_snow % iantrct         )
DEALLOCATE ( soil_snow % pudsto          )
DEALLOCATE ( soil_snow % pudsmx          )
DEALLOCATE ( soil_snow % cls             )
DEALLOCATE ( soil_snow % dfn_dtg         )
DEALLOCATE ( soil_snow % dfh_dtg         )
DEALLOCATE ( soil_snow % dfe_ddq         )
DEALLOCATE ( soil_snow % ddq_dtg         )
DEALLOCATE ( soil_snow % dfe_dtg         )
DEALLOCATE ( soil_snow % evapsn          )
DEALLOCATE ( soil_snow % fwtop           )
DEALLOCATE ( soil_snow % fwtop1          )
DEALLOCATE ( soil_snow % fwtop2          )
DEALLOCATE ( soil_snow % fwtop3          )
DEALLOCATE ( soil_snow % osnowd          )
DEALLOCATE ( soil_snow % potev           )
DEALLOCATE ( soil_snow % runoff          )
DEALLOCATE ( soil_snow % rnof1           )
DEALLOCATE ( soil_snow % rnof2           )
DEALLOCATE ( soil_snow % rtsoil          )
DEALLOCATE ( soil_snow % wbtot1          )
DEALLOCATE ( soil_snow % wbtot2          )
DEALLOCATE ( soil_snow % wb_lake         )
DEALLOCATE ( soil_snow % totwblake       )
DEALLOCATE ( soil_snow % sinfil          )
DEALLOCATE ( soil_snow % qstss           )
DEALLOCATE ( soil_snow % wetfac          )
DEALLOCATE ( soil_snow % owetfac         )
DEALLOCATE ( soil_snow % t_snwlr         )
DEALLOCATE ( soil_snow % tggav           )
DEALLOCATE ( soil_snow % otgg            )
DEALLOCATE ( soil_snow % otss            )
DEALLOCATE ( soil_snow % tprecip         )
DEALLOCATE ( soil_snow % tevap           )
DEALLOCATE ( soil_snow % trnoff          )
DEALLOCATE ( soil_snow % totenbal        )
DEALLOCATE ( soil_snow % totenbal2       )
DEALLOCATE ( soil_snow % fland           )
DEALLOCATE ( soil_snow % ifland          )
DEALLOCATE ( soil_snow % qasrf           )
DEALLOCATE ( soil_snow % qfsrf           )
DEALLOCATE ( soil_snow % qssrf           )
DEALLOCATE ( soil_snow % snage           )
DEALLOCATE ( soil_snow % snowd           )
DEALLOCATE ( soil_snow % smelt           )
DEALLOCATE ( soil_snow % ssdnn           )
DEALLOCATE ( soil_snow % tss             )
DEALLOCATE ( soil_snow % tss_p           )
DEALLOCATE ( soil_snow % deltss          )
DEALLOCATE ( soil_snow % owb1            )
DEALLOCATE ( soil_snow % sconds          )
DEALLOCATE ( soil_snow % sdepth          )
DEALLOCATE ( soil_snow % smass           )
DEALLOCATE ( soil_snow % ssdn            )
DEALLOCATE ( soil_snow % tgg             )
DEALLOCATE ( soil_snow % tggsn           )
DEALLOCATE ( soil_snow % dtmlt           )
DEALLOCATE ( soil_snow % albsoilsn       )
DEALLOCATE ( soil_snow % evapfbl         )
DEALLOCATE ( soil_snow % tilefrac        )
DEALLOCATE ( soil_snow % wbtot           )
DEALLOCATE ( soil_snow % gammzz          )
DEALLOCATE ( soil_snow % wb              )
DEALLOCATE ( soil_snow % wbice           )
DEALLOCATE ( soil_snow % wblf            )
DEALLOCATE ( soil_snow % wbfice          )
DEALLOCATE ( soil_snow % GWwb            )
DEALLOCATE ( soil_snow % GWhk            )
DEALLOCATE ( soil_snow % GWdhkdw         )
DEALLOCATE ( soil_snow % GWdsmpdw        )
DEALLOCATE ( soil_snow % wtd             )
DEALLOCATE ( soil_snow % GWsmp           )
DEALLOCATE ( soil_snow % GWwbeq          )
DEALLOCATE ( soil_snow % GWzq            )
DEALLOCATE ( soil_snow % qhz             )
DEALLOCATE ( soil_snow % satfrac         )
DEALLOCATE ( soil_snow % Qrecharge       )
DEALLOCATE ( soil_snow % rh_srf          )
DEALLOCATE ( soil_snow % rtevap_sat      )
DEALLOCATE ( soil_snow % rtevap_unsat    )
DEALLOCATE ( soil_snow % rt_qh_sublayer  )
DEALLOCATE ( soil_snow % wbeq            )
DEALLOCATE ( soil_snow % zq              )
DEALLOCATE ( soil_snow % icefrac         )
DEALLOCATE ( soil_snow % fracice         )
DEALLOCATE ( soil_snow % hk              )
DEALLOCATE ( soil_snow % smp             )
DEALLOCATE ( soil_snow % dhkdw           )
DEALLOCATE ( soil_snow % dsmpdw          )
DEALLOCATE ( soil_snow % wbliq           )
DEALLOCATE ( soil_snow % wmliq           )
DEALLOCATE ( soil_snow % wmice           )
DEALLOCATE ( soil_snow % wmtot           )
DEALLOCATE ( soil_snow % qhlev           )
DEALLOCATE ( soil_snow % S               )
DEALLOCATE ( soil_snow % Tsoil           )
DEALLOCATE ( soil_snow % SL              )
DEALLOCATE ( soil_snow % TL              )
DEALLOCATE ( soil_snow % h0              )
DEALLOCATE ( soil_snow % rex             )
DEALLOCATE ( soil_snow % wflux           )
DEALLOCATE ( soil_snow % delwcol         )
DEALLOCATE ( soil_snow % zdelta          )
DEALLOCATE ( soil_snow % kth             )
DEALLOCATE ( soil_snow % Tsurface        )
DEALLOCATE ( soil_snow % lE              )
DEALLOCATE ( soil_snow % evap            )
DEALLOCATE ( soil_snow % ciso            )
DEALLOCATE ( soil_snow % cisoL           )
DEALLOCATE ( soil_snow % rlitt           )
DEALLOCATE ( soil_snow % thetai          )
DEALLOCATE ( soil_snow % snowliq         )
DEALLOCATE ( soil_snow % nsteps          )
DEALLOCATE ( soil_snow % TsurfaceFR      )
DEALLOCATE ( soil_snow % Ta_daily        )
DEALLOCATE ( soil_snow % nsnow           )
DEALLOCATE ( soil_snow % Qadv_daily      )
DEALLOCATE ( soil_snow % G0_daily        )
DEALLOCATE ( soil_snow % Qevap_daily     )
DEALLOCATE ( soil_snow % Qprec_daily     )
DEALLOCATE ( soil_snow % Qprec_snow_daily)

RETURN
END SUBROUTINE dealloc_soil_snow_type

SUBROUTINE assoc_soil_snow_type(soil_snow, soil_snow_data )

! Description:
!   Associate the CABLE work pointers in the derived type structure

IMPLICIT NONE

!Arguments
TYPE(soil_snow_type),      INTENT(IN OUT)         :: soil_snow
TYPE(soil_snow_data_type), INTENT(IN OUT), TARGET :: soil_snow_data

CHARACTER(LEN=*), PARAMETER :: RoutineName=''
!End of header

CALL nullify_soil_snow_cbl(soil_snow)

soil_snow% isflag           => soil_snow_data% isflag           
soil_snow% iantrct          => soil_snow_data% iantrct          
soil_snow% pudsto           => soil_snow_data% pudsto           
soil_snow% pudsmx           => soil_snow_data% pudsmx           
soil_snow% cls              => soil_snow_data% cls              
soil_snow% dfn_dtg          => soil_snow_data% dfn_dtg          
soil_snow% dfh_dtg          => soil_snow_data% dfh_dtg          
soil_snow% dfe_ddq          => soil_snow_data% dfe_ddq          
soil_snow% ddq_dtg          => soil_snow_data% ddq_dtg          
soil_snow% dfe_dtg          => soil_snow_data% dfe_dtg          
soil_snow% evapsn           => soil_snow_data% evapsn           
soil_snow% fwtop            => soil_snow_data% fwtop            
soil_snow% fwtop1           => soil_snow_data% fwtop1           
soil_snow% fwtop2           => soil_snow_data% fwtop2           
soil_snow% fwtop3           => soil_snow_data% fwtop3           
soil_snow% osnowd           => soil_snow_data% osnowd           
soil_snow% potev            => soil_snow_data% potev            
soil_snow% runoff           => soil_snow_data% runoff           
soil_snow% rnof1            => soil_snow_data% rnof1            
soil_snow% rnof2            => soil_snow_data% rnof2            
soil_snow% rtsoil           => soil_snow_data% rtsoil           
soil_snow% wbtot1           => soil_snow_data% wbtot1           
soil_snow% wbtot2           => soil_snow_data% wbtot2           
soil_snow% wb_lake          => soil_snow_data% wb_lake          
soil_snow% totwblake        => soil_snow_data% totwblake        
soil_snow% sinfil           => soil_snow_data% sinfil           
soil_snow% qstss            => soil_snow_data% qstss            
soil_snow% wetfac           => soil_snow_data% wetfac           
soil_snow% owetfac          => soil_snow_data% owetfac          
soil_snow% t_snwlr          => soil_snow_data% t_snwlr          
soil_snow% tggav            => soil_snow_data% tggav            
soil_snow% otgg             => soil_snow_data% otgg             
soil_snow% otss             => soil_snow_data% otss             
soil_snow% tprecip          => soil_snow_data% tprecip          
soil_snow% tevap            => soil_snow_data% tevap            
soil_snow% trnoff           => soil_snow_data% trnoff           
soil_snow% totenbal         => soil_snow_data% totenbal         
soil_snow% totenbal2        => soil_snow_data% totenbal2        
soil_snow% fland            => soil_snow_data% fland            
soil_snow% ifland           => soil_snow_data% ifland           
soil_snow% qasrf            => soil_snow_data% qasrf            
soil_snow% qfsrf            => soil_snow_data% qfsrf            
soil_snow% qssrf            => soil_snow_data% qssrf            
soil_snow% snage            => soil_snow_data% snage            
soil_snow% snowd            => soil_snow_data% snowd            
soil_snow% smelt            => soil_snow_data% smelt            
soil_snow% ssdnn            => soil_snow_data% ssdnn            
soil_snow% tss              => soil_snow_data% tss              
soil_snow% tss_p            => soil_snow_data% tss_p            
soil_snow% deltss           => soil_snow_data% deltss           
soil_snow% owb1             => soil_snow_data% owb1             
soil_snow% sconds           => soil_snow_data% sconds           
soil_snow% sdepth           => soil_snow_data% sdepth           
soil_snow% smass            => soil_snow_data% smass            
soil_snow% ssdn             => soil_snow_data% ssdn             
soil_snow% tgg              => soil_snow_data% tgg              
soil_snow% tggsn            => soil_snow_data% tggsn            
soil_snow% dtmlt            => soil_snow_data% dtmlt            
soil_snow% albsoilsn        => soil_snow_data% albsoilsn        
soil_snow% evapfbl          => soil_snow_data% evapfbl          
soil_snow% tilefrac         => soil_snow_data% tilefrac         
soil_snow% wbtot            => soil_snow_data% wbtot            
soil_snow% gammzz           => soil_snow_data% gammzz           
soil_snow% wb               => soil_snow_data% wb               
soil_snow% wbice            => soil_snow_data% wbice            
soil_snow% wblf             => soil_snow_data% wblf             
soil_snow% wbfice           => soil_snow_data% wbfice           
soil_snow% GWwb             => soil_snow_data% GWwb             
soil_snow% GWhk             => soil_snow_data% GWhk             
soil_snow% GWdhkdw          => soil_snow_data% GWdhkdw          
soil_snow% GWdsmpdw         => soil_snow_data% GWdsmpdw         
soil_snow% wtd              => soil_snow_data% wtd              
soil_snow% GWsmp            => soil_snow_data% GWsmp            
soil_snow% GWwbeq           => soil_snow_data% GWwbeq           
soil_snow% GWzq             => soil_snow_data% GWzq             
soil_snow% qhz              => soil_snow_data% qhz              
soil_snow% satfrac          => soil_snow_data% satfrac          
soil_snow% Qrecharge        => soil_snow_data% Qrecharge        
soil_snow% rh_srf           => soil_snow_data% rh_srf           
soil_snow% rtevap_sat       => soil_snow_data% rtevap_sat       
soil_snow% rtevap_unsat     => soil_snow_data% rtevap_unsat     
soil_snow% rt_qh_sublayer   => soil_snow_data% rt_qh_sublayer   
soil_snow% wbeq             => soil_snow_data% wbeq             
soil_snow% zq               => soil_snow_data% zq               
soil_snow% icefrac          => soil_snow_data% icefrac          
soil_snow% fracice          => soil_snow_data% fracice          
soil_snow% hk               => soil_snow_data% hk               
soil_snow% smp              => soil_snow_data% smp              
soil_snow% dhkdw            => soil_snow_data% dhkdw            
soil_snow% dsmpdw           => soil_snow_data% dsmpdw           
soil_snow% wbliq            => soil_snow_data% wbliq            
soil_snow% wmliq            => soil_snow_data% wmliq            
soil_snow% wmice            => soil_snow_data% wmice            
soil_snow% wmtot            => soil_snow_data% wmtot            
soil_snow% qhlev            => soil_snow_data% qhlev            
soil_snow% S                => soil_snow_data% S                
soil_snow% Tsoil            => soil_snow_data% Tsoil            
soil_snow% SL               => soil_snow_data% SL               
soil_snow% TL               => soil_snow_data% TL               
soil_snow% h0               => soil_snow_data% h0               
soil_snow% rex              => soil_snow_data% rex              
soil_snow% wflux            => soil_snow_data% wflux            
soil_snow% delwcol          => soil_snow_data% delwcol          
soil_snow% zdelta           => soil_snow_data% zdelta           
soil_snow% kth              => soil_snow_data% kth              
soil_snow% Tsurface         => soil_snow_data% Tsurface         
soil_snow% lE               => soil_snow_data% lE               
soil_snow% evap             => soil_snow_data% evap             
soil_snow% ciso             => soil_snow_data% ciso             
soil_snow% cisoL            => soil_snow_data% cisoL            
soil_snow% rlitt            => soil_snow_data% rlitt            
soil_snow% thetai           => soil_snow_data% thetai           
soil_snow% snowliq          => soil_snow_data% snowliq          
soil_snow% nsteps           => soil_snow_data% nsteps           
soil_snow% TsurfaceFR       => soil_snow_data% TsurfaceFR       
soil_snow% Ta_daily         => soil_snow_data% Ta_daily         
soil_snow% nsnow            => soil_snow_data% nsnow            
soil_snow% Qadv_daily       => soil_snow_data% Qadv_daily       
soil_snow% G0_daily         => soil_snow_data% G0_daily         
soil_snow% Qevap_daily      => soil_snow_data% Qevap_daily      
soil_snow% Qprec_daily      => soil_snow_data% Qprec_daily      
soil_snow% Qprec_snow_daily => soil_snow_data% Qprec_snow_daily 

RETURN
END SUBROUTINE assoc_soil_snow_type

SUBROUTINE nullify_soil_snow_cbl( soil_snow )

! Description:
!   Nullify the CABLE work pointers in the derived type structure

IMPLICIT NONE

!Arguments
TYPE(soil_snow_type), INTENT(IN OUT) :: soil_snow 

CHARACTER(LEN=*), PARAMETER :: RoutineName='NULLIFY_ASSOC_CBL_TYPES'
!End of header

NULLIFY( soil_snow % isflag           )
NULLIFY( soil_snow % iantrct          )
NULLIFY( soil_snow % pudsto           )
NULLIFY( soil_snow % pudsmx           )
NULLIFY( soil_snow % cls              )
NULLIFY( soil_snow % dfn_dtg          )
NULLIFY( soil_snow % dfh_dtg          )
NULLIFY( soil_snow % dfe_ddq          )
NULLIFY( soil_snow % ddq_dtg          )
NULLIFY( soil_snow % dfe_dtg          )
NULLIFY( soil_snow % evapsn           )
NULLIFY( soil_snow % fwtop            )
NULLIFY( soil_snow % fwtop1           )
NULLIFY( soil_snow % fwtop2           )
NULLIFY( soil_snow % fwtop3           )
NULLIFY( soil_snow % osnowd           )
NULLIFY( soil_snow % potev            )
NULLIFY( soil_snow % runoff           )
NULLIFY( soil_snow % rnof1            )
NULLIFY( soil_snow % rnof2            )
NULLIFY( soil_snow % rtsoil           )
NULLIFY( soil_snow % wbtot1           )
NULLIFY( soil_snow % wbtot2           )
NULLIFY( soil_snow % wb_lake          )
NULLIFY( soil_snow % totwblake        )
NULLIFY( soil_snow % sinfil           )
NULLIFY( soil_snow % qstss            )
NULLIFY( soil_snow % wetfac           )
NULLIFY( soil_snow % owetfac          )
NULLIFY( soil_snow % t_snwlr          )
NULLIFY( soil_snow % tggav            )
NULLIFY( soil_snow % otgg             )
NULLIFY( soil_snow % otss             )
NULLIFY( soil_snow % tprecip          )
NULLIFY( soil_snow % tevap            )
NULLIFY( soil_snow % trnoff           )
NULLIFY( soil_snow % totenbal         )
NULLIFY( soil_snow % totenbal2        )
NULLIFY( soil_snow % fland            )
NULLIFY( soil_snow % ifland           )
NULLIFY( soil_snow % qasrf            )
NULLIFY( soil_snow % qfsrf            )
NULLIFY( soil_snow % qssrf            )
NULLIFY( soil_snow % snage            )
NULLIFY( soil_snow % snowd            )
NULLIFY( soil_snow % smelt            )
NULLIFY( soil_snow % ssdnn            )
NULLIFY( soil_snow % tss              )
NULLIFY( soil_snow % tss_p            )
NULLIFY( soil_snow % deltss           )
NULLIFY( soil_snow % owb1             )
NULLIFY( soil_snow % sconds           )
NULLIFY( soil_snow % sdepth           )
NULLIFY( soil_snow % smass            )
NULLIFY( soil_snow % ssdn             )
NULLIFY( soil_snow % tgg              )
NULLIFY( soil_snow % tggsn            )
NULLIFY( soil_snow % dtmlt            )
NULLIFY( soil_snow % albsoilsn        )
NULLIFY( soil_snow % evapfbl          )
NULLIFY( soil_snow % tilefrac         )
NULLIFY( soil_snow % wbtot            )
NULLIFY( soil_snow % gammzz           )
NULLIFY( soil_snow % wb               )
NULLIFY( soil_snow % wbice            )
NULLIFY( soil_snow % wblf             )
NULLIFY( soil_snow % wbfice           )
NULLIFY( soil_snow % GWwb             )
NULLIFY( soil_snow % GWhk             )
NULLIFY( soil_snow % GWdhkdw          )
NULLIFY( soil_snow % GWdsmpdw         )
NULLIFY( soil_snow % wtd              )
NULLIFY( soil_snow % GWsmp            )
NULLIFY( soil_snow % GWwbeq           )
NULLIFY( soil_snow % GWzq             )
NULLIFY( soil_snow % qhz              )
NULLIFY( soil_snow % satfrac          )
NULLIFY( soil_snow % Qrecharge        )
NULLIFY( soil_snow % rh_srf           )
NULLIFY( soil_snow % rtevap_sat       )
NULLIFY( soil_snow % rtevap_unsat     )
NULLIFY( soil_snow % rt_qh_sublayer   )
NULLIFY( soil_snow % wbeq             )
NULLIFY( soil_snow % zq               )
NULLIFY( soil_snow % icefrac          )
NULLIFY( soil_snow % fracice          )
NULLIFY( soil_snow % hk               )
NULLIFY( soil_snow % smp              )
NULLIFY( soil_snow % dhkdw            )
NULLIFY( soil_snow % dsmpdw           )
NULLIFY( soil_snow % wbliq            )
NULLIFY( soil_snow % wmliq            )
NULLIFY( soil_snow % wmice            )
NULLIFY( soil_snow % wmtot            )
NULLIFY( soil_snow % qhlev            )
NULLIFY( soil_snow % S                )
NULLIFY( soil_snow % Tsoil            )
NULLIFY( soil_snow % SL               )
NULLIFY( soil_snow % TL               )
NULLIFY( soil_snow % h0               )
NULLIFY( soil_snow % rex              )
NULLIFY( soil_snow % wflux            )
NULLIFY( soil_snow % delwcol          )
NULLIFY( soil_snow % zdelta           )
NULLIFY( soil_snow % kth              )
NULLIFY( soil_snow % Tsurface         )
NULLIFY( soil_snow % lE               )
NULLIFY( soil_snow % evap             )
NULLIFY( soil_snow % ciso             )
NULLIFY( soil_snow % cisoL            )
NULLIFY( soil_snow % rlitt            )
NULLIFY( soil_snow % thetai           )
NULLIFY( soil_snow % snowliq          )
NULLIFY( soil_snow % nsteps           )
NULLIFY( soil_snow % TsurfaceFR       )
NULLIFY( soil_snow % Ta_daily         )
NULLIFY( soil_snow % nsnow            )
NULLIFY( soil_snow % Qadv_daily       )
NULLIFY( soil_snow % G0_daily         )
NULLIFY( soil_snow % Qevap_daily      )
NULLIFY( soil_snow % Qprec_daily      )
NULLIFY( soil_snow % Qprec_snow_daily )

RETURN

END SUBROUTINE nullify_soil_snow_cbl

END MODULE cable_soil_snow_type_mod
