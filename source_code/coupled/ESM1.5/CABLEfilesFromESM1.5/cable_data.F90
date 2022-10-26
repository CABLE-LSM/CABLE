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
! Purpose: Defines constants for CABLE
!
! Contact: Jhan.Srbinovsky@csiro.au
!
! History: Combines cable_*_constants from earlier versions
!          Will include define_types in future version.
!
!
! ==============================================================================

module cable_data_module 
   implicit none
   
   public 

   ! definition of major types of constants
    
   TYPE physical_constants
      real ::                                                                  & 
      capp   = 1004.64, & ! air spec. heat (J/kg/K)
      hl = 2.5014e6, & ! air spec. heat (J/kg/K)
      hlf = 0.334e6, & ! latent heat of fusion
      !hl = 2.5104e6, & ! air spec. heat (J/kg/K)
      !hlf = 0.335e6, & ! latent heat of fusion
      dheat  = 21.5E-6, & ! molecular diffusivity for heat
      !grav   = 9.80, & ! gravity acceleration (m/s2)
      grav   = 9.8086, & ! gravity acceleration (m/s2)
      rgas   = 8.3143, & ! universal gas const  (J/mol/K)
      rmair  = 0.02897, & ! molecular wt: dry air (kg/mol)
      rmh2o  = 0.018016, & ! molecular wt: water	(kg/mol)
      sboltz = 5.67e-8, & ! Stefan-Boltz. constant (W/m2/K4)
      tfrz   = 273.16, & ! Temp (K) corresp. to 0 C
      
      ! Teten coefficients
      tetena = 6.106, & ! ??? refs?
      tetenb = 17.27, &
      tetenc = 237.3, &
      
      ! Aerodynamic parameters, diffusivities, water density:
      vonk   = 0.40, & ! von Karman constant
      a33    = 1.25, & ! inertial sublayer sw/us
      csw    = 0.50, & ! canopy sw decay (Weil theory)
      ctl    = 0.40, & ! Wagga wheat (RDD 1992, Challenges)
      apol   = 0.70, & ! Polhausen coeff: single-sided plate
      prandt = 0.71, & ! Prandtl number: visc/diffh
      schmid = 0.60, & ! Schmidt number: visc/diffw
      diffwc = 1.60, & ! diffw/diffc = H2O/CO2 diffusivity
      rhow   = 1000.0, & ! liquid water density   [kg/m3]
      emleaf = 1.0, & ! leaf emissivity
      emsoil = 1.0, & ! soil emissivity
      crd = 0.3,    & ! element drag coefficient
      csd = 0.003,  & ! substrate drag coefficient
      
      !jhan:hardwire for now. note beta2 = crd/csd
      beta2 = 0.3/0.003, & ! ratio cr/cs
      ccd   = 15.0,  & ! constant in d/h equation
      ccw_c = 2.0,   & ! ccw=(zw-d)/(h-d)
      usuhm = 0.3,   & ! (max of us/uh)
      
      ! Turbulence parameters:
      zetmul = 0.4,  & ! if niter=2, final zeta=zetmul*zetar(2)
                       ! NB> niter currently=4 see cable_define_types.F90
      zeta0  = 0.0,  & ! initial value of za/L
      zetneg = -15.0, & ! negative limit on za/L when niter>=3
      zetpos = 1.0,  & ! positive limit on za/L when niter>=3
      zdlin  = 1.0,  & ! height frac of d below which TL linear
      umin   = 0.01
       
   END TYPE physical_constants




   type math_constants
      real :: pi_c = 3.1415927
      !jhan:hardwire for now. note pi180= pi_c/180
      real :: pi180 = 3.1415927/ 180.0 ! radians / degree
   end type math_constants

   type other_constants
      !where 3 = no. radiation bands (nrb in define types)
      real, DIMENSION(3) :: gauss_w=(/0.308,0.514,0.178/) ! Gaussian integ. weights
      !--- jhan: can make these trigger of #defines/namelist
      real:: RAD_THRESH = 0.001 
      real:: LAI_THRESH = 0.001 
   end type other_constants

   type photosynthetic_constants
      integer:: maxiter=20 ! max # interations for leaf temperature
      real :: a1c3 = 9.0
      real :: a1c4 = 4.0
      real :: alpha3 = 0.200
      real :: alpha4  = 0.05
      real :: cfrd3  = 0.010
      real :: cfrd4  = 0.025
      real :: conkc0 = 302.e-6  !mol mol^-1
      real :: conko0 = 256.e-3  !mol mol^-1
      real :: convx3 = 1.0E-2
      real :: convx4 = 0.8
      real :: d0c3 = 1500.0
      real :: d0c4 = 1500.0
      real :: ekc = 59430.0  !J mol^-1
      real :: eko = 36000.0  !J mol^-1
      real :: gam0 = 28.0E-6  !mol mol^-1 @ 20C = 36.9 @ 25C
      real :: gam1 = 0.0509
      real :: gam2 = 0.0010
      real :: gsw03  = 0.01
      real :: gsw04  = 0.04
      real :: rgbwc  = 1.32
      real :: rgswc  = 1.57
      real :: tmaxj  = 45.0
      real :: tmaxv  = 45.0
      real :: tminj  = -5.0
      real :: tminv  = -5.0
      real :: toptj  = 20.0
      real :: toptv  = 20.0
      real :: trefk= 298.2  !reference temperature K
   end type photosynthetic_constants



   ! instantiate major types of constants 
   type( physical_constants ), TARGET :: phys
   type( math_constants ), TARGET  :: math
   type( other_constants ), TARGET  :: other
   type( photosynthetic_constants ), TARGET :: photo
   

   ! TYPEs of local pointers to global constants defined above 
   
   TYPE driver_type
      REAL, POINTER ::                                                         &
         ! physical constants
         TFRZ, EMSOIL, EMLEAF, SBOLTZ
   END TYPE driver_type


   TYPE icbm_type
      REAL, POINTER ::                                                         &
         ! physical constants
         GRAV, CAPP 
   END TYPE icbm_type


   TYPE iair_type
      REAL, POINTER ::                                                         &
         ! physical constants
         TFRZ, RMAIR, RGAS,                                                    &
         TETENA, TETENB, TETENC,                                               &      
         CAPP, RMH2O, HL
   END TYPE iair_type



   TYPE ialbedo_type
      ! local pointers to global constants defined above
      REAL, POINTER ::                                                         &
         ! physical constants
         TFRZ,                                                                 &
      ! other constants
         LAI_THRESH, RAD_THRESH 
   END TYPE ialbedo_type



   TYPE icanopy_type

      REAL, POINTER ::                                                         &
         ! physical constants
         TFRZ, RMAIR, RGAS, DHEAT, ZETNEG,                                     &
         ZETMUL, ZETPOS, GRAV, UMIN, TETENA,                                   &
         TETENB, TETENC, RHOW, CTL, CSW,                                       &      
         EMLEAF, EMSOIL, SBOLTZ, PRANDT, CAPP,                                 &
         RMH2O, APOL, A33, VONK, ZETA0,                                        &
         ! photosynthetic constants
         RGSWC, GAM0, GAM1, GAM2,CONKO0, CONKC0,                               &
         ALPHA3, ALPHA4, D0C3, D0C4, RGBWC,                                    &
         CONVX3, CONVX4, GSW03, GSW04,                                         &
         EKC, EKO, TREFK, A1C3, A1C4, CFRD3, CFRD4,                            &
         ! math constants
         PI_C,                                                                 &
         ! other constants
         LAI_THRESH

      INTEGER, POINTER :: MAXITER
 
   END TYPE icanopy_type



   TYPE icarbon_type
      REAL, POINTER ::                                                         &
         ! physical constants
         TFRZ
   END TYPE icarbon_type



   TYPE irad_type
      REAL, POINTER ::                                                         &
         ! physical constants
         TFRZ, EMSOIL, EMLEAF, SBOLTZ,                                         &
         CAPP, & 
         ! other constants
         LAI_THRESH, RAD_THRESH,                                               & 
         ! math constants
         PI180, PI_C
      REAL, POINTER, DIMENSION(:) ::                                           &
         GAUSS_W
   END TYPE irad_type


   TYPE irough_type
      REAL, POINTER ::                                                         &
         ! physical constants
         CSD, CRD, CCD, CCW_C, USUHM, VONK,                                    &
         A33, CTL,  ZDLIN, CSW, GRAV   
   END TYPE irough_type



   TYPE issnow_type
      REAL, POINTER ::                                                         &
         ! physical constants
         CAPP, TFRZ, HL, HLF
   END TYPE issnow_type



   

 
   INTERFACE point2constants
      MODULE PROCEDURE driver_type_ptr, cbm_type_ptr, air_type_ptr,            &
                       albedo_type_ptr, canopy_type_ptr, carbon_type_ptr,      &
                       rad_type_ptr, rough_type_ptr, ssnow_type_ptr  
   END INTERFACE 

CONTAINS   
   
   ! SUBRs associating local pointers to global constants defined above 
   ! given passed TYPE which is locally declared  

SUBROUTINE driver_type_ptr(C)    
   TYPE(driver_type) :: C
   ! physical constants
   C%TFRZ  => PHYS%TFRZ
   C%EMLEAF => PHYS%EMLEAF
   C%EMSOIL => PHYS%EMSOIL
   C%SBOLTZ => PHYS%SBOLTZ
END SUBROUTINE driver_type_ptr 

! ------------------------------------------------------------------------------

SUBROUTINE cbm_type_ptr(C)    
   TYPE(icbm_type) :: C
   ! physical constants
   C%GRAV  => PHYS%GRAV
   C%CAPP  => PHYS%CAPP
END SUBROUTINE cbm_type_ptr 

! ------------------------------------------------------------------------------

SUBROUTINE air_type_ptr(C)

   TYPE(iair_type) :: C
      
   C%TFRZ  => PHYS%TFRZ
   C%RMAIR => PHYS%RMAIR
   C%RGAS  => PHYS%RGAS            
   C%TETENA => PHYS%TETENA 
   C%TETENB => PHYS%TETENB 
   C%TETENC => PHYS%TETENC 
   C%CAPP  => PHYS%CAPP
   C%RMH2O => PHYS%RMH2O
   C%HL     => PHYS%HL 

END SUBROUTINE air_type_ptr

! ------------------------------------------------------------------------------

SUBROUTINE albedo_type_ptr(C) 
   TYPE(ialbedo_type) :: C
   ! physical constants
   C%TFRZ  => PHYS%TFRZ
   ! other constants
   C%LAI_THRESH => OTHER%LAI_THRESH 
   C%RAD_THRESH => OTHER%RAD_THRESH 
END SUBROUTINE albedo_type_ptr

! ------------------------------------------------------------------------------
 
SUBROUTINE canopy_type_ptr(C)    
   
   TYPE(icanopy_type) :: C

   ! physical constants
   C%TFRZ  => PHYS%TFRZ
   C%RMAIR => PHYS%RMAIR
   C%RGAS  => PHYS%RGAS            
   C%DHEAT => PHYS%DHEAT
   C%ZETNEG => PHYS%ZETNEG
   C%ZETMUL => PHYS%ZETMUL
   C%ZETPOS => PHYS%ZETPOS
   C%GRAV  => PHYS%GRAV            
   C%UMIN  => PHYS%UMIN
   C%TETENA => PHYS%TETENA 
   C%TETENB => PHYS%TETENB 
   C%TETENC => PHYS%TETENC 
   C%RHOW  => PHYS%RHOW
   C%CTL   => PHYS%CTL
   C%CSW   => PHYS%CSW
   C%EMLEAF => PHYS%EMLEAF
   C%EMSOIL => PHYS%EMSOIL
   C%SBOLTZ => PHYS%SBOLTZ
   C%PRANDT => PHYS%PRANDT
   C%CAPP  => PHYS%CAPP
   C%RMH2O => PHYS%RMH2O
   C%APOL  => PHYS%APOL
   C%A33   => PHYS%A33
   C%VONK  => PHYS%VONK
   C%ZETA0 => PHYS%ZETA0
      
   C%MAXITER  => PHOTO%MAXITER ! only integer here

   !photosynthetic constants
   C%RGSWC => PHOTO%RGSWC               
   C%GAM0  => PHOTO%GAM0
   C%GAM2  => PHOTO%GAM2
   C%CONKC0  => PHOTO%CONKC0
   C%CONKO0  => PHOTO%CONKO0
   C%ALPHA3  => PHOTO%ALPHA3
   C%ALPHA4  => PHOTO%ALPHA4
   C%GSW03 => PHOTO%GSW03
   C%CONVX4  => PHOTO%CONVX4
   C%CONVX3  => PHOTO%CONVX3
   C%D0C3  => PHOTO%D0C3
   C%D0C4  => PHOTO%D0C4
   C%RGBWC  => PHOTO%RGBWC
   C%CFRD3  => PHOTO%CFRD3
   C%GSW04  => PHOTO%GSW04
   C%GAM1  => PHOTO%GAM1
   C%EKO  => PHOTO%EKO
   C%EKC  => PHOTO%EKC
   C%TREFK => PHOTO%TREFK
   C%A1C3 => PHOTO%A1C3
   C%A1C4  => PHOTO%A1C4
   C%CFRD4  => PHOTO%CFRD4

   ! math constants
   C%PI_C  => MATH%PI_C
   
   ! other constants
   C%LAI_THRESH  => OTHER%LAI_THRESH

END SUBROUTINE canopy_type_ptr

! ------------------------------------------------------------------------------

SUBROUTINE carbon_type_ptr(C)    
   TYPE(icarbon_type) :: C
   ! physical constants
   C%TFRZ  => PHYS%TFRZ
END SUBROUTINE carbon_type_ptr 

! ------------------------------------------------------------------------------

SUBROUTINE rad_type_ptr(C)    
   TYPE(irad_type) :: C
   
   ! other constants
   C%LAI_THRESH => OTHER%LAI_THRESH 
   C%RAD_THRESH => OTHER%RAD_THRESH 
   C%GAUSS_W  => OTHER%GAUSS_W
   
   ! math constants
   C%PI180  => MATH%PI180
   C%PI_C  => MATH%PI_C
   
   ! physical constants
   C%TFRZ  => PHYS%TFRZ
   C%EMLEAF => PHYS%EMLEAF
   C%EMSOIL => PHYS%EMSOIL
   C%SBOLTZ => PHYS%SBOLTZ
   C%CAPP  => PHYS%CAPP

END SUBROUTINE rad_type_ptr

! ------------------------------------------------------------------------------

SUBROUTINE rough_type_ptr(C)    
   TYPE(irough_type) :: C
   ! physical constants
         C%CSD   => PHYS%CSD                                                     
         C%CRD   => PHYS%CRD                                                      
         C%CCD   => PHYS%CCD                                                      
         C%CSW   => PHYS%CSW                                                      
         C%CCW_C => PHYS%CCW_C                                                      
         C%USUHM => PHYS%USUHM                                                  
         C%VONK  => PHYS%VONK                                                   
         C%A33   => PHYS%A33                                                      
         C%CTL   => PHYS%CTL                                                    
         C%ZDLIN => PHYS%ZDLIN                                                    
         C%GRAV  => PHYS%GRAV                                                   
END SUBROUTINE rough_type_ptr 

! ------------------------------------------------------------------------------

SUBROUTINE ssnow_type_ptr(C)    
   TYPE(issnow_type) :: C
   ! physical constants
   C%CAPP  => PHYS%CAPP
   C%TFRZ  => PHYS%TFRZ
   C%HL    => PHYS%HL
   C%HLF   => PHYS%HLF
   !C% => PHYS%
END SUBROUTINE ssnow_type_ptr 

END MODULE cable_data_module 

