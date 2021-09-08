MODULE cbl_snow_albedo_module

  IMPLICIT NONE

  PUBLIC surface_albedosn
  PRIVATE

CONTAINS

SUBROUTINE surface_albedosn( AlbSnow, AlbSoil, mp, nrb, jls_radiation, surface_type, soil_type, &
                            SnowDepth, SnowODepth, SnowFlag_3L, SnowDensity, &
                            SoilTemp, SnowTemp, SnowAge, & 
                            metTk, coszen )
!H!jhan:Eliminate these USE data statements  
USE cable_common_module, ONLY : kwidth_gl
use cable_phys_constants_mod, ONLY : CTFRZ => TFRZ

implicit none

!re-decl input args
integer :: mp
integer :: nrb
LOGICAL :: jls_radiation            !runtime switch def. in cable_*main routines 
REAL :: AlbSnow(mp,nrb) 
REAL :: AlbSoil(mp,nrb) 
REAL :: MetTk(mp) 
REAL :: coszen(mp) 
REAL :: SnowDepth(mp)
REAL :: SnowODepth(mp)
REAL :: SnowDensity(mp)
REAL :: SoilTemp(mp)
REAL :: SnowTemp(mp)
REAL :: SnowAge(mp)
integer:: SnowFlag_3L(mp)
integer:: surface_type(mp) 
integer:: soil_type(mp) 


    REAL, DIMENSION(mp) ::                                                      &
         alv,     &  ! Snow albedo for visible
         alir,    &  ! Snow albedo for near infra-red
         ar1,     &  ! crystal growth  (-ve)
         ar2,     &  ! freezing of melt water
         ar3,     &  !
         dnsnow,  &  ! new snow albedo
         dtau,    &  !
         fage,    &  ! age factor
         fzenm,   &  !
         sfact,   &  !
         snr,     &  !
         snrat,   &  !
         talb,    &  ! snow albedo
         tmp         ! temporary value

    REAL, PARAMETER ::                                                          &
         alvo  = 0.95,  &  ! albedo for vis. on a new snow
         aliro = 0.70      ! albedo for near-infr. on a new snow

! local vars    
REAL :: SoilAlbsoilf(mp) 

   SoilAlbsoilF = Albsoil(:,1)

    ! lakes: hard-wired number to be removed in future
    WHERE( surface_type == 16 )                                                     &
      soilalbsoilf = -0.022*( MIN( 275., MAX( 260., SoilTemp) ) - 260. ) + 0.45

   WHERE(SnowDepth > 1. .and. surface_type == 16 ) SoilAlbsoilF = 0.85

    sfact = 0.68

   WHERE (SoilAlbsoilF <= 0.14)
       sfact = 0.5
   ELSEWHERE (SoilAlbsoilF > 0.14 .and. SoilAlbsoilF <= 0.20)
       sfact = 0.62
    END WHERE

   AlbSnow(:,2) = 2. * SoilAlbsoilF / (1. + sfact)
   AlbSnow(:,1) = sfact * AlbSnow(:,2)

    ! calc soil albedo based on colour - Ticket #27
   !H!IF (calcsoilalbedo) THEN
   !H!   CALL soilcol_albedo(ssnow, soil)
   !H!END IF

    snrat=0.
    alir =0.
    alv  =0.

   WHERE ( SnowDepth > 1. .AND. .NOT. jls_radiation )

       ! new snow (cm H2O)
      dnsnow = MIN ( 1., .1 * MAX( 0., SnowDepth - SnowODepth ) )

       ! Snow age depends on snow crystal growth, freezing of melt water,
       ! accumulation of dirt and amount of new snow.
      tmp = SnowFlag_3L * SnowTemp + ( 1 - SnowFlag_3L ) * SoilTemp
       tmp = MIN( tmp, CTFRZ )
       ar1 = 5000. * (1. / (CTFRZ-0.01) - 1. / tmp) ! crystal growth  (-ve)
       ar2 = 10. * ar1 ! freezing of melt water
      snr = SnowDepth / max (SnowDensity, 200.)

       WHERE (soil_type == 9)
          ! permanent ice: hard-wired number to be removed in future version
          ar3 = .0000001
          !  NB. dsnow =1,assumes pristine snow; ignores soot etc. ALTERNATIVELY,
          !dnsnow = max (dnsnow, .5) !increase refreshing of snow in Antarctic
          dnsnow = 1.0
          snrat = 1.

       ELSEWHERE

          ! accumulation of dirt
          ar3 = .1
          ! snow covered fraction of the grid
          snrat = MIN (1., snr / (snr + .1) )

       END WHERE

       dtau = 1.e-6 * (EXP( ar1 ) + EXP( ar2 ) + ar3 ) * kwidth_gl

      WHERE (SnowDepth <= 1.0)
         SnowAge = 0.
       ELSEWHERE
         SnowAge = max (0.,(SnowAge+dtau)*(1.-dnsnow))
       END WHERE

      fage = 1. - 1. / (1. + SnowAge ) !age factor

      tmp = MAX( .17365, coszen )
       fzenm = MAX( 0.0, MERGE( 0.0,                                           &
            ( 1. + 1./2. ) / ( 1. + 2.*2. * tmp ) - 1./2., tmp > 0.5 ) )

       tmp = alvo * (1.0 - 0.2 * fage)
       alv = .4 * fzenm * (1. - tmp) + tmp
       tmp = aliro * (1. - .5 * fage)

       ! use dry snow albedo for pernament land ice: hard-wired no to be removed
       WHERE (soil_type == 9)

          tmp = 0.95 * (1.0 - 0.2 * fage)
          alv = .4 * fzenm * (1. - tmp) + tmp
          tmp = 0.75 * (1. - .5 * fage)

       END WHERE

       alir = .4 * fzenm * (1.0 - tmp) + tmp
       talb = .5 * (alv + alir) ! snow albedo

    ENDWHERE        ! snowd > 0

   !H!! when it is called from cable_rad_driver (UM)
   !H!! no need to recalculate snage
WHERE (SnowDepth > 1 .and. jls_radiation )

   snr = SnowDepth / MAX (SnowDensity, 200.)

       WHERE (soil_type == 9)
          ! permanent ice: hard-wired number to be removed
          snrat = 1.
       ELSEWHERE
          snrat = MIN (1., snr / (snr + .1) )
       END WHERE

   fage = 1. - 1. / (1. + SnowAge ) !age factor
   tmp = MAX (.17365, coszen )
       fzenm = MAX( 0., MERGE( 0.0,                                             &
            ( 1. + 1./2. ) / ( 1. + 2. * 2. * tmp ) - 1./2., tmp > 0.5 ) )

       tmp = alvo * (1.0 - 0.2 * fage)
       alv = .4 * fzenm * (1. - tmp) + tmp
       tmp = aliro * (1. - .5 * fage)

       ! use dry snow albedo
       WHERE (soil_type == 9)
          ! permanent ice: hard-wired number to be removed

          tmp = 0.95 * (1.0 - 0.2 * fage)
          alv = .4 * fzenm * (1. - tmp) + tmp
          tmp = 0.75 * (1. - .5 * fage)

       END WHERE

       alir = .4 * fzenm * (1.0 - tmp) + tmp
       talb = .5 * (alv + alir) ! snow albedo

ENDWHERE        ! snowd > 0
!H!jhan:SLI currently not available
    !H!IF(cable_user%SOIL_STRUC=='sli') THEN
    !H!   WHERE (SnowDepth.GT.1.0)
    !H!      snrat = 1.0   ! using default parameterisation, albedo is too low,
    !H!      ! inhibiting snowpack initiation
    !H!   ENDWHERE
    !H!ENDIF
   AlbSnow(:,2) = MIN( aliro,                                          &
                          ( 1. - snrat ) * AlbSnow(:,2) + snrat * alir)

   AlbSnow(:,1) = MIN( alvo,                                           &
                          ( 1. - snrat ) * AlbSnow(:,1) + snrat * alv )

    WHERE (soil_type == 9) ! use dry snow albedo: 1=vis, 2=nir
     AlbSnow(:,1) = alvo - 0.05 ! al*o = albedo appropriate for new snow 
     AlbSnow(:,2) = aliro - 0.05 ! => here al*o LESS arbitrary aging 0.05
    END WHERE

return

END SUBROUTINE surface_albedosn

END MODULE cbl_snow_albedo_module
