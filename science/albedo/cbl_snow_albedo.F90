MODULE cbl_snow_albedo_module

  IMPLICIT NONE

  PUBLIC surface_albedosn
  PRIVATE

CONTAINS

  SUBROUTINE surface_albedosn( AlbSnow, AlbSoil, mp, nrb, surface_type, soil_type, &
                            SnowDepth, SnowDensity, SoilTemp, SnowAge, coszen )

    IMPLICIT NONE

    !re-decl input args
    INTEGER, INTENT(IN) :: mp
    INTEGER, INTENT(IN) :: nrb
    REAL, INTENT(OUT)   :: AlbSnow(mp,nrb) 
    REAL, INTENT(IN)    :: AlbSoil(mp,nrb)           !NB to become IN OUT because of soil colour parameterization
    REAL, INTENT(IN)    :: coszen(mp) 
    REAL, INTENT(IN)    :: SnowDepth(mp)
    REAL, INTENT(IN)    :: SnowDensity(mp)
    REAL, INTENT(IN)    :: SoilTemp(mp)                
    REAL, INTENT(IN)    :: SnowAge(mp)
    INTEGER, INTENT(IN) :: surface_type(mp)          
    INTEGER, INTENT(IN) :: soil_type(mp) 

    !working variables  
    REAL, DIMENSION(mp) ::                                                         &
         alv,     &  ! Snow albedo for visible
         alir,    &  ! Snow albedo for near infra-red
         fage,    &  ! age factor
         fzenm,   &  ! zenith factor
         sfact,   &  ! soil factor
         snrat,   &  ! (1-) fraction of soil 'seen' when evaluating surface albedo
         tmp,     &  ! temporary value
         SoilAlbsoilF

    REAL, PARAMETER ::                                                             &
         alvo  = 0.95,  &  ! albedo for vis. on a new snow
         aliro = 0.70      ! albedo for near-infr. on a new snow

    !hard wired indexes to be substituted with arg list or USEd from module
    INTEGER, PARAMETER :: perm_ice = 9
    INTEGER, PARAMETER :: lake = 16
    !model parameter shared across subroutines -> cable_phys_constants
    REAL, PARAMETER :: snow_depth_thresh = 1.0

    INTEGER :: i    !looping variable

    !initialise to the no-snow value for albedo for all land points
    SoilAlbsoilF = Albsoil(:,1)

    ! lakes - with/without snow cover
    WHERE( surface_type == lake )                                                     
       SoilAlbsoilF = -0.022*( MIN( 275.0, MAX( 260.0, SoilTemp) ) - 260.0 ) + 0.45
    END WHERE
    WHERE(SnowDepth > snow_depth_thresh .and. surface_type == lake )
       SoilAlbsoilF = 0.85
    END WHERE

    sfact(:) = 0.68
    WHERE (SoilAlbsoilF <= 0.14)
       sfact = 0.5
    ELSEWHERE (SoilAlbsoilF > 0.14 .and. SoilAlbsoilF <= 0.20)
       sfact = 0.62
    END WHERE

    !first estimate of snow-affected surface albedos
    AlbSnow(:,2) = 2.0 * SoilAlbsoilF / (1.0 + sfact)
    AlbSnow(:,1) = sfact * AlbSnow(:,2)

    ! calc soil albedo based on colour - Ticket #27
    !H!IF (calcsoilalbedo) THEN
    !H!   CALL soilcol_albedo(ssnow, soil)
    !H!END IF

    !no snow values for working variables.
    snrat(:)=0.0
    alir(:) =0.0
    alv(:)  =0.0

    !Ticket 331 - snow age evaluation moved to cbm module
    !removed permanent ice special conditions as overwritten later
    WHERE (SnowDepth > snow_depth_thresh)

       !snrat is how little (as fraction) of the underlying soil 'seen'
       tmp = SnowDepth / MAX (SnowDensity, 200.0)
       snrat = MIN(1.0, tmp/ (tmp + 0.1))

       !snow age and zenith angle factors
       fage = 1.0 - 1.0 / (1.0 + SnowAge )
       tmp = MAX (0.17365, coszen )
       fzenm = MAX(0.0, MERGE(0.0, 1.5/(1.0+4.0*tmp) - 0.5,tmp>0.5) )

       !alv and alir: aged-snow albedo
       tmp = alvo * (1.0 - 0.2 * fage)
       alv = 0.4 * fzenm * (1.0 - tmp) + tmp
       tmp = aliro * (1.0 - 0.5 * fage )
       alir = 0.4 * fzenm * (1.0 - tmp) + tmp

    END WHERE
    
    !H!jhan:SLI currently not available
    !H!IF(cable_user%SOIL_STRUC=='sli') THEN
    !H!   WHERE (SnowDepth.GT.snow_depth_thresh)
    !H!      snrat = 1.0   ! using default parameterisation, albedo is too low,
    !H!      ! inhibiting snowpack initiation
    !H!   ENDWHERE
    !H!ENDIF

    !final values of soil-snow albedos - 1=vis, 2=nir
    AlbSnow(:,2) = MIN( aliro,                                          &
                          ( 1.0 - snrat ) * AlbSnow(:,2) + snrat * alir)

    AlbSnow(:,1) = MIN( alvo,                                           &
                          ( 1.0 - snrat ) * AlbSnow(:,1) + snrat * alv )

    !except for ice regions
    WHERE (soil_type == perm_ice) ! use dry snow albedo: 1=vis, 2=nir
       AlbSnow(:,1) = alvo - 0.05 ! al*o = albedo appropriate for new snow 
       AlbSnow(:,2) = aliro - 0.05 ! => here al*o LESS arbitrary aging 0.05
    END WHERE

    RETURN

  END SUBROUTINE surface_albedosn

END MODULE cbl_snow_albedo_module
