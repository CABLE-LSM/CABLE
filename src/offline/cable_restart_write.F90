module cable_restart_write_mod
  use cable_restart_mod, only: cable_restart_write_time
  use cable_restart_mod, only: cable_restart_variable_write
  use cable_restart_mod, only: cable_restart_variable_write_darray

  use cable_common_module, only: cable_user

  use cable_def_types_mod, only: met_type
  use cable_def_types_mod, only: soil_parameter_type
  use cable_def_types_mod, only: veg_parameter_type
  use cable_def_types_mod, only: soil_snow_type
  use cable_def_types_mod, only: bgc_pool_type
  use cable_def_types_mod, only: canopy_type
  use cable_def_types_mod, only: roughness_type
  use cable_def_types_mod, only: radiation_type
  use cable_def_types_mod, only: balances_type
  use cable_def_types_mod, only: mvtype, mstype

  use cable_io_vars_module, only: latitude, longitude
  use cable_io_vars_module, only: landpt_global
  use cable_io_vars_module, only: patch

  use cable_netcdf_mod, only: CABLE_NETCDF_INT
  use cable_netcdf_mod, only: CABLE_NETCDF_FLOAT

  implicit none
  private

  public :: cable_restart_write

contains

  subroutine cable_restart_write(current_time, soil, veg, ssnow, canopy, rough, rad, bgc, bal, met)
    real, intent(in) :: current_time !! Current simulation time
    type(met_type), intent(in) :: met !! Meteorological data
    type(soil_parameter_type), intent(in) :: soil !! Soil parameters
    type(veg_parameter_type), intent(in) :: veg !! Vegetation parameters
    type(soil_snow_type), intent(in) :: ssnow !! Soil and snow variables
    type(bgc_pool_type), intent(in) :: bgc !! Carbon pool variables
    type(canopy_type), intent(in) :: canopy !! Vegetation variables
    type(roughness_type), intent(in) :: rough !! Roughness variables
    type(radiation_type), intent(in) :: rad !! Radiation variables
    type(balances_type), intent(in) :: bal !! Energy and water balance variables

    call cable_restart_write_time(current_time)

    call cable_restart_variable_write( &
      var_name="longitude", &
      var_dims=["mland"], &
      data=longitude, &
      var_type=CABLE_NETCDF_FLOAT, &
      long_name="", &
      units="degrees_east" &
    )

    call cable_restart_variable_write( &
      var_name="latitude", &
      var_dims=["mland"], &
      data=latitude, &
      var_type=CABLE_NETCDF_FLOAT, &
      long_name="", &
      units="degrees_north" &
    )

    call cable_restart_variable_write( &
      var_name="nap", &
      var_dims=["mland"], &
      data=landpt_global(:)%nap, &
      var_type=CABLE_NETCDF_FLOAT, &
      long_name="Number of active patches", &
      units="" &
    )

    call cable_restart_variable_write_darray( &
      var_name="patchfrac", &
      var_dims=["mp"], &
      data=patch(:)%frac, &
      var_type=CABLE_NETCDF_FLOAT, &
      long_name="Fraction of vegetated grid cell area occupied by a vegetation/soil patch", &
      units="" &
    )

    call cable_restart_variable_write( &
      var_name="mvtype", &
      data=[mvtype], &
      var_type=CABLE_NETCDF_INT, &
      long_name="Number of vegetation types", &
      units="" &
    )

    call cable_restart_variable_write( &
      var_name="mstype", &
      data=[mstype], &
      var_type=CABLE_NETCDF_INT, &
      long_name="Number of soil types", &
      units="" &
    )

    call cable_restart_variable_write_darray( &
      var_name="tgg", &
      var_dims=["mp","soil"], &
      data=ssnow%tgg, &
      var_type=CABLE_NETCDF_FLOAT, &
      long_name="Average layer soil temperature", &
      units="K" &
    )

    call cable_restart_variable_write_darray( &
      var_name="wb", &
      var_dims=["mp", "soil"], &
      data=ssnow%wb, &
      var_type=CABLE_NETCDF_FLOAT, &
      long_name="Average layer volumetric soil moisture", &
      units="vol/vol" &
    )

    call cable_restart_variable_write_darray( &
      var_name="wbice", &
      var_dims=["mp","soil"], &
      data=ssnow%wbice, &
      var_type=CABLE_NETCDF_FLOAT, &
      long_name="Average layer volumetric soil ice", &
      units="vol/vol" &
    )

    call cable_restart_variable_write_darray( &
      var_name="tss", &
      var_dims=["mp"], &
      data=ssnow%tss, &
      var_type=CABLE_NETCDF_FLOAT, &
      long_name="Combined soil/snow temperature", &
      units="K" &
    )

    call cable_restart_variable_write_darray( &
      var_name="albsoilsn", &
      var_dims=["mp", "rad"], &
      data=ssnow%albsoilsn, &
      var_type=CABLE_NETCDF_FLOAT, &
      long_name="Combined soil/snow albedo", &
      units="-" &
    )

    call cable_restart_variable_write_darray( &
      var_name="rtsoil", &
      var_dims=["mp"], &
      data=ssnow%rtsoil, &
      var_type=CABLE_NETCDF_FLOAT, &
      long_name="Turbulent resistance for soil", &
      units="??" &
    )

    call cable_restart_variable_write_darray( &
      var_name="gammzz", &
      var_dims=["mp","soil"], &
      data=ssnow%gammzz, &
      var_type=CABLE_NETCDF_FLOAT, &
      long_name="Heat capacity for each soil layer", &
      units="J/kg/C" &
    )

    call cable_restart_variable_write_darray( &
      var_name="runoff", &
      var_dims=["mp"], &
      data=ssnow%runoff, &
      var_type=CABLE_NETCDF_FLOAT, &
      long_name="Total runoff", &
      units="mm/timestep" &
    )

    call cable_restart_variable_write_darray( &
      var_name="rnof1", &
      var_dims=["mp"], &
      data=ssnow%rnof1, &
      var_type=CABLE_NETCDF_FLOAT, &
      long_name="Surface runoff", &
      units="mm/timestep" &
    )

    call cable_restart_variable_write_darray( &
      var_name="rnof2", &
      var_dims=["mp"], &
      data=ssnow%rnof2, &
      var_type=CABLE_NETCDF_FLOAT, &
      long_name="Subsurface runoff", &
      units="mm/timestep" &
    )

    call cable_restart_variable_write_darray( &
      var_name="tggsn", &
      var_dims=["mp","snow"], &
      data=ssnow%tggsn, &
      var_type=CABLE_NETCDF_FLOAT, &
      long_name="Average layer snow temperature", &
      units="K" &
    )

    call cable_restart_variable_write_darray( &
      var_name="ssdnn", &
      var_dims=["mp"], &
      data=ssnow%ssdnn, &
      var_type=CABLE_NETCDF_FLOAT, &
      long_name="Average snow density", &
      units="kg/m^3" &
    )

    call cable_restart_variable_write_darray( &
      var_name="ssdn", &
      var_dims=["mp","snow"], &
      data=ssnow%ssdn, &
      var_type=CABLE_NETCDF_FLOAT, &
      long_name="Average layer snow density", &
      units="kg/m^3" &
    )

    call cable_restart_variable_write_darray( &
      var_name="snowd", &
      var_dims=["mp"], &
      data=ssnow%snowd, &
      var_type=CABLE_NETCDF_FLOAT, &
      long_name="Liquid water equivalent snow depth", &
      units="mm" &
    )

    call cable_restart_variable_write_darray( &
      var_name="snage", &
      var_dims=["mp"], &
      data=ssnow%snage, &
      var_type=CABLE_NETCDF_FLOAT, &
      long_name="Snow age", &
      units="??" &
    )

    call cable_restart_variable_write_darray( &
      var_name="smass", &
      var_dims=["mp","snow"], &
      data=ssnow%smass, &
      var_type=CABLE_NETCDF_FLOAT, &
      long_name="Average layer snow mass", &
      units="kg/m^2" &
    )

    call cable_restart_variable_write_darray( &
      var_name="sdepth", &
      var_dims=["mp", "snow"], &
      data=ssnow%sdepth, &
      var_type=CABLE_NETCDF_FLOAT, &
      long_name="Snow layer depth", &
      units="m" &
    )

    call cable_restart_variable_write_darray( &
      var_name="osnowd", &
      var_dims=["mp"], &
      data=ssnow%osnowd, &
      var_type=CABLE_NETCDF_FLOAT, &
      long_name="Previous time step snow depth in water equivalent", &
      units="mm" &
    )

    call cable_restart_variable_write_darray( &
      var_name="isflag", &
      var_dims=["mp"], &
      data=ssnow%isflag, &
      var_type=CABLE_NETCDF_INT, &
      long_name="Snow layer scheme flag", &
      units="-" &
    )

    call cable_restart_variable_write_darray( &
      var_name="cansto", &
      var_dims=["mp"], &
      data=canopy%cansto, &
      var_type=CABLE_NETCDF_FLOAT, &
      long_name="Canopy surface water storage", &
      units="mm" &
    )

    call cable_restart_variable_write_darray( &
      var_name="ghflux", &
      var_dims=["mp"], &
      data=canopy%ghflux, &
      var_type=CABLE_NETCDF_FLOAT, &
      long_name="????", &
      units="W/m^2?" &
    )

    call cable_restart_variable_write_darray( &
      var_name="sghflux", &
      var_dims=["mp"], &
      data=canopy%sghflux, &
      var_type=CABLE_NETCDF_FLOAT, &
      long_name="????", &
      units="W/m^2?" &
    )

    call cable_restart_variable_write_darray( &
      var_name="ga", &
      var_dims=["mp"], &
      data=canopy%ga, &
      var_type=CABLE_NETCDF_FLOAT, &
      long_name="Ground heat flux", &
      units="W/m^2" &
    )

    call cable_restart_variable_write_darray( &
      var_name="dgdtg", &
      var_dims=["mp"], &
      data=canopy%dgdtg, &
      var_type=CABLE_NETCDF_FLOAT, &
      long_name="Derivative of ground heat flux wrt soil temperature", &
      units="W/m^2/K" &
    )

    call cable_restart_variable_write_darray( &
      var_name="fev", &
      var_dims=["mp"], &
      data=canopy%fev, &
      var_type=CABLE_NETCDF_FLOAT, &
      long_name="Latent heat flux from vegetation", &
      units="W/m^2" &
    )

    call cable_restart_variable_write_darray( &
      var_name="fes", &
      var_dims=["mp"], &
      data=canopy%fes, &
      var_type=CABLE_NETCDF_FLOAT, &
      long_name="Latent heat flux from soil", &
      units="W/m^2" &
    )

    call cable_restart_variable_write_darray( &
      var_name="fhs", &
      var_dims=["mp"], &
      data=canopy%fhs, &
      var_type=CABLE_NETCDF_FLOAT, &
      long_name="Sensible heat flux from soil", &
      units="W/m^2" &
    )

    call cable_restart_variable_write_darray( &
      var_name="cplant", &
      var_dims=["mp", "plant_carbon_pools"], &
      data=bgc%cplant, &
      var_type=CABLE_NETCDF_FLOAT, &
      long_name="Plant carbon stores", &
      units="gC/m^2" &
    )

    call cable_restart_variable_write_darray( &
      var_name="csoil", &
      var_dims=["mp", "soil_carbon_pools"], &
      data=bgc%csoil, &
      var_type=CABLE_NETCDF_FLOAT, &
      long_name="Soil carbon stores", &
      units="gC/m^2" &
    )

    call cable_restart_variable_write_darray( &
      var_name="wbtot0", &
      var_dims=["mp"], &
      data=bal%wbtot0, &
      var_type=CABLE_NETCDF_FLOAT, &
      long_name="Initial time step soil water total", &
      units="mm" &
    )

    call cable_restart_variable_write_darray( &
      var_name="osnowd0", &
      var_dims=["mp"], &
      data=bal%osnowd0, &
      var_type=CABLE_NETCDF_FLOAT, &
      long_name="Initial time step snow water total", &
      units="mm" &
    )

    call cable_restart_variable_write_darray( &
      var_name="albedo", &
      var_dims=["mp","rad"], &
      data=rad%albedo, &
      var_type=CABLE_NETCDF_FLOAT, &
      long_name="Albedo for shortwave and NIR radiation", &
      units="-" &
    )

    call cable_restart_variable_write_darray( &
      var_name="trad", &
      var_dims=["mp"], &
      data=rad%trad, &
      var_type=CABLE_NETCDF_FLOAT, &
      long_name="Surface radiative temperature (soil/snow/veg inclusive)", &
      units="K" &
    )

    call cable_restart_variable_write_darray( &
      var_name="iveg", &
      var_dims=["mp"], &
      data=veg%iveg, &
      var_type=CABLE_NETCDF_INT, &
      long_name="Vegetation type", &
      units="-" &
    )

    call cable_restart_variable_write_darray( &
      var_name="isoil", &
      var_dims=["mp"], &
      data=soil%isoilm, &
      var_type=CABLE_NETCDF_INT, &
      long_name="Soil type", &
      units="-" &
    )

    call cable_restart_variable_write( &
      var_name="zse", &
      var_dims=["soil"], &
      data=soil%zse, &
      var_type=CABLE_NETCDF_FLOAT, &
      long_name="Depth of each soil layer", &
      units="m" &
    )

    call cable_restart_variable_write_darray( &
      var_name="albsoil", &
      var_dims=["mp","rad"], &
      data=soil%albsoil, &
      var_type=CABLE_NETCDF_FLOAT, &
      long_name="Soil reflectance", &
      units="-" &
    )

    call cable_restart_variable_write_darray( &
      var_name="GWwb", &
      var_dims=["mp"], &
      data=ssnow%GWwb, &
      var_type=CABLE_NETCDF_FLOAT, &
      long_name="Groundwater water content", &
      units="mm3/mm3" &
    )

    if (cable_user%soil_struc == 'sli' .or. cable_user%fwsoil_switch == 'Haverd2013') then

      call cable_restart_variable_write_darray( &
        var_name="gamma", &
        var_dims=["mp"], &
        data=veg%gamma, &
        var_type=CABLE_NETCDF_FLOAT, &
        long_name="Parameter in root efficiency function (Lai and Katul 2000)", &
        units="-" &
      )

    end if

    if (cable_user%soil_struc == 'sli') then

      call cable_restart_variable_write_darray( &
        var_name="nhorizons", &
        var_dims=["mp"], &
        data=soil%nhorizons, &
        var_type=CABLE_NETCDF_INT, &
        long_name="Number of soil horizons", &
        units="-" &
      )

      call cable_restart_variable_write_darray( &
        var_name="ishorizon", &
        var_dims=["mp"], &
        data=soil%ishorizon, &
        var_type=CABLE_NETCDF_INT, &
        long_name="Horizon number", &
        units="-" &
      )

      call cable_restart_variable_write_darray( &
        var_name="clitt", &
        var_dims=["mp"], &
        data=veg%clitt, &
        var_type=CABLE_NETCDF_FLOAT, &
        long_name="Litter layer carbon content", &
        units="tC/ha" &
      )

      call cable_restart_variable_write_darray( &
        var_name="ZR", &
        var_dims=["mp"], &
        data=veg%ZR, &
        var_type=CABLE_NETCDF_FLOAT, &
        long_name="Maximum rooting depth", &
        units="cm" &
      )

      call cable_restart_variable_write_darray( &
        var_name="F10", &
        var_dims=["mp"], &
        data=veg%F10, &
        var_type=CABLE_NETCDF_FLOAT, &
        long_name="Fraction of roots in top 10 cm", &
        units="-" &
      )

      call cable_restart_variable_write_darray( &
        var_name="S", &
        var_dims=["mp"], &
        data=ssnow%S, &
        var_type=CABLE_NETCDF_FLOAT, &
        long_name="Fractional soil moisture content relative to saturated value", &
        units="-" &
      )

      call cable_restart_variable_write_darray( &
        var_name="Tsoil", &
        var_dims=["mp"], &
        data=ssnow%Tsoil, &
        var_type=CABLE_NETCDF_FLOAT, &
        long_name="Soil temperature", &
        units="degC" &
      )

      call cable_restart_variable_write_darray( &
        var_name="snowliq", &
        var_dims=["mp"], &
        data=ssnow%snowliq, &
        var_type=CABLE_NETCDF_FLOAT, &
        long_name="Liquid water content of snowpack", &
        units="mm" &
      )

      call cable_restart_variable_write_darray( &
        var_name="sconds", &
        var_dims=["mp"], &
        data=ssnow%sconds, &
        var_type=CABLE_NETCDF_FLOAT, &
        long_name="Thermal conductivity of snowpack", &
        units="W/m/K" &
      )

      call cable_restart_variable_write_darray( &
        var_name="h0", &
        var_dims=["mp"], &
        data=ssnow%h0, &
        var_type=CABLE_NETCDF_FLOAT, &
        long_name="Pond height above soil", &
        units="m" &
      )

      call cable_restart_variable_write_darray( &
        var_name="nsnow", &
        var_dims=["mp"], &
        data=ssnow%nsnow, &
        var_type=CABLE_NETCDF_INT, &
        long_name="Number of snow layers", &
        units="-" &
      )

      call cable_restart_variable_write_darray( &
        var_name="Tsurface", &
        var_dims=["mp"], &
        data=ssnow%Tsurface, &
        var_type=CABLE_NETCDF_FLOAT, &
        long_name="Soil or snow surface temperature", &
        units="degC" &
      )

    end if

  end subroutine cable_restart_write

end module
