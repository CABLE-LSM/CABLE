# Meteorological forcing for offline CABLE only

The meteorological forcing is a collection of meteorological variables that need to be read into CABLE's `met` arrays when CABLE is run in offline (uncoupled) mode.

The filename is set using the CABLE namelist variable `filename%met`.
The global offline simulations use multiple met files, one for each meteorological variable recorded in compressed grid method.
The data must be in NetCDF format with the correct units, broadly conforming to the [ALMA](https://web.lmd.jussieu.fr/~polcher/ALMA/convention_output_3.html) standard.

An example input file is provided with the CABLE-AUX release or can can be found at the [NCI THREDDS server](https://geonetwork.nci.org.au/geonetwork/srv/eng/catalog.search#/metadata/f7075_4625_2374_0846).

## Model configuration and grid

Site or regional/global simulations can use either an x-y grid input file or a land only compressed grid-type input file.

Meteorological variables that conform to the ALMA standard can have:

- an x-y grid with 3 dimensions (x,y,t) or 4 dimensions (x,y,z,t), or
- a compressed land-only grid of 2 dimensions (land, t)

For the x-y grid, an "`x`" and a "`y`" dimension variable must be present, even if the simulation is only a single site/gridpoint. Additionally, single precision variables named "`latitude`" and "`longitude`" (or "`nav_lat`" and "`nav_lon`" if using ALMA formatting), both dependent on the x and y dimensions only, must be present. Both sea and land points may be included by using an integer mask(x,y) variable; a value of 1 implies a land gridpoint, anything else is assumed to be ocean.

For the single dimension land-only “compression by gathering” grid (see [here](http://www.lmd.jussieu.fr/~polcher/ALMA/dataformats.html)), a single spatial dimension is used.
An example of the NetCDF header from such a file is shown below:

```
dimensions:
    tstep = UNLIMITED ; // (1461 currently)
    land = 15238 ;
    y = 180 ;
    x = 360 ;
variables:
    float SWdown(tstep, land) ;
        SWdown:axis = "TYX" ;
        SWdown:units = "W/m^2" ;
        SWdown:long_name="Surface incident shortwave radiation" ;
        SWdown:associate = "time (nav_lat nav_lon)" ;
        SWdown:missing_value = 1.e+20f ;
    int land(land) ;
        land:compress = "y x" ;
    float nav_lat(y, x) ;
        nav_lat:units = "degrees_north" ;
        nav_lat:valid_min = -90.f ;
        nav_lat:valid_max = 90.f ;
        nav_lat:long_name = "Latitude" ;
    float nav_lon(y, x) ;
        nav_lon:units = "degrees_east" ;
        nav_lon:valid_min = -180.f ;
        nav_lon:valid_max = 180.f ;
        nav_lon:long_name = "Longitude" ;
    float time(tstep) ;
        time:units = "seconds since 1949-01-01 00:00:00" ;
        time:title = "Time" ;
        time:long_name = "Time axis" ;
        time:time_origin = " 1949-JAN-01 00:00:00" ;
```

Most of the variables will be structured as SWdown is above. The relationship between “land” and “x” and “y” inside the CABLE netcdf driver is:

```fortran
y = INT((landGrid(j)-1)/xdimsize)
x = landGrid(j) - y*xdimsize
y = y + 1
```

Note that the number of points simulated by offline CABLE is not specified in the code, but determined here using the number of land points found in the met file. If a patch dimension exists, the `patchfrac` variable gives the fraction of each vegetation patch in each land grid point. CABLE counts the number of "active" patches (patches with non-zero `patchfrac`) as the number of points for simulation.

## Variables

### Required variables

| Name     | Description                               | Units                  |
|----------|-------------------------------------------|------------------------|
| `time`   | time of each time step (double precision) | \( s \)                |
| `SWdown` | Surface incident shortwave radiation      | \( W \cdot m^{-2} \)   |
| `Tair`   | Surface incident shortwave radiation      | \( K \)                |
| `Qair`   | Near surface specific humidity            | \( kg \cdot kg^{-1} \) |
| `Rainf`  | Rainfall rate                             | \( mm \cdot s^{-1} \)  |
| `Wind` or `Wind_E` and `Wind_E` | Scalar wind speed  | \( m \cdot s^{-1} \)   |

### Optional variables

| Name        | Description | Units |
|-------------|-----------------------------------------------------------------------------------------------------|---------------------|
| `LWdown`    | Surface incident longwave radiation (can be synthesised from `Tair`)                                | \( W \dot m^{-2} \) |
| `PSurf`     | Surface air pressure (can be estimated at a fixed value based on Tair and an "`elevation`" variable | \( Pa \)            |
| `elevation` | Surface elevation (required if PSurf is not present)                                                | \( m \)             |
| `Snowf`     | Snowfall (assumed to be included in Rainf if not present)                                           | \( mm s^{-1} \)     |
| `CO2air`    | CO$_2$ concentration (assumes a fixed value, determined in the `cable.nml` namelist file)           | \( ppm \)           |

### Site specific parameters

Site-specific parameters recorded in the met file will overwrite the default values obtained during initialisation and have top priority over values specified in other files.
For example, the default value of `za` (reference height or measurement height) is 40 m; it will be overwritten by the value read from the meteorological file.
The site-specific parameters include:

| Name    | Description       | Units   |
|---------|-------------------|---------|
| `hc`    | Vegetation height | \( m \) | 
| `za`    | Reference height  | \( m \) |
| `iveg`  | Vegetation type   | \( - \) |
| `frac4` | C4 fraction       | \( - \) |

## Metadata

All variables must have a "units" string attribute and the data must be in the correct units specified in the tables above.

### Details on time and CABLE's time stepping

Meteorological input data must be continuous in time and have regular intervals.

CABLE's running period is deduced from the time period of the meteorolocal forcing data.
The time variables units attribute is of the form "`seconds since <reference_time>`", where the reference time can be any date on or before the starting date (e.g. "`2001-02-22 00:00:00`").
For example, the first value of the "`time`" variable may be 86400, in which case the actual start time might be 2001-02-22 00:00:00 + 86400 seconds; i.e. 2001-02-23 00:00:00.

The "`time`" values for regional or global simulations is assumed to be GMT, while the value for single site/grid cell simulations is assumed to be "local".
Single site simulation "`time`" values will be read as GMT if and only if a "`coordinate`" attribute is present for the "`time`" variable and set to be "GMT".
This time coordinate system will be reported in the log file.

CABLE’s time step size is calculated from the first two values of the time variable, and the run length of the simulation is decided by the length of this same variable.

## Checking ranges

It is possible to check that the ranges of the meteorological variables are physically possible using `check%ranges=.TRUE.` in the `cable.nml` namelist variable.

