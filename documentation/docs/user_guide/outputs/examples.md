# Output configuration examples

These configurations show how the settings of the [configuration file](configuration_file.md) combine. Each one is a complete file that CABLE accepts. File names are relative to the directory the model is run from.

## One stream, a few variables

The smallest useful file: daily means of two variables.

```yaml
streams:
  1:
    file_name: cable_daily.nc
    frequency: daily

variables:
  - name: GPP
    stream: 1
    aggregation: mean
    reduction: grid_cell_average

  - name: Qle
    stream: 1
    aggregation: mean
    reduction: grid_cell_average
```

## Several streams

To keep physical and biochemical variables apart, send each module to its own stream. This writes the monthly means of both to two files.

```yaml
streams:
  1:
    file_name: cable_biogeophysics_monthly.nc
    frequency: monthly
  2:
    file_name: cable_biogeochemistry_monthly.nc
    frequency: monthly

modules:
  - name: biogeophysics
    stream: 1
    aggregation: mean
    reduction: grid_cell_average

  - name: biogeochemistry
    stream: 2
    aggregation: mean
    reduction: grid_cell_average
```

Variables the model cannot provide in the configuration being run (for example the groundwater variables when the groundwater model is off) are left out of a module without any message.

## Adding more output to a stream

Entries are additive. This sends the daily means of two groups to one file.

```yaml
streams:
  1:
    file_name: cable_water_energy.nc
    frequency: daily

groups:
  - name: flux
    stream: 1
    aggregation: mean
    reduction: grid_cell_average

  - name: balances
    stream: 1
    aggregation: mean
    reduction: grid_cell_average
```

## Different frequencies for different variables

Monthly means of the fluxes, and in a second file the three-hourly instantaneous values of two of them, to look at the diurnal cycle.

```yaml
streams:
  1:
    file_name: cable_flux_monthly.nc
    frequency: monthly
  2:
    file_name: cable_flux_diurnal.nc
    frequency: 3hrly

groups:
  - name: flux
    stream: 1
    aggregation: mean
    reduction: grid_cell_average

variables:
  - name: Qle
    stream: 2
    aggregation: instant
    reduction: grid_cell_average

  - name: Qh
    stream: 2
    aggregation: instant
    reduction: grid_cell_average
```

## Different aggregations in the same stream

To write both the three-hourly sums and the instantaneous values of some fluxes to one file, the variables must have different NetCDF names (rule 1). Here the group gets a name template and the variables get explicit names. `Qle` and `Qh` belong to the `flux` group, so each is written twice, once as a sum and once as an instantaneous value.

```yaml
streams:
  1:
    file_name: cable_flux_3hrly.nc
    frequency: 3hrly

groups:
  - name: flux
    stream: 1
    aggregation: sum
    reduction: grid_cell_average
    netcdf_name: "{field_name}_{aggregation}"

variables:
  - name: Qle
    stream: 1
    aggregation: instant
    reduction: grid_cell_average
    netcdf_name: Qle_instant

  - name: Qh
    stream: 1
    aggregation: instant
    reduction: grid_cell_average
    netcdf_name: Qh_instant
```

If the group had no `netcdf_name`, the file would be rejected at start-up: it would contain `Qle` twice (rule 1). If the variables were listed with `aggregation: sum` and a different name, it would be rejected for writing the same variable twice with the same settings (rule 2).

## Choosing how tiles are combined

The same variable can be written with different reductions in one file. Here the area-weighted mean of the air temperature, and its value on the largest tile.

```yaml
streams:
  1:
    file_name: cable_tair_reductions.nc
    frequency: monthly

variables:
  - name: Tair
    stream: 1
    aggregation: mean
    reduction: grid_cell_average
    netcdf_name: Tair_cell_mean

  - name: Tair
    stream: 1
    aggregation: mean
    reduction: dominant_tile
    netcdf_name: Tair_dominant_tile

  - name: Tair
    stream: 1
    aggregation: mean
    reduction: none
    netcdf_name: Tair_per_tile
```

## Compressing a stream

For data kept for a long time, raise the compression. This applies the shuffle filter and level 4 deflate. The file is NetCDF-4.

```yaml
streams:
  1:
    file_name: cable_compressed.nc
    frequency: daily
    shuffle: true
    compression_level: 4

modules:
  - name: biogeophysics
    stream: 1
    aggregation: mean
    reduction: grid_cell_average
```

## Splitting a stream into separate files

Many intercomparison projects want one variable per file. This writes each variable of the `flux` group to its own file, named after the NetCDF name. The `file_name` is not written to, but its directory is where the files go, so the files below are written to `separate/`, for example `separate/Qle_mean.nc`. The directory must exist.

```yaml
streams:
  1:
    file_name: separate/unused.nc
    frequency: daily
    separate_file_per_variable: true
    netcdf_name: "{field_name}_{aggregation}"
    compression_level: 0

groups:
  - name: flux
    stream: 1
    aggregation: mean
    reduction: grid_cell_average
```

## Global attributes and dates in names

Attribute values and names can use the substitutions described in the reference. Here the experiment name and the run dates are recorded in the file, and the variable names carry the aggregation.

```yaml
streams:
  1:
    file_name: cable_soil_monthly.nc
    frequency: monthly
    netcdf_name: "{field_name}_{aggregation}"
    metadata:
      model: CABLE
      experiment: "example run {start_date} to {end_date}"

groups:
  - name: soil
    stream: 1
    aggregation: mean
    reduction: grid_cell_average
    metadata:
      comment: "{aggregation} of {field_name}, written {frequency}"
```

## Restart files and parameters

A restart file is written at the end of a run if `output%restart = .TRUE.` in the namelist, whatever the configuration file selects. To also record the vegetation and soil parameters used by the run, add the `parameters` module. Parameters do not change in time, so they are written once, and they are not averaged. An aggregation must still be given, and `reduction: none` keeps the tiles.

```yaml
streams:
  1:
    file_name: cable_parameters.nc
    frequency: monthly

modules:
  - name: parameters
    stream: 1
    aggregation: instant
    reduction: none
```
