# Output files

CABLE writes its output according to a **YAML output configuration file**, named in the `&cable` namelist by `filename%output_config`. The file describes *streams* (an output file with a write frequency) and the *variables* sent to them. Each variable can be averaged, summed or sampled over the write period, and can have its tiles reduced to one value per grid cell. Variables can be listed one by one or selected in bulk by group or module. The whole file is checked when the model starts, and every problem found is reported together.

```yaml
streams:
  1:
    file_name: cable_output.nc
    frequency: monthly

modules:
  - name: biogeophysics
    stream: 1
    aggregation: mean
    reduction: grid_cell_average
```

This writes the monthly mean, per grid cell, of every biogeophysical variable the current model configuration can provide, to `cable_output.nc`.

## What is on these pages

| Page | Contents |
|---|---|
| [Configuration file reference](configuration_file.md) | Every setting, the rules the file must follow, and how names and attributes are built. |
| [Examples](examples.md) | Worked configurations: several files, mixed frequencies, compression, one file per variable. |
| [Variable catalogue](variable_catalogue.md) | Every variable that can be written, its metadata, and the model variable it comes from. Generated from the source. |

## Namelist settings

Three entries of the `&cable` namelist concern output. See [cable.nml](../inputs/cable_nml.md).

| Setting | Meaning |
|---|---|
| `filename%output_config` | The YAML output configuration file. Required: a run without it stops at start-up. |
| `output%restart` | Write a restart file at the end of the run, to `filename%restart_out`. This does not depend on the configuration file. |
| `output%grid` | Layout of the output files: `'default'` follows the meteorological forcing, `'land'` writes compressed land points, `'mask'` writes a latitude/longitude grid. |

## Changes from the earlier `output%` switches

The `output%met`, `output%flux`, ... group switches, the per-variable `output%<name>` switches, `output%averaging` and all `patchout%` switches have been removed. A namelist that still sets them stops at start-up with a message pointing here. To convert one:

| Earlier setting | Now |
|---|---|
| `output%met` | the `met` group (also in the `forcing` module) |
| `output%flux`, `output%soil`, `output%snow`, `output%radiation`, `output%veg`, `output%balances` | the group of the same name (all in the `biogeophysics` module) |
| `output%carbon`, `output%casa` | the `carbon` and `casa` groups (in the `biogeochemistry` module) |
| `output%params` | the `params` group (the `parameters` module) |
| `output%<name>`, for example `output%GPP` | list the variable by name under `variables` |
| `output%averaging = 'all'` | `frequency: timestep` |
| `output%averaging = 'daily'` or `'monthly'` | `frequency: daily` or `frequency: monthly` |
| `output%averaging = 'user6'` and similar | not available. Use `frequency: 3hrly`, or another listed frequency. |
| `patchout%<name>`, `output%patch` | `reduction: none`, which keeps the individual tiles |
| the default behaviour of a variable | Mostly `aggregation: mean` with `reduction: grid_cell_average`. Parameters were `instant` with the first tile of each cell. The catalogue page lists each variable. |

A group or module now has one aggregation and one reduction for all its members, where the old switches gave each variable its own default. To keep a particular default, list that variable separately. Because the same variable may not appear twice in a stream with the same settings, give the group-level entry a `netcdf_name` template if you also list members individually (see the rules in the reference).

`output%balances` no longer exists. The energy and water balance totals printed at the end of a run now depend on `verbose` only.
