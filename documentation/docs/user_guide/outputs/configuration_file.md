# Output configuration file reference

The output configuration file is a YAML file with up to four sections: `streams`, `variables`, `groups` and `modules`. A **stream** is one output file with one write frequency. Every **variable** is sent to one stream, with an aggregation (how it is combined over the write period) and optionally a reduction (how tiles are combined into a grid cell value). **Groups** and **modules** send many variables at once with the same settings. Names can be built from templates, for example `{field_name}_{aggregation}`. The two rules at the end of this page prevent two variables clashing in one file.

## Streams

```yaml
streams:
  1:
    file_name: cable_output.nc
    frequency: monthly
```

Each stream is introduced by an integer, which is how variables refer to it. The numbers need not be consecutive.

| Setting | Required | Default | Meaning |
|---|---|---|---|
| `file_name` | yes | | File the stream is written to. |
| `frequency` | yes | | How often the file is written, and so the period over which variables are aggregated: `timestep`, `3hrly`, `daily`, `monthly` or `yearly`. |
| `netcdf_name` | no | `"{field_name}"` | Template for the variable names in the file. Individual entries can override it. See [string substitution](#string-substitution). |
| `shuffle` | no | `true` | Apply the NetCDF shuffle filter (only has an effect together with compression). |
| `compression_level` | no | `1` | Deflate compression level from 1 (fastest) to 9 (smallest). `0` turns compression off. |
| `separate_file_per_variable` | no | `false` | Write each variable to a file of its own. See [separate files](#separate-files). |
| `metadata` | no | none | Global attributes for the file. Values can use `{frequency}`, `{start_date}` and `{end_date}`. |

A stream that no variable is sent to does not produce a file.

### File format and compression

With `compression_level` above 0 (the default), the file is NetCDF-4 and each variable is compressed. With `0` the file is the classic NetCDF format, which cannot be compressed. Compression is not yet supported when writing with ParallelIO. Set `compression_level: 0` there.

### Separate files

With `separate_file_per_variable: true`, each variable goes to its own file, named after its NetCDF name with `.nc` added, in the directory part of `file_name`. The `file_name` itself is not used for a file. For example `file_name: out/unused.nc` with the NetCDF name template `{field_name}_{aggregation}` writes `out/GPP_mean.nc` and so on. NetCDF names used this way cannot contain `/`.

## Variables

```yaml
variables:
  - name: GPP
    stream: 1
    aggregation: mean
    reduction: grid_cell_average
```

| Setting | Required | Default | Meaning |
|---|---|---|---|
| `name` | yes | | The variable, as listed in the [variable catalogue](variable_catalogue.md). |
| `stream` | yes | | The stream to write it to. A variable entry goes to one stream; list it again to write it to another. |
| `aggregation` | yes | | How the variable is combined over the stream's period: `mean`, `sum`, `max`, `min` or `instant` (the value at the end of the period). |
| `reduction` | no | `none` | How tiles are combined into a grid cell value: see below. |
| `netcdf_name` | no | the stream's | Name of the variable in the file, as a template. |
| `metadata` | no | none | Attributes for the variable. Values are templates. |

Asking for a variable the current model configuration cannot provide (for example a groundwater variable with the groundwater model off) is an error that says so.

### Aggregation

The model value is sampled at each step the model updates it (every time step for most variables; a few are updated once a day, see the catalogue), and combined over the stream's period: `mean` is the average of the samples, `sum` their total, `max` and `min` the extremes, and `instant` the latest value. `instant` with `frequency: timestep` writes the model value as it is.

A variable the model updates daily cannot be written to a stream with a finer frequency than daily.

### Reduction

CABLE stores most variables for each tile (patch) of each grid cell. The reduction chooses what is written:

| Reduction | Written |
|---|---|
| `none` | One value for every tile, so the output has a tile dimension. |
| `grid_cell_average` | The average over the tiles, weighted by tile area. Not defined for integer variables. |
| `first_tile_on_cell` | The value of the first tile of the cell. Suited to values that are the same on every tile, such as forcing. |
| `dominant_tile` | The value of the tile with the largest area. If tiles are equal, the first of them. |

Only variables defined per tile can be reduced. Scalars and variables on other dimensions (soil layers only, for example) must use `none`.

### Parameters

Variables that do not change in time (the catalogue marks them as parameters) are written once, near the start of the run, and have no time axis. The stream's frequency and the aggregation do not affect them, although an aggregation must still be given.

## Groups and modules

```yaml
groups:
  - name: flux
    stream: 1
    aggregation: mean
    reduction: grid_cell_average

modules:
  - name: biogeophysics
    stream: 1
    aggregation: sum
```

Groups and modules take the same settings as a variable entry, and apply them to every variable they contain. A module contains groups, and the catalogue lists which group and module each variable belongs to. An unavailable member is simply left out of a group or module. Entries are additive: a variable that is reached by a group and also listed by name is written twice, so the two entries must be told apart (see the rules).

Modules are expanded first, then groups, then variables, and variables appear in the file in that order, and within each in catalogue order.

## String substitution

Templates for `netcdf_name` and for the values of `metadata` can contain these names in braces, which are replaced when the file is read:

| Name | Replaced by |
|---|---|
| `{field_name}` | The variable's name, as in the catalogue. |
| `{frequency}` | The stream's frequency. |
| `{aggregation}` | The aggregation of the entry. |
| `{reduction}` | The reduction of the entry. |
| `{start_date}`, `{end_date}` | First and last date of the run, as `YYYY-MM-DD`. |

The metadata of a stream cannot use `{field_name}`, `{aggregation}` or `{reduction}`, because a stream has no single value for them. An unknown name in braces is an error.

## Attributes in the file

Each variable carries `units` and `long_name` from the catalogue, and a `cell_methods` attribute worked out from its aggregation and reduction, for example `area: mean time: mean`. These cannot be set in `metadata`. Attributes you give in `metadata` are added after them. `time` values are the middle of each write period, or the time of the step itself for `frequency: timestep`.

## Rules

Different entries are always additive: every entry adds variables to its stream. Two rules keep the result consistent, and a file that breaks them stops the run at start-up with a message for each.

1. **A stream cannot contain two variables with the same NetCDF name** (after substitution). Writing a variable twice with different aggregations needs a `netcdf_name` that tells them apart, for example `"{field_name}_{aggregation}"`.
2. **A stream cannot contain two identical variables**: the same variable with the same aggregation and reduction is a duplicate even if the NetCDF names differ.

Other mistakes are also found when the model starts and reported together: unknown keys, unknown variables, groups, aggregations, reductions or frequencies, an undefined stream number, two streams writing the same file, a `compression_level` outside 0 to 9, a variable written more often than the model updates it, and attributes that clash with the ones set automatically.
