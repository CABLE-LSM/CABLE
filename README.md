# CABLE Crop branch

## Setup

Before running, set up a python environment if your system does not already have one.

```bash
virtualenv ~/cable-env
source ~/cable-env/bin/activate
pip install numpy xarray pandas netCDF4 matplotlib
```

Then edit the `scripts/run_cable_casa_crop_array.slurm` script to source this virtualenv

```bash
source ~/cable-env/bin/activate
```

### Compiling

```bash
cd offline
./build.ksh
```

## Running

```bash
./scripts/CABLE_CROP_site_runs.ksh
```

This submits `scripts/run_cable_casa_crop_array.slurm` to the queue (petrichor) which in turn
executes `scripts/run_cable_site_CNP_meta.py` which in turn executes `cable`.
