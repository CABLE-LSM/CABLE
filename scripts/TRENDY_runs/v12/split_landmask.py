#!/usr/bin/env python3
"""
Program to split the global landmask into submasks

The split ensures we pick land points from various regions
to avoid having all land points frozen if picking all
of Greenland, for example.

History
    * Written Mar 2024 by Lachlan Whyborn
    * PEP 8 compliance, May 2024, Matthias Cuntz
    * Added file docstring with history, May 2024, Matthias Cuntz
    * More explanatory help messages, May 2024, Matthias Cuntz
    * Removed bug if extent given, May 2024, Matthias Cuntz
    * Output short datatype instead of byte, May 2024, Matthias Cuntz

"""
import argparse
import numpy as np
import xarray as xr
import os


def generate_parser():
    """
    Generate the argument parser

    """
    args = argparse.ArgumentParser(
        description=("Generate land masks for CABLE to split a global run"
                     " in a series of runs."))
    args.add_argument("global_landmask_file",
                      help="Initial global landmask file")
    args.add_argument("nmasks", type=int,
                      help="Number of output sub-landmasks")
    args.add_argument("outpath",
                      help=("Root path for output files. Files will be:"
                            "{outpath}/run{num}/landmask/landmask{num}.nc"))
    args.add_argument("-e", "--extent", type=str, default="global",
                      help=('Extent of latitude and longitude in output file:'
                            ' "global" (default) or'
                            ' "lon_min,lon_max,lat_min,lat_max"'))

    return args


def main():

    parser = generate_parser()
    args = parser.parse_args()

    # Check if extent is global or regional
    bbox = [-180.0, 180.0, -90.0, 90.0]
    if args.extent != "global":
        bbox = [ float(coord) for coord in args.extent.split(",") ]

    # Read in the landmask
    landmask_in = xr.open_dataset(args.global_landmask_file)["land"]

    # Cut the data to the extent required. We need to sort
    # the latitude in ascending order first and then revert if needed
    is_ascending = landmask_in["latitude"][0] < landmask_in["latitude"][1]
    landmask_cut = landmask_in.sortby("latitude").sel(
        latitude=slice(bbox[2], bbox[3]), longitude=slice(bbox[0], bbox[1]))
    landmask_cut = landmask_cut.sortby("latitude", ascending=is_ascending)

    # Create composite mask for all the runs.
    run_mask = create_run_mask(landmask_cut, args.nmasks)

    # Reindex the local run mask to be the same as the global land
    # mask. Set the fill value to 0 i.e. not active.
    run_mask = run_mask.reindex_like(landmask_in, fill_value=0)

    # Save mask to file to allow for an easy check.
    os.makedirs(f"{args.outpath}", exist_ok=True)
    run_mask.to_netcdf(f"{args.outpath}/check_landmask.nc")

    # Write the mask for each run to a separate file
    for i in range(args.nmasks):
        # Write files for each mask. Each mask has 1 over the land
        # points for that mask and 0 everywhere else.
        # Fill missing with 0. and save as bytes.
        landmask_per_run = landmask_in.where(run_mask==i + 1, 0)
        os.makedirs(f"{args.outpath}/run{i+1}/landmask", exist_ok=True)
        landmask_per_run.to_netcdf(
            f"{args.outpath}/run{i+1}/landmask/landmask{i+1}.nc",
            encoding={landmask_per_run.name:
                      {"dtype": np.short, "_FillValue": 0}})


def create_run_mask(
    landmask_in: xr.DataArray,
    nmasks: int) -> xr.DataArray:
    """
    Create the full mask per run. The mask has a value N for points
    done by run N.

    Parameters
    ----------
    landmask_in : xarray.DataArray
        Initial land mask variable. It must have:
           - 1 for land
           - missing value for ocean
    nmasks : int
        Number of runs we want to do.

    Returns
    -------
    xarray.DataArray
        New mask created with value N for a land point for the run N.

    """
    # Number of land points total and
    # number of land points per processor
    ntot = landmask_in.sum().astype(np.int32).item()
    npermask = np.ceil(ntot / nmasks).astype(np.int32).item()

    # Create the repetitive pattern to select the land points
    # for each processor
    run_mask = np.arange(nmasks)
    run_mask = np.tile(run_mask, npermask)[:ntot]

    # New land mask
    mask_land_flat = landmask_in.values.flat
    landmask_in.values.flat[~np.isnan(mask_land_flat)] += run_mask

    return landmask_in


if __name__ == "__main__":
    main()
