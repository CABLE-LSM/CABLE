#!/usr/bin/env python3
"""
Create a land mask file on given latitudes and longitudes

.. moduleauthor:: Matthias Cuntz

History
    * Written May 2024 by Matthias Cuntz (mc (at) macu (dot) de)
      from split_landmask.py and extract_latlon_restart.py

"""
import argparse
import os
import numpy as np
import xarray as xr
import cablepop as cp


def generate_parser():
    """Generate the argument parser"""
    args = argparse.ArgumentParser(
        description=("Generate land mask for CABLE-POP with only latitude,longitude"
                     " set on given latitudes,longitudes."))
    args.add_argument("global_landmask_file", type=str,
                      help="Initial global landmask file")
    args.add_argument("latlon", type=str,
                      help=('latitude,longitude of selected grid point,'
                            ' or file with several latitude,longitude'
                            ' coordinates.'))
    args.add_argument("outfile", type=str,
                      help="Filename of output file")
    args.add_argument("-l", "--land", action='store_true', default=False,
                      help=('Set only latitudes/longitudes that are on land'
                            ' in the global_landmask_file '
                            ' (default: sets also if in ocean)'))

    return args


def create_run_mask(ilandmask, lats, lons):
    """
    Create the sub-land mask on given latitudes and longitudes.

    Parameters
    ----------
    ilandmask : xarray.DataArray
        Initial land mask variable. It must have:
           - 1 for land
           - missing value for ocean
    lats : numpy.array
        Latitudes
    lons : numpy.array
        Longitudes

    Returns
    -------
    xarray.DataArray
        New mask created with value 1 for a land point, missing value otherwise.

    """
    # get lat/lon indexes
    inlats = ilandmask['latitude'].data
    inlons = ilandmask['longitude'].data
    iinlats = []
    iinlons = []
    nlatlon = lats.size
    for ii in range(nlatlon):
        iinlats.append(cp.closest(inlats, lats[ii]))
        iinlons.append(cp.closest(inlons, lons[ii]))

    # set lat/lon indexes to 1, otherwise 0
    olandmask = ilandmask.copy()
    olandmask.data[:] = 0
    olandmask.data[iinlats, iinlons] = 1

    return olandmask


def main():

    parser = generate_parser()
    args = parser.parse_args()

    # Read in the global landmask
    landmask_in = xr.open_dataset(args.global_landmask_file)["land"]

    # Read latitudes/longitudes
    if os.path.exists(args.latlon):
        dat = np.genfromtxt(args.latlon, delimiter=',')
        if len(dat) > 1:
            lats = dat[:, 0]
            lons = dat[:, 1]
        else:
            lats = np.array([dat[0]])
            lons = np.array([dat[1]])
    else:
        lats, lons = args.latlon.split(',')
        lats = np.array([float(lats)])
        lons = np.array([float(lons)])

    # Check that latitudes/longitudes are within global extent
    bbox = [-90., 90., -180., 180]
    bbox[0] = landmask_in['latitude'].data.min()
    bbox[1] = landmask_in['latitude'].data.max()
    bbox[2] = landmask_in['longitude'].data.min()
    bbox[3] = landmask_in['longitude'].data.max()
    for ll in lats:
        assert (ll >= bbox[0]) and (ll <= bbox[1]), (
            f'Latitude {ll} not within extent [{bbox[0]}, {bbox[1]}].')
    for ll in lons:
        assert (ll >= bbox[2]) and (ll <= bbox[3]), (
            f'Longitude {ll} not within extent [{bbox[2]}, {bbox[3]}].')

    # Create sub-land mask
    run_mask = create_run_mask(landmask_in, lats, lons)

    # Check if points on land
    if args.land:
        run_mask.data = run_mask.data == landmask_in.data

    # Write the mask
    # Mask has 1 over land points and 0 everywhere else.
    run_mask.to_netcdf(args.outfile,
                       encoding={run_mask.name:
                                 {"dtype": np.short, "_FillValue": 0}})

    return


if __name__ == "__main__":
    main()
