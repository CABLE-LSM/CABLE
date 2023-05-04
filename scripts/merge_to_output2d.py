#!/usr/bin/env python
"""
usage: merge_to_output2d.py [-h] [-o output_netcdf] [-v] [-z] [files ...]

Merge files of pseudo-parallelised Cable or Casa output in land format to
combined 2D arrays with latitude/longitude.

positional arguments:
  files                 input netcdf files.

options:
  -h, --help            show this help message and exit
  -o output_netcdf, --outfile output_netcdf
                        output netcdf file name
                        (default: first_input-merged.nc).
  -v, --verbose         Feedback during copy (default: no feedback).
  -z, --zip             Use netCDF4 variable compression
                        (default: same format as input file).


Example
-------
  python merge_to_output2d.py -v -z -o cru_out_casa_2009_2011.nc \
      run*/outputs/cru_out_casa_2009_2011.nc


History
-------
Written  Matthias Cuntz, May 2020
             - from unpack_to_output2d.py

"""
import argparse
import sys
import numpy as np
import netCDF4 as nc
import cablepop as cp
import pyjams as pj
import time as ptime

# -------------------------------------------------------------------------
# Command line
#

ofile   = None
verbose = False
izip    = False
parser  = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=('Merge files of pseudo-parallelised Cable or Casa output'
                 ' in land format to combined 2D arrays with'
                 ' latitude/longitude.') )
parser.add_argument('-o', '--outfile', action='store',
                    default=ofile, dest='ofile', metavar='output_netcdf',
                    help=('output netcdf file name (default:'
                          ' first_input-merged.nc).'))
parser.add_argument('-v', '--verbose', action='store_true', default=verbose,
                    dest='verbose',
                    help='Feedback during copy (default: no feedback).')
parser.add_argument('-z', '--zip', action='store_true', default=izip,
                    dest='izip',
                    help=('Use netCDF4 variable compression (default:'
                          ' same format as input file).'))
parser.add_argument('ifiles', nargs='*', default=None, metavar='files',
                    help='input netcdf files.')
args    = parser.parse_args()
ofile   = args.ofile
verbose = args.verbose
izip    = args.izip
ifiles  = args.ifiles
del parser, args

if len(ifiles) == 0:
    raise IOError('Input files must be given.')

if verbose:
    tstart = ptime.time()

# -------------------------------------------------------------------------
# Copy data
#

# Input file
ifile = ifiles[0]
fi = nc.Dataset(ifile, 'r')
if verbose:
    print('Check first input file:', ifile)
ncvars = list(fi.variables.keys())
ntime = fi.dimensions['time'].size

# Output file
if ofile is None:  # Default output filename
    ofile = pj.ncio.set_output_filename(ifile, '-merged')
if verbose:
    print('Output Create output file:', ofile)
if izip:
    fo = nc.Dataset(ofile, 'w', format='NETCDF4')
else:
    if 'file_format' in dir(fi):
        fo = nc.Dataset(ofile, 'w', format=fi.file_format)
    else:
        fo = nc.Dataset(ofile, 'w', format='NETCDF3_64BIT_OFFSET')

# get latitude/longitude
if 'local_lat' in fi.variables:
    ilats = fi.variables['local_lat'][:]
    ilons = fi.variables['local_lon'][:]
else:
    ilats = fi.variables['latitude'][:]
    ilons = fi.variables['longitude'][:]
nland = ilats.size
if 'x' in fi.dimensions:  # use existing grid
    olat = fi.variables['y'][:]
    olon = fi.variables['x'][:]
    nlat = olat.size
    nlon = olon.size
else:                     # create new global grid -60 to +90 latitudes
    # 0.5, 1 degree
    dlon = np.abs(np.diff(np.unique(np.sort(ilons)))).min()
    if np.all(ilats == ilats[0]):
        # if only few cells on one latitude band
        dlat = dlon
    else:
        # 0.5, 1 degree
        dlat = np.abs(np.diff(np.unique(np.sort(ilats)))).min()
    nlat = np.rint(150. / dlat).astype(int)  # 300, 150
    nlon = np.rint(360. / dlon).astype(int)  # 720, 360
    clat = ilats.min() % 1.  # 0.0 or 0.25, 0.0 or 0.5
    clon = ilons.min() % 1.  # 0.0 or 0.25, 0.0 or 0.5
    # new lats
    olat = -60.  + clat + np.arange(nlat) / float(nlat - 1) * (150. - dlat)
    olat = olat[::-1]
    # new lons
    olon = -180. + clon + np.arange(nlon) / float(nlon - 1) * (360. - dlon)
olon2d, olat2d = np.meshgrid(olon, olat)  # new lats, lons in 2D

# Copy global attributes, adding this script
pj.ncio.copy_global_attributes(fi, fo,
                               add={'history': ptime.asctime() + ': ' +
                                    ' '.join(sys.argv)})

# Copy dimensions
pj.ncio.copy_dimensions(fi, fo,
                        removedim=['land', 'ntile'],
                        adddim={'x': nlon, 'y': nlat})

# Create static variables (independent of time)
# if 'local_lat' in fi.variables:
#     renvar = {'latitude': 'nav_lat', 'longitude': 'nav_lon'}
# else:
#     renvar = {}
renvar = {}
pj.ncio.create_variables(fi, fo, time=False, izip=izip, fill=True,
                         chunksizes=False, renamevar=renvar,
                         replacedim={'land': ('y', 'x'),
                                     'ntile': ('y', 'x')})
# create dynamic variables (time dependent)
pj.ncio.create_variables(fi, fo, time=True, izip=izip, fill=True,
                         chunksizes=False, renamevar=renvar,
                         replacedim={'land': ('y', 'x'),
                                     'ntile': ('y', 'x')})

if 'x' not in fi.variables:
    print('Create x')
    nvar = {'name': 'x',
            'dtype': ilons.dtype,
            'dimensions': ('x'),
            'units': 'degrees_east'}
    ovar = pj.ncio.create_new_variable(nvar, fo)
    ovar[:] = olon

if 'y' not in fi.variables:
    print('Create y')
    nvar = {'name': 'y',
            'dtype': ilats.dtype,
            'dimensions': ('y'),
            'units': 'degrees_north'}
    ovar = pj.ncio.create_new_variable(nvar, fo)
    ovar[:] = olat

fi.close()

if verbose:
    print('Get all indexes')
nfiles = len(ifiles)
fis = []
iidl = []
oidx = []
oidy = []
for nfile, ifile in enumerate(ifiles):
    fi = nc.Dataset(ifile, 'r')
    fis.append(fi)
    # Check time
    ntime1 = fi.dimensions['time'].size
    if ntime1 != ntime:
        for fi in fis:
            fi.close()
        fo.close()
        raise ValueError(f'Time not the same in {ifiles[0]} and in {ifile}')
    # get latitude/longitude indices
    if 'local_lat' in fi.variables:
        print('local_lat')
        ilats = fi.variables['local_lat'][:]
        ilons = fi.variables['local_lon'][:]
    else:
        ilats = fi.variables['latitude'][:]
        ilons = fi.variables['longitude'][:]
    nland = ilats.size
    if 'x' in fi.dimensions:  # use existing grid
        olat = fi.variables['y'][:]
        olon = fi.variables['x'][:]
        nlat = olat.size
        nlon = olon.size
    else:                     # create new global grid -60 to +90 latitudes
        # 0.5, 1 degree
        dlon = np.abs(np.diff(np.unique(np.sort(ilons)))).min()
        if np.all(ilats == ilats[0]):
            # if only few cells on one latitude band
            dlat = dlon
        else:
            # 0.5, 1 degree
            dlat = np.abs(np.diff(np.unique(np.sort(ilats)))).min()
        nlat = np.rint(150. / dlat).astype(int)  # 300, 150
        nlon = np.rint(360. / dlon).astype(int)  # 720, 360
        clat = ilats.min() % 1.  # 0.0 or 0.25, 0.0 or 0.5
        clon = ilons.min() % 1.  # 0.0 or 0.25, 0.0 or 0.5
        # new lats
        olat = -60.  + clat + np.arange(nlat) / float(nlat - 1) * (150. - dlat)
        olat = olat[::-1]
        # new lons
        olon = -180. + clon + np.arange(nlon) / float(nlon - 1) * (360. - dlon)
    olon2d, olat2d = np.meshgrid(olon, olat)  # new lats, lons in 2D
    lltree = cp.llKDTree(olat2d, olon2d)  # KD-tree
    iidl.append(np.arange(nland, dtype=int))  # indices of land in input grid
    fi_oidx   = np.empty(nland, dtype=int)    # indices of lon in output grid
    fi_oidy   = np.empty(nland, dtype=int)    # indices of lat in output grid
    for i in range(nland):
        iy, ix = lltree.query(ilats[i], ilons[i])
        fi_oidx[i] = ix
        fi_oidy[i] = iy
    oidx.append(fi_oidx)
    oidy.append(fi_oidy)

#
# Copy variables from in to out expanding the land dimension to y, x
#

# copy static and dynamic variables
if verbose:
    print('Copy input to output')
n = 0
for ncvar in ncvars:
    if verbose:
        tstartvar = ptime.time()
        print(f'    {ncvar}')
    n += 1
    ivar0 = fis[0].variables[ncvar]
    ovar = fo.variables[ncvar]
    if ncvar == 'longitude':
        ovar[:] = olon2d
    elif ncvar == 'latitude':
        ovar[:] = olat2d
    elif (('land' not in ivar0.dimensions) and
          ('ntile' not in ivar0.dimensions)):
        # should not be masked and all the same: check
        for f, fi in enumerate(fis):
            ivar = fi.variables[ncvar]
            if not np.all(ivar0[:] == ivar[:]):
                print('ivar0', ivar0[:])
                print('ivar', ivar[:])
                for fi in fis:
                    fi.close()
                fo.close()
                raise ValueError(f'variable {ncvar} not equal in file'
                                 f' {ifiles[0]} and file {ifiles[f]}')
        ovar[:] = ivar0[:]
    else:
        outvar = np.full(ovar.shape,
                         pj.ncio.get_fill_value_for_dtype(ivar0.dtype))
        for f, fi in enumerate(fis):
            ivar = fi.variables[ncvar]
            # read whole field, otherwise times increasing sharply
            invar = ivar[:]
            # fill in memory, then write to disk in one go
            if len(invar.shape) == 1:
                outvar[oidy, oidx] = invar[iidl]
            else:
                outvar[..., oidy, oidx] = invar[..., iidl]
        # write to disk in one go
        ovar[:] = outvar

fi.close()
fo.close()

# -------------------------------------------------------------------------
# Finish

if verbose:
    tstop = ptime.time()
    print('Finished in [s]: {:.2f}'.format(tstop - tstart))
