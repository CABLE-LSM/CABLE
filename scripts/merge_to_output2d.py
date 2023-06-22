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

Remember:
https://chase-seibert.github.io/blog/2013/08/03/diagnosing-memory-leaks-python.html
https://gist.github.com/schlamar/2311116
https://stackoverflow.com/questions/15455048/releasing-memory-in-python

"""
import argparse
import sys
import numpy as np
import netCDF4 as nc
import cablepop as cp
import pyjams as pj
import time as ptime
import psutil
import gc

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

gb = 1073741824.  # (1024 * 1024 * 1024)

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
    print('Create output file:', ofile)
if izip:
    oformat = 'NETCDF4'
else:
    if 'file_format' in dir(fi):
        oformat = fi.file_format
    else:
        oformat = 'NETCDF3_64BIT_OFFSET'
fo = nc.Dataset(ofile, 'w', format=oformat)

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

# create x and y for cdo, etc.
if 'x' not in fi.variables:
    if verbose:
        print('Create x')
    nvar = {'name': 'x',
            'dtype': ilons.dtype,
            'dimensions': ('x'),
            'units': 'degrees_east'}
    ovar = pj.ncio.create_new_variable(nvar, fo)
    ovar[:] = olon

if 'y' not in fi.variables:
    if verbose:
        print('Create y')
    nvar = {'name': 'y',
            'dtype': ilats.dtype,
            'dimensions': ('y'),
            'units': 'degrees_north'}
    ovar = pj.ncio.create_new_variable(nvar, fo)
    ovar[:] = olat

# write time for correct output shape
if verbose:
    print('Write time')
ivar = fi.variables['time']
ovar = fo.variables['time']
ovar[:] = ivar[:]

fi.close()
fo.close()

if verbose:
    print('Get all indexes')
nfiles = len(ifiles)
iidl = []
oidx = []
oidy = []
for nfile, ifile in enumerate(ifiles):
    fi = nc.Dataset(ifile, 'r')
    # Check time
    ntime1 = fi.dimensions['time'].size
    if ntime1 != ntime:
        fi.close()
        raise ValueError(f'Time not the same in {ifiles[0]} and in {ifile}')
    # get latitude/longitude indices
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
    fi.close()

#
# Copy variables from in to out expanding the land dimension to y, x
#

# copy static and dynamic variables
if verbose:
    print('Copy input to output')
n = 0
for ncvar in ncvars:
    if ncvar == 'time':
        continue
    if verbose:
        tstartvar = ptime.time()
        print(f'    {ncvar}')
    n += 1
    ifile0 = ifiles[0]
    fi0 = nc.Dataset(ifile0, 'r')
    fo = nc.Dataset(ofile, 'a', format=oformat)
    ivar0 = fi0.variables[ncvar]
    ivar0_dtype = ivar0.dtype
    ovar = fo.variables[ncvar]
    if ncvar == 'longitude':
        fi0.close()
        ovar[:] = olon2d
    elif ncvar == 'latitude':
        fi0.close()
        ovar[:] = olat2d
    elif (('land' not in ivar0.dimensions) and
          ('ntile' not in ivar0.dimensions)):
        # should not be masked and all the same: check
        ivar00 = ivar0[:]
        fi0.close()
        for ifile in ifiles:
            fi = nc.Dataset(ifile, 'r')
            ivar = fi.variables[ncvar][:]
            if not np.all(ivar00 == ivar):
                print('ivar0', ivar00)
                print('ivar', ivar)
                fi.close()
                fo.close()
                raise ValueError(f'variable {ncvar} not equal in file'
                                 f' {ifile0} and file {ifile}')
            fi.close()
            del ivar
        ovar[:] = ivar00
        del ivar00
    elif ('time' not in ivar0.dimensions):
        fi0.close()
        outvar = np.full(ovar.shape,
                         pj.ncio.get_fill_value_for_dtype(ivar0_dtype))
        for f, ifile in enumerate(ifiles):
            fi = nc.Dataset(ifile, 'r')
            ivar = fi.variables[ncvar]
            # read whole field
            invar = ivar[:]
            # fill in memory
            if len(ivar.shape) == 1:
                outvar[oidy[f], oidx[f]] = invar[iidl[f]]
            else:
                outvar[..., oidy[f], oidx[f]] = invar[..., iidl[f]]
            fi.close()
            del invar, ivar
        # write to disk in one go
        ovar[:] = outvar
        del outvar
    else:  # has time and land/ntile
        fi0.close()
        print(f'        {ovar.shape}')
        nt = np.ceil(np.prod(ovar.shape) * 8 / gb / 2).astype(int)
        tindexes = np.linspace(0, ntime, nt+1, dtype=int)
        for nn in range(nt):
            oshape = list(ovar.shape)
            i1 = tindexes[nn]
            i2 = tindexes[nn + 1]
            oshape[0] = i2 - i1
            if nt > 1:
                print(f'        {oshape} {i1} {i2}')
            outvar = np.full(oshape,
                             pj.ncio.get_fill_value_for_dtype(ivar0_dtype))
            mem = psutil.Process().memory_info()
            if verbose:
                tstartread = ptime.time()
            for f, ifile in enumerate(ifiles):
                fi = nc.Dataset(ifile, 'r')
                ivar = fi.variables[ncvar]
                    # read time steps
                invar = ivar[i1:i2, ...]
                # fill in memory
                outvar[..., oidy[f], oidx[f]] = invar[..., iidl[f]]
                fi.close()
                del invar, ivar
            # write to disk in one go
            if verbose:
                tstopread = ptime.time()
                print('        Read  {:.2f} s'.format(tstopread - tstartread))
                tstartwrite = tstopread
            ovar[i1:i2, ...] = outvar
            if verbose:
                tstopwrite = ptime.time()
                print('        Wrote {:.2f} s'.format(tstopwrite - tstartwrite))
                mem = psutil.Process().memory_info()
                print(f'        Memory physical [GB]: {mem.rss / gb:.2f} virtual: {mem.vms / gb:.2f}')
            del outvar
    fo.close()
    del ivar0, ovar
    gc.collect()
    if verbose:
        tstopvar = ptime.time()
        print('        Total {:.2f} s'.format(tstopvar - tstartvar))

# -------------------------------------------------------------------------
# Finish

if verbose:
    tstop = ptime.time()
    print('Finished in [s]: {:.2f}'.format(tstop - tstart))
