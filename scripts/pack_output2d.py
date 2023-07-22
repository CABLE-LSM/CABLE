#!/usr/bin/env python3
"""
usage: pack_output2d.py [-h] [-o output_netcdf] [-v] [-z] [file]

Pack Cable 2D output files by removing masked ocean and land points.

positional arguments:
  file                  input netcdf file.

options:
  -h, --help            show this help message and exit
  -o output_netcdf, --outfile output_netcdf
                        output netcdf file name (default: input-packed.nc).
  -v, --verbose         Feedback during merge (default: no feedback).
  -z, --zip             Use netCDF4 variable compression
                        (default: same format as input file).


Example
-------
  python pack_output2d.py -o cru_out_cable_1700_1900-packed.nc \
      run11/outputs/cru_out_cable_1700_1900.nc


History
-------
Written  Matthias Cuntz, May 2023

"""
import sys
import argparse
import numpy as np
import netCDF4 as nc
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
    description=('Pack Cable 2D output files by removing masked ocean and'
                 ' land points.'))
parser.add_argument('-o', '--outfile', action='store',
                    default=ofile, dest='ofile', metavar='output_netcdf',
                    help=('output netcdf file name (default:'
                          ' input-packed.nc).'))
parser.add_argument('-v', '--verbose', action='store_true', default=verbose,
                    dest='verbose',
                    help='Feedback during merge (default: no feedback).')
parser.add_argument('-z', '--zip', action='store_true', default=izip,
                    dest='izip',
                    help=('Use netCDF4 variable compression (default:'
                          ' same format as input file).'))
parser.add_argument('ifile', nargs='?', default=None, metavar='file',
                    help='input netcdf file.')
args    = parser.parse_args()
ofile   = args.ofile
verbose = args.verbose
izip    = args.izip
ifile   = args.ifile
del parser, args

if not ifile:
    raise IOError('Input file must be given.')

if verbose:
    tstart = ptime.time()

# -------------------------------------------------------------------------
# Copy data
#

#
# Check input file and create output file
#
fi = nc.Dataset(ifile, 'r')
if verbose:
    print('Check first input file:  ', ifile)
ncvars = list(fi.variables.keys())
ntime = fi.dimensions['time'].size

# Use patchfrac to determine mask
if 'patchfrac' not in fi.variables:
    fi.close()
    raise ValueError(f'No variable patchfrac in input file {ifile}')
# Get first patchfrac
patchfrac = fi.variables['patchfrac']
if 'time' in patchfrac.dimensions:
    pfrac = patchfrac[0, 0, ...]
else:
    pfrac = patchfrac[0, ...]
mask = ~pfrac.mask
if not np.iterable(mask):
    fi.close()
    raise ValueError(f'patchfrac has no valid mask in file {ifile}')
ii = np.where(mask)
iiy = ii[0]
iix = ii[1]
nii = ii[1].size

# Output file
if ofile is None:  # Default output filename
    ofile = pj.ncio.set_output_filename(ifile, '-packed')
if verbose:
    print('Create output file: ', ofile)
if izip:
    fo = nc.Dataset(ofile, 'w', format='NETCDF4')
else:
    if 'file_format' in dir(fi):
        fo = nc.Dataset(ofile, 'w', format=fi.file_format)
    else:
        fo = nc.Dataset(ofile, 'w', format='NETCDF3_64BIT_OFFSET')

# Copy global attributes, adding script
pj.ncio.copy_global_attributes(
    fi, fo, add={'history': (ptime.asctime() + ': ' + ' '.join(sys.argv))})

# Copy dimensions
pj.ncio.copy_dimensions(fi, fo, removedim=['x', 'y'], adddim={'land': nii})

# create static variables (independent of time)
pj.ncio.create_variables(fi, fo, time=False, izip=izip, fill=True,
                         chunksizes=False, removedim=['x'],
                         replacedim={'y': 'land'},
                         removevar=['x', 'y'])

# create dynamic variables (dependent of time)
pj.ncio.create_variables(fi, fo, time=True, izip=izip, fill=True,
                         chunksizes=False, removedim=['x'],
                         replacedim={'y': 'land'})

# copy static and dynamic variables
if verbose:
    print('Copy variables')
nvar = len(fi.variables.keys())
if nvar < 20:
    n10 = 1
else:
    n10 = nvar // 10
n = 0
for ivar in fi.variables.values():
    if (ivar.name == 'x') or (ivar.name == 'y'):
        continue
    if verbose:
        print(f'    {ivar.name}')
    n += 1
    ovar = fo.variables[ivar.name]
    # read whole field, otherwise execution time increases sharpely
    invar = ivar[:]
    if ( ('x' not in ivar.dimensions) or
         ('y' not in ivar.dimensions) ):
        outvar = invar
    else:
        outvar = pj.pack(invar, mask)
    ovar[:] = outvar

fi.close()
fo.close()

# -------------------------------------------------------------------------
# Finish

if verbose:
    tstop = ptime.time()
    print('Finished in [s]: {:.2f}'.format(tstop - tstart))
