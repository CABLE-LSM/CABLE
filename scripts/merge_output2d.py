#!/usr/bin/env python
"""


Example
-------
  python merge_output2d.py -o cru_out_cable_1700_1900.nc \
      run*/outputs/cru_out_cable_1700_1900.nc


History
-------
Written  Matthias Cuntz, May 2023
             - from sum_patchcasa2d.py

"""
import sys
import argparse
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
    description=('Merge Cable 2D output files from pseudo-parallelised runs.'))
parser.add_argument('-o', '--outfile', action='store',
                    default=ofile, dest='ofile', metavar='output_netcdf',
                    help=('output netcdf file name (default:'
                          ' first_input-merged.nc).'))
parser.add_argument('-v', '--verbose', action='store_true', default=verbose,
                    dest='verbose',
                    help='Feedback during merge (default: no feedback).')
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

#
# Check first input file and make output file
#
ifile = ifiles[0]
fi = nc.Dataset(ifile, 'r')
if verbose:
    print('Check first input file:  ', ifile)
ncvars = list(fi.variables.keys())
ntime = fi.dimensions['time'].size

# Output file
if ofile is None:  # Default output filename
    ofile = cp.set_output_filename(ifile, '-merge')
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
cp.set_global_attributes(
    fi, fo, add={'history': (ptime.asctime() + ': ' + ' '.join(sys.argv))})

# Copy dimensions
cp.create_dimensions(fi, fo)

# create static variables (independent of time)
cp.create_variables(fi, fo, time=False, izip=izip, fill=True,
                    chunksizes=False)

# create dynamic variables (dependent of time)
cp.create_variables(fi, fo, time=True, izip=izip, fill=True,
                    chunksizes=False)

fi.close()

if verbose:
    print('Make masks')
nfiles = len(ifiles)
fis = []
iix = []
iiy = []
nii = []
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
    # Use patchfrac to determine mask
    if 'patchfrac' not in fi.variables:
        for fi in fis:
            fi.close()
        fo.close()
        raise ValueError(f'No variable patchfrac in input file {ifile}')
    # Get first patchfrac
    patchfrac = fi.variables['patchfrac']
    if 'time' in patchfrac.dimensions:
        pfrac = patchfrac[0, 0, ...]
    else:
        pfrac = patchfrac[0, ...]
    mask = ~pfrac.mask
    if not np.iterable(mask):
        for fi in fis:
            fi.close()
        fo.close()
        raise ValueError(f'patchfrac has no valid mask in file {ifile}')
    ii = np.where(mask)
    iiy.append(ii[0])
    iix.append(ii[1])
    nii.append(ii[1].size)

# Iterate over input files
if verbose:
    print('Copy input to output')
for ncvar in ncvars:
    if verbose:
        tstartvar = ptime.time()
        print(f'    {ncvar}')
    ivar0 = fis[0].variables[ncvar]
    ovar = fo.variables[ncvar]
    if ( (ncvar == 'longitude') or
         (ncvar == 'latitude') or
         ('x' not in ivar0.dimensions) or
         ('y' not in ivar0.dimensions) ):
        # should not be masked and all the same: check
        for f, fi in enumerate(fis):
            ivar = fi.variables[ncvar]
            if not np.all(ivar0[:] == ivar[:]):
                for fi in fis:
                    fi.close()
                fo.close()
                print('ivar0', ivar0[:])
                print('ivar', ivar[:])
                raise ValueError(f'variable {ncvar} not equal in file'
                                 f' {ifiles[0]} and file {ifiles[f]}')
        ovar[:] = ivar0[:]
    else:
        if ('time' not in ivar0.dimensions):
            # write in temporary variable
            invar0 = ivar0[:]
            outvar = np.full_like(
                invar0, pj.ncio.get_fill_value_for_dtype(invar0.dtype))
            for f, fi in enumerate(fis):
                ivar = fi.variables[ncvar]
                invar = ivar[:]
                if invar.ndim == 2:
                    outvar[iiy[f], iix[f]] = (
                        invar[iiy[f], iix[f]] )
                else:
                    outvar[..., iiy[f], iix[f]] = (
                        invar[..., iiy[f], iix[f]] )
            # write whole field at once into file
            ovar[:] = outvar
        else:
            invar0 = ivar0[0, ...]
            for tt in range(ntime):
                # tstartread = ptime.time()
                # write in temporary variable
                outvar = np.full_like(
                    invar0, pj.ncio.get_fill_value_for_dtype(invar0.dtype))
                for f, fi in enumerate(fis):
                    ivar = fi.variables[ncvar]
                    invar = ivar[tt, ...]
                    if invar.ndim == 2:
                        outvar[iiy[f], iix[f]] = (
                            invar[iiy[f], iix[f]] )
                    else:
                        outvar[..., iiy[f], iix[f]] = (
                            invar[..., iiy[f], iix[f]] )
                # tstopread = ptime.time()
                # tstartwrite = tstopread
                # write whole field at once into file
                ovar[tt, ...] = outvar
                # tstopwrite = ptime.time()
                # print(f'  r/w [s]: {tstopread - tstartread}'
                #       f' {tstopwrite - tstartwrite}')
    tstopvar = ptime.time()
    print(f'        [s]: {tstopvar - tstartvar}')

for fi in fis:
    fi.close()
fo.close()

# -------------------------------------------------------------------------
# Finish

if verbose:
    tstop = ptime.time()
    print('Finished in [s]: {:.2f}'.format(tstop - tstart))
