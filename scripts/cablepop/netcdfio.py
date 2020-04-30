#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
"""
    NetCDF4 helper tools
"""
import numpy as np


__all__ = ['create_dimensions', 'create_variables', 'get_variable_definition',
           'set_global_attributes', 'set_output_filename']


def create_dimensions(fi, fo, change={}, exclude=[]):
    if exclude: # assert list
        if isinstance(exclude, str): exclude = [exclude]
    for d in fi.dimensions.values():
        if d.name not in exclude:
            if d.name in change.keys():
                nd = change[d.name]
            elif d.isunlimited():
                nd = None
            else:
                nd = len(d)
            fo.createDimension(d.name, nd)


def create_variables(fi, fo, time=None, izip=False, **kwargs):
    for ivar in fi.variables.values():
        if time is None:
            itime = True
        else:
            if time:
                itime = 'time' in ivar.dimensions
            else:
                itime = 'time' not in ivar.dimensions
        if itime:
            invardef = get_variable_definition(ivar, **kwargs)
            if izip: invardef.update({'zlib':True})
            ovar = fo.createVariable(invardef.pop("name"), invardef.pop("dtype"), **invardef)
            for k in ivar.ncattrs():
                iattr = ivar.getncattr(k)
                if (k != 'missing_value') and (k != '_FillValue'):
                    ovar.setncattr(k, iattr)


def get_variable_definition(ncvar, fill=None, chunksizes=True, removedim=[], renamedim={}):
    if removedim: # assert list
        if isinstance(removedim, str): removedim = [removedim]
    out  = ncvar.filters() if ncvar.filters() else {}
    dims = list(ncvar.dimensions)
    for dd in renamedim.keys():
        if dd in dims:
            dims[dims.index(dd)] = renamedim[dd]
    if chunksizes:
        chunks = ncvar.chunking() if not isinstance(ncvar.chunking(), str) else None
    else:
        chunks = None
    for dd in removedim:
        if dd in ncvar.dimensions:
            ip = dims.index(dd)
            dims.remove(dd)
            if chunks:
                _ = chunks.pop(ip)        
    if "missing_value" in dir(ncvar):
        ifill = ncvar.missing_value
    elif "_FillValue" in dir(ncvar):
        ifill = ncvar._FillValue
    else:
        if fill:
            if ncvar.dtype in (np.dtype(np.int8), np.dtype(np.int16), np.dtype(np.int32), np.dtype(np.int64)):
                ifill = None
            else:
                ifill = fill
        else:
            ifill = None
    out.update({
        "name"       : ncvar.name,
        "dtype"      : ncvar.dtype,
        "dimensions" : dims,
        "fill_value" : ifill,
    })
    if chunksizes:
        out.update({"chunksizes" : chunks})
    return out


def set_global_attributes(fi, fo, add={}):
    for k in fi.ncattrs():
        iattr = fi.getncattr(k)
        if k in add.keys():
            iattr += '\n'+add[k]
        fo.setncattr(k, iattr)
    # if global attribute does not exist yet, add it
    for k in add.keys():
        if k not in fi.ncattrs():
            fo.setncattr(k, add[k])


def set_output_filename(ifile, ext):
    sifile = ifile.split('.')
    sifile[-2] = sifile[-2]+ext
    ofile = '.'.join(sifile)
    return ofile


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
