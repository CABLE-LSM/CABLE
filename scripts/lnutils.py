# -*- coding: utf-8 -*-
"""
Created on Tue Sep  2 16:46:29 2014

@author: nie06a
"""
import numpy as np
import sys

def is_leapyear( year ):
    if ( np.mod(year,4) == 0 and np.mod(year,100) != 0 ) or \
         np.mod(year,400) == 0:
        return True
    else:
        return False

def leap_day( year, month=0 ):
    """
    1 if it is a leap_year and month is set 0 (to determine e.g. length-of-year
    ) or 2 (feb), 0 else
    """
    if is_leapyear( year ) and  ( month == 2 or month == 0 ):
        return 1
    else:
        return 0

def length_of_month( year, month ):
    lom = [31,28,31,30,31,30,31,31,30,31,30,31]
    if is_leapyear( year ) and month == 2:
        return 29
    else:
        return lom[month-1]

def handle_err( ex, msg='An unspecified error occured!', kill=-1 ):
    if ex != 0:
        print(msg + " "+ str(ex))
        if kill != 0:
            sys.exit(kill)

#def latlon2ixjy (LAT,LON,nx,ny,mtype='array'):
def latlon2ixjy(LAT,LON,latmin,latmax,lonmin,lonmax,nx,ny,mtype='array'):
    ex = 0
    chk = np.array(LAT,ndmin=1)
    if len(chk[chk<-90.])>0 or len(chk[chk>90.])>0:
        print("Excess LAT in latlon2ixjy")
        ex = 1
    chk = np.array(LON,ndmin=1)
    if len(chk[chk<-180.])>0 or len(chk[chk>180.])>0:
        print("Excess LON in latlon2ixjy")
        ex = 1
    if ex == 1:
        sys.exit(-1)

    valx = np.array((LON + 180.) / 360. * nx)
    # valy = np.array((LAT +  90.) / 180. * ny)
    valy = np.array((LAT -  latmin) / (latmax-latmin) * ny)
    XY = np.copy(valx)
    XY = np.append([XY], [valy], axis=0).reshape(2,-1)
    if mtype == 'array':
        XY = np.fix(XY).astype(int)
        return XY
    elif mtype == 'mask':
        MASK = np.zeros((nx,ny))
        MASK[np.fix(XY[0]).astype(np.int),np.fix(XY[1]).astype(np.int)] = 1
        return MASK.T
    else:
        print("Wrong mtype in latlon2ixjy. Either 'array' or 'mask'!")
        sys.exit(-1)

def cell_area ( LAT, LATRES, LONRES, unit='m2' ):
    """
    Compute Area of a grid cell for a lat-lon grid
    Following http://eos-earthdata.sr.unh.edu/data/dataGuides/global_model_dg.pdf
    LN 02/2014
    Input:
    LAT, LATRES, LONRES : Latitude, Lat/Lon resolution; all [decimal deg]
    unit                : desired output unit=< m2|ha|km2 >; default="m2"
    Output:
    CELL_AREA           : area in units of input "unit"
    """
    from math import pi

    R  = 6371221.3   # Earth radius [m]

    chk = np.array(LAT,ndmin=1)
    if len(chk[chk<-90.])>0 or len(chk[chk>90.])>0:
        print("Wrong LAT in cell_area")
        sys.exit(-1)

    rads = np.array( 90. - LAT + .5 * LATRES ) * pi / 180.
    coss = np.cos(rads) - np.cos(rads + ( LONRES * pi / 180. ))
    CELL_AREA = ( pi * R * R * coss / 360.)

    if unit == 'km2':
        CELL_AREA *= 1./1000000.
    elif unit == 'ha':
        CELL_AREA *= 1./10000.
    elif unit != 'm2':
        print("invalid unit '"+unit+"' in call to cell_area")
        sys.exit(-1)

    return CELL_AREA

def find_closest_match(arr,val,chk_bound=1, is_sorted=1):
    if len(arr.shape) != 1:
        print("Array must have dim = 1")
        exit(-1)
    if len(arr) < 2:
        print("Array is not an array" )
        exit(-1)
    if chk_bound:
        if val < np.min(arr) or val > np.max(arr):
            print("Value is out of bounds")
            print("min "+str(np.min(arr))+" val "+str(val)+" max "+str(np.max(arr)))
            exit(-1)

    sdif = 9.e+12
    idx = 0
    for i in range(len(arr)):
        if np.abs(arr[i] - val) < sdif:
            sdif = np.abs(arr[i] - val)
            idx = i
        elif is_sorted:
            return idx,arr[idx]

    return idx,arr[idx]

def running_mean(arr, n, typ="central", edge='valid'):
    """
    arr: 1-dim input array
    n  : averaging period
    typ: 'central' or 'past'
    edge: how to act on edges
         "valid": leave out, i.e. NaNs
         "edge" : include virtual zeros
         "flow" : in/decrease n on edges
         "fill" : fill with "edge" values
    """
    if n > len(arr):
        exit("Arr to short for period of "+str(n)+" in 'running_mean'")


    if edge == "valid":
        darr = np.zeros((len(arr)))
        darr[:] = np.nan
    elif edge == "edge":
        darr = np.zeros((len(arr)))
    elif edge == "fill":
        darr = np.copy(arr)
        darr[:] = arr[0]
    else:
        print("wrong type 'edge' in call to 'running_mean'")
        exit(-1)

    ns = n - np.int((n+1) / 2.)
    ne = n - ns

    nn = len(arr)

    for i in range(nn):
        if typ == 'central':
            if edge == 'valid':
                if i >= ns and i <= nn-ne:
                    darr[i] = np.average(arr[i-ns:i+ne])
            elif edge == 'edge':
                if i< ns:
                    darr[i] = np.average(arr[:i+ne])
                elif i>nn-ne:
                    darr[i] = np.average(arr[i-ns:])
                else:
                    darr[i] = np.average(arr[i-ns:i+ne])
                print("edge",i,darr[i])
            else:
                exit(" NA")
        elif typ == 'past':
            exit("'past' not yet implemented in lnutils.py:running_mean")
            if edge == 'valid':
                if i >= ns and i <= nn-ne:
                    darr[i] = np.average(arr[i-ns:i+ne])
            elif edge == 'edge':
                if i< ns:
                    darr[i] = np.average(arr[:i+ne])
                elif i>nn-ne:
                    darr[i] = np.average(arr[i-ns:])
                else:
                    darr[i] = np.average(arr[i-ns:i+ne])
                print("edge",i,darr[i])
            else:
                exit(" NA")


    return darr

def get_rel_col( AR1, AR2, facs, cols ):
    zidx = len(facs) + 1
    Xout = np.array([cols[zidx]]*len(AR1))
    Xout[:] = cols[zidx]
    for i in range(len(facs)):
        f = facs[i]
        Xout[np.where(AR1>f*AR2)] = cols[zidx + i]
        f = 1./facs[i]
        Xout[np.where(AR1<f*AR2)] = cols[zidx - i]
    return Xout

def shiftColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    import matplotlib
    import matplotlib.pyplot as plt

    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower ofset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax/(vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highets point in the colormap's range.
          Defaults to 1.0 (no upper ofset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False),
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap

def rank_array( arr, eliminate0='no' ):
    """
    Rank the values of an array equidistantly between 0 and 1
    Input
    arr:   1-dimensional integer or real array
    eliminate0: if 'yes' zeros will be taken outof ranking
           NaNs are taken outanyways
    Output
    ranks: 1-dimensional array with ranks
    """
    from math import isnan

    tmp        = arr.argsort()
    ranks      = np.empty(len(arr),float)
    ranks[tmp] = np.arange(len(arr))

    for i in range(len(arr)):
        if (arr[i] == 0 and eliminate0 != 'no') \
                or isnan(arr[i]):
            ranks[i] = np.nan

    return ranks


