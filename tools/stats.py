#! /usr/bin/env python

"""
Some useful functions
  percentiles(xbins,xarray,yarray,per): obtains percentiles of yarray in xbins

NOTE: this module requires the numpy and scipy libraries to be
      available for import!

"""
import sys
import numpy as np

def percentiles(val,data,weights=None):
    if (val <0 or val >1):
        sys.exit('STOP percentiles: 0<val<1')

    if (weights is None):
        ws = np.zeros(shape=(len(data))) ; ws.fill(1.)
    else:
        ws = weights

    data = np.array(data) ; ws = np.array(ws)
    ind_sorted = np.argsort(data)  # Median calculation from wquantiles
    sorted_data = data[ind_sorted] ; sorted_weights = ws[ind_sorted]
    
    num = np.cumsum(sorted_weights) - 0.5*sorted_weights 
    den = np.sum(sorted_weights) 
    if (den!=0): 
        pn = num/den   
        percentiles = np.interp(val, pn, sorted_data)  
    else:
        sys.exit('STOP percentiles: problem with weights')
    return percentiles

def perc_2arrays(xbins,xarray,yarray,weights,nmin,val):
    """ Returns percentiles of yarray over xbins"""
    xlen = len(xbins)-1
    perc_2arrays = np.zeros(shape=(xlen)) ; perc_2arrays.fill(-999.)

    if len(xarray) != len(yarray):
        sys.exit('ERROR @ perc_2arrays: The lenght of the input arrays should be equal.')

    for i in range(xlen):
        ind = np.where((xarray >= xbins[i]) & (xarray < xbins[i+1]))
        # We require at least nmin points per bin
        if (np.shape(ind)[1] > nmin): 
            data = yarray[ind] ; ws = weights[ind]
            perc_2arrays[i] = percentiles(val,data,weights=ws)

    return perc_2arrays
