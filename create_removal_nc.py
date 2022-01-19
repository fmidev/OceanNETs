#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 16 16:10:54 2021

@author: bergmant
"""

import xarray as xr
#matplotlib inline
import matplotlib.pyplot as plt
import numpy as np
def open_inputfile(fname):
    ds=xr.open_dataset(fname)
    return ds

def remove_extra_data(ds,vars_2_remove=['isf_draft','Bathymetry_isf']):
    ds=ds.drop(vars_2_remove)
    return ds
def change_name(ds,name_dict={"Bathymetry":"co2r"}):
    ds=ds.copy()
    ds=ds.rename(name_dict)
    return ds
def modify_values(ds,value=1):
    ds=ds.where(ds.co2r<1,1)
    return ds
def select_lats(ds,value=20):
    new=ds.where(abs(ds.nav_lat)<value,0)    
    return new
def select_region(ds,lats,lons):
    if lats[1]<lats[0]:
        lats=np.swap(lats)
    if lons[1]<lons[0]:
        lons=np.swap(lons)
    masklon= (ds.nav_lon>lons[0]) & (ds.nav_lon<lons[1])
    masklat= (ds.nav_lat>lats[0]) & (ds.nav_lat<lats[1])
    new=ds.where(masklon & masklat,0)
    return new
def find_coast(ds):
    pass
def calculate_distance(ds,coord):
    pass

def main():
    ds=open_inputfile('bathy_meter.nc')
    ds=remove_extra_data(ds)
    ds=change_name(ds)
    ds=modify_values(ds)
    ds.co2r.plot()
    new2=select_region(ds, [20,40], [-80,-10])
    new1=select_lats(ds,10)    
    plt.figure()
    new1.co2r.plot()
    plt.figure()
    new2.co2r.plot()
    new1.to_netcdf('co2removal.nc')
    plt.show()
if __name__ == "__main__":
    main()