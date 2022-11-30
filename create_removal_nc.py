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
import scipy.spatial as spatial

class KDTreeIndex():

    """ A KD-tree implementation for fast point lookup on a 2D grid

    Keyword arguments:
    dataset -- a xarray DataArray containing lat/lon coordinates
               (named 'lat' and 'lon' respectively)

    """

    def transform_coordinates(self, coords):
        """ Transform coordinates from geodetic to cartesian

        Keyword arguments:
        coords - a set of lan/lon coordinates (e.g. a tuple or
                 an array of tuples)
        """
        # WGS 84 reference coordinate system parameters
        A = 6378.137 # major axis [km]
        E2 = 6.69437999014e-3 # eccentricity squared

        coords = np.asarray(coords).astype(np.float64)

        # is coords a tuple? Convert it to an one-element array of tuples
        if coords.ndim == 1:
            coords = np.array([coords])

        # convert to radiants
        lat_rad = np.radians(coords[:,0])
        lon_rad = np.radians(coords[:,1])

        # convert to cartesian coordinates
        r_n = A / (np.sqrt(1 - E2 * (np.sin(lat_rad) ** 2)))
        x = r_n * np.cos(lat_rad) * np.cos(lon_rad)
        y = r_n * np.cos(lat_rad) * np.sin(lon_rad)
        z = r_n * (1 - E2) * np.sin(lat_rad)

        return np.column_stack((x, y, z))

    def __init__(self, dataset):
        # store original dataset shape
        print(dataset.nav_lat[1:3,1])
        self.shape = dataset.shape
        print(dataset)
        # reshape and stack coordinates
        coords = np.column_stack((dataset.nav_lat.values.ravel(),
                                  dataset.nav_lon.values.ravel()))

        # construct KD-tree
        self.tree = spatial.cKDTree(self.transform_coordinates(coords))

    def query(self, point):
        """ Query the kd-tree for nearest neighbour.

        Keyword arguments:
        point -- a (lat, lon) tuple or array of tuples
        """
        _, index = self.tree.query(self.transform_coordinates(point))

        # regrid to 2D grid
        index = np.unravel_index(index, self.shape)

        # return DataArray indexers
        return xr.DataArray(index[0], dims='pixel'), \
               xr.DataArray(index[1], dims='pixel')

    def query_ball_point(self, point, radius):
        """ Query the kd-tree for all point within distance
        radius of point(s) x

        Keyword arguments:
        point -- a (lat, lon) tuple or array of tuples
        radius -- the search radius (km)
        """

        index = self.tree.query_ball_point(self.transform_coordinates(point),
                                           radius)

        # regrid to 2D grid
        index = np.unravel_index(index[0], self.shape)

        # return DataArray indexers
        return xr.DataArray(index[0], dims='pixel'), \
               xr.DataArray(index[1], dims='pixel')

class Create_mask_removal:
    """
    Create mask for removal of CO2 or liming experiments
    1. read in bathymety file for input
    2. drop extra variables
    3. change name
    4. change values 1 for ocean, 0 for land

    5a. select region or
    5b. latitude range

    """
    def __init__(self,fname):
        self.fname=fname
        self.data_set=None
        self.outputdataset={}
    def open_inputfile(self):
        self.data_set=xr.open_dataset(self.fname)

    def remove_extra_data(self,vars_2_remove=['isf_draft','Bathymetry_isf']):
        self.data_set=self.data_set.drop(vars_2_remove)

    def change_name(self,name_dict={"Bathymetry":"co2r"}):
        self.data_set=self.data_set.rename(name_dict)

    def set_sea_to_1(self,value=1):
        self.data_set=self.data_set.where(self.data_set.co2r<1,1)
        self.data_set.co2r.attrs['units']='1'
        print(self.data_set.co2r)
    def add_months(self,value=12):
        copy=self.data_set.expand_dims({'time':value})
        self.data_set=copy
        for outset in self.outputdataset:
            copy=self.outputdataset[outset].expand_dims({'time':value})
            self.outputdataset[outset]=copy
        # print(copy)
        # test =self.data_set.reindex({'time':value})
        # test.ffill("time")
        # print (test)
        # dasf
        # for t in range(value):
        #     copy.co2r.data[t,:,:]=copy.co2r.data[0,:,:]
        #     print(t)
        #self.data_set.where(self.data_set.co2r<1,1)

    def select_lats(self,value=20):
        new=self.data_set.where(abs(self.data_set.nav_lat)<value,0)
        self.outputdataset[value]=new
        return new

    def select_region(self,lats,lons):
        if lats[1]<lats[0]:
            lats=np.swap(lats)
        if lons[1]<lons[0]:
            lons=np.swap(lons)
        masklon= (self.data_set.nav_lon>lons[0]) & (self.data_set.nav_lon<lons[1])
        masklat= (self.data_set.nav_lat>lats[0]) & (self.data_set.nav_lat<lats[1])
        new=self.data_set.where(masklon & masklat,0)
        self.outputdataset['region']=(new)
        return new

    def shipping_lane_mask(self,lons=(-100,20),lats=(15,60)):
        fullset=xr.open_dataset('../data/co2emi-2014-orca-ym.nc')
        print(fullset)
        ships=fullset.isel(sector=7).isel(time=0)
        print(ships)
        name_dict={'CO2_em_anthro':'lanemask'}
        maxdata=ships.CO2_em_anthro.mean()
        shipsnorm=ships.rename(name_dict)
        #shipsnorm['lanemask']=shipsnorm['lanemask']/maxdata
        mask=shipsnorm['lanemask']<1e-12
        shipsnorm['lanemask']=xr.where(mask,0,shipsnorm['lanemask'])
        mask=shipsnorm['lanemask']>1e-9
        shipsnorm['lanemask']=xr.where(mask,1e-9,shipsnorm['lanemask'])
        maxdata=shipsnorm.lanemask.max()
        shipsnorm['lanemask']=shipsnorm['lanemask']/maxdata
        masklon= (shipsnorm.nav_lon>lons[0]) & (shipsnorm.nav_lon<lons[1])
        masklat= (shipsnorm.nav_lat>lats[0]) & (shipsnorm.nav_lat<lats[1])
        shipsnorm=shipsnorm.where(masklon & masklat,0)
        print(shipsnorm)
        return shipsnorm
    def find_coast(self):
        #mask=self.data_set.copy()

        pass



    def calculate_distance(data_set,coord):
        pass

    def write_data_set(self,which='all'):
        if which=='all':
            for outdata_set in self.outputdataset:
                self.outputdataset[outdata_set].to_netcdf('co2removal_'+str(outdata_set)+'.nc')
        else:
            print ('not implemented yet')

def main():
    # read file to base the mask
    test=Create_mask_removal('../data/bathy_meter.nc')
    test.open_inputfile()
    newmask=test.shipping_lane_mask()
    #clean up
    test.remove_extra_data()
    # variable name change
    test.change_name()
    # create 1/0 mask
    test.set_sea_to_1()

    # show us! all data
    test.data_set.co2r.plot()

    new2=test.select_region([20,40], [-80,-10])
    new1=test.select_lats(10)
    plt.figure()
    new1.co2r.plot()
    plt.figure()
    test.outputdataset[10].co2r.plot()
    plt.figure()
    new2.co2r.plot()
    plt.figure()
    test.add_months(120)
    print(test)
    for i in range(12):
        plt.figure()
        test.outputdataset[10].co2r.isel(time=i).plot()
    plt.figure()
    test.data_set.co2r.isel(time=10).plot()
    test.write_data_set()
    # output
    plt.figure()
    newmask.lanemask.plot()
    plt.show()

# =============================================================================
#     print(test.data_set)
#     puu=KDTreeIndex(test.data_set.co2r)
#     rome = (-30, 12.4964)
#     rome_index = puu.query(rome)
#     print(test.data_set.co2r[rome_index])
#     d200=puu.query_ball_point(rome, 200)
#     d100=puu.query_ball_point(rome, 100)
#     print('heå',d200,np.shape(d200))
#     for i,j in np.ndindex(np.shape(d200)):
#         print(i,j)
#         #print(d200[i,j])
#         #if (i,j) in d100:
#         #    print('hep',(i,j))
#             #d200.pop(i)
#     print(test.data_set.co2r[d200])
#     print(test.data_set.co2r[d100])
# =============================================================================


if __name__ == "__main__":
    main()