#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 09:38:10 2022

@author: bergmant
"""
import numpy as np
import xarray as xr
from itertools import product
import pandas as pd
import matplotlib.pyplot as plt

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import regionmask
import csv
from datetime import datetime, timedelta
import calendar
from mask_aux import create_eez_mask, gridarea, masking_boxes, read_deployment_data, year2mon,  select_areas
#import ESMpy
# import xesmf as xe


start_years=['2030', '2035', '2040']
volumes=['low', 'high']
scenarios=list()

for start_year, volume in product(start_years,volumes):
    scenarios.append(start_year+'-'+volume)

dep_data=read_deployment_data()
dep_data_monthly=dict()
dep_data_monthly_mol_s=dict()
dep_data_monthly_mol_s_m2=dict()


for scenario in scenarios:



    dep_data_monthly[scenario], dep_data_monthly_mol_s[scenario]=year2mon(dep_data[scenario])
    
    
    grid_nc=gridarea()
    
    
    # mlons
    lons=np.linspace(-179.75,179.75,720)
    #lon_bounds=np.linspace(-180,180,721)
    lats=np.linspace(-89.75,89.75,360)
    
    a=regionmask.defined_regions.natural_earth_v5_0_0.countries_50
    b=a.mask(lons,lats)
    china_index=regionmask.defined_regions.natural_earth_v5_0_0.countries_50.map_keys('China')
    US_index=regionmask.defined_regions.natural_earth_v5_0_0.countries_50.map_keys('US')
    EU_index=regionmask.defined_regions.natural_earth_v5_0_0.countries_50.map_keys(['France','Spain','Netherlands','Belgium','Germany','Norway','Sweden',
                      'Italy','Greece','Estonia','Lithuania','Latvia','Denmark','Finland','Poland','Romania','Bulgaria','Croatia','Malta',
                      'United Kingdom','Ireland','Portugal','Slovenia'])
    temp=grid_nc.copy()
    china_m2=temp.m2.where(b==china_index)
    
    temp=grid_nc.copy()
    US_m2=temp.m2.where(b==US_index)
    temp=grid_nc.copy()
    EU_m2=temp.copy()
    EU_m2.data[:,:]=np.nan
    for eu_i in EU_index:
        print(eu_i,(b.data==eu_i).sum())
        print()
        EU_m2.data=np.where(b.data==eu_i,temp.data,EU_m2.data)#EU_m2.data)
    
    plt.figure()
    EU_m2.plot()
    plt.figure()
    china_m2.plot()
    
    
    ## exclusive economical zones as geopandas
    eez, maskeez=create_eez_mask()
    
    
    included_eez_dict=select_areas(eez)
    
    time_dates = np.arange(np.datetime64("2015-01-01"), np.datetime64("2101-01-01"),  np.timedelta64(1, 'M'), dtype='datetime64[M]')
    
    
    
    # prepare data for mask
    data=maskeez.copy()
    
    
    # create two copies with same dimensions
    data2=data.copy()
    data2.data[:,:]=np.nan
    
    for i in included_eez_dict:
        print ('incl',i)
        if included_eez_dict[i]=='China':
            print('Computing China')
            data2.data=np.where((data.data==float(i)),10,data2.data)
    
        elif included_eez_dict[i]=='United States of America':
            print('Computing USA')
            data2.data=np.where((data.data==float(i)),20,data2.data)
        else:
            print('Computing EU')
            data2.data=np.where((data.data==float(i)),30,data2.data)
    
    data2.data=np.where(~np.isnan(data2.data),data2.data,np.nan)
    data2=masking_boxes(data2)
    #data_penalty.data=np.where(~np.isnan(data_penalty.data),data_penalty.data,np.nan)
    print(data2)
    # calculate the area of sea used for different parts
    area_china=float(grid_nc.where(data2.data==10).sum())
    area_us=float(grid_nc.where(data2.data==20).sum())
    area_eu=float(grid_nc.where(data2.data==30).sum())
    print(f"CHINA: {area_china:.2f}")
    print(f"US: {area_us:.2f}")
    print(f"EU: {area_eu:.2f}")
    
    
    # expand time dimension
    data2=data2.expand_dims({'time':time_dates}).copy()
    
    
    
    # Change to liming data as mol m-2 s-1
    dep_data_monthly_mol_s_m2[scenario]=dep_data_monthly_mol_s[scenario].copy()
    dep_data_monthly_mol_s_m2[scenario]['dep-China']=dep_data_monthly_mol_s_m2[scenario]['dep-China']/area_china
    dep_data_monthly_mol_s_m2[scenario]['dep-US']=dep_data_monthly_mol_s_m2[scenario]['dep-US']/area_us
    dep_data_monthly_mol_s_m2[scenario]['dep-Europe']=dep_data_monthly_mol_s_m2[scenario]['dep-Europe']/area_eu
    
    
    
    for i in range(data2.data.shape[0]):
    
        # print(i,time_dates[i],china_lime_yearly[i], penalty_co2_yearly[i])
        imonth=time_dates[i].astype(object).month -1
        iyear=time_dates[i].astype(object).year
        isleap=calendar.isleap(iyear)
    
        
        for region, region_id in zip(['China','US','Europe'],[10,20,30]):
            data2.data[i,:,:]=np.where((data2.data[i,:,:]==region_id),dep_data_monthly_mol_s_m2[scenario]['dep-'+region].iloc[i],data2.data[i,:,:])
    
        
    
    fig,ax=plt.subplots(ncols=1,subplot_kw={'projection': ccrs.Robinson()})
    
    ax.coastlines(color='white')
    
    data3=data2.copy()#.plot(ax=ax)
    data3.data=np.where(data2.data[:,:]>0,1,0)#.plot(ax=ax)
    data3[-1,:,:].plot.pcolormesh(ax=ax,transform=ccrs.PlateCarree(), add_colorbar=False)
    ax.coastlines(color='white')
    plt.savefig('maski.png')
    plt.show()
    
    mask=xr.Dataset(
        data_vars=dict(
            lime_mask=(['time','lat','lon'],data2.data,{'units':'mol m-2 s-1'})
            #gridarea=(['lat','lon'],grid*lat_size,{'units':'m-2'})
            #mask=(['time','lat','lon'],data2.data,{'units':'country number'})
            #data_vars=(['lon','lat'],data,{'units':'mol/s'})
            ),
        coords=dict(
            time=('time',time_dates),
            lon=('lon',lons,{'units':'degrees_east'}),
            lat=('lat',lats,{'units':'degrees_north'})
            ),
        attrs=dict(description='mask for ocean alkalinization for WP4 of OceanNETs project.\
                   Based on the world exclusive economic zones (marineregions.org) and work done by Spyros Foteinis and Phil Renforth \
                       for annual lime production numbers.',
                   author='Tommi Bergman (FMI)'
            )
    )
    
    # write to disk, this includes
    mask.to_netcdf(f'../data/lime_mask_cao_{scenario}.nc')
    
    plt.figure()
    plt.title('eez')
    eez.plot()
    plt.figure()
    
    maskeez.plot()
    plt.figure()
    mask.lime_mask[:,:].plot()
    #mask.mask[0,:,:].plot()
    print(mask.lat)
    mask2=mask
# =============================================================================

# =============================================================================
# print ('maski',grid_nc.where(mask2>0).isel(time=1000).sum())
# print((mask2*grid_nc).isel(time=1000).sum())
# plt.figure()
# grid_nc.where(mask2>0).m2.plot()
# plt.figure()
# #mask2.mask[0,:,:].plot()
# fig,ax=plt.subplots(ncols=1,subplot_kw={'projection': ccrs.Robinson()})
# #mask2.lime_mask[:,:].plot(transform=ccrs.PlateCarree(),add_labels=True)
# ax.coastlines(color='white')
# mask2.to_netcdf('../data/limemask_v2.nc',engine='scipy')

# # regrid to ocean grid
# regrid_to_model(mask2)
# # reference grid
# ingrid=xr.open_dataset('../../data/bathy_meter.nc')
# # makesure lon and lat is proper
# ingrid.rename({'nav_lon':'lon', 'nav_lat':'lat'})
# # regridder creation
# regridder = xe.Regridder(mask2, ingrid, 'bilinear', periodic=True) #since this is global we need to pass periodic
# # regrid
# mask_orca1=regridder(mask2)

# # read reference grid arera
# orcam2=xr.open_dataset('../../data/areacello.nc')
# plt.figure()
# plt.plot(mask_orca1.time,(mask_orca1*orcam2.areacello).lime_mask.sum(dim=('x','y'))*n_ca_oh_2*365*3600*24,label='orca1')
# plt.plot(mask2.time,(mask2*grid_nc).lime_mask.sum(dim=('lon','lat'))*n_ca_oh_2*365*3600*24,label='lonlat')
# plt.legend()
# plt.title('diff')

# # calculate the error from remgridding
# error=((mask2*grid_nc).lime_mask.sum(dim=('lon','lat'))*n_ca_oh_2*365*3600*24/((mask_orca1*orcam2.areacello).lime_mask.sum(dim=('x','y'))*n_ca_oh_2*365*3600*24)).isel(time=range(1020,1031)).mean()
# error=((mask2*grid_nc).lime_mask.sum(dim=('lon','lat'))*n_ca_oh_2*365*3600*24/((mask_orca1*orcam2.areacello).lime_mask.sum(dim=('x','y'))*n_ca_oh_2*365*3600*24)).mean()
# # correct the grid calues with the error
# mask_orca1=mask_orca1*error
# plt.figure()
# plt.plot(mask_orca1.time,(mask_orca1*orcam2.areacello).lime_mask.sum(dim=('x','y'))*n_ca_oh_2*365*3600*24,label='orca')
# #plt.plot(mask2.time,(mask2*grid_nc).lime_mask.sum(dim=('lon','lat'))*n_ca_oh_2*365*3600*24,label='latlon')
# plt.legend()
# plt.title('diff')
# plt.figure()
# plt.plot((mask_orca1*orcam2.areacello).lime_mask.sum(dim=('x','y'))*n_ca_oh_2*365*3600*24/((mask2*grid_nc).lime_mask.sum(dim=('lon','lat'))*n_ca_oh_2*365*3600*24))

# print(((mask2*grid_nc).lime_mask.sum(dim=('lon','lat'))*n_ca_oh_2*365*3600*24/((mask_orca1*orcam2.areacello).lime_mask.sum(dim=('x','y'))*n_ca_oh_2*365*3600*24)).isel(time=range(1020,1031)).mean())

# mask_orca1.to_netcdf('../../data/limemask_v2.orca1.cao.nc',engine='scipy')
# plt.show()#