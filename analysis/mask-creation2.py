#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 09:38:10 2022

@author: bergmant
"""
import numpy as np
import xarray as xr
import geopandas as gpd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import regionmask
import csv
from datetime import datetime, timedelta
import xesmf as xe
# Grid size calculation from Stackoverflow
# https://stackoverflow.com/questions/61338665/how-to-calculate-size-in-m2-of-each-lat-long-grid-square
def gridsize(lat1):
   #https://en.wikipedia.org/wiki/Haversine_formula
   #https://stackoverflow.com/questions/639695/how-to-convert-latitude-or-longitude-to-meters/11172685#11172685
   lon1=200
   import math
   lat2=lat1
   lon2=lon1+1

   R = 6378.137 # // Radius of earth in km
   dLat = lat2 * np.pi / 180 - lat1 * np.pi / 180
   dLon = lon2 * np.pi / 180 - lon1 * np.pi / 180
   a = np.sin(dLat/2) * np.sin(dLat/2) + np.cos(lat1 * np.pi / 180) * np.cos(lat2 * np.pi / 180) * np.sin(dLon/2) * np.sin(dLon/2)
   c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a))
   d = R * c
   return d * 1000 #; // meters
def regrid_to_model(indata,ifile='/Users/bergmant/Documents/projects/H2020-OceanNET/data/bathy_meter.nc'):

    ingrid=xr.open_dataset(ifile)
    if 'nav_lon' in ingrid:
        ingrid.rename({'nav_lon':'lon', 'nav_lat':'lat'})
    regridder = xe.Regridder(indata, ingrid, 'bilinear', periodic=True) #since this is global we need to pass periodic
    mask_orca1=regridder(indata)
    mask_orca1.to_netcdf('/Users/bergmant/Documents/projects/H2020-OceanNET/scripts/limemask.orca1.nc',engine='scipy')
    return mask_orca1

def eez():
    ## exclusive economical zones as geopandas
    eez=gpd.read_file('/Users/bergmant/Documents/projects/H2020-OceanNET/data/World_EEZ_v11_20191118_HR_0_360/eez_v11_0_360.shp')
    print(eez.head())
    lon = np.arange(-180, 180)
    lat = np.arange(-90, 90)
    # print(lons)
    # print(lon)
    # print(lats)
    # print(lat)
    #maskeez = regionmask.mask_geopandas(eez, lon, lat,wrap_lon=True)
    maskeez = regionmask.mask_geopandas(eez, lons, lats,wrap_lon=True)
    return maskeez

def emissions_per_area():
    n_ca_oh_2=74.1 #g/mol
    china_lime=np.array([0,0,0,0.14,0.267718,0.413041,0.571405,0.738835,0.91324,
                         1.09045,1.266208,1.438202,1.606129,1.76801,1.921216,2.0641,2.197527,2.323329])*1e6*1e3*1e3 #g/year

    lime_year=np.linspace(2015,2100,china_lime.size)
    all_years=np.linspace(2015,2100,(2100-2015+1)*12)
    china_lime_yearly=np.interp(all_years, lime_year, china_lime)
    eu_lime_yearly=china_lime_yearly.copy()
    us_lime_yearly=china_lime_yearly.copy()

    print(china_lime_yearly)
    return china_lime_yearly,eu_lime_yearly,us_lime_yearly
def gridarea():
    reso=0.5 #degrees
    mlons,mlats=np.meshgrid(np.linspace(-179.75,179.75,720),np.linspace(-89.75,89.75,360))
    mlons_b,mlats_b=np.meshgrid(np.linspace(-180,180,721),np.linspace(-90,90,361))
    #grid=gridsize(mlats)
    grid=np.cos(np.radians(abs(mlats)))*((2*np.pi*6378.137*reso)*(2*np.pi*6378.137*reso)*1000*1000)
    print(grid)
    grid_nc=xr.DataArray(grid,coords={'lat':mlats[:,1],'lon':mlons[1,:]},dims=['lat','lon'])
    lat_size=110567 #in m
    grid_nc['m2']=grid_nc*lat_size
    grid_nc=grid_nc['m2']
    grid_nc.to_netcdf('earth_m2.nc')
    plt.pcolormesh(mlons[1,:],mlats[:,1],grid_nc)
    plt.colorbar()
    #plt.show()
    return grid,grid_nc
def select_areas():
    ## Zones to be excluded from alkalinity simulations
    excluded_eez=[]
    included_eez=[]
    included_eez_dict={}
    included_region=[]
    included_country=['United States of America','China','France','Spain','Netherlands','Belgium','Germany','Norway','Sweden',
                      'Italy','Greece','Israel','Türkiye','Estonia','Lithuania','Latvia','Denmark','Finland','Poland','Romania','Bulgaria','Croatia','Malta',
                      'United Kingdom of Great Britain and Northern Ireland','Ireland','Portugal']
    with open('/Users/bergmant/Documents/projects/H2020-OceanNET/data/UNSD — Methodology - country numbers for EEZ.csv') as csvfile:
        reader= csv.DictReader(csvfile,delimiter=';')
        for row in reader:
            #if row['Sub-region Name']=='Micronesia':
            if row['Region Name'] not in included_region:
                for area in included_country:
                    #if row['Country or Area'] not in included_country:
                    #print(area)
                    if row['Country or Area'] == area:
                        #print(float(row['M49 Code']), eez['UN_SOV1'])
                        #print(row['Region Name'],row['Country or Area'])
                        if eez['UN_SOV1'][eez['UN_SOV1']==float(row['M49 Code'])].size>0:
                            for i in eez['UN_SOV1'][eez['UN_SOV1']==float(row['M49 Code'])].index:
                                if row['Country or Area']!='Greenland':
                                    print(i,row['Country or Area'])
                                    included_eez.append(i)
                                    included_eez_dict[i]=row['Country or Area']
                                else:
                                    print('greenland')
                        else:
                            print('nocoast',eez['UN_SOV1'][eez['UN_SOV1']==float(row['M49 Code'])].size,float(row['M49 Code']))
                    # else:
                    #     if eez['UN_SOV1'][eez['UN_SOV1']==float(row['M49 Code'])].size>0:
                    #         print('x',row['M49 Code'],eez['UN_SOV1'][eez['UN_SOV1']==float(row['M49 Code'])].index[0],row['Country or Area'])
            # for area in included_region:
            #     #if row['Country or Area'] not in included_country:
            #     if row['Region Name'] == area:
            #         #print(float(row['M49 Code']), eez['UN_SOV1'])
            #         print( row['Region Name'],row['Country or Area'])
            #         if eez['UN_SOV1'][eez['UN_SOV1']==float(row['M49 Code'])].size>0:
            #             for i in eez['UN_SOV1'][eez['UN_SOV1']==float(row['M49 Code'])].index:
            #                 print ('ttt',i,row['M49 Code'],row['Country or Area'],row['Sub-region Name'],row['ISO-alpha3 Code'])# not in [246]:
            #                 if row['Country or Area']!='Greenland':
            #                     included_eez.append(i)
            #         else:
            #             print('nocoast',row['Country or Area'],eez['UN_SOV1'][eez['UN_SOV1']==float(row['M49 Code'])].size,float(row['M49 Code']))
            #     else:
            #        if eez['UN_SOV1'][eez['UN_SOV1']==float(row['M49 Code'])].size>0:
            #            print('x',row['M49 Code'],eez['UN_SOV1'][eez['UN_SOV1']==float(row['M49 Code'])].index[0],row['Country or Area'])### time definitions for the mask

def masking_boxes(mask):
    # southern hemisphere
    mask2=mask.where(mask.lat>0,np.nan)
    #Pacific equatorial
    masklon= ( mask2.lon<-150)
    masklat= (mask2.lat>-40) & (mask2.lat<40)
    mask2=mask2.where( ~(masklon&masklat),np.nan)
    #Pacific north
    masklon= ( mask2.lon>130)
    masklat= ( mask2.lat<40)
    mask2=mask2.where(~(masklon&masklat),np.nan)
    # pacific1
    masklon= ( np.bitwise_and(mask2.lon>-100,mask2.lon<0))
    masklat= ( mask2.lat<21)
    mask2=mask2.where(~np.bitwise_and(masklon,masklat),np.nan)
    # atlantic 1
    masklon= ( np.bitwise_and(mask2.lon>-71 , mask2.lon<-15))
    masklat= ( mask2.lat>55)
    mask2=mask2.where(~np.bitwise_and(masklon,masklat),np.nan)
    # Caribbean
    masklon= ( np.bitwise_and(mask2.lon<-90 , mask2.lon>-125))
    masklat= np.bitwise_and( mask2.lat>0,mask2.lat<20)
    mask2=mask2.where(~np.bitwise_and(masklon,masklat),np.nan)
    #Atlantic
    masklon= ( np.bitwise_and(mask2.lon>-71 , mask2.lon<-12))
    masklat= np.bitwise_and( mask2.lat>0,mask2.lat<35)
    mask2=mask2.where(~np.bitwise_and(masklon,masklat),np.nan)
    # Atlantic
    masklon= ( np.bitwise_and(mask2.lon>-51 , mask2.lon<-20))
    masklat= np.bitwise_and( mask2.lat>30,mask2.lat<55)
    mask2=mask2.where(~np.bitwise_and(masklon,masklat),np.nan)
    # Atlantic
    masklon= ( np.bitwise_and(mask2.lon>-51 , mask2.lon<-14))
    masklat= np.bitwise_and( mask2.lat>20,mask2.lat<40)
    mask2=mask2.where(~np.bitwise_and(masklon,masklat),np.nan)
    # Caribbean
    masklon= ( np.bitwise_and(mask2.lon>-67 , mask2.lon<-50))
    masklat= np.bitwise_and( mask2.lat>20,mask2.lat<40)
    mask2=mask2.where(~np.bitwise_and(masklon,masklat),np.nan)
    #Caribbean
    masklon= ( np.bitwise_and(mask2.lon>-79 , mask2.lon<-50))
    masklat= np.bitwise_and( mask2.lat>18,mask2.lat<25)
    mask2=mask2.where(~np.bitwise_and(masklon,masklat),np.nan)
    #
    masklon= ( np.bitwise_and(mask2.lon>-79 , mask2.lon<-50))
    masklat= np.bitwise_and( mask2.lat>18,mask2.lat<25)
    mask2=mask2.where(~np.bitwise_and(masklon,masklat),np.nan)

    # Baltic sea
    masklon= ( np.bitwise_and(mask2.lon>9.5 , mask2.lon<16))
    masklat= np.bitwise_and( mask2.lat>50,mask2.lat<60)
    mask2=mask2.where(~np.bitwise_and(masklon,masklat),np.nan)
    # Baltic sea 2
    masklon= ( np.bitwise_and(mask2.lon>16 , mask2.lon<50))
    masklat= np.bitwise_and( mask2.lat>50,mask2.lat<80)
    mask2=mask2.where(~np.bitwise_and(masklon,masklat),np.nan)

    # mediterranean sea
    masklon= ( np.bitwise_and(mask2.lon>0 , mask2.lon<60))
    masklat= np.bitwise_and( mask2.lat>15,mask2.lat<48)
    mask2=mask2.where(~np.bitwise_and(masklon,masklat),np.nan)
    # mediterranean sea
    masklon= ( np.bitwise_and(mask2.lon>-5.5 , mask2.lon<0))
    masklat= np.bitwise_and( mask2.lat>15,mask2.lat<40)
    mask2=mask2.where(~np.bitwise_and(masklon,masklat),np.nan)

    # Northpole
    #masklon= ( np.bitwise_and(mask2.lon>-70 , mask2.lon<10))
    masklat= ( mask2.lat>67.5)
    mask2=mask2.where(~masklat,np.nan)
    # masklon= ( mask2.lon>-50 & mask2.lon<0)
    # masklat= ( mask2.lat>50)
    # mask2=mask2.where(~(masklon&masklat),np.nan)
    print(mask2)
    mask2=mask2.where(~mask2.isnull(),0)
    return mask2
# molar mass of Ca(OH)2
n_ca_oh_2=74.1 #g/mol
# Evolution of excess lime
china_lime=np.array([0,0,0,0.14,0.267718,0.413041,0.571405,0.738835,0.91324,
                     1.09045,1.266208,1.438202,1.606129,1.76801,1.921216,2.0641,2.197527,2.323329])*1e6*1e3*1e3 #g/year
# Liming year dimension
lime_year=np.linspace(2015,2100,china_lime.size)
# liming monthly dimension
all_years=np.linspace(2015,2100,(2100-2015+1)*12)
# interpolate to monthly
china_lime_yearly=np.interp(all_years, lime_year, china_lime)

# For now copy china to EU and USA
eu_lime_yearly=china_lime_yearly.copy()
us_lime_yearly=china_lime_yearly.copy()
# check
print(china_lime_yearly)
plt.plot(all_years,china_lime_yearly)
plt.plot(lime_year,china_lime)
# resolution of initial latlon-mask
reso=0.5 #degrees
mlons,mlats=np.meshgrid(np.linspace(-179.75,179.75,720),np.linspace(-89.75,89.75,360))
mlons_b,mlats_b=np.meshgrid(np.linspace(-180,180,721),np.linspace(-90,90,361))
#grid=gridsize(mlats)
# create grid data
grid=np.cos(np.radians(abs(mlats)))*(111100*reso*111100*reso)
print(grid)
# create data array
grid_nc=xr.DataArray(grid,coords={'lat':mlats[:,1],'lon':mlons[1,:]},dims=['lat','lon'])
#lat_size=110567*reso #in m
#create dataset
grid_nc['m2']=grid_nc#*lat_size
grid_nc=grid_nc['m2']
# write to disk
grid_nc.to_netcdf('earth_m2.nc')
#plot
plt.pcolormesh(mlons[1,:],mlats[:,1],grid_nc)
plt.colorbar()
#plt.show()

mlons
lons=np.linspace(-179.75,179.75,720)
#lon_bounds=np.linspace(-180,180,721)
lats=np.linspace(-89.75,89.75,360)
print(lons.shape)
print(lats.shape)

## exclusive economical zones as geopandas
eez=gpd.read_file('/Users/bergmant/Documents/projects/H2020-OceanNET/data/World_EEZ_v11_20191118_HR_0_360/eez_v11_0_360.shp')
print(eez.head())
lon = np.arange(-180, 180)
lat = np.arange(-90, 90)
# print(lons)
# print(lon)
# print(lats)
# print(lat)
#maskeez = regionmask.mask_geopandas(eez, lon, lat,wrap_lon=True)
maskeez = regionmask.mask_geopandas(eez, lons, lats,wrap_lon=True)


## Zones to be excluded from alkalinity simulations
excluded_eez=[]
## Zones to be included
included_eez=[]
included_eez_dict={}
included_region=[]
included_country=['United States of America','China','France','Spain','Netherlands','Belgium','Germany','Norway','Sweden',
                  'Italy','Greece','Israel','Türkiye','Estonia','Lithuania','Latvia','Denmark','Finland','Poland','Romania','Bulgaria','Croatia','Malta',
                  'United Kingdom of Great Britain and Northern Ireland','Ireland','Portugal']
with open('/Users/bergmant/Documents/projects/H2020-OceanNET/data/UNSD — Methodology - country numbers for EEZ.csv') as csvfile:
    reader= csv.DictReader(csvfile,delimiter=';')
    for row in reader:
        #if row['Sub-region Name']=='Micronesia':
        if row['Region Name'] not in included_region:
            for area in included_country:
                #if row['Country or Area'] not in included_country:
                #print(area)
                if row['Country or Area'] == area:
                    #print(float(row['M49 Code']), eez['UN_SOV1'])
                    #print(row['Region Name'],row['Country or Area'])
                    if eez['UN_SOV1'][eez['UN_SOV1']==float(row['M49 Code'])].size>0:
                        for i in eez['UN_SOV1'][eez['UN_SOV1']==float(row['M49 Code'])].index:
                            if row['Country or Area']!='Greenland':
                                print(i,row['Country or Area'])
                                included_eez.append(i)
                                included_eez_dict[i]=row['Country or Area']
                            else:
                                print('greenland')
                    else:
                        print('nocoast',eez['UN_SOV1'][eez['UN_SOV1']==float(row['M49 Code'])].size,float(row['M49 Code']))
                # else:
                #     if eez['UN_SOV1'][eez['UN_SOV1']==float(row['M49 Code'])].size>0:
                #         print('x',row['M49 Code'],eez['UN_SOV1'][eez['UN_SOV1']==float(row['M49 Code'])].index[0],row['Country or Area'])
        # for area in included_region:
        #     #if row['Country or Area'] not in included_country:
        #     if row['Region Name'] == area:
        #         #print(float(row['M49 Code']), eez['UN_SOV1'])
        #         print( row['Region Name'],row['Country or Area'])
        #         if eez['UN_SOV1'][eez['UN_SOV1']==float(row['M49 Code'])].size>0:
        #             for i in eez['UN_SOV1'][eez['UN_SOV1']==float(row['M49 Code'])].index:
        #                 print ('ttt',i,row['M49 Code'],row['Country or Area'],row['Sub-region Name'],row['ISO-alpha3 Code'])# not in [246]:
        #                 if row['Country or Area']!='Greenland':
        #                     included_eez.append(i)
        #         else:
        #             print('nocoast',row['Country or Area'],eez['UN_SOV1'][eez['UN_SOV1']==float(row['M49 Code'])].size,float(row['M49 Code']))
        #     else:
        #        if eez['UN_SOV1'][eez['UN_SOV1']==float(row['M49 Code'])].size>0:
        #            print('x',row['M49 Code'],eez['UN_SOV1'][eez['UN_SOV1']==float(row['M49 Code'])].index[0],row['Country or Area'])### time definitions for the mask
nyears=1
time=np.linspace(0,11,12*nyears)
rampup_time=48 #months
#time_dates = np.arange(datetime(2015,1,1), datetime(2100,1,1), timedelta(months=1)).astype(datetime)
# dates for dataset
time_dates = np.arange(np.datetime64("2015-01-01"), np.datetime64("2100-12-31"),  np.timedelta64(1, 'M'), dtype='datetime64[M]')
#print(excluded_eez)
# prepare data for mask
#data=np.zeros((nyears*12,360,720))
data=maskeez.copy()#expand_dims({'time':time_dates.size})
#data=np.where(data.data>0,1,0)
#for i in excluded_eez:
# create copy with same dimensions
data2=data.copy()
data2.data[:,:]=np.nan
#data2.data[:,:,:]=np.nan
#print (data2)
#for i in range(0,150):
for i in included_eez_dict:
    print ('incl',i)
    if included_eez_dict[i]=='China':
        print('ttt China')
        data2.data=np.where((data.data==float(i)),10,data2.data)
    elif included_eez_dict[i]=='United States of America':
        print('USA')
        data2.data=np.where((data.data==float(i)),20,data2.data)
    else:
        data2.data=np.where((data.data==float(i)),30,data2.data)
data2.data=np.where(~np.isnan(data2.data),data2.data,np.nan)
print(data2)
# calculate the area of sea used for different parts
area_china=grid_nc.where(data2.data==10).sum()
area_us=grid_nc.where(data2.data==20).sum()
area_eu=grid_nc.where(data2.data==30).sum()
print(area_china,grid_nc.sum())
print(area_us,grid_nc.sum())
print(area_eu,grid_nc.sum())
plt.figure()
data2.plot()
data2=masking_boxes(data2)
plt.figure()
data2.plot()
area_china=grid_nc.where(data2.data==10).sum()
area_us=grid_nc.where(data2.data==20).sum()
area_eu=grid_nc.where(data2.data==30).sum()
print(area_china,grid_nc.sum())
print(area_us,grid_nc.sum())
print(area_eu,grid_nc.sum())
# unit change
yearly_kg_2_mol_s=1/n_ca_oh_2/(365*3600*24)
#data2=data2.expand_dims({'time':time_dates.size}).copy()
# expand time dimension
data2=data2.expand_dims({'time':time_dates}).copy()

# Change to liming data as mol m-2 s-1
for i in range(data2.data.shape[0]):
    print(i)
    data2.data[i,:,:]=np.where((data2.data[i,:,:]==10),china_lime_yearly[i]/area_china*yearly_kg_2_mol_s,data2.data[i,:,:])
    data2.data[i,:,:]=np.where((data2.data[i,:,:]==20),us_lime_yearly[i]/area_us*yearly_kg_2_mol_s,data2.data[i,:,:])
    data2.data[i,:,:]=np.where((data2.data[i,:,:]==30),eu_lime_yearly[i]/area_eu*yearly_kg_2_mol_s,data2.data[i,:,:])

#data=np.where(data.data>0,1,0)
#print (data2)
#print (data2.data)
#print (data2.flags)
#china_limec
#data2.data=np.where(data2.data==30,china_lime_yaerly,data2.data)
# time_ramp=time_dates[:rampup_time]
# for i,itime in enumerate(time_ramp):
#     print (i,itime)
#     data2.data[i,:,:]=data2.data[i,:,:]*float((i+1)/rampup_time)
#     #np.where((data.data==float(i)),data.data,data2.data)
# mask dataset
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
               Based on the world exclusive economic zones (marineregions.org) and work done WPX \
                   for annual lime production numbers.',
               author='Tommi Bergman (FMI)'
        )
)
# for i in range(100,150):
#     print (i)
#     mask.mask.data=np.where(np.abs(mask.mask.data-i)<0.5,1,mask.mask.data)
# write to disk, this includes
mask.to_netcdf('/Users/bergmant/Documents/projects/H2020-OceanNET/scripts/test.nc')
plt.figure()
eez.plot()
plt.figure()

maskeez.plot()
plt.figure()
mask.lime_mask[:,:].plot()
#mask.mask[0,:,:].plot()
print(mask.lat)
mask2=mask
# =============================================================================
# # southern hemisphere
# mask2=mask.where(mask.lat>0,np.nan)
# #Pacific equatorial
# masklon= ( mask2.lon<-150)
# masklat= (mask2.lat>-40) & (mask2.lat<40)
# mask_1=(masklon & masklon)
# mask2=mask2.where( ~(masklon&masklat),np.nan)
# #Pacific north
# masklon= ( mask2.lon>130)
# masklat= ( mask2.lat<40)
# mask2=mask2.where(~(masklon&masklat),np.nan)
# # pacific1
# masklon= ( np.bitwise_and(mask2.lon>-100,mask2.lon<0))
# masklat= ( mask2.lat<21)
# mask2=mask2.where(~np.bitwise_and(masklon,masklat),np.nan)
# # atlantic 1
# masklon= ( np.bitwise_and(mask2.lon>-71 , mask2.lon<-15))
# masklat= ( mask2.lat>55)
# mask2=mask2.where(~np.bitwise_and(masklon,masklat),np.nan)
# # Caribbean
# masklon= ( np.bitwise_and(mask2.lon<-90 , mask2.lon>-125))
# masklat= np.bitwise_and( mask2.lat>0,mask2.lat<20)
# mask2=mask2.where(~np.bitwise_and(masklon,masklat),np.nan)
# #Atlantic
# masklon= ( np.bitwise_and(mask2.lon>-71 , mask2.lon<-12))
# masklat= np.bitwise_and( mask2.lat>0,mask2.lat<35)
# mask2=mask2.where(~np.bitwise_and(masklon,masklat),np.nan)
# # Atlantic
# masklon= ( np.bitwise_and(mask2.lon>-51 , mask2.lon<-20))
# masklat= np.bitwise_and( mask2.lat>30,mask2.lat<55)
# mask2=mask2.where(~np.bitwise_and(masklon,masklat),np.nan)
# # Atlantic
# masklon= ( np.bitwise_and(mask2.lon>-51 , mask2.lon<-14))
# masklat= np.bitwise_and( mask2.lat>20,mask2.lat<40)
# mask2=mask2.where(~np.bitwise_and(masklon,masklat),np.nan)
# # Caribbean
# masklon= ( np.bitwise_and(mask2.lon>-67 , mask2.lon<-50))
# masklat= np.bitwise_and( mask2.lat>20,mask2.lat<40)
# mask2=mask2.where(~np.bitwise_and(masklon,masklat),np.nan)
# #Caribbean
# masklon= ( np.bitwise_and(mask2.lon>-79 , mask2.lon<-50))
# masklat= np.bitwise_and( mask2.lat>18,mask2.lat<25)
# mask2=mask2.where(~np.bitwise_and(masklon,masklat),np.nan)
# #
# masklon= ( np.bitwise_and(mask2.lon>-79 , mask2.lon<-50))
# masklat= np.bitwise_and( mask2.lat>18,mask2.lat<25)
# mask2=mask2.where(~np.bitwise_and(masklon,masklat),np.nan)
#
# # Baltic sea
# masklon= ( np.bitwise_and(mask2.lon>9.5 , mask2.lon<16))
# masklat= np.bitwise_and( mask2.lat>50,mask2.lat<60)
# mask2=mask2.where(~np.bitwise_and(masklon,masklat),np.nan)
# # Baltic sea 2
# masklon= ( np.bitwise_and(mask2.lon>16 , mask2.lon<50))
# masklat= np.bitwise_and( mask2.lat>50,mask2.lat<80)
# mask2=mask2.where(~np.bitwise_and(masklon,masklat),np.nan)
#
# # mediterranean sea
# masklon= ( np.bitwise_and(mask2.lon>0 , mask2.lon<60))
# masklat= np.bitwise_and( mask2.lat>15,mask2.lat<48)
# mask2=mask2.where(~np.bitwise_and(masklon,masklat),np.nan)
# # mediterranean sea
# masklon= ( np.bitwise_and(mask2.lon>-5.5 , mask2.lon<0))
# masklat= np.bitwise_and( mask2.lat>15,mask2.lat<40)
# mask2=mask2.where(~np.bitwise_and(masklon,masklat),np.nan)
#
# # Northpole
# #masklon= ( np.bitwise_and(mask2.lon>-70 , mask2.lon<10))
# masklat= ( mask2.lat>67.5)
# mask2=mask2.where(~masklat,np.nan)
# print ('maski',grid_nc.m2.where(~mask2.lime_mask.isnull()).sum())
# # masklon= ( mask2.lon>-50 & mask2.lon<0)
# # masklat= ( mask2.lat>50)
# # mask2=mask2.where(~(masklon&masklat),np.nan)
# mask2=mask2.where(~mask2.lime_mask.isnull(),0)
# =============================================================================
print ('maski',grid_nc.where(mask2>0).isel(time=1000).sum())
print((mask2*grid_nc).isel(time=1000).sum())
plt.figure()
grid_nc.where(mask2>0).m2.plot()
plt.figure()
#mask2.mask[0,:,:].plot()
fig,ax=plt.subplots(ncols=1,subplot_kw={'projection': ccrs.Robinson()})
#mask2.lime_mask[:,:].plot(transform=ccrs.PlateCarree(),add_labels=True)
ax.coastlines(color='white')
mask2.to_netcdf('/Users/bergmant/Documents/projects/H2020-OceanNET/scripts/limemask.nc',engine='scipy')

# regrid to ocean grid
regrid_to_model(mask2)
# reference grid
ingrid=xr.open_dataset('../../data/bathy_meter.nc')
# makesure lon and lat is proper
ingrid.rename({'nav_lon':'lon', 'nav_lat':'lat'})
# regridder creation
regridder = xe.Regridder(mask2, ingrid, 'bilinear', periodic=True) #since this is global we need to pass periodic
# regrid
mask_orca1=regridder(mask2)

# read reference grid arera
orcam2=xr.open_dataset('areacello.nc')
plt.figure()
plt.plot((mask_orca1*orcam2.areacello).lime_mask.sum(dim=('x','y'))*n_ca_oh_2*365*3600*24)
plt.plot((mask2*grid_nc).lime_mask.sum(dim=('lon','lat'))*n_ca_oh_2*365*3600*24)

# calculate the error from remgridding
error=((mask2*grid_nc).lime_mask.sum(dim=('lon','lat'))*n_ca_oh_2*365*3600*24/((mask_orca1*orcam2.areacello).lime_mask.sum(dim=('x','y'))*n_ca_oh_2*365*3600*24)).isel(time=range(1020,1031)).mean()
error=((mask2*grid_nc).lime_mask.sum(dim=('lon','lat'))*n_ca_oh_2*365*3600*24/((mask_orca1*orcam2.areacello).lime_mask.sum(dim=('x','y'))*n_ca_oh_2*365*3600*24)).mean()
# correct the grid calues with the error
mask_orca1=mask_orca1*error
plt.figure()
plt.plot((mask_orca1*orcam2.areacello).lime_mask.sum(dim=('x','y'))*n_ca_oh_2*365*3600*24,label='orca')
plt.plot((mask2*grid_nc).lime_mask.sum(dim=('lon','lat'))*n_ca_oh_2*365*3600*24,label='latlon')
plt.legend()
plt.figure()
plt.plot((mask_orca1*orcam2.areacello).lime_mask.sum(dim=('x','y'))*n_ca_oh_2*365*3600*24/((mask2*grid_nc).lime_mask.sum(dim=('lon','lat'))*n_ca_oh_2*365*3600*24))

print(((mask2*grid_nc).lime_mask.sum(dim=('lon','lat'))*n_ca_oh_2*365*3600*24/((mask_orca1*orcam2.areacello).lime_mask.sum(dim=('x','y'))*n_ca_oh_2*365*3600*24)).isel(time=range(1020,1031)).mean())

mask_orca1.to_netcdf('/Users/bergmant/Documents/projects/H2020-OceanNET/scripts/limemask.orca1.nc',engine='scipy')