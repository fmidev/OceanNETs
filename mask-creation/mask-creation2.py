#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 09:38:10 2022

@author: bergmant
"""
import numpy as np
import xarray as xr
import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import regionmask
import csv
from datetime import datetime, timedelta
import calendar
#import ESMpy
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
    mask_orca1.to_netcdf('/Users/bergmant/Documents/projects/H2020-OceanNET/scripts/limemask_v2.orca1.nc',engine='scipy')
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
    # Atlantic Saint-Pierre-et-Miquelon
    masklon= ( np.bitwise_and(mask2.lon>-58 , mask2.lon<-55))
    masklat= np.bitwise_and( mask2.lat>42,mask2.lat<48)
    mask2=mask2.where(~np.bitwise_and(masklon,masklat),np.nan)
    # Caribbean
    masklon= ( np.bitwise_and(mask2.lon>-67 , mask2.lon<-49))
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
    # masklon= ( np.bitwise_and(mask2.lon>0 , mask2.lon<60))
    #masklat= np.bitwise_and( mask2.lat>15,mask2.lat<48)
    #mask2=mask2.where(~np.bitwise_and(masklon,masklat),np.nan)
    # mediterranean sea
    # masklon= ( np.bitwise_and(mask2.lon>-5.5 , mask2.lon<0))
    # masklat= np.bitwise_and( mask2.lat>15,mask2.lat<40)
    # mask2=mask2.where(~np.bitwise_and(masklon,masklat),np.nan)

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
def read_china_data():
    #with open('../../data/China_spare_capacity.csv', newline='') as csvfile:
    #    datareader = csv.reader(csvfile, delimiter= ',')
    #    for row in datareader:
    #        print(', '.join(row).replace(';','.'))
    ifile='../../data/China_spare_capacity.xlsx'
    ifile='../../data/China_cement_spare_capacity_optimum_scenario.xlsx'

    start_years=[2030, 2035, 2040]
    data={}
    #Open data
    for start_year in start_years:
        data[start_year]=pd.read_excel(ifile, sheet_name=str(start_year), index_col=0)
    print(data)
    return data

def read_eu_us_data(sheet):
    # read data from eu+us excel sheet for spacer capacity for ocean liming
    # name of the sheet as input
    ifile='../../data/Europe_and_US_spare_capacities_for_ocean_liming.xlsx'

    #Open data
    data=pd.read_excel(ifile, sheet_name=sheet, index_col=0)

    # new column names
    data.rename(columns = {'Mt/yr':'2030 low', 'Mt/yr.1':'2030 high',
                             'Unnamed: 3':'2035 low', 'Unnamed: 4':'2035 high',
                             'Unnamed: 5':'2040 low','Unnamed: 6':'2040 high'},
                  inplace=True)

    # add zeros for 2015-2030
    for i in range(2015,2030):
        data.loc[i]=[0,0,0,0,0,0]

    # insert zeros into the dataframe for 2030-2040 when using later start dates
    for i in range(2030,2035):
        data.at[i,'2035 high']=0.0
        data.at[i,'2035 low']=0.0
    for i in range(2030,2040):
        data.at[i,'2040 high']=0.0
        data.at[i,'2040 low']=0.0

    # sort it to get inserted years in the correct location
    data=data.sort_index()

    # drop extra line with no information
    data=data.drop(index=np.nan,axis=0)

    # change objects to floats
    for i in data.keys():
        data[i]=data[i].astype(float)
    # change to gt/yer
    data=data/1000.0
    return data
def read_eu_us_co2penalty():
    # read data from eu+us excel sheet for spacer capacity for ocean liming
    # name of the sheet as input
    ifile='../../data/Europe_and_US_spare_capacities_for_ocean_liming.xlsx'
    sheet='Additonal Lime Production for O'
    #Open data
    data=pd.read_excel(ifile, sheet_name=sheet, index_col=0)
    print(data)
    return data
eu_us_co2penalty_pct=read_eu_us_co2penalty()
usdata=read_eu_us_data('USA (lime and cement)')
eudata=read_eu_us_data('Europe (lime and cement)')
cementdata=read_china_data()
fig,axes=plt.subplots(1,2, figsize=(10,10))
for i, ax in enumerate(fig.axes):
    for k in range(0,3):
        spare_capacity_column=usdata.columns[2*k+i]
        usdata[spare_capacity_column].plot(ax=ax)
ax.legend(loc='upper center', bbox_to_anchor=(-0.1, -0.1),
           ncol=3, frameon=False)
# Quick plots
fig, axes = plt.subplots(2,2, figsize=(10,10))
start_years=[2030, 2035, 2040]
for i, ax in enumerate(fig.axes):
    print(i)
    spare_capacity_column=cementdata[2035].columns[i]
    for start_year in start_years:
        cementdata[start_year][spare_capacity_column].plot(ax=ax, label=start_year)
    ax.set_title(spare_capacity_column, fontsize=7)
ax.legend(loc='upper center', bbox_to_anchor=(-0.1, -0.1),
           ncol=3, frameon=False)
for ycase in start_years:
    for i in range(2015,2030):
        cementdata[ycase].loc[i]=[0,0,0,0]
    cementdata[ycase]=cementdata[ycase].sort_index()
print(cementdata)


# molar mass of Ca(OH)2
n_ca_oh_2=74.1 #g/mol
# molar mass of Ca(OH)2
n_ca_o=56.1 #g/mol
molar_mass={}
molar_mass['CaO']=n_ca_o
molar_mass['Ca(OH)2']=n_ca_oh_2


# hack it now
n_ca_oh_2=n_ca_o

# Evolution of excess lime
china_lime=np.array([0,0,0,0.14,0.267718,0.413041,0.571405,0.738835,0.91324,
                     1.09045,1.266208,1.438202,1.606129,1.76801,1.921216,2.0641,2.197527,2.323329])*1e6*1e3*1e3 #g/year
plt.figure()

# Liming year dimension
lime_year=np.linspace(2015,2100,china_lime.size)
# liming monthly dimension
all_years=np.linspace(2015,2100,(2100-2015+1)*12)

#dates in file start from 2030, but we need dates from 2015 for simulations
months2015_2030=np.linspace(2015,2030,(2030-2015)*12+1) + 1/24 # add half a month to put month descriptor middle month
print(months2015_2030[:-1])

all_months=np.linspace(start_years[0],2100,(2100-start_years[0]+1)*12) + 1/24 # add half a month to put month descriptor middle month
data_year_months=np.linspace(start_years[0],2100,(2100-start_years[0]+1)*12) + 1/24 # add half a month to put month descriptor middle month

#print(all_months)
#test=np.concatenate((months2015_2030[:-1],all_months))
all_months=np.concatenate((months2015_2030[:-1],all_months))
#print(test)
# interpolate to monthly
#china_lime_yearly=np.interp(all_years, lime_year, china_lime)
# add 0.5 to x-values to create annual sum correctly
#spare_capacity_column='Spare capacity'
spare_capacity_header='Total spare capacity (Gt)'
co2_penalty_header='Carbon footprint or carbon penatly (tCO2 per t of cement)'

def year2mon(inputdata,header,start_year=2030,linterp=False):
    print(inputdata[header])
    monlen=[31,28,31,30,31,31,30,31,30,31,30,31]
    monlen_leap=[31,29,31,30,31,31,30,31,30,31,30,31]
    # Liming year dimension
    data_year_months=np.linspace(start_year,2100,(2100-start_year+1)*12) + 1/24 # add half a month to put month descriptor middle month
    if (linterp):
        china_lime_yearly=np.interp(data_year_months, inputdata.index.values ,inputdata[header].values )
    else:
        china_lime_yearly=np.repeat(inputdata[header].values/12,12)
        j=0
        for iyear in range(2015,2101):
            if calendar.isleap(iyear):
                yearlen=366
                for imon in range(12):
                    test=inputdata[header][iyear]/yearlen*monlen_leap[imon]
                    #test=cementdata[start_years[0]][spare_capacity_header][iyear]/yearlen*monlen_leap[imon]
                    china_lime_yearly[j]=test
                    j=j+1
                    print(iyear,test)
            else:
                yearlen=365
                for imon in range(12):
                    test=inputdata[header][iyear]/yearlen*monlen[imon]
                    #test=cementdata[start_years[0]][spare_capacity_header][iyear]/yearlen*monlen[imon]
                    china_lime_yearly[j]=test
                    j=j+1
                    print(iyear,test)

    return china_lime_yearly

print( cementdata[start_years[0]][spare_capacity_header])
#china_lime_yearly=np.interp(all_months, cementdata[start_years[0]].index.values ,cementdata[start_years[0]][spare_capacity_column].values, )
china_lime_yearly=np.interp(data_year_months, cementdata[start_years[0]].index.values ,cementdata[start_years[0]][spare_capacity_header].values/12, ) # divide annual values by 12
#test=year2mon(cementdata[start_years[0]],spare_capacity_header,2015,False)
us_lime_yearly=year2mon(usdata,'2030 high',2015,False)
eu_lime_yearly=year2mon(eudata,'2030 high',2015,False)
print (china_lime_yearly)
#print (test-china_lime_yearly)
print (us_lime_yearly)
zero_months=np.zeros(len(months2015_2030))
#china_lime_yearly=np.interp(all_months, cementdata[start_years[0]].index.values ,cementdata[start_years[0]]['Spare capacity'].values, )
print(china_lime_yearly)
print(sum(china_lime_yearly[:13])/12,cementdata[start_years[0]][spare_capacity_header].values[15:19])
#print(china_lime_yearly2)
print(type(zero_months),type(china_lime_yearly))

# add 2015-2029 zero values in the beginning
china_lime_yearly=np.concatenate((zero_months[:-1],china_lime_yearly))
print(len(np.concatenate((zero_months[:-1],china_lime_yearly))))

china_lime_yearly_interp=china_lime_yearly.copy()
china_lime_yearly_step=np.repeat(cementdata[start_years[0]][spare_capacity_header].values/12,12)
monlen=     [31,28,31,30,31,31,30,31,30,31,30,31]
monlen_leap=[31,29,31,30,31,31,30,31,30,31,30,31]
j=0
for iyear in range(2015,2101):
    if calendar.isleap(iyear):
        yearlen=366
        checki=0
        for imon in range(12):
            print(imon)
            test=cementdata[start_years[0]][spare_capacity_header][iyear]/yearlen*monlen_leap[imon]
            china_lime_yearly_step[j]=test
            j=j+1
            print(iyear,test)
            checki=checki+test
        print(checki)
    else:
        yearlen=365
        checki=0
        for imon in range(12):
            test=cementdata[start_years[0]][spare_capacity_header][iyear]/yearlen*monlen[imon]
            china_lime_yearly_step[j]=test
            j=j+1
            print(iyear,test)
            checki=checki+test
        print(checki,cementdata[start_years[0]][spare_capacity_header][iyear])
china_lime_yearly=china_lime_yearly_step.copy()
plt.plot(china_lime_yearly_step[-50:])
#plt.show()

print(china_lime_yearly[:200],cementdata[start_years[0]][spare_capacity_header].values)
#penalty_co2_yearly=np.interp(all_months, cementdata[start_years[0]].index.values ,cementdata[start_years[0]][co2_penalty_header].values, )

plt.plot((china_lime_yearly_step+eu_lime_yearly+us_lime_yearly))


plt.plot(china_lime_yearly_step)
plt.plot(china_lime_yearly_interp)
plt.plot(china_lime_yearly)
penalty_co2_yearly=np.interp(data_year_months, cementdata[start_years[0]].index.values ,cementdata[start_years[0]][co2_penalty_header].values)
penalty_co2_yearly=np.concatenate((zero_months[:-1],penalty_co2_yearly))

eu_us_penalty_co2=read_eu_us_co2penalty()
eu_us_penalty_co2_yearly=np.interp(data_year_months, eu_us_penalty_co2.index.values ,eu_us_penalty_co2['Carbon penalty (%)'].values)
eu_us_penalty_co2_yearly=np.concatenate((zero_months[:-1],eu_us_penalty_co2_yearly))
plt.figure()
plt.plot(all_months,eu_us_penalty_co2_yearly,lw=3)
#plt.figure()
plt.plot(eu_us_penalty_co2.index.values,eu_us_penalty_co2['Carbon penalty (%)'].values,'r',lw=1)

plt.figure()
plt.plot(all_months,china_lime_yearly)
plt.title('testi')
#plt.figure()
#plt.plot(all_months,test)
plt.title('testi_repeat')
plt.figure()
plt.plot(all_months)
plt.title('testi_repeat')
#print(all_months[160:200])
print(china_lime_yearly[-40:])
#print(test[-40:])
#china_lime_yearly=test
# For now copy china to EU and USA
#eu_lime_yearly=china_lime_yearly.copy()
#us_lime_yearly=china_lime_yearly.copy()
# check
print(china_lime_yearly)
plt.figure()
plt.plot(all_months,china_lime_yearly,'o-')
plt.figure()
plt.plot(lime_year,china_lime,'o-')
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

china_mask_co2=regionmask.defined_regions.natural_earth_v5_0_0.countries_50['China']
us_mask_co2=regionmask.defined_regions.natural_earth_v5_0_0.countries_50['US']
eu_mask_co2=regionmask.defined_regions.natural_earth_v5_0_0.countries_50['China']
print(china_mask_co2)
a=regionmask.defined_regions.natural_earth_v5_0_0.countries_50
b=a.mask(lons,lats)
china_index=regionmask.defined_regions.natural_earth_v5_0_0.countries_50.map_keys('China')
US_index=regionmask.defined_regions.natural_earth_v5_0_0.countries_50.map_keys('US')
EU_index=regionmask.defined_regions.natural_earth_v5_0_0.countries_50.map_keys(['France','Spain','Netherlands','Belgium','Germany','Norway','Sweden',
                  'Italy','Greece','Estonia','Lithuania','Latvia','Denmark','Finland','Poland','Romania','Bulgaria','Croatia','Malta',
                  'United Kingdom','Ireland','Portugal','Slovenia'])
temp=grid_nc.copy()
#temp.m2.data[:,:]=np.nan
print(temp.m2)
plt.figure()
temp.m2.plot()
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
#        data2.data=np.where((data.data==float(i)),10,data2.data)
plt.figure()
EU_m2.plot()
plt.figure()
china_m2.plot()
## exclusive economical zones as geopandas
eez=gpd.read_file('../../data//World_EEZ_v11_20191118_HR_0_360/eez_v11_0_360.shp')
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
                  'Italy','Greece','Estonia','Lithuania','Latvia','Denmark','Finland','Poland','Romania','Bulgaria','Croatia','Malta',
                  'United Kingdom of Great Britain and Northern Ireland','Ireland','Portugal','Slovenia ']
Not_EU=['Turkiye','Israel']
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
                                print (eez['GEONAME'][i])
                                for iii in eez['GEONAME'][i]:
                                    if 'Overlapping claim' in iii:
                                        print(iii)
                                if 'Overlapping claim South China Sea' in eez['GEONAME'][i]:
                                    continue
                                #  print ('hep',eez['GEONAME'])
                                #if eez['UN_SOV1'][eez['UN_SOV1']==
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

time_dates = np.arange(np.datetime64("2015-01-01"), np.datetime64("2101-01-01"),  np.timedelta64(1, 'M'), dtype='datetime64[M]')
print(len(time_dates))
print(len(all_months))
print(time_dates[:40])
print(all_months[:40])
print(time_dates[-40:])
print(all_months[-40:])

#print(excluded_eez)
# prepare data for mask
#data=np.zeros((nyears*12,360,720))
data=maskeez.copy()#expand_dims({'time':time_dates.size})
#data=np.where(data.data>0,1,0)
#for i in excluded_eez:
# create copy with same dimensions
data2=data.copy()
data2.data[:,:]=np.nan
# penalty in CO2
china_penalty=china_m2.copy()
china_penalty.data[:,:]=np.nan
US_penalty=US_m2.copy()
US_penalty.data[:,:]=np.nan
EU_penalty=EU_m2.copy()
EU_penalty[:,:]=np.nan
#data2.data[:,:,:]=np.nan
#print (data2)
#for i in range(0,150):
for i in included_eez_dict:
    print ('incl',i)
    if included_eez_dict[i]=='China':
        print('Computing China')
        data2.data=np.where((data.data==float(i)),10,data2.data)
        china_penalty.data=np.where((~np.isnan(china_m2.data)),10,china_penalty.data)

    elif included_eez_dict[i]=='United States of America':
        print('Computing USA')
        data2.data=np.where((data.data==float(i)),20,data2.data)
        US_penalty.data=np.where((~np.isnan(US_m2.data)),20,US_penalty.data)
    else:
        print('Computing EU')
        data2.data=np.where((data.data==float(i)),30,data2.data)
        EU_penalty.data=np.where((~np.isnan(EU_m2.data)),30,EU_penalty.data)
data2.data=np.where(~np.isnan(data2.data),data2.data,np.nan)
#data_penalty.data=np.where(~np.isnan(data_penalty.data),data_penalty.data,np.nan)
print(data2)
# calculate the area of sea used for different parts
area_china=grid_nc.where(data2.data==10).sum()
area_us=grid_nc.where(data2.data==20).sum()
area_eu=grid_nc.where(data2.data==30).sum()
print(f"CHINA: {area_china:.2f}")
print(f"US: {area_us:.2f}")
print(f"EU: {area_eu:.2f}")
plt.figure()
china_penalty.plot()
plt.figure()
US_penalty.plot()
plt.figure()
EU_penalty.plot()
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
land_area_china=china_m2.sum()
land_area_us=US_m2.sum()
land_area_eu=EU_m2.sum()
# unit change
yearly_kg_2_mol_s=1/n_ca_oh_2/(365*3600*24)*100
monthly_gt_2_mol_s=np.zeros(12)
monthly_gt_2_mol_s_leap=np.zeros(12)
monlen=[31,28,31,30,31,31,30,31,30,31,30,31]
monlen_leap=[31,29,31,30,31,31,30,31,30,31,30,31]
gt_2_g=1e6*1e9
for i in range(12):
    #monthly_gt_2_mol_s[i]=1/n_ca_oh_2/(monlen[i]/365*3600*24)*gt_2_g
    #monthly_gt_2_mol_s_leap[i]=1/n_ca_oh_2/(monlen[i]/366*3600*24)*gt_2_g
    monthly_gt_2_mol_s[i]=1/n_ca_oh_2/(monlen[i]*3600*24)*gt_2_g
    monthly_gt_2_mol_s_leap[i]=1/n_ca_oh_2/(monlen_leap[i]*3600*24)*gt_2_g
    monthly_gt_2_mol_s[i]=1/n_ca_oh_2/(monlen[i]*3600*24)*gt_2_g
    monthly_gt_2_mol_s_leap[i]=1/n_ca_oh_2/(monlen_leap[i]*3600*24)*gt_2_g
print(monthly_gt_2_mol_s)

#penalty kgs-1
penalty_t_2_kg_s=1000/(365*3600*24)

#data2=data2.expand_dims({'time':time_dates.size}).copy()
# expand time dimension
fig,ax=plt.subplots(ncols=1,subplot_kw={'projection': ccrs.Robinson()})
#mask2.lime_mask[:,:].plot(transform=ccrs.PlateCarree(),add_labels=True)
ax.coastlines(color='white')

data3=data2.copy()#.plot(ax=ax)
data3.data=np.where(data2.data[:,:]>0,1,0)#.plot(ax=ax)
data3.plot.pcolormesh(ax=ax,transform=ccrs.PlateCarree(), add_colorbar=False)
ax.coastlines(color='white')
plt.savefig('maski.png')
plt.show()
ff
data2=data2.expand_dims({'time':time_dates}).copy()
china_penalty=china_penalty.expand_dims({'time':time_dates}).copy()
US_penalty=US_penalty.expand_dims({'time':time_dates}).copy()
EU_penalty=EU_penalty.expand_dims({'time':time_dates}).copy()

# free up memory
#del eez
#del data
print(len(us_lime_yearly))

# Change to liming data as mol m-2 s-1
for i in range(data2.data.shape[0]):

    print(i,time_dates[i],china_lime_yearly[i], penalty_co2_yearly[i])
    imonth=time_dates[i].astype(object).month -1
    iyear=time_dates[i].astype(object).year
    isleap=calendar.isleap(iyear)
    if  calendar.isleap(iyear):
        print('leap',time_dates[i],china_lime_yearly[i]/area_china*monthly_gt_2_mol_s_leap[imonth])
    else:
        print(time_dates[i],china_lime_yearly[i]/area_china*monthly_gt_2_mol_s[imonth])

    print(imonth)
    #input()
    if china_lime_yearly[i]<1e-40:
        data2.data[i,:,:]=0.0
    else:
        print(imonth,i,china_lime_yearly[i]/area_china*monthly_gt_2_mol_s[imonth])
        if calendar.isleap(iyear):
            data2.data[i,:,:]=np.where((data2.data[i,:,:]==10),china_lime_yearly[i]/area_china*monthly_gt_2_mol_s_leap[imonth],data2.data[i,:,:])
        else:
            data2.data[i,:,:]=np.where((data2.data[i,:,:]==10),china_lime_yearly[i]/area_china*monthly_gt_2_mol_s[imonth],data2.data[i,:,:])
    if us_lime_yearly[i]>1e-40:
    #        data2.data[i,:,:]=0.0
    #else:
        if calendar.isleap(iyear):
            data2.data[i,:,:]=np.where((data2.data[i,:,:]==20),us_lime_yearly[i]/area_us*monthly_gt_2_mol_s_leap[imonth],data2.data[i,:,:])
        else:
            data2.data[i,:,:]=np.where((data2.data[i,:,:]==20),us_lime_yearly[i]/area_us*monthly_gt_2_mol_s[imonth],data2.data[i,:,:])
    if eu_lime_yearly[i]>1e-40:
    #    data2.data[i,:,:]=0.0
    #else:
        if calendar.isleap(iyear):
            data2.data[i,:,:]=np.where((data2.data[i,:,:]==30),eu_lime_yearly[i]/area_eu*monthly_gt_2_mol_s_leap[imonth],data2.data[i,:,:])
        else:
            data2.data[i,:,:]=np.where((data2.data[i,:,:]==30),eu_lime_yearly[i]/area_eu*monthly_gt_2_mol_s[imonth],data2.data[i,:,:])

# =============================================================================
#     #CO2 pnealty
#     if penalty_co2_yearly[i]<1e-40:
#         china_penalty.data[i,:,:]=0.0
#     else:
#         china_penalty.data[i,:,:]=np.where((china_penalty.data[i,:,:]==10),penalty_co2_yearly[i]*china_lime_yearly[i]/land_area_china*penalty_t_2_kg_s,china_penalty.data[i,:,:])
#     if eu_us_penalty_co2_yearly[i]>1e-40:
#     #    china_penalty.data[i,:,:]=0.0
#     #    US_penalty.data[i,:,:]=0.0
#     #else:
#         china_penalty.data[i,:,:]=np.where((US_penalty.data[i,:,:]==20),eu_us_penalty_co2_yearly[i]*us_lime_yearly[i]/land_area_us*penalty_t_2_kg_s,china_penalty.data[i,:,:])
#         US_penalty.data[i,:,:]=np.where((US_penalty.data[i,:,:]==20),eu_us_penalty_co2_yearly[i]*us_lime_yearly[i]/land_area_us*penalty_t_2_kg_s,china_penalty.data[i,:,:])
#     if eu_us_penalty_co2_yearly[i]>1e-40:
#     #    china_penalty.data[i,:,:]=0.0
#     #    EU_penalty.data[i,:,:]=0.0
#     #else:
#         china_penalty.data[i,:,:]=np.where((EU_penalty.data[i,:,:]==30),eu_us_penalty_co2_yearly[i]*eu_lime_yearly[i]/land_area_eu*penalty_t_2_kg_s,china_penalty.data[i,:,:])
#         EU_penalty.data[i,:,:]=np.where((EU_penalty.data[i,:,:]==30),eu_us_penalty_co2_yearly[i]*eu_lime_yearly[i]/land_area_eu*penalty_t_2_kg_s,china_penalty.data[i,:,:])
#     #data_penalty.data[i,:,:]=np.where((data2.data[i,:,:]==20),us_lime_yearly[i]/area_us*penalty_t_2_kg_s,data_penalty.data[i,:,:])
#     #data_penalty.data[i,:,:]=np.where((data2.data[i,:,:]==30),eu_lime_yearly[i]/area_eu*penalty_t_2_kg_s,data_penalty.data[i,:,:])
# =============================================================================
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
               Based on the world exclusive economic zones (marineregions.org) and work done by Spyros Foteinis and Phil Renforth \
                   for annual lime production numbers.',
               author='Tommi Bergman (FMI)'
        )
)
data_penalty=china_penalty+US_penalty+EU_penalty
co2penalty=xr.Dataset(
    data_vars=dict(
        co2penalty=(['time','lat','lon'],china_penalty.data,{'units':'kg m-2 s-1'})
        #gridarea=(['lat','lon'],grid*lat_size,{'units':'m-2'})
        #mask=(['time','lat','lon'],data2.data,{'units':'country number'})
        #data_vars=(['lon','lat'],data,{'units':'mol/s'})
        ),
    coords=dict(
        time=('time',time_dates),
        lon=('lon',lons,{'units':'degrees_east'}),
        lat=('lat',lats,{'units':'degrees_north'})
        ),
    attrs=dict(description='CO2 penalty for CCS for ocean alkalinization for WP4 of OceanNETs project.\
               Based on the world exclusive economic zones (marineregions.org) and work done by Spyros Foteinis and Phil Renforth \
                   for annual lime production numbers.',
               author='Tommi Bergman (FMI)'
        )
)
# for i in range(100,150):
#     print (i)
#     mask.mask.data=np.where(np.abs(mask.mask.data-i)<0.5,1,mask.mask.data)
# write to disk, this includes
mask.to_netcdf('../../data/lime_mask_v2.cao.nc')
co2penalty.to_netcdf('../../data/co2penalty.nc')
mask.plot()
plt.figure()
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
mask2.to_netcdf('../../data/limemask_v2.nc',engine='scipy')

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
orcam2=xr.open_dataset('../../data/areacello.nc')
plt.figure()
plt.plot(mask_orca1.time,(mask_orca1*orcam2.areacello).lime_mask.sum(dim=('x','y'))*n_ca_oh_2*365*3600*24,label='orca1')
plt.plot(mask2.time,(mask2*grid_nc).lime_mask.sum(dim=('lon','lat'))*n_ca_oh_2*365*3600*24,label='lonlat')
plt.legend()
plt.title('diff')

# calculate the error from remgridding
error=((mask2*grid_nc).lime_mask.sum(dim=('lon','lat'))*n_ca_oh_2*365*3600*24/((mask_orca1*orcam2.areacello).lime_mask.sum(dim=('x','y'))*n_ca_oh_2*365*3600*24)).isel(time=range(1020,1031)).mean()
error=((mask2*grid_nc).lime_mask.sum(dim=('lon','lat'))*n_ca_oh_2*365*3600*24/((mask_orca1*orcam2.areacello).lime_mask.sum(dim=('x','y'))*n_ca_oh_2*365*3600*24)).mean()
# correct the grid calues with the error
mask_orca1=mask_orca1*error
plt.figure()
plt.plot(mask_orca1.time,(mask_orca1*orcam2.areacello).lime_mask.sum(dim=('x','y'))*n_ca_oh_2*365*3600*24,label='orca')
#plt.plot(mask2.time,(mask2*grid_nc).lime_mask.sum(dim=('lon','lat'))*n_ca_oh_2*365*3600*24,label='latlon')
plt.legend()
plt.title('diff')
plt.figure()
plt.plot((mask_orca1*orcam2.areacello).lime_mask.sum(dim=('x','y'))*n_ca_oh_2*365*3600*24/((mask2*grid_nc).lime_mask.sum(dim=('lon','lat'))*n_ca_oh_2*365*3600*24))

print(((mask2*grid_nc).lime_mask.sum(dim=('lon','lat'))*n_ca_oh_2*365*3600*24/((mask_orca1*orcam2.areacello).lime_mask.sum(dim=('x','y'))*n_ca_oh_2*365*3600*24)).isel(time=range(1020,1031)).mean())

mask_orca1.to_netcdf('../../data/limemask_v2.orca1.cao.nc',engine='scipy')
plt.show()#