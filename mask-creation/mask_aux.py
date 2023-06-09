#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 09:02:50 2023

@author: Antti-Ilari Partanen (antti-ilari.partanen@fmi.fi)
"""
import numpy as np
import xarray as xr
import pandas as pd
import geopandas as gpd
import regionmask
from itertools import product
import calendar
import csv

def read_deployment_data():
    ifile='../data/Case study 1 - deployment rates.xlsx'
    
    start_years=['2030', '2035', '2040']
    volumes=['low', 'high']
    scenarios=list()

    for start_year, volume in product(start_years,volumes):
        scenarios.append(start_year+'-'+volume)
        
        
    data=dict()    

    for scenario in scenarios:
        data[scenario]=pd.read_excel(ifile,sheet_name=scenario,index_col=0)
        
        #Replace NaNs with 0's
        data[scenario]=data[scenario].fillna(0.)
        
    return data


def year2mon(inputdata):
    '''
    Divide yearly to months so that each month has same daily deployment rate
    (longer months have higher monthly rate). The time coordinate is also
    expanded to cover months 1/2015 – 12/2100.

    Parameters
    ----------
    inputdata : TYPE
        DESCRIPTION.

    Returns
    -------
    data_monthly : TYPE
        DESCRIPTION.

    '''
    monlen=[31,28,31,30,31,31,30,31,30,31,30,31]
    monlen_leap=[31,29,31,30,31,31,30,31,30,31,30,31]
    # Liming year dimension
    # data_year_months=np.linspace(start_year,2100,(2100-start_year+1)*12) + 1/24 # add half a month to put month descriptor middle month
    all_months=np.linspace(2015,2100+11/12,(2100-2015+1)*12)+ 1/24
    
    data_monthly=pd.DataFrame(0, index=all_months, columns=inputdata.columns) # Unit: Gt CaO / month
    data_monthly_mol_s=pd.DataFrame(0, index=all_months, columns=inputdata.columns) # mol CaO / second

    
    gt_2_g=1e6*1e9 # Gt to g
    n_cao=56.1 #Molar mass of CaO [g/mol]
    for column in data_monthly.columns:
        data_monthly.loc[2030:2101, column]=np.repeat(inputdata[column].values/12,12)
        
    
    j=0
    for iyear in range(2015,2101):
        if calendar.isleap(iyear):
            yearlen=366
            monlen_year=monlen_leap
            
        else:
            yearlen=365
            monlen_year=monlen
        for imon in range(12):
            data_monthly[column].iloc[j]=data_monthly[column].iloc[j]/yearlen*monlen_year[imon]
            
            # Calculate the deployment rates in units of mol CaO / s
            if iyear >= 2030:
                data_monthly_mol_s[column].iloc[j]=inputdata.loc[iyear, column]/n_cao/(365*3600*24)*gt_2_g

            j=j+1
    
    return data_monthly, data_monthly_mol_s
    
   
def select_areas(eez):
   # Zones to be excluded from alkalinity simulations
    excluded_eez=[]
    ## Zones to be included
    included_eez=[]
    included_eez_dict={}
    included_region=[]
    included_country=['United States of America','China','France','Spain','Netherlands','Belgium','Germany','Norway','Sweden',
                      'Italy','Greece','Estonia','Lithuania','Latvia','Denmark','Finland','Poland','Romania','Bulgaria','Croatia','Malta',
                      'United Kingdom of Great Britain and Northern Ireland','Ireland','Portugal','Slovenia ']
    Not_EU=['Turkiye','Israel']
    with open('../data/UNSD — Methodology - country numbers for EEZ.csv') as csvfile:
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
   
    
    return included_eez_dict
    
    
  

def create_eez_mask():
    ## exclusive economical zones as geopandas
    eez=gpd.read_file('../data/World_EEZ_v11_20191118_HR_0_360/eez_v11_0_360.shp')
    print(eez.head())
    
    lons=np.linspace(-179.75,179.75,720)
    #lon_bounds=np.linspace(-180,180,721)
    lats=np.linspace(-89.75,89.75,360)
    
    # lon = np.arange(-180, 180)
    # lat = np.arange(-90, 90)
    # print(lons)
    # print(lon)
    # print(lats)
    # print(lat)
    #maskeez = regionmask.mask_geopandas(eez, lon, lat,wrap_lon=True)
    maskeez = regionmask.mask_geopandas(eez, lons, lats,wrap_lon=True)
    return eez, maskeez

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
    # plt.pcolormesh(mlons[1,:],mlats[:,1],grid_nc)
    # plt.colorbar()
    #plt.show()
    return grid_nc

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
