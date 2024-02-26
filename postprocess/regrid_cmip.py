#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 
"""
regrid cmip output to 1x1 degrees
"""
import numpy as np
import xarray as xr
from itertools import product
import pandas as pd
#import matplotlib.pyplot as plt
import os
#import cartopy.crs as ccrs
#import cartopy.feature as cfeature
#import regionmask
#import csv
from datetime import datetime, timedelta
#import calendar
#from mask_aux import create_eez_mask, gridarea, masking_boxes, read_deployment_data, year2mon,  select_areas
#import ESMpy
#import xesmf as xe
import glob




vars={
'co2':'xco2',
'fgco2':'fgco2',
'nbp':'nbp',
'alkalinization':'fgoae',
'fco2removal2D':'fgdcr',
'spco2':'spco2',
'dissicos':'dissicos',
'talkos':'talkos',
'phos':'pH',
#'omegaa',
'omegac':'omegac',
'gpp':'gpp',
'cVeg':'cveg',
'cSoil':'csoil',
'tas':'tas'
}
 
n_c=12.011
unit_conversions = {'cSoil':1000/n_c,  #kgm-2 -> molCm-2
                    'cVeg':1000/n_c,   #kgm-2 -> molCm-2
                    'gpp':1000/n_c,    #kgm-2s-1 -> molCm-2s-1
                    'nbp':1000/n_c,    #kgm-2s-1 -> molCm-2s-1
                    'fgco2':1000/n_c,    #kgm-2s-1 -> molCm-2s-1
                    'dissicos':1,    #molm-3
                    'talkos':1,    #molm-3
                    'spco2':1e6/101325,    # Pa -> uAtm
                    'co2':1e6,    # mol/mol -> umol/mol or ppm  
                    'phos':1, # no conversion unit of 1
                    'omeagac':1, #no conversion unit of 1
                    'fco2removal2D':1000/n_c, # XXXX
                    'alkalinization':1, # xxxxx
                    'tas':-273.15,
                    'tos':-273.15
                    }  
units = {'csoil':'molC m-2',
                    'cveg':'molC m-2',
                    'gpp':'molC m-2 s-1',
                    'nbp': 'molC m-2 s-1',
                    'fgco2': 'molC m-2 s-1',
                    'dissicos': 'mol m-3',
                    'talkos':'mol m-3',
                    'spco2': 'uAtm',
                    'xco2': 'ppm',
                    'pH':'-',
                    'omeagac':'-',
                    'fgdcr':'molC m-2 s-1',
                    'fgoae': 'mol m-2 s-1',
                    'tas':'degC'}  


expnames={
'C5AF':'esm-ssp534-over',
'D30H':'esm-ssp534-over-2030high',
'D40H':'esm-ssp534-over-2040high',
'E40H':'esm-ssp534-over-2040high-term',
'C40H':'esm-ssp534-over-2040high-dcr'
}

model_id='EC-Earth'
now = datetime.now() # current date and time
year = now.strftime("%Y")
print("year:", year)

month = now.strftime("%m")
print("month:", month)

day = now.strftime("%d")
print("day:", day)
hour = now.strftime("%H:%M")
print("hour:", hour)

version_id=f'v{year}{month}{day}'

time_range='2015-2100'
#syear='2050'
#exp='D30H'
#vars=['fgco2']
expit=['C40H']
#expit=['E40H', 'D40H']
#expit=['C5AF','D30H']
for exp in expit:#expnames.keys():
    experiment_id=expnames[exp]
    path=f'/Volumes/ONETs-SSD/{exp}/1x1/'
    outpath=f'/Volumes/ONETs-SSD/{experiment_id}/2D'
    if not os.path.exists(outpath):
        # Create the folder if it doesn't exist
        os.makedirs(outpath)
        print(f"Folder '{outpath}' created.")
    else:
        print(f"Folder '{outpath}' already exists.")
        
    for var in vars.keys():
        for syear in range (2015,2101): #'2050'
            outputname=f'{vars[var]}_{model_id}_{experiment_id}_{time_range}_{version_id}.{syear}nc'
            print(f"{path}/{vars[var]}_*{syear}*nc")
            files=glob.glob(f"{path}/{var}_*mon*{syear}*nc")
            print(files)
            for fname in files:
                ds=xr.open_dataset(fname,engine='h5netcdf')
                if var=='tas' or var=='tos':
                    ds[var].data=ds[var].data + unit_conversions[var]
                else:
                    ds[var].data=ds[var].data * unit_conversions[var]
                orig_unit=ds[var].attrs['units']
                ds.attrs={'Source':'EC-Earth3-CC', 'Author':'Tommi Bergman','Institute':'Finnish Meteorological Institute',  'project':'OceanNETs','Experiment':experiment_id, 'created':f'{year}-{month}-{day} {hour}',
                        'original_variable':f'{var}','original unit':orig_unit, 'unit_conversion':f'{unit_conversions[var]} (*area when needed)'}
                ds[var].attrs['units']=units[vars[var]]
                if var=='co2':
                    print(ds)
                    ds[var]=ds[var].isel(lev=1)
                    print(ds[var])#ds[var]=ds[var].mean(dim='time')
                #attr_dict=ds[var].attrs
                # Only change in reference time is needed. But use encoding instead of attrs
                #https://github.com/pydata/xarray/issues/3739
                # Change the reference date to "days since 2022-01-01"
                new_reference_time = pd.to_datetime('2015-01-01')
                # Update the units of the time coordinate
                ds['time'].encoding['units'] = f"days since {new_reference_time}"
                #ds['time_bnds'].encoding['units'] = f"days since {new_reference_time}"
                # Display the updated time coordinate
                print("\nUpdated Time Coordinate:")
                print(ds)
                ds=ds.rename({var:vars[var]})
                print (ds)
                
                ds.to_netcdf(outpath+'/'+outputname)
                # if exp == 'C5AF':
                #     print('copy initila files for other exp')
                #     for exp1 in expnames.keys():
                #         print(exp1)
                #         print(exp1)
                #         if exp == 'E30H' and int(syear) < 2081:
                #             outpath=f'/Volumes/ONETs-SSD/{expnames[exp1]}/2D'
                #             #outpath=f'/Volumes/ONETs-SSD/{expnames[exp1]}'
                #             #outputname=f'{vars[var]}_{model_id}_{experiment_id}_{time_range}_{version_id}.{syear}nc'
                #             outputname=f'{vars[var]}_{model_id}_{expnames[exp1]}_{time_range}_{version_id}.{syear}.nc'
                #             ds.attrs={'Source':'EC-Earth3-CC', 'Author':'Tommi Bergman','Institute':'Finnish Meteorological Institute', 'Project':'OceanNETs','Experiment':f'{expnames[exp1]}', 'Created':f'{year}-{month}-{day} {hour}',
                #                     'original_variable':f'{var}','original unit':orig_unit, 'unit_conversion':f'{unit_conversions[var]}'}
                #             ds.to_netcdf(f'{outpath}/{outputname}',engine='h5netcdf')
                #         elif  (exp == 'D40H' or exp == 'C40H') and int(syear) < 2040:
                #             print(exp1)
                #             outpath=f'/Volumes/ONETs-SSD/{expnames[exp1]}/2D'
    
                #             #outpath=f'/Volumes/ONETs-SSD/{expnames[exp1]}'
                #             print(outpath)
                #             outputname=f'{vars[var]}_{model_id}_{expnames[exp1]}_{time_range}_{version_id}.{syear}.nc'
                #             #outputname=f'{vars[var]}_{model_id}_{experiment_id}_{time_range}_{version_id}.{syear}nc'
                #             print(outputname)
                #             ds.attrs={'Source':'EC-Earth3-CC', 'Author':'Tommi Bergman','Institute':'Finnish Meteorological Institute', 'Project':'OceanNETs','Experiment':f'{expnames[exp1]}', 'Created':f'{year}-{month}-{day} {hour}',
                #                     'original_variable':f'{var}','original unit':orig_unit, 'unit_conversion':f'{unit_conversions[var]}'}
                #             ds.to_netcdf(f'{outpath}/{outputname}',engine='h5netcdf')

                #         elif exp == 'D30H'  and int(syear) < 2030:
        
                #             print(exp1)
                #             outpath=f'/Volumes/ONETs-SSD/{expnames[exp1]}/2D'
    
                #             #outpath=f'/Volumes/ONETs-SSD/{expnames[exp1]}'
                #             print(outpath)
                #             outputname=f'{vars[var]}_{model_id}_{expnames[exp1]}_{time_range}_{version_id}.{syear}.nc'
                #             #outputname=f'{vars[var]}_{model_id}_{experiment_id}_{time_range}_{version_id}.{syear}nc'
                #             print(outputname)
                #             ds.attrs={'Source':'EC-Earth3-CC', 'Author':'Tommi Bergman','Institute':'Finnish Meteorological Institute', 'Project':'OceanNETs','Experiment':f'{expnames[exp1]}', 'Created':f'{year}-{month}-{day} {hour}',
                #                     'original_variable':f'{var}','original unit':orig_unit, 'unit_conversion':f'{unit_conversions[var]}'}
                #             ds.to_netcdf(f'{outpath}/{outputname}',engine='h5netcdf')

                # else:
                #     #ds.to_netcdf(f'{outpath}/{vars[var]}_testi_{syear}.nc',engine='h5netcdf')
                #     ds.to_netcdf(f'{outpath}/{outputname}',engine='h5netcdf')
                ds.close()
# ds_in = xe.util.grid_2d(
#     0.5, 359.5, 1.0, -89.5, 89.5, 1.0  # longitude range and resolution
# )  # latitude range and resolution

# # regrid to ocean grid
# # regrid_to_model(mask2)
# # reference grid
# #ingrid=xr.open_dataset('/Users/bergmant/Documents/projects/H2020-OceanNET//data/bathy_meter.nc')
# #ingrid=xr.open_dataset('/Users/bergmant/Documents/projects/H2020-OceanNET/scripts/postprocess/areacello_Ofx_EC-Earth3-CC_esm-pi-cdr-pulse_r1i1p1f1_gr.nc')
# # makesure lon and lat is proper
# #ingrid.rename({'nav_lon':'lon', 'nav_lat':'lat'})
# # regridder creation
# regridder = xe.Regridder(ds, ds_in, 'nearest_s2d', periodic=True) #since this is global we need to pass periodic
# # regrid
# ds_1x1=regridder(ds)
# #print(list(outdict.keys())[0])
# ds_1x1.to_netcdf(f'fgco2_test.nc',engine='h5netcdf')
