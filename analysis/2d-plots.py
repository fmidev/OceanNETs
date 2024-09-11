import numpy as np
import xarray as xr
from itertools import product
import pandas as pd
import matplotlib.pyplot as plt
import os
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import regionmask
import csv
from datetime import datetime, timedelta
import calendar
#import xesmf as xe
import glob 
import cftime

#fig,ax1=plt.subplots()
nrows=2
ncols=2
fig, axs = plt.subplots(nrows=nrows,ncols=ncols,
                        subplot_kw={'projection': ccrs.Robinson()},
                        figsize=(11,8.5))
fig2, axs2 = plt.subplots(nrows=nrows,ncols=ncols,
                        subplot_kw={'projection': ccrs.Robinson()},
                        figsize=(11,8.5))
#axs.flatten()
cbar_kwargs = {'orientation':'horizontal', 'shrink':0.9, 'aspect':40, 'label':'Ocean alkalinity addition regions'}
countrymask = xr.open_dataset('countrymask.1x1.nc')
countrymask['mask']=xr.where(countrymask['mask']>1,countrymask['mask'],np.nan)
countrymask['mask']=xr.where(countrymask['mask']>1,10,countrymask['mask'])
countrymask['mask'].plot(ax=axs[0,0], transform=ccrs.PlateCarree(),add_labels=True, cbar_kwargs=cbar_kwargs)
axs[0,0].coastlines()
countrymask['mask'].plot(ax=axs2[0,0], transform=ccrs.PlateCarree(),add_labels=True, cbar_kwargs=cbar_kwargs)
axs2[0,0].coastlines()
#plt.show()
ctrl='esm-ssp534-over'
datapath='/Volumes/ONETs-SSD/Oceannets-all-models/data/'
vars=['talkos','pH','omegaa']#,'omegac_reg']#,'fgoae_glob']

colors={'EC-Earth':'b','NorESM2-LM':'r','AWI-ESM':'g'}
lintypes={'esm-ssp534-over-2040high-dcr':'--','esm-ssp534-over-2040high':'-','esm-ssp534-over-2040high-term':':'}
models=['NorESM2-LM','AWI-ESM','EC-Earth','FOCI']
labels={'talkos':'Total Alkalinity at \n Ocean Surface','pH':'pH at Ocean Surface','omegaa':'OmegaA Aragonite Saturation State'}
experiments=['esm-ssp534-over-2040high']#,'esm-ssp534-over-2040high']
ctrl_data={}
exp_data={}
diff_data={}
fig.suptitle(experiments[0])
fig2.suptitle(experiments[0])
for var in vars:
    diff={}
    for mod in models:
        if mod=='FOCI':
            if var != 'pH' and var !='omegaa':
                continue
            filename=var+'_'+mod+'*.1x1.nc'
        else:
            filename=var+'_'+mod+'*'
        if var=='omegaa' and (mod=='EC-Earth' or mod=='AWI-ESM'):
            continue
            
        for d in (glob.glob(datapath+'/'+ctrl+'/'+filename)):
            data_a=xr.open_dataset(d)
            a=data_a.copy()
            if data_a.lon[-1]>181:
                b=a.assign_coords({"lon":(((a.lon + 180) % 360 ) -180 )}).sortby('lon')
            else:
                b=a
            ctrl_data[mod]=b
        exp=experiments[0]
        for d in (glob.glob(datapath+'/'+exp+'/'+filename)):
            print(d)
            data_a=xr.open_dataset(d)
            a=data_a.copy()
            if data_a.lon[-1]>181:
                b=a.assign_coords({"lon":(((a.lon + 180) % 360 ) -180 )}).sortby('lon')
            else:
                b=a
            exp_data[mod]=b
        if not d == None:
            diff[mod]=exp_data[mod]-ctrl_data[mod]
    diff_data[var] =diff
model_N=[]
N_plot=1
for var in vars:
    if N_plot==1:
        i=0
        j=1
    elif N_plot==2:
        i=1
        j=0
    elif N_plot==3:
        i=1
        j=1

    all_data=None
    all_N=0
    all_data=0
    for mod in models:
        print(var,mod)
        if mod=="FOCI" and (var != 'pH' and var != 'omegaa'):
            continue
        elif mod=='AWI-ESM' and var=='omegaa':
            continue
        elif mod=='EC-Earth' and var=='omegaa':
            continue
        
        print('----------------------------------------------------------------')
        if all_N==0:
            all_data=all_data+ diff_data[var][mod][var].isel(time=[-10,-9,-8,-7,-6,-5,-4,-3,-2,-1]).mean(dim="time").data
        else:
            all_data+=diff_data[var][mod][var].isel(time=[-10,-9,-8,-7,-6,-5,-4,-3,-2,-1]).mean(dim="time").data
        all_N+=1
    model_N.append(all_N)
    all_data=all_data/all_N
    meandata=xr.zeros_like(diff_data[var]['NorESM2-LM'][var].isel(time=[-10,-9,-8,-7,-6,-5,-4,-3,-2,-1]).mean(dim="time"))
    meandata.data=all_data 
    cbar_kwargs = {'orientation':'horizontal', 'shrink':0.9, 'aspect':40, 'label':labels[var]}
    meandata.plot(ax=axs[i,j], transform=ccrs.PlateCarree(),add_labels=True, cbar_kwargs=cbar_kwargs)
    axs[i,j].coastlines()
    if var=='pH' or var=='talkos':
        datamin=0
        datamax=0.2
    else:
        datamin=0
        datamax=2
    meandata.plot(ax=axs2[i,j], transform=ccrs.PlateCarree(),add_labels=True,  cmap='Reds',vmin=datamin, vmax=datamax, cbar_kwargs=cbar_kwargs)
    axs[i,j].annotate(f'N models:{all_N}',xy=(1, 0), xycoords='axes fraction',
            xytext=(-20, -5), textcoords='offset pixels',
            horizontalalignment='right',
            verticalalignment='bottom')
    axs2[i,j].annotate(f'N models:{all_N}',xy=(1, 0), xycoords='axes fraction',
            xytext=(-20, -5), textcoords='offset pixels',
            horizontalalignment='right',
            verticalalignment='bottom')
    axs2[i,j].coastlines()
    N_plot+=1
fig.tight_layout()
fig2.tight_layout()
fig.savefig(f'fig2_{experiments[0]}_fullscales.png',dpi=150)
fig2.savefig(f'fig2_{experiments[0]}_limitscales.png',dpi=150)
plt.show()

                
                
                