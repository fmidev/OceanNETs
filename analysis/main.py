#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 13:50:25 2022

@author: bergmant
"""

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature


def plot_var(data,ivar):
    plt.figure()
    fig,ax=plt.subplots(ncols=3,subplot_kw={'projection': ccrs.Robinson()})
    #ax=plt.axes(projection=ccrs.Robinson())
    for i,iexp in enumerate(data):
        print(i,ivar,data[iexp][ivar][ivar])
        cbar_kwargs = {'orientation':'horizontal', 'shrink':0.6, 'aspect':40,'label':data[iexp][ivar][ivar].attrs['standard_name']}
        if 'lat' in data[exps[0]][ivar][ivar].coords:
            xlat='lat'
            xlon='lon'
        else:
            xlat='latitude'
            xlon='longitude'
        if 'level' in data[exps[i]][ivar].coords:

            (data[iexp][ivar][ivar]).sel(level=1).mean(dim='time').plot(ax=ax[i],x=xlon,y=xlat, transform=ccrs.PlateCarree(),add_labels=True,cbar_kwargs=cbar_kwargs)
            ax[i].set_title(data[iexp][ivar][ivar].long_name)

        else:
            print(data[iexp][ivar][ivar])
            (data[iexp][ivar][ivar]).mean(dim='time').plot(ax=ax[i],x=xlon,y=xlat, transform=ccrs.PlateCarree(),add_labels=True,cbar_kwargs=cbar_kwargs)
            ax[i].set_title(data[iexp][ivar][ivar].long_name)

        gl = ax[i].gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=0.5, color='gray', alpha=0.5, linestyle='--')

def plot_diffvar(data,ivar):
    plt.figure()
    fig,ax=plt.subplots(ncols=1,subplot_kw={'projection': ccrs.Robinson()})
    # find experiment names from data dict
    print (ivar)
    exps=[]
    for i in data:
        exps.append(i)
    # colorbar settings
    cbar_kwargs = {'orientation':'horizontal', 'shrink':0.6, 'aspect':40, 'label':ivar}
    #print (data[exps[0]])
    if 'lat' in data[exps[0]][ivar][ivar].coords:
        xlat='lat'
        xlon='lon'
    else:
        xlat='latitude'
        xlon='longitude'
    if 'level' in data[exps[0]][ivar].coords:
        (data[exps[1]][ivar][ivar] - data[exps[0]][ivar][ivar]).sel(level=1).mean(dim='time').plot(ax=ax,x=xlon,y=xlat, transform=ccrs.PlateCarree(),add_labels=True,cbar_kwargs=cbar_kwargs)
        ax.set_title(data[exps[0]][ivar][ivar].long_name)
    else:

        # plot
        (data[exps[1]][ivar][ivar] - data[exps[0]][ivar][ivar]).mean(dim='time').plot(ax=ax,x=xlon,y=xlat, transform=ccrs.PlateCarree(),add_labels=True,cbar_kwargs=cbar_kwargs)
        ax.set_title(data[exps[0]][ivar][ivar].long_name)

    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    #ax.set_title(data[exp][var][var].long_name)
def plot_timeseries(data,exp,var,weight,ax=None,operator='sum',type='ocean'):
    if type=='ocean':
        dims=('j','i')
    elif type=='land':
        dims=('lat','lon')
    elif type=='atm':
        dims=('lat','lon')
    else:
        print('unknown type %s'%type)
    if operator=='sum':
        annual_mean=data[exp][var][var].weighted(weight).sum(dim=dims)
    elif operator=='mean':
        annual_mean=data[exp][var][var].weighted(weight).mean(dim=dims)
    else:
        print('unknown operator for timeseries plot')
        return False

    if ax==None:
        fig,ax=plt.subplots(111)
    else:
        print('Using provided axes for timeseries plotting')
    annual_mean.plot(ax=ax,label=exp+': '+var)
    ax.set_title(data[exp][var][var].long_name)
def read_area(inputfile='/Users/bergmant/Documents/projects/H2020-OceanNET/data/ON01/areacello_Ofx_EC-Earth3-CC_esm-pi-cdr-pulse_r1i1p1f1_gr.nc'):
    area=xr.open_dataset(inputfile)
    area2=area.areacello.fillna(0)#/area.areacello.sum()
    area2.name="weights"
    return area2

def read_area_land(inputfile='/Users/bergmant/Documents/projects/H2020-OceanNET/data/areacella_fx_EC-Earth3_ssp245_r13i1p1f2_gr.nc'):
    area=xr.open_dataset(inputfile)
    area2=area.areacella.fillna(0)#/area.areacello.sum()
    area2.name="weights"
    return area2

#exp=['C534','ON01']
#ocean_varlist=['intdic','dissicnatos','talkos','talknatos','fgco2','co3satcalcos','co3os','calcos','zoocos','sos','ppos','sios',
#               'po4os','dissocos','dissicos','phycos','detocos','chldiatos','chlmiscos','zmicroos','zmesoos']
ocean_varlist=['talkos','co3satcalcos','fgco2','co3os','intdic','calcos','phos','dpco2','chldiatos']
               #,'co3satcalcos','co3os','calcos','zoocos','sos','ppos','sios',
               #'po4os','dissocos','dissicos','phycos','detocos','chldiatos','chlmiscos','zmicroos','zmesoos']
#extra_varlist=[]
#land_varlist=['npp','gpp','nbp','evspsblveg','evspsblsoi','lai','cLeaf','cVeg','cLitter','cProduct','cRoot','ra','rh']
land_varlist=['npp','gpp','nbp']#,'evspsblveg','evspsblsoi','lai','cLeaf','cVeg','cLitter','cProduct','cRoot','ra','rh']
atm_1d_varlist=['co2mass','tas']
extra_varlist=['netAtmosLandCO2Flux']#,'nppLut']
varlist={}
varlist['Omon']=ocean_varlist
varlist['Emon']=extra_varlist
varlist['Lmon']=land_varlist
varlist['Amon']=atm_1d_varlist
datapath='/Users/bergmant/Documents/projects/H2020-OceanNET/data/'
#datapath='/Users/bergmant/Documents/projects/H2020-OceanNET/data/oceannets-output'
exps=['C534','C5SI']#r,'ON01']
#exp=exp[0]
year='*2016*'
index=0
expdict={}
for iexp in exps:
    data={}
    for table in varlist:
        for ivar in varlist[table]:
            if table=='Lmon':
                grid='gr'
            else:
                grid='gr'
            #fname=datapath+'/'+iexp+'/'+ivar+'_'+table+'_EC-Earth3-CC_esm-pi-cdr-pulse_r1i1p1f1_'+grid+'_'+year+'01-'+year+'12.nc'
            fname=datapath+'/'+iexp+'/'+ivar+'_'+table+'_EC-Earth3-CC_esm-pi-cdr-pulse_r1i1p1f1_'+grid+'_'+year+'.nc'
            print(fname)
            data[ivar]=xr.open_mfdataset(fname)
            #print(data[ivar][ivar])
    expdict[iexp]=data

for table in varlist:
    for ivar in varlist[table]:
        if ivar=='co2mass':
            continue
        print(ivar)
        plot_diffvar(expdict,ivar)
        plot_var(expdict,ivar)
area2=read_area()
arealand=read_area_land()
#hepuli1=expdict[exps[0]]['fgco2']['fgco2'].mean(dim='time').weighted(area2)
#hepuli2=expdict[exps[1]]['fgco2']['fgco2'].mean(dim='time').weighted(area2)
#print (hepuli1)
#hepulisum1=hepuli1.sum(dim=('j','i'))*3600*24*365
#hepulisum2=hepuli2.sum(dim=('j','i'))*3600*24*365
#print(hepulisum1)
#print(hepulisum2)

fig,ax=plt.subplots(1)
plot_timeseries(expdict, exps[0], 'fgco2', area2,ax,'sum')
plot_timeseries(expdict, exps[1], 'fgco2', area2,ax,'sum')
#plot_timeseries(expdict, exps[2], 'fgco2', area2,ax,'sum')
ax.legend()
fig,ax=plt.subplots(1)
plot_timeseries(expdict, exps[0], 'talkos', area2,ax,'mean')
plot_timeseries(expdict, exps[1], 'talkos', area2,ax,'mean')
#plot_timeseries(expdict, exps[2], 'talkos', area2,ax,'mean')
#plot_timeseries(expdict, exps[1], 'talkos', area2,ax,'mean')
ax.legend()
# fig,ax=plt.subplots(1)
# plot_timeseries(expdict, exps[0], 'intdic', area2,ax,'mean')
# plot_timeseries(expdict, exps[1], 'intdic', area2,ax,'mean')
# ax.legend()
# fig,ax=plt.subplots(1)
# plot_timeseries(expdict, exps[0], 'dissicnatos', area2,ax,'mean')
# plot_timeseries(expdict, exps[1], 'dissicnatos', area2,ax,'mean')
# ax.legend()
for var in ocean_varlist:
    fig,ax=plt.subplots(1)
    plot_timeseries(expdict, exps[0], var, area2,ax,'mean')
    plot_timeseries(expdict, exps[1], var, area2,ax,'mean')
    #plot_timeseries(expdict, exps[2], var, area2,ax,'mean')
    ax.legend()
for var in land_varlist:
    print(var)
    fig,ax=plt.subplots(1)
    plot_timeseries(expdict, exps[0], var, arealand,ax,'mean','land')
    plot_timeseries(expdict, exps[1], var, arealand,ax,'mean','land')
    #plot_timeseries(expdict, exps[2], var, arealand,ax,'mean','land')
    ax.legend()
for var in extra_varlist:
    print(var)
    fig,ax=plt.subplots(1)
    plot_timeseries(expdict, exps[0], var, arealand,ax,'mean','land')
    plot_timeseries(expdict, exps[1], var, arealand,ax,'mean','land')
    #plot_timeseries(expdict, exps[2], var, arealand,ax,'mean','land')
    ax.legend()
for var in atm_1d_varlist:
    print(var)
    fig,ax=plt.subplots(1)
    plot_timeseries(expdict, exps[0], var, arealand,ax,'mean','land')
    plot_timeseries(expdict, exps[1], var, arealand,ax,'mean','land')
    #plot_timeseries(expdict, exps[2], var, arealand,ax,'mean','land')
    ax.legend()
#for var in atm_1d_varlist:
fig,ax=plt.subplots(1)
t1=expdict[exps[0]]['co2mass']#/5.1480e18*(28.9/44)
t2=expdict[exps[1]]['co2mass']#/5.1480e18*(28.9/44)
t1.co2mass.data=t1.co2mass.data/5.1480e18*(28.9/44)*1e6
t2.co2mass.data=t2.co2mass.data/5.1480e18*(28.9/44)*1e6
t1.co2mass.plot(ax=ax,label=exps[0]+':co2mass')
t2.co2mass.plot(ax=ax,label=exps[1]+':co2mass')
#expdict[exps[0]]['co2mass'].co2mass/5.1480e18*(28.9/44).plot(ax=ax,label=exps[0]+':co2mass')
#expdict[exps[1]]['co2mass'].co2mass/5.1480e18*(28.9/44).plot(ax=ax,label=exps[1]+':co2mass')
ax.legend()
# available fixed variables

# layer thickness
# thkcello
# ocean mass
# masscello
# grid area
#areacello