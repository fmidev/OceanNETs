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

'''
Figure 3: Time-series of OAE carbon fluxes globally, and for the three regions (and the fractions of the global carbon fluxes for the three regions in a separate panel)
'''

ctrl='esm-ssp534-over'
datapath='/Volumes/ONETs-SSD/Oceannets-all-models/data/'
#datapath='/Users/bergmant/Documents/projects/H2020-OceanNET/Oceannets-all-models/'
#vars=['talkos_reg','dissicos_reg','fgco2_glob','tasga','xco2atmga','spco2_reg','pH_reg']#,'omegac_reg']#,'fgoae_glob']
#vars=['talkos_reg']#,'dissicos_reg','fgco2_glob','tasga']#,'xco2atmga','fgoae_glob']
#vars=['dissicos_reg','dissicos_reg','fgco2_glob']#,'pH_reg','fgoae_glob']
vars=['fgco2_glob','pH_reg']#,'dissicos_reg','fgco2_glob','fgoae_glob']
#vars=['dissic_glob']#,'pH_reg','dissicos_reg','fgco2_glob','fgoae_glob']
#vars=['xco2atmga']#'spco2_reg']#,'fgoae_glob']
mod='NorESM2-LM'
#mod='AWI-ESM',
colors={'EC-Earth':'b','NorESM2-LM':'r','AWI-ESM':'g','FOCI':'y'}
lintypes={'esm-ssp534-over-2040high-dcr':'--','esm-ssp534-over-2040high':'-','esm-ssp534-over-2040high-term':':'}
#mod='EC-Earth'
models=['NorESM2-LM','AWI-ESM','EC-Earth','FOCI']
#models=['NorESM2-LM','EC-Earth']
#models=['NorESM2-LM']#,'AWI-ESM']#,'EC-Earth']
#models=['EC-Earth']
#models=['FOCI']
#models=['AWI-ESM']#,'EC-Earth']
#var='talkos_glob'
#var='dissicos_reg_global'
#var='global'


varname={}
# varname['talkos_reg','AWI-ESM']=['talkos_reg_global']#,'US_EEZ','CHINA_EEZ','EU_EEZ']
# varname['talkos_reg','EC-Earth']=['global','US_EEZ']#,'CHINA_EEZ','EU_EEZ']
# varname['talkos_reg','NorESM2-LM']=['talkos_glob']#,'US_EEZ','CHINA_EEZ','EU_EEZ']

#varname['talkos_glob','NorESM2-LM']='talkos_glob'
varname['pH_reg','AWI-ESM']='pH_reg_eu_eez'
varname['pH_reg','EC-Earth']='EU_EEZ'
varname['pH_reg','NorESM2-LM']='pH_EU'
varname['pH_reg','FOCI']='pH_3'
# varname['dissicos_reg','AWI-ESM']='dissicos_reg_global'
# varname['dissicos_reg','EC-Earth']='global'
# varname['dissicos_reg','NorESM2-LM']='dissicos_glob'
varname['fgco2_glob','AWI-ESM']='fgco2_glob'
varname['fgco2_glob','EC-Earth']='global'
varname['fgco2_glob','NorESM2-LM']='fgco2_glob'
varname['fgco2_glob','FOCI']='fgco2_glob'
# varname['dissic_glob','AWI-ESM']='dissic_glob'
# varname['dissic_glob','EC-Earth']='global'
# varname['dissic_glob','NorESM2-LM']='dissic_glob'
# varname['xco2atmga','NorESM2-LM']='xco2atmga'
# varname['xco2atmga','EC-Earth']='global'
# varname['xco2atmga','AWI-ESM']='pco2atmga'
# varname['spco2_reg','NorESM2-LM']='spco2_glob'
# varname['spco2_reg','EC-Earth']='global'
# varname['spco2_reg','AWI-ESM']='spco2_reg_global'
# varname['tasga','NorESM2-LM']='tasga'
# varname['tasga','EC-Earth']='global'
# varname['tasga','AWI-ESM']='tasga'
# varname['omegac_reg','NorESM2-LM']='omegac_glob'
# varname['omegac_reg','EC-Earth']='global'
# varname['omegac_reg','AWI-ESM']='omegac_reg_global'

labels={'talkos_reg':'Total Alkalinity at \n Ocean Surface','pH_reg':'pH at Ocean Surface','dissicos_reg':'Dissolved inorganic carbon at\n Ocean surface','fgco2_glob':'CO2 flux to ocean','xco2atmga':'CO2 concentration \n at surface','tasga':'Air temperature at \n the surface','spco2_reg':'CO2 concentration in \n ocean surface','omegac_reg':'OmegaC'}
#f3,ax3=plt.subplots()
# for j in vars:
#     data=None
#     for d in (glob.glob(datapath+'/'+ctrl+'/'+j+'_'+mod+'*')):
#         print(d)
#         data_a=xr.open_dataset(d)
#         print(data_a)
#data_a[varname[vars[0],mod]].plot()
#datapath='/Volumes/ONETs-SSD/Oceannets-all-models/data/'
#mod='EC-Earth'
sumdata={}
ctrl='esm-ssp534-over'
exp='esm-ssp534-over-2040high-term'
exp='esm-ssp534-over-2040high-dcr'
exp='esm-ssp534-over-2040high'
#experiments=['esm-ssp534-over-2040high-term','esm-ssp534-over-2040high-dcr','esm-ssp534-over-2040high']
experiments=['esm-ssp534-over-2040high-dcr','esm-ssp534-over-2040high']
#experiments=['esm-ssp534-over-2040high-dcr' ]
#experiments=['esm-ssp534-over-2040high']
exp_label={'esm-ssp534-over-2040high':'ALK 2040-2100','esm-ssp534-over-2040high-term':'TERM','esm-ssp534-over-2040high-dcr':'DCR 2040-2100','esm-ssp534-over':'SSP534-OS'}
ctrl_data={}
exp_data={}
for j in vars:
    f1,ax1=plt.subplots(2,2,figsize=(9,6))
    f1b,ax1b=plt.subplots(figsize=(9,6))
    f1c,ax1c=plt.subplots(nrows=2,figsize=(9,6))
    f2,ax2=plt.subplots(figsize=(9,6))
    f3,ax3=plt.subplots(figsize=(9,6))
    for mod in models:
        data_a=None
        data_b=None
        #if mod=='EC-Earth':
        #    if j == 'xco2atmga':
        #        continue
        if mod=='FOCI':
            filename='timeseries_*.nc'
        else:
            filename=j+'_'+mod+'*'
        for d in (glob.glob(datapath+'/'+ctrl+'/'+filename)):
            print(d)
            
            data_a=xr.open_dataset(d)
            if mod=='NorESM2-LM' :
                datetimeindex = data_a.indexes['time'].to_datetimeindex()
                data_a['time']=datetimeindex                
            # elif mod=='AWI-ESM':# and data_b.time.year < 2029:
            #     datetimeindex = data_a.indexes['time'].to_datetimeindex()
            #     data_a['time']=datetimeindex                
            #    if exp=='ssp534-over-2040high-dcr':
            #        continue
            elif mod=='EC-Earth':
                print(d)
                print(data_a)
                
            # elif mod=='EC-Earth':
            #     init_a=data_a.sel(time=slice('2015-01-01','2028-12-31'))
            #print(data_a)
            ctrl_data[mod]=data_a
        #print(data_a)
        print(ctrl_data.keys())
        #print(ctrl_data['EC-Earth'])
    #for j in vars:
        #data=None
        # if data_a==None:
        #     continue
        #print(j,mod,data_a)
        if varname[j,mod] in (data_a):
            p2,=(data_a)[varname[j,mod]].rolling(time=60).mean().plot(ax=ax2,label=mod+' '+exp_label[ctrl],lw=3,color=colors[mod],linestyle=lintypes[exp])
        for exp in experiments:
            for d in (glob.glob(datapath+'/'+exp+'/'+filename)):
                print(d)
                # if mod=='AWI-ESM':# and data_b.time.year < 2029:
                #     if exp=='esm-ssp534-over-2040high-dcr':
                #         continue
                print(exp)
                data_b=xr.open_dataset(d)#,decode_times=False)
                if mod=='NorESM2-LM':
                    datetimeindex = data_b.indexes['time'].to_datetimeindex()
                    data_b['time']=datetimeindex                
                elif mod=='EC-Earth':# and data_b.time.year < 2029:
                    if pd.DatetimeIndex(data_b.time).year[0]==2040:
                        print (pd.DatetimeIndex(data_b.time).year[0])
                        print (data_a)
                        print (pd.DatetimeIndex(data_a.time).year[0])
                        
                        #init_a=data_a.sel(time=slice('2015-01-01','2039-12-31'))
                        print(data_a.sel(time=slice('2015-01-01','2039-12-31')).time[-1])
                        print(data_b.time[0])
                        data_c=xr.concat([data_a.sel(time=slice('2015-01-01','2039-12-31')),data_b],dim='time')
                        data_b=data_c
                    elif pd.DatetimeIndex(data_b.time).year[0]==2030:
                    
                        print (pd.DatetimeIndex(data_b.time).year[0])
                        
                        init_a=data_a.sel(time=slice('2015-01-01','2028-12-31'))
                        print(data_a.sel(time=slice('2015-01-01','2028-12-31')).time[-1])
                        print(data_b.time[0])
                        data_c=xr.concat([data_a.sel(time=slice('2015-01-01','2029-12-31')),data_b],dim='time')
                        data_b=data_c.copy()
                #elif mod=='AWI-ESM':# and data_b.time.year < 2029:
                #    if exp=='ssp534-over-2040high-dcr':
                #        continue
                #     dfa
                #     data_c=xr.concat([data_b,data_a.sel(time=slice('2015-01-01','2028-12-31'))],dim='time')
                #     data_b=data_
                print('hep')
                #print(data_b)
            if data_b != None:
                exp_data[mod,exp]=data_b
            #data_b[varname[vars[0],mod]].plot()
            #(data_b-data_a)['CHINA_EEZ'].plot()
            #(data_b-data_a)['EU_EEZ'].plot()
            #(data_b-data_a)['US_EEZ'].plot()
            # if data_b==None:
            #     #print(data_b)
            #     print('b not read')
            #     continue    
            print(varname[j,mod])
            #print('a',data_a)
            #print('b',data_b)
            #print((data_b-data_a).sum().variables)
            print(mod,exp)
            #print(data_a)
            #print(data_b)
            #if exp=='esm-ssp534-over-2040high-dcr' and mod=='AWI-ESM':
            #    continue
            #elif mod=='AWI-ESM' and j=='xco2atmga':
            #    continue
            print(varname[j,mod])
            #print(data_b)
            #print(data_a)
            if varname[j,mod] in (data_b-data_a).sum():
                sumdata[mod,j]=(data_b-data_a).sum()[varname[j,mod]].values
                #print((data_b-data_a).sum()[varname[j,mod]].values)
                #plt.figure()
                #plt.plot((data_b-data_a)[varname[j,mod]].time)
                #print ((data_b-data_a)[varname[j,mod]].time)
                #print(data_b,varname[j,mod],j,mod)
                #print(data_a)
                print((data_b-data_a)[varname[j,mod]].rolling(time=60).mean())
                p1,=(data_b-data_a)[varname[j,mod]].rolling(time=60).mean().plot(ax=ax1[0,0],label=mod+' '+exp_label[exp],lw=3,color=colors[mod],linestyle=lintypes[exp])
                ax1[0,0].set_ylabel('Change in '+labels[j]+' ['+data_b[varname[j,mod]].units+']')
                p1,=(data_b-data_a)[varname[j,mod]].rolling(time=60).mean().plot(ax=ax1b,label=mod+' '+exp_label[exp],lw=3,color=colors[mod],linestyle=lintypes[exp])
                ax1b.set_ylabel('Change in '+labels[j]+' ['+data_b[varname[j,mod]].units+']')
                if j=='pH_reg':
                    if exp=='esm-ssp534-over-2040high':
                        p1,=(data_b-data_a)[varname[j,mod]].rolling(time=60).mean().plot(ax=ax1c[0],label=mod+' '+exp_label[exp],lw=3,color=colors[mod],linestyle=lintypes[exp])
                        ax1c[0].set_ylabel('Change in '+labels[j]+' ['+data_b[varname[j,mod]].units+']')
                    elif exp=='esm-ssp534-over-2040high-dcr':
                        p1,=(data_b-data_a)[varname[j,mod]].rolling(time=60).mean().plot(ax=ax1c[1],label=mod+' '+exp_label[exp],lw=3,color=colors[mod],linestyle=lintypes[exp])
                        ax1c[1].set_ylabel('Change in '+labels[j]+' ['+data_b[varname[j,mod]].units+']')
                else:
                
                    if exp=='esm-ssp534-over-2040high':
                        p1,=(data_b-data_a)[varname[j,mod]].cumsum().plot(ax=ax1c[0],label=mod+' '+exp_label[exp],lw=3,color=colors[mod],linestyle=lintypes[exp])
                        ax1c[0].set_ylabel('Change in '+labels[j]+' ['+data_b[varname[j,mod]].units+']')
                    elif exp=='esm-ssp534-over-2040high-dcr':
                        p1,=(data_b-data_a)[varname[j,mod]].cumsum().plot(ax=ax1c[1],label=mod+' '+exp_label[exp],lw=3,color=colors[mod],linestyle=lintypes[exp])
                        ax1c[1].set_ylabel('Change in '+labels[j]+' ['+data_b[varname[j,mod]].units+']')
                    
            #p1,=(data_b-data_a)[varname[j,mod]].plot(ax=ax1,label=mod)
            #p1,=(data_b-data_a)[varname[j,mod]].plot(ax=ax1,label=mod)
            # if mod=='NorESM2-LM':
            #     plt.figure()
            #     plt.plot(data_a.time)
            #     plt.show()
            #     fsag
            #plt.figure()
            #print(j,mod,varname[j,mod])
            if varname[j,mod] in (data_a):
                p2,=(data_a)[varname[j,mod]].plot(ax=ax2,label=mod+' '+exp_label[ctrl],lw=3,color=colors[mod],linestyle=lintypes[exp])
            if varname[j,mod] in (data_b):
                p3,=(data_b)[varname[j,mod]].rolling(time=60).mean().plot(ax=ax2,label=mod+' '+exp_label[exp],lw=3,color=colors[mod],linestyle=lintypes[exp])
                ax2.set_ylabel(labels[j]+' ['+data_b[varname[j,mod]].units+']',fontsize=20)
                ax2.xaxis.label.set_size(16)
            if varname[j,mod] in (data_b-data_a).sum():
                p4,=(data_b-data_a)[varname[j,mod]].rolling(time=60).mean().plot(ax=ax3,label=mod+' '+exp_label[exp],lw=3,color=colors[mod],linestyle=lintypes[exp])
                ax3.set_ylabel('Change in '+labels[j]+' ['+data_b[varname[j,mod]].units+']',fontsize=20)
                ax3.xaxis.label.set_size(16)
                ax3.yaxis.label.set_size(16)
                ax3.tick_params(labelsize=14)
            #ax3.plot(data_a.time.values)
            #print(data_a.time)
        print (j, exp_data.keys()) 
    ax1[0,0].legend(loc='upper left',fontsize=14)
    ax1b.legend(loc='upper left',fontsize=14)
    ax1c[0].legend(loc='upper left',fontsize=14)
    ax1c[1].legend(loc='upper left',fontsize=14)
    ax2.legend(loc='upper left',fontsize=14)
    ax3.legend(loc='upper left',fontsize=14)
    f1.savefig('fig3_'+labels[j]+'__rolling60mon.png')    
    f1b.savefig('fig3_'+labels[j]+'_rolling60mon.png')    
    f1c.savefig('fig3_'+labels[j]+'_sumulative.png')    
    f2.savefig('draft/'+labels[j]+'_v2.png')    
    f3.savefig('draft/delta_'+labels[j]+'_rolling_v2.png')    
print(sumdata)
plt.show()
