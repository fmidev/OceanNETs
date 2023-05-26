#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 16 14:16:06 2023

@author: Antti-Ilari Partanen (antti-ilari.partanen@fmi.fi)
"""

from itertools import product
import pandas as pd
import matplotlib.pyplot as plt


ifile='../data/Case study 1 - deployment rates.xlsx'

start_years=['2030', '2035', '2040']
volumes=['low', 'high']
scenarios=list()

for start_year, volume in product(start_years,volumes):
    scenarios.append(start_year+'-'+volume)
    
    
data=dict()    

for scenario in scenarios:
    data[scenario]=pd.read_excel(ifile,sheet_name=scenario,index_col=0)
    
    
fig,ax=plt.subplots(2,3,figsize=(15,10))

for i, start_year in enumerate(start_years):
    for j, volume in enumerate(volumes):
        scenario=start_year+'-'+volume
    
        for region in ['China','Europe','US']:
            data[scenario]['dep-'+region].plot(ax=ax[j,i],label=region)
        
        # Plot global total
        data[scenario].iloc[:,0:3].sum(axis=1).plot(ax=ax[j,i],label='Global total')
        ax[j,i].set_title(scenario)
        ax[j,i].set_ylim([0,5.5])
        ax[j,i].set_ylabel('Gt CaO yr$^-1$')
        
ax[j,i].legend(bbox_to_anchor=(0, -0.2),ncol=4)
        
fig.savefig('deployment_rates.png',dpi=150)
        
        