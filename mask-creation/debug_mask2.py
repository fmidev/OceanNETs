#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 11:37:05 2023

@author: Antti-Ilari Partanen (antti-ilari.partanen@fmi.fi)
"""

import xarray as xr
import numpy as np
from mask_aux import read_deployment_data
import matplotlib.pyplot as plt
import pandas as pd

data=xr.open_dataset('../data/lime_mask_cao_2040-high.nc')
data_tb=xr.open_dataset('../data/lime_mask_v2.cao_tommi.nc')


scenario='2040-high'
dep_data=read_deployment_data()

gt_2_g=1e6*1e9 # Gt to g
n_cao=56.1 

area_eu=4426626029151.475
area_china=970136511529.54
area_us=5669494357484.002

seconds_in_year=60**2*24*365

diff_rel=(data-data_tb)/data_tb*100

#Find indeces where difference between datasets is largest (most negative)
indices=np.where((diff_rel.lime_mask.values)==np.nanmin(diff_rel.lime_mask.values))

#Conversion from mol / (s m2) to Gt / year
f_eu=seconds_in_year*area_eu*n_cao/gt_2_g
f_china=seconds_in_year*area_china*n_cao/gt_2_g
f_us=seconds_in_year*area_us*n_cao/gt_2_g

value_new=data.lime_mask[indices[0][0],indices[1][0],indices[2][0]].values*f_eu
value_tb=data_tb.lime_mask[indices[0][0],indices[1][0],indices[2][0]].values*f_eu




print('Annual value in new file',value_new)
print('Annual value in Tommi''s file',value_tb)
print('Value from Excel', dep_data[scenario].loc[2030,'dep-Europe'])


#Pick points in all three regions and plot time series
fig,ax=plt.subplots(3,1)
(data.lime_mask.sel(lon=23.75,lat=33.75)*f_eu).plot(ax=ax[0], label='New')
(data_tb.lime_mask.sel(lon=23.75,lat=33.75)*f_eu).plot(ax=ax[0], label='TB')
ax[0].plot(pd.to_datetime(dep_data[scenario].index, format='%Y'), dep_data[scenario]['dep-Europe'], linestyle='dotted',label='Original')
ax[0].legend()

(data.lime_mask.sel(lon=-123.8,lat=37.25, method='nearest')*f_us).plot(ax=ax[1], label='New')
(data_tb.lime_mask.sel(lon=-123.8,lat=37.25, method='nearest')*f_us).plot(ax=ax[1], label='TB')
ax[1].plot(pd.to_datetime(dep_data[scenario].index, format='%Y'), dep_data[scenario]['dep-US'], linestyle='dotted',label='Original')
ax[1].legend()

(data.lime_mask.sel(lon=123.8,lat=28.6, method='nearest')*f_china).plot(ax=ax[2], label='New')
(data_tb.lime_mask.sel(lon=123.8,lat=28.6, method='nearest')*f_china).plot(ax=ax[2], label='TB')
ax[2].plot(pd.to_datetime(dep_data[scenario].index, format='%Y'), dep_data[scenario]['dep-China'], linestyle='dotted',label='Original')
ax[2].legend()

