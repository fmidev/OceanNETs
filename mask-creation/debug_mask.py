#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 11:37:05 2023

@author: Antti-Ilari Partanen (antti-ilari.partanen@fmi.fi)
"""

import xarray as xr
import numpy as np

data=xr.open_dataset('../data/lime_mask_v2.cao.nc')
data_tb=xr.open_dataset('../data/lime_mask_v2.cao_tommi.nc')


gt_2_g=1e6*1e9 # Gt to g
n_cao=56.1 

area_eu=25395163896767.07

seconds_in_year=60**2*24*365

diff_rel=(data-data_tb)/data_tb*100

#Find indeces where difference between datasets is largest (most negative)
indices=np.where((diff_rel.lime_mask.values)==np.nanmin(diff_rel.lime_mask.values))

#Conversion from mol / (s m2) to Gt / year
f=seconds_in_year*area_eu*n_cao/gt_2_g*area_eu

value_new=data.lime_mask[indices[0][0],indices[1][0],indices[2][0]].values*f
value_tb=data_tb.lime_mask[indices[0][0],indices[1][0],indices[2][0]].f




print('Annual value in new file',value_new*f)
print('Annual value in Tommi''s file',value_tb*f)

