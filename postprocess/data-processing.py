# Processing for data portal

import numpy as np
import xarray as xr
#import paramiko
import os
#import netCDF4 as nc
import glob
import calendar
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import pandas as pd 


vars={
    'talk':'talk_glob',
    #global means
    'fgco2': 'fgco2_glob',
    'alkalinization': 'fgoae_glob',
    'fco2removal': 'fgdcr_glob',
    'tos': 'tosga',
    #'dissic': 'dissic_glob',
    'co2': 'xco2atmga',
    'tas': 'tasga',
    'nbp': 'nbp_glob',
    'gpp': 'gpp_glob',
    'cVeg': 'cveg_glob',
    'cSoil': 'csoil_glob',
    #region means:
    'spco2': 'spco2_reg',
    'dissicos': 'dissicos_reg',
    'talkos': 'talkos_reg',
    'omegac': 'omegac_reg',
    'phos': 'pH_reg',
#    'omegaa': 'omegaa_reg', #not available
    #3d  integrated
    'intdic':'dissic_glob',

    }
units = {
        'fgco2_glob': 'PgC yr-1',#monthly????
        'fgoae_glob': 'Pmol TA yr-1', #monthly
        'fgdcr_glob':'Pg C yr-1',
        'tasga':'degC',
        'tosga':'degC',
        'dissic_glob': 'mol m-3',
        'xco2atmga': 'ppm',
        'nbp_glob': 'Pg C yr-1',
        'gpp_glob':'Pg C yr-1',
        'csoil_glob':'PgC ',
        'cveg_glob':'PgC',
        'spco2_reg': 'uAtm',
        'dissicos_reg': 'mol m-3',
        'talkos_reg':'mol m-3',
        'pH_reg':'-',
        'omegac_reg':'-',
        'dissic_glob':'',
        'talk_glob':'P mol',

}
expnames={
'C5AF':'esm-ssp534-over',
'D30H':'esm-ssp534-over-2030high',
'D40H':'esm-ssp534-over-2040high',
'E40H':'esm-ssp534-over-2040high-term',
'C40H':'esm-ssp534-over-2040high-dcr'}

n_c=12.011
unit_conversions = {'cSoil':[1e-12,True],  #kgm-2 -> PgC m-2
                'cVeg':[1e-12,True],   #kgm-2 -> PgC m-2
                'gpp':[1e-12*3600*24,True],    #kgm-2s-1 -> PgCm-2day-1
                'nbp':[1e-12*3600*24,True],    #kgm-2s-1 -> PgCm-2day-1
                'fgco2':[1e-12*3600*24,True],    #kgm-2s-1 -> PgCm-2day-1
                'dissicos':[1,False],    #molm-3
                'intdic':[1e-12,False],    #molm-3
                'talkos':[1,False],    #molm-3
                'talk':[1e-15,False],    #molm-3
                'spco2':[1e6/101325,False],    # Pa -> uAtm
                'co2':[1e6,False],    # mol/mol -> umol/mol or ppm  
                'ph':[1,False], # no conversion unit of 1
                'omegac':[1,False], #no conversion unit of 1
                'fco2removal':[1e-12*3600*24,True], # XXXX
                'alkalinization':[1e-15*3600*24,True], # xxxxx->Pmol m-3 day-1
                'phos':[1,False], # xxxxx->Pmol m-3 day-1
                'tas':[-273.15,False], # K -> C 
                'tos':[-273.15,False] # K -> C 
                }  
def convert_units(data,var):
    monlen=[31,28,31,30,31,31,30,31,30,31,30,31]
    monlen_leap=[31,29,31,30,31,31,30,31,30,31,30,31]

    for i in data:
        print(data[i],var)
        print(data[i].values[0])
        tmp=data[i].values[0]
        print (unit_conversions[var][1])
        if var=='alkalinization' or var == 'omegac' or var=='fco2removal2D':
            print(data[i].time_counter[0].dt.year)
            iyear = data[i].time_counter[0].dt.year
        else:
            print(var)
            iyear = data[i].time[0].dt.year
        
        print(iyear)
        print(vars[var])
        if vars[var][1]:
            if unit_conversions[var][1]:
                # calculate either new value based on monlenght
                # or year length (s-1 -> mon-1 or year-1)
                if calendar.isleap(iyear):
                    if 'yr-1' in units[vars[var]]:
                        days=366
                    else:
                        days=monlen_leap
                    data[i] = data[i]*days#monlen_leap
                else:   
                    if 'yr-1' in units[vars[var]]:
                        days=365
                    else:
                        days=monlen
                    data[i] = data[i]*days#monlen
        if var == 'tas' or var =='tos':  
            data[i].values = data[i].values + unit_conversions[var][0]
        else:
            print (data[i],unit_conversions[var][0])
            data[i].values = data[i].values * unit_conversions[var][0]
        print('heå',data[i].values[0]/tmp,unit_conversions[var][0]*31.0)
        print(tmp,data[i].values)
        #print(data[i].values[0]*unit_conversions[var][0])
    
    
def read_and_save_netcdf(local_path,exp):
    # Read data from the downloaded NetCDF file
    #matching_files = [f for f in local_files if glob.fnmatch.fnmatch(f, fname)]
    #fname='D40H_*pisces*_2D*'
    #var='fgco2'
    # extract area from areacello
    area=xr.open_dataset('areacello_Ofx_EC-Earth3-CC_esm-pi-cdr-pulse_r1i1p1f1_gr.nc',engine='h5netcdf')
    areab=area.areacello.fillna(0)
    area_ifs=xr.open_dataset('areacella_ifs.nc',engine='h5netcdf')
    area_atm=area_ifs.cell_area#.fillna)0)
    area_tm5=xr.open_dataset('areacella_tm5.nc',engine='h5netcdf')
    area_catm=area_tm5.cell_area#.fillna)0)
    area_pisces=xr.open_dataset('pisces_ga.nc',engine='h5netcdf')
    area_pis=area_pisces.areacello.fillna(0)

    # Download each matching file
    #testout='talkos_'
    #var=fname.split('_')[0]
    #print(var,fname,fname.strip('_'))
    #'fgco2_'
    country_mask_pis=xr.open_dataset('countrymask.b.nc')
    country_mask=xr.open_dataset('countrymask_v2.b.nc')
    #print(country_mask)
    
    # mask_EU=country_mask.where(country_mask.mask==30)
    # mask_US=country_mask.where(country_mask.mask==20)
    # mask_CHINA=country_mask.where(country_mask.mask==10)
    # print(mask_US) 
    
#    print(local_path)
#    print(fname)
    cyear,cmonth,cday,chour=date_string()
   
    version_id=f'v{cyear}{cmonth}{cday}'
    
    experiment_id=expnames[exp]
    ff=xr.open_dataset ('/Volumes/ONETs-SSD//D30H/tas_Amon_EC-Earth3-CC_esm-pi-cdr-pulse_r1i1p1f1_gr_205001-205012.nc',engine='h5netcdf')
    print(ff)
    
    regions=['glob','EU','US','CHINA']
    #vars={'co2':'xco2atmga'}
    if exp=='C5AF':
        startyear=2015
    else:
        byear=exp[1:-1]
        startyear=int(f'20{byear}')
    print (startyear)
    #print(f'20{startyear}')
    endyear=2101
    
    for var in vars:
        print(f'running post processing: {startyear} - {endyear}')
        print(var)
        #for syear in range(startyear,endyear):
        for syear in range(2016,2101):
            print(syear)
            model_id='EC-Earth'
            outpath=f'/Volumes/ONETs-SSD/{experiment_id}/'
            #outpath=f'testout/'
            isExist = os.path.exists(outpath)
            time_range='2015-2100'
            #syear='2050'
            if not isExist:
                # Create a new directory because it does not exist
                os.makedirs(outpath)
            testout=var+'_'
            outputname=f'{vars[var]}_{model_id}_{experiment_id}_{time_range}_{version_id}.{syear}.nc'
            if var=='co2':
                fname=f'tm5/{var}_*mon*{syear}*nc'   
            elif var=='alkalinization' or var== 'fco2removal2D' or var=='omegac':
                fname=f'{var}_*1m*{syear}*nc'   
            else:
                fname=f'{var}_*mon*{syear}*nc'   
            local_files = glob.glob(os.path.join(local_path, fname))
            print(local_files)
            print(os.path.join(local_path, fname))
            for local_file in local_files:
                print(local_file)
                year=local_file[-9:][:4]
                outti=testout+year+'.nc'
                # with nc.Dataset(os.path.join(local_path, local_file), 'r') as nc_file:
                #     # Do something with the NetCDF data
                #     # For example, print information about the variables and dimensionsema
                #     print(nc_file.variables)
                #     print(nc_file.dimensions)
                syear=local_file.split('_')[6][:4]
                xrds=xr.open_dataset(local_file,engine='h5netcdf')
                #plt.figure()
                data={}
                print(vars[var],vars[var][-2:])
                if vars[var][-3:]=='reg':# and vars[var][-2:]!='ga':
                    for reg in regions:
                        if reg =='CHINA':
                            match=10
                        elif reg =='US':
                            match=20
                        elif reg =='EU':
                            match=30
                        else:
                            match=-999
                        if  var=='omegac':
                            if reg== 'glob':
                                print(reg)
                                data[reg]=(xrds[var]).weighted(area_pis).mean(dim=('y','x'), skipna=True)
                            else:
                                print(reg)
                                data[reg]=xrds[var].where(country_mask_pis.mask==match).weighted(area_pis).mean(dim=('y','x'), skipna=True)                    
                        else:
                            if reg== 'glob':
                                data[reg]=xrds[var].weighted(areab).mean(dim=('i','j'), skipna=True)

                            else:
                                data[reg]=xrds[var].where(country_mask.mask==match).weighted(areab).mean(dim=('i','j'), skipna=True)
                else:
                    reg='glob'
                    if var=='tas':
                        data[reg]=(xrds[var]).weighted(area_atm).mean(dim=('lat','lon'), skipna=True)
                        print(f'tas {data[reg]}')
                        
                    elif var=='tos':
                        data[reg]=(xrds[var]).weighted(areab).mean(dim=('i','j'), skipna=True)
                        print(f'tas {data[reg]}')
                        
                    elif var in ['nbp','gpp','cSoil','cVeg']:
                        if unit_conversions[var][1]:
                            data[reg]=(xrds[var]*area_atm).sum(dim=('lat','lon'), skipna=True)
                        else:
                            data[reg]=(xrds[var]).weighted(area_atm).mean(dim=('lat','lon'), skipna=True)
                    elif var=='co2':
                        data[reg]=(xrds[var]).isel(lev=1).weighted(area_catm).mean(dim=('lat','lon'), skipna=True)
                        print(f'co2 {data[reg]}')
                        print(xrds[var])
                    elif var=='talk':
                        thk=xr.open_dataset(f'{local_path}/thkcello_Omon_EC-Earth3-CC_esm-pi-cdr-pulse_r1i1p1f1_gr_{syear}01-{syear}12.nc',engine='h5netcdf')
                        thkb=thk.thkcello.fillna(0)
                        data[reg]=(xrds[var]*thkb*areab).sum(dim=('i','j','lev'), skipna=True)
                    elif var=='alkalinization'  or var=='fco2removal2D':
                        data[reg]=(xrds[var]*area_pis).sum(dim=('y','x'), skipna=True)
                
                    
                    else:
                        print(var)
                        if unit_conversions[var][1]:
                            data[reg]=(xrds[var]*areab).sum(dim=('i','j'), skipna=True)
                        else:
                            data[reg]=(xrds[var]).weighted(areab).mean(dim=('i','j'), skipna=True)
                    print(data[reg])
                    print('--------------------------------')
                convert_units(data,var)
                if var=='alkalinization' or var == 'omegac' or var=='fco2removal2D':
                    timevalues=xrds.time_counter.values
                else:
                    timevalues=xrds.time.values
                print(timevalues)
                print(data['glob'].values)
                orig_unit=xrds[var].attrs['units']
                if vars[var][-3:]=='reg':
                    # create Dataset for given variables
                    outds = xr.Dataset(
                        {
                            "global": (["time"], data['glob'].values,{'units':units[vars[var]]}),
                            "US_EEZ": (["time"], data['US'].values,{'units':units[vars[var]]}),
                            "EU_EEZ": (["time"], data['EU'].values,{'units':units[vars[var]]}),
                            "CHINA_EEZ": (["time"], data['CHINA'].values,{'units':units[vars[var]]}),
                            },
                        coords={
                            "time": timevalues,
                            #"reference_time": pd.Timestamp("2015-01-01"),
                            }
                    )
                else:
                    outds = xr.Dataset(
                        {
                            "global": (["time"], data['glob'].values,{'units':units[vars[var]]}),
                            },
                        coords={
                            "time": timevalues,
                            #"reference_time": pd.Timestamp("2015-01-01"),
                            }
                    )
                
                reftime=pd.Timestamp("2015-01-01")
                outds['time'].encoding['dtype']=float
                outds['time'].encoding['units']=f"days since {reftime}"
                outds.attrs={'Source':'EC-Earth3-CC', 'Author':'Tommi Bergman','Institute':'Finnish Meteorological Institute', 'Project':'OceanNETs','Experiment':experiment_id, 'Created':f'{cyear}-{cmonth}-{cday} {chour}',
                            'original_variable':f'{var}','original unit':orig_unit, 'unit_conversion':f'{unit_conversions[var]}'}
                print(outds)
                print('write to disk')
                print(f'{outpath}/{outputname}')
                
                #outds.to_netcdf(f'{outpath}/{vars[var]}_testi_{syear}.nc',engine='h5netcdf')
                outds.to_netcdf(f'{outpath}/{outputname}',engine='h5netcdf')
                print(f'make a copy {exp},{var}')
                print(exp)
                if exp == 'C5AF':
                    print('copy initila files for other exp')
                    for exp1 in expnames.keys():
                        print(exp1)
                        print(exp1)
                        if exp == 'E40H' and int(syear) < 2080:
                            outpath=f'/Volumes/ONETs-SSD/{expnames[exp1]}'
                            outputname=f'{vars[var]}_{model_id}_{expnames[exp1]}_{time_range}_{version_id}.{syear}.nc'
                            outds.attrs={'Source':'EC-Earth3-CC', 'Author':'Tommi Bergman','Institute':'Finnish Meteorological Institute', 'Project':'OceanNETs','Experiment':f'{expnames[exp1]}', 'Created':f'{cyear}-{cmonth}-{cday} {chour}',
                                    'original_variable':f'{var}','original unit':orig_unit, 'unit_conversion':f'{unit_conversions[var]}'}
                            outds.to_netcdf(f'{outpath}/{outputname}',engine='h5netcdf')
                        elif  (exp == 'D40H' or exp == 'C40H') and int(syear) < 2040:
        
                            print(exp1)
                            outpath=f'/Volumes/ONETs-SSD/{expnames[exp1]}'
                            print(outpath)
                            outputname=f'{vars[var]}_{model_id}_{expnames[exp1]}_{time_range}_{version_id}.{syear}.nc'
                            print(outputname)
                            outds.attrs={'Source':'EC-Earth3-CC', 'Author':'Tommi Bergman','Institute':'Finnish Meteorological Institute', 'Project':'OceanNETs','Experiment':f'{expnames[exp1]}', 'Created':f'{cyear}-{cmonth}-{cday} {chour}',
                                    'original_variable':f'{var}','original unit':orig_unit, 'unit_conversion':f'{unit_conversions[var]}'}
                            outds.to_netcdf(f'{outpath}/{outputname}',engine='h5netcdf')
                        elif  exp == 'D30H'  and int(syear) < 2030:
        
                            print(exp1)
                            outpath=f'/Volumes/ONETs-SSD/{expnames[exp1]}'
                            print(outpath)
                            outputname=f'{vars[var]}_{model_id}_{expnames[exp1]}_{time_range}_{version_id}.{syear}.nc'
                            print(outputname)
                            outds.attrs={'Source':'EC-Earth3-CC', 'Author':'Tommi Bergman','Institute':'Finnish Meteorological Institute', 'Project':'OceanNETs','Experiment':f'{expnames[exp1]}', 'Created':f'{cyear}-{cmonth}-{cday} {chour}',
                                    'original_variable':f'{var}','original unit':orig_unit, 'unit_conversion':f'{unit_conversions[var]}'}
                            outds.to_netcdf(f'{outpath}/{outputname}',engine='h5netcdf')


                else:
                    #outds.to_netcdf(f'{outpath}/{vars[var]}_testi_{syear}.nc',engine='h5netcdf')
                    outds.to_netcdf(f'{outpath}/{outputname}',engine='h5netcdf')
                outds.close()
                xrds.close()
                if var=='talk':
                    thk.close()
    
def change_origin_Date(ds,date_string):
    # Extract the time values
    time_values = data.values

    # Change the reference date to "days since 2022-01-01"
    new_reference_date = pd.to_datetime('2015-01-01')
    time_values += (new_reference_date - pd.to_datetime('2000-01-01')).days

    # Update the time coordinate values
    #data['time'] = time_values
    #ds = xr.decode_cf(data)
    #ds.attrs['calendar'] = 'standard'
    #ds.attrs['units'] = f'days since 2015-01-01'
    return ds    
def date_string():
    now = datetime.now() # current date and time
    year = now.strftime("%Y")
    print("year:", year)

    month = now.strftime("%m")
    print("month:", month)

    day = now.strftime("%d")
    print("day:", day)
    hour = now.strftime("%H:%M")
    print("hour:", hour)
    return year,month,day,hour

if __name__ == "__main__":
    outputnames={'phos':'phos'}
    # Paths
    #exp=['C5AF','D30H','D40H','C40H','E40H']
    #exp=['C5AF']
    #exp=['D30H']
    #exp=['D40H']
    exp=['C40H']
    #exp=['E40H']

    #mask_eu_us_china=xr.open_dataset('../mask-creation/mask-eu-us-china.nc')
    #download_data_via_ssh(hostname, username, password, remote_path, local_path, 'areacello*nc')
    # Download data via SSH
    outpath=f'testout/'+exp[0]+'/'
    #for prefix in vars.keys():
    #    variable_name=prefix
    #prefix=var
    #filname = prefix + "_*mon*gr_20[3]*.nc"
    #print(local_path,filname)
    #download_data_via_ssh(hostname, username, password, remote_path, local_path, filname)

    # Read and save data in NetCDF format
    for ex in exp:
        localdatapath='/Volumes/ONETs-SSD/'
        local_path = localdatapath+'/'+ex+'/'
        isExist = os.path.exists(local_path)
        if not isExist:
            # Create a new directory because it does not exist
            os.makedirs(local_path)
        read_and_save_netcdf(local_path,ex)


   
    
