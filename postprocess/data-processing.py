# Processing for data portal

import numpy as np
import xarray as xr
#import paramiko
import os
#import netCDF4 as nc
import glob
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import pandas as pd 
# def download_data_via_ssh(hostname, username, password, remote_path, local_path, file_pattern):
#     # Connect to the remote server
#     ssh_client = paramiko.SSHClient()
#     ssh_client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
#     print(ssh_client)
#     try:
#         if password != None:
#             ssh_client.connect(hostname, username=username, password=password)
#         else:
#             ssh_client.connect(hostname, username=username)
        
#         # SFTP connection to download the file
#         sftp = ssh_client.open_sftp()
#         #remote_file = os.path.join(remote_path, fname)
#         #local_file = os.path.join(local_path, fname)
#         # List files matching the pattern on the remote server
#         print(sftp)
#         remote_files = sftp.listdir(remote_path)
#         print(remote_path)
#         #print(remote_files)
#         matching_files = [f for f in remote_files if glob.fnmatch.fnmatch(f, file_pattern)]

#         #sftp.get(remote_file_path,'areacello_Ofx_EC-Earth3-CC_esm-pi-cdr-pulse_r1i1p1f1_gr.nc')
#         # Download each matching file
#         for remote_file in matching_files:
#             remote_file_path = os.path.join(remote_path, remote_file)
#             local_file_path = os.path.join(local_path, remote_file)
#             sftp.get(remote_file_path, local_file_path)
#             print(f'Downloaded: {remote_file_path} to {local_file_path}')

#         sftp.close()
        
#         print(f'Data downloaded from {remote_file} to {local_file}')
        
#     except Exception as e:
#         print(f'Error: {e}')
#     finally:
#         ssh_client.close()

vars={
    #global means
    'fgco2': 'fgco2_glob',
    'fgoae': 'fgoae_glob',
    'fgdcr': 'fgdcr_glob',
    'tosga': 'tosga',
    'dissic': 'dissic_glob',
    #'talk': 'talk_glob',
    'xco2atm': 'xco2atmga',
    'tas': 'tasga',
    'nbp': 'nbp_glob',
    'gpp': 'gpp_glob',
    'cVeg': 'cveg_glob',
    'cSoil': 'csoil_glob',
    #region means:
    'spco2': 'spco2_reg',
    'dissicos': 'dissicos_reg',
    'talkos': 'talkos_reg',
    'phos': 'pH_reg',
#    'omegaa': 'omegaa_reg', #not available
    'omegac': 'omegac_reg',
    #3d  integrated
    'intdic':'dissic_glob',
    'talk':'talk_glob'

    }
units = {
        'fgco2_glob': 'molC m-2 s-1',
        'fgoae': 'Pmol TA yr-1', #monthly
        'fgdcr':'Pg C yr-1',
        'tasga':'degC',
        'tas':'degC',
        'dissic_glob': 'mol m-3',
        'xco2': 'ppm',
        'nbp_glob': 'Pg C yr-1',
        'gpp_glob':'Pg C yr-1',
        'cSoil_glob':'PgC ',
        'cVeg_glob':'PgC',
        'spco2': 'uAtm',
        'dissicos_reg': 'mol m-3',
        'talkos':'mol m-3',
        'ph':'-',
        'omeagac':'-',
        'dissic_glob':'',
        'talk':'P mol',

}
expnames={
'C5AF':'esm-ssp534-over',
'D30H':'esm-ssp534-over-2030high',
'D40H':'esm-ssp534-over-2040high',
'E40H':'esm-ssp534-over-2040high-term',
'C40H':'esm-ssp534-over-2040high-dcr'}

n_c=12.011
unit_conversions = {'cSoil':1e-12,  #kgm-2 -> PgC m-2
                'cVeg':1e-12,   #kgm-2 -> PgC m-2
                'gpp':1e-12*3600*24,    #kgm-2s-1 -> PgCm-2day-1
                'nbp':1e-12*3600*24,    #kgm-2s-1 -> PgCm-2day-1
                'fgco2':1e-12*3600*24,    #kgm-2s-1 -> PgCm-2day-1
                'dissicos':1,    #molm-3
                'talkos':1,    #molm-3
                'spco2':1e6/101325,    # Pa -> uAtm
                'xco2':1e6,    # mol/mol -> umol/mol or ppm  
                'ph':1, # no conversion unit of 1
                'omeagac':1, #no conversion unit of 1
                'fgdcr':1e-12*3600*24, # XXXX
                'fgoae':1e-15*3600*24, # xxxxx->Pmol m-3 day-1
                'tas':-273.15 # K -> C 
                }  
def convert_units(data):
    for i in data:
        print(i)

def read_and_save_netcdf(local_path,fname):
    # Read data from the downloaded NetCDF file
    #matching_files = [f for f in local_files if glob.fnmatch.fnmatch(f, fname)]
    #fname='D40H_*pisces*_2D*'
    model_id='EC-Earth'
    year,month,day,hour=date_string()
   
    version_id=f'v{year}{month}{day}'
    
    experiment_id=expnames[exp[0]]
    outpath=f'/Volumes/ONETs-SSD/{experiment_id}/2D'
    isExist = os.path.exists(outpath)
    time_range='2030-2100'
    syear='2050'
    if not isExist:
        # Create a new directory because it does not exist
        os.makedirs(outpath)
    
    local_files = glob.glob(os.path.join(local_path, fname))
    area=xr.open_dataset('areacello_Ofx_EC-Earth3-CC_esm-pi-cdr-pulse_r1i1p1f1_gr.nc',engine='h5netcdf')
    # Download each matching file
    #testout='talkos_'
    var=fname.split('_')[0]
    print(var,fname,fname.strip('_'))
    testout=var+'_'
    #'fgco2_'
    country_mask=xr.open_dataset('countrymask.b.nc')
    country_mask=xr.open_dataset('countrymask_v2.b.nc')
    #print(country_mask)
    
    # mask_EU=country_mask.where(country_mask.mask==30)
    # mask_US=country_mask.where(country_mask.mask==20)
    # mask_CHINA=country_mask.where(country_mask.mask==10)
    # print(mask_US) 
    outputname=f'{variable_name}_{model_id}_{experiment_id}_{time_range}_{version_id}.nc'

    areab=area.areacello.fillna(0)
    print(local_files)
    print(local_path)
    print(fname)
    regions=['glob','EU','US','CHINA']
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
        plt.figure()
        data={}
        for reg in regions:
            if reg =='CHINA':
                match=10
            elif reg =='US':
                match=20
            elif reg =='EU':
                match=30
            else:
                match=-999
            if reg== 'glob':
                data[reg]=xrds[var].mean(dim=('i','j'), skipna=True)

            else:
                data[reg]=xrds[var].where(country_mask.mask==match).mean(dim=('i','j'), skipna=True)
        # create DS
        timevalues=xrds.time.values
        print(timevalues)
        print(data['glob'].values)
        convert_units(data)
        outds = xr.Dataset(
            {
                "global": (["time"], data['glob'].values),
                "US_EEZ": (["time"], data['US'].values),
                "EU_EEZ": (["time"], data['EU'].values),
                "CHINA_EEZ": (["time"], data['CHINA'].values),
                },
            coords={
                "time": timevalues,
                #"reference_time": pd.Timestamp("2015-01-01"),
                }
        )
        reftime=pd.Timestamp("2015-01-01")
        outds['time'].encoding['units']=f"days since {reftime}"
        outds.attrs={'Source':'EC-Earth3-CC', 'Author':'Tommi Bergman','Institute':'Finnish Meteorological Institute', 'Project':'OceanNETs','Experiment':experiment_id, 'Created':f'{year}-{month}-{day} {hour}'}
        print(outds)
        outds.to_netcdf(f'{outpath}/{var}_testi_{syear}.nc')
        #xrds.dpco2.where(country_mask.mask==10).mean(dim=('y','x'), skipna=True).plot()
        #plt.figure()
        #xrds.dpco2.where(country_mask.mask==20).mean(dim=('y','x'), skipna=True).plot()
        #xrds.phos.where(country_mask.mask==20).mean(dim=('i','j'), skipna=True).plot()
        #plt.figure()
        #xrds.dpco2.where(country_mask.mask==30).mean(dim=('y','x'), skipna=True).plot()
        #xrds.phos.where(country_mask.mask==30).mean(dim=('i','j'), skipna=True).plot()
        #print(areatest)
        #plt.show()
        #(((xrds.fgco2*areab).sum(dim=('i','j')))*3600*24*30).to_netcdf(local_path+'/'+outti )
   
    #plt.show()
    
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
    exp=['D30H','C40H']
    localdatapath='.'
    local_path = localdatapath+'/'+exp[0]+'/'
    isExist = os.path.exists(local_path)
    if not isExist:
        # Create a new directory because it does not exist
        os.makedirs(local_path)
    
    #mask_eu_us_china=xr.open_dataset('../mask-creation/mask-eu-us-china.nc')
    #download_data_via_ssh(hostname, username, password, remote_path, local_path, 'areacello*nc')
    # Download data via SSH
    outpath=f'testout/'
    for prefix in vars.keys():
        variable_name=prefix
        filname = prefix + "_*mon*gr_20[3]*.nc"
        print(local_path,filname)
        #download_data_via_ssh(hostname, username, password, remote_path, local_path, filname)

        # Read and save data in NetCDF format
        read_and_save_netcdf(local_path,filname)


   
    
