#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 15:16:04 2024

@author: vpremier
"""
import pandas as pd
import xarray as xr
import glob
import os
import sys
import numpy as np
import geopandas as gpd
import time
import matplotlib.pyplot as plt

from tqdm import tqdm  # Progress bar (optional)

from correction import *
from utils import *
from loading import *

from kriging import get_kriging

    
    
def corection(sca_slice, df, HR_dates, status, w=5):
    # convert the time window in pd.Timedelta if working with a netcdf
    w=pd.Timedelta(str(w)+'D')
    # LOAD THE XARRAY DATASETS
    sca_slice = sca_slice.load()
    # status = status.assign_coords(lat=sca_slice.lat, lon=sca_slice.lon)
  
           
    # initialize t_min_melt and t_min_acc
    t_min_melt = 0
    t_min_acc = 0
    
    for date in df.index[1:]:

        ix = df.index.get_loc(date)
        print(ix,date)

    
        start_time = time.time()
        
        # t_min = np.minimum(t_min_melt, t_min_acc)
        
        HR_prev = sca_slice.sel(time=date-pd.Timedelta('1D'))['SCA']
        HR_curr = sca_slice.sel(time=date)['SCA']   
        
        # valid mask: snow in previous and current date
        valid_HR = (HR_prev <= 100) & (HR_curr <= 100)
        
        # accumulation: transition from snow free to snow
        acc_HR = valid_HR & np.logical_and(HR_prev==0, HR_curr==100)
       
        # melting: transition from snow to snow free
        melt_HR = valid_HR & np.logical_and(HR_prev==100, HR_curr==0)
    
    
        # get index of the previous accumulation and melting
        status_slice = status.sel(time=slice(df.index[0],date-pd.Timedelta('1D')))
        
        # Reverse the time axis to start from the end
        status_slice_rev = status_slice[::-1]
        
        #previous accumulation
        ix_prev_acc = (status_slice_rev == 1).argmax(dim='time', skipna=True)
        ix_prev_acc = len(status_slice.time) - ix_prev_acc - 1
        ix_prev_acc = xr.where((status_slice_rev == 1).sum(dim='time') > 0, ix_prev_acc, 0)
    
        t_min_acc = ix_prev_acc.min().compute().item()
            
            
        #previous melting
        ix_prev_melt = (status_slice_rev == 0).argmax(dim='time', skipna=True)
        ix_prev_melt = len(status_slice.time) - ix_prev_melt - 1
        ix_prev_melt =  xr.where((status_slice_rev == 0).sum(dim='time') > 0, ix_prev_melt, 0)
    
        t_min_melt = ix_prev_melt.min().compute().item()
    
        print("--- %s seconds ---" % (time.time() - start_time))
          
       
        # how much time from the last accumulation?
        delta_last_acc = ix + 1 - ix_prev_acc 
    
        
        sca_slice = correct_rec_acc(date, sca_slice, status, melt_HR, ix_prev_melt,
                            delta_last_acc, w, t_min_melt)
        
        sca_slice = correct_old_acc(date, sca_slice, status, HR_dates, melt_HR, 
                                    ix_prev_acc, ix_prev_melt, delta_last_acc, w, t_min_melt)
    
        sca_slice = correct_rec_melt(date, sca_slice, status, acc_HR, ix_prev_acc,
                              delta_last_acc, w, t_min_acc)
        
        sca_slice = correct_old_melt(date, sca_slice, status, HR_dates, acc_HR, 
                                    ix_prev_acc, delta_last_acc, w, t_min_acc)
        
        print("--- %s seconds ---" % (time.time() - start_time))

  


    del status, status_slice, status_slice_rev, ix_prev_acc, ix_prev_melt
    
    return sca_slice
    



def run_harmonization(config_path):
    
    config = load_config(config_path)

    dirname = config['dirname']
    outdir = config['outdir']
    os.makedirs(outdir, exist_ok=True)
    
    dem_path = config['DEM_path']
    temp_dir = config['temp_dir']
    era5_dir = config['era5_dir']


    hy_xxxx = config['hy_xxxx']
    basin = config['catchment']
    pixel_ratio = config["pixel_ratio"]

    subbasin = None
    
    # load SCA
    sca_path = glob.glob(dirname + os.sep + '*' + hy_xxxx + '*.nc')[0]
    csv_path = sca_path.replace('.nc','.csv')
    
    outname = os.path.join(outdir, os.path.basename(sca_path).replace('.nc','_harm.nc'))

    if os.path.exists(outname):
        print('File %s has been already created' %outname)
    
    else:
        
        DEM, mask_shape_cut, info = get_mask_info(dem_path, subbasin, 
                                                  resType='HR', 
                                                  pixel_ratio=pixel_ratio)
        
        stack_HR, df, epsg_code = upload_sca(sca_path, dem_path, subbasin)
        stack_HR = stack_HR.rio.write_crs(epsg_code['projection'], inplace=True)  

        ta = load_micromet(temp_dir, hy_xxxx)
        
        era5 = load_era5land(era5_dir, hy_xxxx)
        pr = era5.tp
        pr = pr.rio.write_crs("EPSG:4326", inplace=True)  

        pr_reprojected = pr.rio.reproject_match(stack_HR)


        # dates with a HR acquisition
        HR_dates = (df.filter(like="HR_fileName") != 'M').any(axis=1)
        
        # status
        status = (ta['t2m'] < 274.15) & (pr_reprojected > 1)
        # status = buffer(status, n = 10)
        # status = set_all_to_one_per_timestep(status)
   
        status = status.load()

        stack_HR = stack_HR.load()
        sca_corr = corection(stack_HR, df, HR_dates, status)
        
                 
        array =  sca_corr.transpose("y", "x","time").SCA.values

        save_nc(outname, array, info, df, 'SCA', 'percentage', 
                scale=1, dtype = 'int32', complevel = 9)
        
        # plot 
        sca_harm_ts = sca_corr.where(sca_corr<=100).mean(dim=['x','y']).SCA.values       
        df = pd.read_csv(csv_path, parse_dates=[0], index_col=0)
        
        ta_mean = ta.mean(dim=['x','y']).t2m.values
        pr_mean = pr_reprojected.where(pr_reprojected>0.1).mean(dim=['x','y']).values
        status_ts = status.mean(dim=['x','y']).values
    
    
    
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 6), sharex=True)
        
        # === Upper subplot: SCA final and harmonized ===
        # Create the twin axis for status bars
        ax1_twinx = ax1.twinx()
        
        # Draw status bars FIRST, with low zorder and transparent alpha
        ax1_twinx.bar(
            df.index, status_ts*100, width=1,
            color='cyan', alpha=0.4, label='% area with snowfall',
            zorder=0
        )
        
        ax1_twinx.set_ylim([-5, 105])
        # Now draw the SCA lines with higher zorder so they appear on top
        ax1.plot(df.index, df['SCA final'], color='black', label='SCA Final', zorder=3)
        ax1.plot(df.index, sca_harm_ts, color='tab:orange', label='SCA Harmonized', zorder=3)
        
        ax1.set_ylabel('SCA [%]')
        ax1.set_title('SCA Time Series')
        # Merge legends from both y-axes
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax1_twinx.get_legend_handles_labels()
        ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right')
        
        ax1.grid(True, linestyle='--', alpha=0.3, zorder=1)
        ax1.set_ylim([-5, 105])
    
    
        
    
        
        # === Lower subplot: Temperature and Precipitation ===
        # Left y-axis → Temperature
        ax2.plot(df.index, ta_mean, color='tab:red', label='Temperature (K)')
        ax2.set_ylabel('Mean Temperature (K)', color='tab:red')
        ax2.tick_params(axis='y', labelcolor='tab:red')
        ax2.axhline(274.15, color='gray', linestyle='--', label='1°C (274.15 K)')  # Zero degree line
        
        # Right y-axis → Precipitation (bar)
        ax3 = ax2.twinx()
        ax3.bar(df.index, pr_mean, width=1, color='tab:blue', alpha=0.5, label='Precipitation (mm)')
        ax3.set_ylabel('Mean Precipitation (mm)', color='tab:blue')
        ax3.tick_params(axis='y', labelcolor='tab:blue')
        
        # Merge legends from both y-axes
        lines1, labels1 = ax2.get_legend_handles_labels()
        lines2, labels2 = ax3.get_legend_handles_labels()
        ax2.legend(lines1 + lines2, labels1 + labels2, loc='upper right')
        
        # Formatting
        ax2.grid(True, linestyle='--', alpha=0.4)
        ax3.grid(False)
        ax3.set_xlabel('Date')
        
        plt.tight_layout()
        plt.savefig(outname.replace('.nc','.png'))
        plt.savefig(outname.replace('.nc','.svg'))





if __name__ == "__main__":   
    
    if len(sys.argv) != 2:
        print("Usage: python main.py path_to_config.json")
    else:
        config_path = sys.argv[1]
        config_path = r'/home/vpremier/Documents/git/dailySCA_harmonization/config.json'
        start_time = time.time()
    
        run_harmonization(config_path)
    
        end_time = time.time()
        elapsed = end_time - start_time
        elapsed_min = int(elapsed // 60)
        elapsed_sec = int(elapsed % 60)
    
        config = load_config(config_path)
        
        print("\nThe harmonization workflow run succefully.")
        print(f"Execution time: {elapsed_min} minutes and {elapsed_sec} seconds")
        
        
     

                





