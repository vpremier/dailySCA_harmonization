#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  7 15:09:40 2024

@author: vpremier
"""
import pandas as pd
import xarray as xr
import numpy as np
import time
import os



# accumulation ok, correct false melting
def correct_rec_acc(date, stack_HR, status, 
                   melt_HR, ix_prev_melt,
                   delta_last_acc, w, t_min, old_acc=10):
    
    # changes: exclude date and ceil for accumulation

    # mask the recent accumulation
    rec_acc = ((status.sel(time=date)==1) & (delta_last_acc<old_acc))

    # check the dates within a temporal window (of +- 5 days)
    stack_HR_window = stack_HR.sel(time=slice(date-w,date+w))['SCA']
    stack_HR_window = stack_HR_window.sel(time=stack_HR_window['time'] != date)
    
    # ensure that the two datasets have the same dimensions
    rec_acc_3d = rec_acc.broadcast_like(stack_HR_window)
    
    # mask the snow array: pixels that have experienced a recent accumulation
    snow_array = stack_HR_window.where(rec_acc_3d)
    
    # compute the most frequent label in the given time window
    snow_occ = (snow_array.where(melt_HR) == 100).sum(dim='time')/snow_array.sizes['time']
    


   
    # correct melting
    HR_curr = stack_HR.sel(time=date)['SCA'].copy()
    
    
    # replace pixels for the date with the most frequent label   
    condition = (melt_HR & rec_acc) & np.isfinite(snow_occ.round())
    stack_HR.loc[dict(time=date)] = xr.where(condition, 
                                              np.ceil(snow_occ).astype(int) * 100, 
                                              HR_curr)
    
    # replace the pixels from the date of the last ablation until 
    # date with snow free (when the most frequent label is snow free)
    
    # mask current date is snow free
    mask_sf = stack_HR.sel(time=date)['SCA'] == 0
         
    # Use apply_ufunc to replace each value with its time index
    stack_HR_slice = stack_HR.sel(time=slice(stack_HR.time[t_min],date))
    
    time_indices = t_min + ((
                        stack_HR_slice.time - stack_HR.time[t_min]
                                             )/(60*60*24*10**9)).astype(int)
    
    
    # condition to be verified: date is after the previous melting event
    # and pixels that are snow free for the current date
    condition = ((time_indices >= ix_prev_melt) & (melt_HR & mask_sf & rec_acc
                                                   ).broadcast_like(stack_HR_slice))
    
    stack_HR.loc[dict(time=slice(stack_HR.time[t_min], date))] = xr.where(~condition, 
                              stack_HR_slice,
                              0)
    return stack_HR




       
def correct_old_acc(date, stack_HR, status, HR_dates,
               melt_HR, ix_prev_acc, ix_prev_melt, 
               delta_last_acc, w, t_min, old_acc=10):
    
    # mask old accumulation
    old_acc = ((status.sel(time=date)==1) & (delta_last_acc>=old_acc)) 
    
    # select only the indices that correspond to a HR acquisition
    snow_array = stack_HR.sel(time=HR_dates[HR_dates].index)  
    
    nearest_date = snow_array.time.sel(time=date, method="nearest").values
    target_idx = snow_array.time.get_index("time").get_loc(nearest_date)
    
    
    # filter up to current date (t)
    # snow_array = snow_array.sel(time=slice(stack_HR.time[0],date))
    snow_array = snow_array.isel(time=slice(max(target_idx-5,0),target_idx))
    
    # consider up to the last 5 HR images
    # snow_array = snow_array.isel(time=slice(-5, None))
    
    # filter dates after the previous accumulation
    time_indices = t_min + ((
                            snow_array.time - stack_HR.time[t_min]
                            )/(60*60*24*10**9)).astype(int)
    
    
    condition = (time_indices>=ix_prev_acc) & (old_acc & melt_HR).broadcast_like(snow_array)
    snow_array_fltd = snow_array.where(condition)
    # occurences
    snow_occ = (snow_array_fltd == 100).sum(dim='time')/(snow_array_fltd.count(dim='time'))  
    
    
    # condition to be satisfied: the date is after the last previous accumulation
    # snow_array_fltd = snow_array.where(time_indices>=ix_prev_acc)
    
    # old_acc_3d = old_acc.broadcast_like(snow_array_fltd)
    
    
    # snow_array_fltd = snow_array_fltd.where(old_acc_3d)

    # occurences
    # snow_occ = (snow_array_fltd.where(melt_HR) == 100).sum(
    #                     dim='time')/(snow_array_fltd.where(melt_HR).sum(dim='time')/100)

    #correct melting
    HR_curr = stack_HR.sel(time=date)['SCA'].copy()
    condition = (melt_HR & old_acc) & np.isfinite(snow_occ.round())
    stack_HR.loc[dict(time=date)] = xr.where(condition, 
                                             snow_occ.round() * 100, 
                                             HR_curr)
    
 
    mask_sf = stack_HR.sel(time=date)['SCA'] == 0
    
    
    # new xarray with the time indices
    stack_HR_slice = stack_HR.sel(time=slice(stack_HR.time[t_min], date))
    
    time_indices = t_min + ((stack_HR_slice.time - stack_HR.time[t_min])/(60*60*24*10**9)).astype(int)
    
    condition = ((time_indices >= ix_prev_melt) & (melt_HR & mask_sf & old_acc).broadcast_like(stack_HR_slice))

    stack_HR.loc[dict(time=slice(stack_HR.time[t_min], date))] = xr.where(~condition, 
                                                                          stack_HR_slice,
                                                                          0)
                                                    
    return stack_HR
 



# melting ok, correct false accumulation
def correct_rec_melt(date, stack_HR, status, 
                     acc_HR, ix_prev_acc,
                     delta_last_acc, w, t_min, old_acc=10):
    
    # recent accumulation, full season. Check also what happens later
    rec_melt = ((status.sel(time=date)==0) & (delta_last_acc<old_acc)) 
 
    # check the dates within a temporal window (of +- 5 days)
    stack_HR_window = stack_HR.sel(time=slice(date-w,date+w))['SCA']
    stack_HR_window = stack_HR_window.sel(time=stack_HR_window['time'] != date)
    
    # ensure that the two datasets have the same dimensions
    rec_melt_3d = rec_melt.broadcast_like(stack_HR_window)
    
    snow_array = stack_HR_window.where(rec_melt_3d)
 
    snow_occ = (snow_array.where(acc_HR) == 100).sum(dim='time')/snow_array.sizes['time']
   
    # very recent accumulation
    snow_occ = snow_occ.where(~(delta_last_acc < 5), 1)
    
    #correct melting
    HR_curr = stack_HR.sel(time=date)['SCA'].copy()
    
    # replace pixels for the date with the most frequent label   
    condition = (acc_HR & rec_melt) & np.isfinite(snow_occ.round())
    stack_HR.loc[dict(time=date)] = xr.where(condition, 
                                             snow_occ.round() * 100, 
                                             HR_curr)

    mask_snow = stack_HR.sel(time=date)['SCA'] == 100
    
    # Use apply_ufunc to replace each value with its time index
    stack_HR_slice = stack_HR.sel(time=slice(stack_HR.time[t_min],date))
    
    time_indices = t_min + ((stack_HR_slice.time - stack_HR.time[t_min])/(60*60*24*10**9)).astype(int)
    
    
    # condition to be verified: date is after the previous melting event
    # and pixels that are snow free for the current date
    condition = ((time_indices >= ix_prev_acc) & (acc_HR &  mask_snow & rec_melt).broadcast_like(stack_HR_slice))
    
    stack_HR.loc[dict(time=slice(stack_HR.time[t_min], date))] = xr.where(~condition, 
                              stack_HR_slice,
                              100)
    return stack_HR
 

def correct_old_melt(date, stack_HR, status, HR_dates, acc_HR, 
                     ix_prev_acc, delta_last_acc, w, t_min, old_acc=10):
    # old accumulation
    old_melt = ((status.sel(time=date)==0) & (delta_last_acc>=old_acc)) #.astype(int) 
    
    
    # select only the indices that correspond to a HR acquisition
    snow_array = stack_HR.sel(time=HR_dates[HR_dates].index)
    
    nearest_date = snow_array.time.sel(time=date, method="nearest").values
    target_idx = snow_array.time.get_index("time").get_loc(nearest_date)
    snow_array = snow_array.isel(time=slice(max(target_idx-5,0),target_idx))

    # 
    # filter up to current date (t)
    # snow_array = snow_array.sel(time=slice(stack_HR.time[0],date))
    # consider up to the last 5 HR images
    # snow_array = snow_array.isel(time=slice(-5, None))
    # filter from the previous accumulation
 
    # new xarray with the time indices
    time_indices = t_min + ((snow_array.time - stack_HR.time[t_min])/(60*60*24*10**9)).astype(int)
    

    condition = (time_indices>=ix_prev_acc) & (old_melt & acc_HR).broadcast_like(snow_array)
    snow_array_fltd = snow_array.where(condition)
    # occurences
    snow_occ = (snow_array_fltd == 100).sum(dim='time')/(snow_array_fltd.count(dim='time'))  
    
    # old_melt_3d = old_melt.broadcast_like(snow_array_fltd)
    
    # snow_array_fltd = snow_array_fltd.where(old_melt_3d)

    

                            
    # snow_occ = (snow_array_fltd.where(acc_HR) == 100).sum(
    #                     dim='time')/(snow_array_fltd.where(acc_HR).sum(
    #                     dim='time')/100)   
 
    #correct melting
    HR_curr = stack_HR.sel(time=date)['SCA'].copy()
    condition = (acc_HR & old_melt) & np.isfinite(snow_occ.round())
    stack_HR.loc[dict(time=date)] = xr.where(condition, 
                                             snow_occ.round() * 100, 
                                             HR_curr)
 
    mask_snow = stack_HR.sel(time=date)['SCA'] == 100
    
    # Use apply_ufunc to replace each value with its time index
    stack_HR_slice = stack_HR.sel(time=slice(stack_HR.time[t_min],date))
    
    time_indices = t_min +  ((stack_HR_slice.time - stack_HR.time[t_min])/(60*60*24*10**9)).astype(int)


    # condition to be verified: date is after the previous melting event
    # and pixels that are snow free for the current date
    condition = ((time_indices >= ix_prev_acc) & (acc_HR & mask_snow & old_melt).broadcast_like(stack_HR_slice))
    
    stack_HR.loc[dict(time=slice(stack_HR.time[t_min], date))] = xr.where(~condition, 
                              stack_HR_slice,
                              100)       

    return stack_HR
    





