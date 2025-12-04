#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 17:16:54 2025

@author: vpremier
"""

import xarray as xr
import geopandas as gpd
import pandas as pd
from shapely.geometry import Point, box
import numpy as np
from pykrige.ok3d import OrdinaryKriging3D
from pykrige.uk3d import UniversalKriging3D
import rasterio
from utils import *
from rasterio.warp import reproject, Resampling



def create_grid(info):
    # get the information about the extent and resolution
    pixel_size = info['geotransform'][1]
    x_min = info['extent'][0] 
    x_max = info['extent'][2]
    y_min = info['extent'][1]
    y_max = info['extent'][3]
    eps = 0.0001

    # get the coordinates of the center of each cell.
    # An epsilon is added to consider always the last cell, that may be excluded 
    # due to floating errors
    xcoords = np.arange(x_min + pixel_size/2, 
                        x_max - pixel_size/2 + eps, pixel_size)
    ycoords = np.arange(y_min + pixel_size/2, 
                        y_max - pixel_size/2 + eps, pixel_size)

    gridx, gridy = np.meshgrid(xcoords, ycoords)
    
    return gridx, gridy 





def ds_to_df(ds, var_name):
    #Convert clipped data to DataFrames
    lon, lat = ds['x'].values, ds['y'].values
    values = ds[var_name].values
    
    # Create meshgrid of coordinates
    xx, yy = np.meshgrid(lon, lat)
    df = pd.DataFrame({
        'lon': xx.flatten(),
        'lat': yy.flatten(),
        var_name: values.flatten()
    })
    return df


def get_kriging(var, varname, date, stack_HR, info, DEM, mask, dem_path, sample = False):

    # select a given time step
    var_sel = var.sel(time=date)

    # Get the bounding box of stack_HR in its native CRS
    stack_hr_bounds = stack_HR.rio.bounds()
    stack_hr_geom = box(*stack_hr_bounds)

    # Convert to GeoDataFrame for reprojection
    stack_hr_gdf = gpd.GeoDataFrame(geometry=[stack_hr_geom], crs=info['projection'])
    # stack_hr_gdf["geometry"] = stack_hr_gdf.buffer(6000)  # Apply buffer of 5 km

    # Reproject the bounding box to match the CRS of pr and ta (EPSG:4326)
    stack_hr_gdf_4326 = stack_hr_gdf.to_crs("EPSG:4326")

    #Clip datasets to the extent of stack_HR
    var_clipped = var_sel.rio.clip(stack_hr_gdf_4326.geometry, stack_hr_gdf_4326.crs, all_touched =True)

    # get the point as a dataframe
    var_df = ds_to_df(var_clipped, varname)

    #Filter only points within stack_HR area
    # Convert merged DataFrame to GeoDataFrame
    gdf_points = gpd.GeoDataFrame(
        var_df, 
        geometry=gpd.points_from_xy(var_df.lon, var_df.lat),
        crs="EPSG:4326"
    )

    # Reproject gdf_points to your reference system
    gdf_points_utm = gdf_points.to_crs(info['projection'])



    # Open your reference xarray Dataset (e.g., era5land)
    ref = var_clipped  # your xarray dataset
    ref_crs = ref.rio.crs
    ref_res = ref.rio.resolution()
    ref_transform = ref.rio.transform()
    ref_width = ref.rio.width
    ref_height = ref.rio.height

    # Read the DEM and reproject/resample to match reference grid
    with rasterio.open(dem_path) as src:
        # Prepare output array
        dst_array = np.empty((ref_height, ref_width), dtype=np.float32)
    
        # Perform resampling using mean aggregation
        rppj_dem = reproject(
                    source=rasterio.band(src, 1),
                    destination=dst_array,
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=ref_transform,
                    dst_crs=ref_crs,
                    dst_resolution=ref_res,
                    resampling=Resampling.average  # This is key for mean aggregation
                )


        # Sample elevation values at these coordinates
        elevations = list(rppj_dem[0].flatten())

    # # Open DEM with rasterio
    # with rasterio.open(dem_path) as dem:
    #     # Extract (x, y) coordinates as a list of tuples
    #     coords = [(x, y) for x, y in zip(gdf_points_utm.geometry.x, gdf_points_utm.geometry.y)]
        
    #     # Sample elevation values at these coordinates
    #     elevations = [val[0] for val in dem.sample(coords)]
        


    # Add elevation values to GeoDataFrame
    gdf_points_utm["elevation"] = elevations
    
    if sample:
        gdf_points_utm = gdf_points_utm.sample(n=100, random_state=42)

    parameter = gdf_points_utm[varname]
    X = gdf_points_utm.geometry.x
    Y = gdf_points_utm.geometry.y
    Z = gdf_points_utm['elevation'].values
        

    gridx, gridy = create_grid(info)
    
    gridy = np.flip(gridy, axis=0)

    func = lambda x, y, z: z

    krig = UniversalKriging3D(
        X[~np.isnan(parameter)],
        Y[~np.isnan(parameter)],
        Z[~np.isnan(parameter)],
        parameter[~np.isnan(parameter)],
        variogram_model="linear",
        drift_terms=["functional"],
        functional_drift=[func])
    
    z, ss = krig.execute("points", gridx.flatten(), gridy.flatten(), DEM.flatten())

    T_interp = np.reshape(z,np.shape(DEM))
    T_interp[~mask] = np.nan
    
    return T_interp
    

