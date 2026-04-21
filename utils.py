#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  4 20:36:30 2021

@author: vpremier
"""

from osgeo import gdal, osr, ogr
import numpy as np
import os
import re
import xarray as xr
import rioxarray as rxr
from datetime import datetime
import pickle
import geopandas as gpd
import shapely.vectorized
from affine import Affine
from rasterio.features import geometry_mask
import json

def load_config(config_path):
    with open(config_path, 'r') as f:
        config = json.load(f)
    return config

def save_tif(outdir, fname, info, array):
    """Save a tiff
    
    Parameters
    ----------
    outdir : str
        path to the output folder
    fname : str
        output file name
    info : dict
        dictionary containing image metadata  (see open_image)
    array : np.array
        array with the image
 
    """

    fileName_output = outdir + os.sep + fname + '.tif'
    #Create the map
    out = gdal.GetDriverByName('GTiff').Create(fileName_output,info['X_Y_raster_size'][0],
                              info['X_Y_raster_size'][1],1,gdal.GDT_Float32)     
    outband = out.GetRasterBand(1)
    
    # Set the geographic information
    out.SetGeoTransform(info['geotransform'])
    out.SetProjection(info['projection'])
    outband.WriteArray(array)
        
    outband.FlushCache()
    out = None
    
    
def reproj_point(x, y, srIn, srOut):
    """Reproject a point into a defined coordinate system.
    
    Parameters
    ----------
    x : float
        x-coordinate
    y : float
        y-coordinate
    srIn : osgeo.osr.SpatialReference
        spatial reference of the input coordinate system
    srOut : osgeo.osr.SpatialReference
        spatial reference of the output coordinate system
        
    Returns
    -------
    (x, y) : tuple
        the transformed coordinates
    """
    
    epsgList = ['4326', '31287']
    # create a geometry from coordinates
    point = ogr.Geometry(ogr.wkbPoint)
    
    if int(gdal.__version__[0]) >= 3 and srIn.GetAttrValue("AUTHORITY", 1) in epsgList:
        # GDAL 3 changes axis order: https://github.com/OSGeo/gdal/issues/1546
        point.AddPoint(y, x)
    else:
        point.AddPoint(x, y)
    
    coordTransform = osr.CoordinateTransformation(srIn,srOut)
    
    # transform point
    point.Transform(coordTransform)
    
    if int(gdal.__version__[0]) >= 3 and srOut.GetAttrValue("AUTHORITY", 1) in epsgList:
        (y, x) = point.GetX(), point.GetY()
    else:
        (x, y) = point.GetX(), point.GetY()
        
    return (x,y)


def load_ts(path):
    """
    Get the extension of the data   
    """
    
    ext = os.path.basename(path).split('.')[-1]
    if ext == 'nc':
        data = xr.open_dataset(path)
        data = data.load()
    elif ext == 'pickle':
        data = pickle.load(open(path, 'rb'))
 
    return data


def open_image(image_path):
    """Opens an image and reads its metadata.
    
    Parameters
    ----------
    image_path : str
        path to an image
    
    Returns
    -------
    image : osgeo.gdal.Dataset
        the opened image
    information : dict
        dictionary containing image metadata    
    """
    
    image = gdal.Open(image_path)
    
    if image is None:
        print('could not open ' + image_path)
        return
        
    cols = image.RasterXSize
    rows = image.RasterYSize
    geotransform = image.GetGeoTransform()
    proj = image.GetProjection()
    minx = geotransform[0]
    maxy = geotransform[3]
    maxx = minx + geotransform[1] * cols
    miny = maxy + geotransform[5] * rows
    X_Y_raster_size = [cols, rows]
    extent = [minx, miny, maxx, maxy]
    information = {}
    information['geotransform'] = geotransform
    information['extent'] = extent
    information['X_Y_raster_size'] = X_Y_raster_size
    information['projection'] = proj

    return image, information



def save_nc(outname, array, info, df, varname, unit, scale=1, dtype = 'int32',
            complevel = 9):
   

    time = df.index
    reference_time = df.index[0]
    
        
 
    if type(array) is xr.core.dataset.Dataset:
        da = array
       
    else:
        # in this way, the coordinate is the upper left corner
        nx = np.shape(array)[1]; ny = np.shape(array)[0]
        
        x = info['geotransform'][0] + info['geotransform'][1]/2 + \
            info['geotransform'][1]*np.arange(nx)
        y = info['geotransform'][3] + info['geotransform'][5]/2 + \
            info['geotransform'][5]*np.arange(ny)
        
 
    
        da = xr.DataArray(
            name = varname,
            data=array,
            dims=["y","x",  "time"],
            coords=dict(
                x=(["x"], x),
                y=(["y"], y),
                time=time,
                reference_time=reference_time,
            ),
            attrs=dict(
                units=unit,
            ),
        )
    
        da = da.transpose("time", "y", "x")

    if os.path.exists(outname):
        print("The output netcdf already exists")
        da_0 = load_ts(outname)
        da = xr.combine_by_coords([da_0, da])
        os.remove(outname)
        
    #get EPSG
    srs = osr.SpatialReference()
    srs.ImportFromWkt(info['projection'])
    
    # da = da.assign_coords({"crs": info['projection']})



    da.rio.write_crs("epsg:" + srs.GetAttrValue('AUTHORITY',1), 
                     inplace=True).rio.set_spatial_dims(
                         x_dim="x",
                         y_dim="y",inplace=True).rio.write_coordinate_system(inplace=True)
                         
    da.rio.write_coordinate_system("epsg:" + srs.GetAttrValue('AUTHORITY',1))
    
    encode = {
        varname: {
            'zlib': True,
            'complevel': complevel,
            'dtype': dtype
        }
    }
    
    # Only add scale_factor if different from 1
    if scale != 1:
        encode[varname]['scale_factor'] = scale
    
    da.to_netcdf(outname, encoding=encode)
    
    

def get_scene_extent(i1, i2):
    """Calculates the extent of one image in another one, eventually with
    different resolutions. This function is used to retrieve the extent of the
    high resolution Landsat or Sentinel 2 scenes in the lower resolution
    reference products AVHRR, MODIS, etc.
    
    Parameters
    ----------
    i1 : str
        path to a .tif image
    i2 : str
        path to a .tif image

    Returns
    -------
    x_tl : float
        the x-coordinate of the top left corner point of i1 in i2
    y_br : float
        the y-coordinate of the bottom right corner point of i1 in i2
    x_br : float
        the x-coordinate of the bottom right corner point of i1 in i2
    y_tl : float
        the y-coordinate of the top left corner point of i1 in i2
    """
    
    # read the two images
    i1_ds, i1_info = open_image(i1)
    i2_ds, i2_info = open_image(i2)
    
    i1_ds = None
    i2_ds = None
    
    # difference in i2 grid points of the x-coordinate of the top
    # left corner between i1 and i2
    diff_x = ((i1_info['geotransform'][0] - i2_info['geotransform'][0]) /
            i2_info['geotransform'][1])
    
    # x-coordinate of the top left corner of i1 relative to i2
    x_tl = int(diff_x)
    
    # number of i1 grid points i1 is shifted in x-direction with respect to the
    # nearest i2 grid point
    x_shift = int((round((diff_x % 1), 5) * i2_info['geotransform'][1]) /
               i1_info['geotransform'][1])
    
    # difference in i2 grid points of the y-coordinate of the top
    # left corner between i1 and i2
    diff_y = ((i1_info['geotransform'][3] - i2_info['geotransform'][3]) /
            i2_info['geotransform'][5])
    
    # y-coordinate of the top left corner of i1 relative to i2
    y_tl = int(diff_y)
    
    # number of i1 grid points i1 is shifted in y-direction with respect to the
    # nearest i2 grid point
    y_shift = int((round((diff_y % 1), 5) * i2_info['geotransform'][5]) /
               i1_info['geotransform'][5])

    # get the indices of the bottom right corner of i1 relative to i2
    x_br = int(np.ceil(x_tl + ((x_shift + i1_info['X_Y_raster_size'][0]) *
                               i1_info['geotransform'][1]) /
               i2_info['geotransform'][1]))
    y_br = int(np.ceil(y_tl + ((y_shift + i1_info['X_Y_raster_size'][1]) *
                               i1_info['geotransform'][5]) /
               i2_info['geotransform'][5]))
    
    # check if the extent of i1 is contained in i2
    if x_br > i2_ds.RasterXSize:
        raise IndexError('{} domain exceeds {} domain. '
                         ''.format(os.path.basename(os.path.normpath(i1)),
                                   os.path.basename(os.path.normpath(i2))))
        return
    
    return x_tl, y_br, x_br, y_tl, x_shift, y_shift




def dateFromFileName(string):
    """Return a date from a file name
    
    Parameters
    ----------
    string : str
        path to an image
    
    Returns
    -------
    date : datetime
        the date  
    """
    match = re.search(r'\d{8}', string)
    date = datetime.strptime(match.group(), '%Y%m%d').date()
    
    return(date)
 
    
    
def get_new_path(df, newdir_S2, newdir_L8):
    """
    The function looks for the correct subfolder when changing the main
    directory that contains the data and creates a new dataframe with the 
    correct paths.
    """
    df_new = df.copy()
    
    def find_LC08(inpath):
#        l8_id = r'LC08_L1TP_\d{6}_\d{8}_\d{8}_01_T1'
        l8_id = r'L[CE]0[78]_L1TP_\d{6}_\d{8}_\d{8}_\d{2}_T[12]'
        match = re.search(l8_id, inpath)
        return match.group(0)
    
    def find_S2(inpath):

        s2_id = r'S2[AB]_MSIL1C_\d{8}T\d{6}_N\d{4}_R\d{3}_.*_\d{8}T\d{6}'
        match = re.search(s2_id, inpath)
        try: 
            return match.group(0)
        except:
            s2_id = r'S2[AB]_OPER_MSI_L1C_TL_[A-Z][A-Z][A-Z]__V\d{8}T\d{6}_A\d{6}_T\d{2}[A-Z][A-Z][A-Z]'
            match = re.search(s2_id, inpath)   
        return match.group(0)
    
    if "S2A_fileName" in df:
        df_new["S2A_fileName"] = df_new["S2A_fileName"].apply(
                                    lambda x: x if x == 'M' else newdir_S2)
        # df_new["S2A_fileName"] = df_new["S2A_fileName"].apply(
        #                             lambda x: x if x == 'M' else 
        #                             os.path.abspath(os.path.join(x ,"../..")))       
    if "LC8_fileName" in df:
        df_new["LC8_fileName"] = df_new["LC8_fileName"].apply(
                                    lambda x: x if x == 'M' else newdir_L8)

        # df_new["LC8_fileName"] = df_new["LC8_fileName"].apply(
        #                             lambda x: x if x == 'M' else 
        #                             os.path.abspath(os.path.join(x ,"../..")))           
    return df_new



def compareLR(df, mask, date, countLR, pixel_ratio, cut_nan = False):
    """
    This function combines two different LR images.
    RULE: take an average vale for the SCF 
    """    
    # get the indices of valid pixels
    (j0, jend, i0, iend) = get_mask_indices(mask, cut_nan)

    
    i0_lr = int(i0/pixel_ratio)
    iend_lr = int(iend/pixel_ratio)+1
    j0_lr = int(j0/pixel_ratio)
    jend_lr = int(jend/pixel_ratio)+1
    
    open_first = True
    n_stack=[]
    for col in countLR:

        if open_first:
            img_LR = gdal.Open(df[col].loc[date])
            scfMap = img_LR.GetRasterBand(1).ReadAsArray(i0_lr, j0_lr, 
                                             iend_lr-i0_lr, jend_lr-j0_lr)
            img_LR = None
            
            open_first = False
            
        else:
            
            img_LR = gdal.Open(df[col].loc[date])
            scfMap_tmp = img_LR.GetRasterBand(1).ReadAsArray(i0_lr, j0_lr, 
                                             iend_lr-i0_lr, jend_lr-j0_lr)
            img_LR = None
            
            c_diff = scfMap != scfMap_tmp
            c_nodata = np.logical_xor(scfMap>100,scfMap_tmp>100)
            
            valid_1 = np.logical_and(scfMap>=0,scfMap<=100)
            valid_2 = np.logical_and(scfMap_tmp>=0,scfMap_tmp<=100)
            c_valid = np.logical_and(valid_1,valid_2)
            
            valid_1 = None
            valid_2 = None
            
            n_stack.append(c_valid)
            
            nValid = np.sum(np.array(n_stack),axis=0)
        
            # mean SCF if there are more valid values
            averageMap = (scfMap[c_valid]*nValid[c_valid] + \
                          scfMap_tmp[c_valid])/(nValid[c_valid]+1)
                       
            # if they are different but one is valid, keep the valid value 
            map1_tmp = scfMap[np.logical_and(c_diff,c_nodata)]
            map2_tmp = scfMap_tmp[np.logical_and(c_diff,c_nodata)] 
            map1_tmp[map1_tmp>100] = 0
            map2_tmp[map2_tmp>100] = 0
            
            scfMap[c_diff] = 205
            scfMap[c_valid] = averageMap   
            scfMap[np.logical_and(c_diff,c_nodata)] = map1_tmp + map2_tmp
         
            averageMap = None
            map1_tmp = None
            map2_tmp = None
            
    return scfMap


def scf2binary(snowMap, snowT):
    """
    Convert a HR map from scf to binary map.
    """
    snowMap[np.logical_and(snowMap>snowT,snowMap<=100)] = 100
    snowMap[snowMap <= snowT] = 0
    snowMap = snowMap.astype(np.int16)
    
    return snowMap


   
def compare2HR(path_scene1, path_scene2, mask, snowT, cut_nan = False):
    """
    This function combines two different HR images.
    RULE: if one is no data keep the information of the other,
    if same information keep it, otherwise if different info
    substitute with no data value  
    """
    
    # get the indices of valid pixels
    (j0, jend, i0, iend) = get_mask_indices(mask, cut_nan)
        
    
    
    img_1= gdal.Open(path_scene1)
    snowMap_1 = img_1.GetRasterBand(1).ReadAsArray(i0, j0, iend-i0+1, jend-j0+1)
    img_1 = None
    snowMap_1 = scf2binary(snowMap_1, snowT)        
            
    img_2 = gdal.Open(path_scene2)
    snowMap_2 = img_2.GetRasterBand(1).ReadAsArray(i0, j0, iend-i0+1, jend-j0+1)
    img_2 = None
    snowMap_2 = scf2binary(snowMap_2, snowT)

    snowMap = np.zeros((np.shape(snowMap_2)[0],
                        np.shape(snowMap_2)[1]), dtype=np.int16)
    snowMap[:] = 255
    
    # logical conditions
    c_diff = snowMap_1 != snowMap_2
    c_nodata = np.logical_xor(snowMap_1>100,snowMap_2>100)
    
    # if they have the same value, keep it
    snowMap[snowMap_1 == snowMap_2] = snowMap_1[snowMap_1 == snowMap_2] 
    
    # if they are different, replace with no data
    snowMap[c_diff] = 205
    
    # if they are different but one is valid, keep the valid value 
    map1_tmp = snowMap_1[np.logical_and(c_diff,c_nodata)]
    map2_tmp = snowMap_2[np.logical_and(c_diff,c_nodata)] 
    map1_tmp[map1_tmp>100] = 0
    map2_tmp[map2_tmp>100] = 0
    
    snowMap[np.logical_and(c_diff,c_nodata)] = map1_tmp + map2_tmp
    
    map1_tmp = None
    map2_tmp = None
    
    return snowMap


def normalized_difference(img_path_1, img_path_2, mask, cut_nan = False):
    """
    The function computes the normalized diffference between
    two images (bands)
    """
    # get the indices of valid pixels
    (j0, jend, i0, iend) = get_mask_indices(mask, cut_nan)
    
    # open the images
    img_1, info_1 = open_image(img_path_1)
    img_2, info_2 = open_image(img_path_2)
    
#    if info_1['X_Y_raster_size'] != info_2['X_Y_raster_size']:
#        # get the correct extent if the images are not the same
#        extent = get_scene_extent(dem_path, img_path_1)
#        
#        if (extent[0] == 0) and (extent[3] == 0):
#            # there is a shift
#            band_1 = img_1.ReadAsArray()
#            band_2 = img_2.ReadAsArray()
#            
#        else:
#            band_1 = img_1.ReadAsArray(extent[0],extent[3],
#                                       extent[2]-extent[0],
#                                       extent[1]-extent[3])
#            band_2 = img_2.ReadAsArray(extent[0],extent[3],
#                                       extent[2]-extent[0],
#                                       extent[1]-extent[3])  
#    else:
        
    band_1 = img_1.ReadAsArray(i0, j0, iend-i0+1, jend-j0+1)
    band_2 = img_2.ReadAsArray(i0, j0, iend-i0+1, jend-j0+1)   

    img_1 = None
    img_2 = None
    
    norm_diff = (band_1 - band_2) / (band_1 + band_2)
    
    band_1 = None
    band_2 = None
    
    return norm_diff



def mask_nv(array, mask, codes=[205,210,255,254]):
    
    nonvalid = np.isin(array, codes)
    array[nonvalid] = 205
    
    array[mask] = 255
    
    return array



def get_LR_mask(mask, pixel_ratio=20, nv_thres=40):
    
    # rows and columns of the new mask
    nrows = int(np.shape(mask)[0]/pixel_ratio)
    ncols = int(np.shape(mask)[1]/pixel_ratio)
    
    #initialize array
    mask_LR = np.zeros(shape=(nrows, ncols), dtype=bool)
    
    # iterate over rows
    y_j = 0
    x_i = 0
    for j in range(0, nrows):                                    
        # reset column counter
        x_i = 0             
        # iterate over columns
        for i in range(0, ncols):

            # read the slice of the scene matching the current
            # estimator pixel
            data_ij = mask[y_j:y_j + pixel_ratio,x_i: x_i + pixel_ratio]

            # check how many pixels are not valid
            nv_sum = np.sum(data_ij)

            if nv_sum >= nv_thres:
                # if the number of masked pixels exceed the threshold, set as False
                mask_LR[j, i] = True

            # advance column counter by number of high resolution pixels
            # contained in one low resoution pixels
            x_i += pixel_ratio
        
        # advance row counter by number of high resolution pixels
        # contained in one low resoution pixels
        y_j += pixel_ratio
                
    return mask_LR



def get_LR_dem(dem, pixel_ratio=20):
    
    # rows and columns of the new mask
    nrows = int(np.shape(dem)[0]/pixel_ratio)
    ncols = int(np.shape(dem)[1]/pixel_ratio)
    
    #initialize array
    dem_LR = np.zeros(shape=(nrows, ncols))
    
    # iterate over rows
    y_j = 0
    x_i = 0
    for j in range(0, nrows):                                    
        # reset column counter
        x_i = 0             
        # iterate over columns
        for i in range(0, ncols):

            # read the slice of the scene matching the current
            # estimator pixel
            dem_LR[j, i] = np.nanmean(dem[y_j:y_j + pixel_ratio,x_i: x_i + pixel_ratio])

            # advance column counter by number of high resolution pixels
            # contained in one low resoution pixels
            x_i += pixel_ratio
        
        # advance row counter by number of high resolution pixels
        # contained in one low resoution pixels
        y_j += pixel_ratio
                
    return dem_LR


def get_max(array, pixel_ratio=20):
    
    # rows and columns of the new mask
    nrows = int(np.shape(array)[0]/pixel_ratio)
    ncols = int(np.shape(array)[1]/pixel_ratio)
    
    #initialize array
    output = np.zeros(shape=(nrows, ncols))
    
    # iterate over rows
    y_j = 0
    x_i = 0
    for j in range(0, nrows):                                    
        # reset column counter
        x_i = 0             
        # iterate over columns
        for i in range(0, ncols):

            # read the slice of the scene matching the current
            # estimator pixel
            output[j, i] = np.nanmax(array[y_j:y_j + pixel_ratio,x_i: x_i + pixel_ratio])

            # advance column counter by number of high resolution pixels
            # contained in one low resoution pixels
            x_i += pixel_ratio
        
        # advance row counter by number of high resolution pixels
        # contained in one low resoution pixels
        y_j += pixel_ratio
                
    return output


def get_CP(SCFMap, scf_range_dic, CP_scf_dic, pixel_ratio=20):
    
    SCFMapHR = SCFMap.repeat(pixel_ratio, axis=0).repeat(pixel_ratio, axis=1) 
    
    for key in scf_range_dic:
        cp = CP_scf_dic[key]
        
        if key=='0_20':
            c1 = np.logical_and( SCFMapHR >= scf_range_dic[key][0], 
                                SCFMapHR <= scf_range_dic[key][1])
        else:
                
            c1 = np.logical_and(SCFMapHR > scf_range_dic[key][0], 
                                SCFMapHR <= scf_range_dic[key][1])
        cpMap = cp.copy()
        cpMap[~c1] = np.nan
        
    cpMap[SCFMapHR==0] = 0
    cpMap[SCFMapHR==100] = 1
    
    SCFMapHR = None
    cp = None
    
    return cpMap



def interpolate_scf(stack_SCF, stack_SCF_min, stack_SCF_max, df):
    # convert to xarray to perform the linear interpolation
    SCF = xr.DataArray(stack_SCF, dims=('x', 'y','time'), 
                       coords={'time': df.index})
    SCF = SCF.where(SCF!=205)
    SCF_interp = SCF.interpolate_na(dim='time',method='linear')
    
    del SCF
    
    SCF_min = xr.DataArray(stack_SCF_min, dims=('x', 'y','time'), 
                           coords={'time': df.index})
    SCF_min = SCF_min.where(SCF_min!=205)
    SCF_min_interp = SCF_min.interpolate_na(dim='time',method='linear')
    
    del SCF_min
    
    SCF_max = xr.DataArray(stack_SCF_max, dims=('x', 'y','time'), 
                           coords={'time': df.index})
    SCF_max = SCF_max.where(SCF_max!=205)
    SCF_max_interp = SCF_max.interpolate_na(dim='time',method='linear')
    
    del SCF_max
    
    SCF_interp = np.array(SCF_interp)
    SCF_min_interp = np.array(SCF_min_interp)
    SCF_max_interp = np.array(SCF_max_interp)
    
    valid_mask = np.logical_and(SCF_interp>=0,SCF_interp<=100)
            
    c_min = np.logical_and(np.less(SCF_interp,SCF_min_interp),valid_mask)  
    c_max = np.logical_and(np.greater(SCF_interp,SCF_max_interp),valid_mask)
    
    SCF_interp[c_min] = SCF_min_interp[c_min]
    SCF_interp[c_max] = SCF_max_interp[c_max]
    
    return SCF_interp, SCF_min_interp, SCF_max_interp  


def get_window(mask, ncol, nrow, pixel_ratio):
    
    rows = np.shape(mask)[0]
    cols = np.shape(mask)[1]
    
    rows_lr = int(rows/pixel_ratio)
    cols_lr = int(cols/pixel_ratio)
    
    
    row_slice = int(np.ceil(rows_lr/nrow))
    col_slice = int(np.ceil(cols_lr/ncol))
    
    lr_mask = np.ones(shape=(rows_lr, cols_lr)) * 0
    hr_mask = np.ones(shape=(rows, cols)) * 0

    nwindow = 1
    for j in range(0,rows_lr,row_slice):  
        for i in range(0, cols_lr, col_slice):
 
            lr_mask[j:j + row_slice,i: i + col_slice]+=nwindow
            hr_mask[j*pixel_ratio:j*pixel_ratio + pixel_ratio*row_slice,
                    i*pixel_ratio: i*pixel_ratio + pixel_ratio*col_slice]+=nwindow
            nwindow+=1

    return hr_mask, nwindow
  

def get_new_info(mask, info):
    
    # Find indices where mask is False (valid region)
    indices = np.where(~mask)
    
    # Ensure there are valid indices before accessing
    if indices[0].size > 0 and indices[1].size > 0:
        imin = indices[1].min()
        imax = indices[1].max()
        
        jmin = indices[0].min()
        jmax = indices[0].max()
    else:
        # Handle the case where all values are masked (e.g., set defaults)
        imin, imax, jmin, jmax = None, None, None, None
    
    
    xmin_new = info['geotransform'][0] + imin*info['geotransform'][1] 
    xmax_new = info['geotransform'][0] + (imax+1)*info['geotransform'][1] 

    ymax_new = info['geotransform'][3] + jmin*info['geotransform'][5] 
    ymin_new = info['geotransform'][3] + (jmax+1)*info['geotransform'][5] 
    
    info_new = info.copy()
    
    info_new['geotransform']= (xmin_new, info['geotransform'][1], 0,
                                   ymax_new, 0, info['geotransform'][5])
    info_new['extent'] = [xmin_new, ymin_new, xmax_new, ymax_new]
    info_new['X_Y_raster_size']= [imax-imin+1, jmax-jmin+1]
    
    return info_new


def get_mask_info(dem_path, subbasin, resType='HR', pixel_ratio=10):    
    # read the dem
    DEM, info = open_image(dem_path)
    DEM = DEM.ReadAsArray()
    mask = np.logical_or(np.isnan(DEM),(DEM <= 0.001))


    indices, mask_shape =  get_mask_indices(mask, info, subbasin, 
                                            resType=resType, pixel_ratio=pixel_ratio)
    (i0, iend, j0, jend) = indices
    mask_shape_cut = mask_shape[j0:jend+1,i0:iend+1]
    
    info_new = get_new_info(~mask_shape, info)
    DEM_new = DEM[j0:jend+1,i0:iend+1]
    
    return DEM_new, mask_shape_cut, info_new



def extract_ds_info(dataset, epsg):

    # Extract coordinate and raster size information
    x_coords = dataset.coords['lon']  # Replace 'lon' with the name of the longitude coordinate
    y_coords = dataset.coords['lat']  # Replace 'lat' with the name of the latitude coordinate
    
    # Extract raw numerical values for geotransform calculation
    x_min = x_coords[0].item()
    x_res = ((x_coords[-1] - x_coords[0]) / (len(x_coords) - 1)).item()
    y_max = y_coords[0].item()
    y_res = ((y_coords[-1] - y_coords[0]) / (len(y_coords) - 1)).item()
    
    geotransform = (x_min - x_res/2,
                    x_res, 
                    0.0, 
                    y_max - y_res/2, 
                    0.0, 
                    y_res)
    
    # Extract extent
    extent = [
        x_coords[0].item() - x_res/2,  # xmin
        y_coords[-1].item() + y_res/2,  # ymin (last because latitude is usually descending)
        x_coords[-1].item() - x_res/2,  # xmax
        y_coords[0].item() - y_res/2   # ymax
    ]

    # Determine raster size
    X_Y_raster_size = [len(x_coords), len(y_coords)]

    # Define the projection
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(epsg)  # Assuming WGS84
    projection = srs.ExportToWkt()

    # Create the dictionary
    result = {
        'geotransform': geotransform,
        'extent': extent,
        'X_Y_raster_size': X_Y_raster_size,
        'projection': projection
    }
    
    return result


def open_ds(dirname, basin, varname, dsname, date_start = None, date_end = None):
    
    dir_mask = dirname + basin + os.sep + 'forcings'
    
    if varname in ['pr', 'ta']:
        dir_var = dir_mask
    else:
        dir_var = dir_mask.replace('forcings', 'parameters')
        
    if basin == 'Alpenrhein':
        mask224 = xr.open_mfdataset(os.path.join(
                        dir_mask + os.sep + '224', 'my_mask.nc'))
        mask227 = xr.open_mfdataset(os.path.join(
                        dir_mask + os.sep + '227', 'my_mask.nc'))
        
        var224 = xr.open_mfdataset(os.path.join(
                        dir_var + os.sep + '224', varname + '*.nc'))
        var227 = xr.open_mfdataset(os.path.join(
                        dir_var + os.sep + '227', varname + '*.nc'))
    
        var224 = var224[dsname].where(mask224.area.values == 1)
        var227 = var227[dsname].where(mask227.area.values == 1)
 
        var_xr = var224.combine_first(var227)
        
    else:
        
        mask = xr.open_mfdataset(os.path.join(dir_mask, 'my_mask.nc'))
        if varname == 'elv':      
            var = xr.open_mfdataset(os.path.join(dir_var, varname + '.nc'))
        else:
            var = xr.open_mfdataset(os.path.join(dir_var, varname + '*.nc'))
        var_xr = var[dsname].where(mask.area.values == 1)
    
    if varname in ['pr', 'ta']:
        outvar = var_xr.sortby('time')
        outvar = outvar.resample(time='D', closed='right', label='left').mean()
        outvar = outvar.sel(time=slice(date_start, date_end))
    else:
        outvar = var_xr
            
    return outvar


def get_basin_mask(shp_path, ds_nc, factor = 1):
    
    catchment_outline = gpd.read_file(shp_path) 
    bbox = catchment_outline.bounds.iloc[0]
    
    
    final_pixel_size_x = (ds_nc.coords['lon'][-1] - ds_nc.coords['lon'][-2]).values.item()/factor
    final_pixel_size_y = -(ds_nc.coords['lat'][-1] - ds_nc.coords['lat'][-2]).values.item()/factor


    # new spatial coordinates
    xmin = bbox.minx + final_pixel_size_x/2
    xmax = bbox.maxx - final_pixel_size_x/2
    # segno invertito perchè res y è negativa
    ymin = bbox.miny - final_pixel_size_y/2
    ymax = bbox.maxy + final_pixel_size_y/2
    
    # get the coordinates
    coord_x = np.linspace(start=xmin, 
                          stop=xmax,
                          num=len(ds_nc.coords['lon'])*factor)
    coord_y = np.linspace(start=ymin, 
                          stop=ymax,
                          num=len(ds_nc.coords['lat'])*factor)
    

    # mask of the basin
    xx, yy = np.meshgrid(coord_x, coord_y)
    basin_mask = shapely.vectorized.contains(catchment_outline.dissolve().geometry.item(), xx, yy)
    basin_mask = np.flip(basin_mask)
    basin_mask = np.flip(basin_mask, axis=1)

    return basin_mask



def get_mask_indices(mask, info, shape, resType='HR', pixel_ratio=20):
    
    if shape:
        # Create the affine transformation from the geotransform
        affine_transform = Affine(info['geotransform'][1],
                                  info['geotransform'][2],
                                  info['geotransform'][0],
                                  info['geotransform'][4],
                                  info['geotransform'][5],
                                  info['geotransform'][3])
        dims = info['X_Y_raster_size'][1],info['X_Y_raster_size'][0]
        mask = geometry_mask(
                        shape, 
                        transform=affine_transform, 
                        invert=True,  # True means inside the shape is True
                        out_shape=dims
                    )
        
        # convert to LR-> at least 1 valid to be valid 
        mask_LR = get_LR_mask(mask, pixel_ratio=pixel_ratio, nv_thres=1)
        
        # get back to HR
        mask_HR = np.repeat(mask_LR, pixel_ratio, axis=0)
        mask_HR = np.repeat(mask_HR, pixel_ratio, axis=1)
        
        if resType=='HR':
            mask_shape = mask_HR
        else:
            mask_shape = mask_LR
            
        j0 = int(np.where(mask_shape.any(axis=1))[0][0])
        jend = int(np.where(mask_shape.any(axis=1))[0][-1])
        i0 = int(np.where(mask_shape.any(axis=0))[0][0])
        iend = int(np.where(mask_shape.any(axis=0))[0][-1])
    else:
        j0 = 0
        jend = np.shape(mask)[0]-1
        i0 = 0
        iend = np.shape(mask)[1]-1

    indices = (i0, iend, j0, jend)
    mask_shape = ~mask
    return indices, mask_shape

