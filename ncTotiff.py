import netCDF4 as nc
from osgeo import gdal, osr, ogr
import os
import glob
import datetime
import pyresample
import matplotlib as mpl
import matplotlib.pyplot as plt
import proplot as pplt
from pylab import imread
import cartopy
import cartopy.crs as ccrs
import cmaps
from numpy import meshgrid, deg2rad, gradient, cos
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from geopandas import *
import rioxarray # for CRS management in Xarray
import rasterio as rio
import copy
from skimage import exposure
from skimage import morphology
def main_algorithm(nc_file):
    # 读取反射率波段数据
    geo_grp_file = nc_file.groups["geophysical_data"]
    geophysical_variables = geo_grp_file.variables
    Rrs_412 = geophysical_variables['Rrs_412'][:]
    Rrs_443 = geophysical_variables['Rrs_443'][:]
    Rrs_469 = geophysical_variables['Rrs_469'][:]
    Rrs_488 = geophysical_variables['Rrs_488'][:]
    Rrs_531 = geophysical_variables['Rrs_531'][:]
    Rrs_547 = geophysical_variables['Rrs_547'][:]
    Rrs_555 = geophysical_variables['Rrs_555'][:]
    Rrs_645 = geophysical_variables['Rrs_645'][:]
    Rrs_667 = geophysical_variables['Rrs_667'][:]
    Rrs_678 = geophysical_variables['Rrs_678'][:]
    Rrs_748 = geophysical_variables['Rrs_748'][:]
    Rrs_859 = geophysical_variables['Rrs_859'][:]
    # 对无效值进行掩膜
    Rrs_412 = Rrs_412.filled(np.nan) #-0.01
    Rrs_443 = Rrs_443.filled(np.nan)
    Rrs_469 = Rrs_469.filled(np.nan)
    Rrs_488 = Rrs_488.filled(np.nan)
    Rrs_531 = Rrs_531.filled(np.nan)
    Rrs_547 = Rrs_547.filled(np.nan)
    Rrs_555 = Rrs_555.filled(np.nan)
    Rrs_645 = Rrs_645.filled(np.nan)
    Rrs_667 = Rrs_667.filled(np.nan)
    Rrs_678 = Rrs_678.filled(np.nan)
    Rrs_748 = Rrs_748.filled(np.nan)
    Rrs_859 = Rrs_859.filled(np.nan)
    # 读取mask波段信息
    l2_flags = geophysical_variables['l2_flags'][:]
    checkLAND = True
    checkHISATZEN = True
    checkHIGLINT = True
    checkHISOLZEN = True
    checkCLDICE = True
    checkSTRAYLIGHT = True
    useFlags = 0
    useFlags += checkLAND * 2 ** (2 - 1)
    useFlags += checkHISATZEN * 2 ** (6 - 1)
    useFlags += checkHIGLINT * 2 ** (4 - 1)
    useFlags += checkHISOLZEN * 2 ** (13 - 1)
    useFlags += checkCLDICE * 2 ** (10 - 1)
    useFlags += checkSTRAYLIGHT * 2 ** (9 - 1)
    invalid1 = np.bitwise_and(l2_flags, useFlags) > 0
    # invalid2=((Rrs_412 < 0) | (Rrs_443 < 0) | (Rrs_469 < 0 ) | (Rrs_488 < 0) | (Rrs_531 < 0) | (Rrs_547 < 0)
    # | (Rrs_555 < 0) | (Rrs_645 < 0) | (Rrs_678 < 0) | (Rrs_748 < 0) | (Rrs_859 < 0))
    # bad_idx=invalid1 | invalid2
    Rrs_412[invalid1] = np.nan
    Rrs_443[invalid1] = np.nan
    Rrs_469[invalid1] = np.nan
    Rrs_488[invalid1] = np.nan
    Rrs_531[invalid1] = np.nan
    Rrs_547[invalid1] = np.nan
    Rrs_555[invalid1] = np.nan
    Rrs_645[invalid1] = np.nan
    Rrs_667[invalid1] = np.nan
    Rrs_678[invalid1] = np.nan
    Rrs_748[invalid1] = np.nan
    Rrs_859[invalid1] = np.nan
def save_to_tif(geo_extent, out_var, out_folder, out_var_name, band_num):
    maxlon = geo_extent['maxLon']
    minlon = geo_extent['minLon']
    maxlat = geo_extent['maxLat']
    minlat = geo_extent['minLat']
    # 分辨率计算
    N_Lat = out_var.shape[0]
    N_Lon = out_var.shape[1]
    Lon_Res = (maxlon - minlon) / (float(N_Lon) - 1)
    Lat_Res = (maxlat - minlat) / (float(N_Lat) - 1)
    # 创建.tif文件
    driver = gdal.GetDriverByName('GTiff')
    out_tif_name = os.path.join(out_folder, out_var_name + '.tif')
    out_tif = driver.Create(out_tif_name, N_Lon, N_Lat, band_num, gdal.GDT_Byte, options=["TILED=YES", "COMPRESS=LZW"]) #gdal.GDT_Byte , options=["TILED=YES", "COMPRESS=LZW"]
    # 设置影像的显示范围
    # -Lat_Res一定要是-的
    geotransform = (minlon, Lon_Res, 0, maxlat, 0, -Lat_Res)
    out_tif.SetGeoTransform(geotransform)
    out_tif.SetProjection('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')  # 给新建图层赋予投影信息
    # 数据写出
    if band_num == 1:
        out_tif.GetRasterBand(1).WriteArray(out_var)
    else:
        for i in range(band_num):
            out_tif.GetRasterBand(i + 1).WriteArray(out_var[:, :, i])
    out_tif.FlushCache() #将数据写入硬盘
    out_tif.BuildOverviews('average', [2, 4, 8, 16, 32]) #建立金字塔
    del out_tif #注意必须关闭tif文件
    print("L2 Process Finished -- %s " % out_var_name)

def swath_resampling(src_data: np.ma.array, src_lon: np.array, src_lat: np.array,
                     trg_lon: np.array, trg_lat: np.array, search_radius: float, fill_value: None):
    if len(trg_lon.shape) == 1:
        grid_def = pyresample.geometry.SwathDefinition(*np.meshgrid(trg_lon, trg_lat))
    else:
        grid_def = pyresample.geometry.SwathDefinition(lons=trg_lon, lats=trg_lat)

    swath_def = pyresample.geometry.SwathDefinition(lons=src_lon, lats=src_lat)
    """Use pyresample to resample all files to a uniform grid, then Average each grid cell"""
    result = pyresample.kd_tree.resample_nearest(swath_def, src_data, grid_def, epsilon=0.5,
                                                 fill_value=fill_value, radius_of_influence=search_radius)
    return result, grid_def
def process_coccolith(in_folder, out_folder):
    nc_files = glob.glob(os.path.join(in_folder, "*.nc"))
    for ncfile in nc_files:
        dirname, filename = os.path.split(ncfile)
        out_var_name = filename[0:14]
        nc_file = nc.Dataset(ncfile)
        navigation_variables = nc_file.groups["navigation_data"].variables
        lon = navigation_variables['longitude'][:]
        lat = navigation_variables['latitude'][:]

        if lon.mask.all() | isinstance(np.nanmin(lon), np.ma.core.MaskedConstant):
            print('L2 Process Failed   -- %s ************************ Masks all around' % out_var_name)
        else:
            classdata, rdata, gdata, bdata = main_algorithm(nc_file)
            minlat = np.min(lat)
            minlon = np.min(lon)
            maxlat = np.max(lat)
            maxlon = np.max(lon)
            x = np.arange(minlon, maxlon, 0.01)  # 1 km grid,
            y = np.arange(maxlat, minlat, -0.01)

            """left mask 对应 0< lon <180"""
            split_flag = (lon >= 0)
            left_lon = lon[split_flag]
            left_lat = lat[split_flag]
            left_var = classdata[split_flag]
            left_r = rdata[split_flag]
            left_g = gdata[split_flag]
            left_b = bdata[split_flag]

            right_lon = lon[~split_flag]
            right_lat = lat[~split_flag]
            right_var = classdata[~split_flag]
            right_r = rdata[~split_flag]
            right_g = gdata[~split_flag]
            right_b = bdata[~split_flag]
            #len_left_lon = len(left_lon)
            #len_right_lon = len(right_lon)

            if (len(left_lon) == 0) | (len(right_lon) == 0):
                var_re, grid = swath_resampling(classdata, lon, lat, x, y, 3000, np.nan)  # 1 km grid
                r_re, grid = swath_resampling(rdata, lon, lat, x, y, 3000, np.nan)
                g_re, grid = swath_resampling(gdata, lon, lat, x, y, 3000, np.nan)
                b_re, grid = swath_resampling(bdata, lon, lat, x, y, 3000, np.nan)
                rgbArray = np.dstack((r_re, g_re, b_re))
                extent_json = {
                    'maxLon': maxlon,
                    'minLon': minlon,
                    'maxLat': maxlat,
                    'minLat': minlat
                }
                save_to_tif(extent_json, var_re, out_folder, out_var_name, 1)
                #save_to_tif(extent_json, rgbArray, out_folder, out_var_name + '_RGB', 3)
                #extent_json['tifname'] = out_var_name
            else:
                left_maxlon = np.max(left_lon)
                left_minlon = np.min(left_lon)
                left_maxlat = np.max(left_lat)
                left_minlat = np.min(left_lat)
                left_x = np.arange(left_minlon, left_maxlon, 0.01)  # 1 km grid,
                left_y = np.arange(left_maxlat, left_minlat, -0.01)

                right_maxlon = np.max(right_lon)
                right_minlon = np.min(right_lon)
                right_maxlat = np.max(right_lat)
                right_minlat = np.min(right_lat)
                right_x = np.arange(right_minlon, right_maxlon, 0.01)  # 1 km grid,
                right_y = np.arange(right_maxlat, right_minlat, -0.01)

                if ((left_minlon - right_maxlon) > 300) | ((left_maxlon - right_minlon) > 300):
                    left_var_re, left_grid = swath_resampling(left_var, left_lon.filled(), left_lat.filled(), left_x,
                                                              left_y, 3000, np.nan)  # 1 km grid
                    left_r_re, left_grid = swath_resampling(left_r, left_lon.filled(), left_lat.filled(), left_x,
                                                              left_y, 3000, np.nan)
                    left_g_re, left_grid = swath_resampling(left_g, left_lon.filled(), left_lat.filled(), left_x,
                                                              left_y, 3000, np.nan)
                    left_b_re, left_grid = swath_resampling(left_b, left_lon.filled(), left_lat.filled(), left_x,
                                                              left_y, 3000, np.nan)
                    left_rgbArray = np.dstack((left_r_re, left_g_re, left_b_re))
                    left_extent_json = {
                        'maxLon': left_maxlon,
                        'minLon': left_minlon,
                        'maxLat': left_maxlat,
                        'minLat': left_minlat
                    }

                    right_var_re, right_grid = swath_resampling(right_var, right_lon.filled(), right_lat.filled(),
                                                                right_x, right_y, 3000, np.nan)  # 1 km grid
                    right_r_re, right_grid = swath_resampling(right_r, right_lon.filled(), right_lat.filled(),
                                                                right_x, right_y, 3000, np.nan)
                    right_g_re, right_grid = swath_resampling(right_g, right_lon.filled(), right_lat.filled(),
                                                                right_x, right_y, 3000, np.nan)
                    right_b_re, right_grid = swath_resampling(right_b, right_lon.filled(), right_lat.filled(),
                                                                right_x, right_y, 3000, np.nan)
                    right_rgbArray = np.dstack((right_r_re, right_g_re, right_b_re))
                    right_extent_json = {
                        'maxLon': right_maxlon,
                        'minLon': right_minlon,
                        'maxLat': right_maxlat,
                        'minLat': right_minlat
                    }

                    save_to_tif(left_extent_json, left_var_re, out_folder, out_var_name + '_1', 1)
                    save_to_tif(right_extent_json, right_var_re, out_folder, out_var_name + '_2', 1)
                    #save_to_tif(left_extent_json, left_rgbArray, out_folder, out_var_name + '_RGB_1', 3)
                    #save_to_tif(right_extent_json, right_rgbArray, out_folder, out_var_name + '_RGB_2', 3)
                    # left_extent_json['tifname'] = out_var_name + '_1'
                    # right_extent_json['tifname'] = out_var_name + '_2'

                else:
                    var_re, grid = swath_resampling(classdata, lon, lat, x, y, 3000, np.nan)  # 1 km grid
                    r_re, grid = swath_resampling(rdata, lon, lat, x, y, 3000, np.nan)
                    g_re, grid = swath_resampling(gdata, lon, lat, x, y, 3000, np.nan)
                    b_re, grid = swath_resampling(bdata, lon, lat, x, y, 3000, np.nan)
                    rgbArray = np.dstack((r_re, g_re, b_re))
                    extent_json = {
                        'maxLon': maxlon,
                        'minLon': minlon,
                        'maxLat': maxlat,
                        'minLat': minlat
                    }
                    #extent_json['tifname'] = out_var_name
                    save_to_tif(extent_json, var_re, out_folder, out_var_name, 1)
                    #save_to_tif(extent_json, rgbArray, out_folder, out_var_name + '_RGB', 3)
