#!/usr/bin/env python
import time
import traceback
import glob
import os, random, shutil
import sys
import argparse
import pandas as pd

import numpy as np
import netCDF4 as nc
from scipy import ndimage
import pyresample
from osgeo import gdal, osr, gdal_array
from skimage import morphology

import multiprocessing as mp
import matplotlib.pyplot as plt


def main_algorithm(nc_file):
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
    Rrs_1640 = geophysical_variables['Rrs_1640'][:]

    l2_flags = geophysical_variables['l2_flags'][:]

    checkLAND = True
    checkHISATZEN = True
    checkHIGLINT = True
    checkHISOLZEN = True
    checkCLDICE = True
    checkSTRAYLIGHT = True
    checkCOCCOLITH = True

    useFlags = 0
    useFlags += checkLAND * 2 ** (2 - 1)
    useFlags += checkHISATZEN * 2 ** (6 - 1)
    useFlags += checkHIGLINT * 2 ** (4 - 1)
    useFlags += checkHISOLZEN * 2 ** (13 - 1)
    useFlags += checkCLDICE * 2 ** (10 - 1)
    # useFlags += checkCOCCOLITH * 2 ** (11 - 1)
    useFlags += checkSTRAYLIGHT * 2 ** (9 - 1)

    invalid1 = np.bitwise_and(l2_flags, useFlags) > 0

    "mask掉rrc小于0的像元 饱和像元"
    invalid2 = ((Rrs_412 < 0.0015) | (Rrs_443 < 0) | (Rrs_469 < 0) | (Rrs_488 < 0) \
                     | (Rrs_531 < 0) | (Rrs_547 < 0) | (Rrs_555 < 0) | (Rrs_645 < 0) \
                     | (Rrs_667 < 0) | (Rrs_678 < 0) | (Rrs_748 < 0))

    "计算云阴影及噪声指数"
    nRrs443 = Rrs_443 / Rrs_488
    nRrs488 = Rrs_488 / Rrs_488
    nRrs555 = Rrs_555 / Rrs_488
    nCSI = nRrs488 - nRrs443 - (nRrs555 - nRrs443) * 0.5
    invalid3 = (nCSI >= 0.12)

    "计算泥沙浑浊水体指数"
    nRrs645 = Rrs_645 / Rrs_667
    nRrs667 = Rrs_667 / Rrs_667
    nRrs678 = Rrs_678 / Rrs_667
    nSI = nRrs667 - nRrs645 - (nRrs678 - nRrs645) * 0.5
    # nSI = nRrs667 - nRrs645 - (nRrs678 - nRrs645) * (667 - 645) / (678 - 645)
    invalid4 = (nSI <= -0.4)

    nRrs469 = Rrs_469 / (Rrs_469 + Rrs_555 + Rrs_645)
    nRrs555 = Rrs_555 / (Rrs_469 + Rrs_555 + Rrs_645)
    nRrs645 = Rrs_645 / (Rrs_469 + Rrs_555 + Rrs_645)
    normal = nRrs555 - (nRrs469 + (nRrs645 - nRrs469) * (555 - 469) / (645 - 469))
    invalid5 = (normal <= 0.005)

    "计算cie坐标xy"
    X = 2.7689 * Rrs_748 + 1.7517 * Rrs_678 + 1.1302 * Rrs_667
    Y = 1.0000 * Rrs_748 + 4.5907 * Rrs_678 + 0.0601 * Rrs_667
    Z = 0.0000 * Rrs_748 + 0.0565 * Rrs_678 + 5.5943 * Rrs_667
    xx = X / (X + Y + Z)
    yy = Y / (X + Y + Z)
    "模型提取"
    Y1_diff = yy - (4.809313 * np.power(xx, 2) - 3.095812 * xx + 0.835671)  # final_2_area1
    Y2_diff = yy - (4.904038 * np.power(xx, 2) - 3.575923 * xx + 0.986154)  # final_2_area2

    bad_idx2_ = invalid2 | invalid3 | invalid4  ## 1km
    bad_idx2 = ndimage.binary_dilation(bad_idx2_, structure=np.ones((3, 3)))

    Y1_diff[bad_idx2] = -999
    Y2_diff[bad_idx2] = -999

    bad_idx1 = invalid1 | invalid5  ## 1km
    Y1_diff[bad_idx1] = -999
    Y2_diff[bad_idx1] = -999

    bloom_idx1 = (xx <= 0.333333) & (Y1_diff >= 0)
    bloom_idx2 = (xx > 0.333333) & (Y2_diff >= 0)

    """np.int8 (-128, 127) np.unit8 (0, 255)"""
    classdata = np.zeros(Rrs_645.shape, dtype='uint8')

    classdata[bloom_idx1] = 255
    classdata[bloom_idx2] = 255
    classdata[invalid1] = 100

    if np.sum(invalid1) > 0:
        Rrs_443[invalid1] = np.nan
        Rrs_469[invalid1] = np.nan
        Rrs_488[invalid1] = np.nan
        Rrs_531[invalid1] = np.nan
        Rrs_555[invalid1] = np.nan
        Rrs_645[invalid1] = np.nan
        Rrs_667[invalid1] = np.nan
        Rrs_678[invalid1] = np.nan
        Rrs_748[invalid1] = np.nan
        Rrs_1640[invalid1] = np.nan

    if np.sum(invalid2) > 0:
        Rrs_443[invalid2] = np.nan
        Rrs_469[invalid2] = np.nan
        Rrs_488[invalid2] = np.nan
        Rrs_531[invalid2] = np.nan
        Rrs_555[invalid2] = np.nan
        Rrs_645[invalid2] = np.nan
        Rrs_667[invalid2] = np.nan
        Rrs_678[invalid2] = np.nan
        Rrs_748[invalid2] = np.nan
        Rrs_1640[invalid2] = np.nan

    nan_idx = np.where(np.isnan(Rrs_555))
    valid_idx = np.where(~np.isnan(Rrs_555))

    lgRrs555 = np.zeros_like(Rrs_555)
    lgRrs488 = np.zeros_like(Rrs_488)
    lgRrs443 = np.zeros_like(Rrs_443)

    Rrs_555 = np.ma.masked_where(Rrs_555 == 0, Rrs_555)
    Rrs_488 = np.ma.masked_where(Rrs_488 == 0, Rrs_488)
    Rrs_443 = np.ma.masked_where(Rrs_443 == 0, Rrs_443)
    lgRrs555[valid_idx] = np.log(Rrs_555[valid_idx]/0.01)/np.log(1/0.01)
    lgRrs488[valid_idx] = np.log(Rrs_488[valid_idx]/0.01)/np.log(1/0.01)
    lgRrs443[valid_idx] = np.log(Rrs_443[valid_idx]/0.01)/np.log(1/0.01)
    lgRrs555[nan_idx] = np.nan
    lgRrs488[nan_idx] = np.nan
    lgRrs443[nan_idx] = np.nan

    rgbArray = get_rgb(lgRrs555, lgRrs488, lgRrs443)

    return classdata, rgbArray


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

def normalize(arr):
    arr_min = np.nanmin(arr)
    arr_max = np.nanmax(arr)
    norm = (arr - arr_min) / (arr_max - arr_min)
    return norm

def stretch_data(data, num_stddev):
    mean = np.mean(data)
    std_range = np.std(data) * num_stddev
    new_min = max(mean - std_range, np.min(data))
    new_max = min(mean + std_range, np.max(data))
    clipped_data = np.clip(data, new_min, new_max)
    return clipped_data / (new_max - new_min)

def get_rgb(r, g, b):
    rgb = np.dstack((normalize(r), normalize(g), normalize((b))))
    # rgb = np.dstack((stretch_data(r, 2), stretch_data(g, 2), stretch_data(b, 2)))
    return rgb

def save_to_png(rgbArray, out_folder, out_var_name):
    fig, ax = plt.subplots()
    rgbArray = np.rot90(rgbArray, k=2)
    ax.imshow(rgbArray, aspect='equal')
    plt.axis('off')

    # 去除图像周围白边
    height, width, channels = rgbArray.shape
    # 如果dpi=300，那么图像大小 height*width
    fig.set_size_inches(width/100.0/3.0, height/100.0/3.0)
    plt.gca().xaxis.set_major_locator(plt.NullLocator())
    plt.gca().yaxis.set_major_locator(plt.NullLocator())
    plt.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0)
    plt.margins(0, 0)

    plt.savefig(os.path.join(out_folder, out_var_name + ".png"), dpi=300)
    plt.close('all')


def save_to_tif(geo_extent, out_var, out_folder, out_var_name, non_value):
    maxlon = geo_extent['maxLon']
    minlon = geo_extent['minLon']
    maxlat = geo_extent['maxLat']
    minlat = geo_extent['minLat']

    # 分辨率计算
    N_Lat = out_var.shape[0]
    N_Lon = out_var.shape[1]
    Lon_Res = (maxlon - minlon) / (float(N_Lon) - 1)
    Lat_Res = (maxlat - minlat) / (float(N_Lat) - 1)

    #去除单个像元
    out_va = np.ma.masked_array(out_var, mask=(out_var == non_value) | (out_var == 100))
    binary_image = np.ma.where(out_va == 255, 1, 0)
    cleaned = morphology.remove_small_objects(binary_image.astype(bool), min_size=2, connectivity=1).astype(int)
    # black out pixels
    mask_x = np.ma.where(cleaned == 0)
    out_va[mask_x] = 0

    # 创建.tif文件
    driver = gdal.GetDriverByName('GTiff')

    out_tif_name = os.path.join(out_folder, out_var_name + '.tif')
    out_tif = driver.Create(out_tif_name, N_Lon, N_Lat, 1, gdal.GDT_Byte, options=["TILED=YES", "COMPRESS=LZW"])

    # 设置影像的显示范围
    # -Lat_Res一定要是-的
    geotransform = (minlon, Lon_Res, 0, maxlat, 0, -Lat_Res)
    out_tif.SetGeoTransform(geotransform)

    out_tif.SetProjection('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')  # 给新建图层赋予投影信息

    """字节仅接受8位值，因此不支持NaN """
    # 数据写出
    out_tif.GetRasterBand(1).WriteArray(out_va)  # 将数据写入内存，此时没有写入硬盘
    out_tif.GetRasterBand(1).SetNoDataValue(non_value)

    out_tif.FlushCache()  # 将数据写入硬盘
    out_tif = None  # 注意必须关闭tif文件
    print("L2 Process Finished -- %s " % out_var_name)

def save_to_csv(tmp, out_folder, out_var_name):
    with open(os.path.join(out_folder, "geo_" + out_var_name[1:] + ".csv"), 'w') as f:
        [f.write('{0},{1}\n'.format(key, value)) for key, value in tmp.items()]


def process_redtide(in_path, *args):
    dirname, filename = os.path.split(in_path)
    out_var_name = str(filename[:14])

    out_folder = args[0]
    out_rgbfolder = args[1]
    out_geofolder = args[2]
    non_value = 99
    with nc.Dataset(in_path, 'r') as nc_file:

        navigation_variables = nc_file.groups["navigation_data"].variables
        lon = navigation_variables['longitude'][:]
        lat = navigation_variables['latitude'][:]

        if lon.mask.all() | isinstance(np.nanmin(lon), np.ma.core.MaskedConstant):
            print('L2 Process Failed   -- %s ************************ Masks all around'% out_var_name)
            return out_var_name
        else:
            classdata, rgbarray = main_algorithm(nc_file)
            # minlat = float(nc_file.getncattr('geospatial_lat_min'))
            # minlon = float(nc_file.getncattr('geospatial_lon_min'))
            # maxlat = float(nc_file.getncattr('geospatial_lat_max'))
            # maxlon = float(nc_file.getncattr('geospatial_lon_max'))
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

            right_lon = lon[~split_flag]
            right_lat = lat[~split_flag]
            right_var = classdata[~split_flag]

            if (len(left_lon) == 0) | (len(right_lon) == 0):
                var_re, grid = swath_resampling(classdata, lon, lat, x, y, 3000, non_value)  # 1 km grid
                extent_json = {
                    'maxLon': maxlon,
                    'minLon': minlon,
                    'maxLat': maxlat,
                    'minLat': minlat
                }
                save_to_tif(extent_json, var_re, out_folder, out_var_name, non_value)
                extent_json['tifname'] = out_var_name
                save_to_csv(extent_json, out_geofolder, out_var_name)
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

                if (left_minlon - right_maxlon) > 300:
                    left_var_re, left_grid = swath_resampling(left_var, left_lon.filled(), left_lat.filled(), left_x,
                                                              left_y, 3000, non_value)  # 1 km grid
                    left_extent_json = {
                        'maxLon': left_maxlon,
                        'minLon': left_minlon,
                        'maxLat': left_maxlat,
                        'minLat': left_minlat
                    }

                    right_var_re, right_grid = swath_resampling(right_var, right_lon.filled(), right_lat.filled(), right_x,
                                                                right_y, 3000, non_value)  # 1 km grid
                    right_extent_json = {
                        'maxLon': right_maxlon,
                        'minLon': right_minlon,
                        'maxLat': right_maxlat,
                        'minLat': right_minlat
                    }

                    save_to_tif(left_extent_json, left_var_re, out_folder, out_var_name + '_1', non_value)
                    save_to_tif(right_extent_json, right_var_re, out_folder, out_var_name + '_2', non_value)
                    left_extent_json['tifname'] = out_var_name + '_1'
                    right_extent_json['tifname'] = out_var_name + '_2'

                    save_to_csv(left_extent_json, out_geofolder, out_var_name + '_1')
                    save_to_csv(right_extent_json, out_geofolder, out_var_name + '_2')
                else:
                    var_re, grid = swath_resampling(classdata, lon, lat, x, y, 3000, non_value)  # 1 km grid
                    extent_json = {
                        'maxLon': maxlon,
                        'minLon': minlon,
                        'maxLat': maxlat,
                        'minLat': minlat
                    }
                    extent_json['tifname'] = out_var_name

                    save_to_tif(extent_json, var_re, out_folder, out_var_name, non_value)
                    save_to_csv(extent_json, out_geofolder, out_var_name)
            save_to_png(rgbarray, out_rgbfolder, out_var_name)


def batch_extract_gridded_redtide(in_folder, out_folder, out_rgbfolder, out_logfolder, out_geofolder):
    invalid_list = []
    def results(val):
        if (val != None) & isinstance(val, str):
           invalid_list.append(val)

    kssj = int(time.time() * 1000)  # 定义开始时间到毫秒，因此*1000
    print('开始时间数据（毫秒）：', kssj)
    pool = mp.Pool(processes=50)  # 定义一个进程池Pool，并定义CPU核数量
    files = [i for i in os.listdir(in_folder) if (("L2_LAC_OC" in i) & (i.endswith(".nc")))]
    args = (out_folder, out_rgbfolder, out_geofolder, )  # 类型异常：元组只能与元组进行连接, 因为不加逗号的话，(out_filename)并不是包含一个元素的元组
    for i in range(len(files)):
        pool.apply_async(process_redtide, args=(os.path.join(in_folder, files[i]), ) + args, callback=results)

    pool.close()
    pool.join()
    jssj = int(time.time() * 1000)
    ms = jssj - kssj
    print('结束时间数据（毫秒）：', jssj)
    print('时间差额（毫秒）:', ms)
    # 保留两位小数，但若ms太小，h就会显示为0。
    s = round(ms / 1000, 2)
    m = round(s / 60, 2)
    h = round(m / 60, 2)
    print('{0}换算后等于{1}秒，等于{2}分钟，等于{3}小时'.format(ms, s, m, h))
    df = pd.DataFrame(invalid_list, columns=['invalid'])
    df.to_csv(os.path.join(out_logfolder, files[0][1:5] + '_invalid.csv'), index=None)


def main():
    parser = argparse.ArgumentParser(description='extract gridded redtide')
    parser.add_argument('--in_folder', required=True, default=None)
    parser.add_argument('--out_folder', required=True, default=None)
    parser.add_argument('--out_rgbfolder', required=True, default=None)
    parser.add_argument('--out_logfolder', required=True, default=None)
    parser.add_argument('--out_geofolder', required=True, default=None)
    args = parser.parse_args()
    batch_extract_gridded_redtide(**vars(args))


if __name__ == '__main__':
    # batch_extract_gridded_redtide(
    #     r'E:\l2_process\test1',
    #     r'E:\l2_process\test_result',
    #     r'E:\l2_process\test_rgb_result',
    #     r'E:\l2_process\test_log_result',
    #     r'E:\l2_process\test_geo_result',
    # )
    main()

