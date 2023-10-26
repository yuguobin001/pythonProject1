import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import matplotlib as mpl
import os
import matplotlib.path as mpath
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.ticker as tk
from matplotlib import gridspec
import xarray as xr
import cartopy
# 导入Cartopy专门提供的经纬度的Formatter
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import rioxarray as rxr
land = cfeature.LAND.with_scale('110m')
# 读取并打开 TIFF 栅格图像
input_tif_path = "E:/coccolith/result_dai/correct/fre/average/slope_map10km.tif"

def export_Fig2a(in_slopetif, out_folder):
    font = {
        'family': 'Arial',
        'color': 'black',
        'weight': 'normal',
        'size': 6
    }
    title_font = {
        'family': 'Arial',
        'color': 'black',
        'weight': 'bold',
        'size': 8
    }

    fig = plt.figure(num=1, figsize=(7, 2.5))  # figsize: w,h letter size
    titles = ['N', 'S']
    gs = gridspec.GridSpec(1, 2)  # the middle row is to place the colorbar axes
    gs.update(bottom=0.1, top=0.98, left=0, right=1,wspace=0.2)

    #  子图1，北极地区投影
    ax1 = plt.subplot(gs[0, 0], projection=ccrs.NorthPolarStereo(central_longitude=0))
    ax1.text(-185, 80, titles[0], fontdict=title_font)
    left_extent = [-180, 180, 45, 90]  # 子图1的范围
    ax1.set_extent(left_extent, crs=ccrs.PlateCarree())  # 设置地图范围
    ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='grey', linestyle='--', zorder=3)
    ax1.add_feature(land, facecolor='0.8', linewidth=1, zorder=2)
    # 绘制极坐标轴
    theta = np.linspace(0, 2 * np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)
    ax1.set_boundary(circle, transform=ax1.transAxes)

    # 子图2，南极地区投影
    ax2 = plt.subplot(gs[0, 1], projection=ccrs.SouthPolarStereo(central_longitude=0))  # 使用南极投影
    ax2.text(-185, -80, titles[1], fontdict=title_font)  # 添加标题
    right_extent = [-180, 180, -45, -90]  # 子图2的范围
    ax2.set_extent(right_extent, crs=ccrs.PlateCarree())  # 设置地图范围
    ax2.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='grey', linestyle='--', zorder=3)
    ax2.add_feature(land, facecolor='0.8', linewidth=1, zorder=2)
    # 绘制极坐标轴
    ax2.set_boundary(circle, transform=ax2.transAxes)

    export_slope = rxr.open_rasterio(in_slopetif)
    lats = np.arange(900, -900, -1)
    lons = np.arange(-1800, 1800, 1)
    lons, lats = np.meshgrid(lons, lats)

    rgb = (
        [8, 87, 156], [33, 113, 181], [66, 146, 199], [90, 160, 205], [120, 191, 214], [170, 220, 230], [219, 245, 255],
        [255, 255, 255], [255, 224, 224],
        [252, 187, 170], [252, 146, 114], [251, 106, 74], [240, 60, 43], [204, 24, 30], [166, 15, 20])
    rgb = np.array(rgb) / 255.0
    cmap = mcolors.ListedColormap(rgb, name='my_color')
    # cmap_color = icmap

    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmap = LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N)
    bounds_tick1 = np.array([-1, -0.5, 0, 0.5, 1])
    norm1 = mpl.colors.BoundaryNorm(bounds_tick1, cmap.N)
    slope_vals = export_slope.values[0, :, :]

    # 在子图1中绘制北极地区的栅格数据
    cs1 = ax1.pcolormesh(lons, lats, slope_vals*100,
                        transform=ccrs.PlateCarree(), cmap=cmap,
                        vmax=1, vmin=-1,
                        zorder=2)
    # 在子图2中绘制南极地区的栅格数据
    cs2 = ax2.pcolormesh(lons, lats, slope_vals * 100,
                         transform=ccrs.PlateCarree(), cmap=cmap,
                         vmax=1, vmin=-1,
                         zorder=2)

    # 在子图1中绘制图例
    cbar_ax = ax1.inset_axes([0.62, 0.07, 0.4, 0.03])
    cb = fig.colorbar(mappable=cs1, cax=cbar_ax, orientation='horizontal',
                      extend='both',
                      pad=0.01, shrink=0.5, fraction=0.6)

    tick_label_list = [str(i) for i in bounds_tick1]
    cb.set_ticks(bounds_tick1)
    cb.ax.set_xticklabels(['-1.0', '-0.5', '0', '0.5', '1.0'], fontdict=font)
    cb.ax.tick_params(width=0.7, length=3, which='major')
    cb.ax.tick_params(width=0.5, length=2, which='minor')
    cb.set_label(r'Slope (10$^2$ yr$^{-1}$)', labelpad=-29, x=0.5)
    cb.ax.xaxis.set_minor_locator(tk.MultipleLocator(0.25))

    plt.rcParams['pdf.fonttype'] = 42
    fig.savefig(os.path.join(out_folder, 'export_figNS_test.pdf'), format='pdf', dpi=300, bbox_inches='tight',
                pad_inches=0.1)

export_Fig2a(r"E:/coccolith/result_dai/correct/fre/average/slope_map10km.tif",
             r"E:/coccolith/result_dai/correct/fre/average/map")
