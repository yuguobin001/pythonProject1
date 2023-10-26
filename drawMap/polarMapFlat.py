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
from osgeo import gdal


def export_Fig2a(in_slopetif, out_folder):
    font = {
        'family': 'Arial',
        'color': 'black',
        'weight': 'normal',
        'size': 7
    }
    title_font = {
        'family': 'Arial',
        'color': 'black',
        'weight': 'bold',
        'size': 8
    }
    mpl.rcParams["font.family"] = 'Arial'  # 默认字体类型
    mpl.rcParams["mathtext.fontset"] = 'cm'  # 数学文字字体
    mpl.rcParams["font.size"] = 12  # 字体大小
    mpl.rcParams["axes.linewidth"] = 1

    fig = plt.figure(num=1, figsize=(8, 4))  # figsize: w,h letter size
    titles = ['a', 'b']
    gs = gridspec.GridSpec(1, 1)  # the middle row is to place the colorbar axes
    gs.update(bottom=0.1, top=0.98, left=0, right=1,wspace=0.2)

    extent = [-180, 180, -90, 90]
    ax1 = plt.subplot(gs[0, 0], projection=ccrs.PlateCarree(central_longitude=0))
    ax1.set_extent(extent, crs=ccrs.PlateCarree())
    ax1.add_feature(cartopy.feature.NaturalEarthFeature('physical', 'land', '50m', facecolor='#D3D3D3', zorder=2))
    ax1.axis('off')


    ds = gdal.Open(in_slopetif)
    band = ds.GetRasterBand(1)
    data = band.ReadAsArray()
    export_slope = data * 100


    rgb = (
        [8, 87, 156], [33, 113, 181], [66, 146, 199], [90, 160, 205], [120, 191, 214], [170, 220, 230], [219, 245, 255],
        [255, 255, 255], [255, 224, 224],[252, 187, 170], [252, 146, 114], [251, 106, 74], [240, 60, 43], [204, 24, 30], [166, 15, 20])
    rgb = np.array(rgb) / 255.0
    cmap = mcolors.ListedColormap(rgb, name='my_color')
    # cmap_color = icmap

    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmap = LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N)
    bounds_tick1 = np.array([-1, -0.5, 0, 0.5, 1])
    #norm1 = mpl.colors.BoundaryNorm(bounds_tick1, cmap.N)

    cs = ax1.imshow(export_slope,
                        transform=ccrs.PlateCarree(), cmap=cmap,
                        extent= extent,
                        vmax=1, vmin=-1,
                        zorder=1)

    cbar_ax = ax1.inset_axes([0.62, 0.07, 0.2, 0.03])
    cb = fig.colorbar(mappable=cs, cax=cbar_ax, orientation='horizontal',
                      extend='both',
                      pad=0.01, shrink=0.5, fraction=0.6)

    tick_label_list = [str(i) for i in bounds_tick1]
    cb.set_ticks(bounds_tick1)
    cb.ax.set_xticklabels(['-1.0', '-0.5', '0', '0.5', '1.0'], fontdict=font)
    cb.ax.tick_params(width=0.7, length=3, which='major')
    cb.ax.tick_params(width=0.5, length=2, which='minor')
    cb.set_label(r'Slope (10$^2$ yr$^{-1}$)', labelpad=-32, x=0.5)
    cb.ax.xaxis.set_minor_locator(tk.MultipleLocator(0.25))
    background_color = '#D3D3D3'  # 设置背景颜色
    cb.ax.add_patch(plt.Rectangle((-1.2, -1.2), 3, 3, color=background_color, zorder=0))

    plt.rcParams['pdf.fonttype'] = 42
    fig.savefig(os.path.join(out_folder, 'export_fig10m.pdf'), format='pdf', dpi=300, bbox_inches='tight',
                pad_inches=0.1)
    plt.show()

export_Fig2a(r"E:/coccolith/result_dai/correct/fre/average/slope_map10km.tif",
             r"E:/coccolith/result_dai/correct/fre/average/map")
