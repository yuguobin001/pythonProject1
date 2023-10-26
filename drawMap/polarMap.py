import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import matplotlib as mpl
from osgeo import gdal
import matplotlib.path as mpath
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.ticker as tk
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
# 读取并打开 TIFF 栅格图像
input_tif_path = "E:/coccolith/result_dai/correct/fre/average/slope_map100km.tif"  # 请替换为你的输入文件路径
ds = gdal.Open(input_tif_path)
raster_data = ds.ReadAsArray()*100

# 定义拉伸范围
vmin,vmax = -1,1
# vmin = np.min(raster_data)
# vmax = np.max(raster_data)
# 创建图形和坐标系
mpl.rcParams["font.family"] = 'Arial'  #默认字体类型
mpl.rcParams["mathtext.fontset"] = 'cm' #数学文字字体
mpl.rcParams["font.size"] = 12   #字体大小
mpl.rcParams["axes.linewidth"] = 1
# 设置左边图形的范围和位置
left_extent = [-180, 180, 45, 90]
left_position = [0.05, 0.1, 0.4, 0.8]
# 设置右边图形的范围和位置
right_extent = [-180, 180, -45, -90]
right_position = [0.55, 0.1, 0.4, 0.8]
land = cfeature.LAND.with_scale('110m')

# 创建自定义 coolwarm 色带
rgb = (
        [8, 87, 156], [33, 113, 181], [66, 146, 199],
        [90, 160, 205], [120, 191, 214],
        [170, 220, 230], [219, 245, 255],
        [255, 255, 255], [255, 224, 224],
        [252, 187, 170], [252, 146, 114],
        [251, 106, 74], [240, 60, 43], [204, 24, 30], [166, 15, 20])
rgb = np.array(rgb) / 255.0
cmap = mcolors.ListedColormap(rgb, name='my_color')
# cmap_color = icmap
cmaplist = [cmap(i) for i in range(cmap.N)]
cmap = LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N)
bounds_tick1 = np.array([-1, -0.5, 0, 0.5, 1])
norm1 = mpl.colors.BoundaryNorm(bounds_tick1, cmap.N)

plt.figure(figsize=(7,2.5))#设置画布大小
fig, axes = plt.subplots(nrows=1, ncols=2)
# 去掉子图的框线
for ax in axes:
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)

#绘制地图要素
axes[0] = fig.add_axes(left_position,projection = ccrs.NorthPolarStereo(central_longitude=0))#绘制地图位置
axes[0].gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='grey',linestyle='--',zorder=3)
axes[0].set_extent(left_extent, ccrs.PlateCarree())
axes[0].add_feature(land,facecolor='0.8', linewidth=1,zorder=2)
#axes[0].add_feature(cfeature.COASTLINE.with_scale('110m'))

# 绘制极坐标轴
theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)
axes[0].set_boundary(circle, transform=axes[0].transAxes)


axes[1] = fig.add_axes(right_position,projection = ccrs.SouthPolarStereo(central_longitude=0))#绘制地图位置
axes[1].gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='grey',linestyle='--',zorder=3)
axes[1].set_extent(right_extent, ccrs.PlateCarree())
axes[1].add_feature(land,facecolor='0.8', linewidth=1,zorder=2)
#axes[0].add_feature(cfeature.COASTLINE.with_scale('110m'))
# 绘制极坐标轴
theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)
axes[1].set_boundary(circle, transform=axes[1].transAxes)


# 在地图上绘制北极数据
cs = axes[0].imshow(
    raster_data,
    origin="upper",
    cmap=cmap,  # 可根据需要选择色带
    transform=ccrs.PlateCarree(),
    vmin=vmin,
    vmax=vmax,
)
# 添加色标

# 在地图上绘制南极数据
cs2 = axes[1].imshow(
    raster_data,
    origin="upper",
    cmap=cmap,  # 可根据需要选择色带
    transform=ccrs.PlateCarree(),
    vmin=vmin,
    vmax=vmax,
)

# 显示图例
cbar_ax = axes[1].inset_axes([0.62, 0.07, 0.4, 0.03])
# cbar.set_label('Pixel Values')
cb = fig.colorbar(mappable=cs2, cax=cbar_ax, orientation='horizontal',
                  extend='both',
                  pad=0.01, shrink=0.5, fraction=0.6)

cb.set_ticks(bounds_tick1)
cb.set_ticklabels('-1','-0.5','0','0.5','1')
cb.ax.tick_params(width=1, length=3, which='major')
cb.ax.tick_params(width=0.5, length=2, which='minor')
cb.set_label('Slope (10$^2$ yr$^{-1}$)', labelpad=-35, x=0.5, fontdict=font)
cb.ax.xaxis.set_minor_locator(tk.MultipleLocator(0.25))
plt.rcParams['pdf.fonttype'] = 42
plt.show()
