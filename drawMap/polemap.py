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
# 读取并打开 TIFF 栅格图像
input_tif_path = "E:/coccolith/result_dai/correct/fre/average/slope_map100km.tif"  # 请替换为你的输入文件路径
ds = gdal.Open(input_tif_path)
raster_data = ds.ReadAsArray()

# 定义拉伸范围
vmin, vmax = -0.01, 0.01

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
cmap1 = cm.RdBu_r
norm1 = mcolors.TwoSlopeNorm(vmin=-0.01, vcenter=0, vmax=0.01)
bins = np.array([-0.01, -0.009, -0.008, -0.007,-0.006,-0.005,-0.004,-0.003,-0.002,-0.001,0, 0.001, 0.002, 0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01])
#bins = np.array([-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1])
nbin = len(bins)-1
n_negative = np.count_nonzero(bins < 0)
n_positive = np.count_nonzero(bins > 0)
colors = np.vstack((
    cmap1(np.linspace(0, 0.5, n_negative))[:-1],
    cmap1(np.linspace(0.5, 1, n_positive))
))  # 根据bins的区间数新建colormap.
cmap2 = mcolors.ListedColormap(colors)
norm2 = mcolors.BoundaryNorm(bins, nbin)
levels = np.linspace(bins.min(), bins.max(), 16)

#coolwarm_custom = ListedColormap(['#08579c',(123,178,212),(187,220,236),(235,244,226),(254,238,168),(251,191,113),(243,124,72),(216,54,43),(108,13,18)])
#coolwarm_custom = ListedColormap([(69,117,177),(123,178,212),(187,220,236),(235,244,226),(254,238,168),(251,191,113),(243,124,72),(216,54,43),(108,13,18)])
plt.figure(figsize=(40,20))#设置画布大小
fig, axes = plt.subplots(nrows=1, ncols=2)
# 去掉子图的框线
for ax in axes:
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
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
im = axes[0].imshow(
    raster_data,
    origin="upper",
    cmap=cmap1,  # 可根据需要选择色带
    norm=norm1,
    transform=ccrs.PlateCarree(),
    vmin=vmin,
    vmax=vmax,
)
# 添加色标
cbar = plt.colorbar(im, ax=axes[0], ticks=levels[1::2])
# cbar.set_label('Pixel Values')

# 在地图上绘制南极数据
imS = axes[1].imshow(
    raster_data,
    origin="upper",
    cmap=cmap2,  # 可根据需要选择色带
    norm=norm2,
    transform=ccrs.PlateCarree(),
    vmin=vmin,
    vmax=vmax,
)
# 添加色标
cbar = plt.colorbar(imS, ax=axes[1], ticks=bins)
cbar.set_label('Pixel Values')
# 显示图像
plt.show()
