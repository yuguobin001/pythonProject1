import rasterio
from rasterio.plot import show
import geopandas as gpd
import matplotlib.pyplot as plt

# 打开tif格式栅格影像
input_raster_path = "E:/coccolith/result_dai/correct/fre/average/average.tif"
src = rasterio.open(input_raster_path)

# 读取栅格数据
raster_data = src.read(1)

# 设置色带范围（这里使用默认的0-255，可以根据影像实际情况调整）
vmin, vmax = 0, 0.3

# 打开Shapefile文件
shapefile_path = "E:/coccolith/World_Continents/World_Continents.shp"
gdf = gpd.read_file(shapefile_path)

# 使用matplotlib绘制影像
plt.figure(figsize=(10, 10))
ax = show(raster_data, cmap='viridis', vmin=vmin, vmax=vmax)

# 在同一图上绘制Shapefile文件
gdf.plot(ax=ax, color='red', edgecolor='black')

# 添加色带图例
cbar = plt.colorbar(ax.get_images()[0], ax=ax)
cbar.set_label('Pixel Value')

# 显示影像
plt.title('Raster Image with Shapefile Overlay')
plt.show()
