from osgeo import gdal
import numpy as np
import matplotlib.pyplot as plt

# 打开栅格影像文件
file_path = 'E:/coccolith/result_dai/correct/fre/average/slope_map10km.tif'
dataset = gdal.Open(file_path)

if dataset is None:
    print("无法打开栅格影像文件")
    exit()

# 读取栅格数据
band = dataset.GetRasterBand(1)
raster_data = band.ReadAsArray()

# 关闭数据集
dataset = None

# 计算直方图
histogram = np.histogram(raster_data, bins=17, range=(-0.05, 0.05))

# 提取直方图数据
hist_values, bin_edges = histogram

# 绘制直方图
plt.figure(figsize=(10, 5))
plt.bar(bin_edges[:-1], hist_values, width=1)
plt.title("栅格影像直方图")
plt.xlabel("像元值")
plt.ylabel("像元数量")
plt.show()
