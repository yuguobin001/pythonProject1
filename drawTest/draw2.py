# from osgeo import gdal
# import matplotlib.pyplot as plt
#
# dataset = gdal.Open("E:/coccolith/trueFre/ture001/cooSum2003_sum2003.tif", gdal.GA_ReadOnly)
# band = dataset.GetRasterBand(1) # 波段序号从1开始，而不是0
# plt.figure(figsize=(10, 10))
# plt.imshow(band.ReadAsArray())
# plt.show()
# !/usr/bin/python

str = "this is string example....wow!!!";

suffix = "wow";
print(str.endswith(suffix, -6, -3))