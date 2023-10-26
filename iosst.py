#Import libraries
import xarray as xr
#import cmocean as cm
import cartopy.crs as ccrs
import pylab as plt
import numpy as np
# Inline plotting
from xarrayMannKendall import Mann_Kendall_test
from utils import area,ccrs_land,add_patches
import datetime as datetime
from dask.distributed import Client
c = Client()