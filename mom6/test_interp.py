from netCDF4 import Dataset
import numpy as np
from stripack import trmesh
import sys, time
import _pickle as cPickle

if len(sys.argv) > 1:
    res = sys.argv[1]
else:
    raise SystemExit("input mom6 resolution (025 for 0.25 deg)")

# path to mom6 test data
mom6testdata = './ocn_2016_01_01_01.nc'

# perform interpolation of orography using saved triangulation.

# first read data from each of the 6 tile files, append into 1d array
nc = Dataset(mom6testdata)
data = nc['temp'][0,:,:]
nc.close()
data = (np.asarray(data,dtype=np.float64)).ravel()

# then load trmesh object from pickle (pre-generated using generate_mom6mesh.py)

picklefile = 'mom6_mx%s_grid.pickle' % res
tri = cPickle.load(open(picklefile,'rb'))

# do interpolation (nearest neighbor).
# 1/4 deg output grid.
lons2d,lats2d,latlon_data_nn = tri.interp_latlon(1440,data,order=0) 
latlon_data_nn = np.ma.masked_less(latlon_data_nn,-1.e10)
# do interpolation (bilinear)
lons2d,lats2d,latlon_data_lin = tri.interp_latlon(1440,data,order=1) 
latlon_data_lin = np.ma.masked_less(latlon_data_lin,-1.e10)
# combine (use un-masked bilinear points, with nearest neighbor mask)
latlon_data = latlon_data_nn
ix = latlon_data_lin.mask == False
latlon_data[ix] = latlon_data_lin[ix]

# make a plot.
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
fig = plt.figure()
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
ax.contourf(lons2d,lats2d,latlon_data,15)
plt.savefig('test_interp.png')
