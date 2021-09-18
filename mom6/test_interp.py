from __future__ import print_function
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

# do interpolation.
lons2d,lats2d,latlon_data = tri.interp_latlon(1440,data,order=0) # 1/4 deg output grid
latlon_data = np.ma.masked_less(latlon_data,-1.e30)

# make a plot.
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
fig = plt.figure()
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
ax.contourf(lons2d,lats2d,latlon_data,15)
plt.savefig('test_interp.png')
