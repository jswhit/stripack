from __future__ import print_function
from netCDF4 import Dataset
import numpy as np
from stripack import trmesh
import sys, time
import _pickle as cPickle

if len(sys.argv) > 1:
    res = int(sys.argv[1])
else:
    raise SystemExit("input cubed sphere resolution")

# path to FV3 fix files (to read orography and test interpolation)
fixfv3 = '/scratch1/NCEPDEV/global/glopara/fix_NEW/fix_fv3_gmted2010/'

# perform interpolation of orography using saved triangulation.

# first read data from each of the 6 tile files, append into 1d array
data = []
for ntile in range(1,7,1):
    gridfile = '%s/C%s/C%s_oro_data.tile%s.nc'% (fixfv3,res,res,ntile)
    nc = Dataset(gridfile)
    data.append(nc['orog_filt'][:,:])
    nc.close()
data = (np.array(data,dtype=np.float64)).ravel()

# then load trmesh object from pickle (pre-generated using generate_fv3mesh.py)

picklefile = 'C%s_grid.pickle' % res
tri = cPickle.load(open(picklefile,'rb'))

# do interpolation.
lons2d,lats2d,latlon_data = tri.interp_latlon(1440,data) # 1/4 deg output grid

# make a plot.
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
fig = plt.figure()
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
ax.contourf(lons2d,lats2d,latlon_data,15)
plt.savefig('test_interp.png')
