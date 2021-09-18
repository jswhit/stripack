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
# read data from each of the 6 tile files, append into 1d array
data = []
for ntile in range(1,7,1):
    gridfile = '%s/C%s/C%s_oro_data.tile%s.nc'% (fixfv3,res,res,ntile)
    nc = Dataset(gridfile)
    data.append(nc['orog_filt'][:,:])
    nc.close()
data = (np.array(data,dtype=np.float64)).ravel()

# load trmesh object from pickle (pre-generated using generate_mesh.py)

picklefile = 'C%s_grid.pickle' % res
tri = cPickle.load(open(picklefile,'rb'))

# generate regular 2d 0.25deg lat/lon grid
nlons = 1440; nlats = nlons/2 
olons = (360./nlons)*np.arange(nlons)
olats = -90 + 0.5*(360./nlons) + (360./nlons)*np.arange(nlats)
olonsd, olatsd = np.meshgrid(olons, olats) # degrees
olons = np.radians(olons);  olats = np.radians(olats) # radians
olons, olats = np.meshgrid(olons, olats)

# interpolate linearly to the reg grid
# interp_cubic and interp_nn (nearest neighbor) also available
t1 = time.time()
if tri._shuffle:
    data = data[tri._ix] # points were randomly shuffled to make triangulation faster
latlon_data = tri.interp_linear(olons,olats,data) # expects radians
print('time to interpolation =',time.time()-t1)
print(latlon_data.shape, latlon_data.min(), latlon_data.max())

# make a plot
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
fig = plt.figure()
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
ax.contourf(olonsd,olatsd,latlon_data,15)
plt.savefig('test_interp.png')
