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

# then load trmesh object from pickle (pre-generated using generate_mesh.py)

picklefile = 'C%s_grid.pickle' % res
tri = cPickle.load(open(picklefile,'rb'))

# convenience function to generate output grid and do interpolation.
def interp_latlon(nlons,data,tri,order=1):
    '''
    nlons:  number of longitudes on output grid
    data: 1d array of cubed sphere data
    tri: stripack triangulation object for cubed sphere grid
    order: order of interpolation (0 for nearest neighbor, 
           1 for linear, 3 for cubic - default is 1).

    returns lons2d,lats2d,latlon_data where
    lons2d: 2d array of longitudes on output grid (in degrees)
    lats2d: 2d array of longitudes on output grid (in degrees)
    latlon_data:  2d array of interpolated data on output grid.
    '''
    # generate regular 2d lat/lon grid (including poles,
    # but not wrap-around longitude).
    nlats = nlons/2 
    olons = (360./nlons)*np.arange(nlons)
    olats = -90 + 0.5*(360./nlons) + (360./nlons)*np.arange(nlats)
    olonsd, olatsd = np.meshgrid(olons, olats) # degrees
    # interpolate to the reg lat/lon grid
    t1 = time.time()
    if tri._shuffle:
        data = data[tri._ix] # points were randomly shuffled to make triangulation faster
    latlon_data = tri.interp(np.radians(olonsd),np.radians(olatsd),data,order=order) # expects radians
    print('time to interpolation =',time.time()-t1)
    print(latlon_data.shape, latlon_data.min(), latlon_data.max())
    return olonsd,olatsd,latlon_data

# do interpolation.
lons2d,lats2d,latlon_data = interp_latlon(1440,data,tri) # 1/4 deg output grid

# make a plot.
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
fig = plt.figure()
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
ax.contourf(lons2d,lats2d,latlon_data,15)
plt.savefig('test_interp.png')
