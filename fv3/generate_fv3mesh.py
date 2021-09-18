from __future__ import print_function
from netCDF4 import Dataset
import numpy as np
import time, sys
from stripack import trmesh
import _pickle as cPickle

# generate FV3 mesh triangulation objects, save to pickle files.

if len(sys.argv) > 1:
    res = int(sys.argv[1])
else:
    raise SystemExit("input cubed sphere resolution")

shuffle = True
if len(sys.argv) > 2:
    shuffle = bool(int(sys.argv[1]))
print('generating triangulation for C%s' % res)
if shuffle: print('grid points randomly shuffled')

# path to FV3 fix files 
fixfv3 = '/scratch1/NCEPDEV/global/glopara/fix_NEW/fix_fv3_gmted2010/'

# read in lat/lon of native tile points from oro_data files.
lons = []; lats = []
for ntile in range(1,7,1):
    gridfile = '%s/C%s/C%s_oro_data.tile%s.nc'% (fixfv3,res,res,ntile)
    nc = Dataset(gridfile)
    lons2 = nc['geolon'][:,:]
    lats2 = nc['geolat'][:,:]
    lons.append(lons2); lats.append(lats2)
    nc.close()
lons = np.radians(np.array(lons,dtype=np.float64)).ravel()
lats = np.radians(np.array(lats,dtype=np.float64)).ravel()

# randomly shuffle points (may speed up triangulation)
ix = np.arange(len(lons))
if shuffle:
    np.random.shuffle(ix)
    lons = lons[ix]; lats = lats[ix]

# generate mesh triangulation object
tri = trmesh(lons, lats)
tri._shuffle = shuffle
tri._ix = ix

# pickle it and save it.
picklefile = 'C%s_grid.pickle' % res
cPickle.dump(tri,open(picklefile,'wb'),-1)
