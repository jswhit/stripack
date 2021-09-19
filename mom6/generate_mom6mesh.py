from netCDF4 import Dataset
import numpy as np
import sys
from stripack import trmesh
import _pickle as cPickle

# generate mom6 mesh triangulation objects, save to pickle files.

if len(sys.argv) > 1:
    res = sys.argv[1]
else:
    raise SystemExit("input mom6 resolution (025 for 0.25 deg)")

shuffle = True
if len(sys.argv) > 2:
    shuffle = bool(int(sys.argv[1]))
print('generating triangulation for mx%s' % res)
if shuffle: print('grid points randomly shuffled')

# path to MOM6 geometry file 
mom6geom = './ocean_geometry_mx%s.nc' % res

# read in lat/lon of native tile points from oro_data files.
nc = Dataset(mom6geom)
lons = np.radians(nc['geolon'][:,:].astype(np.float64)).ravel()
lats = np.radians(nc['geolat'][:,:].astype(np.float64)).ravel()
nc.close()

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
picklefile = 'mom6_mx%s_grid.pickle' % res
cPickle.dump(tri,open(picklefile,'wb'),-1)
