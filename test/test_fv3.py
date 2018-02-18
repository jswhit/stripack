from netCDF4 import Dataset
import numpy as np
import time
from stripack import trmesh
import cPickle

# test fv3 interpolation from native history files to random points.

res = 128 
fixfv3 = '/scratch4/NCEPDEV/global/save/glopara/svn/fv3gfs/fix/fix_fv3_gmted2010'
#fixfv3='/lustre/f1/unswept/Jeffrey.S.Whitaker/fv3_reanl/fv3gfs/global_shared.v15.0.0/fix/fix_fv3_gmted2010'
datapath = '/scratch3/BMC/gsienkf/whitaker/fv3_ics64'
#datapath = '/lustre/f1/unswept/Jeffrey.S.Whitaker/fv3_reanl/ics'
date = '2016100100'
var = 'pressfc'
ntime = 0


# perform triangulation.
lons = []; lats = []
for ntile in range(1,7,1):
    gridfile = '%s/C%s/C%s_grid.tile%s.nc'% (fixfv3,res,res,ntile)
    nc = Dataset(gridfile)
    lonsmid = nc['x'][1::2,1::2]
    latsmid = nc['y'][1::2,1::2]
    lons.append(lonsmid); lats.append(latsmid)
    nc.close()
lons = np.radians(np.array(lons,dtype=np.float64)).ravel()
lats = np.radians(np.array(lats,dtype=np.float64)).ravel()
t1 = time.clock()
print 'triangulation of', len(lons),' points'
tri = trmesh(lons, lats)
print 'triangulation took',time.clock()-t1,' secs'

# pickle it.
picklefile = 'C%s_grid.pickle' % res
cPickle.dump(tri,open(picklefile,'wb'),-1)
# load from pickle
t1 = time.clock()
tri2 = cPickle.load(open(picklefile,'rb'))
print 'load from pickle file took',time.clock()-t1,' secs'

# generate random points on a sphere
# so that every small area on the sphere is expected
# to have the same number of points.
# http://mathworld.wolfram.com/SpherePointPicking.html
npts = 10000
np.random.seed(42) # set random seed for reproducibility
u = np.random.uniform(0.,1.,size=npts)
v = np.random.uniform(0.,1.,size=npts)
olons = np.radians(360.*u)
olats = np.radians((180./np.pi)*np.arccos(2*v-1) - 90.)

# read data from history files.
t1 = time.clock()
data = []
for ntile in range(1,7,1):
    datafile = '%s/C%s_%s/mem001/fv3_history2d.tile%s.nc'% (datapath,res,date,ntile)
    nc = Dataset(datafile)
    data.append(nc[var][ntime,...])
    nc.close()
data = (np.array(data,dtype=np.float64)).ravel()
print 'read %s data from history files took' % var,time.clock()-t1,' secs'
    
# interpolate to random points.
t1 = time.clock()
data_interp = tri2.interp_linear(olons,olats,data)
print 'interpolation %s points took' % npts,time.clock()-t1,' secs'
print 'min/max interpolated field:',data_interp.min(), data_interp.max()
