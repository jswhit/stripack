from netCDF4 import Dataset
import numpy as np
import time
from stripack import trmesh
import cPickle

# test fv3 interpolation from native history files to random points.

res = 96  
fixfv3 = '/scratch4/NCEPDEV/global/save/glopara/svn/fv3gfs/fix_fv3'

gridspecfiles = []
for ntile in range(1,7,1):
    gridspecfiles.append('%s/C%s/C%s_grid_spec.tile%s.nc'% (fixfv3,res,res,ntile))

# test computation of cell centers given cell corners.
#from spherical_geometry import vector, great_circle_arc
#for gridspecfile in gridspecfiles:
#    print 'reading ',gridspecfile
#    nc = Dataset(gridspecfile)
#    lonscorner = nc['grid_lon'][:]
#    latscorner = nc['grid_lat'][:]
#    latsA = latscorner[:-1,:-1].ravel(); lonsA = lonscorner[:-1,:-1].ravel() # lower left
#    latsC = latscorner[1:,:-1].ravel(); lonsC = lonscorner[1:,:-1].ravel() # upper left
#    latsB = latscorner[1:,1:].ravel(); lonsB = lonscorner[1:,1:].ravel() # upper right
#    latsD = latscorner[:-1,1:].ravel(); lonsD = lonscorner[:-1,1:].ravel() # lower right
#    nxp1,nyp1 = lonscorner.shape; nx = nxp1-1; ny = nyp1-1
#    npts = len(latsA)
#    A = np.empty((npts,3),dtype=np.float)
#    B = np.empty((npts,3),dtype=np.float)
#    C = np.empty((npts,3),dtype=np.float)
#    D = np.empty((npts,3),dtype=np.float)
#    A[:,0],A[:,1],A[:,2] = vector.radec_to_vector(lonsA, latsA, degrees=True)
#    B[:,0],B[:,1],B[:,2] = vector.radec_to_vector(lonsB, latsB, degrees=True)
#    C[:,0],C[:,1],C[:,2] = vector.radec_to_vector(lonsC, latsC, degrees=True)
#    D[:,0],D[:,1],D[:,2] = vector.radec_to_vector(lonsD, latsD, degrees=True)
#    # find intersection of great circle arcs from each corner of grid box
#    T = great_circle_arc.intersection(A, B, C, D)
#    lonsmid, latsmid =\
#    vector.vector_to_radec(T[:,0],T[:,1],T[:,2],degrees=True)
#    lonsmid = lonsmid.reshape((ny,nx))
#    latsmid = latsmid.reshape((ny,nx))
#    lonsmid2 = nc['grid_lont'][:]
#    latsmid2 = nc['grid_latt'][:]
#    difflon = lonsmid-lonsmid2; difflat = latsmid-latsmid2
#    print 'lon error min/max',difflon.min(), difflon.max()
#    print 'lat error min/max',difflat.min(), difflat.max()

# perform triangulation.
lons = []; lats = []
for gridspecfile in gridspecfiles:
    nc = Dataset(gridspecfile)
    lonsmid = nc['grid_lont'][:]
    latsmid = nc['grid_latt'][:]
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
datapath = '/scratch3/BMC/gsienkf/whitaker/fv3_ics64'
date = '2016100100'
var = 'pressfc'
ntime = 0
datafiles = []
for ntile in range(1,7,1):
    datafiles.append('%s/C%s_%s/mem001/fv3_history2d.tile%s.nc'% (datapath,res,date,ntile))
t1 = time.clock()
data = []
for datafile in datafiles:
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
