# make cell-fill plot on native FV3 mesh
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import sys, os
from netCDF4 import Dataset
try:
    from spherical_geometry import vector, great_circle_arc
except ImportError:
    msg='requires spherical_geometry module from https://github.com/spacetelescope/sphere'
    raise ImportError(msg)

# define variable name to plot, bounds of plot region, min/max of color range.
varname = 'orog'
lon_0_fv3 = -120.
lat_0_fv3 = 40.
vmin = 0; vmax = 3200.
mapLeft   = lon_0_fv3-20.   # longitude
mapRight  = lon_0_fv3+20.   # longitude
mapBottom = lat_0_fv3-20.   # latitude
mapTop    = lat_0_fv3+20.   # latitude

fig = plt.figure()
ax = fig.add_subplot(1,1,1)

for ntile in range(1,7):
    tile = 'tile%s' % ntile
    nc = Dataset(os.path.join('grid_spec.%s.nc' % tile))
    lonVertex = nc['grid_lon'][:]
    latVertex = nc['grid_lat'][:]
    latsA = latVertex[:-1,:-1].ravel(); lonsA = lonVertex[:-1,:-1].ravel() # lower left
    latsC = latVertex[1:,:-1].ravel(); lonsC = lonVertex[1:,:-1].ravel() # upper left
    latsB = latVertex[1:,1:].ravel(); lonsB = lonVertex[1:,1:].ravel() # upper right
    latsD = latVertex[:-1,1:].ravel(); lonsD = lonVertex[:-1,1:].ravel() # lower right
    nxp1,nyp1 = lonVertex.shape; nx = nxp1-1; ny = nyp1-1
    npts = len(latsA)
    A = np.empty((npts,3),dtype=np.float)
    B = np.empty((npts,3),dtype=np.float)
    C = np.empty((npts,3),dtype=np.float)
    D = np.empty((npts,3),dtype=np.float)
    A[:,0],A[:,1],A[:,2] = vector.radec_to_vector(lonsA, latsA, degrees=True)
    B[:,0],B[:,1],B[:,2] = vector.radec_to_vector(lonsB, latsB, degrees=True)
    C[:,0],C[:,1],C[:,2] = vector.radec_to_vector(lonsC, latsC, degrees=True)
    D[:,0],D[:,1],D[:,2] = vector.radec_to_vector(lonsD, latsD, degrees=True)
    # find intersection of great circle arcs from each corner of grid box
    T = great_circle_arc.intersection(A, B, C, D)
    lonCell, latCell =\
    vector.vector_to_radec(T[:,0],T[:,1],T[:,2],degrees=True)
    nCells = len(lonCell)
    latVertices = np.empty((nCells,4), np.float)
    lonVertices = np.empty((nCells,4), np.float)
    print tile
    print latsA.min(), latsA.max()
    print lonsA.min(), lonsA.max()
    latVertices[:,0]=latsD
    latVertices[:,1]=latsA
    latVertices[:,2]=latsC
    latVertices[:,3]=latsB
    lonVertices[:,0]=lonsD
    lonVertices[:,1]=lonsA
    lonVertices[:,2]=lonsC
    lonVertices[:,3]=lonsB
    lonCell = np.where(lonCell > 180., lonCell-360., lonCell)
    print lonCell.min(), lonCell.max()
    lonVertices = np.where(lonVertices > 180., lonVertices-360., lonVertices)
    nc.close()

    filename = 'fv3_history.%s.nc' % tile
    nc = Dataset(filename)
    datin = nc[varname][:].ravel()
    print datin.min(), datin.max(), lonCell.shape, datin.shape
    nc.close()

    patches = []; data = []
    for iCell in range(nCells):
        if  latCell[iCell] >= mapBottom and \
            latCell[iCell] <= mapTop    and \
            lonCell[iCell] >= mapLeft   and \
            lonCell[iCell] <= mapRight:
            xypoly = np.empty((4,2), np.float)
            xypoly[:,0] = lonVertices[iCell,:]
            xypoly[:,1] = latVertices[iCell,:]
            polygon = Polygon(xypoly, closed=True)
            patches.append(polygon)
            data.append(datin[iCell])
    print len(patches), ' patches'
    p = PatchCollection(patches, cmap=matplotlib.cm.hot_r, edgecolors='none')
    data = np.array(data)
    p.set_array(data)
    p.set_clim(vmin=vmin,vmax=vmax)
    p.set_clip_on(True)
    ax.add_collection(p)

plt.xlim(mapLeft, mapRight)
plt.ylim(mapBottom, mapTop)
plt.grid(True)
plt.colorbar(p)
plt.title('FV3 cell-fill orography')
plt.tight_layout()

plt.show()
