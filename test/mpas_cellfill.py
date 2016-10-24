# make cell-fill plot on native MPAS mesh
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
varname = 'oro_uf'
lon_0 = -120.
lat_0 = 40.
vmin = 0; vmax = 3200.
mapLeft   = lon_0-20.   # longitude
mapRight  = lon_0+20.   # longitude
mapBottom = lat_0-20.   # latitude
mapTop    = lat_0+20.   # latitude

fig = plt.figure()
ax = fig.add_subplot(1,1,1)

filename = 'mpas_restart.nc'
nc = Dataset(filename)
latCell = np.degrees(nc['latCell'][:])
lonCell = np.degrees(nc['lonCell'][:])
latVertex = np.degrees(nc['latVertex'][:])
lonVertex = np.degrees(nc['lonVertex'][:])
nEdgesOnCell = nc['nEdgesOnCell'][:]
verticesOnCell = nc['verticesOnCell'][:]
nCells = len(latCell)
print latCell.shape, verticesOnCell.shape, nEdgesOnCell.shape
maxEdges = verticesOnCell.shape[1]

datin = nc[varname][:]
print datin.min(), datin.max(), datin.shape
nc.close()
lonCell = np.where(lonCell > 180., lonCell-360., lonCell)
lonVertex = np.where(lonVertex > 180., lonVertex-360., lonVertex)
patches = []; data = []
for iCell in range(nCells):
    if  latCell[iCell] >= mapBottom and \
        latCell[iCell] <= mapTop    and \
        lonCell[iCell] >= mapLeft   and \
        lonCell[iCell] <= mapRight:
        xypoly = np.empty((nEdgesOnCell[iCell],2), np.float)
        for i in range(nEdgesOnCell[iCell]):
            xypoly[i,0] = lonVertex[verticesOnCell[iCell,i]-1]
            xypoly[i,1] = latVertex[verticesOnCell[iCell,i]-1]
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
plt.title('MPAS cell-fill orography')
plt.tight_layout()

plt.show()
