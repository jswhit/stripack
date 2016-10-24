from netCDF4 import Dataset
import numpy as np
import time
from stripack import trmesh
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import pygrib
from regrid_mod import fv3_regrid

datafiles = []
gridspecfiles = []
for ntile in range(1,7,1):
    datafiles.append('fv3_history.tile%s.nc'%ntile)
    gridspecfiles.append('grid_spec.tile%s.nc'%ntile)

nlons = 768; nlats = nlons/2 # T382 gaussian grid
olons = (360./nlons)*np.arange(nlons)
olats = pygrib.gaulats(nlats)

varname = 'psfc'
slice_out = np.s_[-1,:,:]
tri = None
latlon_data, tri = fv3_regrid(datafiles, gridspecfiles, olons, olats, varname,\
        slice_out,tri=tri)
latlon_data, tri = fv3_regrid(datafiles, gridspecfiles, olons, olats, varname,\
       slice_out,tri=tri)

# make plot on output mesh
m = Basemap(lon_0=180)
m.drawcoastlines()
m.drawmapboundary()
olons, olats = np.meshgrid(olons, olats)
m.contourf(olons,olats,latlon_data,15)
m.colorbar()
plt.show()
