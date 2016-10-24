from netCDF4 import Dataset
import numpy as np
import time
from stripack import trmesh
import matplotlib.pyplot as plt
import pygrib
from mpl_toolkits.basemap import Basemap
from regrid_mod import mpas_regrid

filename='mpas_restart.nc'

nlons = 768; nlats = nlons/2 # T382 gaussian grid
olons = (360./nlons)*np.arange(nlons,dtype=np.float)
olats = pygrib.gaulats(nlats)

varname = 'surface_pressure'
slice_out = np.s_[-1,...]
tri = None
latlon_data, tri = mpas_regrid(filename, olons, olats, varname,\
        slice_out,tri=tri)
latlon_data, tri = mpas_regrid(filename, olons, olats, varname,\
       slice_out,tri=tri)

# make plot on output mesh
m = Basemap(lon_0=180)
m.drawcoastlines()
m.drawmapboundary()
olons, olats = np.meshgrid(olons, olats)
m.contourf(olons,olats,latlon_data,15)
m.colorbar()
plt.show()
