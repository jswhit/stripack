from netCDF4 import Dataset
import numpy as np
import time
from stripack import trmesh
from spherical_geometry import vector, great_circle_arc

def _regrid(lons, lats, data, olons, olats, shuffle=False, tri=None):
    """
    performs triangulation and linear interpolation.

    lons,lats: 1d arrays of mesh lat/lon values in radians.
    data: 1d array of mesh data values.
    olons, olats: 1d arrays describing 2d output mesh (degrees)
    shuffle : if True, randomly shuffle mesh points
    tri : if not None, use existing triangulation.

    returns 2d data array on regular lat/lon mesh"""
    if tri is None:
        if shuffle:
            # randomly shuffle points (may speed up triangulation)
            ix = np.arange(len(lons))
            np.random.shuffle(ix)
            lons = lons[ix]; lats = lats[ix]; data = data[ix]
        t1 = time.clock()
        print 'triangulation of', len(lons),' points'
        tri = trmesh(lons, lats)
        if shuffle:
            tri._shuffle = shuffle
            tri._ix = ix
        print 'triangulation took',time.clock()-t1,' secs'
    olons = np.radians(olons);  olats = np.radians(olats)
    olons, olats = np.meshgrid(olons, olats)
    t1 = time.clock()
    latlon_data = tri.interp_linear(olons,olats,data)
    print 'interpolation took',time.clock()-t1,' secs'
    print 'min/max field:',latlon_data.min(), latlon_data.max()
    return latlon_data, tri

def mpas_regrid(filename, olons, olats, varname, slice_out, tri=None):
    """
    Regrid MPAS icos data to reg lat/lon mesh, using
    triangulation + linear interpolation.

    filename : data filename.
    olons : 1d array of lons (in degrees) for output mesh
    olats : 1d array of lats (in degrees) for output mesh
    varname : var name to read from data files
    slice_out : python slice object for data slice
    tri : existing triangulation instance.  If None, it is computed.

    returns 2d data array on regular lat/lon mesh, triangulation instance."""
    print 'reading ',filename
    nc = Dataset(filename)
    if tri is None:
        lons = nc.variables['lonCell'][:].squeeze()
        lats = nc.variables['latCell'][:].squeeze()
    else:
        lons = tri.lats
        lats = tri.lons
    data = nc.variables[varname][slice_out].squeeze()
    nc.close()
    if getattr(tri, '_shuffle',False): data = data[tri._ix]
    return _regrid(lons, lats, data, olons, olats, shuffle=False, tri=tri)

def fv3_regrid(datafiles,gridspecfiles,olons,olats,varname,slice_out,tri=None):
    """
    Regrid FV3 cubed sphere data to reg lat/lon mesh using
    triangulation + linear interpolation.

    datafiles : list of data filenames for each tile
    gridspecfiles : list of grid_spec filenames for each tile
    olons : 1d array of lons (in degrees) for output mesh
    olats : 1d array of lats (in degrees) for output mesh
    varname : var name to read from data files
    slice_out : python slice object for data slice
    tri : existing triangulation instance.  If None, it is computed.

    returns 2d data array on regular lat/lon mesh, triangulation instance."""
    lons = []; lats = []; data = []; A = None
    for datafile, gridspecfile in zip(datafiles,gridspecfiles):
        print 'reading ',datafile
        nc = Dataset(datafile)
        arr = nc.variables[varname][slice_out]
        data.append(arr)
        nc.close()
        if tri is None:
            print 'reading ',gridspecfile
            nc = Dataset(gridspecfile)
            lons_tmp = nc.variables['grid_lon'][:]
            lats_tmp = nc.variables['grid_lat'][:]
            latsA = lats_tmp[:-1,:-1].ravel(); lonsA = lons_tmp[:-1,:-1].ravel() # lower left
            latsC = lats_tmp[1:,:-1].ravel(); lonsC = lons_tmp[1:,:-1].ravel() # upper left
            latsB = lats_tmp[1:,1:].ravel(); lonsB = lons_tmp[1:,1:].ravel() # upper right
            latsD = lats_tmp[:-1,1:].ravel(); lonsD = lons_tmp[:-1,1:].ravel() # lower right
            if A is None:
                nxp1,nyp1 = lons_tmp.shape; nx = nxp1-1; ny = nyp1-1
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
            lonsmid, latsmid =\
            vector.vector_to_radec(T[:,0],T[:,1],T[:,2],degrees=True)
            lonsmid = lonsmid.reshape((ny,nx))
            latsmid = latsmid.reshape((ny,nx))
            lons.append(lonsmid); lats.append(latsmid)
            nc.close()
    if tri is None:
        lons = np.radians(np.array(lons,dtype=np.float64)).ravel()
        lats = np.radians(np.array(lats,dtype=np.float64)).ravel()
    else:
        lons = tri.lons
        lats = tri.lats
    data = (np.array(data,dtype=np.float64)).ravel()
    if getattr(tri, '_shuffle',False): data = data[tri._ix]
    return _regrid(lons, lats, data, olons, olats, shuffle=True, tri=tri)
