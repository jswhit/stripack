import  _stripack
import time
import numpy as np

__version__ = "1.2"

class trmesh(object):
    def __init__(self, lons, lats):
        """
given mesh points (lons,lats in radians) define triangulation.
n is size of input mesh (length of 1-d arrays lons and lats).

Algorithm:
 R. J. Renka, "ALGORITHM 772: STRIPACK: Delaunay triangulation
 and Voronoi diagram on the surface of a sphere"
 ACM Trans. Math. Software, Volume 23 Issue 3, Sept. 1997
 pp 416-434.

Instance variables:

lats, lons: lat/lon values of N nodes (in radians)
npts:  number of nodes (N).
x, y, z: cartesian coordinates corresponding to lats, lons.

lst: 6*(N-2) nodal indices, which along with
with lptr and lend, define the triangulation as a set of N
adjacency lists; counterclockwise-ordered sequences of neighboring nodes
such that the first and last neighbors of a boundary node are boundary
nodes (the first neighbor of an interior node is arbitrary).  In order to
distinguish between interior and boundary nodes, the last neighbor of
each boundary node is represented by the negative of its index.
The indices are 1-based (as in Fortran), not zero based (as in python).

lptr:  Set of 6*(N-2) pointers (indices)
in one-to-one correspondence with the elements of lst.
lst(lptr(i)) indexes the node which follows lst(i) in cyclical
counterclockwise order (the first neighbor follows the last neighbor).
The indices are 1-based (as in Fortran), not zero based (as in python).

lend: N pointers to adjacency lists.
lend(k) points to the last neighbor of node K.  lst(lend(K)) < 0 if and
only if K is a boundary node.
The indices are 1-based (as in Fortran), not zero based (as in python).
 """
        if len(lons.shape) != 1 or len(lats.shape) != 1:
            raise ValueError('lons and lats must be 1d')
        lats = lats.astype(np.float64,copy=False)
        lons = lons.astype(np.float64,copy=False)
        lons = lons.clip(-2.*np.pi,2.*np.pi)
        if (np.abs(lons)).max() > 2.*np.pi:
            msg="lons must be in radians (-2*pi <= lon <= 2*pi)"
            raise ValueError(msg)
        if (np.abs(lats)).max() > 0.5*np.pi:
            msg="lats must be in radians (-pi/2 <= lat <= pi/2)"
            raise ValueError(msg)
        npts = len(lons)
        if len(lats) != npts:
            raise ValueError('lons and lats must have same length')
        # compute cartesian coords on unit sphere.
        x,y,z = _stripack.trans(lats,lons,npts)
        lst,lptr,lend,ierr = _stripack.trmesh(x,y,z,npts)
        if ierr != 0:
            raise ValueError('ierr = %s in trmesh' % ierr)
        self.lons = lons; self.lats = lats; self.npts = npts
        self.x = x; self.y = y; self.z = z
        self.lptr = lptr; self.lst = lst; self.lend = lend
        self._shuffle=False; self._ix=None
    def interp(self,olons,olats,data,order=1):
        """
given a triangulation, perform interpolation on
output mesh defined by olons,olats (in radians), return result in data.
olons, olats can be 1d or 2d (output data array has same shape as olats,lons).
order of interpolation specified by 'order' kwarg, can be 0 (nearest neighbor),
1 (linear), or 3 (hermite cubic).

Algorithms:
 R. J. Renka, "ALGORITHM 623:  Interpolation on the Surface of a
 Sphere", ACM Trans. Math. Software, Vol. 10, No. 4, December 1984,
 pp. 437-439.
"""
        shapeout = olons.shape
        if len(shapeout) not in [1,2]:
            raise ValueError('olons,olats must be 1d or 2d')
        olons1 = (olons.astype(np.float64,copy=False)).ravel()
        olats1 = (olats.astype(np.float64,copy=False)).ravel()
        nptso = len(olons1)
        if (np.abs(olons1)).max() > 2.*np.pi:
            msg="lons must be in radians (-2*pi <= lon <= 2*pi)"
            raise ValueError(msg)
        if (np.abs(olats1)).max() > 0.5*np.pi:
            msg="lats must be in radians (-pi/2 <= lat <= pi/2)"
            raise ValueError(msg)
        if len(olats1) != nptso:
            raise ValueError('lons and lats must have same length')
        if len(data) != self.npts:
            raise ValueError('input data wrong size')
        if order not in [0,1,3]:
            raise ValueError('order must be 0,1 or 3')
        else:
            if self._shuffle:
                data = data[self._ix] # points were randomly shuffled to make triangulation faster
            odata,ierr = \
            _stripack.interp_n(order, olats1, olons1,\
            self.x, self.y, self.z, data.astype(np.float64),\
            self.lst,self.lptr,self.lend,self.npts,nptso)
        if ierr != 0:
            raise ValueError('ierr = %s in intrpc0_n' % ierr)
        return odata.reshape(shapeout)
    def interp_nn(self,olons,olats,data):
        """
same as interp(olons,olats,data,order=0)"""
        return self.interp(olons,olats,data,order=0)
    def interp_linear(self,olons,olats,data):
        """
same as interp(olons,olats,data,order=1)"""
        return self.interp(olons,olats,data,order=1)
    def interp_cubic(self,olons,olats,data):
        """
same as interp(olons,olats,data,order=3)"""
        return self.interp(olons,olats,data,order=3)
    def interp_latlon(self,nlons,data,order=1,quiet=False):
        """
convenience method to generate a reg lat/lon output grid and do interpolation.

nlons:  number of longitudes on output grid
        (a reg lat/lon grid with nlons longitudes and
         nlons/2 latitudes, not including poles or wrap-
         around longitude).
data: 1d array of cubed sphere data
order: order of interpolation (0 for nearest neighbor, 
       1 for linear, 3 for cubic - default is 1).
quiet: Suppress diagnostic messages (default False)

returns lons2d,lats2d,latlon_data where
lons2d: 2d array of longitudes on output grid (in degrees)
lats2d: 2d array of latiitudes on output grid (in degrees)
latlon_data:  2d array of interpolated data on output grid.
        """
        # generate regular 2d lat/lon grid (not including poles,
        # or wrap-around longitude).
        nlats = nlons/2; delta = 360./nlons
        olons = delta*np.arange(nlons)
        olats = -90. + delta*(0.5 + np.arange(nlats))
        olonsd, olatsd = np.meshgrid(olons, olats) # degrees
        # interpolate to the reg lat/lon grid
        if not quiet: t1 = time.time()
        latlon_data = self.interp(np.radians(olonsd),np.radians(olatsd),data,order=order) # expects radians
        if not quiet:
            t = time.time()-t1
            print('time to interpolate to %s by %s lat/lon grid = %s' %
                  (nlons,nlats,t))
            print('shape/min/max',latlon_data.shape, latlon_data.min(), latlon_data.max())
        return olonsd,olatsd,latlon_data

if __name__ == "__main__":
    # test function.
    def _test(npts=14000,nlons=360,tol=[0.2,1.e-2,None,3.e-3]):
        def fibonacci_pts(npts):
            # return lats and lons of N=npts fibonacci grid on a sphere.
            pi = np.pi
            inc = pi * (3.0 - np.sqrt(5.0))
            off = 2. / npts
            lats = []; lons = []
            for k in range(npts):
               y = k*off - 1. + 0.5*off
               r = np.sqrt(1 - y**2)
               phi = k * inc
               x = np.cos(phi)*r
               z = np.sin(phi)*r
               theta = np.arctan2(np.sqrt(x**2+y**2),z)
               phi = np.arctan2(y,x)
               lats.append( 0.5*pi-theta )
               if phi < 0.: phi = 2.*pi+phi
               lons.append( phi )
            return np.array(lats), np.array(lons)
        # fake test data.
        def test_func(lon, lat):
            nexp = 8
            return np.cos(nexp*lon)*np.sin(0.5*lon)**nexp*np.cos(lat)**nexp+np.sin(lat)**nexp
        # function to check error
        def check_err(latlon_data, latlon_datax, order):
            err = (np.abs(latlon_datax-latlon_data)).max()
            print('order = %s: max abs error %s specified tolerance %s' %\
                    (order,err,tol[order]))
            assert(err < tol[order])
        # input mesh (fibonacci spiral points)
        lats, lons = fibonacci_pts(npts)
        icos_data = test_func(lons,lats)
        # triangulation
        tri = trmesh(lons, lats)
        assert( tri.npts == npts )
        assert( np.array_equal(tri.lons, lons) )
        assert( np.array_equal(tri.lats, lats) )
        # output mesh
        nlats = nlons//2; delta = 360./nlons
        olons = delta*np.arange(nlons)
        olats = -90. + delta*(0.5 + np.arange(nlats))
        olons = np.radians(olons);  olats = np.radians(olats)
        olons, olats = np.meshgrid(olons, olats)
        # nearest neighbor interpolation
        order = 0 # can be 0 (nearest neighbor) or 1 (linear) or 3 (cubic)
        latlon_data = tri.interp(olons,olats,icos_data,order=order)
        # check error
        latlon_datax = test_func(olons,olats)
        check_err(latlon_data,test_func(olons,olats),order)
        # test interp_nn alias.
        latlon_data = tri.interp_nn(olons,olats,icos_data)
        latlon_datax = test_func(olons,olats)
        check_err(latlon_data,test_func(olons,olats),order)
        # linear interpolation
        order = 1 # can be 0 (nearest neighbor) or 1 (linear) or 3 (cubic)
        latlon_data = tri.interp(olons,olats,icos_data,order=order)
        # check error
        check_err(latlon_data,test_func(olons,olats),order)
        # test interp_linear alias
        latlon_data = tri.interp_linear(olons,olats,icos_data)
        latlon_datax = test_func(olons,olats)
        check_err(latlon_data,test_func(olons,olats),order)
        # cubic interpolation
        order = 3 # can be 0 (nearest neighbor) or 1 (linear) or 3 (cubic)
        latlon_data = tri.interp(olons,olats,icos_data,order=order)
        # check error
        latlon_datax = test_func(olons,olats)
        check_err(latlon_data,latlon_datax,order)
        # test interp_cubic alias
        latlon_data = tri.interp_cubic(olons,olats,icos_data)
        check_err(latlon_data,latlon_datax,order)
        # check convenience method
        lons2d, lats2d, latlon_data =\
        tri.interp_latlon(nlons,icos_data,quiet=True)
        check_err(latlon_data,latlon_datax,1)
    import unittest
    class RegridTest(unittest.TestCase):
        def test(self):
            _test()
    unittest.main()
