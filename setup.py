from numpy.distutils.core  import setup, Extension

# interface for Renka's algorithm 772 fortran code
ext = Extension(name  = '_stripack',
                sources       = ['_stripack.pyf','stripack.f90'])

if __name__ == "__main__":
    setup(name = 'stripack',
          author            = "Jeff Whitaker",
          author_email      = "jeffrey.s.whitaker@noaa.gov",
          url               = "https://github.com/jswhit/stripack",
          download_url      = "https://github.com/jswhit/stripack/releases",
          version           = "1.2",
          description       = "Python interface to STRIPACK fortran code for triangulation/interpolation on a sphere",
          ext_modules       = [ext],
          packages          = ['stripack'],
          )
