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
          version           = "1.0",
          description       = "Python interface to TOMS 772 (STRIPACK) fortran code",
          ext_modules       = [ext],
          packages          = ['stripack'],
          )
