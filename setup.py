from numpy.distutils.core  import setup, Extension

# interface for Renka's algorithm 772 fortran code
ext = Extension(name  = '_stripack',
                sources       = ['_stripack.pyf','stripack.f90'])

if __name__ == "__main__":
    setup(name = 'stripack',
          version           = "1.0",
          description       = "Python interface to TOMS 772 (STRIPACK) fortran code",
          author            = "Jeff Whitaker",
          author_email      = "jeffrey.s.whitaker@noaa.gov",
          ext_modules       = [ext],
          packages          = ['stripack'],
          )
