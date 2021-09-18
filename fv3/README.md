tools for interpolating native FV3 cubed sphere tile grids to lat/lon

* generate_fv3mesh.py:  pre-generate stripack mesh triangulation objects and store in pickle files.
                        edit to change path to FV3 C###_oro_data.tile#.nc files before running.
                        takes resolution as command line arg ('python generate_fv3mesh.py 768' generates
                        C786_grid.pickle).
* generate_fv3mesh.sh:  batch script to generated pickled mesh objects on hera.
* test_interp.py:       reads orography from C###_oro_data.tile#.nc, interpolates to 1/4 deg
                        grid using pickled triangulation object generated by generate_fv3mesh.py
                        and makes a plot. May need to edit to change path to C###_oro_data.tile#.nc files.