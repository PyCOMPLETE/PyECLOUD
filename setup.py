from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

setup(
    ext_modules = cythonize([Extension("boris_cython", ["boris_cython.pyx", 'boris_c_function.c']), 
								Extension("geom_impact_poly_cython", ["geom_impact_poly_cython.pyx"])]
								, annotate=True)
)
