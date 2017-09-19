from __future__ import division, print_function
import numpy
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

setup(
 ext_modules = cythonize([Extension("boris_cython", ["boris_cython.pyx", 'boris_c_function.c'],
                                include_dirs=[numpy.get_include()]),
                             Extension("geom_impact_poly_cython", ["geom_impact_poly_cython.pyx"],
                                include_dirs=[numpy.get_include()], annotate=True)]))

