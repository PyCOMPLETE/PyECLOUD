import numpy
from setuptools import setup, Extension
from Cython.Build import cythonize

setup(
 ext_modules=cythonize([Extension("boris_cython", ["boris_cython.pyx", 'boris_c_function.c'],
                                  include_dirs=[numpy.get_include()]),
                        Extension("geom_impact_poly_cython", ["geom_impact_poly_cython.pyx"],
                                  include_dirs=[numpy.get_include()], annotate=True)]))

