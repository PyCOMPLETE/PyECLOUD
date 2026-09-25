"""Native extension builds; package metadata lives in pyproject.toml."""

from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import numpy

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext
from setuptools.errors import CompileError


class F2PyBuildExt(build_ext):
    """Let setuptools manage placement and editable installs of F2PY modules."""

    def build_extension(self, ext):
        if not all(Path(source).suffix == '.f' for source in ext.sources):
            return super().build_extension(ext)
        module_name = ext.name.rsplit('.', 1)[-1]
        sources = [str(Path(source).resolve()) for source in ext.sources]
        destination = Path(self.get_ext_fullpath(ext.name)).resolve()
        build_dir = Path(self.build_temp).resolve()
        build_dir.mkdir(parents=True, exist_ok=True)
        # A fresh directory prevents an old binary from masking a failed build.
        with tempfile.TemporaryDirectory(prefix=ext.name + '-', dir=build_dir) as tmp:
            command = [sys.executable, '-m', 'numpy.f2py', '-c',
                       '--backend', 'meson', '-m', module_name, *sources]
            try:
                subprocess.run(command, cwd=tmp, check=True)
            except subprocess.CalledProcessError as exc:
                raise CompileError(
                    f'Failed to build {ext.name}. A C and Fortran compiler '
                    'must be available; see README.md for installation details.'
                ) from exc
            binary = Path(tmp) / self.get_ext_filename(module_name)
            if not binary.is_file():
                raise CompileError(f'F2PY did not produce the expected file: {binary}')
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(binary, destination)


extensions = [
    Extension('PyECLOUD.' + name, ['PyECLOUD/fortran/' + source])
    for name, source in [
        ('rhocompute', 'compute_rho.f'),
        ('int_field_for', 'interp_field_for.f'),
        ('hist_for', 'compute_hist.f'),
        ('seg_impact', 'update_seg_impact.f'),
        ('boris_step', 'boris_step.f'),
        ('vectsum', 'vectsum.f'),
        ('errffor', 'errfff.f'),
    ]
]

extensions += [
    Extension('PyECLOUD.boris_cython',
              ['PyECLOUD/boris_cython.pyx', 'PyECLOUD/boris_c_function.c'],
              depends=['PyECLOUD/boris_cython.pxd', 'PyECLOUD/boris_c_function.h'],
              include_dirs=[numpy.get_include(), 'PyECLOUD']),
    Extension('PyECLOUD.geom_impact_poly_cython',
              ['PyECLOUD/geom_impact_poly_cython.pyx'],
              include_dirs=[numpy.get_include()]),
]

setup(ext_modules=extensions, cmdclass={'build_ext': F2PyBuildExt})
