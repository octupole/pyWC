# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
"""pyWC: Willard–Chandler surface analysis toolkit"""

# To use a consistent encoding
import codecs
import os
import sys
import shutil

# Always prefer setuptools over distutils
try:
    from setuptools import find_packages
    from Cython.Distutils import build_ext
    import numpy
    import pybind11
except ImportError as mod_error:
    mod_name = mod_error.message.split()[3]
    sys.stderr.write("Error : " + mod_name + " is not installed\n"
                     "Use pip install " + mod_name + "\n")
    exit(100)

from setuptools import setup
from distutils.extension import Extension

pywc_dbscan = Extension(
    "pywc_dbscan", ["pywc/dbscan_inner.pyx"],
    language="c++",
    include_dirs=[numpy.get_include()])

compile_args = ["-O3", "-std=c++17"]
link_args = []
if sys.platform.startswith("win"):
    compile_args = ["/O2", "/std:c++17", "/openmp"]
elif sys.platform == "darwin":
    compile_args += ["-Xpreprocessor", "-fopenmp"]
    link_args += ["-lomp"]
else:
    compile_args += ["-fopenmp"]
    link_args += ["-fopenmp"]

wc_kde = Extension(
    "pywc._wc_kde",
    ["pywc/_wc_kde.cpp"],
    language="c++",
    include_dirs=[numpy.get_include(), pybind11.get_include()],
    extra_compile_args=compile_args,
    extra_link_args=link_args,
)

center_fast = Extension(
    "pywc.center_fast",
    ["pywc/center_fast.cpp"],
    language="c++",
    include_dirs=[numpy.get_include(), pybind11.get_include()],
    extra_compile_args=compile_args,
    extra_link_args=link_args,
)

center_fast_full = Extension(
    "pywc.center_fast_full",
    ["pywc/center_fast_full.cpp"],
    language="c++",
    include_dirs=[numpy.get_include(), pybind11.get_include()],
    extra_compile_args=compile_args,
    extra_link_args=link_args,
)

here = os.path.abspath(os.path.dirname(__file__))

# Get the long description from the README file
with codecs.open(os.path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

# This fixes the default architecture flags of Apple's python
if sys.platform == 'darwin' and os.path.exists('/usr/bin/xcodebuild'):
    os.environ['ARCHFLAGS'] = ''

# Get version from the file version.py
version = {}
with open("pywc/version.py") as fp:
    exec(fp.read(), version)

BASE_DEPENDENCIES = [
    "numpy>=2.1.3",
    "scipy>=1.11.3",
    "scikit-image>=0.24.0",
    "MDAnalysis>=2.8.0",
    "packaging>=23.0",
]

EXTRAS_REQUIRE = {
    "gpu": ["cupy>=12.0"],
    "dev": ["nose>=1.3.7", "coverage"],
}


def _cuda_available() -> bool:
    """Return True if a CUDA toolkit/driver seems to be present."""
    cuda_env = os.environ.get("CUDA_HOME") or os.environ.get("CUDA_PATH")
    if cuda_env and os.path.isdir(cuda_env):
        return True

    candidate_bins = ["nvcc", "nvidia-smi"]
    if any(shutil.which(cmd) for cmd in candidate_bins):
        return True

    candidate_dirs = ("/usr/local/cuda", "/opt/cuda")
    return any(os.path.isdir(path) for path in candidate_dirs)


def _should_require_cupy() -> bool:
    """Determine whether CuPy should be auto-required."""
    skip_flag = os.environ.get("PYWC_SKIP_CUPY", "").lower()
    if skip_flag in {"1", "true", "yes"}:
        return False
    return _cuda_available()


install_requires = BASE_DEPENDENCIES.copy()
if _should_require_cupy():
    cupy_pkg = os.environ.get("PYWC_CUPY_PACKAGE")
    gpu_requires = [cupy_pkg] if cupy_pkg else list(EXTRAS_REQUIRE["gpu"])
    EXTRAS_REQUIRE["gpu"] = gpu_requires
    install_requires.extend(gpu_requires)

setup(
    name='pywc',
    ext_modules=[pywc_dbscan, wc_kde, center_fast, center_fast_full],
    cmdclass={
        'build_ext': build_ext,
    },
    
    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=find_packages(),
    install_requires=install_requires,
    extras_require=EXTRAS_REQUIRE,

    # List additional groups of dependencies here (e.g. development
    # dependencies). You can install these using the following syntax,
    # for example:
    # If there are data files included in your packages that need to be
    # installed, specify them here.  If using Python 2.6 or less, then these
    # have to be included in MANIFEST.in as well.
    package_data={
        'pywc': ['data/*'],
    },

)
