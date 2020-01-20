# http://docs.cython.org/src/userguide/tutorial.html

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("Metropolis", ["Metropolis.pyx"],
                   extra_compile_args=["-O3",])]
)

# To compile
# python setup.py build_ext --inplace
