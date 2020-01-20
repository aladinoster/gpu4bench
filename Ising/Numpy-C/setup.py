from distutils.core import setup, Extension
import numpy

# define the extension module
array_module_np = Extension('array_module_np', sources=['array_module_np.c'],include_dirs=[numpy.get_include()])

# run the setup
setup(ext_modules=[array_module_np])
