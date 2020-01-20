import numpy
import sys

npzOCLLG=numpy.load("Ising2D_OCLL_256_16777216_GPU_atlas.npz")
npzOCLLC=numpy.load("Ising2D_OCLL_256_16777216_CPU_atlas.npz")
npzCUDAL=numpy.load("Ising2D_CUDAL_256_16777216_GPU_atlas.npz")

npzOCLLG2=numpy.load("Ising2D_OCLL_256_16777216_GPU_o745-2.npz")
npzCUDAL2=numpy.load("Ising2D_CUDAL_256_16777216_GPU_o745-2.npz")


T=npzOCLLG['arr_0'][0]
EOCLLG=npzOCLLG['arr_0'][1]
POCLLG=npzOCLLG['arr_0'][2]
EOCLLC=npzOCLLC['arr_0'][1]
POCLLC=npzOCLLC['arr_0'][2]
ECUDAL=npzCUDAL['arr_0'][1]
PCUDAL=npzCUDAL['arr_0'][2]

EOCLLG2=npzOCLLG2['arr_0'][1]
POCLLG2=npzOCLLG2['arr_0'][2]
ECUDAL2=npzCUDAL2['arr_0'][1]
PCUDAL2=npzCUDAL2['arr_0'][2]

dPOCLLG=numpy.diff(POCLLG)
dPOCLLC=numpy.diff(POCLLC)
dPCUDAL=numpy.diff(PCUDAL)

dPOCLLG2=numpy.diff(POCLLG2)
dPCUDAL2=numpy.diff(PCUDAL2)

dPOCLLG=numpy.insert(dPOCLLG,0,0)
dPOCLLC=numpy.insert(dPOCLLC,0,0)
dPCUDAL=numpy.insert(dPCUDAL,0,0)

dPOCLLG2=numpy.insert(dPOCLLG2,0,0)
dPCUDAL2=numpy.insert(dPCUDAL2,0,0)

from matplotlib.pyplot import *

plot(T,EOCLLG,T,EOCLLC,T,ECUDAL,T,dPOCLLG,T,dPOCLLC,T,dPCUDAL,
     T,EOCLLG2,T,ECUDAL2,T,dPOCLLG2,T,dPCUDAL2)
legend((r'OCL GPU GTX560',r'OCL CPU Xeon4c',r'Cuda GTX560',
        r'OCL GPU GT8400',r'Cuda GT8400'))
xlabel('Temperature',{'fontsize':20})
ylabel('Energy',{'fontsize':20})

matplotlib.pyplot.show()

