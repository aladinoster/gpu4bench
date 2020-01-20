import numpy
import sys

npzSerial=numpy.load("Ising2D_Serial_256_16777216.npz")
npzMP=numpy.load("Ising2D_MP_256_16777216.npz")
npzMPI=numpy.load("Ising2D_MPI_256_16777216.npz")
npzOCLG=numpy.load("Ising2D_OCLG_256_16777216.npz")
npzOCLL=numpy.load("Ising2D_OCLL_256_16777216.npz")
npzCUDAG=numpy.load("Ising2D_CUDAG_256_16777216.npz")
npzCUDAL=numpy.load("Ising2D_CUDAL_256_16777216.npz")

T=npzSerial['arr_0'][0]
ES=npzSerial['arr_0'][1]
EMP=npzMP['arr_0'].transpose()[1]
EMPI=npzMPI['arr_0'].transpose()[1]
EOCLG=npzOCLG['arr_0'][1]
EOCLL=npzOCLL['arr_0'][1]
ECUDAG=npzCUDAG['arr_0'][1]
ECUDAL=npzCUDAL['arr_0'][1]

import matplotlib.pyplot
matplotlib.pyplot.plot(T,ES,T,EMP,T,EMPI,T,EOCLG,T,EOCLL,T,ECUDAG,T,ECUDAL)
matplotlib.pyplot.show()

