import numpy
import sys

npzSerial=numpy.load("Ising2D_Serial_256_16777216.npz")
npzMP=numpy.load("Ising2D_MP_256_16777216.npz")
npzMPI=numpy.load("Ising2D_MPI_256_16777216.npz")

T=npzSerial['arr_0'][0]
ES=npzSerial['arr_0'][1]
EMP=npzMP['arr_0'].transpose()[1]
EMPI=npzMPI['arr_0'].transpose()[1]

import matplotlib.pyplot
matplotlib.pyplot.plot(T,ES,T,EMP,T,EMPI)
matplotlib.pyplot.show()

