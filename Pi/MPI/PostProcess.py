import numpy
npzfile=numpy.load("Ising2D_MPI_256_16777216.npz")
npzfile['arr_0']
J=npzfile['arr_0'].transpose()[0]
E=npzfile['arr_0'].transpose()[1]
M=npzfile['arr_0'].transpose()[2]
import matplotlib.pyplot 
matplotlib.pyplot.plot(J,E)
matplotlib.pyplot.show()
