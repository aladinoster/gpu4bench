import numpy
npzfile=numpy.load("Ising2D_MP_256_16777216.npz")
npzfile['arr_0']
T=npzfile['arr_0'].transpose()[0]
E=npzfile['arr_0'].transpose()[1]
M=npzfile['arr_0'].transpose()[2]
print T
print E
print M
import matplotlib.pyplot 
matplotlib.pyplot.plot(T,E)
matplotlib.pyplot.show()
