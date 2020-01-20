import numpy
import sys
npzfile=numpy.load(sys.argv[1])
T=npzfile['arr_0'][0]
E=npzfile['arr_0'][1]
ED=npzfile['arr_0'][2]
import matplotlib.pyplot
matplotlib.pyplot.plot(T,E)
#matplotlib.pyplot.plot(T,E,T,ED)
matplotlib.pyplot.show()
