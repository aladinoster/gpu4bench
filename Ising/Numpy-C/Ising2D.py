#!/usr/bin/env python
#
# Ising2D model in serial mode using external C library
#
# CC BY-NC-SA 2011 : <emmanuel.quemener@ens-lyon.fr> 

import sys
import numpy
from PIL import Image
from math import exp
from random import random
import time
import getopt
import matplotlib.pyplot as plt
import array_module_np
from numpy.random import randint as nprnd

def ImageOutput(sigma,prefix):
    Max=sigma.max()
    Min=sigma.min()
    
    # Normalize value as 8bits Integer
    SigmaInt=(255*(sigma-Min)/(Max-Min)).astype('uint8')
    image = Image.fromarray(SigmaInt)
    image.save("%s.jpg" % prefix)
    
def Metropolis(sigma,J,B,T,iterations): 
    start=time.time()

    SizeX,SizeY=sigma.shape
    
    for p in xrange(0,iterations):
        # Random access coordonate
        X,Y=numpy.random.randint(SizeX),numpy.random.randint(SizeY)
        
        DeltaE=J*sigma[X,Y]*(2*(sigma[X,(Y+1)%SizeY]+
                                sigma[X,(Y-1)%SizeY]+
                                sigma[(X-1)%SizeX,Y]+
                                sigma[(X+1)%SizeX,Y])+B)
        
        if DeltaE < 0. or random() < exp(-DeltaE/T):
            sigma[X,Y]=-sigma[X,Y]
    duration=time.time()-start
    return(duration)

def Magnetization(sigma,M):
    return(numpy.sum(sigma)/(sigma.shape[0]*sigma.shape[1]*1.0))

def Energy(sigma,J):
    # Copier et caster 
    E=numpy.copy(sigma).astype(numpy.float32)
    
    # Appel par slice
    E[1:-1,1:-1]=-J*E[1:-1,1:-1]*(E[:-2,1:-1]+E[2:,1:-1]+
                                  E[1:-1,:-2]+E[1:-1,2:])
    
    # Bien nettoyer la peripherie
    E[:,0]=0
    E[:,-1]=0
    E[0,:]=0
    E[-1,:]=0
    
    Energy=numpy.sum(E)

    return(Energy/(E.shape[0]*E.shape[1]*1.0))

def DisplayCurves(T,E,M,J,B):

    plt.xlabel("Temperature")
    plt.ylabel("Energy")

    Experience,=plt.plot(T,E,label="Energy") 

    plt.legend()
    plt.show()


if __name__=='__main__':

    # Set defaults values
    # Coupling factor
    J=1.
    # Magnetic Field
    B=0.
    # Size of Lattice
    Size=256
    # Default Temperatures (start, end, step)
    Tmin=0.1
    Tmax=5
    Tstep=0.1
    # Default Number of Iterations
    Iterations=Size*Size

    # Curves is True to print the curves
    Curves=False

    try:
        opts, args = getopt.getopt(sys.argv[1:],"hcj:b:z:i:s:e:p:",["coupling=","magneticfield=","size=","iterations=","tempstart=","tempend=","tempstep="])
    except getopt.GetoptError:
        print '%s -j <Coupling Factor> -b <Magnetic Field> -z <Size of Lattice> -i <Iterations> -s <Minimum Temperature> -e <Maximum Temperature> -p <steP Temperature> -c (Print Curves)' % sys.argv[0]
        sys.exit(2)
    
 
    for opt, arg in opts:
        if opt == '-h':
            print '%s -j <Coupling Factor> -b <Magnetic Field> -z <Size of Lattice> -i <Iterations> -s <Minimum Temperature> -e <Maximum Temperature> -p <steP Temperature> -c (Print Curves)' % sys.argv[0]
            sys.exit()
        elif opt == '-c':
            Curves=True
        elif opt in ("-j", "--coupling"):
            J = float(arg)
        elif opt in ("-b", "--magneticfield"):
            B = float(arg)
        elif opt in ("-s", "--tempmin"):
            Tmin = float(arg)
        elif opt in ("-e", "--tempmax"):
            Tmax = float(arg)
        elif opt in ("-p", "--tempstep"):
            Tstep = float(arg)
        elif opt in ("-i", "--iterations"):
            Iterations = int(arg)
        elif opt in ("-z", "--size"):
            Size = int(arg)
      
    print "Coupling Factor J : %s" % J
    print "Magnetic Field B :  %s" % B
    print "Size of lattice : %s" % Size
    print "Iterations : %s" % Iterations
    print "Temperature on start : %s" % Tmin
    print "Temperature on end : %s" % Tmax 
    print "Temperature step : %s" % Tstep
    Trange=numpy.arange(Tmin,Tmax+Tstep,Tstep)
    print "Temperatures explored : %s" % Trange

    LAPIMAGE=False

    sigmaIn=numpy.where(numpy.random.randn(Size,Size)>0,1,-1).astype(numpy.int32)

    ImageOutput(sigmaIn,"Ising2D_Serial_%i_Initial" % (Size))


    E=[]
    M=[]

    for T in Trange:
        # Indispensable d'utiliser copy : [:] ne fonctionne pas avec numpy !
        sigma=numpy.copy(sigmaIn)
        # duration=Metropolis(sigma,J,B,T,Iterations)
        SeedW,SeedZ=numpy.int32(nprnd(2**31-1)),numpy.int32(nprnd(2**31-1))
        start=time.time()
        array_module_np.array_metropolis_np(sigma,J,B,T,Iterations,SeedW,SeedZ)
        duration=time.time()-start
        E=numpy.append(E,Energy(sigma,J))
        M=numpy.append(M,Magnetization(sigma,B))
        ImageOutput(sigma,"Ising2D_Serial_%i_%1.1f_Final" % (Size,T))
    
        print "CPU Time : %f" % (duration)
        print "Total Energy at Temperature %f : %f" % (T,E[-1])
        print "Total Magnetization at Temperature %f : %f" % (T,M[-1])

    if Curves:
        DisplayCurves(Trange,E,M,J,B)

    # Save output
    numpy.savez("Ising2D_Serial_%i_%.8i" % (Size,Iterations),(Trange,E,M))
    
