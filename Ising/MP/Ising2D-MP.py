#!/usr/bin/env python
#
# Ising2D model in serial mode
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

from multiprocessing import Pool

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
        
        DeltaE=sigma[X,Y]*(2*J*(sigma[X,(Y+1)%SizeY]+
                                sigma[X,(Y-1)%SizeY]+
                                sigma[(X-1)%SizeX,Y]+
                                sigma[(X+1)%SizeX,Y])+B)
        
        if DeltaE < 0. or random() < exp(-DeltaE/T):
            sigma[X,Y]=-sigma[X,Y]
    duration=time.time()-start
    return(duration)

def Magnetization(sigma,M):
    return(numpy.sum(sigma)/(sigma.shape[0]*sigma.shape[1]*1.0))

def Energy(sigma,J,B):
    # Copier et caster 
    E=numpy.copy(sigma).astype(numpy.float32)
    
    # Appel par slice
    E[1:-1,1:-1]=E[1:-1,1:-1]*(2.*J*(E[:-2,1:-1]+E[2:,1:-1]+
                                     E[1:-1,:-2]+E[1:-1,2:])+B)
    
    # Bien nettoyer la peripherie
    E[:,0]=0
    E[:,-1]=0
    E[0,:]=0
    E[-1,:]=0
    
    Energy=numpy.sum(E)

    return(Energy/(E.shape[0]*E.shape[1]*1.0))

def CriticalT(T,E):

    Epoly=numpy.poly1d(numpy.polyfit(T,E,T.size/3))
    dEpoly=numpy.diff(Epoly(T))
    dEpoly=numpy.insert(dEpoly,0,0)
    return(T[numpy.argmin(dEpoly)])

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
    Iterations=Size*Size*Size
    # Default Number of Procs Used
    Procs=4

    # Curves is True to print the curves
    Curves=False

    try:
        opts, args = getopt.getopt(sys.argv[1:],"hcj:b:z:i:s:e:p:u:",["coupling=","magneticfield=","size=","iterations=","tempstart=","tempend=","tempstep=","units"])
    except getopt.GetoptError:
        print '%s -j <Coupling Factor> -b <Magnetic Field> -z <Size of Square Lattice> -i <Iterations> -s <Minimum Temperature> -e <Maximum Temperature> -p <steP Temperature> -u <Units of processing> -c (Print Curves)' % sys.argv[0]
        sys.exit(2)
    
    for opt, arg in opts:
        if opt == '-h':
            print '%s -j <Coupling Factor> -b <Magnetic Field> -z <Size of Square Lattice> -i <Iterations> -s <Minimum Temperature> -e <Maximum Temperature> -p <steP Temperature> -u <Units of Processing> -c (Print Curves)' % sys.argv[0]
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
            Tmax = arg
        elif opt in ("-p", "--tempstep"):
            Tstep = numpy.uint32(arg)
        elif opt in ("-i", "--iterations"):
            Iterations = int(arg)
        elif opt in ("-z", "--size"):
            Size = int(arg)
        elif opt in ("-u", "--units"):
            Procs = int(arg)
      
    print "Process Units : %s" % Procs
    print "Coupling Factor J : %s" % J
    print "Magnetic Field B :  %s" % B
    print "Size of lattice : %s" % Size
    print "Iterations : %s" % Iterations
    print "Temperature on start : %s" % Tmin
    print "Temperature on end : %s" % Tmax
    print "Temperature step : %s" % Tstep

    LAPIMAGE=False

    sigmaIn=numpy.where(numpy.random.randn(Size,Size)>0,1,-1).astype(numpy.int8)

    ImageOutput(sigmaIn,"Ising2D_Serial_%i_Initial" % (Size))

    Trange=numpy.arange(Tmin,Tmax+Tstep,Tstep)

    def MetropolisStrip(T):
        # Indispensable d'utiliser copy : [:] ne fonctionne pas avec numpy !
        sigma=numpy.copy(sigmaIn)
        duration=Metropolis(sigma,J,B,T,Iterations)
        ImageOutput(sigma,"Ising2D_MP_%i_%1.1f_Final" % (Size,T))
        print "CPU Time : %f" % (duration)
        indice=numpy.where(Trange==T)[0][0]
        E,M=Energy(sigma,J,B),Magnetization(sigma,B)
        print "Total Energy at Temperature %f : %f" % (T,E)
        print "Total Magnetization at Temperature %f : %f" % (T,M)
        return([T,E,M])

    pool=Pool(processes=Procs)
    print "Start on %i processors..." % Procs
    # Apply MAP to POOL
    # MAP: distribution of usecases T to be applied to MetropolisStrip 
    Results=pool.map(MetropolisStrip,Trange)

    E=numpy.array(Results)[:,1]
    M=numpy.array(Results)[:,2]
    
    # Save output
    numpy.savez("Ising2D_MP_%i_%.8i" % (Size,Iterations),(Trange,E,M))
      
    # Estimate Critical temperature
    print "The critical temperature on %ix%i lattice with J=%f, B=%f is %f " % (Size,Size,J,B,CriticalT(Trange,E))
  
    if Curves:
        DisplayCurves(Trange,E,M,J,B)

