#!/usr/bin/env python
#!/usr/bin/env python
#
# Ising2D model using mpi4py MPI implementation for Python
#
# CC BY-NC-SA 2011 : <emmanuel.quemener@ens-lyon.fr> 
#
# Thanks to Lisandro Dalcin for MPI4PY :
# http://mpi4py.scipy.org/

import sys
import numpy
#import pylab
from PIL import Image
from math import exp
from random import random
import time
# MPI librairie call
from mpi4py import MPI
import getopt
import matplotlib.pyplot as plt

LAPIMAGE=False

def partition ( lst, n ):
    return [ lst[i::n] for i in xrange(n) ]

def ImageOutput(sigma,prefix):
    Max=sigma.max()
    Min=sigma.min()

    # Normalize value as 8bits Integer
    SigmaInt=(255*(sigma-Min)/(Max-Min)).astype('uint8')
    image = Image.fromarray(SigmaInt)
    image.save("%s.jpg" % prefix)

def Metropolis(sigma,J,B,T,Iterations): 
    start=time.time()

    SizeX,SizeY=sigma.shape
    
    for p in xrange(0,Iterations):
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

if __name__=='__main__':
	
    ToSave=[]
	
    # MPI Init
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
	
    # Define number of Nodes on with computing is performed (exclude 0)
    NODES=comm.Get_size()-1

    # pass explicit MPI datatypes
    if rank == 0:

        # Set defaults values
        # Coupling factor
        J=1.
        # Magnetic Field
        B=0.
        # Size of Lattice
        Size=256
        # Default Temperatures (start, end, step)
        Tmin=0.1
        Tmax=5.
        Tstep=0.1
        # Default Number of Iterations
        Iterations=Size*Size*Size

        # Curves is True to print the curves
        Curves=False

        try:
            opts, args = getopt.getopt(sys.argv[1:],"hcj:b:z:i:s:e:p:",["coupling=","magneticfield=","size=","Iterations=","tempstart=","tempend=","tempstep=","units"])
        except getopt.GetoptError:
            print '%s -j <Coupling Factor> -b <Magnetic Field> -z <Size of Square Lattice> -i <Iterations> -s <Minimum Temperature> -e <Maximum Temperature> -p <steP Temperature> -c (Print Curves)' % sys.argv[0]
            sys.exit(2)
    
        for opt, arg in opts:
            if opt == '-h':
                   print '%s -j <Coupling Factor> -b <Magnetic Field> -z <Size of Square Lattice> -i <Iterations> -s <Minimum Temperature> -e <Maximum Temperature> -p <steP Temperature> -c (Print Curves)' % sys.argv[0]
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
                                
        print "Coupling Factor J : %s" % J
        print "Magnetic Field B :  %s" % B
        print "Size of lattice : %s" % Size
        print "Iterations : %s" % Iterations
        print "Temperature on start : %s" % Tmin
        print "Temperature on end : %s" % Tmax
        print "Temperature step : %s" % Tstep
        
        LAPIMAGE=False

        sigmaIn=numpy.where(numpy.random.randn(Size,Size)>0,1,-1).astype(numpy.int8)
                
        ImageOutput(sigmaIn,"Ising2D_MPI_%i_Initial" % (Size))

        Trange=numpy.arange(Tmin,Tmax+Tstep,Tstep)

        sigmaIn=numpy.where(numpy.random.randn(Size,Size)>0
                            ,1,-1).astype(numpy.int8)
		
        numpy.random.seed(int(time.time()))

        # Master control distribution of computing
        print "Distributing work to %i node(s)..." % NODES
        Distribution=numpy.array_split(Trange,NODES)
        
        for i in range(NODES):
                
            Input=Distribution[i]                    
            print "Send from 0 to %i %s" % (i+1,Input)
            ToSend=sigmaIn,J,B,Iterations,Input
            # Send MPI call to each node
            comm.send(ToSend, dest=i+1, tag=11)
                
        print "Retreive results..."

        Results=[]
        for i in range(NODES):
            # Receive MPI call from each node
            Output=comm.recv(source=i+1,tag=11)
            print "Result from %i: %s" % (i+1,Output)
            Results+=Output

        E=numpy.array(Results)[:,1]
        M=numpy.array(Results)[:,2]
            
        numpy.savez("Ising2D_MPI_%i_%.8i" % (Size,Iterations),(Trange,E,M))

        # Estimate Critical temperature
        print "The critical temperature on %ix%i lattice with J=%f, B=%f is %f " % (Size,Size,J,B,CriticalT(Trange,E))

        if Curves:
            DisplayCurves(Trange,E,M,J,B)
    else:
        numpy.random.seed(int(time.time()/rank))
        # Slave applies simulation to set provided by master
        # Receive MPI call with Input set
        ToSplit=comm.recv(source=0, tag=11)
        sigmaIn,J,B,Iterations,Input=ToSplit
        print "Rank %i receive with %ix%i lattice at J=%.2f, B=%.2f with %i iterations and T=%s" % (rank,sigmaIn.shape[0],sigmaIn.shape[1],J,B,Iterations,Input) 
        Output=[]
        # Launch simulations on the set, one by one
        for T in Input:
            print "Processing T=%.2f on rank %i" % (T,rank)
            # Reinitialize to original
            sigma=numpy.copy(sigmaIn)
            duration=Metropolis(sigma,J,B,T,Iterations)
            print "CPU Time : %f" % (duration)
            E=Energy(sigma,J,B)
            M=Magnetization(sigma,0.)
            ImageOutput(sigma,"Ising2D_MPI_%i_%1.1f_Final" % 
                        (sigmaIn.shape[0],T))
            Output.append([T,E,M])
        comm.send(Output, dest=0, tag=11)
