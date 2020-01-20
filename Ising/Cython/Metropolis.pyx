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
import Metropolis

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
