#!/usr/bin/env python
#
# Ising2D model in serial mode
#
# CC BY-NC-SA 2011 : <emmanuel.quemener@ens-lyon.fr> 
#
# Thanks to George Marsaglia : http://en.wikipedia.org/wiki/Multiply-with-carry

import sys
import numpy
#import pylab
from PIL import Image
from math import exp
from random import random
import time

z=0
w=0

def ImageOutput(sigma,prefix):
        Max=sigma.max()
        Min=sigma.min()

        # Normalize value as 8bits Integer
        SigmaInt=(255*(sigma-Min)/(Max-Min)).astype('uint8')
        image = Image.fromarray(SigmaInt)
        image.save("%s.jpg" % prefix)

SIZE=256

def Metropolis(sigma,J,T,iterations,seed_w,seed_z):
	Energy=0
	
	global z
	global w

	z=seed_z
	w=seed_w
	
	for p in range(0,iterations):
		# Random access coordonate
		#X,Y=numpy.random.randint(SIZE),numpy.random.randint(SIZE)
		
		# def znew():
		# 	global z
		# 	z=numpy.uint32(36969*(z&65535)+(z>>16))
		# 	return(z)
		# def wnew():
		# 	global w
		# 	w=numpy.uint32(18000*(w&65535)+(w>>16))
		# 	return(w)
		# def MWC(): return(numpy.uint32((znew()<<16)+wnew()))
		# def MWCfp(): return((MWC() + 1.0) * 2.328306435454494e-10)
		
		def MWC(): 
			global w
			global z
			z=numpy.uint32(36969*(z&65535)+(z>>16))
			w=numpy.uint32(18000*(w&65535)+(w>>16))
			return(numpy.uint32((z<<16)+w))

		def MWCfp(): 
			global w
			global z
			z=numpy.uint32(36969*(z&65535)+(z>>16))
			w=numpy.uint32(18000*(w&65535)+(w>>16))
			return(((numpy.uint32((z<<16)+w))+1.0)*
			       2.328306435454494e-10)

		X=numpy.uint32(SIZE*MWCfp())
		Y=numpy.uint32(SIZE*MWCfp())

		#print X,Y
		
		DeltaE=2.*J*sigma[X,Y]*(sigma[X,(Y+1)%SIZE]+
					sigma[X,(Y-1)%SIZE]+
					sigma[(X-1)%SIZE,Y]+
					sigma[(X+1)%SIZE,Y])
		
		if DeltaE < 0. or random() < exp(-DeltaE/T):
			sigma[X,Y]=-sigma[X,Y]
			Energy+=DeltaE
			
	return(Energy)

if __name__=='__main__':

    J=1.
    T=0.5

    iterations=numpy.uint32(SIZE*SIZE*SIZE)
    
    sigma=numpy.where(numpy.random.randn(SIZE,SIZE)>0,1,-1).astype(numpy.int8)

    ImageOutput(sigma,"Ising2D_OpenCL_%i_Initial" % (SIZE))

    seed_w=numpy.uint32(37)
    seed_z=numpy.uint32(91)
	
    i=0
    step=262144
    while (step*i < iterations):
        start=time.time()
        Metropolis(sigma,J,T,step,seed_w*(i+1),seed_z*((i+1)%SIZE))
        stop=time.time()
        elapsed = (stop-start)
        print "Iteration %i in %f: " % (i,elapsed)
        ImageOutput(sigma,"Ising2D_OpenCL_%i_%.3i_Final" % (SIZE,i))
        i=i+1
    
