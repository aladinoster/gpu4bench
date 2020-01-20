#!/usr/bin/env python
#
# PiMC model using mpi4py MPI implementation for Python
#
# CC BY-NC-SA 2011 : <emmanuel.quemener@ens-lyon.fr> 
#
# Thanks to Lisandro Dalcin for MPI4PY :
# http://mpi4py.scipy.org/

import sys
from math import exp
from random import random
import time
# MPI librairie call
import mpi4py
from mpi4py import MPI

def MainLoop(iterations):

    total=0
    for i in xrange(iterations):
        # Random access coordonate
        x,y=random(),random()

	if ((x*x+y*y) < 1.0):
		total+=1

    return(total)

if __name__=='__main__':
	
   # MPI Init
   comm = MPI.COMM_WORLD
   rank = comm.Get_rank()
	
   # Define number of Nodes on with computing is performed (exclude 0)
   NODES=comm.Get_size()

   # Au dessus de 4 10^7, le MPI ne veut plus se lancer...
   Iterations=1000000000

   # pass explicit MPI datatypes
   total=0
   if rank == 0:
	   # Define iterations to send to each node
	   if Iterations%NODES==0:
		   iterations=Iterations/NODES
	   else:
		   iterations=Iterations/NODES+1
	   print "%i iterations will be send to each node" % iterations
	   
	   for i in range(1,NODES):
		   
		   print "Send from 0 to node %i" % i
		   ToSend=iterations
		   # Send MPI call to each node
		   comm.send(ToSend, dest=i, tag=11)

	   # Master does part of job !
	   total=MainLoop(iterations)
	   print "Partial Result from master %i: %i" % (rank,total)

	   print "Retreive results..."
	   for i in range(1,NODES):
		   # Receive MPI call from each node
		   Output=comm.recv(source=i,tag=11)
		   print "Partial Result from %i: %s" % (i,Output)
		   total+=Output

	   print "Global Result: %i" % (total)
   else:
	   # Slave applies simulation to set provided by master
	   # Receive MPI call with Input set
	   ToReceive=comm.recv(source=0, tag=11)
	   iterations=ToReceive
	   print "Rank %i receives with job with %i" % (rank,iterations) 
	   Output=MainLoop(iterations)

	   comm.send(Output, dest=0, tag=11)
