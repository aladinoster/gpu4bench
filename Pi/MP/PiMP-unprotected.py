#!/usr/bin/env python

# Version using Multiprocessing module
#
# Pi-by-MC
#
# CC BY-NC-SA 2013 : <emmanuel.quemener@ens-lyon.fr> 
#

import getopt
import numpy

from random import random
# Multithread library call
from multiprocessing import Pool

# Common tools
import sys
import getopt
import time
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit
from socket import gethostname

# Predicted Amdahl Law (Reduced with s=1-p)  
def AmdahlR(N, T1, p):
    return (T1*(1-p+p/N))

# Predicted Amdahl Law
def Amdahl(N, T1, s, p):
    return (T1*(s+p/N))

# Predicted Mylq Law with first order
def Mylq(N, T1,s,c,p):
    return (T1*(s+c*N+p/N))

# Predicted Mylq Law with second order
def Mylq2(N, T1,s,c1,c2,p):
    return (T1*(s+c1*N+c2*N*N+p/N))

def MainLoop(iterations):
    
    total=0
    # Absulute necessary to use xrange instead of range to avoid out of memory !
    for i in xrange(iterations):
        # Random access coordonate
        x,y=random(),random()
        
        if ((x*x+y*y) < 1.0):
            total+=1
            
    return(total)
        
def MetropolisMP(circle,iterations,steps,jobs):

    MyPi=numpy.zeros(steps)
    MyDuration=numpy.array([])
    
    # Define iterations to send to each node
    if iterations%jobs==0:
        iterationsMP=iterations/jobs
    else:
        iterationsMP=iterations/jobs+1
    print "%i iterations will be send to each core" % iterationsMP
    
    WorkItems=[]
    for i in xrange(0,jobs):
        WorkItems.append(iterationsMP)

    for i in xrange(steps):
        start=time.time()
        pool=Pool(processes=jobs)
        # Define the number of processes to be launched at a time
        # pool=Pool(processes=CORES)
        print "Start on %i processors..." % jobs
        # MAP: distribution of usecases T to be applied to MetropolisStrip 
        # WorkLaunches=[ pool.apply_async(MainLoop,(wi,)) for wi in WorkItems ]
        # circle=[wl.get() for wl in WorkLaunches]
        circle=pool.map(MainLoop,WorkItems)
        MyDuration=numpy.append(MyDuration,time.time()-start)
        MyPi[i]=4.*float(numpy.sum(circle))/float(iterationsMP)/float(jobs)

    return(numpy.mean(MyDuration),numpy.median(MyDuration),numpy.std(MyDuration))

def FitAndPrint(N,D,Curves):

    try:
      coeffs_Amdahl, matcov_Amdahl = curve_fit(Amdahl, N, D)

      D_Amdahl=Amdahl(N,coeffs_Amdahl[0],coeffs_Amdahl[1],coeffs_Amdahl[2])
      coeffs_Amdahl[1]=coeffs_Amdahl[1]*coeffs_Amdahl[0]/D[0]
      coeffs_Amdahl[2]=coeffs_Amdahl[2]*coeffs_Amdahl[0]/D[0]
      coeffs_Amdahl[0]=D[0]
      print "Amdahl Normalized: T=%.2f(%.6f+%.6f/N)" % \
            (coeffs_Amdahl[0],coeffs_Amdahl[1],coeffs_Amdahl[2])
    except:
        print "Impossible to fit for Amdahl law : only %i elements" % len(D) 

    try:
        coeffs_AmdahlR, matcov_AmdahlR = curve_fit(AmdahlR, N, D)

        D_AmdahlR=AmdahlR(N,coeffs_AmdahlR[0],coeffs_AmdahlR[1])
        coeffs_AmdahlR[1]=coeffs_AmdahlR[1]*coeffs_AmdahlR[0]/D[0]
        coeffs_AmdahlR[0]=D[0]
        print "Amdahl Reduced Normalized: T=%.2f(%.6f+%.6f/N)" % \
              (coeffs_AmdahlR[0],1-coeffs_AmdahlR[1],coeffs_AmdahlR[1])
    except:
        print "Impossible to fit for Reduced Amdahl law : only %i elements" % len(D) 

    try:
        coeffs_Mylq, matcov_Mylq = curve_fit(Mylq, N, D)

        coeffs_Mylq[1]=coeffs_Mylq[1]*coeffs_Mylq[0]/D[0]
        coeffs_Mylq[2]=coeffs_Mylq[2]*coeffs_Mylq[0]/D[0]
        coeffs_Mylq[3]=coeffs_Mylq[3]*coeffs_Mylq[0]/D[0]
        coeffs_Mylq[0]=D[0]
        print "Mylq Normalized : T=%.2f(%.6f+%.6f*N+%.6f/N)" % (coeffs_Mylq[0],
                                                                coeffs_Mylq[1],
                                                                coeffs_Mylq[2],
                                                                coeffs_Mylq[3])
        D_Mylq=Mylq(N,coeffs_Mylq[0],coeffs_Mylq[1],coeffs_Mylq[2],
                    coeffs_Mylq[3])
    except:
        print "Impossible to fit for Mylq law : only %i elements" % len(D) 

    if Curves:
        plt.xlabel("Number of Threads/work Items")
        plt.ylabel("Total Elapsed Time")

        Experience,=plt.plot(N,D,'ro') 
        try:
            pAmdahl,=plt.plot(N,D_Amdahl,label="Loi de Amdahl")    
            pMylq,=plt.plot(N,D_Mylq,label="Loi de Mylq")
        except:
            print "Fit curves seem not to be available"

        plt.legend()
        plt.show()

if __name__=='__main__':

    # Set defaults values
    # Iterations is integer
    Iterations=1000000
    # JobStart in first number of Jobs to explore
    JobStart=1
    # JobEnd is last number of Jobs to explore
    JobEnd=1
    # Redo is the times to redo the test to improve metrology
    Redo=1
    # OutMetrology is method for duration estimation : False is GPU inside
    OutMetrology=False
    # Curves is True to print the curves
    Curves=False

    try:
        opts, args = getopt.getopt(sys.argv[1:],"hoci:s:e:r:",["alu=","gpustyle=","parastyle=","iterations=","jobstart=","jobend=","redo=","device="])
    except getopt.GetoptError:
        print '%s -o (Out of Core Metrology) -c (Print Curves) -i <Iterations> -s <JobStart> -e <JobEnd> -r <RedoToImproveStats>' % sys.argv[0]
        sys.exit(2)
    
    for opt, arg in opts:
        if opt == '-h':
            print '%s -o (Out of Core Metrology) -c (Print Curves) -i <Iterations> -s <JobStart> -e <JobEnd> -r <RedoToImproveStats>' % sys.argv[0]
            sys.exit()
        elif opt == '-o':
            OutMetrology=True
        elif opt == '-c':
            Curves=True
        elif opt in ("-i", "--iterations"):
            Iterations = numpy.uint32(arg)
        elif opt in ("-s", "--jobstart"):
            JobStart = int(arg)
        elif opt in ("-e", "--jobend"):
            JobEnd = int(arg)
        elif opt in ("-r", "--redo"):
            Redo = int(arg)

    print "Iterations : %s" % Iterations
    print "Number of threads on start : %s" % JobStart
    print "Number of threads on end : %s" % JobEnd
    print "Number of redo : %s" % Redo
    print "Metrology done out : %r" % OutMetrology

    average=numpy.array([]).astype(numpy.float32)
    median=numpy.array([]).astype(numpy.float32)
    stddev=numpy.array([]).astype(numpy.float32)

    ExploredJobs=numpy.array([]).astype(numpy.uint32)

    Jobs=JobStart

    while Jobs <= JobEnd:
        avg,med,std=0,0,0
        ExploredJobs=numpy.append(ExploredJobs,Jobs)
        circle=numpy.zeros(Jobs).astype(numpy.uint32)
        
        if OutMetrology:
            duration=numpy.array([]).astype(numpy.float32)
            for i in range(Redo):
                start=time.time()
                MetropolisMP(circle,Iterations,Redo,Jobs)
                duration=numpy.append(duration,time.time()-start)
            avg=numpy.mean(duration)
            med=numpy.median(duration)
            std=numpy.std(duration)
        else:
            avg,med,std=MetropolisMP(circle,Iterations,Redo,Jobs)

        if (avg,med,std) != (0,0,0):
            print "avg,med,std",avg,med,std
            average=numpy.append(average,avg)
            median=numpy.append(median,med)
            stddev=numpy.append(stddev,std)
        else:
            print "Values seem to be wrong..."
        # THREADS*=2
        numpy.savez("Pi_%s_%i_%.8i_%s" % (JobStart,JobEnd,Iterations,gethostname()),(ExploredJobs,average,median,stddev))
        Jobs+=1

    FitAndPrint(ExploredJobs,median,Curves)
