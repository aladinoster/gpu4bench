#!/usr/bin/env python3

#
# Pi-by-MonteCarlo using PyCUDA/PyOpenCL
#
# CC BY-NC-SA 2011 : Emmanuel QUEMENER <emmanuel.quemener@gmail.com> 
# Cecill v2 : Emmanuel QUEMENER <emmanuel.quemener@gmail.com>
#
# Thanks to Andreas Klockner for PyCUDA:
# http://mathema.tician.de/software/pycuda
# Thanks to Andreas Klockner for PyOpenCL:
# http://mathema.tician.de/software/pyopencl
# 

# 2013-01-01 : problems with launch timeout
# http://stackoverflow.com/questions/497685/how-do-you-get-around-the-maximum-cuda-run-time
# Option "Interactive" "0" in /etc/X11/xorg.conf

# Common tools
import numpy
from numpy.random import randint as nprnd
import sys
import getopt
import time
import math
import itertools
from socket import gethostname

import mpi4py
from mpi4py import MPI

from PiXPU import *

if __name__=='__main__':

    # MPI Init
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
	
    # Define number of Nodes on with computing is performed (exclude 0)
    RankSize=comm.Get_size()

    if rank == 0:
    
        # Set defaults values
  
        # Id of Device : 1 is for first find !
        Device=1
        # GPU style can be Cuda (Nvidia implementation) or OpenCL
        GpuStyle='OpenCL'
        # Iterations is integer
        Iterations=10000000
        # BlocksBlocks in first number of Blocks to explore
        BlocksBegin=1
        # BlocksEnd is last number of Blocks to explore
        BlocksEnd=16
        # BlocksStep is the step of Blocks to explore
        BlocksStep=1
        # ThreadsBlocks in first number of Blocks to explore
        ThreadsBegin=1
        # ThreadsEnd is last number of Blocks to explore
        ThreadsEnd=1
        # ThreadsStep is the step of Blocks to explore
        ThreadsStep=1
        # Redo is the times to redo the test to improve metrology
        Redo=1
        # OutMetrology is method for duration estimation : False is GPU inside
        OutMetrology=False
        Metrology='InMetro'
        # Curves is True to print the curves
        Curves=False
        # Fit is True to print the curves
        Fit=False
        # Marsaglia RNG
        RNG='MWC'
        # Value type : INT32, INT64, FP32, FP64
        ValueType='FP32'
        # Inside based on If
        IfThen=False
        
        HowToUse='%s -c (Print Curves) -k (Case On IfThen) -d <DeviceId> -g <CUDA/OpenCL> -i <Iterations> -b <BlocksBegin> -e <BlocksEnd> -s <BlocksStep> -f <ThreadsFirst> -l <ThreadsLast> -t <ThreadssTep> -r <RedoToImproveStats> -m <SHR3/CONG/MWC/KISS> -v <INT32/INT64/FP32/FP64>'
    
        try:
            opts, args = getopt.getopt(sys.argv[1:],"hckg:i:b:e:s:f:l:t:r:d:m:v:",["gpustyle=","iterations=","blocksBegin=","blocksEnd=","blocksStep=","threadsFirst=","threadsLast=","threadssTep=","redo=","device=","marsaglia=","valuetype="])
        except getopt.GetoptError:
            print(HowToUse % sys.argv[0])
            sys.exit(2)

        # List of Devices
        Devices=[]
        Alu={}
        
        for opt, arg in opts:
            if opt == '-h':
                print(HowToUse % sys.argv[0])

                print("\nInformations about devices detected under OpenCL:")
                # For PyOpenCL import
                try:
                    import pyopencl as cl
                    Id=0
                    for platform in cl.get_platforms():
                        for device in platform.get_devices():
                            #deviceType=cl.device_type.to_string(device.type)
                            deviceType="xPU"
                            print("Device #%i from %s of type %s : %s" % (Id,platform.vendor.lstrip(),deviceType,device.name.lstrip()))
                            Id=Id+1

                    print
                except:
                    print("Your platform does not seem to support OpenCL")

                print("\nInformations about devices detected under CUDA API:")
                # For PyCUDA import
                try:
                    import pycuda.driver as cuda
                    cuda.init()
                    for Id in range(cuda.Device.count()):
                        device=cuda.Device(Id)
                        print("Device #%i of type GPU : %s" % (Id,device.name()))
                    print
                except:
                    print("Your platform does not seem to support CUDA")
        
                sys.exit()
                    
            elif opt == '-c':
                Curves=True
            elif opt == '-k':
                IfThen=True
            elif opt in ("-d", "--device"):
                Devices.append(int(arg))
            elif opt in ("-g", "--gpustyle"):
                GpuStyle = arg
            elif opt in ("-m", "--marsaglia"):
                RNG = arg
            elif opt in ("-v", "--valuetype"):
                ValueType = arg
            elif opt in ("-i", "--iterations"):
                Iterations = numpy.uint64(arg)
            elif opt in ("-b", "--blocksbegin"):
                BlocksBegin = int(arg)
                BlocksEnd = BlocksBegin
            elif opt in ("-e", "--blocksend"):
                BlocksEnd = int(arg)
            elif opt in ("-s", "--blocksstep"):
                BlocksStep = int(arg)
            elif opt in ("-f", "--threadsfirst"):
                ThreadsBegin = int(arg)
                ThreadsEnd = ThreadsBegin
            elif opt in ("-l", "--threadslast"):
                ThreadsEnd = int(arg)
            elif opt in ("-t", "--threadsstep"):
                ThreadsStep = int(arg)
            elif opt in ("-r", "--redo"):
                Redo = int(arg)

        print("Devices Identification : %s" % Devices)
        print("GpuStyle used : %s" % GpuStyle)
        print("Iterations : %s" % Iterations)
        print("Number of Blocks on begin : %s" % BlocksBegin)
        print("Number of Blocks on end : %s" % BlocksEnd)
        print("Step on Blocks : %s" % BlocksStep)
        print("Number of Threads on begin : %s" % ThreadsBegin)
        print("Number of Threads on end : %s" % ThreadsEnd)
        print("Step on Threads : %s" % ThreadsStep)
        print("Number of redo : %s" % Redo)
        print("Metrology done out of XPU : %r" % OutMetrology)
        print("Type of Marsaglia RNG used : %s" % RNG)
        print("Type of variable : %s" % ValueType)

        if GpuStyle=='CUDA':
            try:
                # For PyCUDA import
                import pycuda.driver as cuda

                cuda.init()
                for Id in range(cuda.Device.count()):
                    device=cuda.Device(Id)
                    print("Device #%i of type GPU : %s" % (Id,device.name()))
                    if Id in Devices:
                        Alu[Id]='GPU'
            except ImportError:
                print("Platform does not seem to support CUDA")

        if GpuStyle=='OpenCL':
            try:
                # For PyOpenCL import
                import pyopencl as cl
                Id=0
                for platform in cl.get_platforms():
                    for device in platform.get_devices():
                        #deviceType=cl.device_type.to_string(device.type)
                        deviceType="*PU"
                        print("Device #%i from %s of type %s : %s" % (Id,platform.vendor.lstrip().rstrip(),deviceType,device.name.lstrip().rstrip()))

                        if Id in Devices:
                            # Set the Alu as detected Device Type
                            Alu[Id]=deviceType
                        Id=Id+1
            except ImportError:
                print("Platform does not seem to support OpenCL")

        print(Devices,Alu)
            
        BlocksList=range(BlocksBegin,BlocksEnd+BlocksStep,BlocksStep)
        ThreadsList=range(ThreadsBegin,ThreadsEnd+ThreadsStep,ThreadsStep)
    
        ExploredJobs=numpy.array([]).astype(numpy.uint32)
        ExploredBlocks=numpy.array([]).astype(numpy.uint32)
        ExploredThreads=numpy.array([]).astype(numpy.uint32)
        avgD=numpy.array([]).astype(numpy.float32)
        medD=numpy.array([]).astype(numpy.float32)
        stdD=numpy.array([]).astype(numpy.float32)
        minD=numpy.array([]).astype(numpy.float32)
        maxD=numpy.array([]).astype(numpy.float32)
        avgR=numpy.array([]).astype(numpy.float32)
        medR=numpy.array([]).astype(numpy.float32)
        stdR=numpy.array([]).astype(numpy.float32)
        minR=numpy.array([]).astype(numpy.float32)
        maxR=numpy.array([]).astype(numpy.float32)

        IterationsMPI=numpy.uint64(Iterations/len(Devices))
        if Iterations%len(Devices)!=0:
            IterationsMPI+=1

        for Blocks,Threads in itertools.product(BlocksList,ThreadsList):
        
            ExploredJobs=numpy.append(ExploredJobs,Blocks*Threads)
            ExploredBlocks=numpy.append(ExploredBlocks,Blocks)
            ExploredThreads=numpy.append(ExploredThreads,Threads)

            DurationItem=numpy.array([]).astype(numpy.float32)
            Duration=numpy.array([]).astype(numpy.float32)
            Rate=numpy.array([]).astype(numpy.float32)
            for i in range(Redo):
                time_start=time.time()

                r=1
                # Distribution of Devices over nodes
                InputCL={}
                InputCL['Iterations']=IterationsMPI
                InputCL['Steps']=1
                InputCL['Blocks']=Blocks
                InputCL['Threads']=Threads
                InputCL['RNG']=RNG
                InputCL['ValueType']=ValueType
                InputCL['GpuStyle']=GpuStyle
                InputCL['IfThen']=IfThen

                for Device in Devices[1:]:
                    print("Send to device %i on rank %i" % (Device,r))
                    InputCL['Device']=Device
                    comm.send('CONTINUE',dest=r,tag=11)
                    comm.send(InputCL,dest=r,tag=11)
                    r+=1

                # Compute on rank 0
                print("Compute on rank 0")
                InputCL['Device']=Devices[0]

                if GpuStyle=='CUDA':
                    try:
                        OutputCL=MetropolisCuda(InputCL)
                    except:
                        print("Problem with (%i,%i) // computations on Cuda" % (Blocks,Threads))
                elif GpuStyle=='OpenCL':
                    try:
                        OutputCL=MetropolisOpenCL(InputCL)
                    except:
                        print("Problem with (%i,%i) // computations on OpenCL" %  (Blocks,Threads))

                Inside=OutputCL['Inside']
                NewIterations=OutputCL['NewIterations']

                for slave in range(1,len(Devices)):
                    print("Get OutputCL from %i" % slave)
                    OutputCL=comm.recv(source=slave,tag=11)
                    print(OutputCL)
                    NewIterations+=OutputCL['NewIterations']
                    Inside+=OutputCL['Inside']

                print("Pi estimation %.8f" % (4./NewIterations*Inside))
                
                Duration=numpy.append(Duration,time.time()-time_start)
                Rate=numpy.append(Rate,NewIterations/Duration[-1])
                                            
            avgD=numpy.append(avgD,numpy.average(Duration))
            medD=numpy.append(medD,numpy.median(Duration))
            stdD=numpy.append(stdD,numpy.std(Duration))
            minD=numpy.append(minD,numpy.min(Duration))
            maxD=numpy.append(maxD,numpy.max(Duration))
            avgR=numpy.append(avgR,numpy.average(Rate))
            medR=numpy.append(medR,numpy.median(Rate))
            stdR=numpy.append(stdR,numpy.std(Rate))
            minR=numpy.append(minR,numpy.min(Rate))
            maxR=numpy.append(maxR,numpy.max(Rate))

            print("%.2f %.2f %.2f %.2f %.2f %i %i %i %i %i" % (avgD[-1],medD[-1],stdD[-1],minD[-1],maxD[-1],avgR[-1],medR[-1],stdR[-1],minR[-1],maxR[-1]))
        
            numpy.savez("PiMPI_%s_%s_%s_%s_%s_%s_%s_%s_%.8i_Device%i_%s_%s" % (ValueType,RNG,Alu[Devices[0]],GpuStyle,BlocksBegin,BlocksEnd,ThreadsBegin,ThreadsEnd,Iterations,Devices[0],Metrology,gethostname()),(ExploredBlocks,ExploredThreads,avgD,medD,stdD,minD,maxD,avgR,medR,stdR,minR,maxR))
            ToSave=[ ExploredBlocks,ExploredThreads,avgD,medD,stdD,minD,maxD,avgR,medR,stdR,minR,maxR ]
            numpy.savetxt("PiMPI_%s_%s_%s_%s_%s_%s_%s_%i_%.8i_Device%i_%s_%s" % (ValueType,RNG,Alu[Devices[0]],GpuStyle,BlocksBegin,BlocksEnd,ThreadsBegin,ThreadsEnd,Iterations,Devices[0],Metrology,gethostname()),numpy.transpose(ToSave),fmt='%i %i %e %e %e %e %e %i %i %i %i %i')

            if Fit:
                FitAndPrint(ExploredJobs,median,Curves)
        # Send MPI exit tag
        for slave in range(1,RankSize):
            comm.send('BREAK',dest=slave,tag=11)
            
    else:
        while True:
            Signal=comm.recv(source=0,tag=11)
            if Signal=='CONTINUE':
                # Receive information from Master
                InputCL=comm.recv(source=0,tag=11)
                print("Parameters retreive for rank %s of %s on %s from master:" % (rank,RankSize,gethostname()))
                print("Input CL:" % InputCL)
                # Execute on slave

                if InputCL['GpuStyle']=='CUDA':
                    try:
                        OutputCL=MetropolisCuda(InputCL)
                    except:
                        print("Problem with (%i,%i) // computations on Cuda" % (InputCL['Blocks'],InputCL['Threads']))
                elif InputCL['GpuStyle']=='OpenCL':
                    try:
                        OutputCL=MetropolisOpenCL(InputCL)
                    except:
                        print("Problem with (%i,%i) // computations on OpenCL" %  (InputCL['Blocks'],InputCL['Threads']))

                print("Output CL:" % OutputCL)
                # Send information to Master
                comm.send(OutputCL,dest=0,tag=11)
                print("Data sent to master")
            else:
                print('Exit signal from Master')
                break
