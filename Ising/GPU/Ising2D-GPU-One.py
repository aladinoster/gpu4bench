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
from numpy.random import randint as nprnd

KERNEL_CODE_OPENCL="""

// Marsaglia RNG very simple implementation
#define znew  ((z=36969*(z&65535)+(z>>16))<<16)
#define wnew  ((w=18000*(w&65535)+(w>>16))&65535)
#define MWC   (znew+wnew)
#define SHR3  (jsr=(jsr=(jsr=jsr^(jsr<<17))^(jsr>>13))^(jsr<<5))
#define CONG  (jcong=69069*jcong+1234567)
#define KISS  ((MWC^CONG)+SHR3)

#define MWCfp MWC * 2.328306435454494e-10f
#define KISSfp KISS * 2.328306435454494e-10f

__kernel void MainLoopOne(__global char *s,float T,float J,float B,
                          uint sizex,uint sizey,
                          uint iterations,uint seed_w,uint seed_z)

{
   uint z=seed_z;
   uint w=seed_w;

   for (uint i=0;i<iterations;i++) {

      uint x=(uint)(MWC%sizex) ;
      uint y=(uint)(MWC%sizey) ;

      int p=s[x+sizex*y];

      int d=s[x+sizex*((y+1)%sizey)];
      int u=s[x+sizex*((y-1)%sizey)];
      int l=s[((x-1)%sizex)+sizex*y];
      int r=s[((x+1)%sizex)+sizex*y];

      float DeltaE=2.0f*p*(J*(u+d+l+r)+B);

      int factor=((DeltaE < 0.0f) || (MWCfp < exp(-DeltaE/T))) ? -1:1;
      s[x%sizex+sizex*(y%sizey)] = (char)factor*p;
   }
   barrier(CLK_GLOBAL_MEM_FENCE);
   
}
"""

KERNEL_CODE_CUDA="""

// Marsaglia RNG very simple implementation
#define znew  ((z=36969*(z&65535)+(z>>16))<<16)
#define wnew  ((w=18000*(w&65535)+(w>>16))&65535)
#define MWC   (znew+wnew)
#define SHR3  (jsr=(jsr=(jsr=jsr^(jsr<<17))^(jsr>>13))^(jsr<<5))
#define CONG  (jcong=69069*jcong+1234567)
#define KISS  ((MWC^CONG)+SHR3)

#define MWCfp MWC * 2.328306435454494e-10f
#define KISSfp KISS * 2.328306435454494e-10f

__global__ void MainLoopOne(char *s,float T,float J,float B,
                            uint sizex,uint sizey,
                            uint iterations,uint seed_w,uint seed_z)

{
   uint z=seed_z;
   uint w=seed_w;

   for (uint i=0;i<iterations;i++) {

      uint x=(uint)(MWC%sizex) ;
      uint y=(uint)(MWC%sizey) ;

      int p=s[x+sizex*y];

      int d=s[x+sizex*((y+1)%sizey)];
      int u=s[x+sizex*((y-1)%sizey)];
      int l=s[((x-1)%sizex)+sizex*y];
      int r=s[((x+1)%sizex)+sizex*y];

      float DeltaE=2.0f*p*(J*(u+d+l+r)+B);

      int factor=((DeltaE < 0.0f) || (MWCfp < exp(-DeltaE/T))) ? -1:1;
      s[x%sizex+sizex*(y%sizey)] = (char)factor*p;
   }
   __syncthreads();
   
}
"""

def ImageOutput(sigma,prefix):
    Max=sigma.max()
    Min=sigma.min()
    
    # Normalize value as 8bits Integer
    SigmaInt=(255*(sigma-Min)/(Max-Min)).astype('uint8')
    image = Image.fromarray(SigmaInt)
    image.save("%s.jpg" % prefix)
    
def Metropolis(sigma,T,J,B,iterations): 
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

def MetropolisCuda(sigma,T,J,B,iterations,ParaStyle,Alu,Device):
    
    # Avec PyCUDA autoinit, rien a faire !

    sigmaCU=cuda.InOut(sigma)

    mod = SourceModule(KERNEL_CODE_CUDA)
    
    MetropolisCU=mod.get_function("MainLoopOne")

    start = pycuda.driver.Event()
    stop = pycuda.driver.Event()

    SizeX,SizeY=sigma.shape

    start.record()
    start.synchronize()
    MetropolisCU(sigmaCU,
                 numpy.float32(T), 
                 numpy.float32(J),
                 numpy.float32(B),
                 numpy.uint32(SizeX),
                 numpy.uint32(SizeY),
                 numpy.uint32(iterations),
                 numpy.uint32(nprnd(2**31-1)),
                 numpy.uint32(nprnd(2**31-1)),
                 grid=(1,1),
                 block=(1,1,1))
    
    print "%s with %i %s done" % (Alu,1,ParaStyle)

    stop.record()
    stop.synchronize()
                
    #elapsed = stop.time_since(start)*1e-3
    elapsed = start.time_till(stop)*1e-3

    return(elapsed)


def MetropolisOpenCL(sigma,T,J,B,iterations,ParaStyle,Alu,Device):

    # Initialisation des variables en les CASTant correctement
    
    # Je detecte un peripherique GPU dans la liste des peripheriques
    # for platform in cl.get_platforms():
    #     for device in platform.get_devices():
    #             if cl.device_type.to_string(device.type)=='GPU':
    #                     GPU=device
    #print "GPU detected: ",device.name
    
    HasGPU=False
    Id=1
    # Primary Device selection based on Device Id
    for platform in cl.get_platforms():
        for device in platform.get_devices():
            #deviceType=cl.device_type.to_string(device.type)
            deviceType="xPU"
            if Id==Device and not HasGPU:
                GPU=device
                print "CPU/GPU selected: ",device.name
                HasGPU=True
            Id=Id+1
            
    # Je cree le contexte et la queue pour son execution
    # ctx = cl.create_some_context()
    ctx = cl.Context([GPU])
    queue = cl.CommandQueue(ctx,
                            properties=cl.command_queue_properties.PROFILING_ENABLE)
    
    # Je recupere les flag possibles pour les buffers
    mf = cl.mem_flags

    # Attention au CAST ! C'est un int8 soit un char en OpenCL !
    sigmaCL = cl.Buffer(ctx, mf.WRITE_ONLY|mf.COPY_HOST_PTR,hostbuf=sigma)
   
    MetropolisCL = cl.Program(ctx,KERNEL_CODE_OPENCL).build( \
        options = "-cl-mad-enable -cl-fast-relaxed-math")

    SizeX,SizeY=sigma.shape

    if ParaStyle=='Blocks':
        # Call OpenCL kernel
        # (1,) is Global work size (only 1 work size)
        # (1,) is local work size
        CLLaunch=MetropolisCL.MainLoopOne(queue,(1,),None,
                                          sigmaCL,
                                          numpy.float32(T), 
                                          numpy.float32(J),
                                          numpy.float32(B),
                                          numpy.uint32(SizeX),
                                          numpy.uint32(SizeY),
                                          numpy.uint32(iterations),
                                          numpy.uint32(nprnd(2**31-1)),
                                          numpy.uint32(nprnd(2**31-1)))
        print "%s with %i %s done" % (Alu,1,ParaStyle)
    else:
        # en OpenCL, necessaire de mettre un Global_id identique au local_id
        CLLaunch=MetropolisCL.MainLoopOne(queue,(1,),(1,),
                                          sigmaCL,
                                          numpy.float32(T),
                                          numpy.float32(J),
                                          numpy.float32(B),
                                          numpy.uint32(SizeX),
                                          numpy.uint32(SizeY),
                                          numpy.uint32(iterations),
                                          numpy.uint32(nprnd(2**31-1)),
                                          numpy.uint32(nprnd(2**31-1)))
        print "%s with %i %s done" % (Alu,1,ParaStyle)
        
    CLLaunch.wait()
    cl.enqueue_copy(queue, sigma, sigmaCL).wait()
    elapsed = 1e-9*(CLLaunch.profile.end - CLLaunch.profile.start)
    sigmaCL.release()
        
    return(elapsed)

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
    # Alu can be CPU or GPU
    Alu='CPU'
    # Id of GPU : 0 will use the first find !
    Device=0
    # GPU style can be Cuda (Nvidia implementation) or OpenCL
    GpuStyle='OpenCL'
    # Parallel distribution can be on Threads or Blocks
    ParaStyle='Blocks'
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
        opts, args = getopt.getopt(sys.argv[1:],"hcj:b:z:i:s:e:p:a:d:g:t:",["coupling=","magneticfield=","size=","iterations=","tempstart=","tempend=","tempstep=","alu=","gpustyle=","parastyle="])
    except getopt.GetoptError:
        print '%s -j <Coupling Factor> -b <Magnetic Field> -z <Size of Lattice> -i <Iterations> -s <Minimum Temperature> -e <Maximum Temperature> -p <steP Temperature> -c (Print Curves) -a <CPU/GPU> -d <DeviceId> -g <CUDA/OpenCL> -p <Threads/Blocks> -t <ParaStyle>' % sys.argv[0]
        sys.exit(2)
    
 
    for opt, arg in opts:
        if opt == '-h':
            print '%s -j <Coupling Factor> -b <Magnetic Field> -z <Size of Lattice> -i <Iterations> -s <Minimum Temperature> -e <Maximum Temperature> -p <steP Temperature> -c (Print Curves) -a <CPU/GPU> -d <DeviceId> -g <CUDA/OpenCL> -p <Threads/Blocks> -t <ParaStyle>' % sys.argv[0]
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
        elif opt in ("-a", "--alu"):
            Alu = arg
        elif opt in ("-d", "--device"):
            Device = int(arg)
        elif opt in ("-g", "--gpustyle"):
            GpuStyle = arg
        elif opt in ("-t", "--parastyle"):
            ParaStyle = arg

    if Alu=='CPU' and GpuStyle=='CUDA':
        print "Alu can't be CPU for CUDA, set Alu to GPU"
        Alu='GPU'

    if ParaStyle not in ('Blocks','Threads','Hybrid'):
        print "%s not exists, ParaStyle set as Threads !" % ParaStyle
        ParaStyle='Blocks'
   
    print "Compute unit : %s" % Alu
    print "Device Identification : %s" % Device
    print "GpuStyle used : %s" % GpuStyle
    print "Parallel Style used : %s" % ParaStyle
    print "Coupling Factor : %s" % J
    print "Magnetic Field :  %s" % B
    print "Size of lattice : %s" % Size
    print "Iterations : %s" % Iterations
    print "Temperature on start : %s" % Tmin
    print "Temperature on end : %s" % Tmax
    print "Temperature step : %s" % Tstep

    if GpuStyle=='CUDA':
        # For PyCUDA import
        import pycuda.driver as cuda
        import pycuda.gpuarray as gpuarray
        import pycuda.autoinit
        from pycuda.compiler import SourceModule

    if GpuStyle=='OpenCL':
        # For PyOpenCL import
        import pyopencl as cl
        Id=1
        for platform in cl.get_platforms():
            for device in platform.get_devices():
                #deviceType=cl.device_type.to_string(device.type)
                deviceType="xPU"
                print "Device #%i of type %s : %s" % (Id,deviceType,device.name)
                Id=Id+1

    LAPIMAGE=False

    sigmaIn=numpy.where(numpy.random.randn(Size,Size)>0,1,-1).astype(numpy.int8)

    ImageOutput(sigmaIn,"Ising2D_Serial_%i_Initial" % (Size))

    Trange=numpy.arange(Tmin,Tmax+Tstep,Tstep)

    E=[]
    M=[]

    for T in Trange:
        sigma=numpy.copy(sigmaIn)
        if GpuStyle=='CUDA':
            duration=MetropolisCuda(sigma,T,J,B,Iterations,ParaStyle,Alu,Device)
        else:
            duration=MetropolisOpenCL(sigma,T,J,B,Iterations,ParaStyle,Alu,Device)
            
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
    
