#!/usr/bin/env python

#
# Splutter-by-MonteCarlo using PyCUDA/PyOpenCL
#
# CC BY-NC-SA 2014 : <emmanuel.quemener@ens-lyon.fr> 
#
# Thanks to Andreas Klockner for PyCUDA:
# http://mathema.tician.de/software/pycuda
# http://mathema.tician.de/software/pyopencl
# 

# 2013-01-01 : problems with launch timeout
# http://stackoverflow.com/questions/497685/how-do-you-get-around-the-maximum-cuda-run-time
# Option "Interactive" "0" in /etc/X11/xorg.conf

# Marsaglia elements about RNG 

# Common tools
import numpy
from numpy.random import randint as nprnd
import sys
import getopt
import time
import math
from socket import gethostname

# find prime factors of a number
# Get for WWW :
# http://pythonism.wordpress.com/2008/05/17/looking-at-factorisation-in-python/
def PrimeFactors(x):
  factorlist=numpy.array([]).astype('uint32')
  loop=2
  while loop<=x:
    if x%loop==0:
      x/=loop
      factorlist=numpy.append(factorlist,[loop])
    else:
      loop+=1
  return factorlist
  
# Try to find the best thread number in Hybrid approach (Blocks&Threads)
# output is thread number
def BestThreadsNumber(jobs):
  factors=PrimeFactors(jobs)
  matrix=numpy.append([factors],[factors[::-1]],axis=0)
  threads=1
  for factor in matrix.transpose().ravel():
    threads=threads*factor
    if threads*threads>jobs:
      break
  return(long(threads))

# Predicted Amdahl Law (Reduced with s=1-p)  
def AmdahlR(N, T1, p):
  return (T1*(1-p+p/N))

# Predicted Amdahl Law
def Amdahl(N, T1, s, p):
  return (T1*(s+p/N))

# Predicted Mylq Law with first order
def Mylq(N, T1,s,c,p):
  return (T1*(s+p/N)+c*N)

# Predicted Mylq Law with second order
def Mylq2(N, T1,s,c1,c2,p):
  return (T1*(s+p/N)+c1*N+c2*N*N)

KERNEL_CODE_CUDA="""

// Marsaglia RNG very simple implementation

#define znew  ((z=36969*(z&65535)+(z>>16))<<16)
#define wnew  ((w=18000*(w&65535)+(w>>16))&65535)

#define MWC   (znew+wnew)
#define SHR3  (jsr=(jsr=(jsr=jsr^(jsr<<17))^(jsr>>13))^(jsr<<5))
#define CONG  (jcong=69069*jcong+1234567)
#define KISS  ((MWC^CONG)+SHR3)

#define CONGfp CONG * 2.328306435454494e-10f
#define SHR3fp SHR3 * 2.328306435454494e-10f
#define MWCfp MWC * 2.328306435454494e-10f
#define KISSfp KISS * 2.328306435454494e-10f

#define MAX (ulong)4294967296
#define UMAX (uint)2147483648

__global__ void SplutterGlobal(uint *s,const uint space,const ulong iterations,const uint seed_w,const uint seed_z)
{
    const ulong id=(ulong)(blockIdx.x);
   
    uint z=seed_z-(uint)id;
    uint w=seed_w+(uint)id;

    uint jsr=seed_z;
    uint jcong=seed_w;

   for ( ulong i=0;i<iterations;i++) {

      // All version
      uint position=(uint)( ((ulong)MWC*(ulong)space)/MAX );

      // UMAX is set to avoid round over overflow
      atomicInc(&s[position],UMAX);
   }

   __syncthreads();
}

__global__ void SplutterGlobalDense(uint *s,const uint space,const ulong iterations,const uint seed_w,const uint seed_z)
{
    const ulong id=(ulong)(threadIdx.x+blockIdx.x*blockDim.x);
    const ulong size=(ulong)(gridDim.x*blockDim.x);
    const ulong block=(ulong)space/(ulong)size;
   
    uint z=seed_z-(uint)id;
    uint w=seed_w+(uint)id;

    uint jsr=seed_z;
    uint jcong=seed_w;

   for ( ulong i=0;i<iterations;i++) {

      // Dense version
       uint position=(uint)( ((ulong)MWC+id*MAX)*block/MAX );

      s[position]++;
   }

   __syncthreads();
}

__global__ void SplutterGlobalSparse(uint *s,const uint space,const ulong iterations,const uint seed_w,const uint seed_z)
{ 
    const ulong id=(ulong)(threadIdx.x+blockIdx.x*blockDim.x);
    const ulong size=(ulong)(gridDim.x*blockDim.x);
    const ulong block=(ulong)space/(ulong)size;
   
    uint z=seed_z-(uint)id;
    uint w=seed_w+(uint)id;

    uint jsr=seed_z;
    uint jcong=seed_w;

   for ( ulong i=0;i<iterations;i++) {

      // Sparse version
       uint position=(uint)( (ulong)MWC*block/MAX*size+id );

      s[position]++;
   }

   __syncthreads();
}

__global__ void SplutterLocalDense(uint *s,const uint space,const ulong iterations,const uint seed_w,const uint seed_z)
{
    const ulong id=(ulong)(threadIdx.x);
    const ulong size=(ulong)(blockDim.x);
    const ulong block=(ulong)space/(ulong)size;
   
    uint z=seed_z-(uint)id;
    uint w=seed_w+(uint)id;

    uint jsr=seed_z;
    uint jcong=seed_w;

   for ( ulong i=0;i<iterations;i++) {

      // Dense version
       size_t position=(size_t)( ((ulong)MWC+id*MAX)*block/MAX );

      s[position]++;
   }


   __syncthreads();

}

__global__ void SplutterLocalSparse(uint *s,const uint space,const ulong iterations,const uint seed_w,const uint seed_z)
{
    const ulong id=(ulong)threadIdx.x;
    const ulong size=(ulong)blockDim.x;
    const ulong block=(ulong)space/(ulong)size;
   
    uint z=seed_z-(uint)id;
    uint w=seed_w+(uint)id;

    uint jsr=seed_z;
    uint jcong=seed_w;

   for ( ulong i=0;i<iterations;i++) {

      // Sparse version
       size_t position=(size_t)( (ulong)MWC*block/MAX*size+id );

      s[position]++;
   }

   __syncthreads();

}

__global__ void SplutterHybridDense(uint *s,const uint space,const ulong iterations,const uint seed_w,const uint seed_z)
{
    const ulong id=(ulong)(blockIdx.x);
    const ulong size=(ulong)(gridDim.x);
    const ulong block=(ulong)space/(ulong)size;
   
    uint z=seed_z-(uint)id;
    uint w=seed_w+(uint)id;

    uint jsr=seed_z;
    uint jcong=seed_w;

   for ( ulong i=0;i<iterations;i++) {

      // Dense version
      size_t position=(size_t)( ((ulong)MWC+id*MAX)*block/MAX );

      s[position]++;
   }
      
   __syncthreads();
}

__global__ void SplutterHybridSparse(uint *s,const uint space,const ulong iterations,const uint seed_w,const uint seed_z)
{
    const ulong id=(ulong)(blockIdx.x);
    const ulong size=(ulong)(gridDim.x);
    const ulong block=(ulong)space/(ulong)size;
   
    uint z=seed_z-(uint)id;
    uint w=seed_w+(uint)id;

    uint jsr=seed_z;
    uint jcong=seed_w;

   for ( ulong i=0;i<iterations;i++) {

      // Sparse version
      size_t position=(size_t)( (((ulong)MWC*block)/MAX)*size+id );

      s[position]++;

   }

   //s[blockIdx.x]=blockIdx.x;
   __syncthreads();
}

"""

KERNEL_CODE_OPENCL="""
#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable

// Marsaglia RNG very simple implementation
#define znew  ((z=36969*(z&65535)+(z>>16))<<16)
#define wnew  ((w=18000*(w&65535)+(w>>16))&65535)

#define MWC   (znew+wnew)
#define SHR3  (jsr=(jsr=(jsr=jsr^(jsr<<17))^(jsr>>13))^(jsr<<5))
#define CONG  (jcong=69069*jcong+1234567)
#define KISS  ((MWC^CONG)+SHR3)

#define CONGfp CONG * 2.328306435454494e-10f
#define SHR3fp SHR3 * 2.328306435454494e-10f
#define MWCfp MWC * 2.328306435454494e-10f
#define KISSfp KISS * 2.328306435454494e-10f

#define MAX (ulong)4294967296

uint rotl(uint value, int shift) {
    return (value << shift) | (value >> (sizeof(value) * CHAR_BIT - shift));
}
 
uint rotr(uint value, int shift) {
    return (value >> shift) | (value << (sizeof(value) * CHAR_BIT - shift));
}

__kernel void SplutterGlobal(__global uint *s,const uint space,const ulong iterations,const uint seed_w,const uint seed_z)
{
   __private const ulong id=(ulong)get_global_id(0);
   
   __private uint z=seed_z-(uint)id;
   __private uint w=seed_w+(uint)id;

   __private uint jsr=seed_z;
   __private uint jcong=seed_w;

   for (__private ulong i=0;i<iterations;i++) {

      // Dense version
      __private size_t position=(size_t)( MWC%space );

      atomic_inc(&s[position]);
   }

   barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);

}

__kernel void SplutterLocal(__global uint *s,const uint space,const ulong iterations,const uint seed_w,const uint seed_z)
{
   __private const ulong id=(ulong)get_local_id(0);
   
   __private uint z=seed_z-(uint)id;
   __private uint w=seed_w+(uint)id;

   __private uint jsr=seed_z;
   __private uint jcong=seed_w;

   for (__private ulong i=0;i<iterations;i++) {

      // Dense version
      //__private size_t position=(size_t)( (MWC+id*block)%space );
      __private size_t position=(size_t)( MWC%space );

      atomic_inc(&s[position]);
   }

   barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);

}

__kernel void SplutterHybrid(__global uint *s,const uint space,const ulong iterations,const uint seed_w,const uint seed_z)
{
   __private const ulong id=(ulong)(get_global_id(0)+get_local_id(0));
   
   __private uint z=seed_z-(uint)id;
   __private uint w=seed_w+(uint)id;

   __private uint jsr=seed_z;
   __private uint jcong=seed_w;

   for (__private ulong i=0;i<iterations;i++) {

      // Dense version
      __private size_t position=(size_t)( MWC%space );

      atomic_inc(&s[position]);
   }
      
}

"""

def MetropolisCuda(circle,iterations,steps,jobs,ParaStyle,Density,Memory):

  # Avec PyCUDA autoinit, rien a faire !

  circleCU = cuda.InOut(circle)
  
  mod = SourceModule(KERNEL_CODE_CUDA)

  if Density=='Dense':
    MetropolisBlocksCU=mod.get_function("SplutterGlobalDense")
    MetropolisThreadsCU=mod.get_function("SplutterLocalDense")
    MetropolisHybridCU=mod.get_function("SplutterHybridDense")
  elif Density=='Sparse':
    MetropolisBlocksCU=mod.get_function("SplutterGlobalSparse")
    MetropolisThreadsCU=mod.get_function("SplutterLocalSparse")
    MetropolisHybridCU=mod.get_function("SplutterHybridSparse")
  else:
    MetropolisBlocksCU=mod.get_function("SplutterGlobal")
    
  start = pycuda.driver.Event()
  stop = pycuda.driver.Event()
  
  MySplutter=numpy.zeros(steps)
  MyDuration=numpy.zeros(steps)

  if iterations%jobs==0:
    iterationsCL=numpy.uint64(iterations/jobs)
  else:
    iterationsCL=numpy.uint64(iterations/jobs+1)
    
  iterationsNew=iterationsCL*jobs

  Splutter=numpy.zeros(jobs*16).astype(numpy.uint32)

  for i in range(steps):

    start_time=time.time()
    Splutter[:]=0
    
    print Splutter,len(Splutter)

    SplutterCU = cuda.InOut(Splutter)

    start.record()
    start.synchronize()
    if ParaStyle=='Blocks':
      MetropolisBlocksCU(SplutterCU,
                         numpy.uint32(len(Splutter)),
                         numpy.uint64(iterationsCL),
                         numpy.uint32(nprnd(2**30/jobs)),
                         numpy.uint32(nprnd(2**30/jobs)),
                         grid=(jobs,1),
                         block=(1,1,1))
        
      print "%s with (WorkItems/Threads)=(%i,%i) %s method done" % \
            (Alu,jobs,1,ParaStyle)      
    elif ParaStyle=='Hybrid':
      threads=BestThreadsNumber(jobs)
      MetropolisHybridCU(SplutterCU,
                         numpy.uint32(len(Splutter)),
                         numpy.uint64(iterationsCL),
                         numpy.uint32(nprnd(2**30/jobs)),
                         numpy.uint32(nprnd(2**30/jobs)),
                         grid=(jobs,1),
                         block=(threads,1,1))
      print "%s with (WorkItems/Threads)=(%i,%i) %s method done" % \
            (Alu,jobs/threads,threads,ParaStyle)
    else:
      MetropolisThreadsCU(SplutterCU,
                       numpy.uint32(len(Splutter)),
                       numpy.uint64(iterationsCL),
                       numpy.uint32(nprnd(2**30/jobs)),
                       numpy.uint32(nprnd(2**30/jobs)),
                       grid=(1,1),
                       block=(jobs,1,1))
      print "%s with (WorkItems/Threads)=(%i,%i) %s method done" % \
            (Alu,1,jobs,ParaStyle)
    stop.record()
    stop.synchronize()
                
#    elapsed = start.time_till(stop)*1e-3
    elapsed = time.time()-start_time

    print Splutter,sum(Splutter)
    MySplutter[i]=numpy.median(Splutter)
    print numpy.mean(Splutter),MySplutter[i],numpy.std(Splutter)

    MyDuration[i]=elapsed

    #AllPi=4./numpy.float32(iterationsCL)*circle.astype(numpy.float32)
    #MyPi[i]=numpy.median(AllPi)
    #print MyPi[i],numpy.std(AllPi),MyDuration[i]


  print jobs,numpy.mean(MyDuration),numpy.median(MyDuration),numpy.std(MyDuration)

  return(numpy.mean(MyDuration),numpy.median(MyDuration),numpy.std(MyDuration))


def MetropolisOpenCL(circle,iterations,steps,jobs,
                     ParaStyle,Alu,Device,Memory):
	
  # Initialisation des variables en les CASTant correctement

  MaxMemoryXPU=0
  MinMemoryXPU=0

  if Device==0:
    print "Enter XPU selector based on ALU type: first selected"
    HasXPU=False
    # Default Device selection based on ALU Type
    for platform in cl.get_platforms():
      for device in platform.get_devices():
        #deviceType=cl.device_type.to_string(device.type)
        deviceMemory=device.max_mem_alloc_size
        if deviceMemory>MaxMemoryXPU:
          MaxMemoryXPU=deviceMemory
        if deviceMemory<MinMemoryXPU or MinMemoryXPU==0:
          MinMemoryXPU=deviceMemory
        if not HasXPU:        
          XPU=device
          print "XPU selected with Allocable Memory %i: %s" % (deviceMemory,device.name)
          HasXPU=True
          MemoryXPU=deviceMemory
          
  else:
    print "Enter XPU selector based on device number & ALU type"
    Id=1
    HasXPU=False
    # Primary Device selection based on Device Id
    for platform in cl.get_platforms():
      for device in platform.get_devices():
        #deviceType=cl.device_type.to_string(device.type)
        deviceMemory=device.max_mem_alloc_size
        if deviceMemory>MaxMemoryXPU:
          MaxMemoryXPU=deviceMemory
        if deviceMemory<MinMemoryXPU or MinMemoryXPU==0:
          MinMemoryXPU=deviceMemory
        if Id==Device  and HasXPU==False:
          XPU=device
          print "CPU/GPU selected with Allocable Memory %i: %s" % (deviceMemory,device.name)
          HasXPU=True
          MemoryXPU=deviceMemory
        Id=Id+1
    if HasXPU==False:
      print "No XPU #%i of type %s found in all of %i devices, sorry..." % \
          (Device,Alu,Id-1)
      return(0,0,0)

  print "Allocable Memory is %i, between %i and %i " % (MemoryXPU,MinMemoryXPU,MaxMemoryXPU)

  # Je cree le contexte et la queue pour son execution
  ctx = cl.Context([XPU])
  queue = cl.CommandQueue(ctx,properties=cl.command_queue_properties.PROFILING_ENABLE)
  
  # Je recupere les flag possibles pour les buffers
  mf = cl.mem_flags

  MetropolisCL = cl.Program(ctx,KERNEL_CODE_OPENCL).build(options = "-cl-mad-enable -cl-fast-relaxed-math")
      
  MyDuration=numpy.zeros(steps)
  
  if iterations%jobs==0:
    iterationsCL=numpy.uint64(iterations/jobs)
  else:
    iterationsCL=numpy.uint64(iterations/jobs+1)
    
  iterationsNew=numpy.uint64(iterationsCL*jobs)

  MySplutter=numpy.zeros(steps)

  MaxWorks=2**(int)(numpy.log2(MinMemoryXPU/4))
  print MaxWorks,2**(int)(numpy.log2(MemoryXPU))
  
  #Splutter=numpy.zeros((MaxWorks/jobs)*jobs).astype(numpy.uint32)
  #Splutter=numpy.zeros(jobs*16).astype(numpy.uint32)
  Splutter=numpy.zeros(Memory).astype(numpy.uint32)

  for i in range(steps):
		
    #Splutter=numpy.zeros(2**(int)(numpy.log2(MemoryXPU/4))).astype(numpy.uint32)
    #Splutter=numpy.zeros(1024).astype(numpy.uint32)
 
    #Splutter=numpy.zeros(jobs).astype(numpy.uint32)

    Splutter[:]=0

    print Splutter,len(Splutter)

    h2d_time=time.time()
    SplutterCL = cl.Buffer(ctx, mf.WRITE_ONLY|mf.COPY_HOST_PTR,hostbuf=Splutter)
    print('From Host to Device time %f' % (time.time()-h2d_time))

    start_time=time.time()
    if ParaStyle=='Blocks':
      # Call OpenCL kernel
      # (1,) is Global work size (only 1 work size)
      # (1,) is local work size
      # circleCL is lattice translated in CL format
      # SeedZCL is lattice translated in CL format
      # SeedWCL is lattice translated in CL format
      # step is number of iterations
      # CLLaunch=MetropolisCL.MainLoopGlobal(queue,(jobs,),None,
      #                                      SplutterCL,
      #                                      numpy.uint32(len(Splutter)),
      #                                      numpy.uint64(iterationsCL),
      #                                      numpy.uint32(nprnd(2**30/jobs)),
      #                                      numpy.uint32(nprnd(2**30/jobs)))
      CLLaunch=MetropolisCL.SplutterGlobal(queue,(jobs,),None,
                                           SplutterCL,
                                           numpy.uint32(len(Splutter)),
                                           numpy.uint64(iterationsCL),
                                           numpy.uint32(nprnd(2**30/jobs)),
                                           numpy.uint32(nprnd(2**30/jobs)))
        
      print "%s with (WorkItems/Threads)=(%i,%i) %s method done" % \
            (Alu,jobs,1,ParaStyle)
    elif ParaStyle=='Hybrid':
      #threads=BestThreadsNumber(jobs)
      threads=BestThreadsNumber(256)
      print "print",threads      
      # en OpenCL, necessaire de mettre un Global_id identique au local_id
      CLLaunch=MetropolisCL.SplutterHybrid(queue,(jobs,),(threads,),
                                           SplutterCL,
                                           numpy.uint32(len(Splutter)),
                                           numpy.uint64(iterationsCL),
                                           numpy.uint32(nprnd(2**30/jobs)),
                                           numpy.uint32(nprnd(2**30/jobs)))
        
      print "%s with (WorkItems/Threads)=(%i,%i) %s method done" % \
            (Alu,jobs/threads,threads,ParaStyle)
    else:
      # en OpenCL, necessaire de mettre un global_id identique au local_id
      CLLaunch=MetropolisCL.SplutterLocal(queue,(jobs,),(jobs,),
                                          SplutterCL,
                                          numpy.uint32(len(Splutter)),
                                          numpy.uint64(iterationsCL),
                                          numpy.uint32(nprnd(2**30/jobs)),
                                          numpy.uint32(nprnd(2**30/jobs)))
        
        
      print "%s with %i %s done" % (Alu,jobs,ParaStyle)

    CLLaunch.wait()
    d2h_time=time.time()
    cl.enqueue_copy(queue, Splutter, SplutterCL).wait()
    print('From Device to Host %f' % (time.time()-d2h_time))
    
#    elapsed = 1e-9*(CLLaunch.profile.end - CLLaunch.profile.start)
    elapsed = time.time()-start_time
    print('Elapsed compute time %f' % elapsed)

    MyDuration[i]=elapsed
    #print Splutter,sum(Splutter)
    #MySplutter[i]=numpy.median(Splutter)
    #print numpy.mean(Splutter)*len(Splutter),MySplutter[i]*len(Splutter),numpy.std(Splutter)
    
  SplutterCL.release()

  print jobs,numpy.mean(MyDuration),numpy.median(MyDuration),numpy.std(MyDuration)
	
  return(numpy.mean(MyDuration),numpy.median(MyDuration),numpy.std(MyDuration))


def FitAndPrint(N,D,Curves):

  from scipy.optimize import curve_fit
  import matplotlib.pyplot as plt

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
    # coeffs_Mylq[2]=coeffs_Mylq[2]*coeffs_Mylq[0]/D[0]
    coeffs_Mylq[3]=coeffs_Mylq[3]*coeffs_Mylq[0]/D[0]
    coeffs_Mylq[0]=D[0]
    print "Mylq Normalized : T=%.2f(%.6f+%.6f/N)+%.6f*N" % (coeffs_Mylq[0],
                                                            coeffs_Mylq[1],
                                                            coeffs_Mylq[3],
                                                            coeffs_Mylq[2])
    D_Mylq=Mylq(N,coeffs_Mylq[0],coeffs_Mylq[1],coeffs_Mylq[2],
                coeffs_Mylq[3])
  except:
    print "Impossible to fit for Mylq law : only %i elements" % len(D) 

  try:
    coeffs_Mylq2, matcov_Mylq2 = curve_fit(Mylq2, N, D)

    coeffs_Mylq2[1]=coeffs_Mylq2[1]*coeffs_Mylq2[0]/D[0]
    # coeffs_Mylq2[2]=coeffs_Mylq2[2]*coeffs_Mylq2[0]/D[0]
    # coeffs_Mylq2[3]=coeffs_Mylq2[3]*coeffs_Mylq2[0]/D[0]
    coeffs_Mylq2[4]=coeffs_Mylq2[4]*coeffs_Mylq2[0]/D[0]
    coeffs_Mylq2[0]=D[0]
    print "Mylq 2nd order Normalized: T=%.2f(%.6f+%.6f/N)+%.6f*N+%.6f*N^2" % \
        (coeffs_Mylq2[0],coeffs_Mylq2[1],
         coeffs_Mylq2[4],coeffs_Mylq2[2],coeffs_Mylq2[3])

  except:
    print "Impossible to fit for 2nd order Mylq law : only %i elements" % len(D) 

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
  
  # Alu can be CPU, GPU or ACCELERATOR
  Alu='CPU'
  # Id of GPU : 1 is for first find !
  Device=0
  # GPU style can be Cuda (Nvidia implementation) or OpenCL
  GpuStyle='OpenCL'
  # Parallel distribution can be on Threads or Blocks
  ParaStyle='Blocks'
  # Iterations is integer
  Iterations=10000000
  # JobStart in first number of Jobs to explore
  JobStart=1
  # JobEnd is last number of Jobs to explore
  JobEnd=16
  # JobStep is the step of Jobs to explore
  JobStep=1
  # Redo is the times to redo the test to improve metrology
  Redo=1
  # OutMetrology is method for duration estimation : False is GPU inside
  OutMetrology=False
  Metrology='InMetro'
  # Curves is True to print the curves
  Curves=False
  # Fit is True to print the curves
  Fit=False
  # Memory of vector explored
  Memory=1024

  try:
    opts, args = getopt.getopt(sys.argv[1:],"hocfa:g:p:i:s:e:t:r:d:m:",["alu=","gpustyle=","parastyle=","iterations=","jobstart=","jobend=","jobstep=","redo=","device="])
  except getopt.GetoptError:
    print '%s -o (Out of Core Metrology) -c (Print Curves) -f (Fit to Amdahl Law) -a <CPU/GPU/ACCELERATOR> -d <DeviceId> -g <CUDA/OpenCL> -p <Threads/Hybrid/Blocks> -i <Iterations> -s <JobStart> -e <JobEnd> -t <JobStep> -r <RedoToImproveStats> -m <MemoryRaw>' % sys.argv[0]
    sys.exit(2)
    
  for opt, arg in opts:
    if opt == '-h':
      print '%s -o (Out of Core Metrology) -c (Print Curves) -f (Fit to Amdahl Law) -a <CPU/GPU/ACCELERATOR> -d <DeviceId> -g <CUDA/OpenCL> -p <Threads/Hybrid/Blocks> -i <Iterations> -s <JobStart> -e <JobEnd> -t <JobStep> -r <RedoToImproveStats> -m <MemoryRaw>' % sys.argv[0]

      print "\nInformations about devices detected under OpenCL:"
      # For PyOpenCL import
      try:
        import pyopencl as cl
        Id=1
        for platform in cl.get_platforms():
          for device in platform.get_devices():
            #deviceType=cl.device_type.to_string(device.type)
            deviceMemory=device.max_mem_alloc_size
            print "Device #%i from %s with memory %i : %s" % (Id,platform.vendor,deviceMemory,device.name.lstrip())
            Id=Id+1

        print
        sys.exit()
      except ImportError:
        print "Your platform does not seem to support OpenCL"
        
    elif opt == '-o':
      OutMetrology=True
      Metrology='OutMetro'
    elif opt == '-c':
      Curves=True
    elif opt == '-f':
      Fit=True
    elif opt in ("-a", "--alu"):
      Alu = arg
    elif opt in ("-d", "--device"):
      Device = int(arg)
    elif opt in ("-g", "--gpustyle"):
      GpuStyle = arg
    elif opt in ("-p", "--parastyle"):
      ParaStyle = arg
    elif opt in ("-i", "--iterations"):
      Iterations = numpy.uint64(arg)
    elif opt in ("-s", "--jobstart"):
      JobStart = int(arg)
    elif opt in ("-e", "--jobend"):
      JobEnd = int(arg)
    elif opt in ("-t", "--jobstep"):
      JobStep = int(arg)
    elif opt in ("-r", "--redo"):
      Redo = int(arg)
    elif opt in ("-m", "--memory"):
      Memory = int(arg)

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
  print "Iterations : %s" % Iterations
  print "Number of threads on start : %s" % JobStart
  print "Number of threads on end : %s" % JobEnd
  print "Number of redo : %s" % Redo
  print "Memory  : %s" % Memory
  print "Metrology done out of CPU/GPU : %r" % OutMetrology

  if GpuStyle=='CUDA':
    try:
      # For PyCUDA import
      import pycuda.driver as cuda
      import pycuda.gpuarray as gpuarray
      import pycuda.autoinit
      from pycuda.compiler import SourceModule
    except ImportError:
      print "Platform does not seem to support CUDA"

  if GpuStyle=='OpenCL':
    try:
      # For PyOpenCL import
      import pyopencl as cl
      Id=1
      for platform in cl.get_platforms():
        for device in platform.get_devices():
          #deviceType=cl.device_type.to_string(device.type)
          print "Device #%i : %s" % (Id,device.name)
          if Id == Device:
            # Set the Alu as detected Device Type
            Alu='xPU'
          Id=Id+1
    except ImportError:
      print "Platform does not seem to support CUDA"
      
  average=numpy.array([]).astype(numpy.float32)
  median=numpy.array([]).astype(numpy.float32)
  stddev=numpy.array([]).astype(numpy.float32)

  ExploredJobs=numpy.array([]).astype(numpy.uint32)

  Jobs=JobStart

  while Jobs <= JobEnd:
    avg,med,std=0,0,0
    ExploredJobs=numpy.append(ExploredJobs,Jobs)
    circle=numpy.zeros(Jobs).astype(numpy.uint64)

    if OutMetrology:
      duration=numpy.array([]).astype(numpy.float32)
      for i in range(Redo):
        start=time.time()
        if GpuStyle=='CUDA':
          try:
            a,m,s=MetropolisCuda(circle,Iterations,1,Jobs,ParaStyle,
                                 Memory)
          except:
            print "Problem with %i // computations on Cuda" % Jobs
        elif GpuStyle=='OpenCL':
          try:
            a,m,s=MetropolisOpenCL(circle,Iterations,1,Jobs,ParaStyle,
                                   Alu,Device,Memory)
          except:
            print "Problem with %i // computations on OpenCL" % Jobs            
        duration=numpy.append(duration,time.time()-start)
      if (a,m,s) != (0,0,0):
        avg=numpy.mean(duration)
        med=numpy.median(duration)
        std=numpy.std(duration)
      else:
        print "Values seem to be wrong..."
    else:
      if GpuStyle=='CUDA':
        try:
          avg,med,std=MetropolisCuda(circle,Iterations,Redo,
                                     Jobs,ParaStyle,Memory)
        except:
          print "Problem with %i // computations on Cuda" % Jobs
      elif GpuStyle=='OpenCL':
        try:
          avg,med,std=MetropolisOpenCL(circle,Iterations,Redo,Jobs,
                                       ParaStyle,Alu,Device,Memory)
        except:
          print "Problem with %i // computations on OpenCL" % Jobs            

    if (avg,med,std) != (0,0,0):
      print "jobs,avg,med,std",Jobs,avg,med,std
      average=numpy.append(average,avg)
      median=numpy.append(median,med)
      stddev=numpy.append(stddev,std)
    else:
      print "Values seem to be wrong..."
    #THREADS*=2
    if len(average)!=0:
      numpy.savez("Splutter_%s_%s_%s_%i_%i_%.8i_Device%i_%s_%s" % (Alu,GpuStyle,ParaStyle,JobStart,JobEnd,Iterations,Device,Metrology,gethostname()),(ExploredJobs,average,median,stddev))
      ToSave=[ ExploredJobs,average,median,stddev ]
      numpy.savetxt("Splutter_%s_%s_%s_%i_%i_%.8i_Device%i_%s_%s" % (Alu,GpuStyle,ParaStyle,JobStart,JobEnd,Iterations,Device,Metrology,gethostname()),numpy.transpose(ToSave))
    Jobs+=JobStep

  if Fit:
    FitAndPrint(ExploredJobs,median,Curves)

