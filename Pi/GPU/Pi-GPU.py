#!/usr/bin/env python

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
from socket import gethostname

Marsaglia={'CONG':0,'SHR3':1,'MWC':2,'KISS':3}
Computing={'INT32':0,'INT64':1,'FP32':2,'FP64':3}

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
    if threads*threads>jobs or threads>512:
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
#define TCONG 0
#define TSHR3 1
#define TMWC 2
#define TKISS 3

#define TINT32 0
#define TINT64 1
#define TFP32 2
#define TFP64 3

// Marsaglia RNG very simple implementation

#define znew  ((z=36969*(z&65535)+(z>>16))<<16)
#define wnew  ((w=18000*(w&65535)+(w>>16))&65535)
#define MWC   (znew+wnew)
#define SHR3  (jsr=(jsr=(jsr=jsr^(jsr<<17))^(jsr>>13))^(jsr<<5))
#define CONG  (jcong=69069*jcong+1234567)
#define KISS  ((MWC^CONG)+SHR3)

#define MWCfp MWC * 2.328306435454494e-10f
#define KISSfp KISS * 2.328306435454494e-10f
#define SHR3fp SHR3 * 2.328306435454494e-10f
#define CONGfp CONG * 2.328306435454494e-10f

__device__ ulong MainLoop(ulong iterations,uint seed_w,uint seed_z,size_t work)
{

#if TRNG == TCONG
   uint jcong=seed_z+work;
#elif TRNG == TSHR3
   uint jsr=seed_w+work;
#elif TRNG == TMWC
   uint z=seed_z+work;
   uint w=seed_w+work;
#elif TRNG == TKISS
   uint jcong=seed_z+work;
   uint jsr=seed_w+work;
   uint z=seed_z-work;
   uint w=seed_w-work;
#endif

   ulong total=0;

   for (ulong i=0;i<iterations;i++) {

#if TYPE == TINT32
    #define THEONE 1073741824
    #if TRNG == TCONG
        uint x=CONG>>17 ;
        uint y=CONG>>17 ;
    #elif TRNG == TSHR3
        uint x=SHR3>>17 ;
        uint y=SHR3>>17 ;
    #elif TRNG == TMWC
        uint x=MWC>>17 ;
        uint y=MWC>>17 ;
    #elif TRNG == TKISS
        uint x=KISS>>17 ;
        uint y=KISS>>17 ;
    #endif
#elif TYPE == TINT64
    #define THEONE 4611686018427387904
    #if TRNG == TCONG
        ulong x=(ulong)(CONG>>1) ;
        ulong y=(ulong)(CONG>>1) ;
    #elif TRNG == TSHR3
        ulong x=(ulong)(SHR3>>1) ;
        ulong y=(ulong)(SHR3>>1) ;
    #elif TRNG == TMWC
        ulong x=(ulong)(MWC>>1) ;
        ulong y=(ulong)(MWC>>1) ;
    #elif TRNG == TKISS
        ulong x=(ulong)(KISS>>1) ;
        ulong y=(ulong)(KISS>>1) ;
    #endif
#elif TYPE == TFP32
    #define THEONE 1.0f
    #if TRNG == TCONG
        float x=CONGfp ;
        float y=CONGfp ;
    #elif TRNG == TSHR3
        float x=SHR3fp ;
        float y=SHR3fp ;
    #elif TRNG == TMWC
        float x=MWCfp ;
        float y=MWCfp ;
    #elif TRNG == TKISS
      float x=KISSfp ;
      float y=KISSfp ;
    #endif
#elif TYPE == TFP64
    #define THEONE 1.0f
    #if TRNG == TCONG
        double x=(double)CONGfp ;
        double y=(double)CONGfp ;
    #elif TRNG == TSHR3
        double x=(double)SHR3fp ;
        double y=(double)SHR3fp ;
    #elif TRNG == TMWC
        double x=(double)MWCfp ;
        double y=(double)MWCfp ;
    #elif TRNG == TKISS
        double x=(double)KISSfp ;
        double y=(double)KISSfp ;
    #endif
#endif

      ulong inside=((x*x+y*y) <= THEONE) ? 1:0;
      total+=inside;
   }

   return(total);
}

__global__ void MainLoopBlocks(ulong *s,ulong iterations,uint seed_w,uint seed_z)
{
   ulong total=MainLoop(iterations,seed_z,seed_w,blockIdx.x);
   s[blockIdx.x]=total;
   __syncthreads();

}

__global__ void MainLoopThreads(ulong *s,ulong iterations,uint seed_w,uint seed_z)
{
   ulong total=MainLoop(iterations,seed_z,seed_w,threadIdx.x);
   s[threadIdx.x]=total;
   __syncthreads();

}

__global__ void MainLoopHybrid(ulong *s,ulong iterations,uint seed_w,uint seed_z)
{
   ulong total=MainLoop(iterations,seed_z,seed_w,blockDim.x*blockIdx.x+threadIdx.x);
   s[blockDim.x*blockIdx.x+threadIdx.x]=total;
   __syncthreads();
}

"""

KERNEL_CODE_OPENCL="""
#define TCONG 0
#define TSHR3 1
#define TMWC 2
#define TKISS 3

#define TINT32 0
#define TINT64 1
#define TFP32 2
#define TFP64 3

// Marsaglia RNG very simple implementation
#define znew  ((z=36969*(z&65535)+(z>>16))<<16)
#define wnew  ((w=18000*(w&65535)+(w>>16))&65535)

#define MWC   (znew+wnew)
#define SHR3  (jsr=(jsr=(jsr=jsr^(jsr<<17))^(jsr>>13))^(jsr<<5))
#define CONG  (jcong=69069*jcong+1234567)
#define KISS  ((MWC^CONG)+SHR3)

#define MWCfp MWC * 2.328306435454494e-10f
#define KISSfp KISS * 2.328306435454494e-10f
#define CONGfp CONG * 2.328306435454494e-10f
#define SHR3fp SHR3 * 2.328306435454494e-10f

ulong MainLoop(ulong iterations,uint seed_z,uint seed_w,size_t work)
{

#if TRNG == TCONG
   uint jcong=seed_z+work;
#elif TRNG == TSHR3
   uint jsr=seed_w+work;
#elif TRNG == TMWC
   uint z=seed_z+work;
   uint w=seed_w+work;
#elif TRNG == TKISS
   uint jcong=seed_z+work;
   uint jsr=seed_w+work;
   uint z=seed_z-work;
   uint w=seed_w-work;
#endif

   ulong total=0;

   for (ulong i=0;i<iterations;i++) {

#if TYPE == TINT32
    #define THEONE 1073741824
    #if TRNG == TCONG
        uint x=CONG>>17 ;
        uint y=CONG>>17 ;
    #elif TRNG == TSHR3
        uint x=SHR3>>17 ;
        uint y=SHR3>>17 ;
    #elif TRNG == TMWC
        uint x=MWC>>17 ;
        uint y=MWC>>17 ;
    #elif TRNG == TKISS
        uint x=KISS>>17 ;
        uint y=KISS>>17 ;
    #endif
#elif TYPE == TINT64
    #define THEONE 4611686018427387904
    #if TRNG == TCONG
        ulong x=(ulong)(CONG>>1) ;
        ulong y=(ulong)(CONG>>1) ;
    #elif TRNG == TSHR3
        ulong x=(ulong)(SHR3>>1) ;
        ulong y=(ulong)(SHR3>>1) ;
    #elif TRNG == TMWC
        ulong x=(ulong)(MWC>>1) ;
        ulong y=(ulong)(MWC>>1) ;
    #elif TRNG == TKISS
        ulong x=(ulong)(KISS>>1) ;
        ulong y=(ulong)(KISS>>1) ;
    #endif
#elif TYPE == TFP32
    #define THEONE 1.0f
    #if TRNG == TCONG
        float x=CONGfp ;
        float y=CONGfp ;
    #elif TRNG == TSHR3
        float x=SHR3fp ;
        float y=SHR3fp ;
    #elif TRNG == TMWC
        float x=MWCfp ;
        float y=MWCfp ;
    #elif TRNG == TKISS
      float x=KISSfp ;
      float y=KISSfp ;
    #endif
#elif TYPE == TFP64
#pragma OPENCL EXTENSION cl_khr_fp64: enable
    #define THEONE 1.0f
    #if TRNG == TCONG
        double x=(double)CONGfp ;
        double y=(double)CONGfp ;
    #elif TRNG == TSHR3
        double x=(double)SHR3fp ;
        double y=(double)SHR3fp ;
    #elif TRNG == TMWC
        double x=(double)MWCfp ;
        double y=(double)MWCfp ;
    #elif TRNG == TKISS
        double x=(double)KISSfp ;
        double y=(double)KISSfp ;
    #endif
#endif

      ulong inside=((x*x+y*y) <= THEONE) ? 1:0;
      total+=inside;
   }
   
   return(total);
}

__kernel void MainLoopGlobal(__global ulong *s,ulong iterations,uint seed_w,uint seed_z)
{
   ulong total=MainLoop(iterations,seed_z,seed_w,get_global_id(0));
   barrier(CLK_GLOBAL_MEM_FENCE);
   s[get_global_id(0)]=total;     
}


__kernel void MainLoopLocal(__global ulong *s,ulong iterations,uint seed_w,uint seed_z)
{
   ulong total=MainLoop(iterations,seed_z,seed_w,get_local_id(0));
   barrier(CLK_LOCAL_MEM_FENCE);
   s[get_local_id(0)]=total;
}

__kernel void MainLoopHybrid(__global ulong *s,ulong iterations,uint seed_w,uint seed_z)
{
   ulong total=MainLoop(iterations,seed_z,seed_w,get_global_id(0));
   barrier(CLK_GLOBAL_MEM_FENCE || CLK_LOCAL_MEM_FENCE);
   s[get_global_id(0)]=total;
}

"""

def MetropolisCuda(circle,iterations,steps,jobs,ParaStyle,RNG,ValueType):

  # Avec PyCUDA autoinit, rien a faire !
  
  circleCU = cuda.InOut(circle)

  try:
    mod = SourceModule(KERNEL_CODE_CUDA,options=['--compiler-options','-Wall -DTRNG=%i -DTYPE=%s' % (Marsaglia[RNG],Computing[ValueType])])
  except:
    print "Compilation seems to brake"
  
  MetropolisBlocksCU=mod.get_function("MainLoopBlocks")
  MetropolisJobsCU=mod.get_function("MainLoopThreads")
  MetropolisHybridCU=mod.get_function("MainLoopHybrid")
  
  start = pycuda.driver.Event()
  stop = pycuda.driver.Event()
  
  MyPi=numpy.zeros(steps)
  MyDuration=numpy.zeros(steps)

  if iterations%jobs==0:
    iterationsCL=numpy.uint64(iterations/jobs)
    iterationsNew=iterationsCL*jobs
  else:
    iterationsCL=numpy.uint64(iterations/jobs+1)
    iterationsNew=iterations

  for i in range(steps):
    start.record()
    start.synchronize()
    if ParaStyle=='Blocks':
      MetropolisBlocksCU(circleCU,
                         numpy.uint64(iterationsCL),
                         numpy.uint32(nprnd(2**30/jobs)),
                         numpy.uint32(nprnd(2**30/jobs)),
                         grid=(jobs,1),
                         block=(1,1,1))
      print "%s with (WorkItems/Threads)=(%i,%i) %s method done" % \
            (Alu,jobs,1,ParaStyle)      
    elif ParaStyle=='Hybrid':
      threads=BestThreadsNumber(jobs)
      MetropolisHybridCU(circleCU,
                         numpy.uint64(iterationsCL),
                         numpy.uint32(nprnd(2**30/jobs)),
                         numpy.uint32(nprnd(2**30/jobs)),
                         grid=(jobs,1),
                         block=(threads,1,1))
      print "%s with (WorkItems/Threads)=(%i,%i) %s method done" % \
            (Alu,jobs/threads,threads,ParaStyle)
    else:
      MetropolisJobsCU(circleCU,
                       numpy.uint64(iterationsCL),
                       numpy.uint32(nprnd(2**30/jobs)),
                       numpy.uint32(nprnd(2**30/jobs)),
                       grid=(1,1),
                       block=(jobs,1,1))
      print "%s with (WorkItems/Threads)=(%i,%i) %s method done" % \
            (Alu,jobs,1,ParaStyle)
    stop.record()
    stop.synchronize()
                
    elapsed = start.time_till(stop)*1e-3

    MyDuration[i]=elapsed
    AllPi=4./numpy.float32(iterationsCL)*circle.astype(numpy.float32)
    MyPi[i]=numpy.median(AllPi)
    print MyPi[i],numpy.std(AllPi),MyDuration[i]


  print jobs,numpy.mean(MyDuration),numpy.median(MyDuration),numpy.std(MyDuration),numpy.mean(Iterations/MyDuration),numpy.median(Iterations/MyDuration),numpy.std(Iterations/MyDuration)

  return(numpy.mean(MyDuration),numpy.median(MyDuration),numpy.std(MyDuration),numpy.mean(Iterations/MyDuration),numpy.median(Iterations/MyDuration),numpy.std(Iterations/MyDuration))


def MetropolisOpenCL(circle,iterations,steps,jobs,ParaStyle,Alu,Device,
                     RNG,ValueType):
	
  # Initialisation des variables en les CASTant correctement
    
  if Device==0:
    print "Enter XPU selector based on ALU type: first selected"
    HasXPU=False
    # Default Device selection based on ALU Type
    for platform in cl.get_platforms():
      for device in platform.get_devices():
        deviceType=cl.device_type.to_string(device.type)
        if deviceType=="GPU" and Alu=="GPU" and not HasXPU:
          XPU=device
          print "GPU selected: ",device.name
          HasXPU=True
        if deviceType=="CPU" and Alu=="CPU" and not HasXPU:        
          XPU=device
          print "CPU selected: ",device.name
          HasXPU=True
  else:
    print "Enter XPU selector based on device number & ALU type"
    Id=1
    HasXPU=False
    # Primary Device selection based on Device Id
    for platform in cl.get_platforms():
      for device in platform.get_devices():
        deviceType=cl.device_type.to_string(device.type)
        if Id==Device and Alu==deviceType and HasXPU==False:
          XPU=device
          print "CPU/GPU selected: ",device.name.lstrip()
          HasXPU=True
        Id=Id+1
    if HasXPU==False:
      print "No XPU #%i of type %s found in all of %i devices, sorry..." % \
          (Device,Alu,Id-1)
      return(0,0,0)
				
  # Je cree le contexte et la queue pour son execution
  ctx = cl.Context([XPU])
  queue = cl.CommandQueue(ctx,
                          properties=cl.command_queue_properties.PROFILING_ENABLE)

  # Je recupere les flag possibles pour les buffers
  mf = cl.mem_flags
	
  circleCL = cl.Buffer(ctx, mf.WRITE_ONLY|mf.COPY_HOST_PTR,hostbuf=circle)

  
  MetropolisCL = cl.Program(ctx,KERNEL_CODE_OPENCL).build( \
    options = "-cl-mad-enable -cl-fast-relaxed-math -DTRNG=%i -DTYPE=%s" % (Marsaglia[RNG],Computing[ValueType]))

  i=0

  MyPi=numpy.zeros(steps)
  MyDuration=numpy.zeros(steps)
  
  if iterations%jobs==0:
    iterationsCL=numpy.uint64(iterations/jobs)
    iterationsNew=numpy.uint64(iterationsCL*jobs)
  else:
    iterationsCL=numpy.uint64(iterations/jobs+1)
    iterationsNew=numpy.uint64(iterations)

  for i in range(steps):
		
    if ParaStyle=='Blocks':
      # Call OpenCL kernel
      # (1,) is Global work size (only 1 work size)
      # (1,) is local work size
      # circleCL is lattice translated in CL format
      # SeedZCL is lattice translated in CL format
      # SeedWCL is lattice translated in CL format
      # step is number of iterations
      CLLaunch=MetropolisCL.MainLoopGlobal(queue,(jobs,),None,
                                           circleCL,
                                           numpy.uint64(iterationsCL),
                                           numpy.uint32(nprnd(2**30/jobs)),
                                           numpy.uint32(nprnd(2**30/jobs)))
      print "%s with (WorkItems/Threads)=(%i,%i) %s method done" % \
            (Alu,jobs,1,ParaStyle)
    elif ParaStyle=='Hybrid':
      threads=BestThreadsNumber(jobs)
      # en OpenCL, necessaire de mettre un Global_id identique au local_id
      CLLaunch=MetropolisCL.MainLoopHybrid(queue,(jobs,),(threads,),
                                          circleCL,
                                          numpy.uint64(iterationsCL),
                                          numpy.uint32(nprnd(2**30/jobs)),
                                          numpy.uint32(nprnd(2**30/jobs)))
        
      print "%s with (WorkItems/Threads)=(%i,%i) %s method done" % \
            (Alu,jobs/threads,threads,ParaStyle)
    else:
      # en OpenCL, necessaire de mettre un Global_id identique au local_id
      CLLaunch=MetropolisCL.MainLoopLocal(queue,(jobs,),(jobs,),
                                          circleCL,
                                          numpy.uint64(iterationsCL),
                                          numpy.uint32(nprnd(2**30/jobs)),
                                          numpy.uint32(nprnd(2**30/jobs)))
      print "%s with %i %s done" % (Alu,jobs,ParaStyle)

    CLLaunch.wait()
    cl.enqueue_copy(queue, circle, circleCL).wait()

    elapsed = 1e-9*(CLLaunch.profile.end - CLLaunch.profile.start)

    print circle,numpy.mean(circle),numpy.median(circle),numpy.std(circle)
    MyDuration[i]=elapsed
    AllPi=4./numpy.float32(iterationsCL)*circle.astype(numpy.float32)
    MyPi[i]=numpy.median(AllPi)
    print MyPi[i],numpy.std(AllPi),MyDuration[i]

  circleCL.release()

  print jobs,numpy.mean(MyDuration),numpy.median(MyDuration),numpy.std(MyDuration),numpy.mean(Iterations/MyDuration),numpy.median(Iterations/MyDuration),numpy.std(Iterations/MyDuration)
	
  return(numpy.mean(MyDuration),numpy.median(MyDuration),numpy.std(MyDuration),numpy.mean(Iterations/MyDuration),numpy.median(Iterations/MyDuration),numpy.std(Iterations/MyDuration))


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
  Iterations=100000000
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
  # Marsaglia RNG
  RNG='MWC'
  # Value type : INT32, INT64, FP32, FP64
  ValueType='FP32'
  
  try:
    opts, args = getopt.getopt(sys.argv[1:],"hocfa:g:p:i:s:e:t:r:d:m:v:",["alu=","gpustyle=","parastyle=","iterations=","jobstart=","jobend=","jobstep=","redo=","device=","marsaglia=","valuetype="])
  except getopt.GetoptError:
    print '%s -o (Out of Core Metrology) -c (Print Curves) -f (Fit to Amdahl Law) -a <CPU/GPU/ACCELERATOR> -d <DeviceId> -g <CUDA/OpenCL> -p <Threads/Hybrid/Blocks> -i <Iterations> -s <JobStart> -e <JobEnd> -t <JobStep> -r <RedoToImproveStats> -m <SHR3/CONG/MWC/KISS> -v <INT32/INT64/FP32/FP64> ' % sys.argv[0]
    sys.exit(2)
    
  for opt, arg in opts:
    if opt == '-h':
      print '%s -o (Out of Core Metrology) -c (Print Curves) -f (Fit to Amdahl Law) -a <CPU/GPU/ACCELERATOR> -d <DeviceId> -g <CUDA/OpenCL> -p <Threads/Hybrid/Blocks> -i <Iterations> -s <JobStart> -e <JobEnd> -t <JobStep> -r <RedoToImproveStats> -m <SHR3/CONG/MWC/KISS> -v <INT32/INT64/FP32/FP64>' % sys.argv[0]

      print "\nInformations about devices detected under OpenCL:"
      # For PyOpenCL import
      try:
        import pyopencl as cl
        Id=1
        for platform in cl.get_platforms():
          for device in platform.get_devices():
            deviceType=cl.device_type.to_string(device.type)
            print "Device #%i from %s of type %s : %s" % (Id,platform.vendor.lstrip(),deviceType,device.name.lstrip())
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
    elif opt in ("-m", "--marsaglia"):
      RNG = arg
    elif opt in ("-v", "--valuetype"):
      ValueType = arg
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

  if Alu=='CPU' and GpuStyle=='CUDA':
    print "Alu can't be CPU for CUDA, set Alu to GPU"
    Alu='GPU'

  if ParaStyle not in ('Blocks','Threads','Hybrid'):
    print "%s not exists, ParaStyle set as Threads !" % ParaStyle
    ParaStyle='Threads'

  print "Compute unit : %s" % Alu
  print "Device Identification : %s" % Device
  print "GpuStyle used : %s" % GpuStyle
  print "Parallel Style used : %s" % ParaStyle
  print "Iterations : %s" % Iterations
  print "Number of threads on start : %s" % JobStart
  print "Number of threads on end : %s" % JobEnd
  print "Number of redo : %s" % Redo
  print "Metrology done out of CPU/GPU : %r" % OutMetrology
  print "Type of Marsaglia RNG used : %s" % RNG
  print "Type of variable : %s" % ValueType

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
          deviceType=cl.device_type.to_string(device.type)
          print "Device #%i from %s of type %s : %s" % (Id,platform.vendor.lstrip(),deviceType,device.name.lstrip())

          if Id == Device:
            # Set the Alu as detected Device Type
            Alu=deviceType
          Id=Id+1
    except ImportError:
      print "Platform does not seem to support CUDA"
      
  average=numpy.array([]).astype(numpy.float32)
  median=numpy.array([]).astype(numpy.float32)
  stddev=numpy.array([]).astype(numpy.float32)
  averageRate=numpy.array([]).astype(numpy.float32)
  medianRate=numpy.array([]).astype(numpy.float32)
  stddevRate=numpy.array([]).astype(numpy.float32)

  ExploredJobs=numpy.array([]).astype(numpy.uint32)

  Jobs=JobStart

  while Jobs <= JobEnd:
    avg,med,std=0,0,0
    ExploredJobs=numpy.append(ExploredJobs,Jobs)
    circle=numpy.zeros(Jobs).astype(numpy.uint64)

    if OutMetrology: 
      duration=numpy.array([]).astype(numpy.float32)
      rate=numpy.array([]).astype(numpy.float32)
      for i in range(Redo):
        start=time.time()
        if GpuStyle=='CUDA':
          try:
            a,m,s,aR,mR,sR=MetropolisCuda(circle,Iterations,1,Jobs,ParaStyle,RNG,ValueType)
          except:
            print "Problem with %i // computations on Cuda" % Jobs
        elif GpuStyle=='OpenCL':
          try:
            a,m,s,aR,mR,sR=MetropolisOpenCL(circle,Iterations,1,Jobs,ParaStyle,Alu,Device,RNG,ValueType)
          except:
            print "Problem with %i // computations on OpenCL" % Jobs            
        duration=numpy.append(duration,time.time()-start)
        rate=numpy.append(rate,Iterations/(time.time()-start))
      if (a,m,s) != (0,0,0): 
        avg=numpy.mean(duration)
        med=numpy.median(duration)
        std=numpy.std(duration)
        avgR=numpy.mean(Iterations/duration)
        medR=numpy.median(Iterations/duration)
        stdR=numpy.std(Iterations/duration)
      else:
        print "Values seem to be wrong..."
    else:
      if GpuStyle=='CUDA':
        try:
          avg,med,std,avgR,medR,stdR=MetropolisCuda(circle,Iterations,Redo,Jobs,ParaStyle,RNG,ValueType)
        except:
          print "Problem with %i // computations on Cuda" % Jobs
      elif GpuStyle=='OpenCL':
        try:
          avg,med,std,avgR,medR,stdR=MetropolisOpenCL(circle,Iterations,Redo,Jobs,ParaStyle,Alu,Device,RNG,ValueType)
        except:
          print "Problem with %i // computations on OpenCL" % Jobs            

    if (avg,med,std) != (0,0,0):
      print "jobs,avg,med,std",Jobs,avg,med,std
      average=numpy.append(average,avg)
      median=numpy.append(median,med)
      stddev=numpy.append(stddev,std)
      averageRate=numpy.append(averageRate,avgR)
      medianRate=numpy.append(medianRate,medR)
      stddevRate=numpy.append(stddevRate,stdR)
    else:
      print "Values seem to be wrong..."
    #THREADS*=2
    if len(average)!=0:
      averageRate=averageRate.astype(int)
      medianRate=medianRate.astype(int)
      stddevRate=stddevRate.astype(int)
      numpy.savez("Pi_%s_%s_%s_%s_%s_%s_%i_%.8i_Device%i_%s_%s" % (ValueType,RNG,Alu,GpuStyle,ParaStyle,JobStart,JobEnd,Iterations,Device,Metrology,gethostname()),(ExploredJobs,average,median,stddev,averageRate,medianRate,stddevRate))
      ToSave=[ ExploredJobs,average,median,stddev,averageRate,medianRate,stddevRate ]
      numpy.savetxt("Pi_%s_%s_%s_%s_%s_%s_%i_%.8i_Device%i_%s_%s" % (ValueType,RNG,Alu,GpuStyle,ParaStyle,JobStart,JobEnd,Iterations,Device,Metrology,gethostname()),numpy.transpose(ToSave),fmt='%i %e %e %e %i %i %i')
    Jobs+=JobStep

  if Fit:
    FitAndPrint(ExploredJobs,median,Curves)
