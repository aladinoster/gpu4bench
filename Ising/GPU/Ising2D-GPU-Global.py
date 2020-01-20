#!/usr/bin/env python
#
# Ising2D model using PyOpenCL 
#
# CC BY-NC-SA 2011 : <emmanuel.quemener@ens-lyon.fr> 
#
# Thanks to Andreas Klockner for PyOpenCL:
# http://mathema.tician.de/software/pyopencl
# 
# Interesting links:
# http://viennacl.sourceforge.net/viennacl-documentation.html
# http://enja.org/2011/02/22/adventures-in-pyopencl-part-1-getting-started-with-python/

import pyopencl as cl
import numpy
from numpy.random import randint as nprnd
from PIL import Image
import time
import sys
import getopt
import matplotlib.pyplot as plt

# 2097152 on HD5850 (with 1GB of RAM)
#  262144 on GT218
#STEP=262144
#STEP=1048576
#STEP=2097152
#STEP=4194304
#STEP=8388608
STEP=16777216

KERNEL_CODE="""

// Marsaglia RNG very simple implementation

#define znew (z=36969*(z&65535)+(z>>16))
#define wnew (w=18000*(w&65535)+(w>>16))
#define MWC ((znew<<16)+wnew )
#define MWCfp MWC * 2.328306435454494e-10f

__kernel void MainLoop(__global int *s,float J,float B,float T,uint size,
                       uint iterations,uint seed_w,uint seed_z)
{
   uint z=seed_z;
   uint w=seed_w;

   for (uint i=0;i<iterations;i++) {

      uint x=(uint)(MWC%size) ;
      uint y=(uint)(MWC%size) ;

      int p=s[x+size*y];

      int d=s[x+size*((y+1)%size)];
      int u=s[x+size*((y-1)%size)];
      int l=s[((x-1)%size)+size*y];
      int r=s[((x+1)%size)+size*y];

      float DeltaE=p*(2.0f*J*(u+d+l+r)+B);

      int factor=((DeltaE < 0.0f) || (MWCfp < exp(-DeltaE/T))) ? -1:1;
      s[x%size+size*(y%size)] = factor*p;
      // barrier(CLK_GLOBAL_MEM_FENCE);
   }
   barrier(CLK_GLOBAL_MEM_FENCE);
   
}
"""

def ImageOutput(sigma,prefix):
	
	Max=sigma.max()
	Min=sigma.min()
	
	# Normalize value as 8bits Integer
	SigmaInt=(255*(sigma-Min)/(Max-Min)).astype('uint8')
	image = Image.fromarray(SigmaInt)
	image.save("%s.jpg" % prefix)
	
def Metropolis(sigma,J,B,T,iterations,Device):
		
	# Initialisation des variables en les CASTant correctement
    
	# Je detecte un peripherique GPU dans la liste des peripheriques
    Id=1
    HasXPU=False
    for platform in cl.get_platforms():
        for device in platform.get_devices():
            if Id==Device:
                XPU=device
                print "CPU/GPU selected: ",device.name.lstrip()
                HasXPU=True
            Id+=1

    if HasXPU==False:
        print "No XPU #%i found in all of %i devices, sorry..." % (Device,Id-1)
        sys.exit()
     
    # Je cree le contexte et la queue pour son execution
    ctx = cl.Context([XPU])
    queue = cl.CommandQueue(ctx,
        properties=cl.command_queue_properties.PROFILING_ENABLE)

	# Je recupere les flag possibles pour les buffers
    mf = cl.mem_flags
	
    sigmaCL = cl.Buffer(ctx, mf.WRITE_ONLY|mf.COPY_HOST_PTR,hostbuf=sigma)

    MetropolisCL = cl.Program(ctx,KERNEL_CODE).build(
        options = "-cl-mad-enable -cl-fast-relaxed-math")

    i=0
    step=STEP
    duration=0.
	
    while (step*i < iterations):
        # Call OpenCL kernel
        # (1,) is global work size (only 1 work size)
        # (1,) is local work size
        # sigmaCL is lattice translated in CL format
        # step is number of iterations
        
        start_time=time.time()
        CLLaunch=MetropolisCL.MainLoop(queue,(1,),None,
					        sigmaCL,
					        numpy.float32(J),numpy.float32(B),
					        numpy.float32(T),
					        numpy.uint32(sigma.shape[0]),
					        numpy.uint32(step),
					        numpy.uint32(nprnd(2**32)),
					        numpy.uint32(nprnd(2**32)))
        CLLaunch.wait()
        # Does not seem to work under AMD/ATI
        # elapsed = 1e-9*(CLLaunch.profile.end - CLLaunch.profile.start)
        elapsed = time.time()-start_time
        print "Iteration %i with T=%f and %i iterations in %f: " % (i,T,step,elapsed)
        if LAPIMAGE:
            cl.enqueue_copy(queue, sigma, sigmaCL).wait()
            ImageOutput(sigma,"Ising2D_GPU_Global_%i_%1.1f_%.3i_Lap" %  
                (SIZE,T,i))
        i=i+1
        duration=duration+elapsed

	cl.enqueue_copy(queue, sigma, sigmaCL).wait()
	sigmaCL.release()
	
	return(duration)

def Magnetization(sigma,M):
    return(numpy.sum(sigma)/(sigma.shape[0]*sigma.shape[1]*1.0))

def Energy(sigma,J,B):
    # Copier et caster 
    E=numpy.copy(sigma).astype(numpy.float32)

    # Appel par slice
    E[1:-1,1:-1]=E[1:-1,1:-1]*(2.0*J*(E[:-2,1:-1]+E[2:,1:-1]+
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
    # External Magnetic Field is null
    B=0.
    # Size of Lattice
    Size=256
    # Default Temperatures (start, end, step)
    Tmin=0.1
    Tmax=5
    Tstep=0.1
    # Default Number of Iterations
    Iterations=Size*Size*Size
    # Default Device is first
    Device=1

    # Curves is True to print the curves
    Curves=False

    OCL_vendor={}
    OCL_type={}
    OCL_description={}

    try:
        import pyopencl as cl
 
        print "\nHere are available OpenCL devices:"
        Id=1
        for platform in cl.get_platforms():
            for device in platform.get_devices():
                OCL_vendor[Id]=platform.vendor.lstrip().rstrip()
                #OCL_type[Id]=cl.device_type.to_string(device.type)
                OCL_type[Id]="xPU"
                OCL_description[Id]=device.name.lstrip().rstrip()
                print "* Device #%i from %s of type %s : %s" % (Id,OCL_vendor[Id],OCL_type[Id],OCL_description[Id])
                Id=Id+1
        OCL_MaxDevice=Id-1
        print
        
    except ImportError:
        print "Your platform does not seem to support OpenCL"
        sys.exit(0)   
    
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hcj:b:z:i:s:e:p:d:",["coupling=","magneticfield=","size=","iterations=","tempstart=","tempend=","tempstep=","units",'device='])
    except getopt.GetoptError:
        print '%s -d <Device Id> -j <Coupling Factor> -b <Magnetic Field> -z <Size of Square Lattice> -i <Iterations> -s <Minimum Temperature> -e <Maximum Temperature> -p <steP Temperature> -c (Print Curves)' % sys.argv[0]
        sys.exit(2)
    
    for opt, arg in opts:
        if opt == '-h':
            print '%s -d <Device Id> -j <Coupling Factor> -b <Magnetic Field> -z <Size of Square Lattice> -i <Iterations> -s <Minimum Temperature> -e <Maximum Temperature> -p <steP Temperature> -c (Print Curves)' % sys.argv[0]
            sys.exit()
        elif opt == '-c':
            Curves=True
        elif opt in ("-d", "--device"):
            Device = int(arg)
            if Device>OCL_MaxDevice:
                "Device #%s seems not to be available !"
                sys.exit()
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

    print "Here are parameters of simulation:"
    print "* Device selected #%s: %s of type %s from %s" % (Device,OCL_description[Device],OCL_type[Device],OCL_vendor[Device])
    print "* Coupling Factor J : %s" % J
    print "* Magnetic Field B :  %s" % B
    print "* Size of lattice : %sx%s" % (Size,Size)
    print "* Iterations : %s" % Iterations
    print "* Temperatures from %s to %s by %s" % (Tmin,Tmax,Tstep)

    LAPIMAGE=False
	
    if Iterations<STEP:
        STEP=Iterations
    
    sigmaIn=numpy.where(numpy.random.randn(Size,Size)>0,1,-1).astype	(numpy.int32)
	
    ImageOutput(sigmaIn,"Ising2D_GPU_Global_%i_Initial" % (Size))

    Trange=numpy.arange(Tmin,Tmax+Tstep,Tstep)

    E=[]
    M=[]

    for T in Trange:
 	    sigma=numpy.copy(sigmaIn)
	    duration=Metropolis(sigma,J,B,T,Iterations,Device)
	    E=numpy.append(E,Energy(sigma,J,B))
	    M=numpy.append(M,Magnetization(sigma,B))	    
	    ImageOutput(sigma,"Ising2D_GPU_Global_%i_%1.1f_Final" % (Size,T))
	    print "GPU Time : %f" % (duration)
	    print "Total Energy at Temperature %f : %f" % (T,E[-1])
	    print "Total Magnetization at Temperature %f : %f" % (T,M[-1])

    # Save output
    numpy.savez("Ising2D_GPU_Global_%i_%.8i" % (Size,Iterations),(Trange,E,M))
    
    # Estimate Critical temperature
    print "The critical temperature on %ix%i lattice with J=%f, B=%f is %f " % (Size,Size,J,B,CriticalT(Trange,E))

    if Curves:
        DisplayCurves(Trange,E,M,J,B)

