#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Demonstrateur OpenCL pour l'ANR Epidevomath

Emmanuel QUEMENER <emmanuel.quemener@ens-lyon.fr> CeCILLv2
"""
import getopt
import sys
import time
import numpy as np
import pyopencl as cl
import pyopencl.array as cl_array
from numpy.random import randint as nprnd


BlobOpenCL= """
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

#define LENGTH 1.

#define PI 3.14159265359

#define SMALL_NUM 0.000000001

__kernel void SplutterPoints(__global float8* clData, float box,
                               uint seed_z,uint seed_w)
{
    int gid = get_global_id(0);
    uint z=seed_z+(uint)gid;
    uint w=seed_w-(uint)gid;

    clData[gid].s01234567 = (float8) (box*MWCfp,box*MWCfp,box*MWCfp,0.,0.,0.,0.,0.);
}

__kernel void ExtendSegment(__global float8* clData, float length,
                               uint seed_z,uint seed_w)
{
    int gid = get_global_id(0);
    uint z=seed_z+(uint)gid;
    uint w=seed_w-(uint)gid;

    float theta=MWCfp*PI;
    float phi=MWCfp*PI*2.;
    float sinTheta=sin(theta);
    clData[gid].s4=clData[gid].s0+length*sinTheta*cos(phi);
    clData[gid].s5=clData[gid].s1+length*sinTheta*sin(phi);
    clData[gid].s6=clData[gid].s2+length*cos(theta);

}

__kernel void EstimateLength(__global float8* clData,__global float* clSize)
{
    int gid = get_global_id(0);
    
    clSize[gid]=distance(clData[gid].lo,clData[gid].hi);
}

// Get from http://geomalgorithms.com/a07-_distance.html

__kernel void ShortestDistance(__global float8* clData,__global float* clDistance)
{
    int gidx = get_global_id(0);
    int ggsz = get_global_size(0);
    int gidy = get_global_id(1);
    
    float4   u = clData[gidx].hi - clData[gidx].lo;
    float4   v = clData[gidy].hi - clData[gidy].lo;
    float4   w = clData[gidx].lo - clData[gidy].lo;     

    float    a = dot(u,u);         // always >= 0
    float    b = dot(u,v);
    float    c = dot(v,v);         // always >= 0
    float    d = dot(u,w);
    float    e = dot(v,w);
   
    float    D = a*c - b*b;        // always >= 0
    float    sc, sN, sD = D;       // sc = sN / sD, default sD = D >= 0
    float    tc, tN, tD = D;       // tc = tN / tD, default tD = D >= 0

    // compute the line parameters of the two closest points
    if (D < SMALL_NUM) { // the lines are almost parallel
        sN = 0.0;         // force using point P0 on segment S1
        sD = 1.0;         // to prevent possible division by 0.0 later
        tN = e;
        tD = c;
    }
    else {                 // get the closest points on the infinite lines
        sN = (b*e - c*d);
        tN = (a*e - b*d);
        if (sN < 0.0) {        // sc < 0 => the s=0 edge is visible
            sN = 0.0;
            tN = e;
            tD = c;
        }
        else if (sN > sD) {  // sc > 1  => the s=1 edge is visible
            sN = sD;
            tN = e + b;
            tD = c;
        }
    }

    if (tN < 0.0) {            // tc < 0 => the t=0 edge is visible
        tN = 0.0;
        // recompute sc for this edge
        if (-d < 0.0)
            sN = 0.0;
        else if (-d > a)
            sN = sD;
        else {
            sN = -d;
            sD = a;
        }
    }
    else if (tN > tD) {      // tc > 1  => the t=1 edge is visible
        tN = tD;
        // recompute sc for this edge
        if ((-d + b) < 0.0)
            sN = 0;
        else if ((-d + b) > a)
            sN = sD;
        else {
            sN = (-d +  b);
            sD = a;
        }
    }
    // finally do the division to get sc and tc
    sc = (fabs(sN) < SMALL_NUM ? 0.0 : sN / sD);
    tc = (fabs(tN) < SMALL_NUM ? 0.0 : tN / tD);

    // get the difference of the two closest points
    float4   dP = w + (sc * u) - (tc * v);  // =  S1(sc) - S2(tc)

    clDistance[ggsz*gidy+gidx]=length(dP);   // return the closest distance
}

 """

if __name__=='__main__':
    
    # Set defaults values
  
    # Id of Device : 1 is for first find !
    Device=1
    # Iterations is integer
    Number=16384
    # Size of box
    SizeOfBox=1000.
    # Size of segment
    LengthOfSegment=1.
    # Redo the last process
    Redo=1

    HowToUse='%s -d <DeviceId> -n <NumberOfSegments> -s <SizeOfBox> -l <LengthOfSegment>'

    try:
        opts, args = getopt.getopt(sys.argv[1:],"hd:n:s:l:r:",["device=","number=","size=","length=","redo="])
    except getopt.GetoptError:
        print HowToUse % sys.argv[0]
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            print HowToUse % sys.argv[0]

            print "\nInformations about devices detected under OpenCL:"
            try:
                Id=1
                for platform in cl.get_platforms():
                    for device in platform.get_devices():
                        deviceType=cl.device_type.to_string(device.type)
                        print "Device #%i from %s of type %s : %s" % (Id,platform.vendor.lstrip(),deviceType,device.name.lstrip())
                        Id=Id+1
                sys.exit()
            except ImportError:
                print "Your platform does not seem to support OpenCL"
                sys.exit()

        elif opt in ("-d", "--device"):
            Device=int(arg)
        elif opt in ("-n", "--number"):
            Number=int(arg)
        elif opt in ("-s", "--size"):
            SizeOfBox=np.float32(arg)
        elif opt in ("-l", "--length"):
            LengthOfSegment=np.float32(arg)
        elif opt in ("-r", "--redo"):
            Redo=int(arg)
            
    print "Device choosed : %s" % Device
    print "Number of segments : %s" % Number
    print "Size of Box : %s" % SizeOfBox
    print "Length of Segment % s" % LengthOfSegment
    print "Redo the last process % s" % Redo
    
    MyData = np.zeros(Number, dtype=cl_array.vec.float8)

    Id=1
    HasXPU=False
    for platform in cl.get_platforms():
        for device in platform.get_devices():
            if Id==Device:
                PlatForm=platform
                XPU=device
                print "CPU/GPU selected: ",device.name.lstrip()
                HasXPU=True
            Id+=1

    if HasXPU==False:
        print "No XPU #%i found in all of %i devices, sorry..." % (Device,Id-1)
        sys.exit()      

    # Je cree le contexte et la queue pour son execution
    try:
        ctx = cl.Context([XPU])
        queue = cl.CommandQueue(ctx,properties=cl.command_queue_properties.PROFILING_ENABLE)
    except:
        print "Crash during context creation"
   

    MyRoutines = cl.Program(ctx, BlobOpenCL).build()

    mf = cl.mem_flags
    clData = cl.Buffer(ctx, mf.READ_WRITE, MyData.nbytes)

    print 'Tous au meme endroit',MyData

    MyRoutines.SplutterPoints(queue,(Number,1),None,clData,np.float32(SizeOfBox-LengthOfSegment),np.uint32(nprnd(2**32)),np.uint32(nprnd(2**32)))

    cl.enqueue_copy(queue, MyData, clData)

    print 'Tous distribues',MyData

    MyRoutines.ExtendSegment(queue,(Number,1),None,clData,np.float32(LengthOfSegment),np.uint32(nprnd(2**32)),np.uint32(nprnd(2**32)))

    cl.enqueue_copy(queue, MyData, clData)

    print 'Tous avec leur extremite',MyData

    MySize = np.zeros(len(MyData), dtype=np.float32)
    clSize = cl.Buffer(ctx, mf.READ_WRITE, MySize.nbytes)

    MyRoutines.EstimateLength(queue, (Number,1), None, clData, clSize)
    cl.enqueue_copy(queue, MySize, clSize)

    print 'La distance de chacun avec son extremite',MySize

    MyDistance = np.zeros(len(MyData)*len(MyData), dtype=np.float32)
    clDistance = cl.Buffer(ctx, mf.READ_WRITE, MyDistance.nbytes)

    time_start=time.time()
    for i in xrange(Redo):
        CLLaunch=MyRoutines.ShortestDistance(queue, (Number,Number), None, clData, clDistance)
        sys.stdout.write('.')
        CLLaunch.wait()
    print "\nDuration on %s for each %s" % (Device,(time.time()-time_start)/Redo)
    cl.enqueue_copy(queue, MyDistance, clDistance)

    MyDistance=np.reshape(MyDistance,(Number,Number))
    clDistance.release()

    print 'La distance de chacun',MyDistance

    clData.release()