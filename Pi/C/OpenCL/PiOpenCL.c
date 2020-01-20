// Pi Dart Dash in OpenCL in C, illustrative example
//
// Emmanuel Quemener <emmanuel.quemener@gmail.com>
//
// CC BY-NC-SA 2011 : Emmanuel QUEMENER <emmanuel.quemener@gmail.com> 
// Copyleft Cecill v2
//
// -h : print the documentation and detect devices as (platform,device)
//
// classical use:
// #1 OpenCL Plateform ID: get this information with -h option
// #2 OpenCL Device ID: get this information with -h option
// #3 Minimal number of iterations: 
// #4 Parallel Rate: scattering global work in parts executed //ly
// #5 Loops (to improve statistics)
// #6 Type of variables INT32, INT64, FP32, FP64
// 
// To compile : gcc -o PiOpenCL PiOpenCL.c -lOpenCL

#define CL_USE_DEPRECATED_OPENCL_1_2_APIS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <CL/cl.h>
#include <stdint.h>
#include <sys/time.h>

#define TINT32 0
#define TINT64 1
#define TFP32 2
#define TFP64 3

int DetectOpenCLDevices(void) 
{
  int i, j;
  char* value;
  size_t valueSize;
  cl_uint platformCount;
  cl_platform_id* platforms;
  cl_uint deviceCount;
  cl_device_id* devices;
  cl_uint maxComputeUnits;
  cl_int maxWorkGroupSize;
  cl_int maxWorkItemSizes;
  cl_device_type dev_type;

  // get all platforms
  clGetPlatformIDs(0, NULL, &platformCount);
  platforms = (cl_platform_id*) malloc(sizeof(cl_platform_id) * platformCount);
  clGetPlatformIDs(platformCount, platforms, NULL);

  printf("OpenCL statistics: %d platform(s) detected\n\n",platformCount);

  for (i = 0; i < platformCount; i++) {

    // get all devices
    clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, 0, NULL, &deviceCount);
    devices = (cl_device_id*) malloc(sizeof(cl_device_id) * deviceCount);
    clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, deviceCount, devices, NULL);

    // for each device print critical attributes
    for (j = 0; j < deviceCount; j++) {
      
      // print device name
      clGetDeviceInfo(devices[j], CL_DEVICE_NAME, 0, NULL, &valueSize);
      value = (char*) malloc(valueSize);
      clGetDeviceInfo(devices[j], CL_DEVICE_NAME, valueSize, value, NULL);
      printf("Device (%d,%d): %s\n",i, j, value);
      free(value);

      // print type device CPU/GPU/ACCELERATOR
      clGetDeviceInfo(devices[j], CL_DEVICE_TYPE, sizeof(dev_type), &dev_type, NULL);
      printf("\tDevice Type: ");
      if(dev_type & CL_DEVICE_TYPE_GPU)
	printf("CL_DEVICE_TYPE_GPU ");
      if(dev_type & CL_DEVICE_TYPE_CPU)
	printf("CL_DEVICE_TYPE_CPU ");
      if(dev_type & CL_DEVICE_TYPE_ACCELERATOR)
	printf("CL_DEVICE_TYPE_ACCELERATOR ");
      if(dev_type & CL_DEVICE_TYPE_DEFAULT)
	printf("CL_DEVICE_TYPE_DEFAULT ");
      printf("\n");

      // print device vendor
      clGetDeviceInfo(devices[j], CL_DEVICE_VENDOR, 0, NULL, &valueSize);
      value = (char*) malloc(valueSize);
      clGetDeviceInfo(devices[j], CL_DEVICE_VENDOR, valueSize, value, NULL);
      printf("\tDevice vendor: %s\n", value);
      free(value);

      // print hardware device version
      clGetDeviceInfo(devices[j], CL_DEVICE_VERSION, 0, NULL, &valueSize);
      value = (char*) malloc(valueSize);
      clGetDeviceInfo(devices[j], CL_DEVICE_VERSION, valueSize, value, NULL);
      printf("\tHardware version: %s\n", value);
      free(value);

      // print software driver version
      clGetDeviceInfo(devices[j], CL_DRIVER_VERSION, 0, NULL, &valueSize);
      value = (char*) malloc(valueSize);
      clGetDeviceInfo(devices[j], CL_DRIVER_VERSION, valueSize, value, NULL);
      printf("\tSoftware version: %s\n", value);
      free(value);
      
      // print c version supported by compiler for device
      clGetDeviceInfo(devices[j], CL_DEVICE_OPENCL_C_VERSION, 0, NULL, &valueSize);
      value = (char*) malloc(valueSize);
      clGetDeviceInfo(devices[j], CL_DEVICE_OPENCL_C_VERSION, valueSize, value, NULL);
      printf("\tOpenCL C version: %s\n", value);
      free(value);

      // print parallel compute units
      clGetDeviceInfo(devices[j], CL_DEVICE_MAX_COMPUTE_UNITS,
		      sizeof(maxComputeUnits), &maxComputeUnits, NULL);
      printf("\tParallel compute units: %d\n", maxComputeUnits);
      
      // print max work group size
      clGetDeviceInfo(devices[j], CL_DEVICE_MAX_WORK_GROUP_SIZE,
		      sizeof(maxWorkGroupSize), &maxWorkGroupSize, NULL);
      printf("\tMaximum Work Group Size: %d\n", maxWorkGroupSize);
      
      // print max work items size
      clGetDeviceInfo(devices[j], CL_DEVICE_MAX_WORK_ITEM_SIZES,
		      sizeof(maxWorkItemSizes), &maxWorkItemSizes, NULL);
      printf("\tMaximum Work Item Sizes: %d\n", maxWorkItemSizes);
      
    }
    printf("\n");
    free(devices);
  }

  free(platforms);
  return 0;

}

const char* OpenCLSource[] = {
  "#pragma OPENCL EXTENSION cl_khr_fp64: enable \n",
  "// Marsaglia RNG very simple implementation \n",
  "#define znew  ((z=36969*(z&65535)+(z>>16))<<16) \n",
  "#define wnew  ((w=18000*(w&65535)+(w>>16))&65535) \n",
  "#define MWC   (znew+wnew) \n",
  "#define SHR3  (jsr=(jsr=(jsr=jsr^(jsr<<17))^(jsr>>13))^(jsr<<5)) \n",
  "#define CONG  (jcong=69069*jcong+1234567) \n",
  "#define KISS  ((MWC^CONG)+SHR3) \n",
  "#define MWCfp MWC * 2.328306435454494e-10f \n",
  "#define KISSfp KISS * 2.328306435454494e-10f \n",
  "#define CONGfp CONG * 2.328306435454494e-10f \n",
  "#define SHR3fp SHR3 * 2.328306435454494e-10f \n",
  "#define TINT32 0 \n",
  "#define TINT64 1 \n",
  "#define TFP32 2 \n",
  "#define TFP64 3 \n",
  "#define THEONE32I 1073741824 \n",
  "#define THEONE32F 1.e0f \n",
  "#define THEONE64I 4611686018427387904 \n",
  "#define THEONE64F (double)1.e0f \n",
  "ulong MainLoop32I(ulong iterations,uint seed_z,uint seed_w,size_t work)",
  "{",
  "   uint z=seed_z+work;",
  "   uint w=seed_w+work;",
  "   ulong total=0;",
  "   for (ulong i=0;i<iterations;i++)",
  "   {",
  "      uint x= MWC>>17;",
  "      uint y= MWC>>17;",
  "      ulong inside=((x*x+y*y) <= THEONE32I) ? 1:0;",
  "      total+=inside;",
  "   }",
  "   return(total);",
  "}",
  "ulong MainLoop32F(ulong iterations,uint seed_z,uint seed_w,size_t work)",
  "{",
  "   uint z=seed_z+work;",
  "   uint w=seed_w+work;",
  "   ulong total=0;",
  "   for (ulong i=0;i<iterations;i++)",
  "   {",
  "      float x=(float)MWCfp ;",
  "      float y=(float)MWCfp ;",
  "      ulong inside=((x*x+y*y) <= THEONE32F) ? 1:0;",
  "      total+=inside;",
  "   }",
  "   return(total);",
  "}",
  "ulong MainLoop64I(ulong iterations,uint seed_z,uint seed_w,size_t work)",
  "{",
  "   uint z=seed_z+work;",
  "   uint w=seed_w+work;",
  "   ulong total=0;",
  "   for (ulong i=0;i<iterations;i++)",
  "   {",
  "      ulong x=(ulong)(MWC>>1);",
  "      ulong y=(ulong)(MWC>>1);",
  "      ulong inside=((x*x+y*y) <= THEONE64I) ? 1:0;",
  "      total+=inside;",
  "   }",
  "   return(total);",
  "}",
  "ulong MainLoop64F(ulong iterations,uint seed_z,uint seed_w,size_t work)",
  "{",
  "   uint z=seed_z+work;",
  "   uint w=seed_w+work;",
  "   ulong total=0;",
  "   for (ulong i=0;i<iterations;i++)",
  "{",
  "        double x=(double)MWCfp ;",
  "        double y=(double)MWCfp ;",
  "      ulong inside=((x*x+y*y) <= THEONE64F) ? 1:0;",
  "      total+=inside;",
  "}",
  "   return(total);",
  "}",
  "__kernel void MainLoopGlobal(__global ulong *s,ulong iterations,uint seed_w,uint seed_z,uint MyType)",
  "{",
  "   ulong total;",
  "   if (MyType==TFP32) {",
  "      total=(ulong)MainLoop32F(iterations,seed_z,seed_w,get_global_id(0));",
  "   }",
  "   else if (MyType==TFP64) {",
  "      total=(ulong)MainLoop64F(iterations,seed_z,seed_w,get_global_id(0));",
  "   }",  
  "   else if (MyType==TINT32) {",
  "      total=(ulong)MainLoop32I(iterations,seed_z,seed_w,get_global_id(0));",
  "   }",  
  "   else if (MyType==TINT64) {",
  "      total=(ulong)MainLoop64I(iterations,seed_z,seed_w,get_global_id(0));",
  "   }",  
  "   barrier(CLK_GLOBAL_MEM_FENCE);",
  "   s[get_global_id(0)]=(ulong)total;",
  "}",
  "__kernel void MainLoopLocal(__global ulong *s,ulong iterations,uint seed_w,uint seed_z,uint MyType)",
  "{",
  "   ulong total;",
  "   if (MyType==TFP32) {",
  "      total=(ulong)MainLoop32F(iterations,seed_z,seed_w,get_local_id(0));",
  "   }",
  "   else if (MyType==TFP64) {",
  "      total=(ulong)MainLoop64F(iterations,seed_z,seed_w,get_local_id(0));",
  "   }",  
  "   else if (MyType==TINT32) {",
  "      total=(ulong)MainLoop32I(iterations,seed_z,seed_w,get_local_id(0));",
  "   }",  
  "   else if (MyType==TINT64) {",
  "      total=(ulong)MainLoop64I(iterations,seed_z,seed_w,get_local_id(0));",
  "   }",  
  "   barrier(CLK_LOCAL_MEM_FENCE);",
  "   s[get_local_id(0)]=(ulong)total;",
  "}",
  "__kernel void MainLoopHybrid(__global ulong *s,ulong iterations,uint seed_w,uint seed_z,uint MyType)",
  "{",
  "   ulong total;",
  "   if (MyType==TFP32) {",
  "      total=(ulong)MainLoop32F(iterations,seed_z,seed_w,get_global_id(0));",
  "   }",
  "   else if (MyType==TFP64) {",
  "      total=(ulong)MainLoop64F(iterations,seed_z,seed_w,get_global_id(0));",
  "   }",  
  "   else if (MyType==TINT32) {",
  "      total=(ulong)MainLoop32I(iterations,seed_z,seed_w,get_global_id(0));",
  "   }",  
  "   else if (MyType==TINT64) {",
  "      total=(ulong)MainLoop64I(iterations,seed_z,seed_w,get_global_id(0));",
  "   }",  
  "   barrier(CLK_GLOBAL_MEM_FENCE || CLK_LOCAL_MEM_FENCE);",
  "   s[get_global_id(0)]=total;",
  "}"
};

int main(int argc, char **argv)
{
  if ((argc==1)||
      (strcmp(argv[1],"-h")==0)||
      (strcmp(argv[1],"--help")==0))
    {
      printf("\nPerforms a Pi estimation by Dart Dash:\n\n"
	     "\t#1 OpenCL Plateform ID (default 0)\n"
	     "\t#2 OpenCL Device ID (default 0)\n"
	     "\t#3 Minimal number of iterations (default 1000000)\n"
	     "\t#4 Parallel Rate (default 1024)\n"
	     "\t#5 Loops (default 1)\n"
	     "\t#6 Type of variable: INT32, INT64, FP32, FP64 (default FP32)\n\n");
      DetectOpenCLDevices();
    }
  else
    {
      
      int MyPlatform=atoi(argv[1]);
      int MyDevice=atoi(argv[2]);

      struct timeval tv1,tv2;
      struct timezone tz;
      
      uint64_t Iterations=1000000;
      if (argc>3) {
	Iterations=(uint64_t)atoll(argv[3]);
      }
      
      uint32_t ParallelRate=64;
      if (argc>4) {
	ParallelRate=(uint32_t)atoi(argv[4]);
      }
      
      uint32_t Loops=1;
      if (argc>5) {
	Loops=(uint32_t)atoi(argv[5]);
      }
      
      uint32_t MyType=TFP32;
      if (argc>6) {
	if (strcmp(argv[6],"INT32")==0) {
	  MyType=(uint32_t)TINT32;
	}
	else if (strcmp(argv[6],"INT64")==0) {
	  MyType=(uint32_t)TINT64;
	}
	else if (strcmp(argv[6],"FP32")==0) {
	  MyType=(uint32_t)TFP32;
	}
	else if (strcmp(argv[6],"FP64")==0) {
	  MyType=(uint32_t)TFP64;
	}
      }

      printf("MyType %d\n",MyType);
      
      cl_int err;
      cl_context_properties props[3] = { CL_CONTEXT_PLATFORM, 0, 0 };
      
      // Detect, scan, get & initialize platform and device
      cl_uint platformCount;
      cl_platform_id* platforms;
      cl_uint deviceCount;
      cl_device_id* devices;      
      size_t valueSize;
      
      /* Setup OpenCL environment. */
     
      // Get all platforms
      err = clGetPlatformIDs(0, NULL, &platformCount);
      platforms = (cl_platform_id*) malloc(sizeof(cl_platform_id) * platformCount);
      err = clGetPlatformIDs(platformCount, platforms, NULL);

      // Get Device defined
      err = clGetDeviceIDs(platforms[MyPlatform], CL_DEVICE_TYPE_ALL, 0, NULL, &deviceCount);
      devices = (cl_device_id*) malloc(sizeof(cl_device_id) * deviceCount);
      err = clGetDeviceIDs(platforms[MyPlatform], CL_DEVICE_TYPE_ALL, deviceCount, devices, NULL);  

      // print device name
      err = clGetDeviceInfo(devices[MyDevice], CL_DEVICE_NAME, 0, NULL, &valueSize);
      char* deviceName=(char*)malloc(valueSize);
      err = clGetDeviceInfo(devices[MyDevice], CL_DEVICE_NAME, valueSize, deviceName, NULL);
      err = clGetDeviceInfo(devices[MyDevice], CL_DEVICE_VENDOR, 0, NULL, &valueSize);
      char* vendorName=(char*)malloc(valueSize);      
      err = clGetDeviceInfo(devices[MyDevice], CL_DEVICE_VENDOR, valueSize, vendorName, NULL);
      printf("\nDevice (%d,%d):\n\t- vendor: %s\n\t- device: %s\n\n",MyPlatform,MyDevice, vendorName,deviceName);
      free(deviceName);
      free(vendorName);
      
      props[1] = (cl_context_properties)platforms[MyPlatform];
      
      cl_context GPUContext = clCreateContext(props, 1, &devices[MyDevice], NULL, NULL, &err);
      cl_command_queue cqCommandQueue = clCreateCommandQueue(GPUContext,devices[MyDevice], 0, &err);

      cl_mem GPUInside = clCreateBuffer(GPUContext, CL_MEM_WRITE_ONLY,
					      sizeof(uint64_t) * ParallelRate, NULL, NULL);
      
      // 51 is the number of line for OpenCL code
      // 66, sans test
      cl_program OpenCLProgram = clCreateProgramWithSource(GPUContext, 130 ,OpenCLSource,NULL,NULL);
      clBuildProgram(OpenCLProgram, 0, NULL, NULL, NULL, NULL);
      cl_kernel OpenCLMainLoopGlobal = clCreateKernel(OpenCLProgram, "MainLoopGlobal", NULL);
      cl_kernel OpenCLMainLoopLocal = clCreateKernel(OpenCLProgram, "MainLoopLocal", NULL);
      cl_kernel OpenCLMainLoopHybrid = clCreateKernel(OpenCLProgram, "MainLoopHybrid", NULL);

      // Divide the total number of iterations by the parallel rate
      // Add +1 to the number of per work iterations if division not integer
      uint64_t IterationsEach=((Iterations%ParallelRate)==0)?Iterations/ParallelRate:Iterations/ParallelRate+1;
      // Initialize seeds for MWC RNG generator from Marsaglia
      uint32_t seed_w=110271;
      uint32_t seed_z=101008;

      // Set the values of arguments for OpenCL function call
      clSetKernelArg(OpenCLMainLoopGlobal, 0, sizeof(cl_mem),&GPUInside);
      clSetKernelArg(OpenCLMainLoopGlobal, 1, sizeof(uint64_t),&IterationsEach);
      clSetKernelArg(OpenCLMainLoopGlobal, 2, sizeof(uint32_t),&seed_w);
      clSetKernelArg(OpenCLMainLoopGlobal, 3, sizeof(uint32_t),&seed_z);
      clSetKernelArg(OpenCLMainLoopGlobal, 4, sizeof(uint32_t),&MyType);
      
      size_t WorkSize[1] = {ParallelRate}; // one dimensional Range

      uint64_t HostInside[ParallelRate];

      for (uint32_t loop=0;loop<Loops;loop++) {
	// Set start timer
	gettimeofday(&tv1, &tz);
	
    	// Execute the OpenCL kernel with datas
	clEnqueueNDRangeKernel(cqCommandQueue, OpenCLMainLoopGlobal, 1, NULL,
			       WorkSize, NULL, 0, NULL, NULL);
	// Copy each result for each PR from Device to Host
	clEnqueueReadBuffer(cqCommandQueue, GPUInside, CL_TRUE, 0,
			    ParallelRate * sizeof(uint64_t), HostInside, 0, NULL, NULL);
	uint64_t inside=0;

	for (int i= 0; i < ParallelRate; i++) {
	  inside+=HostInside[i];
	}
	  
	// Set stop timer
	gettimeofday(&tv2, &tz);

	double elapsed=(double)((tv2.tv_sec-tv1.tv_sec) * 1000000L +
				(tv2.tv_usec-tv1.tv_usec))/1000000;  

	double itops=(double)(ParallelRate*IterationsEach)/elapsed;
      
	printf("Inside/Total %ld %ld\nParallelRate %i\nElapsed Time %.2f\nItops %.0f\nPi estimation %f\n\n",inside,ParallelRate*IterationsEach,ParallelRate,elapsed,itops,(4.*(float)inside/((float)(ParallelRate)*(float)(IterationsEach))));
      }
      printf("\n\n");
      
      clReleaseKernel(OpenCLMainLoopGlobal);
      clReleaseProgram(OpenCLProgram);
      clReleaseCommandQueue(cqCommandQueue);
      clReleaseContext(GPUContext);
      clReleaseMemObject(GPUInside);

      
      return 0;
    }
}

