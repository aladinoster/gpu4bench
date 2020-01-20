/* 
   Performs matrix multiply on several BLAS implementations 
   Copyleft Emmanuel QUEMENER <emmanuel.quemener@gmail.com> under GPLv3

   2014-03-14 : Add clBLAS implementation

   Thanks for help from aurel32@debian.org
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>
#include <string.h>

#ifdef CLBLAS
#include <clBLAS.h>
/* Precise here to avoid new specific bench function */
int MyPlatform;
int MyDevice;
#elif CUBLAS
#include <cublas.h>
#define CUBLAS_WRAPPER_ERROR_NOERR      0
#define CUBLAS_WRAPPER_ERROR_ALLOC      1
#define CUBLAS_WRAPPER_ERROR_SET        2
#define CUBLAS_WRAPPER_ERROR_GET        3
#define CUBLAS_WRAPPER_ERROR_STUB       4
#elif THUNKING
#include <cublas.h>
#include "fortran_common.h"
#include "fortran_thunking.h"
#elif FBLAS
#include <cblas.h>
#include <cblas_f77.h>
#elif GSL
#include <gsl_cblas.h>
#elif ACML
#include <acml.h>
#else
#include <cblas.h>
#include <blaswrap.h>
#endif

#ifdef CLBLAS

#ifdef DOUBLE
#define LENGTH cl_double
#else
#define LENGTH cl_float
#endif

#else

#ifdef DOUBLE
#define LENGTH double
#else
#define LENGTH float
#endif

#endif

#ifdef FBLAS

#ifdef DOUBLE

void F77_dgemm(FCHAR, FCHAR, FINT, FINT, FINT, const double *, const double *, FINT, 
	       const double *, FINT, const double *, double *, FINT);

#else

void F77_sgemm(FCHAR, FCHAR, FINT, FINT, FINT, const float *, const float *, FINT, 
	       const float *, FINT, const float *, float *, FINT);

#endif
#endif

/* Matrix with only defined triangular terms */
/* Even if there are 0 in matrix, must be defined at all ! */

/* Get from fortran.c */

#ifdef CUBLAS
static char *errMsg[5] = 
{
    "no error",
    "allocation error",
    "setVector/setMatrix error",
    "getVector/getMatrix error",
    "not implemented"
};

static void wrapperError (const char *funcName, int error)
{
    printf ("cublas%s wrapper: %s\n", funcName, errMsg[error]);
    fflush (stdout);
}
#endif

int printVector(const int dimVector,const LENGTH *dataVector,
		char *nameVector,char *mesgVector)
{
#ifndef QUIET

  int i;
  printf("\n%s of %s, size %i:\n",mesgVector,nameVector,dimVector);
  for (i=0;i<dimVector;i++)
    {
      printf("%s[%i]=%2.10e\n",nameVector,i,dataVector[i]);
    }
#endif

  return 0;
}
  
int printResults(const int dimVector,const LENGTH *dataVector,
		 char *nameVector,char *mesgVector)
{
#ifdef RESULTS
  int i;

  printf("\n%s of %s, size %i:\n",mesgVector,nameVector,dimVector);
  for (i=0;i<dimVector;i++)
    {
      printf("%s[%i]=%2.10e\n",nameVector,i,dataVector[i]);
    }
#endif
  return 0;
}
  
#ifdef CUBLAS
int printVectorGPU(const int dimVector,const LENGTH *dataVector,
		   char *nameVector,char *mesgVector)
{
#ifndef QUIET
  int i;
  cublasStatus stat;
  LENGTH *P=0;
  int incx=1;

  P=malloc(dimVector*sizeof(LENGTH));
  
  stat=cublasGetVector(dimVector,sizeof(P[0]),dataVector,incx,P,incx);

  if (stat != CUBLAS_STATUS_SUCCESS) {
    wrapperError ("ToGet", CUBLAS_WRAPPER_ERROR_GET);
  }  

  printf("\n%s of %s, size %i:\n",mesgVector,nameVector,dimVector);
  for (i=0;i<dimVector;i++)
    {
      printf("%s[%i]=%2.10e\n",nameVector,i,P[i]);
    }

  free(P);  
#endif

  return 0;
}
#endif

int bench(int dim,int RUNS)
{
  /*
  int dim=1000;
  int RUNS=100;
  int incx=1;
  */
#ifdef PRINT
  LENGTH factor=1.;
#endif

  LENGTH alpha=1.,beta=0.;
  LENGTH *A,*B,*C,*D;

  /* checkBefore checkAfter checks */
  LENGTH *checksA,*checksB;

  int i=0, j=0;

  double duration;

  struct timeval tv1,tv2;
  struct timezone tz;

  /* Create 4 Matrix of dimension dim by dim  */

  A=malloc(dim*dim*sizeof(LENGTH));
  B=malloc(dim*dim*sizeof(LENGTH));
  C=malloc(dim*dim*sizeof(LENGTH));
  D=malloc(dim*dim*sizeof(LENGTH));

  /* Create 2 vectors for checker Before and After */

  checksA=malloc(RUNS*sizeof(LENGTH));
  checksB=malloc(RUNS*sizeof(LENGTH));

  /* Initialize elements with random numbers */
  /* Initialize the seed for rand() */
  /* srand(time()); */

  for (i=0; i<dim; i++) {
    for (j=0; j<dim; j++) {
	A[i*dim+j]=(LENGTH)rand()/(RAND_MAX+1.)
	  *(LENGTH)(i+1.)/(LENGTH)(j+1.); 
	B[i*dim+j]=(LENGTH)rand()/(RAND_MAX+1.)
	  *(LENGTH)(i+1.)/(LENGTH)(j+1.);
	C[i*dim+j]=0.;
	D[i*dim+j]=0.;
    }
  }
  /*
  A[0]=1;
  A[1]=2;
  A[2]=3;
  A[3]=4;
  
  B[0]=5;
  B[1]=6;
  B[2]=7;
  B[3]=8;
  */

  /* Print the matrix */
 
#ifdef QUIET
#else
  for (i=0; i<dim; i++) {
    for (j=0; j<dim; j++) printf("A[%i,%i]=%1.5f ", i,j,A[i*dim+j]);
    putchar('\n');
  }
  putchar('\n');
  for (i=0; i<dim; i++) {
    for (j=0; j<dim; j++) printf("B[%i,%i]=%1.5f ", i,j,B[i*dim+j]);
    putchar('\n');
  }
  putchar('\n');
#endif

 /* Get first timer before launching */
  gettimeofday(&tv1, &tz);

  /* Compute with CLBLAS library  */
#ifdef CLBLAS

  cl_uint platformCount;
  cl_platform_id* platforms;
  cl_uint deviceCount;
  cl_device_id* devices;

  cl_int err,errA,errB,errC,errD;
  cl_context_properties props[3] = { CL_CONTEXT_PLATFORM, 0, 0 };
  cl_context ctx = 0;
  cl_command_queue queue = 0;
  cl_mem bufA, bufB, bufC, bufD;
  cl_event event = NULL;

  char* value;
  size_t valueSize;

  // tv3 Put on Device: Allocate & Write buffer
  // tv4 Compute
  struct timeval tv3,tv4;

  printf("Using CLBLAS: %i iterations for %ix%i matrix on (%d,%d)\n",
	 RUNS,dim,dim,MyPlatform,MyDevice);

  /* Setup OpenCL environment. */
  /* - get all platforms and select MyPlatform */
  /* - get all devices from MyPlatform and select MyDevice */

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
  value = (char*) malloc(valueSize);
  err = clGetDeviceInfo(devices[MyDevice], CL_DEVICE_NAME, valueSize, value, NULL);
  printf("Device (%d,%d): %s\n",MyPlatform,MyDevice, value);
  free(value);

  props[1] = (cl_context_properties)platforms[MyPlatform];

  /* Initialize Context */
  ctx = clCreateContext( props, 1, &devices[MyDevice], NULL, NULL, &err );
  queue = clCreateCommandQueue( ctx, devices[MyDevice], 0, &err );

  /* Setup clBLAS */
  err = clblasSetup( );

  /* Prepare OpenCL memory objects and place matrices inside them. */
  bufA = clCreateBuffer(ctx,CL_MEM_READ_ONLY,dim*dim*sizeof(*A),NULL,&errA );
  bufB = clCreateBuffer(ctx,CL_MEM_READ_ONLY,dim*dim*sizeof(*B),NULL,&errB );
  bufC = clCreateBuffer(ctx,CL_MEM_READ_WRITE,dim*dim*sizeof(*C),NULL,&errC );
  bufD = clCreateBuffer(ctx,CL_MEM_READ_WRITE,dim*dim*sizeof(*D),NULL,&errD );

  errA = clEnqueueWriteBuffer( queue,bufA,CL_TRUE,0,
			       dim*dim*sizeof(*A),A,0,NULL,NULL );
  errB = clEnqueueWriteBuffer( queue, bufB, CL_TRUE,0,
			       dim*dim*sizeof(*B),B,0,NULL,NULL );
  errC = clEnqueueWriteBuffer( queue, bufC, CL_TRUE,0,
  			       dim*dim*sizeof(*C),C,0,NULL,NULL );
  errD = clEnqueueWriteBuffer( queue, bufD, CL_TRUE,0,
  			       dim*dim*sizeof(*D),D,0,NULL,NULL );
  
  /* Get third timer after memory operation */
  gettimeofday(&tv3, &tz);

#ifdef DOUBLE

  for (i=0;i<RUNS;i++)
    {
      err = clblasDgemm( clblasRowMajor,clblasNoTrans,clblasNoTrans, 
			 dim,dim,dim,alpha,bufA,0,dim,bufB,0,dim,beta,
			 bufC,0,dim,1,&queue,0,NULL,&event );

      err = clblasDgemm( clblasRowMajor,clblasTrans,clblasTrans, 
			 dim,dim,dim,alpha,bufB,0,dim,bufA,0,dim,beta,
			 bufD,0,dim,1,&queue,0,NULL,&event );

    }
  
  if (err != CL_SUCCESS) {
    printf("clblasDgemm() failed with %d\n", err);
  }

#else

  for (i=0;i<RUNS;i++)
    {

      err = clblasSgemm( clblasRowMajor,clblasNoTrans,clblasNoTrans, 
			 dim,dim,dim,alpha,bufA,0,dim,bufB,0,dim,beta,
			 bufC,0,dim,1,&queue,0,NULL,&event );

      err = clblasSgemm( clblasRowMajor,clblasTrans,clblasTrans, 
			 dim,dim,dim,alpha,bufB,0,dim,bufA,0,dim,beta,
			 bufD,0,dim,1,&queue,0,NULL,&event );
    }
  
  if (err != CL_SUCCESS) {
    printf("clblasSgemm() failed with %d\n", err);
  }

#endif

  /* Wait for calculations to be finished. */
  err = clWaitForEvents( 1, &event );
  
  /* Get fourth timer after memory free */
  gettimeofday(&tv4, &tz);

  /* Fetch results of calculations from GPU memory. */
  errC = clEnqueueReadBuffer( queue,bufC,CL_TRUE,0,dim*dim * sizeof(*C),
			     C,0,NULL,NULL );

  /* Fetch results of calculations from GPU memory. */
  errD = clEnqueueReadBuffer( queue,bufD,CL_TRUE,0,dim*dim*sizeof(*D),
			     D,0,NULL,NULL );

  /* Release OpenCL memory objects. */
  clReleaseMemObject( bufD );
  clReleaseMemObject( bufC );
  clReleaseMemObject( bufB );
  clReleaseMemObject( bufA );

  /* Finalize work with clBLAS */
  clblasTeardown( );

  /* Release OpenCL working objects. */
  clReleaseCommandQueue( queue );
  clReleaseContext( ctx );
  

  /* Compute with CuBLAS library  */
#elif CUBLAS
  LENGTH *devPtrA=0, *devPtrB=0, *devPtrC=0, *devPtrD=0;
  cublasStatus stat1, stat2, stat3, stat4;
  struct timeval tv3,tv4;

  /* Order is Row */
  /* Have to swap uplo and trans */
  char transa='N',transb='T';

  printf("Using CuBLAS: %i iterations for %ix%i matrix\n",
	 RUNS,dim,dim);

  stat1=cublasAlloc(dim*dim,sizeof(devPtrA[0]),(void**)&devPtrA);
  stat2=cublasAlloc(dim*dim,sizeof(devPtrB[0]),(void**)&devPtrB);
  stat3=cublasAlloc(dim*dim,sizeof(devPtrC[0]),(void**)&devPtrC);
  stat4=cublasAlloc(dim*dim,sizeof(devPtrD[0]),(void**)&devPtrD);

  if ((stat1 != CUBLAS_STATUS_SUCCESS) || 
      (stat2 != CUBLAS_STATUS_SUCCESS) ||
      (stat3 != CUBLAS_STATUS_SUCCESS) ||
      (stat4 != CUBLAS_STATUS_SUCCESS) ) {
    wrapperError ("xGEMM", CUBLAS_WRAPPER_ERROR_ALLOC);
    cublasFree (devPtrA);
    cublasFree (devPtrB);
    cublasFree (devPtrC);
    cublasFree (devPtrD);
    return 1;
  }

  stat1=cublasSetMatrix(dim,dim,sizeof(A[0]),A,dim,devPtrA,dim);
  stat2=cublasSetMatrix(dim,dim,sizeof(B[0]),B,dim,devPtrB,dim);
  stat3=cublasSetMatrix(dim,dim,sizeof(C[0]),C,dim,devPtrC,dim);
  stat4=cublasSetMatrix(dim,dim,sizeof(D[0]),D,dim,devPtrD,dim);
  
  if ((stat1 != CUBLAS_STATUS_SUCCESS) ||
      (stat2 != CUBLAS_STATUS_SUCCESS) ||
      (stat3 != CUBLAS_STATUS_SUCCESS) ||
      (stat4 != CUBLAS_STATUS_SUCCESS) ) {
    wrapperError ("xGEMM", CUBLAS_WRAPPER_ERROR_SET);
    cublasFree (devPtrA);
    cublasFree (devPtrB);
    cublasFree (devPtrC);
    cublasFree (devPtrD);
    return 1;
  }

  /* Get third timer after memory operation */
  gettimeofday(&tv3, &tz);

#ifdef DOUBLE

  for (i=0;i<RUNS;i++)
    {
      cublasDgemm(transa,transa,dim,dim,dim,alpha,devPtrB,dim,
		  devPtrA,dim,beta,devPtrC,dim);
      cublasDgemm(transb,transb,dim,dim,dim,alpha,devPtrA,dim,
		  devPtrB,dim,beta,devPtrD,dim);
    }
  
#else

  for (i=0;i<RUNS;i++)
    {
      cublasSgemm(transa,transa,dim,dim,dim,alpha,devPtrB,dim,
		  devPtrA,dim,beta,devPtrC,dim);
      cublasSgemm(transb,transb,dim,dim,dim,alpha,devPtrA,dim,
		  devPtrB,dim,beta,devPtrD,dim);
    }
  
#endif

  stat3=cublasGetMatrix(dim,dim,sizeof(C[0]),devPtrC,dim,C,dim);
  stat4=cublasGetMatrix(dim,dim,sizeof(D[0]),devPtrD,dim,D,dim);
  
  /* Get fourth timer before memory free */
  gettimeofday(&tv4, &tz);

  cublasFree (devPtrA);
  cublasFree (devPtrB);
  cublasFree (devPtrC);
  cublasFree (devPtrD);
  
  if ((stat1 != CUBLAS_STATUS_SUCCESS) ) {
    wrapperError ("xGEMM", CUBLAS_WRAPPER_ERROR_GET);
  }
  

#elif THUNKING
  
  /* Order is Row : Have to swap uplo='U' and trans='N' */
  char transa='N',transb='T';
  printf("Using CuBLAS/Thunking: %i iterations for %ix%i matrix\n",
	 RUNS,dim,dim);

#ifdef DOUBLE

  for (i=0;i<RUNS;i++)
    {      
      CUBLAS_DGEMM(&transa,&transa,
      		   &dim,&dim,&dim,&alpha,B,&dim,A,&dim,&beta,C,&dim);
      CUBLAS_DGEMM(&transb,&transb,
      		   &dim,&dim,&dim,&alpha,A,&dim,B,&dim,&beta,D,&dim);
    }

#else

  for (i=0;i<RUNS;i++)
    {      
      CUBLAS_SGEMM(&transa,&transa,
      		   &dim,&dim,&dim,&alpha,B,&dim,A,&dim,&beta,C,&dim);
      CUBLAS_SGEMM(&transb,&transb,
      		   &dim,&dim,&dim,&alpha,A,&dim,B,&dim,&beta,D,&dim);
    }
  
#endif

#elif FBLAS
  
  /* Order is Row : Have to swap uplo='U' and trans='N' */
      char transa='N',transb='T';
  
  printf("Using FBLAS: %i iterations for %ix%i matrix\n",
	 RUNS,dim,dim);
  
#ifdef DOUBLE

  for (i=0;i<RUNS;i++)
    {    
      F77_dgemm(&transa,&transa,&dim,&dim,&dim,&alpha,B,&dim,A,&dim,&beta,C,&dim);
      F77_dgemm(&transb,&transb,&dim,&dim,&dim,&alpha,A,&dim,B,&dim,&beta,D,&dim);
    }

#else

  for (i=0;i<RUNS;i++)
    {    
      F77_sgemm(&transa,&transa,&dim,&dim,&dim,&alpha,B,&dim,A,&dim,&beta,C,&dim);
      F77_sgemm(&transb,&transb,&dim,&dim,&dim,&alpha,A,&dim,B,&dim,&beta,D,&dim);
    }

#endif

#elif ACML
  
  /* Order is Row : Have to swap uplo='U' and trans='N' */
      char transa='N',transb='T';
  
  printf("Using ACML: %i iterations for %ix%i matrix\n",
	 RUNS,dim,dim);
  
#ifdef DOUBLE

  for (i=0;i<RUNS;i++)
    {    
      dgemm(transa,transa,dim,dim,dim,alpha,B,dim,A,dim,beta,C,dim);
      dgemm(transb,transb,dim,dim,dim,alpha,A,dim,B,dim,beta,D,dim);
    }

#else

  for (i=0;i<RUNS;i++)
    {    
      sgemm(transa,transa,dim,dim,dim,alpha,B,dim,A,dim,beta,C,dim);
      sgemm(transb,transb,dim,dim,dim,alpha,A,dim,B,dim,beta,D,dim);
    }

#endif

#elif GSL

  printf("Using GSL: %i iterations for %ix%i matrix\n",RUNS,dim,dim);

  /* 
     RowMajor : Matrix is read row by row
     Upper : the no null elements are on top
     NoTrans : no transposition before estimation
     NonUnit : Matrix is not unit
   */

#ifdef DOUBLE

  for (i=0;i<RUNS;i++)
    {  
      cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
		  dim,dim,dim,alpha,A,dim,B,dim,beta,C,dim);
      cblas_dgemm(CblasRowMajor,CblasTrans,CblasTrans,
		  dim,dim,dim,alpha,B,dim,A,dim,beta,D,dim);
    }
  
#else

  for (i=0;i<RUNS;i++)
    {  
      cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
		  dim,dim,dim,alpha,A,dim,B,dim,beta,C,dim);
      cblas_sgemm(CblasRowMajor,CblasTrans,CblasTrans,
		  dim,dim,dim,alpha,B,dim,A,dim,beta,D,dim);
    }
  
#endif
      
#else

  printf("Using CBLAS: %i iterations for %ix%i matrix\n",RUNS,dim,dim);

  /* 
     RowMajor : Matrix is read row bu row
     Upper : the no null elements are on top
     NoTrans : no transposition before estimation
     NonUnit : Matrix is not unit
   */

#ifdef DOUBLE

  for (i=0;i<RUNS;i++)
    {  
      cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
		  dim,dim,dim,alpha,A,dim,B,dim,beta,C,dim);
      cblas_dgemm(CblasRowMajor,CblasTrans,CblasTrans,
		  dim,dim,dim,alpha,B,dim,A,dim,beta,D,dim);
    }
  
#else

  for (i=0;i<RUNS;i++)
    {  
      cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
		  dim,dim,dim,alpha,A,dim,B,dim,beta,C,dim);
      cblas_sgemm(CblasRowMajor,CblasTrans,CblasTrans,
		  dim,dim,dim,alpha,B,dim,A,dim,beta,D,dim);
    }
  
#endif

#endif

  /* Get second timer after launching */
  gettimeofday(&tv2, &tz);

  /* Store the checker of errors */
  checksA[0]=0.;
  
  for (i=0; i<dim; i++) {
    for (j=0; j<dim; j++) {
      checksA[0]=checksA[0]+fabs(D[i*dim+j]-C[j*dim+i]);
    }
  }

  /* Print the matrix */
 
#ifdef QUIET
#else
  for (i=0; i<dim; i++) {
    for (j=0; j<dim; j++) printf("C[%i,%i]=%1.5f ", i,j,C[i*dim+j]);
    putchar('\n');
  }
  putchar('\n');
  for (i=0; i<dim; i++) {
    for (j=0; j<dim; j++) printf("D[%i,%i]=%1.5f ", i,j,D[i*dim+j]);
    putchar('\n');
  }
  putchar('\n');
#endif

  /* Free 1 Matrix and 2 Vectors of dimension dim  */

  free(A);
  free(B);
  free(C);
  free(D);

  putchar('\n');

#ifdef CLBLAS
  double memoryIn,memoryOut;

  memoryIn=(double)((tv3.tv_sec-tv1.tv_sec) * 1000000L +	\
		    (tv3.tv_usec-tv1.tv_usec))/1000000.;  

  memoryOut=(double)((tv2.tv_sec-tv4.tv_sec) * 1000000L +	\
		    (tv2.tv_usec-tv4.tv_usec))/1000000.;  

  duration=(double)((tv4.tv_sec-tv3.tv_sec) * 1000000L +	\
		    (tv4.tv_usec-tv3.tv_usec))/1000000./RUNS;  

  printf("Duration of memory allocation : %2.10f s\n",memoryIn);
  printf("Duration of memory free : %2.10f s\n",memoryOut);
#elif CUBLAS
  double memoryIn,memoryOut;

  memoryIn=(double)((tv3.tv_sec-tv1.tv_sec) * 1000000L +	\
		    (tv3.tv_usec-tv1.tv_usec))/1000000.;  

  memoryOut=(double)((tv2.tv_sec-tv4.tv_sec) * 1000000L +	\
		    (tv2.tv_usec-tv4.tv_usec))/1000000.;  

  duration=(double)((tv4.tv_sec-tv3.tv_sec) * 1000000L +	\
		    (tv4.tv_usec-tv3.tv_usec))/1000000./RUNS;  

  printf("Duration of memory allocation : %2.10f s\n",memoryIn);
  printf("Duration of memory free : %2.10f s\n",memoryOut);
#else
  duration=(double)((tv2.tv_sec-tv1.tv_sec) * 1000000L +	\
		    (tv2.tv_usec-tv1.tv_usec))/1000000./RUNS;  

#endif

  printf("Duration of each cycle : %2.10f s\n",duration);

  printf("Number of GFlops : %2.3f \n",
	 dim*dim*2.*(2.*dim-1)/duration/1000000000.);

  printf("Error %1.10f\n",checksA[0]);
  printResults(RUNS,checksA,"C","Errors cumulated");

  putchar('\n');

  /* Free 2 vectors for checker Before and After */

  free(checksA);
  free(checksB);

  return 0;
}

#ifdef CLBLAS

int DelectOpenCLDevices() 
{
  /* */
  /* Not needed to import CL.h, already done in CLBLAS.h */

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
#endif

int main(int argc,char **argv)
{
  if ((argc==1)||
      (strcmp(argv[1],"-h")==0)||
      (strcmp(argv[1],"--help")==0))
    {
#ifdef CLBLAS
      printf("\nPerforms a bench using BLAS library implementation:\n\n"
	     "\t#1 Size of square matrices \n"
	     "\t#2 Number of iterations \n"
	     "\t#3 OpenCL Plateform ID\n"
	     "\t#4 OpenCL Device ID\n\n");
      DelectOpenCLDevices();
#else
      printf("\nPerforms a bench using BLAS library implementation:\n\n"
	     "\t#1 Size of square matrices \n"
	     "\t#2 Number of iterations\n\n");
#endif
    }
  else if ((atoi(argv[1])>=2)&&
	   (atoi(argv[2])>=1))
    {
#ifdef CLBLAS
      MyPlatform=atoi(argv[3]);
      MyDevice=atoi(argv[4]);
#endif
      bench(atoi(argv[1]),atoi(argv[2]));
    }

  return 0;
}
