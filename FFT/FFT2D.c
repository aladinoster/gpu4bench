/* 
   Performs matrix multiply

   Thanks for help from aurel32@debian.org
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>
#include <string.h>

#ifdef CUFFT
#include <complex.h>
#include <cufft.h>
#else
#include <fftw3.h>
#endif

#ifdef DOUBLE
#define LENGTH double
#else
#define LENGTH float
#endif

// kind of plan estimations for FFTW3
#define ESTIMATE 1
#define MEASURE 2
#define PATIENT 3
#define EXHAUSTIVE 4

// Default values
#define NTHREADS 1
#define NLOOPS 10
#define DIMENSION 1024

// Redefine FFTW function calls 
#ifdef FLOAT
#define FFTW_COMPLEX fftwf_complex
#define FFTW_PLAN fftwf_plan
#define FFTW_INIT_THREADS fftwf_init_threads
#define FFTW_PLAN_WITH_NTHREADS fftwf_plan_with_nthreads
#define FFTW_MALLOC fftwf_malloc
#define	FFTW_PLAN_DFT_2D fftwf_plan_dft_2d
#define FFTW_EXECUTE fftwf_execute
#define FFTW_DESTROY_PLAN fftwf_destroy_plan
#define FFTW_FREE fftwf_free

#define CUFFTCOMPLEX cufftComplex
#define CUFFT_X2X CUFFT_C2C
#define CUFFTEXECX2X cufftExecC2C

#elif DOUBLE
#define FFTW_COMPLEX fftw_complex
#define FFTW_PLAN fftw_plan
#define FFTW_INIT_THREADS fftw_init_threads
#define FFTW_PLAN_WITH_NTHREADS fftw_plan_with_nthreads
#define FFTW_MALLOC fftw_malloc
#define	FFTW_PLAN_DFT_2D fftw_plan_dft_2d
#define FFTW_EXECUTE fftw_execute
#define FFTW_DESTROY_PLAN fftw_destroy_plan
#define FFTW_FREE fftw_free

#define CUFFTCOMPLEX cufftDoubleComplex
#define CUFFT_X2X CUFFT_Z2Z
#define CUFFTEXECX2X cufftExecZ2Z

#endif

int main(int argc,char **argv)
{
  int i;
  int nloops=NLOOPS,nthreads=NTHREADS,size;
  int type_plan=0;
  double dplan=0.,dcompute=0.;
  
  struct timeval tv1,tv2;
  struct timezone tz;
    
  if ((argc==1)||
      (strcmp(argv[1],"-h")==0)||
      (strcmp(argv[1],"--help")==0))
    {
      printf("\n\tCalculate fake FFTW images\n\n"
	     "\tParameters to give :\n\n"
	     "\t1> Plan ESTIMATE MEASURE PATIENT EXHAUSTIVE\n"
	     "\t2> Size of Image (2^n)\n"
	     "\t3> Number of Loops (integer)\n"
	     "\t4> Number of Threads (integer)\n\n");
      
      return(0);
    }

  // initialization of FFTW threads
  
  if (atoi(argv[4])<1) {
    nthreads=1;
  }
  else {
    nthreads=atoi(argv[4]);
  }
  
  if (atoi(argv[3])<1) {
    nloops=1;
  }
  else {
    nloops=atoi(argv[3]);
  }
  
  if (atoi(argv[2])<1) {
    size=2;
  }
  else {
    size=atoi(argv[2]);
  }
    
  printf("%i %i\n",size,nthreads);
        
#ifdef CUFFT
  cufftHandle plan;
  CUFFTCOMPLEX *in, *devin;

  size_t arraySize = sizeof(CUFFTCOMPLEX)*size*size;
  cudaMallocHost((void**) &in, arraySize);
  cudaMalloc((void**) &devin, arraySize);
    
  // initialisation of arrays
  for (i = 0; i < size*size; i++) {
    in[i].x=1.;
    in[i].y=0.;
  }
  in[0].x=1.;
  in[0].y=0.;

  // First timer to start plan
  gettimeofday(&tv1, &tz);
  
  // Plan & copy to device
  cufftPlan2d(&plan, size, size, CUFFT_X2X);
  cudaMemcpy(devin, in, arraySize, cudaMemcpyHostToDevice);

  // Second timer to end plan      
  gettimeofday(&tv2, &tz);	
  dplan=(double)((tv2.tv_sec-tv1.tv_sec) * 1000000L +		\
		     (tv2.tv_usec-tv1.tv_usec))/1000000.;  

  // First timer to start process
  gettimeofday(&tv1, &tz);

  // Process
  for (i=0;i<nloops;i++) {
    CUFFTEXECX2X(plan, devin, devin, CUFFT_FORWARD);
  }
    
  cudaMemcpy(in, devin, arraySize, cudaMemcpyDeviceToHost);

  printf("%i=%f,%f\n",0,in[0].x,in[0].y);

  // Second timer to end process
  gettimeofday(&tv2, &tz);

  cufftDestroy(plan);
  cudaFreeHost(in);
  cudaFree(devin);

#elif FFTW3      

  FFTW_COMPLEX *in;
  FFTW_PLAN p;

  FFTW_INIT_THREADS();
  FFTW_PLAN_WITH_NTHREADS(nthreads);

  in = (FFTW_COMPLEX*) FFTW_MALLOC(sizeof(FFTW_COMPLEX)*size*size);

  // First timer to start plan
  gettimeofday(&tv1, &tz);
  
  if (strcmp(argv[1],"MEASURE")==0) {
    p = FFTW_PLAN_DFT_2D(size,size,in,in,FFTW_FORWARD,FFTW_MEASURE);
    type_plan=MEASURE;
  }
  
  if (strcmp(argv[1],"PATIENT")==0) {
    p = FFTW_PLAN_DFT_2D(size,size,in,in,FFTW_FORWARD,FFTW_PATIENT);
    type_plan=PATIENT;
  }
  
  if (strcmp(argv[1],"EXHAUSTIVE")==0) {
    p = FFTW_PLAN_DFT_2D(size,size,in,in,FFTW_FORWARD,FFTW_EXHAUSTIVE);
    type_plan=EXHAUSTIVE;
  }
  
  if ((type_plan==0)||(strcmp(argv[1],"ESTIMATE")==0)) {
    p = FFTW_PLAN_DFT_2D(size,size,in,in,FFTW_FORWARD,FFTW_ESTIMATE);
    type_plan=ESTIMATE;
  }
  // Second timer to end plan      
  gettimeofday(&tv2, &tz);	
  dplan=(double)((tv2.tv_sec-tv1.tv_sec) * 1000000L +		\
		 (tv2.tv_usec-tv1.tv_usec))/1000000.;  
  
  // initialisation of arrays
  for (i = 0; i < size*size; i++) {
    (in[i])[0]=0.;
    (in[i])[1]=0.;
  }
  (in[0])[0]=1.;
  (in[0])[1]=0.;
  
  // First timer to start process
  gettimeofday(&tv1, &tz);
  for (i=0;i<nloops;i++) {
    FFTW_EXECUTE(p);
  }
  // Second timer to end process
  gettimeofday(&tv2, &tz);
  FFTW_DESTROY_PLAN(p);
  FFTW_FREE(in); 
#endif

  dcompute=(double)((tv2.tv_sec-tv1.tv_sec) * 1000000L +		\
		    (tv2.tv_usec-tv1.tv_usec))/1000000./(double)nloops;  
  
  printf("%s,%i,%i,%i,%i,%lf,%lf\n",argv[1],
	 size,size,nthreads,nloops,dplan,dcompute);
  
  /* if ((file=fopen(argv[7],"a"))==NULL) */
  /* 	{ */
  /* 	  printf("Impossible to append result in %s file\n", */
  /* 		 argv[7]); */
  /* 	  return(0); */
  /* 	} */
  
  /* fprintf(file,"%s,%i,%i,%i,%i,%lf,%lf\n", */
  /* 	  argv[1],size,size,nthreads,nloops,dplan,dcompute); */
  
  /* fclose(file); */
  
  return 0;
}
