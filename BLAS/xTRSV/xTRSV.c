/* 
   Performs a linear system solving of random generated system
   Estimates a test

   Matrix is triangular
   
   Thanks for help from aurel32@debian.org
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>
#include <string.h>

#ifdef CUBLAS
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

#ifdef DOUBLE
#define LENGTH double
#else
#define LENGTH float
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
  */
  int incx=1;
#ifdef PRINT
  LENGTH factor=1.;
#endif

  LENGTH alpha=1.,beta=0.,beta2=-1.;
  LENGTH *A,*X,*Y;

  /* checkBefore checkAfter checks */
  LENGTH *checksA,*checksB;

  int i=0, j=0;

  double duration;

  struct timeval tv1,tv2;
  struct timezone tz;

  /* Create 1 Matrix and 2 Vectors of dimension dim  */

  A=malloc(dim*dim*sizeof(LENGTH));
  X=malloc(dim*sizeof(LENGTH));
  Y=malloc(dim*sizeof(LENGTH));

  /* Create 2 vectors for checker Before and After */

  checksA=malloc(RUNS*sizeof(double));
  checksB=malloc(RUNS*sizeof(double));

  /* Initialize elements with random numbers */
  /* Initialize the seed for rand() */
  /* srand(time()); */

#ifdef UNIT
  /* Fill the matrix and vector with random numbers */
  for (i=0; i<dim; i++) {
    for (j=0; j<dim; j++) 
      if (j>=i)
	{
	  /* Normalization is necessary to avoid problems */
	  A[i*dim+j]=1.;
	}
      else
	{
	   A[i*dim+j]=0.;
	}
    X[i]=1;
  }
#else
  for (i=0; i<dim; i++) {
    for (j=0; j<dim; j++) 
      if (j>i)
	{
	  /* Normalization is necessary to avoid problems */
	  /* A[i*dim+j]=(LENGTH)rand()/(RAND_MAX+1.) */
	  /*   *(LENGTH)(i+1.)/(LENGTH)(j+1.); */
	  A[i*dim+j]=(LENGTH)rand()/(RAND_MAX+1.)
	    /(LENGTH)(dim-j);
	}
      else if (j==i)
	{
	   A[i*dim+j]=1.;
	}
      else
	{
	   A[i*dim+j]=0.;
	}
    X[i]=(LENGTH)rand()/(RAND_MAX+1.);
  }
#endif

  /* Print the matrix */

#ifdef QUIET
#else
  for (i=0; i<dim; i++) {
    for (j=0; j<dim; j++) printf("A[%i,%i]=%1.5f ", i,j,A[i*dim+j]);
    printf("\tX[%i]=%1.5f ", i,X[i]);
    putchar('\n');
  }
  putchar('\n');
#endif

  /* Get first timer before launching */
  gettimeofday(&tv1, &tz);

  /* Compute with CuBLAS library  */
#ifdef CUBLAS
  LENGTH *devPtrA=0, *devPtrX=0, *devPtrY=0;
  cublasStatus stat1, stat2, stat3;
  struct timeval tv3,tv4;

  /* Order is Row */
  /* Have to swap uplo and trans */
  char uplo='L',trans='T',diag='N';

  printf("Using CuBLAS: %i iterations for %ix%i matrix\n",
	 RUNS,dim,dim);

  stat1=cublasAlloc(dim*dim,sizeof(devPtrA[0]),(void**)&devPtrA);
  stat2=cublasAlloc(dim,sizeof(devPtrX[0]),(void**)&devPtrX);
  stat3=cublasAlloc(dim,sizeof(devPtrY[0]),(void**)&devPtrY);

  if ((stat1 != CUBLAS_STATUS_SUCCESS) || 
      (stat2 != CUBLAS_STATUS_SUCCESS) ||
      (stat3 != CUBLAS_STATUS_SUCCESS)) {
    wrapperError ("Dtrsv", CUBLAS_WRAPPER_ERROR_ALLOC);
    cublasFree (devPtrA);
    cublasFree (devPtrX);
    cublasFree (devPtrY);
    return 1;
  }

  stat1=cublasSetMatrix(dim,dim,sizeof(A[0]),A,dim,devPtrA,dim);
  stat2=cublasSetVector(dim,sizeof(X[0]),X,incx,devPtrX,incx);
  stat3=cublasSetVector(dim,sizeof(Y[0]),Y,incx,devPtrY,incx);
  
  if ((stat1 != CUBLAS_STATUS_SUCCESS) ||
      (stat2 != CUBLAS_STATUS_SUCCESS) ||
      (stat3 != CUBLAS_STATUS_SUCCESS)) {
    wrapperError ("Dtrsv", CUBLAS_WRAPPER_ERROR_SET);
    cublasFree (devPtrA);
    cublasFree (devPtrX);
    cublasFree (devPtrY);
    return 1;
  }

  /* Get third timer after memory operation */
  gettimeofday(&tv3, &tz);

  for (i=0;i<RUNS;i++)
    {
#ifdef DOUBLE

      printVectorGPU(dim,devPtrX,"X","Roots");

      /* Multiply Y <- A.X */
      cublasDgemv(trans,dim,dim,alpha,devPtrA,dim,
		  devPtrX,incx,beta,devPtrY,incx);

      printVectorGPU(dim,devPtrY,"Y","Results");

      /* Solve linear system A.X=Y : Y <- A-1.Y */
      cublasDtrsv(uplo,trans,diag,dim,devPtrA,dim,devPtrY,incx);

      printVectorGPU(dim,devPtrY,"Y","Solutions");

      /* Estimate the difference between X and Y : Y <- -Y+X */
      cublasDaxpy(dim,beta2,devPtrY,incx,devPtrX,incx);

      printVectorGPU(dim,devPtrX,"X","Errors");

      /* Estimate the second checker */
      checksA[i]=(double)cublasDnrm2(dim,devPtrX,incx);

      /* Swap vector X and Y */
      cublasDswap(dim,devPtrX,incx,devPtrY,incx);

#else

      printVectorGPU(dim,devPtrX,"X","Roots");

      /* Multiply Y <- A.X */
      cublasSgemv(trans,dim,dim,alpha,devPtrA,dim,
		  devPtrX,incx,beta,devPtrY,incx);

      printVectorGPU(dim,devPtrY,"Y","Results");

      /* Solve linear system Y <- A-1.Y */
      cublasStrsv(uplo,trans,diag,dim,devPtrA,dim,devPtrY,incx);

      printVectorGPU(dim,devPtrY,"Y","Solutions");

      /* Add vectors X and -Y */
      cublasSaxpy(dim,beta2,devPtrY,incx,devPtrX,incx);

      printVectorGPU(dim,devPtrX,"X","Errors");

      /* Estimate the second checker */
      checksA[i]=(double)cublasSnrm2(dim,devPtrX,incx);

      /* Swap vector X and Y */
      cublasSswap(dim,devPtrX,incx,devPtrY,incx);

#endif
  
    }

  stat1=cublasGetMatrix(dim,dim,sizeof(A[0]),devPtrA,dim,A,dim);
  stat2=cublasGetVector(dim,sizeof(X[0]),devPtrX,incx,X,incx);
  stat3=cublasGetVector(dim,sizeof(Y[0]),devPtrY,incx,Y,incx);
  
  cublasFree (devPtrA);
  cublasFree (devPtrX);
  cublasFree (devPtrY);
  
  if ((stat1 != CUBLAS_STATUS_SUCCESS) ||
      (stat2 != CUBLAS_STATUS_SUCCESS) ||
      (stat3 != CUBLAS_STATUS_SUCCESS)) {
    wrapperError ("LinearSystem", CUBLAS_WRAPPER_ERROR_GET);
  }
  
  /* Get fourth timer after memory free */
  gettimeofday(&tv4, &tz);

#elif THUNKING
  
  /* Order is Row : Have to swap uplo='U' and trans='N' */
  char uplo='L',trans='T',diag='N';
  printf("Using CuBLAS/Thunking: %i iterations for %ix%i matrix\n",
	 RUNS,dim,dim);

  for (i=0;i<RUNS;i++)
    {
#ifdef DOUBLE
      
      printVector(dim,X,"X","Roots");
      
      /* Multiply A by X as Y <- A.X */
      CUBLAS_DGEMV(&trans,&dim,&dim,&alpha,A,&dim,X,&incx,&beta,Y,&incx);
      
      printVector(dim,Y,"Y","Results");

      /* Solve linear system */
      CUBLAS_DTRSV(&uplo,&trans,&diag,&dim,A,&dim,Y,&incx);
      
      printVector(dim,Y,"Y","Solutions");

      /* Compare the roots X and Y */
      CUBLAS_DAXPY(&dim,&beta2,Y,&incx,X,&incx);

      printVector(dim,X,"X","Errors");

      /* Store the checker of errors */
      checksA[i]=(double)CUBLAS_DNRM2(&dim,X,&incx);

      /* Swap vector X and Y */
      CUBLAS_DSWAP(&dim,X,&incx,Y,&incx);
#else

      printVector(dim,X,"X","Roots");
      
      /* Multiply A by X as Y <- A.X */
      CUBLAS_SGEMV(&trans,&dim,&dim,&alpha,A,&dim,X,&incx,&beta,Y,&incx);
      
      printVector(dim,Y,"Y","Results");

      /* Solve linear system */
      CUBLAS_STRSV(&uplo,&trans,&diag,&dim,A,&dim,Y,&incx);
      
      printVector(dim,Y,"Y","Solutions");

      /* Compare the roots X and Y */
      CUBLAS_SAXPY(&dim,&beta2,Y,&incx,X,&incx);

      printVector(dim,X,"X","Errors");

      /* Store the checker of errors */
      checksA[i]=(double)CUBLAS_SNRM2(&dim,X,&incx);

      /* Swap vector X and Y */
      CUBLAS_SSWAP(&dim,X,&incx,Y,&incx);
#endif

#ifdef PRINT
      printf("Iteration %i, checker is %2.5f and error is %2.10f\n",
	     i,checksA[i],fabs(checksB[i]-checksA[i])/factor);
#endif
    }

#elif FBLAS
  
  /* Order is Row : Have to swap uplo='U' and trans='N' */
  char uplo='L',trans='T',diag='N';
  
  printf("Using FBLAS: %i iterations for %ix%i matrix\n",
	 RUNS,dim,dim);
  
  for (i=0;i<RUNS;i++)
    {
#ifdef DOUBLE
      
      printVector(dim,X,"X","Initial roots");
      
      /* Multiply A by X as Y <- A.X */
      dgemv_(&trans,&dim,&dim,&alpha,A,&dim,X,&incx,&beta,Y,&incx);
      
      printVector(dim,Y,"Y<-A.X","Estimated results");
      
      /* Solve linear system */
      dtrsv_(&uplo,&trans,&diag,&dim,A,&dim,Y,&incx);
      
      printVector(dim,Y,"X","Solutions from A.X=Y");
      
      /* Compare the roots X and Y */
      daxpy_(&dim,&beta2,Y,&incx,X,&incx);
      
      printVector(dim,X,"X","Differences initial and estimated");
      
      /* Store the checker of errors */
      dnrm2_(&dim,X,&incx,&checksA[i]);
            
      /* Swap vector X and Y */
      dswap_(&dim,X,&incx,Y,&incx);

#else

      printVector(dim,X,"X","Roots");
      
      /* Multiply A by X as Y <- A.X */
      sgemv_(&trans,&dim,&dim,&alpha,A,&dim,X,&incx,&beta,Y,&incx);
      
      printVector(dim,Y,"Y","Results");

      /* Solve linear system */
      strsv_(&uplo,&trans,&diag,&dim,A,&dim,Y,&incx);
      
      printVector(dim,Y,"Y","Solutions");

      /* Compare the roots X and Y */
      saxpy_(&dim,&beta2,Y,&incx,X,&incx);

      printVector(dim,X,"X","Errors");

      /* Store the checker of errors */
      snrm2_(&dim,X,&incx,&checksA[i]);

      /* Swap vector X and Y */
      sswap_(&dim,X,&incx,Y,&incx);
#endif

    }

#elif ACML
  
  /* Order is Row : Have to swap uplo='U' and trans='N' */
  char uplo='L',trans='T',diag='N';
  
  printf("Using ACML: %i iterations for %ix%i matrix\n",
	 RUNS,dim,dim);
  
  for (i=0;i<RUNS;i++)
    {
#ifdef DOUBLE
      
      printVector(dim,X,"X","Roots");
      
      /* Multiply A by X as Y <- A.X */
      dgemv(trans,dim,dim,alpha,A,dim,X,incx,beta,Y,incx);
      
      printVector(dim,Y,"Y","Results");
      
      /* Solve linear system */
      dtrsv(uplo,trans,diag,dim,A,dim,Y,incx);
      
      printVector(dim,Y,"Y","Solutions");
      
      /* Compare the roots X and Y */
      daxpy(dim,beta2,Y,incx,X,incx);
      
      printVector(dim,X,"X","Errors");
      
      /* Store the checker of errors */
      dnrm2_(&dim,X,&incx,&checksA[i]);
            
      /* Swap vector X and Y */
      dswap(dim,X,incx,Y,incx);

#else

      printVector(dim,X,"X","Roots");
      
      /* Multiply A by X as Y <- A.X */
      sgemv(trans,dim,dim,alpha,A,dim,X,incx,beta,Y,incx);
      
      printVector(dim,Y,"Y","Results");

      /* Solve linear system */
      strsv(uplo,trans,diag,dim,A,dim,Y,incx);
      
      printVector(dim,Y,"Y","Solutions");

      /* Compare the roots X and Y */
      saxpy(dim,beta2,Y,incx,X,incx);

      printVector(dim,X,"X","Errors");

      /* Store the checker of errors */
      snrm2_(&dim,X,&incx,&checksA[i]);

      /* Swap vector X and Y */
      sswap(dim,X,incx,Y,incx);
#endif

    }

#elif GSL

  printf("Using GSL: %i iterations for %ix%i matrix\n",RUNS,dim,dim);

  /* 
     RowMajor : Matrix is read row by row
     Upper : the no null elements are on top
     NoTrans : no transposition before estimation
     NonUnit : Matrix is not unit
   */

  for (i=0;i<RUNS;i++)
    {  

#ifdef DOUBLE

      printVector(dim,X,"X","Roots");

      /* Multiply A by X as Y <- A.X */
      cblas_dgemv(CblasRowMajor,CblasNoTrans,
		  dim,dim,alpha,A,dim,X,incx,beta,Y,incx);

      printVector(dim,Y,"Y","Results");

      /* Solve linear system : Y <- A-1.Y */
      cblas_dtrsv(CblasRowMajor,CblasUpper,CblasNoTrans,CblasNonUnit,
		  dim,A,dim,Y,incx);

      printVector(dim,Y,"Y","Solutions");
      
      cblas_daxpy(dim,beta2,Y,incx,X,incx);

      printVector(dim,X,"X","Errors");

      /* Store the checker of errors */
      checksA[i]=(double)cblas_dnrm2(dim,X,incx);

      cblas_dswap(dim,X,incx,Y,incx);
      
#else

      printVector(dim,X,"X","Roots");

      /* Multiply A by X as Y <- A.X */
      cblas_sgemv(CblasRowMajor,CblasNoTrans,
		  dim,dim,alpha,A,dim,X,incx,beta,Y,incx);

      printVector(dim,Y,"Y","Results");

      /* Solve linear system : Y <- A-1.Y */
      cblas_strsv(CblasRowMajor,CblasUpper,CblasNoTrans,CblasNonUnit,
		  dim,A,dim,Y,incx);

      printVector(dim,Y,"Y","Solutions");
      
      cblas_saxpy(dim,beta2,Y,incx,X,incx);

      printVector(dim,X,"X","Errors");

      /* Store the checker of errors */
      checksA[i]=(double)cblas_snrm2(dim,X,incx);

      cblas_sswap(dim,X,incx,Y,incx);
      
#endif
      
    }
#else

  printf("Using CBLAS: %i iterations for %ix%i matrix\n",RUNS,dim,dim);

  /* 
     RowMajor : Matrix is read row bu row
     Upper : the no null elements are on top
     NoTrans : no transposition before estimation
     NonUnit : Matrix is not unit
   */

  for (i=0;i<RUNS;i++)
    {  

#ifdef DOUBLE

      printVector(dim,X,"X","Roots");

      /* Multiply A by X as Y <- A.X */
      cblas_dgemv(CblasRowMajor,CblasNoTrans,
		  dim,dim,alpha,A,dim,X,incx,beta,Y,incx);

      printVector(dim,Y,"Y","Results");

      /* Solve linear system : Y <- A-1.Y */
      cblas_dtrsv(CblasRowMajor,CblasUpper,CblasNoTrans,CblasNonUnit,
		  dim,A,dim,Y,incx);

      printVector(dim,Y,"Y","Solutions");
      
      cblas_daxpy(dim,beta2,Y,incx,X,incx);

      printVector(dim,X,"X","Errors");

      /* Store the checker of errors */
      checksA[i]=(double)cblas_dnrm2(dim,X,incx);

      cblas_dswap(dim,X,incx,Y,incx);
      
#else

      printVector(dim,X,"X","Roots");

      /* Multiply A by X as Y <- A.X */
      cblas_sgemv(CblasRowMajor,CblasNoTrans,
		  dim,dim,alpha,A,dim,X,incx,beta,Y,incx);

      printVector(dim,Y,"Y","Results");

      /* Solve linear system : Y <- A-1.Y */
      cblas_strsv(CblasRowMajor,CblasUpper,CblasNoTrans,CblasNonUnit,
		  dim,A,dim,Y,incx);

      printVector(dim,Y,"Y","Solutions");
      
      cblas_saxpy(dim,beta2,Y,incx,X,incx);

      printVector(dim,X,"X","Errors");

      /* Store the checker of errors */
      checksA[i]=(double)cblas_snrm2(dim,X,incx);
      
      cblas_sswap(dim,X,incx,Y,incx);
      
#endif

    }
#endif
  putchar('\n');

  /* Get second timer after launching */
  gettimeofday(&tv2, &tz);

#ifdef CUBLAS
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

  printResults(RUNS,checksA,"C","Errors cumulated");

  putchar('\n');

  /*
#ifdef PRINT
  for (i=0; i<dim; i++) {
    for (j=0; j<dim; j++) printf("A[%i,%i]=%1.5f ", i,j,A[i*dim+j]);
    putchar('\n');
  }

  for (i=0; i<dim; i++) {
    printf("X[%i]=%2.5f",i,X[i]);
    putchar('\n');
  }
  putchar('\n');
  for (i=0; i<dim; i++) {
    printf("Y[%i]=%2.5f",i,Y[i]);
    putchar('\n');
  }
#endif
  */

  return 0;
}

int main(int argc,char **argv)
{
  if ((argc==1)||
      (strcmp(argv[1],"-h")==0)||
      (strcmp(argv[1],"--help")==0))
    {
      printf("\nPerforms a bench using BLAS library implementation:\n\n"
	     "\t#1 Size on triangular system\n"
	     "\t#2 Number of iterations\n\n");
    }
  else if ((atoi(argv[1])>=2)&&
	   (atoi(argv[2])>=1))
    {
      bench(atoi(argv[1]),atoi(argv[2]));
    }

  return 0;
}
