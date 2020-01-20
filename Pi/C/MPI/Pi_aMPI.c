//
// Estimation of Pi using Monte Carlo exploration process
// gcc -std=c99 -O3 -o Pi_aMPI Pi_aMPI.c -lm 
// Emmanuel Quemener <emmanuel.quemener@ens-lyon.fr>

// Needed for gethostname
#define _BSD_SOURCE
#include <sys/unistd.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <mpi.h>
#include <stddef.h>

#include <sys/time.h>

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

#define ITERATIONS 1000000000

#ifdef LONG
#define LENGTH long long
#else
#define LENGTH int
#endif

typedef struct compute_node {
        LENGTH inside;
        long int useconds;
} result;

unsigned int rotl(unsigned int value, int shift) {
    return (value << shift) | (value >> (sizeof(value) * CHAR_BIT - shift));
}
 
unsigned int rotr(unsigned int value, int shift) {
    return (value >> shift) | (value << (sizeof(value) * CHAR_BIT - shift));
}

LENGTH MainLoopGlobal(LENGTH iterations,unsigned int seed_w,unsigned int seed_z)
{
#if defined TCONG
   unsigned int jcong=seed_z;
#elif defined TSHR3
   unsigned int jsr=seed_w;
#elif defined TMWC
   unsigned int z=seed_z;
   unsigned int w=seed_w;
#elif defined TKISS
   unsigned int jcong=seed_z;
   unsigned int jsr=seed_w;
   unsigned int z=seed_z;
   unsigned int w=seed_w;
#endif

  LENGTH total=0;

   for (LENGTH i=0;i<iterations;i++) {

#if defined TINT32
    #define THEONE 1073741824
    #if defined TCONG
        unsigned int x=CONG>>17 ;
        unsigned int y=CONG>>17 ;
    #elif defined TSHR3
        unsigned int x=SHR3>>17 ;
        unsigned int y=SHR3>>17 ;
    #elif defined TMWC
        unsigned int x=MWC>>17 ;
        unsigned int y=MWC>>17 ;
    #elif defined TKISS
        unsigned int x=KISS>>17 ;
        unsigned int y=KISS>>17 ;
    #endif
#elif defined TINT64
    #define THEONE 4611686018427387904
    #if defined TCONG
        unsigned long x=(unsigned long)(CONG>>1) ;
        unsigned long y=(unsigned long)(CONG>>1) ;
    #elif defined TSHR3
        unsigned long x=(unsigned long)(SHR3>>1) ;
        unsigned long y=(unsigned long)(SHR3>>1) ;
    #elif defined TMWC
        unsigned long x=(unsigned long)(MWC>>1) ;
        unsigned long y=(unsigned long)(MWC>>1) ;
    #elif defined TKISS
        unsigned long x=(unsigned long)(KISS>>1) ;
        unsigned long y=(unsigned long)(KISS>>1) ;
    #endif
#elif defined TFP32
    #define THEONE 1.0f
    #if defined TCONG
        float x=CONGfp ;
        float y=CONGfp ;
    #elif defined TSHR3
        float x=SHR3fp ;
        float y=SHR3fp ;
    #elif defined TMWC
        float x=MWCfp ;
        float y=MWCfp ;
    #elif defined TKISS
      float x=KISSfp ;
      float y=KISSfp ;
    #endif
#elif defined TFP64
    #define THEONE 1.0f
    #if defined TCONG
        double x=(double)CONGfp ;
        double y=(double)CONGfp ;
    #elif defined TSHR3
        double x=(double)SHR3fp ;
        double y=(double)SHR3fp ;
    #elif defined TMWC
        double x=(double)MWCfp ;
        double y=(double)MWCfp ;
    #elif defined TKISS
        double x=(double)KISSfp ;
        double y=(double)KISSfp ;
    #endif
#endif

      // Matching test
      unsigned long inside=((x*x+y*y) < THEONE) ? 1:0;
      total+=inside;
   }
  
   return(total);

}

int main(int argc, char *argv[]) {

  unsigned int seed_z=362436069,seed_w=52128862;
  LENGTH iterations=ITERATIONS,inside[8192],insides,part_inside,part_iterations;
  int numtasks,rank,rc,tag=1,i;
  float pi;

  char hostname[128];

  gethostname(hostname, sizeof hostname);

#ifdef TIME
  struct timeval start,end;
  long int useconds;
#endif

  MPI_Status Stat;
  MPI_Request RequestSend, RequestRecv, RequestSend2, RequestRecv2;
     
  rc = MPI_Init(&argc,&argv);
  if (rc != MPI_SUCCESS) {
    printf ("Error starting MPI program. Terminating.\n");
    MPI_Abort(MPI_COMM_WORLD, rc);
  }

  MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  
  MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  const int nitems=2;
  int blocklengths[2] = {1,1};

#ifdef LONG
  MPI_Datatype types[2] = {MPI_LONG, MPI_LONG};
#else
  MPI_Datatype types[2] = {MPI_INT, MPI_LONG};
#endif

  MPI_Datatype mpi_result_type;
  MPI_Aint     offsets[2];

  offsets[0] = offsetof(result, inside);
  offsets[1] = offsetof(result, useconds);
  
  MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_result_type);
  MPI_Type_commit(&mpi_result_type);

  if (rank==0) {
    
    if (argc > 1) {
      iterations=(LENGTH)atoll(argv[1]);
    }
    else {
      printf("\n\tPi : Estimate Pi with Monte Carlo exploration\n\n");
      printf("\t\t#1 : number of iterations (default 1 billion)\n\n");
    }
    
    printf ("\tSizeof int = %lld bytes.\n", (long long)sizeof(int));
    printf ("\tSizeof long = %lld bytes.\n", (long long)sizeof(long));
    printf ("\tSizeof long long = %lld bytes.\n", (long long)sizeof(long long));
    
    printf ("\tMax int = %u\n", INT_MAX);
    printf ("\tMax long = %ld\n", LONG_MAX);
    printf ("\tMax long long = %lld\n", LLONG_MAX);
    
    part_iterations=((iterations%numtasks) == 0) ? iterations/numtasks:iterations/numtasks+1 ;
    
    // Split part of code
    for (i=1;i<numtasks;i++) {
      
#ifdef LONG
      rc = MPI_Isend(&part_iterations, 1, MPI_LONG_LONG, i, tag, 
		     MPI_COMM_WORLD,&RequestSend);
#else
      rc = MPI_Isend(&part_iterations, 1, MPI_INT, i, tag, 
		     MPI_COMM_WORLD,&RequestSend);
#endif      
    }
    if (numtasks>1) {
      MPI_Wait(&RequestSend, &Stat);
    }
    
    gettimeofday(&start,(struct timezone *)0);
    
    insides=MainLoopGlobal(part_iterations,seed_w,seed_z);
    
    gettimeofday(&end,(struct timezone *)0);
    useconds=(end.tv_sec-start.tv_sec)*1000000+end.tv_usec-start.tv_usec;
    
    printf("\tOn %s with rank %i, find %lld inside in %lu useconds.\n",
	     hostname,rank,(long long)insides,useconds);
      
    // Join part of code
      for (i=1;i<numtasks;i++) {

	result recv;
	
	rc = MPI_Irecv(&recv, 1, mpi_result_type, i, tag, 
		      MPI_COMM_WORLD, &RequestRecv2);

	if (numtasks>1) {
	    MPI_Wait(&RequestRecv2, &Stat);
	}
	
	inside[i]=recv.inside;
	useconds=recv.useconds;
	
	printf("\tReceive from %i, find %lld inside in %lu useconds\n",i,(long long)inside[i],useconds);
	
	insides+=inside[i];
      }
      
      pi=4.*(float)insides/(float)(part_iterations*numtasks);
      
      printf("\n\tPi=%.40f\n\twith error %.40f\n\twith %lld iterations\n\n",pi,
	     fabs(pi-4*atan(1.))/pi,(long long)iterations);

  }
  else
    {
      // Receive information from master
      
#ifdef LONG
      rc = MPI_Irecv(&part_iterations, 1, MPI_LONG_LONG, 0, tag, 
                    MPI_COMM_WORLD, &RequestRecv);
#else
      rc = MPI_Irecv(&part_iterations, 1, MPI_INT, 0, tag, 
                    MPI_COMM_WORLD, &RequestRecv);
#endif
      MPI_Wait(&RequestRecv, &Stat);
      
      /* Not such a good idea to print information... 
	 printf("\tOn %s with rank %i, receive from master %lld\n",
	 hostname,rank,(long long)part_iterations); */
      
      gettimeofday(&start,(struct timezone *)0);

      part_inside=MainLoopGlobal(part_iterations,rotr(seed_w,rank),rotl(seed_z,rank));
      
      gettimeofday(&end,(struct timezone *)0);
      useconds=(end.tv_sec-start.tv_sec)*1000000+end.tv_usec-start.tv_usec;

      /*
      printf("\tOn %s with rank %i, find %lld inside in %lu useconds.\n",
	     hostname,rank,(long long)part_inside,useconds);
      */

      result send;
      send.inside=part_inside;
      send.useconds=useconds;      
      
      rc = MPI_Isend(&send, 1, mpi_result_type, 0, tag, MPI_COMM_WORLD,&RequestSend2);

      MPI_Wait(&RequestSend2, &Stat);

    }

  MPI_Type_free(&mpi_result_type);
  MPI_Finalize();
  
}
