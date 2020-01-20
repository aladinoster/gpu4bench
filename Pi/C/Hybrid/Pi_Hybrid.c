//
// Estimation of Pi using Monte Carlo exploration process
// gcc -std=c99 -O3 -o Pi Pi.c -lm 
// Emmanuel Quemener <emmanuel.quemener@ens-lyon.fr>
// Cecill v2

// Needed for gethostname
#define _BSD_SOURCE
#include <sys/unistd.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <mpi.h>
#include <stddef.h>

#ifdef TIME
#include <sys/time.h>
#endif

// Marsaglia RNG very simple implementation
#define znew  ((z=36969*(z&65535)+(z>>16))<<16)
#define wnew  ((w=18000*(w&65535)+(w>>16))&65535)
#define MWC   (znew+wnew)
#define SHR3  (jsr=(jsr=(jsr=jsr^(jsr<<17))^(jsr>>13))^(jsr<<5))
#define CONG  (jcong=69069*jcong+1234567)
#define KISS  ((MWC^CONG)+SHR3)

#define ITERATIONS 1000000000

#define MWCfp MWC * 2.328306435454494e-10f
#define KISSfp KISS * 2.328306435454494e-10f
#define SHR3fp SHR3 * 2.328306435454494e-10f
#define CONGfp CONG * 2.328306435454494e-10f

#define PROCESS 1

#ifdef LONG
#define LENGTH long long
#else
#define LENGTH int
#endif

typedef struct compute_node {
  LENGTH iterations;
  int process;
} node;

typedef struct compute_result {
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
  
  unsigned int seed_z=362436069,seed_w=52128862,process=PROCESS;
  // Number of NP or OpenMP processes <1024
  LENGTH iterations=ITERATIONS,insideMPI[8192],insideOpenMP[1024],
    part_inside=0,part_iterations,insides=0;
  int numtasks,rank,rc,tag=1,i;
  float pi;
  
  // Hostname supposed to be <128 characters
  char hostname[128];
  
  gethostname(hostname, sizeof hostname);
  
  struct timeval start,end;
  long int useconds;
  
  MPI_Status Stat;
  
  rc = MPI_Init(&argc,&argv);
  if (rc != MPI_SUCCESS) {
    printf ("Error starting MPI program. Terminating.\n");
    MPI_Abort(MPI_COMM_WORLD, rc);
  }

  MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  
  const int nitems=2;
  int blocklengths[2] = {1,1};
  
#ifdef LONG
  MPI_Datatype types_node[2] = {MPI_LONG, MPI_INT};
  MPI_Datatype types_result[2] = {MPI_LONG, MPI_LONG};
#else
  MPI_Datatype types_node[2] = {MPI_INT, MPI_INT};
  MPI_Datatype types_result[2] = {MPI_INT, MPI_LONG};
#endif
  
  MPI_Datatype mpi_node_type,mpi_result_type;
  MPI_Aint     offsets[2],offsetsr[2];
  
  offsets[0] = offsetof(node, iterations);
  offsets[1] = offsetof(node, process);
  
  MPI_Type_create_struct(nitems, blocklengths, offsets, types_node, &mpi_node_type);
  MPI_Type_commit(&mpi_node_type);
  
  offsetsr[0] = offsetof(result, inside);
  offsetsr[1] = offsetof(result, useconds);
  
  MPI_Type_create_struct(nitems, blocklengths, offsetsr, types_result, &mpi_result_type);
  MPI_Type_commit(&mpi_result_type);
  
  if (rank==0) {
    
    if (argc > 1) {
      iterations=(LENGTH)atoll(argv[1]);
      process=atoi(argv[2]);
    }
    else {
      printf("\n\tPi : Estimate Pi with Monte Carlo exploration\n\n");
      printf("\t\t#1 : number of iterations (default 1 billion)\n");
      printf("\t\t#2 : number of OpenMP processes (default 1)\n\n");
    }
    
    printf ("\n\tInformation about architecture:\n\n");
    
    printf ("\tSizeof int = %lld bytes.\n", (long long)sizeof(int));
    printf ("\tSizeof long = %lld bytes.\n", (long long)sizeof(long));
    printf ("\tSizeof long long = %lld bytes.\n", (long long)sizeof(long long));
    
    printf ("\tMax int = %u\n", INT_MAX);
    printf ("\tMax long = %ld\n", LONG_MAX);
    printf ("\tMax long long = %lld\n\n", LLONG_MAX);
    
    part_iterations=(((iterations%numtasks)%process) == 0) ? iterations/numtasks/process:iterations/numtasks/process+1 ;
    
    node send;
    send.iterations=part_iterations;
    send.process=process;

    // Split part of code
    for (i=1;i<numtasks;i++) {      
      rc = MPI_Send(&send, 1, mpi_node_type, i, tag, MPI_COMM_WORLD);
    }
    
    gettimeofday(&start,(struct timezone *)0);
    
#pragma omp parallel for
    for (int i=0 ; i<process; i++) {
      insideOpenMP[i]=MainLoopGlobal(part_iterations,
				     rotr(seed_w,i),
				     rotl(seed_z,i));
      /*
	printf("\t(%s,%i) found %lld for process %i\n",hostname,0,
	(long long)insideOpenMP[i],i); */
    }
    /*
      printf("\n");
    */
    
    insides=0;
    for (int i=0 ; i<process; i++) {
      insides+=insideOpenMP[i];
    }
    
    gettimeofday(&end,(struct timezone *)0);
    useconds=(end.tv_sec-start.tv_sec)*1000000+end.tv_usec-start.tv_usec;
    
    printf("\tOn %s with rank #%i find %lld inside in %lu useconds.\n",
	   hostname,rank,(long long)insides,useconds);
    
    // Join part of code
    for (i=1;i<numtasks;i++) {
      
      result recv;
      
      rc = MPI_Recv(&recv, 1, mpi_result_type, i, tag, MPI_COMM_WORLD,&Stat);
      
      insideMPI[i]=recv.inside;
      useconds=recv.useconds;
      
      printf("\tReceive from rank #%i, find %lld inside in %lu useconds\n",i,(long long)insideMPI[i],useconds);
    
      insides+=insideMPI[i];
    }
    
    pi=4.*(float)insides/(float)(part_iterations*numtasks*process);
    
    printf("\n\tPi=%.40f\n\twith error %.40f\n\twith %lld iterations\n\n",pi,
	   fabs(pi-4*atan(1.))/pi,(long long)(part_iterations*numtasks*process));
    
  }
  else
    {
      // Receive information from master
      
      node recv;
      
      rc = MPI_Recv(&recv, 1, mpi_node_type, 0, tag, MPI_COMM_WORLD,&Stat);
      /*   
      printf("\t(%s,%i) receive from master %lld with %i process\n",
      hostname,rank,(long long)recv.iterations,recv.process);
      */
      
      gettimeofday(&start,(struct timezone *)0);
      
#pragma omp parallel for
      for (int i=0 ; i<recv.process; i++) {
	insideOpenMP[i]=MainLoopGlobal(recv.iterations,rotr(seed_w,rank+i),rotl(seed_z,rank-i));
	/*
	  printf("\t(%s,%i) found %lld for process %i\n",hostname,rank,
	  (long long)insideOpenMP[i],i);
	*/
      }
      
      /* printf("\n"); */
      
      for (int i=0 ; i<recv.process; i++) {
	part_inside+=insideOpenMP[i];
      }
      
      gettimeofday(&end,(struct timezone *)0);
      useconds=(end.tv_sec-start.tv_sec)*1000000+end.tv_usec-start.tv_usec;
      /*
	printf("\tOn %s rank %i find %lld inside in %lu useconds.\n",
	hostname,rank,(long long)part_inside,useconds);
      */
      result send;
      send.inside=part_inside;
      send.useconds=useconds;
      
      rc = MPI_Send(&send, 1, mpi_result_type, 0, tag, MPI_COMM_WORLD);
      
    }
  
  MPI_Type_free(&mpi_node_type);
  MPI_Type_free(&mpi_result_type);
  
  MPI_Finalize();
  
}
