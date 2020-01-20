//
// Estimation of Pi using Monte Carlo exploration process
// Cecill v2 Emmanuel QUEMENER <emmanuel.quemener@gmail.com>
// Exploit OpenMP on Xeon Phi
// source /opt/intel/bin/compilervars.sh intel64
// icpc -std=c99 -O3 -o Pi Pi.c -lm 
//

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <limits.h>
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

#define PARALLELRATE 1024

#ifdef LONG
#define LENGTH long long
#else
#define LENGTH int
#endif

#pragma omp declare target
LENGTH splitter(LENGTH,int,int,int);

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

LENGTH splitter(LENGTH iterations,int seed_w,int seed_z,int ParallelRate) {

  LENGTH *inside,insides=0;
  int i;
  struct timeval tv1,tv2;
  struct timezone tz;
  LENGTH IterationsEach=((iterations%ParallelRate)==0)?iterations/ParallelRate:iterations/ParallelRate+1;

  inside=(LENGTH*)malloc(sizeof(LENGTH)*ParallelRate);
  
  gettimeofday(&tv1, &tz);
  
#pragma omp target device(0)
#pragma omp teams num_teams(60) thread_limit(4)
// #pragma omp parallel for
#pragma omp distribute
  for (int i=0 ; i<ParallelRate; i++) {
    inside[i]=MainLoopGlobal(IterationsEach,seed_w+i,seed_z+i);
  }

  for (int i=0 ; i<ParallelRate; i++) {
    insides+=inside[i];
  }

  gettimeofday(&tv2, &tz);

  for (int i=0 ; i<ParallelRate; i++) {
    printf("\tFound %lld for process %i\n",(long long)inside[i],i);
  }
  printf("\n");

  double elapsed=(double)((tv2.tv_sec-tv1.tv_sec) * 1000000L +
			  (tv2.tv_usec-tv1.tv_usec))/1000000;
  
  double itops=(double)(ParallelRate*IterationsEach)/elapsed;

  printf("ParallelRate %i\nElapsed Time %.2f\nItops %.0f\n",ParallelRate,elapsed,itops);

  free(inside);
  return(insides);
}
 
int main(int argc, char *argv[]) {

  unsigned int seed_w=110271,seed_z=101008,ParallelRate=PARALLELRATE;
  LENGTH iterations=ITERATIONS;
  LENGTH insides=0;

  if (argc > 1) {
    iterations=(LENGTH)atoll(argv[1]);
    ParallelRate=atoi(argv[2]);
  }
  else {
    printf("\n\tPi : Estimate Pi with Monte Carlo exploration\n\n");
    printf("\t\t#1 : number of iterations (default 1 billion)\n");
    printf("\t\t#2 : ParallelRate (default 1024)\n\n");
  }

  printf ("\n\tInformation about architecture:\n\n");

  printf ("\tSizeof int = %lld bytes.\n", (long long)sizeof(int));
  printf ("\tSizeof long = %lld bytes.\n", (long long)sizeof(long));
  printf ("\tSizeof long long = %lld bytes.\n\n", (long long)sizeof(long long));

  printf ("\tMax int = %u\n", INT_MAX);
  printf ("\tMax long = %ld\n", LONG_MAX);
  printf ("\tMax long long = %lld\n\n", LLONG_MAX);

  insides=splitter(iterations,seed_w,seed_z,ParallelRate);

  LENGTH total=((iterations%ParallelRate)==0)?iterations:(iterations/ParallelRate+1)*ParallelRate;

  printf("Inside/Total %ld %ld\nPi estimation %f\n\n",insides,total,(4.*(float)insides/total));
  
}
