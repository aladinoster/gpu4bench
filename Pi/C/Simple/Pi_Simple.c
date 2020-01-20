//
// Estimation of Pi using Monte Carlo exploration process
// Cecill v2 Emmanuel QUEMENER <emmanuel.quemener@gmail.com>
// gcc -O3 -o Pi_Simple Pi_Simple.c -lm 
//

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <sys/time.h>

// Marsaglia RNG very simple implementation
#define znew  ((z=36969*(z&65535)+(z>>16))<<16)
#define wnew  ((w=18000*(w&65535)+(w>>16))&65535)
#define MWC   (znew+wnew)
#define SHR3  (jsr=(jsr=(jsr=jsr^(jsr<<17))^(jsr>>13))^(jsr<<5))
#define CONG  (jcong=69069*jcong+1234567)
#define KISS  ((MWC^CONG)+SHR3)

#define MWCfp MWC*2.328306435454494e-10f
#define KISSfp KISS*2.328306435454494e-10f
#define SHR3fp SHR3*2.328306435454494e-10f
#define CONGfp CONG*2.328306435454494e-10f

#define ITERATIONS 1000000000

#ifdef LONG
#define LENGTH long long
#else
#define LENGTH int
#endif

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

   LENGTH total=0,i;
   unsigned long inside;
     
   for (i=0;i<iterations;i++) {

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
#define THEONE (float)1.0f
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
#define THEONE (double)1.0f
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

      inside=((x*x+y*y) < THEONE) ? 1:0;
      total+=inside;
   }

   return(total);
}

int main(int argc, char *argv[]) {

  unsigned int seed_w=110271,seed_z=101008;
  LENGTH iterations=ITERATIONS,inside=0;
  struct timeval tv1,tv2;
  struct timezone tz;

  if (argc > 1) {
    iterations=(LENGTH)atoll(argv[1]);
  }
  else {
    printf("\n\tPi : Estimate Pi with Monte Carlo exploration\n\n");
    printf("\t\t#1 : number of iterations (default 1 billion)\n");
  }

  printf ("\n\tInformation about architecture:\n\n");

  printf ("\tSizeof int = %lld bytes.\n", (long long)sizeof(int));
  printf ("\tSizeof long = %lld bytes.\n", (long long)sizeof(long));
  printf ("\tSizeof long long = %lld bytes.\n\n", (long long)sizeof(long long));

  printf ("\tMax int = %u\n", INT_MAX);
  printf ("\tMax long = %ld\n", LONG_MAX);
  printf ("\tMax long long = %lld\n\n", LLONG_MAX);

  // Set start timer
  gettimeofday(&tv1, &tz);

  inside=MainLoopGlobal(iterations,seed_w,seed_z);

  // Set stop timer
  gettimeofday(&tv2, &tz);
  
  double elapsed=(double)((tv2.tv_sec-tv1.tv_sec) * 1000000L +
			  (tv2.tv_usec-tv1.tv_usec))/1000000;  
  
  double itops=(double)(iterations)/elapsed;
  
  printf("Inside/Total %ld %ld\nElapsed Time %.2f\nItops %.0f\nPi estimation %f\n\n",inside,iterations,elapsed,itops,(4.*(float)inside/(float)(iterations)));

  printf("\n\n");

  return 0;
  
}
