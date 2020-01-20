//
// Estimation of Pi using Monte Carlo exploration process
// gcc -std=c99 -O3 -o Pi Pi.c -lm 
//

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

// Marsaglia RNG very simple implementation
#define znew  ((z=36969*(z&65535)+(z>>16))<<16)
#define wnew  ((w=18000*(w&65535)+(w>>16))&65535)
#define MWC   (znew+wnew)
#define SHR3  (jsr=(jsr=(jsr=jsr^(jsr<<17))^(jsr>>13))^(jsr<<5))
#define CONG  (jcong=69069*jcong+1234567)
#define KISS  ((MWC^CONG)+SHR3)

#define CONGfp CONG * 2.328306435454494e-10f
#define MWCfp MWC * 2.328306435454494e-10f
#define KISSfp KISS * 2.328306435454494e-10f

#define ITERATIONS 1000000000

#ifdef LONG
#define LENGTH long long
#else
#define LENGTH int
#endif

LENGTH MainLoopGlobal(LENGTH iterations,unsigned int seed_w,unsigned int seed_z)
{
   unsigned int z=seed_z;
   unsigned int w=seed_w;

   LENGTH total=0;

   for (LENGTH i=0;i<iterations;i++) {

      float x=MWCfp ;
      float y=MWCfp ;

      // Matching test

#ifdef IFSQRT
      if ( sqrt(x*x+y*y) < 1.0f ) {
	total+=1;
      }
#elif IFWOSQRT
      if ( (x*x+y*y) < 1.0f ) {
	total+=1;
      }
#else
      int inside=((x*x+y*y) < 1.0f) ? 1:0;
      total+=inside;
#endif

   }

   return(total);

}

int main(int argc, char *argv[]) {

  unsigned int seed_w=10,seed_z=10;
  LENGTH iterations=ITERATIONS;

  if (argc > 1) {
    iterations=(LENGTH)atoll(argv[1]);
  }
  else {
    printf("\n\tPi : Estimate Pi with Monte Carlo exploration\n\n\t\t#1 : number of iterations (default 1 billion)\n\n");
  }

  printf ("\n\tInformation about architecture:\n\n");

  printf ("\tSizeof int = %lld bytes.\n", (long long)sizeof(int));
  printf ("\tSizeof long = %lld bytes.\n", (long long)sizeof(long));
  printf ("\tSizeof long long = %lld bytes.\n\n", (long long)sizeof(long long));

  printf ("\tMax int = %u\n", INT_MAX);
  printf ("\tMax long = %ld\n", LONG_MAX);
  printf ("\tMax long long = %lld\n\n", LLONG_MAX);

  float pi=(float)MainLoopGlobal(iterations,seed_w,seed_z)/(float)iterations*4;

  printf("\tPi=%.40f\n\twith error %.40f\n\twith %lld iterations\n\n",pi,
         fabs(pi-4*atan(1))/pi,(long long)iterations);
  
}
