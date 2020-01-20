//
// Estimation of Pi using Monte Carlo exploration process
// gcc -std=c99 -O3 -o Pi Pi.c -lm 
//
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// Marsaglia RNG very simple implementation
#define znew  ((z=36969*(z&65535)+(z>>16))<<16)
#define wnew  ((w=18000*(w&65535)+(w>>16))&65535)
#define MWC   (znew+wnew)
#define SHR3  (jsr=(jsr=(jsr=jsr^(jsr<<17))^(jsr>>13))^(jsr<<5))
#define CONG  (jcong=69069*jcong+1234567)
#define KISS  ((MWC^CONG)+SHR3)

#define MWCfp MWC * 2.328306435454494e-10f
#define KISSfp KISS * 2.328306435454494e-10f

#define ITERATIONS 1000000000

unsigned int MainLoopGlobal(unsigned int iterations,unsigned int seed_w,unsigned int seed_z)
{
   unsigned int z=seed_z;
   unsigned int w=seed_w;

   unsigned int total=0;

   for (unsigned int i=0;i<iterations;i++) {

      float x=MWCfp ;
      float y=MWCfp ;

      // Matching test
      int inside=((x*x+y*y) < 1.0f) ? 1:0;
      total+=inside;
   }

   return(total);

}

int main(int argc, char *argv[]) {

  unsigned int seed_w=10,seed_z=10;
  unsigned int iterations=ITERATIONS;

  if (argc > 1) {
    iterations=(unsigned int)atol(argv[1]);
  }
  else {
    printf("\n\tPi : Estimate Pi with Monte Carlo exploration\n\n\t\t#1 : number of iterations (default 1 billion)\n");
  }

  float pi=(float)MainLoopGlobal(iterations,seed_w,seed_z)/(float)iterations*4;

  printf("Pi=%f with %u iterations\n",pi,iterations);
  
}
