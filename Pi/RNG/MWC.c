#include <stdio.h>
#include <math.h>

#define znew (z=36969*(z&65535)+(z>>16))
#define wnew (w=18000*(w&65535)+(w>>16))
#define MWC ((znew<<16)+wnew)
#define MWCfp (MWC + 1.0f) * 2.328306435454494e-10f

int main(void) {
  int i,z,w;

  z=37;
  w=91;

  for (i=1;i<100;i++) {

    printf("%i %i %i %i\n",i,znew,wnew,MWC);
  
  }

}
