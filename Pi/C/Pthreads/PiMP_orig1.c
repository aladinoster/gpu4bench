/******************************************************************************
* FILE: hello.c
* DESCRIPTION:
*   A "hello world" Pthreads program.  Demonstrates thread creation and
*   termination.
* AUTHOR: Blaise Barney
* LAST REVISED: 08/09/11
******************************************************************************/
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#define NUM_THREADS	4

#include <math.h>
#include <stdio.h>

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

int MainLoopGlobal(unsigned int iterations,unsigned int seed_w,unsigned int seed_z)
{
   unsigned int z=seed_z;
   unsigned int w=seed_w;

   int total=0;

   for (unsigned int i=0;i<iterations;i++) {

      float x=MWCfp ;
      float y=MWCfp ;

      // Matching test
      int inside=((x*x+y*y) < 1.0f) ? 1:0;
      total+=inside;
   }

   return(total);

}

void *PrintHello(void *threadid)
{
   long tid;
   tid = (long)threadid;
   printf("Hello World! It's me, thread #%ld!\n", tid);

   pthread_exit(NULL);
}

void *MainLoopThread(void *threadid)
{
   long tid;
   tid = (long)threadid;
   printf("Hello World! It's me, thread #%ld!\n", tid);

   /* unsigned int z=seed_z; */
   /* unsigned int w=seed_w; */
   unsigned int z=16;
   unsigned int w=64;

   int total=0;

   for (unsigned int i=0;i<ITERATIONS;i++) {

      float x=MWCfp ;
      float y=MWCfp ;

      // Matching test
      int inside=((x*x+y*y) < 1.0f) ? 1:0;
      total+=inside;
   }

   printf("In #%ld, found %i inside\n",tid,total);
   
   pthread_exit(NULL);
}

int main(int argc, char *argv[])
{
   pthread_t threads[NUM_THREADS];
   int rc;
   long t;
   for(t=0;t<NUM_THREADS;t++){
     printf("In main: creating thread %ld\n", t);
     rc = pthread_create(&threads[t], NULL, MainLoopThread, (void *)t);
     if (rc){
       printf("ERROR; return code from pthread_create() is %d\n", rc);
       exit(-1);
       }
     }

   /* Last thing that main() should do */
   pthread_exit(NULL);
}
