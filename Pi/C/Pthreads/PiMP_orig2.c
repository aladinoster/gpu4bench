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
#include <math.h>
#include <stdio.h>

#define NUM_THREADS	4

struct thread_data
{
  int thread_id;
  unsigned int seed_w;
  unsigned int seed_z;
  unsigned int iterations;
  unsigned int inside;
};

struct thread_data thread_data_array[NUM_THREADS];

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

void *MainLoopThread(void *threadarg)
{
  int taskid;
  unsigned int iterations,z,w,total=0;

  struct thread_data *my_data;

  my_data=(struct thread_data *) threadarg;

  taskid = my_data->thread_id;
  iterations = my_data->iterations;
  z = my_data->seed_z;
  w = my_data->seed_w;
  iterations = my_data->iterations;

  printf("Hello World! It's me, thread #%ld, with seeds (%i,%i) and %i !\n", taskid,z,w,iterations);

   for (unsigned int i=0;i<iterations;i++) {

      float x=MWCfp ;
      float y=MWCfp ;

      // Matching test
      int inside=((x*x+y*y) < 1.0f) ? 1:0;
      total+=inside;
   }

   printf("In #%ld, found %i inside\n",taskid,total);
   
   my_data->inside=total;

   pthread_exit(NULL);
}

int main(int argc, char *argv[])
{

pthread_t threads[NUM_THREADS];
int *taskids[NUM_THREADS];
int rc, t, sum;

sum=0;

for(t=0;t<NUM_THREADS;t++) {

  thread_data_array[t].thread_id = t;
  thread_data_array[t].iterations = ITERATIONS/NUM_THREADS;
  thread_data_array[t].seed_w = (t+1)<<4;
  thread_data_array[t].seed_z = (1048576*(t+1))>>4;
  thread_data_array[t].inside = 0;

  printf("Creating thread %d\n", t);
  rc = pthread_create(&threads[t], NULL, MainLoopThread, (void *) 
       &thread_data_array[t]);
  printf("%i inside\n", thread_data_array[t].inside);
  if (rc) {
    printf("ERROR; return code from pthread_create() is %d\n", rc);
    exit(-1);
    }
  }
pthread_exit(NULL);


}
