/******************************************************************************
* FILE: Pi_Threads
* Pthreads based on Hello from Blaise Barney, LLNL Tutorial
* DESCRIPTION:
*   A Pi by Monte Carlo Pthreads program.  Demonstrates thread creation and
*   termination.
* AUTHOR: Emmanuel Quemener from Blaise Barney
* LAST REVISED: 2013-03-30
******************************************************************************/

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <limits.h>

#define NUM_THREADS	1024

#define ITERATIONS 1000000000

#ifdef LONG
#define LENGTH long long
#else
#define LENGTH int
#endif

struct thread_data
{
  int thread_id;
  unsigned int seed_w;
  unsigned int seed_z;
  LENGTH iterations;
  LENGTH inside;
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

LENGTH MainLoopGlobal(LENGTH iterations,unsigned int seed_w,unsigned int seed_z)
{
  unsigned int z=seed_z;
  unsigned int w=seed_w;
  
  LENGTH total=0;
  
  for (LENGTH i=0;i<iterations;i++) {
    
    float x=MWCfp ;
    float y=MWCfp ;
    
    // Matching test
    LENGTH inside=((x*x+y*y) < 1.0f) ? 1:0;
    total+=inside;
   }
  
  return(total);
  
}

void *MainLoopThread(void *threadarg)
{
  int taskid;
  LENGTH iterations,total=0;
  unsigned int z,w;
  
  struct thread_data *my_data;
  
  my_data=(struct thread_data *) threadarg;
  
  taskid = my_data->thread_id;
  iterations = my_data->iterations;
  z = my_data->seed_z;
  w = my_data->seed_w;

  printf("\tThread #%i, with seeds (%u,%u) and %lld !\n",
         taskid,z,w,(long long)iterations);
  
  for (LENGTH i=0;i<iterations;i++) {
    
    float x=MWCfp ;
    float y=MWCfp ;
    
    // Matching test
    LENGTH inside=((x*x+y*y) < 1.0f) ? 1:0;
    total+=inside;
  }
  
  my_data->inside=total;
  
  pthread_exit((void*) my_data);
}

int main(int argc, char *argv[])
{
  pthread_t threads[NUM_THREADS];
  pthread_attr_t attr;
  int rc, t;
  unsigned int num_threads;
  LENGTH iterations,inside=0;
  void *status;
  float pi;
  
  if (argc > 1) {
    iterations=(LENGTH)atoll(argv[1]);
    num_threads=atoi(argv[2]);
  }
  else {
    iterations=ITERATIONS;
    num_threads=1;
    printf("\n\tPi : Estimate Pi with Monte Carlo exploration\n\n");
    printf("\t\t#1 : number of iterations (default 1 billion)\n");
    printf("\t\t#2 : number of threads (default 1)\n\n");
  }

  printf ("\n\tInformation about architecture:\n\n");

  printf ("\tSizeof int = %lld bytes.\n", (long long)sizeof(int));
  printf ("\tSizeof long = %lld bytes.\n", (long long)sizeof(long));
  printf ("\tSizeof long long = %lld bytes.\n\n", (long long)sizeof(long long));

  printf ("\tMax int = %u\n", INT_MAX);
  printf ("\tMax long = %ld\n", LONG_MAX);
  printf ("\tMax long long = %lld\n\n", LLONG_MAX);

  printf("\tNumber of threads defined to %u\n",num_threads);
  printf("\tNumber of iterations defined to %lld\n\n",(long long)iterations);

  /* Initialize and set thread detached attribute */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for(t=0;t<num_threads;t++) {
    
    thread_data_array[t].thread_id = t;
    thread_data_array[t].iterations = iterations/(long long)num_threads;
    thread_data_array[t].seed_w = (t+1)<<4;
    thread_data_array[t].seed_z = (1048576*(t+1))>>4;
    thread_data_array[t].inside = 0;
    
    printf("\tCreating thread %d\n", t);
    rc = pthread_create(&threads[t], &attr, MainLoopThread, (void *) 
                        &thread_data_array[t]);

    if (rc) {
      printf("\tERROR; return code from pthread_create() is %d\n", rc);
      exit(-1);
    }
  }
  
  /* Free attribute and wait for the other threads */
  pthread_attr_destroy(&attr);

  for(t=0; t<num_threads; t++) {
    rc = pthread_join(threads[t], &status);
    if (rc) {
      printf("\tERROR; return code from pthread_join() is %d\n", rc);
      exit(-1);
    }
    printf("\tMain: completed join with thread %i having a status of %ld\n",
           t,(long)status);
  }

  for(t=0;t<num_threads;t++) {
    printf("\tReturn to main with %i thread, %lld inside\n", 
           t,(long long)thread_data_array[t].inside);
    inside+=thread_data_array[t].inside;
  }

  printf("\tMain: program completed. Exiting.\n\n");

  pi=4.*(float)inside/(float)iterations;

  printf("\tPi=%f with error %f and %lld iterations\n\n",pi,
         fabs(pi-4*atan(1))/pi,(long long)iterations);

  pthread_exit(NULL);
}
