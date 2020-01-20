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

#define NUM_THREADS	1024

struct thread_data
{
  int thread_id;
  unsigned int seed_w;
  unsigned int seed_z;
  unsigned long iterations;
  unsigned long inside;
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

int MainLoopGlobal(unsigned long iterations,unsigned int seed_w,unsigned int seed_z)
{
  unsigned int z=seed_z;
  unsigned int w=seed_w;
  
  unsigned long total=0;
  
  for (unsigned long i=0;i<iterations;i++) {
    
    float x=MWCfp ;
    float y=MWCfp ;
    
    // Matching test
    unsigned long inside=((x*x+y*y) < 1.0f) ? 1:0;
    total+=inside;
   }
  
  return(total);
  
}

void *MainLoopThread(void *threadarg)
{
  int taskid;
  unsigned long iterations,total=0;
  unsigned int z,w;
  
  struct thread_data *my_data;
  
  my_data=(struct thread_data *) threadarg;
  
  taskid = my_data->thread_id;
  iterations = my_data->iterations;
  z = my_data->seed_z;
  w = my_data->seed_w;
  
  printf("Hello World! It's me, thread #%ld, with seeds (%ld,%ld) and %ld !\n", taskid,z,w,iterations);
  
  for (unsigned long i=0;i<iterations;i++) {
    
    float x=MWCfp ;
    float y=MWCfp ;
    
    // Matching test
    unsigned long inside=((x*x+y*y) < 1.0f) ? 1:0;
    total+=inside;
  }
  
  my_data->inside=total;
  
  printf("Thread %ld done, found %ld inside.\n",taskid,total);
  
  pthread_exit((void*) my_data);
}

int main(int argc, char *argv[])
{
  pthread_t threads[NUM_THREADS];
  pthread_attr_t attr;
  int *taskids[NUM_THREADS];
  int rc, t, sum;
  unsigned long num_threads,iterations;
  sum=0;
  void *status;
  
  iterations=atol(argv[1]);
  num_threads=atoi(argv[2]);
  
  printf("Number of threads defined to %ld\n",num_threads);
  printf("Number of iterations defined to %ld\n",iterations);

  /* Initialize and set thread detached attribute */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for(t=0;t<num_threads;t++) {
    
    thread_data_array[t].thread_id = t;
    thread_data_array[t].iterations = iterations/num_threads;
    thread_data_array[t].seed_w = (t+1)<<4;
    thread_data_array[t].seed_z = (1048576*(t+1))>>4;
    thread_data_array[t].inside = 0;
    
    printf("Creating thread %d\n", t);
    rc = pthread_create(&threads[t], &attr, MainLoopThread, (void *) 
                        &thread_data_array[t]);

    if (rc) {
      printf("ERROR; return code from pthread_create() is %d\n", rc);
      exit(-1);
    }
  }
  
  /* Free attribute and wait for the other threads */
  pthread_attr_destroy(&attr);

  for(t=0; t<num_threads; t++) {
    rc = pthread_join(threads[t], &status);
    if (rc) {
      printf("ERROR; return code from pthread_join() is %d\n", rc);
      exit(-1);
    }
    printf("Main: completed join with thread %ld having a status of %ld\n",
           t,(long)status);
  }

  for(t=0;t<num_threads;t++) {
    printf("Return to main with %ld thread, %ld inside\n", 
           t,thread_data_array[t].inside);
  }

  printf("Main: program completed. Exiting.\n");
  pthread_exit(NULL);

}
