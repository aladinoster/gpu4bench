/*
//@HEADER
// ************************************************************************
// 
//                        Kokkos v. 2.0
//              Copyright (2014) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
// 
// ************************************************************************
//@HEADER
*/

#include <Kokkos_Core.hpp>
#include <cstdio>
#include <typeinfo>
#include <math.h>
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

// On Kokkos, vector, matrix and + are "View"s
typedef Kokkos::View<LENGTH*> view;

struct splitter {

  view Inside;
  unsigned int seed_w;
  unsigned int seed_z;
  LENGTH iterations;
  
  splitter(view Inside_,LENGTH iterations,unsigned int seed_w_,unsigned int seed_z_) :
    Inside (Inside_),iterations (iterations),seed_w (seed_w_),seed_z (seed_z_)
  {}
  
  KOKKOS_INLINE_FUNCTION
  void operator() (int i) const {

    // MainLoopGlobal totally copied inside operator()
#if defined TCONG
   unsigned int jcong=seed_z+i;
#elif defined TSHR3
   unsigned int jsr=seed_w+i;
#elif defined TMWC
   unsigned int z=seed_z+i;
   unsigned int w=seed_w-i;
#elif defined TKISS
   unsigned int jcong=seed_z+i;
   unsigned int jsr=seed_w-i;
   unsigned int z=seed_z+i;
   unsigned int w=seed_w-i;
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

	unsigned long inside=((x*x+y*y) < THEONE) ? 1:0;
	total+=inside;
   }
    Inside(i)=total;
  }
};

struct print {

  view Inside;

  print(view Inside_) :
    Inside (Inside_)
  {}
  
  KOKKOS_INLINE_FUNCTION
  void operator() (const int i) const {
    printf ("Inside of %i = %lld\n", i,Inside(i));
  }
};

// Reduction functor that reads the View given to its constructor.
struct ReduceFunctor {
  view Inside;

  ReduceFunctor (view Inside_) : Inside (Inside_) {}

  typedef LENGTH value_type;

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, LENGTH &lsum) const {
    lsum += Inside(i);
  }
};

int main (int argc, char* argv[]) {

  unsigned int seed_w=110271,seed_z=101008,ParallelRate=PARALLELRATE;
  LENGTH iterations=ITERATIONS,insides=0;
  struct timeval tv1,tv2;
  struct timezone tz;

  if (argc > 1) {
    iterations=(LENGTH)atoll(argv[1]);
    ParallelRate=atoi(argv[2]);
  }
  else {
    printf("\n\tPi : Estimate Pi with Monte Carlo exploration\n\n");
    printf("\t\t#1 : number of iterations (default 1 billion)\n");
    printf("\t\t#2 : number of ParallelRate (default 1024)\n\n");
  }

  printf ("\n\tInformation about architecture:\n\n");

  printf ("\tSizeof int = %lld bytes.\n", (long long)sizeof(int));
  printf ("\tSizeof long = %lld bytes.\n", (long long)sizeof(long));
  printf ("\tSizeof long long = %lld bytes.\n\n", (long long)sizeof(long long));

  printf ("\tMax int = %u\n", INT_MAX);
  printf ("\tMax long = %ld\n", LONG_MAX);
  printf ("\tMax long long = %lld\n\n", LLONG_MAX);
  
  Kokkos::initialize (argc, argv);

  printf ("Pi Dart Dash on Kokkos execution space %s\n",
          typeid (Kokkos::DefaultExecutionSpace).name ());

  view Inside("Inside",ParallelRate);
  LENGTH IterationsEach=((iterations%ParallelRate)==0)?iterations/ParallelRate:iterations/ParallelRate+1;

  gettimeofday(&tv1, &tz);


  // Core of Kokkos : parallel_for & parallel_reduce
  Kokkos::parallel_for (ParallelRate,splitter(Inside,IterationsEach,seed_w,seed_z));
  //  Kokkos::parallel_for (ParallelRate,print(Inside));
  Kokkos::parallel_reduce (ParallelRate, ReduceFunctor (Inside), insides);

  gettimeofday(&tv2, &tz);
	
  Kokkos::finalize ();

  double elapsed=(double)((tv2.tv_sec-tv1.tv_sec) * 1000000L +
			  (tv2.tv_usec-tv1.tv_usec))/1000000;  

  double itops=(double)(ParallelRate*IterationsEach)/elapsed;
  
  printf("\n");

  printf("Inside/Total %lld %lld\nParallelRate %i\nElapsed Time %.2f\nItops %.0f\nPi estimation %f\n\n",(long long)insides,(long long)ParallelRate*IterationsEach,ParallelRate,elapsed,itops,(4.*(float)insides/((float)(ParallelRate)*(float)(IterationsEach))));
  
}

