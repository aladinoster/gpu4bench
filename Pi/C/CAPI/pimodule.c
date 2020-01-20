#include <Python.h>
#include <math.h>

// Marsaglia RNG very simple implementation
#define znew  ((z=36969*(z&65535)+(z>>16))<<16)
#define wnew  ((w=18000*(w&65535)+(w>>16))&65535)
#define MWC   (znew+wnew)
#define SHR3  (jsr=(jsr=(jsr=jsr^(jsr<<17))^(jsr>>13))^(jsr<<5))
#define CONG  (jcong=69069*jcong+1234567)
#define KISS  ((MWC^CONG)+SHR3)

#define MWCfp MWC * 2.328306435454494e-10f
#define KISSfp KISS * 2.328306435454494e-10f

int _InsideCircle(unsigned int iterations,
		  unsigned int seed_w,unsigned int seed_z)
{
   unsigned int z=seed_z;
   unsigned int w=seed_w;
   unsigned int i;

   int total=0;

   for (i=0;i<iterations;i++) {

      float x=MWCfp ;
      float y=MWCfp ;

      // Matching test
      int inside=((x*x+y*y) < 1.0f) ? 1:0;
      total+=inside;
   }

   return(total);
}

static PyObject* InsideCircle(PyObject* self, PyObject* args)
{
  int i,z,w;
 
  if (!PyArg_ParseTuple(args, "iii", &i,&z,&w))
    return NULL;
 
  return Py_BuildValue("i", _InsideCircle(i,z,w));
}
 
static PyMethodDef MonteCarloMethods[] = {
  {"InsideCircle", InsideCircle, METH_VARARGS, "Inside Circle Randomly Check"},
  {NULL, NULL, 0, NULL}
};
 
PyMODINIT_FUNC
initMonteCarlo(void)
{
  (void) Py_InitModule("MonteCarlo", MonteCarloMethods);
}
