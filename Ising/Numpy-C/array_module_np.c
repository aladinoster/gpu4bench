/*  Example of change values in array in the Numpy-C-API. */
/* http://scipy-lectures.github.io/advanced/interfacing_with_c/interfacing_with_c.html */
/* http://dsnra.jpl.nasa.gov/software/Python/numpydoc/numpy-13.html */
/* http://enzo.googlecode.com/.../Grid_ConvertToNumpy.C */
/* http://docs.scipy.org/doc/numpy/reference/c-api.array.html#data-access */

/* 
   What is absolutely necessary to know:
   - 
*/

#include <Python.h>
#include <numpy/arrayobject.h>
#include <math.h>
#include <stdio.h>
#include <sys/time.h>

/* All RNG numbers of Marsaglia are Unsigned int32 ! */
#define znew  ((z=36969*(z&65535)+(z>>16))<<16)
#define wnew  ((w=18000*(w&65535)+(w>>16))&65535)
#define MWC   (znew+wnew)
#define SHR3  (jsr=(jsr=(jsr=jsr^(jsr<<17))^(jsr>>13))^(jsr<<5))
#define CONG  (jcong=69069*jcong+1234567)
#define KISS  ((MWC^CONG)+SHR3)
#define LFIB4 (t[c]=t[c]+t[c+58]+t[c+119]+t[++c+178])
#define SWB   (t[c+237]=(x=t[c+15])-(y=t[++c]+(x<y)))
#define UNI   (KISS*2.328306e-10)
#define VNI   ((long) KISS)*4.656613e-10

#define MWCfp (float)MWC * 2.328306435454494e-10f
#define KISSfp (float)KISS * 2.328306435454494e-10f
#define SHR3fp (float)SHR3 * 2.328306435454494e-10f

// static PyObject* array_metropolis_np(PyObject* self, PyObject *args)
static PyObject* array_metropolis_np(PyObject* self, PyObject *args)
{

  PyObject *in_array;
  PyArrayObject *local;
  unsigned int z,w;
  long iterations=1;
  int seed_w=521288629,seed_z=362436069,jcong=380116160,jsr=123456789;
  int factor;
  float B=0,J=0,T=0;
  double elapsed=0;
  unsigned long i;
  struct timeval tv1,tv2;
  struct timezone tz;

  /*  parse single numpy array argument */
  /* http://docs.python.org/release/1.5.2p2/ext/parseTuple.html */
  /* O for in_array */
  /* fff  for 3 float values (Coupling, Magnetic Field, Temperature)*/
  /* l for 1 long int value (iterations) */
  /* ii for 2 int values (2 seeds) */
  if (!PyArg_ParseTuple(args, "Offflii", &in_array,
			&J,&B,&T,&iterations,&seed_w,&seed_z))
    {
      printf("Not null argument provided !\n");
      return NULL;
    }
  
  /* display the extra parameters */
  printf("Parameters for physics: J=%f\nB=%f\nT=%f\n",J,B,T);
  printf("Parameters for simulation: iterations=%ld\nseed_w=%i\nseed_z=%i\n",
	 iterations,seed_w,seed_z);
  
  local = (PyArrayObject *) PyArray_ContiguousFromAny(in_array, NPY_INT32, 2, 2);
  

  if (PyArray_NDIM(local)!=2) {
    printf("This simulation is only for 2D arrays !\n");
    return NULL;
  }
  
  npy_intp size_x=*(npy_intp *)(PyArray_DIMS(local));
  npy_intp size_y=*(npy_intp *)(PyArray_DIMS(local)+1);
  /* Cast necessary, PyArg_ParseTuple() does not support any unsigned type ! */
  w=(unsigned int)seed_w;
  z=(unsigned int)seed_z;
 
  gettimeofday(&tv1, &tz);
  printf("Simulation started: ...");

  for (i=0;i<(unsigned long)iterations;i++) {
    
    npy_intp x=MWC%size_x ;
    npy_intp y=MWC%size_y ;

    int p=*(int *)PyArray_GETPTR2(local,x,y);
    
    int d=*(int *)PyArray_GETPTR2(local,x,(y+1)%size_y);
    int u=*(int *)PyArray_GETPTR2(local,x,(y-1)%size_y);
    int l=*(int *)PyArray_GETPTR2(local,(x-1)%size_x,y);
    int r=*(int *)PyArray_GETPTR2(local,(x+1)%size_x,y);
    
    float DeltaE=J*(float)p*(2.*(float)(u+d+l+r)+B);
    
    factor=((DeltaE < 0.0f) || (MWCfp < exp(-DeltaE/T))) ? -1:1;

    *(int *)PyArray_GETPTR2(local,x,y)=(factor*p);
  }
  gettimeofday(&tv2, &tz);    
  elapsed=(double)((tv2.tv_sec-tv1.tv_sec) * 1000000L +	\
                   (tv2.tv_usec-tv1.tv_usec))/1000000.; 
  printf("done in %.2f seconds !\n",elapsed);
  
  return PyArray_Return(local);
  
}

static PyObject* array_display_np(PyObject* self, PyObject *args)
{
  
  PyObject *in_array;
  PyArrayObject      *local;
  unsigned int i,x,y;

  /*  parse single numpy array argument */
  /* http://docs.python.org/release/1.5.2p2/ext/parseTuple.html */
  /* O for in_array */
  if (!PyArg_ParseTuple(args, "O", &in_array)) {
    printf("Not null argument provided !\n");
    return NULL;
  }
  
  local = (PyArrayObject *) PyArray_ContiguousFromObject(in_array, PyArray_NOTYPE, 2, 2);
  
  /* display the size of array */
  printf("Array has %i dimension(s) : size of ",PyArray_NDIM(local));
  if (PyArray_NDIM(local)==1) {
      printf("%i\n",(int)*(PyArray_DIMS(local)));
  }
  else {
    for (i=0;i<PyArray_NDIM(local)-1;i++) {
      printf("%ix",(int)*(PyArray_DIMS(local)+i));
    }
    printf("%i\n",(int)*(PyArray_DIMS(local)+PyArray_NDIM(local)-1));
  }
  
  /* display the array */
  printf("[");
  for (x=0;x<(unsigned int)(int)*(PyArray_DIMS(local));x++) {
    printf("[");
    for (y=0;y<(unsigned int)(int)*(PyArray_DIMS(local)+1);y++)
      printf("%2i ",* (int *)PyArray_GETPTR2(local,x,y));
    printf("]\n ");
  }
  printf("]\n");
  
  // Py_INCREF(local);

  /*  construct the output array, like the input array */
  return PyArray_Return(local);

}


/*  define functions in module */
static PyMethodDef ArrayMethods[] =
{
     {"array_metropolis_np", array_metropolis_np, METH_VARARGS,
         "evaluate a metropolis simulation on a numpy array"},
     {"array_display_np", array_display_np, METH_VARARGS,
         "evaluate a metropolis simulation on a numpy array"},
     {NULL, NULL, 0, NULL}
};

/* module initialization */
PyMODINIT_FUNC

initarray_module_np(void)
{
     (void) Py_InitModule("array_module_np", ArrayMethods);
     /* IMPORTANT: this must be called */
     import_array();
}
