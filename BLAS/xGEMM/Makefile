SOURCE=xGEMM.c

CC=gcc
FC=gfortran
CFLAGS=-Wall -O3 -ffast-math
LDFLAGS=-lm
CUDASRC=/usr/share/doc/nvidia-cuda-doc/examples/
THUNKING=fortran_thunking.c
CUDASRCINC=fortran_common.h
CUDAINC=/usr/include
CUDALIB=/usr/lib/x86_64-linux-gnu/
CLBLASINC=/usr/include
CLBLASLIB=/usr/lib
# CLBLASINC=/opt/clBLAS/src
# CLBLASLIB=/opt/clBLAS/src/library/

PATCHTHUNKING=patch_thunking.h

GSLINC=/usr/include/gsl

GOTO2=/opt/GotoBLAS2

#OPENBLAS=/usr/lib/openblas-base
OPENBLAS=/usr/lib

ACML=/opt/acml
ACMLINC=$(ACML)/gfortran64_mp/include
ACMLLIB=$(ACML)/gfortran64_mp/lib

ATLAS=/usr/lib/atlas-base

#EXECUTABLE=cblas fblas gsl cublas thunking openblas acml clblas
#EXECUTABLE=cblas fblas gsl cublas thunking openblas clblas
#EXECUTABLE=gsl openblas clblas cublas thunking
#EXECUTABLE=openblas clblas cublas thunking
EXECUTABLE=acml gsl cblas fblas openblas clblas cublas thunking

FORMAT=DOUBLE
#FORMAT=FLOAT

#DIRECTIVES=-D$(FORMAT) -DPRINT -DUNIT
#DIRECTIVES=-D$(FORMAT) -DUNIT -DRESULTS -DQUIET
DIRECTIVES=-DUNIT -DQUIET

all: $(EXECUTABLE)

cblas: $(SOURCE)

	$(CC) -L$(ATLAS) $(CFLAGS) $(DIRECTIVES) -DFLOAT -DCBLAS $(LDFLAGS) \
		$(SOURCE) -latlas -lcblas -o $(SOURCE:.c=)_SP_$@

	$(CC) -L$(ATLAS) $(CFLAGS) $(DIRECTIVES) -DDOUBLE -DCBLAS $(LDFLAGS) \
		$(SOURCE) -latlas -lcblas -o $(SOURCE:.c=)_DP_$@

gotoblas: $(SOURCE)

	$(CC) $(CFLAGS) $(DIRECTIVES) -DFLOAT -DCBLAS $(LDFLAGS) \
		$(SOURCE) $(GOTO2)/libgoto2.a -lpthread -o $(SOURCE:.c=)_SP_$@

	$(CC) $(CFLAGS) $(DIRECTIVES) -DDOUBLE -DCBLAS $(LDFLAGS) \
		$(SOURCE) $(GOTO2)/libgoto2.a -lpthread -o $(SOURCE:.c=)_DP_$@

openblas: $(SOURCE)

	$(CC) $(CFLAGS) $(DIRECTIVES) -L/usr/include/openblas -fopenmp -DFLOAT -DCBLAS $(LDFLAGS) \
		$(SOURCE) $(OPENBLAS)/libopenblas.a -lpthread -o $(SOURCE:.c=)_SP_$@

	$(CC) $(CFLAGS) $(DIRECTIVES) -L/usr/include/openblas -fopenmp -DDOUBLE -DCBLAS $(LDFLAGS) \
		$(SOURCE) $(OPENBLAS)/libopenblas.a -lpthread -o $(SOURCE:.c=)_DP_$@

acml: $(SOURCE)

	$(CC) -I$(ACMLINC) $(CFLAGS) $(DIRECTIVES) -DFLOAT -DACML $(LDFLAGS) \
		$(SOURCE) -L$(ACMLLIB) -lacml_mp \
		-lgomp -lgfortran -lpthread -o $(SOURCE:.c=)_SP_$@

	$(CC) -I$(ACMLINC) $(CFLAGS) $(DIRECTIVES) -DDOUBLE -DACML $(LDFLAGS) \
		$(SOURCE) -L$(ACMLLIB) -lacml_mp \
		-lgomp -lgfortran -lpthread -o $(SOURCE:.c=)_DP_$@

fblas: $(SOURCE)

	$(CC) -L$(ATLAS) $(CFLAGS) $(DIRECTIVES) -DFLOAT -DFBLAS $(LDFLAGS) \
		$(SOURCE) -latlas -lf77blas -o $(SOURCE:.c=)_SP_$@

	$(CC) -L$(ATLAS) $(CFLAGS) $(DIRECTIVES) -DDOUBLE -DFBLAS $(LDFLAGS) \
		$(SOURCE) -latlas -lf77blas -o $(SOURCE:.c=)_DP_$@


gsl: $(SOURCE)

	$(CC) -I$(GSLINC) $(CFLAGS) $(DIRECTIVES) -DFLOAT \
		-DGSL $(LDFLAGS) \
		$(SOURCE) -lgslcblas -o $(SOURCE:.c=)_SP_$@

	$(CC) -I$(GSLINC) $(CFLAGS) $(DIRECTIVES) -DDOUBLE \
		-DGSL $(LDFLAGS) \
		$(SOURCE) -lgslcblas -o $(SOURCE:.c=)_DP_$@

clblas: $(SOURCE)

	$(CC) -I$(CLBLASINC) -L$(CLBLASLIB) $(CFLAGS) -DFLOAT \
		-DCLBLAS $(LDFLAGS) \
		$(DIRECTIVES) $(SOURCE) -lclBLAS -lOpenCL -o $(SOURCE:.c=)_SP_$@

	$(CC) -I$(CLBLASINC) -L$(CLBLASLIB) $(CFLAGS) -DDOUBLE \
		-DCLBLAS $(LDFLAGS) \
		$(DIRECTIVES) $(SOURCE) -lclBLAS -lOpenCL -o $(SOURCE:.c=)_DP_$@

cublas: $(SOURCE)

	$(CC) -I$(CUDAINC) -L$(CUDALIB) $(CFLAGS) -DFLOAT \
		-DCUBLAS $(LDFLAGS) \
		$(DIRECTIVES) $(SOURCE) -lcublas -o $(SOURCE:.c=)_SP_$@

	$(CC) -I$(CUDAINC) -L$(CUDALIB) $(CFLAGS) -DDOUBLE \
		-DCUBLAS $(LDFLAGS) \
		$(DIRECTIVES) $(SOURCE) -lcublas -o $(SOURCE:.c=)_DP_$@

thunking: $(SOURCE)

	# Copy of source for thunking CUBLAS approach
	cp $(CUDASRC)/$(THUNKING) $(CUDASRC)/$(THUNKING:.c=.h) \
		$(CUDASRC)/$(CUDASRCINC) .
	# Patch Thunking prototypes to compile on Debian Lenny
	$(CC) -DCUBLAS_GFORTRAN -I$(CUDAINC) $(CFLAGS) -c $(THUNKING)

	$(CC) -DCUBLAS_GFORTRAN \
		-I$(CUDAINC) -L$(CUDALIB) $(CFLAGS) -DFLOAT -DTHUNKING \
		$(LDFLAGS) $(DIRECTIVES) $(SOURCE) \
		$(THUNKING:.c=.o) -lcublas -o $(SOURCE:.c=)_SP_$@

	$(CC) -DCUBLAS_GFORTRAN \
		-I$(CUDAINC) -L$(CUDALIB) $(CFLAGS) -DDOUBLE -DTHUNKING \
		$(LDFLAGS) $(DIRECTIVES) $(SOURCE) \
		$(THUNKING:.c=.o) -lcublas -o $(SOURCE:.c=)_DP_$@

clean: $(SOURCE)
	find . -name "$(SOURCE:.c=)_*" -exec rm {} \;
	find . -name "$(THUNKING:.c=)*" -exec rm {} \;
	find . -name "*~" -exec rm {} \;
	find . -name "$(CUDASRCINC)" -exec rm {} \;

check: $(EXECUTABLE)

	$(SOURCE:.c=)_SP_$(EXECUTABLE) 2 1

