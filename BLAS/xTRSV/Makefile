SOURCE=xTRSV.c

CC=gcc
FC=gfortran
CFLAGS=-Wall -O3
LDFLAGS=-lm
CUDASRC=/usr/share/doc/nvidia-cuda-doc/examples/
THUNKING=fortran_thunking.c
CUDASRCINC=fortran_common.h
CUDAINC=/usr/include
CUDALIB=/usr/lib/x86_64-linux-gnu/

GSLINC=/usr/include/gsl

GOTO2=/opt/GotoBLAS2

OPENBLAS=/usr/lib

ATLAS=/usr/lib/atlas-base

ACML=/opt/acml
ACMLINC=$(ACML)/gfortran64_mp/include
ACMLLIB=$(ACML)/gfortran64_mp/lib

#EXECUTABLE=cblas fblas gsl cublas thunking gotoblas acml openblas
EXECUTABLE=cblas fblas gsl cublas thunking openblas

#FORMAT=DOUBLE
FORMAT=FLOAT

#DIRECTIVES=-D$(FORMAT) -DPRINT -DUNIT
#DIRECTIVES=-D$(FORMAT) -DUNIT -DRESULTS -DQUIET
#DIRECTIVES=-DQUIET -DUNIT
DIRECTIVES=-DRESULTS -DQUIET

all: $(EXECUTABLE)

cblas: $(SOURCE)

	$(CC) $(CFLAGS) $(DIRECTIVES) -DFLOAT -DCBLAS $(LDFLAGS) \
		$(SOURCE) $(ATLAS)/libcblas.a $(ATLAS)/libatlas.a -o $(SOURCE:.c=)_SP_$@

	$(CC) $(CFLAGS) $(DIRECTIVES) -DDOUBLE -DCBLAS $(LDFLAGS) \
		$(SOURCE) $(ATLAS)/libcblas.a $(ATLAS)/libatlas.a -o $(SOURCE:.c=)_DP_$@

openblas: $(SOURCE)

	$(CC) $(CFLAGS) $(DIRECTIVES) -DFLOAT -DCBLAS $(LDFLAGS) \
		$(SOURCE) $(OPENBLAS)/libopenblas.a -lpthread -o $(SOURCE:.c=)_SP_$@

	$(CC) $(CFLAGS) $(DIRECTIVES) -DDOUBLE -DCBLAS $(LDFLAGS) \
		$(SOURCE) $(OPENBLAS)/libopenblas.a -lpthread -o $(SOURCE:.c=)_DP_$@

acml: $(SOURCE)

	$(CC) -I$(ACMLINC) $(CFLAGS) $(DIRECTIVES) -DFLOAT -DACML $(LDFLAGS) \
		$(SOURCE) -L$(ACMLLIB) -lacml_mp -lacml_mv \
		-lgomp -lgfortran -lpthread -o $(SOURCE:.c=)_SP_$@

	$(CC) -I$(ACMLINC) $(CFLAGS) $(DIRECTIVES) -DDOUBLE -DACML $(LDFLAGS) \
		$(SOURCE) -L$(ACMLLIB) -lacml_mp -lacml_mv \
		-lgomp -lgfortran -lpthread -o $(SOURCE:.c=)_DP_$@

fblas: $(SOURCE)

	$(CC) $(CFLAGS) $(DIRECTIVES) -DFLOAT -DFBLAS $(LDFLAGS) \
		$(SOURCE) $(ATLAS)/libf77blas.a $(ATLAS)/libatlas.a -lgfortran -o $(SOURCE:.c=)_SP_$@

	$(CC) $(CFLAGS) $(DIRECTIVES) -DDOUBLE -DFBLAS $(LDFLAGS) \
		$(SOURCE) $(ATLAS)/libf77blas.a $(ATLAS)/libatlas.a -lgfortran -o $(SOURCE:.c=)_DP_$@


gsl: $(SOURCE)

	$(CC) -I$(GSLINC) $(CFLAGS) $(DIRECTIVES) -DFLOAT -DGSL $(LDFLAGS) \
		$(SOURCE) -lgslcblas -o $(SOURCE:.c=)_SP_$@

	$(CC) -I$(GSLINC) $(CFLAGS) $(DIRECTIVES) -DDOUBLE -DGSL $(LDFLAGS) \
		$(SOURCE) -lgslcblas -o $(SOURCE:.c=)_DP_$@


cublas: $(SOURCE)

	$(CC) -I$(CUDAINC) -L$(CUDALIB) $(CFLAGS) -DCUBLAS -DFLOAT $(LDFLAGS) \
		$(DIRECTIVES) $(SOURCE) -lcublas -o $(SOURCE:.c=)_SP_$@

	$(CC) -I$(CUDAINC) -L$(CUDALIB) $(CFLAGS) -DCUBLAS -DDOUBLE $(LDFLAGS) \
		$(DIRECTIVES) $(SOURCE) -lcublas -o $(SOURCE:.c=)_DP_$@

thunking: $(SOURCE)

# Copy of source for thunking CUBLAS approach
	cp $(CUDASRC)/$(THUNKING) $(CUDASRC)/$(THUNKING:.c=.h) \
		$(CUDASRC)/$(CUDASRCINC) .

	$(CC) -DCUBLAS_GFORTRAN -I$(CUDAINC) $(CFLAGS) -c $(THUNKING)

	$(CC) -DCUBLAS_GFORTRAN -I$(CUDAINC) -L$(CUDALIB) $(CFLAGS) -DTHUNKING \
		$(LDFLAGS) $(DIRECTIVES) $(SOURCE) -DFLOAT \
		$(THUNKING:.c=.o) -lcublas -o $(SOURCE:.c=)_SP_$@

	$(CC) -DCUBLAS_GFORTRAN -I$(CUDAINC) -L$(CUDALIB) $(CFLAGS) -DTHUNKING \
		$(LDFLAGS) $(DIRECTIVES) $(SOURCE) -DDOUBLE \
		$(THUNKING:.c=.o) -lcublas -o $(SOURCE:.c=)_DP_$@

clean: $(SOURCE)
	find . -name "$(SOURCE:.c=)_*" -exec rm {} \;
	find . -name "$(THUNKING:.c=)*" -exec rm {} \;
	find . -name "*~" -exec rm {} \;
	find . -name "$(CUDASRCINC)" -exec rm {} \;
