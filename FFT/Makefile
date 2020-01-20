SOURCE=FFT2D.c

CC=gcc
CFLAGS=-Wall -O3
LDFLAGS=-lm
CUDADIR=/opt/cuda
CUDASRC=$(CUDADIR)/src
CUDAINC=$(CUDADIR)/include
CUDALIB=$(CUDADIR)/lib64

ACML=/opt/acml
ACMLINC=$(ACML)/gfortran64_mp/include
ACMLLIB=$(ACML)/gfortran64_mp/lib

EXECUTABLE=fftw3 cufft

FORMAT=DOUBLE
#FORMAT=FLOAT

#DIRECTIVES=-D$(FORMAT) -DPRINT -DUNIT
#DIRECTIVES=-D$(FORMAT) -DUNIT -DRESULTS -DQUIET
DIRECTIVES=-DUNIT -DQUIET

all: $(EXECUTABLE)

fftw3: $(SOURCE)

	$(CC) $(CFLAGS) $(DIRECTIVES) -DFLOAT -DFFTW3 $(LDFLAGS) \
		$(SOURCE) \
		-lm -lfftw3f_threads -lfftw3f -lpthread \
		-o $(SOURCE:.c=)_SP_$@

	$(CC) $(CFLAGS) $(DIRECTIVES) -DDOUBLE -DFFTW3 $(LDFLAGS) \
		$(SOURCE) \
		-lm -lfftw3_threads -lfftw3 -lpthread \
		-o $(SOURCE:.c=)_DP_$@

cufft: $(SOURCE)

	$(CC) -I$(CUDAINC) -L$(CUDALIB) $(CFLAGS) -DFLOAT \
		-DCUFFT $(LDFLAGS) \
		$(DIRECTIVES) $(SOURCE) -lcufft -o $(SOURCE:.c=)_SP_$@

	$(CC) -I$(CUDAINC) -L$(CUDALIB) $(CFLAGS) -DDOUBLE \
		-DCUFFT $(LDFLAGS) \
		$(DIRECTIVES) $(SOURCE) -lcufft -o $(SOURCE:.c=)_DP_$@

clean: $(SOURCE)
	find . -name "$(SOURCE:.c=)_*" -exec rm {} \;
	find . -name "*~" -exec rm {} \;

check: $(EXECUTABLE)

	$(SOURCE:.c=)_SP_$(EXECUTABLE) 2 1

