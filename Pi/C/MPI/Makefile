SOURCE=Pi_aMPI.c Pi_MPI.c

COMPUTING=INT32 INT64 FP32 FP64

MARSAGLIA=SHR3 CONG MWC KISS

CC=mpicc.openmpi
CFLAGS=-Wall -Wno-unused-variable -Wno-sequence-point -O3 -std=c99
LIBRARY=-lm

all: $(SOURCE)

	$(foreach SRC,$(SOURCE),$(foreach TVAR,$(COMPUTING),$(foreach TRND,$(MARSAGLIA),$(CC) $(CFLAGS) -DT$(TVAR) -DT$(TRND) -DLONG -DTIME -o $(SRC:.c=)_$(TVAR)_$(TRND) $(SRC) $(LIBRARY); ) ) )

.PHONY: clean check mrproper

mrproper: 
	rm -rf $(foreach SRC,$(SOURCE),$(foreach TVAR,$(COMPUTING),$(foreach TRND,$(MARSAGLIA),$(SRC:.c=)_$(TVAR)_$(TRND) ) ) )
	find . -name "*~" -exec rm {} \;

clean:
	find . -name "*~" -exec rm {} \;
