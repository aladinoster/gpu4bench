COMPUTING=INT32 INT64 FP32 FP64

MARSAGLIA=SHR3 CONG MWC KISS

#SOURCE=Pi_Hybrid.c Pi_aHybrid.c
SOURCE=Pi_Hybrid.c

CC=mpicc.openmpi
CFLAGS=-Wall -Wno-sequence-point -O3 -std=c99 -fopenmp 
LIBRARY=-lm

all: $(SOURCE)
	$(foreach SRC,$(SOURCE),$(foreach TVAR,$(COMPUTING),$(foreach TRND,$(MARSAGLIA),$(CC) $(CFLAGS) -DT$(TVAR) -DT$(TRND) -DLONG -DTIME -o $(SRC:.c=)_$(TVAR)_$(TRND) $(SRC) $(LIBRARY); ) ) )

.PHONY: clean check mrproper

mrproper: $(SOURCE)
	$(foreach SRC,$(SOURCE),$(foreach TVAR,$(COMPUTING),$(foreach TRND,$(MARSAGLIA),rm $(SRC:.c=)_$(TVAR)_$(TRND); ) ) )
	find . -name "*~" -exec rm {} \;

clean:
	find . -name "*~" -exec rm {} \;

check: $(EXECUTABLE)

	@echo "To be Defined"
