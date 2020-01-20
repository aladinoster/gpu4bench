# Compile all version of Pi_OpenACC using Nvidia GPU
# Cecill v2 Emmanuel QUEMENER <emmanuel.quemener@gmail.com>

SOURCE=Pi_OpenACC.c

COMPUTING=INT32 INT64 FP32 FP64
MARSAGLIA=SHR3 CONG MWC KISS

CC=pgcc
CFLAGS=-O3 
LIBRARY=
ACCFLAGS=-acc -ta=tesla:cc60 -Minfo=accel

all: $(SOURCE)

	$(foreach TVAR,$(COMPUTING),$(foreach TRND,$(MARSAGLIA),$(CC) $(CFLAGS) $(ACCFLAGS) -DT$(TVAR) -DT$(TRND) -DLONG -DTIME -o $(<:.c=)_$(TVAR)_$(TRND) $< $(LIBRARY); ) )

.PHONY: clean check mrproper

mrproper:
	rm -rf $(foreach SRC,$(SOURCE),$(foreach TVAR,$(COMPUTING),$(foreach TRND,$(MARSAGLIA),$(SRC:.c=)_$(TVAR)_$(TRND) ) ) )
	find . -name "*~" -exec rm {} \;

clean:
	find . -name "*~" -exec rm {} \;

check: $(EXECUTABLE)

	@echo "To be Defined"
