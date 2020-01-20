# Compile all version of Pi_Simple
# Cecill v2 Emmanuel QUEMENER <emmanuel.quemener@gmail.com>

SOURCE=Pi_Simple.c

COMPUTING=INT32 INT64 FP32 FP64
MARSAGLIA=SHR3 CONG MWC KISS

CC=gcc
CFLAGS=-Wall -O3 -g
LIBRARY=-lm 

all: $(SOURCE)

	$(foreach TVAR,$(COMPUTING),$(foreach TRND,$(MARSAGLIA),$(CC) $(CFLAGS) -DT$(TVAR) -DT$(TRND) -DLONG -DTIME -o $(<:.c=)_$(TVAR)_$(TRND) $< $(LIBRARY); ) )

.PHONY: clean check mrproper

mrproper: $(SOURCE)
	rm -rf $(foreach TVAR,$(COMPUTING),$(foreach TRND,$(MARSAGLIA),$(<:.c=)_$(TVAR)_$(TRND) ) )

	find . -name "*~" -exec rm {} \;

clean:
	find . -name "*~" -exec rm {} \;

