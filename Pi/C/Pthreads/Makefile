EXECUTABLE=Pi_Pthreads_LONG Pi_Pthreads_INT

SOURCE=Pi_Pthreads.c

CC=gcc
CFLAGS=-Wall -O3 -std=c99
LIBRARY=-lm -lpthread

all: $(EXECUTABLE)

Pi_Pthreads_LONG: $(SOURCE)

	$(CC) $(CFLAGS) $(DIRECTIVES) -DLONG -o $@ $< $(LIBRARY)

Pi_Pthreads_INT: $(SOURCE)

	$(CC) $(CFLAGS) $(DIRECTIVES) -DINTEGER -o $@ $< $(LIBRARY)

.PHONY: clean check mrproper

mrproper: 
	rm -rf $(EXECUTABLE)
	find . -name "*~" -exec rm {} \;

clean:
	find . -name "*~" -exec rm {} \;

check: $(EXECUTABLE)

	@echo "To Be Defined"
