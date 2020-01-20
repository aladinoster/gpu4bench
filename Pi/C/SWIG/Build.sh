#!/bin/bash
# Delete files created by Swig compilation process
rm _MonteCarlo.so MonteCarlo_wrap.o MonteCarlo_wrap.c MonteCarlo.py MonteCarlo.pyc MCmodule.o
# MonteCarlo.i holds prototype (like .h) of InsideCircle function
# Swig creates MonteCarlo_wrap.c and MonteCarlo.py
swig -python MonteCarlo.i
# Compilation of MCmodule.c as wrapper
gcc -O3 -fpic -c MonteCarlo_wrap.c MCmodule.c -I/usr/include/python2.7/
# Link to create library
gcc -O3 -shared MCmodule.o MonteCarlo_wrap.o -o _MonteCarlo.so
# Execution
python MC.py