rm -r build MonteCarlo.so 
python setup.py build
cp build/lib.linux-x86_64-2.7/MonteCarlo.so .
python MC.py
