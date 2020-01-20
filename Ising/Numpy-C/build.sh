python setup.py build
[ -e array_module_np.so ] && rm array_module_np.so
ln -s build/lib.linux-x86_64-2.7/array_module_np.so array_module_np.so