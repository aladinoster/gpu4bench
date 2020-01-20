import array_module_np
import numpy as np
import pylab

# x = np.arange(0, 2 * np.pi, 0.1)
# y = array_module_np.array_cos_np(x)
# z = array_module_np.array_sin_np(x)
# pylab.plot(x, y,x,z)
# pylab.show()

# x = np.arange(0, 2 * np.pi, 0.1)
# print array_module_np.array_cos_arg_np(x,1)

# x = np.arange(0,64,1).reshape(8,8).astype(np.int64)

# print x

# print array_module_np.array_operation_np(x,1,1,1,1000000000,2008,1010)

SIZE=8

x=np.where(np.random.randn(SIZE,SIZE)>0,1,-1).astype('int32')

print x

y=array_module_np.array_metropolis_np(x,1,0,0.1,10000000,128,256)

array_module_np.array_display_np(x)


print y,y.dtype

print x,x.dtype
