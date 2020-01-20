#!/usr/bin/env python

import numpy

SIZE=256

if __name__=='__main__':

    z=numpy.uint32(37)
    w=numpy.uint32(91)

    for i in range(1000):

        z=numpy.uint32(36969*(z&65535)+(z>>16))
        w=numpy.uint32(18000*(w&65535)+(w>>16))
        MWC= lambda : numpy.uint32((z<<16)+w)
        MWCfp= lambda: (MWC() + 1.0) * 2.328306435454494e-10

        print i,MWC(),numpy.uint32(MWCfp()*SIZE)



