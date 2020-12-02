import os
import numpy as np
from fparam import *
from default_functions_proc import *

class offset_class:

    def __init__(self): # {{{

        self.vx = None
        self.vy = None
        self.snr = None
    # }}}
    def load(self,filename,par_file=None,snr_file=None,init=False): # {{{

        if par_file is None:
            par_file = os.path.splitext(filename)[0]+'.par'
    
        p = off_param()
        p.load(par_file)

        if check_if_binary_is_good(filename,par_file,datatype='complex',file_type='off'):
            dtype = np.complex64
        elif check_if_binary_is_good(filename,par_file,datatype='float',file_type='off'):
            dtype = np.float32
        else:
            dtype = np.complex64

        if dtype == np.complex64:

            off=np.fromfile(filename,dtype=dtype).reshape(p.nrec,p.npix).byteswap()
            self.vx = off.real
            self.vy = off.imag

        elif dtype == np.float32:

            vx_filename = os.path.splitext(filename)[0]+'.vx'
            self.vx = np.fromfile(vx_filename,dtype=dtype).reshape(p.nrec,p.npix).byteswap()
            vy_filename = os.path.splitext(filename)[0]+'.vy'
            self.vy = np.fromfile(vy_filename,dtype=dtype).reshape(p.nrec,p.npix).byteswap()

        #TODO add SNR file
   
    # }}}
#  def plot(self,minimum=None,maximum=None):
#
#        import matplotlib.pyplot as plt
#        v = np.sqrt( self.vx**2 + self.vy**2)
#
#        if minimum is None:
#            minimum = np.min(v)
#        if maximum is None:
#            minimum = np.max(v)

#        v = np.clip(v,minimum,maximum)

#        fig = plt.imshow(v)
#        plt.show(block=False)

