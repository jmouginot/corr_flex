#!/usr/bin/env python
# coding: utf-8

import os
filename = os.environ.get('PYTHONSTARTUP')
if filename and os.path.isfile(filename):
    with open(filename) as fobj:
        startup_file = fobj.read()
        exec(startup_file)

import numpy as np
from corr_class import corr_class
from fparam import isp_param

id1 = '20180618_3s'
id2 = '20180630_3s'
p1 = isp_param()
p1.load('../data/'+id1+'.par')
p2 = isp_param()
p2.load('../data/'+id2+'.par')

# initialize corr_class
corr = corr_class()

# input file definitions
corr.im1file = '../data/'+id1+'.slc'
corr.im2file = '../data/'+id2+'.slc'
corr.output_file = id1+'-'+id2+'.offmap_10'
corr.nx_im1 = p1.npix
corr.ny_im1 = p1.nrec
corr.nx_im2 = p2.npix
corr.ny_im2 = p2.nrec

# size of the search windows
corr.nx_search = 64
corr.ny_search = 32

# initial shift between images
corr.x0_shift = -30
corr.y0_shift = -2

# area to be processed in image
corr.x_start = 400
corr.x_end = 600
corr.y_start = 400
corr.y_end = 600

# flag to activate functionalities or debug
corr.off.IsUsed = False
corr.verbose = True
corr.progress = True
corr.debug = True

# output spacing
corr.nx_step = 64
corr.ny_step = 32
# image correlation window size
corr.nx_win = 128
corr.ny_win = 64

corr.p_output_off()

#process data
corr.ampcor()


