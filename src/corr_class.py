import os
import traceback
import numpy as np
from fparam import off_param
from offset_class import offset_class

class corr_class:

    def __init__(self): # {{{

        self.im1file = ''
        self.im2file = ''
        self.output_file = ''
        self.nx_im1 = 0
        self.nx_im2 = 0
        self.ny_im1 = 0
        self.ny_im2 = 0
        self.datatype=np.complex64
        self.nx_step = 10
        self.ny_step = 10
        self.nx_search = 32
        self.ny_search = 32
        self.nx_win = 32
        self.ny_win = 32

        self.x0_shift = 0
        self.y0_shift = 0

        self.scaling = 1
        self.use_scaling = False
        
        self.x0_polyshift = None
        self.y0_polyshift = None
        self.use_polyshift = False
        
        self.x_start = None
        self.x_end = None
        self.y_start = None
        self.y_end = None
        
        self.mask_file = None # filename for mask
        self.mask = None # mask
        self.use_mask = False

        self.init_shiftmap = None # filename for initial shift map
        self.shiftmap = off_param() # shift map and associated parameters 
        self.shiftmap.m_xshift = np.zeros( (2,2) ,dtype=np.float32)
        self.shiftmap.m_yshift = np.zeros( (2,2) ,dtype=np.float32)
        self.use_shiftmap = False
        
        self.init_searchmap = None
        self.searchmap = off_param()
        self.searchmap.npix = 2
        self.searchmap.nrec = 2
        self.searchmap.m_nx_search = np.zeros( (2,2) ,dtype=np.float32)
        self.searchmap.m_ny_search = np.zeros( (2,2) ,dtype=np.float32)
        self.use_searchmap = False

        self.restart = False
        self.progress = False
        self.verbose = False
        self.debug = False
        
        self.off = off_param()
        self.off.IsUsed=False
        self.off.data = offset_class()
# }}}
    def describe(self): # {{{

        print('Input parameters: ')
        print('self.im1file',self.im1file)
        print('self.im2file',self.im2file)
        print('self.output_file',self.output_file)
        print('self.nx_im1',self.nx_im1)
        print('self.nx_im2',self.nx_im2)
        print('self.ny_im1',self.ny_im1)
        print('self.ny_im2',self.ny_im2)
        print('self.datatype',self.datatype)
        print('\n')
        print('Output parameters: ')
        print('self.nx_step',self.nx_step)
        print('self.ny_step',self.ny_step)
        print('self.nx_search',self.nx_search)
        print('self.ny_search',self.ny_search)
        print('self.nx_win',self.nx_win)
        print('self.ny_win',self.ny_win)
        print('self.x0_shift',self.x0_shift)
        print('self.y0_shift',self.y0_shift)
        print('self.x0_polyshift',self.x0_polyshift)
        print('self.y0_polyshift',self.y0_polyshift)
        print('self.x_start',self.x_start)
        print('self.x_end',self.x_end)
        print('self.y_start',self.y_start)
        print('self.y_end',self.y_end)
        print('\n')
        print('Processing options :')
        print('self.use_shiftmap',self.use_shiftmap)
        print('self.use_searchmap',self.use_searchmap)
        print('self.use_polyshift',self.use_polyshift)
        print('self.use_mask',self.use_mask)
        print('self.restart',self.restart)
        print('\n')
        print('self.restart',self.restart)
        print('self.progress',self.progress)
        print('self.verbose',self.verbose)
        print('self.debug',self.debug)
        print('self.off',self.off)
        print('self.off.IsUsed',self.off.IsUsed)
# }}}
    def p_output_off(self): # {{{
        # prepare/initialize output offmap 
        ny_off = np.int32( (self.y_end-self.y_start)/self.ny_step )
        nx_off = np.int32( (self.x_end-self.x_start)/self.nx_step )

        self.off.data.vx = np.zeros([ny_off,nx_off])
        self.off.data.vy = np.zeros([ny_off,nx_off])
        self.off.data.snr = np.zeros([ny_off,nx_off])
        self.off.r0 = self.x0_shift
        self.off.z0 = self.y0_shift
        self.off.npix = nx_off
        self.off.nrec = ny_off
        self.off.x_start = self.x_start
        self.off.x_end = self.x_start + self.nx_step*nx_off
        self.off.rgsp = self.nx_step
        self.off.y_start = self.y_start
        self.off.y_end = self.y_start + self.ny_step*ny_off
        self.off.azsp = self.ny_step
        self.off.ofw_w = self.nx_win
        self.off.ofw_h = self.ny_win

    #}}}
    def ampcor(self): # {{{
        import ampcor_flex
        if self.progress:
            self.describe()
        # default values for x_start,end y_start,end {{{
        if self.x_start == None:
            self.x_start = self.nx_win/2
        if self.x_end == None:
            self.x_end = self.nx_im1-self.nx_win/2
        if self.y_start == None:
            self.y_start = self.ny_win/2
        if self.y_end == None:
            self.y_end = self.ny_im1-self.ny_win/2
        # }}}
        # writing initial file for offset map{{{
        inif = open(self.output_file+'.in','w')
        inif.write(self.im1file+'\n') #input 1
        inif.write(self.im2file+'\n') #input 2
        inif.write(self.output_file+'\n') # output
        inif.write('%s %s\n'%(self.nx_im1,self.nx_im2))
        inif.write('%s %s %s \n'%(self.y_start,self.y_end,self.ny_step)) #ymin ymax stepy
        inif.write('%s %s %s \n'%(self.x_start,self.x_end,self.nx_step)) #xmin xmax stepx
        inif.write('%s %s \n'%(self.nx_win,self.ny_win)) # window size x y
        inif.write('%s %s \n'%(self.nx_search,self.ny_search)) # search in x y
        inif.write('1 1 \n') # oversampling factor x y
        inif.write('%s %s \n'%(self.x0_shift,self.y0_shift)) # initial offset x y
        inif.write('0. 1.e10\n') # snr threshold, covariance threshold (as set we take everything and do not filter) 
        inif.write('f f\n') # Debug and Display Flags T/F
        inif.close()

        ny_off = np.int32( (self.y_end-self.y_start)/self.ny_step )
        nx_off = np.int32( (self.x_end-self.x_start)/self.nx_step ) # x from self.nx_win/2 to self.nx_im1-self.nx_win/2 with step = nx_step

        if self.off.IsUsed and self.off.data.vx is None:
            self.p_offmap_out()
        #}}}
        if self.restart: # {{{
            if os.path.exists(self.output_file):
                statinfo = os.stat(self.output_file)
                if statinfo.st_size != 0:
                    try:
                        x, offx, y, offy, snr0, tmp1, tmp2, tmp3 = np.loadtxt(self.output_file,unpack=True)
                        if max(y) >= self.y_end:
                            return
                        if max(y) > self.y_start:
                            self.y_start = max(y)-self.ny_step # new y start from previous run
                    except:
                        self.restart=False
                else:
                    self.restart=False
            else:
                self.restart=False
        # }}}
        # read initial offsetmap if provided {{{
        # should be 2 files : init_shiftmap+'
        if self.init_shiftmap is not None and self.use_shiftmap:
            self.shiftmap = off_param()
            self.shiftmap.load(self.init_shiftmap+'.par')

            f=open(self.init_shiftmap+'.vx','rb')
            vx_shiftmap = np.fromfile(f,dtype=np.float32).reshape(self.shiftmap.nrec,self.shiftmap.npix).byteswap()
            f.close

            f=open(self.init_shiftmap+'.vy','rb')
            vy_shiftmap = np.fromfile(f,dtype=np.float32).reshape(self.shiftmap.nrec,self.shiftmap.npix).byteswap()
            f.close

            #shiftmap = vx_shiftmap + 1j * vy_shiftmap
            self.shiftmap.m_xshift = vx_shiftmap
            self.shiftmap.m_yshift = vy_shiftmap
            
            del vx_shiftmap
            del vy_shiftmap
        # }}}    
        # read initial search-distance map if provided {{{
        # should be 2 files : binary + par'
        if self.init_searchmap is not None and self.use_searchmap:
            self.searchmap = off_param()
            self.searchmap.load(self.init_searchmap+'.par')

            f=open(self.init_searchmap+'.vx','rb')
            vx_searchmap = np.fromfile(f,dtype=np.float32).reshape(self.searchmap.nrec,self.searchmap.npix).byteswap()
            f.close

            f=open(init_searchmap+'.vy','rb')
            vy_searchmap = np.fromfile(f,dtype=np.float32).reshape(self.searchmap.nrec,self.searchmap.npix).byteswap()
            f.close

            #searchmap = vx_searchmap + 1j * vy_searchmap
            self.searchmap.m_nx_search = vx_searchmap
            self.searchmap.m_ny_search = vy_searchmap
            del vx_searchmap
            del vy_searchmap
        # }}}
        # read initial mask map if provided {{{
        if self.mask_file is not None and self.use_mask:
            self.mask = off_param()
            self.load(self.mask_file+'.par')
        
            f=open(self.mask_file,'rb')
            self.mask.data = np.fromfile(f,dtype=np.uint8).reshape(self.mask.nrec,self.mask.npix).byteswap()
        # }}} 
        # initialization input, output {{{
        im1 = np.zeros([self.ny_win,self.nx_im1],dtype=self.datatype)
    
        # read master and slave files (may cause memory problem if the file is too large)
        
        f1 = open(self.im1file,'rb')
        
        jstart_master = np.int32(self.y_start-self.ny_win/2).clip(0,self.ny_im1-self.ny_win)
        f1.seek(self.nx_im1*np.dtype(self.datatype).itemsize*jstart_master)
        ny_master = np.int32( self.y_end - self.y_start + self.ny_win + 1 ) # y from self.ny_win/2 to self.ny_im1-self.ny_win/2 with step = ny_step
        if (jstart_master + ny_master) > self.ny_im1:
            ny_master = self.ny_im1 - jstart_master
        if ny_master < self.ny_win:
            return
        
        try:
            master = np.fromfile(f1,count=ny_master*self.nx_im1,dtype=self.datatype).byteswap().reshape(ny_master,self.nx_im1)
            f1.close()
        except Exception:
            print('ny_master:',ny_master)
            print('self.nx_im1:',self.nx_im1)
            print('datatype:',datatype)
    
            print(traceback.format_exc())
            return
        
        if self.debug:
            import matplotlib.pyplot as plt
            plt.imshow(np.abs(master)**0.35)
            plt.show()

        min_y0_shift = self.y0_shift
        max_y0_shift = self.y0_shift

        if self.y0_polyshift is not None: # {{{
            min_y0_shift = min( [ min( np.int32( self.y0_polyshift[0] + np.arange(self.nx_im1)*self.y0_polyshift[1] + self.y_start*self.y0_polyshift[2])), \
                    min(np.int32(self.y0_polyshift[0] + np.arange(self.nx_im1)*self.y0_polyshift[1] + self.y_end*self.y0_polyshift[2])) ])
            max_y0_shift = max( [ max(np.int32( self.y0_polyshift[0] + np.arange(self.nx_im1)*self.y0_polyshift[1] + self.y_start*self.y0_polyshift[2])), \
                    max(np.int32( self.y0_polyshift[0] + np.arange(self.nx_im1) * self.y0_polyshift[1] + self.y_end*self.y0_polyshift[2]))])
        # }}}

        if self.use_shiftmap: # {{{
            # find
            j_shift_start = np.int32((self.y_start - self.shiftmap.y_start) / self.shiftmap.azsp).clip(0,self.shiftmap.nrec-1)
            j_shift_end = np.int32(( ny_master + self.y_start - self.shiftmap.y_start) / self.shiftmap.azsp).clip( 0, self.shiftmap.nrec-1)
            min_y0_shift = np.int32( np.min(self.shiftmap.m_yshift[j_shift_start:j_shift_end,:]) + min_y0_shift)
            max_y0_shift = np.int32( np.max(self.shiftmap.m_yshift[j_shift_start:j_shift_end,:]) + max_y0_shift)
        #}}}

        if self.use_searchmap: # {{{
            j_search_start = np.int32((self.y_start - self.searchmap.y_start) / self.searchmap.azsp).clip(0,self.searchmap.nrec-1)
            j_search_end = np.int32(( ny_master + self.y_start - self.searchmap.y_start) / self.searchmap.azsp).clip( 0, self.searchmap.nrec-1)
            max_ny_search = np.int32( np.max(self.searchmap.m_ny_search[j_search_start:j_search_end,:]) )
        else:
            max_ny_search = self.ny_search
        #}}}

        jmin = np.int32( self.y_start + min_y0_shift - self.ny_win/2 - 2*max_ny_search).clip(0 , self.ny_im2 - self.ny_win - max_ny_search)
        jmax = np.int32( ny_master + self.y_start + max_y0_shift + self.ny_win/2 + 2*max_ny_search).clip( 0, self.ny_im2 - self.ny_win - max_ny_search) #clip between 0 and number of lines in im2

        ny_slave = jmax - jmin + 1
        
        if ny_slave < (self.ny_win+2*max_ny_search):
            return
        
        jstart_slave = jmin
       
        try:
            f2 = open(self.im2file,'rb')
            f2.seek(self.nx_im2*np.dtype(self.datatype).itemsize*jmin)
            slave = np.fromfile(f2,count=ny_slave*self.nx_im2,dtype=self.datatype).byteswap().reshape(ny_slave,self.nx_im2)
            f2.close()
        except Exception:
            print('ny_slave:',ny_slave)
            print('nx_im2:',self.nx_im2)
            print('datatype:',self.datatype)
            print(traceback.format_exc())
            return
    
        ny_off = np.int32( (self.y_end-self.y_start)/self.ny_step)
        nx_off = np.int32( (self.x_end-self.x_start)/self.nx_step) # x from self.nx_win/2 to self.nx_im1-self.nx_win/2 with step = nx_step
        
        jmin_previous = 0
        jmax_previous = 0
        jmin_slave_previous = 0
        jmax_slave_previous = 0
   
        if not self.off.IsUsed:
            if self.restart:
                f = open(self.output_file,'a')
            else:
                f = open(self.output_file,'w')
        # }}}
        for j_off in range(0,ny_off):# loop over y-direction (azimuth for SAR) {{{
    
            if self.verbose:
                print('\n',j_off,'/',ny_off)

            j_im = j_off*self.ny_step+self.y_start
    
            if (j_im+self.y0_shift > 2*self.ny_search+self.ny_win/2) and \
                    (j_im+self.y0_shift < self.ny_im2-self.ny_win/2-2*self.ny_search) and \
                    (j_im > self.ny_win/2) and \
                    (j_im < self.ny_im1-self.ny_win/2):
                
                # read mask if provided {{{
                # if entire line is masked, process_this_line is set to False
                if self.mask is not None and self.use_mask:
                    j_mask = np.int32( (j_im - self.mask.y_start) / self.mask.azsp) # conversion from im y-coord to mask y-coord 
                    if np.max(self.mask.data[j_mask,:]) == 1:
                        process_this_line = True
                    else:
                        process_this_line = False
                else:
                    process_this_line = True
                # }}}
     
                jmin_im1 = np.int32(j_im-self.ny_win/2) # starting line to be read
                jmax_im1 = np.int32(j_im+self.ny_win/2) # ending line to be read
    
                # find lines to read based on init shift {{{
                min_y0_shift = self.y0_shift
                max_y0_shift = self.y0_shift

                if self.y0_polyshift is not None and self.use_polyshift:
                    # not using self.y0_shift
                    min_y0_shift = min(np.int32(self.y0_polyshift[0] + np.arange(self.nx_im1)*self.y0_polyshift[1] + j_im*self.y0_polyshift[2]))
                    max_y0_shift = max(np.int32(self.y0_polyshift[0] + np.arange(self.nx_im1)*self.y0_polyshift[1] + j_im*self.y0_polyshift[2]))

                if self.shiftmap is not None and self.use_shiftmap:
                    # shift is self.y0_shift + self.shiftmap.m_yshift
                    j_shift = np.int32((j_im - self.shiftmap.y_start) / self.shiftmap.azsp) # conversion from image y-coord to offsetmap y-corrd
                    min_y0_shift = np.int32( min( self.shiftmap.m_yshift[j_shift,:]) + min_y0_shift)
                    max_y0_shift = np.int32( max( self.shiftmap.m_yshift[j_shift,:]) + max_y0_shift)
                # }}}
                if self.searchmap is not None and  self.use_searchmap:
                    j_search = np.int32((j_im - self.searchmap.y_start) / self.searchmap.azsp) # conversion from image y-coord to offsetmap y-corrd
                    max_ny_search = np.int32( max( self.searchmap.m_ny_search[j_search,:]) )
                    if self.verbose:
                        print('max_ny_search (search_map):',max_ny_search)
                else:
                    max_ny_search = self.ny_search

                jmin_im2 = np.int32( j_im + min_y0_shift - self.ny_win/2 - 2*max_ny_search)
                jmax_im2 = np.int32( j_im + max_y0_shift + self.ny_win/2 + 2*max_ny_search)
                
                if (jmin_im2 < 0) or (jmax_im2 > self.ny_im2):
                    if self.verbose:
                        print('line not processed: outside slave image, jmin_im2=',jmin_im2)
                    process_this_line=False
    
                if process_this_line:
                    # read master {{{
                         
                    im1 = master[ jmin_im1 - jstart_master : jmax_im1 - jstart_master , : ]
                    im2 = slave[ jmin_im2 - jstart_slave : jmax_im2 - jstart_slave , : ]
                    # }}}
                
                    
                    for i_off in range(0,nx_off):# loop over x-direction (range for SAR) {{{
                        if self.progress and i_off == 0:
                            print("\r>> "+str(int((i_off+j_off*nx_off)/nx_off/ny_off*100.))+' % line='+str(j_im), end='')
                        if self.verbose:
                            print('\n',i_off,'/',nx_off,'- j=',j_off,'/',ny_off)
                            print('Lines to be processed : ', j_im)
                            print('Lines extract for ref : ', jmin_im1, jmax_im1)
                            print('Lines extract for sla : ', jmin_im2, jmax_im2)
    
                        i_im = np.int32(i_off*self.nx_step+self.x_start) # from output offset x-coord to im x-coord
    
                        # read mask if provided {{{
                        if self.mask is not None and self.use_mask:
                            j_mask = np.int32((j_im - self.mask.y_start) / self.mask.azsp) # conversion from im y-coord to mask y-coord
                            i_mask = np.int32((i_im - self.mask.x_start) / self.mask.rgsp) # conversion from im x-coord to mask x-coord
                            if self.mask.data[j_mask,i_mask] == 1:
                                process_this_point = True
                            else:
                                process_this_point = False
                        else:
                            process_this_point = True
                        # }}}
                        
                        
                        if process_this_point:
                            # read init offset if provided {{{
                            x_shift = self.x0_shift
                            y_shift = self.y0_shift
                            if self.verbose:
                                print('Initial offset :',x_shift,y_shift)

                            if self.x0_polyshift  is not None and self.use_polyshift:
                                # use polynom instead of constant as initial shift
                                x_shift = np.int32(self.x0_polyshift[0] + i_im*self.x0_polyshift[1] + j_im*self.x0_polyshift[2])
                                y_shift = np.int32(self.y0_polyshift[0] + i_im*self.y0_polyshift[1] + j_im*self.y0_polyshift[2])
                             
                            if self.shiftmap is not None and self.use_shiftmap:
                                # adding shiftmap to initial offset
                                if self.verbose:
                                    print('adding shiftmap to initial offset')
                                j_shift = np.int32((j_im - self.shiftmap.y_start) / self.shiftmap.azsp)
                                i_shift = np.int32((i_im - self.shiftmap.x_start) / self.shiftmap.rgsp)
                                if self.verbose:
                                    print('j_im,i_im',j_im,i_im)
                                    print('j_shift,i_shift',j_shift,i_shift)
                                x_shift = np.int32(self.shiftmap.m_xshift[j_shift,i_shift] + x_shift)
                                y_shift = np.int32(self.shiftmap.m_yshift[j_shift,i_shift] + y_shift)
    
                            # }}}
                            if self.searchmap is not None and self.use_searchmap:
                                j_search = np.int32((j_im - self.searchmap.y_start) / self.searchmap.azsp)
                                i_search = np.int32((i_im - self.searchmap.x_start) / self.searchmap.rgsp)
                                nx_search = np.int32(self.searchmap.m_nx_search[j_search,i_search])
                                ny_search = np.int32(self.searchmap.m_ny_search[j_search,i_search])
                            else:
                                nx_search = self.nx_search
                                ny_search = self.ny_search
                            

                            if self.verbose:
                                print('Initial offset: ',x_shift,y_shift)
    
                            if (i_im+x_shift > 2*nx_search + self.nx_win/2) and \
                                    ( i_im + x_shift < self.nx_im2 - self.nx_win / 2 - 2 * nx_search) and \
                                    ( i_im > self.nx_win / 2) and \
                                    ( i_im < self.nx_im1 - self.nx_win / 2): # if inside im2
                                
                                # extract ref and srch sub-images from lines previously read {{{
                                ref = im1[ : , i_im - np.int32(self.nx_win/2) : i_im + np.int32(self.nx_win/2) ]
                                # IF ny_win and nx_win are not even then you slave and master sub-images might be shifted by 1 pixel
                                srch = im2[ np.int32(j_im + y_shift - self.ny_win/2 - 2*ny_search - jmin_im2): \
                                        np.int32(j_im + y_shift + self.ny_win/2 + 2*ny_search - jmin_im2), \
                                        i_im + np.int32(x_shift - self.nx_win/2 - 2*nx_search) : i_im + np.int32(x_shift + self.nx_win/2 + 2*nx_search) ]
                                # }}}
                                if self.verbose:
                                    print('Columns extract for ref : ', i_im - np.int32(self.nx_win/2), i_im + np.int32(self.nx_win/2))
                                    print('Columns extract for sla : ', i_im + np.int32(x_shift - self.nx_win/2 - 2*nx_search), i_im + np.int32(x_shift + self.nx_win/2 + 2*nx_search))
                                # initialize output values {{{
                                r_shftxosc = np.float64(0)
                                r_shftyosc = np.float64(0)
                                r_snr = np.float64(0)
                                r_cov = np.zeros(3, dtype=np.float64)
                                # }}}
                                
                                if self.scaling is not None and self.use_scaling: # {{{
                                    
                                    if self.verbose:
                                        print(np.size(srch,0),np.size(ref,0),4*ny_search,np.size(ref,0) + 4*ny_search)
                                        
                                    ref2 = np.zeros( ( int(ref.shape[0]/self.scaling), int(ref.shape[1]/self.scaling)) ,dtype=np.complex64)
                                    srch2 = np.zeros( ( int(srch.shape[0]/self.scaling), int(srch.shape[1]/self.scaling)),dtype=np.complex64 )
                                    
                                    for i_scal in range(0,self.scaling):
                                        for j_scal in range(0,self.scaling):
                                            ref2 += ref[j_scal::self.scaling,i_scal::self.scaling][0:ref2.shape[0],0:ref2.shape[1]]
                                            srch2 += srch[j_scal::self.scaling,i_scal::self.scaling][0:srch2.shape[0],0:srch2.shape[1]]
                                            
                                    ref2 = ref2 / self.scaling / self.scaling # normalize
                                    ref = ref2
                                    del ref2
                                    
                                    srch2 = srch2 / self.scaling / self.scaling # normalize
                                    srch = srch2
                                    del srch2
                                    
                                    ny_search = ny_search/self.scaling
                                    nx_search = nx_search/self.scaling
                                    
                                    if self.verbose:
                                        print(np.size(srch,0),np.size(ref,0),4*ny_search,np.size(ref,0) + 4*ny_search)
                                # }}}
                                
                                if np.size(srch,0) == (np.size(ref,0) + 4*ny_search):
                                    if self.verbose:
                                        print('Running fortran..')
                                    r_shftxosc, r_shftyosc, r_snr, r_cov = \
                                            ampcor_flex.ampcor_flex(ref.T,srch.T,nx_search,ny_search) # .T transpose and put in Fortran order
    
                                    if self.scaling is not None and self.use_scaling: #{{{
                                        r_shftxosc = r_shftxosc*self.scaling
                                        r_shftyosc = r_shftyosc*self.scaling
                                    # }}}
                                    if self.verbose:
                                        print(r_shftxosc, r_shftyosc, r_snr, r_cov)

                                    #r_shftxosc = \
                                    #        ampcor_flex_ifort.ampcor_flex(ref.T) # .T transpose and put in Fortran order
    		   
                                del ref
                                del srch

                                
                                if r_snr !=0:
                                    #format(1x,i7,1x,f9.3,1x,i7,1x,f11.3,1x,f10.5,1x,f10.6,1x,f10.6,1x,f10.6) in fortran code
                                    if self.off.IsUsed:
                                        xoff = np.int32((i_im - self.off.x_start) / self.off.rgsp)
                                        yoff = np.int32( (j_im - self.off.y_start) / self.off.azsp)
                                        #if r_snr > self.off.data.snr[ yoff , xoff ]:
                                        self.off.data.vx[ yoff , xoff ] = r_shftxosc+x_shift
                                        self.off.data.vy[ yoff , xoff ] = r_shftyosc+y_shift
                                        self.off.data.snr[ yoff , xoff ] = r_snr
                                    else:
                                        f.write(' '+'{:7d}'.format(np.int32(i_im))+' '+'{:9.3f}'.format(r_shftxosc+x_shift)+\
                                                ' '+'{:7d}'.format(np.int32(j_im))+' '+'{:11.3f}'.format(r_shftyosc+y_shift)+\
                                                ' '+'{:10.5f}'.format(r_snr)+ \
                                                ' '+'{:10.6f}'.format(r_cov[0])+ \
                                                ' '+'{:10.6f}'.format(r_cov[1])+ \
                                                ' '+'{:10.6f}'.format(r_cov[2])+'\n' )
    			
    		                    #r_shftxosc_arr[j_off,i_off]=r_shftxosc+x0_shift
    		                    #r_shftyosc_arr[j_off,i_off]=r_shftyosc+y0_shift
    		                    #r_snr_arr[j_off,i_off]=r_snr
    		                    #r_cov_arr[:,j_off,i_off]=r_cov
    		
    		                    #print(i_off,nx_off,'-',j_off,ny_off)
    	                    #}}}
    	                #}}}
        if not self.off.IsUsed:
            f.close
    # }}}    
    def ampcor2(self): # {{{
       
        import ampcor_tflex
        if self.progress:
            self.describe()
        # default values for x_start,end y_start,end {{{
        if self.x_start == None:
            self.x_start = self.nx_win/2
        if self.x_end == None:
            self.x_end = self.nx_im1-self.nx_win/2
        if self.y_start == None:
            self.y_start = self.ny_win/2
        if self.y_end == None:
            self.y_end = self.ny_im1-self.ny_win/2
        # }}}
        # writing initial file for offset map{{{
        inif = open(self.output_file+'.in','w')
        inif.write(self.im1file+'\n') #input 1
        inif.write(self.im2file+'\n') #input 2
        inif.write(self.output_file+'\n') # output
        inif.write('%s %s\n'%(self.nx_im1,self.nx_im2))
        inif.write('%s %s %s \n'%(self.y_start,self.y_end,self.ny_step)) #ymin ymax stepy
        inif.write('%s %s %s \n'%(self.x_start,self.x_end,self.nx_step)) #xmin xmax stepx
        inif.write('%s %s \n'%(self.nx_win,self.ny_win)) # window size x y
        inif.write('%s %s \n'%(self.nx_search,self.ny_search)) # search in x y
        inif.write('1 1 \n') # oversampling factor x y
        inif.write('%s %s \n'%(self.x0_shift,self.y0_shift)) # initial offset x y
        inif.write('0. 1.e10\n') # snr threshold, covariance threshold (as set we take everything and do not filter) 
        inif.write('f f\n') # Debug and Display Flags T/F
        inif.close()

        ny_off = np.int32( (self.y_end-self.y_start)/self.ny_step )
        nx_off = np.int32( (self.x_end-self.x_start)/self.nx_step ) # x from self.nx_win/2 to self.nx_im1-self.nx_win/2 with step = nx_step

        if self.off.IsUsed and self.off.data.vx is None:
            self.p_offmap_out()
        #}}}
        if self.restart: # {{{
            if os.path.exists(self.output_file):
                statinfo = os.stat(self.output_file)
                if statinfo.st_size != 0:
                    try:
                        x, offx, y, offy, snr0, tmp1, tmp2, tmp3 = np.loadtxt(self.output_file,unpack=True)
                        if max(y) >= self.y_end:
                            return
                        if max(y) > self.y_start:
                            self.y_start = max(y)-self.ny_step # new y start from previous run
                    except:
                        self.restart=False
                else:
                    self.restart=False
            else:
                self.restart=False
        # }}}
        # read initial offsetmap if provided {{{
        # should be 2 files : init_shiftmap+'
        if self.init_shiftmap is not None and self.use_shiftmap:
            self.shiftmap = off_param()
            self.shiftmap.load(self.init_shiftmap+'.par')

            f=open(self.init_shiftmap+'.vx','rb')
            vx_shiftmap = np.fromfile(f,dtype=np.float32).reshape(self.shiftmap.nrec,self.shiftmap.npix).byteswap()
            f.close

            f=open(self.init_shiftmap+'.vy','rb')
            vy_shiftmap = np.fromfile(f,dtype=np.float32).reshape(self.shiftmap.nrec,self.shiftmap.npix).byteswap()
            f.close

            #shiftmap = vx_shiftmap + 1j * vy_shiftmap
            self.shiftmap.m_xshift = vx_shiftmap
            self.shiftmap.m_yshift = vy_shiftmap
            
            del vx_shiftmap
            del vy_shiftmap
        # }}}    

        # initialization input, output {{{
        im1 = np.zeros([self.ny_win,self.nx_im1],dtype=self.datatype)
    
        # read master and slave files (may cause memory problem if the file is too large)
        
        f1 = open(self.im1file,'rb')
        
        jstart_master = np.int32(self.y_start-self.ny_win/2).clip(0,self.ny_im1-self.ny_win)
        f1.seek(self.nx_im1*np.dtype(self.datatype).itemsize*jstart_master)
        ny_master = np.int32( self.y_end - self.y_start + self.ny_win + 1 ) # y from self.ny_win/2 to self.ny_im1-self.ny_win/2 with step = ny_step
        if (jstart_master + ny_master) > self.ny_im1:
            ny_master = self.ny_im1 - jstart_master
        if ny_master < self.ny_win:
            return
        
        try:
            master = np.fromfile(f1,count=ny_master*self.nx_im1,dtype=self.datatype).byteswap().reshape(ny_master,self.nx_im1)
            f1.close()
        except Exception:
            print('ny_master:',ny_master)
            print('self.nx_im1:',self.nx_im1)
            print('datatype:',datatype)
    
            print(traceback.format_exc())
            return
        
        if self.debug:
            import matplotlib.pyplot as plt
            plt.imshow(np.abs(master)**0.35)
            plt.show()

        min_y0_shift = self.y0_shift
        max_y0_shift = self.y0_shift

        if self.y0_polyshift is not None: # {{{
            min_y0_shift = min( [ min( np.int32( self.y0_polyshift[0] + np.arange(self.nx_im1)*self.y0_polyshift[1] + self.y_start*self.y0_polyshift[2])), \
                    min(np.int32(self.y0_polyshift[0] + np.arange(self.nx_im1)*self.y0_polyshift[1] + self.y_end*self.y0_polyshift[2])) ])
            max_y0_shift = max( [ max(np.int32( self.y0_polyshift[0] + np.arange(self.nx_im1)*self.y0_polyshift[1] + self.y_start*self.y0_polyshift[2])), \
                    max(np.int32( self.y0_polyshift[0] + np.arange(self.nx_im1) * self.y0_polyshift[1] + self.y_end*self.y0_polyshift[2]))])
        # }}}

        if self.use_shiftmap: # {{{
            # find
            j_shift_start = np.int32((self.y_start - self.shiftmap.y_start) / self.shiftmap.azsp).clip(0,self.shiftmap.nrec-1)
            j_shift_end = np.int32(( ny_master + self.y_start - self.shiftmap.y_start) / self.shiftmap.azsp).clip( 0, self.shiftmap.nrec-1)
            min_y0_shift = np.int32( np.min(self.shiftmap.m_yshift[j_shift_start:j_shift_end,:].imag) + min_y0_shift)
            max_y0_shift = np.int32( np.max(self.shiftmap.m_yshift[j_shift_start:j_shift_end,:].imag) + max_y0_shift)
        #}}}

        if self.use_searchmap: # {{{
            j_search_start = np.int32((self.y_start - self.searchmap.y_start) / self.searchmap.azsp).clip(0,self.searchmap.nrec-1)
            j_search_end = np.int32(( ny_master + self.y_start - self.searchmap.y_start) / self.searchmap.azsp).clip( 0, self.searchmap.nrec-1)
            max_ny_search = np.int32( np.max(self.searchmap.m_ny_search[j_search_start:j_search_end,:]) )
        else:
            max_ny_search = self.ny_search
        #}}}

        jmin = np.int32( self.y_start + min_y0_shift - self.ny_win/2 - 2*max_ny_search).clip(0 , self.ny_im2 - self.ny_win - max_ny_search)
        jmax = np.int32( ny_master + self.y_start + max_y0_shift + self.ny_win/2 + 2*max_ny_search).clip( 0, self.ny_im2 - self.ny_win - max_ny_search) #clip between 0 and number of lines in im2

        ny_slave = jmax - jmin + 1
        
        if ny_slave < (self.ny_win+2*max_ny_search):
            return
        
        jstart_slave = jmin
       
        try:
            f2 = open(self.im2file,'rb')
            f2.seek(self.nx_im2*np.dtype(self.datatype).itemsize*jmin)
            slave = np.fromfile(f2,count=ny_slave*self.nx_im2,dtype=self.datatype).byteswap().reshape(ny_slave,self.nx_im2)
            f2.close()
        except Exception:
            print('ny_slave:',ny_slave)
            print('nx_im2:',self.nx_im2)
            print('datatype:',self.datatype)
            print(traceback.format_exc())
            return
    
        ny_off = np.int32( (self.y_end-self.y_start)/self.ny_step)
        nx_off = np.int32( (self.x_end-self.x_start)/self.nx_step) # x from self.nx_win/2 to self.nx_im1-self.nx_win/2 with step = nx_step
        
        jmin_previous = 0
        jmax_previous = 0
        jmin_slave_previous = 0
        jmax_slave_previous = 0
   
        if not self.off.IsUsed:
            if self.restart:
                f = open(self.output_file,'a')
            else:
                f = open(self.output_file,'w')

        mask_fortran = np.zeros( (self.off.nrec,self.off.npix) ,dtype=np.int8) #Byte (-128 to 127)
        nx_search_fortran = np.zeros( (self.off.nrec,self.off.npix) ,dtype=np.int16)
        ny_search_fortran = np.zeros( (self.off.nrec,self.off.npix) ,dtype=np.int16)
        j0_master_fortran = np.zeros( (self.off.nrec,self.off.npix) ,dtype=np.int32) #Integer (-2147483648 to 2147483647)
        jsize_master_fortran = np.zeros( (self.off.nrec,self.off.npix) ,dtype=np.int16) #Integer (-32768 to 32767)
        i0_master_fortran = np.zeros( (self.off.nrec,self.off.npix) ,dtype=np.int32)
        isize_master_fortran = np.zeros( (self.off.nrec,self.off.npix) ,dtype=np.int16)
        j0_slave_fortran = np.zeros( (self.off.nrec,self.off.npix) ,dtype=np.int32) #Integer (-2147483648 to 2147483647)
        jsize_slave_fortran = np.zeros( (self.off.nrec,self.off.npix) ,dtype=np.int16) #Integer (-32768 to 32767)
        i0_slave_fortran = np.zeros( (self.off.nrec,self.off.npix) ,dtype=np.int32)
        isize_slave_fortran = np.zeros( (self.off.nrec,self.off.npix) ,dtype=np.int16)
        # }}}
        for j_off in range(0,ny_off):# loop over y-direction (azimuth for SAR) {{{
    
            if self.verbose:
                print('\n',j_off,'/',ny_off)

            j_im = j_off*self.ny_step+self.y_start
    
            if (j_im+self.y0_shift > 2*self.ny_search+self.ny_win/2) and \
                    (j_im+self.y0_shift < self.ny_im2-self.ny_win/2-2*self.ny_search) and \
                    (j_im > self.ny_win/2) and \
                    (j_im < self.ny_im1-self.ny_win/2):
                
                # read mask if provided {{{
                # if entire line is masked, process_this_line is set to False
                if self.mask is not None and self.use_mask:
                    j_mask = np.int32( (j_im - self.mask.y_start) / self.mask.azsp) # conversion from im y-coord to mask y-coord 
                    if np.max(self.mask.data[j_mask,:]) == 1:
                        process_this_line = True
                    else:
                        process_this_line = False
                else:
                    process_this_line = True
                # }}}
     
                jmin_im1 = np.int32(j_im-self.ny_win/2) # starting line to be read
                jmax_im1 = np.int32(j_im+self.ny_win/2) # ending line to be read
    
                # find lines to read based on init shift {{{
                min_y0_shift = self.y0_shift
                max_y0_shift = self.y0_shift

                if self.y0_polyshift is not None and self.use_polyshift:
                    # not using self.y0_shift
                    min_y0_shift = min(np.int32(self.y0_polyshift[0] + np.arange(self.nx_im1)*self.y0_polyshift[1] + j_im*self.y0_polyshift[2]))
                    max_y0_shift = max(np.int32(self.y0_polyshift[0] + np.arange(self.nx_im1)*self.y0_polyshift[1] + j_im*self.y0_polyshift[2]))

                if self.shiftmap is not None and self.use_shiftmap:
                    # shift is self.y0_shift + self.shiftmap.m_yshift
                    j_shift = np.int32((j_im - self.shiftmap.y_start) / self.shiftmap.azsp) # conversion from image y-coord to offsetmap y-corrd
                    min_y0_shift = np.int32( min( self.shiftmap.m_yshift[j_shift,:].imag) + min_y0_shift)
                    max_y0_shift = np.int32( max( self.shiftmap.m_yshift[j_shift,:].imag) + max_y0_shift)
                # }}}
                if self.searchmap is not None and  self.use_searchmap:
                    j_search = np.int32((j_im - self.searchmap.y_start) / self.searchmap.azsp) # conversion from image y-coord to offsetmap y-corrd
                    max_ny_search = np.int32( max( self.searchmap.m_ny_search[j_search,:]) )
                    if self.verbose:
                        print('max_ny_search (search_map):',max_ny_search)
                else:
                    max_ny_search = self.ny_search

                jmin_im2 = np.int32( j_im + min_y0_shift - self.ny_win/2 - 2*max_ny_search)
                jmax_im2 = np.int32( j_im + max_y0_shift + self.ny_win/2 + 2*max_ny_search)
                
                if (jmin_im2 < 0) or (jmax_im2 > self.ny_im2):
                    process_this_line=False
    
                if process_this_line:
                     
                    for i_off in range(0,nx_off):# loop over x-direction (range for SAR) {{{
                        
                        if self.verbose:
                            print('\n',i_off,'/',nx_off,'- j=',j_off,'/',ny_off)

                        i_im = np.int32(i_off*self.nx_step+self.x_start) # from output offset x-coord to im x-coord
    
                        # read mask if provided {{{
                        if self.mask is not None and self.use_mask:
                            j_mask = np.int32((j_im - self.mask.y_start) / self.mask.azsp) # conversion from im y-coord to mask y-coord
                            i_mask = np.int32((i_im - self.mask.x_start) / self.mask.rgsp) # conversion from im x-coord to mask x-coord
                            if self.mask.data[j_mask,i_mask] == 1:
                                process_this_point = True
                            else:
                                process_this_point = False
                        else:
                            process_this_point = True
                        # }}}
                        
                        if process_this_point:
                            # read init offset if provided {{{
                            x_shift = self.x0_shift
                            y_shift = self.y0_shift

                            if self.x0_polyshift  is not None and self.use_polyshift:
                                # use polynom instead of constant as initial shift
                                x_shift = np.int32(self.x0_polyshift[0] + i_im*self.x0_polyshift[1] + j_im*self.x0_polyshift[2])
                                y_shift = np.int32(self.y0_polyshift[0] + i_im*self.y0_polyshift[1] + j_im*self.y0_polyshift[2])
                             
                            if self.shiftmap is not None and self.use_shiftmap:
                                # adding shiftmap to initial offset
                                j_shift = np.int32((j_im - self.shiftmap.y_start) / self.shiftmap.azsp)
                                i_shift = np.int32((i_im - self.shiftmap.x_start) / self.shiftmap.rgsp)
                                x_shift = np.int32(self.shiftmap.m_xshift[j_shift,i_shift] + x_shift)
                                y_shift = np.int32(self.shiftmap.m_yshift[j_shift,i_shift].imag + y_shift)
    
                            # }}}
                            if self.searchmap is not None and self.use_searchmap:
                                j_search = np.int32((j_im - self.searchmap.y_start) / self.searchmap.azsp)
                                i_search = np.int32((i_im - self.searchmap.x_start) / self.searchmap.rgsp)
                                nx_search = np.int32(self.searchmap.m_nx_search[j_search,i_search])
                                ny_search = np.int32(self.searchmap.m_ny_search[j_search,i_search])
                            else:
                                nx_search = self.nx_search
                                ny_search = self.ny_search
                            # }}}

                            if (i_im+x_shift > 2*nx_search + self.nx_win/2) and \
                                    ( i_im + x_shift < self.nx_im2 - self.nx_win / 2 - 2 * nx_search) and \
                                    ( i_im > self.nx_win / 2) and \
                                    ( i_im < self.nx_im1 - self.nx_win / 2): # if inside im2
                                
                                # extract ref and srch sub-images from lines previously read {{{
                                #ref = im1[ : , i_im - np.int32(self.nx_win/2) : i_im + np.int32(self.nx_win/2) ]
                                xoff = np.int32((i_im - self.off.x_start) / self.off.rgsp)
                                yoff = np.int32( (j_im - self.off.y_start) / self.off.azsp)
                                
                                nx_search_fortran[yoff,xoff] = nx_search
                                ny_search_fortran[yoff,xoff] = ny_search

                                j0_master_fortran[yoff,xoff] = np.int32(j_im-self.ny_win/2)
                                jsize_master_fortran[yoff,xoff] = np.int32(self.ny_win)
                                i0_master_fortran[yoff,xoff] = i_im - np.int32(self.nx_win/2)
                                isize_master_fortran[yoff,xoff] = np.int32(self.nx_win)

                                j0_slave_fortran[yoff,xoff] = np.int32(j_im + y_shift - self.ny_win/2 - 2*ny_search ) 
                                jsize_slave_fortran[yoff,xoff] = np.int32(self.ny_win + 4*ny_search)
                                i0_slave_fortran[yoff,xoff] = i_im + np.int32(x_shift - self.nx_win/2 - 2*nx_search)
                                isize_slave_fortran[yoff,xoff] = np.int32(self.nx_win + 4*nx_search)

                                mask_fortran[yoff,xoff] = 1

                                if self.verbose:
                                    print('nx_search, ny_search', nx_search_fortran[yoff,xoff],ny_search_fortran[yoff,xoff])
                                    print('j0_master_fortran, i0_master_fortran',j0_master_fortran[yoff,xoff],i0_master_fortran[yoff,xoff])
                                    print('jsize_master_fortran,isize_master_fortran',jsize_master_fortran[yoff,xoff],isize_master_fortran[yoff,xoff])
                                    print('j0_slave_fortran,i0_slave_fortran',j0_slave_fortran[yoff,xoff],i0_slave_fortran[yoff,xoff])
                                    print('jsize_slave_fortran,isize_slave_fortran',jsize_slave_fortran[yoff,xoff],isize_slave_fortran[yoff,xoff])
    	                    #}}}
    	                #}}}

        
        xoff_out_map, yoff_out_map, r_snr_map, r_cov_map = ampcor_tflex.ampcor_tflex(master.T,slave.T, \
                nx_search_fortran.T, ny_search_fortran.T, mask_fortran.T, \
                i0_master_fortran.T, j0_master_fortran.T, isize_master_fortran.T, jsize_master_fortran.T, \
                i0_slave_fortran.T, j0_slave_fortran.T, isize_slave_fortran.T, jsize_slave_fortran.T)

        del master, slave, ny_search_fortran, nx_search_fortran, mask_fortran, i0_master_fortran, j0_master_fortran, isize_master_fortran
        del i0_slave_fortran, j0_slave_fortran, isize_slave_fortran, jsize_slave_fortran

        for j_off in range(0,ny_off):# loops to store results in self.off {{{
            j_im = j_off*self.ny_step+self.y_start
    
            if (j_im+self.y0_shift > 2*self.ny_search+self.ny_win/2) and \
                    (j_im+self.y0_shift < self.ny_im2-self.ny_win/2-2*self.ny_search) and \
                    (j_im > self.ny_win/2) and \
                    (j_im < self.ny_im1-self.ny_win/2):
                
                # read mask if provided {{{
                # if entire line is masked, process_this_line is set to False
                if self.mask is not None and self.use_mask:
                    j_mask = np.int32( (j_im - self.mask.y_start) / self.mask.azsp) # conversion from im y-coord to mask y-coord 
                    if np.max(self.mask.data[j_mask,:]) == 1:
                        process_this_line = True
                    else:
                        process_this_line = False
                else:
                    process_this_line = True
                # }}}
     
                jmin_im1 = np.int32(j_im-self.ny_win/2) # starting line to be read
                jmax_im1 = np.int32(j_im+self.ny_win/2) # ending line to be read
    
                # find lines to read based on init shift {{{
                min_y0_shift = self.y0_shift
                max_y0_shift = self.y0_shift

                if self.y0_polyshift is not None and self.use_polyshift:
                    # not using self.y0_shift
                    min_y0_shift = min(np.int32(self.y0_polyshift[0] + np.arange(self.nx_im1)*self.y0_polyshift[1] + j_im*self.y0_polyshift[2]))
                    max_y0_shift = max(np.int32(self.y0_polyshift[0] + np.arange(self.nx_im1)*self.y0_polyshift[1] + j_im*self.y0_polyshift[2]))

                if self.shiftmap is not None and self.use_shiftmap:
                    # shift is self.y0_shift + self.shiftmap.m_yshift
                    j_shift = np.int32((j_im - self.shiftmap.y_start) / self.shiftmap.azsp) # conversion from image y-coord to offsetmap y-corrd
                    min_y0_shift = np.int32( min( self.shiftmap.m_yshift[j_shift,:].imag) + min_y0_shift)
                    max_y0_shift = np.int32( max( self.shiftmap.m_yshift[j_shift,:].imag) + max_y0_shift)
                # }}}
                if self.searchmap is not None and  self.use_searchmap:
                    j_search = np.int32((j_im - self.searchmap.y_start) / self.searchmap.azsp) # conversion from image y-coord to offsetmap y-corrd
                    max_ny_search = np.int32( max( self.searchmap.m_ny_search[j_search,:]) )
                else:
                    max_ny_search = self.ny_search

                jmin_im2 = np.int32( j_im + min_y0_shift - self.ny_win/2 - 2*max_ny_search)
                jmax_im2 = np.int32( j_im + max_y0_shift + self.ny_win/2 + 2*max_ny_search)
                
                if (jmin_im2 < 0) or (jmax_im2 > self.ny_im2):
                    if self.verbose:
                        print('line not processed: outside slave image, jmin_im2=',jmin_im2)
                    process_this_line=False
    
                if process_this_line:
                     
                    for i_off in range(0,nx_off):# loop over x-direction (range for SAR) {{{
    
                        i_im = np.int32(i_off*self.nx_step+self.x_start) # from output offset x-coord to im x-coord
    
                        # read mask if provided {{{
                        if self.mask is not None and self.use_mask:
                            j_mask = np.int32((j_im - self.mask.y_start) / self.mask.azsp) # conversion from im y-coord to mask y-coord
                            i_mask = np.int32((i_im - self.mask.x_start) / self.mask.rgsp) # conversion from im x-coord to mask x-coord
                            if self.mask.data[j_mask,i_mask] == 1:
                                process_this_point = True
                            else:
                                process_this_point = False
                        else:
                            process_this_point = True
                        # }}}
                        
                        if process_this_point:
                            # read init offset if provided {{{
                            x_shift = self.x0_shift
                            y_shift = self.y0_shift

                            if self.x0_polyshift  is not None and self.use_polyshift:
                                # use polynom instead of constant as initial shift
                                x_shift = np.int32(self.x0_polyshift[0] + i_im*self.x0_polyshift[1] + j_im*self.x0_polyshift[2])
                                y_shift = np.int32(self.y0_polyshift[0] + i_im*self.y0_polyshift[1] + j_im*self.y0_polyshift[2])
                             
                            if self.shiftmap is not None and self.use_shiftmap:
                                # adding shiftmap to initial offset
                                if self.verbose:
                                    print('adding shiftmap to initial offset')
                                j_shift = np.int32((j_im - self.shiftmap.y_start) / self.shiftmap.azsp)
                                i_shift = np.int32((i_im - self.shiftmap.x_start) / self.shiftmap.rgsp)
                                if self.verbose:
                                    print('j_im,i_im',j_im,i_im)
                                    print('j_shift,i_shift',j_shift,i_shift)
                                x_shift = np.int32(self.shiftmap.m_xshift[j_shift,i_shift] + x_shift)
                                y_shift = np.int32(self.shiftmap.m_yshift[j_shift,i_shift] + y_shift)
    
                            # }}}
    
                            if (i_im+x_shift > 2*nx_search + self.nx_win/2) and \
                                    ( i_im + x_shift < self.nx_im2 - self.nx_win / 2 - 2 * nx_search) and \
                                    ( i_im > self.nx_win / 2) and \
                                    ( i_im < self.nx_im1 - self.nx_win / 2): # if inside im2
                                
                                if r_snr_map[xoff, yoff] !=0:
                                    #format(1x,i7,1x,f9.3,1x,i7,1x,f11.3,1x,f10.5,1x,f10.6,1x,f10.6,1x,f10.6) in fortran code
                                    if self.off.IsUsed:
                                        xoff = np.int32((i_im - self.off.x_start) / self.off.rgsp)
                                        yoff = np.int32( (j_im - self.off.y_start) / self.off.azsp)
                                        #if r_snr > self.off.data.snr[ yoff , xoff ]:
                                        self.off.data.vx[ yoff , xoff ] = xoff_out_map[xoff, yoff]+x_shift
                                        self.off.data.vy[ yoff , xoff ] = yoff_out_map[xoff, yoff]+y_shift
                                        self.off.data.snr[ yoff , xoff ] = r_snr_map[xoff, yoff]
                                    else:
                                        f.write(' '+'{:7d}'.format(np.int32(i_im))+' '+'{:9.3f}'.format(r_shftxosc+x_shift)+\
                                                ' '+'{:7d}'.format(np.int32(j_im))+' '+'{:11.3f}'.format(r_shftyosc+y_shift)+\
                                                ' '+'{:10.5f}'.format(r_snr)+ \
                                                ' '+'{:10.6f}'.format(r_cov[0])+ \
                                                ' '+'{:10.6f}'.format(r_cov[1])+ \
                                                ' '+'{:10.6f}'.format(r_cov[2])+'\n' )

    	                    #}}}
    	                #}}}


        if not self.off.IsUsed:
            f.close
    # }}}    

    def off_write(self,binary=False): # {{{
        
        # default values for x_start,end y_start,end {{{
        if self.x_start == None:
            self.x_start = self.nx_win/2
        if self.x_end == None:
            self.x_end = self.nx_im1-self.nx_win/2
        if self.y_start == None:
            self.y_start = self.ny_win/2
        if self.y_end == None:
            self.y_end = self.ny_im1-self.ny_win/2
        # }}}
        # writing initial file for offset map{{{
        inif = open(self.output_file+'.in','w')
        inif.write(self.im1file+'\n') #input 1
        inif.write(self.im2file+'\n') #input 2
        inif.write(self.output_file+'\n') # output
        inif.write('%s %s\n'%(self.nx_im1,self.nx_im2))
        inif.write('%s %s %s \n'%(self.y_start,self.y_end,self.ny_step)) #ymin ymax stepy
        inif.write('%s %s %s \n'%(self.x_start,self.x_end,self.nx_step)) #xmin xmax stepx
        inif.write('%s %s \n'%(self.nx_win,self.ny_win)) # window size x y
        inif.write('%s %s \n'%(self.nx_search,self.ny_search)) # search in x y
        inif.write('1 1 \n') # oversampling factor x y
        inif.write('%s %s \n'%(self.x0_shift,self.y0_shift)) # initial offset x y
        inif.write('0. 1.e10\n') # snr threshold, covariance threshold (as set we take everything and do not filter) 
        inif.write('f f\n') # Debug and Display Flags T/F
        inif.close()

        ny_off = np.int32( (self.y_end-self.y_start)/self.ny_step )
        nx_off = np.int32( (self.x_end-self.x_start)/self.nx_step ) # x from self.nx_win/2 to self.nx_im1-self.nx_win/2 with step = nx_step

        if self.off.IsUsed and self.off.data.vx is None:
            self.p_offmap_out()
        #}}}
        # initialization input, output {{{
        ny_off = np.int32( (self.y_end-self.y_start)/self.ny_step)
        nx_off = np.int32( (self.x_end-self.x_start)/self.nx_step) # x from self.nx_win/2 to self.nx_im1-self.nx_win/2 with step = nx_step
        
        f = open(self.output_file,'w')
        # }}}
        for j_off in range(0,ny_off):# loop over y-direction (azimuth for SAR) {{{
            
            j_im = j_off*self.ny_step+self.y_start

            for i_off in range(0,nx_off): 
                 
                i_im = np.int32(i_off*self.nx_step+self.x_start) # from output offset x-coord to im x-coord 

                #xoff = np.int32((i_im - self.off.x_start) / self.off.rgsp)
                #yoff = np.int32( (j_im - self.off.y_start) / self.off.azsp)

                if self.off.data.snr[ j_off , i_off ] !=0:
                    f.write(' '+'{:7d}'.format(np.int32(i_im))+' '+'{:9.3f}'.format(self.off.data.vx[ j_off , i_off ])+\
                            ' '+'{:7d}'.format(np.int32(j_im))+' '+'{:11.3f}'.format(self.off.data.vy[ j_off , i_off ])+\
                            ' '+'{:10.5f}'.format(self.off.data.snr[ j_off , i_off ])+ \
                            ' '+'{:10.6f}'.format(0)+ \
                            ' '+'{:10.6f}'.format(0)+ \
                            ' '+'{:10.6f}'.format(0)+'\n' )
    			
        f.close

        if binary:
            self.off.write(self.output_file+'.par')
            off = np.float32(self.off.data.vx) + 1j * np.float32(self.off.data.vy)
            off.byteswap().tofile(self.output_file+'.off')

    # }}}   
    #}}}
    def off4shiftmap(self,size=9,thre=1): # {{{
        from median_filter_off import median_filter_off
        from scipy.interpolate import griddata
        from scipy.linalg import lstsq
        from scipy.signal import medfilt

        self.shiftmap = off_param()
        self.shiftmap.r0 = self.off.r0
        self.shiftmap.z0 = self.off.z0
        self.shiftmap.x_start = self.off.x_start
        self.shiftmap.x_end = self.off.x_end
        self.shiftmap.npix = self.off.npix
        self.shiftmap.rgsp = self.off.rgsp
        self.shiftmap.y_start = self.off.y_start
        self.shiftmap.y_end = self.off.y_end
        self.shiftmap.nrec = self.off.nrec
        self.shiftmap.azsp = self.off.azsp
        self.shiftmap.ofw_w = self.off.ofw_w
        self.shiftmap.ofw_h = self.off.ofw_h
        self.shiftmap.ofw_thr = self.off.ofw_thr
        
        self.use_shiftmap = True
        self.use_polyshift = True
        self.use_mask = True

        off = self.off.data.vx + (self.off.data.vy)*1j
        mask = median_filter_off(off,size=size,thre=thre)
        off[mask]= 0 + 1j * 0
        mask = mask | (medfilt(off.real,3) == 0) | (medfilt(off.imag,3) == 0)
        off[mask]= 0 + 1j * 0
        off.real = medfilt(off.real,3)
        off.imag = medfilt(off.imag,3)
        mask = mask | (off.real == 0) | (off.imag ==0)
        # mask is used to indicate where to process (1=process, 0=no process) 
        self.mask = off_param()
        
        self.mask.data = mask.astype(np.uint8)

        #import matplotlib.pyplot as plt
        #f, axarr = plt.subplots(2,2)
        #axarr[0,0].imshow((off.real == 0))
        #axarr[0,1].imshow((off.imag == 0))
        #axarr[1,1].imshow(self.off.data.vy)
        #axarr[1,0].imshow(mask)
        #plt.show()

        self.mask.r0 = self.off.r0
        self.mask.z0 = self.off.z0
        self.mask.x_start = self.off.x_start
        self.mask.x_end = self.off.x_end
        self.mask.npix = self.off.npix
        self.mask.rgsp = self.off.rgsp
        self.mask.y_start = self.off.y_start
        self.mask.y_end = self.off.y_end
        self.mask.nrec = self.off.nrec
        self.mask.azsp = self.off.azsp
        self.mask.ofw_w = self.off.ofw_w
        self.mask.ofw_h = self.off.ofw_h
        self.mask.ofw_thr = self.off.ofw_thr
        
        vx_masked = np.ma.array( off.real, mask=mask)
        vy_masked = np.ma.array( off.imag, mask=mask)
         
        x, y = np.mgrid[0:vx_masked.shape[0], 0:vx_masked.shape[1]]
        x1 = x[~vx_masked.mask]
        y1 = y[~vx_masked.mask]
        vx1 = off.real[~vx_masked.mask]
        vy1 = off.imag[~vx_masked.mask]
        off = None
        vx_masked = None
        vy_masked = None
        data = np.c_[x1,y1,vx1]
        A = np.c_[data[:,0], data[:,1], np.ones(data.shape[0])]
        C,_,_,_ = lstsq(A, data[:,2])
        #if self.verbose:
        print('POLYNOM X: ',C[0],'*X+',C[1],'*Y+',C[2])
    
        self.x0_polyshift = [C[2]-C[1]*self.off.x_start/self.off.rgsp-C[0]*self.off.y_start/self.off.azsp,C[1]/self.off.rgsp,C[0]/self.off.azsp] # Y-X are inverted, scale poly from off pixels to slc pixels

        vx1 = vx1 - (C[0]*x1 + C[1]*y1 + C[2])
        
        data = np.c_[x1,y1,vy1]
        A = np.c_[data[:,0], data[:,1], np.ones(data.shape[0])]
        C,_,_,_ = lstsq(A, data[:,2])
        #if self.verbose:
        print('POLYNOM Y: ',C[0],'*X+',C[1],'*Y+',C[2])
        self.y0_polyshift = [C[2]-C[1]*self.off.x_start/self.off.rgsp-C[0]*self.off.y_start/self.off.azsp,C[1]/self.off.rgsp,C[0]/self.off.azsp] # Y-X are inverted, scale poly from off to slc
        vy1 = vy1 - (C[0]*x1 + C[1]*y1 + C[2])
        
        grid_vx = griddata((x1,y1), vx1, (x, y), method='linear')
        grid_vy = griddata((x1,y1), vy1, (x, y), method='linear')
        import matplotlib.pyplot as plt
        f, axarr = plt.subplots(3)
        axarr[0].imshow(mask)
        axarr[1].imshow(self.off.data.vx)
        axarr[2].imshow(self.off.data.vy)
        self.shiftmap.m_xshift = grid_vx
        self.shiftmap.m_yshift = grid_vy
    # }}}
    def off4polyfit(self): # {{{
        from scipy.linalg import lstsq
        
        mask = (self.off.data.vx == 0)
        
        x, y = np.mgrid[0:self.off.data.vx.shape[0], 0:self.off.data.vy.shape[1]]
        x1 = x[~mask]
        y1 = y[~mask]
        vx1 = self.off.data.vx[~mask]
        vy1 = self.off.data.vy[~mask]
        snr1 = self.off.data.snr[~mask]

        index = np.argsort(snr1)
        index = index[ int(len(index)/2): ] # keep second of best snr
        x1 = x1[index]
        y1 = y1[index]
        vx1 = vx1[index]
        vy1 = vy1[index]
            
        off = None
        vx_masked = None
        vy_masked = None
        
        data = np.c_[x1,y1,vx1]
        A = np.c_[data[:,0], data[:,1], np.ones(data.shape[0])]
        C,_,_,_ = lstsq(A, data[:,2])
        #if self.verbose:
        print('POLYNOM X: ',C[0],'*X+',C[1],'*Y+',C[2])
    
        self.x0_polyshift = [C[2]-C[1]*self.off.x_start/self.off.rgsp-C[0]*self.off.y_start/self.off.azsp,C[1]/self.off.rgsp,C[0]/self.off.azsp] # Y-X are inverted, scale poly from off pixels to slc pixels

        vx1 = vx1 - (C[0]*x1 + C[1]*y1 + C[2])
        
        data = np.c_[x1,y1,vy1]
        A = np.c_[data[:,0], data[:,1], np.ones(data.shape[0])]
        C,_,_,_ = lstsq(A, data[:,2])
        #if self.verbose:
        print('POLYNOM Y: ',C[0],'*X+',C[1],'*Y+',C[2])
        self.y0_polyshift = [C[2]-C[1]*self.off.x_start/self.off.rgsp-C[0]*self.off.y_start/self.off.azsp,C[1]/self.off.rgsp,C[0]/self.off.azsp] # Y-X are inverted, scale poly from off to slc
        vy1 = vy1 - (C[0]*x1 + C[1]*y1 + C[2])
        
    # }}}
    def display_offmap(self): # {{{

        import matplotlib.pyplot as plt

        y, x = np.mgrid[0:self.off.data.vx.shape[0], 0:self.off.data.vx.shape[1]]
        x = x * self.off.rgsp + self.off.x_start
        y = y * self.off.azsp + self.off.y_end

        vx = self.off.data.vx #- self.x0_polyshift[0] + self.x0_polyshift[1]*x + self.x0_polyshift[2]*y
        vy = self.off.data.vy #- self.y0_polyshift[0] + self.y0_polyshift[1]*x + self.y0_polyshift[2]*y
        print('X POLYNOM:')
        print(self.x0_polyshift)
        print('Y POLYNOM:')
        print(self.y0_polyshift)
        vx[self.off.data.vx ==0]=0
        vy[self.off.data.vx ==0]=0
        v = np.sqrt( vx**2 + vy**2)
        mask = (v == 0)
        vx[mask]=np.nan
        vy[mask]=np.nan
        v[mask]=np.nan

        import matplotlib.pyplot as plt
        f, axarr = plt.subplots(2)
        im1 = axarr[0].imshow(vx,cmap='gray')
        f.colorbar(im1,ax=axarr[0],orientation='horizontal')

        im2 = axarr[1].imshow(vy,cmap='gray')
        f.colorbar(im2, ax=axarr[1], orientation='horizontal')
        plt.show()

    # }}}
    def display_init(self): # {{{

        import matplotlib.pyplot as plt

        plt.imshow(self.shiftmap.m_xshift)
        plt.show()
        plt.imshow(self.mask.data)
        plt.show()
    # }}}
