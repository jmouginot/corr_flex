from fparam import off_param, isp_param, geo_param
from geocode_im_gdal import geocode_im_gdal
import numpy as np
import scipy.misc
import os

def r_off_optical(id1,id2):

    offin = id1+'-'+id2+'.offmap_1.in'

    f = open(offin,'r')
    tmp = f.readlines()
    f.close

    yposting = np.int32((tmp[4].split()[2]))
    xposting = np.int32((tmp[5].split()[2]))
    r0 = np.int32((tmp[9].split()[0]))
    z0 = np.int32((tmp[9].split()[1]))
    ofw_w = np.int32((tmp[6].split()[0]))
    ofw_h = np.int32((tmp[6].split()[1]))

    print('Offmap parameters')
    print('yposting :',yposting)
    print('xposting :',xposting)
    print('r0 :',r0)
    print('z0 :',z0)

    os.system('rm off.off')

    os.system("ls "+id1+"-"+id2+".offmap_?  | xargs cat | grep -v '*' | awk 'length($0)>80' >> off.off")
    os.system("ls "+id1+"-"+id2+".offmap_??  | xargs cat | grep -v '*' | awk 'length($0)>80' >> off.off")
    
    statinfo = os.stat('off.off')
    if statinfo.st_size < 830:
        print('Too few lines in off.off. File is (almost) empty..')
        return False
    

    #if True:
    try:
        x, offx, y, offy, snr0, tmp1, tmp2, tmp3 = np.loadtxt('off.off',unpack=True)
    except:
        print('Something wrong with off.off (offsetmap)')
        return False

    print(len(x))
    tmp1 = None
    tmp2 = None
    tmp3 = None

    xmin = min(x)
    xmax = max(x)
    ymin = min(y)
    ymax = max(y)

    print('xmin, xmax, ymin, ymax')
    print(xmin,xmax,ymin,ymax)

    npix = np.int32((xmax-xmin)/xposting) + 1
    nrec = np.int32((ymax-ymin)/yposting) + 1

    print('Offset map dimensions: ',npix,nrec)
    
    ## x and y axis are transposed compared to IDL
    off = np.zeros([nrec,npix],dtype=np.complex64)

    x = np.int32((x - xmin) / xposting)
    y = np.int32((y - ymin) / yposting)

    #w2, = np.where(offx == 0.)
    #if np.size(w2) != 0:
    #    offx = offx + 0.00001

    #w2, = np.where(offy == 0.)
    #if np.size(w2) != 0:
    #    offy = offy + 0.00001
    #w2 = None

    off[y,x] = offx + offy*1j
    offx = None
    offy = None

    snr = np.zeros([nrec,npix],dtype=np.float32)
    snr[y,x] = snr0
    x = None
    y = None
    snr0 = None

    poff = off_param()
    poff.r0 = np.int32(r0)
    poff.z0 = np.int32(z0)
    poff.x_start = np.int32(xmin)
    poff.x_end = np.int32(xmin+(npix-1)*xposting)
    poff.npix = np.int32(npix)
    poff.rgsp =  np.int32(xposting)
    poff.y_start = np.int32(ymin)
    poff.y_end = np.int32(ymin+(nrec-1)*yposting)
    poff.nrec = np.int32(nrec)
    poff.azsp = np.int32(yposting)
    poff.ofw_w = np.int32(ofw_w)
    poff.ofw_h = np.int32(ofw_h)
    poff.ofw_thr = 3
    poff.xoff[0] = r0
    poff.yoff[0] = z0

    # xposting, yposting is in pixel not meters
    p1 = isp_param()
    p1.load(id1+'.par')
    poff.nrec_i = np.int32(p1.nrec/yposting)
    poff.npix_i = np.int32(p1.npix/xposting)
    poff.xnlook = xposting
    poff.ynlook = yposting
    poff.rgsp_i = xposting*p1.rgsp
    poff.azsp_i = yposting*p1.azsp

    # load gc_par
    #modification by SJ at 06/20/2019 for compatibility with SENTINEL-1 AMPCOR results
    geo1=geo_param()
    try:
        geo1.load(f=id1+'.DEM_gc_par')
    except:
        print('WARNING: Cannot find {}.DEM_gc_par. Skip loading.'.format(id1))
    
    geo2=geo_param()
    try:
        geo2.load(f=id2+'.DEM_gc_par')
    except:
        print('WARNING: Cannot find {}.DEM_gc_par. Skip loading.'.format(id2))
    

    # adapt to mli
    if geo1.projection == 'UTM':
        geomli1=geo_param(utm=True)
        geomli1.projection_zone=geo1.projection_zone
    else:
        geomli1=geo_param()
        
    ymin = geo1.ymax + (geo1.nrec-1.)*geo1.yposting
    pmli1=isp_param()
   
    # test if id1.mli.par exists {{{
    if not os.path.exists(id1+'.mli.par'):
        pmli1.sensor = p1.sensor
        pmli1.title = p1.title
        pmli1.date = p1.date
        pmli1.npix = np.int32(p1.npix/10.)
        pmli1.nrec = np.int32(p1.nrec/10.)
        pmli1.rgsp = p1.rgsp*10.
        pmli1.azsp = p1.azsp*10.
        pmli1.write(id1+'.mli.par')
        #os.system("multi_look_MLI "+id1+".slf "+id1+".par "+id1+".mli "+id1+".mli.par 10 10")
   
    if os.path.exists(id1+'.mli.par'):
        statinfo = os.stat(id1+'.mli.par')
        if statinfo.st_size < 10:
            pmli1.sensor = p1.sensor
            pmli1.title = p1.title
            pmli1.date = p1.date
            pmli1.npix = np.int32(p1.npix/10.)
            pmli1.nrec = np.int32(p1.nrec/10.)
            pmli1.rgsp = p1.rgsp*10.
            pmli1.azsp = p1.azsp*10.
            pmli1.write(id1+'.mli.par')
            #os.system("multi_look_MLI "+id1+".slf "+id1+".par "+id1+".mli "+id1+".mli.par 10 10")
    # }}}
    pmli1.load(id1+'.mli.par')
    geomli1.npix = pmli1.npix
    geomli1.nrec = pmli1.nrec
    geomli1.posting = pmli1.rgsp
    geomli1.xposting = pmli1.rgsp
    geomli1.yposting = -pmli1.azsp
    geomli1.xmin = geo1.xmin
    geomli1.ymax = ymin - (geomli1.nrec-1.)*geomli1.yposting

    if geo1.projection == 'UTM':
        geomli2=geo_param(utm=True)
        geomli2.projection_zone=geo2.projection_zone
    else:
        geomli2=geo_param()
    ymin = geo2.ymax + (geo2.nrec-1.)*geo2.yposting
    pmli2=isp_param()
    # test if id2.mli.par exists {{{
    if not os.path.exists(id2+'.mli.par'):
        p2 = isp_param()
        p2.load(id2+'.par')
        pmli2.sensor = p2.sensor
        pmli2.title = p2.title
        pmli2.date = p2.date
        pmli2.npix = np.int32(p2.npix/10.)
        pmli2.nrec = np.int32(p2.nrec/10.)
        pmli2.rgsp = p2.rgsp*10.
        pmli2.azsp = p2.azsp*10.
        pmli2.write(id2+'.mli.par')
        #os.system("multi_look_MLI "+id2+".slf "+id2+".par "+id2+".mli "+id2+".mli.par 10 10")

    if os.path.exists(id2+'.mli.par'):
        statinfo = os.stat(id2+'.mli.par')
        if statinfo.st_size < 10:
            p2 = isp_param()
            p2.load(id2+'.par')
            pmli2.sensor = p2.sensor
            pmli2.title = p2.title
            pmli2.date = p2.date
            pmli2.npix = np.int32(p2.npix/10.)
            pmli2.nrec = np.int32(p2.nrec/10.)
            pmli2.rgsp = p2.rgsp*10.
            pmli2.azsp = p2.azsp*10.
            pmli2.write(id2+'.mli.par')
            #os.system("multi_look_MLI "+id2+".slf "+id2+".par "+id2+".mli "+id2+".mli.par 10 10")
    # }}} 
    pmli2.load(id2+'.mli.par')
    geomli2.npix = pmli2.npix
    geomli2.nrec = pmli2.nrec
    geomli2.posting = pmli1.rgsp
    geomli2.xposting = pmli1.rgsp
    geomli2.yposting = -pmli1.azsp
    geomli2.xmin = geo2.xmin
    geomli2.ymax = ymin - (geomli2.nrec-1.)*geomli2.yposting


    p1 = isp_param()
    p1.load(id1+'.par')
    geo1.npix = np.int32(p1.npix/poff.xnlook)
    geo1.nrec = np.int32(p1.nrec/poff.ynlook)
    geo1.posting = p1.rgsp*poff.xnlook
    geo1.xposting = p1.rgsp*poff.xnlook
    geo1.yposting = -p1.azsp*poff.ynlook

    #geomli1.describe()
    print(geomli1.osr_spatial_reference())
    #geo1.describe()
    print(geo1.osr_spatial_reference())

    pwr1 = geocode_im_gdal(id1+'.mli',geolocal=geo1, \
            geoim=geomli1,datatype=np.float32, \
            missing=0.,verbose=True)
    print(np.shape(pwr1))
    print(pwr1.dtype)
    fout = open(id1+'.pwr1','w')
    pwr1.byteswap().tofile(fout)
    fout.close()
    pwr1=None

    pwr2 = geocode_im_gdal(id2+'.mli',geolocal=geo1, \
            geoim=geomli2,datatype=np.float32, \
            missing=0.,verbose=True)
    fout = open(id2+'.pwr2','w')
    pwr2.byteswap().tofile(fout)
    fout.close()
    pwr2=None

    poff.write(id1+'-'+id2+'.offmap.par')
    
    os.system('ln -s '+id1+'-'+id2+'.offmap.par '+id1+'-'+id2+'.offmap.par.interp')
    
    print('Save off ..')
    print('dtype off:',off.dtype)
    fout = open(id1+'-'+id2+'.offmap.off','w')
    off.byteswap().tofile(fout)
    fout.close()
    
    print('Save snr ..')
    fout = open(id1+'-'+id2+'.offmap.snr','w')
    snr.byteswap().tofile(fout)
    fout.close()
    
    print('Save DEM_gc_par.'+id1+'-'+id2+'.offmap.par.interp ..')
    geo1.write(f='DEM_gc_par.'+id1+'-'+id2+'.offmap.par.interp')
   
    # the lines below uses GAMMA software and have been commented
    # files are not needed and just create an visual output (bmp).
    # (JM 5Jun2020)

    #cmd='offset_sub '+id1+'-'+id2+'.offmap.off '+id1+'-'+id2+'.offmap.par '+id1+'-'+id2+'.offmap.off.new'
    #print(cmd)
    #os.system(cmd)

    #cmd = 'rasmph '+id1+'-'+id2+'.offmap.off.new '+str(poff.npix)+' - - - - - - - '+id1+'-'+id2+'.offmap.off.new.bmp' 
    #print(cmd)
    #os.system(cmd)

    #cmd = 'rm '+id1+'-'+id2+'.offmap.off.new'
    #print(cmd)
    #os.system(cmd)

    return True


    

