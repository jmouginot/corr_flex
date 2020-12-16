from fparam import off_param, isp_param, geo_param
import numpy as np
import scipy.misc
import os

def r_off(id1,id2):

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

    poff.write(id1+'-'+id2+'.offmap.par')
    
    print('Save off ..')
    print('dtype off:',off.dtype)
    fout = open(id1+'-'+id2+'.offmap.off','w')
    off.byteswap().tofile(fout)
    fout.close()
    
    print('Save snr ..')
    fout = open(id1+'-'+id2+'.offmap.snr','w')
    snr.byteswap().tofile(fout)
    fout.close()
    
    return True


    

