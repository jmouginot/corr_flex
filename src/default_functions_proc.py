import os
import json
import numpy as np

def get_defaults(verbose=False):

    server = os.uname()[1]
    if verbose: print('HOSTNAME: ', server
    )
    if server == 'luke':
        path='/bettik/jmougino/'
        sql_host='ucimac_j'
    elif server == 'f-dahu':
        path='/bettik/jmougino/'
        sql_host='ucimac_j'
    elif server == 'f-dahu2':
        path='/bettik/jmougino/'
        sql_host='ucimac_j'
    elif server== 'stallo':
        path='/global/work/pha051/'
        sql_host='ucimac_j'
    elif server== 'astrolabe':
        path='/mnt/data/mouginot/'
        sql_host='ucimac_j'

    return server, path, sql_host


def create_default(region_dir,year,sensor,cycle,server='luke',fast=False,region='mountains'):#{{{

    pwd = os.getcwd()
    if server == 'luke':
        path='/bettik/jmougino/'
        sql_host='ucimac_j'
    elif server == 'dahu':
        path='/bettik/jmougino/'
        sql_host='ucimac_j'
    elif server== 'stallo':
        path='/global/work/pha051/'
        sql_host='ucimac_j'
        


    if sensor == 'SENTINEL-2':
        sensor_d = 'sentinel2'
        type_d   = 'OPTICAL'
        repo     = path+'SENTINEL2_REPOSITORY/'

    elif sensor == 'LANDSAT-8':
        sensor_d = 'landsat'
        type_d   = 'OPTICAL'
        repo     = path+'LANDSAT_REPOSITORY/'

    elif sensor == 'venus':
        sensor_d = 'venus'
        type_d   = 'OPTICAL'
        if region_dir == 'HIMALAYA':
            repo='/bettik/millanr/VENUS_REPOSITORY/KHUMBU/'
        if region_dir == 'NEW_ZEALAND':
            repo='/bettik/millanr/VENUS_REPOSITORY/TASMAN/'
    else:
        sensor_d=''
        type_d='OPTICAL'
    
   
    if (np.uint32(cycle)<=35):
        nx_search=4
        ny_search=4
    elif (np.uint32(cycle)>35 and np.int32(cycle)<=80):
        nx_search=8
        ny_search=8
    elif (np.uint32(cycle)>80 and np.int32(cycle)<=160):
        nx_search=16
        ny_search=16
    elif (np.uint32(cycle)>160 and np.int32(cycle)<=240):
        nx_search=24
        ny_search=24
    elif (np.uint32(cycle)>240 and np.int32(cycle)<=320):
        nx_search=32
        ny_search=32
    elif (np.uint32(cycle)>320):
        nx_search=48
        ny_search=48
    
    #Fast moving glacier{{{
    if fast:
        if (np.uint32(cycle)<=35):
            nx_search=16
            ny_search=16
        elif (np.uint32(cycle)>35 and np.int32(cycle)<=80):
            nx_search=32
            ny_search=32
        elif (np.uint32(cycle)>80 and np.int32(cycle)<=160):
            nx_search=64
            ny_search=64
        elif (np.uint32(cycle)>160):
            nx_search=128
            ny_search=128
    #}}}

    nx_win = 16
    ny_win = 16

    if sensor == 'LANDSAT-8':
        nx_step = 3
        ny_step = 3
    
    else:
        nx_step = 5
        ny_step = 5

    if region == 'greenland' or region == 'antarctic': # {{{
        if (np.uint32(cycle)<=35):
            nx_search=4
            ny_search=4
        elif (np.uint32(cycle)>35 and np.int32(cycle)<=80):
            nx_search=8
            ny_search=8
        elif (np.uint32(cycle)>80 and np.int32(cycle)<=160):
            nx_search=16
            ny_search=16
        elif (np.uint32(cycle)>160):
            nx_search=32
            ny_search=32
        
        nx_win = 32
        ny_win = 32

        nx_step = 10
        ny_step = 10
    #}}}
    def_param_j = { \
            'sensor':sensor_d,\
            'type':type_d,\
            'region':region_dir,\
            'year':year,\
            'path':{\
                'local':pwd,\
                'repo':repo,\
                'datab':pwd+'/MOSAIC/'},\
                'repeat_cycle':cycle,\
                'nx_search':nx_search,\
                "ny_search":ny_search, \
            'nx_win':nx_win, \
            'ny_win':ny_win, \
            'nx_step':nx_step, \
            'ny_step':ny_step, \
            'sql_host':sql_host, \
            'RegionShp':path+region_dir+'/'+region_dir+'_RoughPolygon.shp' }

    lun = open('defaults_STvar.json','w')
    print(path+region_dir+'/'+region_dir+'_RoughPolygon.shp') 
    json.dump(def_param_j,lun,indent=4,sort_keys=True)
    lun.close()
#}}}
def init_var_proc(filename='defaults_STvar.json'): #{{{
    
    
    if os.path.exists(filename):
        print()
        print(os.getcwd())
        print()
        out = open(filename)
        defaults = json.load(out)
        
        
        #a bit of trimming on the strings
        if defaults['path']['local'][-1]!='/':
            defaults['path']['local']+='/'

        if defaults['path']['repo'][-1]!='/':
            defaults['path']['repo']+='/'

        if defaults['path']['datab'][-1]!='/':
            defaults['path']['datab']+='/'

    else:
        print('Cannot find {}. Making blank defaults JSON.')
        defaults = { \
                'sensor':'',\
                'type':'',\
                'region_dir':'',\
                'year':'',\
                'path':{\
                    'local':'',\
                    'repo':'',\
                    'datab':''},\
                'repeat_cycle':0}

    return defaults
#}}}
def find_in_list(listjson): # {{{

    # RETURNS: list elements with all parameters (ID1, ID2 etc ...) and associated DEFAULTS processing parameters

    lun=open(listjson,'r')
    lst=json.load(lun)
    lun.close()
    pwd=os.getcwd()+'/'

    found_lste=False
    for i in range(0,np.size(lst)):
        if lst[i]['directory']==pwd:
            lste=lst[i]
            ind = i
            found_lste=True

    if not found_lste:
        print(pwd,' not found ..')
        return False, 0, 0, 0

    if os.path.exists('../../defaults_STvar.json'):
        defaults = init_var_proc('../../defaults_STvar.json')
    elif os.path.exists('../defaults_STvar.json'):
        defaults = init_var_proc('../defaults_STvar.json')
    elif os.path.exists('../../../defaults_STvar.json'):
        defaults =init_var_proc('../../../defaults_STvar.json')
    else:
        print('Could not find defaults_STvar.json ..')
        return False, 0, 0, 0

    return True, lste, ind, defaults
# }}}
def copy_slc2irods(lste,id12): # {{{

    irod_dir = '/cigri/home/jmougino'+lste['directory']
    print('irod_dir:',irod_dir)
    print('imkdir -p '+irod_dir)
    status=os.system('imkdir -p '+irod_dir)

    print('iput -f '+id12+'_sobel.slc '+irod_dir)
    status=status+os.system('iput -f '+id12+'.slc '+irod_dir)

    if status == 0 : #status = 0 means all commands worked
        os.system('type ils > irods.is.in.use')
        print('Deleting :',id12+'.slc')
        os.remove(id12+'.slc')

# }}}
def check_if_binary_is_good(binary_file,par_file,datatype='complex',file_type='isp'):

    from fparam import isp_param, off_param, geo_param

    if datatype == 'complex':
        size_data=4*2
    elif datatype == 'float':
        size_data=4
    elif datatype == 'float32':
        size_data=4
    elif datatype == 'integer':
        size_data=2
    elif datatype == 'byte':
        size_data=1
    else:
        size_data=4*2

    slc1_ready=True

    if os.path.exists(par_file):

        if os.path.getsize(par_file) == 0:
            return False

        if file_type == 'isp':
            p1=isp_param()
        elif file_type == 'off':
            p1=off_param()
        elif file_type == 'geo':
            p1=geo_param()
        else:
            p1=isp_param()

        p1.load(par_file)

    else:

        print('ERROR1: '+par_file+' does not exist.')
        slc1_ready=False

    if os.path.exists(binary_file):
        if not os.path.getsize(binary_file) == p1.nrec*p1.npix*size_data:
            print('ERROR2: '+binary_file+' has not the good size.')
            slc1_ready=False

    else:

        print('ERROR2: '+binary_file+' does not exist.')
        slc1_ready=False

    return slc1_ready
