import os
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
import numpy as np
import h5py
import math

def divisorGenerator(n):
    large_divisors = []
    for i in range(1, int(math.sqrt(n) + 1)):
        if n % i == 0:
            yield i
            if i*i != n:
                large_divisors.append(n / i)

    for divisor in reversed(large_divisors):
        yield divisor

def quickplot(filename):

    sbpltdim =  list(divisorGenerator(len(filename)))

    if len(sbpltdim) > 2:
        print(int(sbpltdim[1]), int(len(filename)/sbpltdim[1]))
        fig, axs = plt.subplots(int(sbpltdim[1]), int(len(filename)/sbpltdim[1]))
    # elif len(sbpltdim) == 2:
    #     print(int(sbpltdim[0]), int(len(filename)/sbpltdim[0]))
    #     fig, axs = plt.subplots(int(sbpltdim[1]), int(len(filename)/sbpltdim[1]))
    else:
        fig, axs = plt.subplots(1, len(filename))

    axind = 0
    for axi in axs.ravel():
        model =  getmodel(filename[axind])

        if "xobs" in model.keys():
            axi.plot(model['xobs'], model['flux_norm'])
            axi.invert_xaxis()
            # set lables 
            axi.set_title(f'{axind+1}')
            plt.xlabel(r'$\Delta v (v/v_{th})$')
            plt.ylabel(r'$F_{tot}/F_{cont}$')
        elif 'x' in model.keys():
            axi.pcolor(model['p'], model['zeta'], model['icont_surface'])
            # set lables 
            axi.set_title(f'{axind+1}')
            plt.xlabel(r'$P$')
            plt.ylabel(r'$\Zeta$')
        else:
            raise NotImplementedError("function note yet mplemented")
        
        axind += 1
        
    plt.show()




def getmodel(filename):

    if not os.path.exists(filename): 
        raise FileExistsError()

    if ".dat" in filename:
        model = get_flux(filename)
    elif "int" in filename:
        model = get_int2d(filename)
    elif "surf" in filename:
        model = get_surface(filename)
    else:
        raise NameError('Invalide file name')
    
    return model


def get_int2d(filename):

    if 'surf' in filename:
        raise Exception('Unsupported file type')

    int2d= {'x' : [], \
            'z' : [], \
            'params' : {'alpha' : None, \
                        'gamma' : None, \
                        'vth_fiducial' : None, \
                        'xic1' : None, \
                        'xnue0' : None, \
                        'xobs' : None},
            'iabs2d' : [], \
            'iemi2d' : [], \
            'int2d' : [], \
            'tau2d' : [], \
            'vn2d' : []}
    
    with h5py.File(filename, "r") as fio:
        a_ray = fio.get('along_ray')
        
        int2d['iabs2d'] = np.array(a_ray.get('iabs2d'))
        int2d['iemi2d'] = np.array(a_ray.get('iemi2d'))
        int2d['int2d'] = np.array(a_ray.get('int2d'))
        int2d['tau2d'] = np.array(a_ray.get('tau2d'))
        int2d['vn2d'] = np.array(a_ray.get('vn2d'))

        Cords = fio.get('coordinates')
    
        int2d['x'] = np.array(Cords.get('x'))
        int2d['z'] = np.array(Cords.get('z'))

        params = fio.get('input_parameters')

        int2d['params']['alpha'] = params.attrs['alpha']
        int2d['params']['gamma'] = params.attrs['gamma']
        int2d['params']['vth_fiducial'] = params.attrs['vth_fiducial']
        int2d['params']['xic1'] = params.attrs['xic1']
        int2d['params']['xnue0'] = params.attrs['xnue0']
        int2d['params']['xobs'] = params.attrs['xobs']


    return int2d

def get_surface(filename):

    if 'int' in filename:
        raise Exception('Unsupported file type')

    surface =  {'p' : [], \
                'zeta' : [], \
                'params' : {'alpha' : None, \
                            'gamma' : None, \
                            'vth_fiducial' : None, \
                            'xic1' : None, \
                            'xnue0' : None, \
                            'xobs' : None},
                'iabs_surface' : [], \
                'icont_surface' : [], \
                'iem_surface' : [], \
                'iemi_surface' : []}
    
    with h5py.File(filename, "r") as fio:

        cord = fio.get('coordinates')

        surface['p'] = np.array(cord.get('p'))
        surface['zeta'] = np.array(cord.get('zeta'))


        surf = fio.get('surfb')

        surface['iabs_surface'] = np.array(surf.get('iabs_surface'))
        surface['icont_surface'] = np.array(surf.get('icont_surface'))
        surface['iem_surface'] = np.array(surf.get('iem_surface'))
        surface['iemi_surface'] = np.array(surf.get('iemi_surface'))

        params = fio.get('input_parameters')

        surface['params']['alpha'] = params.attrs['alpha']
        surface['params']['gamma'] = params.attrs['gamma']
        surface['params']['vth_fiducial'] = params.attrs['vth_fiducial']
        surface['params']['xic1'] = params.attrs['xic1']
        surface['params']['xnue0'] = params.attrs['xnue0']
        surface['params']['xobs'] = params.attrs['xobs']

    return surface

def get_flux(filename):

    flux = {'alpha' : None, \
            'gamma' : None, \
            'xobs' : [], \
            'flux_tot' : [], \
            'flux_cont' : [], \
            'flux_norm' : []} 

    with open(filename, 'r') as fio:
        # read some parameters 
        fio.readline()

        flux['alpha'] = fio.readline()
        flux['gamma'] = fio.readline()
        # skip the column header
        fio.readline()
        # read data from file 
        dataio = np.loadtxt(fio, dtype="float")

    
    flux['xobs'] = dataio[:,0]
    flux['flux_tot'] = dataio[:,1]
    flux['flux_cont'] = dataio[:,2]

    # plot normalised fuls vs xobs 
    if "DEBUG" in filename :
        flux['flux_norm'] =  np.divide(dataio[:,1],dataio[:,2])
    else:
        flux['flux_norm']  = dataio[:,3]

    return flux