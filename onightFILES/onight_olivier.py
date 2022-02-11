import sys, os
import numpy as np
import const_cgs
from lev_model_laws import *


def calc_orbit_phase(mass1, mass2, eccentricity, period, phase):

    #on input: masses in g
    #          period in seconds
    #          phase in rad

    #calculates 2-body orbit as a function of phase of 2-body problem with eccentric orbits

    #calculate  semi-latus rectum
    psep = (1.-eccentricity**2)*(const_cgs.cgs_grav*(mass1+mass2)*period**2. /4./np.pi**2)**(1./3.)

    #separation of both objects
    rsep = psep/(1.+eccentricity*np.cos(phase))

    
    #derivatives
    drdphi = rsep*eccentricity*np.sin(phase) / (1.+eccentricity*np.cos(phase))
    dphidt = np.sqrt(const_cgs.cgs_grav * (mass1+mass2) * psep)/rsep**2

    #calculate actual position in reduced frame
    x00 = np.cos(phase)*rsep
    y00 = np.sin(phase)*rsep

    #calculate the velocities in reduced frame
    vx00 = drdphi*dphidt * np.cos(phase) - rsep*dphidt*np.sin(phase)
    vy00 = drdphi*dphidt * np.sin(phase) + rsep*dphidt*np.cos(phase)

    #calculate the positions in the rest-frame (center-of-mass-frame)
    x01=mass2/(mass2+mass1)*x00
    y01=mass2/(mass2+mass1)*y00
    x02=-mass1/(mass2+mass1)*x00
    y02=-mass1/(mass2+mass1)*y00

    #calculate the velocities in the rest-frame (center-of-mass-frame)
    vx01=mass2/(mass2+mass1)*vx00
    vy01=mass2/(mass2+mass1)*vy00
    vx02=-mass1/(mass2+mass1)*vx00
    vy02=-mass1/(mass2+mass1)*vy00
    
    return {'pos1': np.array([x01, y01]), 'vel1': np.array([vx01, vy01]), 'pos2': np.array([x02, y02]), 'vel2': np.array([vx02, vy02])}


def get_mass2(mass1, massf, incl):
    #calculate the mass of the secondary
    #mass 1 and massf are mass of primary and mass-function
    #incl is inclination

    tol=1.e-7
    #try fix-point iteration
    acoeff = np.sin(inclination)
    acoeff = acoeff**3.
    bcoeff = -massf
    ccoeff = -2.*mass1*massf
    dcoeff = -mass1**2 * massf

    #find mass2 by secant method
    mass2_lower = 0.0001
    mass2_upper = 10.
    fct_lower = acoeff*mass2_lower**3 + bcoeff*mass2_lower**2 + ccoeff*mass2_lower + dcoeff
    fct_upper = acoeff*mass2_upper**3 + bcoeff*mass2_upper**2 + ccoeff*mass2_upper + dcoeff

    for i in range(0,1000):
        mass2_new = mass2_lower - (mass2_upper-mass2_lower)/(fct_upper-fct_lower)*fct_lower
        fct_new = acoeff*mass2_new**3 + bcoeff*mass2_new**2 + ccoeff*mass2_new + dcoeff
        if fct_new*fct_lower > 0:
            mass2_lower=mass2_new
            fct_lower=fct_new
        elif fct_new*fct_upper > 0:
            mass2_upper=mass2_new
            fct_upper=fct_new
        else:
            print('error in secant method')
            exit()

        if np.abs(fct_new) < tol:
            return mass2_new


def number_to_string(number, dtype):
    
    #sign
    str_sign=''
    if number < 0.: str_sign='-'

    number=np.abs(number)

    if dtype == 'dp':
        if number == 0.:
            str_all = '0.d0'            
        elif number > 1.:
            iexp=0
            fdum=1.
            for i in np.arange(1,100):
                if number/fdum < 10.: break
                fdum = fdum*10.
                iexp = iexp+1
            number = number/fdum
            str_all = str_sign + '{number:}'.format(number=number) + 'd{iexp:}'.format(iexp=iexp)
        elif number < 1.:
            iexp = 0
            fdum=1.
            for i in np.arange(1,100):
                if number*fdum > 0.1: break
                fdum = fdum*10.
                iexp = iexp-1
            number = number*fdum*10.  #this weird formalism only to have a nice representation of numbers after the comma
            iexp=iexp-1
            str_all = str_sign + '{number:}'.format(number=number) + 'd{iexp:}'.format(iexp=iexp)
        else:
            number = '1.'
            iexp = 0
            str_all = str_sign + '1.d0'

    elif dtype == 'int':
        str_all = str_sign + '{number:}'.format(number=number)
    else:
        exit('error in number_to_string: dtype not properly specified')
        
    return str_all

def numbers_to_string(numbers, dtype):
    str_all=[]
    for number in numbers:
        number=number_to_string(number,dtype)
        str_all.append(number)
    return str_all

def get_model(model,
              inclination=None,
              mass1=None,
              massf=None,
              period=None,
              eccentricity=None,
              phase=None,
              kline=1.):

    
    #calculate mass of secondary
    mass2 = get_mass2(mass1, massf, inclination)

    #masses in cgs
    mass1 = mass1*const_cgs.cgs_msu
    mass2 = mass2*const_cgs.cgs_msu
    massf = massf*const_cgs.cgs_msu

    #calculate orbit at a given phase
    data = calc_orbit_phase(mass1, mass2, eccentricity, period, phase)

    p_object01 = [data['pos1'][0]/model['unit_length']/const_cgs.cgs_rsu, data['pos1'][1]/model['unit_length']/const_cgs.cgs_rsu, 0.]
    p_object02 = [data['pos2'][0]/model['unit_length']/const_cgs.cgs_rsu, data['pos2'][1]/model['unit_length']/const_cgs.cgs_rsu, 0.]

    v_object01 = [data['vel1'][0]/1.e5, data['vel1'][1]/1.e5, 0.]
    v_object02 = [data['vel2'][0]/1.e5, data['vel2'][1]/1.e5, 0.]
    

    model['p_object01']=p_object01
    model['v_object01']=v_object01
    model['p_object02']=p_object02
    model['v_object02']=v_object02

    model['kline'] = kline

    return model
#
#--------------------------------------------------------------------------
#
def indat_modelspec(model, dir_base=None, fname='indat_modelspec_olivier.nml'):

    rstar1 = number_to_string(model['rstar1'], 'dp')
    rstar2 = number_to_string(model['rstar2'], 'dp')
    rmin1 = number_to_string(model['rmin1'], 'dp')
    rmin2 = number_to_string(model['rmin2'], 'dp')
    rmax1 = number_to_string(model['rmax1'], 'dp')
    rmax2 = number_to_string(model['rmax2'], 'dp')
    teff1 = number_to_string(model['teff1'], 'dp')
    teff2 = number_to_string(model['teff2'], 'dp')
    trad1 = number_to_string(model['trad1'], 'dp')
    trad2 = number_to_string(model['trad2'], 'dp')
    logg1 = number_to_string(model['logg1'], 'dp')
    logg2 = number_to_string(model['logg2'], 'dp')
    yhe1 = number_to_string(model['yhe1'], 'dp')
    yhe2 = number_to_string(model['yhe2'], 'dp')
    fehe1 = number_to_string(model['fehe1'], 'dp')
    fehe2 = number_to_string(model['fehe2'], 'dp')
    aenh1 = number_to_string(model['aenh1'], 'dp')
    aenh2 = number_to_string(model['aenh2'], 'dp')       
    vrot1 = number_to_string(model['vrot1'], 'dp')
    vrot2 = number_to_string(model['vrot2'], 'dp')
    vmicro1 = number_to_string(model['vmicro1'], 'dp')
    vmicro2 = number_to_string(model['vmicro2'], 'dp')

    kline = number_to_string(model['kline'], 'dp')    

    p_object01 = numbers_to_string(model['p_object01'], 'dp')
    p_object02 = numbers_to_string(model['p_object02'], 'dp')
    v_object01 = numbers_to_string(model['v_object01'], 'dp')        
    v_object02 = numbers_to_string(model['v_object02'], 'dp')


    ex01 = numbers_to_string(model['ex01'], 'dp')
    ey01 = numbers_to_string(model['ey01'], 'dp')
    ez01 = numbers_to_string(model['ez01'], 'dp')        
    rot_axis01 = numbers_to_string(model['rot_axis01'], 'dp')
    ex02 = numbers_to_string(model['ex02'], 'dp')
    ey02 = numbers_to_string(model['ey02'], 'dp')
    ez02 = numbers_to_string(model['ez02'], 'dp')        
    rot_axis02 = numbers_to_string(model['rot_axis02'], 'dp')

    unit_length = number_to_string(model['unit_length'], 'dp')
    vth_fiducial = number_to_string(model['vth_fiducial'], 'dp')


    model_unitl = number_to_string(model['model_unitl'], 'dp')
    model_unitd = number_to_string(model['model_unitd'], 'dp')
    model_unitv = number_to_string(model['model_unitv'], 'dp')
    model_unitt = number_to_string(model['model_unitt'], 'dp')
    cs2_nphi = number_to_string(model['cs2_nphi'], 'int')
    tjet = number_to_string(model['tjet'], 'dp')    

    
    file = open(fname, "w")
    file.write("&input_options\n")
    file.write("input_file = ''\n")
    file.write("input_file2 = ''\n")
    if not dir_base:
        file.write("output_file = './outputFILES/models_olivier/modspec_model00.h5'\n")
    else:
        command00 = 'mkdir ' + dir_base
        os.system(command00)
        file.write("output_file = '" + dir_base + "modspec_model00.h5'\n")
        
    file.write("input_mod = 9 /\n")
    file.write("\n")


    file.write("&input_model1\n")
    file.write("rstar1 = " + rstar1 + "\n")
    file.write("rmin1 = " + rmin1 + "\n")
    file.write("rmax1 = " + rmax1 + "\n")
    file.write("teff1 = " + teff1 + "\n")
    file.write("trad1 = " + trad1 + "\n")
    file.write("logg1 = " + logg1 + "\n")
    file.write("yhe1  = " + yhe1 + "\n")
    file.write("fehe1 = " + fehe1 + "\n")
    file.write("aenh1 = " + aenh1 + "\n")
    file.write("vrot1 = " + vrot1 + "\n")
    file.write("vmicro1 = " + vmicro1 + "\n")
    file.write("p_object01 = " + p_object01[0] + ", " + p_object01[1] + ", " + p_object01[2] + "\n")
    file.write("v_object01 = " + v_object01[0] + ", " + v_object01[1] + ", " + v_object01[2] + "\n")
    file.write("ex01 = " + ex01[0] + ", " + ex01[1] + ", " + ex01[2] + "\n")
    file.write("ey01 = " + ey01[0] + ", " + ey01[1] + ", " + ey01[2] + "\n")
    file.write("ez01 = " + ez01[0] + ", " + ez01[1] + ", " + ez01[2] + "\n")
    file.write("rot_axis01 = " + rot_axis01[0] + ", " + rot_axis01[1] + ", " + rot_axis01[2] + "\n")        
    file.write("/\n")
    file.write("\n")
    
    file.write("&input_model2\n")
    file.write("rstar2 = " + rstar2 + "\n")
    file.write("rmin2 = " + rmin2 + "\n")
    file.write("rmax2 = " + rmax2 + "\n")
    file.write("teff2 = " + teff2 + "\n")
    file.write("trad2 = " + trad2 + "\n")
    file.write("logg2 = " + logg2 + "\n")
    file.write("yhe2  = " + yhe2 + "\n")
    file.write("fehe2 = " + fehe2 + "\n")
    file.write("aenh2 = " + aenh2 + "\n")
    file.write("vrot2 = " + vrot2 + "\n")
    file.write("vmicro2 = " + vmicro2 + "\n")
    file.write("p_object02 = " + p_object02[0] + ", " + p_object02[1] + ", " + p_object02[2] + "\n")
    file.write("v_object02 = " + v_object02[0] + ", " + v_object02[1] + ", " + v_object02[2] + "\n")
    file.write("ex02 = " + ex02[0] + ", " + ex02[1] + ", " + ex02[2] + "\n")
    file.write("ey02 = " + ey02[0] + ", " + ey02[1] + ", " + ey02[2] + "\n")
    file.write("ez02 = " + ez02[0] + ", " + ez02[1] + ", " + ez02[2] + "\n")
    file.write("rot_axis02 = " + rot_axis02[0] + ", " + rot_axis02[1] + ", " + rot_axis02[2] + "\n")        
    file.write("/\n")
    file.write("\n")
    
    file.write("&input_line\n")
    file.write("iline=1\n")
    file.write("eps_line = 0.d0\n")
    file.write("kline = " + kline + "\n")
    file.write("/\n")
    
    file.write("iline=1: halpha\n")
    file.write("iline=2: hbeta\n")
    file.write("iline=10: civ resonance line\n")
    file.write("\n")
    file.write("\n")
    
    file.write("&input_units\n")
    file.write("unit_length = " + unit_length + "\n")
    file.write("vth_fiducial = " + vth_fiducial + "\n")    
    file.write("/\n")
    file.write("\n")
    
    file.write("&input_usr\n")
    file.write("fname_model = '" + model['fname_model'] + "'\n")    
    file.write("model_unitl = " + model_unitl + "\n")
    file.write("model_unitd = " + model_unitd + "\n")
    file.write("model_unitv = " + model_unitv + "\n")
    file.write("model_unitt = " + model_unitt + "\n")
    file.write("cs2_nphi = " +   cs2_nphi + "\n")
    file.write("tjet = " + tjet + "\n")        
    file.write("/\n")
    file.write("\n")
    
    file.write("&dum\n")
    file.write("/")
    
    file.close
#
#---------------------------------------------------------------------
#
def indat_spec(fname='indat_spec_olivier.nml', dir_base=None):


    file = open(fname, "w")
    file.write("&input_options\n")
    file.write("input_mod=2\n")
    if not dir_base:
        file.write("input_file='outputFILES/models_olivier/modspec_model00.h5'\n")
        file.write("output_dir='outputFILES/models_olivier'\n")
    else:
        file.write("input_file = '" + dir_base + "modspec_model00.h5'\n")
        file.write("output_dir = '" + dir_base[:-1] + "'\n")
        

    file.write("opt_photprof1=0 \n")
    file.write("opt_photprof2=0 \n")
    file.write("opt_obsdir_read=t\n")
    file.write("opt_surface=f\n")
    file.write("opt_int2d=f\n")
    file.write("opt_incl_gdark1=f\n")
    file.write("opt_incl_sdist1=f\n")
    file.write("opt_incl_gdark2=f\n")
    file.write("opt_incl_sdist2=f\n")
    file.write("opt_pgrid01 = 'log'\n")
    file.write("opt_rgrid01 = 'log'\n")
    file.write("opt_pgrid02 = 'lin'\n")
    file.write("opt_rgrid02 = 'lin'\n")
    file.write("nalpha=1\n")
    file.write("ngamma=1 /\n")
    file.write("\n")

    file.write("input_mod=0   1D input model\n")
    file.write("input_mod=1   3d input model (cartesian coordinates)\n")
    file.write("input_mod=2   3d input model (spherical coordinates)\n")
    file.write("\n")

    file.write("&input_model\n")
    file.write("vth_fiducial = 1.d2\n")
    file.write("/\n")
    file.write("\n")
    file.write("\n")

    file.write("&input_surface\n")
    file.write("alpha_surface=1.57d0\n")
    file.write("gamma_surface=3.35103 !3.14d0\n")
    file.write("xobs_surface=-2.d0 /\n")
    file.write("0.785d0\n")
    file.write("\n")
    file.write("\n")

    file.write("&dum\n")
    file.write("/")
    
    file.close
#
#---------------------------------------------------------------------
#
def copy_to_directory(dir_base,
                      indx=0):
    
    dir_extension =   'phase{indx:03d}'.format(indx=indx) 

    dir_out = dir_base+dir_extension

    command01 = 'mkdir ' + dir_out
    command02 = 'mv ' + dir_base + 'FLUXEM*.dat ' + dir_out
    command03 = 'mv ' + dir_base + 'photprof*.dat ' + dir_out
    command04 = 'mv ' + dir_base + 'spec_points*.dat ' + dir_out
    command05 = 'mv ' + dir_base + 'spec_triangles*.dat ' + dir_out
    command06 = 'mv ' + dir_base + 'modspec_model00.h5 ' + dir_out
    command07 = 'mv ./indatFILES/indat_modelspec_olivier.nml ' + dir_out
    command08 = 'mv ./indatFILES/indat_spec_olivier.nml ' + dir_out
    command09 = 'mv ./output_modelspec.log ' + dir_out
    command10 = 'mv ./output_spec.log ' + dir_out

    os.system(command01)
    os.system(command02)
    os.system(command03)
    os.system(command04)
    os.system(command05)
    os.system(command06)
    os.system(command07)
    os.system(command08)
    os.system(command09)
    os.system(command10)
#
#---------------------------------------------------------------------
#
def calc_all(dir_base, k, nd):

    print('--------------------running model {k:} / {nd:}---------------------'.format(k=k,nd=nd))
    print('')

    print('-----------------running modelspec_vbin.eo-------------------')
    print('')

    os.system("./modelspec_vbin.eo < in_modelspec_vbin > output_modelspec.log")


    print('--------------------running spec_vbin.eo---------------------')
    print('')

    os.system("./spec_vbin.eo < in_spec_vbin > output_spec.log")

    print('------------------copy everything to directory---------------')
    print('')

    copy_to_directory(dir_base, k)
#
#---------------------------------------------------------------------
#----------------------model p05s10-----------------------------------
#---------------------------------------------------------------------
#
inclination = 0.942477796077
mass1 = 0.6
massf = 0.274
period = 1288.6*24.0*3600.0
eccentricity = 0.23

model={'rstar1': 0.59 * const_cgs.cgs_au/const_cgs.cgs_rsu,
       'rmin1': 1.0,
       'rmax1': 1.1,
       'teff1': 6000.,
       'trad1': 6000.,
       'logg1': 1.,
       'yhe1': 0.1,
       'fehe1': -1.,
       'aenh1': 0.,
       'vrot1': 10./np.sin(inclination),
       'vmicro1': 20.,
       'p_object01': None,
       'v_object01': None,
       'ex01': [1., 0., 0.],
       'ey01': [0., 1., 0.],
       'ez01': [0., 0., 1.],
       'rot_axis01': [0., 0., 1.],
       'rstar2': 2.8754468919386458,
       'rmin2': 1.,
       'rmax2': 430.,
       'teff2': 100.,
       'trad2': 100.,
       'logg2': 3.,
       'yhe2': 0.,
       'fehe2': -1.,
       'aenh2': 0.,
       'vrot2': 0.,
       'vmicro2': 20.,
       'p_object02': None,
       'v_object02': None,
       'ex02': [1., 0., 0.],
       'ey02': [0., 1., 0.],
       'ez02': [0., 0., 1.],
       'rot_axis02': [0., 0., 1.],       
       'unit_length': 0.59 * const_cgs.cgs_au/const_cgs.cgs_rsu,
       'vth_fiducial': 100.,
       'fname_model': 'models/olivier/model1143.h5',
       'model_unitl': 1.49598e13,
       'model_unitd': 1.,
       'model_unitv': 1.e5,
       'model_unitt': 1.,
       'cs2_nphi': 101,
       'tjet': 5000.}

#Argument of periapsis [deg].
omega = 27.4 #297.4

#to get transit, set actual phases:
phases = np.linspace(omega+180.-90., omega+180.+90., 17)
phases = np.linspace(0.,360., 33)
phases = phases*np.pi/180.
#phases = np.array([omega-180.])*np.pi/180.

#phases = np.linspace(np.pi/2,3./2.*np.pi, 17)
nphases = len(phases)
indx = 1

dir_base = '/lhome/levin/Postdoc/line3D/outputFILES/models_olivier/'

kline=1.e0
dir_base = '/lhome/levin/Postdoc/line3D/outputFILES/models_olivier/photprof0_sline0_vmicro20_kline1d0_test/'



for phase in phases:

    #calculate the orbital parameters for given phase
    model = get_model(model,
                      inclination=inclination,
                      mass1=mass1,
                      massf=massf,
                      period=period,
                      eccentricity=eccentricity,
                      phase=phase,
                      kline=kline)
#    model['v_object01'] = np.zeros(3)
#    model['v_object02'] = np.zeros(3)          
    print(phase*180./np.pi, model['p_object01'], model['p_object02'], model['v_object02'])

    #set the indat file for modelspec
    indat_modelspec(model,
                    fname='indatFILES/indat_modelspec_olivier.nml',
                    dir_base=dir_base)


    #set the indat file for spec
    indat_spec(fname='indatFILES/indat_spec_olivier.nml', dir_base=dir_base)

    #calculate everything and copy to directories
    calc_all(dir_base, indx, nphases)
    indx=indx+1
