import sys, os
import numpy as np
from const_cgs import *
from lev_model_laws import *

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
#
#---------------------------------------------------------------------
#
def indat_modelspec(fname='indat_modelspec_lb1.nml',
                    rstar1=None,
                    rstar2=None,
                    rmax1=None,
                    rmax2=None,
                    teff1=None,
                    teff2=None,
                    trad1=None,
                    trad2=None,
                    p_object01=None,
                    p_object02=None,
                    v_object01=None,
                    v_object02=None,              
                    opt_temperature1=None,
                    opt_temperature2=None,
                    opt_vrot1=None,
                    opt_vrot2=None,
                    vmicro1=None,
                    vmicro2=None,
                    vrot1=None,
                    vrot2=None,
                    logg1=None,
                    logg2=None,
                    tdisc1=None,
                    tdisc2=None,
                    mdisc1=None,
                    mdisc2=None,
                    rdisc1_max=None,
                    rdisc2_max=None,
                    slope1=None,
                    slope2=None):

    lerr=False
    if rstar1 == None: lerr=True
    if rstar2 == None: lerr=True
    if rmax1 == None: lerr=True
    if rmax2 == None: lerr=True
    if teff1 == None: lerr=True
    if teff2 == None: lerr=True
    if trad1 == None: lerr=True
    if trad2 == None: lerr=True
    if p_object01 == None: lerr=True
    if p_object02 == None: lerr=True
    if v_object01 == None: lerr=True
    if v_object02 == None: lerr=True
    if opt_temperature1 == None: lerr=True
    if opt_temperature2 == None: lerr=True
    if opt_vrot1 == None: lerr=True
    if opt_vrot2 == None: lerr=True
    if vmicro1 == None: lerr=True
    if vmicro2 == None: lerr=True
    if vrot1 == None: lerr=True
    if vrot2 == None: lerr=True    
    if logg1 == None: lerr=True
    if logg2 == None: lerr=True    
    if tdisc1 == None: lerr=True
    if tdisc2 == None: lerr=True
    if mdisc1 == None: lerr=True
    if mdisc2 == None: lerr=True
    if rdisc1_max == None: lerr=True
    if rdisc2_max == None: lerr=True
    if slope1 == None: lerr=True
    if slope2 == None: lerr=True
    if lerr: exit('error in indat_modelspec: argument not specified')


    rstar1 = number_to_string(rstar1, 'dp')
    rstar2 = number_to_string(rstar2, 'dp')
    rmax1 = number_to_string(rmax1, 'dp')
    rmax2 = number_to_string(rmax2, 'dp')
    teff1 = number_to_string(teff1, 'dp')
    teff2 = number_to_string(teff2, 'dp')
    trad1 = number_to_string(trad1, 'dp')
    trad2 = number_to_string(trad2, 'dp')

    p_object01 = numbers_to_string(p_object01, 'dp')
    p_object02 = numbers_to_string(p_object02, 'dp')
    v_object01 = numbers_to_string(v_object01, 'dp')        
    v_object02 = numbers_to_string(v_object02, 'dp')

    opt_temperature1 = number_to_string(opt_temperature1, 'int')
    opt_temperature2 = number_to_string(opt_temperature2, 'int')
    opt_vrot1 = number_to_string(opt_vrot1, 'int')
    opt_vrot2 = number_to_string(opt_vrot2, 'int')    
    vmicro1 = number_to_string(vmicro1, 'dp')
    vmicro2 = number_to_string(vmicro2, 'dp')
    vrot1 = number_to_string(vrot1, 'dp')
    vrot2 = number_to_string(vrot2, 'dp')
    logg1 = number_to_string(logg1, 'dp')
    logg2 = number_to_string(logg2, 'dp')        
    tdisc1 = number_to_string(tdisc1, 'dp')  
    tdisc2 = number_to_string(tdisc2, 'dp')    
    mdisc1 = number_to_string(mdisc1, 'dp')
    mdisc2 = number_to_string(mdisc2, 'dp')
    rdisc1_max = number_to_string(rdisc1_max, 'dp')
    rdisc2_max = number_to_string(rdisc2_max, 'dp')
    slope1 = number_to_string(slope1, 'dp')
    slope2 = number_to_string(slope2, 'dp')
    
    file = open(fname, "w")
    file.write("&input_options\n")
    file.write("input_file = './outputFILES/output_model00.h5'\n")
    file.write("input_file2 = './inputFILES/model2d.h5'\n")
    file.write("output_file = './outputFILES/modspec_model00.h5'\n")
    file.write("input_mod = 5 /\n")
    file.write("\n")


    file.write("&input_model1\n")
    file.write("rstar1 = " + rstar1 + "\n")
    file.write("rmin1 = 1.d0\n")
    file.write("rmax1 = " + rmax1 + "\n")
    file.write("teff1 = " + teff1 + "\n")
    file.write("trad1 = " + trad1 + "\n")
    file.write("p_object01 = " + p_object01[0] + ", " + p_object01[1] + ", " + p_object01[2] + "\n")
    file.write("v_object01 = " + v_object01[0] + ", " + v_object01[1] + ", " + v_object01[2] + "\n")
    file.write("ex01 = 1.d0, 0.d0, 0.d0\n")
    file.write("ey01 = 0.d0, 1.d0, 0.d0\n")
    file.write("ez01 = 0.d0, 0.d0, 1.d0\n")
    file.write("rot_axis01 = 0.d0, 0.d0, 1.d0\n")
    file.write("/\n")
    
    file.write("teff = 18.d3\n")
    file.write("vy01 = \n")
    file.write("rmax1 = 10.223429999999997d0\n")
    file.write("\n")
    
    file.write("&input_model2\n")
    file.write("rstar2 = " + rstar2 + "\n")
    file.write("rmin2 = 1.d0\n")
    file.write("rmax2 = " + rmax2 + "\n")
    file.write("teff2 = " + teff2 + "\n")
    file.write("trad2 = " + trad2 + "\n")
    file.write("p_object02 = " + p_object02[0] + ", " + p_object02[1] + ", " + p_object02[2] + "\n")
    file.write("v_object02 = " + v_object02[0] + ", " + v_object02[1] + ", " + v_object02[2] + "\n")
    file.write("ex02 = 1.d0, 0.d0, 0.d0\n")
    file.write("ey02 = 0.d0, 1.d0, 0.d0\n")
    file.write("ez02 = 0.d0, 0.d0, 1.d0\n")
    file.write("rot_axis02 = 0.d0, 0.d0, 1.d0\n")
    file.write("/\n")
    
    file.write("rmax2 = 10.d0\n")
    file.write("\n")
    file.write("\n")
    
    file.write("&input_line\n")
    file.write("iline=1\n")
    file.write("eps_line = 0.d0\n")
    file.write("kline = 1.d0\n")
    file.write("/\n")
    
    file.write("iline=1: halpha\n")
    file.write("iline=2: hbeta\n")
    file.write("iline=10: civ resonance line\n")
    file.write("\n")
    file.write("\n")
    
    file.write("&input_units\n")
    file.write("unit_length = 3.7d0\n")
    file.write("/\n")
    file.write("\n")
    
    file.write("&input_usr\n")
    file.write("opt_temperature1 = " + opt_temperature1 + "\n")
    file.write("opt_temperature2 = " + opt_temperature2 + "\n")
    file.write("opt_vrot1 = " + opt_vrot1 + "\n")
    file.write("opt_vrot2 = " + opt_vrot2 + "\n")
    file.write("kcont=1.d0\n")
    file.write("kline=1.d0\n")
    file.write("vmicro1 = " + vmicro1 + "\n")
    file.write("vrot1 = " + vrot1 + "\n")
    file.write("logg1 = " + logg1 + "   !4.77\n")
    file.write("yhe1 = 0.1d0\n")
    file.write("hi1 = 1.d0\n")
    file.write("hei1 = 2.d0\n")
    file.write("lstar1 = 1.d6\n")
    file.write("tdisc1 = " + tdisc1 + "\n")
    file.write("tdisc1a = 2.d3\n")
    file.write("tdisc1b = 2.d3\n")
    file.write("mdisc1 = " + mdisc1 + "\n")
    file.write("rdisc1_max = " + rdisc1_max + "\n")
    file.write("dtheta_disc1 = 45.d0\n")
    file.write("slope1 = " + slope1 + "\n")
    file.write("vmicro2 = " + vmicro2 + "\n")
    file.write("vrot2 = " + vrot2 + "  !11.1d0 !21.4d0\n")
    file.write("logg2 = " + logg2 + "\n")
    file.write("yhe2 = 0.2d0\n")
    file.write("hi2 = 1.d0\n")
    file.write("hei2 = 2.d0\n")
    file.write("lstar2 = 1.d6\n")
    file.write("tdisc2 = " + tdisc2 + "\n")
    file.write("tdisc2a = 12.7d3\n")
    file.write("tdisc2b = 5.d3\n")
    file.write("mdisc2 = " + mdisc2 + "\n")
    file.write("rdisc2_max = " + rdisc2_max + "\n")
    file.write("dtheta_disc2 = 45.d0\n")
    file.write("slope2 = " + slope2 + "\n")
    file.write("vrot2fac=1.d0\n")
    file.write("vth_fiducial=1.d2\n")
    file.write("/\n")
    
    file.write("mdisc2 = 1.78d-12\n")
    file.write("rdisc2_max = 13.5d0\n")
    file.write("\n")
    
    file.write("&dum\n")
    file.write("/")
    
    file.close
#
#---------------------------------------------------------------------
#
def indat_spec(fname='indat_spec_lb1.nml',
               opt_photprof1=None,
               opt_photprof2=None):

    lerr=False
    if opt_photprof1 == None: lerr=True
    if opt_photprof2 == None: lerr=True
    if lerr: exit('error in indat_modelspec: argument not specified')

    opt_photprof1 = number_to_string(opt_photprof1, 'int')
    opt_photprof2 = number_to_string(opt_photprof2, 'int')
    
    file = open(fname, "w")
    file.write("&input_options\n")
    file.write("input_mod=2\n")
    file.write("input_file='outputFILES/modspec_model00.h5'\n")
    file.write("output_dir='outputFILES'\n")
    file.write("opt_photprof1=" + opt_photprof1 + "\n")
    file.write("opt_photprof2=" + opt_photprof2 + "\n")
    file.write("opt_obsdir_read=t\n")
    file.write("opt_surface=f\n")
    file.write("opt_int2d=f\n")
    file.write("opt_incl_gdark1=f\n")
    file.write("opt_incl_sdist1=f\n")
    file.write("opt_incl_gdark2=f\n")
    file.write("opt_incl_sdist2=f\n")
    file.write("opt_pgrid01 = 'log'\n")
    file.write("opt_rgrid01 = 'log'\n")
    file.write("opt_pgrid02 = 'log'\n")
    file.write("opt_rgrid02 = 'log'\n")
    file.write("nalpha=1\n")
    file.write("ngamma=17 /\n")
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
                      vmicro1=None,
                      vmicro2=None,
                      mdisc1=None,
                      mdisc2=None,
                      tdisc1=None,
                      tdisc2=None,
                      slope1=None,
                      slope2=None,
                      file_out=None):
    
    dir_extension =   'mdisc1v{mdisc1:.2e}'.format(mdisc1=mdisc1) + '_' \
                    + 'mdisc2v{mdisc2:.2e}'.format(mdisc2=mdisc2) + '_' \
                    + 'tdisc1v{tdisc1:.0f}'.format(tdisc1=tdisc1) + '_' \
                    + 'tdisc2v{tdisc2:.0f}'.format(tdisc2=tdisc2) + '_' \
                    + 'vmicro1v{vmicro1:.1f}'.format(vmicro1=vmicro1) + '_' \
                    + 'vmicro2v{vmicro2:.1f}'.format(vmicro2=vmicro2) + '_' \
                    + 'slope1v{slope1:.1f}'.format(slope1=slope1) + '_' \
                    + 'slope2v{slope2:.1f}'.format(slope2=slope2)

    dir_out = dir_base+'_'+dir_extension
    if file_out==None:
        print(dir_out)
    else:
        file_out.write(dir_out)
        file_out.write('\n \n')


    command01 = 'mkdir ' + dir_out
    command02 = 'mv ./outputFILES/FLUXEM*.dat ' + dir_out
    command03 = 'mv ./outputFILES/photprof*.dat ' + dir_out
    command04 = 'mv ./outputFILES/spec_points*.dat ' + dir_out
    command05 = 'mv ./outputFILES/spec_triangles*.dat ' + dir_out
    command06 = 'mv ./outputFILES/modspec_model00.h5 ' + dir_out
    command07 = 'mv ./indat_modelspec_lb1.nml ' + dir_out
    command08 = 'mv ./indat_spec_lb1.nml ' + dir_out
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
def calc_all_p05s10(value_arr, value_name=None, file_out=None):

    dir_base='/lhome/levin/Postdoc/papers/paperIII/models/lb1/parameter_study/model_p05s10'
#default values
    rstar1 = 3.7
    rmax1 = 20.
    teff1 = 18.e3
    trad1 = 18.e3
    p_object01 = [6.638, 0., 0.]
    v_object01 = [0., 15.737, 0.]

    rstar2 = 5.4
    rmax2 = 12.
    teff2 = 12.7e3
    trad2 = 12.7e3
    p_object02 = [-35.745, 0., 0.]
    v_object02 = [0., -84.736, 0.]


    opt_temperature1 = 4
    opt_temperature2 = 4
    opt_vrot1 = 1
    opt_vrot2 = 0
    vrot1 = 0.
    logg1 = 4.15
    yhe1 = 0.1
    rdisc1_max = 13.5
    vrot2 = 11.1
    logg2 = 3.1
    yhe2 = 0.2
    rdisc2_max = 14.5

    opt_photprof1=3
    opt_photprof2=3


    tdisc1 = 10.e3
    mdisc1 = 6.3e-13
    slope1 = 1.5
    vmicro1 = 100.

    tdisc2 = 10.e3
    mdisc2 = 9.2e-12
    slope2 = 3.5
    vmicro2 = 30.

    nd = np.size(value_arr)
    k = 1
    
    for value in value_arr:

        if file_out == None:
            print('--------------------running model {k:} / {nd:}---------------------'.format(k=k,nd=nd))
            print('')
        else:
            file_out.write('--------------------running model {k:} / {nd:}---------------------'.format(k=k,nd=nd))
            file_out.write('\n \n')
            file_out.flush()
        k=k+1

        #reset to default values
        vmicro1_dum = np.copy(vmicro1)
        vmicro2_dum = np.copy(vmicro2)
        mdisc1_dum = np.copy(mdisc1)
        mdisc2_dum = np.copy(mdisc2)
        tdisc1_dum = np.copy(tdisc1)
        tdisc2_dum = np.copy(tdisc2)
        slope1_dum = np.copy(slope1)
        slope2_dum = np.copy(slope2)

        #only overwrite current value variation
        if value_name == 'vmicro1':
            vmicro1_dum = value
        elif value_name == 'vmicro2':
            vmicro2_dum = value
        elif value_name == 'mdisc1':
            mdisc1_dum = value
        elif value_name == 'mdisc2':
            mdisc2_dum = value
        elif value_name == 'tdisc1':
            tdisc1_dum = value
        elif value_name == 'tdisc2':
            tdisc2_dum = value
        elif value_name == 'slope1':
            slope1_dum = value
        elif value_name == 'slope2':
            slope2_dum = value
        else:
            if file_out == None:
                print('error in calc_all_p05s10: value_name not properly specified')
            else:
                file_out.write('error in calc_all_p05s10: value_name not properly specified')
            exit()

        indat_modelspec(rstar1=rstar1,
                        rstar2=rstar2,
                        rmax1=rmax1,
                        rmax2=rmax2,
                        teff1=teff1,
                        teff2=teff2,
                        trad1=trad1,
                        trad2=trad2,
                        p_object01=p_object01,
                        p_object02=p_object02,
                        v_object01=v_object01,
                        v_object02=v_object02, 
                        opt_temperature1=opt_temperature1,
                        opt_temperature2=opt_temperature2,
                        opt_vrot1=opt_vrot1,
                        opt_vrot2=opt_vrot2,
                        vmicro1=vmicro1_dum,
                        vmicro2=vmicro2_dum,
                        vrot1=vrot1,
                        vrot2=vrot2,
                        logg1=logg1,
                        logg2=logg2,
                        tdisc1=tdisc1_dum,
                        tdisc2=tdisc2_dum,
                        mdisc1=mdisc1_dum,
                        mdisc2=mdisc2_dum,
                        rdisc1_max=rdisc1_max,
                        rdisc2_max=rdisc2_max,
                        slope1=slope1_dum,
                        slope2=slope2_dum)

        indat_spec(opt_photprof1=opt_photprof1,
                   opt_photprof2=opt_photprof2)

        if file_out == None:
            print('-----------------running modelspec_vbin.eo-------------------')
            print('')
        else:
            file_out.write('-----------------running modelspec_vbin.eo-------------------')
            file_out.write('\n \n')
            file_out.flush()
            
        os.system("./modelspec_vbin.eo < in_modelspec_vbin > output_modelspec.log")
        #
        #
        #
        if file_out == None:
            print('--------------------running spec_vbin.eo---------------------')
            print('')
        else:
            file_out.write('--------------------running spec_vbin.eo---------------------')
            file_out.write('\n \n')
            file_out.flush()
        #
        os.system("./spec_vbin.eo < in_spec_vbin > output_spec.log")
        #
        #
        #
        if file_out == None:
            print('------------------copy everything to directory---------------')
            print('')
        else:
            file_out.write('------------------copy everything to directory---------------')
            file_out.write('\n \n')
            file_out.flush()
            
        copy_to_directory(dir_base,
                          mdisc1=mdisc1_dum,
                          mdisc2=mdisc2_dum,
                          vmicro1=vmicro1_dum,
                          vmicro2=vmicro2_dum,
                          tdisc1=tdisc1_dum,
                          tdisc2=tdisc2_dum,
                          slope1=slope1_dum,
                          slope2=slope2_dum,
                          file_out=file_out)
#
#---------------------------------------------------------------------
#
def calc_all_p02s12(value_arr, value_name=None, file_out=None):

    dir_base='/lhome/levin/Postdoc/papers/paperIII/models/lb1/parameter_study/model_p02s12'
#default values
    rstar1 = 1.
    rmax1 = 180.
    teff1 = 18.e1
    trad1 = 18.e1
    p_object01 = [8.34, 0., 0.]
    v_object01 = [0., 19.78, 0.]

    rstar2 = 5.3
    rmax2 = 2.
    teff2 = 14.e3
    trad2 = 14.e3
    p_object02 = [-59.64, 0., 0.]
    v_object02 = [0., -141.3, 0.]


    opt_temperature1 = 0
    opt_temperature2 = 4
    opt_vrot1 = 1
    opt_vrot2 = 0
    vrot1 = 0.
    logg1 = 5.914
    yhe1 = 0.1
    rdisc1_max = 45.
    vrot2 = 21.4
    logg2 = 3.5
    yhe2 = 0.2
    rdisc2_max = 0.

    opt_photprof1=0
    opt_photprof2=1


    tdisc1 = 8.e3
    mdisc1 = 1.062e-14
    slope1 = 1.5
    vmicro1 = 100.

    tdisc2 = 10.e3
    mdisc2 = 1.78e-22
    slope2 = 2.5
    vmicro2 = 50.

    nd = np.size(value_arr)
    k = 1
    
    for value in value_arr:


        if file_out == None:
            print('--------------------running model {k:} / {nd:}---------------------'.format(k=k,nd=nd))
            print('')
        else:
            file_out.write('--------------------running model {k:} / {nd:}---------------------'.format(k=k,nd=nd))
            file_out.write('\n \n')
            file_out.flush()
            
        k=k+1

        #reset to default values
        vmicro1_dum = np.copy(vmicro1)
        vmicro2_dum = np.copy(vmicro2)
        mdisc1_dum = np.copy(mdisc1)
        mdisc2_dum = np.copy(mdisc2)
        tdisc1_dum = np.copy(tdisc1)
        tdisc2_dum = np.copy(tdisc2)
        slope1_dum = np.copy(slope1)
        slope2_dum = np.copy(slope2)

        #calculate rho0
        if value_name == 'tdisc1':
            tdisc1_dum = value['tdisc1']
            rho0_disc1 = value['rho0_disc1']

            #adapt mdisc1 input parameter to get correct rho0_disc1
            mmw1 = mean_molecular_weight(ih=1.,ihe=2.,yhe=yhe1)
            csound = vsound(tdisc1_dum,mmw1)
            sr1 = rstar1*cgs_rsu
            mstar1_cgs = sr1**2 * 10.**logg1/cgs_grav
#            rho0_disc1 = mdisc1**cgs_msu*np.sqrt(cgs_grav*mstar1_cgs)*np.log(10.)/2./np.sqrt(np.pi**3)/csound/sr1**3.5
            mdisc1_dum = rho0_disc1/(cgs_msu*np.sqrt(cgs_grav*mstar1_cgs)*np.log(10.)/2./np.sqrt(np.pi**3)/csound/sr1**3.5)
            print('adapting mdisc1 to', mdisc1_dum)

        #only overwrite current value variation
        elif value_name == 'vmicro1':
            vmicro1_dum = value
        elif value_name == 'vmicro2':
            vmicro2_dum = value
        elif value_name == 'mdisc1':
            mdisc1_dum = value
        elif value_name == 'mdisc2':
            mdisc2_dum = value
        elif value_name == 'tdisc2':
            tdisc2_dum = value
        elif value_name == 'slope1':
            slope1_dum = value
        elif value_name == 'slope2':
            slope2_dum = value
        else:
            if file_out == None:
                print('error in calc_all_p05s10: value_name not properly specified')
            else:
                file_out.write('error in calc_all_p05s10: value_name not properly specified')
            exit()

        indat_modelspec(rstar1=rstar1,
                        rstar2=rstar2,
                        rmax1=rmax1,
                        rmax2=rmax2,
                        teff1=teff1,
                        teff2=teff2,
                        trad1=trad1,
                        trad2=trad2,
                        p_object01=p_object01,
                        p_object02=p_object02,
                        v_object01=v_object01,
                        v_object02=v_object02, 
                        opt_temperature1=opt_temperature1,
                        opt_temperature2=opt_temperature2,
                        opt_vrot1=opt_vrot1,
                        opt_vrot2=opt_vrot2,
                        vmicro1=vmicro1_dum,
                        vmicro2=vmicro2_dum,
                        vrot1=vrot1,
                        vrot2=vrot2,
                        logg1=logg1,
                        logg2=logg2,
                        tdisc1=tdisc1_dum,
                        tdisc2=tdisc2_dum,
                        mdisc1=mdisc1_dum,
                        mdisc2=mdisc2_dum,
                        rdisc1_max=rdisc1_max,
                        rdisc2_max=rdisc2_max,
                        slope1=slope1_dum,
                        slope2=slope2_dum)

        indat_spec(opt_photprof1=opt_photprof1,
                   opt_photprof2=opt_photprof2)

        if file_out == None:
            print('-----------------running modelspec_vbin.eo-------------------')
            print('')
        else:
            file_out.write('-----------------running modelspec_vbin.eo-------------------')
            file_out.write('\n \n')
            file_out.flush()
            
        os.system("./modelspec_vbin.eo < in_modelspec_vbin > output_modelspec.log")
        #
        #
        #
        if file_out == None:
            print('--------------------running spec_vbin.eo---------------------')
            print('')
        else:
            file_out.write('--------------------running spec_vbin.eo---------------------')
            file_out.write('\n \n')
            file_out.flush()
        #
        os.system("./spec_vbin.eo < in_spec_vbin > output_spec.log")
        #
        #
        #
        if file_out == None:
            print('------------------copy everything to directory---------------')
            print('')
        else:
            file_out.write('------------------copy everything to directory---------------')
            file_out.write('\n \n')
            file_out.flush()
            
        copy_to_directory(dir_base,
                          mdisc1=mdisc1_dum,
                          mdisc2=mdisc2_dum,
                          vmicro1=vmicro1_dum,
                          vmicro2=vmicro2_dum,
                          tdisc1=tdisc1_dum,
                          tdisc2=tdisc2_dum,
                          slope1=slope1_dum,
                          slope2=slope2_dum)
#
#---------------------------------------------------------------------
#
def calc_all_p02s12_chi2(value1_arr, value2_arr, value_names=None, file_out=None):

    dir_base='/lhome/levin/Postdoc/papers/paperIII/models/lb1/parameter_study/model_p02s12'
#default values
    rstar1 = 1.
    rmax1 = 180.
    teff1 = 18.e1
    trad1 = 18.e1
    p_object01 = [8.34, 0., 0.]
    v_object01 = [0., 19.78, 0.]

    rstar2 = 5.3
    rmax2 = 2.
    teff2 = 14.e3
    trad2 = 14.e3
    p_object02 = [-59.64, 0., 0.]
    v_object02 = [0., -141.3, 0.]


    opt_temperature1 = 0
    opt_temperature2 = 4
    opt_vrot1 = 1
    opt_vrot2 = 0
    vrot1 = 0.
    logg1 = 5.914
    yhe1 = 0.1
    rdisc1_max = 45.
    vrot2 = 21.4
    logg2 = 3.5
    yhe2 = 0.2
    rdisc2_max = 0.

    opt_photprof1=0
    opt_photprof2=1


    tdisc1 = 8.e3
    mdisc1 = 1.062e-14
    slope1 = 1.5
    vmicro1 = 100.

    tdisc2 = 10.e3
    mdisc2 = 1.78e-22
    slope2 = 2.5
    vmicro2 = 50.

    nd = np.size(value1_arr)*np.size(value2_arr)
    k = 1

    value1_name=value_names[0]
    value2_name=value_names[1]    
    
    for value1 in value1_arr:
        for value2 in value2_arr:

            if file_out == None:
                print('--------------------running model {k:} / {nd:}---------------------'.format(k=k,nd=nd))
                print('')
            else:
                file_out.write('--------------------running model {k:} / {nd:}---------------------'.format(k=k,nd=nd))
                file_out.write('\n \n')
                file_out.flush()
                
            k=k+1

            #reset to default values
            vmicro1_dum = np.copy(vmicro1)
            vmicro2_dum = np.copy(vmicro2)
            mdisc1_dum = np.copy(mdisc1)
            mdisc2_dum = np.copy(mdisc2)
            tdisc1_dum = np.copy(tdisc1)
            tdisc2_dum = np.copy(tdisc2)
            slope1_dum = np.copy(slope1)
            slope2_dum = np.copy(slope2)

            #calculate rho0
            if value1_name == 'tdisc1':
                if value2_name == 'mdisc1':
                    #calculate rho0 for the standard model with 8 kK
                    tdisc1_dum = 8.e3                    
                    mmw1 = mean_molecular_weight(ih=1.,ihe=2.,yhe=yhe1)
                    csound = vsound(tdisc1_dum,mmw1)
                    sr1 = rstar1*cgs_rsu
                    mstar1_cgs = sr1**2 * 10.**logg1/cgs_grav
                    rho0_disc1 = value2*cgs_msu*np.sqrt(cgs_grav*mstar1_cgs)*np.log(10.)/2./np.sqrt(np.pi**3)/csound/sr1**3.5
                else:
                    rho0_disc1 = value1['rho0_disc1']

                tdisc1_dum = value1['tdisc1']
                #adapt mdisc1 input parameter to get correct rho0_disc1
                mmw1 = mean_molecular_weight(ih=1.,ihe=2.,yhe=yhe1)
                csound = vsound(tdisc1_dum,mmw1)
                sr1 = rstar1*cgs_rsu
                mstar1_cgs = sr1**2 * 10.**logg1/cgs_grav
                #rho0_disc1 = mdisc1**cgs_msu*np.sqrt(cgs_grav*mstar1_cgs)*np.log(10.)/2./np.sqrt(np.pi**3)/csound/sr1**3.5
                mdisc1_dum = rho0_disc1/(cgs_msu*np.sqrt(cgs_grav*mstar1_cgs)*np.log(10.)/2./np.sqrt(np.pi**3)/csound/sr1**3.5)
            elif value1_name == 'vmicro1':
                vmicro1_dum = value1
            elif value1_name == 'vmicro2':
                vmicro2_dum = value1
            elif value1_name == 'mdisc1':
                mdisc1_dum = value1
            elif value1_name == 'mdisc2':
                mdisc2_dum = value1
            elif value1_name == 'tdisc2':
                tdisc2_dum = value1
            elif value1_name == 'slope1':
                slope1_dum = value1
            elif value1_name == 'slope2':
                slope2_dum = value1
            else:
                if file_out == None:
                    print('error in calc_all_p02s12: value_name not properly specified')
                else:
                    file_out.write('error in calc_all_p02s12: value_name not properly specified')
                exit()


            if value2_name == 'tdisc1':
                tdisc1_dum = value2['tdisc1']
                rho0_disc1 = value2['rho0_disc1']
                #adapt mdisc1 input parameter to get correct rho0_disc1
                mmw1 = mean_molecular_weight(ih=1.,ihe=2.,yhe=yhe1)
                csound = vsound(tdisc1_dum,mmw1)
                sr1 = rstar1*cgs_rsu
                mstar1_cgs = sr1**2 * 10.**logg1/cgs_grav
                #rho0_disc1 = mdisc1**cgs_msu*np.sqrt(cgs_grav*mstar1_cgs)*np.log(10.)/2./np.sqrt(np.pi**3)/csound/sr1**3.5
                mdisc1_dum = rho0_disc1/(cgs_msu*np.sqrt(cgs_grav*mstar1_cgs)*np.log(10.)/2./np.sqrt(np.pi**3)/csound/sr1**3.5)
                #print('adapting mdisc1 to', mdisc1_dum)
                #only overwrite current value variation
            elif value2_name == 'vmicro1':
                vmicro1_dum = value2
            elif value2_name == 'vmicro2':
                vmicro2_dum = value2
            elif value2_name == 'mdisc1':
                mdisc1_dum = value2
            elif value2_name == 'mdisc2':
                mdisc2_dum = value2
            elif value2_name == 'tdisc2':
                tdisc2_dum = value2
            elif value2_name == 'slope1':
                slope1_dum = value2
            elif value2_name == 'slope2':
                slope2_dum = value2
            else:
                if file_out == None:
                    print('error in calc_all_p02s12: value_name not properly specified')
                else:
                    file_out.write('error in calc_all_p02s12: value_name not properly specified')
                exit()

            indat_modelspec(rstar1=rstar1,
                            rstar2=rstar2,
                            rmax1=rmax1,
                            rmax2=rmax2,
                            teff1=teff1,
                            teff2=teff2,
                            trad1=trad1,
                            trad2=trad2,
                            p_object01=p_object01,
                            p_object02=p_object02,
                            v_object01=v_object01,
                            v_object02=v_object02, 
                            opt_temperature1=opt_temperature1,
                            opt_temperature2=opt_temperature2,
                            opt_vrot1=opt_vrot1,
                            opt_vrot2=opt_vrot2,
                            vmicro1=vmicro1_dum,
                            vmicro2=vmicro2_dum,
                            vrot1=vrot1,
                            vrot2=vrot2,
                            logg1=logg1,
                            logg2=logg2,
                            tdisc1=tdisc1_dum,
                            tdisc2=tdisc2_dum,
                            mdisc1=mdisc1_dum,
                            mdisc2=mdisc2_dum,
                            rdisc1_max=rdisc1_max,
                            rdisc2_max=rdisc2_max,
                            slope1=slope1_dum,
                            slope2=slope2_dum)

            indat_spec(opt_photprof1=opt_photprof1,
                       opt_photprof2=opt_photprof2)

            if file_out == None:
                print('-----------------running modelspec_vbin.eo-------------------')
                print('')
            else:
                file_out.write('-----------------running modelspec_vbin.eo-------------------')
                file_out.write('\n \n')
                file_out.flush()
            #
            os.system("./modelspec_vbin.eo < in_modelspec_vbin > output_modelspec.log")
            #
            #
            #
            if file_out == None:
                print('--------------------running spec_vbin.eo---------------------')
                print('')
            else:
                file_out.write('--------------------running spec_vbin.eo---------------------')
                file_out.write('\n \n')
                file_out.flush()
            #
            os.system("./spec_vbin.eo < in_spec_vbin > output_spec.log")
            #
            #
            #
            if file_out == None:
                print('------------------copy everything to directory---------------')
                print('')
            else:
                file_out.write('------------------copy everything to directory---------------')
                file_out.write('\n \n')
                file_out.flush()

            copy_to_directory(dir_base,
                              mdisc1=mdisc1_dum,
                              mdisc2=mdisc2_dum,
                              vmicro1=vmicro1_dum,
                              vmicro2=vmicro2_dum,
                              tdisc1=tdisc1_dum,
                              tdisc2=tdisc2_dum,
                              slope1=slope1_dum,
                              slope2=slope2_dum,
                              file_out=file_out)
#
#---------------------------------------------------------------------
#----------------------model p05s10-----------------------------------
#---------------------------------------------------------------------
#
mdisc1_arr = [6.3e-15, 6.3e-14, 1.12e-13, 1.99e-13, 3.54e-13, 6.3e-13, 1.12e-12, 1.99e-12, 3.54e-12, 6.3e-12, 6.3e-11]
mdisc2_arr = [9.2e-14, 9.2e-13, 1.64e-12, 2.91e-12, 5.17e-12, 9.2e-12, 1.64e-11, 2.91e-11, 5.17e-11, 9.2e-11, 1.9e-10]
slope1_arr = [1.5, 2., 2.5, 3., 3.5, 4.]
slope2_arr = [1.5, 2., 2.5, 3., 3.5, 4.]
vmicro1_arr = [15.,30.,45.,60.,75.,90.,105.]
vmicro2_arr = [15.,30.,45.,60.,75.,90.,105.]

#calc_all_p05s10(mdisc1_arr, value_name='mdisc1')
#calc_all_p05s10(mdisc2_arr, value_name='mdisc2')
#calc_all_p05s10(slope1_arr, value_name='slope1')
#calc_all_p05s10(slope2_arr, value_name='slope2')
#calc_all_p05s10(vmicro1_arr, value_name='vmicro1')
#calc_all_p05s10(vmicro2_arr, value_name='vmicro2')
#
#---------------------------------------------------------------------
#----------------------model p02s12-----------------------------------
#---------------------------------------------------------------------
#
mdisc1_arr = [1.062e-16, 1.062e-15, 1.8885327e-15, 3.3583389e-15, 5.9720649e-15, 1.062e-14, 1.8885327e-14, 3.3583389e-14, 5.9720649e-14, 1.0620000e-13, 1.062e-12]
slope1_arr = [1.5, 2., 2.5, 3., 3.5, 4.]
vmicro1_arr = [15.,30.,45.,60.,75.,90.,105.]
#tdisc1_arr = [6.e3,7.e3,8.e3,9.e3,10.e3,11.e3,12.e3,13.e3,14.e3,15.e3]
tdisc1_arr = [{'tdisc1': 6.e3, 'rho0_disc1': 2.97e-12},
              {'tdisc1': 7.e3, 'rho0_disc1': 2.97e-12},
              {'tdisc1': 8.e3, 'rho0_disc1': 2.97e-12},              
              {'tdisc1': 9.e3, 'rho0_disc1': 2.97e-12},
              {'tdisc1': 10.e3, 'rho0_disc1': 2.97e-12},
              {'tdisc1': 11.e3, 'rho0_disc1': 2.97e-12},
              {'tdisc1': 12.e3, 'rho0_disc1': 2.97e-12},
              {'tdisc1': 13.e3, 'rho0_disc1': 2.97e-12},
              {'tdisc1': 14.e3, 'rho0_disc1': 2.97e-12},
              {'tdisc1': 15.e3, 'rho0_disc1': 2.97e-12}]
#calc_all_p02s12(mdisc1_arr, value_name='mdisc1')
#calc_all_p02s12(slope1_arr, value_name='slope1')
#calc_all_p02s12(vmicro1_arr, value_name='vmicro1')
#calc_all_p02s12(tdisc1_arr, value_name='tdisc1')

file_out = open("nohup.log", "w")
#calc_all_p02s12_chi2(tdisc1_arr, mdisc1_arr, value_names=['tdisc1','mdisc1'], file_out=file_out)
#calc_all_p02s12_chi2(mdisc1_arr, slope1_arr, value_names=['mdisc1','slope1'], file_out=file_out)
#calc_all_p02s12_chi2(mdisc1_arr, vmicro1_arr, value_names=['mdisc1','vmicro1'], file_out=file_out)
calc_all_p02s12_chi2(tdisc1_arr, slope1_arr, value_names=['tdisc1','slope1'], file_out=file_out)
calc_all_p02s12_chi2(tdisc1_arr, vmicro1_arr, value_names=['tdisc1','vmicro1'], file_out=file_out)
calc_all_p02s12_chi2(slope1_arr, vmicro1_arr, value_names=['slope1','vmicro1'], file_out=file_out)

file_out.close()

