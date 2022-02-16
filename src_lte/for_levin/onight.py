import sys, os
import numpy as np
from const_cgs import *
#
#---------------------------------------------------------------------
#
def indat_file(dir_out='output',
               lgTmin  = 4.,
               lgTmax  = 5.5,
               N_lgT   = 40,
               lgDmin  = -15.,
               lgDmax  = -7.,
               N_lgD   = 40,
               X_mass = 0.,
               Z_mass = 0.02,
               SpecZ   = 6,
               SpecI   = 3,
               SpecL   = 9,
               ver     = '.FALSE.'):
    
    fname='indat.par'
    
    file = open(fname, "w")
    file.write("&init_param\n")
    file.write("DIR = '" + dir_out + "' \n")      
    file.write("lgTmin = {x:} \n".format(x=lgTmin))
    file.write("lgTmax = {x:} \n".format(x=lgTmax))
    file.write("N_lgT = {x:} \n".format(x=N_lgT))    
    file.write("lgDmin = {x:} \n".format(x=lgDmin))
    file.write("lgDmax = {x:} \n".format(x=lgDmax))    
    file.write("N_lgD = {x:} \n".format(x=N_lgD))  
    file.write("X_mass = {x:} \n".format(x=X_mass))  
    file.write("Z_mass = {x:} \n".format(x=Z_mass))
    file.write("SpecZ = {x:} \n".format(x=SpecZ))  
    file.write("SpecI = {x:} \n".format(x=SpecI))  
    file.write("SpecL = {x:} \n".format(x=SpecL))
    file.write("ver = .false. \n")
    file.write("/\n")
   
    file.close()
#
#---------------------------------------------------------------------
#
def copy_to_directory(dir_base,
                      y_mass=0.98, SpecZ=None, SpecI=None, SpecL=None):
    dir_ext = 'error'
    if y_mass == 0.0999: dir_ext='Y00999'
    if y_mass == 0.1: dir_ext='Y01000'
    if y_mass == 0.28: dir_ext='Y02800'
    if y_mass == 0.59: dir_ext='Y05900'
    if y_mass == 0.61: dir_ext='Y06100'
    if y_mass == 0.98: dir_ext='Y09800'
    if y_mass == 0.999: dir_ext='Y09990'
    if y_mass == 1.: dir_ext='Y10000'
    
    dir_out = dir_base+'/'+dir_ext

    log_out = dir_out + '/output_{SpecZ:02d}_{SpecI:02d}_{SpecL:02d}.log'.format(SpecZ=SpecZ,SpecL=SpecL,SpecI=SpecI)

    command01 = 'mkdir ' + dir_out
    command02 = 'mv indat.par ' + dir_out
    command03 = 'mv ./output/*.txt ' + dir_out
    command04 = 'mv output.log ' + log_out

#    print(command01)
#    print(command02)
#    print(command03)
#    print(command04)    
    os.system(command01)
    os.system(command02)
    os.system(command03)
#    os.system(command04)
#
#---------------------------------------------------------------------
#
def calc_all(y_mass=1.):

    dir_base='../../lte_tables'

    z_mass=0.02
    #same as for OPAL tables
    if y_mass == 0.0999: z_mass=0.0001
    if y_mass == 0.1: z_mass=0.
    if y_mass == 0.28: z_mass=0.02
    if y_mass == 0.59: z_mass=0.06
    if y_mass == 0.61: z_mass=0.04
    if y_mass == 0.98: z_mass=0.02
    if y_mass == 0.999: z_mass=0.001
    if y_mass == 1.: z_mass=0.
    X_mass = 1.-y_mass-z_mass

    #up until Z=26 (iron)
    for SpecZ in np.arange(1,27):
        print('calculating element species', SpecZ)
        for SpecI in np.arange(1,SpecZ+1):
           print('calculating ionization stage', SpecI)
           for SpecL in np.arange(1,17):
               print('calculating transitions to level', SpecL)
               indat_file(dir_out='output',
                          lgTmin  = 4.,
                          lgTmax  = 6.,
                          N_lgT   = 40,
                          lgDmin  = -16.,
                          lgDmax  = -6.,
                          N_lgD   = 40,
                          X_mass = X_mass,
                          Z_mass = z_mass,
                          SpecZ   = SpecZ,
                          SpecI   = SpecI,
                          SpecL   = SpecL)

               os.system("./run -i indat.par > output.log")
               copy_to_directory(dir_base, y_mass=y_mass, SpecZ=SpecZ, SpecI=SpecI, SpecL=SpecL)
#

#calc_all(y_mass=0.98)
#calc_all(y_mass=1.)
#calc_all(y_mass=0.999)
#calc_all(y_mass=0.61)

#calc_all(y_mass=0.28)
#calc_all(y_mass=0.59)
#calc_all(y_mass=0.1)
calc_all(y_mass=0.0999)

