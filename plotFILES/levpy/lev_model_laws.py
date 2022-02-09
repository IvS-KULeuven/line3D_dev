import numpy as np
from const_cgs import *
#
#-----------------------------------------------------------------------
#----------------mean molecular weight----------------------------------
#-----------------------------------------------------------------------
#
def mean_molecular_weight(ih=None,ihe=None,yhe=None):
#
#     calculates mean molecular weight (if contribution of metals besides helium
#              is negligible for both the electron density and the mass
#
#input: ih:   number free electrons per hydrogen atom
#       ihe:  number free electrons per helium atom
#       yhe:  helium abundance by number (N_He/N_H)
#
    mmw = (1. + 4.*yhe)/(1.+ih + (1.+ihe)*yhe)

    return mmw
#
#-----------------------------------------------------------------------
#----------------isothermal sound speed---------------------------------
#-----------------------------------------------------------------------
#
def vsound(temp, mmw):
#
#     calculates isothermal sound speed (sqrt(p/rho)) for
#   given temperature (temp) and mean molecular weight (mmw)
#
#input: temp:   temperature in Kelvin
#       mmw:    mean molecular weight
#  
    vsound = np.sqrt(cgs_kb*temp/mmw/cgs_mp)
    
    return vsound
