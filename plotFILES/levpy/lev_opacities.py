from const_cgs import *

def opac_thomson(yhe=None, ihe=None, rho=None, kcont=None):
    #!calculate thomson opacity, assuming fully ionized hydrogen
    #!
    #!  INPUT: yhe:     helium abundance by number
    #!         hei:     number of free electrons per helium atom (=2 for full ionization)
    #!         rho:     density in cgs
    #!         kcont:   multiplication factor
    #!  OUTPUT: thomson opacity in 1/cm

    c1=(1.+4.*yhe)*cgs_mp
    c2=(1.+ihe*yhe)/c1
    ne=c2*rho
    opac_thomson=cgs_sigmae*ne*kcont

    return opac_thomson
