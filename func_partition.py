import numpy     as np
from   func_flux import *

def f_xH2O(x, phi, Mmantle, Mtot, rp, g, mu_avg):
    fatm = 4*np.pi*rp**2/g * (x/(phi*AH2O))**(1./0.7) * muv/mu_avg
    fmat = x*Mmantle
    
    return fatm + fmat - Mtot

def calc_xH2O(phi, Mres, Mmantle, xH2O_init, rp, g, mu_avg):
    xv0 = xH2O_init*1e-4
    xv1 = xH2O_init*10
    f0  = f_xH2O(xv0, phi, Mmantle, Mres, rp, g, mu_avg)
    f1  = f_xH2O(xv1, phi, Mmantle, Mres, rp, g, mu_avg)
    
    eps = 1e-12
    while (np.log(xv1/xv0)>eps):
        # new point
        xvA = (xv0+xv1)*0.5
        fA  = f_xH2O(xvA, phi, Mmantle, Mres, rp, g, mu_avg)
        
        if (f0*fA < 0):
            xv1 = xvA
            f1  = fA
        else:
            xv0 = xvA
            f0  = fA

    #print(xv0, xv1, xvA, f0, f1, fA)
    #print(fA)
    return (xv0+xv1)*0.5
