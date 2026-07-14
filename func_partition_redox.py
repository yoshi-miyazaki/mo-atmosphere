import numpy as np


from func_real_redox import molar_mass

# import molar mass
muH2, muv, muCO, muCO2 = molar_mass[0], molar_mass[1], molar_mass[2], molar_mass[3]

# constants
AH2O = 6.8e-8

def calc_atmos_composition(xH2_init, xH2O_init, xCO_init, xCO2_init, MN2, M_mantle, rp, g):

    # solve for partitioning
    # average molar mass is needed, but make a guess first.
    mu_avg = muCO2

    for i in range(5):
        PCO2      = M_mantle * xCO2_init / (4*np.pi*rp*rp/g * muCO2/mu_avg + M_mantle*4.4e-12)
        #print(PCO2/1e5, " bar ", PCO2*4.4e-12*100, " wt%")

        PCO       = M_mantle * xCO_init  / (4*np.pi*rp*rp/g *  muCO/mu_avg + M_mantle*0.55e-12)
        PH2       = M_mantle * xH2_init  / (4*np.pi*rp*rp/g *  muH2/mu_avg + M_mantle*5.e-12)
        
        # for Figure 1
        # PCO2 = 200e5

        xH2O      = calc_xH2O(1., xH2O_init*M_mantle, M_mantle, rp, g, mu_avg)
        PH2O      = (xH2O/AH2O)**(1./0.7) 
        #print(PH2O/1e5, " bar")

        MCO2 = PCO2 * 4*np.pi*rp*rp/g * muCO2/mu_avg
        MCO  = PCO  * 4*np.pi*rp*rp/g * muCO /mu_avg
        MH2O = PH2O * 4*np.pi*rp*rp/g * muv  /mu_avg
        MH2  = PH2  * 4*np.pi*rp*rp/g * muH2 /mu_avg

        r_atmos = np.array([MH2, MH2O, MCO, MCO2, MN2])
        x_atmos = r_atmos / molar_mass
        x_atmos = x_atmos / x_atmos.sum()
        print(x_atmos, "\t", mu_avg - np.sum(x_atmos * molar_mass))
        
        mu_avg = np.sum(x_atmos * molar_mass)
        PN2 = MN2 * g / (4*np.pi*rp*rp) * mu_avg/molar_mass[2]
        

    
    # create composition array
    P_atmos = np.array([PH2, PH2O, PCO, PCO2, PN2])
    
    print(PH2/1e5,  " bar ", PH2*5e-12*100, " wt%")
    print(PH2O/1e5, " bar ")
    print(PCO/1e5,  " bar ", PCO*0.55e-12*100, " wt%")
    print(PCO2/1e5, " bar ", PCO2*4.4e-12*100, " wt%")

    return P_atmos


def f_xH2O(xH2O, phi, MH2O_tot, Mmantle, rp, g, mu_avg):
    fatm = 4*np.pi*rp**2/g * (xH2O/(phi*AH2O))**(1./0.7) * muv/mu_avg
    fmat = xH2O*Mmantle
    
    return fatm + fmat - MH2O_tot

def calc_xH2O(phi, MH2O_tot, Mmantle, rp, g, mu_avg):
    xH2O_max = MH2O_tot/Mmantle
    xv0      = xH2O_max*1e-4
    xv1      = xH2O_max
    
    f0  = f_xH2O(xv0, phi, MH2O_tot, Mmantle, rp, g, mu_avg)
    f1  = f_xH2O(xv1, phi, MH2O_tot, Mmantle, rp, g, mu_avg)
    
    eps = 1e-12
    #print(f"f0, f1 = {f0:.3e}\t{f1:.3e} at {xv1:.5e}")
    while ((xv1/xv0)>1.0000001):
        # new point
        xvA = (xv0+xv1)*0.5
        fA  = f_xH2O(xvA, phi, MH2O_tot, Mmantle, rp, g, mu_avg)
        
        if (f0*fA < 0):
            xv1 = xvA
            f1  = fA
        else:
            xv0 = xvA
            f0  = fA

    #print(f"f0, f1 = {f0:.3e}\t{f1:.3e} at {xv1:.5e}")
    
    return (xv0+xv1)*0.5


