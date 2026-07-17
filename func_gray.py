import numpy as np
import sys
from numba import njit

from const import *


# each array contains information in the order of H2O, CO2, N2
kappa         = np.array([1.e-2, 1.3e-4,  0])         # m2/kg
molar_mass    = np.array([18,    44,     28]) * 1e-3  # kg/mol
heat_capacity = np.array([4.,    3.5,   3.5]) * R     # J/(K mol)

P_ref = 1e5

# kpaCO2 = kappa[1]

''' water saturation '''
pA,lH2O = 1.4e11, 4.3655e4

@njit
def P_watersaturation(T):
    return pA * np.exp(- lH2O / (R*T))

@njit
def dP_watersaturation(T):
    return pA * np.exp(- lH2O / (R*T)) * (lH2O/(R*T*T))


''' convective T gradient '''
@njit
def moist_dTdP(dP,T,P,q,Cp):
    # Use second order Runge-Kutta to calculate moist adiabat
    # 
    # q:  water mixing ratio (= xH2O / (1.-xH2O))
    # Cp: avearge heat capacity

    def calc_moist(T,P,q,Cp):
        numer = (R*T)/(P*Cp) * (1.+ (q*lH2O)/(R*T))
        denom = 1. + (q*lH2O*lH2O)/(R*T*T*Cp)
    
        return numer/denom

    # estimate T at P = P+dP*0.5
    Gamma1 = calc_moist(T,P,q,Cp)
    T1     = T + Gamma1 * dP*.5

    # estimate T at P = P
    Gamma2 = calc_moist(T1,P+dP*.5,q,Cp)
    # print(Gamma1, Gamma2, "\t", T, T1)
    
    return Gamma2
    

''' determine tau, T at the bottom of a uniform stratosphere '''
def tau_saturation(rH, x_atmos, Fnet, g):
    # return (tau, T) at the bottom of a uniform stratosphere
    #
    # In a compositionally uniform stratosphere, two conditions apply:
    #   1. Optical depth:         tau = kappa * P / g
    #   2. Radiative equilibrium: σ * T^4 = (3/4) * Fnet * (tau + 2/3)
    #
    # At the base of the stratosphere, water vapor reaches saturation and the
    # atmospheric composition begins to change. This conditions determines
    # pressure, and thus tau, T.
    
    
    T0, T1  = 10, tau_upper_search_limit(rH, x_atmos, Fnet, g)
    f0, tau = dtau(T0, rH, x_atmos, Fnet, g)
    f1, tau = dtau(T1, rH, x_atmos, Fnet, g)

    if (f0 > 0 and f1 > 0):
        # print(f"flux over the limit - dtau {f0:.3e}, {f1:.3e} > 0 \t at T = {T1:.3f} for {Fnet:.0f}, {x_atmos[0]:.3e}")
        return -1, -1
        
    elif (f0*f1 > 0):
        print(f"bisection error - dtau {f0:.3e}, {f1:.3e} \t {T1:.3f}", x_atmos)
        sys.exit()

    # print("search between", T0, T1, " f =", f0, f1)
    while (f0*f1 < 0 and np.abs(T1-T0) > 1e-3):
        TA = (T0+T1)*.5
        fA, tau = dtau(TA, rH, x_atmos, Fnet, g)

        if (f0*fA < 0):
            T1, f1 = TA, fA
        else:
            T0, f0 = TA, fA            

    # print("tau/T at saturation depth", tau, TA)

    return tau, TA

@njit
def dtau(T, rH, x_atmos, Fnet, g):
    # tau_saturation findes T where dtau = 0
    # dtau = (tau calculated from Eq. 1) - (tau calculated from Eq.2)
    
    # set opacity
    kpa = np.sum(kappa * x_atmos * molar_mass) / np.sum(x_atmos * molar_mass)
    
    # (humidity) * (water saturation pressure) = (total pressure) * (water ratio)
    P = rH * P_watersaturation(T) / x_atmos[0]
        
    # compare tau derived from two equations
    # (optical depth) = (opacity) * (surface density)
    tauo = kpa * P / g
    
    # radiative equilibrium
    taur = 4.*sB*(T*T*T*T) / (3*Fnet) - 2./3 

    return tauo - taur, tauo

def tau_upper_search_limit(rH, x_atmos, Fnet, g):
    # find the temperature where dtau(T) reaches its minimum.
    #
    # Since dtau(T) has either 0 or 2 roots, choosing a valid search
    # range is important. When 2 roots exist, the larger root corresponds
    # to an unstable solution, so we take the smaller one. The minimum of
    # dtau(T) always lies between the two roots, making it a useful target.
    #
    # The derivative is:
    #   d(dtau)/dT = (kappa_H2O / g) * rH * dPsat/dT  –  (16/3) * σ * T^3 / Fnet
    
    def ddtau(T, kpa, rH, x_atmos, Fnet, g):
        dPsatdT = dP_watersaturation(T)
        return kpa*rH*dPsatdT/(x_atmos[0]*g) - 16.*sB*(T*T*T)/(3.*Fnet)

    # set opacity
    kpa = np.sum(kappa * x_atmos * molar_mass)/np.sum(x_atmos * molar_mass)

    #for T in range(10, 300, 10):
    #print(T, "\t", dtau(T, rH, x_atmos, Fnet, g), ddtau(T, kpa, rH, x_atmos, Fnet, g))
    
    T0, T1 = 10, 400
    f0, f1 = ddtau(T0, kpa, rH, x_atmos, Fnet, g), ddtau(T1, kpa, rH, x_atmos, Fnet, g)

    if (f0*f1 > 0):
        print("bisection error - ddtau", f0, f1)
        sys.exit()
                    
    while (f0*f1 < 0 and np.abs(T1-T0) > 1e-3):
        TA = (T0+T1)*.5
        fA = ddtau(TA, kpa, rH, x_atmos, Fnet, g)
        
        if (f0*fA < 0):
            T1, f1 = TA, fA
        else:
            T0, f0 = TA, fA        

    # print("dtau minimum at", TA, f0, f1)
    return TA

''' '''
def determine_min_xH2O(Mnc, rH, x_atmos_nc, Fnet, rp, g):
    # Compute the minimum allowed x_H2O_top used by `calc_Tprofile`
    #
    # If x_H2O_top is too small, the column mass above the saturation pressure (P_sat)
    # becomes excessively large because # P_total ≈ P_sat / x_H2O_top.
    # therefore x_H2O_top has a lower bound — it cannot be arbitrarily smal.
    
    def calc_Mupnc(xH2O_top):
        # set composition w/ water
        xH2O_array = np.zeros(len(x_atmos_nc))
        xH2O_array[0] = xH2O_top
        x_atmos = xH2O_array + (1-xH2O_top) * x_atmos_nc
        
        tau_sat, T_sat = tau_saturation(rH, x_atmos, Fnet, g)
        if (tau_sat < 0):
            return Mnc * 100
        P_sat = rH * P_watersaturation(T_sat) / x_atmos[0]
        
        mu_avg = np.sum(molar_mass * x_atmos)
        M_atmos = P_sat * (4*np.pi*rp*rp) / g * (molar_mass * x_atmos) / mu_avg

        return M_atmos[1:].sum()
    
    x0, x1 = 1e-20, 0.999
    f0, f1 = calc_Mupnc(x0) - Mnc, calc_Mupnc(x1) - Mnc
    # print(f"dM = {f0:.3e} \t {f1:.3e}")
    
    if (f0 < 0):
        return x0

    count = 0
    while (f0 > 0 and x1/x0 > 1.00000001):
        xA = np.sqrt(x0*x1)
        fA = calc_Mupnc(xA) - Mnc

        if (fA < 0):
            x1, f1 = xA, fA
        else:
            x0, f0 = xA, fA
        # print(count, "\t", x1/x0, fA/Mnc, "\t", f0, f1)

        count += 1
        if (count > 40):
            print(f"error: determine_min_xH2O. Search for miminum xH2O_top is not converging.")
            sys.exit()

    # check whether saturated tropopause could form
    #
    xH2O_array = np.zeros(len(x_atmos_nc))
    xH2O_array[0] = x1
    x_atmos = xH2O_array + (1-x1) * x_atmos_nc
    T_dtaumin = tau_upper_search_limit(rH, x_atmos, Fnet, g)
    dt, tau = dtau(T_dtaumin, rH, x_atmos, Fnet, g)

    # FLAG
    # print(f"mass check: {x0:.3e} -> {f1 + Mnc}. Set atmos mass: {f1+Mnc:.3e}")
    # print(f"dtau check: {dt:.3e}")
    
    return x1

        

''' guess T profile '''

''' heat flux '''
@njit
def compute_fluxes(T, tau):
    # compute upward fluxes
    # using 2 different formulations
    
    D   = 1.5          # Eddington approximation
    N   = len(T)
    Fup = np.zeros(N)

    # 
    Fnet  = 2*sB*T[0]**4
    Fsurf = Fnet*(1.+0.75*tau[-1]) # sB * T[-1]**4 #
    if (Fsurf > sB*T[-1]**4 or 1): # FLAG
        Fsurf = sB*T[-1]**4 
    
    for i in range(N):
        Fup[i]  = sB * T[i]**4  + (Fsurf - sB * T[-1]**4)*np.exp(D*(-tau[-1]+tau[i]))
        # Fup[i] = sB *T[-1]**4 * np.exp(-tau[-1] + tau[i])

        for j in range(i+1, N):
            dT   = T[j]   - T[j-1]
            dtau = tau[j] - tau[j-1]
            
            Fup[i]  += 4*sB* (T[j]**3      * np.exp(D*(-tau[j]  +tau[i]))
                              + T[j-1]**3  * np.exp(D*(-tau[j-1]+tau[i]))) * 0.5 * dT
            #Fup2[i] +=   sB* (T[j]**4      * np.exp(-tau[j]  +tau[i])
                              # + T[j-1]**4  * np.exp(-tau[j-1]+tau[i])) * 0.5 * dtau

    return Fup


@njit
def compute_fluxes_down(T, tau):
    # compute downward fluxes
    
    D   = 1.5          # Eddington approximation
    N   = len(T)
    Fdown = np.zeros(N)
    
    for i in range(N):
        Fdown[i]  = sB*T[i]**4 + (0. - sB*T[0]**4) * np.exp(-D*tau[i])

        for j in range(1, i):
            dT   = T[j]   - T[j-1]
            dtau = tau[j] - tau[j-1]
            
            Fdown[i] -= 4*sB* (T[j]**3      * np.exp(D*(-tau[i] +tau[j]))
                               + T[j-1]**3  * np.exp(D*(-tau[i] +tau[j-1]))) * 0.5 * dT
            #Fup2[i] +=   sB* (T[j]**4      * np.exp(-tau[j]  +tau[i])
                              # + T[j-1]**4  * np.exp(-tau[j-1]+tau[i])) * 0.5 * dtau

    return Fdown
