import numpy              as     np
import matplotlib.pyplot  as     plt
from func_const           import *
from petitRADTRANS.radtrans        import Radtrans


def calc_dF(Ts, Tp, d, Ptot, xCO2, xH2O, radiation):
    
    Ra = alpha*rhom*g*(Tp-Ts)*(d**3)/(kappa * nu)
    F  = 0.089*kT*(Tp-Ts)/d * (Ra**(1./3))

    FO, radiation = calc_flux(Ts, F, Ptot, xCO2, xH2O, radiation)
    
    #print(Ts, "\t", F, "\t OLR: ", FO)

    return F-FO, F, radiation


''' 
calculate the outgoing longwave radiation
'''
def calc_flux(Tsurf, Fout, Pmax, xCO2_i, xH2O_i, radiation):
    # calculate flux
    N = 100
    T, P, tau, Fup, mu, kpa, Cp, xCO2, xH2O = np.zeros(N), np.zeros(N), np.zeros(N), np.zeros(N), np.zeros(N), np.zeros(N), np.zeros(N), np.zeros(N), np.zeros(N)

    # set-up grid
    Pmin, Tmin = 1e0, 200   # P in Pa, T in K
    
    # 
    for i in range(N):
        # set up pressure
        P[i] = Pmin * (Pmax/Pmin)**(i/(N-1))
        
        # initial guess for species
        xCO2[i] = xCO2_i
        xH2O[i] = xH2O_i        
        mu[i]   = mun*(1-xCO2[i]-xH2O[i]) + muCO2*xCO2[i] + muv*xH2O[i]
        Cp[i]   = Cpn*(1-xCO2[i]-xH2O[i]) + CpCO2*xCO2[i] + Cpv*xH2O[i]


    # lapse rate
    Th = Tsurf
    T[-1] = Tsurf
    for i in range(N-1, 0, -1):
        rho = P[i]/(R*Th)*mu[i]
    
        Gamma = mu[i]/(rho*Cp[i])  # g/kpa[i] removed
                        
        Th  -= Gamma * (P[i]-P[i-1])  # originall (tau[i] - tau[i-1])
        T[i-1] = max(Tmin, Th)

        # set up opacity
        mu[i]  = mun*(1-xCO2[i]-xH2O[i]) + muCO2*xCO2[i] + muv*xH2O[i]
        kpa[i] = (kpav*muv*xH2O[i] + kpaCO2*muCO2*xCO2[i])/mu[i]

    # optical depth
    for i in range(N-1):
        tau[i+1] = tau[i] + kpa[i]/g*(P[i+1]-P[i])
    
    # simple flux
    for i in range(N):
        Fup[i] = sb*(T[i]**4)
        for j in range(N-1, i, -1):
            Fup[i] -= np.exp(-1.5*(tau[j]-tau[i])) * sb * (T[j-1]**4 - T[j]**4)

    '''
    if (i>0):
    print(Ts1, "\t", -np.trapz(tmprd.flux,tmprd.freq)*1e-7/1e-4)
    i = 1
    print(i, tau[i], T[i], "\t", Fup[i], " W/m2")
    
    fig, ax = plt.subplots(1,2)
    ax[0].semilogy(T,P)
    ax[0].set_ylim([np.max(P),np.min(P)])

    ax[1].semilogy(T,tau)
    ax[1].set_ylim([np.max(tau),np.min(tau)])

    plt.tight_layout()
    plt.savefig("./fig_atmos.pdf")
    '''

    # Now, use petitRADTRANS
    # pRT reads pressures in bar! (and equidistant in pressure log-space)
    #radiation.setup_opa_structure(P*1e-5)
    
    # store species info
    # this is in mass fraction -- NOT number fraction...
    mass_fractions = {}
    mass_fractions['H2O'] = 0 #xH2O * muv  /mu
    mass_fractions['CO2'] = 0.01 #xCO2 * muCO2/mu #_all_iso_HITEMP
    #mass_fractions['H2'] = 0.99

    for i in range(N):
        mu[i]  = 2e-3*0.99 + muCO2*0.01
    #mass_fractions['CO_all_iso_HITEMP'] = xCO2 * muCO2/mu
    
    # calc flux
    w, flux, _ = radiation.calculate_flux(temperatures = T,                   # again, flip the profile!
                                          mass_fractions =mass_fractions,
                                          reference_gravity = g*100,               # WARNING! pRT reads gravity in cgs!
                                          mean_molar_masses = mu)                  # again, flip the profile!

    
    # calc temp
    erg2J  = 1e-7 # ergs to Joules conversion factor
    cm22m2 = 1e-4 # cm² to m² converion factor
    OLR    = np.trapz(flux,(w*1e-2))*erg2J/cm22m2
    print(Tsurf, 'OLR =',OLR,'W/m2, \t', Fup[0], " W/m2", "\t", Fup[0]/OLR)

    return OLR, radiation
