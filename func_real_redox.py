import os
import numpy as np
from   petitRADTRANS.radtrans import Radtrans
from   petitRADTRANS.config   import petitradtrans_config_parser

# replace this path with the actual path to your input_data directory
# print(os.environ["pRT_INPUT_DATA_PATH"])
petitradtrans_config_parser.set_input_data_path('/Users/yoshi/petitRADTRANS/input_data')

from func_gray import *

# each array contains information in the order of H2, H2O, CO, CO2, N2
kappa         = np.array([1.3e-4, 1.e-2,  1.3e-4, 1.3e-4, 0]) # m2/kg
molar_mass = np.array([     2,   18,    28,    44,    28]) * 1e-3  # kg/mol
heat_capacity = np.array([3.5,   4.,   3.5,   3.5,   3.5]) * R     # J/(K mol)


def compute_pRT_fluxes(T0, P0, xv0, x_atmos_nc, g):
    
    radtrans = Radtrans(
        pressures = P0*1e-5,                   # enter pressure in bar
        line_species=['H2','H2O','CO','CO2'],
        gas_continuum_contributors=['N2--N2', 'H2--H2'],
        wavelength_boundaries=[0.1, 100]        # wavelength in micron (mu-m)
    )
    
    
    # mean molar mass is in g/mol. multiple by 1000 for kg -> g conversion
    kg2g = 1000
    mean_molar_masses = np.zeros(len(T0))
    for i in range(len(T0)):
        # prepare an array of atmospheric composition w/ water
        xH2O_array    = np.zeros(len(x_atmos_nc))
        xH2O_array[1] = xv0[i]
        x_atmos       = xH2O_array + (1-xv0[i]) * x_atmos_nc

        # calc average molar mass - take average by molar ratio
        mu_avg = np.sum(molar_mass * x_atmos)     
        mean_molar_masses[i] = mu_avg * kg2g

        # calc mass ratio
        r_atmos = molar_mass * x_atmos / mu_avg

        rH2  = r_atmos[0]
        rH2O = r_atmos[1]
        rCO  = r_atmos[2]
        rCO2 = r_atmos[3]
        rN2  = r_atmos[4]
        # print(i, "\t", xv0[i], mean_molar_masses[i])

    
    # set mass fraction
    mass_fractions = {
        'H2O': rH2O,
        'H2':  rH2,
        'CO':  rCO,
        'CO2': rCO2,
        'N2':  rN2,
    }
    
    # calculate flux using petitRADTRANS
    lmb, flux, ex = radtrans.calculate_flux(
        temperatures=T0,
        mass_fractions=mass_fractions,
        mean_molar_masses = mean_molar_masses,
        reference_gravity = g*100,
        return_opacities = True
    )
    
    erg2J, cm2m = 1e-7, 1e-2
    OLR = np.trapz(flux, lmb) * erg2J / cm2m / cm2m
    
    ''' test using Blackbody radiation 
    Teff = 252.36910891723454
    
    lmb_SI = lmb*0.01
    fi     = 2*h*c*c / (lmb_SI**5 * (np.exp(h*c/(lmb_SI*kB*Teff)) - 1.)) # in SI unit
    OLRi = np.trapz(fi, lmb_SI)
    # print(f"{sB * Teff**4:.1f}", "\t OLR:", OLRi * np.pi)

    # plot in micron - erg/(cm2*s) / cm
    #ax.semilogx(lmb*1e4, flux)
    #ax.semilogx(lmb*1e4, fi / erg2J * cm2m**3)
    '''

    # k = radtrans.k.shape()
    
    return lmb, flux, OLR

