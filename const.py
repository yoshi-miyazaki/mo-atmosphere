#!/usr/bin/env python
# coding: utf-8

# Universal constants
sigma = 5.670374419e-8
h = 6.626e-34
sB = 5.670374419e-8  #Named different in new code you use, eiter sigma or sB works.

G   = 6.674e-11 #m^3/(kg*s^2)
kB  = 1.380649e-23    # Boltzmann constant [J/K]
R   = 8.314           # Gas constant [J/(mol*K)]
N_A = 6.02214076e23   # Avogrado number
c   = 299792458

grav = 9.81 #m/s^2
lH2O= 2.5e6 #J/(kg*K)
lH2O_mol=43655 #J/mol
Racr=1000
T_sun=5780 #K
R_sun = 6.957e8 #m
d_sun_earth = 1.496e11#m
M_earth=5.97219e24 #kg
R_earth=6.371e6 #m
# Temp w/o atmosphere
F_sun=sigma*T_sun**4*(R_sun/d_sun_earth)**2 #Solar flux on Earth

mH2O = 18.01528e-3 #kg/mol
mCO2 = 44.009e-3 #kg/mol
mOTHER = 28.97e-3 #kg/mol
gammaH2O=1.3
gammaCO2=1.3
gammaOTHER=1.4
cpH2O = 4 * R/mH2O  #J/(kg*K)
cpCO2 = 3.5 * R/mCO2  #J/(kg*K)
cpOTHER = 3.5 * R/mOTHER #J/(kg*K)
CpH2O=4*R #J/(mol*K)
CpCO2=3.5*R #J/(mol*K)
CpOTHER=3.5*R #J/(mol*K)
kpaH2O = 1e-2 #unitless
kpaCO2 = 1.3e-4 #unitless


# In[ ]:




