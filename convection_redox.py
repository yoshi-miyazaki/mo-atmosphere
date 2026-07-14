import numpy as np
import matplotlib.pyplot as plt
from   matplotlib import rc
from   matplotlib.font_manager import FontProperties
from   matplotlib.ticker       import MultipleLocator, FormatStrFormatter

plt.rcParams['font.family']='sans-serif'
plt.rcParams['axes.linewidth']= 1.0
plt.rcParams['font.size']     = 10
plt.rcParams['figure.figsize']= 4*1.618*4, 4.

from func_real_redox      import *
from func_partition_redox import *


''' set planet-related parameters here '''
# planet size
ME = 5.972e24        # planet mass
rE = 6371e3          # planet radius (m)

Mp = ME * 1.
rp = rE * (Mp/ME)**0.26
g  = G*Mp/(rp*rp)    # gravity       (m/s2)

 

# atmopshperic mass & composition
# M_atmos  = 1e20     # atmospheric mass (kg)
xCO2     = 9e-1     # CO2 fraction     (molar ratio)
rH       = 1.

x_atmos_nc = np.array([0., xCO2, 1.-xCO2])


#
def calc_Tconv_profile(T_surf, M_atmos, x_atmos, rp, g):
    
    # heat capacity, opacity
    mu_avg = np.sum(molar_mass * x_atmos)
    Cp     = np.sum(heat_capacity * x_atmos)  # take average by molar ratio
    kpa    = np.sum(kappa * x_atmos * molar_mass) / mu_avg
    
    gamma = R/Cp

    print(kappa, x_atmos)
    
    # Top & surface pressures
    P_surf = M_atmos * g / (4 * np.pi * rp*rp)
    P_top  = 1 #XS0.01 * g/kpa
    
    N = 200
    T_vals, P_vals, tau_vals = np.zeros(N), np.zeros(N), np.zeros(N)
    for i in range(N):
        P_vals[i] = P_top * (P_surf/P_top)**(i/(N-1))
        T_vals[i] = T_surf * (P_vals[i]/P_surf)**gamma

    tau_vals[0] = 1e-6
    for i in range(N-1):
        tau_vals[i+1] = tau_vals[i] + (P_vals[i+1]-P_vals[i]) * kpa/g * P_vals[i]/P_ref

    return T_vals, P_vals, tau_vals

def calc_Tprofile(Tsurf, M_atmos, x_atmos, rp, g):

    # first, assume  adiabatic structure
    T0, P0, tau0 = calc_Tconv_profile(Tsurf, M_atmos, x_atmos, rp, g)    

    # calculate flux (gray)
    Fgray_up = compute_fluxes(T0, tau0)
    Fgray_dn = compute_fluxes_down(T0, tau0)
    
    
    # check where the net-IR heating rate turns negative for a gray atmosphere
    dF = Fgray_up - Fgray_dn
    
    tropopause_index = 0
    for i in range(len(dF)-1):
        IR_net = dF[i+1] - dF[i]
        #print(i, f"{P0[i]:.2e}\t{tau0[i]:.2e}\t{IR_net:.2f}")

        if (IR_net < 0 and abs(IR_net) > dF[i]*0.01):
            tropopause_index = i
            
            for j in range(i, 0, -1):
                IR_net = dF[j+1] - dF[j]
                print(j, f"* {P0[j]:.2e}\t{tau0[j]:.2e}\t{IR_net:.2f}")

                if (IR_net > 0):
                    tropopause_index = j+1
                    break
                
            break

    # return a troposphere
    Tt   = T0[j+1:]
    Pt   = P0[j+1:]
    taut = tau0[j+1:]

    return Tt, Pt, taut, tropopause_index


'''
# estimate modern N2 amount
M_atmos_now = 101325 * (4*np.pi*rp*rp) / g
x_atmos_now = np.array([.2, 0., .8])
r_atmos_now = x_atmos_now * molar_mass
r_atmos_now = r_atmos_now / r_atmos_now.sum()

MN2 = M_atmos_now * r_atmos_now[2]
'''

# set parameters
#rp = rE

# set atmospheric compsition
#
# For Fig. 1
# xCO2_init, xH2O_init = 2.49e-4, 1.047e-3  # for CHILI

# IW+4
xC_base, xH_base = 2.e-4, 4.e-4
xCO_init, xCO2_init = xC_base * 28/44 * 0., xC_base * 1.
xH2_init, xH2O_init = xH_base  * 2/18 * 0., xH_base * 1.

# IW+1
#xC_base, xH_base = 2.e-4, 1.e-3
#xCO_init, xCO2_init = xC_base * 28/44 * 0.77, xC_base * 0.23
#xH2_init, xH2O_init = xH_base * 2/18 * 0.4, xH_base * 0.6

MN2 = 4.54349737285796e+18 # modern N2 amount

M_mantle = 4.01e24 * Mp/ME

# calcualte?
P_atmos = calc_atmos_composition(xH2_init, xH2O_init, xCO_init, xCO2_init, MN2, M_mantle, rp, g)
x_atmos = P_atmos / P_atmos.sum()
M_atmos = P_atmos.sum() * (4*np.pi*rp*rp) / g

# weight ratio
r_atmos = x_atmos * molar_mass
r_atmos = r_atmos / r_atmos.sum()
print(f"Atmospheric mass =", M_atmos * r_atmos)
print(f"Atmospheric P    =", P_atmos / 1e5)

#sys.exit()

# (non-gray)
x_atmos_nc = np.copy(x_atmos)
x_atmos_nc[1] = 0
x_atmos_nc = x_atmos_nc / x_atmos_nc.sum()

# sys.exit()

# calculate T profile
Tsurf_vals = np.linspace(1600, 1000, 61) # 401
Ftop_vals  = np.zeros(len(Tsurf_vals))
Fgray_vals = np.zeros(len(Tsurf_vals))

for i, Tsurf in enumerate(Tsurf_vals):
    T0, P0, tau0, tp0= calc_Tprofile(Tsurf, M_atmos, x_atmos, rp, g)

    # print(T0, P0, tau0, tp0)
    
    xv = np.ones_like(T0) * x_atmos[1]
    lmb, flux_top, Fr_conv = compute_pRT_fluxes(T0, P0, xv, x_atmos_nc, g)

    Fgray_up = compute_fluxes(T0, tau0)
    Fgray_dn = compute_fluxes_down(T0, tau0)

    Ftop_vals[i]  = Fr_conv
    Fgray_vals[i] = Fgray_up[0]

data_output = np.vstack([Tsurf_vals, Ftop_vals, Fgray_vals]).T
np.savetxt(f"./data_flux/fluxM{int(Mp/ME*10)}_C{int(xCO2_init*1e4)}-H{int(xH2O_init*1e4)}.txt", data_output, fmt="%.6e", delimiter="\t")
    
fig, ax = plt.subplots(1,4)

# plot thermal structure
ax[0].semilogy(T0, P0, color="b")
ax[0].set(ylabel = "pressure")
ax[0].set_ylim([np.max(P0), np.min(P0)])

ax[1].semilogy(T0, tau0, color="b")
ax[1].set(ylabel = "optical depth")
ax[1].set_ylim([np.max(tau0), np.min(tau0)])


ax[2].semilogy(Fgray_up - Fgray_dn, tau0, color="k", alpha = .5, ls = "--")
ax[2].set_xlim([0, 1.5*(Fgray_up[0]-Fgray_dn[0])])
ax[2].set_ylim([np.max(tau0), np.min(tau0)])

#ax[2].set_xscale('log')


# spectrum
ax[3].semilogx(lmb*1e4, flux_top)
ax[3].set(xlabel="wavelength [$\mu$m]")

# reference blackbody spectrum
for Tref in [T0[0], T0[7]]:
    lmb_SI = lmb*0.01
    fi     = np.pi * 2*h*c*c / (lmb_SI**5 * (np.exp(h*c/(lmb_SI*kB*Tref)) - 1.)) # in SI unit
    OLRi   = np.trapz(fi, lmb_SI)
        
    # plot in micron - erg/(cm2*s) / cm
    erg2J, cm2m = 1e-7, 1e-2
    ax[3].semilogx(lmb*1e4, fi / erg2J * cm2m**3)
    ax[3].set_ylim([0, np.max(fi / erg2J * cm2m**3)*1.1])

# ax[3].set_xscale('linear')
plt.tight_layout()
plt.savefig("./fig_conv.pdf")
