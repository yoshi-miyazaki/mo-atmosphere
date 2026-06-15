import numpy as np
import matplotlib.cm      as cm
import matplotlib.pyplot  as plt
import matplotlib.patches as ptch
from   matplotlib import rc
from   matplotlib.font_manager import FontProperties
from   matplotlib.ticker       import MultipleLocator, FormatStrFormatter
import sys

plt.rcParams['font.family']='sans-serif'
plt.rcParams['axes.linewidth']= 1.0
plt.rcParams['font.size']     = 14
plt.rcParams['figure.figsize']= 5*2, 4
plt.rcParams['lines.linewidth'] = 1.3
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text',usetex=True)
rc('text.latex', preamble=r'\usepackage{sfmath}')

# run
# -- export pRT_input_data_path="/Users/yoshi/moai/MOAI/data/input_data"
# first!

from func_flux      import *
from func_partition import *

N = 100
P = np.zeros(N)

'''
for i in range(N):
    Ts = 3000 + i * 100
    dF = calc_dF(Ts, Tp, d, radiation)
'''
# CO2 partitioning
xCO2_i = 1e-4
PCO2_i = MEM * xCO2_i / (4*np.pi*rE*rE/g + MEM*4.4e-12)
print(PCO2_i/1e5, " bar ", PCO2_i*4.4e-12*100, " wt%")

xH2O_i = calc_xH2O(1., 0.5*MEO, MEM, 0.1, rE, g, muv)
PH2O_i = (xH2O_i/AH2O)**(1./0.7) 
print(PH2O_i/1e5, " bar")
# 3.5 - 10, 16 - 100

muavg = (PCO2_i*muCO2 + PH2O_i*muv)/(PCO2_i + PH2O_i) 

MatCO2_i = (4*np.pi*rE*rE / g) * (PCO2_i*muCO2)/muavg 
MatH2O_i = (4*np.pi*rE*rE / g) * (PH2O_i*muv  )/muavg
Mat_i    = MatCO2_i + MatH2O_i

# set initial atmosphere
Psurf = Mat_i * g / (4*np.pi*rE*rE)
xCO2, xH2O =  0.01, 0.99 #PCO2_i/Psurf, PH2O_i/Psurf

# set-up grid
Pmin, Pmax = 1e0, Psurf
for i in range(N):
    # set up pressure
    P[i] = Pmin * (Pmax/Pmin)**(i/(N-1))


radiation = Radtrans(line_species = ['H2O','CO2'], #_all_iso_HITEMP'],
                     wavelength_boundaries = [0.1, 251.],
                     rayleigh_species=['H2'],
                     gas_continuum_contributors=['H2-H2'],
                     pressures = P*1e-5)

dt = 10.*yr2sec

M = 300
t, v_Tp, v_Ts, v_d, v_Fconv = np.zeros(M), np.zeros(M), np.zeros(M), np.zeros(M), np.zeros(M)

# set initial condition
v_Tp[0] = 4000
v_d[0]  = 2900e3
for i in range(M-1):
    if (i==0):
        t[i] = 0
    else:
        t[i] = t[i-1]+dt
        
        
    # calculate surface temperature
    # -- using bisection search (F_conv = F_rad from the top of the atmosphere)
    Ts0, Ts1 = v_Tp[i]-500, v_Tp[i]
    (dF0, tmp, tmprd) = calc_dF(Ts0, v_Tp[i], v_d[i], Psurf, xCO2, xH2O, radiation)
    (dF1, tmp, tmprd) = calc_dF(Ts1, v_Tp[i], v_d[i], Psurf, xCO2, xH2O, radiation)

    eps = 0.01
    while (Ts1-Ts0 > eps):
        TsA = (Ts0 + Ts1)*0.5
        dFA, Fconv, rad = calc_dF(TsA, v_Tp[i], v_d[i], Psurf, xCO2, xH2O, radiation)

        if (dFA*dF0 < 0):
            Ts1, dF1 = TsA, dFA
        else:
            Ts0, dF0 = TsA, dFA

    v_Ts[i]    = TsA
    v_Fconv[i] = Fconv
    
    print("bisection check \t", dF0, dF1, "\t", TsA, "\t", Fconv)
    mu = clt/rad.freq/1e-4
    fl = rad.flux/1e-6    
    
    # thermal evolution
    dEdt = Fconv * (4*np.pi*rE*rE)
    dTdt = dEdt/7e27  # heat capacity
    
    v_Tp[i+1] = v_Tp[i] - dTdt * (dt)
    #v_d[i+1]  = v_d[i]
    v_d[i+1] = 2900e3 - max(0,3000-v_Tp[i+1])/1200*2900e3
    # save results
    output = np.array([t/yr2sec, v_Ts, v_Tp, v_Fconv])
    np.savetxt("./d_timeevo-H2O-200bar_c.txt", np.transpose(output), delimiter=",")

    if (i > 10):
        dt = 5/dTdt

v_Ts[-1] = v_Ts[-2]


# plot
fig, ax = plt.subplots(1,2)
ax[0].plot(v_Tp)
ax[0].plot(v_Ts, ls="--")
#ax[0].set_ylim([np.max(P),np.min(P)])

ax[1].semilogx(mu, fl)
ax[1].set_ylim([0., 2.5])

plt.tight_layout()
plt.savefig("./fig_MO-cooling.pdf")

