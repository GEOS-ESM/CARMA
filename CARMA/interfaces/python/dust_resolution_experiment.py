import numpy as np
import matplotlib.pyplot as plt
import carma_tools.carmabins as carmabins
from carma_tools.dust_tools import kok_size_distribution
from parma_column_dust import run_column

colors = ['#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02']
markers = ['x','d','^','o']
verbose = False

# Time and grid parameters
nz       = 30 # number of (1km) layers
dt       = 900 # timestep (900 seconds (15 minutes))
nt       = 3 #960 # number of timesteps (10 days)
nt_carma = 960 # number carma timesteps per parma "timestep"

# CARMA bin parameters
nbins_list = [24,47,93,185] # number of bins

rmrat_su_list = [3.75125201 ** x for x in [1., 0.5, 0.25, 0.125]]
rmin_su = 2.6686863e-8 # cm

rmrat_mx_list = [2.2587828 ** x for x in [1., 0.5, 0.25, 0.125]] # mass ratio between bins
rmin_mx = 5.e-06 # cm

rhop_su = 1.923 # g cm-3
rhop_du = 2.65 # g cm-3

constant_h2so4 = False # Keep H2SO4 at initial value?
h2o   = np.zeros(nz) + 1e-7 # H2O mmr (kg/kg)
h2so4 = np.zeros(nz) + 4.71e-13 # H2SO4 mmr (kg/kg)


if verbose is True:
    fig, axs = plt.subplots(1,3)

i = 0
for nbins, rmrat_su, rmrat_mx in zip(nbins_list, rmrat_su_list, rmrat_mx_list):
    print('Running: nbins = ' + str(nbins))

    rmass_su, rmassup_su, r_su, rup_su, dr_su, rlow_su, masspart_su = \
        carmabins.carmabins(nbins, rmrat_su, rmin_su, rhop_su)
    rmass_mx, rmassup_mx, r_mx, rup_mx, dr_mx, rlow_mx, masspart_mx = \
        carmabins.carmabins(nbins, rmrat_mx, rmin_mx, rhop_du)
    if nbins == 24:
        r_su_24 = r_su.copy()

    dlogr_su = np.log(rup_su) - np.log(rlow_su)
    dlogr_mx = np.log(rup_mx) - np.log(rlow_mx)

    kok_dm = kok_size_distribution(nbins, r_mx, rlow_mx, rup_mx, rhop_du, rhop_du)

    rhoa = 0.118

    su = np.zeros((nz,nbins)) # Pure sulfate aerosol mmr (kg/kg)
    for i in range(nz):
        su[i,:] = np.interp(r_su, r_su_24,
                            [4.2579942e-28, 1.5634370e-27, 1.1281235e-26,
                             7.5493089e-24, 4.6745146e-22, 1.1954208e-20,
                             1.5612010e-19, 1.3430505e-18, 2.1747372e-17,
                             1.5991898e-15, 6.5651301e-14, 1.1735656e-12,
                             1.3214742e-11, 7.1381352e-11, 2.0359722e-10,
                             3.1819103e-10, 1.3236401e-10, 1.2854880e-11,
                             3.6080232e-13, 2.5853792e-15, 4.4601091e-18,
                             2.1512019e-21, 2.8390088e-25, 5.5137755e-29])
    mxdu = np.zeros((nz,nbins)) + 1e-6 * kok_dm*dlogr_mx # Mixed group dust aerosol mmr (kg/kg)
    mxsu = np.zeros((nz,nbins)) + mxdu + 1e-20*dlogr_mx # Mixed group sulfate aerosol mmr (kg/kg)
    su = su * dlogr_su

    if verbose is True:
        axs[0].loglog(r_su*1e-2, su[0], marker=markers[i], color=colors[0])
        axs[1].loglog(r_mx*1e-2, mxdu[0], marker=markers[i], color=colors[1])
        axs[2].loglog(r_mx*1e-2, mxsu[0], marker=markers[i], color=colors[2])

    i += 1

    run_column(nz, dt, nt, nt_carma, nbins, rmrat_su, rmin_su, rmrat_mx, \
                rmin_mx, rhop_su, rhop_du, h2o, h2so4, su, mxsu, mxdu, \
                constant_h2so4 = constant_h2so4, exp_name = str(nbins))

if verbose is True:
    plt.show()
