"""Example script for PARMA, a python-wrapped CARMA box model.

This is an example of how to run the CARMA box model interface using carma_box,
compiled using f2py. For more information on how to run PARMA, see the README
in the interfaces directory.

Author: Parker Case
Version: v0.1 (2023/10/12)
"""
import numpy as np
import matplotlib.pyplot as plt
from time import perf_counter
from carma_column import carma_column
import carma_tools.carmabins as carmabins
from carma_tools.wetr import grow_v75

# Simulation parameters
dt    = 900 # timestep (900 seconds (15 minutes))
nt    = int(768/8) # number of timesteps (8 days)
nt_carma = 1 # number carma timesteps per timestep
t     = np.zeros(30) + 220 # temperature (kelvin)
p     = np.zeros(30) + 2300 # pressure (bar)
h2o   = np.zeros(30) + 1e-7 # H2O mmr (kg/kg)
h2o   = h2o[::-1]
h2o[24:26] = 125.1e-6
h2so4 = np.zeros(30) + 4.71e-13 # H2SO4 mmr (kg/kg)
#h2so4 = np.linspace(0, 4.7e-11, 30) + 4.7e-10
h2so4[24:26] = 4.71e-10
print(h2so4)
print(h2o)
constant_h2so4 = False # Keep H2SO4 at initial value?
su = np.zeros((30,24))
for i in range(30):
    su[i,:] = [4.2579942e-28, 1.5634370e-27, 1.1281235e-26, 
               7.5493089e-24, 4.6745146e-22, 1.1954208e-20, 
               1.5612010e-19, 1.3430505e-18, 2.1747372e-17,
               1.5991898e-15, 6.5651301e-14, 1.1735656e-12, 
               1.3214742e-11, 7.1381352e-11, 2.0359722e-10,
               3.1819103e-10, 1.3236401e-10, 1.2854880e-11, 
               3.6080232e-13, 2.5853792e-15, 4.4601091e-18, 
               2.1512019e-21, 2.8390088e-25, 5.5137755e-29]
su_0 = su.copy()

# CARMA bin structure
nbins = 24
rmrat = 3.75125201
rmin = 2.6686863e-8 # cm
rhop = 1.923 # g cm-3

# Run the model!
su_lib = np.zeros((nt, 30, nbins))
time_0 = perf_counter()
for i in range(nt):
    su_out, t_out, p_out, h2so4_out, h2o_out = carma_column(rmrat, rmin, rhop, t, p, h2so4, h2o, su, dt, 1, constant_h2so4, nbins)
    su_lib[i,:,:] = su_out
    su = su_out
time_1 = perf_counter()
print('CARMA took ' + str(time_1 - time_0) + ' seconds to run '
      + str(nt) + ' timesteps.')
print('Time per timestep = ' + str((time_1-time_0)/nt) + ' seconds')

# Convert to mKs units for analysis:
rhop = rhop * 1000 # kg m-3
rmin = rmin * 1e-2 # m

# Build CARMA bin structure for analysis
airdensity = 0.188 # kg m-3
rmass, rmassup, r, rup, dr, rlow, masspart = carmabins.carmabins(nbins, rmrat,
                                                                 rmin, rhop)
dndr_lib = su_lib / dr / (4/3 * np.pi * r**3 * rhop) * airdensity

fig, axs = plt.subplots(3,1)

time = np.linspace(0, nt*dt, nt)/60/60/24
levels = np.arange(-3,22.1,0.5)
cf1 = axs[0].contourf(time, r, np.log10(dndr_lib[:,25,:].T), levels, cmap='magma', extend='both')
axs[0].set_yscale('log')
plt.colorbar(cf1, ax=axs[0])

cf2 = axs[1].contourf(time, r, np.log10(dndr_lib[:,15,:].T), levels, cmap='magma', extend='both')
axs[1].set_yscale('log')
plt.colorbar(cf2, ax=axs[1])

cf3 = axs[2].contourf(time, r, np.log10(dndr_lib[:,5,:].T), levels, cmap='magma', extend='both')
axs[2].set_yscale('log')
plt.colorbar(cf3, ax=axs[2])
plt.show()
plt.cla()

fig, axs = plt.subplots(1,2)
z = np.linspace(0,30,30)
axs[0].contourf(r, z, np.log10(dndr_lib[-1,:,:]), levels, cmap='magma', extend='both')
axs[0].set_xscale('log')

axs[1].plot(t_out, z, label='Temperature')
axs[1].plot(p_out, z, label='Pressure')
axs[1].plot(h2o_out, z, label='$H_2O$')
axs[1].plot(h2so4_out, z, label='$H_2SO_4$')
plt.show()

