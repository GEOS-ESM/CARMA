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
from carma_box import carma_box
import carma_tools.carmabins as carmabins
from carma_tools.wetr import grow_v75

# Simulation parameters
dt    = 900 # timestep (900 seconds (15 minutes))
nt    = 96  # number of timesteps (1 days)
t     = np.asarray((283,)) # temperature (kelvin)
p     = np.asarray((100000,))   # pressure (bar)
h2o   = np.asarray((1e-2,))  # H2O mmr (kg/kg)
h2so4 = np.asarray((2.5e-14,))  # H2SO4 mmr (kg/kg)
substeps = np.asarray((0.,))
supersat = np.asarray((0.,))
constant_h2so4 = False # Keep H2SO4 at initial value?
su = np.asarray(
    [1.00928194e-17, 3.58698844e-15, 2.40942977e-14, 3.54669229e-14,
     5.91490968e-14, 1.42859235e-13, 3.65422536e-13, 9.66915133e-13,
     3.87790624e-12, 1.93728714e-11, 8.03272934e-11, 2.93220032e-10,
     5.02852038e-10, 3.62854635e-10, 1.61361605e-10, 5.27720402e-11,
     1.20563759e-11, 2.17942027e-12, 4.07592471e-13, 6.44141587e-14,
     4.42018723e-15, 1.32088171e-16, 2.83900880e-25, 5.51377550e-29]
     ) # sulfate aerosol mmr (kg/kg)
su_0 = su.copy()

# CARMA bin structure
nbins = 24
rmrat = 3.75125201
rmin = 2.6686863e-8 # cm
rhop = 1.923 # g cm-3

# History
substeps_history = []
h2so4_history = []
supersat_history = []
temperature_history = []

# Run the model!
time_0 = perf_counter()
for ti in range(nt):
    carma_box(rmrat, rmin, rhop, t, p, h2so4, h2o, su, dt, 1, constant_h2so4, 6, substeps, supersat)
    h2so4_history.append(h2so4[0])
    substeps_history.append(substeps[0])
    supersat_history.append(supersat[0])
    temperature_history.append(t[0])
time_1 = perf_counter()
print('CARMA took ' + str(time_1 - time_0) + ' seconds to run '
      + str(nt) + ' timesteps.')
print('Time per timestep = ' + str((time_1-time_0)/nt) + ' seconds')

hours = [dt/60/60 * n for n in range(nt)]
print(supersat_history)

fig, axs = plt.subplots(3,1)
ax = axs[0]
ax.plot(hours,substeps_history, marker='x',color='black')
ax.set_ylabel('substeps')
ax.set_xlabel('hours')
ax.grid()

ax2 = axs[1]
ax2.plot(hours,h2so4_history, marker='x', color='red')
ax2.set_ylabel('$H_2SO_4$ concentration', color='red')
ax2.set_yscale('log')
ax2.set_xlabel('hours')
ax2.grid()

ax3 = axs[2]
ax3.plot(hours,supersat_history, marker='x', color='blue')
ax3.set_ylabel('$H_2SO_4$ supersaturation', color='blue')
ax3.axhline(1, ls='--', color='deepskyblue')
ax3.set_xlabel('hours')
ax3.grid()

plt.tight_layout()
plt.show()


# Convert to mKs units for analysis:
rhop = rhop * 1000 # kg m-3
rmin = rmin * 1e-2 # m

# Build CARMA bin structure for analysis
airdensity = 0.188 # kg m-3
rmass, rmassup, r, rup, dr, rlow, masspart = carmabins.carmabins(nbins, rmrat,
                                                                 rmin, rhop)
dndr_in  = su_0 / dr / (4/3 * np.pi * r**3 * rhop) * airdensity
dndr_out = su / dr / (4/3 * np.pi * r**3 * rhop) * airdensity

# Plot size distributions
fig, axs = plt.subplots(3,1, figsize=(8,8))

axs[0].loglog(r, dndr_in, '+', color='grey', label='t = 0 hours')
axs[0].loglog(r, dndr_out, 'x', color='black',
              label='t = ' + str(nt*dt/3600) + ' hours')
axs[0].legend()
axs[0].set_title('Number Distribution (dN/dr)')
axs[0].set_xlabel('Radius ($m$)')
axs[0].set_ylabel('# / dr ($m^{-1}$)')
axs[0].grid()

axs[1].loglog(r, 4*np.pi*r**2*dndr_in, '+', color='grey', label='t = 0 hours')
axs[1].loglog(r, 4*np.pi*r**2*dndr_out, 'x', color='black',
              label='t = ' + str(nt*dt/3600) + ' hours')
axs[1].legend()
axs[1].set_title('Surface Area Distribution (dA/dr)')
axs[1].set_xlabel('Radius ($m$)')
axs[1].set_ylabel('Surface Area / dr ($m^2\ m^{-1}$)')
axs[1].grid()

axs[2].loglog(r, 4/3*np.pi*r**3*dndr_in, '+', color='grey', label='t = 0 hours')
axs[2].loglog(r, 4/3*np.pi*r**3*dndr_out, 'x', color='black',
              label='t = ' + str(nt*dt/3600) + ' hours')
axs[2].legend()
axs[2].set_title('Volume Distribution (dV/dr)')
axs[2].set_xlabel('Radius ($m$)')
axs[2].set_ylabel('Volume / dr ($m^3\ m^{-1}$)')
axs[2].grid()

plt.tight_layout()
plt.show()
