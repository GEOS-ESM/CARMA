"""Example script for analysis of the PARMA column sulfate model

This script should read in the output of parma_column_sulfate.py and plot up
some variables.

Author: Parker Case
Version: v0.1 (2024/08/19) First commit
"""
import numpy as np
import matplotlib.pyplot as plt
import netCDF4

################################################################################
# Analysis parameters
################################################################################
alt_idx = 25 # layer to analyze for distribution plots

################################################################################
# Fetch simulation output
################################################################################
data = netCDF4.Dataset('parma_column_sulfate.nc4')

alt = data.variables['alt'][:]/1000
time = data.variables['time'][:]/60/60/24

m    = data.variables['M'][:]
n    = data.variables['N'][:]
sa   = data.variables['SA'][:]

dndlogr   = data.variables['dNdlogr'][:]
dadlogr   = data.variables['dAdlogr'][:]
dmdlogr   = data.variables['dMdlogr'][:]
r = data.variables['r'][:]
r *= 1e6


################################################################################
# Create plot
################################################################################
fig, axs = plt.subplots(2,3, figsize=(12,8))

# Total number contour plots
cf = axs[0,0].contourf(time, alt, n.T, cmap='magma_r')
axs[0,0].set_title('Number: SU')
plt.colorbar(cf, ax=axs[0,0], label='$m^{-3}$')

# Total surface area contour plots
cf = axs[0,1].contourf(time, alt, sa.T, cmap='magma_r')
axs[0,1].set_title('SA: SU')
plt.colorbar(cf, ax=axs[0,1], label='$m^2\ m^{-3}$')

# Total mass contour plots
cf = axs[0,2].contourf(time, alt, m.T, cmap='magma_r')
axs[0,2].set_title('Mass: SU')
plt.colorbar(cf, ax=axs[0,2], label='$kg\ m^{-3}$')

[axs[0,i].set_ylabel('km') for i in range(3)]
[axs[0,i].set_xlabel('days') for i in range(3)]

# Time dependent number distribution
axs[1,0].loglog(r, dndlogr[0,alt_idx,:], color='blue', label='SU Initial')
axs[1,0].loglog(r, dndlogr[-1,alt_idx,:], '--', color='blue', \
        label='SU t = ' + str(time[-1]))

upper_lim = np.max(dndlogr)
axs[1,0].set_ylim((1e-3,upper_lim))
axs[1,0].set_title('Number distribution')
axs[1,0].set_ylabel('dN/dlog(r)')
axs[1,0].set_xlabel('r ($\mu m$)')
axs[1,0].grid()

# Time dependent surface area distribution
axs[1,1].loglog(r, dadlogr[0,alt_idx,:], color='blue', label='SU Initial')
axs[1,1].loglog(r, dadlogr[-1,alt_idx,:], '--', color='blue', \
        label='SU t = ' + str(time[-1]))

upper_lim = np.max(dadlogr)
axs[1,1].set_ylim((1e-10,upper_lim))
axs[1,1].set_title('Surface area distribution')
axs[1,1].set_ylabel('dA/dlog(r)')
axs[1,1].set_xlabel('r ($\mu m$)')
axs[1,1].grid()

# Time dependent mass distribution
axs[1,2].loglog(r, dmdlogr[0,alt_idx,:], color='blue', label='SU Initial')
axs[1,2].loglog(r, dmdlogr[-1,alt_idx,:], '--', color='blue', \
        label='SU t = ' + str(time[-1]))

upper_lim = np.max(dmdlogr)
axs[1,2].set_ylim((1e-15,upper_lim))
axs[1,2].set_title('Mass distribution')
axs[1,2].set_ylabel('dM/dlog(r)')
axs[1,2].set_xlabel('r ($\mu m$)')
axs[1,2].grid()

axs[1,2].legend(bbox_to_anchor=(1.05,1), loc='upper left', borderaxespad=0.)

plt.tight_layout()
plt.show()
