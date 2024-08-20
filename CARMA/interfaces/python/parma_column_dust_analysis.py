"""Example script for analysis of the PARMA column dust model

This script should read in the output of parma_column_dust.py and plot up
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
data = netCDF4.Dataset('parma_column_dust.nc4')

alt = data.variables['alt'][:]/1000
time = data.variables['time'][:]/60/60/24

m_mxdu  = data.variables['M_mxdu'][:]
m_mxsu  = data.variables['M_mxsu'][:]
m_su    = data.variables['M_su'][:]
n_mxdu  = data.variables['N_mxdu'][:]
n_mxsu  = data.variables['N_mxsu'][:]
n_su    = data.variables['N_su'][:]
sa_mxdu = data.variables['SA_mxdu'][:]
sa_mxsu = data.variables['SA_mxsu'][:]
sa_su   = data.variables['SA_su'][:]

dndlogr_mxdu = data.variables['dNdlogr_mxdu'][:]
dndlogr_mxsu = data.variables['dNdlogr_mxsu'][:]
dndlogr_su   = data.variables['dNdlogr_su'][:]
dadlogr_mxdu = data.variables['dAdlogr_mxdu'][:]
dadlogr_mxsu = data.variables['dAdlogr_mxsu'][:]
dadlogr_su   = data.variables['dAdlogr_su'][:]
dmdlogr_mxdu = data.variables['dMdlogr_mxdu'][:]
dmdlogr_mxsu = data.variables['dMdlogr_mxsu'][:]
dmdlogr_su   = data.variables['dMdlogr_su'][:]
r_mx = data.variables['r_mx'][:]
r_su = data.variables['r_su'][:]
r_mx *= 1e6
r_su *= 1e6


################################################################################
# Create plot
################################################################################
fig, axs = plt.subplots(4,3, figsize=(12,8))

# Total number contour plots
cf = axs[0,0].contourf(time, alt, n_mxdu.T, cmap='magma_r')
axs[0,0].set_title('Number: MXDU')
plt.colorbar(cf, ax=axs[0,0], label='$m^{-3}$')
cf = axs[0,1].contourf(time, alt, n_mxsu.T, cmap='magma_r')
axs[0,1].set_title('Number: MXSU')
plt.colorbar(cf, ax=axs[0,1], label='$m^{-3}$')
cf = axs[0,2].contourf(time, alt, n_su.T, cmap='magma_r')
axs[0,2].set_title('Number: SU')
plt.colorbar(cf, ax=axs[0,2], label='$m^{-3}$')

[axs[0,i].set_ylabel('km') for i in range(3)]
[axs[0,i].set_xlabel('days') for i in range(3)]

# Total surface area contour plots
cf = axs[1,0].contourf(time, alt, sa_mxdu.T, cmap='magma_r')
axs[1,0].set_title('SA: MXDU')
plt.colorbar(cf, ax=axs[1,0], label='$m^2\ m^{-3}$')
cf = axs[1,1].contourf(time, alt, sa_mxsu.T, cmap='magma_r')
axs[1,1].set_title('SA: MXSU')
plt.colorbar(cf, ax=axs[1,1], label='$m^2\ m^{-3}$')
cf = axs[1,2].contourf(time, alt, sa_su.T, cmap='magma_r')
axs[1,2].set_title('SA: SU')
plt.colorbar(cf, ax=axs[1,2], label='$m^2\ m^{-3}$')

[axs[1,i].set_ylabel('km') for i in range(3)]
[axs[1,i].set_xlabel('days') for i in range(3)]

# Total mass contour plots
cf = axs[2,0].contourf(time, alt, m_mxdu.T, cmap='magma_r')
axs[2,0].set_title('Mass: MXDU')
plt.colorbar(cf, ax=axs[2,0], label='$kg\ m^{-3}$')
cf = axs[2,1].contourf(time, alt, m_mxsu.T, cmap='magma_r')
axs[2,1].set_title('Mass: MXSU')
plt.colorbar(cf, ax=axs[2,1], label='$kg\ m^{-3}$')
cf = axs[2,2].contourf(time, alt, m_su.T, cmap='magma_r')
axs[2,2].set_title('Mass: SU')
plt.colorbar(cf, ax=axs[2,2], label='$kg\ m^{-3}$')

[axs[2,i].set_ylabel('km') for i in range(3)]
[axs[2,i].set_xlabel('days') for i in range(3)]

# Time dependent number distribution
axs[3,0].loglog(r_su, dndlogr_su[0,alt_idx,:], color='blue', label='SU Initial')
axs[3,0].loglog(r_su, dndlogr_su[-1,alt_idx,:], '--', color='blue', \
        label='SU t = ' + str(time[-1]))
axs[3,0].loglog(r_mx[0,alt_idx,:], dndlogr_mxdu[0,alt_idx,:], color='black', \
        label='MXDU Initial')
axs[3,0].loglog(r_mx[-1,alt_idx,:], dndlogr_mxdu[-1,alt_idx,:], '--', \
        color='black', label='MXDU t = ' + str(time[-1]))
axs[3,0].loglog(r_mx[0,alt_idx,:], dndlogr_mxsu[0,alt_idx,:], color='green', \
        label='MXSU Initial')
axs[3,0].loglog(r_mx[-1,alt_idx,:], dndlogr_mxsu[-1,alt_idx,:], '--', \
        color='green', label='MXSU t = ' + str(time[-1]))

upper_lim = np.max((dndlogr_su, dndlogr_mxsu, dndlogr_mxdu))
axs[3,0].set_ylim((1e-3,upper_lim))
axs[3,0].set_title('Number distribution')
axs[3,0].set_ylabel('dN/dlog(r)')
axs[3,0].set_xlabel('r ($\mu m$)')
axs[3,0].grid()

# Time dependent surface area distribution
axs[3,1].loglog(r_su, dadlogr_su[0,alt_idx,:], color='blue', label='SU Initial')
axs[3,1].loglog(r_su, dadlogr_su[-1,alt_idx,:], '--', color='blue', \
        label='SU t = ' + str(time[-1]))
axs[3,1].loglog(r_mx[0,alt_idx,:], dadlogr_mxdu[0,alt_idx,:], color='black', \
        label='MXDU Initial')
axs[3,1].loglog(r_mx[-1,alt_idx,:], dadlogr_mxdu[-1,alt_idx,:], '--', \
        color='black', label='MXDU t = ' + str(time[-1]))
axs[3,1].loglog(r_mx[0,alt_idx,:], dadlogr_mxsu[0,alt_idx,:], color='green', \
        label='MXSU Initial')
axs[3,1].loglog(r_mx[-1,alt_idx,:], dadlogr_mxsu[-1,alt_idx,:], '--', \
        color='green', label='MXSU t = ' + str(time[-1]))

upper_lim = np.max((dadlogr_su, dadlogr_mxsu, dadlogr_mxdu))
axs[3,1].set_ylim((1e-10,upper_lim))
axs[3,1].set_title('Surface area distribution')
axs[3,1].set_ylabel('dA/dlog(r)')
axs[3,1].set_xlabel('r ($\mu m$)')
axs[3,1].grid()

# Time dependent mass distribution
axs[3,2].loglog(r_su, dmdlogr_su[0,alt_idx,:], color='blue', label='SU Initial')
axs[3,2].loglog(r_su, dmdlogr_su[-1,alt_idx,:], '--', color='blue', \
        label='SU t = ' + str(time[-1]))
axs[3,2].loglog(r_mx[0,alt_idx,:], dmdlogr_mxdu[0,alt_idx,:], color='black', \
        label='MXDU Initial')
axs[3,2].loglog(r_mx[-1,alt_idx,:], dmdlogr_mxdu[-1,alt_idx,:], '--', \
        color='black', label='MXDU t = ' + str(time[-1]))
axs[3,2].loglog(r_mx[0,alt_idx,:], dmdlogr_mxsu[0,alt_idx,:], color='green', \
        label='MXSU Initial')
axs[3,2].loglog(r_mx[-1,alt_idx,:], dmdlogr_mxsu[-1,alt_idx,:], '--', \
        color='green', label='MXSU t = ' + str(time[-1]))

upper_lim = np.max((dmdlogr_su, dmdlogr_mxsu, dmdlogr_mxdu))
axs[3,2].set_ylim((1e-15,upper_lim))
axs[3,2].set_title('Mass distribution')
axs[3,2].set_ylabel('dM/dlog(r)')
axs[3,2].set_xlabel('r ($\mu m$)')
axs[3,2].grid()

axs[3,2].legend(bbox_to_anchor=(1.05,1), loc='upper left', borderaxespad=0.)

plt.tight_layout()
plt.show()
