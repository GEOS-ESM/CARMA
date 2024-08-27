"""Example script for analysis of the PARMA column dust model

This script should read in the output of parma_column_dust.py and plot up
some variables.

Author: Parker Case
Version: v0.1 (2024/08/19) First commit
"""
import numpy as np
import matplotlib.pyplot as plt
import netCDF4
colors = ['#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02']
markers = ['x','d','^','o']

################################################################################
# Analysis parameters
################################################################################
alt_idx = 25 # layer to analyze for distribution plots

################################################################################
# Fetch simulation output
################################################################################
data = netCDF4.Dataset('parma_column_dust_185.nc4')
m_mxdu_baseline = data.variables['M_mxdu'][:,25]
fig, axs = plt.subplots(2,2, figsize=(7,5))
gs = axs[0,0].get_gridspec()
for ax in axs[0,:]:
    ax.remove()
axbig = fig.add_subplot(gs[0,:])
i = 0
for nbins in ['24','47','93','185']:
    data = netCDF4.Dataset('parma_column_dust_'+str(nbins)+'.nc4')

    alt = data.variables['alt'][:]/1000
    time = data.variables['time'][:]*960/60/60/24

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
    cf = axbig.plot(time, m_mxdu[:,25]/m_mxdu_baseline, color=colors[i], label='nbins = ' + str(nbins))

    cf = axs[1,0].semilogx(r_mx[0,25,:], dmdlogr_mxdu[0,25,:], color=colors[i], markersize=2, marker=markers[i])

    cf = axs[1,1].semilogx(r_mx[0,25,:], dmdlogr_mxdu[-1,25,:], color=colors[i], markersize=2, marker=markers[i])

    #print(r_su)
    #print(r_su[::2])
    #print(r_su[::4])
    #print(r_su[::8])

    i += 1

axbig.legend()
axbig.set_xlabel('Days')
axbig.set_ylabel('Normalized mass')
axbig.set_title('Dust Mass as a fraction of 185 bin simulation')
axs[1,0].set_xlabel('Radius ($\mu m$)')
axs[1,0].set_ylabel('dM/dlogr')
axs[1,0].set_title('t = 0.0 days')
axs[1,0].set_xlim((1e-1,1e1))
axs[1,1].set_xlabel('Radius ($\mu m$)')
axs[1,1].set_ylabel('dM/dlogr')
axs[1,1].set_title('t = ' + str(time[-1]) + ' days')
axs[1,1].set_xlim((1e-1,1e1))
plt.tight_layout()
plt.show()
