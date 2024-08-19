"""Script for making CARMA bin size distributions

Author: Parker A. Case
Changelog: Created 2/1/2019
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
import carmabins
plt.style.use('pcase13')

# CARMA bin parameters:
nbins = 24
rmrat = 3.7515201
rmin = 2.6686863e-10 # m
rhop = 1923. # kg m-3

# CARMA bins
rmass, rmassup, r, rup, dr, rlow, masspart = carmabins.carmabins(nbins, rmrat, rmin, rhop)

# Lognormal parameters
rm = 0.086e-6 # median radius in meters
sigma = 1.59 # lognormal spread parameter
N0 = 1.06183 #lognormal number parameter

# Calculate log normal in CARMA bins
dNdr = N0/((2*np.pi)**0.5*r*np.log(sigma)) * \
        np.exp(-1* (np.log(r) - np.log(rm))**2 / (2*np.log(sigma)**2))

# Print bin values in # m-3
print("Effective radius (m):")
print((np.nansum(dNdr*dr*r**3))/(np.nansum(dNdr*dr*r**2)))
print("CARMA bin radii (m):")
print(r)
print("---------------------------------------------------------------------")
print("CARMA bin values (# m-3):")
print(dNdr*dr)
print("---------------------------------------------------------------------")
print("Total number concentration (# m-3)")
print(np.sum(dNdr*dr))

