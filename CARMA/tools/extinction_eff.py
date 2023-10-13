"""Script for making CARMA bin size distributions

Author: Parker A. Case
Changelog: Created 2/1/2019
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
import geos_io
import geos_utils
import carmabins
from wetr import grow_v75
from netCDF4 import Dataset
plt.style.use('pcase13')

# CARMA bin parameters
nbins = 24
rmrat = 3.7515201
rmin = 2.6686863e-10
rhop = 1923. # kg m-3

# CARMA bins
rmass, rmassup, r, rup, dr, rlow, masspart = carmabins.carmabins(nbins, rmrat, rmin, rhop)

# Read in extinction
ncfile = Dataset('carma_opticsBands_SU.v6.nbin=24.nc', 'r')
qext = ncfile.variables['qext'][:]

# Plot extinction efficiency
plt.figure(figsize=(10,10))
plt.semilogx(r, qext[:,0,0], color='#5D6C89', lw=3)
plt.legend()
plt.xlim(1e-10, 1e-4)
plt.xlabel('Radius ($m$)')
plt.ylabel('Extinction Efficiency')
plt.title('Extinction Efficiency')
plt.tight_layout()
plt.show()
plt.cla()

