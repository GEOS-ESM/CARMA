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
plt.style.use('pcase13')

# CARMA bin parameters:
nbins = 100
rmrat = 3.0000000
rmin = 1.0e-12 # m
rhop = 1923. # kg m-3

# CARMA bin parameters
nbins = 24
rmrat = 3.7515201
rmin = 2.6686863e-10
rhop = 1923. # kg m-3

# CARMA bins
rmass, rmassup, r, rup, dr, rlow, masspart = carmabins.carmabins(nbins, rmrat, rmin, rhop)

# Geospatial parameters
lat_min = [12, 12, 12, 12, 12, 12]
lat_max = [17, 17, 17, 17, 17, 17]
lon_min = [90, 90, 75, 65, 40, 20]
lon_max = [100, 100, 85, 75, 50, 30]

# Read in CARMA data
variables = carmabins.get_variables('su', nbins, 3)
files = ['/discover/nobackup/pcase1/c90Fc_I33pac_full00/holding/tavg3d_carma_v/199108/c90Fc_I33pac_full00.tavg3d_carma_v.19910815_0900z.nc4',]
files.append('/discover/nobackup/pcase1/c90Fc_I33pac_full00/holding/tavg3d_carma_v/199108/c90Fc_I33pac_full00.tavg3d_carma_v.19910816_0900z.nc4')
files.append('/discover/nobackup/pcase1/c90Fc_I33pac_full00/holding/tavg3d_carma_v/199108/c90Fc_I33pac_full00.tavg3d_carma_v.19910817_0900z.nc4')
files.append('/discover/nobackup/pcase1/c90Fc_I33pac_full00/holding/tavg3d_carma_v/199108/c90Fc_I33pac_full00.tavg3d_carma_v.19910818_0900z.nc4')
files.append('/discover/nobackup/pcase1/c90Fc_I33pac_full00/holding/tavg3d_carma_v/199108/c90Fc_I33pac_full00.tavg3d_carma_v.19910819_0900z.nc4')
files.append('/discover/nobackup/pcase1/c90Fc_I33pac_full00/holding/tavg3d_carma_v/199108/c90Fc_I33pac_full00.tavg3d_carma_v.19910820_0900z.nc4')
k = 0
for file in files:
    print(file)
    airdensity = geos_io.fetch_var(file, 'AIRDENS')
    rh = geos_io.fetch_var(file, 'RH')
    data = geos_io.fetch_vars(file, variables)
    data = np.asarray(data)[:,0,:,:,:] # Day, Bin, vertical, latitude, longitude
    airdensity = np.asarray(airdensity)[0,:,:,:]
    rh = np.asarray(rh)[0,:,:,:]

    lat = np.load('/discover/nobackup/pcase1/working_data/lats.npy')
    lon = np.load('/discover/nobackup/pcase1/working_data/lons.npy')
    alt = np.load('/discover/nobackup/pcase1/working_data/72alts.npy')

    r_global = np.repeat(r[:,np.newaxis], rh.shape[0], axis=1)
    r_global = np.repeat(r_global[:,:,np.newaxis], rh.shape[1], axis=2)
    r_global = np.repeat(r_global[:,:,:,np.newaxis], rh.shape[2], axis=3)
    rwet = grow_v75(rh, r_global*100)
    rwet = np.swapaxes(rwet, 0,1)*0.01

    rwet = rwet[:, :, (lat > lat_min[k]) & (lat < lat_max[k]), :]
    rwet = rwet[:, :, :, (lon > lon_min[k]) & (lon < lon_max[k])]

    data = data[:, :, (lat > lat_min[k]) & (lat < lat_max[k]), :]
    data = data[:, :, :, (lon > lon_min[k]) & (lon < lon_max[k])]

# RH bins
    rh = ['00', '50', '70', '80', '90', '95', '99']

# Lognormal parameters
    rm = [0.089e-6, 0.118e-6, 0.129e-6, 0.138e-6, 0.155e-6, 0.178e-6, 0.251e-6] # median radius in meters
    sigma = np.e**(0.60) # lognormal spread parameter
    N0 = 1.00000 #lognormal number parameter

# Lognormal parameters 2
    rm2 = [0.058e-6, 0.068e-6, 0.078e-6, 0.088e-6, 0.108e-6, 0.136e-6, 0.229e-6] # median radius in meters
    sigma2 = np.e**(1.60) # lognormal spread parameter
    N02 = 1.00000 #lognormal number parameter

# Lognormal parameters 3
    rm3 = [0.35e-6, 1.4*0.35e-6, 1.5*0.35e-6, 1.6*0.35e-6, 1.8*0.35e-6, 1.9*0.35e-6, 4.6*0.35e-6] # median radius in meters
    sigma3 = np.e**(1.25) # lognormal spread parameter
    N03 = 1.00000 #lognormal number parameter

# Lognormal parameters 4
    rm4 = [0.0695e-6, 1.4*0.0695e-6, 1.5*0.0695e-6, 1.6*0.0695e-6, 1.8*0.0695e-6, 1.9*0.0695e-6, 2.2*0.0695e-6] # median radius in meters
    sigma4 = np.e**(2.03) # lognormal spread parameter
    N04 = 1.00000 #lognormal number parameter

# Calculate log normal in CARMA bins
    dNdr = np.zeros((len(rh), nbins))
    dNdr2 = np.zeros((len(rh), nbins))
    dNdr3 = np.zeros((len(rh), nbins))
    dNdr4 = np.zeros((len(rh), nbins))
    for i in range(len(rh)):
        dNdr[i,:] = N0/((2*np.pi)**0.5*r*np.log(sigma)) * \
                np.exp(-1* (np.log(r) - np.log(rm[i]))**2 / (2*np.log(sigma)**2))
        dNdr2[i,:] = N02/((2*np.pi)**0.5*r*np.log(sigma2)) * \
                np.exp(-1* (np.log(r) - np.log(rm2[i]))**2 / (2*np.log(sigma2)**2))
        dNdr3[i,:] = N03/((2*np.pi)**0.5*r*np.log(sigma3)) * \
                np.exp(-1* (np.log(r) - np.log(rm3[i]))**2 / (2*np.log(sigma3)**2))
        dNdr4[i,:] = N04/((2*np.pi)**0.5*r*np.log(sigma4)) * \
                np.exp(-1* (np.log(r) - np.log(rm4[i]))**2 / (2*np.log(sigma4)**2))

# Print bin values in # m-3
    print("Effective radius (m):")
    print((np.nansum(dNdr[0]*dr*r**3))/(np.nansum(dNdr[0]*dr*r**2)))
    print("CARMA bin radii (m):")
    print(r)
    print("---------------------------------------------------------------------")
    print("CARMA bin values (# m-3):")
    print(dNdr[0]*dr)
    print("---------------------------------------------------------------------")
    print("Total number concentration (# m-3)")
    print(np.sum(dNdr[0]*dr))

    plt.figure(figsize=(10,10))
    for i in range(len(rh)):
        if i == 0:
            plt.loglog(r, dNdr[i,:]/np.max(dNdr[i,:]), color='#5D6C89', label='GC Strat Dry', lw=3)
            plt.loglog(r, dNdr3[i,:]/np.max(dNdr3[i,:]), color='#CB7D36', label='GOCART Strat Dry', lw=3)
        #    plt.loglog(r, dNdr2[i,:]/np.max(dNdr2[i,:]), color='black', label='GC Sulfate Dry', lw=3)
        #    plt.loglog(r, dNdr4[i,:]/np.max(dNdr4[i,:]), color='#36D6E7', label='GOCART Sulfate Dry', lw=3)
        #elif i == len(rh)-1:
        #    plt.loglog(r, dNdr[i,:]/np.max(dNdr[i,:]), '--', color='#5D6C89', lw=3)
        #    plt.loglog(r, dNdr3[i,:]/np.max(dNdr3[i,:]), '--', color='#CB7D36', lw=3)
        #    plt.loglog(r, dNdr2[i,:]/np.max(dNdr2[i,:]), '--', color='black', lw=3)
        #    plt.loglog(r, dNdr4[i,:]/np.max(dNdr4[i,:]), '--', color='#36D6E7', lw=3)
#    else:
#        plt.loglog(r, dNdr[i,:]/np.max(dNdr[i,:]), '--', color='grey')
#        plt.loglog(r, dNdr2[i,:]/np.max(dNdr2[i,:]), color='grey')
    print(rwet.shape)
    print(data.shape)
    print(np.nanmean(rwet, axis=(0,2,3)).shape)
    print(np.nanmean(data, axis=(1,2,3)).shape)
    plt.loglog(np.nanmean(rwet[32,:,:,:]/3, axis=(1,2)), np.nanmean(data[:,32,:,:], axis=(1,2))/np.max(np.nanmean(data[:,32,:,:], axis=(1,2))), color='#36D6E7', label='CARMA', lw=3)
    #plt.loglog(np.nanmean(rwet[32,:,:,:], axis=(1,2)), color='#36D6E7', label='CARMA', lw=3)
    plt.legend()
    plt.ylim(1e-4, 1e1)
    plt.xlim(1e-10, 1e-4)
    plt.xlabel('Radius ($m$)')
    plt.ylabel('dN/dr ($m^{-3}$)')
    plt.title('Normalized Number Distribution Function')
    plt.tight_layout()
    plt.savefig('plots/old_' + file[-13:-10] + '_asd.png')
    #plt.savefig('background.png')
    plt.cla()
    k+=1

