"""Example script for a PARMA column model

This is an example of how to run the CARMA column model interface using carma_column,
compiled using f2py. For more information on how to run PARMA, see the README
in the interfaces directory.

TODO: Add RH swelling and wet radius output
      Make a non-sulfate (no condensing gas) version

Author: Parker Case
Version: v0.2 (2024/08/14) Readability, RH output, air density output,
                            now using netcdf4
         v0.1 (2023/10/12) First commit
"""
import numpy as np
from time import perf_counter
import netCDF4
from stdatmosphere import get_standard_atmosphere_1d
from carma_column_dust import carma_column_dust
import carma_tools.carmabins as carmabins
from carma_tools.wetr import grow_v75

################################################################
# Simulation parameters
################################################################
# Time and grid parameters
nz       = 30 # number of (1km) layers
dt       = 900 # timestep (900 seconds (15 minutes))
nt       = 768 # number of timesteps (8 days)
nt_carma = 1 # number carma timesteps per parma "timestep"

# CARMA bin parameters
nbins = 24 # number of bins

rmrat_su = 3.75125201 # mass ratio between bins
rmin_su = 2.6686863e-8 # cm

rmrat_mx = 2.2587828 # mass ratio between bins
rmin_mx = 5.e-06 # cm

rhop_su = 1.923 # g cm-3
rhop_du = 2.65 # g cm-3

# Background atmospheric parameters
h2o   = np.zeros(nz) + 1e-7 # H2O mmr (kg/kg)
h2so4 = np.zeros(nz) + 4.71e-13 # H2SO4 mmr (kg/kg)
su = np.zeros((nz,nbins)) # Pure sulfate aerosol mmr (kg/kg)
for i in range(nz):
    su[i,:] = [4.2579942e-28, 1.5634370e-27, 1.1281235e-26,
               7.5493089e-24, 4.6745146e-22, 1.1954208e-20,
               1.5612010e-19, 1.3430505e-18, 2.1747372e-17,
               1.5991898e-15, 6.5651301e-14, 1.1735656e-12,
               1.3214742e-11, 7.1381352e-11, 2.0359722e-10,
               3.1819103e-10, 1.3236401e-10, 1.2854880e-11,
               3.6080232e-13, 2.5853792e-15, 4.4601091e-18,
               2.1512019e-21, 2.8390088e-25, 5.5137755e-29]
mxsu = np.zeros((nz,nbins)) + 1e-27 # Mixed group sulfate aerosol mmr (kg/kg)
mxdu = np.zeros((nz,nbins)) + 1e-27 # Mixed group dust aerosol mmr (kg/kg)

# Add plume 24-26 km
h2o[24:26] = 125.1e-6
h2so4[24:26] = 4.71e-10
su[24:26, :] = [4.2579942e-28, 1.5634370e-27, 1.1281235e-26,
               7.5493089e-24, 4.6745146e-22, 1.1954208e-20,
               1.5612010e-19, 1.3430505e-18, 2.1747372e-17,
               1.5991898e-15, 6.5651301e-14, 1.1735656e-12,
               1.3214742e-11, 7.1381352e-11, 2.0359722e-10,
               3.1819103e-10, 1.3236401e-10, 1.2854880e-11,
               3.6080232e-13, 2.5853792e-15, 4.4601091e-18,
               2.1512019e-21, 2.8390088e-10, 5.5137755e-29]

# Prep the rest of the variables for CARMA
alt = np.linspace(0,nz*1000+.01,nz)
p, t = get_standard_atmosphere_1d(alt*1000)
constant_h2so4 = False # Keep H2SO4 at initial value?
su_0 = su.copy()
mxsu_0 = mxsu.copy()
mxdu_0 = mxdu.copy()
h2so4_0 = h2so4.copy()
h2o_0 = h2o.copy()


################################################################
# Run the model!
################################################################
su_lib = np.zeros((nt, nz, nbins))
mxsu_lib = np.zeros((nt, nz, nbins))
mxdu_lib = np.zeros((nt, nz, nbins))
h2so4_lib = np.zeros((nt, nz))
rh_lib = np.zeros((nt, nz))
rhoa_lib = np.zeros((nt, nz))
time_0 = perf_counter()
for i in range(nt):
    su_out, mxsu_out, mxdu_out, t_out, p_out, h2so4_out, h2o_out, rhoa_out, rh_out = carma_column_dust(rmrat_su, rmrat_mx, rmin_su, rmin_mx, rhop_su, rhop_du, t, p, h2so4, h2o, su, mxsu, mxdu, dt, 1, constant_h2so4, nbins, nz)
    su_lib[i,:,:] = su_out
    mxsu_lib[i,:,:] = mxsu_out
    mxdu_lib[i,:,:] = mxdu_out
    h2so4_lib[i,:] = h2so4_out
    rh_lib[i,:] = rh_out
    rhoa_lib[i,:] = rhoa_out
    su = su_out
    mxsu = mxsu_out
    mxdu = mxdu_out
    h2so4 = h2so4_out
    h2o = h2o_out
    t = t_out
    p = p_out
time_1 = perf_counter()
print('CARMA took ' + str(time_1 - time_0) + ' seconds to run '
      + str(nt) + ' timesteps.')
print('Time per timestep = ' + str((time_1-time_0)/nt) + ' seconds')

################################################################
# Convert to mKs units for analysis:
################################################################
rhop_su = rhop_su * 1000 # kg m-3
rhop_du = rhop_du * 1000 # kg m-3
rmin_su = rmin_su * 1e-2 # m
rmin_mx = rmin_mx * 1e-2 # m
rhoa_lib = rhoa_lib * 1000 # kg m-3

################################################################
# Build CARMA bin structure for analysis
################################################################
rmass_su, rmassup_su, r_su, rup_su, dr_su, rlow_su, masspart_su = \
    carmabins.carmabins(nbins, rmrat_su, rmin_su, rhop_su)
rmass_mx, rmassup_mx, r_mx, rup_mx, dr_mx, rlow_mx, masspart_mx = \
    carmabins.carmabins(nbins, rmrat_mx, rmin_mx, rhop_du)

dlogr_su = np.log10(rup_su) - np.log10(rlow_su)
dndlogr_su_lib = su_lib / dlogr_su / (4/3 * np.pi * r_su**3 * rhop_su) * np.repeat(rhoa_lib[:,:,np.newaxis], nbins, axis=2)
dndlogr_su_0   = su_0   / dlogr_su / (4/3 * np.pi * r_su**3 * rhop_su) * np.repeat(rhoa_lib[0,:,np.newaxis], nbins, axis=1)
dadlogr_su_lib = dndlogr_su_lib * 4 * np.pi * r_su**2
dadlogr_su_0   = dndlogr_su_0   * 4 * np.pi * r_su**2
dmdlogr_su_lib = dndlogr_su_lib * (4/3) * np.pi * r_su**3 * rhop_su
dmdlogr_su_0   = dndlogr_su_0   * (4/3) * np.pi * r_su**3 * rhop_su

dlogr_mx = np.log10(rup_mx) - np.log10(rlow_mx)

dndlogr_mxsu_lib = mxsu_lib / dlogr_mx / (4/3 * np.pi * r_mx**3 * rhop_su) * np.repeat(rhoa_lib[:,:,np.newaxis], nbins, axis=2)
dndlogr_mxsu_0   = mxsu_0   / dlogr_mx / (4/3 * np.pi * r_mx**3 * rhop_su) * np.repeat(rhoa_lib[0,:,np.newaxis], nbins, axis=1)
dadlogr_mxsu_lib = dndlogr_mxsu_lib * 4 * np.pi * r_mx**2
dadlogr_mxsu_0   = dndlogr_mxsu_0   * 4 * np.pi * r_mx**2
dmdlogr_mxsu_lib = dndlogr_mxsu_lib * (4/3) * np.pi * r_mx**3 * rhop_su
dmdlogr_mxsu_0   = dndlogr_mxsu_0   * (4/3) * np.pi * r_mx**3 * rhop_su

dndlogr_mxdu_lib = mxdu_lib / dlogr_mx / (4/3 * np.pi * r_mx**3 * rhop_du) * np.repeat(rhoa_lib[:,:,np.newaxis], nbins, axis=2)
dndlogr_mxdu_0   = mxdu_0   / dlogr_mx / (4/3 * np.pi * r_mx**3 * rhop_du) * np.repeat(rhoa_lib[0,:,np.newaxis], nbins, axis=1)
dadlogr_mxdu_lib = dndlogr_mxdu_lib * 4 * np.pi * r_mx**2
dadlogr_mxdu_0   = dndlogr_mxdu_0   * 4 * np.pi * r_mx**2
dmdlogr_mxdu_lib = dndlogr_mxdu_lib * (4/3) * np.pi * r_mx**3 * rhop_du
dmdlogr_mxdu_0   = dndlogr_mxdu_0   * (4/3) * np.pi * r_mx**3 * rhop_du

################################################################
# Create .nc4 file
################################################################
ncfile = netCDF4.Dataset('parma_column_dust.nc4',mode='w')
ncfile.createDimension('layer', nz)
ncfile.createDimension('bin', nbins)
ncfile.createDimension('time', nt+1)
ncfile.title='PARMA Column Output'
ncfile.subtitle='Created by parma_column_dust.py'

ncr = ncfile.createVariable('r_su', np.float64, ('bin'))
ncr.long_name = 'Sulfate bin (geometric) center radius'
ncr.units = 'm'
ncr[:] = r_su

ncr = ncfile.createVariable('r_mx', np.float64, ('bin'))
ncr.long_name = 'Mixed bin (geometric) center radius'
ncr.units = 'm'
ncr[:] = r_mx

ncalt = ncfile.createVariable('alt', np.float64, ('layer'))
ncalt.long_name = 'Altitude'
ncalt.units = 'm'
ncalt[:] = np.linspace(0,(nz-1)*1000,nz)+500

nctime = ncfile.createVariable('time', np.float64, ('time'))
nctime.long_name = 'Simulated time elapsed'
nctime.units = 's'
nctime[:] = np.linspace(0,dt*nt,nt+1)

ncdndlogr = ncfile.createVariable('dNdlogr_su', np.float64, ('time', 'layer', 'bin'))
ncdndlogr.long_name = 'Number concentration as a function of log10(r)'
ncdndlogr.units = 'm-3'
ncdndlogr[0,:,:] = dndlogr_su_0
ncdndlogr[1:,:,:] = dndlogr_su_lib

ncdadlogr = ncfile.createVariable('dAdlogr_su', np.float64, ('time', 'layer', 'bin'))
ncdadlogr.long_name = 'Surface area concentration as a function of log10(r)'
ncdadlogr.units = 'm2 m-3'
ncdadlogr[0,:,:] = dadlogr_su_0
ncdadlogr[1:,:,:] = dadlogr_su_lib

ncdmdlogr = ncfile.createVariable('dMdlogr_su', np.float64, ('time', 'layer', 'bin'))
ncdmdlogr.long_name = 'Mass concentration as a function of log10(r)'
ncdmdlogr.units = 'kg m-3'
ncdmdlogr[0,:,:] = dmdlogr_su_0
ncdmdlogr[1:,:,:] = dmdlogr_su_lib

ncmmr = ncfile.createVariable('MMR_su', np.float64, ('time', 'layer', 'bin'))
ncmmr.long_name = 'Aerosol mass mixing ratio'
ncmmr.units = 'kg kg-1'
ncmmr[0,:,:] = su_0
ncmmr[1:,:,:] = su_lib

ncmass = ncfile.createVariable('M_su', np.float64, ('time', 'layer'))
ncmass.long_name = 'Total aerosol mass concentration'
ncmass.units = 'kg m-3'
ncmass[0,:] = np.sum(dmdlogr_su_0 * dlogr_su, axis=1)
ncmass[1:,:] = np.sum(dmdlogr_su_lib * dlogr_su, axis=2)

ncarea = ncfile.createVariable('SA_su', np.float64, ('time', 'layer'))
ncarea.long_name = 'Total aerosol surface area concentration'
ncarea.units = 'm2 m-3'
ncarea[0,:] = np.sum(dadlogr_su_0 * dlogr_su, axis=1)
ncarea[1:,:] = np.sum(dadlogr_su_lib * dlogr_su, axis=2)

ncn = ncfile.createVariable('N_su', np.float64, ('time', 'layer'))
ncn.long_name = 'Total aerosol number concentration'
ncn.units = '# m-3'
ncn[0,:] = np.sum(dndlogr_su_0 * dlogr_su, axis=1)
ncn[1:,:] = np.sum(dndlogr_su_lib * dlogr_su, axis=2)

ncdndlogr = ncfile.createVariable('dNdlogr_mxsu', np.float64, ('time', 'layer', 'bin'))
ncdndlogr.long_name = 'Number concentration as a function of log10(r)'
ncdndlogr.units = 'm-3'
ncdndlogr[0,:,:] = dndlogr_mxsu_0
ncdndlogr[1:,:,:] = dndlogr_mxsu_lib

ncdadlogr = ncfile.createVariable('dAdlogr_mxsu', np.float64, ('time', 'layer', 'bin'))
ncdadlogr.long_name = 'Surface area concentration as a function of log10(r)'
ncdadlogr.units = 'm2 m-3'
ncdadlogr[0,:,:] = dadlogr_mxsu_0
ncdadlogr[1:,:,:] = dadlogr_mxsu_lib

ncdmdlogr = ncfile.createVariable('dMdlogr_mxsu', np.float64, ('time', 'layer', 'bin'))
ncdmdlogr.long_name = 'Mass concentration as a function of log10(r)'
ncdmdlogr.units = 'kg m-3'
ncdmdlogr[0,:,:] = dmdlogr_mxsu_0
ncdmdlogr[1:,:,:] = dmdlogr_mxsu_lib

ncmmr = ncfile.createVariable('MMR_mxsu', np.float64, ('time', 'layer', 'bin'))
ncmmr.long_name = 'Aerosol mass mixing ratio'
ncmmr.units = 'kg kg-1'
ncmmr[0,:,:] = mxsu_0
ncmmr[1:,:,:] = mxsu_lib

ncmass = ncfile.createVariable('M_mxsu', np.float64, ('time', 'layer'))
ncmass.long_name = 'Total aerosol mass concentration'
ncmass.units = 'kg m-3'
ncmass[0,:] = np.sum(dmdlogr_mxsu_0 * dlogr_mx, axis=1)
ncmass[1:,:] = np.sum(dmdlogr_mxsu_lib * dlogr_mx, axis=2)

ncarea = ncfile.createVariable('SA_mxsu', np.float64, ('time', 'layer'))
ncarea.long_name = 'Total aerosol surface area concentration'
ncarea.units = 'm2 m-3'
ncarea[0,:] = np.sum(dadlogr_mxsu_0 * dlogr_mx, axis=1)
ncarea[1:,:] = np.sum(dadlogr_mxsu_lib * dlogr_mx, axis=2)

ncn = ncfile.createVariable('N_mxsu', np.float64, ('time', 'layer'))
ncn.long_name = 'Total aerosol number concentration'
ncn.units = '# m-3'
ncn[0,:] = np.sum(dndlogr_mxsu_0 * dlogr_mx, axis=1)
ncn[1:,:] = np.sum(dndlogr_mxsu_lib * dlogr_mx, axis=2)

ncdndlogr = ncfile.createVariable('dNdlogr_mxdu', np.float64, ('time', 'layer', 'bin'))
ncdndlogr.long_name = 'Number concentration as a function of log10(r)'
ncdndlogr.units = 'm-3'
ncdndlogr[0,:,:] = dndlogr_mxdu_0
ncdndlogr[1:,:,:] = dndlogr_mxdu_lib

ncdadlogr = ncfile.createVariable('dAdlogr_mxdu', np.float64, ('time', 'layer', 'bin'))
ncdadlogr.long_name = 'Surface area concentration as a function of log10(r)'
ncdadlogr.units = 'm2 m-3'
ncdadlogr[0,:,:] = dadlogr_mxdu_0
ncdadlogr[1:,:,:] = dadlogr_mxdu_lib

ncdmdlogr = ncfile.createVariable('dMdlogr_mxdu', np.float64, ('time', 'layer', 'bin'))
ncdmdlogr.long_name = 'Mass concentration as a function of log10(r)'
ncdmdlogr.units = 'kg m-3'
ncdmdlogr[0,:,:] = dmdlogr_mxdu_0
ncdmdlogr[1:,:,:] = dmdlogr_mxdu_lib

ncmmr = ncfile.createVariable('MMR_mxdu', np.float64, ('time', 'layer', 'bin'))
ncmmr.long_name = 'Aerosol mass mixing ratio'
ncmmr.units = 'kg kg-1'
ncmmr[0,:,:] = mxdu_0
ncmmr[1:,:,:] = mxdu_lib

ncmass = ncfile.createVariable('M_mxdu', np.float64, ('time', 'layer'))
ncmass.long_name = 'Total aerosol mass concentration'
ncmass.units = 'kg m-3'
ncmass[0,:] = np.sum(dmdlogr_mxdu_0 * dlogr_mx, axis=1)
ncmass[1:,:] = np.sum(dmdlogr_mxdu_lib * dlogr_mx, axis=2)

ncarea = ncfile.createVariable('SA_mxdu', np.float64, ('time', 'layer'))
ncarea.long_name = 'Total aerosol surface area concentration'
ncarea.units = 'm2 m-3'
ncarea[0,:] = np.sum(dadlogr_mxdu_0 * dlogr_mx, axis=1)
ncarea[1:,:] = np.sum(dadlogr_mxdu_lib * dlogr_mx, axis=2)

ncn = ncfile.createVariable('N_mxdu', np.float64, ('time', 'layer'))
ncn.long_name = 'Total aerosol number concentration'
ncn.units = '# m-3'
ncn[0,:] = np.sum(dndlogr_mxdu_0 * dlogr_mx, axis=1)
ncn[1:,:] = np.sum(dndlogr_mxdu_lib * dlogr_mx, axis=2)

nch2so4 = ncfile.createVariable('H2SO4', np.float64, ('time', 'layer'))
nch2so4.long_name = 'Sulfuric acid vapor concentration'
nch2so4.units = 'kg m-3'
nch2so4[0,:] = h2so4_0 * rhoa_lib[0,:]
nch2so4[1:,:] = h2so4_lib * rhoa_lib

ncrh = ncfile.createVariable('RH', np.float64, ('time', 'layer'))
ncrh.long_name = 'Relative humidity'
ncrh.units = '%'
ncrh[0,:] = rh_lib[0,:]
ncrh[1:,:] = rh_lib

ncfile.close()

