"""Primary interface for running PARMA column sulfate only setup

This is an interface for running PARMA with a pure sulfate group. The
run_column() function can be used to run these simulations. The main
method at the bottom of this file gives an example of using the
run_column() function.

For more information on how to run PARMA, see the README in the interfaces
directory.

Author: Parker Case
Version: v1.0 (2024/08/20) Implemented function format from dust column
                            and other features implemented there
         v0.2 (2024/08/14) Readability, RH output, air density output,
                            now using netcdf4
         v0.1 (2023/10/12) First commit
"""
import numpy as np
from time import perf_counter
import netCDF4
from stdatmosphere import get_standard_atmosphere_1d
from carma_column_sulfate import carma_column_sulfate
import carma_tools.carmabins as carmabins
from carma_tools.wetr import grow_v75

def run_column(nz, dt, nt, nt_carma, nbins, rmrat, rmin, rhop, h2o, h2so4, \
                su, constant_h2so4 = False, exp_name = '0'):
    """Runs a CARMA column with the specified parameters

    Arguments:
    nz -- number of layers to simulate in the column
    dt -- length of timestep (seconds)
    nt -- number of PARMA timesteps to simulate
    nt_carma -- number of CARMA timesteps per PARMA timestep
    nbins -- number of CARMA bins
    rmrat -- ratio between bins in the pure sulfate group
    rmin -- minimum radius in the pure sulfate group (cm)
    rhop -- density of sulfate (g cm-3)
    h2o -- array (nz,) of water vapor mass mixing ratios (g g-1)
    h2so4 -- array (nz,) of h2so4 vapor mass mixing ratios (g g-1)
    su -- array (nz, nbins) of sulfate mass mixing ratios (g g-1)
    constant_h2so4 -- boolean indicating whether CARMA should deplete h2so4
    """
    # Prep the rest of the variables for CARMA
    alt = np.linspace(0,nz*1000+.01,nz)
    p, t = get_standard_atmosphere_1d(alt*1000)
    su_0 = su.copy()
    h2so4_0 = h2so4.copy()
    h2o_0 = h2o.copy()


    ###########################################################################
    # Run the model!
    ###########################################################################
    su_lib = np.zeros((nt, nz, nbins))
    h2so4_lib = np.zeros((nt, nz))
    rh_lib = np.zeros((nt, nz))
    rhoa_lib = np.zeros((nt, nz))
    time_0 = perf_counter()
    for i in range(nt):
        su_out, t_out, p_out, h2so4_out, h2o_out, rhoa_out, rh_out = \
                carma_column_sulfate(rmrat, rmin, rhop, t, p, h2so4, h2o, su, \
                dt, 1, constant_h2so4, nbins, nz)
        su_lib[i,:,:] = su_out
        h2so4_lib[i,:] = h2so4_out
        rh_lib[i,:] = rh_out
        rhoa_lib[i,:] = rhoa_out
        su = su_out
        h2so4 = h2so4_out
        h2o = h2o_out
        t = t_out
        p = p_out
    time_1 = perf_counter()
    print('CARMA took ' + str(time_1 - time_0) + ' seconds to run '
          + str(nt) + ' timesteps.')
    print('Time per timestep = ' + str((time_1-time_0)/nt) + ' seconds')

    ###########################################################################
    # Convert to mKs units for analysis:
    ###########################################################################
    rhop = rhop * 1000 # kg m-3
    rmin = rmin * 1e-2 # m
    rhoa_lib = rhoa_lib * 1000 # kg m-3

    ###########################################################################
    # Build CARMA bin structure for analysis
    ###########################################################################
    rmass, rmassup, r, rup, dr, rlow, masspart = carmabins.carmabins(nbins, rmrat,
                                                                     rmin, rhop)
    dlogr = np.log10(rup) - np.log10(rlow)
    dndlogr_lib = su_lib / dlogr / (4/3 * np.pi * r**3 * rhop) * np.repeat(rhoa_lib[:,:,np.newaxis], nbins, axis=2)
    dndlogr_0   = su_0   / dlogr / (4/3 * np.pi * r**3 * rhop) * np.repeat(rhoa_lib[0,:,np.newaxis], nbins, axis=1)
    dadlogr_lib = dndlogr_lib * 4 * np.pi * r**2
    dadlogr_0   = dndlogr_0   * 4 * np.pi * r**2
    dmdlogr_lib = dndlogr_lib * (4/3) * np.pi * r**3 * rhop
    dmdlogr_0   = dndlogr_0   * (4/3) * np.pi * r**3 * rhop

    ###########################################################################
    # Create .nc4 file
    ###########################################################################
    ncfile = netCDF4.Dataset('parma_column_sulfate_' + exp_name + '.nc4', mode='w')
    ncfile.createDimension('layer', nz)
    ncfile.createDimension('bin', nbins)
    ncfile.createDimension('time', nt+1)
    ncfile.title='PARMA Column Output'
    ncfile.subtitle='Created by parma_column.py'

    ncr = ncfile.createVariable('r', np.float64, ('bin'))
    ncr.long_name = 'Bin (geometric) center radius'
    ncr.units = 'm'
    ncr[:] = r

    ncalt = ncfile.createVariable('alt', np.float64, ('layer'))
    ncalt.long_name = 'Altitude'
    ncalt.units = 'm'
    ncalt[:] = np.linspace(0,(nz-1)*1000,nz)+500

    nctime = ncfile.createVariable('time', np.float64, ('time'))
    nctime.long_name = 'Simulated time elapsed'
    nctime.units = 's'
    nctime[:] = np.linspace(0,dt*nt,nt+1)

    ncdndlogr = ncfile.createVariable('dNdlogr', np.float64, ('time', 'layer', 'bin'))
    ncdndlogr.long_name = 'Number concentration as a function of log10(r)'
    ncdndlogr.units = 'm-3'
    ncdndlogr[0,:,:] = dndlogr_0
    ncdndlogr[1:,:,:] = dndlogr_lib

    ncdadlogr = ncfile.createVariable('dAdlogr', np.float64, ('time', 'layer', 'bin'))
    ncdadlogr.long_name = 'Surface area concentration as a function of log10(r)'
    ncdadlogr.units = 'm2 m-3'
    ncdadlogr[0,:,:] = dadlogr_0
    ncdadlogr[1:,:,:] = dadlogr_lib

    ncdmdlogr = ncfile.createVariable('dMdlogr', np.float64, ('time', 'layer', 'bin'))
    ncdmdlogr.long_name = 'Mass concentration as a function of log10(r)'
    ncdmdlogr.units = 'kg m-3'
    ncdmdlogr[0,:,:] = dmdlogr_0
    ncdmdlogr[1:,:,:] = dmdlogr_lib

    ncmmr = ncfile.createVariable('MMR', np.float64, ('time', 'layer', 'bin'))
    ncmmr.long_name = 'Aerosol mass mixing ratio'
    ncmmr.units = 'kg kg-1'
    ncmmr[0,:,:] = su_0
    ncmmr[1:,:,:] = su_lib

    ncmass = ncfile.createVariable('M', np.float64, ('time', 'layer'))
    ncmass.long_name = 'Total aerosol mass concentration'
    ncmass.units = 'kg m-3'
    ncmass[0,:] = np.sum(dmdlogr_0 * dlogr, axis=1)
    ncmass[1:,:] = np.sum(dmdlogr_lib * dlogr, axis=2)

    ncarea = ncfile.createVariable('SA', np.float64, ('time', 'layer'))
    ncarea.long_name = 'Total aerosol surface area concentration'
    ncarea.units = 'm2 m-3'
    ncarea[0,:] = np.sum(dadlogr_0 * dlogr, axis=1)
    ncarea[1:,:] = np.sum(dadlogr_lib * dlogr, axis=2)

    ncn = ncfile.createVariable('N', np.float64, ('time', 'layer'))
    ncn.long_name = 'Total aerosol number concentration'
    ncn.units = '# m-3'
    ncn[0,:] = np.sum(dndlogr_0 * dlogr, axis=1)
    ncn[1:,:] = np.sum(dndlogr_lib * dlogr, axis=2)

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


if __name__ == '__main__':
    ############################################################################
    # Example simulation parameters
    ############################################################################
    # Time and grid parameters
    nz    = 30
    dt    = 900 # timestep (900 seconds (15 minutes))
    nt    = 768 # number of timesteps (8 days)
    nt_carma = 1 # number carma timesteps per parma "timestep"

    # CARMA bin parameters
    nbins = 24 # number of bins
    rmrat = 3.75125201 # mass ratio between bins
    rmin = 2.6686863e-8 # cm
    rhop = 1.923 # g cm-3

    # Background atmospheric parameters
    constant_h2so4 = False
    h2o   = np.zeros(nz) + 1e-7 # H2O mmr (kg/kg)
    h2so4 = np.zeros(nz) + 4.71e-13 # H2SO4 mmr (kg/kg)
    su = np.zeros((nz,nbins)) # Sulfate aerosol mmr (kg/kg)
    for i in range(nz):
        su[i,:] = [4.2579942e-28, 1.5634370e-27, 1.1281235e-26,
                   7.5493089e-24, 4.6745146e-22, 1.1954208e-20,
                   1.5612010e-19, 1.3430505e-18, 2.1747372e-17,
                   1.5991898e-15, 6.5651301e-14, 1.1735656e-12,
                   1.3214742e-11, 7.1381352e-11, 2.0359722e-10,
                   3.1819103e-10, 1.3236401e-10, 1.2854880e-11,
                   3.6080232e-13, 2.5853792e-15, 4.4601091e-18,
                   2.1512019e-21, 2.8390088e-25, 5.5137755e-29]

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

    run_column(nz, dt, nt, nt_carma, nbins, rmrat, rmin, rhop, h2o, h2so4, \
                su, constant_h2so4 = constant_h2so4)
