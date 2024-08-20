"""Primary interface for running PARMA box sulfate only setup

This is an interface for running PARMA with a pure sulfate group. The
run_box() function can be used to run these simulations. The main
method at the bottom of this file gives an example of using the
run_box() function.

For more information on how to run PARMA, see the README in the interfaces
directory.

Author: Parker Case
Version: v1.0 (2024/08/20) Added column-like structure, readability, and
                            broke it out into a function
         v0.1 (2023/10/12)
"""
import numpy as np
import matplotlib.pyplot as plt
import netCDF4
from time import perf_counter
from carma_box_sulfate import carma_box_sulfate
import carma_tools.carmabins as carmabins
from carma_tools.wetr import grow_v75

def run_box(dt, nt, nt_carma, nbins, rmrat, rmin, rhop, h2o, h2so4, \
                su, t, p, constant_h2so4 = False):
    """Runs a CARMA column with the specified parameters

    Arguments:
    dt -- length of timestep (seconds)
    nt -- number of PARMA timesteps to simulate
    nt_carma -- number of CARMA timesteps per PARMA timestep
    nbins -- number of CARMA bins
    rmrat -- ratio between bins in the pure sulfate group
    rmin -- minimum radius in the pure sulfate group (cm)
    rhop -- density of sulfate (g cm-3)
    h2o -- water vapor mass mixing ratio (g g-1)
    h2so4 -- h2so4 vapor mass mixing ratio (g g-1)
    su -- array (nbins) of sulfate mass mixing ratio (g g-1)
    t -- temperature (K)
    p -- pressure (Pa)
    constant_h2so4 -- boolean indicating whether CARMA should deplete h2so4
    """
    # Run the model!
    su_0 = su.copy()
    h2so4_0 = h2so4
    su_lib = np.zeros((nt, nbins))
    h2so4_lib = np.zeros((nt,))
    rh_lib = np.zeros((nt,))
    rhoa_lib = np.zeros((nt,))
    time_0 = perf_counter()
    for i in range(nt):
        su_out, t_out, p_out, h2so4_out, h2o_out, rhoa_out, rh_out = \
                carma_box_sulfate(rmrat, rmin, rhop, t, p, h2so4, h2o, su, \
                dt, 1, constant_h2so4, nbins)
        su_lib[i,:] = su_out
        h2so4_lib[i] = h2so4_out
        rh_lib[i] = rh_out
        rhoa_lib[i] = rhoa_out
        su = su_out
        h2so4 = h2so4_out
        h2o = h2o_out
        t = t_out
        p = p_out
    time_1 = perf_counter()
    print('CARMA took ' + str(time_1 - time_0) + ' seconds to run '
          + str(nt) + ' timesteps.')
    print('Time per timestep = ' + str((time_1-time_0)/nt) + ' seconds')

    # Convert to mKs units for analysis:
    rhop = rhop * 1000 # kg m-3
    rmin = rmin * 1e-2 # m

    # Build CARMA bin structure for analysis
    rmass, rmassup, r, rup, dr, rlow, masspart = carmabins.carmabins(nbins, rmrat,
                                                                     rmin, rhop)
    dlogr = np.log10(rup) - np.log10(rlow)
    dndlogr_lib = su_lib / dlogr / (4/3 * np.pi * r**3 * rhop) * np.repeat(rhoa_lib[:,np.newaxis], nbins, axis=1)
    dndlogr_0  = su_0 / dlogr / (4/3 * np.pi * r**3 * rhop) * rhoa_lib[0]
    dadlogr_lib = dndlogr_lib * 4 * np.pi * r**2
    dadlogr_0   = dndlogr_0   * 4 * np.pi * r**2
    dmdlogr_lib = dndlogr_lib * (4/3) * np.pi * r**3 * rhop
    dmdlogr_0   = dndlogr_0   * (4/3) * np.pi * r**3 * rhop

    ncfile = netCDF4.Dataset('parma_box_sulfate.nc4',mode='w')
    ncfile.createDimension('bin', nbins)
    ncfile.createDimension('time', nt+1)
    ncfile.title='PARMA Box Output'
    ncfile.subtitle='Created by parma_box.py'

    ncr = ncfile.createVariable('r', np.float64, ('bin'))
    ncr.long_name = 'Bin (geometric) center radius'
    ncr.units = 'm'
    ncr[:] = r

    nctime = ncfile.createVariable('time', np.float64, ('time'))
    nctime.long_name = 'Simulated time elapsed'
    nctime.units = 's'
    nctime[:] = np.linspace(0,dt*nt,nt+1)

    ncdndlogr = ncfile.createVariable('dNdlogr', np.float64, ('time', 'bin'))
    ncdndlogr.long_name = 'Number concentration as a function of log10(r)'
    ncdndlogr.units = 'm-3'
    ncdndlogr[0,:] = dndlogr_0
    ncdndlogr[1:,:] = dndlogr_lib

    ncdadlogr = ncfile.createVariable('dAdlogr', np.float64, ('time', 'bin'))
    ncdadlogr.long_name = 'Surface area concentration as a function of log10(r)'
    ncdadlogr.units = 'm2 m-3'
    ncdadlogr[0,:] = dadlogr_0
    ncdadlogr[1:,:] = dadlogr_lib

    ncdmdlogr = ncfile.createVariable('dMdlogr', np.float64, ('time', 'bin'))
    ncdmdlogr.long_name = 'Mass concentration as a function of log10(r)'
    ncdmdlogr.units = 'kg m-3'
    ncdmdlogr[0,:] = dmdlogr_0
    ncdmdlogr[1:,:] = dmdlogr_lib

    ncmmr = ncfile.createVariable('MMR', np.float64, ('time', 'bin'))
    ncmmr.long_name = 'Aerosol mass mixing ratio'
    ncmmr.units = 'kg kg-1'
    ncmmr[0,:] = su_0
    ncmmr[1:,:] = su_lib

    ncmass = ncfile.createVariable('M', np.float64, ('time'))
    ncmass.long_name = 'Total aerosol mass concentration'
    ncmass.units = 'kg m-3'
    ncmass[0] = np.sum(dmdlogr_0 * dlogr, axis=0)
    ncmass[1:] = np.sum(dmdlogr_lib * dlogr, axis=1)

    ncarea = ncfile.createVariable('SA', np.float64, ('time'))
    ncarea.long_name = 'Total aerosol surface area concentration'
    ncarea.units = 'm2 m-3'
    ncarea[0] = np.sum(dadlogr_0 * dlogr, axis=0)
    ncarea[1:] = np.sum(dadlogr_lib * dlogr, axis=1)

    ncn = ncfile.createVariable('N', np.float64, ('time'))
    ncn.long_name = 'Total aerosol number concentration'
    ncn.units = '# m-3'
    ncn[0] = np.sum(dndlogr_0 * dlogr, axis=0)
    ncn[1:] = np.sum(dndlogr_lib * dlogr, axis=1)

    nch2so4 = ncfile.createVariable('H2SO4', np.float64, ('time'))
    nch2so4.long_name = 'Sulfuric acid vapor concentration'
    nch2so4.units = 'kg m-3'
    nch2so4[0] = h2so4_0 * rhoa_lib[0]
    nch2so4[1:] = h2so4_lib * rhoa_lib

    ncrh = ncfile.createVariable('RH', np.float64, ('time'))
    ncrh.long_name = 'Relative humidity'
    ncrh.units = '%'
    ncrh[0] = rh_lib[0]
    ncrh[1:] = rh_lib

    ncfile.close()

if __name__ == '__main__':
    ############################################################################
    # Example simulation parameters
    ############################################################################
    # Simulation parameters
    dt    = 900 # timestep (900 seconds (15 minutes))
    nt    = 768 # number of timesteps (8 days)
    nt_carma = 1
    t     = 220 # temperature (kelvin)
    p     = 2400 # pressure (bar)
    h2o   = 125e-6 # H2O mmr (kg/kg)
    h2so4 = 4.7e-10 # H2SO4 mmr (kg/kg)
    constant_h2so4 = False # Keep H2SO4 at initial value?
    su = np.asarray(
         [4.2579942e-28, 1.5634370e-27, 1.1281235e-26, 7.5493089e-24, 4.6745146e-22,
          1.1954208e-20, 1.5612010e-19, 1.3430505e-18, 2.1747372e-17, 1.5991898e-15,
          6.5651301e-14, 1.1735656e-12, 1.3214742e-11, 7.1381352e-11, 2.0359722e-10,
          3.1819103e-10, 1.3236401e-10, 1.2854880e-11, 3.6080232e-13, 2.5853792e-15,
          4.4601091e-18, 2.1512019e-21, 2.8390088e-25, 5.5137755e-29]
         ) # sulfate aerosol mmr (kg/kg)

    # CARMA bin structure
    nbins = 24
    rmrat = 3.75125201
    rmin = 2.6686863e-8 # cm
    rhop = 1.923 # g cm-3
    run_box(dt, nt, nt_carma, nbins, rmrat, rmin, rhop, h2o, h2so4, \
                su, t, p, constant_h2so4 = constant_h2so4)
