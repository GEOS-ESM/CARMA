"""Primary interface for running PARMA column with dust + sulfate setup

This is an interface for running PARMA with a pure sulfate group and a mixed
group with two elements: sulfate and dust. The run_column() function can be
used to run these simulations. The main method at the bottom of this file
gives an example of using the run_column() function.

For more information on how to run PARMA, see the README in the interfaces
directory.

Author: Parker Case
Version: v1.0 (2024/08/20) Broke out simulation and netCDF writing into
                            function.
         v0.4 (2024/08/19) Implemented kok size distribution and checks
                            for negative
         v0.3 (2024/08/16) Implemented variable density for mixed group
         v0.2 (2024/08/14) Readability, RH output, air density output,
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
from carma_tools.dust_tools import kok_size_distribution

def run_column(nz, dt, nt, nt_carma, nbins, rmrat_su, rmin_su, rmrat_mx, \
                rmin_mx, rhop_su, rhop_du, h2o, h2so4, su, mxsu, mxdu, \
                constant_h2so4 = False, exp_name = '0'):
    """Runs a CARMA column with the specified parameters

    Arguments:
    nz -- number of layers to simulate in the column
    dt -- length of timestep (seconds)
    nt -- number of PARMA timesteps to simulate
    nt_carma -- number of CARMA timesteps per PARMA timestep
    nbins -- number of CARMA bins
    rmrat_su -- ratio between bins in the pure sulfate group
    rmin_su -- minimum radius in the pure sulfate group (cm)
    rmrat_mx -- ratio between bins in the mixed group
    rmin_mx -- minimum radius in the mixed group (cm)
    rhop_su -- density of sulfate (g cm-3)
    rhop_du -- density of dust (g cm-3)
    h2o -- array (nz,) of water vapor mass mixing ratios (g g-1)
    h2so4 -- array (nz,) of h2so4 vapor mass mixing ratios (g g-1)
    su -- array (nz, nbins) of sulfate mass mixing ratios (g g-1)
    mxsu -- array (nz, nbins) of mixed group sulfate mass mixing ratios (g g-1)
    mxdu -- array (nz, nbins) of mixed group dust mass mixing ratios (g g-1)
    constant_h2so4 -- boolean indicating whether CARMA should deplete h2so4
    """

    # Prep the rest of the variables for CARMA
    alt = np.linspace(0,nz*1000+.01,nz)
    p, t = get_standard_atmosphere_1d(alt*1000)
    su_0 = su.copy()
    mxsu_0 = mxsu.copy()
    mxdu_0 = mxdu.copy()
    h2so4_0 = h2so4.copy()
    h2o_0 = h2o.copy()
    mxsu_0 = mxsu_0 - mxdu_0


    ############################################################################
    # Run the model!
    ############################################################################
    su_lib = np.zeros((nt, nz, nbins))
    mxsu_lib = np.zeros((nt, nz, nbins))
    mxdu_lib = np.zeros((nt, nz, nbins))
    h2so4_lib = np.zeros((nt, nz))
    rh_lib = np.zeros((nt, nz))
    rhoa_lib = np.zeros((nt, nz))
    time_0 = perf_counter()
    for i in range(nt):
        su_out, mxsu_out, mxdu_out, t_out, p_out, h2so4_out, h2o_out, \
                rhoa_out, rh_out = carma_column_dust(rmrat_su, rmrat_mx, \
                rmin_su, rmin_mx, rhop_su, rhop_du, t, p, h2so4, h2o, su, \
                mxsu, mxdu, dt, nt_carma, constant_h2so4, nbins, nz)
        su_lib[i,:,:] = su_out
        mxsu_lib[i,:,:] = mxsu_out - mxdu_out
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

    ############################################################################
    # Check for negative values produced by post processing subtraction
    ############################################################################
    mxsu_lib[mxsu_lib < 0.] = 1e-28

    ############################################################################
    # Convert to mKs units for analysis:
    ############################################################################
    rhop_su = rhop_su * 1000 # kg m-3
    rhop_du = rhop_du * 1000 # kg m-3
    rmin_su = rmin_su * 1e-2 # m
    rmin_mx = rmin_mx * 1e-2 # m
    rhoa_lib = rhoa_lib * 1000 # kg m-3

    ############################################################################
    # Build CARMA bin structure for analysis
    ############################################################################
    rhomx_lib = (mxsu_lib * rhop_su + mxdu_lib * rhop_du) / \
            (mxsu_lib + mxdu_lib)
    rmass_su, rmassup_su, r_su, rup_su, dr_su, rlow_su, masspart_su = \
        carmabins.carmabins(nbins, rmrat_su, rmin_su, rhop_su)
    rmass_mx = np.zeros((nt,nz,nbins))
    rmassup_mx = np.zeros((nt,nz,nbins))
    r_mx = np.zeros((nt,nz,nbins))
    rup_mx = np.zeros((nt,nz,nbins))
    dr_mx = np.zeros((nt,nz,nbins))
    rlow_mx = np.zeros((nt,nz,nbins))
    masspart_mx = np.zeros((nt,nz,nbins))
    for i in range(nt):
        for j in range(nz):
            rmass_mx[i,j,:], rmassup_mx[i,j,:], r_mx[i,j,:], rup_mx[i,j,:], \
                    dr_mx[i,j,:], rlow_mx[i,j,:], masspart_mx[i,j,:] = \
                    carmabins.carmabins_mx(nbins, rmrat_mx, rmin_mx, \
                    rhomx_lib[i,j,:])

    dlogr_su = np.log10(rup_su) - np.log10(rlow_su)
    dndlogr_su_lib = su_lib / dlogr_su / (4/3 * np.pi * r_su**3 * rhop_su) * \
            np.repeat(rhoa_lib[:,:,np.newaxis], nbins, axis=2)
    dndlogr_su_0   = su_0   / dlogr_su / (4/3 * np.pi * r_su**3 * rhop_su) * \
            np.repeat(rhoa_lib[0,:,np.newaxis], nbins, axis=1)
    dadlogr_su_lib = dndlogr_su_lib * 4 * np.pi * r_su**2
    dadlogr_su_0   = dndlogr_su_0   * 4 * np.pi * r_su**2
    dmdlogr_su_lib = dndlogr_su_lib * (4/3) * np.pi * r_su**3 * rhop_su
    dmdlogr_su_0   = dndlogr_su_0   * (4/3) * np.pi * r_su**3 * rhop_su

    dlogr_mx = np.log10(rup_mx) - np.log10(rlow_mx)

    dndlogr_mxsu_lib = mxsu_lib / dlogr_mx / (4/3 * np.pi * r_mx**3 * rhomx_lib) * \
            np.repeat(rhoa_lib[:,:,np.newaxis], nbins, axis=2)
    dndlogr_mxsu_0   = mxsu_0   / dlogr_mx[0,:,:] / (4/3 * np.pi * r_mx[0,:,:]**3 \
            * rhomx_lib[0,:,:]) * np.repeat(rhoa_lib[0,:,np.newaxis], nbins, axis=1)
    dadlogr_mxsu_lib = dndlogr_mxsu_lib * 4 * np.pi * r_mx**2
    dadlogr_mxsu_0   = dndlogr_mxsu_0   * 4 * np.pi * r_mx[0,:,:]**2
    dmdlogr_mxsu_lib = dndlogr_mxsu_lib * (4/3) * np.pi * r_mx**3 * rhomx_lib
    dmdlogr_mxsu_0   = dndlogr_mxsu_0   * (4/3) * np.pi * r_mx[0,:,:]**3 * \
            rhomx_lib[0,:,:]

    dndlogr_mxdu_lib = mxdu_lib / dlogr_mx / (4/3 * np.pi * r_mx**3 * rhomx_lib) * \
            np.repeat(rhoa_lib[:,:,np.newaxis], nbins, axis=2)
    dndlogr_mxdu_0   = mxdu_0   / dlogr_mx[0,:,:] / (4/3 * np.pi * r_mx[0,:,:]**3 \
            * rhomx_lib[0,:,:]) * np.repeat(rhoa_lib[0,:,np.newaxis], nbins, axis=1)
    dadlogr_mxdu_lib = dndlogr_mxdu_lib * 4 * np.pi * r_mx**2
    dadlogr_mxdu_0   = dndlogr_mxdu_0   * 4 * np.pi * r_mx[0,:,:]**2
    dmdlogr_mxdu_lib = dndlogr_mxdu_lib * (4/3) * np.pi * r_mx**3 * rhomx_lib
    dmdlogr_mxdu_0   = dndlogr_mxdu_0   * (4/3) * np.pi * r_mx[0,:,:]**3 * \
            rhomx_lib[0,:,:]

    ############################################################################
    # Create .nc4 file
    ############################################################################
    ncfile = netCDF4.Dataset('parma_column_dust_' + exp_name + '.nc4',mode='w')
    ncfile.createDimension('layer', nz)
    ncfile.createDimension('bin', nbins)
    ncfile.createDimension('time', nt+1)
    ncfile.title='PARMA Column Output'
    ncfile.subtitle='Created by parma_column_dust.py'

    ncr = ncfile.createVariable('r_su', np.float64, ('bin'))
    ncr.long_name = 'Sulfate bin (geometric) center radius'
    ncr.units = 'm'
    ncr[:] = r_su

    ncrmx = ncfile.createVariable('r_mx', np.float64, ('time', 'layer', 'bin'))
    ncrmx.long_name = 'Mixed bin (geometric) center radius'
    ncrmx.units = 'm'
    ncrmx[0,:,:] = r_mx[0,:,:]
    ncrmx[1:,:,:] = r_mx

    ncalt = ncfile.createVariable('alt', np.float64, ('layer'))
    ncalt.long_name = 'Altitude'
    ncalt.units = 'm'
    ncalt[:] = np.linspace(0,(nz-1)*1000,nz)+500

    nctime = ncfile.createVariable('time', np.float64, ('time'))
    nctime.long_name = 'Simulated time elapsed'
    nctime.units = 's'
    nctime[:] = np.linspace(0,dt*nt,nt+1)

    ncdndlogr = ncfile.createVariable('dNdlogr_su', np.float64, \
            ('time', 'layer', 'bin'))
    ncdndlogr.long_name = 'Number concentration as a function of log10(r)'
    ncdndlogr.units = 'm-3'
    ncdndlogr[0,:,:] = dndlogr_su_0
    ncdndlogr[1:,:,:] = dndlogr_su_lib

    ncdadlogr = ncfile.createVariable('dAdlogr_su', np.float64, \
            ('time', 'layer', 'bin'))
    ncdadlogr.long_name = 'Surface area concentration as a function of log10(r)'
    ncdadlogr.units = 'm2 m-3'
    ncdadlogr[0,:,:] = dadlogr_su_0
    ncdadlogr[1:,:,:] = dadlogr_su_lib

    ncdmdlogr = ncfile.createVariable('dMdlogr_su', np.float64, \
            ('time', 'layer', 'bin'))
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

    ncdndlogrmxsu = ncfile.createVariable('dNdlogr_mxsu', np.float64, \
            ('time', 'layer', 'bin'))
    ncdndlogrmxsu.long_name = 'Number concentration as a function of log10(r)'
    ncdndlogrmxsu.units = 'm-3'
    ncdndlogrmxsu[0,:,:] = dndlogr_mxsu_0
    ncdndlogrmxsu[1:,:,:] = dndlogr_mxsu_lib

    ncdadlogrmxsu = ncfile.createVariable('dAdlogr_mxsu', np.float64, \
            ('time', 'layer', 'bin'))
    ncdadlogrmxsu.long_name = 'Surface area concentration as a function of log10(r)'
    ncdadlogrmxsu.units = 'm2 m-3'
    ncdadlogrmxsu[0,:,:] = dadlogr_mxsu_0
    ncdadlogrmxsu[1:,:,:] = dadlogr_mxsu_lib

    ncdmdlogrmxsu = ncfile.createVariable('dMdlogr_mxsu', np.float64, \
            ('time', 'layer', 'bin'))
    ncdmdlogrmxsu.long_name = 'Mass concentration as a function of log10(r)'
    ncdmdlogrmxsu.units = 'kg m-3'
    ncdmdlogrmxsu[0,:,:] = dmdlogr_mxsu_0
    ncdmdlogrmxsu[1:,:,:] = dmdlogr_mxsu_lib

    ncmmrmxsu = ncfile.createVariable('MMR_mxsu', np.float64, \
            ('time', 'layer', 'bin'))
    ncmmrmxsu.long_name = 'Aerosol mass mixing ratio'
    ncmmrmxsu.units = 'kg kg-1'
    ncmmrmxsu[0,:,:] = mxsu_0
    ncmmrmxsu[1:,:,:] = mxsu_lib

    ncmassmxsu = ncfile.createVariable('M_mxsu', np.float64, ('time', 'layer'))
    ncmassmxsu.long_name = 'Total aerosol mass concentration'
    ncmassmxsu.units = 'kg m-3'
    ncmassmxsu[0,:] = np.sum(dmdlogr_mxsu_0 * dlogr_mx[0,:,:], axis=1)
    ncmassmxsu[1:,:] = np.sum(dmdlogr_mxsu_lib * dlogr_mx, axis=2)

    ncareamxsu = ncfile.createVariable('SA_mxsu', np.float64, ('time', 'layer'))
    ncareamxsu.long_name = 'Total aerosol surface area concentration'
    ncareamxsu.units = 'm2 m-3'
    ncareamxsu[0,:] = np.sum(dadlogr_mxsu_0 * dlogr_mx[0,:,:], axis=1)
    ncareamxsu[1:,:] = np.sum(dadlogr_mxsu_lib * dlogr_mx, axis=2)

    ncnmxsu = ncfile.createVariable('N_mxsu', np.float64, ('time', 'layer'))
    ncnmxsu.long_name = 'Total aerosol number concentration'
    ncnmxsu.units = '# m-3'
    ncnmxsu[0,:] = np.sum(dndlogr_mxsu_0 * dlogr_mx[0,:,:], axis=1)
    ncnmxsu[1:,:] = np.sum(dndlogr_mxsu_lib * dlogr_mx, axis=2)

    ncdndlogrmxdu = ncfile.createVariable('dNdlogr_mxdu', np.float64, \
            ('time', 'layer', 'bin'))
    ncdndlogrmxdu.long_name = 'Number concentration as a function of log10(r)'
    ncdndlogrmxdu.units = 'm-3'
    ncdndlogrmxdu[0,:,:] = dndlogr_mxdu_0
    ncdndlogrmxdu[1:,:,:] = dndlogr_mxdu_lib

    ncdadlogrmxdu = ncfile.createVariable('dAdlogr_mxdu', np.float64, \
            ('time', 'layer', 'bin'))
    ncdadlogrmxdu.long_name = 'Surface area concentration as a function of log10(r)'
    ncdadlogrmxdu.units = 'm2 m-3'
    ncdadlogrmxdu[0,:,:] = dadlogr_mxdu_0
    ncdadlogrmxdu[1:,:,:] = dadlogr_mxdu_lib

    ncdmdlogrmxdu = ncfile.createVariable('dMdlogr_mxdu', np.float64, \
            ('time', 'layer', 'bin'))
    ncdmdlogrmxdu.long_name = 'Mass concentration as a function of log10(r)'
    ncdmdlogrmxdu.units = 'kg m-3'
    ncdmdlogrmxdu[0,:,:] = dmdlogr_mxdu_0
    ncdmdlogrmxdu[1:,:,:] = dmdlogr_mxdu_lib

    ncmmrmxdu = ncfile.createVariable('MMR_mxdu', np.float64, \
            ('time', 'layer', 'bin'))
    ncmmrmxdu.long_name = 'Aerosol mass mixing ratio'
    ncmmrmxdu.units = 'kg kg-1'
    ncmmrmxdu[0,:,:] = mxdu_0
    ncmmrmxdu[1:,:,:] = mxdu_lib

    ncmassmxdu = ncfile.createVariable('M_mxdu', np.float64, ('time', 'layer'))
    ncmassmxdu.long_name = 'Total aerosol mass concentration'
    ncmassmxdu.units = 'kg m-3'
    ncmassmxdu[0,:] = np.sum(dmdlogr_mxdu_0 * dlogr_mx[0,:,:], axis=1)
    ncmassmxdu[1:,:] = np.sum(dmdlogr_mxdu_lib * dlogr_mx, axis=2)

    ncareamxdu = ncfile.createVariable('SA_mxdu', np.float64, ('time', 'layer'))
    ncareamxdu.long_name = 'Total aerosol surface area concentration'
    ncareamxdu.units = 'm2 m-3'
    ncareamxdu[0,:] = np.sum(dadlogr_mxdu_0 * dlogr_mx[0,:,:], axis=1)
    ncareamxdu[1:,:] = np.sum(dadlogr_mxdu_lib * dlogr_mx, axis=2)

    ncnmxdu = ncfile.createVariable('N_mxdu', np.float64, ('time', 'layer'))
    ncnmxdu.long_name = 'Total aerosol number concentration'
    ncnmxdu.units = '# m-3'
    ncnmxdu[0,:] = np.sum(dndlogr_mxdu_0 * dlogr_mx[0,:,:], axis=1)
    ncnmxdu[1:,:] = np.sum(dndlogr_mxdu_lib * dlogr_mx, axis=2)

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

    rmass_su, rmassup_su, r_su, rup_su, dr_su, rlow_su, masspart_su = \
        carmabins.carmabins(nbins, rmrat_su, rmin_su, rhop_su)
    rmass_mx, rmassup_mx, r_mx, rup_mx, dr_mx, rlow_mx, masspart_mx = \
        carmabins.carmabins(nbins, rmrat_mx, rmin_mx, rhop_du)
    kok_dm = kok_size_distribution(nbins, r_mx, rlow_mx, rup_mx, rhop_du, rhop_du)
    kok_dm = kok_dm / np.sum(kok_dm)

    # Background atmospheric parameters
    constant_h2so4 = False # Keep H2SO4 at initial value?
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
    mxdu = np.zeros((nz,nbins)) + 1e-10 * kok_dm # Mixed group dust aerosol mmr (kg/kg)
    mxsu = np.zeros((nz,nbins)) + mxdu + 1e-20 # Mixed group sulfate aerosol mmr (kg/kg)

    run_column(nz, dt, nt, nt_carma, nbins, rmrat_su, rmin_su, rmrat_mx, \
                rmin_mx, rhop_su, rhop_du, h2o, h2so4, su, mxsu, mxdu, \
                constant_h2so4 = constant_h2so4)

