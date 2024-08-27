"""Module to create CARMA like bins for size distribution analysis
Translated from Dr. Peter Colarco's IDL program

Author: Parker A. Case
Changelog: Created 6/21/2017
Added get_variables and dndr calculation 1/31/2019
"""
import numpy as np

def get_variables(prefix, nbin, fill):
    """Create list of GEOS-5 CARMA variable names

    Arguments:
    prefix -- the name of the group ex. su
    nbin -- number of bins in the group
    fill -- number of digits in variable names (3 usually)
    """
    binvarnumbers = [str(i).zfill(fill) for i in range(1, nbin+1)]
    variables = [prefix + binvarnumber for binvarnumber in binvarnumbers]
    return variables

def carmabins(nbin, rmrat, rmin, rhop):
    """Create CARMA like size distribution bins

    Arguments:
    nbin -- number of bins
    rmrat -- mass ratio
    rmin -- minimum radius
    rhop -- density of particles
    """
    cpi = 4./3. * np.pi
    # Convert minimum radius to mass
    rmassmin = cpi*rhop*rmin**3.
    # Volume to radius factor
    vrfact = ((3./2. / np.pi / (rmrat+1))**(1./3.))*(rmrat**(1./3.) - 1.)

    rmass = np.zeros(nbin)
    rmassup = np.zeros(nbin)
    r = np.zeros(nbin)
    rup = np.zeros(nbin)
    dr = np.zeros(nbin)
    rlow = np.zeros(nbin)

    for ibin in range(0, nbin):
        rmass[ibin]   = rmassmin*rmrat**ibin # Bin median mass
        rmassup[ibin] = 2.*rmrat/(rmrat+1.)*rmass[ibin] # Bin maximum radius
        r[ibin]       = (rmass[ibin]/rhop/cpi)**(1./3.) # Bin median radius
        rup[ibin]     = (rmassup[ibin]/rhop/cpi)**(1./3.) # Bin maximum radius
        dr[ibin]      = vrfact*(rmass[ibin]/rhop)**(1./3.) # Bin width in radius
        rlow[ibin]    = rup[ibin] - dr[ibin] # Bin minimum radius

    masspart = 4./3. * np.pi * r**3. * rhop # Mass of median radius

    return rmass, rmassup, r, rup, dr, rlow, masspart

def carmabins_mx(nbin, rmrat, rmin, rhop):
    """Create CARMA like size distribution bins

    Arguments:
    nbin -- number of bins
    rmrat -- mass ratio
    rmin -- minimum radius
    rhop -- density of particles
    """
    cpi = 4./3. * np.pi
    # Convert minimum radius to mass
    rmassmin = cpi*rhop[0]*rmin**3.
    # Volume to radius factor
    vrfact = ((3./2. / np.pi / (rmrat+1))**(1./3.))*(rmrat**(1./3.) - 1.)

    rmass = np.zeros(nbin)
    rmassup = np.zeros(nbin)
    r = np.zeros(nbin)
    rup = np.zeros(nbin)
    dr = np.zeros(nbin)
    rlow = np.zeros(nbin)

    for ibin in range(0, nbin):
        rmass[ibin]   = rmassmin*rmrat**ibin # Bin median mass
        rmassup[ibin] = 2.*rmrat/(rmrat+1.)*rmass[ibin] # Bin maximum radius
        r[ibin]       = (rmass[ibin]/rhop[ibin]/cpi)**(1./3.) # Bin median radius
        rup[ibin]     = (rmassup[ibin]/rhop[ibin]/cpi)**(1./3.) # Bin maximum radius
        dr[ibin]      = vrfact*(rmass[ibin]/rhop[ibin])**(1./3.) # Bin width in radius
        rlow[ibin]    = rup[ibin] - dr[ibin] # Bin minimum radius

    masspart = 4./3. * np.pi * r**3. * rhop[ibin] # Mass of median radius

    return rmass, rmassup, r, rup, dr, rlow, masspart

def dndr(mmr, dr, r, rhop, airdensity):
    """Calculate dndr from mass mixing ratio that comes from CARMA output

    Arguments:
    mmr -- 3 dimensional CARMA data (bin, altitude, latitude)
    dr -- spacing between bins from carmabins()
    r -- radius of bin centers from carmabins()
    rhop -- density of particles
    airdensity -- 2 dimensional GEOS-5 airdensity
    """
    dndr = np.zeros(mmr.shape)
    for j in range(mmr.shape[0]):
        dndr[j,:,] = mmr[j,:] / dr[j] / (4/3 * np.pi *
                r[j]**3 * rhop) * airdensity[:]
    return dndr

def dndr_0d(mmr, dr, r, rhop, airdensity):
    """Calculate dndr from mass mixing ratio that comes from CARMA output,
        for a single ASD. Use dndr for higher dimension data.

    Arguments:
    mmr -- 3 dimensional CARMA data (bin)
    dr -- spacing between bins from carmabins()
    r -- radius of bin centers from carmabins()
    rhop -- density of particles
    airdensity -- GEOS airdensity
    """
    dndr = mmr / dr / (4/3 * np.pi * r**3 * rhop) * airdensity
    return dndr

def cdf(dndr, rwet, dr, bin_mins):
    """Calculate cumulitive distribution function with custom bins
    Developed for comparing with Deshler balloon data

    Arguments:
    dndr -- dndr as calculated by dndr()
    rwet -- wet center radii of each bin
    bin_mins -- radii to sum above for each bin
    """
    cdf = np.zeros(len(bin_mins))
    dn = dndr * dr
    if len(dr) != len(bin_mins):
        for min, i in zip(bin_mins, range(len(bin_mins))):
            cdf[i] = np.nansum(dn[rwet > min])
    else:
        for i in range(len(dr)):
            cdf[i] = np.nansum(dn[i:])
    return np.asarray(cdf)

def lognormal(r, N, mu, sigma):
    """Create a lognormal distribution, sampled at CARMA bins. Useful for
    initializing CARMA or comparing to data.

    Arguments:
    r -- center radius bins, like those from carmabins()
    N -- scaling variable, total number concentration
    mu -- median radius, same units as r
    sigma -- distribution width parameter
    """
    return N/(np.sqrt(2*np.pi)*sigma*r)*np.exp(-((np.log(r)-np.log(mu))**2)/(2*sigma**2))

