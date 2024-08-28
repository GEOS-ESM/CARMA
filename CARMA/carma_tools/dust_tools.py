"""Utilities related to CARMA dust

Useful functions for initializing and analyzing dust in CARMA

Author: Parker Case
Version: v0.1 (2024/08/19) First commit
"""
import numpy as np
from scipy.special import erf
import carma_tools.carmabins as carmabins

def kok_size_distribution(nbin, r, rlow, rup, rhod, rhog):
    """Generate a mass distribution following brittle fragmentation theory

    See Kok PNAS (2011)

    Arguments:
    nbin -- number of bins
    r -- center radius of each bin (cm)
    rlow -- lower limit radius of each bin (cm)
    rup -- upper limit radius of each bin (cm)
    rhod -- density of dust (g cm-3)
    rhog -- density of mixed group (g cm-3)
    """
    ds = 3.4e-6
    sigma = 3.
    lam = 12.e-6
    cv = 12.62e-6
    cn = 0.9539e-6

    nbin_ = 1000
    rho_dust = rhod*1000
    rho_grp  = rhog*1000
    r    = r * 1e-2
    rlow = rlow * 1e-2
    rup  = rup * 1e-2

    dm = np.zeros((nbin,))
    for ibin in range(nbin):
        rmrat_ = (rup[ibin]/rlow[ibin])**(3./nbin_)
        rmin_  = rlow[ibin]*((1.+rmrat_)/2.)**(1./3.)

        rmass_, rmassup_, r_, rup_, dr_, rlow_, masspart_ = \
                carmabins.carmabins(nbin, rmrat_, rmin_, rho_dust)

        r_ = r_*(rho_grp/rho_dust)**(1./3.)
        dr_ = dr_*(rho_grp/rho_dust)**(1./3.)

        dm_ = (dr_/r_) * (2.*r_/cv) * \
                (1. + erf(np.log(2.* r_/ ds) / np.sqrt(2.) / np.log(sigma))) * \
                np.exp(-1*(2.*r_/lam)**3)
        dm[ibin] = np.sum(dm_)

    return dm


