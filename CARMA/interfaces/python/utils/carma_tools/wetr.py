"""Module for calculating wet radius of CARMA bins in model output
Translated from Dr. Peter Colarco's IDL program
Method originally from Tabazadeh

Author: Parker A. Case
Changelog: Created 2/9/2019
"""
import numpy as np

def wtpct_sulf(relhum, T=220):
    rhopdry = 1.923
    if T is None:
        T = 220 # 220 Kelvin as default temperature value

    activ = relhum
    atab1 = np.zeros(activ.shape)
    btab1 = np.zeros(activ.shape)
    ctab1 = np.zeros(activ.shape)
    dtab1 = np.zeros(activ.shape)
    atab2 = np.zeros(activ.shape)
    btab2 = np.zeros(activ.shape)
    ctab2 = np.zeros(activ.shape)
    dtab2 = np.zeros(activ.shape)

    idx_lt05 = activ < 0.05
    idx_gt05lt85 = (activ >= 0.05) & (activ <=0.85)
    idx_gt85 = activ > 0.85

    activ[activ < 1.e-6] = 1.e-6
    activ[activ > 1.] = 1.

    atab1[idx_lt05] = 12.37208932
    btab1[idx_lt05] = -0.16125516114
    ctab1[idx_lt05] = -30.490657554
    dtab1[idx_lt05] = -2.1133114241
    atab2[idx_lt05] = 13.455394705
    btab2[idx_lt05] = -0.1921312255
    ctab2[idx_lt05] = -34.285174607
    dtab2[idx_lt05] = -1.7620073078

    atab1[idx_gt05lt85] = 11.820654354
    btab1[idx_gt05lt85] = -0.20786404244
    ctab1[idx_gt05lt85] = -4.807306373
    dtab1[idx_gt05lt85] = -5.1727540348
    atab2[idx_gt05lt85] = 12.891938068
    btab2[idx_gt05lt85] = -0.23233847708
    ctab2[idx_gt05lt85] = -6.4261237757
    dtab2[idx_gt05lt85] = -4.9005471319

    atab1[idx_gt85] = -180.06541028
    btab1[idx_gt85] = -0.38601102592
    ctab1[idx_gt85] = -93.317846778
    dtab1[idx_gt85] = 273.88132245
    atab2[idx_gt85] = -176.95814097
    btab2[idx_gt85] = -0.36257048154
    ctab2[idx_gt85] = -90.469744201
    dtab2[idx_gt85] = 267.45509988

    contl = atab1*(activ**btab1)+ctab1*activ+dtab1
    conth = atab2*(activ**btab2)+ctab2*activ+dtab2

    contt = contl + (conth-contl) * ((T - 190.)/70.)
    conwtp = (contt*98.) + 1000.

    wtpct_tabaz = (100.*contt*98.)/conwtp
    # Restrict values to between 1 and 100
    wtpct_tabaz[wtpct_tabaz > 100.] = 100.
    wtpct_tabaz[wtpct_tabaz < 1.] = 1.

    return wtpct_tabaz

def dens(relhum, T=220):
    wtp = wtpct_sulf(relhum, T=T)

    dnwtp =[ 0., 1., 5., 10., 20., 25., 30., 35., 40., \
       41., 45., 50., 53., 55., 56., 60., 65., 66., 70., \
       72., 73., 74., 75., 76., 78., 79., 80., 81., 82., \
       83., 84., 85., 86., 87., 88., 89., 90., 91., 92., \
       93., 94., 95., 96., 97., 98., 100. ]

    dnc0 =[ 1., 1.13185, 1.17171, 1.22164, 1.3219, 1.37209,       \
       1.42185, 1.4705, 1.51767, 1.52731, 1.56584, 1.61834, 1.65191, \
       1.6752, 1.68708, 1.7356, 1.7997, 1.81271, 1.86696, 1.89491,   \
       1.9092, 1.92395, 1.93904, 1.95438, 1.98574, 2.00151, 2.01703, \
       2.03234, 2.04716, 2.06082, 2.07363, 2.08461, 2.09386, 2.10143,\
       2.10764, 2.11283, 2.11671, 2.11938, 2.12125, 2.1219, 2.12723, \
       2.12654, 2.12621, 2.12561, 2.12494, 2.12093 ]

    dnc1 = [ 0.,  -0.000435022, -0.000479481, -0.000531558, -0.000622448,\
       -0.000660866, -0.000693492, -0.000718251, -0.000732869, -0.000735755, \
       -0.000744294, -0.000761493, -0.000774238, -0.00078392, -0.000788939,  \
       -0.00080946, -0.000839848, -0.000845825, -0.000874337, -0.000890074,  \
       -0.00089873, -0.000908778, -0.000920012, -0.000932184, -0.000959514,  \
       -0.000974043, -0.000988264, -0.00100258, -0.00101634, -0.00102762,    \
       -0.00103757, -0.00104337, -0.00104563, -0.00104458, -0.00104144,      \
       -0.00103719, -0.00103089, -0.00102262, -0.00101355, -0.00100249,      \
       -0.00100934, -0.000998299, -0.000990961, -0.000985845, -0.000984529,  \
       -0.000989315 ]

    idx = np.zeros((wtp.shape))
    den1 = np.zeros((wtp.shape))
    den2 = np.zeros((wtp.shape))
    frac = np.zeros((wtp.shape))
    for i in range(len(dnwtp)):
        idx[wtp > dnwtp[i]] = i
        den1[wtp > dnwtp[i]] = dnc0[i-1] + dnc1[i-1] * T
        den2[wtp > dnwtp[i]] = dnc0[i] + dnc1[i] * T
        frac[wtp > dnwtp[i]] = (dnwtp[i] - wtp[wtp > dnwtp[i]]) / (dnwtp[i] - dnwtp[i-1])
    dens = den1*frac + den2*(1.-frac)
    return dens

def grow_v75(relhum, rdry, T=273.16):
    # rdry should be in cm here
    rhopdry = 1.923 # g cm-3
    mw_h2so4 = 98. # molecular weight of H2SO4 g mol-1
    rgas = 8.31447e7 # erg mol-1 K-1

    # First calculate the mass concentration of water given rh and T
    # saturation vapor pressure (Curry & Webster, 4.31)
    llv = 2.501e6 # J kg-1 from C&W
    rv = 461. # J K-1 kg-1 from C&W
    ttr = 273.16 # K Triple Point
    estr = 611. # Pa
    es = estr*np.exp(llv/rv*(1./ttr - 1/T))
    # Ideal gas law, pv = nkT, rearrange to n/v = p/(kT)
    k = 1.38e-23 # Boltzman constant mks
    n_v = relhum * es / (k * T) # number m-3
    navogad = 6.022e23 # mole-1
    mw_h2o = 0.018 # kg mole-1
    # Finally mass concentration, and note units to g cm-3
    h2o_mass = n_v / navogad * mw_h2o * 1000. / 1.e6

    # Kelvin effect adjustment
    wtpkelv = 80.
    den1 = 2.00151 - 0.000974043 * T  # density at 79 wt %
    den2 = 2.01703 - 0.000988264 * T  # density at 80 wt %
    drho_dwt = den2-den1                   # change in density for change in 1 wt %

    sig1 = 79.3556 - 0.0267212 * T   # surface tension at 79.432 wt %
    sig2 = 75.608  - 0.0269204 * T   # surface tension at 85.9195 wt %
    dsigma_dwt = (sig2-sig1) / (85.9195 - 79.432) # change in density for change in 1 wt %
    sigkelv = sig1 + dsigma_dwt * (80.0 - 79.432)

    rwet = rdry * (100. * rhopdry / wtpkelv / den2)**(1. / 3.)

    rkelvinH2O_b = 1. + wtpkelv * drho_dwt / den2 - 3. * wtpkelv \
            * dsigma_dwt / (2.*sigkelv)

    rkelvinH2O_a = 2. * mw_h2so4 * sigkelv / (den1 * rgas * T * rwet)

    rkelvinH2O = np.exp(rkelvinH2O_a*rkelvinH2O_b)

    h2o_kelv = h2o_mass / rkelvinH2O

    # wtpct just wants a relative humidity, but it is recalculated here
    # based on the calculated terms from above
    k_cgs = 1.3807e-16 # cm2 g s-2 K-1
    mw_h2o_cgs = 18. # g mol-1
    relhum_ = h2o_kelv*navogad/mw_h2o_cgs*k_cgs*T / (es*10.) # es now dyne cm-2
    rhopwet = dens(relhum_,T=T)
    rwet    = rdry * (100. * rhopdry / wtpct_sulf(relhum_) / rhopwet)**(1. / 3.)

    return rwet
