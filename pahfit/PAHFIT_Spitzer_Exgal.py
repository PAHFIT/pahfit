from __future__ import (absolute_import, print_function, division)

import numpy as np

import matplotlib.pyplot as plt
import matplotlib

from astropy.io import fits
from astropy.modeling.fitting import LevMarLSQFitter

from PAHFITBase import PAHFITBase


if __name__ == '__main__':
    # define the instrument resolution
    inst_res = 100.

    # the blackbodies for the continuum
    bb_temps = [3000., 2000., 1000., 800., 600., 450., 300., 200., 135.,
                90., 65., 50., 40., 35.]
    bb_amps = np.full((len(bb_temps)), 1.0e-10)
    bb_amps_limits = [(0.0, None) for k in range(len(bb_temps))]
    bb_info = (bb_amps, bb_temps, bb_amps_limits)

    # the dust features (mainly the PAH/aromatic features)
    df_cwave = np.array([5.27, 5.70, 6.22, 6.69,
                        7.42, 7.60, 7.85, 8.33,
                        8.61, 10.68, 11.23, 11.33,
                        11.99, 12.62, 12.69, 13.48,
                        14.04, 14.19, 15.9, 16.45,
                        17.04, 17.375, 17.87, 18.92,
                        33.1])
    df_frac_fwhm = np.array([0.034, 0.035, 0.030, 0.07,
                             0.126, 0.044, 0.053, 0.05,
                             0.039, 0.02, 0.012, 0.032,
                             0.045, 0.042, 0.013, 0.04,
                             0.016, 0.025, 0.02, 0.014,
                             0.065, 0.012, 0.016, 0.019,
                             0.05])
    df_fwhm = df_frac_fwhm*df_cwave
    n_df = len(df_cwave)
    df_amps = np.full((n_df), 100.)
    df_amps_limits = [(0.0, None) for k in range(n_df)]
    df_cwave_limits = [(cwave - 0.1, cwave + 0.1) for cwave in df_cwave]
    df_fwhm_limits = [(cfwhm*0.9, cfwhm*1.1) for cfwhm in df_fwhm]
    dust_features = (df_amps, df_cwave, df_fwhm,
                     df_amps_limits, df_cwave_limits, df_fwhm_limits)

    # define the H2 features
    h2_cwaves = np.array([5.5115, 6.1088, 6.9091, 8.0258,
                          9.6649, 12.2785, 17.0346, 28.2207])
    h2_fwhm = h2_cwaves/inst_res
    h2_names = np.array(["H2 S(7)", "H2 S(6)", "H2 S(5)", "H2 S(4)",
                         "H2 S(3)", "H2 S(2)", "H2 S(1)", "H2 S(0)"])
    n_h2 = len(h2_cwaves)
    h2_amps = np.full((n_h2), 100.)
    h2_amps_limits = [(0.0, None) for k in range(n_h2)]
    h2_cwaves_limits = [(cwave - 0.05, cwave + 0.05) for cwave in h2_cwaves]
    h2_fwhm_limits = [(cfwhm*0.5, cfwhm*1.5) for cfwhm in h2_fwhm]
    h2_features = (h2_amps, h2_cwaves, h2_fwhm,
                   h2_amps_limits, h2_cwaves_limits, h2_fwhm_limits,
                   h2_names)

    # define the ionic features
    ion_cwaves = np.array([6.985274, 8.99138, 10.5105, 12.813,
                           15.555, 18.713, 25.91, 25.989,
                           33.480, 34.8152, 35.349])
    ion_fwhm = ion_cwaves/inst_res
    ion_names = np.array(["[ArII]", "[ArIII]", "[SIV]", "[NeII]",
                          "[NeIII]", "[SIII] 18", "[OIV]", "[FeII]",
                          "[SIII] 33", "[SiII]", "[FeII]"])
    n_ion = len(ion_cwaves)
    ion_amps = np.full((n_ion), 100.)
    ion_amps_limits = [(0.0, None) for k in range(n_ion)]
    ion_cwaves_limits = [(cwave - 0.05, cwave + 0.05) for cwave in ion_cwaves]
    ion_fwhm_limits = [(cfwhm*0.5, cfwhm*1.5) for cfwhm in ion_fwhm]
    ion_features = (ion_amps, ion_cwaves, ion_fwhm,
                    ion_amps_limits, ion_cwaves_limits, ion_fwhm_limits,
                    ion_names)

    pmodel = PAHFITBase(bb_info=bb_info,
                        dust_features=dust_features,
                        h2_features=h2_features,
                        ion_features=ion_features)

    # print(pmodel.model.amplitude_0.bounds)

    # read in an example spectrum (from M101)
    # obs = Table()
    # obs.read('data/NGC5471_irs.fits')
    # print(obs)
    hdul = fits.open('data/Nucleus_irs.fits')
    obs_x = hdul[1].data['WAVELENGTH']
    obs_y = hdul[1].data['FLUX']
    obs_unc = hdul[1].data['SIGMA']
    obs_npts = hdul[1].data['NPTS']
    hdul.close()

    # fit the model to the data
    fit = LevMarLSQFitter()
    obs_fit = fit(pmodel.model, obs_x, obs_y, weights=1./obs_unc,
                  maxiter=5000)
    # print(fit.fit_info)
    print(obs_fit.name)
    print(obs_fit.param_names)
    print(obs_fit.parameters)

    # plot result
    fontsize = 18
    font = {'size': fontsize}
    matplotlib.rc('font', **font)
    matplotlib.rc('lines', linewidth=2)
    matplotlib.rc('axes', linewidth=2)
    matplotlib.rc('xtick.major', width=2)
    matplotlib.rc('ytick.major', width=2)

    fig, ax = plt.subplots(figsize=(15, 10))

    ax.plot(obs_x, pmodel.model(obs_x))
    ax.plot(obs_x, obs_fit(obs_x))
    ax.plot(obs_x, obs_y)
    ax.set_yscale('log')
    ax.set_xscale('log')

    plt.show()
