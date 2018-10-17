from __future__ import (absolute_import, print_function, division)

import numpy as np

import matplotlib.pyplot as plt
import matplotlib

from astropy.modeling.functional_models import Lorentz1D


if __name__ == '__main__':
    dust_features_cwave = np.array([5.27, 5.70, 6.22, 6.69,
                                    7.42, 7.60, 7.85, 8.33,
                                    8.61, 10.68, 11.23, 11.33,
                                    11.99, 12.62, 12.69, 13.48,
                                    14.04, 14.19, 15.9, 16.45,
                                    17.04, 17.375, 17.87, 18.92,
                                    33.1])

    dust_features_fwhm = np.array([0.034, 0.035, 0.030, 0.07,
                                   0.126, 0.044, 0.053, 0.05,
                                   0.039, 0.02, 0.012, 0.032,
                                   0.045, 0.042, 0.013, 0.04,
                                   0.016, 0.025, 0.02, 0.014,
                                   0.065, 0.012, 0.016, 0.019,
                                   0.05])
    n_df = len(dust_features_cwave)
    dust_features_amps = np.full((n_df), 100.)

    h2_cwaves = np.array([5.5115, 6.1088, 6.9091, 8.0258,
                          9.6649, 12.2785, 17.0346, 28.2207])
    h2_names = np.array(["H2 S(7)", "H2 S(6)", "H2 S(5)", "H2 S(4)",
                         "H2 S(3)", "H2 S(2)", "H2 S(1)", "H2 S(0)"])
    n_h2 = len(h2_cwaves)
    h2_amps = np.full((n_h2), 100.)

    ion_cwaves = np.array([6.985274, 8.99138, 10.5105, 12.813,
                           15.555, 18.713, 25.91, 25.989,
                           33.480, 34.8152, 35.349])
    ion_names = np.array(["[ArII]", "[ArIII]", "[SIV]", "[NeII]",
                          "[NeIII]", "[SIII] 18", "[OIV]", "[FeII]",
                          "[SIII] 33", "[SiII]", "[FeII]"])
    n_ion = len(ion_cwaves)
    ion_amps = np.full((n_ion), 100.)

    pmodel = Lorentz1D(amplitude=dust_features_amps[0],
                       x_0=dust_features_cwave[0],
                       fwhm=dust_features_fwhm[0])
    for k in range(1, n_df):
        pmodel = pmodel + Lorentz1D(amplitude=dust_features_amps[k],
                                    x_0=dust_features_cwave[k],
                                    fwhm=dust_features_fwhm[k])
    for k in range(n_h2):
        pmodel = pmodel + Guass

    print(pmodel)

    x = np.arange(5.0, 40.0, 0.1)

    # plot result
    fontsize = 18
    font = {'size': fontsize}
    matplotlib.rc('font', **font)
    matplotlib.rc('lines', linewidth=2)
    matplotlib.rc('axes', linewidth=2)
    matplotlib.rc('xtick.major', width=2)
    matplotlib.rc('ytick.major', width=2)

    fig, many_ax = plt.subplots(figsize=(15, 10))

    plt.plot(x, pmodel(x))

    plt.show()
