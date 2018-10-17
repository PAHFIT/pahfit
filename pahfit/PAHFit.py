from __future__ import (absolute_import, print_function, division)

import numpy as np

import matplotlib.pyplot as plt
import matplotlib

# import astropy.units as u
from astropy.modeling import Fittable1DModel, Parameter
from astropy.modeling.functional_models import (Lorentz1D,
                                                Gaussian1D)
# from astropy.modeling.blackbody import BlackBody1D


class BlackBody1D(Fittable1DModel):
    """
    Current astropy BlackBody1D does not play well with Lorentz1D and Gauss1D
    maybe...need to check again...possibly a units issue
    """
    inputs = ('x',)
    outputs = ('y',)

    amplitude = Parameter()
    temperature = Parameter()

    @staticmethod
    def evaluate(x, amplitude, temperature):
        """
        """
        return (amplitude*((9.7/x)**2)*3.97289e13/x**3
                / (np.exp(1.4387752e4/x/temperature)-1.))


def PAHFit_Spitzer_Exgal(bb_info,
                         dust_features,
                         h2_features,
                         ion_features):
    """
    PAHFit "Classic"

    For Spitzer IRS spectra (low-res) and optimized for extragalactic objects
    """
    amps, temps = bb_info
    n_bb = len(amps)
    cont_model = BlackBody1D(temperature=temps[0],
                             amplitude=amps[0],
                             fixed={'temperature': True})
    cont_model.amplitude.bounds = (0.0, None)
    for k in range(1, n_bb):
        new_model = BlackBody1D(temperature=temps[k],
                                amplitude=amps[k],
                                fixed={'temperature': True})
        new_model.amplitude.bounds = (0.0, None)
        cont_model = cont_model + new_model

    n_df = len(dust_features[0])
    amps, cwaves, fwhms = dust_features
    df_model = Lorentz1D(amplitude=amps[0],
                         x_0=cwaves[0],
                         fwhm=fwhms[0],
                         bounds={'amplitude': (0, 1e6)})
    for k in range(1, n_df):
        df_model = df_model + Lorentz1D(amplitude=amps[k],
                                        x_0=cwaves[k],
                                        fwhm=fwhms[k])

    n_h2 = len(h2_features[0])
    amps, cwaves, fwhms, names = h2_features
    h2_model = Gaussian1D(amplitude=amps[0],
                          mean=cwaves[0],
                          stddev=fwhms[0])
    for k in range(n_h2):
        h2_model = h2_model + Gaussian1D(amplitude=amps[k],
                                         mean=cwaves[k],
                                         stddev=fwhms[k]/2.355)

    n_ion = len(ion_features[0])
    amps, cwaves, fwhms, names = ion_features
    ion_model = Gaussian1D(amplitude=amps[0],
                           mean=cwaves[0],
                           stddev=fwhms[0])
    for k in range(n_ion):
        ion_model = ion_model + Gaussian1D(amplitude=amps[k],
                                           mean=cwaves[k],
                                           stddev=fwhms[k]/2.355)

    return cont_model + df_model + h2_model + ion_model


if __name__ == '__main__':
    # define the instrument resolution
    inst_res = 100.

    # the blackbodies for the continuum
    bb_temps = [3000., 2000., 1000., 800., 600., 450., 300., 200., 135.,
                90., 65., 50., 40., 35.]
    bb_amps = np.full((len(bb_temps)), 1.0e-10)
    bb_info = (bb_amps, bb_temps)

    # the dust features (mainly the PAH/aromatic features)
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
    dust_features = (dust_features_amps, dust_features_cwave,
                     dust_features_fwhm)

    # define the H2 features
    h2_cwaves = np.array([5.5115, 6.1088, 6.9091, 8.0258,
                          9.6649, 12.2785, 17.0346, 28.2207])
    h2_fwhm = h2_cwaves/inst_res
    h2_names = np.array(["H2 S(7)", "H2 S(6)", "H2 S(5)", "H2 S(4)",
                         "H2 S(3)", "H2 S(2)", "H2 S(1)", "H2 S(0)"])
    n_h2 = len(h2_cwaves)
    h2_amps = np.full((n_h2), 100.)
    h2_features = (h2_amps, h2_cwaves, h2_fwhm, h2_names)

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
    ion_features = (ion_amps, ion_cwaves, ion_fwhm, ion_names)

    pmodel = PAHFit_Spitzer_Exgal(bb_info,
                                  dust_features,
                                  h2_features,
                                  ion_features)

    print(pmodel)
    print(pmodel.amplitude_0.bounds)

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
