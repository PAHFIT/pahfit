from __future__ import (absolute_import, print_function, division)

import numpy as np

from astropy.modeling import Fittable1DModel, Parameter
from astropy.modeling.functional_models import (Lorentz1D,
                                                Gaussian1D)

__all__ = ['PAHFITBase']


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


class PAHFITBase():
    """
    Base class for PAHFIT variants.  Each varient nomially specifies the valid
    wavelength range, instrument, and type of astronomnical objects.

    For example, the original IDL version of PAHFIT was valid for
    Spitzer/IRS spectra (5-38 micron) and observations of parts or all of
    external galaxies.
    """

    def __init__(self,
                 bb_info,
                 dust_features,
                 h2_features,
                 ion_features):
        """
        Setup a varient based on inputs.  Generates an astropy.modeling
        compound model.

        Parameters
        ----------
        bb_info : type
            Description of parameter `bb_info`.
        dust_features : type
            Description of parameter `dust_features`.
        h2_features : type
            Description of parameter `h2_features`.
        ion_features : type
            Description of parameter `ion_features`.

        Returns
        -------
        astropy.modeling compond model
            Model that includes all the different components including
            blackbodies for the continuum, lorentizians for the dust
            emission features, and Gaussians for the gas emission features.

        Note
        ----
        Would be great to rename the parameters such that they uniquely
        identify the component (dust, gas, specific line, etc.).  This is
        possible - say a discussion on the stsci slack channel - James Davies?
        """
        model_comps = []
        self.bb_info = bb_info
        if bb_info is not None:
            amps, temps, amps_limits = bb_info
            n_bb = len(amps)
            cont_model = BlackBody1D(temperature=temps[0],
                                     amplitude=amps[0],
                                     fixed={'temperature': True})
            cont_model.amplitude.bounds = amps_limits[0]
            for k in range(1, n_bb):
                new_model = BlackBody1D(temperature=temps[k],
                                        amplitude=amps[k],
                                        fixed={'temperature': True})
                new_model.amplitude.bounds = amps_limits[k]
                cont_model = cont_model + new_model

            self.cont_model = cont_model
            model_comps.append(cont_model)

        # dust features should be a Drude profile
        # using Lorentizian for prototyping
        # need to define the appropriate Drude1D astropy.model
        self.dust_features = dust_features
        if dust_features is not None:
            amps, cwaves, fwhms, amps_limits, \
                x_0_limits, fwhms_limits = dust_features
            n_df = len(amps)
            df_model = Lorentz1D(amplitude=amps[0],
                                 x_0=cwaves[0],
                                 fwhm=fwhms[0],
                                 bounds={'amplitude': amps_limits[0],
                                         'x_0': x_0_limits[0],
                                         'fwhm': fwhms_limits[0]})
            for k in range(1, n_df):
                df_model = df_model + Lorentz1D(
                    amplitude=amps[k],
                    x_0=cwaves[k],
                    fwhm=fwhms[k],
                    bounds={'amplitude': amps_limits[k],
                            'x_0': x_0_limits[k],
                            'fwhm': fwhms_limits[k]})

            self.df_model = df_model
            model_comps.append(df_model)

        self.h2_features = h2_features
        if h2_features is not None:
            amps, cwaves, fwhms, \
                amps_limits, cwaves_limits, fwhms_limits, \
                names = h2_features
            n_h2 = len(amps)
            h2_model = Gaussian1D(
                amplitude=amps[0],
                mean=cwaves[0],
                stddev=fwhms[0]/2.355,
                bounds={'amplitude': amps_limits[0],
                        'x_0': x_0_limits[0],
                        'stddev': (fwhms_limits[0][0]/2.355,
                                   fwhms_limits[0][1]/2.355)})
            for k in range(n_h2):
                h2_model = h2_model + Gaussian1D(
                    amplitude=amps[k],
                    mean=cwaves[k],
                    stddev=fwhms[k]/2.355,
                    bounds={'amplitude': amps_limits[k],
                            'x_0': x_0_limits[k],
                            'stddev': (fwhms_limits[k][0]/2.355,
                                       fwhms_limits[k][1]/2.355)})

            self.h2_model = h2_model
            model_comps.append(h2_model)

        self.ion_features = ion_features
        if ion_features is not None:
            n_ion = len(ion_features[0])
            amps, cwaves, fwhms, \
                amps_limits, cwaves_limits, fwhms_limits, \
                names = ion_features
            ion_model = Gaussian1D(
                name=names[0],
                amplitude=amps[0],
                mean=cwaves[0],
                stddev=fwhms[0]/2.355,
                bounds={'amplitude': amps_limits[0],
                        'x_0': x_0_limits[0],
                        'stddev': (fwhms_limits[0][0]/2.355,
                                   fwhms_limits[0][1]/2.355)})
            for k in range(n_ion):
                ion_model = ion_model + Gaussian1D(
                    name=names[k],
                    amplitude=amps[k],
                    mean=cwaves[k],
                    stddev=fwhms[k]/2.355,
                    bounds={'amplitude': amps_limits[k],
                            'x_0': x_0_limits[k],
                            'stddev': (fwhms_limits[k][0]/2.355,
                                       fwhms_limits[k][1]/2.355)})
            self.ion_model = ion_model
            model_comps.append(ion_model)

        self.model = model_comps[0]
        for cmodel in model_comps[1:]:
            self.model += cmodel

    def plot(self, ax):
        """
        Plot model using axis object.
        """
        pass

    def save(self, filename):
        """
        Save the model parameters to a file.
        Format TBD
        """
        pass

    def read(self, filename):
        """
        Read the model parameters from a file.
        Format TBD
        This could be the location of how the data files giving the
        model packs are read.  Could also include an optional name of
        such a file for the init function and make all the inputs
        optional.  Probably cleaner that way.
        """
        pass
