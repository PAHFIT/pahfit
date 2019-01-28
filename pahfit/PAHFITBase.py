from __future__ import (absolute_import, print_function, division)

from astropy.modeling.functional_models import Gaussian1D

from component_models import (BlackBody1D, Drude1D,
                                     S07_attenuation)

from astropy.table import Table

import numpy as np

__all__ = ['PAHFITBase']


class PAHFITBase():
    """
    Base class for PAHFIT variants. Each variant nominally specifies the valid
    wavelength range, instrument, and type of astronomical objects.

    For example, the original IDL version of PAHFIT was valid for
    Spitzer/IRS spectra (5-38 micron) and observations of parts or all of
    external galaxies.

    Mainly sets up the astropy.modeling compound model
    that includes all the different components including
    blackbodies for the continuum, lorentizians for the dust
    emission features, and Gaussians for the gas emission features.

    Parameters
    ----------
    bb_info : dict
        dict with {'amps', 'temps', 'amps_limits'}, each a vector
    dust_features,
    h2_features,
    ion_features : dict
        dict with {amps, x_0, fwhm,
                   amps_limits, x_0_limits, fwhms_limits}, each a vector
    """

    def __init__(self,
                 bb_info,
                 dust_features,
                 h2_features,
                 ion_features):
        """
        Setup a variant based on inputs.  Generates an astropy.modeling
        compound model.

        Note
        ----
        Would be great to rename the parameters such that they uniquely
        identify the component (dust, gas, specific line, etc.).  This is
        possible - say a discussion on the stsci slack channel - James Davies?
        """
        model_comps = []
        self.bb_info = bb_info
        if bb_info is not None:
            amps = bb_info['amps']
            temps = bb_info['temps']
            amps_limits = bb_info['amps_limits']
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
            amps = dust_features['amps']
            x_0 = dust_features['x_0']
            fwhms = dust_features['fwhms']
            amps_limits = dust_features['amps_limits']
            x_0_limits = dust_features['x_0_limits']
            fwhms_limits = dust_features['fwhms_limits']
            n_df = len(amps)
            df_model = Drude1D(amplitude=amps[0],
                               x_0=x_0[0],
                               fwhm=fwhms[0],
                               bounds={'amplitude': amps_limits[0],
                                       'x_0': x_0_limits[0],
                                       'fwhm': fwhms_limits[0]})
            for k in range(1, n_df):
                df_model = df_model + Drude1D(
                    amplitude=amps[k],
                    x_0=x_0[k],
                    fwhm=fwhms[k],
                    bounds={'amplitude': amps_limits[k],
                            'x_0': x_0_limits[k],
                            'fwhm': fwhms_limits[k]})

            self.df_model = df_model
            model_comps.append(df_model)

        self.h2_features = h2_features
        if h2_features is not None:
            amps = h2_features['amps']
            x_0 = h2_features['x_0']
            fwhms = h2_features['fwhms']
            amps_limits = h2_features['amps_limits']
            x_0_limits = h2_features['x_0_limits']
            fwhms_limits = h2_features['fwhms_limits']
            names = h2_features['names']
            n_h2 = len(amps)
            h2_model = Gaussian1D(
                name=names[0],
                amplitude=amps[0],
                mean=x_0[0],
                stddev=fwhms[0]/2.355,
                bounds={'amplitude': amps_limits[0],
                        'x_0': x_0_limits[0],
                        'stddev': (fwhms_limits[0][0]/2.355,
                                   fwhms_limits[0][1]/2.355)})
            for k in range(1, n_h2):
                h2_model = h2_model + Gaussian1D(
                    name=names[k],
                    amplitude=amps[k],
                    mean=x_0[k],
                    stddev=fwhms[k]/2.355,
                    bounds={'amplitude': amps_limits[k],
                            'x_0': x_0_limits[k],
                            'stddev': (fwhms_limits[k][0]/2.355,
                                       fwhms_limits[k][1]/2.355)})

            self.h2_model = h2_model
            model_comps.append(h2_model)

        self.ion_features = ion_features
        if ion_features is not None:
            amps = ion_features['amps']
            x_0 = ion_features['x_0']
            fwhms = ion_features['fwhms']
            amps_limits = ion_features['amps_limits']
            x_0_limits = ion_features['x_0_limits']
            fwhms_limits = ion_features['fwhms_limits']
            names = ion_features['names']
            n_ion = len(amps)
            ion_model = Gaussian1D(
                name=names[0],
                amplitude=amps[0],
                mean=x_0[0],
                stddev=fwhms[0]/2.355,
                bounds={'amplitude': amps_limits[0],
                        'x_0': x_0_limits[0],
                        'stddev': (fwhms_limits[0][0]/2.355,
                                   fwhms_limits[0][1]/2.355)})
            for k in range(1, n_ion):
                ion_model = ion_model + Gaussian1D(
                    name=names[k],
                    amplitude=amps[k],
                    mean=x_0[k],
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

        # need to make the type of attenuation model a passed variable
        self.model *= S07_attenuation()

    def plot(self, ax, x, y, model):
        """
        Plot model using axis object.

        Parameters
        ----------
        ax : matplotlib.axis object
            where to put the plot
        x : floats
            wavelength points
        y : floats
            observed spectrum
        model : PAHFITBase model
            model giving all the components and parameters
        """
        ax.plot(x, model(x)/x, 'g-')
        ax.plot(x, y/x, 'ks', fillstyle='none')

        # get the extinction model (probably a better way to do this)
        for cmodel in model:
            if isinstance(cmodel, S07_attenuation):
                ax.plot(x, cmodel(x)*max(y/x), 'k--')
                ext_model = cmodel(x)

        # create the continum compound model (base for plotting lines)
        cont_components = []
        for cmodel in model:
            if isinstance(cmodel, BlackBody1D):
                cont_components.append(cmodel)
                # plot as we go
                ax.plot(x, cmodel(x)*ext_model/x, 'r-')
        cont_model = cont_components[0]
        for cmodel in cont_components[1:]:
            cont_model += cmodel
        cont_y = cont_model(x)

        # now plot the dust and gas lines
        for cmodel in model:
            if isinstance(cmodel, Gaussian1D):
                ax.plot(x, (cont_y + cmodel(x))*ext_model/x,
                        color='tab:purple')
            if isinstance(cmodel, Drude1D):
                ax.plot(x, (cont_y + cmodel(x))*ext_model/x,
                        color='tab:blue')

        ax.plot(x, cont_y*ext_model/x, 'k-')

        ax.set_xlabel(r'$\lambda$ [$\mu m$]')
        ax.set_ylabel(r'$\nu F_{\nu}$')

        # print(obs_fit.name)
        # print(obs_fit.param_names)
        # print(obs_fit.parameters)

    def save(self, obs_fit, filename, outform):
        """
        Save the model parameters to a user defined file format.

        Parameters
        ----------
        obs_fit : PAHFITBase model
            Model giving all the components and parameters.
        filename : string
            String used to name the output file.
            Currently using the input data file name.
        outform : string
            Sets the output file format (ascii, fits, csv, etc.).

        """
        # Instantiating lists
        Name, Form, Fixed = ([] for i in range(3))
        amp, amp_min, amp_max = ([] for i in range(3))
        x_0, x_0_min, x_0_max = ([] for i in range(3))
        fwhm, fwhm_min, fwhm_max = ([] for i in range(3))
        mean, stddev, stddev_min, stddev_max = ([] for i in range(4))

        # Instantiating mask lists
        x_0_mask, x_0_min_mask, x_0_max_mask = ([] for i in range(3))
        fwhm_mask, fwhm_min_mask, fwhm_max_mask = ([] for i in range(3))
        mean_mask, stddev_mask, stddev_min_mask, stddev_max_mask = ([] for i in range(4))

        # Dust feature component index
        DFi = 0

        for component in obs_fit:
            # Getting object name
            comp_name = (component.__class__.__name__)

            if comp_name == 'BlackBody1D':
                Name.append('BB{}'.format(int(component.temperature.value)))
                Form.append(comp_name)
                Fixed.append(component.temperature.fixed)
                amp.append(component.amplitude.value)
                amp_min.append(component.amplitude.bounds[0])
                amp_max.append(component.amplitude.bounds[1])
                x_0.append(np.nan)
                x_0_min.append(np.nan)
                x_0_max.append(np.nan)
                x_0_mask.append(True)
                x_0_min_mask.append(True)
                x_0_max_mask.append(True)
                fwhm.append(np.nan)
                fwhm_min.append(np.nan)
                fwhm_max.append(np.nan)
                fwhm_mask.append(True)
                fwhm_min_mask.append(True)
                fwhm_max_mask.append(True)
                mean.append(np.nan)
                mean_mask.append(True)
                stddev.append(np.nan)
                stddev_min.append(np.nan)
                stddev_max.append(np.nan)
                stddev_mask.append(True)
                stddev_min_mask.append(True)
                stddev_max_mask.append(True)

            elif comp_name == 'Drude1D':
                DFi += 1
                Name.append('DF{}'.format(DFi))
                Form.append(comp_name)
                Fixed.append(False)
                amp.append(component.amplitude.value)
                amp_min.append(component.amplitude.bounds[0])
                amp_max.append(component.amplitude.bounds[1])
                x_0.append(component.x_0.value)
                x_0_min.append(component.x_0.bounds[0])
                x_0_max.append(component.x_0.bounds[1])
                x_0_mask.append(False)
                x_0_min_mask.append(False)
                x_0_max_mask.append(False)
                fwhm.append(component.fwhm.value)
                fwhm_min.append(component.fwhm.bounds[0])
                fwhm_max.append(component.fwhm.bounds[1])
                fwhm_mask.append(False)
                fwhm_min_mask.append(False)
                fwhm_max_mask.append(False)
                mean.append(np.nan)
                mean_mask.append(True)
                stddev.append(np.nan)
                stddev_min.append(np.nan)
                stddev_max.append(np.nan)
                stddev_mask.append(True)
                stddev_min_mask.append(True)
                stddev_max_mask.append(True)

            elif comp_name == 'Gaussian1D':
                Name.append(component.name)
                Form.append(comp_name)
                Fixed.append(False)
                amp.append(component.amplitude.value)
                amp_min.append(component.amplitude.bounds[0])
                amp_max.append(component.amplitude.bounds[1])
                x_0.append(np.nan)
                x_0_min.append(np.nan)
                x_0_max.append(np.nan)
                x_0_mask.append(True)
                x_0_min_mask.append(True)
                x_0_max_mask.append(True)
                fwhm.append(component.fwhm)
                fwhm_min.append(np.nan)
                fwhm_max.append(np.nan)
                fwhm_mask.append(False)
                fwhm_min_mask.append(True)
                fwhm_max_mask.append(True)
                mean.append(component.mean.value)
                mean_mask.append(False)
                stddev.append(component.stddev.value)
                stddev_min.append(component.stddev.bounds[0])
                stddev_max.append(component.stddev.bounds[1])
                stddev_mask.append(False)
                stddev_min_mask.append(False)
                stddev_max_mask.append(False)

            elif comp_name == 'S07_attenuation':
                Name.append('tau_si')
                Form.append(comp_name)
                Fixed.append(False)
                amp.append(component.tau_si.value)
                amp_min.append(component.tau_si.bounds[0])
                amp_max.append(component.tau_si.bounds[1])
                x_0.append(np.nan)
                x_0_min.append(np.nan)
                x_0_max.append(np.nan)
                x_0_mask.append(True)
                x_0_min_mask.append(True)
                x_0_max_mask.append(True)
                fwhm.append(np.nan)
                fwhm_min.append(np.nan)
                fwhm_max.append(np.nan)
                fwhm_mask.append(True)
                fwhm_min_mask.append(True)
                fwhm_max_mask.append(True)
                mean.append(np.nan)
                mean_mask.append(True)
                stddev.append(np.nan)
                stddev_min.append(np.nan)
                stddev_max.append(np.nan)
                stddev_mask.append(True)
                stddev_min_mask.append(True)
                stddev_max_mask.append(True)

        # Table column names
        colnames = ['Name', 'Form', 'Fixed',
                    'amp', 'amp_min', 'amp_max',
                    'x_0', 'x_0_min', 'x_0_max',
                    'fwhm', 'fwhm_min', 'fwhm_max',
                    'mean', 'stddev', 'stddev_min', 'stddev_max']

        # Creating Table
        t = Table([Name, Form, Fixed,
                   amp, amp_min, amp_max,
                   x_0, x_0_min, x_0_max,
                   fwhm, fwhm_min, fwhm_max,
                   mean, stddev, stddev_min, stddev_max],
                  dtype=('U25', 'U25', 'U25',
                         'f8', 'f8', 'f8',
                         'f8', 'f8', 'f8',
                         'f8', 'f8', 'f8',
                         'f8', 'f8', 'f8', 'f8'),
                  names=colnames, masked=True)

        # Assigning masks
        t['x_0'].mask = x_0_mask
        t['x_0_min'].mask = x_0_min_mask
        t['x_0_max'].mask = x_0_max_mask

        t['fwhm'].mask = fwhm_mask
        t['fwhm_min'].mask = fwhm_min_mask
        t['fwhm_max'].mask = fwhm_max_mask

        t['mean'].mask = mean_mask

        t['stddev'].mask = stddev_mask
        t['stddev_min'].mask = stddev_min_mask
        t['stddev_max'].mask = stddev_max_mask

        # Writing output table
        t.write('{}_output.{}'.format(filename, outform), format=outform, overwrite=True)

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
