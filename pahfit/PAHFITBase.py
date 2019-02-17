from astropy.modeling.functional_models import Gaussian1D

from pahfit.component_models import (BlackBody1D, Drude1D,
                                     S07_attenuation)

from astropy.table import Table, vstack

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
    blackbodies for the continuum, drudes for the dust
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
                 param_info=None,
                 filename=None):
        """
        Setup a variant based on inputs.  Generates an astropy.modeling
        compound model.

        Notes
        -----
        Would be great to rename the parameters such that they uniquely
        identify the component (dust, gas, specific line, etc.).  This is
        possible - say a discussion on the stsci slack channel - James Davies?
        """
        # check that param_info or filename is set
        if filename is None and param_info is None:
            raise ValueError('Either param_info or filename need to be set \
                             when initializing a PAHFITBase object')

        # read in the parameter info from a file
        if filename is not None:
            param_info = self.read(filename)

        bb_info = param_info[0]
        dust_features = param_info[1]
        h2_features = param_info[2]
        ion_features = param_info[3]
        # setup the model
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
                        'mean': x_0_limits[0],
                        'stddev': (fwhms_limits[0][0]/2.355,
                                   fwhms_limits[0][1]/2.355)})
            for k in range(1, n_h2):
                h2_model = h2_model + Gaussian1D(
                    name=names[k],
                    amplitude=amps[k],
                    mean=x_0[k],
                    stddev=fwhms[k]/2.355,
                    bounds={'amplitude': amps_limits[k],
                            'mean': x_0_limits[k],
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
        self.model *= S07_attenuation(name='tau_sil')

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
        # setup the tables for the different components
        bb_table = Table(names=('Name', 'Form',
                                'temp', 'temp_min', 'temp_max', 'temp_fixed',
                                'amp', 'amp_min', 'amp_max', 'amp_fixed'),
                         dtype=('U25', 'U25',
                                'float64', 'float64', 'float64', 'bool',
                                'float64', 'float64', 'float64', 'bool'))
        line_table = Table(names=('Name', 'Form',
                                  'x_0', 'x_0_min', 'x_0_max', 'x_0_fixed',
                                  'amp', 'amp_min', 'amp_max', 'amp_fixed',
                                  'fwhm', 'fwhm_min', 'fwhm_max',
                                  'fwhm_fixed'),
                           dtype=('U25', 'U25',
                                  'float64', 'float64', 'float64', 'bool',
                                  'float64', 'float64', 'float64', 'bool',
                                  'float64', 'float64', 'float64', 'bool'))
        att_table = Table(names=('Name', 'Form',
                                 'amp', 'amp_min', 'amp_max', 'amp_fixed'),
                          dtype=('U25', 'U25',
                                 'float64', 'float64', 'float64', 'bool'))

        for component in obs_fit:
            comp_type = (component.__class__.__name__)

            if comp_type == 'BlackBody1D':
                bb_table.add_row([component.name, comp_type,
                                  component.temperature.value,
                                  component.temperature.bounds[0],
                                  component.temperature.bounds[1],
                                  component.temperature.fixed,
                                  component.amplitude.value,
                                  component.amplitude.bounds[0],
                                  component.amplitude.bounds[1],
                                  component.amplitude.fixed])
            elif comp_type == 'Drude1D':
                line_table.add_row([component.name, comp_type,
                                    component.x_0.value,
                                    component.x_0.bounds[0],
                                    component.x_0.bounds[1],
                                    component.x_0.fixed,
                                    component.amplitude.value,
                                    component.amplitude.bounds[0],
                                    component.amplitude.bounds[1],
                                    component.amplitude.fixed,
                                    component.fwhm.value,
                                    component.fwhm.bounds[0],
                                    component.fwhm.bounds[1],
                                    component.fwhm.fixed])
            elif comp_type == 'Gaussian1D':
                line_table.add_row([component.name, comp_type,
                                    component.mean.value,
                                    component.mean.bounds[0],
                                    component.mean.bounds[1],
                                    component.mean.fixed,
                                    component.amplitude.value,
                                    component.amplitude.bounds[0],
                                    component.amplitude.bounds[1],
                                    component.amplitude.fixed,
                                    2.355*component.stddev.value,
                                    2.355*component.stddev.bounds[0],
                                    2.355*component.stddev.bounds[1],
                                    component.stddev.fixed])
            elif comp_type == 'S07_attenuation':
                att_table.add_row([component.name, comp_type,
                                   component.tau_sil.value,
                                   component.tau_sil.bounds[0],
                                   component.tau_sil.bounds[1],
                                   component.tau_sil.fixed])

        # stack the tables (handles missing columns between tables)
        out_table = vstack([bb_table, line_table, att_table])

        # Writing output table
        out_table.write('{}_output.{}'.format(filename, outform),
                        format=outform, overwrite=True)

    def read(self, filename):
        """
        Read the model parameters from a file.

        Parameters
        ----------
        filename : string
            The name of the input file containing fit results.

        Returns
        -------
        readout : tuple
            Tuple containing dictionaries of all components from
            the input file.
        """
        # Getting file extension
        ext = filename.split('.')[1]

        # Reading the input file as table
        t = Table.read(filename, format=ext)

        # Getting indices for the different components
        bb_ind = np.concatenate(np.argwhere(t['Form'] == 'BlackBody1D'))
        df_ind = np.concatenate(np.argwhere(t['Form'] == 'Drude1D'))
        ga_ind = np.concatenate(np.argwhere(t['Form'] == 'Gaussian1D'))
        at_ind = np.concatenate(np.argwhere(t['Form'] == 'S07_attenuation'))

        # now split the gas emission lines between H2 and ions
        names = [str(i) for i in np.take(t['Name'], ga_ind)]
        h2_temp = np.concatenate(np.where(np.char.find(names, 'H2') >= 0))
        ion_temp = np.concatenate(np.where(np.char.find(names, 'H2') == -1))
        h2_ind = np.take(ga_ind, h2_temp)
        ion_ind = np.take(ga_ind, ion_temp)

        # Obtaining the blackbody components
        bb_names = np.take(t['Name'], bb_ind)
        bb_amps = np.take(t['amp'], bb_ind)
        bb_amp_min = np.take(t['amp_min'], bb_ind)
        bb_amp_max = np.take(t['amp_max'], bb_ind)
        bb_amp_fixed = np.take(t['amp_fixed'], bb_ind)
        bb_amps_limits = [(i, j) for i, j in zip(bb_amp_min, bb_amp_max)]

        # Creating the blackbody dict
        bb_info = {'names': np.array(t['Name'][bb_ind].data),
                   'temps': np.array(t['temp'][bb_ind].data),
                   'amps': np.array(t['amp'][bb_ind].data),
                   'amps_limits': list(zip(t['amp_min'][bb_ind].data,
                                           t['amp_max'][bb_ind].data)),
                   'amps_fixed': np.array(t['amp_fixed'][bb_ind].data)}
        print(bb_info)
        exit()

        # Obtaining the dust features components
        df_amps = np.take(t['amp'], df_ind)
        df_amp_min = np.take(t['amp_min'], df_ind)
        df_amp_max = np.take(t['amp_max'], df_ind)
        df_amps_limits = [(i, j) for i, j in zip(df_amp_min, df_amp_max)]
        df_cwave = np.take(t['x_0'], df_ind)
        df_cwave_min = np.take(t['x_0_min'], df_ind)
        df_cwave_max = np.take(t['x_0_max'], df_ind)
        df_cwave_limits = [(i, j) for i, j in zip(df_cwave_min, df_cwave_max)]
        df_fwhm = np.take(t['fwhm'], df_ind)
        df_fwhm_min = np.take(t['fwhm_min'], df_ind)
        df_fwhm_max = np.take(t['fwhm_max'], df_ind)
        df_fwhm_limits = [(i, j) for i, j in zip(df_fwhm_min, df_fwhm_max)]

        # Creating the dust_features dict
        dust_features = {'amps': df_amps,
                         'x_0': df_cwave,
                         'fwhms': df_fwhm,
                         'amps_limits': df_amps_limits,
                         'x_0_limits': df_cwave_limits,
                         'fwhms_limits': df_fwhm_limits}

        # Obtaining the H2 components
        h2_amps = np.take(t['amp'], h2_ind)
        h2_amp_min = np.take(t['amp_min'], h2_ind)
        h2_amp_max = np.take(t['amp_max'], h2_ind)
        h2_amps_limits = [(i, j) for i, j in zip(h2_amp_min, h2_amp_max)]
        h2_cwave = np.take(t['x_0'], h2_ind)
        h2_cwave_min = np.take(t['x_0_min'], h2_ind)
        h2_cwave_max = np.take(t['x_0_max'], h2_ind)
        h2_cwave_limits = [(i, j) for i, j in zip(h2_cwave_min, h2_cwave_max)]
        h2_fwhm = np.take(t['fwhm'], h2_ind)
        h2_fwhm_min = np.take(t['fwhm_min'], h2_ind)
        h2_fwhm_max = np.take(t['fwhm_max'], h2_ind)
        h2_fwhm_limits = [(i, j) for i, j in zip(h2_fwhm_min, h2_fwhm_max)]
        h2_names = np.take(t['Name'], h2_ind)

        # Creating the H2 dict
        h2_features = {'amps': h2_amps,
                       'x_0': h2_cwave,
                       'fwhms': h2_fwhm,
                       'amps_limits': h2_amps_limits,
                       'x_0_limits': h2_cwave_limits,
                       'fwhms_limits': h2_fwhm_limits,
                       'names': h2_names}

        # Obtaining the ion components
        ion_amps = np.take(t['amp'], ion_ind)
        ion_amp_min = np.take(t['amp_min'], ion_ind)
        ion_amp_max = np.take(t['amp_max'], ion_ind)
        ion_amps_limits = [(i, j) for i, j in zip(ion_amp_min, ion_amp_max)]
        ion_cwave = np.take(t['x_0'], ion_ind)
        ion_cwave_min = np.take(t['x_0_min'], ion_ind)
        ion_cwave_max = np.take(t['x_0_max'], ion_ind)
        ion_cwave_limits = [(i, j)
                            for i, j in zip(ion_cwave_min, ion_cwave_max)]
        ion_fwhm = np.take(t['fwhm'], ion_ind)
        ion_fwhm_min = np.take(t['fwhm_min'], ion_ind)
        ion_fwhm_max = np.take(t['fwhm_max'], ion_ind)
        ion_fwhm_limits = [(i, j) for i, j in zip(ion_fwhm_min, ion_fwhm_max)]
        ion_names = np.take(t['Name'], ion_ind)

        # Creating the ion dict
        ion_features = {'amps': ion_amps,
                        'x_0': ion_cwave,
                        'fwhms': ion_fwhm,
                        'amps_limits': ion_amps_limits,
                        'x_0_limits': ion_cwave_limits,
                        'fwhms_limits': ion_fwhm_limits,
                        'names': ion_names}

        # Create the attenuation dict
        TBD

        # Create output tuple
        readout = (bb_info, dust_features, h2_features, ion_features, at_ind)

        return readout
