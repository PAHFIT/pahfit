from __future__ import (absolute_import, print_function, division)

from astropy.modeling.functional_models import Gaussian1D

from component_models import (BlackBody1D, Drude1D,
                                     S07_attenuation)
from astropy.io import fits

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
            for k in range(n_h2):
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
            for k in range(n_ion):
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

    def save(self, obs_fit, filename):
        """
        Save the model parameters to a file.
        Format TBD
        """
        # Instantiating lists
        bb_temp, bb_amp = ([] for i in range(2))
        pah_amp, pah_x_0, pah_fwhm = ([] for i in range(3))
        line_name, line_amp, line_mean, line_std = ([] for i in range(4))
        att = []
        # TODO: find the variable that sets the components number (currently 55).
        for i in range(55):
            # Getting object name
            comp_name = (obs_fit[i].__class__.__name__)

            if comp_name == 'BlackBody1D':
                bb_temp.append(obs_fit[i].temperature.value)
                bb_amp.append(obs_fit[i].amplitude.value)

            elif comp_name == 'Drude1D':
                pah_amp.append(obs_fit[i].amplitude.value)
                pah_x_0.append(obs_fit[i].x_0.value)
                pah_fwhm.append(obs_fit[i].fwhm.value)

            elif comp_name == 'Gaussian1D':
                line_name.append(obs_fit[i].name)
                line_amp.append(obs_fit[i].amplitude.value)
                line_mean.append(obs_fit[i].mean.value)
                line_std.append(obs_fit[i].stddev.value)

            elif comp_name == 'S07_attenuation':
                att.append(obs_fit[i].tau_si.value)

        # Creating columns for tables
        bb_col1 = fits.Column(name='temperature', format='E', array=bb_temp)
        bb_col2 = fits.Column(name='amplitude', format='E', array=bb_amp)

        pah_col1 = fits.Column(name='amplitude', format='E', array=pah_amp)
        pah_col2 = fits.Column(name='x_0', format='E', array=pah_x_0)
        pah_col3 = fits.Column(name='fwhm', format='E', array=pah_fwhm)

        line_col1 = fits.Column(name='name', format='A10', array=line_name)
        line_col2 = fits.Column(name='amplitude', format='E', array=line_amp)
        line_col3 = fits.Column(name='mean', format='E', array=line_mean)
        line_col4 = fits.Column(name='stddev', format='E', array=line_std)

        tau_si_col = fits.Column(name='tau_si', format='E', array=att)

        # Creating (empty) Primary HDU
        # TODO: populate Primary HDU
        hdu_pr = fits.PrimaryHDU()

        # Creating BinTableHDU for BlackBody1D
        coldefs1 = fits.ColDefs([bb_col1, bb_col2])
        hdu_bb = fits.BinTableHDU.from_columns(coldefs1)
        hdu_bb.header.set('EXTNAME', 'BlackBody1D')

        # Creating BinTableHDU for Drude1D
        coldefs2 = fits.ColDefs([pah_col1, pah_col2, pah_col3])
        hdu_pah = fits.BinTableHDU.from_columns(coldefs2)
        hdu_pah.header.set('EXTNAME', 'Drude1D')

        # Creating BinTableHDU for Gaussian1D
        coldefs3 = fits.ColDefs([line_col1, line_col2, line_col3, line_col4])
        hdu_line = fits.BinTableHDU.from_columns(coldefs3)
        hdu_line.header.set('EXTNAME', 'Gaussian1D')

        # Creating BinTableHDU for S07_attenuation
        coldefs4 = fits.ColDefs([tau_si_col])
        hdu_att = fits.BinTableHDU.from_columns(coldefs4)
        hdu_att.header.set('EXTNAME', 'S07_attenuation')

        # Joining extensions to a list of HDUs.
        hdul = fits.HDUList([hdu_pr, hdu_bb, hdu_pah, hdu_line, hdu_att])

        # Writing fits to disk.
        hdul.writeto('{}_output.fits'.format(filename), overwrite=True)


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
