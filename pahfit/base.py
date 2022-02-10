from astropy.modeling.functional_models import Gaussian1D

from pahfit.component_models import BlackBody1D, S07_attenuation

from astropy.table import Table, vstack
from astropy.modeling.physical_models import Drude1D

from scipy import interpolate

import numpy as np

import matplotlib as mpl

from pahfit.feature_strengths import pah_feature_strength, line_strength, featcombine, eqws

__all__ = ["PAHFITBase"]


def _ingest_limits(min_vals, max_vals):
    """
    Ingest the limits read from a file and generate the appropriate
    internal format (list of tuples).  Needed to handle the case
    where a limit is not desired as numpy arrays cannot have elements
    of None, instead a value of nan is used.

    Limits that are not set are designated as 'nan' in files and
    these are changed to the python None to be compatible with
    the astropy.modeling convention.

    Parameters
    ----------
    min_vals,
    max_vals : numpy.array
        min/max values of the limits for a parameter
        nan designates no limit

    Returns
    -------
    plimits : list of tuples
        tuples give the min/max limits for a parameter
    """
    plimits = []
    for cmin, cmax in zip(min_vals, max_vals):
        if np.isnan(cmin):
            cmin = None
        if np.isnan(cmax):
            cmax = None
        plimits.append((cmin, cmax))

    return plimits


def _ingest_fixed(fixed_vals):
    """
    Ingest the fixed value read from a file and generate the appropriate
    internal format (list of booleans).  Needed as booleans but
    represented in files as True/False strings.

    Parameters
    ----------
    min_vals : numpy.array (string)
        fixed designations

    Returns
    -------
    pfixed : numpy.array (boolean)
        True/False designation for parameters
    """
    pfixed = []
    for cfixed in fixed_vals:
        if cfixed == "True":
            cfixed = True
        if cfixed == "False":
            cfixed = False
        pfixed.append(cfixed)

    return pfixed


class PAHFITBase:
    """
    Base class for PAHFIT variants. Each variant nominally specifies the valid
    wavelength range, instrument, and type of astronomical objects.

    For example, the original IDL version of PAHFIT was valid for
    Spitzer/IRS SL/LL spectra (5-38 micron) and observations of parts or all
    of external galaxies.

    Mainly sets up the astropy.modeling compound model
    that includes all the different components including
    blackbodies for the continuum, drudes for the dust
    emission features, and Gaussians for the gas emission features.

    Parameters
    ----------
    obs_x and obs_y: np.array
        the input spectrum

    filename: string
        filename giving the pack that contains all the
        info described for param_info

    tformat: string
        table format of filename (compatible with astropy Table.read)

    param_info: tuple of dics
        The dictonaries contain info for each type of component.  Each
        component of the dictonaries is a vector.
        bb_info -
        dict with {name, temps, temps_limits, temps_fixed,
        amps, amps_limits, amps_fixed}, each a vector
        dust_features, h2_features, ion_features -
        dict with {name amps, amps_limits, amps_fixed,
        x_0, x_0_limits, x_0_fixed, fwhms, fwhms_limits, fwhm_fixed}.
    """

    def __init__(self, obs_x, obs_y, estimate_start=False, param_info=None, filename=None, tformat=None):
        """
        Setup a variant based on inputs.  Generates an astropy.modeling
        compound model.
        """
        # check that param_info or filename is set
        if filename is None and param_info is None:
            raise ValueError(
                "Either param_info or filename need to be set \
                             when initializing a PAHFITBase object"
            )

        # read in the parameter info from a file
        if filename is not None:
            param_info = self.read(filename, tformat=tformat)

        if estimate_start:
            # guess values and update starting point (if not set fixed) based on the input spectrum
            param_info = self.estimate_init(obs_x, obs_y, param_info)

        self.param_info = param_info

        bb_info = param_info[0]
        dust_features = param_info[1]
        h2_features = param_info[2]
        ion_features = param_info[3]
        att_info = param_info[4]

        # setup the model
        self.bb_info = bb_info
        if bb_info is not None:
            # 1st component defines the overall model variable
            self.model = BlackBody1D(
                name=bb_info["names"][0],
                temperature=bb_info["temps"][0],
                amplitude=bb_info["amps"][0],
                bounds={
                    "temperature": bb_info["temps_limits"][0],
                    "amplitude": bb_info["amps_limits"][0],
                },
                fixed={
                    "temperature": bb_info["temps_fixed"][0],
                    "amplitude": bb_info["amps_fixed"][0],
                },
            )
            for k in range(1, len(bb_info["names"])):
                self.model += BlackBody1D(
                    name=bb_info["names"][k],
                    temperature=bb_info["temps"][k],
                    amplitude=bb_info["amps"][k],
                    bounds={
                        "temperature": bb_info["temps_limits"][k],
                        "amplitude": bb_info["amps_limits"][k],
                    },
                    fixed={
                        "temperature": bb_info["temps_fixed"][k],
                        "amplitude": bb_info["amps_fixed"][k],
                    },
                )

        self.dust_features = dust_features
        if dust_features is not None:
            for k in range(len(dust_features["names"])):
                self.model += Drude1D(
                    name=dust_features["names"][k],
                    amplitude=dust_features["amps"][k],
                    x_0=dust_features["x_0"][k],
                    fwhm=dust_features["fwhms"][k],
                    bounds={
                        "amplitude": dust_features["amps_limits"][k],
                        "x_0": dust_features["x_0_limits"][k],
                        "fwhm": dust_features["fwhms_limits"][k],
                    },
                    fixed={
                        "amplitude": dust_features["amps_fixed"][k],
                        "x_0": dust_features["x_0_fixed"][k],
                        "fwhm": dust_features["fwhms_fixed"][k],
                    },
                )

        self.h2_features = h2_features
        if h2_features is not None:
            for k in range(len(h2_features["names"])):
                self.model += Gaussian1D(
                    name=h2_features["names"][k],
                    amplitude=h2_features["amps"][k],
                    mean=h2_features["x_0"][k],
                    stddev=h2_features["fwhms"][k] / 2.355,
                    bounds={
                        "amplitude": h2_features["amps_limits"][k],
                        "mean": h2_features["x_0_limits"][k],
                        "stddev": (
                            h2_features["fwhms_limits"][k][0] / 2.355,
                            h2_features["fwhms_limits"][k][1] / 2.355,
                        ),
                    },
                    fixed={
                        "amplitude": h2_features["amps_fixed"][k],
                        "mean": h2_features["x_0_fixed"][k],
                        "stddev": h2_features["fwhms_fixed"][k],
                    },
                )

        self.ion_features = ion_features
        if ion_features is not None:
            for k in range(len(ion_features["names"])):
                self.model += Gaussian1D(
                    name=ion_features["names"][k],
                    amplitude=ion_features["amps"][k],
                    mean=ion_features["x_0"][k],
                    stddev=ion_features["fwhms"][k] / 2.355,
                    bounds={
                        "amplitude": ion_features["amps_limits"][k],
                        "mean": ion_features["x_0_limits"][k],
                        "stddev": (
                            ion_features["fwhms_limits"][k][0] / 2.355,
                            ion_features["fwhms_limits"][k][1] / 2.355,
                        ),
                    },
                    fixed={
                        "amplitude": ion_features["amps_fixed"][k],
                        "mean": ion_features["x_0_fixed"][k],
                        "stddev": ion_features["fwhms_fixed"][k],
                    },
                )

        # apply the attenuation to *all* the components
        self.model *= S07_attenuation(
            name=att_info["names"][0],
            tau_sil=att_info["amps"][0],
            bounds={"tau_sil": att_info["amps_limits"][0]},
            fixed={"tau_sil": att_info["amps_fixed"][0]},
        )

    def plot(self, axs, x, y, yerr, model, scalefac_resid=2):
        """
        Plot model using axis object.

        Parameters
        ----------
        axs : matplotlib.axis objects
            where to put the plot
        x : floats
            wavelength points
        y : floats
            observed spectrum
        yerr: floats
            observed spectrum uncertainties
        model : PAHFITBase model (astropy modeling CompoundModel)
            model giving all the components and parameters
        scalefac_resid : float
            Factor multiplying the standard deviation of the residuals to adjust plot limits
        """
        # remove units if they are present
        if hasattr(x, "value"):
            x = x.value
        if hasattr(y, "value"):
            y = y.value
        if hasattr(yerr, "value"):
            yerr = yerr.value

        # spectrum and best fit model

        ax = axs[0]
        ax.set_yscale("linear")
        ax.set_xscale("log")
        ax.minorticks_on()
        ax.tick_params(axis="both", which='major', top="on", right="on", direction='in', length=10)
        ax.tick_params(axis="both", which='minor', top="on", right="on", direction='in', length=5)

        ax_att = ax.twinx()  # axis for plotting the extinction curve
        ax_att.tick_params(which='minor', direction="in", length=5)
        ax_att.tick_params(which='major', direction='in', length=10)
        ax_att.minorticks_on()
        # get the extinction model (probably a better way to do this)
        for cmodel in model:
            if isinstance(cmodel, S07_attenuation):
                ax_att.plot(x, cmodel(x), "k--", alpha=0.5)
                ext_model = cmodel(x)
        ax_att.set_ylabel("Attenuation")
        ax_att.set_ylim(0, 1.1)

        # Define legend lines
        Leg_lines = [mpl.lines.Line2D([0], [0], color="k", linestyle="--", lw=2),
                     mpl.lines.Line2D([0], [0], color="#FE6100", lw=2),
                     mpl.lines.Line2D([0], [0], color="#648FFF", lw=2, alpha=0.5),
                     mpl.lines.Line2D([0], [0], color="#DC267F", lw=2, alpha=0.5),
                     mpl.lines.Line2D([0], [0], color="#785EF0", lw=2, alpha=1),
                     mpl.lines.Line2D([0], [0], color="#FFB000", lw=2, alpha=0.5)]

        # create the continum compound model (base for plotting lines)
        cont_components = []

        for cmodel in model:
            if isinstance(cmodel, BlackBody1D):
                cont_components.append(cmodel)
                # plot as we go
                ax.plot(x, cmodel(x) * ext_model / x, "#FFB000", alpha=0.5)
        cont_model = cont_components[0]
        for cmodel in cont_components[1:]:
            cont_model += cmodel
        cont_y = cont_model(x)

        # now plot the dust and gas lines
        for cmodel in model:
            if isinstance(cmodel, Gaussian1D):
                ax.plot(x, (cont_y + cmodel(x)) * ext_model / x, "#DC267F", alpha=0.5)
            if isinstance(cmodel, Drude1D):
                ax.plot(x, (cont_y + cmodel(x)) * ext_model / x, "#648FFF", alpha=0.5)

        ax.plot(x, cont_y * ext_model / x, "#785EF0", alpha=1)

        ax.plot(x, model(x) / x, "#FE6100", alpha=1)
        ax.errorbar(x, y / x, yerr=yerr / x,
                    fmt='o', markeredgecolor='k', markerfacecolor='none',
                    ecolor='k', elinewidth=0.2, capsize=0.5, markersize=6)

        ax.set_ylim(0)
        ax.set_ylabel(r"$\nu F_{\nu}$")

        ax.legend(Leg_lines, ["S07_attenuation",
                              "Spectrum Fit",
                              "Dust Features",
                              r"Atomic and $H_2$ Lines",
                              "Total Continuum Emissions",
                              "Continuum Components"], prop={'size': 10}, loc="best", facecolor="white", framealpha=1,
                  ncol=3)
        # residuals
        res = (y - model(x)) / x
        std = np.std(res)
        ax = axs[1]

        ax.set_yscale("linear")
        ax.set_xscale("log")
        ax.tick_params(axis="both", which='major', top="on", right="on", direction='in', length=10)
        ax.tick_params(axis="both", which='minor', top="on", right="on", direction='in', length=5)
        ax.minorticks_on()

        # Custom X axis ticks
        ax.xaxis.set_ticks([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 20, 25, 30, 40])

        ax.axhline(0, linestyle='--', color='gray', zorder=0)
        ax.plot(x, res, 'ko-', fillstyle='none', zorder=1)
        ax.set_ylim(-scalefac_resid * std, scalefac_resid * std)
        ax.set_xlim(np.floor(np.amin(x)), np.ceil(np.amax(x)))
        ax.set_xlabel(r"$\lambda$ [$\mu m$]")
        ax.set_ylabel('Residuals [%]')

        # non scientific x-axis
        ax.xaxis.set_minor_formatter(mpl.ticker.ScalarFormatter())
        ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())

    @staticmethod
    def save(obs_fit, filename, outform):
        """
        Save the model parameters to a user defined file format.

        Parameters
        ----------
        obs_fit : PAHFITBase model
            Model giving all the components and parameters.
        filename : string
            String used to name the output file.
            Currently using the input data file name (without the file's extension).
        outform : string
            Sets the output file format (ascii, fits, csv, etc.).
        """
        # setup the tables for the different components
        bb_table = Table(
            names=(
                "Name",
                "Form",
                "temp",
                "temp_min",
                "temp_max",
                "temp_fixed",
                "amp",
                "amp_min",
                "amp_max",
                "amp_fixed",
            ),
            dtype=(
                "U25",
                "U25",
                "float64",
                "float64",
                "float64",
                "bool",
                "float64",
                "float64",
                "float64",
                "bool",
            ),
        )
        line_table = Table(
            names=(
                "Name",
                "Form",
                "x_0",
                "x_0_min",
                "x_0_max",
                "x_0_fixed",
                "amp",
                "amp_min",
                "amp_max",
                "amp_fixed",
                "fwhm",
                "fwhm_min",
                "fwhm_max",
                "fwhm_fixed",
                "strength",
                "strength_unc",
                "eqw"
            ),
            dtype=(
                "U25",
                "U25",
                "float64",
                "float64",
                "float64",
                "bool",
                "float64",
                "float64",
                "float64",
                "bool",
                "float64",
                "float64",
                "float64",
                "bool",
                "float64",
                "float64",
                "float64"
            ),
        )
        att_table = Table(
            names=("Name", "Form", "amp", "amp_min", "amp_max", "amp_fixed"),
            dtype=("U25", "U25", "float64", "float64", "float64", "bool"),
        )

        for component in obs_fit:
            comp_type = component.__class__.__name__

            if comp_type == "BlackBody1D":
                bb_table.add_row(
                    [
                        component.name,
                        comp_type,
                        component.temperature.value,
                        component.temperature.bounds[0],
                        component.temperature.bounds[1],
                        component.temperature.fixed,
                        component.amplitude.value,
                        component.amplitude.bounds[0],
                        component.amplitude.bounds[1],
                        component.amplitude.fixed,
                    ]
                )
            elif comp_type == "Drude1D":

                # Calculate feature strength.
                strength = pah_feature_strength(component.amplitude.value,
                                                component.fwhm.value,
                                                component.x_0.value)

                strength_unc = None

                # Calculate feature EQW.
                if strength != 0.:
                    eqw = eqws(comp_type,
                               component.x_0.value,
                               component.amplitude.value,
                               component.fwhm,
                               obs_fit)
                else:
                    eqw = 0.

                line_table.add_row(
                    [
                        component.name,
                        comp_type,
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
                        component.fwhm.fixed,
                        strength,
                        strength_unc,
                        eqw
                    ]
                )
            elif comp_type == "Gaussian1D":

                # Calculate feature strength.
                strength = line_strength(component.amplitude.value,
                                         component.mean.value,
                                         component.stddev.value)

                strength_unc = None

                # Calculate feature EQW.
                if strength != 0.:
                    eqw = eqws(comp_type,
                               component.mean.value,
                               component.amplitude.value,
                               component.stddev.value,
                               obs_fit)
                else:
                    eqw = 0.

                line_table.add_row(
                    [
                        component.name,
                        comp_type,
                        component.mean.value,
                        component.mean.bounds[0],
                        component.mean.bounds[1],
                        component.mean.fixed,
                        component.amplitude.value,
                        component.amplitude.bounds[0],
                        component.amplitude.bounds[1],
                        component.amplitude.fixed,
                        2.355 * component.stddev.value,
                        2.355 * component.stddev.bounds[0],
                        2.355 * component.stddev.bounds[1],
                        component.stddev.fixed,
                        strength,
                        strength_unc,
                        eqw
                    ]
                )
            elif comp_type == "S07_attenuation":
                att_table.add_row(
                    [
                        component.name,
                        comp_type,
                        component.tau_sil.value,
                        component.tau_sil.bounds[0],
                        component.tau_sil.bounds[1],
                        component.tau_sil.fixed,
                    ]
                )

        # Call featcombine to calculate combined dust feature strengths.
        cftable = featcombine(line_table)

        # stack the tables (handles missing columns between tables)
        out_table = vstack([bb_table, line_table, att_table, cftable])

        # Writing output table
        out_table.write(
            "{}_output.{}".format(filename, outform), format=outform, overwrite=True
        )

    def read(self, filename, tformat=None):
        """
        Read the model parameters from a file.

        Parameters
        ----------
        filename : string
            The name of the input file containing fit results.

        tformat: string
            table format of filename (compatible with astropy Table.read)

        Returns
        -------
        readout : tuple
            Tuple containing dictionaries of all components from
            the input file.
        """
        # get the table format
        if tformat is None:
            tformat = filename.split(".")[-1]

        # Reading the input file as table
        t = Table.read(filename, format=tformat)

        # Getting indices for the different components
        bb_ind = np.concatenate(np.argwhere(t["Form"] == "BlackBody1D"))
        df_ind = np.concatenate(np.argwhere(t["Form"] == "Drude1D"))
        ga_ind = np.concatenate(np.argwhere(t["Form"] == "Gaussian1D"))
        at_ind = np.concatenate(np.argwhere(t["Form"] == "S07_attenuation"))

        # now split the gas emission lines between H2 and ions
        names = [str(i) for i in np.take(t["Name"], ga_ind)]
        h2_temp = np.concatenate(np.where(np.char.find(names, "H2") >= 0))
        ion_temp = np.concatenate(np.where(np.char.find(names, "H2") == -1))
        h2_ind = np.take(ga_ind, h2_temp)
        ion_ind = np.take(ga_ind, ion_temp)

        # Creating the blackbody dict
        bb_info = {
            "names": np.array(t["Name"][bb_ind].data),
            "temps": np.array(t["temp"][bb_ind].data),
            "temps_limits": _ingest_limits(
                t["temp_min"][bb_ind].data, t["temp_max"][bb_ind].data
            ),
            "temps_fixed": _ingest_fixed(t["temp_fixed"][bb_ind].data),
            "amps": np.array(t["amp"][bb_ind].data),
            "amps_limits": _ingest_limits(
                t["amp_min"][bb_ind].data, t["amp_max"][bb_ind].data
            ),
            "amps_fixed": _ingest_fixed(t["amp_fixed"][bb_ind].data),
        }

        # Creating the dust_features dict
        df_info = {
            "names": np.array(t["Name"][df_ind].data),
            "x_0": np.array(t["x_0"][df_ind].data),
            "x_0_limits": _ingest_limits(
                t["x_0_min"][df_ind].data, t["x_0_max"][df_ind].data
            ),
            "x_0_fixed": _ingest_fixed(t["x_0_fixed"][df_ind].data),
            "amps": np.array(t["amp"][df_ind].data),
            "amps_limits": _ingest_limits(
                t["amp_min"][df_ind].data, t["amp_max"][df_ind].data
            ),
            "amps_fixed": _ingest_fixed(t["amp_fixed"][df_ind].data),
            "fwhms": np.array(t["fwhm"][df_ind].data),
            "fwhms_limits": _ingest_limits(
                t["fwhm_min"][df_ind].data, t["fwhm_max"][df_ind].data
            ),
            "fwhms_fixed": _ingest_fixed(t["fwhm_fixed"][df_ind].data),
        }

        # Creating the H2 dict
        h2_info = {
            "names": np.array(t["Name"][h2_ind].data),
            "x_0": np.array(t["x_0"][h2_ind].data),
            "x_0_limits": _ingest_limits(
                t["x_0_min"][h2_ind].data, t["x_0_max"][h2_ind].data
            ),
            "x_0_fixed": _ingest_fixed(t["x_0_fixed"][h2_ind].data),
            "amps": np.array(t["amp"][h2_ind].data),
            "amps_limits": _ingest_limits(
                t["amp_min"][h2_ind].data, t["amp_max"][h2_ind].data
            ),
            "amps_fixed": _ingest_fixed(t["amp_fixed"][h2_ind].data),
            "fwhms": np.array(t["fwhm"][h2_ind].data),
            "fwhms_limits": _ingest_limits(
                t["fwhm_min"][h2_ind].data, t["fwhm_max"][h2_ind].data
            ),
            "fwhms_fixed": _ingest_fixed(t["fwhm_fixed"][h2_ind].data),
        }

        # Creating the ion dict
        ion_info = {
            "names": np.array(t["Name"][ion_ind].data),
            "x_0": np.array(t["x_0"][ion_ind].data),
            "x_0_limits": _ingest_limits(
                t["x_0_min"][ion_ind].data, t["x_0_max"][ion_ind].data
            ),
            "x_0_fixed": _ingest_fixed(t["x_0_fixed"][ion_ind].data),
            "amps": np.array(t["amp"][ion_ind].data),
            "amps_limits": _ingest_limits(
                t["amp_min"][ion_ind].data, t["amp_max"][ion_ind].data
            ),
            "amps_fixed": _ingest_fixed(t["amp_fixed"][ion_ind].data),
            "fwhms": np.array(t["fwhm"][ion_ind].data),
            "fwhms_limits": _ingest_limits(
                t["fwhm_min"][ion_ind].data, t["fwhm_max"][ion_ind].data
            ),
            "fwhms_fixed": _ingest_fixed(t["fwhm_fixed"][ion_ind].data),
        }

        # Create the attenuation dict
        att_info = {
            "names": np.array(t["Name"][at_ind].data),
            "amps": np.array(t["amp"][at_ind].data),
            "amps_limits": _ingest_limits(
                t["amp_min"][at_ind].data, t["amp_max"][at_ind].data
            ),
            "amps_fixed": _ingest_fixed(t["amp_fixed"][at_ind].data),
        }

        # Create output tuple
        readout = (bb_info, df_info, h2_info, ion_info, att_info)

        return readout

    @staticmethod
    def estimate_init(obs_x, obs_y, param_info):
        """
        Estimate and return updated starting point in param_info based on the input spectrum.

        Parameters
        ----------
        obs_x and obs_y: np.array
            input spectrum

        param_info: tuple of dics
            The dictonaries contain info for each type of component.
        """

        # simple linear interpolation function for spectrum
        sp = interpolate.interp1d(obs_x, obs_y)

        # guess starting point of bb
        for i, (fix, temp) in enumerate(zip(param_info[0]['amps_fixed'], param_info[0]['temps'])):

            if (fix is False) & (temp >= 2500):  # stellar comoponent is defined by BB that has T>=2500 K
                bb = BlackBody1D(1, temp)
                if min(obs_x) < 5:
                    lam = min(obs_x) + 0.1  # the wavelength used to compare
                else:  # if min(obs_x) > 5, use 5.5 um
                    lam = 5.5
                amp_guess = sp(lam) / bb(lam)
                param_info[0]['amps'][i] = amp_guess

            elif fix is False:
                fmax_lam = 2898. / temp
                bb = BlackBody1D(1, temp)
                if (fmax_lam >= min(obs_x)) & (fmax_lam <= max(obs_x)):
                    lam = fmax_lam
                    amp_guess = sp(lam) / bb(lam) * 0.2
                elif (fmax_lam > max(obs_x)):
                    lam = max(obs_x)
                    amp_guess = obs_y[np.argmax(obs_x)] / bb(lam) * 0.2
                else:
                    lam = min(obs_x)
                    amp_guess = obs_y[np.argmin(obs_x)] / bb(lam) * 0.2
                param_info[0]['amps'][i] = amp_guess
            else:
                pass

        # guess starting point of dust features and lines
        # set to half of the median (non-negative) intensity of the entire input spectrum

        # dust
        for i, fix in enumerate(param_info[1]['amps_fixed']):

            if fix is False:
                amp_guess = 0.5 * np.median(obs_y)
                param_info[1]['amps'][i] = amp_guess

        # h2
        for i, fix in enumerate(param_info[2]['amps_fixed']):

            if fix is False:
                amp_guess = 0.5 * np.median(obs_y)
                param_info[2]['amps'][i] = amp_guess

        # ion
        for i, fix in enumerate(param_info[3]['amps_fixed']):

            if fix is False:
                amp_guess = 0.5 * np.median(obs_y)
                param_info[3]['amps'][i] = amp_guess

        return param_info
