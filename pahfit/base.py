from astropy.modeling.functional_models import Gaussian1D

from pahfit.component_models import (
    BlackBody1D,
    ModifiedBlackBody1D,
    S07_attenuation,
    att_Drude1D,
)
from astropy.table import Table, vstack
from astropy.modeling.physical_models import Drude1D

from scipy import interpolate

import numpy as np

import matplotlib as mpl

from pahfit.instrument import within_segment, fwhm_recommendation
from pahfit.feature_strengths import (
    pah_feature_strength,
    line_strength,
    featcombine,
    eqws,
)

from pahfit.features import Features

__all__ = ["PAHFITBase"]


def _ingest_limits(min_vals, max_vals):
    """
    Ingest the limits read from yaml file and generate the appropriate
    internal format (list of tuples).  Needed to handle the case
    where a limit is not desired as numpy arrays cannot have elements
    of None, instead a value of nan is used.

    Limits that are not set are designated as 'nan' in files and
    these are changed to the python None to be compatible with
    the astropy.modeling convention.

    Parameters
    ----------
    min_vals,
    max_vals : numpy.array (masked arrays)
        min/max values of the limits for a parameter
        nan designates no limit

    Returns
    -------
    plimits : list of tuples
        tuples give the min/max limits for a parameter
    """
    plimits = []
    mask_min = min_vals.mask
    data_min = min_vals.data
    mask_max = max_vals.mask
    data_max = max_vals.data

    mask_min_ind = np.where(mask_min == False)[0]
    mask_max_ind = np.where(mask_max == False)[0]

    min_vals = np.zeros(len(mask_min))
    min_vals[mask_min_ind] = data_min[mask_min_ind]

    max_vals = np.zeros(len(mask_max))
    max_vals[mask_max_ind] = data_max[mask_max_ind]

    plimits = []
    for cmin, cmax in zip(min_vals, max_vals):
        if np.isinf(cmin):
            cmin = None
        if np.isinf(cmax):
            cmax = None
        plimits.append((cmin, cmax))

    return plimits


def _ingest_fixed(fixed_vals):
    """
    Ingest the fixed value read from a file and generate the appropriate
    internal format (list of booleans). Since this information is indirectly
    hidden in the parameter of a feature, this function is needed to 
    extract that. 

    Parameters
    ----------
    min_vals : numpy.array (masked array)
        fixed designations

    Returns
    -------
    pfixed : list (boolean)
        True/False designation for parameters
    """

    check_mask = fixed_vals.mask
    mask_false_ind = np.where(check_mask == False)[0]
    fixed_vals = ["True"] * len(check_mask)
    for i in range(0, len(mask_false_ind)):
        ind = mask_false_ind[i]
        fixed_vals[ind] = "False"

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

    instrumentname: string
        the instrument which which the input spectrum
        was observed.    

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
    def __init__(
        self,
        obs_x,
        obs_y,
        instrumentname,
        estimate_start=False,
        param_info=None,
        filename=None,
        tformat=None,
    ):
        """
        Setup a variant based on inputs.  Generates an astropy.modeling
        compound model.
        """
        # check that param_info or filename is set
        if filename is None and param_info is None:
            raise ValueError("Either param_info or filename need to be set \
                             when initializing a PAHFITBase object")

        # read in the parameter info from a file
        if filename is not None:
            param_info = self.read(filename, instrumentname, tformat=tformat)

        if estimate_start:
            # guess values and update starting point (if not set fixed) based on the input spectrum
            param_info = self.estimate_init(obs_x, obs_y, param_info)

        if not param_info:
            raise ValueError("No parameter information set.")

        self.param_info = param_info

        bb_info = param_info[0]
        dust_features = param_info[1]
        h2_features = param_info[2]
        ion_features = param_info[3]
        att_info = param_info[4]

        # setup the model
        self.model = None
        self.bb_info = bb_info
        if bb_info is not None:
            bbs = []
            for k in range(len(bb_info["names"])):
                BBClass = ModifiedBlackBody1D if bb_info["modified"][
                    k] else BlackBody1D
                bbs.append(
                    BBClass(
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
                    ))
            self.model = sum(bbs[1:], bbs[0])

        self.dust_features = dust_features
        if dust_features is not None:
            df = []
            for k in range(len(dust_features["names"])):
                df.append(
                    Drude1D(
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
                    ))

            df = sum(df[1:], df[0])
            if self.model:
                self.model += df
            else:
                self.model = df

        self.h2_features = h2_features
        if h2_features is not None:
            h2 = []
            for k in range(len(h2_features["names"])):
                h2.append(
                    Gaussian1D(
                        name=h2_features["names"][k],
                        amplitude=h2_features["amps"][k],
                        mean=h2_features["x_0"][k],
                        stddev=h2_features["fwhms"][k] / 2.355,
                        bounds={
                            "amplitude":
                            h2_features["amps_limits"][k],
                            "mean":
                            h2_features["x_0_limits"][k],
                            "stddev": (
                                h2_features["fwhms"][k] * 0.9 / 2.355,
                                h2_features["fwhms"][k] * 1.1 / 2.355,
                            ),
                        },
                        fixed={
                            "amplitude": h2_features["amps_fixed"][k],
                            "mean": h2_features["x_0_fixed"][k],
                            "stddev": h2_features["fwhms_fixed"][k],
                        },
                    ))
            h2 = sum(h2[1:], h2[0])
            if self.model:
                self.model += h2
            else:
                self.model = h2

        self.ion_features = ion_features
        if ion_features is not None:
            ions = []
            for k in range(len(ion_features["names"])):
                ions.append(
                    Gaussian1D(
                        name=ion_features["names"][k],
                        amplitude=ion_features["amps"][k],
                        mean=ion_features["x_0"][k],
                        stddev=ion_features["fwhms"][k] / 2.355,
                        bounds={
                            "amplitude":
                            ion_features["amps_limits"][k],
                            "mean":
                            ion_features["x_0_limits"][k],
                            "stddev": (
                                ion_features["fwhms"][k] * 0.9 / 2.355,
                                ion_features["fwhms"][k] * 1.1 / 2.355,
                            ),
                        },
                        fixed={
                            "amplitude": ion_features["amps_fixed"][k],
                            "mean": ion_features["x_0_fixed"][k],
                            "stddev": ion_features["fwhms_fixed"][k],
                        },
                    ))
            ions = sum(ions[1:], ions[0])
            if self.model:
                self.model += ions
            else:
                self.model = ions

        # add additional att components to the model if necessary
        if not self.model:
            raise ValueError("No model components found")

        self.att_info = att_info
        if att_info is not None:
            for k in range(len(att_info["names"])):
                if (
                        att_info["names"][k] == "S07_att"
                ):  # Only loop through att components that can be parameterized
                    self.model *= S07_attenuation(
                        name=att_info["names"][k],
                        tau_sil=att_info["amps"][k],
                        bounds={"tau_sil": att_info["amps_limits"][k]},
                        fixed={"tau_sil": att_info["amps_fixed"][k]},
                    )
                else:
                    self.model *= att_Drude1D(
                        name=att_info["names"][k],
                        tau=att_info["amps"][k],
                        x_0=att_info["x_0"][k],
                        fwhm=att_info["fwhms"][k],
                        bounds={
                            "tau": att_info["amps_limits"][k],
                            "fwhm": att_info["fwhms_limits"][k],
                        },
                        fixed={"x_0": att_info["x_0_fixed"][k]},
                    )

    @staticmethod
    def plot(axs, x, y, yerr, model, model_samples=1000, scalefac_resid=2):
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
        model_samples : int
            Total number of wavelength points to allocate to the model display
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

        # Fine x samples for model fit
        x_mod = np.logspace(np.log10(min(x)), np.log10(max(x)), model_samples)

        # spectrum and best fit model
        ax = axs[0]
        ax.set_yscale("linear")
        ax.set_xscale("log")
        ax.minorticks_on()
        ax.tick_params(axis="both",
                       which="major",
                       top="on",
                       right="on",
                       direction="in",
                       length=10)
        ax.tick_params(axis="both",
                       which="minor",
                       top="on",
                       right="on",
                       direction="in",
                       length=5)

        ax_att = ax.twinx()  # axis for plotting the extinction curve
        ax_att.tick_params(which="minor", direction="in", length=5)
        ax_att.tick_params(which="major", direction="in", length=10)
        ax_att.minorticks_on()

        # get the extinction model (probably a better way to do this)
        ext_model = None
        for cmodel in model:
            if isinstance(cmodel, S07_attenuation):
                ext_model = cmodel(x_mod)

        # get additional extinction components that can be
        # characterized by functional forms (Drude profile in this case)
        for cmodel in model:
            if isinstance(cmodel, att_Drude1D):
                if ext_model is not None:
                    ext_model *= cmodel(x_mod)
                else:
                    ext_model = cmodel(x_mod)
        ax_att.plot(x_mod, ext_model, "k--", alpha=0.5)
        ax_att.set_ylabel("Attenuation")
        ax_att.set_ylim(0, 1.1)

        # Define legend lines
        Leg_lines = [
            mpl.lines.Line2D([0], [0], color="k", linestyle="--", lw=2),
            mpl.lines.Line2D([0], [0], color="#FE6100", lw=2),
            mpl.lines.Line2D([0], [0], color="#648FFF", lw=2, alpha=0.5),
            mpl.lines.Line2D([0], [0], color="#DC267F", lw=2, alpha=0.5),
            mpl.lines.Line2D([0], [0], color="#785EF0", lw=2, alpha=1),
            mpl.lines.Line2D([0], [0], color="#FFB000", lw=2, alpha=0.5),
        ]

        # create the continum compound model (base for plotting lines)
        cont_components = []

        for cmodel in model:
            if isinstance(cmodel, BlackBody1D):
                cont_components.append(cmodel)
                # plot as we go
                ax.plot(x_mod,
                        cmodel(x_mod) * ext_model / x_mod,
                        "#FFB000",
                        alpha=0.5)
        cont_model = cont_components[0]
        for cmodel in cont_components[1:]:
            cont_model += cmodel
        cont_y = cont_model(x_mod)

        # now plot the dust bands and lines
        for cmodel in model:
            if isinstance(cmodel, Gaussian1D):
                ax.plot(
                    x_mod,
                    (cont_y + cmodel(x_mod)) * ext_model / x_mod,
                    "#DC267F",
                    alpha=0.5,
                )
            if isinstance(cmodel, Drude1D):
                ax.plot(
                    x_mod,
                    (cont_y + cmodel(x_mod)) * ext_model / x_mod,
                    "#648FFF",
                    alpha=0.5,
                )

        ax.plot(x_mod, cont_y * ext_model / x_mod, "#785EF0", alpha=1)

        ax.plot(x_mod, model(x_mod) / x_mod, "#FE6100", alpha=1)
        ax.errorbar(
            x,
            y / x,
            yerr=yerr / x,
            fmt="o",
            markeredgecolor="k",
            markerfacecolor="none",
            ecolor="k",
            elinewidth=0.2,
            capsize=0.5,
            markersize=6,
        )

        ax.set_ylim(0)
        ax.set_ylabel(r"$\nu F_{\nu}$")

        ax.legend(
            Leg_lines,
            [
                "S07_attenuation",
                "Spectrum Fit",
                "Dust Features",
                r"Atomic and $H_2$ Lines",
                "Total Continuum Emissions",
                "Continuum Components",
            ],
            prop={"size": 10},
            loc="best",
            facecolor="white",
            framealpha=1,
            ncol=3,
        )

        # residuals, lower sub-figure
        res = (y - model(x)) / x
        std = np.std(res)
        ax = axs[1]

        ax.set_yscale("linear")
        ax.set_xscale("log")
        ax.tick_params(axis="both",
                       which="major",
                       top="on",
                       right="on",
                       direction="in",
                       length=10)
        ax.tick_params(axis="both",
                       which="minor",
                       top="on",
                       right="on",
                       direction="in",
                       length=5)
        ax.minorticks_on()

        # Custom X axis ticks
        ax.xaxis.set_ticks(
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 20, 25, 30, 40])

        ax.axhline(0, linestyle="--", color="gray", zorder=0)
        ax.plot(x, res, "ko-", fillstyle="none", zorder=1)
        ax.set_ylim(-scalefac_resid * std, scalefac_resid * std)
        ax.set_xlim(np.floor(np.amin(x)), np.ceil(np.amax(x)))
        ax.set_xlabel(r"$\lambda$ [$\mu m$]")
        ax.set_ylabel("Residuals [%]")

        # scalar x-axis marks
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
                "eqw",
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
                "float64",
            ),
        )
        att_table = Table(
            names=("Name", "Form", "amp", "amp_min", "amp_max", "amp_fixed"),
            dtype=("U25", "U25", "float64", "float64", "float64", "bool"),
        )

        # attenuation components that can be characterized by a functional form
        att_funct_table = Table(
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
            ),
        )

        for component in obs_fit:
            comp_type = component.__class__.__name__

            if isinstance(component, BlackBody1D):
                bb_table.add_row([
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
                ])
            elif isinstance(component, Drude1D):

                # Calculate feature strength.
                strength = pah_feature_strength(component.amplitude.value,
                                                component.fwhm.value,
                                                component.x_0.value)

                strength_unc = None

                # Calculate feature EQW.
                if strength != 0.0:
                    eqw = eqws(
                        comp_type,
                        component.x_0.value,
                        component.amplitude.value,
                        component.fwhm,
                        obs_fit,
                    )
                else:
                    eqw = 0.0

                line_table.add_row([
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
                    eqw,
                ])
            elif isinstance(component, Gaussian1D):

                # Calculate feature strength.
                strength = line_strength(
                    component.amplitude.value,
                    component.mean.value,
                    component.stddev.value,
                )

                strength_unc = None

                # Calculate feature EQW.
                if strength != 0.0:
                    eqw = eqws(
                        comp_type,
                        component.mean.value,
                        component.amplitude.value,
                        component.stddev.value,
                        obs_fit,
                    )
                else:
                    eqw = 0.0

                line_table.add_row([
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
                    eqw,
                ])
            elif isinstance(component, S07_attenuation):
                att_table.add_row([
                    component.name,
                    comp_type,
                    component.tau_sil.value,
                    component.tau_sil.bounds[0],
                    component.tau_sil.bounds[1],
                    component.tau_sil.fixed,
                ])

            elif isinstance(component, att_Drude1D):
                att_funct_table.add_row([
                    component.name,
                    comp_type,
                    component.x_0.value,
                    component.x_0.bounds[0],
                    component.x_0.bounds[1],
                    component.x_0.fixed,
                    component.tau.value,
                    component.tau.bounds[0],
                    component.tau.bounds[1],
                    component.tau.fixed,
                    component.fwhm.value,
                    component.fwhm.bounds[0],
                    component.fwhm.bounds[1],
                    component.fwhm.fixed,
                ])

        # Call featcombine to calculate combined dust feature strengths.
        cftable = featcombine(line_table)

        # stack the tables (handles missing columns between tables)
        out_table = vstack(
            [bb_table, line_table, att_table, att_funct_table, cftable])

        # Writing output table
        out_table.write("{}_output.{}".format(filename, outform),
                        format=outform,
                        overwrite=True)

    @staticmethod
    def update_dictionary(feature_dict, instrumentname, update_fwhms=False):
        """
        Update parameter dictionary based on the instrument name.
        Based on the instrument name, this function removes the 
        features outside of the wavelength range and
        updates the FWHMs of the lines. 


        Parameters
        ----------
        feature_dict : dictionary
            Dictionary created by reading in a science pack.
            
        instrumentname : string
            Name of the instrument with which the input spectrum
            is observed.
            
        update_fwhms = Boolean
            True for h2_info and ion_info
            False for df_info

        Returns
        -------
        updated feature_dict
        """

        ind = np.nonzero(within_segment(feature_dict["x_0"], instrumentname))[0]

        # select the valid entries in these arrays
        array_keys = ("x_0", "amps", "fwhms", "names")
        new_values_1 = {key: feature_dict[key][ind] for key in array_keys}

        # these are lists instead
        list_keys = (
            "amps_fixed",
            "fwhms_fixed",
            "x_0_fixed",
            "x_0_limits",
            "amps_limits",
            "fwhms_limits",
        )
        new_values_2 = {
            key: [feature_dict[key][i] for i in ind] for key in list_keys
        }

        feature_dict.update(new_values_1)
        feature_dict.update(new_values_2)

        if update_fwhms:
            waves = feature_dict["x_0"]
            values, fixed, mins, maxes = fwhm_recommendation(instrumentname, waves)
            feature_dict.update(
                {
                    "fwhms": values,
                    "fwhms_fixed": fixed,
                    "fwhms_limits": list(zip(mins, maxes)),
                }
            )

        return feature_dict

    @staticmethod
    def parse_table(pack_table):
        """
        Load the model parameters from a Table

        Parameters
        ----------
        pack_table : Table
            Table created by reading in a science pack.

        Returns
        -------
        readout : tuple
            Tuple containing dictionaries of all components from the
            input Table. Can be used to create PAHFITBase instance using
            param_info argument. Dictionary in tuple is None if no
            components of that type were specified.
        """
        # Getting indices for the different components
        bb_ind = np.where((pack_table["kind"] == "starlight")
                          | (pack_table["kind"] == "dust_continuum"))[0]
        df_ind = np.where(pack_table["kind"] == "dust_feature")[0]
        ga_ind = np.where(pack_table["kind"] == "line")[0]
        at_ind = np.where((pack_table["kind"] == "attenuation")
                          | (pack_table["kind"] == "absorption"))[0]

        # now split the gas emission lines between H2 and ions
        names = [str(i) for i in pack_table["name"][ga_ind]]
        if len(names) > 0:
            # this has trouble with empty list
            h2_temp = np.char.find(names, "H2") >= 0
            ion_temp = np.char.find(names, "H2") == -1
            h2_ind = ga_ind[h2_temp]
            ion_ind = ga_ind[ion_temp]
        else:
            h2_ind = []
            ion_ind = []
        # the rest works fine with empty list

        # Creating the blackbody dict
        bb_info = None
        if len(bb_ind) > 0:
            bb_info = {
                "names":
                np.array(pack_table["name"][bb_ind].data),
                "temps":
                np.array(pack_table["temperature"][:, 0][bb_ind].data),
                "temps_limits":
                _ingest_limits(
                    pack_table["temperature"][:, 1][bb_ind],
                    pack_table["temperature"][:, 2][bb_ind],
                ),
                "temps_fixed":
                _ingest_fixed(pack_table["temperature"][:, 1][bb_ind]),
                "amps":
                np.array(pack_table["tau"][:, 0][bb_ind].data),
                "amps_limits":
                _ingest_limits(
                    pack_table["tau"][:, 1][bb_ind],
                    pack_table["tau"][:, 2][bb_ind],
                ),
                "amps_fixed":
                _ingest_fixed(pack_table["tau"][:, 1][bb_ind]),
                "modified":
                np.array(pack_table["kind"][bb_ind] == "dust_continuum"),
            }

        # Creating the dust_features dict
        df_info = None
        if len(df_ind) > 0:
            df_info = {
                "names":
                np.array(pack_table["name"][df_ind].data),
                "x_0":
                np.array(pack_table["wavelength"][:, 0][df_ind].data),
                "x_0_limits":
                _ingest_limits(
                    pack_table["wavelength"][:, 1][df_ind],
                    pack_table["wavelength"][:, 2][df_ind],
                ),
                "x_0_fixed":
                _ingest_fixed(pack_table["wavelength"][:, 1][df_ind]),
                "amps":
                np.array(pack_table["power"][:, 0][df_ind].data),
                "amps_limits":
                _ingest_limits(
                    pack_table["power"][:, 1][df_ind],
                    pack_table["power"][:, 2][df_ind],
                ),
                "amps_fixed":
                _ingest_fixed(pack_table["power"][:, 1][df_ind]),
                "fwhms":
                np.array(pack_table["fwhm"][:, 0][df_ind].data),
                "fwhms_limits":
                _ingest_limits(
                    pack_table["fwhm"][:, 1][df_ind],
                    pack_table["fwhm"][:, 2][df_ind],
                ),
                "fwhms_fixed":
                _ingest_fixed(pack_table["fwhm"][:, 1][df_ind]),
            }

        # Creating the H2 dict
        h2_info = None
        if len(h2_ind) > 0:
            h2_info = {
                "names":
                np.array(pack_table["name"][h2_ind].data),
                "x_0":
                np.array(pack_table["wavelength"][:, 0][h2_ind].data),
                "x_0_limits":
                _ingest_limits(
                    pack_table["wavelength"][:, 1][h2_ind],
                    pack_table["wavelength"][:, 2][h2_ind],
                ),
                "x_0_fixed":
                _ingest_fixed(pack_table["wavelength"][:, 1][h2_ind]),
                "amps":
                np.array(pack_table["power"][:, 0][h2_ind].data),
                "amps_limits":
                _ingest_limits(
                    pack_table["power"][:, 1][h2_ind],
                    pack_table["power"][:, 2][h2_ind],
                ),
                "amps_fixed":
                _ingest_fixed(pack_table["power"][:, 1][h2_ind]),
                "fwhms":
                np.array(pack_table["fwhm"][:, 0][h2_ind].data),
                "fwhms_limits":
                _ingest_limits(
                    pack_table["fwhm"][:, 1][h2_ind],
                    pack_table["fwhm"][:, 2][h2_ind],
                ),
                "fwhms_fixed":
                _ingest_fixed(pack_table["fwhm"][:, 1][h2_ind]),
            }

        # Creating the ion dict
        ion_info = None
        if len(ion_ind) > 0:
            ion_info = {
                "names":
                np.array(pack_table["name"][ion_ind].data),
                "x_0":
                np.array(pack_table["wavelength"][:, 0][ion_ind].data),
                "x_0_limits":
                _ingest_limits(
                    pack_table["wavelength"][:, 1][ion_ind],
                    pack_table["wavelength"][:, 2][ion_ind],
                ),
                "x_0_fixed":
                _ingest_fixed(pack_table["wavelength"][:, 1][ion_ind]),
                "amps":
                np.array(pack_table["power"][:, 0][ion_ind].data),
                "amps_limits":
                _ingest_limits(
                    pack_table["power"][:, 1][ion_ind],
                    pack_table["power"][:, 2][ion_ind],
                ),
                "amps_fixed":
                _ingest_fixed(pack_table["power"][:, 1][ion_ind]),
                "fwhms":
                np.array(pack_table["fwhm"][:, 0][ion_ind].data),
                "fwhms_limits":
                _ingest_limits(
                    pack_table["fwhm"][:, 1][ion_ind].data,
                    pack_table["fwhm"][:, 2][ion_ind].data,
                ),
                "fwhms_fixed":
                _ingest_fixed(pack_table["fwhm"][:, 1][ion_ind].data),
            }

        # Create the attenuation dict
        att_info = None
        if len(at_ind) > 0:
            att_info = {
                "names":
                np.array(pack_table["name"][at_ind].data),
                "x_0":
                np.array(pack_table["wavelength"][:, 0][at_ind].data),
                "x_0_limits":
                _ingest_limits(
                    pack_table["wavelength"][:, 1][at_ind],
                    pack_table["wavelength"][:, 2][at_ind],
                ),
                "x_0_fixed":
                _ingest_fixed(pack_table["wavelength"][:, 1][at_ind]),
                "amps":
                np.array(pack_table["tau"][:, 0][at_ind].data),
                "amps_limits":
                _ingest_limits(
                    pack_table["tau"][:, 0][at_ind],
                    pack_table["tau"][:, 1][at_ind],
                ),
                "amps_fixed":
                _ingest_fixed(pack_table["tau"][:, 1][at_ind]),
                "fwhms":
                np.array(pack_table["fwhm"][:, 0][at_ind].data),
                "fwhms_limits":
                _ingest_limits(
                    pack_table["fwhm"][:, 1][at_ind],
                    pack_table["fwhm"][:, 2][at_ind],
                ),
                "fwhms_fixed":
                _ingest_fixed(pack_table["fwhm"][:, 1][at_ind]),
            }

        if att_info["names"] == 'silicate':
            name_array = np.array(['S07_att'])
            att_info.update({'names': name_array})

        return (bb_info, df_info, h2_info, ion_info, att_info)

    @staticmethod
    def read(filename, instrumentname, tformat=None):
        """
        Create model by reading the parameters from a file.

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
        t = Features.read(filename)
        bb_info, df_info, h2_info, ion_info, att_info = PAHFITBase.parse_table(
            t)
        df_info = PAHFITBase.update_dictionary(df_info, instrumentname)
        h2_info = PAHFITBase.update_dictionary(h2_info,
                                               instrumentname,
                                               update_fwhms=True)
        ion_info = PAHFITBase.update_dictionary(ion_info,
                                                instrumentname,
                                                update_fwhms=True)

        return (bb_info, df_info, h2_info, ion_info, att_info)

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
        for i, (fix, temp) in enumerate(
                zip(param_info[0]["amps_fixed"], param_info[0]["temps"])):

            if (fix is False) & (
                    temp >= 2500
            ):  # stellar comoponent is defined by BB that has T>=2500 K
                bb = BlackBody1D(1, temp)
                if min(obs_x) < 5:
                    lam = min(obs_x) + 0.1  # the wavelength used to compare
                else:  # if min(obs_x) > 5, use 5.5 um
                    lam = 5.5
                amp_guess = sp(lam) / bb(lam)
                param_info[0]["amps"][i] = amp_guess

            elif fix is False:
                fmax_lam = 2898.0 / temp
                bb = BlackBody1D(1, temp)
                if (fmax_lam >= min(obs_x)) & (fmax_lam <= max(obs_x)):
                    lam = fmax_lam
                    amp_guess = sp(lam) / bb(lam) * 0.2
                elif fmax_lam > max(obs_x):
                    lam = max(obs_x)
                    amp_guess = obs_y[np.argmax(obs_x)] / bb(lam) * 0.2
                else:
                    lam = min(obs_x)
                    amp_guess = obs_y[np.argmin(obs_x)] / bb(lam) * 0.2
                param_info[0]["amps"][i] = amp_guess
            else:
                pass

        # guess starting point of dust features and lines
        # set to half of the median (non-negative) intensity of the entire input spectrum

        # dust (1), h2 (2), and ion (3)
        for j in range(1, 4):
            if param_info[j] is not None:
                for i, fix in enumerate(param_info[j]["amps_fixed"]):
                    if fix is False:
                        amp_guess = 0.5 * np.median(obs_y)
                        param_info[j]["amps"][i] = amp_guess

        return param_info
