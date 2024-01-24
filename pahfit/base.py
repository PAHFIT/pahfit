import numpy as np
import matplotlib as mpl

from pahfit.instrument import within_segment, fwhm
from pahfit.errors import PAHFITModelError
from pahfit.component_models import (
    BlackBody1D,
    ModifiedBlackBody1D,
    S07_attenuation,
    att_Drude1D,
)
from astropy.modeling.physical_models import Drude1D
from astropy.modeling.functional_models import Gaussian1D

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

    mask_min_ind = np.where(np.logical_not(mask_min))[0]
    mask_max_ind = np.where(np.logical_not(mask_max))[0]

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
    mask_false_ind = np.where(np.logical_not(fixed_vals.mask))[0]
    fixed_vals = ["True"] * len(fixed_vals.mask)
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
    Old implementation. Some functions are still used by the new Model
    class. The unused functionality has been removed.

    Construct that is still used for now

    param_info: tuple of dicts called (bb_info, df_info, h2_info, ion_info, abs_info, att_info)
        The dictionaries contain info for each type of component. Each
        component of the dictonaries is a vector.
        bb_info -
        dict with {name, temps, temps_limits, temps_fixed,
        amps, amps_limits, amps_fixed}, each a vector
        dust_features, h2_features, ion_features -
        dict with {name amps, amps_limits, amps_fixed,
        x_0, x_0_limits, x_0_fixed, fwhms, fwhms_limits, fwhm_fixed}.
    """
    @staticmethod
    def model_from_param_info(param_info):
        # setup the model
        bb_info = param_info[0]
        dust_features = param_info[1]
        h2_features = param_info[2]
        ion_features = param_info[3]
        abs_info = param_info[4]
        att_info = param_info[5]

        model = None
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
                    )
                )
            model = sum(bbs[1:], bbs[0])

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
            if model:
                model += df
            else:
                model = df

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
            if model:
                model += h2
            else:
                model = h2

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
            if model:
                model += ions
            else:
                model = ions

        # add additional att components to the model if necessary
        if not model:
            raise PAHFITModelError("No model components found")

        if abs_info is not None:
            for k in range(len(abs_info["names"])):
                model *= att_Drude1D(
                    name=abs_info["names"][k],
                    tau=abs_info["amps"][k],
                    x_0=abs_info["x_0"][k],
                    fwhm=abs_info["fwhms"][k],
                    bounds={
                        "tau": abs_info["amps_limits"][k],
                        "fwhm": abs_info["fwhms_limits"][k],
                    },
                    fixed={"x_0": abs_info["x_0_fixed"][k]},
                )

        if att_info is not None:
            model *= S07_attenuation(
                name=att_info["name"],
                tau_sil=att_info["tau_sil"],
                bounds={"tau_sil": att_info["tau_sil_limits"]},
                fixed={"tau_sil": att_info["tau_sil_fixed"]},
            )

        return model

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
    def update_dictionary(feature_dict, instrumentname, update_fwhms=False, redshift=0):
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
        if feature_dict is None:
            return None

        # convert from physical feature, to observed wavelength
        def redshifted_waves():
            return feature_dict["x_0"] * (1 + redshift)

        ind = np.nonzero(within_segment(redshifted_waves(), instrumentname))[0]

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

        if len(feature_dict['names']) == 0:
            # if we removed all the things, be careful here. Setting to
            # None should make the model construction function behave
            # normally.
            feature_dict = None
            return feature_dict

        if update_fwhms:
            # observe the lines at the redshifted wavelength
            fwhm_min_max = fwhm(instrumentname, redshifted_waves(), as_bounded=True)
            # shift the observed fwhm back to the rest frame (where the
            # observed data will be moved, and its features will become
            # narrower)
            fwhm_min_max /= (1 + redshift)
            # For astropy a numpy.bool does not work for the 'fixed'
            # parameter. It needs to be a regular bool. Doing tolist()
            # instead of using the array mask directly solves this.
            feature_dict.update(
                {
                    "fwhms": fwhm_min_max[:, 0],
                    # masked means there is no min/max, i.e. they need to be fixed
                    "fwhms_fixed": fwhm_min_max[:, 1].mask.tolist(),
                    "fwhms_limits": fwhm_min_max[:, 1:].tolist(),
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
        at_ind = np.where(pack_table["kind"] == "attenuation")[0]
        ab_ind = np.where(pack_table["kind"] == "absorption")[0]

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

        # Create the attenuation dict (could be absorption drudes
        # and S07 model)
        abs_info = None
        if len(ab_ind) > 0:
            abs_info = {
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

        att_info = None
        if len(at_ind) > 1:
            raise NotImplementedError("More than one attenuation component not supported")
        elif len(at_ind) == 1:
            i = at_ind[0]
            att_info = {"name": pack_table["name"][i],
                        "tau_sil": pack_table["tau"][i][0],
                        "tau_sil_limits": pack_table["tau"][i][1:],
                        "tau_sil_fixed": True if pack_table["tau"][i].mask[1] else False}

        return [bb_info, df_info, h2_info, ion_info, abs_info, att_info]
