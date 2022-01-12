import os
import pkg_resources

import astropy.units as u
from astropy.table import Table

from astropy.modeling.fitting import LevMarLSQFitter

from pahfit.base import PAHFITBase

from pahfit.component_models import BlackBody1D, S07_attenuation
from astropy.modeling.physical_models import Drude1D
from astropy.modeling.functional_models import Gaussian1D


__all__ = ["read_spectrum", "initialize_model", "fit_spectrum", "calculate_compounds"]


def read_spectrum(specfile, colnames=["wavelength", "flux", "sigma"]):
    """
    Read in a spectrum and convert intput units to the expected internal PAHFIT units.

    Parameters
    ----------
    specfile : string
        file with the spectrum to be fit

    colnames : list of strings
        list giving the column names of the wavelength, flux, and flux uncertainty
        in the spectrum file with default =  ["wavelength", "flux", "sigma"]

    Returns
    -------
    obsdata : dict
        x is wavelength in microns and y/unc are the spectrum/unc in units
        of Jy
    """
    # read in the observed spectrum
    # assumed to be astropy table compatibile and include units
    if not os.path.isfile(specfile):
        pack_path = pkg_resources.resource_filename("pahfit", "data/")
        test_specfile = "{}/{}".format(pack_path, specfile)
        if os.path.isfile(test_specfile):
            specfile = test_specfile
        else:
            raise ValueError("Input spectrumfile {} not found".format(specfile))

    # get the table format (from extension of filename)
    tformat = specfile.split(".")[-1]
    if tformat == "ecsv":
        tformat = "ascii.ecsv"
    obs_spectrum = Table.read(specfile, format=tformat)
    obsdata = {}
    obsdata["x"] = obs_spectrum[colnames[0]].to(u.micron, equivalencies=u.spectral())
    obsdata["y"] = obs_spectrum[colnames[1]].to(
        u.Jy, equivalencies=u.spectral_density(obsdata["x"])
    )
    obsdata["unc"] = obs_spectrum[colnames[2]].to(
        u.Jy, equivalencies=u.spectral_density(obsdata["x"])
    )

    return obsdata


def initialize_model(packfile, obsdata, estimate_start=False):
    """
    Initialize a model based on the packfile

    Parameters
    ----------
    packfile : string
        file with the PAHFIT pack information

    obsdata : dict
        observed data where x = wavelength, y = SED, and unc = uncertainties

    estimate_start : boolean
        estimate the starting parameters based on the observed data

    Returns
    -------
    pmodel : PAHFITBase model
        PAHFIT model
    """

    # read in the pack file
    if not os.path.isfile(packfile):
        pack_path = pkg_resources.resource_filename("pahfit", "packs/")
        test_packfile = "{}/{}".format(pack_path, packfile)
        if os.path.isfile(test_packfile):
            packfile = test_packfile
        else:
            raise ValueError("Input packfile {} not found".format(packfile))

    pmodel = PAHFITBase(
        obsdata["x"].value,
        obsdata["y"].value,
        estimate_start=estimate_start,
        filename=packfile,
    )

    return pmodel


def fit_spectrum(obsdata, pmodel, maxiter=1000, verbose=True):
    """
    Fit the observed data using the input PAHFIT model.

    Parameters
    ----------
    obsdata : dict
        observed data where x = wavelength, y = SED, and unc = uncertainties

    pmodel : PAHFITBase model
        PAHFIT model

    maxiter : int
        maximum number of fitting iterations

    verbose : boolean
        set to provide screen output

    Returns
    -------
    obsfit : PAHFITBase model (astropy modeling CompoundModel)
        PAHFIT model with best fit parameters
    """

    # pick the fitter
    fit = LevMarLSQFitter()

    # fit
    obs_fit = fit(
        pmodel.model,
        obsdata["x"].value,
        obsdata["y"].value,
        weights=1.0 / obsdata["unc"].value,
        maxiter=maxiter,
        epsilon=1e-10,
        acc=1e-10,
    )
    if verbose:
        print(fit.fit_info["message"])

    return obs_fit


def calculate_compounds(obsdata, pmodel):
    """
    Determine model compounds for total continuum, stellar continuum,
    total dust continuum, combined dust features,
    combined atomic and H2 lines, combined H2 lines,
    combined atomic lines, and extinction model

    Parameters
    ----------
    obsdata : dict
        observed data where x = wavelength, y = SED, and unc = uncertainties

    pmodel : PAHFITBase model
        model giving all the components and parameters

    Returns
    -------
    compounds : dict
        x = wavelength in microns;
        tot_cont = total continuum;
        stellar_cont = stellar continuum;
        dust_cont = total dust continuum;
        dust_features = combined dust features;
        tot_lines = combined atomic and H2 emission lines;
        h2_lines = combined H2 lines;
        atomic_lines = combined atomic lines;
        extinction_model = extinction model
    """

    # get wavelength array
    x = obsdata["x"].value

    # calculate total dust continuum and total continuum (including stellar continuum)
    # v2.0: first BlackBody1D is stellar continuum
    cont_components = []

    for cmodel in pmodel.model:
        if isinstance(cmodel, BlackBody1D):
            cont_components.append(cmodel)
    stellar_cont_model = cont_components[0]
    dust_cont_model = cont_components[1]
    for cmodel in cont_components[2:]:
        dust_cont_model += cmodel
    totcont = dust_cont_model(x) + stellar_cont_model(x)

    # calculate total dust features
    dust_features = []

    for cmodel in pmodel.model:
        if isinstance(cmodel, Drude1D):
            dust_features.append(cmodel)
    dust_features_model = dust_features[0]
    for cmodel in dust_features[1:]:
        dust_features_model += cmodel

    # calculate H2 spectrum
    h2_features = []

    for cmodel in pmodel.model:
        if isinstance(cmodel, Gaussian1D):
            if cmodel.name[0:2] == 'H2':
                h2_features.append(cmodel)
    h2_features_model = h2_features[0]
    for cmodel in h2_features[1:]:
        h2_features_model += cmodel

    # calculate atomic line spectrum
    atomic_features = []

    for cmodel in pmodel.model:
        if isinstance(cmodel, Gaussian1D):
            if cmodel.name[0:2] != 'H2':
                atomic_features.append(cmodel)
    atomic_features_model = atomic_features[0]
    for cmodel in atomic_features[1:]:
        atomic_features_model += cmodel

    # all atomic and H2 lines
    totlines = h2_features_model(x) + atomic_features_model(x)

    # get extinction model
    for cmodel in pmodel.model:
        if isinstance(cmodel, S07_attenuation):
            ext_model = cmodel(x)

    # save compounds in dictionary
    compounds = {}
    compounds["x"] = x
    compounds["tot_cont"] = totcont
    compounds["stellar_cont"] = stellar_cont_model(x)
    compounds["dust_cont"] = dust_cont_model(x)
    compounds["dust_features"] = dust_features_model(x)
    compounds["tot_lines"] = totlines
    compounds["h2_lines"] = h2_features_model(x)
    compounds["atomic_lines"] = atomic_features_model(x)
    compounds["extinction_model"] = ext_model

    return compounds
