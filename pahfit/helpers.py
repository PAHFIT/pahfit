import os
import pkg_resources

import astropy.units as u
from astropy.table import Table

from astropy.modeling.fitting import LevMarLSQFitter

from pahfit.base import PAHFITBase

from pahfit.component_models import BlackBody1D
from astropy.modeling.physical_models import Drude1D
from astropy.modeling.functional_models import Gaussian1D


__all__ = ["read_spectrum", "initialize_model", "fit_spectrum"]


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
    obsfit : PAHFITBase model
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


def get_compounds(obsdata, pmodel):
    """
    Determine model compounds for total continuum, dust continuum, dust feature, H2 lines

    Parameters
    ----------
    obsdata : dict
        observed data where x = wavelength, y = SED, and unc = uncertainties

    obs_fit : CompoundModel
        PAHFIT model (output from fit_spectrum)

    Returns
    -------
    compounds : dict
        x = wavelength in microns; 
        tot_cont = total continuum; 
        dust_cont = total dust continuum; 
        dust_feat = combined dust features; 
        h2_feat = combined H2 lines
    """

    # get wavelength array
    x = obsdata["x"].value

    # calculate total dust continuum and total continuum (oncluding stellar continuum)
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
    H2_features = []

    for cmodel in pmodel.model:
        if isinstance(cmodel, Gaussian1D):
            if cmodel.name[0:2] == 'H2':
                H2_features.append(cmodel)
    H2_features_model = H2_features[0]
    for cmodel in H2_features[1:]:
        H2_features_model += cmodel

    # save compounds in dictionary
    compounds = {}
    compounds["x"] = x
    compounds["tot_cont"] = totcont
    compounds["dust_cont"] = dust_cont_model(x)
    compounds["dust_feat"] = dust_features_model(x)
    compounds["h2_feat"] = H2_features_model(x)

    return compounds
