import os
import pkg_resources

import astropy.units as u
from astropy.table import Table

from astropy.modeling.fitting import LevMarLSQFitter

from pahfit.base import PAHFITBase

__all__ = ["read_spectrum", "initialize_model", "fit_spectrum"]


def read_spectrum(specfile):
    """
    Read in a spectrum and convert intput units to the expected internal PAHFIT units.

    Parameters
    ----------
    specfile : string
        file with the spectrum to be fit

    Returns
    -------
    obsdata : dict of ndarrays
        x is wavelength in microns and y/unc are the spectrum/unc in units
        of Jy
    """

    # read in the observed spectrum
    # assumed to be astropy table compatibile and include units
    specfile = specfile
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
    obsdata["x"] = obs_spectrum["wavelength"].to(u.micron, equivalencies=u.spectral())
    obsdata["y"] = obs_spectrum["flux"].to(u.Jy, equivalencies=u.spectral_density(obsdata["x"]))
    obsdata["unc"] = obs_spectrum["sigma"].to(u.Jy, equivalencies=u.spectral_density(obsdata["x"]))

    return obsdata


def initialize_model(packfile, obsdata, estimate_start):
    """
    Initialize a model based on the packfile

    Parameters
    ----------

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
        obsdata["x"].value, obsdata["y"].value, estimate_start=estimate_start, filename=packfile
    )

    return pmodel


def fit_spectrum(obsdata, pmodel):
    """
    Read a spectrum from disk and fit it

    Parameters
    ----------

    Returns
    -------
    obsfit
    """

    # pick the fitter
    fit = LevMarLSQFitter()

    # fit
    obs_fit = fit(
        pmodel.model,
        obsdata["x"].value,
        obsdata["y"].value,
        weights=1.0 / obsdata["unc"].value,
        maxiter=200,
        epsilon=1e-10,
        acc=1e-10,
    )
    print(fit.fit_info["message"])

    return obs_fit
