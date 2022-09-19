import os
import pkg_resources

import astropy.units as u
from astropy.table import Table

from specutils import Spectrum1D
from astropy.nddata import StdDevUncertainty

from pahfit.component_models import BlackBody1D, S07_attenuation
from astropy.modeling.physical_models import Drude1D
from astropy.modeling.functional_models import Gaussian1D

__all__ = ["read_spectrum", "calculate_compounds"]


def find_packfile(packfile):
    """Determine packfile path.

    If packfile not in current directory, try to find one in the packs
    directory that comes with PAHFIT. If nothing is found, throw an
    error.

    Parameters
    ----------
    packfile : str
        Name or path of the pack file.

    Returns
    -------
    packfile_found : str
        Full path of the pack file. Will be a file in pahfit/packs when
        the given file name was not found.

    """
    if os.path.isfile(packfile):
        packfile_found = packfile
    else:
        pack_path = pkg_resources.resource_filename("pahfit", "packs/science/")
        test_packfile = "{}/{}".format(pack_path, packfile)
        if os.path.isfile(test_packfile):
            packfile_found = test_packfile
        else:
            raise ValueError("Input packfile {} not found".format(packfile))

    return packfile_found


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
    spec1d : Spectrum1D
        spectral_axis in microns, flux and uncertainties in units of Jy
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
    x = obs_spectrum[colnames[0]].to(u.micron, equivalencies=u.spectral())
    y = obs_spectrum[colnames[1]].to(u.Jy, equivalencies=u.spectral_density(x))
    unc = obs_spectrum[colnames[2]].to(u.Jy, equivalencies=u.spectral_density(x))

    return Spectrum1D(spectral_axis=x, flux=y, uncertainty=StdDevUncertainty(unc))


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
            if cmodel.name[0:2] == "H2":
                h2_features.append(cmodel)
    h2_features_model = h2_features[0]
    for cmodel in h2_features[1:]:
        h2_features_model += cmodel

    # calculate atomic line spectrum
    atomic_features = []

    for cmodel in pmodel.model:
        if isinstance(cmodel, Gaussian1D):
            if cmodel.name[0:2] != "H2":
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
