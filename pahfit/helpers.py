import os
from importlib import resources

from pahfit import units

from specutils import Spectrum1D
from astropy import units as u
from astropy.nddata import StdDevUncertainty

__all__ = ["find_packfile", "read_spectrum"]


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
        test_packfile = resources.files("pahfit") / "packs/science" / packfile
        if os.path.isfile(test_packfile):
            packfile_found = test_packfile
        else:
            raise ValueError("Input packfile {} not found".format(packfile))

    return str(packfile_found)


def read_spectrum(specfile, format=None):
    """
    Find included spectrum file and read it in as a Spectrum1D object.

    Parameters
    ----------
    specfile : string
        File name. Will resolve to a path relative to the working
        directory, or if not found, a path relative to the PAHFIT
        included data directory (pahfit/data).

    format : string
        Format option to pass to Spectum1D.read

    Returns
    -------
    spec1d : Spectrum1D
        spectral_axis in microns, flux and uncertainties in units of Jy
    """
    # resolve filename
    if not os.path.isfile(specfile):
        test_specfile = resources.files("pahfit") / "data" / specfile
        if os.path.isfile(test_specfile):
            specfile = test_specfile
        else:
            raise ValueError("Input spectrumfile {} not found".format(specfile))

    # File assumed to be compatible with specutils.Spectrum1D.read
    # Default option is to auto-identify format
    tformat = None
    # process user-specified or filename extension based format
    if format is None:
        suffix = specfile.name.split(".")[-1].lower()
        if suffix == "ecsv":
            tformat = "ECSV"
        elif suffix == "ipac":
            tformat = "IPAC"
    else:
        tformat = format

    s = Spectrum1D.read(specfile, format=tformat)

    # Convert to intensity units by assuming an arbitrary solid angle
    # for now. To be removed when dual unit support (intensity and flux
    # density) is supported.
    if s.flux.unit.is_equivalent(units.flux_density):
        solid_angle = (3 * u.arcsec) ** 2
        alt_flux = (s.flux / solid_angle).to(units.intensity)
        alt_unc_array = (s.uncertainty.array * s.flux.unit / solid_angle).to(
            units.intensity
        )
        s = Spectrum1D(
            alt_flux, s.spectral_axis, uncertainty=StdDevUncertainty(alt_unc_array)
        )

    return s
