"""pahfit.packs.instrument

Simple interface to instrument packs, which are responsible for
defining wavelength-dependent unresolved-line resolution, and range.

Note that individual spectral ranges are called 'segments', and PAHFIT
accepts a pre-defined list of segments, as fully-qualified
heirarchical names, e.g. 'iso.sws.speed0.3a'.  For the full list of
supported telescopes, instruments, and instrument modes, see TBD.
"""

import glob, os
from numpy.polynomial import Polynomial
from astropy.io.misc import yaml
from pkg_resources import resource_filename
from pahfit.errors import PAHFITPackError

packs = {}

def read_instrument_packs():
    """Read all instrument packs into the 'packs' variable."""
    for pack in glob.glob(resource_filename("pahfit", "packs/instrument/*.yaml")):
        try:
            with open(pack) as fd:
                p = yaml.load(fd)
        except IOError as e:
            raise PAHFITPackError("Error reading instrument pack file\n"
                                  f"\t{pack}\n\t{repr(e)}")
        else:
            telescope = os.path.basename(pack).rstrip('.yaml')
            packs[telescope] = p
            

def pack_element(segment):
    """Return the pack element for the given segment name.

    Arguments:
    ----------
      segment: The fully qualified segment name, as a string.

    Returns:
    --------

    The dictionary structure including the keys 'range' and
    'coefficients' of the named segment.

    Raises:
    -------

    PAHFITPackError: If the named segment cannot be identified.
    """
    global packs
    d = packs
    for key in segment.split('.'):
        d = d.get(key)
        if d is None:
            raise PAHFITPackError(f"Could not locate instrument segment {key} of {segment}")
    if d.get('polynomial') is None:
        try:
            d['polynomial'] = Polynomial(d['coefficients'])
        except KeyError:
            raise PAHFITPackError(f"Incomplete segment name {segment}")

    return d


def ins_name(path, tree):
    if isinstance(tree, dict) and not tree.get("range"):
        _dot = "." if len(path) > 0 else ""
        for key, sub in tree.items():
            yield from ins_name(f"{path}{_dot}{key}", sub)
    else:
        yield path

def instruments():
    """Return the fully qualified set of supported instrument names."""
    if len(packs) == 0:
        read_instrument_packs()
    return ins_name('', packs)
        
def resolution(segment, wave_micron):
    p = pack_element(segment)['polynomial']
    return p(wave_micron)


def fwhm(segment, wave_micron):
    """Return the FWHM for SEGMENT at one more more wavelengths.

    Arguments:
    ----------
      segment: The fully qualified segment name, as a string.

      wave_micron: The observed-frame (instrument-relative) wavelength
      in microns, as a scalar or numpy array.

    Returns:
    --------

    The full-width at half maximum of an unresolved line spread
    function, at the relevant wavelength(s).
    """
    return wave_micron/resolution(segment, wave_micron)

def wave_range(segment):
    """Return the segment wavelength range (a two element list) in microns."""
    return pack_element(segment)['range']

read_instrument_packs()
