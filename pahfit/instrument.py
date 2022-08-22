"""pahfit.packs.instrument

Simple interface to instrument packs, which are responsible for
defining wavelength-dependent unresolved-line resolution, and range.

Note that individual spectral ranges are called 'segments', and PAHFIT
accepts a pre-defined list of segments, as fully-qualified
heirarchical names, e.g. 'iso.sws.speed0.3a'.  For the full list of
supported telescopes, instruments, and instrument modes, see TBD.
"""

import os
from pathlib import Path
import glob
import numpy as np
from numpy.polynomial import Polynomial
from astropy.io.misc import yaml
from pkg_resources import resource_filename
from pahfit.errors import PAHFITPackError

packs = {}


def read_instrument_packs():
    """Read all instrument packs into the 'packs' variable."""
    global packs
    for pack in glob.glob(resource_filename("pahfit", "packs/instrument/*.yaml")):
        try:
            with open(pack) as fd:
                p = yaml.load(fd)
        except IOError as e:
            raise PAHFITPackError(
                "Error reading instrument pack file\n" f"\t{pack}\n\t{repr(e)}"
            )
        else:
            telescope = os.path.basename(pack).rstrip(".yaml")
            packs[telescope] = p
    packs = dict(_ins_items("", packs))  # Flatten list


def pack_element(segments):
    """Return the pack element(s) for the given segment name(s).

    Arguments:
    ----------

      segments: The fully qualified segment name, as a string. Can be a
        scalar or a list, for returning more than one pack.  Each
        segment name can optionally contain glob characters recognized
        by Pathlib.Path.match (e.g. `*' or [123]).  This can be used
        to match multiple segments.

    Returns:
    --------

    A list of dictionary structures including the keys 'range' and
    'coefficients' of the named segment(s).

    Raises:
    -------

    PAHFITPackError: If none of the named segments can be identified.

    """
    global packs
    ret = []

    if not isinstance(segments, (tuple, list)):
        segments = (segments,)

    for segment in segments:
        sm = [p for p in packs.keys() if Path(p).match(segment)]
        if not sm:
            raise PAHFITPackError(f"Could not locate instrument segment {segment}")
        for s in sm:
            if packs[s].get("polynomial") is None:
                try:
                    packs[s]["polynomial"] = Polynomial(packs[s]["coefficients"])
                except KeyError:
                    raise PAHFITPackError(f"Invalid instrument pack {s}")
            ret.append(packs[s])

    return ret


def _ins_items(path, tree):
    """Generator for traversing packs."""
    if isinstance(tree, dict) and not tree.get("range"):
        _dot = "." if len(path) > 0 else ""
        for key, sub in tree.items():
            yield from _ins_items(f"{path}{_dot}{key}", sub)
    else:
        yield path, tree


def instruments(match=None):
    """Returns fully qualified set of supported instrument names.

    Arguments:
    ----------

      match (Optional): A string or iterable of strings to match
       against the full instrument list. If provide, only the list of
       matching instruments is returned.  Useful for checking a
       glob-style instrument name results in the desired instrument
       set.

    Returns:
    --------

      The list of matching instruments, if any.

    """
    if len(packs) == 0:
        read_instrument_packs()
    ins = packs.keys()
    ret = []
    if match:
        if not isinstance(match, (tuple, list)):
            match = (match,)
        for m in match:
            ret.extend([p for p in ins if Path(p).match(m)])
        return ret
    else:
        return list(ins)


def resolution(segment, wave_micron):
    """
    Parameters
    ----------
    segment: The fully qualified segment name, as a string.

    wave_micron: The observed-frame (instrument-relative) wavelength
        in microns, as a scalar or numpy array of size N.

    Returns
    -------
    r: 1D array shape (N,) or 2D masked_array shape (N, 3)
        The spectral resolution at the given wavelengths. In case of
        multiple segments: r[:, 0] is mean resolution, r[:, 1] is min
        resolution or masked., r[:, 2] is max resolution or masked. The
        values of the latter two arrays are masked if only one segment
        covers this wavelength.
    """
    _packs = pack_element(segment)
    npk = len(_packs)
    if npk == 1:
        return _packs[0]["polynomial"](wave_micron)

    wave_micron = np.atleast_1d(wave_micron)

    res = np.ma.empty((npk,) + wave_micron.shape)
    for i, p in enumerate(_packs):
        inside = (wave_micron >= p["range"][0]) & (wave_micron <= p["range"][1])
        res[i, inside] = p["polynomial"](wave_micron[inside])
        res[i, ~inside] = np.ma.masked

    out = np.ma.empty_like(wave_micron, shape=wave_micron.shape + (3,))
    out[:] = np.ma.masked  # mask all by default
    out[..., 0] = np.mean(res, 0)

    multi = np.count_nonzero(res, 0) > 1  # waves matching more than one segment
    if np.count_nonzero(multi) > 0:
        out[multi, 1] = np.min(res[:, multi], 0)
        out[multi, 2] = np.max(res[:, multi], 0)

    return out


def fwhm(segment, wave_micron):
    """Return the FWHM for SEGMENT at one more more wavelengths.

    Arguments:
    ----------
      segment: The fully qualified segment name, as a string.

      wave_micron: The observed-frame (instrument-relative) wavelength
        in microns, as a scalar or numpy array of size N.

    Returns:
    --------

    The full-width at half maxima of unresolved line spread functions
    at the relevant wavelength(s), either as a scalar (if only one
    segment applied) or as a (N, 3) shape array of:

      [[value, low bound, high bound], ...] for every wavelength

    """

    r = resolution(segment, wave_micron)
    if np.ma.isMaskedArray(r):
        ret = np.expand_dims(wave_micron, -1) / r
        return ret[..., (0, 2, 1)]  # swap lower with upper (inverse!)
    else:
        return wave_micron / r


def fwhm_recommendation(segment, wave_micron):
    """Returns recommended parameters for the fwhm.

    This function checks the shape of the output of fhwm(), and returns consistent values.
    - value as a starting point
    - upper and lower, or masked where there is no overlap
    - 'fixed' flag

    When a wavelength is covered by only one segment, the recommendation
    is to fix the fwhm. In case of multiple, it should be variable,
    between the given upper and lower bounds.

    Parameters
    ----------
    segment: The fully qualified segment name, as a string.

    wave_micron: The observed-frame (instrument-relative) wavelength
        in microns, as a scalar or numpy array of size N.

    Returns
    -------
    Tuple of lists / arrays. Each one of size len(wave_micron)
        (array float, list bool, array float, array float)
        (fwhm, bool fixed, low bound, high bound)

    """
    N = len(wave_micron)
    fwhm_output = fwhm(segment, wave_micron)
    if len(fwhm_output.shape) == 1:
        return fwhm_output, [True] * N, [0.] * N, [0.] * N
    else:
        # We need to be careful here, because for astropy a numpy.bool
        # does not work for the 'fixed' parameter. It needs to be a
        # regular bool. Tolist() solves this.
        fixed = fwhm_output[: , 1].mask.tolist()
        return fwhm_output[:, 0], fixed, fwhm_output[:, 1], fwhm_output[:, 2]


def wave_range(segment):
    """Return the segment wavelength range(s) in microns.  If more
    than one segment is specified (either directly, or via glob-style
    matching), the return value is a list of two element lists [min,
    max].
    """
    ret = [x["range"] for x in pack_element(segment)]
    if len(ret) == 1:
        return ret[0]
    else:
        return ret


def test_wave_minmax(wave_micron, segments):
    """Return True if the instrument range encompasses the given wavelengths

    For now, it is just based on the mininum and maximum, but it could
    check all the wavelengths in the future, to make sure that all
    observational data is covered by a valid instrument model.

    Parameters
    ----------
    wave_micron: 1D array-like (not quantity)
        list of wavelengths in micron

    segments: list of str
        list of instrument segment names, as defined by the pack_element function

    """
    ranges = [x["range"] for x in pack_element(segments)]
    wmin = min(wave_micron)
    wmax = max(wave_micron)
    rmin = min(r[0] for r in ranges)
    rmax = max(r[1] for r in ranges)
    return rmin <= wmin and wmax <= rmax


def test_waves_in_any_segment(wave_micron, segments):
    """Test if each given wavelength is covered by at least one instrument

    Parameters
    ----------
    wave_micron: dimensionless numpy array of wavelengths in micron

    segments: list of segments

    Returns
    -------
    np.array of bool, True where wave_micron[i] is is the range of any of the segments."""

    # for each instrument range, test all the wavelengths (not necessarily sorted)
    ranges = [x["range"] for x in pack_element(segments)]
    in_range_per_segment = [
        (rmin <= wave_micron) & (wave_micron <= rmax) for (rmin, rmax) in ranges
    ]
    in_range_any_segment = np.any(np.array(in_range_per_segment), axis=0)
    return in_range_any_segment


read_instrument_packs()
