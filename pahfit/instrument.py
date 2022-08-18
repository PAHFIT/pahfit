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
from pahfit.errors import PAHFITPackError, PAHFITWarning
from warnings import warn

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


def _ins_items(path, tree):
    """Generator for traversing packs."""
    if isinstance(tree, dict) and not tree.get("range"):
        _dot = "." if len(path) > 0 else ""
        for key, sub in tree.items():
            yield from _ins_items(f"{path}{_dot}{key}", sub)
    else:
        yield path, tree


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
                    p = Polynomial(packs[s]['coefficients'])
                    packs[s]['polynomial'] = p
                    packs[s]['range_fwhm'] = [x / p(x) for x in packs[s]['range']]
                except KeyError:
                    raise PAHFITPackError(f"Invalid instrument pack {s}")
            ret.append(packs[s])

    return ret


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


def resolution(segment, wave_micron, fwhm_near=0.0, as_bounded=False):
    """Return the resolution for `segment' at one more more
    wavelengths.

    Arguments:
    ----------

      segment: The segment name or list of names, potentially
        including glob chars.  See pack_element.

      wave_micron: The observed-frame (instrument-relative) wavelength
        in microns, as a scalar or numpy array of any shape.

      fwhm_near (optional, default: 0.0): Compute resolution for any
        wavelength within this many FWHMs of the segment end.  Note
        that, to avoid extrapolation, the FWHM of interest is
        calculated at the segment's endpoints.

      as_bounded (optional, default: False): always return bounded
        variables (i.e. with shape (..., 3)), even when only one segment
        is used.

    Returns:
    --------

    The resolution (wavelength divided by full-width at half maxima)
    of unresolved line spread functions at the relevant wavelength(s),
    either as a scalar or array (if only one segment applied and
    `as_bounded' is not True) or as a masked array of 3 element tuples
    specifying:

      (value, low bound, high bound)

    otherwise.  Note that, when more than one segment matches (and/or
    `as_bounded` is True), only wavelengths inside (or near, if
    `fwhm_near` is non-zero) segments have resolution values returned.
    They are masked otherwise.

    """

    _packs = pack_element(segment)
    npk = len(_packs)

    if npk == 1 and not as_bounded:
        return _packs[0]['polynomial'](wave_micron)

    wave_micron = np.atleast_1d(wave_micron)

    res = np.ma.empty((npk,) + wave_micron.shape)
    for i, p in enumerate(_packs):
        inside = ((wave_micron >= p['range'][0] - p['range_fwhm'][0] * fwhm_near)
                  & (wave_micron <= p['range'][1] + p['range_fwhm'][1] * fwhm_near))
        res[i, inside] = p['polynomial'](wave_micron[inside])
        res[i, ~inside] = np.ma.masked

    out = np.ma.empty_like(wave_micron, shape=wave_micron.shape + (3,))
    out[:] = np.ma.masked  # mask all by default
    out[..., 0] = np.mean(res, 0)

    multi = np.count_nonzero(res, 0) > 1  # waves matching more than one segment
    if np.count_nonzero(multi) > 0:  # add lower/upper bounds
        out[multi, 1] = np.min(res[:, multi], 0)
        out[multi, 2] = np.max(res[:, multi], 0)

    return out


def fwhm(segment, wave_micron, **kwargs):
    """Return the FWHM for `segment' at one more more wavelengths.

    Arguments:
    ----------

      segment: The segment name or list of names, potentially
        including glob chars.  See pack_element.

      wave_micron: The observed-frame (instrument-relative) wavelength
        in microns, as a scalar or numpy array of any shape.

      kwargs: other keyword args as accepted by `resolution'.

    Returns:
    --------

    The fwhm in microns for `wave_micron', in the format specified by
    `resolution'.
    """

    r = resolution(segment, wave_micron, **kwargs)
    if np.ma.isMaskedArray(r):
        ret = np.expand_dims(wave_micron, -1) / r
        return ret[..., (0, 2, 1)]  # swap lower with upper (inverse!)
    else:
        return wave_micron / r


def fwhm_recommendation(segment, wave_micron):
    """Returns recommended parameters for the fwhm.

    This function checks the shape of the output of fhwm(), and returns consistent values.
    - Value as a starting point
    - Min and max where there is segment overlap, masked where there is no overlap
    - 'fixed' flag, False where there is overlap

    Parameters
    ----------
    segment: The fully qualified segment name, as a string.

    wave_micron: The observed-frame (instrument-relative) wavelength
        in microns, as a scalar or numpy array of size N.

    Returns
    -------
    Tuple of lists / arrays. Each one of size len(wave_micron)
        (float, bool, float, float)
        (fwhm, bool fixed, low bound, high bound)

    """
    N = len(wave_micron)
    fwhm_output = fwhm(segment, wave_micron)
    if len(fwhm_output.shape) == 1:
        return fwhm_output, [True] * N, [0.] * N, [0.] * N
    else:
        # when a wavelength is covered by only one segment, the
        # recommendation is to fix the fwhm. In case of multiple, it
        # should be variable, between the given upper and lower bounds
        fixed = fwhm_output[: , 1].mask
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


_frac = np.array([0.03, 0.01, 0.])  # Fractional bounds for error and warnings


def check_range(wave_bounds, segments):
    """Check that the wave_bounds of the passed spectra are consistent
    with the given instrument segment(s).

    Arguments:
    ----------
      wave_bounds: The min and max observed
        wavelength bounds in microns in the observed spectrum or
        spectral segment of interest.

      segments: The segment name or list of names, potentially
        including glob chars.  See `pack_element'.

    Returns:
    --------

    True, if the wave_bounds are inside or near the wavelength range
    given the specified segment(s).  Note that only the ends of the
    segment ranges are checked.  If wave_bounds exceed the specified
    instrument pack segment, a warning is issued.  If they exceed at
    either end by more than 3% of the segment width, an error is
    thrown.

    Raises:
    -------

    PAHFITPackError: If the wave bounds exceed the ends of the
      specified instrument segment(s) by more than 3%.

    """

    els = pack_element(segments)
    low = [s['range'][0] for s in els]
    high = [s['range'][1] for s in els]

    mmpos = (np.argmin(low), np.argmax(high))
    mn, mx = np.min(low), np.max(high)

    full_range = [els[x]['range'][1] - els[x]['range'][0] for x in mmpos]
    bins = [np.searchsorted(mn - _frac * full_range[0], wave_bounds[0], side='right'),
            np.searchsorted(mx + _frac[::-1] * full_range[1], wave_bounds[1])]

    if bins[0] == 3 and bins[1] == 0:  # inside range
        return True

    st = (f":\n\tSegment(s) cover: {mn:.2f}-{mx:.2f}, "
          + f"wave_bounds: {wave_bounds[0]:.2f}-{wave_bounds[1]:.2f}")
    if bins[0] == 0 or bins[1] == 3:
        raise PAHFITPackError(
            f"\nInput wavelength bounds exceed instrument segment(s) by more than 3%{st}")
    else:
        w = ' by more than 1%' if bins[0] == 1 or bins[1] == 2 else ''
        warn(f"\nInput wavelength bounds exceed instrument segment(s){w}{st}",
             PAHFITWarning)
        return True


def within_segment(wave_micron, segments, fwhm_near=None, wave_bounds=None):
    """Return a mask indicating which wavelengths are in (or near) the
    wavelength range of the given segment(s).

    Arguments:
    ----------
      segment: The segment name or list of names, potentially
        including glob chars.  See `pack_element'.

      wave_micron: The observed-frame (instrument-relative) wavelength
        in microns, as a scalar or numpy array of any shape.

      fwhm_near (optional, default: None): If not None, mark a
        wavelength in a segment (or set of segments) if it falls
        within `fwhm_near'*FWHM of any segment's end.  Note that, to
        avoid extrapolation, the FWHM of interest is calculated at the
        segment endpoint.

      wave_bounds (optional, default: None): The min and max observed
        wavelength bounds in microns in the observed spectrum or
        spectral component corresponding to segment.  It is assumed
        that the bounds have been checked against the instrument
        segments using check_range.

    Returns:
    --------

    A boolean mask of identical shape as `wave_micron', indicating
    which wavelength are within (or near to) some segment.

    See Also:
    ---------
    check_range

    """
    res = []

    els = pack_element(segments)
    low = [s['range'][0] for s in els]
    high = [s['range'][1] for s in els]
    mnpos, mxpos = np.argmin(low), np.argmax(high)
    for i, s in enumerate(els):
        rng = s['range'].copy()
        if wave_bounds:  # replace extreme segment bounds with actual limits
            if i == mnpos:
                rng[0] = wave_bounds[0]
            elif i == mxpos:
                rng[1] = wave_bounds[1]
        if fwhm_near:  # Account for wing overlap
            rng[0] -= s['range_fwhm'][0] * fwhm_near
            rng[1] += s['range_fwhm'][1] * fwhm_near
        res.append((wave_micron >= rng[0]) & (wave_micron <= rng[1]))
    return np.any(res, 0)


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
