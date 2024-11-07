"""pahfit.util General pahfit.features utility functions."""
import numpy as np
from pahfit.errors import PAHFITFeatureError


def value_bounds(val, bounds):
    """Compute bounds for a bounded value.

    Arguments:
    ----------
      val: The value to bound.

      bounds: Either None for no relevant bounds (i.e. fixed), or a
        two element iterable specifying (min, max) bounds.  Each of
        min and max can be a numerical value, None (for infinite
        bounds, either negative or positive, as appropriate), or a
        string ending in:in

          #: an absolute offset from the value
          %: a percentage offset from the value

        Offsets are necessarily negative for min bound, positive for
        max bounds.

    Examples:
    ---------

      A bound of ('-1.5%', '0%') would indicate a minimum bound
        1.5% below the value, and a max bound at the value itself.

      A bound of ('-0.1#', None) would indicate a minimum bound 0.1 below
        the value, and infinite maximum bound.

    Returns:
    -------

      A 3 element tuple (value, min, max).
        Any missing bound is replaced with the numpy.nan value.

    Raises:
    -------

      ValueError: if bounds are specified and the value does not fall
        between them.

    """
    if val is None:
        val = np.ma.masked
    if not bounds:
        return (val,) + 2 * (np.nan,)  # (val,nan,nan) indicates fixed
    ret = [val]
    for i, b in enumerate(bounds):
        if isinstance(b, str):
            if b.endswith('%'):
                b = val * (1. + float(b[:-1]) / 100.)
            elif b.endswith('#'):
                b = val + float(b[:-1])
            else:
                raise PAHFITFeatureError(f"Incorrectly formatted bound {b}")
        elif b is None:
            b = np.inf if i else -np.inf
        ret.append(b)

    if (val < ret[1] or val > ret[2]):
        raise PAHFITFeatureError(f"Value <{ret[0]}> is not between bounds: {ret[1:]}")
    return tuple(ret)


def bounded_is_missing(val):
    """Return a mask array indicating which of the bounded values
    are missing.  A missing bounded value has a masked value."""
    return getattr(val['val'], 'mask', None) or np.zeros_like(val['val'], dtype=bool)


def bounded_is_fixed(val):
    """Return a mask array indicating which of the bounded values
    are fixed.  A fixed bounded value has masked bounds."""
    return np.isnan(val['min']) & np.isnan(val['max'])


def bounded_min(val):
    """Return the minimum of each bounded value passed.
    Either the lower bound, or, if no such bound is set, the value itself."""
    lower = val['min']
    return np.where(lower, lower, val['val'])


def bounded_max(val):
    """Return the maximum of each bounded value passed.
    Either the upper bound, or, if no such bound is set, the value itself."""
    upper = val['max']
    return np.where(upper, upper, val['val'])
