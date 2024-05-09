"""pahfit.util General pahfit.features utility functions."""
import numpy as np


def bounded_is_missing(val):
    """Return a mask array indicating which of the bounded values
    are missing.  A missing bounded value has a masked value."""
    return getattr(a['val'], 'mask', None) or np.zeros_like(a['val'], dtype=bool)


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
