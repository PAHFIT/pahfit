"""pahfit.util General pahfit.features utility functions."""
import numpy as np
import numpy.ma as ma


def bounded_is_missing(val):
    """Return a mask array indicating which of the bounded values
    are missing.  A missing bounded value has a masked value."""
    return ma.getmask(val)[..., 0]


def bounded_is_fixed(val):
    """Return a mask array indicating which of the bounded values
    are fixed.  A fixed bounded value has masked bounds."""
    return ma.getmask(val)[..., -2:].all(-1)


def bounded_min(val):
    """Return the minimum of each bounded value passed.
    Either the lower bound, or, if no such bound is set, the value itself."""
    lower = val[..., 1]
    return np.where(lower, lower, val[..., 0])


def bounded_max(val):
    """Return the maximum of each bounded value passed.
    Either the upper bound, or, if no such bound is set, the value itself."""
    upper = val[..., 2]
    return np.where(upper, upper, val[..., 0])
