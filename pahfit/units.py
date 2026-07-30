import numpy as np
import astropy.units as u
from astropy.units import CompositeUnit


# Working/default unitParameter default units: flux density/intensity/power
# These are PAHFITs default science packs parameter and output units
temperature = u.K
wavelength = u.um
flux_density = u.mJy
flux_power = CompositeUnit(1e-22, (u.W, u.m), (1, -2))
solid_angle = u.sr
intensity = u.MJy / u.sr
intensity_power = CompositeUnit(1e-10, (u.W, u.m, u.sr), (1, -2, -1))

# Note: integrated power units of 1e-22 W/m^2 (from flux) corresponds
# to the unit 1e-10 W/m^2/sr (from intensity) if it occurs uniformly
# over a solid angle 0.21" on a side (about a small JWST IFU pixel)

def is_surface_brightness(unit):
    """True if `unit` is a surface brightness (has 1/sr), False if a flux.
    Method suggested by J.D. Smith in issue #28."""

    decomposed = unit.decompose()
    return (u.rad, -2) in zip(decomposed.bases, decomposed.powers)

def get_quantity(features, name, column):
    """Get one feature's parameter value with its unit attached
    (e.g. "5.2 mJy"). Returns None if masked, instead of an
    incorrect 0. Columns with no assigned unit (e.g. tau, which
    is dimensionless) are returned as plain dimensionless
    Quantities, not multiplied against None.
    """
    row = features.loc[name]
    if row[column] is np.ma.masked:
        return None
    unit = features[column].unit
    if unit is None:
        unit = u.dimensionless_unscaled
    return row[column]["val"] * unit


def working_units(is_flux):
    """Return (value_unit, power_unit) for the given track.

    is_flux=True  -> (flux_density, flux_power)       e.g. mJy, 1e-22 W/m^2
    is_flux=False -> (intensity, intensity_power)      e.g. MJy/sr, 1e-10 W/m^2/sr

    Centralizes the flux-vs-surface-brightness unit choice that was
    previously duplicated across _convert_spec_data, guess(), and the
    Power* fitting components.
    """
    if is_flux:
        return flux_density, flux_power
    return intensity, intensity_power
