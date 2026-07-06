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
    """Return True if `unit` is a surface brightness, False if a flux.

    Surface brightness units carry a per-solid-angle component (e.g.
    MJy/sr). We detect this by decomposing the unit to its base units
    and checking for an inverse-solid-angle term. Steradian decomposes
    to rad**2, so "per steradian" shows up as (u.rad, -2).

    Method suggested by J.D. Smith in issue #28.
    """
    decomposed = unit.decompose()
    return (u.rad, -2) in zip(decomposed.bases, decomposed.powers)

def get_quantity(features, name, column):
    """Get one feature's parameter value, with its unit attached.

    Parameters
    ----------
    features : Features table
        The PAHFIT features table (e.g. model.features).
    name : str
        Feature name, e.g. '[NeII]'.
    column : str
        Parameter name, e.g. 'power', 'wavelength', 'fwhm'.

    Returns
    -------
    astropy.units.Quantity
        The fitted value, with its unit attached (e.g. "5.2 mJy").
    """
    row = features.loc[name]
    return row[column]["val"] * features[column].unit