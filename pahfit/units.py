import astropy.units as u
from astropy.units import CompositeUnit

from enum import Enum

# Working/default unitParameter default units: flux density/intensity/power
# These are PAHFITs default science packs parameter and output units
class UNITS(Enum):
    temperature = u.K
    wavelength = u.um
    fwhm = u.um
    flux_density = u.mJy
    flux_power = CompositeUnit(1e-22, (u.W, u.m), (1, -2))
    intensity = u.MJy / u.sr
    intensity_power = CompositeUnit(1e-10, (u.W, u.m, u.sr), (1, -2, -1))

# integrated power units of 1e-22 W/m^2 (from flux) corresponds to the
# unit 1e-10 W/m^2/sr (from intensity) if it occurs uniformly over a
# solid angle 0.21" on a side (about a small JWST IFU pixel)
