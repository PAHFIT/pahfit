import numpy as np
from astropy import constants as const


def pah_feature_strength(ampl, fwhm, x_0):
    """
    Calculate PAH feature strengths.

    Parameters
    ----------
    ampl : float
        Drude1D amplitude
    fwhm : float
        Drude1D fwhm
    x_0 : float
        Drude1D peak position

    Returns
    -------
    strength : float
        PAH feature strength in W/m^2/Hz units.

    """

    # TODO: Use input units to determine conversion and feature strength units.
    strength = (np.pi * const.c.to('micron/s') / 2) * (ampl * fwhm / x_0**2) * 1e-26

    return strength

def line_strength(ampl, mean, stddev):
    """
    Calculate emission line strengths.

    Parameters
    ----------
    ampl : float
        Gaussian1D amplitude
    mean : float
        Mean of the Gaussian
    stddev : float
        Gaussian1D stddev

    Returns
    -------
    strength : float
        Emission line strength in W/m^2/Hz units.

    """

    # TODO: Use input units to determine conversion and line strength units.
    fwhm = 2 * stddev * np.sqrt(2 * np.log(2))
    strength = ((np.sqrt(np.pi / np.log(2))) * const.c.to('micron/s').value / 2) * (ampl * fwhm / mean**2) * 1e-26

    return strength
