from astropy.modeling.functional_models import Gaussian1D
from pahfit.component_models import BlackBody1D

import numpy as np

from astropy import constants as const
from astropy.table import Table
from astropy.modeling.physical_models import Drude1D

from scipy import integrate


def pah_feature_strength(ampl, fwhm, x_0):
    """
    Calculate the integrated area of PAH feature profiles.

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


def featcombine(ftable):
    """
    Combine dust features strengths.
    Parameters
    ----------
    ftable : data table
        Astropy table containing lines.
    Returns
    -------
    cftable : data table
        Astropy table containing combined dust feature strengths.
    """
    # Define combined strength dictionary.
    cfdic = {'PAH_62': {'range': [6.2, 6.3]},
             'PAH_77_C': {'range': [7.3, 7.9]},
             'PAH_83': {'range': [8.3, 8.4]},
             'PAH_86': {'range': [8.6, 8.7]},
             'PAH_113_C': {'range': [11.2, 11.4]},
             'PAH_120': {'range': [11.9, 12.1]},
             'PAH_126_C': {'range': [12.6, 12.7]},
             'PAH_136': {'range': [13.4, 13.6]},
             'PAH_142': {'range': [14.1, 14.2]},
             'PAH_164': {'range': [16.4, 16.5]},
             'PAH_17_C': {'range': [16.4, 17.9]},
             'PAH_174': {'range': [17.35, 17.45]}
             }

    # Create combined strength, unc, and eqw table.
    cftable = Table(
        names=("Name", "range_min", "range_max", "strength", "strength_unc", "eqw"),
        dtype=("U25", "float64", "float64", "float64", "float64", "float64"))

    # Get indices of dust features.
    df_ind = np.concatenate(np.argwhere(ftable["Form"] == "Drude1D"))

    # Get table subset containing only the dust features.
    dfs = ftable[df_ind[0]:df_ind[-1] + 1]

    # Iterate keys and calculate combined strengths.
    for feat in cfdic.keys():
        mask = np.logical_and(dfs['x_0'] >= cfdic[feat]['range'][0], dfs['x_0'] <= cfdic[feat]['range'][-1])
        cftable.add_row([feat,
                         cfdic[feat]['range'][0],
                         cfdic[feat]['range'][-1],
                         np.sum(dfs[mask]['strength']),
                         None,
                         np.sum(dfs[mask]['eqw'])
                         ])

    return cftable


def eqws(comp_type, x_0, amp, fwhm_stddev, obs_fit):
    """
    Calculate the emission features equivalent width
    (integral[(I_nu-I_cont)/I_cont d_lam]) in microns.

    Parameters
    ----------
    comp_type : string
        type of emission component (Drude1D/Gaussian)
    x_0 : float
        central wavelength of the feature.
    amp : float
        central intensity of the feature.
    fwhm_stddev : float
        fwhm or stddev of the feature depending on comp_type.
    Returns
    -------
    eqw : float
        the equivalent width of the feature
    """
    # Check if the emission component is Gaussian and calculate fwhm.
    if comp_type == 'Gaussian1D':
        fwhm = 2 * fwhm_stddev * np.sqrt(2 * np.log(2))
    else:
        fwhm = fwhm_stddev

    # Get range and wavelength region for integration.
    low = x_0 - (fwhm * 6)
    lmin = low if low > 0 else 0.
    lmax = x_0 + (fwhm * 6)
    # lam = np.arange(100) / 99 * (lmax - lmin) + lmin
    lam = np.linspace(lmin, lmax, num=100)

    # Calculate the continuum and feature components in the integration range.
    cont_components = []
    for cmodel in obs_fit:
        if isinstance(cmodel, BlackBody1D):
            cont_components.append(cmodel)
    cont_model = cont_components[0]
    for cmodel in cont_components[1:]:
        cont_model += cmodel
    continuum = np.nan_to_num(cont_model(lam))

    if comp_type == 'Drude1D':
        drude = Drude1D(amplitude=amp,
                        x_0=x_0,
                        fwhm=fwhm)
        lnu = drude(lam)

    elif comp_type == 'Gaussian1D':
        gauss = Gaussian1D(amplitude=amp,
                           mean=x_0,
                           stddev=fwhm_stddev)
        lnu = gauss(lam)

    # Following the EQW calculation in IDL PAHFIT, we define a default broad limit (bl).
    # Set bl to replace the continuum in the EQW calculation
    # by its profile-averaged value for fractional FWHM greater than that limit.
    # Useful when the EQW requires extrapolation to regions where the continuum
    # can vanish, such that f_line/f_continuum diverges.
    # Default value bl = 0.05.
    bl = 0.05

    # Calculate EQW.
    if fwhm / x_0 > bl:
        ilam = integrate.simpson(lnu, lam)
        weighted_cont = integrate.simpson(lnu * continuum, lam) / ilam
        eqw = ilam / weighted_cont
    else:
        eqw = integrate.simpson(lnu / continuum, lam)

    return eqw
