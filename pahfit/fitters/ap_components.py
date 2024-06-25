import numpy as np
from scipy import interpolate
from astropy.modeling.physical_models import Drude1D
from astropy.modeling import Fittable1DModel
from astropy.modeling import Parameter
from astropy import constants
from pahfit import units

__all__ = ["BlackBody1D", "ModifiedBlackBody1D", "S07_attenuation", "att_Drude1D"]


class BlackBody1D(Fittable1DModel):
    """
    A blackbody component.

    Current astropy BlackBody1D does not play well with Lorentz1D and Gauss1D
    maybe, need to check again, possibly a units issue
    """

    amplitude = Parameter()
    temperature = Parameter()

    @staticmethod
    def evaluate(x, amplitude, temperature):
        """
        """
        return (
            amplitude
            * 3.97289e13
            / x ** 3
            / (np.exp(1.4387752e4 / x / temperature) - 1.0)
        )


class ModifiedBlackBody1D(BlackBody1D):
    """
    Modified blackbody with an emissivity propoportional to nu^2
    """

    @staticmethod
    def evaluate(x, amplitude, temperature):
        return BlackBody1D.evaluate(x, amplitude, temperature) * ((9.7 / x) ** 2)


class S07_attenuation(Fittable1DModel):
    """
    Smith, Draine, et al. (2007) kvt attenuation model calculation.
    Calculation is for a fully mixed geometrically model.
    Uses an extinction curve based on the silicate profiles from
    Kemper, Vriend, & Tielens (2004, apJ, 609, 826).
    Constructed as a weighted sum of two components: silicate profiles,
    and an exponent 1.7 power-law.

    Attenuation curve for a mixed case calculated from
    .. math::

        Att(x) = \\frac{1 - e^{-\\tau_{x}}}{\\tau_{x}}

    Parameters
    ----------
    kvt_amp : float
      amplitude
    """

    # Attenuation tau
    tau_sil = Parameter(description="kvt term: amplitude", default=1.0, min=0.0, max=10)

    @staticmethod
    def kvt(in_x):
        """
        Output the kvt extinction curve
        """
        # fmt: off
        kvt_wav = np.array([8.0, 8.2, 8.4, 8.6, 8.8, 9.0, 9.2, 9.4, 9.6, 9.7,
                            9.75, 9.8, 10.0, 10.2, 10.4, 10.6, 10.8, 11.0,
                            11.2, 11.4, 11.6, 11.8, 12.0, 12.2, 12.4, 12.6,
                            12.7])

        kvt_int = np.array([.06, .09, .16, .275, .415, .575, .755, .895, .98,
                            .99, 1.0, .99, .94, .83, .745, .655, .58, .525,
                            .43, .35, .27, .20, .13, .09, .06, .045, .04314])
        # fmt: on

        # Extend kvt profile to shorter wavelengths
        if min(in_x) < min(kvt_wav):
            kvt_wav_short = in_x[in_x < min(kvt_wav)]
            kvt_int_short_tmp = min(kvt_int) * np.exp(2.03 * (kvt_wav_short - min(kvt_wav)))
            # Since kvt_int_shoft_tmp does not reach min(kvt_int),
            # we scale it to stitch it.
            kvt_int_short = kvt_int_short_tmp * (kvt_int[0] / max(kvt_int_short_tmp))

            spline_x = np.concatenate([kvt_wav_short, kvt_wav])
            spline_y = np.concatenate([kvt_int_short, kvt_int])
        else:
            spline_x = kvt_wav
            spline_y = kvt_int

        intfunc = interpolate.interp1d(spline_x, spline_y)
        in_x_spline = in_x[in_x < max(kvt_wav)]
        new_spline_y = intfunc(in_x_spline)

        nf = Drude1D(amplitude=0.4, x_0=18.0, fwhm=0.247 * 18.0)
        in_x_drude = in_x[in_x >= max(kvt_wav)]

        ext = np.concatenate([new_spline_y, nf(in_x_drude)])

        # Extend to ~2 um
        # assuming beta is 0.1
        beta = 0.1
        y = (1.0 - beta) * ext + beta * (9.7 / in_x) ** 1.7

        return y

    def evaluate(self, in_x, tau_si):
        if tau_si == 0.0:
            return np.full((len(in_x)), 1.0)
        else:
            tau_x = tau_si * self.kvt(in_x)
            return (1.0 - np.exp(-1.0 * tau_x)) / tau_x


class att_Drude1D(Fittable1DModel):
    """
    Attenuation components that can be parameterized by Drude profiles.
    """

    tau = Parameter()
    x_0 = Parameter()
    fwhm = Parameter()

    @staticmethod
    def evaluate(x, tau, x_0, fwhm):
        if tau == 0.0:
            return np.full((len(x)), 0.0)
        else:
            profile = Drude1D(amplitude=1.0, fwhm=fwhm, x_0=x_0)
            tau_x = tau * profile(x)
            return (1.0 - np.exp(-1.0 * tau_x)) / tau_x


class PowerDrude1D(Fittable1DModel):
    """Drude profile with amplitude determined by power.

    Special attention is needed for the units. An implementation for
    this is 'unitful' because 'power' is defined as integral over
    frequency, while the profile is formulated as a function of
    wavelength. If we assume the output has unit(flux), then without
    additional conversions the power unit will be unit(flux) * frequency
    = unit(flux) * unit(c) / unit(lam).

    Example: if x_0 and fwhm are in micron, and the flux is in MJy / sr,
    then the unit of the fitted power will be MJy sr-1 * unit(c) /
    micron, which differs by a constant factor from MJy sr-1 Hz,
    depending on the chosen unit of c.

    For efficiency and to prevent ambiguity, we assume that all units
    are the internal pahfit units defined in pahfit.units, and
    precalculate a conversion factor.

    TODO: We need to check if the flux is 'intensity' or 'flux_density',
    and assume the power parameter has 'intensity_power' or
    'flux_density_power' units respectively. For now, only intensity is
    supported.

    """

    power = Parameter(min=0.0)
    x_0 = Parameter(min=0.0)
    fwhm = Parameter(default=1, min=0.0)

    # constant factors in the equation to convert power to amplitude of
    # the profile.
    intensity_amplitude_factor = (
        (2 * units.intensity_power * units.wavelength / (constants.c * np.pi))
        .to(units.intensity)
        .value
    )

    @staticmethod
    def evaluate(x, power, x_0, fwhm):
        """Smith, et al. (2007) dust features model. Calculation is for
        a Drude profile (equation in section 4.1.4).

        The intensity profile as a function of
        wavelength is

        Inu(lambda) = (b * g**2) / ((lambda / x0 - x0 / lambda)**2 + g**2)

        With
        b = amplitude (has same unit as flux)
        g = fwhm / x0
        x0 = central wavelength

        The integrated power (Fnu dnu) of the profile is

        P = (pi * c / 2) * (b * g / x0)

        Which can be solved for the amplitude b.

        b = (P * 2 * x0) / (pi * c * g) = 2P / (pi nu0 g).

        According to the above equations, without additional
        conversions, the resulting amplitude unit will be unit(P) *
        Hz-1. This will result in very small values for for Inu(lambda),
        or very large values for P. To avoid numerical problems with the
        fitting algorithm, we apply conversions so that Inu(lambda) and
        P are in the internal units.

        Parameters
        ----------
        power : float
        fwhm : float
        central intensity (x_0) : float

        """
        # The equation and unit conversion for the amplitude:
        # b = (2 * P * x_0 / (pi * c * g)).to(output_unit).value

        # Use predetermined factor that deals with units and constants.
        # factor = (2 * unit(power) * unit(wavelength) / (pi * c)).to(unit(intensity))

        g = fwhm / x_0
        b = power * x_0 / g * PowerDrude1D.intensity_amplitude_factor
        return b * g**2 / ((x / x_0 - x_0 / x) ** 2 + g**2)


class PowerGaussian1D(Fittable1DModel):
    """Gaussian profile with amplitude derived from power.

    Implementation analogous to PowerDrude1D.

    The amplitude of a gaussian profile given its power P, is P /
    (stddev sqrt(2 pi)). Since stddev is given in wavelength units, this
    equation gives the peak density per wavelength interval, Alambda.
    The profile Flambda is then

    Flambda(lambda) = Alambda * G(lambda; mean, stddev)

    where G is a gaussian profile with amplitude 1.Converting this to
    Fnu units yields

    Fnu(lambda) = lambda**2 / c * Flambda(lambda)

    Approximating the lambda**2 factor as a constant (the central
    wavelength = 'mean'), then yields

    Fnu(lambda) = mean**2 / c * Alambda * G(lambda; mean, stddev)

    In other words, for narrow lines, the per-frequency profile is
    approximately a gaussian with amplitude

    Anu = P * mean**2 / (c * stddev sqrt(2 pi)).

    So the constant factor we can set is
    (unit(power) * unit(wavelength)**2 / (c * unit(wavelength) * sqrt(2 pi))).to(intensity)

    """

    power = Parameter(min=0.0)
    mean = Parameter()
    stddev = Parameter(default=1, min=0.0)

    intensity_amplitude_factor = (
        (
            units.intensity_power
            * (units.wavelength) ** 2
            / (constants.c * units.wavelength * np.sqrt(2 * np.pi))
        )
        .to(units.intensity)
        .value
    )

    @staticmethod
    def evaluate(x, power, mean, stddev):
        """Evaluate F_nu(lambda) given the power.

        See class description for equations and unit notes."""

        # amplitude in intensity units
        Anu = power * mean**2 / stddev * PowerGaussian1D.intensity_amplitude_factor
        return Anu * np.exp(-0.5 * np.square((x - mean) / stddev))
