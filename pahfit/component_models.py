import numpy as np
from scipy import interpolate
from astropy.modeling.physical_models import Drude1D
from astropy.modeling import Fittable1DModel
from astropy.modeling import Parameter


__all__ = ["BlackBody1D", "S07_attenuation"]


class BlackBody1D(Fittable1DModel):
    """
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
            * ((9.7 / x) ** 2)
            * 3.97289e13
            / x ** 3
            / (np.exp(1.4387752e4 / x / temperature) - 1.0)
        )


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
        kvt_wav_short = in_x[in_x < min(kvt_wav)]
        kvt_int_short_tmp = min(kvt_int) * np.exp(2.03 * (kvt_wav_short - min(kvt_wav)))
        # Since kvt_int_shoft_tmp does not reach min(kvt_int),
        # we scale it to stitch it.
        kvt_int_short = kvt_int_short_tmp * (kvt_int[0] / max(kvt_int_short_tmp))

        spline_x = np.concatenate([kvt_wav_short, kvt_wav])
        spline_y = np.concatenate([kvt_int_short, kvt_int])

        # spline_rep = interpolate.splrep(spline_x, spline_y)
        intfunc = interpolate.interp1d(spline_x, spline_y)
        in_x_spline = in_x[in_x < max(kvt_wav)]
        # new_spline_y = interpolate.splev(in_x_spline, spline_rep, der=0)
        new_spline_y = intfunc(in_x_spline)

        nf = Drude1D(amplitude=0.4, x_0=18.0, fwhm=0.247 * 18.0)
        in_x_drude = in_x[in_x >= max(kvt_wav)]

        ext = np.concatenate([new_spline_y, nf(in_x_drude)])

        # Extend to ~2 um
        # assuing beta is 0.1
        beta = 0.1
        y = (1.0 - beta) * ext + beta * (9.7 / in_x) ** 1.7

        return y

    def evaluate(self, in_x, tau_si):
        if tau_si == 0.0:
            return np.full((len(in_x)), 1.0)
        else:
            tau_x = tau_si * self.kvt(in_x)
            return (1.0 - np.exp(-1.0 * tau_x)) / tau_x
