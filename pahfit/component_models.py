from collections import OrderedDict

import numpy as np
from scipy import interpolate

from astropy.modeling import Fittable1DModel
from astropy.modeling import Parameter


__all__ = ['BlackBody1D', 'Drude1D', 'SD07_attenuation']


class BlackBody1D(Fittable1DModel):
    """
    Current astropy BlackBody1D does not play well with Lorentz1D and Gauss1D
    maybe...need to check again...possibly a units issue
    """
    inputs = ('x',)
    outputs = ('y',)

    amplitude = Parameter()
    temperature = Parameter()

    @staticmethod
    def evaluate(x, amplitude, temperature):
        """
        """
        return (amplitude*((9.7/x)**2)*3.97289e13/x**3
                / (np.exp(1.4387752e4/x/temperature)-1.))


class Drude1D(Fittable1DModel):
    """
    One dimensional Drude model.

    Parameters
    ----------
    amplitude : float
        Peak value
    x_0 : float
        Position of the peak
    fwhm : float
        Full width at half maximum

    See Also
    --------
    Gaussian1D, Box1D, MexicanHat1D, Lorentz1D

    Notes
    -----
    Model formula:
    .. math::
        f(x) = \\frac{A \\gamma^{2}}{(x/x_0 - x_0/x)^{2} + \\gamma^{2}}

    Examples
    --------
    .. plot::
        :include-source:
        import numpy as np
        import matplotlib.pyplot as plt
        from models import Drude1D
        plt.figure()
        s1 = Drude1D()
        r = np.arange(-5, 5, .01)
        for factor in range(1, 4):
            s1.amplitude = factor
            plt.plot(r, s1(r), color=str(0.25 * factor), lw=2)
        plt.axis([-5, 5, -1, 4])
        plt.show()
    """

    amplitude = Parameter(default=1)
    x_0 = Parameter(default=0)
    fwhm = Parameter(default=1)

    @staticmethod
    def evaluate(x, amplitude, x_0, fwhm):
        """
        One dimensional Drude model function
        """

        return (amplitude * ((fwhm / x_0) ** 2) / ((x/x_0 - x_0/x) ** 2 +
                                                   (fwhm / x_0) ** 2))

    # @staticmethod
    # def fit_deriv(x, amplitude, x_0, fwhm):
        # """
        # One dimensional Lorentzian model derivative with respect
        # to parameters
        # """

        # d_amplitude = (fwhm/x_0) ** 2 / ((x/x_0 - x_0/x) ** 2 +
        #                                  (fwhm / x_0) ** 2)
        # need updating below (d_amplitude already updated)
        # d_x_0 = (amplitude * d_amplitude * (2 * x - 2 * x_0) /
        #         (fwhm ** 2 + (x - x_0) ** 2))
        # d_fwhm = 2 * amplitude * d_amplitude / fwhm * (1 - d_amplitude)
        # return [d_amplitude, d_x_0, d_fwhm]

    def bounding_box(self, factor=25):
        """Tuple defining the default ``bounding_box`` limits,
        ``(x_low, x_high)``.
        Parameters
        ----------
        factor : float
            The multiple of FWHM used to define the limits.
            Default is chosen to include most (99%) of the
            area under the curve, while still showing the
            central feature of interest.
        """
        x0 = self.x_0
        dx = factor * self.fwhm

        return (x0 - dx, x0 + dx)

    @property
    def input_units(self):
        if self.x_0.unit is None:
            return None
        else:
            return {'x': self.x_0.unit}

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        return OrderedDict([('x_0', inputs_unit['x']),
                            ('fwhm', inputs_unit['x']),
                            ('amplitude', outputs_unit['y'])])


class SD07_attenuation(Fittable1DModel):
    """
    Smith & Draine (2007) kvt attenuation model calculation.
    Calculation is for a fully mixed geometrically model.
    Uses an extinction curve based on the silicate profiles from
    Kemper, Vreind, & Tielens (2004, apJ, 609, 826).
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
    inputs = ('x',)
    outputs = ('y',)

    # Attenuation tau
    tau_si = Parameter(description="kvt term: amplitude",
                       default=1.0, min=0.0, max=10)

    @staticmethod
    def kvt(in_x):
        """
        Output the kvt extinction curve
        """
        kvt_wav = np.array([8.0, 8.2, 8.4, 8.6, 8.8, 9.0, 9.2, 9.4, 9.6, 9.7,
                            9.75, 9.8, 10.0, 10.2, 10.4, 10.6, 10.8, 11.0,
                            11.2, 11.4, 11.6, 11.8, 12.0, 12.2, 12.4, 12.6,
                            12.7])

        kvt_int = np.array([.06, .09, .16, .275, .415, .575, .755, .895, .98,
                            .99, 1.0, .99, .94, .83, .745, .655, .58, .525,
                            .43, .35, .27, .20, .13, .09, .06, .045, .04314])

        # Extend kvt profile to shorter wavelengths
        kvt_wav_short = in_x[in_x < min(kvt_wav)]
        kvt_int_short_tmp = (min(kvt_int)
                             * np.exp(2.03*(kvt_wav_short-min(kvt_wav))))
        # Since kvt_int_shoft_tmp does not reach min(kvt_int),
        # we scale it to stitch it.
        kvt_int_short = kvt_int_short_tmp * (kvt_int[0]/max(kvt_int_short_tmp))

        # Extend kvt profile to longer wavelengths
        kvt_wav_long = in_x[in_x > max(kvt_wav)]
        # Need to input the value instead of using zeros
        kvt_int_long = np.zeros(len(kvt_wav_long))

        spline_x = np.concatenate([kvt_wav_short, kvt_wav, kvt_wav_long])
        spline_y = np.concatenate([kvt_int_short, kvt_int, kvt_int_long])

        spline_rep = interpolate.splrep(spline_x, spline_y)
        new_spline_y = interpolate.splev(in_x, spline_rep, der=0)

        nf = Drude1D(amplitude=0.4,
                     x_0=18.,
                     fwhm=0.247*18.)
        ext = nf(in_x) + new_spline_y

        # Extend to ~2 um
        # assuing beta is 0.1
        beta = 0.1
        y = (1.0-beta)*ext + beta*(9.7/in_x)**1.7

        return y

    def evaluate(self, in_x, tau_si):
        if tau_si == 0.0:
            return np.full((len(in_x)), 0.0)
        else:
            tau_x = tau_si*self.kvt(in_x)
            return (1.0 - np.exp(-1.0*tau_x))/tau_x
