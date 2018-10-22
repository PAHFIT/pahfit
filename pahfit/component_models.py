from collections import OrderedDict

import numpy as np

from astropy.modeling import Fittable1DModel
from astropy.modeling import Parameter


__all__ = ['BlackBody1D', 'Drude1D']


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
        # One dimensional Lorentzian model derivative with respect to parameters
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
