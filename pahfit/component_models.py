import numpy as np
import pkg_resources
from scipy import interpolate
from astropy.modeling.physical_models import Drude1D
from astropy.modeling import Fittable1DModel
from astropy.modeling import Parameter
from astropy.table import Table

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

        kvt_prof = Table.read(pkg_resources.resource_filename("pahfit", "packs/")+'att_profile/kvt_Smith2007.ecsv')

        f = interpolate.interp1d(kvt_prof['lam'], kvt_prof['optical_depth'])

        new_f = f(in_x)

        return new_f
     

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
            return np.full((len(x)), 1.0)
        else:
            profile = Drude1D(amplitude=1.0, fwhm=fwhm, x_0=x_0)
            tau_x = tau * profile(x)
            return (1.0 - np.exp(-1.0 * tau_x)) / tau_x
