from __future__ import (absolute_import, print_function, division)

import numpy as np

import astropy.units as u
from astropy.modeling import (Fittable1DModel,
                              Parameter, InputParameterError)

__all__ = ['PAHFit']


class PAHFit(Fittable1DModel):
    """
    PAHFit model calcultion

    Parameters
    ----------


    Notes
    -----
    FM90 extinction model

    From Fitzpatrick & Massa (1990)

    Only applicable at UV wavelengths

    Example showing a FM90 curve with components identified.

    .. plot::
        :include-source:

        import numpy as np
        import matplotlib.pyplot as plt
        import astropy.units as u

        from pypahfit.pahfit import PAHFit

        fig, ax = plt.subplots()

        # generate the curves and plot them
        x = np.arange(5.0,40.0,0.1)*u.micron

        mir_model = PAHFit()
        ax.plot(x,mir_model(x),label='total')

        ax.legend(loc='best')
        plt.show()
    """

    def __init__(self):

        self._bb_temps = np.array([300., 200., 135., 90., 65.,
                                  50., 40., 35.])

        pnames = ["bb".format(k+1) for k in range(self._bb_temps)]
        self._param_names = tuple(pnames)

    def evaluate(self, x, *bbparams):
        """
        Need docs
        """
        pass
