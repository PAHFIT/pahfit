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
    Base PAHFit model
    proto-type - only blackbodies for continuum right now

    .. plot::
        :include-source:

        import numpy as np
        import matplotlib.pyplot as plt
        import astropy.units as u

        from pypahfit.PAHFit import PAHFit

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
        self._bb_amps_default = np.array([1.0, 1.0, 1.0, 1.0, 1.0,
                                          1.0, 1.0, 1.0])

        pnames = ["bb%i" % (k+1) for k in range(len(self._bb_temps))]
        self._param_names = tuple(pnames)

        super().__init__()

        for k, pname in enumerate(self._param_names):
            param = Parameter(pname, default=0.0, model=self)
            param.__set__(self, self._bb_amps_default[k])

    def __getattr__(self, attr):
        if self._param_names and attr in self._param_names:
            return Parameter(attr, default=0.0, model=self)

        raise AttributeError(attr)

    def __setattr__(self, attr, value):
        # TODO: Support a means of specifying default values for coefficients
        # Check for self._ndim first--if it hasn't been defined then the
        # instance hasn't been initialized yet and self.param_names probably
        # won't work.
        # This has to vaguely duplicate the functionality of
        # Parameter.__set__.
        # TODO: I wonder if there might be a way around that though...
        if attr[0] != '_' and self._param_names and attr in self._param_names:
            param = Parameter(attr, default=0.0, model=self)
            # This is a little hackish, but we can actually reuse the
            # Parameter.__set__ method here
            param.__set__(self, value)
        else:
            super().__setattr__(attr, value)

    def evaluate(self, x, *bbparams):
        """
        Need docs
        """
        pass
