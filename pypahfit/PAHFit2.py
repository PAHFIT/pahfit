# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains models representing polynomials and polynomial series.
"""

from collections import OrderedDict

import numpy as np

from astropy.modeling import FittableModel
from astropy.modeling.parameters import Parameter
from astropy.modeling.utils import poly_map_domain, comb

__all__ = [
    'Polynomial1D'
]


class PolynomialBase(FittableModel):
    """
    Base class for all polynomial-like models with an arbitrary number of
    parameters in the form of coefficients.

    In this case Parameter instances are returned through the class's
    ``__getattr__`` rather than through class descriptors.
    """

    # Default _param_names list; this will be filled in by the implementation's
    # __init__
    _param_names = ()

    linear = True
    col_fit_deriv = False

    @property
    def param_names(self):
        """Coefficient names generated based on the model's polynomial degree
        and number of dimensions.

        Subclasses should implement this to return parameter names in the
        desired format.

        On most `Model` classes this is a class attribute, but for polynomial
        models it is an instance attribute since each polynomial model instance
        can have different parameters depending on the degree of the polynomial
        and the number of dimensions, for example.
        """

        return self._param_names

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


class PolynomialModel(PolynomialBase):
    """
    Base class for polynomial models.

    Its main purpose is to determine how many coefficients are needed
    based on the polynomial order and dimension and to provide their
    default values, names and ordering.
    """

    def __init__(self, degree, n_models=None, model_set_axis=None,
                 name=None, meta=None, **params):
        self._degree = degree

        super().__init__(
            n_models=n_models, model_set_axis=model_set_axis, name=name,
            meta=meta, **params)

    def __repr__(self):
        return self._format_repr([self.degree])

    def __str__(self):
        return self._format_str([('Degree', self.degree)])

    @property
    def degree(self):
        """Degree of polynomial."""

        return self._degree

    def get_num_coeff(self, ndim):
        """
        Return the number of coefficients in one parameter set
        """

        if self.degree < 0:
            raise ValueError("Degree of polynomial must be positive or null")
        # deg+1 is used to account for the difference between iraf using
        # degree and numpy using exact degree
        if ndim != 1:
            nmixed = comb(self.degree, ndim)
        else:
            nmixed = 0
        numc = self.degree * ndim + nmixed + 1
        return numc

    def _invlex(self):
        c = []
        lencoeff = self.degree + 1
        for i in range(lencoeff):
            for j in range(lencoeff):
                if i + j <= self.degree:
                    c.append((j, i))
        return c[::-1]


class PF_BlackBodies(PolynomialModel):
    r"""
    PAHFIT blackbodies component
    """

    inputs = ('x',)
    outputs = ('y',)
    # _separable = True

    def __init__(self,
                 bb_temps=np.array([300., 200., 135., 90.,
                                    65., 50., 40., 35.]),
                 domain=[-1, 1], window=[-1, 1], n_models=None,
                 model_set_axis=None, name=None, meta=None, **params):
        self.domain = domain
        self.window = window

        self._bb_temps = bb_temps

        self._order = len(self._bb_temps)

        self._bb_temps = np.array([300., 200., 135., 90., 65.,
                                  50., 40., 35.])
        # not clear how to get defaults in this setup
        self._bb_amps_default = np.array([1.0, 1.0, 1.0, 1.0, 1.0,
                                          1.0, 1.0, 1.0])

        self._param_names = self._generate_coeff_names()

        super().__init__(
            self._order, n_models=n_models, model_set_axis=model_set_axis,
            name=name, meta=meta, **params)

    def _generate_coeff_names(self):
        names = []
        for cbbtemp in self._bb_temps:
            names.append('bb{:d}'.format(int(cbbtemp)))

        return tuple(names)

    def prepare_inputs(self, x, **kwargs):
        inputs, format_info = super().prepare_inputs(x, **kwargs)

        x = inputs[0]
        return (x,), format_info

    @staticmethod
    def _bb_MJysr(waves, T):
        return 3.97289e13/waves**3/(np.exp(1.4387752e4/waves/T)-1.)

    def evaluate(self, x, *coeffs):
        bbs = np.zeros(len(x))
        for i, ctemp in enumerate(self._bb_temps):
            bbs += coeffs[i]*((9.7/x)**2)*self._bb_MJysr(x, ctemp)

        return bbs

    @property
    def input_units(self):
        if self.degree == 0 or self.c1.unit is None:
            return None
        else:
            return {'x': self.c0.unit / self.c1.unit}

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        mapping = []
        for i in range(self.degree + 1):
            par = getattr(self, 'c{0}'.format(i))
            mapping.append((par.name, outputs_unit['y'] / inputs_unit['x'] ** i))
        return OrderedDict(mapping)
