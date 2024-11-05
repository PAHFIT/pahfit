from abc import ABC, abstractmethod


class Fitter(ABC):
    """Abstract base class for internal Fitter API.

    All shared methods should have the same arguments, enforced by this
    abstract class. Any API-specific options preferably go into the
    constructor of the subclass, although some general-purpose
    dictionaries could also be used if absolutely necessary.

    The main functionalities of a Fitter subclass:

    1. Convert the numbers that are in the Features table to a fittable
       model configuration for a certain framework. The details of the
       fitting framework are hidden behind the respective subclass.

    2. Fit the model to the spectrum without any additional
       assumptions.  The Fitter will fit the given data using the
       given model without needing to be aware of redshift, units, or
       other instrumental effects.
    
    3. Retrieve the fitted quantities, which are the values that were
       passed during step 1.  When fit uncertainties are implemented,
       they will also need to be retrieved through this API.

    4. Access to the evaluation of the underlying model (again with no
       assumptions like in step 2.).

    A few notes on how the above is achieved:

    For the model setup, there is one function per type of component
    supported by PAHFIT, and the arguments of these functions will ask
    for certain standard parameters (in practice, these are the values
    stored in the Features table). This abstract Fitter class ensures
    that the function signatures are the same between different Fitter
    implementations, so that uniform logic can be implemented outside
    of the Fitter (in practice, this is a loop over the Features table
    implemented in :class:`pahfit.model.Model`).

    During the Fitter setup, the initial values, bounds, and "fixed"
    flags are passed using one function call for each component, e.g.
    :meth:`~pahfit.fitters.Fitter.add_feature_line`.  Once all
    components have been added, the
    :meth:`~pahfit.fitters.Fitter.finalize` function should be called;
    some subclasses (e.g. :class:`pahfit.fitters.APFitter`) need to
    consolidate the registered components to prepare the model that
    they manage for fitting. After this,
    :meth:`~pahfit.fitters.Fitter.fit` can be called to apply the
    model and the fitter to the data. The results will then be
    retrievable for one component at a time, by passing the component
    name to get_result().
    """

    @abstractmethod
    def finalize(self):
        """Process the registered features and prepare for fitting.

        The register functions below allow adding individual features.
        The exact implementation of how features are added, and
        finalized in to a single fittable model, depend on the
        underlying implementation.

        """
        pass

    @abstractmethod
    def add_feature_starlight(self, name, temperature, tau):
        """Register a starlight feature.

        The exact representation depends on the implementation, but the
        meaning of the parameters should be equivalent.

        Parameters
        ----------
        name : str
            Unique name. Will be used to allow retrieval of the results
            after the fitting.

        temperature : array of size 3 or scalar
            Blackbody temperature. Given as [value, lower_bound,
            upper_bound] if the parameter should be variable (and
            fitted). Given as scalar if parameter should be fixed.

        tau : array of size 3
            Analogously, used as power.

        """
        pass

    @abstractmethod
    def add_feature_dust_continuum(self, name, temperature, tau):
        """Register a dust continuum feature."""
        pass

    @abstractmethod
    def add_feature_line(self, name, power, wavelength, fwhm):
        """Register an emission line feature.

        Typically a Gaussian profile.

        """
        pass

    @abstractmethod
    def add_feature_dust_feature(self, name, power, wavelength, fwhm):
        """Register a dust feature.

        Typically a Drude profile.

        """
        pass

    @abstractmethod
    def add_feature_attenuation(self, name, tau, model="S07", geometry="screen"):
        """Register the S07 attenuation component.

        Other types of attenuation might be possible in the
        future. Multiplicative.
        """
        pass

    @abstractmethod
    def add_feature_absorption(self, name, tau, wavelength, fwhm, geometry="screen"):
        """Register an absorption feature.

        Modeled by a Drude profile.  Multiplicative.

        """
        pass

    @abstractmethod
    def evaluate(self, lam):
        """Evaluate the fitting function at the given wavelengths.

        Parameters
        ----------
        lam : array
            Rest frame wavelengths in micron

        Returns
        -------
        flux : array
            Rest frame flux in internal units

        """
        pass

    @abstractmethod
    def fit(self, lam, flux, unc, maxiter=1000):
        """Perform the fit using the framework of the subclass.

        :class:`~pahfit.fitters.Fitter` is unit agnostic, and deals
        with the numbers the :class:`~pahfit.model.Model` tells it to
        deal with.  In practice, the input spectrum should be expected
        in internal units (see :mod:`pahfit.units`), and corrected for
        redshift (models operate in the rest frame).

        After the fit, the results can be retrieved via get_result().

        Parameters
        ----------
        lam : array
            Rest frame wavelengths in microns.

        flux : array
            Rest frame flux in internal units.

        unc : array
            Uncertainty in the rest-frame flux.  Same units as flux.
        """
        pass

    @abstractmethod
    def get_result(self, feature_name):
        """Retrieve results from the underlying model after fit.

        Parameters
        ----------
        component_name : str
            One of the names provided to any of the
            :meth:`~pahfit.fitters.Fitter.add_feature` calls made
            during setup.

        Returns
        -------
        feature_info : dict
           feature parameters according to the relevant PAHFIT
           definitions. Key names are the same as the function
           arguments of the relevant register function. Values are in
           the same format as :class:`~pahfit.feature.Features`, and
           can therefore be directly filled in.

           E.g. ``{'name': 'line0', 'power': value, 'fwhm': value,
        'wavelength': value}``.
        """
        pass
