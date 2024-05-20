from abc import ABC, abstractmethod


class Fitter(ABC):
    """Abstract base class for interal Fitter API.

    All shared methods should have the same arguments, enforced by this
    abstract class. Any API-specific options preferably go into the
    constructor of the subclass, although some general-purpose
    dictionaries could also be used if absolutely necessary.

    The main functionalities of a Fitter subclass:
    1. Convert the numbers that are in the Features table to a fittable
       model configuration for a certain framework. The details of the
       fitting framework are hidden behind the respective subclass.
    2. Fit the model to the spectrum without any additional assumptions.
       The Fitter will fit the given data using the given model without
       thinking about redshift, units, instrumental effects.
    3. Retrieve the fitted quantities, which are the values that were
       passed during step 1. When fit result uncertainties are
       implemented, they will also need to be retrieved through this
       API.
    4. Access to the evaluation of the underlying model (again with no
       assumptions like in step 2.).

    A few notes on how the above is achieved:

    For the model setup, there is one function per type of component
    supported by PAHFIT, and the arguments of these functions will ask
    for certain standard PAHFIT quantities (in practice, these are the
    values stored in the Features table). The abstract Fitter class
    ensure that the function signatures are the same between different
    Fitter implementations, so that only a single logic has to be
    implemented to up the Fitter (in practice, this is a loop over the
    Features table implemented in Model).

    During the Fitter setup, the initial values, bounds, and "fixed"
    flags are passed using one function call for each component, e.g.
    add_feature_line()). Once all components have been added, the
    finalize() function should be called; some subclasses (e.g.
    APFitter) need to consolidate the registered components to prepare
    the model that they manage for fitting. After this, fit() can be
    called to apply the model and the astropy fitter to the data. The
    results will then be retrievable for one component at a time, by
    passing the component name to get_result().

    """

    @abstractmethod
    def components(self):
        """Return list of features.

        Only works after finalize(). Will return the names passed
        using the register functions.

        """
        pass

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
    def add_feature_attenuation(self, name, tau, model='S07', geometry='screen'):
        """Register the S07 attenuation component.

        Other types of attenuation might be possible in the future. Is
        multiplicative.

        """
        pass

    @abstractmethod
    def add_feature_absorption(self, name, tau, wavelength, fwhm, geometry='screen'):
        """Register an absorption feature.

        Modeled by a Drude profile. Is multiplicative.

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

        Fitter is unit agnostic, and deals with the numbers the Model
        tells it to deal with. In practice, the input spectrum is
        expected to be in internal units, and corrected for redshift
        (models operate in the rest frame).

        After the fit, the results can be retrieved via get_result().

        Parameters
        ----------
        lam : array
            Rest frame wavelengths in micron

        flux : array
            Rest frame flux in internal units.

        unc : array
            Uncertainty on rest frame flux. Same units as flux.

        """
        pass

    @abstractmethod
    def get_result(self, feature_name):
        """Retrieve results from underlying model after fit.

        Parameters
        ----------
        component_name : str
            One of the names provided to any of the add_feature_() calls
            made during setup. See also Fitter.components().

        Returns
        -------
        dict : parameters according to the PAHFIT definitions. Keys are
        the same as the function signature of the relevant register
        function. Values are in the same format as Features, and can
        therefore be directly filled in.

        e.g. {'name': 'line0', 'power': value, 'fwhm': value, 'wavelength'}

        """
        pass
