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
       thinking about redshift, units, instrumental effects.)
    3. Retrieve the fitted quantities, which are the values that were
       passed during step 1. When fit result uncertainties are
       implemented, they will also need to be retrieved through this
       API.
    4. Access to the evaluation of the underlying model (again with no
       assumptions like in step 2.).

    For the model setup, multiple functions are used, so a few notes are
    provided here. There is one function per type of component supported
    by PAHFIT, and the arguments of these functions will ask for
    different "standard" PAHFIT numbers, i.e. those from the Features
    table. These functions have the same signature between all Fitter
    implementations, so that the Model class can use a single
    implementation to set up the Fitter. The Model has access to the
    Features table and the instrument model, and needs to set up Fitter
    with the correct initial values, bounds, and "fixed" flags (e.g.
    setting a fixed FWHM based on the instrument for the lines). After
    all the components have been added, the finalize_model() function
    can be called to finish setting up the internal astropy model. After
    this has finished, fit() can be called to apply the model and the
    astropy fitter to the data.

    """

    @abstractmethod
    def clear(self):
        """Reset model.

        After reset, register_() and finalize_model() can be used again.

        """
        pass

    @abstractmethod
    def components(self):
        """Return list of features.

        Only works after finalize_model(). Will return the names passed
        using the register functions.

        """
        pass

    @abstractmethod
    def finalize_model(self):
        """Process the registered features and prepare for fitting.

        The register functions below allow adding individual features.
        The exact implementation of how features are added, and
        finalized in to a single fittable model, depend on the
        underlying implementation.

        """
        pass

    @abstractmethod
    def register_starlight(self, name, temperature, tau):
        """Register a starlight feature.

        The exact representation depends on the implementation, but the
        meaning of the parameters should be equivalent.

        Parameters
        ----------
        name : str
            Unique name. Will be used to allow retrieval of the results
            after the fitting.

        temperature : array of size 3
            Blackbody temperature, given as [value, lower_bound,
            upper_bound]. The bounds assume the same system as the
            features table: if they are masked, the parameter will be
            fixed, while +inf or -inf mean unbounded.

        tau : array of size 3
            Analogously, used as power.

        """
        pass

    @abstractmethod
    def register_dust_continuum(self, name, temperature, tau):
        """Register a dust continuum feature."""
        pass

    @abstractmethod
    def register_line(self, name, power, wavelength, fwhm):
        """Register an emission line feature.

        Typically a Gaussian profile.

        """
        pass

    @abstractmethod
    def register_dust_feature(self, name, power, wavelength, fwhm):
        """Register a dust feature.

        Typically a Drude profile.

        """
        pass

    @abstractmethod
    def register_attenuation(self, name, tau):
        """Register the S07 attenuation component.

        Other types of attenuation might be possible in the future. Is
        multiplicative.

        """
        pass

    @abstractmethod
    def register_absorption(self, name, tau, wavelength, fwhm):
        """Register an absorption feature.

        Typically a Drude profile. Is multiplicative.

        """
        pass

    @abstractmethod
    def evaluate_model(self, xz):
        """Evaluate the model at the given wavelengths.

        Parameters
        ----------
        xz : array
            Rest frame wavelengths in micron

        Returns
        -------
        yz : array
            Rest frame flux in internal units

        """
        pass

    @abstractmethod
    def fit(self, xz, yz, uncz, maxiter=1000):
        """Perform the fit using the framework of the subclass.

        Fitter is unit agnostic, and deals with the numbers the Model
        tells it to deal with. In practice, the input spectrum is
        expected to be in internal units, and corrected for redshift
        (models operate in the rest frame).

        After the fit, the results can be retrieved via get_result().

        Parameters
        ----------
        xz : array
            Rest frame wavelengths in micron

        yz : array
            Rest frame flux in internal units.

        uncz : array
            Uncertainty on rest frame flux. Same units as yz.

        """
        pass

    @abstractmethod
    def get_result(self, feature_name):
        """Retrieve results from underlying model after fit.

        Parameters
        ----------
        component_name : str
            One of the names provided to any of the register_() calls
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
