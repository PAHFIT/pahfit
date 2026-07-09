from . import Fitter
from pahfit.errors import PAHFITModelError
from .ap_components import (
    BlackBody1D,
    ModifiedBlackBody1D,
    S07_attenuation,
    att_Drude1D,
    PowerDrude1D,
    PowerGaussian1D)
from astropy.modeling.fitting import LevMarLSQFitter
from scipy.optimize import least_squares
import numpy as np


class APFitter(Fitter):
    """Astropy fitting implementation using
    :class:`pahfit.fitters.Fitter` API.

    Fitter API is still subject to change. This draft was written with
    what is needed, and the Fitter class was set up based on this draft.

    APFitter implements the Fitter methods using actions that involve an
    astropy-based CompoundModel. Some more implementation-specific
    details for these tasks are given in the documentation of the
    functions.

    Multi-segment fitting was not implemented for this subclass, but
    here is a sketch for how it could work:
    1. During set up, an extra argument is passed, indicating the
       segment number. Each segment gets a different list of components
       to use.
    2. During finalize the joint model / something else for joint
       fitting is set up. Alternatively, all the separate models are
       constructed, and then a set of "unique" parameters is gathered.
       An appropriate name would be the "parameter union" (as it is like
       taking the union of the parameters sets for the individual
       segments).
    3. During the fitting, the objective function (chi2) is the sum of
       the objective functions for the individual models. The arguments
       of this objective function, are the unique parameters that were
       derived above. At every fitting step, these parameters will
       change. Inside the objective function, these changes should be
       propagated to the individual models, after which they can be
       evaluated again. All the fitting algorithm needs, is the
       parameter space, and the objective function value.

    Alternative idea, based on concatenated segments (not spectrally
    merged! Just a raw data concatenation).
    - concatenate x
    - concatenate y
    - bookkeep the boundaries between the concatenations
    - when evaluating, fill the concatenated model spectrum using a
      model that depends on the index (according to the boundaries that
      were set up)
    - Use the 'tied' option for parameters that should be equivalent
      (e.g. starlight in segment 1 should have same temperature as
      starlight in segment 2).
    """

    def __init__(self):
        """Construct a new fitter.

        After construction, use the add_feature_() functions to start
        setting up a model, then call finalize().

        """
        self.additive_components = []
        self.multiplicative_components = []
        self.feature_types = {}
        self.model = None
        self.message = None
        self.methods = ["lm", "trf"]

    def finalize(self):
        """Sum the registered components into one CompoundModel.

        To be called after a series of "add_feature_()" calls, and before
        using the other functionality (fitting and evaluating).

        """
        if len(self.additive_components) > 1:
            self.model = sum(self.additive_components[1:], self.additive_components[0])
        elif len(self.additive_components) == 1:
            self.model = self.additive_components[0]
        else:
            raise PAHFITModelError("No components were set up for APFitter!")

        for c in self.multiplicative_components:
            self.model *= c

    def _add_component(self, astropy_model_class, multiplicative=False, **kwargs):
        """Generically add any feature as an astropy model component.

        To be finalized with finalize()

        Parameters
        ----------
        astropy_model_class : class
            Class symbol, e.g. Blackbody1D.

        multiplicative : bool
            During finalize_model(), all components with multiplicative
            == False will be added together, and then the model will be
            multiplied with the components with multiplicative == True.

        kwargs : dict
            Arguments for the astropy model constructor, including a
            unique value for "name". Should be generated with
            self._astropy_model_kwargs; the add_feature_() functions show how
            to do this for each type of feature.

        """
        if multiplicative:
            self.multiplicative_components.append(astropy_model_class(**kwargs))
        else:
            self.additive_components.append(astropy_model_class(**kwargs))

    def add_feature_starlight(self, name, temperature, tau):
        """Register a BlackBody1D.

        Parameters
        ----------
        name : str
            Unique name.

        temperature : array of size 3
            Temperature for the blackbody, given as [value, lower_bound,
            upper_bound]. The bounds assume the same system as the
            features table, where masked == fixed, and +inf or -inf mean
            unbounded.

        tau : analogous. Used as amplitude.

        """
        self.feature_types[name] = "starlight"
        kwargs = self._astropy_model_kwargs(
            name, ["temperature", "amplitude"], [temperature, tau]
        )
        self._add_component(BlackBody1D, **kwargs)

    def add_feature_dust_continuum(self, name, temperature, tau, model="default"):
        """Register a ModifiedBlackBody1D.

        Analogous. Temperature and tau are used as temperature and
        amplitude

        Parameters
        ----------
        model : str (available values: 'default' and 'special')
            Keyword to activate alternate continuum shapes. 'default' is
            a modified blackbody with a lambda**-2 factor. With
            model='special', the modification factor is the MW average
            extinction law for Rv=5.5.

        """
        self.feature_types[name] = "dust_continuum"
        kwargs = self._astropy_model_kwargs(
            name, ["temperature", "amplitude"], [temperature, tau])
        ap_class = ModifiedBlackBody1D
        self._add_component(ap_class, **kwargs)

    def add_feature_line(self, name, power, wavelength, fwhm):
        """Register a PowerGaussian1D

        Analogous. Uses an implementation of the Gaussian profile, that
        directly fits the power based on the internal PAHFIT units.

        """
        self.feature_types[name] = "line"

        kwargs = self._astropy_model_kwargs(
            name,
            ["power", "mean", "stddev"],
            [power, wavelength, fwhm / 2.355],
        )
        self._add_component(PowerGaussian1D, **kwargs)

    def add_feature_dust_feature(self, name, power, wavelength, fwhm):
        """Register a PowerDrude1D.

        Analogous. Uses an implementation of the Drude profile that
        directly fits the power based on the internal PAHFIT units.

        """
        self.feature_types[name] = "dust_feature"
        kwargs = self._astropy_model_kwargs(
            name, ["power", "x_0", "fwhm"], [power, wavelength, fwhm]
        )
        self._add_component(PowerDrude1D, **kwargs)

    def add_feature_attenuation(self, name, tau, model="S07", geometry="screen"):
        """Register the S07 attenuation component.

        Analogous. Uses tau as tau_sil for S07_attenuation. Is
        multiplicative.

        model : string to select attenuation shape. Only 'S07' is supported for now.

        geometry : string to select different geometries. Only 'screen'
        is available for now.

        """
        self.feature_types[name] = "attenuation"
        kwargs = self._astropy_model_kwargs(name, ["tau_sil"], [tau])
        self._add_component(S07_attenuation, multiplicative=True, **kwargs)

    def add_feature_absorption(self, name, tau, wavelength, fwhm, geometry="screen"):
        """Register an absorbing Drude1D component.

        Analogous. Is multiplicative.

        """
        self.feature_types[name] = "absorption"
        kwargs = self._astropy_model_kwargs(name, ["tau", "x_0", "fwhm"], [tau, wavelength, fwhm])
        self._add_component(att_Drude1D, multiplicative=True, **kwargs)


    def evaluate(self, lam):
        """Evaluate internal astropy model with its current parameters.

        Parameters
        ----------
        lam : array
            Rest frame wavelengths in micron

        Returns
        -------
        flux : array
            Rest frame flux in internal units
        """
        return self.model(lam)


    def _free_parameter_state(self):
        free_param_names = []
        x0 = []
        lower_bounds = []
        upper_bounds = []

        fixed = getattr(self.model, "fixed", {}) 
        tied = getattr(self.model, "tied", {}) 
        bounds = getattr(self.model, "bounds", {}) 

        for pname in self.model.param_names:
            if fixed.get(pname, False):      
               continue
            if tied.get(pname, False):
               continue

            lo_bound, hi_bound = bounds.get(pname, (None, None)) 

            if lo_bound is None or np.ma.is_masked(lo_bound):
                lo_bound = -np.inf
            else:
                lo_bound = float(lo_bound)

            if hi_bound is None or np.ma.is_masked(hi_bound):
                hi_bound = np.inf
            else:
                hi_bound = float(hi_bound)

            if np.isfinite(lo_bound) and np.isfinite(hi_bound) and lo_bound == hi_bound:
                getattr(self.model, pname).value = lo_bound
                continue

            value = np.asarray(getattr(self.model, pname).value, dtype=float).item()

            if np.isfinite(lo_bound) and value <= lo_bound:
                step = 1e-12 * max(1.0, abs(lo_bound))
                value = lo_bound + step
                if np.isfinite(hi_bound) and value >= hi_bound:
                   value = 0.5 * (lo_bound + hi_bound)
	
            if np.isfinite(hi_bound) and value >= hi_bound:
                step = 1e-12 * max(1.0, abs(hi_bound))
                value = hi_bound - step
                if np.isfinite(lo_bound) and value <= lo_bound:
                   value = 0.5 * (lo_bound +hi_bound)

            free_param_names.append(pname)
            x0.append(value)
            lower_bounds.append(lo_bound)
            upper_bounds.append(hi_bound)

        return (free_param_names, np.asarray(x0, dtype=float),
            (np.asarray(lower_bounds, dtype=float),
             np.asarray(upper_bounds, dtype=float)))

    def _set_free_parameters(self, free_param_names, values):
        for pname, value in zip(free_param_names, values):
            getattr(self.model, pname).value = value


    def fit_methods_available(self):
        return self.methods

    def fit(self, lam, flux, unc, maxiter=10000, method=None, x_scale="jac", acc=1e-7, epsilon=1e-7):
        """Fit the internal model using.

        "lm" uses Astropy LevMarLSQFitter

        "trf" calls scipy.optimize.least_squares directly because Astropy
        TRFLSQFitter does not expose x_scale. Here, x_scale="jac" is 
        used because the result parameters can differ by orders of magnitude.
        For example, line power and dust-feature powers can have very different
        numerical scales. Without scaling, TRF can take ineffective trust-region
        steps and terminate before small-scale parameters move correctly.

        The fitter class is unit agnostic, and deal with the numbers the
        Model tells it to deal with. Internal renormalizations could be
        good to consider, as long as any values are converted back to
        the original system before returning them. In practice, the
        input spectrum is expected to be in internal units, and orrected
        for redshift (models operate in the rest frame).

        After the fit, the results can be retrieved via get_result().

        Retrieval of uncertainties and fit details is yet to be
        implemented throught fit_info. The LM path uses Astropy fit_info.
        The TRF path stores the scipy OptimizeResult fields needed for
        diagnostics, chi2, redchi, and dof.

        CAVEAT: flux unit (flux) is still ambiguous, since it can be
        flux density or intensity, according to the options defined in
        pahfit.units. After the fit, the return units of "power" in
        get_results depend on the given spectrum (they will be flux unit
        times wavelength unit).

        Parameters
        ----------
        lam : array
            Rest frame wavelengths in micron

        flux : array
            Rest frame flux in internal units.

        unc : array
            Uncertainty on rest frame flux. Same units as flux.

        maxiter : int
            Maximum number of fitting iterations or function evaluations.

        method : str
            Fit method. Available options: "lm" uses Astropy
            LevMarLSQFitter. "trf" uses scipy.optimize.least_squares
            TRFLSQFitter.

        x_scale : str (or array)
            Scaling passed to scipy.optimize.least_squares for the TRF
            path. The defualt is "jac", which updates parameter scales
            using the inverse norms of the Jacobian columns.

        acc : float
            Relative error desired in the approximate solution. This is
            passed as acc to the Astropy LM fitter and as xtol to scipy
            least_squares for TRF.

        epsilon : float
            Step_size control for numerical Jacobian estimation. The TRF
            path follows Astropy's convention and passes diff_step=npsqrt(epsilon)
            to scipy.optimize.least_squares.
        """
        w = 1 / unc
        mask = np.isfinite(lam) & np.isfinite(flux) & np.isfinite(w)

        self.fit_info = {}

        # select method based on given string or pick default if None
        method_str = self.methods[0] if method is None else method
        if method_str not in self.methods:
            raise PAHFITModelError(
                f"Selected method {method} not available for APFitter backend."
            )

        if method_str == "lm":
            print("Running lmlsq fit")

            fit = LevMarLSQFitter(calc_uncertainties=True)
            temp_result = fit(
                self.model,
                lam[mask],
                flux[mask],
                weights=w[mask],
                maxiter=maxiter,
                epsilon=epsilon,
                acc=acc)

            self.model = temp_result
            self.fit_info = fit.fit_info
            self.message = fit.fit_info["message"]

            model_y = self.model(lam[mask])
            resid = (flux[mask] - model_y) * w[mask]
            chi2 = np.nansum(resid**2)

            free_param_names, _, _ = self._free_parameter_state()
            dof = max(np.count_nonzero(mask) - len(free_param_names), 1)
            redchi = chi2 / dof

            self.fit_info["chi2"] = chi2
            self.fit_info["redchi"] = redchi
            self.fit_info["dof"] = dof
            self.fit_info["method"] = "lm"
            self.fit_info["x_scale"] = None
            return

        if method_str == "trf":
            print(f"Running trflsq fit with scipy least_squares x_scale={x_scale}")

            fit_lam = lam[mask]
            fit_flux = flux[mask]
            fit_w = w[mask]

            free_param_names, x0, bounds = self._free_parameter_state()

            if len(free_param_names) == 0:
                model_y = self.model(fit_lam)
                resid = (fit_flux - model_y) * fit_w
                chi2 = np.nansum(resid**2)
                dof = max(np.count_nonzero(mask), 1)

                self.message = "No free parameters."
                self.fit_info = {
                    "message": self.message,
                    "success": True,
                    "status": 0,
                    "method": "trf",
                    "x_scale": x_scale,
                    "nfev": 0,
                    "njev": 0,
                    "cost": 0.5 * chi2,
                    "optimality": np.nan,
                    "active_mask": np.asarray([], dtype=int),
                    "param_names": free_param_names,
                    "x": x0,
                    "chi2": chi2,
                    "redchi": chi2 / dof,
                    "dof": dof}
                return

            def residuals(values):
                self._set_free_parameters(free_param_names, values)
                return (fit_flux - self.model(fit_lam)) * fit_w

            result = least_squares(
                residuals,
                x0,
                bounds=bounds,
                method="trf",
                jac="2-point",
                x_scale=x_scale,
                max_nfev=maxiter,
                diff_step=np.sqrt(epsilon),
                xtol=acc)

            self._set_free_parameters(free_param_names, result.x)

            resid = (fit_flux - self.model(fit_lam)) * fit_w
            chi2 = np.nansum(resid**2)
            dof = max(np.count_nonzero(mask) - len(free_param_names), 1)
            redchi = chi2 / dof

            self.message = result.message
            self.fit_info = {
                "message": result.message,
                "success": result.success,
                "status": result.status,
                "method": "trf",
                "x_scale": x_scale,
                "nfev": result.nfev,
                "njev": result.njev,
                "cost": result.cost,
                "optimality": result.optimality,
                "active_mask": result.active_mask,
                "param_names": free_param_names,
                "x": result.x,
                "fun": result.fun,
                "jac": result.jac,
                "chi2": chi2,
                "redchi": redchi,
                "dof": dof}
            return

    def get_result(self, component_name):
        """Retrieve results from astropy model component.

        Every relevant value inside the astropy model, need to be
        written to the right position in the features table. For some
        cases (amplitude/power, fwhm/stddev), conversions are necessary
        (generally the inverse conversions of what was done in the
        register function).

        Parameters
        ----------
        component_name : str
            One of the names provided to any of the add_feature_*() calls
            made during setup.

        Returns
        -------
        dict with Parameters according to the PAHFIT definitions.

        e.g., for a feature with amplitude, stddev, and mean parameters:
        {'power': converted from amplitude, 'fwhm': converted from
        stddev, 'mean': wavelength}

        """
        if self.model is None:
            raise PAHFITModelError("Model not finalized yet.")

        if hasattr(self.model, "submodel_names"):
            component = self.model[component_name]
        else:
            # deals with edge case with single component, so is not
            # CompoundModel but normal single-component model.
            component = self.model

        c_type = self.feature_types[component_name]
        if c_type == "starlight" or c_type == "dust_continuum":
            return {
                "temperature": component.temperature.value,
                "tau": component.amplitude.value}
        elif c_type == "line":
            return {
                "power": component.power.value,
                "wavelength": component.mean.value,
                "fwhm": component.stddev.value * 2.355}
        elif c_type == "dust_feature":
            return {
                "power": component.power.value,
                "wavelength": component.x_0.value,
                "fwhm": component.fwhm.value}
        elif c_type == "attenuation":
            return {"tau": component.tau_sil.value}
        elif c_type == "absorption":
            return {
                "tau": component.tau.value,
                "wavelength": component.x_0.value,
                "fwhm": component.fwhm.value}
        else:
            raise PAHFITModelError(f"Unsupported component type: {c_type}")

    @staticmethod
    def _astropy_model_kwargs(component_name, param_names, param_values):
        """Create kwargs for an astropy model constructor.

        This is a utility that deduplicates the logic for going from
        (value, min, max) tuples, to astropy model constructor keyword
        arguments as in the following example:

        AstropyModelClass(name="feature name",
                    param1=value1,
                    param2=value2,
                    bounds={param1: (min,max), param2:(min,max)},
                    fixed={param1: True, param2: False})

        The returned arguments are in a dict that looks as follows, and
        can be passed to the appropriate astropy model constructor using
        **kwargs.

        {"name": "feature name"
         param_name: double, ...,
         "bounds": {param_name: array of size 2, ...},
         "fixed": {param_name: True or False, ...}}

        Parameters:
        -----------
        component_name : str
            Unique name for the component. Will later be used for
            indexing the components in the astropy model.

        param_names : list of str
            Names of the parameters for the astropy model, e.g.
            ["dust_feature1", "dust_feature2"]

        param_values : list of (array of size 3 OR scalar)
            One for each param name, each in the format of [value,
            min_bound, max_bound] for variable parameters, or a scalar
            (single float) for fixed parameters.

        Returns
        -------
        dict : kwargs to be used in an astropy model constructor

        """
        # basic format of constructor parameters of astropy model
        kwargs = {"name": component_name, "bounds": {}, "fixed": {}}

        for param_name, tuple_or_scalar in zip(param_names, param_values):
            if np.isscalar(tuple_or_scalar):
                is_fixed = True
                value, lo_bound, up_bound = tuple_or_scalar, 0, 0
            else:
                is_fixed = False
                value, lo_bound, up_bound = tuple_or_scalar

            # For the limits, use 0 if fixed, the raw values if
            # variable, but None if infinite (this is the convention for
            # astropy modeling)
            kwargs[param_name] = value
            kwargs["fixed"][param_name] = is_fixed
            kwargs["bounds"][param_name] = [
                None if np.isinf(x) else x for x in (lo_bound, up_bound)
            ]

        return kwargs
