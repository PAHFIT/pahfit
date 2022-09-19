from specutils import Spectrum1D
from astropy import units as u
import copy
from astropy.modeling.fitting import LevMarLSQFitter
from matplotlib import pyplot as plt

from pahfit.helpers import find_packfile
from pahfit.features import Features
from pahfit.base import PAHFITBase
from pahfit import instrument


class Model:
    """This class acts as the main API for PAHFIT.

    The users deal with model objects, of which the state is modified
    during initalization, initial guessing, and fitting. What the model
    STORES is a description of the physics: what features are there and
    what are their properties, regardless of the instrument with which
    those features are observed. The methods provided by this class,
    form the connection between those physics, and what is observed.
    During fitting and plotting, those physics are converted into a
    model for the observation, by applying instrumental parameters from
    the instrument.py module.

    The main thing that defines a model, is the features table, loaded
    from a YAML file given to the constructor. After construction, the
    Model can be edited by accessing the stored features table directly.
    Changing numbers in this table, is allowed, and the updated numbers
    will be reflected when the next fit or initial guess happens. At the
    end of these actions, the fit or guess results are stored in the
    same table.

    The model can be saved.

    The model can be copied.

    Attributes
    ----------
    features : Features
        Instance of the Features class. Can be edited on-the-fly.
        Non-breaking behavior by the user is expected. Changes will be
        reflected at the next fit, guess, or plot call.

    """

    def __init__(self, features: Features, instrumentname, redshift):
        """
        Parameters
        ----------
        features: Features
            Features table.

        instrumentname : str or list of str
            Qualified instrument name, see instrument.py. This will
            determine what the line widths are, when going from the
            features table to a fittable/plottable model.

        redshift : float
            Redshift used to shift from the physical model, to the
            observed model.

        """
        self.redshift = redshift
        self.instrumentname = instrumentname
        self.features = features

    @classmethod
    def from_yaml(cls, pack_file, instrumentname, redshift):
        """
        Generate feature table from YAML file.

        Parameters
        ----------
        pack_file : str
            Path to YAML file, or name of one of the default YAML files.

        Returns
        -------
        Model instance

        """
        path = find_packfile(pack_file)
        features = Features.read(path)
        return cls(features, instrumentname, redshift)

    @classmethod
    def from_saved(cls, saved_model_file):
        """
        Parameters
        ----------
        saved_model_file : str
           Path to file generated by Model.save()

        Returns
        -------
        Model instance
        """
        # features.read automatically switches to astropy table reader.
        # Maybe needs to be more advanced here in the future.
        features = Features.read(saved_model_file, format="ascii.ecsv")
        return cls(features, features.meta["instrumentname"], features.meta["redshift"])

    def save(self, fn, **write_kwargs):
        """Save the model to disk.

        Only ECSV supported for now. Models saved this way can be read
        back in, with metadata.

        TODO: store details about the fit results somehow. Uncertainties
        (covariance matrix) should be retrievable. Use Table metadata?

        Parameters
        ----------
        fn : file name

        **write_kwargs : kwargs passed to astropy.table.Table.write

        """
        if fn.split(".")[-1] != "ecsv":
            raise NotImplementedError("Only ecsv is supported for now")

        self.features.meta.update(
            {"redshift": self.redshift, "instrumentname": self.instrumentname}
        )
        self.features.write(fn, format="ascii.ecsv", **write_kwargs)

    def guess(self, spec: Spectrum1D):
        """Make an initial guess of the physics, based on the given
        observational data.

        Parameters
        ----------
        spec : Spectrum1D
            1D (not 2D or 3D) spectrum object, containing the
            observational data. (TODO: should support list of spectra,
            for the segment-based joint fitting). Initial guess will be
            based on the flux in this spectrum.

        Returns
        -------
        Nothing, but internal feature table is updated.

        """
        x = spec.spectral_axis.to(u.micron).value
        y = spec.flux.value  # TODO figure out right unit
        xz = x / (1 + self.redshift)

        # remake param_info to make sure we have any feature updates from the user
        param_info = self._kludge_param_info()
        param_info = PAHFITBase.estimate_init(xz, y, param_info)
        self._backport_param_info(param_info)

    def fit(self, spec: Spectrum1D, maxiter=1000, verbose=True):
        """Fit the observed data.

        The model setup is based on the features table and instrument specification.

        The last fit results can accessed through the variable
        model.astropy_result. The results are also stored back to the
        model.features table.

        CAVEAT: any features that do not overlap with the data range
        will not be included in the model, for performance and numerical
        stability. Their values in the features table will be left
        untouched.

        Parameters
        ----------
        spec : Spectrum1D
            Observed spectrum containing wavelengths, flux, and
            uncertainties. Needs to be compatible with the instrument
            specification of this Model object.
            TODO: convert flux to preferred units

        maxiter : int
            maximum number of fitting iterations

        verbose : boolean
            set to provide screen output

        """
        x = spec.spectral_axis.to(u.micron).value
        y = spec.flux.value
        w = 1.0 / spec.uncertainty.array

        # check if spectrum is compatible with instrument model
        instrument.check_range([min(x), max(x)], self.instrumentname)

        # transform observed wavelength to "physical" wavelength
        xz = x / (1 + self.redshift)

        # construct model
        astropy_model = self._construct_astropy_model()

        # pick the fitter
        fit = LevMarLSQFitter()

        # fit
        self.astropy_result = fit(
            astropy_model,
            xz,
            y,
            weights=w,
            maxiter=maxiter,
            epsilon=1e-10,
            acc=1e-10,
        )
        if verbose:
            print(fit.fit_info["message"])

        self._parse_astropy_result(self.astropy_result)

    def info(self):
        """Print out the last fit results."""
        print(self.astropy_result)

    def plot(self, spec=None):
        """Plot model, and optionally compare to observational data.

        Parameters
        ----------
        spec : Spectrum1D
            Observational data. Does not have to be the same data that
            was used for guessing or fitting.
        """
        # copied some stuff from plot_pahfit
        fig, axs = plt.subplots(
            ncols=1,
            nrows=2,
            figsize=(15, 10),
            gridspec_kw={"height_ratios": [3, 1]},
            sharex=True,
        )

        x = spec.spectral_axis.to(u.micron).value
        xz = x / (1 + self.redshift)
        y = spec.flux.value
        unc = spec.uncertainty.array
        astropy_model = self._construct_astropy_model()
        PAHFITBase.plot(axs, xz, y, unc, astropy_model)

        fig.subplots_adjust(hspace=0)

    def copy(self):
        """Copy the model.

        Main use case: use this model as a parent model for more
        fits.

        Currently uses copy.deepcopy. We should do something smarter if
        we run into memory problems or sluggishness.

        Returns
        -------
        model_copy : Model
        """
        # We could do this
        # make new instance
        # model_copy = type(self)(self.pack_file, self.instrumentname, self.redshift)
        # copy over all the variables that might have changed
        # make sure to deep copy the table!

        # But maybe try this first
        return copy.deepcopy(self)

    def _kludge_param_info(self):
        param_info = PAHFITBase.parse_table(self.features)
        # edit line widths and drop lines out of range
        param_info[2] = PAHFITBase.update_dictionary(
            param_info[2],
            self.instrumentname,
            update_fwhms=True,
            redshift=self.redshift,
        )
        param_info[3] = PAHFITBase.update_dictionary(
            param_info[3],
            self.instrumentname,
            update_fwhms=True,
            redshift=self.redshift,
        )
        param_info[4] = PAHFITBase.update_dictionary(
            param_info[4], self.instrumentname, redshift=self.redshift
        )

        return param_info

    def _backport_param_info(self, param_info):
        """Convert param_info to values in features table.

        Temporary hack to make the new system compatible with the old system.

        TODO: if we remove the param_info stuff entirely, we won't need this

        """
        # unfortunately, there is no implementation for this, even in
        # the original code. That one goes straight from astropy model
        # to table... But we can do a kludge here: convert to model
        # first, and then back to table.
        astropy_model = PAHFITBase.model_from_param_info(param_info)
        self._parse_astropy_result(astropy_model)

    def _construct_astropy_model(self):
        """Convert the features table into a fittable model.

        TODO: Make sure the features outside of the data range are
        removed. The instrument-based feature check is done in
        _kludge_param_info(), but the observational data might only
        cover a part of the instrument range."""
        param_info = self._kludge_param_info()
        return PAHFITBase.model_from_param_info(param_info)

    def _parse_astropy_result(self, astropy_model):
        """Store the result of the astropy fit into the features table.

        Every relevant value inside the astropy model, is written to the
        right position in the features table. This way, the astropy
        model and the features table are kept in sync.

        Doing things this way, makes it possible for the user to make
        edits to the features table, and makes it easy to store the
        model (just store the features table)

        """
        # Some translation rules between astropy model components and
        # feature table names and values.

        # these have the same value but different (fixed) names
        param_name_equivalent = {
            "temperature": "temperature",
            "fwhm": "fwhm",
            "x_0": "wavelength",
            "mean": "wavelength",
            "tau_sil": "tau",
        }

        def param_conversion(features_kind, param_name, param_value):
            # default conversion
            if param_name in param_name_equivalent:
                new_name = param_name_equivalent[param_name]
                new_value = param_value
            # certain types of features use tau instead of amplitude
            elif param_name == "amplitude":
                if features_kind in ["starlight", "dust_continuum", "absorption"]:
                    new_name = "tau"
                else:
                    new_name = "power"
                new_value = param_value
            # convert stddev to fwhm
            elif param_name == "stddev":
                new_name = "fwhm"
                new_value = param_value * 2.355
            else:
                raise NotImplementedError(
                    f"no conversion rule for model parameter {param_name}"
                )
            return new_name, new_value

        # now apply these rules to fill in the parameters in the right
        # spots of the table
        for component in astropy_model:
            # find the matching row
            if component.name == "S07_att":
                # Just need to watch out for S07_att. It is named silicate
                # in the features table.
                feature_name = "silicate"
            else:
                feature_name = component.name
            row = self.features.loc[feature_name]

            # write the new values (every element is masked array [value
            # lowerbound upperbound])
            for param_name in component.param_names:
                param_value = getattr(component, param_name).value
                col_name, col_value = param_conversion(
                    row["kind"], param_name, param_value
                )
                # even though values specified as fixed in the table
                # shouldn't change, there is one exception to this:
                # lines in overlapping segments will be forced to have
                # variable widths. Therefore, do an explicit check here
                # to make sure that we're not writing to fixed values
                # TODO: how do we store the variable line widths then?
                table_element = row[col_name]
                fixed = table_element.mask[1]
                if not fixed:
                    table_element = col_value
