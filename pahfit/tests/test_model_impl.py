from pahfit.helpers import read_spectrum
from pahfit.model import Model
import tempfile
import numpy as np
import os
from astropy import units as u


def assert_features_table_equality(features1, features2):
    for string_col in ["name", "group", "kind", "model", "geometry"]:
        assert (features1[string_col] == features2[string_col]).all()
    for param_col in ["temperature", "tau", "wavelength", "power", "fwhm"]:
        for k in ("val", "min", "max"):
            np.testing.assert_allclose(
                features1[param_col][k], features2[param_col][k], rtol=1e-6, atol=1e-6
            )


def default_spec_and_model_fit(fit=True):
    spectrumfile = "M101_Nucleus_irs.ipac"
    spec = read_spectrum(spectrumfile)
    packfile = "classic.yaml"
    # use a spitzer instrument model that covers the required range. SL1, SL2, LL1, LL2 should do
    spec.meta["instrument"] = "spitzer.irs.*.[12]"
    model = Model.from_yaml(packfile)
    if fit:
        model.guess(spec)
        model.fit(spec)
    return spec, model


def test_feature_table_model_conversion():
    spec, model = default_spec_and_model_fit()

    # two versions of the fit result: one is fresh from the fit. The fit
    # results were then written to model.features. If everything went
    # correct, reconstructing the model from model.features should
    # result in the exact same model.

    # test only works for the astropy-based implementation at the moment.
    fit_result = model.fitter
    reconstructed_fit_result = model._set_up_fitter(
        instrumentname=spec.meta["instrument"], redshift=0, use_instrument_fwhm=False
    )
    for name in fit_result.components():
        par_dict1 = fit_result.get_result(name)
        par_dict2 = reconstructed_fit_result.get_result(name)
        for key in par_dict1:
            assert par_dict1[key] == par_dict2[key]


def test_model_edit():
    # make model from default feature list and a copy
    _, model = default_spec_and_model_fit(fit=False)
    model_to_edit = model.copy()

    # remember this random parameter
    feature = "dust_cont00"
    col = "temperature"
    originalT = model.features.loc[feature][col][0]

    # edit the same parameter in the copy
    newT = 123

    i = np.where(model_to_edit.features["name"] == feature)[0]
    model_to_edit.features[col]["val"][i] = newT

    # make sure the original value is still the same
    j = np.where(model.features["name"] == feature)[0]
    assert model.features[col]["val"][j] == originalT

    # construct astropy model with dummy instrument
    fitter_edit = model_to_edit._set_up_fitter(
        instrumentname="spitzer.irs.*", redshift=0
    )

    # Make sure the change is reflected in this astropy model. Very
    # handy that we can access the right component by the feature name!
    assert fitter_edit.get_result(feature)["temperature"] == newT


def test_model_tabulate():
    """This function has several options. The following tests runs will
    make these more maintainable."""
    # do not fit yet
    spec, model = default_spec_and_model_fit(False)

    # test tabulate before fitting, unit needs to be unitless
    df = model.tabulate(
        instrumentname="spitzer.irs.*.[12]",
        feature_mask=model.features["kind"] == "dust_feature",
    )
    assert df.unit == u.dimensionless_unscaled

    # with spec given (should still be dimensionless)
    df = model.tabulate(
        wavelengths=spec,
        instrumentname="spitzer.irs.*.[12]",
        feature_mask=model.features["kind"] == "dust_feature",
    )
    assert df.unit == u.dimensionless_unscaled

    # after fitting
    model.fit(spec)

    # default wavelength grid. Unit should be the same as spec.
    df = model.tabulate(
        instrumentname="spitzer.irs.*.[12]",
        feature_mask=model.features["kind"] == "dust_feature",
    )
    assert df.unit == spec.unit

    # spec wavelength grid. Length should be the same as spec.
    tab_Jy = model.tabulate(
        wavelengths=spec,
        instrumentname="spitzer.irs.*.[12]",
    )
    assert tab_Jy.shape == spec.shape


def test_save_load():
    _, model = default_spec_and_model_fit()

    with tempfile.TemporaryDirectory() as d:
        fn = os.path.join(d, "save_test.ecsv")
        model.save(fn, overwrite=True)
        model_loaded = Model.from_saved(fn)

    assert_features_table_equality(model.features, model_loaded.features)
