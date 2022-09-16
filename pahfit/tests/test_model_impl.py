from pahfit.helpers import read_spectrum
from pahfit.model import Model
import tempfile
import numpy as np


def assert_features_table_equality(features1, features2):
    for string_col in ["name", "group", "kind", "model", "geometry"]:
        assert (features1[string_col] == features2[string_col]).all()
    for param_col in ["temperature", "tau", "wavelength", "power", "fwhm"]:
        np.testing.assert_allclose(
            features1[param_col], features2[param_col], rtol=1e-6, atol=1e-6
        )


def default_spec_and_model_fit(fit=True):
    spectrumfile = "M101_Nucleus_irs.ipac"
    spec = read_spectrum(spectrumfile, spec1d=True)
    packfile = "classic.yaml"
    # use a spitzer instrument model that covers the required range. SL1, SL2, LL1, LL2 should do
    instrumentname = "spitzer.irs.*.[12]"
    model = Model.from_yaml(packfile, instrumentname, 0)
    if fit:
        model.guess(spec)
        model.fit(spec)
    return spec, model


def test_feature_table_model_conversion():
    _, model = default_spec_and_model_fit()

    # two versions of the fit result: one is fresh from the fit. The fit
    # results were then written to model.features. If everything went
    # correct, reconstructing the model from model.features should
    # result in the exact same model.
    fit_result = model.astropy_result
    reconstructed_fit_result = model._construct_astropy_model()
    for p in fit_result.param_names:
        p1 = getattr(fit_result, p)
        p2 = getattr(reconstructed_fit_result, p)
        assert p1 == p2


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
    model_to_edit.features.loc[feature][col][0] = newT

    # make sure the original value is still the same
    assert model.features.loc[feature][col][0] == originalT

    # construct astropy model
    astropy_model_edit = model_to_edit._construct_astropy_model()

    # Make sure the change is reflected in this model. Very handy that
    # we can access the right component by the feature name!
    assert astropy_model_edit[feature].temperature == newT


def test_save_load():
    _, model = default_spec_and_model_fit()
    temp_file = tempfile.NamedTemporaryFile(suffix=".ecsv")
    model.save(temp_file.name, overwrite=True)
    model_loaded = Model.from_saved(temp_file.name)

    # all the things we want to recover
    assert model.redshift == model_loaded.redshift
    assert model.instrumentname == model_loaded.instrumentname
    assert_features_table_equality(model.features, model_loaded.features)
