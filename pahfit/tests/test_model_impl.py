from pahfit.helpers import read_spectrum
from specutils import Spectrum1D


def test_feature_table_model_conversion():
    # do a fit
    spectrumfile = "M101_Nucleus_irs.ipac"
    spec = read_spectrum(spectrumfile, spec1d=True)
    packfile = "classic.yaml"
    # use a spitzer instrument model that covers the required range. SL1, SL2, LL1, LL2 should do
    instrumentname = f"spitzer.irs.*.[12]"

    model = Model.from_yaml(packfile, instrumentname, 0)
    model.guess(spec)
    model.fit(spec)

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


def test_feature_table_manual_edit():
    # test outline:
    # make model
    # construct astropy model
    # copy model
    # edit features table
    # construct astropy model 2
    # change should be reflected in astropy model 2, but not astropy model 1
    pass
