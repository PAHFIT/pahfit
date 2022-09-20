from pahfit.model import Model
from pahfit.helpers import find_packfile
from pahfit.features import Features


def test_feature_parsing():
    """
    Goal
    ----
    Test if the model is built correctly with certain features removed
    from the science pack. For example, when this test was written,
    generating the model without specifying any gaussians will cause the
    code to crash. This test will try to provoke such crashes.

    Desired behavior
    ----------------

    The PAHFITBase instance is generated correctly, without crashing.

    Functions that depend on specific model contents (lines, dust
    features, ...) can deal with those feature not being there.

    """
    # random instrument name
    instrumentname = "spitzer.irs.sl.2"

    # choose any science pack
    packfile = find_packfile("classic.yaml")

    # convert the yaml prescription to a features table
    features = Features.read(packfile)

    def test_parsing(features_edit):
        m = Model(features_edit, instrumentname, 0)
        amodel = m._construct_astropy_model()
        m._parse_astropy_result(amodel)

    # Case 0: the whole table
    test_parsing(features)

    # Cases 1, 2, ... (can it run without any of the following?)
    remove_kinds = [
        "dust_continuum",
        "dust_feature",
        "line",
        "starlight",
        "attenuation",
    ]
    for kind in remove_kinds:
        features_edit = features[features["kind"] != kind]
        test_parsing(features_edit)


if __name__ == "__main__":
    test_feature_parsing()
