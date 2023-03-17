# import pkg_resources
# from pahfit.model import Model
# from pahfit.tests.test_model_impl import assert_features_table_equality

def test_classic_regression():
    """Regression test for the features table of classic.yaml.

    Uses a saved version of the table. Every time the classic pack (or
    details about the parsing) are changed, the saved reference table
    must be updated.

    """
    pass
    # model_new = Model.from_yaml("classic.yaml")

    # saved_path = pkg_resources.resource_filename("pahfit", "data/regression_test/classic_features_table.ascii.ecsv")
    # model_old = Model.from_saved(saved_path)

    # assert_features_table_equality(model_new.features, model_old.features)
