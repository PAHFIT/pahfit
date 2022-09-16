import numpy as np

from pahfit.PAHFIT_Spitzer_Exgal import InstPackSpitzerIRSSLLL, SciPackExGal
from pahfit.model import Model


def test_classic_pack():
    """Compares param_info between one that was constructed from the classic
    input file, and one that has been hardcoded."""
    # first use the static code to generate the feature dictonaries
    # instrument pack
    spitzer_sl_ll_pack = InstPackSpitzerIRSSLLL()
    # science pack
    sci_pack = SciPackExGal(spitzer_sl_ll_pack)
    oparam_info = (
        sci_pack.bb_info,
        sci_pack.dust_features,
        sci_pack.h2_features,
        sci_pack.ion_features,
        sci_pack.att_info,
    )

    # now read in the equivalent info from a file

    # setup the model
    packfile = "classic.yaml"
    # NOTE: Since the new instrument model modifies line widths, this
    # test will fail.
    instrumentname = "spitzer.irs.sl.2"

    model = Model.from_yaml(packfile, instrumentname, 0)
    nparam_info = model._kludge_param_info()

    # check the different dictonaries are equivalent
    for k in range(len(oparam_info)):
        for ckey in oparam_info[k].keys():
            if isinstance(oparam_info[k][ckey][0], tuple):
                # loop over each tuple and test
                # complicated as tuple and has None values
                assert len(oparam_info[k][ckey]) == len(nparam_info[k][ckey])
                for i in range(len(oparam_info[k][ckey])):
                    o1, o2 = oparam_info[k][ckey][i]
                    n1, n2 = nparam_info[k][ckey][i]
                    if (o1 is not None) and (n1 is not None):
                        np.testing.assert_allclose(o1, n1)
                    else:
                        assert o1 == n1
                    if (o2 is not None) and (n2 is not None):
                        np.testing.assert_allclose(o2, n2)
                    else:
                        assert o2 == n2
            elif isinstance(oparam_info[k][ckey][0], float):
                np.testing.assert_allclose(oparam_info[k][ckey], nparam_info[k][ckey])
            else:
                np.testing.assert_equal(oparam_info[k][ckey], nparam_info[k][ckey])
