

def test_fitting_m101():

    # get an observed extinction curve to fit
    g03_model = G03_LMCAvg()

    x = g03_model.obsdata_x
    # convert to E(x-V)/E(B0V)
    y = (g03_model.obsdata_axav - 1.0) * g03_model.Rv
    # only fit the UV portion (FM90 only valid in UV)
    gindxs, = np.where(x > 3.125)

    fm90_init = FM90()
    fit = LevMarLSQFitter()
    g03_fit = fit(fm90_init, x[gindxs], y[gindxs])
    fit_vals = [
        g03_fit.C1.value,
        g03_fit.C2.value,
        g03_fit.C3.value,
        g03_fit.C4.value,
        g03_fit.xo.value,
        g03_fit.gamma.value,
    ]

    good_vals = np.array(
        [
            -0.958016797002,
            1.0109751831,
            2.96430606652,
            0.313137860902,
            4.59996300532,
            0.99000982258,
        ]
    )

    np.testing.assert_allclose(good_vals, fit_vals)
