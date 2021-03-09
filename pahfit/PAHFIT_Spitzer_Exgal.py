import numpy as np


__all__ = ["InstPackSpitzerIRSSLLL", "SciPackExGal"]


class InstPackSpitzerIRSSLLL:
    """
    Spitzer IRS SL and LL Instrument Pack

    Note: probably would be good to define a generic inst pack and have
    this class inherit it.  Would provide defaults.
    """

    def __init__(self):
        # what telescope/instrument/modes are covered by this pack
        self.telescope = "Spitzer"
        self.instrument = "IRS"
        self.modes = ["SL", "LL"]

    @staticmethod
    def fwhm(x):
        """
        Get the Spitzer fwhm starting values based on module origin

        Parameters
        ----------
        x: float
            central wavelengths of the line(s)
            """
        # in micron
        breaks = [0.0, 7.55, 14.6, 20.7, 1e6]
        # in micron
        widths = [0.053, 0.1, 0.14, 0.34]
        fwhms = np.zeros((len(x)))
        for k in range(len(breaks) - 1):
            indxs = np.where(np.logical_and(breaks[k] <= x, x < breaks[k + 1]))
            fwhms[indxs] = np.full((len(indxs)), widths[k])

        return fwhms


class SciPackExGal:
    """
    Science pack for spectra of whole or regions of galaxies

    Note: probably would be good to define a generic inst pack and have
    this class inherit it.  Would provide defaults.
    Also would be good to have some way of limiting the features
    used based on the wavelength range of the spectra to be fit.

    Parameters
    ----------
    inst_pack : InstPack class
        supplies
    fwhm_func : function
        function that takes wavelengths and returns the fwhms based
        on the instrument pack

    Attributes
    ----------
    bb_info : dict with {'amps', 'temps', 'amps_limits'}
    df_features,
    h2_features,
    ion_features: dict with {amps, x_0, fwhm,
                             amps_limits, x_0_limits, fwhms_limits}
    """

    def __init__(self, inst_pack=None):
        """
        Initialize the information for the components and subcomponents.

        Note: Currently given in python code.  Enhancement to have at
        least some (all?) of the values read in from a file needed.  Some
        of the values may be better derived from the file based information
        (e.g., fwhm of gas lines).
        """
        # the blackbodies for the continuum
        # 5000 K BB is for the stellar continuum,
        # rest for for the dust emission
        bb_temps = [5000.0, 300.0, 200.0, 135.0, 90.0, 65.0, 50.0, 40.0, 35.0]
        n_bb = len(bb_temps)
        bb_names = ["BB{}".format(k) for k in range(n_bb)]
        bb_amps = np.full(n_bb, 1.0e-10)
        self.bb_info = {
            "names": bb_names,
            "amps": bb_amps,
            "temps": bb_temps,
            "temps_limits": [(0.0, None) for k in range(n_bb)],
            "amps_limits": [(0.0, None) for k in range(n_bb)],
            "temps_fixed": [True for k in range(n_bb)],
            "amps_fixed": [False for k in range(n_bb)],
        }

        # the dust features (mainly the PAH/aromatic features)
        df_cwave = np.array(
            [
                5.27,
                5.70,
                6.22,
                6.69,
                7.42,
                7.60,
                7.85,
                8.33,
                8.61,
                10.68,
                11.23,
                11.33,
                11.99,
                12.62,
                12.69,
                13.48,
                14.04,
                14.19,
                15.9,
                16.45,
                17.04,
                17.375,
                17.87,
                18.92,
                33.1,
            ]
        )
        df_frac_fwhm = np.array(
            [
                0.034,
                0.035,
                0.030,
                0.07,
                0.126,
                0.044,
                0.053,
                0.05,
                0.039,
                0.02,
                0.012,
                0.032,
                0.045,
                0.042,
                0.013,
                0.04,
                0.016,
                0.025,
                0.02,
                0.014,
                0.065,
                0.012,
                0.016,
                0.019,
                0.05,
            ]
        )
        df_fwhm = df_frac_fwhm * df_cwave
        n_df = len(df_cwave)
        df_names = ["DF{}".format(k) for k in range(n_df)]
        df_amps = np.full((n_df), 100.0)
        df_cwave_limits = [(cwave - 0.1, cwave + 0.1) for cwave in df_cwave]
        df_fwhm_limits = [(cfwhm * 0.9, cfwhm * 1.1) for cfwhm in df_fwhm]
        self.dust_features = {
            "names": df_names,
            "amps": df_amps,
            "x_0": df_cwave,
            "fwhms": df_fwhm,
            "amps_limits": [(0.0, None) for k in range(n_df)],
            "x_0_limits": df_cwave_limits,
            "fwhms_limits": df_fwhm_limits,
            "amps_fixed": [False for k in range(n_df)],
            "x_0_fixed": [True for k in range(n_df)],
            "fwhms_fixed": [True for k in range(n_df)],
        }

        # define the H2 features
        h2_cwaves = np.array(
            [5.5115, 6.1088, 6.9091, 8.0258, 9.6649, 12.2785, 17.0346, 28.2207]
        )
        if inst_pack is None:
            h2_fwhm = np.full((len(h2_cwaves)), 0.1)
        else:
            h2_fwhm = inst_pack.fwhm(h2_cwaves)
        h2_names = np.array(
            [
                "H2 S(7)",
                "H2 S(6)",
                "H2 S(5)",
                "H2 S(4)",
                "H2 S(3)",
                "H2 S(2)",
                "H2 S(1)",
                "H2 S(0)",
            ]
        )
        n_h2 = len(h2_cwaves)
        h2_amps = np.full((n_h2), 100.0)
        h2_cwaves_limits = [(cwave - 0.05, cwave + 0.05) for cwave in h2_cwaves]
        h2_fwhm_limits = [(cfwhm * 0.5, cfwhm * 1.5) for cfwhm in h2_fwhm]
        self.h2_features = {
            "names": h2_names,
            "amps": h2_amps,
            "x_0": h2_cwaves,
            "fwhms": h2_fwhm,
            "amps_limits": [(0.0, None) for k in range(n_h2)],
            "x_0_limits": h2_cwaves_limits,
            "fwhms_limits": h2_fwhm_limits,
            "amps_fixed": [False for k in range(n_h2)],
            "x_0_fixed": [False for k in range(n_h2)],
            "fwhms_fixed": [False for k in range(n_h2)],
        }

        # define the ionic features
        ion_cwaves = np.array(
            [
                6.985274,
                8.99138,
                10.5105,
                12.813,
                15.555,
                18.713,
                25.91,
                25.989,
                33.480,
                34.8152,
            ]
        )
        if inst_pack is None:
            ion_fwhm = np.full((len(ion_cwaves)), 0.1)
        else:
            ion_fwhm = inst_pack.fwhm(ion_cwaves)
        ion_names = np.array(
            [
                "[ArII]",
                "[ArIII]",
                "[SIV]",
                "[NeII]",
                "[NeIII]",
                "[SIII] 18",
                "[OIV]",
                "[FeII]",
                "[SIII] 33",
                "[SiII]",
            ]
        )
        n_ion = len(ion_cwaves)
        ion_amps = np.full((n_ion), 100.0)
        ion_cwaves_limits = [(cwave - 0.05, cwave + 0.05) for cwave in ion_cwaves]
        ion_fwhm_limits = [(cfwhm * 0.5, cfwhm * 1.5) for cfwhm in ion_fwhm]
        self.ion_features = {
            "names": ion_names,
            "amps": ion_amps,
            "x_0": ion_cwaves,
            "fwhms": ion_fwhm,
            "amps_limits": [(0.0, None) for k in range(n_ion)],
            "x_0_limits": ion_cwaves_limits,
            "fwhms_limits": ion_fwhm_limits,
            "amps_fixed": [False for k in range(n_ion)],
            "x_0_fixed": [False for k in range(n_ion)],
            "fwhms_fixed": [False for k in range(n_ion)],
        }

        self.att_info = {
            "names": ["S07_att"],
            "amps": [1.0],
            "amps_limits": [(0.0, 10.0)],
            "amps_fixed": [False],
        }
