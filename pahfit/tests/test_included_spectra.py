from pahfit.helpers import read_spectrum


def test_read_spectrum():
    """Tests to make sure that ALL the included files can be read in
    out of the box. Later, we could add a test that makes sure that all
    included data at least don't make PAHFIT crash (e.g. issue 239
    https://github.com/PAHFIT/pahfit/issues/239)
    """
    files = [
        "1C_snr_gt20.ecsv",
        "ISO_SWS/Orion_Brgamma_ISO-SWS_merged.ipac",
        "Lai2020_1C_akari_spitzer.ecsv",
        "M101_Nucleus_irs.ipac",
        "MIRI_MRS/VV114_MIRI_ch1_LONG_s1d.fits",
        "orion_bar_SWS_with7023shape.ipac",
    ]
    for f in files:
        _ = read_spectrum(f)
