#!/usr/bin/env python

import os
import pkg_resources
import argparse

import matplotlib.pyplot as plt
import matplotlib as mpl

import astropy.units as u
from astropy.table import Table

from pahfit.base import PAHFITBase


def initialize_parser():
    """
    Command line parser for plot_pahfit

    Returns
    -------
    parser : argparse object
    """
    plottypes = [
        "png",
        "jpg",
        "jpeg",
        "pdf",
        "ps",
        "eps",
        "rgba",
        "svg",
        "tiff",
        "tif",
        "pgf",
        "svgz",
        "raw",
    ]
    parser = argparse.ArgumentParser()
    parser.add_argument("spectrumfile", help="name of file with observed spectrum")
    parser.add_argument("fitfilename", help="name of PAHFIT results file")
    parser.add_argument(
        "--savefig",
        action="store",
        default=None,
        choices=plottypes,
        help="Save figure to a file of specified type",
    )
    parser.add_argument(
        "--scalefac_resid",
        action="store",
        type=int,
        default=2,
        help="Factor multiplying the standard deviation of the residuals to adjust plot limits",
    )

    return parser


def main():

    # commandline parser
    parser = initialize_parser()
    args = parser.parse_args()

    # read an observed spectrum
    # read in the observed spectrum
    # assumed to be astropy table compatibile and include units
    specfile = args.spectrumfile
    outputname = specfile.split(".")[0]
    if not os.path.isfile(specfile):
        pack_path = pkg_resources.resource_filename("pahfit", "data/")
        test_specfile = "{}/{}".format(pack_path, specfile)
        if os.path.isfile(test_specfile):
            specfile = test_specfile
        else:
            raise ValueError("Input spectrumfile {} not found".format(specfile))

    # get the table format (from extension of filename)
    tformat = specfile.split(".")[-1]
    if tformat == "ecsv":
        tformat = "ascii.ecsv"
    obs_spectrum = Table.read(specfile, format=tformat)
    obs_x = obs_spectrum["wavelength"].to(u.micron, equivalencies=u.spectral())
    obs_y = obs_spectrum["flux"].to(u.Jy, equivalencies=u.spectral_density(obs_x))
    obs_unc = obs_spectrum["sigma"].to(u.Jy, equivalencies=u.spectral_density(obs_x))

    # strip units as the observed spectrum is in the internal units
    obs_x = obs_x.value
    obs_y = obs_y.value

    # read in the PAHFIT results
    pmodel = PAHFITBase(obs_x, obs_y, filename=args.fitfilename)

    # plot result
    fontsize = 18
    font = {"size": fontsize}
    mpl.rc("font", **font)
    mpl.rc("lines", linewidth=2)
    mpl.rc("axes", linewidth=2)
    mpl.rc("xtick.major", size=5, width=1)
    mpl.rc("ytick.major", size=5, width=1)
    mpl.rc("xtick.minor", size=3, width=1)
    mpl.rc("ytick.minor", size=3, width=1)

    fig, axs = plt.subplots(ncols=1, nrows=2, figsize=(15, 10),
                            gridspec_kw={'height_ratios': [3, 1]},
                            sharex=True)

    pmodel.plot(axs, obs_x, obs_y, obs_unc.value, pmodel.model, scalefac_resid=args.scalefac_resid)

    # use the whitespace better
    fig.subplots_adjust(hspace=0)

    # show or save
    if args.savefig:
        fig.savefig("{}.{}".format(outputname, args.savefig))
    else:
        plt.show()


if __name__ == "__main__":
    main()
