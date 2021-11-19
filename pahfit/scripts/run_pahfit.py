#!/usr/bin/env python

import argparse

import matplotlib.pyplot as plt
import matplotlib as mpl

from pahfit.helpers import read_spectrum, initialize_model, fit_spectrum


def initialize_parser():
    """
    Command line parser for run_pahfit

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
    savetypes = ["fits", "votable", "ipac", "ascii.ecsv"]
    parser = argparse.ArgumentParser()
    parser.add_argument("spectrumfile", help="name of file with observed spectrum")
    parser.add_argument("packfile", help="name of PAHFIT pack file")
    parser.add_argument(
        "--savefig",
        action="store",
        default="pdf",
        choices=plottypes,
        help="Save figure to a file of specified type",
    )
    parser.add_argument(
        "--showplot", action="store_true", help="display plot to the screen"
    )
    parser.add_argument(
        "--saveoutput",
        action="store",
        default="ipac",
        choices=savetypes,
        help="Save fit results to a file of specified type",
    )
    parser.add_argument(
        "--no_starting_estimate",
        action="store_true",
        help="Bypass the estimation of the fit starting point based on the input spectrum.",
    )
    parser.add_argument(
        "--scalefac_resid",
        action="store",
        type=float,
        default=2.0,
        help="Factor multiplying the standard deviation of the residuals to adjust plot limits",
    )
    parser.add_argument(
        "--fit_maxiter",
        default=1000,
        type=int,
        help="Maximum number of interations for the fitting",
    )

    return parser


def main():

    # setup and parse the command line
    parser = initialize_parser()
    args = parser.parse_args()

    # read in the spectrum
    obsdata = read_spectrum(args.spectrumfile)

    # setup the model
    pmodel = initialize_model(args.packfile, obsdata, not args.no_starting_estimate)

    # fit the spectrum
    obsfit = fit_spectrum(obsdata, pmodel, maxiter=args.fit_maxiter)

    # save fit results to file
    outputname = args.spectrumfile.split(".")[0]
    pmodel.save(obsfit, outputname, args.saveoutput)

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

    pmodel.plot(axs, obsdata["x"], obsdata["y"], obsdata["unc"], obsfit, scalefac_resid=args.scalefac_resid)

    # use the whitespace better
    fig.subplots_adjust(hspace=0)

    # show
    if args.showplot:
        plt.show()
    # save (always)
    fig.savefig("{}.{}".format(outputname, args.savefig))


if __name__ == "__main__":
    main()
