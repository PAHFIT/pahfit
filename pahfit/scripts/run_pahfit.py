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
        "--estimate_start",
        action="store_true",
        help="Estimate of starting point based on the input spectrum",
    )

    return parser


def main():

    # setup and parse the command line
    parser = initialize_parser()
    args = parser.parse_args()

    # read in the spectrum
    obsdata = read_spectrum(args.spectrumfile)

    # setup the model
    pmodel = initialize_model(args.packfile, obsdata, args.estimate_start)

    # fit the spectrum
    obsfit = fit_spectrum(obsdata, pmodel)

    # save fit results to file
    outputname = args.spectrumfile.split(".")[0]
    pmodel.save(obsfit, outputname, args.saveoutput)

    # plot result
    fontsize = 18
    font = {"size": fontsize}
    mpl.rc("font", **font)
    mpl.rc("lines", linewidth=2)
    mpl.rc("axes", linewidth=2)
    mpl.rc("xtick.major", width=2)
    mpl.rc("ytick.major", width=2)

    fig, ax = plt.subplots(figsize=(15, 10))

    pmodel.plot(ax, obsdata["x"], obsdata["y"], obsfit)

    ax.set_yscale("linear")
    ax.set_xscale("log")

    # use the whitespace better
    fig.tight_layout()

    # show
    if args.showplot:
        plt.show()
    # save (always)
    fig.savefig("{}.{}".format(outputname, args.savefig))


if __name__ == "__main__":
    main()
