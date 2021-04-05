#!/usr/bin/env python

import argparse

import matplotlib.pyplot as plt
import matplotlib as mpl

from pahfit.helpers import read_spectrum, initialize_model


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
    return parser


def main():

    # commandline parser
    parser = initialize_parser()
    args = parser.parse_args()

    # read in the spectrum
    obsdata = read_spectrum(args.spectrumfile)

    # setup the model
    pmodel = initialize_model(args.fitfilename, obsdata, estimate_start=False)

    # plot result
    fontsize = 18
    font = {"size": fontsize}
    mpl.rc("font", **font)
    mpl.rc("lines", linewidth=2)
    mpl.rc("axes", linewidth=2)
    mpl.rc("xtick.major", width=2)
    mpl.rc("ytick.major", width=2)

    fig, ax = plt.subplots(figsize=(15, 10))

    pmodel.plot(ax, obsdata["x"].value, obsdata["y"].value, pmodel.model)

    ax.set_yscale("linear")
    ax.set_xscale("log")

    # use the whitespace better
    fig.tight_layout()

    # show or save
    outputname = args.spectrumfile.split(".")[0]
    if args.savefig:
        fig.savefig("{}.{}".format(outputname, args.savefig))
    else:
        plt.show()


if __name__ == "__main__":
    main()
