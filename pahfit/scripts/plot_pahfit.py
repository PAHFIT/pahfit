#!/usr/bin/env python

import argparse

import matplotlib.pyplot as plt
import matplotlib as mpl

from pahfit.model import Model
from pahfit.helpers import read_spectrum


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
        type=float,
        default=2.0,
        help="Factor multiplying the standard deviation of the residuals to adjust plot limits",
    )

    return parser


def main():
    # commandline parser
    parser = initialize_parser()
    args = parser.parse_args()

    # read in the spectrum
    spec = read_spectrum(args.spectrumfile)

    # setup the model from saved one
    model = Model.from_saved(args.fitfilename)

    fig = default_layout_plot(spec, model, args.scalefac_resid)

    # show or save
    outputname = args.spectrumfile.split(".")[0]
    if args.savefig:
        fig.savefig("{}.{}".format(outputname, args.savefig))
    else:
        plt.show()


def default_layout_plot(spec, model, scalefac_resid):
    fontsize = 18
    font = {"size": fontsize}
    mpl.rc("font", **font)
    mpl.rc("lines", linewidth=2)
    mpl.rc("axes", linewidth=2)
    mpl.rc("xtick.major", size=5, width=1)
    mpl.rc("ytick.major", size=5, width=1)
    mpl.rc("xtick.minor", size=3, width=1)
    mpl.rc("ytick.minor", size=3, width=1)

    fig = model.plot(spec, scalefac_resid=scalefac_resid)
    fig.subplots_adjust(hspace=0)
    return fig


if __name__ == "__main__":
    main()
