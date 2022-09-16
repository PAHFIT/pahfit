#!/usr/bin/env python

import argparse

import matplotlib.pyplot as plt
import matplotlib as mpl

from pahfit.model import Model
from pahfit.base import PAHFITBase
from pahfit.helpers import read_spectrum

from astropy import units as u


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
    spec = read_spectrum(args.spectrumfile, spec1d=True)

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
    """
    Returns
    -------
    fig : Figure object

    """

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

    fig, axs = plt.subplots(
        ncols=1,
        nrows=2,
        figsize=(15, 10),
        gridspec_kw={"height_ratios": [3, 1]},
        sharex=True,
    )
    PAHFITBase.plot(
        axs,
        spec.wavelength.to(u.micron).value,
        spec.flux,
        spec.uncertainty.array,
        model._construct_astropy_model(),
        scalefac_resid=scalefac_resid,
    )

    # use the whitespace better
    fig.subplots_adjust(hspace=0)
    return fig


if __name__ == "__main__":
    main()
