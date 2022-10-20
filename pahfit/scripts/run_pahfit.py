#!/usr/bin/env python

import argparse

import matplotlib.pyplot as plt

from pahfit.helpers import read_spectrum
from pahfit.model import Model
from pahfit.scripts.plot_pahfit import default_layout_plot


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
    instrumenttypes = [
        "spitzer.irs.sl.1",
        "spitzer.irs.sl.2",
        "spitzer.irs.sl.3",
        "spitzer.irs.ll.1",
        "spitzer.irs.ll.2",
        "spitzer.irs.ll.3",
        "spitzer.irs.sh",
        "spitzer.irs.lh",
        "iso.sws.speed0",
        "iso.sws.speed1",
        "iso.sws.speed2",
        "iso.sws.speed3",
        "iso.sws.speed4",
    ]
    parser = argparse.ArgumentParser()
    parser.add_argument("spectrumfile", help="name of file with observed spectrum")
    parser.add_argument("packfile", help="name of PAHFIT pack file")
    parser.add_argument(
        "instrumentname",
        # choices=instrumenttypes,
        help="Name of the instrument. Available:" + str(instrumenttypes),
    )
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
        default="ascii.ecsv",
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
    spec = read_spectrum(args.spectrumfile)

    # setup the model
    model = Model.from_yaml(args.packfile, args.instrumentname, 0)

    # initial guess if not explicitly disabled
    if not args.no_starting_estimate:
        model.guess(spec)

    # fit the spectrum
    model.fit(spec, maxiter=args.fit_maxiter)

    # save fit results to file
    basename = args.spectrumfile.split(".")[0]
    extension = args.saveoutput
    outputname = f"{basename}_output.{extension}"
    print("Writing result to ", outputname)
    model.save(outputname)

    fig = default_layout_plot(spec, model, args.scalefac_resid)

    # show
    if args.showplot:
        plt.show()
    # save (always)
    fig.savefig("{}.{}".format(outputname, args.savefig))


if __name__ == "__main__":
    main()
