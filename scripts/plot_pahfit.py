#!/usr/bin/env python

import pkg_resources
import argparse

import matplotlib.pyplot as plt
import matplotlib as mpl

from astropy.io import fits

from pahfit.base import PAHFITBase


def initialize_parser():
    """
    Command line parser for plot_PAHFIT

    Returns
    -------
    parser : argparse object
    """
    plottypes = ['png', 'jpg', 'jpeg', 'pdf', 'ps', 'eps', 'rgba',
                 'svg', 'tiff', 'tif', 'pgf', 'svgz', 'raw']
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', help='name of PAHFIT results file')
    parser.add_argument('-f', '--savefig', action='store',
                        default=None, choices=plottypes,
                        help='Save figure to a file of specified type \
                            Must be one \
                            of: "{}"'.format('", "'.join(plottypes))
                        )
    return parser


if __name__ == "__main__":

    # commandline parser
    parser = initialize_parser()
    args = parser.parse_args()

    # read an observed spectrum

    # Example 1
    # read in an example spectrum (from M101)
    data_path = pkg_resources.resource_filename('pahfit',
                                                'data/')
    data = 'Nucleus_irs.fits'
    name = data.split('.')[0]
    hdul = fits.open('{}{}'.format(data_path, data))
    obs_x = hdul[1].data['WAVELENGTH']
    obs_y = hdul[1].data['FLUX']
    obs_unc = hdul[1].data['SIGMA']
    obs_npts = hdul[1].data['NPTS']
    hdul.close()
    weights = 1./obs_unc

    packfilename = pkg_resources.resource_filename(
        'pahfit',
        'packs/scipack_ExGal_SpitzerIRSSLLL.ipac')
    pmodel = PAHFITBase(filename=packfilename)

    # read in the PAHFIT results
    pmodel = PAHFITBase(filename=args.filename)

    # plot result
    fontsize = 18
    font = {'size': fontsize}
    mpl.rc('font', **font)
    mpl.rc('lines', linewidth=2)
    mpl.rc('axes', linewidth=2)
    mpl.rc('xtick.major', width=2)
    mpl.rc('ytick.major', width=2)

    fig, ax = plt.subplots(figsize=(15, 10))

    pmodel.plot(ax, obs_x, obs_y, pmodel.model)

    ax.set_yscale('linear')
    ax.set_xscale('log')

    # show or save
    if args.savefig:
        fig.savefig('{}.{}'.format(name, args.savefig), bbox_inches='tight')
    else:
        plt.show()
