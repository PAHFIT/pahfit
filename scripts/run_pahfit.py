#!/usr/bin/env python

import pkg_resources
import argparse

import matplotlib.pyplot as plt
import matplotlib as mpl

from astropy.io import fits
from astropy.modeling.fitting import LevMarLSQFitter

from pahfit.base import PAHFITBase
from pahfit.PAHFIT_Spitzer_Exgal import (InstPackSpitzerIRSSLLL, SciPackExGal)


def initialize_parser():
    """
    Command line parser for PAHFIT

    Returns
    -------
    parser : argparse object
    """
    plottypes = ['png', 'jpg', 'jpeg', 'pdf', 'ps', 'eps', 'rgba',
                 'svg', 'tiff', 'tif', 'pgf', 'svgz', 'raw']
    savetypes = ['fits', 'votable', 'ipac', 'ascii.ecsv']
    parser = argparse.ArgumentParser()
    parser.add_argument('packfile', help='name of PAHFIT pack file')
    parser.add_argument('-f', '--savefig', action='store',
                        default='pdf', choices=plottypes,
                        help='Save figure to a file of specified type \
                            Must be one \
                            of: "{}"'.format('", "'.join(plottypes))
                        )
    parser.add_argument('-o', '--saveoutput', action='store',
                        default='ipac', choices=savetypes,
                        help='Save fit results to a file of specified type \
                            Must be one \
                            of: "{}"'.format('", "'.join(savetypes))
                        )
    return parser


if __name__ == '__main__':

    parser = initialize_parser()
    args = parser.parse_args()

    # would be good to use argparse to allow for command line input
    # of the observed spectra.  Need to standardize the expected format
    # to allow this.

    # instrument pack
    spitzer_sl_ll_pack = InstPackSpitzerIRSSLLL()

    # science pack
    sci_pack = SciPackExGal(spitzer_sl_ll_pack)

    param_info = (sci_pack.bb_info, sci_pack.dust_features,
                  sci_pack.h2_features, sci_pack.ion_features,
                  sci_pack.att_info)

    opmodel = PAHFITBase(param_info=param_info)
    opmodel.save(opmodel.model, 'Pack_ExGal_SpitzerIRSSLLL', args.saveoutput)

    # pmodel = PAHFITBase(filename=args.packfile)
    pmodel = opmodel

    # pick the fitter
    fit = LevMarLSQFitter()

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

    # Example 2 (from Thomas - not sure the galaxy)
    # import pandas as pd
    # df = pd.read_table(data_path+'1120229'+'.txt',
    #                    header=None,
    #                    names=['w', 'f', 'f_l', 'f_h'])
    # obs_x = df.w.values
    # obs_y = df.f.values
    # weights = None

    obs_fit = fit(pmodel.model, obs_x, obs_y, weights=weights,
                  maxiter=1000)
    print(fit.fit_info['message'])

    # save results to fits file
    pmodel.save(obs_fit, name, args.saveoutput)

    # plot result
    fontsize = 18
    font = {'size': fontsize}
    mpl.rc('font', **font)
    mpl.rc('lines', linewidth=2)
    mpl.rc('axes', linewidth=2)
    mpl.rc('xtick.major', width=2)
    mpl.rc('ytick.major', width=2)

    fig, ax = plt.subplots(figsize=(15, 10))

    pmodel.plot(ax, obs_x, obs_y, obs_fit)

    ax.set_yscale('linear')
    ax.set_xscale('log')

    # use the whitespace better
    fig.tight_layout()

    # show
    plt.show()
    # and save
    fig.savefig('{}.{}'.format(name, args.savefig))
