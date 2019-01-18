from __future__ import (absolute_import, print_function, division)

import pkg_resources

import matplotlib.pyplot as plt
import matplotlib

from astropy.io import fits
from astropy.modeling.fitting import LevMarLSQFitter

from PAHFITBase import PAHFITBase
from PAHFIT_Spitzer_Exgal import (InstPackSpitzerIRSSLLL, SciPackExGal)


if __name__ == '__main__':
    # would be good to use argparse to allow for command line input
    # of the observed spectra.  Need to standardize the expected format
    # to allow this.

    # instrument pack
    spitzer_sl_ll_pack = InstPackSpitzerIRSSLLL()

    # science pack
    sci_pack = SciPackExGal(spitzer_sl_ll_pack)

    pmodel = PAHFITBase(bb_info=sci_pack.bb_info,
                        dust_features=sci_pack.dust_features,
                        h2_features=sci_pack.h2_features,
                        ion_features=sci_pack.ion_features)

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
    # define output format (ascii, fits, csv, etc.)
    outform = 'fits'
    pmodel.save(obs_fit, name, outform)

    # plot result
    fontsize = 18
    font = {'size': fontsize}
    matplotlib.rc('font', **font)
    matplotlib.rc('lines', linewidth=2)
    matplotlib.rc('axes', linewidth=2)
    matplotlib.rc('xtick.major', width=2)
    matplotlib.rc('ytick.major', width=2)

    fig, ax = plt.subplots(figsize=(15, 10))

    pmodel.plot(ax, obs_x, obs_y, obs_fit)

    ax.set_yscale('linear')
    ax.set_xscale('log')

    plt.show()
    plt.savefig('{}.pdf'.format(name), bbox_inches='tight')
