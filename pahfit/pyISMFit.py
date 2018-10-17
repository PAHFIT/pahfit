#!/usr/bin/env python
#
# Python version of PAHFIT
#
# written Nov 2015 by Karl Gordon (kgordon@stsci.edu)
#
from __future__ import print_function

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import argparse
import matplotlib

#print(plt.get_backend())

import scipy.optimize as op
from lmfit import minimize, Parameters

from astropy.table import Table

# global variables
temps = [3000.,2000.,1000.,800.,600.,450.,300.,200.,135.,90.,65.,50.,40.,35.]

dust_features_cwave = np.array([5.27, 5.70, 6.22, 6.69,
                                7.42, 7.60, 7.85, 8.33,
                                8.61, 10.68, 11.23, 11.33,
                                11.99, 12.62, 12.69, 13.48,
                                14.04, 14.19, 15.9, 16.45,
                                17.04, 17.375, 17.87, 18.92,
                                33.1])
                                
dust_features_fwhm = np.array([0.034, 0.035, 0.030, 0.07,
                               0.126, 0.044, 0.053, 0.05,
                               0.039, 0.02, 0.012, 0.032,
                               0.045, 0.042, 0.013, 0.04,
                               0.016, 0.025, 0.02, 0.014,
                               0.065, 0.012, 0.016, 0.019,
                               0.05])

h2_cwaves = np.array([5.5115, 6.1088, 6.9091, 8.0258,
                      9.6649, 12.2785, 17.0346, 28.2207])
h2_names = np.array(["H2 S(7)", "H2 S(6)", "H2 S(5)", "H2 S(4)",
                     "H2 S(3)", "H2 S(2)", "H2 S(1)", "H2 S(0)"])

ion_cwaves = np.array([6.985274, 8.99138, 10.5105, 12.813,
                       15.555, 18.713, 25.91, 25.989,
                       33.480, 34.8152, 35.349])

ion_names = np.array(["[ArII]", "[ArIII]", "[SIV]", "[NeII]",
                      "[NeIII]", "[SIII] 18", "[OIV]", "[FeII]", 
                      "[SIII] 33", "[SiII]", "[FeII]"])

def ism_bb_MJysr(waves, T):
    return 3.97289e13/waves**3/(np.exp(1.4387752e4/waves/T)-1.)

def ism_gaussian(waves, cwave, fwhm):
    return np.exp(-((waves-cwave)**2)*2.7725887222397811/ 
                  (fwhm*cwave)**2)

def ism_lorentzian(waves, cwave, fwhm):
  g = (fwhm*cwave/2.)**2
  return g/((waves-cwave)**2+g)

def ismfunc_vec_params(params):
    p = params.valuesdict()
    a1 = [p['bb' + str(k+1)] for k in range(len(temps))]
    a2 = [p['df' + str(k+1)] for k in range(len(dust_features_cwave))]
    a3 = [p['h2_iten' + str(k+1)] for k in range(len(h2_cwaves))]
    a4 = [p['h2_width' + str(k+1)] for k in range(len(h2_cwaves))]
    a5 = [p['ion_iten' + str(k+1)] for k in range(len(ion_cwaves))]
    a6 = [p['ion_width' + str(k+1)] for k in range(len(ion_cwaves))]
    return np.concatenate([a1,a2,a3,a4,a5,a6])

def ismfunc_cont(a, waves):
    # backbodies
    ismmodel = np.zeros(len(waves))
    for i, ctemp in enumerate(temps):
        ismmodel += a[i]*((9.7/waves)**2)*ism_bb_MJysr(waves, ctemp)
    
    return ismmodel

def ismfunc(a, waves):
    # backbodies
    ismmodel = np.zeros(len(waves))
    for i, ctemp in enumerate(temps):
        ismmodel += a[i]*((9.7/waves)**2)*ism_bb_MJysr(waves, ctemp)
    
    # dust features
    ioffset = len(temps)
    for i, cwave in enumerate(dust_features_cwave):
        ismmodel += a[i+ioffset]*ism_lorentzian(waves,
                                                cwave,
                                                dust_features_fwhm[i])
    # h2 lines
    ioffset += len(dust_features_cwave)
    eoffset = len(h2_cwaves)
    for i, cwave in enumerate(h2_cwaves):
        ismmodel += a[i+ioffset]*ism_gaussian(waves,
                                              cwave,
                                              a[i+ioffset+eoffset])
    # ionized lines
    ioffset += 2*eoffset
    eoffset = len(ion_cwaves)
    for i, cwave in enumerate(ion_cwaves):
        ismmodel += a[i+ioffset]*ism_gaussian(waves,
                                              cwave,
                                              a[i+ioffset+eoffset])

    return ismmodel

def ismfunc_cont_residuals(params, x, y, yerr, **kws):
    a = ismfunc_vec_params(params)
    model = ismfunc_cont(a, x)
    return (y - model)/yerr

def ismfunc_residuals(params, x, y, yerr, **kws):
    a = ismfunc_vec_params(params)
    model = ismfunc(a, x)
    return (y - model)/yerr

def ismfunc_plot(params, niter, resid, *args, **kws):
    if (niter % 50) == 0:
        model = ismfunc(ismfunc_vec_params(params), args[0])
        kws['mod_line'].set_ydata(model)
        kws['mod_resid'].set_ydata(args[1] - model)
        kws['many_ax'][1].relim()
        kws['many_ax'][1].autoscale_view()
        kws['many_ax'][1].set_ylim(-50.,50.)
        print(niter, 0.5*np.sum(resid**2)/len(args[0]))
        #plt.show()        
        plt.draw()
        plt.pause(0.001)

def lnlike(a, x, y, yerr):
    model = ismfunc(a, x)
    ans = -0.5*np.sum((y-model)**2/(yerr**2))
    return ans

if __name__ == "__main__":

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("filename",help="file with ISO SWS spectrum")
    parser.add_argument("--png", help="save figure as a png file",
                        action="store_true")
    parser.add_argument("--eps", help="save figure as an eps file",
                        action="store_true")
    parser.add_argument("-pdf", help="save figure as a pdf file",
                        action="store_true")
    args = parser.parse_args()

    # read in the spectrum to fit
    a = Table.read(args.filename, format="ascii.commented_header")

    indxs = np.where(np.logical_and(5.0 <= a['wave'].data,
                                    a['wave'].data <= 28.8))
    x = a['wave'].data[indxs]
    y = a['flux'].data[indxs]
    yerr = a['unc'].data[indxs] + y*0.05

    # scipy optimze
    #nll = lambda *args: -lnlike(*args)
    #a = np.array([1.0, 1.0])*1e-7
    #result = op.minimize(nll, a, args=(x, y, yerr))
    #print(result["x"])

    # lmfit setup
    params = Parameters()
    for k, ctemp in enumerate(temps):
        if k <= 2:
            val = 1e-12
        elif k <= 5:
            val = 1e-10
        else:
            val = 1e-7
        params.add('bb' + str(k+1),value=val,min=0.0)

    for k, cwave in enumerate(dust_features_cwave):
        params.add('df' + str(k+1),value=100.,min=0.0)

    for k, cwave in enumerate(h2_cwaves):
        params.add('h2_iten' + str(k+1),value=1e2,min=0.0)
        params.add('h2_width' + str(k+1),value=cwave/3000.,
                   min=cwave/6000.,max=cwave/1000.)

    for k, cwave in enumerate(ion_cwaves):
        params.add('ion_iten' + str(k+1),value=1e2,min=0.0)
        params.add('ion_width' + str(k+1),value=cwave/2500.,
                   min=cwave/4500.,max=cwave/1000.)
        
    # plot stuff
    fontsize = 18
    font = {'size'   : fontsize}
    matplotlib.rc('font', **font)
    matplotlib.rc('lines', linewidth=2)
    matplotlib.rc('axes', linewidth=2)
    matplotlib.rc('xtick.major', width=2)
    matplotlib.rc('ytick.major', width=2)

    fig, many_ax = plt.subplots(ncols=1, nrows=2, figsize=(15,10))

    # spectrum and current model
    ax = many_ax[0]
    ax.plot(x,y)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(2,50.)
    ax.set_xlabel(r'$\lambda$ [$\mu$m]')
    ax.set_ylabel('Flux [Jy]')

    model = ismfunc(ismfunc_vec_params(params), x)
    mod_line, = ax.plot(x, model)

    # residuals
    ax = many_ax[1]
    mod_resid, = ax.plot(x,y-model)
    ax.set_xscale('log')
    ax.set_xlim(2,50.)
    ax.set_xlabel(r'$\lambda$ [$\mu$m]')
    ax.set_ylabel('Residuals [Jy]')
    
    plt.tight_layout()    
    #plt.show()
    
    plt.ion()
    plt.draw()
    plt.pause(0.001)

    # run lmfit
    ismfunc_plot_keywords = {'mod_line': mod_line, 'mod_resid': mod_resid,
                             'many_ax': many_ax}
    out = minimize(ismfunc_residuals, params,
                   args=(x, y, yerr),
                   iter_cb=ismfunc_plot,
                   kws=ismfunc_plot_keywords)
    print(out.params)

    # show or save
    basename = args.filename
    basename.replace('.txt','')
    if args.png:
        fig.savefig(basename+'.png')
    elif args.eps:
        fig.savefig(basename+'.eps')
    elif args.pdf:
        fig.savefig(basename+'.pdf')
    else:
        plt.show()
    
    
