#!/usr/bin/env python
#
# Python version of PAHFIT
#
# written Nov 2015 by Karl Gordon (kgordon@stsci.edu)
#
from __future__ import print_function

import numpy as np
import astropy.io.fits as pyfits
import matplotlib as mpl
import matplotlib.pyplot as plt
import argparse
import scipy.optimize as op
from lmfit import minimize, Parameters, fit_report, Model
from astropy.table import Table
from astropy.io import ascii
from astropy import units as u

# define input spectrum file
spectrum = 'NGC7023-NW-BRIGHT.txt'

# read in the spectrum to fit
a = Table.read(spectrum, format="ascii.commented_header")

indxs = np.where(np.logical_and(5.0 <= a['wave'].data,
                                a['wave'].data <= 28.8))

x = a['wave'].data[indxs]
y = a['flux'].data[indxs]
yerr = a['unc'].data[indxs] + y * 0.05

# define PAH feautures profile: ('lorentzian' or 'drude')
fitprof = 'lorentzian'

# global variables
temps = [3000., 2000., 1000., 800., 600., 450., 300., 200., 135., 90., 65., 50., 40., 35.]

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
h2_names = np.array(["H2_S(7)", "H2_S(6)", "H2_S(5)", "H2_S(4)",
                     "H2_S(3)", "H2_S(2)", "H2_S(1)", "H2_S(0)"])

ion_cwaves = np.array([6.985274, 8.99138, 10.5105, 12.813,
                       15.555, 18.713, 25.91, 25.989,
                       33.480, 34.8152, 35.349])

ion_names = np.array(["[ArII]", "[ArIII]", "[SIV]", "[NeII]",
                      "[NeIII]", "[SIII]_18", "[OIV]", "[FeII]",
                      "[SIII]_33", "[SiII]", "[FeII]"])


def ism_bb_MJysr(waves, T):
    return 3.97289e13 / waves ** 3 / (np.exp(1.4387752e4 / waves / T) - 1.)


def ism_gaussian(waves, cwave, fwhm):
    return np.exp(-((waves - cwave) ** 2) * 2.7725887222397811 /
                  (fwhm * cwave) ** 2)


def ism_lorentzian(waves, cwave, fwhm):
    g = (fwhm * cwave / 2.) ** 2
    return g / ((waves - cwave) ** 2 + g)


def ism_drude(waves, cwave, fwhm):
    return fwhm ** 2 / ((waves / cwave - cwave / waves) ** 2 + fwhm ** 2)


def ismfunc_vec_params(params):
    p = params.valuesdict()
    a1 = [p['bb' + str(k + 1)] for k in range(len(temps))]
    a2 = [p['df' + str(k + 1)] for k in range(len(dust_features_cwave))]
    a3 = [p['h2_iten' + str(k + 1)] for k in range(len(h2_cwaves))]
    a4 = [p['h2_width' + str(k + 1)] for k in range(len(h2_cwaves))]
    a5 = [p['ion_iten' + str(k + 1)] for k in range(len(ion_cwaves))]
    a6 = [p['ion_width' + str(k + 1)] for k in range(len(ion_cwaves))]
    return np.concatenate([a1, a2, a3, a4, a5, a6])


def ismfunc_cont(a, waves):
    # blackbodies
    ismmodel = np.zeros(len(waves))
    for i, ctemp in enumerate(temps):
        ismmodel += a[i] * ((9.7 / waves) ** 2) * ism_bb_MJysr(waves, ctemp)

    return ismmodel


def ismfunc(a, waves):
    # blackbodies
    ismmodel = np.zeros(len(waves))
    for i, ctemp in enumerate(temps):
        ismmodel += a[i] * ((9.7 / waves) ** 2) * ism_bb_MJysr(waves, ctemp)

    if fitprof == 'lorentzian':
        # dust features
        ioffset = len(temps)
        for i, cwave in enumerate(dust_features_cwave):
            ismmodel += a[i + ioffset] * ism_lorentzian(waves,
                                                        cwave,
                                                        dust_features_fwhm[i])

    elif fitprof == 'drude':
        # dust features
        ioffset = len(temps)
        for i, cwave in enumerate(dust_features_cwave):
            ismmodel += a[i + ioffset] * ism_drude(waves,
                                                   cwave,
                                                   dust_features_fwhm[i])

    # h2 lines
    ioffset += len(dust_features_cwave)
    eoffset = len(h2_cwaves)
    for i, cwave in enumerate(h2_cwaves):
        ismmodel += a[i + ioffset] * ism_gaussian(waves,
                                                  cwave,
                                                  a[i + ioffset + eoffset])
    # ionized lines
    ioffset += 2 * eoffset
    eoffset = len(ion_cwaves)
    for i, cwave in enumerate(ion_cwaves):
        ismmodel += a[i + ioffset] * ism_gaussian(waves,
                                                  cwave,
                                                  a[i + ioffset + eoffset])

    return ismmodel


# def ismfunc_cont_residuals(params, x, y, yerr, **kws):
#     a = ismfunc_vec_params(params)
#     model = ismfunc_cont(a, x)
#     return (y - model)/yerr

def ismfunc_residuals(params, x, y, yerr, **kws):
    a = ismfunc_vec_params(params)
    model = ismfunc(a, x)
    return (y - model) / yerr


def lnlike(a, x, y, yerr):
    model = ismfunc(a, x)
    ans = -0.5 * np.sum((y - model) ** 2 / (yerr ** 2))
    return ans


# lmfit setup
params = Parameters()
for k, ctemp in enumerate(temps):
    if k <= 2:
        val = 1e-12
    elif k <= 5:
        val = 1e-10
    else:
        val = 1e-7
    params.add('bb' + str(k + 1), value=val, min=0.0)

for k, cwave in enumerate(dust_features_cwave):
    params.add('df' + str(k + 1), value=100., min=0.0)

for k, cwave in enumerate(h2_cwaves):
    params.add('h2_iten' + str(k + 1), value=1e2, min=0.0)
    params.add('h2_width' + str(k + 1), value=cwave / 3000.,
               min=cwave / 6000., max=cwave / 1000.)

for k, cwave in enumerate(ion_cwaves):
    params.add('ion_iten' + str(k + 1), value=1e2, min=0.0)
    params.add('ion_width' + str(k + 1), value=cwave / 2500.,
               min=cwave / 4500., max=cwave / 1000.)

# running lmfit
print('running lmfit, please wait...')

out = minimize(ismfunc_residuals, params, args=(x, y, yerr))

# get the residual array
res = out.residual

# calculate the best fit model
best_model = [k - h for k, h in zip(y, res)]

# open output files
o = open('output.txt', 'w')
bm = open('best_model.txt', 'w')

# print and write fit results
print(fit_report(out))
o.write(str(fit_report(out)))

# create best-fit and residuals file
bm.write('#wave\t\tfit\t\t\tres\n')

for i, line in enumerate(best_model):
    bm.write(str(x[i]) + '\t' + str(best_model[i]) + '\t' + str(res[i]) + '\n')

bm.close()
o.close()
print('fit completed')

# reading output file
bb, bb_val, df, df_val, h2_iten, h2_iten_val, h2_width, h2_width_val, ion_iten, \
ion_iten_val, ion_width, ion_width_val = ([] for i in range(12))

with open('output.txt') as out:
    temp = out.readlines()[9:]
    params = [word.strip() for word in temp]
    for i, line in enumerate(params):
        if 'bb' in line:
            temp_val = line.split()
            bb.append(temp_val[0].strip(':'))
            bb_val.append(float(temp_val[1]))
        if 'df' in line:
            temp_val = line.split()
            df.append(temp_val[0].strip(':'))
            df_val.append(float(temp_val[1]))
        if 'h2_iten' in line:
            temp_val = line.split()
            h2_iten.append(temp_val[0].strip(':'))
            h2_iten_val.append(float(temp_val[1]))
        if 'h2_width' in line:
            temp_val = line.split()
            h2_width.append(temp_val[0].strip(':'))
            h2_width_val.append(float(temp_val[1]))
        if 'ion_iten' in line:
            temp_val = line.split()
            ion_iten.append(temp_val[0].strip(':'))
            ion_iten_val.append(float(temp_val[1]))
        if 'ion_width' in line:
            temp_val = line.split()
            ion_width.append(temp_val[0].strip(':'))
            ion_width_val.append(float(temp_val[1]))

# calculating best-fit H2 gaussians
h2_dict = {}
for i in range(len(h2_cwaves)):
    h2_key = h2_names[i]
    h2_dict_val = [x * h2_iten_val[i] for x in ism_gaussian(x, h2_cwaves[i], h2_width_val[i])]
    h2_dict[h2_key] = h2_dict_val

# calculating best-fit ion gaussians
ion_dict = {}
for i in range(len(ion_cwaves)):
    ion_key = ion_names[i]
    ion_dict_val = [x * ion_iten_val[i] for x in ism_gaussian(x, ion_cwaves[i], ion_width_val[i])]
    ion_dict[ion_key] = ion_dict_val

# calculating best-fit PAH features
if fitprof == 'lorentzian':
    df_dict = {}
    for i in range(len(dust_features_cwave)):
        ind = str(i + 1).zfill(2)
        df_key = 'df' + str(ind)
        df_dict_val = [x * df_val[i] for x in ism_lorentzian(x, dust_features_cwave[i], dust_features_fwhm[i])]
        df_dict[df_key] = df_dict_val

if fitprof == 'drude':
    df_dict = {}
    for i in range(len(dust_features_cwave)):
        ind = i + 1
        df_key = 'df' + str(ind)
        df_dict_val = [x * df_val[i] for x in ism_drude(x, dust_features_cwave[i], dust_features_fwhm[i])]
        df_dict[df_key] = df_dict_val

# calculating BBs
bb_dict = {}
for i in range(len(temps)):
    ind = str(i + 1).zfill(2)
    bb_key = 'bb' + str(ind)
    bb_dict_val = [x * bb_val[i] for x in ism_bb_MJysr(x, temps[i])]
    bb_dict[bb_key] = bb_dict_val

# obtaining the integrated flux of all features
inf = open('integ_flux.txt', 'w')
inf.write('#feature\tc_wave\t\tint_flux\n')

# integrated flux of PAHs
for index, key, in enumerate(sorted(df_dict.items(), key=lambda x: x[0])):
    temp = key[1]
    new_value = []
    for j in range(0, len(temp) - 1):
        temp_flux = (10 ** 6 * temp[j]) * u.Jy
        temp_wave = a['wave'].data[j]
        new_flux = temp_flux.to(u.Watt / u.m ** 2 / u.micron,
                                equivalencies=u.spectral_density(temp_wave * u.micron))
        new_value.append(new_flux.value)
    inf.write(key[0] + '\t\t' + str(dust_features_cwave[index]) + '\t\t' + '%.3e' % np.sum(new_value) + '\n')

# integrated flux of ions
for key, value in sorted(ion_dict.items(), key=lambda x: x[0]):
    ind = np.where(ion_names == key)
    wave = ion_cwaves[ind][0]
    temp = value
    new_value = []
    for j in range(0, len(temp) - 1):
        temp_flux = (10 ** 6 * temp[j]) * u.Jy
        temp_wave = a['wave'].data[j]
        new_flux = temp_flux.to(u.Watt / u.m ** 2 / u.micron,
                                equivalencies=u.spectral_density(temp_wave * u.micron))
        new_value.append(new_flux.value)
    inf.write(key + '\t\t' + str(wave) + '\t\t' + '%.3e' % np.sum(new_value) + '\n')

# integrated flux of H2
for key, value in sorted(h2_dict.items(), key=lambda x: x[0], reverse=True):
    ind = np.where(h2_names == key)
    wave = h2_cwaves[ind]
    temp = value
    new_value = []
    for j in range(0, len(temp) - 1):
        temp_flux = (10 ** 6 * temp[j]) * u.Jy
        temp_wave = a['wave'].data[j]
        new_flux = temp_flux.to(u.Watt / u.m ** 2 / u.micron,
                                equivalencies=u.spectral_density(temp_wave * u.micron))
        new_value.append(new_flux.value)
    inf.write(key + '\t\t' + str(wave[0]) + '\t\t' + '%.3e' % np.sum(new_value) + '\n')

# integrated flux of black bodies
for index, key, in enumerate(sorted(bb_dict.items(), key=lambda x: x[0])):
    temp = key[1]
    new_value = []
    for j in range(0, len(temp) - 1):
        temp_flux = (10 ** 6 * temp[j]) * u.Jy
        temp_wave = a['wave'].data[j]
        new_flux = temp_flux.to(u.Watt / u.m ** 2 / u.micron,
                                equivalencies=u.spectral_density(temp_wave * u.micron))
        new_value.append(new_flux.value)
    inf.write(key[0] + '\t\t' + str(temps[index]) + '\t\t' + '%.3e' % np.sum(new_value) + '\n')

inf.close()

# reading best fit model
bfm = np.genfromtxt('best_model.txt', dtype=None, names=True)
best_model = bfm['fit']

fontsize = 21
font = {'size': fontsize}
mpl.rc('font', **font)
mpl.rc('lines', linewidth=1.7)
mpl.rc('axes', linewidth=2)
mpl.rc('xtick.major', width=2)
mpl.rc('ytick.major', width=2)
mpl.rc('xtick.major', size=5)
mpl.rc('ytick.major', size=5)
mpl.rc('xtick.minor', size=3)
mpl.rc('ytick.minor', size=3)

fig, many_ax = plt.subplots(ncols=1, nrows=2, figsize=(15, 10), sharex=True)

# spectrum and best fit model
ax = many_ax[0]
ax.minorticks_on()
ax.tick_params(which='major', top='on', direction='in')
ax.tick_params(which='minor', top='on', direction='in')
ax.plot(x, best_model, color='darkgreen', label='Model')
ax.errorbar(x, y, yerr=yerr, fmt='o', markeredgecolor='k', markerfacecolor='none', ecolor='k', markersize=6)

# plotting PAHs and generating total PAH spectrum
pah_sum = np.zeros(x.shape)
for i, (key, value) in enumerate(df_dict.items()):
    pah_sum += value
    if i == len(df_dict) - 1:
        ax.plot(x, value, color='royalblue', label='PAHs')
    else:
        ax.plot(x, value, color='royalblue')

# saving total PAH spectrum to file
ascii.write([x, pah_sum], '%s_pahs.txt' % spectrum, names=['wavelength', 'surface brightness'], overwrite=True)

# plotting H2
for i, (key, value) in enumerate(h2_dict.items()):
    if i == len(h2_dict) - 1:
        ax.plot(x, value, color='mediumspringgreen', label='H$_{2}$')
    else:
        ax.plot(x, value, color='mediumspringgreen')

# plotting ions
for i, (key, value) in enumerate(ion_dict.items()):
    if i == len(ion_dict) - 1:
        ax.plot(x, value, color='darkmagenta', label='Ions')
    else:
        ax.plot(x, value, color='darkmagenta')

# plotting BBs
for i, (key, value) in enumerate(bb_dict.items()):
    if i == len(bb_dict) - 1:
        ax.plot(x, value, color='red', label='dust')
    else:
        ax.plot(x, value, color='red')

ax.set_ylabel('Flux [MJy/sr]')

# top panel legend
ax.legend(loc=1)

# residuals
ax = many_ax[1]
ax.minorticks_on()
ax.tick_params(which='major', top='on', direction='in')
ax.tick_params(which='minor', top='on', direction='in')
ax.plot(x, y - best_model, color='k', label='Residuals')
ax.set_xlim(x[0], x[-1])
ax.set_xlabel(r'$\lambda$ [$\mu$m]')
ax.set_ylabel('Residuals [MJy/sr]')

# bottom panel legend
ax.legend(loc=1)

# adjust subplot space
fig.subplots_adjust(hspace=0)

# saving and showing plot
plt.savefig('%s.pdf' % spectrum, bbox_inches='tight')
plt.show()
print('plot created')