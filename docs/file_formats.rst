############
File Formats
############

Model Packs
===========

YAML form
---------

The main way for the user to set up models, is by using (or copying and editing) one of the default science packs provided as YAML files.
These YAML files are internally converted to an astropy Table (see next section `Table form`_).

YAML is a data serialization language.
The model packs are constructed as lists of dictionaries of different feature components (eg., Line Features, Dust Features).
Each list consists of dictionaries which define the feature and the parameters associated with it.
For example ::

 PAH_5.3:
     kind : dust_feature
     wavelength: 5.27
     fwhm: 0.17918

This shows a single PAH feature and its corresponding attributes (wavelength, fwhm).
The parameters specified depend on the kind of the feature and its respective functional form.
Bounds (min and max) can also be set for different parameters as well as for all the parameters of a certain group.

The fit parameters for continua (`starlight` and `dust_continuum`) are temperature and tau (emission strength).
For features like `line` and `dust_feature`, the fit parameters are wavelength, fwhm and power.
For `attenuation`, the fit parameters are tau, model, and geometry, while `absorption` features have tau, wavelength, fwhm and geometry.

See :ref:`Example Input <example_input>` for an example of YAML input model file.
For a detailed account of the YAML based science packs and the relevant syntax, please follow the `Wiki <https://github.com/PAHFIT/pahfit/wiki/PAHFIT-2022>`_.

Table form
----------

When a model is set up from a YAML file, or when the fit results are saved to disk, the model is described by a table compatible with astropy, called the "feature table".
This table has one row for each model component.
Each row has a value for the `name` (e.g., BB3, DF10, H2 S(3), [OIV]), `group`, and `kind` columns.
The name and group columns come from the YAML file, and can be any string, as defined by the user in the YAML file from which this table was generated.
The kind column indicates one of the supported component types, and determines the function by which the component wil be represented (e.g., BlackBody1D for dust_continuum, Drude1D for dust_feature, Gaussian1D for line).
Following these three columns, there are multiple columns for the fit parameters.
These columns contain 3-tuples with the value and the bounds: `[value, min, max]`.
The bounds `min` and `max` are masked if a parameter is fixed.
If a parameter is not relevant for a certain feature, then `value` will be masked.
To reduce the number of columns, some parameters have different meaning depending on the kind of the feature.
For the dust_continuum components, the fit parameters are temperature (`temperature`) and the "intrinsic" or "unattenuated" amplitude (`tau`).
For the dust_feature and (gas) line models, the fit parameters are "intrinsic" or "unattenuated" amplitude (`power`), center wavelength (`wavelength`), and full width at half maximum (`fwhm`).
For the attenuation model, the fit is a measure for the optical depth (`tau`).

The model tables can be stored as ECSV tables as these provide readable (ASCII) files with metadata support.
The ECSV metadata is used to read in a model from file, after a fit result has already been saved.

The astropy version of `Drude1D <https://docs.astropy.org/en/stable/modeling/physical_models.html#drude1d>`_ and `Gaussian1D <https://docs.astropy.org/en/stable/api/astropy.modeling.functional_models.Gaussian1D.html#astropy.modeling.functional_models.Gaussian1D>`_ are used under the hood.
See :ref:`reference API <reference_API>` for description of BlackBody1D and the dust attenuation model.

See :ref:`example fit output <example_fit_output>` for an example of
a feature table/output file.

.. note::
   IDL PAHFIT uses fractional FWHM as input/output for the Drude profiles, python
   PAHFIT uses FWHM as input/output for the Drude profiles.

.. _example_input:

Example Input (YAML)
--------------------
::

 # PAHFIT Classic Model Pack
 # Implements the IDL-based PAHFIT v1.2 model
 # v0.1, May 2022

 #######################
 # Starlight continuum #
 #######################
 starlight:
     kind: starlight_continuum
     temperature: 5000

 ##################
 # Dust continuum #
 ##################
 dust_cont: # Modified Blackbodies
     kind: dust_continuum
     temperature: [300, 200, 135, 90, 65, 50, 40, 35]

 ##############################
 # H_2 Line Emission Features #
 ##############################
 H2_lines:
     kind: line
     wavelength:
         H2_S(7):    5.5115
         H2_S(6):    6.1088
         H2_S(5):    6.9091
         H2_S(4):    8.0258
         H2_S(3):    9.6649
         H2_S(2):   12.2785
         H2_S(1):   17.0346
         H2_S(0):   28.2207

 ################################
 # Ionic Line Emission Features #
 ################################
 ionic_lines:
     kind: line
     wavelength:
         '[ArII]':     6.985274
         '[ArIII]':    8.99138
         '[SIV]':     10.5105
         '[NeII]':    12.813
         '[NeIII]':   15.555
         '[SIII]_18': 18.713
         '[OIV]':     25.91
         '[FeII]':    25.989
         '[SIII]_33': 33.480
         '[SiII]':    34.8152
         '[FeII]_35':    35.349

 #################
 # Dust Features #
 #################
 PAH_5.3:
     kind: dust_feature
     wavelength: 5.27
     fwhm: 0.17918

 PAH_5.7:
     kind: dust_feature
     wavelength: 5.7
     fwhm: 0.1995

 PAH_6.2:
     kind: dust_feature
     wavelength: 6.22
     fwhm: 0.1866

 PAH_6.7:
     kind: dust_feature
     wavelength: 6.69
     fwhm: 0.4683

 PAH_7.7_cmp:
     kind: dust_feature
     features:
         PAH_7.7a:
             wavelength: 7.42
             fwhm: 0.93492
         PAH_7.7b:
             wavelength: 7.6
             fwhm: 0.3344
         PAH_7.7c:
             wavelength: 7.85
             fwhm: 0.41605

 PAH_8.3:
     kind: dust_feature
     wavelength: 8.33
     fwhm: 0.4165

 PAH_8.6:
     kind: dust_feature
     wavelength: 8.61
     fwhm: 0.33579

 PAH_10.7:
     kind: dust_feature
     wavelength: 10.68
     fwhm: 0.2136

 PAH_11.3_cmp:
     kind: dust_feature
     features:
         PAH_11.3a:
             wavelength: 11.23
             fwhm: 0.13476
         PAH_11.3b:
             wavelength: 11.33
             fwhm: 0.36256

 PAH_12:
     kind: dust_feature
     wavelength: 11.99
     fwhm: 0.53955

 PAH_12.6_cmp:
     kind: dust_feature
     features:
         PAH_12.6a:
             wavelength: 12.62
             fwhm: 0.53004
         PAH_12.6b:
             wavelength: 12.69
             fwhm: 0.16497

 PAH_13.48:
     kind: dust_feature
     wavelength: 13.48
     fwhm: 0.5392

 PAH_14.04:
     kind: dust_feature
     wavelength: 14.04
     fwhm: 0.22464

 PAH_14.19:
     kind: dust_feature
     wavelength: 14.19
     fwhm: 0.35475

 PAH_15.9:
     kind: dust_feature
     wavelength: 15.9
     fwhm: 0.318

 PAH_17_cmp:
     kind: dust_feature
     features:
         PAH_17a:
             wavelength: 16.45
             fwhm: 0.2303
         PAH_17b:
             wavelength: 17.04
             fwhm: 1.1076
         PAH_17c:
             wavelength: 17.37
             fwhm: 0.2085
         PAH_17d:
             wavelength: 17.87
             fwhm: 0.28592

 # This dust feature, in the PAHFIT classic model, is attributed to C60 and omitted.
 #PAH_18.92:
 #    kind: dust_feature
 #    wavelength: 18.92
 #    fwhm: 0.35948

 PAH_33.1:
     kind: dust_feature
     wavelength: 33.1
     fwhm: 1.655

 ##########################
 # Attenuation Model      #
 ##########################
 silicate:
     kind: attenuation
     model: S07_attenuation
     geometry: mixed

Instrument packs
================

Rationale
---------

While the physical parameters are provided by the science packs, the observed emission lines will be spectroscopically unresolved for certain instruments.
Therefore, the width of the lines will be determined mainly by instrumental effects.
For this purpose, instrument packs have been developed, which describe the spectral resolution as a function of wavelength.
Each instrument pack (usually pertaining to one observatory) contains multiple instrument models, typically one per spectral segment.
See the example below. One or multiple instrument configurations need to be provided when setting up a fit, either on the command line or when a Model instance is set up (see <fit_spectrum>).

Example file
------------

The user can find the instrument configurations inspecting the contents of each instrument pack.
For example, the Spitzer pack (spitzer.yaml) is shown below.
It has a hierarchical structure to define the possible instrument configurations.
To use a single segment configuration, the user needs to provide a string containing the name of the science pack, followed by the names in the hierarchy separated by dots,
e.g. `"spitzer.irs.sl.1"`.
Because most spectroscopic obervations will combine multiple of those configurations, a list of strings or wildcard string can be passed as well.
For example, if observations were taken with the SL1, SL2, LL1, and LL2 modes, then the instrument configuration can be writen down as `"spitzer.irs.*.[12]"`.
The asterisk will match any name at level 2, while `[12]` will match either '1' or '2' at level 3.

.. code-block:: yaml

  # ---
  # PAHFIT Instrument Pack, see https://github.com/PAHFIT/pahfit/wiki/File-Formats for more information

  # This is the Instrument Pack for Spitzer observations using staring
  # or spectral mapping mode.  For each of the four instrument
  # (Short-Low = SL, Long-Low = LL, Short-High = SH, Long-High = LH,
  # with LL/SL divided into order/subslits denoted 1, 2, and 3) this
  # file first provides the wavelength range covered (in microns).  The
  # resolving power (lambda / delta_lambda, where delta_lamba is the
  # FWHM of the resolution element) is represented by a polynomial of
  # degree 3, and this file lists the coefficients CC of that
  # polynomial, i.e.
  #
  # R = CC[0] + CC[1] * lambda + CC[2] * lambda^2 + CC[3] * lambda^3
  #
  # where lambda is expressed in microns as well.
  #
  # Resolution data from the IRS Instrument Handbook v5.0, Section 4.1.3
  # (see Figures 4.8 and 4.9).  Wave range data from the
  # irs_b[0-3]_WAVECUTvN.tbl file of the CUBISM package.

  irs:
    sl:
      '1':
        range: [7.51, 14.78]
        coefficients: [0.0, 8.2667]
      '2':
        range: [5.24, 7.6]
        coefficients: [0.0, 16.5333]
      '3':
        range: [7.34, 8.7]
        coefficients: [0.0, 8.2667]
    ll:
      '1':
        range: [20.5, 38.5]
        coefficients: [0.0, 2.9524]
      '2':
        range: [14.5, 21.1]
        coefficients: [0.0, 5.9048]
      '3':
        range: [19.4, 21.65]
        coefficients: [0.0, 2.9524]
      sh:
        range: [9.9661, 19.4386]
        coefficients: [600.]
      lh:
        range: [19.1095, 37.1661]
        coefficients: [600.]

Fit Outputs
===========

The fit results are saved in the same format as the feature table.
This enables the fit output file to be used as an input file to reload a model.
This is potentially useful for fitting multiple spectra that are similar (e.g., JWST spectral cubes).
