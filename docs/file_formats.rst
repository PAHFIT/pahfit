############
File Formats
############

Model Packs
===========
YAML
=====
YAML is a data serialization language. 
The model packs are costructed as
lists of dictionaries of different feature components 
(eg., Line Features, Dust Features).
Each list consists of dictionaries which define the feature and the 
parameters associated with it.
For example ::
 
 PAH_5.3:
     kind : dust_feature
     wavelength: 5.27
     fwhm: 0.17918

This shows a single PAH feature and its corresponding
attributes (wavelength, fwhm).
The parameters specified depend on the kind of 
features and their respective functional forms.
Bounds (min and max) can also be set for different 
parameters as well as for 
all the parameters of a certain group.


The fit parameters for `Continua` (starlight and dust continuum) 
are temperature and tau (emmission strength).
For features like line and dust features, the 
fit parameters are wavelength, fwhm and power.
For `attenuation`, the fit parameters are tau, 
model and geometry and for `absorption features`
have tau, wavelength, fwhm and geometry.

See :ref:`Example Input <example_input>` for an 
example of YAML input model file.
For a detailed account of the YAML based science packs 
and the relevant syntax, please follow the
`Wiki <https://github.com/PAHFIT/pahfit/wiki/PAHFIT-2022#an-outline-of-the-syntax-
of-a-pahfits-yaml-science-packs>`_.

IPAC
======
The model packs are saved as IPAC tables as these provide readable (ASCII)
files that include unit information.
The model packs are given in a table with a row for each model component.
For each fit parameter, there are columns giving
the `Name` (e.g., BB3, DF10, H2 S(3), [OIV]),
the functional `Form` (e.g., BlackBody1D, Drude1D, Gaussian1D),
and multiple columns for the fit parameters.
For each fit parameter, the columns are for the actual value,
the bounds (i.e., min and max), and if the parameter is fixed or not.


For the `BlackBody1D` models, the fit parameters are temperature (`temp`)
and  "intrinsic" or "unattenuated" amplitude (`amp`).
For the `Drude1D` and `Gaussian1D` models, the fit parameters are
"intrinsic" or "unattenuated" amplitude (`amp`), center wavelength (`x_0`), and
full width at half maximum (`fwhm`).
For the dust attenuation model (`S07_att`), the fit parameters are
amplitude (`amp`).

The astropy version of `Drude1D
<https://docs.astropy.org/en/stable/modeling/physical_models.html#drude1d>`_
and `Gaussian1D
<https://docs.astropy.org/en/stable/api/astropy.modeling.functional_models.Gaussian1D.
html#astropy.modeling.functional_models.Gaussian1D>`_
are used. See :ref:`reference API <reference_API>` for a
description of BlackBody1D and the dust attenuation model.

See :ref:`example fit output <example_fit_output>` for an example of
a model pack/output file.

.. note::
   IDL PAHFIT uses fractional FWHM as input/output for the Drude profiles, python
   PAHFIT uses FWHM as input/output for the Drude profiles.

Fit Outputs
===========

Same as model packs.  This enables the fit output file to be used as a
model pack file, potentially useful for fitting multiple spectra that are
similar (e.g., JWST spectral cubes).

.. _example_input:


Example Input (YAML)
====================
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
