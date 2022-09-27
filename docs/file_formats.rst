############
File Formats
############

Model Packs
===========

YAML form
---------

The main way for the user to set up models, is by using (or copying and editing) one of the
default science packs provided as YAML files. These YAML files are internally converted to an
astropy Table (see below).

YAML format description here.

Table form
----------

When a model is set up from a YAML file, or when the fit results are saved to disk, the model is
described by a table compatible with astropy, called the "feature table".
This table has one row for each model component, each of which has a name, group, and kind column.
The name and group columns come from the YAML file, and can be any string, as defined by the user in the YAML file from which this table was generated.
The kind column indicates one of the supported component types, and determines the function by which the component wil be represented (e.g., BlackBody1D for dust_continuum, Drude1D for dust_feature, Gaussian1D for line)
Following these three columns, there is one column for each fit parameter (with nan values if the parameter is not relevant for a certain feature).
The columns contain 3-tuples with the value and the bounds: [value, min, max]
The bounds are masked if a parameter is fixed.
To reduce the number of columns, some parameters have different meaning depending on the kind of the feature.
For the dust_continuum components, the fit parameters are
temperature (`temperature`) and the "intrinsic" or "unattenuated" amplitude (`tau`).
For the dust_feature and (gas) line models, the fit parameters are
"intrinsic" or "unattenuated" amplitude (`power`), center wavelength (`wavelength`), and full width at half maximum (`fwhm`).
For the attenuation model, the fit is a measure for the optical depth (`tau`).

The astropy version of `Drude1D
<https://docs.astropy.org/en/stable/modeling/physical_models.html#drude1d>`_
and `Gaussian1D
<https://docs.astropy.org/en/stable/api/astropy.modeling.functional_models.Gaussian1D.html#astropy.modeling.functional_models.Gaussian1D>`_
are used. See :ref:`reference API <reference_API>` for a
description of BlackBody1D and the dust attenuation model.

The model tables are saved as ECSV tables as these provide readable (ASCII)
files with metadata support.
The metadata is used to read in a model from file, after a fit result has already been saved.

See :ref:`example fit output <example_fit_output>` for an example of
a feature table/output file.

.. note::
   IDL PAHFIT uses fractional FWHM as input/output for the Drude profiles, python
   PAHFIT uses FWHM as input/output for the Drude profiles.

Fit Outputs
-----------

Same as feature table.  This enables the fit output file to be used as an
input file to reload a model. This is potentially useful for fitting multiple spectra that are
similar (e.g., JWST spectral cubes).

Instrument packs
================

